// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#include "DesignBoundary.h"

#include <algorithm>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <utility>

#include "NLName.h"
#include "SNLBitNet.h"
#include "SNLBitTerm.h"
#include "SNLBundleTerm.h"
#include "SNLBusTerm.h"
#include "SNLBusTermBit.h"
#include "SNLDesign.h"
#include "SNLInstance.h"
#include "SNLInstTerm.h"
#include "SNLNetComponent.h"

namespace KEPLER_FORMAL {
namespace {

using naja::NL::NLID;
using naja::NL::NLName;
using naja::NL::SNLBitNet;
using naja::NL::SNLBitTerm;
using naja::NL::SNLBusTermBit;
using naja::NL::SNLDesign;
using naja::NL::SNLInstance;
using naja::NL::SNLInstTerm;
using naja::NL::SNLNetComponent;
using naja::NL::SNLTerm;

struct PinSpec {
  BoundaryPort port;
  SNLBitTerm* term = nullptr;
};

struct BoundarySpec {
  size_t pairIndex = 0;
  std::string path;
  std::vector<std::string> components;
  std::vector<PinSpec> pins;
};

std::string quote(const std::string& text) {
  return "`" + text + "`";
}

std::vector<std::string> splitPath(const std::string& path) {
  if (path.empty()) {
    throw std::invalid_argument("boundary instance path must not be empty");
  }
  if (path.front() == '/' || path.back() == '/') {
    throw std::invalid_argument(
        "boundary instance path " + quote(path) +
        " must be top-relative without a leading or trailing slash");
  }

  std::vector<std::string> result;
  size_t begin = 0;
  while (begin < path.size()) {
    const size_t end = path.find('/', begin);
    const size_t length =
        end == std::string::npos ? path.size() - begin : end - begin;
    if (length == 0) {
      throw std::invalid_argument(
          "boundary instance path " + quote(path) +
          " contains an empty component");
    }
    result.emplace_back(path.substr(begin, length));
    if (end == std::string::npos) {
      break;
    }
    begin = end + 1;
  }
  return result;
}

bool isPathPrefix(const std::vector<std::string>& prefix,
                  const std::vector<std::string>& path) {
  return prefix.size() <= path.size() &&
         std::equal(prefix.begin(), prefix.end(), path.begin());
}

std::string encodeName(std::string_view name) {
  std::ostringstream encoded;
  encoded << std::hex << std::setfill('0');
  for (const unsigned char character : name) {
    encoded << std::setw(2) << static_cast<unsigned>(character);
  }
  return encoded.str();
}

std::string bitToken(int32_t bit) {
  if (bit < 0) {
    return "n" + std::to_string(-static_cast<int64_t>(bit));
  }
  return "p" + std::to_string(bit);
}

std::string makeTopTermName(size_t pairIndex,
                            const std::string& pinName,
                            int32_t bit) {
  return "__kepler_boundary_p" + std::to_string(pairIndex) + "_n" +
         encodeName(pinName) + "_b" + bitToken(bit);
}

SNLInstance* resolveInstance(SNLDesign* top,
                             const std::vector<std::string>& components,
                             const std::string& path) {
  SNLDesign* current = top;
  SNLInstance* instance = nullptr;
  std::string resolved;
  for (const auto& component : components) {
    instance = current->getInstance(NLName(component));
    if (instance == nullptr) {
      throw std::invalid_argument(
          "boundary instance path " + quote(path) +
          " does not resolve at " +
          quote(resolved.empty() ? component : resolved + "/" + component));
    }
    if (!resolved.empty()) {
      resolved += '/';
    }
    resolved += component;
    current = instance->getModel();
  }
  return instance;
}

bool isDriver(const SNLNetComponent* component) {
  const auto direction = component->getDirection();
  if (dynamic_cast<const SNLInstTerm*>(component) != nullptr) {
    return direction == SNLTerm::Direction::Output ||
           direction == SNLTerm::Direction::InOut;
  }
  if (dynamic_cast<const SNLBitTerm*>(component) != nullptr) {
    return direction == SNLTerm::Direction::Input ||
           direction == SNLTerm::Direction::InOut;
  }
  return false;
}

void validateInputConnectivity(const std::string& path,
                               const SNLInstTerm* input) {
  SNLBitNet* net = input->getNet();
  if (net == nullptr) {
    throw std::invalid_argument(
        "boundary instance " + quote(path) + " pin " +
        quote(input->getBitTerm()->getName().getString()) +
        " is unconnected");
  }
  if (net->isConstant()) {
    return;
  }

  size_t driverCount = 0;
  for (const auto* component : net->getComponents()) {
    driverCount += isDriver(component) ? 1 : 0;
  }
  if (driverCount != 1) {
    throw std::invalid_argument(
        "boundary input pin " +
        quote(input->getBitTerm()->getName().getString()) +
        " on instance " + quote(path) +
        " must have exactly one driver (found " +
        std::to_string(driverCount) + ")");
  }
}

void validateOutputConnectivity(const std::string& path,
                                const SNLInstTerm* output) {
  SNLBitNet* net = output->getNet();
  if (net == nullptr) {
    // An unused output still supplies a paired environment variable.
    return;
  }
  if (net->isConstant()) {
    throw std::invalid_argument(
        "boundary output pin " + quote(output->getBitTerm()->getName().getString()) +
        " on instance " + quote(path) + " is connected to a constant net");
  }

  size_t driverCount = 0;
  bool targetIsDriver = false;
  for (const auto* component : net->getComponents()) {
    if (!isDriver(component)) {
      continue;
    }
    ++driverCount;
    targetIsDriver = targetIsDriver || component == output;
  }
  if (!targetIsDriver || driverCount != 1) {
    throw std::invalid_argument(
        "boundary output pin " + quote(output->getBitTerm()->getName().getString()) +
        " on instance " + quote(path) +
        " must be the only driver of its net");
  }
}

BoundaryPort describePort(size_t pairIndex, const SNLBitTerm* bitTerm) {
  if (bitTerm->isUnnamed()) {
    throw std::invalid_argument(
        "boundary instances with unnamed pins are not supported");
  }

  BoundaryPort port;
  port.pairIndex = pairIndex;
  port.pinName = bitTerm->getName().getString();
  port.bit = bitTerm->getBit();
  port.isInput = bitTerm->getDirection() == SNLTerm::Direction::Input;

  if (const auto* busBit = dynamic_cast<const SNLBusTermBit*>(bitTerm)) {
    const auto* bus = busBit->getBus();
    if (bus->getBundleOwner() != nullptr) {
      throw std::invalid_argument(
          "bundled boundary pins are not supported: " + quote(port.pinName));
    }
    port.width = static_cast<size_t>(bus->getWidth());
    port.msb = bus->getMSB();
    port.lsb = bus->getLSB();
  } else {
    if (bitTerm->getBundleOwner() != nullptr) {
      throw std::invalid_argument(
          "bundled boundary pins are not supported: " + quote(port.pinName));
    }
    port.width = 1;
    port.msb = 0;
    port.lsb = 0;
  }

  port.topTermName = makeTopTermName(pairIndex, port.pinName, port.bit);
  return port;
}

std::vector<BoundarySpec> preflight(SNLDesign* top,
                                    const BoundaryPairs& pairs,
                                    size_t side) {
  if (top == nullptr) {
    throw std::invalid_argument("boundary top design must not be null");
  }
  if (side > 1) {
    throw std::invalid_argument("boundary side must be 0 or 1");
  }
  if (top->isPrimitive()) {
    throw std::invalid_argument("a primitive design cannot be a boundary top");
  }

  std::vector<BoundarySpec> specs;
  specs.reserve(pairs.size());
  for (size_t pairIndex = 0; pairIndex < pairs.size(); ++pairIndex) {
    if (pairs[pairIndex].first.empty() || pairs[pairIndex].second.empty()) {
      throw std::invalid_argument(
          "boundary pair " + std::to_string(pairIndex) +
          " must name an instance path on both sides");
    }
    BoundarySpec spec;
    spec.pairIndex = pairIndex;
    spec.path = side == 0 ? pairs[pairIndex].first : pairs[pairIndex].second;
    spec.components = splitPath(spec.path);
    specs.push_back(std::move(spec));
  }

  for (size_t i = 0; i < specs.size(); ++i) {
    for (size_t j = i + 1; j < specs.size(); ++j) {
      if (specs[i].components == specs[j].components) {
        throw std::invalid_argument(
            "duplicate boundary instance path " + quote(specs[i].path));
      }
      if (isPathPrefix(specs[i].components, specs[j].components) ||
          isPathPrefix(specs[j].components, specs[i].components)) {
        throw std::invalid_argument(
            "nested boundary instance paths are not supported: " +
            quote(specs[i].path) + " and " + quote(specs[j].path));
      }
    }
  }

  std::set<std::string> generatedTopNames;
  for (auto& spec : specs) {
    SNLInstance* target = resolveInstance(top, spec.components, spec.path);
    if (!target->getModel()->getInstances().empty()) {
      throw std::invalid_argument(
          "boundary instance " + quote(spec.path) +
          " is not a leaf: its model contains child instances");
    }
    for (SNLInstTerm* instTerm : target->getInstTerms()) {
      SNLBitTerm* bitTerm = instTerm->getBitTerm();
      const auto direction = bitTerm->getDirection();
      if (direction == SNLTerm::Direction::InOut ||
          direction == SNLTerm::Direction::Undefined) {
        throw std::invalid_argument(
            "boundary instance " + quote(spec.path) + " pin " +
            quote(bitTerm->getName().getString()) +
            " has unsupported direction " + direction.getString());
      }
      if (direction == SNLTerm::Direction::Input) {
        validateInputConnectivity(spec.path, instTerm);
      }
      if (direction == SNLTerm::Direction::Output) {
        validateOutputConnectivity(spec.path, instTerm);
      }

      PinSpec pin;
      pin.port = describePort(spec.pairIndex, bitTerm);
      pin.term = bitTerm;
      if (!generatedTopNames.insert(pin.port.topTermName).second) {
        throw std::invalid_argument(
            "boundary pins generate duplicate top term " +
            quote(pin.port.topTermName));
      }
      // SEC reports and aligns observed outputs by their display names.
      if (top->getTerm(NLName(pin.port.topTermName)) != nullptr) {
        throw std::invalid_argument(
            "boundary top term collides with existing term " +
            quote(pin.port.topTermName));
      }

      spec.pins.push_back(std::move(pin));
    }

    if (spec.pins.empty()) {
      throw std::invalid_argument(
          "boundary instance " + quote(spec.path) + " has no pins");
    }
    std::sort(
        spec.pins.begin(), spec.pins.end(), [](const PinSpec& lhs, const PinSpec& rhs) {
          return std::tie(lhs.port.pinName, lhs.port.bit, lhs.port.isInput) <
                 std::tie(rhs.port.pinName, rhs.port.bit, rhs.port.isInput);
        });
  }
  return specs;
}

using PortKey = std::tuple<size_t, std::string, int32_t>;

std::map<PortKey, const BoundaryPort*> indexPorts(
    const std::vector<BoundaryPort>& ports,
    const char* sideName) {
  std::map<PortKey, const BoundaryPort*> indexed;
  for (const auto& port : ports) {
    const PortKey key{port.pairIndex, port.pinName, port.bit};
    if (!indexed.emplace(key, &port).second) {
      throw std::invalid_argument(
          std::string("duplicate ") + sideName + " boundary signature for pair " +
          std::to_string(port.pairIndex) + " pin " + quote(port.pinName) +
          " bit " + std::to_string(port.bit));
    }
  }
  return indexed;
}

std::string describePortKey(const PortKey& key) {
  return "pair " + std::to_string(std::get<0>(key)) + " pin " +
         quote(std::get<1>(key)) + " bit " + std::to_string(std::get<2>(key));
}

}  // namespace

BoundarySelection::BoundarySelection(SNLDesign* top,
                                     const BoundaryPairs& pairs, size_t side) {
  for (const auto& spec : preflight(top, pairs, side)) {
    for (const auto& pin : spec.pins) {
      ports_.push_back(pin.port);
    }
  }
}

const BoundaryPort* LeafBoundary::getPort(DNLID id) const {
  const auto found = portIndices_.find(id);
  return found == portIndices_.end() ? nullptr : &ports_[found->second];
}

LeafBoundary::LeafBoundary(const naja::DNL::DNLFull& dnl,
                           const BoundaryPairs& pairs, size_t side) {
  using namespace naja::DNL;
  const auto specs = preflight(
      const_cast<SNLDesign*>(dnl.getTop().getSNLModel()), pairs, side);
  for (const auto& spec : specs) {
    const DNLInstanceFull* instance = &dnl.getTop();
    for (const auto& component : spec.components) {
      instance = &instance->getChildInstance(
          instance->getSNLModel()->getInstance(NLName(component)));
    }
    instances_.insert(instance->getID());
    for (const auto& pin : spec.pins) {
      const auto& term = instance->getTerminalFromBitTerm(pin.term);
      const DNLID id = term.getID();
      // Leaf outputs must remain actual drivers in the original flattened DNL.
      // Internal wire aliases or constant nets cannot be overridden by PI flags.
      // An unconnected output needs no iso: it is an unused free input.
      if (!pin.port.isInput && term.getIsoID() != DNLID_MAX) {
        const auto& iso = dnl.getDNLIsoDB().getIsoFromIsoIDconst(term.getIsoID());
        if (iso.isConstant() || iso.getDrivers().size() != 1 ||
            iso.getDrivers().front() != id) {
          throw std::invalid_argument(
              "boundary leaf output " + quote(spec.path + "/" + pin.port.pinName) +
              " must be the sole nonconstant driver of its DNL iso");
        }
      }
      portIndices_.emplace(id, ports_.size());
      ports_.push_back(pin.port);
      (pin.port.isInput ? outputs_ : inputs_).push_back(id);
    }
  }
}

void validateBoundaryInterfaces(const std::vector<BoundaryPort>& left,
                                const std::vector<BoundaryPort>& right) {
  const auto leftPorts = indexPorts(left, "left");
  const auto rightPorts = indexPorts(right, "right");
  if (leftPorts.size() != rightPorts.size()) {
    throw std::invalid_argument(
        "boundary interfaces have different pin-bit counts (left=" +
        std::to_string(leftPorts.size()) + ", right=" +
        std::to_string(rightPorts.size()) + ")");
  }

  for (const auto& [key, leftPort] : leftPorts) {
    const auto rightIt = rightPorts.find(key);
    if (rightIt == rightPorts.end()) {
      throw std::invalid_argument(
          "right boundary interface is missing " + describePortKey(key));
    }
    const BoundaryPort& rightPort = *rightIt->second;
    if (leftPort->isInput != rightPort.isInput) {
      throw std::invalid_argument(
          "boundary direction mismatch for " + describePortKey(key));
    }
    if (leftPort->width != rightPort.width ||
        leftPort->msb != rightPort.msb || leftPort->lsb != rightPort.lsb) {
      throw std::invalid_argument(
          "boundary bus-shape mismatch for " + describePortKey(key));
    }
    if (leftPort->topTermName != rightPort.topTermName) {
      throw std::invalid_argument(
          "boundary synthetic-name mismatch for " + describePortKey(key));
    }
  }
}

}  // namespace KEPLER_FORMAL
