// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "KeplerFormalUtils.h"

#include <cctype>
#include <fstream>
#include <sstream>

// LCOV_EXCL_START
std::string sanitizeFileToken(const std::string& input) {
  std::string out;
  out.reserve(input.size());
  for (unsigned char ch : input) {
    if (std::isalnum(ch) || ch == '_' || ch == '-' || ch == '.') {
      out.push_back(static_cast<char>(ch));
    } else {
      out.push_back('_');
      // LCOV_EXCL_STOP
    }
  }
  // LCOV_EXCL_START
  if (out.empty()) {
    out = "scope";
  }
  return out;
}
// LCOV_EXCL_STOP

namespace {

// LCOV_EXCL_START
std::string formatStringList(const std::vector<std::string>& values) {
  std::ostringstream oss;
  oss << "[";
  for (size_t i = 0; i < values.size(); ++i) {
    if (i) {
      oss << ", ";
    }
    oss << values[i];
  }
  oss << "]";
  return oss.str();
}
// LCOV_EXCL_STOP

}  // namespace

// LCOV_EXCL_START
void writeBoundaryTermsReport(
// LCOV_EXCL_STOP
    const std::filesystem::path& reportPath,
    const std::vector<KEPLER_FORMAL::SEC::ExtractedBoundaryReportEntry>& reports) {
  // LCOV_EXCL_START
  if (reports.empty()) {
    return;
    // LCOV_EXCL_STOP
  }

  // LCOV_EXCL_START
  std::ofstream report(reportPath, std::ios::trunc);
  report << "# SEC boundary terms report\n";
  report << "# Categories:\n";
  report << "# - top_input / top_output: original top-level interface terms.\n\n";
  for (size_t i = 0; i < reports.size(); ++i) {
    const auto& entry = reports[i];
    report << "- design: " << entry.design << "\n";
    report << "  signal: " << entry.signal << "\n";
    report << "  roles: " << formatStringList(entry.roles) << "\n";
    if (!entry.connectivitySkip.empty()) {
      report << "  connectivity_skip: " << entry.connectivitySkip << "\n";
    }
    if (i + 1 != reports.size()) {
      report << "\n";
    }
  }
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
void writeResetUnanchoredSkippedOutputsReport(
// LCOV_EXCL_STOP
    const std::filesystem::path& reportPath,
    const std::vector<std::string>& skippedOutputs) {
  // LCOV_EXCL_START
  if (skippedOutputs.empty()) {
    return;
    // LCOV_EXCL_STOP
  }

  // LCOV_EXCL_START
  std::ofstream report(reportPath, std::ios::trunc);
  report << "# SEC reset-unanchored skipped observed outputs\n";
  report << "# These top outputs were removed from the proof surface because\n";
  report << "# their cones depend on internal state without an inductive\n";
  report << "# cross-design anchor. SEC does not assume internal flop equality\n";
  report << "# by name; only top-level interface signals are name-aligned.\n\n";
  for (const auto& skippedOutput : skippedOutputs) {
    report << "- " << skippedOutput << "\n";
    // LCOV_EXCL_STOP
  }
// LCOV_EXCL_START
}
// LCOV_EXCL_STOP

// LCOV_EXCL_START
void writeMultiClockDomainSkippedOutputsReport(
// LCOV_EXCL_STOP
    const std::filesystem::path& reportPath,
    const std::vector<std::string>& skippedOutputs) {
  // LCOV_EXCL_START
  if (skippedOutputs.empty()) {
    return;
    // LCOV_EXCL_STOP
  }

  // LCOV_EXCL_START
  std::ofstream report(reportPath, std::ios::trunc);
  report << "# SEC multi-clock-domain skipped observed outputs\n";
  report << "# These top outputs were removed from the proof surface because\n";
  report << "# their cones span more than one extracted clock domain. CDC\n";
  report << "# modeling is intentionally outside this SEC pass, so the result\n";
  report << "# is reported as skipped coverage instead of assumed synchronous.\n\n";
  for (const auto& skippedOutput : skippedOutputs) {
    report << "- " << skippedOutput << "\n";
    // LCOV_EXCL_STOP
  }
// LCOV_EXCL_START
}
// LCOV_EXCL_STOP

void writeOpaqueCellSkippedOutputsReport(
    const std::filesystem::path& reportPath,
    const std::vector<std::string>& skippedOutputs) {
  if (skippedOutputs.empty()) {
    return;
  }

  std::ofstream report(reportPath, std::ios::trunc);
  report << "# SEC opaque-cell skipped top-level outputs\n";
  report << "# These outputs were not verified because backward cone traversal\n";
  report << "# reached an internal cell or pin without usable semantics. No free\n";
  report << "# or shared proof symbol was substituted for the opaque element.\n\n";
  for (const auto& skippedOutput : skippedOutputs) {
    report << "- " << skippedOutput << "\n";
  }
}
