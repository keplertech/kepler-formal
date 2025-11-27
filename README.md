# kepler-formal

## Introduction

Kepler-Formal is a logic equivalence checking (LEC) tool that operates on the naja interchange format(https://github.com/najaeda/naja-if). It is designed for verifying incremental edits produced by the najaeda Python library(https://pypi.org/project/najaeda/) or any workflow that preserves stable indices across modifications so that corresponding items retain the same identifiers. Kepler-Formal focuses today on combinational equivalence checking only — sequential boundary changes are not supported yet and remain planned work.

### Acknowledgement

[<img src="https://nlnet.nl/logo/banner.png" width=100>](https://nlnet.nl/project/Naja)
[<img src="https://nlnet.nl/image/logos/NGI0Entrust_tag.svg" width=100>](https://nlnet.nl/project/Naja)

This project is supported and funded by NLNet through the [NGI0 Entrust](https://nlnet.nl/entrust) Fund.

## Dependencies

On Ubuntu:

```bash
sudo apt-get install g++ libboost-dev python3.9-dev capnproto libcapnp-dev libtbb-dev pkg-config bison flex doxygen
```

On macOS, using [Homebrew](https://brew.sh/):

```bash
brew install cmake doxygen capnp tbb bison flex boost
```

Ensure the versions of `bison` and `flex` installed via Homebrew take precedence over the macOS defaults by modifying your $PATH environment variable as follows:

```bash
export PATH="/opt/homebrew/opt/flex/bin:/opt/homebrew/opt/bison/bin:$PATH"
```
## Build

```bash
git clone --recurse-submodules https://github.com/keplertech/kepler-formal.git
cd kepler-formal
mkdir build
cd build
cmake ..
make
```

## Usage

```bash
"build/src/bin/kepler_formal <-naja_if/-verilog> <netlist1> <netlist2> [<liberty-file>...]"
```

## Contact

contact@keplertech.io
