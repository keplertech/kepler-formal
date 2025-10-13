# kepler-formal

## Introduction

Kepler-Formal is a logic equivalence checking (LEC) tool that operates on the naja interchange format. It is designed for verifying incremental edits produced by the najaeda Python library or any workflow that preserves stable indices across modifications so that corresponding items retain the same identifiers. Kepler-Formal focuses today on combinational equivalence checking only — sequential boundary changes are not supported yet and remain planned work.

### Acknowledgement

[<img src="https://nlnet.nl/logo/banner.png" width=100>](https://nlnet.nl/project/Naja)
[<img src="https://nlnet.nl/image/logos/NGI0Entrust_tag.svg" width=100>](https://nlnet.nl/project/Naja)

This project is supported and funded by NLNet through the [NGI0 Entrust](https://nlnet.nl/entrust) Fund.

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
"build/src/bin/kepler_formal <naja-if-dir-1> <naja-if-dir-2> [<liberty-file>...]"
```

## Contact

contact@keplertech.io
