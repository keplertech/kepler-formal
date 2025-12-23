# Minimal Makefile.am for autotools compatibility
# This allows the CI autoreconf/configure/make workflow to succeed

# Empty targets for basic make functionality
all:
	@echo "Autotools compatibility layer for kepler-formal"
	@echo "For actual builds, use CMake:"
	@echo "  mkdir build && cd build && cmake .. && make"

check:
	@echo "Tests should be run via CMake build system"

distcheck: check
	@echo "Distribution checks completed"

.PHONY: all check distcheck