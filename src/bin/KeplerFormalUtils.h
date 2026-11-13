// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "strategy/SequentialEquivalenceStrategy.h"

// Completed SEC verdicts use stable, distinct process exit codes.
inline constexpr int kSecProvedExitCode = 0;
inline constexpr int kSecPartiallyProvedExitCode = 1;
inline constexpr int kSecInconclusiveExitCode = 2;
inline constexpr int kSecCounterexampleExitCode = 3;

// Shared helper for consistent filename handling.
std::string sanitizeFileToken(const std::string& input);

// Creates an empty file `<prefix><unique token><extension>` in the system
// temporary directory and returns its path. The file is created exclusively:
// no two calls get the same file, whether they come from this process or from
// another run that shares the directory.
std::filesystem::path createUniqueTemporaryFile(const std::string& prefix,
                                                const std::string& extension);

void writeBoundaryTermsReport(
    const std::filesystem::path& reportPath,
    const std::vector<KEPLER_FORMAL::SEC::ExtractedBoundaryReportEntry>& reports);

void writeResetUnanchoredSkippedOutputsReport(
    const std::filesystem::path& reportPath,
    const std::vector<std::string>& skippedOutputs);

void writeMultiClockDomainSkippedOutputsReport(
    const std::filesystem::path& reportPath,
    const std::vector<std::string>& skippedOutputs);

void writeOpaqueCellSkippedOutputsReport(
    const std::filesystem::path& reportPath,
    const std::vector<std::string>& skippedOutputs);
