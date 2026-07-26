// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

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

void writeBoundaryTermsReport(
    const std::filesystem::path& reportPath,
    const std::vector<KEPLER_FORMAL::SEC::ExtractedBoundaryReportEntry>& reports);

void writeResetUnanchoredSkippedOutputsReport(
    const std::filesystem::path& reportPath,
    const std::vector<std::string>& skippedOutputs);

void writeMultiClockDomainSkippedOutputsReport(
    const std::filesystem::path& reportPath,
    const std::vector<std::string>& skippedOutputs);
