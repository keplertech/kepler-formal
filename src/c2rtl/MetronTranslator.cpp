// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#include "MetronTranslator.h"

#include "metrolib/core/Log.h"
#include "metron/main/main.h"

#include <cerrno>
#include <cstring>
#include <mutex>
#include <stdexcept>
#include <string>

#include <fcntl.h>
#include <unistd.h>

Options global_options;

namespace KEPLER_FORMAL::C2RTL {
namespace {

class ScopedFileDescriptorRedirect {
 public:
  ScopedFileDescriptorRedirect(int targetFd,
                               const std::filesystem::path& destination)
      : targetFd_(targetFd), savedFd_(::dup(targetFd)) {
    if (savedFd_ < 0) {
      throw std::runtime_error("failed to duplicate file descriptor: " +
                               std::string(std::strerror(errno)));
    }

    const int logFd = ::open(destination.c_str(), O_CREAT | O_TRUNC | O_WRONLY,
                             0644);
    if (logFd < 0) {
      const auto error = std::string(std::strerror(errno));
      ::close(savedFd_);
      throw std::runtime_error("failed to open C2RTL log `" +
                               destination.string() + "`: " + error);
    }

    if (::dup2(logFd, targetFd_) < 0) {
      const auto error = std::string(std::strerror(errno));
      ::close(logFd);
      ::close(savedFd_);
      throw std::runtime_error("failed to redirect C2RTL log `" +
                               destination.string() + "`: " + error);
    }
    ::close(logFd);
  }

  ScopedFileDescriptorRedirect(const ScopedFileDescriptorRedirect&) = delete;
  ScopedFileDescriptorRedirect& operator=(const ScopedFileDescriptorRedirect&) =
      delete;

  ~ScopedFileDescriptorRedirect() {
    fflush(nullptr);
    if (savedFd_ >= 0) {
      ::dup2(savedFd_, targetFd_);
      ::close(savedFd_);
    }
  }

 private:
  int targetFd_;
  int savedFd_;
};

std::mutex& metronMutex() {
  static std::mutex mutex;
  return mutex;
}

std::vector<std::string> toStringVector(
    const std::vector<std::filesystem::path>& paths) {
  std::vector<std::string> result;
  result.reserve(paths.size());
  for (const auto& path : paths) {
    result.push_back(path.string());
  }
  return result;
}

void addBundledMetronIncludePath(std::vector<std::string>& paths) {
#if defined(KEPLER_METRON_ROOT)
  paths.push_back(KEPLER_METRON_ROOT);
#endif
}

void copyBundledSystemVerilogSupport(const std::filesystem::path& outputPath) {
#if defined(KEPLER_METRON_ROOT)
  const auto source =
      std::filesystem::path(KEPLER_METRON_ROOT) / "metron" / "metron_tools.sv";
  const auto destinationDir = outputPath.parent_path() / "metron";
  const auto destination = destinationDir / "metron_tools.sv";
  std::error_code ec;
  std::filesystem::create_directories(destinationDir, ec);
  ec.clear();
  std::filesystem::copy_file(source, destination,
                             std::filesystem::copy_options::overwrite_existing,
                             ec);
#else
  (void)outputPath;
#endif
}

}  // namespace

void translateWithMetron(const MetronTranslationOptions& options) {
  std::lock_guard lock(metronMutex());

  if (options.inputPath.empty()) {
    throw std::runtime_error("C2RTL translation requires an input path");
  }
  if (options.outputPath.empty()) {
    throw std::runtime_error("C2RTL translation requires an output path");
  }

  std::filesystem::create_directories(options.outputPath.parent_path());
  std::filesystem::create_directories(options.stdoutLog.parent_path());
  std::filesystem::create_directories(options.stderrLog.parent_path());
  copyBundledSystemVerilogSupport(options.outputPath);

  ScopedFileDescriptorRedirect stdoutRedirect(STDOUT_FILENO,
                                              options.stdoutLog);
  ScopedFileDescriptorRedirect stderrRedirect(STDERR_FILENO,
                                              options.stderrLog);

  global_options = Options();
  global_options.quiet = true;
  global_options.monochrome = true;
  global_options.src_name = options.inputPath.string();
  global_options.dst_name = options.outputPath.string();
  global_options.inc_paths = toStringVector(options.includePaths);
  addBundledMetronIncludePath(global_options.inc_paths);

  TinyLog::get().reset();
  TinyLog::get().mute();
  TinyLog::get().mono();

  const int rc = main_new();
  fflush(nullptr);
  if (rc != 0) {
    throw std::runtime_error("Metron C2RTL translation failed for `" +
                             options.inputPath.string() + "`");
  }
}

}  // namespace KEPLER_FORMAL::C2RTL
