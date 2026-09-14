// Copyright 2024-2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#ifdef putc_unlocked
#error "The Windows stdio replacement must remain private to Glucose"
#endif

#include <cstdio>
#include "core/Solver.h"

int main() {
  FILE* output = std::tmpfile();
  if (!output) {
    std::perror("tmpfile");
    return 1;
  }

  Glucose::Solver solver;
  solver.certifiedOutput = output;
  for (unsigned int byte = 0; byte < 256; ++byte) {
    solver.write_char(static_cast<unsigned char>(byte));
  }
  if (std::fflush(output) != 0 || std::fseek(output, 0, SEEK_SET) != 0) {
    std::perror("rewind proof output");
    std::fclose(output);
    return 1;
  }
  for (int byte = 0; byte < 256; ++byte) {
    const int actual = std::fgetc(output);
    if (actual != byte) {
      std::fprintf(stderr, "Proof byte %d was %d\n", byte, actual);
      std::fclose(output);
      return 1;
    }
  }
  const bool exact_size = std::fgetc(output) == EOF && !std::ferror(output);
  const bool closed = std::fclose(output) == 0;
  solver.certifiedOutput = nullptr;
  return exact_size && closed ? 0 : 1;
}
