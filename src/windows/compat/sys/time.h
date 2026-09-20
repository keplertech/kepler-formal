/* Copyright 2024-2026 keplertech.io
 * SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#ifdef _WIN32
/* Winsock supplies Windows' timeval definition; it must precede windows.h.
 * No Winsock functions are called and no socket runtime is required. */
#include <winsock2.h>
#include <windows.h>

static inline int gettimeofday(struct timeval *value, void *timezone) {
  FILETIME now;
  ULARGE_INTEGER ticks;
  (void) timezone;
  GetSystemTimeAsFileTime(&now);
  ticks.LowPart = now.dwLowDateTime;
  ticks.HighPart = now.dwHighDateTime;
  /* FILETIME counts 100-ns intervals since 1601; timeval uses Unix seconds. */
  const unsigned long long unix_ticks = ticks.QuadPart - 116444736000000000ULL;
  value->tv_sec = (long) (unix_ticks / 10000000ULL);
  value->tv_usec = (long) ((unix_ticks % 10000000ULL) / 10ULL);
  return 0;
}
#else
#include_next <sys/time.h>
#endif
