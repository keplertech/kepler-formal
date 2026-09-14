/* Copyright 2024-2026 keplertech.io
 * SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#ifdef _WIN32
#include <io.h>
#include <process.h>
#include <stdio.h>
#include <sys/stat.h>

/* Only the CRT calls used by the embedded solver libraries. Function-like
 * macros avoid rewriting fields such as kissat_file.close. */
#define access(...) _access(__VA_ARGS__)
#define isatty(...) _isatty(__VA_ARGS__)
#define fileno(...) _fileno(__VA_ARGS__)
#define getpid(...) _getpid(__VA_ARGS__)
#define unlink(...) _unlink(__VA_ARGS__)
#define popen(...) _popen(__VA_ARGS__)
#define pclose(...) _pclose(__VA_ARGS__)
#ifndef R_OK
#define R_OK 4
#endif
#ifndef W_OK
#define W_OK 2
#endif
#ifndef S_ISDIR
#define S_ISDIR(mode) (((mode) & _S_IFMT) == _S_IFDIR)
#endif
#ifndef S_ISFIFO
#define S_ISFIFO(mode) (((mode) & _S_IFMT) == _S_IFIFO)
#endif
#else
#include_next <unistd.h>
#endif
