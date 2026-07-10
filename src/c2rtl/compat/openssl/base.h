// Copyright 2026 keplertech.io
// SPDX-License-Identifier: GPL-3.0-only

#ifndef KEPLER_C2RTL_COMPAT_OPENSSL_BASE_H_
#define KEPLER_C2RTL_COMPAT_OPENSSL_BASE_H_

// XLS includes BoringSSL's openssl/base.h. Homebrew OpenSSL provides the BN
// APIs XLS uses here, but not that BoringSSL umbrella header.
#include <openssl/crypto.h>
#include <openssl/opensslv.h>

#endif  // KEPLER_C2RTL_COMPAT_OPENSSL_BASE_H_
