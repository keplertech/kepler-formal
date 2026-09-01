# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

file(REMOVE_RECURSE "${DST}")
file(COPY "${SRC}/" DESTINATION "${DST}")
file(COPY "${NATIVE}" DESTINATION "${DST}")
