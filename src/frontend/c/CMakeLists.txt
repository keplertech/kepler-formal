# Copyright 2024-2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

add_library(kepler_c_frontend
    CFrontend.cpp
)

target_include_directories(kepler_c_frontend PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(kepler_c_frontend
    PRIVATE
        kepler_c2rtl
)
