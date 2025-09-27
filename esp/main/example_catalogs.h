#pragma once

#include "star_solver.h"

#define PATTERN_CATALOG_SIZE 4

static const pattern_entry_t EXT_RAM_BSS_ATTR pattern_catalog[] = {
    {{0, 0, 0, 0}},
    {{0, 0, 0, 0}},
    {{0, 0, 0, 0}},
    {{0, 0, 0, 0}},
};

#define STAR_CATALOG_SIZE 5

static const star_entry_t star_catalog[] = {
    {1.675310f, -0.919713f, -0.063225f, 0.602741f, -0.795428f, -0.720f},
    {3.837066f, -1.061691f, -0.374199f, -0.312298f, -0.873181f, -0.010f},
    {0.426383f, -0.998973f, 0.492715f, 0.223816f, -0.840915f, 0.460f},
    {3.681866f, -1.053711f, -0.423938f, -0.254278f, -0.869264f, 0.610f},
    {3.349799f, -1.041764f, -0.493798f, -0.104324f, -0.863296f, 1.250f}
};

