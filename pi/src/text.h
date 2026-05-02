#pragma once

#include "canvas.h"
#include <string>

namespace text {

// Glyphs are 6x8. `large` doubles each axis to 12x16. The degree symbol (0xB0)
// is mapped to glyph index 95.
void draw(ICanvas &canvas, int x, int y, const std::string &s, bool large = false);

}
