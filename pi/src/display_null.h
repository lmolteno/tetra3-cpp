#pragma once

#include "display.h"

class NullDisplay : public IDisplay {
public:
    void clear() override {}
    void draw_text(int, int, const std::string &, bool) override {}
    void draw_pixel(int, int) override {}
    void update() override {}
};
