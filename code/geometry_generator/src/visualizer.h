#pragma once
#include <string>
#include <copengl.h>

using std::string;

class Visualizer
{
    bool is_running;
public:
    Visualizer(int width, int height, string window_title, bool full_screen, double camera_speed);
    COpenGL *opengl;
    char window_title[1000];
    void begin_draw();
    void end_draw();
    bool tick();
};
