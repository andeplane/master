#pragma once
#include <GL/glfw.h>      // Include OpenGL Framework library
#include <vector>
#include <string>
#include <cbitmap.h>

using std::vector;
using std::string;

class COpenGL;
// internal texture structure
class COpenGLTexture {
public:
    string      category;
    GLuint		id;
    GLsizei		imgHeight, imgWidth;
    GLsizei		txDimensions;
    GLfloat		alpha;
    bool        has_alpha;
};

class CTexture
{
public:
    COpenGL *opengl;
    CBitMap bmp;
    GLuint texture_id;
    vector<COpenGLTexture> textures;
    vector<string> names;

    CTexture(COpenGL *ogl);
    void create_sphere1(string name, int w);
    void load_texture(CBitMap* bmp, COpenGLTexture* texture, bool has_alpha);
    void render_billboards(vector<float> &positions);
};
