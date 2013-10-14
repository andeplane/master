#include <defines.h>

#ifdef OPENGL
#include <visualizer.h>
#endif
#include <cinifile.h>
#include <marchingcubes.h>
#include <cvector.h>
#include <complexgeometry.h>
#include <testshader.h>
#include <camera.h>

using std::cout;
using std::endl;

int main(int argc, char **argv)
{
    cout << "Controls: Use WSAD and the mouse to move around!" << endl;
    CIniFile ini;
    ini.load("settings.ini");
    int screen_width = ini.getint("screen_width");
    int screen_height = ini.getint("screen_height");
    double Lx = ini.getdouble("Lx");
    double Ly = ini.getdouble("Ly");
    double Lz = ini.getdouble("Lz");
    double threshold = ini.getdouble("threshold");

    ComplexGeometry cg;
    cg.create_perlin_geometry(100, 100, 100, 4,3,2,3, threshold, true, 1);
    // cg.create_empty_space(100, 100, 100, true, 1.0);
    cg.save_to_file_2("./perlin4/",CVector(1,1,1));
    return 0;
    // cg.save_to_file("perlin.bin");

//    string text_files_base_filename = ini.getstring("text_files_base_filename");

//    cg.load_text_files(text_files_base_filename,CVector(100, 100, 50), threshold);

    CVector system_length = CVector(Lx, Ly, Lz);
    MarchingCubes c;
    c.create_marching_cubes_from_complex_geometry(cg, system_length, threshold, false);
    // c.create_marching_cubes_from_complex_geometry(cg, CVector(Lx, Ly, Lz), threshold);

    #ifdef OPENGL
    char *window_title = new char[1000];
    sprintf(window_title, "DSMC Geometry Visualizer (DSMCGV) - [%.2f fps]", 60.0);
    Visualizer v(screen_width, screen_height, string(window_title), false, 0.1);
//    for(int i=0; i<marching_cubes.size(); i++) {
//        marching_cubes[i].build_vbo();
//    }
    c.build_vbo(v.opengl);
    TestShader shader;
    shader.COpenGLPointer = v.opengl;
    try {
        shader.Initialize("Balle");
    } catch (string ex) {
        cout << ex << endl;
    }
    
    float t = 0;
    float dt = 0.01;
    float pi = 3.141592653;
    float frequency = 0.01;
    float omega = 2*pi*frequency;
    float A = 10.0;
    while(true) {
        try {
            CVector close_to_me = CVector(0,0,0);//v.opengl->camera->position;
            close_to_me.x += A*cos(omega*t);
            close_to_me.y += A*sin(omega*t);
            // close_to_me.z += A*cos(omega*t)*A*sin(omega*t);
            t += dt;
            shader.lightpos = close_to_me;
            shader.targetdir = v.opengl->camera->target;
            v.render_begin();
            shader.Start();
            c.render_vbo();
            shader.End();
            v.render_end();
            if(!v.is_running) break;
        } catch (string ex) {
            cout << "Error: " << ex << endl;
        }
        
    }
    #endif

    return 0;
}
