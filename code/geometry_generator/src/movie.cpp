#include <defines.h>

#ifdef OPENGL
#include <visualizer.h>
#endif
#include <cinifile.h>
#include <marchingcubes.h>
#include <cvector.h>
#include <complexgeometry.h>

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
    cg.create_perlin_geometry(200, 200, 200, 1,1,1,3,threshold, false);
    CVector system_length = CVector(Lx, Ly, Lz);
    MarchingCubes c;
    c.create_marching_cubes_from_complex_geometry(cg, system_length, threshold);

    #ifdef OPENGL
    char *window_title = new char[1000];
    sprintf(window_title, "DSMC Geometry Visualizer (DSMCGV) - [%.2f fps]", 60.0);
    Visualizer v(screen_width, screen_height, string(window_title), false, 0.1);
    c.build_vbo();

    while(true) {
         v.render_begin();
         c.render_vbo();
         v.render_end();
        if(!v.is_running) break;
    }
    #endif

    return 0;
}
