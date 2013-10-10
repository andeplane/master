#include <defines.h>

#ifdef OPENGL
#include <visualizer.h>
#include <copengl.h>
#endif
#include <cinifile.h>
#include <marchingcubes.h>
#include <cvector.h>
#include <complexgeometry.h>
#include <moviedata.h>

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
    cg.create_perlin_geometry(100, 100, 100, 1,1,1,3, threshold, true);
    // cg.create_sphere(100, 100, 100, 0.8, true, true, 1);
    cg.save_to_file("perlin.bin");
    CVector system_length = CVector(10*Lx, 10*Ly, 10*Lz);
    MarchingCubes c;
    c.create_marching_cubes_from_complex_geometry(cg, system_length, threshold, false);

    #ifdef OPENGL
    char *window_title = new char[1000];
    sprintf(window_title, "DSMC Geometry Visualizer (DSMCGV) - [%.2f fps]", 60.0);
    Visualizer v(screen_width, screen_height, string(window_title), false, 0.1);
    c.build_vbo();

    string state_folder = "/projects/master/code/base_code";
    MovieData movie_data(1,1000);
    movie_data.load_movie_files(state_folder);
    Timestep *timestep = movie_data.first_timestep;

    while(true) {
         v.render_begin();
         if(v.opengl->bool1) c.render_vbo();
         timestep->render();
         v.render_end();
         timestep = timestep->next;
        if(!v.is_running) break;
    }
    #endif

    return 0;
}
