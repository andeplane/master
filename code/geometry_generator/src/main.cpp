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
    cg.create_perlin_geometry(100, 100, 100, 1,1,1,3, threshold, true);
    cg.save_to_file("perlin.bin");

    string text_files_base_filename = ini.getstring("text_files_base_filename");

    cg.load_text_files(text_files_base_filename,CVector(100, 100, 50), threshold);

    CVector system_length = CVector(Lx, Ly, Lz);
    MarchingCubes c;
    c.create_marching_cubes_from_complex_geometry(cg, system_length, threshold, true);
    // c.create_marching_cubes_from_complex_geometry(cg, CVector(Lx, Ly, Lz), threshold);

    #ifdef OPENGL
    char *window_title = new char[1000];
    sprintf(window_title, "DSMC Geometry Visualizer (DSMCGV) - [%.2f fps]", 60.0);
    Visualizer v(screen_width, screen_height, string(window_title), false, 0.1);
//    for(int i=0; i<marching_cubes.size(); i++) {
//        marching_cubes[i].build_vbo();
//    }
    c.build_vbo();

    while(true) {
         v.render_begin();
         c.render_vbo();
//         for(int i=0; i<marching_cubes.size(); i++) {
//             marching_cubes[i].render_vbo();
//         }
//         glBegin(GL_POINTS);
//         float scale = 50.0;
//         float scale_z = 10.0;
//         CVector offset(-25,-25,-5);
//         for(int i=0; i<cg.nx; i++) {
//             for(int j=0; j<cg.ny; j++) {
//                 for(int k=0; k<cg.nz; k++) {
//                     int index = i + j*cg.nx + k*cg.nx*cg.ny;
//                     if(cg.vertices[index] == 1) {
//                         glVertex3f( (scale*i)/cg.nx+offset.x, (scale*j)/cg.ny + offset.y, (scale_z*k)/cg.nz + offset.z);
//                     }
//                 }
//             }

//         }
//         glEnd();
         v.render_end();
        if(!v.is_running) break;
    }
    #endif

    return 0;
}
