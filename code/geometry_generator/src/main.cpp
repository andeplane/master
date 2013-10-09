#include <marchingcubes.h>
#include <cvector.h>
#include <complexgeometry.h>

using std::cout;
using std::endl;
int main(int argc, char **argv)
{
    ComplexGeometry cg;
    cg.create_perlin_geometry(100, 100, 100, 1,1,1,3,0.2);

    MarchingCubes c;
    CVector system_length(0.1, 0.1, 0.1);
    c.create_marching_cubes_from_complex_geometry(cg,system_length,0.2);

    cout << "HORE" << endl;
    return 0;
}
