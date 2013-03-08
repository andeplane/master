#include <grid.h>

#include <fstream>
#include <system.h>

#define min(a,b)                      (((a) < (b)) ? (a) : (b))
#define max(a,b)                      (((a) > (b)) ? (a) : (b))
#define clamp(value, lb, ub)          max( lb, min( ub, value ))

void Grid::read_matrix(string filename) {
    ifstream file (filename.c_str(), ios::in | ios::binary);
    file.read (reinterpret_cast<char*>(&Nx), sizeof(int));
    file.read (reinterpret_cast<char*>(&Ny), sizeof(int));
    file.read (reinterpret_cast<char*>(&Nz), sizeof(int));
    points = Nx*Ny*Nz;

    voxels = new unsigned char[points];
    normals   = new float[3*points];
    tangents1 = new float[3*points];
    tangents2 = new float[3*points];

    file.read (reinterpret_cast<char*>(voxels), points*sizeof(unsigned char));
    file.read (reinterpret_cast<char*>(normals), 3*points*sizeof(float));
    file.read (reinterpret_cast<char*>(tangents1), 3*points*sizeof(float));
    file.read (reinterpret_cast<char*>(tangents2), 3*points*sizeof(float));
    file.close();
}

Grid::Grid(string filename, System *system_)
{
    system = system_;
    read_matrix(filename);
}

unsigned char *Grid::get_voxel(const int &i, const int &j, const int &k) {
    if(i < 0 || i >= Nx || j < 0 || j >= Ny || k < 0 || k >= Nz) {
        return &voxels[((i+Nx)%Nx) + ((j+Ny)%Ny)*Nx+ ((k+Nz)%Nz)*Nx*Ny];
    }

    return &voxels[i + j*Nx + k*Nx*Ny];
}

unsigned char *Grid::get_voxel(const double &x, const double &y, const double &z) {
    int i =  x/system->Lx*Nx;
    int j =  y/system->Ly*Ny;
    int k =  z/system->Lz*Nz;

    return get_voxel(i,j,k);
}

unsigned char *Grid::get_voxel(const vec3 &r) {
    int i =  r(0)/system->Lx*Nx;
    int j =  r(1)/system->Ly*Ny;
    int k =  r(2)/system->Lz*Nz;

    return get_voxel(i,j,k);
}

int Grid::get_index_of_voxel(const vec3 &r) {
    int i =  r(0)/system->Lx*Nx;
    int j =  r(1)/system->Ly*Ny;
    int k =  r(2)/system->Lz*Nz;

    return i + j*Nx + k*Nx*Ny;
}

/*
void Grid::calculate_normals() {
    for(int i=0;i<cols;i++) {
        for(int j=0;j<rows;j++) {
            for(int k=0;k<slices;k++) {
                GridPoint *p = get_grid_point(i,j,k);

                for(int di=-1;di<=1;di++) {
                    for(int dj=-1;dj<=1;dj++) {
                        for(int dk=-1;dk<=1;dk++) {
                            p->normal.x -= get_grid_point(i+di,j+dj, k+dk)->is_wall*di;
                            p->normal.y -= get_grid_point(i+di,j+dj, k+dk)->is_wall*dj;
                            p->normal.z -= get_grid_point(i+di,j+dj, k+dk)->is_wall*dk;
                        }
                    }
                }
                p->normal = p->normal.Normalize();
            }
        }
    }
}

void Grid::calculate_inner_points() {
    for(int i=0;i<cols;i++) {
        for(int j=0;j<rows;j++) {
            for(int k=0;k<slices;k++) {
                GridPoint *p = get_grid_point(i,j,k);

                if(p->is_wall && p->normal.Length() > 0) {
                    p->is_wall_boundary = true;
                }
            }
        }
    }
}

void Grid::calculate_tangents() {
    for(int i=0;i<cols;i++) {
        for(int j=0;j<rows;j++) {
            for(int k=0;k<slices;k++) {
                GridPoint *p = get_grid_point(i,j,k);
                p->tangent1 = CVector(1,1,1);
                p->tangent1.RandomUniform();
                p->tangent1 = p->tangent1 - p->normal*p->normal.Dot(p->tangent1);
                p->tangent1.Normalize();
                p->tangent2 = p->normal.Cross(p->tangent1);
                p->tangent2.Normalize();
            }
        }
    }
}
*/
