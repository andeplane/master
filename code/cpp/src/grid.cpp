#include "grid.h"

#define min(a,b)                      (((a) < (b)) ? (a) : (b))
#define max(a,b)                      (((a) > (b)) ? (a) : (b))
#define clamp(value, lb, ub)          max( lb, min( ub, value ))

Grid::Grid(mat &M,System *system_)
{
    system = system_;
    cols = M.n_cols;
    rows = M.n_rows;
    slices = M.n_cols;

    points.reserve(M.n_elem*slices);

    for(int k=0;k<slices;k++) {
        for(int i=0;i<cols;i++) {
            for(int j=0;j<rows;j++) {
                GridPoint p(M(j,i));

                p.T = system->wall_temperature;
                p.i = i;
                p.j = j;
                p.k = k;
                points.push_back(p);
            }
        }
    }
    calculate_normals();
    calculate_tangents();
    calculate_inner_points();
}

inline GridPoint* Grid::get_grid_point(const int &i, const int &j, const int &k) {
    if(i < 0 || i >= cols || j < 0 || j >= rows || k < 0 || k >= slices) {
        return &points[((j+rows)%rows) + ((i+cols)%cols)*rows + ((k+slices)%slices)*rows*cols];
    }

    return &points[j + i*rows + k*rows*cols];
}

GridPoint* Grid::get_grid_point(const double &x, const double &y, const double &z) {
    int i =  clamp((int)(x/system->Lx*cols),0,cols-1);
    int j =  clamp((int)(y/system->Ly*rows),0,rows-1);
    int k =  clamp((int)(z/system->Lz*slices),0,slices-1);

    return get_grid_point(i,j,k);
}

GridPoint* Grid::get_grid_point(double *r) {
    int i =  clamp((int)(r[0]/system->Lx*cols),0,cols-1);
    int j =  clamp((int)(r[1]/system->Ly*rows),0,rows-1);
    int k =  clamp((int)(r[2]/system->Lz*slices),0,slices-1);

    return get_grid_point(i,j,k);
}

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
                    p->normal = p->normal.Normalize();
                }
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
