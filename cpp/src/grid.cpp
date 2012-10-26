#include "grid.h"

#define min(a,b)                      (((a) < (b)) ? (a) : (b))
#define max(a,b)                      (((a) > (b)) ? (a) : (b))
#define clamp(value, lb, ub)          max( lb, min( ub, value ))

Grid::Grid(mat M,System *_system)
{
    system = _system;
    points.reserve(M.n_elem);
    cols = M.n_cols;
    rows = M.n_rows;

    for(int i=0;i<cols;i++)
        for(int j=0;j<rows;j++) {
            GridPoint p(M(j,i));
            // p.T = system->T + (system->wall_temperature - system->T)*((double)i/(M.n_cols-1));
            p.T = system->wall_temperature;
            p.i = i;
            p.j = j;
            points.push_back(p);
        }

    calculate_normals();
    calculate_tangents();
    calculate_inner_points();
}

GridPoint* Grid::get_grid_point(const int i, const int j) {
    if(i < 0 || i >= cols || j < 0 || j >= rows) {
        return &points[((j+rows)%rows) + ((i+cols)%cols)*rows];
    }

    return &points[j + i*rows];
}

GridPoint* Grid::get_grid_point(const double x, const double y) {
    int i =  clamp((int)(x/system->width*cols),0,cols-1);
    int j =  clamp((int)(y/system->height*rows),0,rows-1);

    return get_grid_point(i,j);
}

GridPoint* Grid::get_grid_point(const vec r, const int idx) {
    int i =  clamp((int)(r(0)/system->width*cols),0,cols-1);
    int j =  clamp((int)(r(1)/system->height*rows),0,rows-1);

    return get_grid_point(i,j);
}

// #define DEBUG

void Grid::calculate_normals() {
    for(int i=0;i<cols;i++)
        for(int j=0;j<rows;j++) {
            GridPoint *p = get_grid_point(i,j);

            for(int k=-1;k<=1;k++)
                for(int l=-1;l<=1;l++) {
                    p->normal.x -= get_grid_point(i+k,j+l)->is_wall*k;
                    p->normal.y -= get_grid_point(i+k,j+l)->is_wall*l;
                }
            p->normal = p->normal.Normalize();
        }
}

void Grid::calculate_inner_points() {
    for(int i=0;i<cols;i++)
        for(int j=0;j<rows;j++) {
            GridPoint *p = get_grid_point(i,j);

            if(p->is_wall && p->normal.Length() > 0)
                p->is_wall_boundary = true;
        }
}

void Grid::calculate_tangents() {
    for(int i=0;i<cols;i++)
        for(int j=0;j<rows;j++) {
            GridPoint *p = get_grid_point(i,j);
            p->tangent.x = p->normal.y;
            p->tangent.y = -p->normal.x;
        }
}
