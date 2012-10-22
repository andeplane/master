#include "grid.h"

#define min(a,b)                      (((a) < (b)) ? (a) : (b))
#define max(a,b)                      (((a) > (b)) ? (a) : (b))
#define clamp(value, lb, ub)          max( lb, min( ub, value ))

Grid::Grid(mat M)
{
    points.reserve(M.n_elem);
    cols = M.n_cols;
    rows = M.n_rows;

    for(int i=0;i<cols;i++)
        for(int j=0;j<rows;j++) {
            GridPoint p(M(j,i));
            points.push_back(p);
        }

    calculate_normals();
}

GridPoint* Grid::get_grid_point(const int i, const int j) {
    if(i < 0 || i >= cols || j < 0 || j >= rows) {
        return new GridPoint(true);
    }

    return &points[j + i*rows];
}

GridPoint* Grid::get_grid_point(const double x, const double y) {
    int i =  clamp((int)(x/system->width*cols),0,cols-1);
    int j =  clamp((int)(y/system->height*rows),0,rows-1);
    return get_grid_point(i,j);
}

void Grid::calculate_normals() {
    for(int i=0;i<cols;i++)
        for(int j=0;j<rows;j++) {
            GridPoint *p = get_grid_point(i,j);
            for(int k=-1;k<=1;k++)
                for(int l=-1;l<=1;l++) {
                    p->normal.x += get_grid_point(i+k,j+l)->is_wall*k;
                    p->normal.y += get_grid_point(i+k,j+l)->is_wall*l;
                }
        }
}
