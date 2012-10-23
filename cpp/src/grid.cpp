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
        return new GridPoint(true);
    }

    return &points[j + i*rows];
}

GridPoint* Grid::get_grid_point(const double x, const double y) {
    int i =  clamp((int)(x/system->width*cols),0,cols-1);
    int j =  clamp((int)(y/system->height*rows),0,rows-1);

    return get_grid_point(i,j);
}

GridPoint* Grid::get_grid_point(const vec r, const int idx) {
    if(idx==594) {
        // cout << "YEAH" << endl;
    }

    int i =  clamp((int)(r(0)/system->width*cols),0,cols-1);
    int j =  clamp((int)(r(1)/system->height*rows),0,rows-1);

    return get_grid_point(i,j);
}

#define DEBUG

void Grid::calculate_normals() {
#ifdef DEBUG
    int t_i = 524;
    int t_j = 293;
    mat m(3,3);
#endif
    for(int i=0;i<cols;i++)
        for(int j=0;j<rows;j++) {
            GridPoint *p = get_grid_point(i,j);

            for(int k=-1;k<=1;k++)
                for(int l=-1;l<=1;l++) {
                    p->normal.x += get_grid_point(i+k,j+l)->is_wall*k;
                    p->normal.y -= get_grid_point(i+k,j+l)->is_wall*l;
#ifdef DEBUG
                    if(i==t_i && j==t_j) {
                        m(l+1,k+1) = get_grid_point(i+k,j+l)->is_wall;
                        // cout << get_grid_point(i+k,j+l)->is_wall*l << endl;
                        cout << p->normal.y << endl;
                    }
#endif
                }
#ifdef DEBUG
            if(i==t_i && j==t_j)
                cout << "Normal vector before norm: " << endl << p->normal << endl;
#endif
            p->normal = p->normal.Normalize();
#ifdef DEBUG
            if(i==t_i && j==t_j) {
                cout << "Normal vector: " << endl << p->normal << endl;
            }
#endif
        }
#ifdef DEBUG
    cout << m << endl;
#endif
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
            p->tangent.y = p->normal.x;
        }
}
