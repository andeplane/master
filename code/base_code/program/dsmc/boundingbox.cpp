#include "boundingbox.h"

BoundingBox::BoundingBox(int particle_index, System *system_)
{
    system = system_;

    corner[0] = system->r[3*particle_index+0] - system->sticky_particle_radius;
    corner[1] = system->r[3*particle_index+1] - system->sticky_particle_radius;
    corner[2] = system->r[3*particle_index+2] - system->sticky_particle_radius;
    size[0] = 2*system->sticky_particle_radius;
    size[1] = 2*system->sticky_particle_radius;
    size[2] = 2*system->sticky_particle_radius;
    particle_indices.push_back(particle_index);
}

void BoundingBox::add_particle(int particle_index) {
    double *r = &system->r[3*particle_index];
    double x_min = min(corner[0], r[0]);
    double x_max = max(corner[0] + size[0], r[0]);

    double y_min = min(corner[1], r[1]);
    double y_max = max(corner[1] + size[1], r[1]);

    double z_min = min(corner[2], r[2]);
    double z_max = max(corner[2] + size[2], r[2]);

    particle_indices.push_back(particle_index);

    corner[0] = x_min;
    corner[1] = y_min;
    corner[2] = z_min;

    size[0]= x_max - x_min;
    size[1]= y_max - y_min;
    size[2]= z_max - z_min;
}

bool BoundingBox::point_is_inside_box(double *r) {
    return (r[2] >= corner[2] && r[2] <= corner[2]+size[2] && r[0] >= corner[0] && r[0] <= corner[0]+size[0] && r[1] >= corner[1] && r[1] <= corner[1]+size[1]);
}

bool BoundingBox::did_collide_with_a_particle(double *r) {
    double sticky_radius_squared = 4*system->sticky_particle_radius*system->sticky_particle_radius;

    for(int n=0; n<particle_indices.size(); n++) {
        int index = particle_indices[n];
        double dx = r[0] - system->r[3*index+0];
        double dy = r[1] - system->r[3*index+1];
        double dz = r[2] - system->r[3*index+2];
        double dr2 = dx*dx + dy*dy + dz*dz;
        if(dr2<sticky_radius_squared) return true;
    }

    return false;
}
