#include <copengl.h>
#include <ctexture.h>
#include <camera.h>

CTexture::CTexture(COpenGL *ogl)
{
    opengl = ogl;
}


void CTexture::create_sphere1(string name, int size) {
    CBitMap bmp;
    COpenGLTexture texture;

    bmp.width = size;
    bmp.height= size;
    bmp.data = new unsigned char[4*size*size];
    for (int i=0;i<size; i++) {
        for (int j=0;j<size; j++) {
            double x = (i-size/2)/(double)size;
            double y = (j-size/2)/(double)size;
            double dr = sqrt(x*x + y*y)*2;
            dr = min(dr,1.0);

            if(dr>=0.99) {
                bmp.data[4*i + 4*size*j  +0] = 0;
                bmp.data[4*i + 4*size*j  +1] = 0;
                bmp.data[4*i + 4*size*j  +2] = 0;
                bmp.data[4*i + 4*size*j  +3] = 0;
            } else {
                bmp.data[4*i + 4*size*j  +0] = (unsigned char)(255.0*(1 - 0.7*dr));
                bmp.data[4*i + 4*size*j  +1] = (unsigned char)(255.0*(1 - 0.7*dr));
                bmp.data[4*i + 4*size*j  +2] = (unsigned char)(255.0*(1 - 0.7*dr));
                bmp.data[4*i + 4*size*j  +3] = 255;
            }
        }
    }

    load_texture(&bmp, &texture, true);
    textures.push_back(texture);
    names.push_back(name);
}

void CTexture::load_texture(CBitMap* bmp, COpenGLTexture* texture, bool has_alpha) {
    glGenTextures(1, &texture->id);

    texture->has_alpha = true;
    glBindTexture(GL_TEXTURE_2D, texture->id);

    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );

    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);

    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

    if(has_alpha) gluBuild2DMipmaps( GL_TEXTURE_2D, GL_RGBA, bmp->width, bmp->height, GL_RGBA, GL_UNSIGNED_BYTE, bmp->data  );
    else gluBuild2DMipmaps( GL_TEXTURE_2D, GL_RGB, bmp->width, bmp->height, GL_RGB, GL_UNSIGNED_BYTE, bmp->data  );
}

void CTexture::render_billboards(vector<float> &positions) {
    Camera *camera = opengl->camera;
    glEnable(GL_LIGHT0);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    COpenGLTexture texture = textures[0];
    GLuint texture_id = texture.id;

    glEnable( GL_TEXTURE_2D );
    glBindTexture(GL_TEXTURE_2D, texture_id);
    glDepthMask(GL_TRUE);

    glEnable(GL_ALPHA_TEST);
    glAlphaFunc(GL_GREATER,0.9);

    glBegin(GL_QUADS);

    CVector left, up, right, direction, v0, v1, v2, v3;
    double cam_x = camera->position.x; double cam_y = camera->position.y; double cam_z = camera->position.z;

    CVector up_on_screen = opengl->coord_to_ray(0,opengl->window_height/2.0);
    direction = camera->target;

    left = direction.Cross(up_on_screen);
    up = (direction.Cross(left)).Normalize();
    right = (direction.Cross(up)).Normalize();

    v0 = (right + up);
    v1 = (right*-1 + up);
    v2 = (right*-1 + up*-1);
    v3 = (right + up*-1);
    float one_over_color_cutoff = 1.0/1000;
    float scale = 0.01;

    glNormal3f(direction.x, direction.y, direction.z);
    int num_particles = positions.size()/3;
    for(int index=0; index<num_particles; index++) {
        float x = positions[3*index+0];
        float y = positions[3*index+1];
        float z = positions[3*index+2];

        double delta_x = x - cam_x;
        double delta_y = y - cam_y;
        double delta_z = z - cam_z;

        double dr2 = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z;
        if(dr2 < 0.1) continue;
        // if(dr2 > dr2_max) continue;

        double cam_target_times_dr = delta_x*direction.x + delta_y*direction.y + delta_z*direction.z;
        if(cam_target_times_dr < 0) continue;

        double color_factor = max( (1-dr2*one_over_color_cutoff),0.3);
        // double color_list[7][3] = {{1,1,1},{230.0/255,230.0/255,0},{0,0,1},{1.0,1.0,1.0},{1,0,0},{9.0/255,92.0/255,0},{95.0/255,216.0/255,250.0/255}};
        glColor4f(color_factor*230.0/255,color_factor*230.0/255,0, 1.0);
        // glColor4f(1.0, 1.0, 1.0, 1.0);

        glTexCoord2f(0,0);
        glVertex3f(v0.x*scale + x,v0.y*scale + y,v0.z*scale + z);
        glTexCoord2f(1,0);
        glVertex3f(v1.x*scale + x,v1.y*scale + y,v1.z*scale + z);
        glTexCoord2f(1,1);
        glVertex3f(v2.x*scale + x,v2.y*scale + y,v2.z*scale + z);
        glTexCoord2f(0,1);
        glVertex3f(v3.x*scale + x,v3.y*scale + y,v3.z*scale + z);
    }

    glEnd();
    glDisable(GL_BLEND);
    // glEnable(GL_CULL_FACE);
    glDisable(GL_ALPHA_TEST);

    glDepthMask(GL_TRUE);
    glDisable( GL_TEXTURE_2D );
    glColor4f(1.0,1.0,1.0,1.0);
}
