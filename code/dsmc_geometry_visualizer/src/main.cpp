

#include <DSMCOpenGL.h>
#include <iostream>
#include <string>
#include <GLUT/glut.h>
#include <Camera.h>       // Include our Camera header so we can work with Camera objects
#include <FpsManager.hpp> // Include our FpsManager class
#include <Vec3.hpp>       // Include our Vec3 class
#include <CIniFile.h>
#include <time.h>
#include <marching.cpp>
#include <mesh.h>

// Specify default namespace for commonly used elements
using std::string;
using std::cout;
using std::endl;

typedef enum {
    voxel_type_empty = 0,
    voxel_type_wall = 1,
    voxel_type_boundary = 2
} voxel_type;

DSMCOpenGL dsmcopengl;
bool paused = false;
long num_voxels = 0;
int Nx, Ny, Nz;

bool draw_normals = false;
bool draw_surface_only = false;
bool marching_cubes = true;
vector<unsigned char> point_type;
vector<vector<float> > points;
vector<vector<float> > normals;
Mesh mesh;
float t = 0;

void draw_marching_cubes() {
    if(!marching_cubes) return;
    // mesh.enable_blend(false);
    mesh.render_vbo();
    // mesh.disable_blend();
}

void draw_points() {
    glBegin(GL_POINTS);
    glColor4f(0.0, 1.0, 0.0, 1.0);
    glPointSize(50.0);
    for(int n=0; n<points.size(); n++) {
        if(draw_surface_only && point_type[n] != voxel_type_boundary) continue;

        float x = points[n][0];
        float y = points[n][1];
        float z = points[n][2];

        glVertex3f(x,y,z);
    }

    glEnd();
    if(draw_normals) {
        glBegin(GL_LINES);
        glColor4f(1.0, 0.0, 0.0, 1.0);
        for(int n=0; n<points.size(); n++) {
            if(point_type[n] == voxel_type_boundary) {
                float x = points[n][0];
                float y = points[n][1];
                float z = points[n][2];

                float nx = 0.1*normals[n][0];
                float ny = 0.1*normals[n][1];
                float nz = 0.1*normals[n][2];

                glVertex3f(x,y,z);
                glVertex3f(x+nx,y+ny,z+nz);
            }
            
        }
        glEnd();
    }
}

// Function to draw our scene
void drawScene()
{
    // Clear the screen and depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 
    // Reset the matrix
    glMatrixMode(GL_MODELVIEW);
    glDrawBuffer(GL_BACK_RIGHT);
    glLoadIdentity();
 
    // Move the camera to our location in space
    glRotatef(dsmcopengl.camera->get_rot_x(), 1.0f, 0.0f, 0.0f); // Rotate our camera on the x-axis (looking up and down)
    glRotatef(dsmcopengl.camera->get_rot_y(), 0.0f, 1.0f, 0.0f); // Rotate our camera on the  y-axis (looking left and right)

    // Translate the ModelView matrix to the position of our camera - everything should now be drawn relative
    // to this position!
    glTranslatef( -dsmcopengl.camera->position.x, -dsmcopengl.camera->position.y, -dsmcopengl.camera->position.z );
    
    // draw_points();
    dsmcopengl.set_standard_light(t);
    // draw_points();
    draw_marching_cubes();
    
    
    CVector up_on_screen = dsmcopengl.coord_to_ray(0,dsmcopengl.window_height/2.0);
    
    // ----- Stop Drawing Stuff! ------ 
    glfwSwapBuffers(); // Swap the buffers to display the scene (so we don't have to watch it being drawn!)
}

// Callback function to handle mouse movements
void handle_mouse_move(int mouse_x, int mouse_y) {
    dsmcopengl.camera->handle_mouse_move(mouse_x, mouse_y);
}

// Callback function to handle keypresses
void handle_keypress(int theKey, int theAction) {
    // If a key is pressed, toggle the relevant key-press flag
    if (theAction == GLFW_PRESS)
    {
        switch (theKey)
        {
        case 'W':
            dsmcopengl.camera->holding_forward = true;
            break;
        case 'S':
            dsmcopengl.camera->holding_backward = true;
            break;
        case 'A':
            dsmcopengl.camera->holding_left_strafe = true;
            break;
        case 'D':
            dsmcopengl.camera->holding_right_strafe = true;
            break;
        case '[':
            dsmcopengl.fps_manager.set_target_fps(dsmcopengl.fps_manager.get_target_fps() - 10);
            break;
        case ']':
            dsmcopengl.fps_manager.set_target_fps(dsmcopengl.fps_manager.get_target_fps() + 10);
            break;
        case 'N':
            draw_normals = !draw_normals;
            break;
        case 'M':
            marching_cubes = !marching_cubes;
            break;
        case 'Q':
            draw_surface_only = !draw_surface_only;
            break;
        case ' ':
            paused = !paused;
            break;
        }
    }
    else // If a key is released, toggle the relevant key-release flag
    {
        switch (theKey)
        {
        case 'W':
            dsmcopengl.camera->holding_forward = false;
            break;
        case 'S':
            dsmcopengl.camera->holding_backward = false;
            break;
        case 'A':
            dsmcopengl.camera->holding_left_strafe = false;
            break;
        case 'D':
            dsmcopengl.camera->holding_right_strafe = false;
            break;
        default:
            // Do nothing...
            break;
        }
    }
}

void read_matrix(string filename) {
    ifstream file (filename.c_str(), ios::in | ios::binary);
    if(!file.is_open()) {
        cout << "Error, could not open file " << filename << endl;
        exit(1);
    }

    file.read (reinterpret_cast<char*>(&Nx), sizeof(int));
    file.read (reinterpret_cast<char*>(&Ny), sizeof(int));
    file.read (reinterpret_cast<char*>(&Nz), sizeof(int));
    num_voxels = Nx*Ny*Nz;

    unsigned char *voxels_ = new unsigned char[num_voxels];
    float *normals_   = new float[3*num_voxels];
    float *tangents1_ = new float[3*num_voxels];
    float *tangents2_ = new float[3*num_voxels];

    file.read (reinterpret_cast<char*>(voxels_), num_voxels*sizeof(unsigned char));
    file.read (reinterpret_cast<char*>(normals_), 3*num_voxels*sizeof(float));
    file.read (reinterpret_cast<char*>(tangents1_), 3*num_voxels*sizeof(float));
    file.read (reinterpret_cast<char*>(tangents2_), 3*num_voxels*sizeof(float));
    file.close();

    for(int i=0; i<Nx; i++) {
        for(int j=0; j<Ny; j++) {
            for(int k=0; k<Nz; k++) {
                int index = ((i+Nx)%Nx) + ((j+Ny)%Ny)*Nx+ ((k+Nz)%Nz)*Nx*Ny;
                float x = 10*(float)i / Nx;
                float y = 10*(float)j / Ny;
                float z = 10*(float)k / Nz;
                if(voxels_[index] > voxel_type_empty) {
                    // We have a wall point
                    vector<float> pos;
                    vector<float> normal;
                    pos.resize(3);
                    normal.resize(3);

                    pos[0] = x; pos[1] = y; pos[2] = z;
                    normal[0] = normals_[3*index + 0];
                    normal[1] = normals_[3*index + 1];
                    normal[2] = normals_[3*index + 2];

                    points.push_back(pos);
                    normals.push_back(normal);
                    point_type.push_back(voxels_[index]);
                }
            }
        }
    }

    delete voxels_;
    delete normals_;
    delete tangents1_;
    delete tangents2_;
    
    cout << "We have " << points.size() << " wall points." << endl;
    int total_triangles = 0;

    int mesh_vertex_counter = 0;
    int mesh_normal_counter = 0;
    int mesh_color_counter = 0;
    CVector direction = CVector(-1,-1,-1).Normalize();
    int max_num_triangles = 0;
    TRIANGLE triangles[100];

    for(int i=0; i<Nx-1; i++) {
        cout << i << " / " << Nx << endl;
        for(int j=0; j<Ny-1; j++) {
            for(int k=0; k<Nz-1; k++) {
                GRIDCELL cell;
                int vertex_count = 0;

                vector<int> x; x.resize(8,0);
                x[0] = 0; x[1] = 1; x[2] = 1; x[3] = 0;
                x[4] = 0; x[5] = 1; x[6] = 1; x[7] = 0;

                vector<int> y; y.resize(8,0);
                y[0] = 0; y[1] = 0; y[2] = 0; y[3] = 0;
                y[4] = 1; y[5] = 1; y[6] = 1; y[7] = 1;

                vector<int> z; z.resize(8,0);
                z[0] = 0; z[1] = 0; z[2] = 1; z[3] = 1;
                z[4] = 0; z[5] = 0; z[6] = 1; z[7] = 1;

                for(int balle=0; balle<8; balle++) {
                    int index = i+x[balle] + (j+y[balle])*Nx + (k+z[balle])*Nx*Ny;
                    float xx = 10*(float)(i-Nx/2.0+x[balle]) / Nx;
                    float yy = 10*(float)(j-Ny/2.0+y[balle]) / Ny;
                    float zz = 10*(float)(k-Nz/2.0+z[balle]) / Nz;
                    
                    cell.p[vertex_count].x = xx;
                    cell.p[vertex_count].y = yy;
                    cell.p[vertex_count].z = zz;
                    cell.val[vertex_count] = (voxels_[index] == voxel_type_boundary);
                    vertex_count++;
                }
                int num = Polygonise(cell, 1, triangles);
                max_num_triangles = max(max_num_triangles, num);

                for(int n=0; n<num; n++) {
                    CVector normal = (triangles[n].p[0] - triangles[n].p[1]).Cross((triangles[n].p[0] - triangles[n].p[2])).Normalize();
                    float d = direction.Dot(normal);
                    if(d<0.3) d=0.3;

                    CVector color(i/(float)Nx, 1 - j/(float)Ny, (k+j)/(2.0*Ny));
                    
                    for(int point = 0; point < 3; point++) {
                        mesh.add_vertex(triangles[n].p[point]);
                        mesh.add_normal(normal);
                        mesh.add_color(color, 0.2);
                    }
                }
                total_triangles += num;
            }
        }
    }
    

    cout << "Marching cube done, we have " << total_triangles << " triangles." << endl;
    cout << "Marching cubes gave max_num_triangles = " << max_num_triangles << endl;
}
 
// Fire it up...
int main(int argc, char **argv)
{
    cout << "Controls: Use WSAD and the mouse to move around!" << endl;
    
    CIniFile ini;
    
    ini.load("dsmc_geometry_visualizer.ini");
    string world_file = ini.getstring("world_file");
    cout << world_file << endl;
    read_matrix(world_file);

    char *window_title = new char[1000];
    sprintf(window_title, "DSMC geometry Visualizer (DSMCGV) - [%.2f fps]", 60.0);
    dsmcopengl.initialize(ini.getint("screen_width"),ini.getint("screen_height"), string(window_title), handle_keypress, handle_mouse_move, false, 1);
    GLenum error = glewInit();

    mesh.build_vbo();

    bool running = true;
    while (running)
    {
        t += 0.01;
        // Calculate our camera movement
        dsmcopengl.camera->move();
 
        // Draw our scene
        drawScene();
 
        // exit if ESC was pressed or window was closed
        running = !glfwGetKey(GLFW_KEY_ESC) && glfwGetWindowParam(GLFW_OPENED);
        // Call our fps manager to limit the FPS
        dsmcopengl.fps_manager.enforce_fps();
        double fps = dsmcopengl.fps_manager.average_fps;
        // Calculate the current time in pico seconds to show in the title bar
        sprintf(window_title, "Molecular Dynamics Visualizer (MDV) - [%.2f fps]",fps);
        dsmcopengl.set_window_title(string(window_title));
    }
 
    // Clean up GLFW and exit
    glfwTerminate();
 
    return 0;
}