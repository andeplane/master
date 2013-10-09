#pragma once
#include <vector>
#include <map>
#include <string>
#include <GL/glfw.h>      // Include OpenGL Framework library

using std::vector;
using std::map;
using std::string;

class CVector;

using std::vector;

class Mesh {
	vector<float> vertices;
	vector<int>   indices;
	vector<float> normals;
	vector<float> colors;
	// VBO
	GLuint vbo_buffers[3];
public:
	Mesh() {
		num_vertices = 0;
		vertices.reserve(20000000);
		normals.reserve(20000000);
		colors.reserve(20000000);
	}
	int num_vertices;
	void create_sequential(int num_vertices_);
	void render_triangles();
	void set_vertex(const int &index, float value);
	void set_color(const int &index, float value);
	void set_normal(const int &index, float value);
	void add_vertex(CVector &v);
	void add_color(CVector &v, float alpha);
	void add_normal(CVector &v);
	void build_vbo();
	void render_vbo();
	void generate_smooth_normals(map<string, vector<int> > &vertex_map );
	void enable_blend(bool inverse);
	void disable_blend();
};