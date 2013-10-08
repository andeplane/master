#include <mesh.h>
#include <CVector.h>


void Mesh::create_sequential(int num_vertices_) {
	num_vertices = num_vertices_;
	vertices.reserve(3*num_vertices);
	normals.reserve(3*num_vertices);
	colors.reserve(4*num_vertices);
}

void Mesh::render_triangles() {
	glBegin(GL_TRIANGLES);
	for(int i=0; i<num_vertices; i++) {
		glNormal3f(normals[3*i+0], normals[3*i+1], normals[3*i+2]);
		glVertex3f(vertices[3*i+0], vertices[3*i+1], vertices[3*i+2]);
		glColor4f(colors[4*i+0], colors[4*i+1], colors[4*i+2], colors[4*i+3]);
	}
	glEnd();
}

void Mesh::render_vbo() {
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_buffers[0]);
	glVertexPointer(3, GL_FLOAT, 0, (char*)NULL);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_buffers[1]);
	glNormalPointer(GL_FLOAT, 0, (char*)NULL);

	glBindBuffer(GL_ARRAY_BUFFER, vbo_buffers[2]);
	glColorPointer(4, GL_FLOAT, 0, (char*)NULL);

	glDrawArrays(GL_TRIANGLES, 0, num_vertices);
}

void Mesh::build_vbo() {
	// generate a buffer for our triangle
    glGenBuffers(3, vbo_buffers);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_buffers[0]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*num_vertices, &vertices[0], GL_STATIC_DRAW);
    
    glBindBuffer(GL_ARRAY_BUFFER, vbo_buffers[1]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*num_vertices, &normals[0], GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, vbo_buffers[2]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*4*num_vertices, &colors[0], GL_STATIC_DRAW);
}

void Mesh::enable_blend(bool inverse) {
	glEnable(GL_BLEND);
	if(inverse) glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	else glBlendFunc(GL_ONE, GL_ONE);
}

void Mesh::disable_blend() {
	glDisable(GL_BLEND);
}

void Mesh::generate_smooth_normals() {
	vector<float> normals_original = normals;

	for(int i=0; i<num_vertices; i++) {
		cout << "Calculating normals for i=" << i << endl;
		CVector v1(vertices[3*i+0], vertices[3*i+1], vertices[3*i+2]);
		CVector n1(normals_original[3*i+0], normals_original[3*i+1], normals_original[3*i+2]);
		int k1 = max(0, i-12500);
		int k2 = min(num_vertices-1, i+12500);
		for(int j=k1; j<k2; j++) {
			if(i != j) {
				CVector v2(vertices[3*j+0], vertices[3*j+1], vertices[3*j+2]);
				if((v1 - v2).Length2() < 1e-5) {
					CVector n2(normals_original[3*j+0], normals_original[3*j+1], normals_original[3*j+2]);
					n1 = n1+n2;
				}
			}
		} // end j

		n1 = n1.Normalize();
		normals[3*i+0] = n1.x;
		normals[3*i+1] = n1.y;
		normals[3*i+2] = n1.z;
	}
}

void Mesh::set_vertex(const int &index, float value) {
	if(vertices.size() <= index) vertices.resize(index+1);
	vertices[index] = value;
}

void Mesh::set_color(const int &index, float value) {
	if(colors.size() <= index) colors.resize(index+1);
	colors[index] = value;
}

void Mesh::set_normal(const int &index, float value) {
	if(normals.size() <= index) normals.resize(index+1);
	normals[index] = value;
}

void Mesh::add_vertex(CVector &v) {
	num_vertices++;
	vertices.push_back(v.x);
	vertices.push_back(v.y);
	vertices.push_back(v.z);
}

void Mesh::add_color(CVector &v, float alpha) {
	colors.push_back(v.x);
	colors.push_back(v.y);
	colors.push_back(v.z);
	colors.push_back(alpha);
}

void Mesh::add_normal(CVector &v) {
	normals.push_back(v.x);
	normals.push_back(v.y);
	normals.push_back(v.z);
}