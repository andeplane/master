#include <mesh.h>
#include <CVector.h>
#include <iomanip>

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

void Mesh::generate_smooth_normals(map<string, vector<int> > &vertex_map ) {
	vector<float> normals_original = normals;
	char vertex_key[100];
	// Loop over atoms
    int old_progress = -1;
	for(int i=0; i<num_vertices; i++) {
		// if(num_vertices >= 100 && !(i%(num_vertices/100))) {
  //           printf("Creating smooth normal vectors: %d%% \n",(100*i)/num_vertices);
  //           fflush(stdout);
  //       }

		if (i % 1000 == 0) { // put couts at the top to avoid getting killed by "continue" further down
            int progress = floor(i/double(num_vertices) / 0.1);
            if (progress > old_progress) {
                cout << "Creating smooth normal vectors ";
                cout << "[";
                cout << setfill('#') << setw(progress) << "";
                cout << setfill('-') << setw(10-progress) << "";
                cout << "]" << endl;
                old_progress = progress;
            }
        }

		CVector v1(vertices[3*i+0], vertices[3*i+1], vertices[3*i+2]);
		CVector n1(normals_original[3*i+0], normals_original[3*i+1], normals_original[3*i+2]);
		sprintf(vertex_key, "%.2f%.2f%.2f",vertices[3*i+0], vertices[3*i+1], vertices[3*i+2]);
		string key = vertex_key;
		std::map<string, vector<int> >::iterator iterator = vertex_map.find(key);
		vector<int> &neighbor_triangle_indices = iterator->second;

		for(int j=0; j<neighbor_triangle_indices.size(); j++) {
			int neighbor_triangle_index = neighbor_triangle_indices[j];
			if(i != neighbor_triangle_index) {
				CVector v2(vertices[3*neighbor_triangle_index+0], vertices[3*neighbor_triangle_index+1], vertices[3*neighbor_triangle_index+2]);
				if((v1 - v2).Length2() < 1e-5) {
					CVector n2(normals_original[3*neighbor_triangle_index+0], normals_original[3*neighbor_triangle_index+1], normals_original[3*neighbor_triangle_index+2]);
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