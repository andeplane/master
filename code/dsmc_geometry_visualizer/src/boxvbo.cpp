#include <boxvbo.h>

void BoxVBO::create_vbo() {
	vertices = new float[6*4*3]; // 6 faces, 4 points each face
    indices  = new int  [4*6];

    glGenBuffers(1, &vertices_id);
    glBindBuffer(GL_ARRAY_BUFFER, vertices_id);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*6*num_indices, NULL, GL_DYNAMIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(GLfloat)*6*num_indices, vertices);

    glGenBuffers(1, &indices_id);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indices_id);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*num_indices, NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, sizeof(GLuint)*num_indices, indices);
}

void BoxVBO::draw() {
	glNormal3f(direction.x, direction.y, direction.z);
    // glNormalPointer(GL_FLOAT, 6*sizeof(GLfloat), (float*)(3*sizeof(GLfloat)));

    glBindBuffer(GL_ARRAY_BUFFER, vertices_id);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*3*num_vertices, NULL, GL_DYNAMIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(GLfloat)*3*num_vertices, vertices);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indices_id);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*num_vertices, NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, sizeof(GLuint)*num_vertices, indices);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);

    glVertexPointer(3, GL_FLOAT, 3*sizeof(GLfloat), NULL);
    glColorPointer (4, GL_FLOAT, 4, NULL);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indices_id);
    glBindBuffer(GL_ARRAY_BUFFER, vertices_id);
    glDrawElements(GL_QUADS, num_vertices, GL_UNSIGNED_INT, 0);

    glDisableClientState(GL_VERTEX_ARRAY);
}