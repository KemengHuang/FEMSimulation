#include "GL\glew.h"
#include "GL\freeglut.h"

#include <iostream>
#include "fem.h"
using namespace std;

//mesh_circle mesh;
mesh2D mesh;

void draw_box(float ox, float oy, float width, float height)
{
    glLineWidth(2.5f);
    glColor3f(0.8f, 0.8f, 0.8f);

    glBegin(GL_LINES);

    glVertex3f(ox, oy, 0);
    glVertex3f(ox + width, oy, 0);

    glVertex3f(ox, oy, 0);
    glVertex3f(ox, oy + height, 0);

    glVertex3f(ox + width, oy, 0);
    glVertex3f(ox + width, oy + height, 0);

    glVertex3f(ox + width, oy + height, 0);
    glVertex3f(ox, oy + height, 0);

    glEnd();
}

void draw_mesh()
{
    glLineWidth(1.5f);
    glColor3f(0.8f, 0.1f, 0.8f);
    glBegin(GL_LINES);
    for (int i = 0; i < mesh.triangleNum; i++) {
        //float a[] = { circle.vertex[circle.index[i][0]][0], circle.vertex[circle.index[i][0]][1], 0.0f };
        //float b[] = { circle.vertex[circle.index[i][1]][0], circle.vertex[circle.index[i][1]][1], 0.0f };
        //float c[] = { circle.vertex[circle.index[i][2]][0], circle.vertex[circle.index[i][2]][1], 0.0f };
        //glVertex2fv(a);
        //glVertex2fv(b);
        //glVertex2fv(c);
        glVertex3f(mesh.vertexes[mesh.triangles[i][0]][0], mesh.vertexes[mesh.triangles[i][0]][1], 0.0f);
        glVertex3f(mesh.vertexes[mesh.triangles[i][1]][0], mesh.vertexes[mesh.triangles[i][1]][1], 0.0f);

        glVertex3f(mesh.vertexes[mesh.triangles[i][0]][0], mesh.vertexes[mesh.triangles[i][0]][1], 0.0f);
        glVertex3f(mesh.vertexes[mesh.triangles[i][2]][0], mesh.vertexes[mesh.triangles[i][2]][1], 0.0f);

        glVertex3f(mesh.vertexes[mesh.triangles[i][1]][0], mesh.vertexes[mesh.triangles[i][1]][1], 0.0f);
        glVertex3f(mesh.vertexes[mesh.triangles[i][2]][0], mesh.vertexes[mesh.triangles[i][2]][1], 0.0f);
    }
    glEnd();
}

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT);

    draw_box(-1, -1, 2, 2);
    //glColor3f(1.0, 1.0, 0.0);
    //glBegin(GL_TRIANGLES);
    //glVertex3f(To.vertexs[0](0) / 5 - 0, To.vertexs[0](1) / 5 - 1, 0.0f);
    //glVertex3f(To.vertexs[1](0) / 5 - 0, To.vertexs[1](1) / 5 - 1, 0.0f);
    //glVertex3f(To.vertexs[2](0) / 5 - 0, To.vertexs[2](1) / 5 - 1, 0.0f);
    //glColor3f(1.0, 1.0, 0.0);
    //glVertex3f(Tc.vertexs[0](0) / 5 + 0, Tc.vertexs[0](1) / 5 + 0, 0.0f);
    //glVertex3f(Tc.vertexs[1](0) / 5 + 0, Tc.vertexs[1](1) / 5 + 0, 0.0f);
    //glVertex3f(Tc.vertexs[2](0) / 5 + 0, Tc.vertexs[2](1) / 5 + 0, 0.0f);
    
    draw_mesh();

    glFlush();
    //fem_explicit();
    //Projected_Newton_2DTest();
    fem_explicit2D(mesh);
    //Projected_Newton2D(mesh);
}

void init(void)
{
    // select clearing color: purple
    glClearColor(0.0, 0.0, 0.0, 0.0);
    initMesh(mesh, 1, 0.2);
    // initialize viewing values
    //glMatrixMode(GL_PROJECTION);
    //glLoadIdentity();
    //glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
}




void idle_func()
{
    glutPostRedisplay();
}

void reshape_func(GLint width, GLint height)
{
    //window_width = width;
    //window_height = height;

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(45.0, (float)width / height, 0.001, 500.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0.0f, 0.0f, -3.0f);
}

void keyboard_func(unsigned char key, int x, int y)
{
    glutPostRedisplay();
}

void special_keyboard_func(int key, int x, int y)
{
    glutPostRedisplay();
}

void mouse_func(int button, int state, int x, int y)
{
    glutPostRedisplay();
}

void motion_func(int x, int y)
{
    glutPostRedisplay();
}

int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("FEM");
    init();
    glutDisplayFunc(display);


    //glutDisplayFunc(display_func);
    glutReshapeFunc(reshape_func);
    //glutKeyboardFunc(keyboard_func);
    //glutSpecialFunc(special_keyboard_func);
    //glutMouseFunc(mouse_func);
    //glutMotionFunc(motion_func);
    glutIdleFunc(idle_func);


    glutMainLoop();
    //return 0;
}

