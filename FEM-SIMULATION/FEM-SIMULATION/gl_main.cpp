#include "GL\glew.h"
#include "GL\freeglut.h"

#include <iostream>
#include "fem2D.h"
#include "fem3D.h"
using namespace std;

int s_dimention = 2;

mesh2D mesh2d;
mesh3D mesh3d;

void draw_box2D(float ox, float oy, float width, float height)
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

void draw_mesh2D()
{
    glLineWidth(1.5f);
    glColor3f(0.8f, 0.1f, 0.8f);
    glBegin(GL_LINES);
    for (int i = 0; i < mesh2d.triangleNum; i++) {
        //float a[] = { circle.vertex[circle.index[i][0]][0], circle.vertex[circle.index[i][0]][1], 0.0f };
        //float b[] = { circle.vertex[circle.index[i][1]][0], circle.vertex[circle.index[i][1]][1], 0.0f };
        //float c[] = { circle.vertex[circle.index[i][2]][0], circle.vertex[circle.index[i][2]][1], 0.0f };
        //glVertex2fv(a);
        //glVertex2fv(b);
        //glVertex2fv(c);
        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][0]][0], mesh2d.vertexes[mesh2d.triangles[i][0]][1], 0.0f);
        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][1]][0], mesh2d.vertexes[mesh2d.triangles[i][1]][1], 0.0f);

        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][0]][0], mesh2d.vertexes[mesh2d.triangles[i][0]][1], 0.0f);
        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][2]][0], mesh2d.vertexes[mesh2d.triangles[i][2]][1], 0.0f);

        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][1]][0], mesh2d.vertexes[mesh2d.triangles[i][1]][1], 0.0f);
        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][2]][0], mesh2d.vertexes[mesh2d.triangles[i][2]][1], 0.0f);
    }
    glEnd();
}

void draw_Scene2D() {
    glClear(GL_COLOR_BUFFER_BIT);
    draw_box2D(-1, -1, 2, 2);
    draw_mesh2D();
    glFlush();
}

void display(void)
{
    

    if (s_dimention == 2) {
        draw_Scene2D();
        fem_implicit2D(mesh2d);
    }
    else if (s_dimention == 3) {

    }
}

void init(void)
{
    // select clearing color: purple
    glClearColor(0.0, 0.0, 0.0, 0.0);

    if (s_dimention == 2) {
        initMesh2D(mesh2d, 1, 0.2);
    }
    else if (s_dimention == 3) {
        initMesh3D(mesh3d, 1, 0.2);
    }
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

