#include "GL\glew.h"
#include "GL\freeglut.h"
#include <fstream>
#include <iostream>
#include "fem2D.h"
#include "fem3D.h"
#include "fem_timer.h"
using namespace std;
int step = 0;
float xRot = 0.0f;
float yRot = 0.f;
float xTrans = 0;
float yTrans = 0;
float zTrans = 0;
int ox;
int oy;
int buttonState;
float xRotLength = 0.0f;
float yRotLength = 0.0f;
float window_width = 600;
float window_height = 600;
int s_dimention = 3;
bool stop = true;
bool screenshot = false;
mesh2D mesh2d;
mesh3D mesh3d;

bool WriteBitmapFile(int width, int height, const std::string& file_name, unsigned char* bitmapData)
{
    BITMAPFILEHEADER bitmapFileHeader;
    memset(&bitmapFileHeader, 0, sizeof(BITMAPFILEHEADER));
    bitmapFileHeader.bfSize = sizeof(BITMAPFILEHEADER);
    bitmapFileHeader.bfType = 0x4d42;   //BM  
    bitmapFileHeader.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);

    BITMAPINFOHEADER bitmapInfoHeader;
    memset(&bitmapInfoHeader, 0, sizeof(BITMAPINFOHEADER));
    bitmapInfoHeader.biSize = sizeof(BITMAPINFOHEADER);
    bitmapInfoHeader.biWidth = width;
    bitmapInfoHeader.biHeight = height;
    bitmapInfoHeader.biPlanes = 1;
    bitmapInfoHeader.biBitCount = 24;
    bitmapInfoHeader.biCompression = BI_RGB;
    bitmapInfoHeader.biSizeImage = width * abs(height) * 3;

    //////////////////////////////////////////////////////////////////////////  
    FILE* filePtr;
    unsigned char tempRGB;
    int imageIdx;

    for (imageIdx = 0; imageIdx < (int)bitmapInfoHeader.biSizeImage; imageIdx += 3)
    {
        tempRGB = bitmapData[imageIdx];
        bitmapData[imageIdx] = bitmapData[imageIdx + 2];
        bitmapData[imageIdx + 2] = tempRGB;
    }

    filePtr = fopen(file_name.c_str(), "wb");
    if (NULL == filePtr)
    {
        return false;
    }

    fwrite(&bitmapFileHeader, sizeof(BITMAPFILEHEADER), 1, filePtr);

    fwrite(&bitmapInfoHeader, sizeof(BITMAPINFOHEADER), 1, filePtr);

    fwrite(bitmapData, bitmapInfoHeader.biSizeImage, 1, filePtr);

    fclose(filePtr);
    return true;
}

void SaveScreenShot(int width, int height, const std::string& file_name)
{
    int data_len = height * width * 3;      // bytes
    void* screen_data = malloc(data_len);
    memset(screen_data, 0, data_len);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, screen_data);

    WriteBitmapFile(width, height, file_name + ".bmp", (unsigned char*)screen_data);

    free(screen_data);
}

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

void draw_box3D(float ox, float oy, float oz, float width, float height, float length)
{
    glLineWidth(2.5f);
    glColor3f(0.8f, 0.8f, 0.8f);

    glBegin(GL_LINES);

    glVertex3f(ox, oy, oz);
    glVertex3f(ox + width, oy, oz);

    glVertex3f(ox, oy, oz);
    glVertex3f(ox, oy + height, oz);

    glVertex3f(ox, oy, oz);
    glVertex3f(ox, oy, oz + length);

    glVertex3f(ox + width, oy, oz);
    glVertex3f(ox + width, oy + height, oz);

    glVertex3f(ox + width, oy + height, oz);
    glVertex3f(ox, oy + height, oz);

    glVertex3f(ox, oy + height, oz + length);
    glVertex3f(ox, oy, oz + length);

    glVertex3f(ox, oy + height, oz + length);
    glVertex3f(ox, oy + height, oz);

    glVertex3f(ox + width, oy, oz);
    glVertex3f(ox + width, oy, oz + length);

    glVertex3f(ox, oy, oz + length);
    glVertex3f(ox + width, oy, oz + length);

    glVertex3f(ox + width, oy + height, oz);
    glVertex3f(ox + width, oy + height, oz + length);

    glVertex3f(ox + width, oy + height, oz + length);
    glVertex3f(ox + width, oy, oz + length);

    glVertex3f(ox, oy + height, oz + length);
    glVertex3f(ox + width, oy + height, oz + length);

    glEnd();
}

void draw_mesh2D()
{
    glEnable(GL_DEPTH_TEST);
    glLineWidth(1.5f);
    glColor3f(0.8f, 0.1f, 0.8f);
    glBegin(GL_LINES);
    for (int i = 0; i < mesh2d.triangleNum; i++) {
        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][0]][0], mesh2d.vertexes[mesh2d.triangles[i][0]][1], 0.0f);
        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][1]][0], mesh2d.vertexes[mesh2d.triangles[i][1]][1], 0.0f);

        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][0]][0], mesh2d.vertexes[mesh2d.triangles[i][0]][1], 0.0f);
        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][2]][0], mesh2d.vertexes[mesh2d.triangles[i][2]][1], 0.0f);

        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][1]][0], mesh2d.vertexes[mesh2d.triangles[i][1]][1], 0.0f);
        glVertex3f(mesh2d.vertexes[mesh2d.triangles[i][2]][0], mesh2d.vertexes[mesh2d.triangles[i][2]][1], 0.0f);
    }
    glEnd();
}

void draw_mesh3D()
{
    glEnable(GL_DEPTH_TEST);
    glLineWidth(1.5f);
    glColor3f(0.8f, 0.1f, 0.8f);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < mesh3d.tetrahedraNum; i++) {
        for (int j = 0; j < 4; j++) {
            glVertex3f(mesh3d.vertexes[mesh3d.tetrahedras[i][j]][0], mesh3d.vertexes[mesh3d.tetrahedras[i][j]][1], mesh3d.vertexes[mesh3d.tetrahedras[i][j]][2]);
            glVertex3f(mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 1) % 4]][0], mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 1) % 4]][1], mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 1) % 4]][2]);
            glVertex3f(mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 2) % 4]][0], mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 2) % 4]][1], mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 2) % 4]][2]);
        }
    }
    glEnd();
    glColor3f(0.8f, 0.8f, 0.1f);
    //glDisable(GL_DEPTH_TEST);
    glLineWidth(1.5f);
    glBegin(GL_LINES);
    double offset = 1.0;
    for (int i = 0; i < mesh3d.tetrahedraNum; i++) {
        for (int j = 0; j < 4; j++) {
            glVertex3f(mesh3d.vertexes[mesh3d.tetrahedras[i][j]][0] * offset, mesh3d.vertexes[mesh3d.tetrahedras[i][j]][1] * offset, mesh3d.vertexes[mesh3d.tetrahedras[i][j]][2] * offset);
            glVertex3f(mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 1) % 4]][0] * offset, mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 1) % 4]][1] * offset, mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 1) % 4]][2] * offset);

            glVertex3f(mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 1) % 4]][0] * offset, mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 1) % 4]][1] * offset, mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 1) % 4]][2] * offset);
            glVertex3f(mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 2) % 4]][0] * offset, mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 2) % 4]][1] * offset, mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 2) % 4]][2] * offset);

            glVertex3f(mesh3d.vertexes[mesh3d.tetrahedras[i][j]][0] * offset, mesh3d.vertexes[mesh3d.tetrahedras[i][j]][1] * offset, mesh3d.vertexes[mesh3d.tetrahedras[i][j]][2] * offset);
            glVertex3f(mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 2) % 4]][0] * offset, mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 2) % 4]][1] * offset, mesh3d.vertexes[mesh3d.tetrahedras[i][(j + 2) % 4]][2] * offset);
        }
    }
    glEnd();
}

void draw_Scene2D() {
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);   //sf ±³¾°ÑÕÉ«
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    draw_box2D(-1, -1, 2, 2);
    draw_mesh2D();
    glutSwapBuffers();
}

void draw_Scene3D() {
    glEnable(GL_DEPTH_TEST);

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);   //sf ±³¾°ÑÕÉ«
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);

    glPushMatrix();
    glTranslatef(xTrans, yTrans, zTrans);
    glRotatef(xRot, 1.0f, 0.0f, 0.0f);
    glRotatef(yRot, 0.0f, 1.0f, 0.0f);
    //glTranslatef(-0.5, -0.5, -0.5);


    draw_box3D(-1, -1, -1, 2, 2, 2);
    draw_mesh3D();

    glPopMatrix();
    glutSwapBuffers();
    //glFlush();
}
double mfsum = 0;
double total_time = 0;
int total_cg_iterations = 0;
int total_newton_iterations = 0;
void display(void)
{
    if (s_dimention == 2) {
        draw_Scene2D();
    }
    else if (s_dimention == 3) {
        draw_Scene3D();
    }

    if (stop) return;
    
    if(screenshot)
    {
        std::stringstream ss;
        ss << "saveScreen/step_";
        ss.fill('0');
        ss.width(5);
        ss << step;
        std::string file_path = ss.str();
        SaveScreenShot(window_width, window_height, file_path);
        step++;
    }

    bool isImplicit = true;
    if (s_dimention == 2) {
        if (isImplicit) {
            fem_implicit2D(mesh2d);
        }
        else {
            fem_explicit2D(mesh2d);
        }
    }
    else if (s_dimention == 3) {     
        if (isImplicit) {
            fem_implicit3D(mesh3d);
            /*double tempS = mfsum;
            HighResolutionTimerForWin timer;
            timer.set_start();
            Projected_Newton3D(mesh3d, mfsum, total_cg_iterations, total_newton_iterations);
            timer.set_end();
            total_time += timer.get_millisecond()/1000;
            cout << "cost time:         "<<total_time << endl;
            if (abs(mfsum - tempS) < 1e-3) {
                ofstream output("simulationData.txt");
                output << "total_time: " << total_time << endl;
                output << "total_cg_iterations: " << total_cg_iterations << endl;
                output << "total_newton_iterations: " << total_newton_iterations << endl;
                output.close();
                exit(0);
            }*/
        }
        else {
            fem_explicit3D(mesh3d);
        }
    }
    
    if (step == 20000) {
        exit(0);
    }
}

void init(void)
{
    // select clearing color: purple
    glClearColor(0.0, 0.0, 0.0, 1.0);

    if (s_dimention == 2) {
        initMesh2D(mesh2d, 1, 0.2);
    }
    else if (s_dimention == 3) {
        initMesh3D(mesh3d, 1, 0.2);
    }

    //glewInit();

    glViewport(0, 0, window_width, window_height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(45.0, (float)window_width / window_height, 10.0f, 500.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0.0f, 0.0f, -3.0f);
    //glTranslatef(0.5f, 0.5f, -3.0f);
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
    //glTranslatef(0.5f, 0.5f, -4.0f);
}

void keyboard_func(unsigned char key, int x, int y)
{
    if (key == 'w')
    {
        zTrans += .06f;
    }

    if (key == 's')
    {
        zTrans -= .06f;
    }

    if (key == 'a')
    {
        xTrans += .06f;
    }

    if (key == 'd')
    {
        xTrans -= .06f;
    }

    if (key == 'q')
    {
        yTrans -= .06f;
    }

    if (key == 'e')
    {
        yTrans += .06f;
    }

    if (key == ' ')
    {
        stop = !stop;
    }

    if (key == '/')
    {
        screenshot = !screenshot;
    }
    glutPostRedisplay();
}

void special_keyboard_func(int key, int x, int y)
{
    glutPostRedisplay();
}

void mouse_func(int button, int state, int x, int y)
{
    if (state == GLUT_DOWN)
    {
        buttonState = 1;
    }
    else if (state == GLUT_UP)
    {
        buttonState = 0;
    }

    ox = x; oy = y;

    glutPostRedisplay();
}

void motion_func(int x, int y)
{
    float dx, dy;
    dx = (float)(x - ox);
    dy = (float)(y - oy);

    if (buttonState == 1)
    {
        xRot += dy / 5.0f;
        yRot += dx / 5.0f;
    }

    ox = x; oy = y;

    glutPostRedisplay();
}

int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(window_width, window_height);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("FEM");
    init();
    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);
    glutDisplayFunc(display);


    //glutDisplayFunc(display_func);
    glutReshapeFunc(reshape_func);
    glutKeyboardFunc(keyboard_func);
    //glutSpecialFunc(special_keyboard_func);
    glutMouseFunc(mouse_func);
    glutMotionFunc(motion_func);
    glutIdleFunc(idle_func);


    glutMainLoop();
    //return 0;
}

