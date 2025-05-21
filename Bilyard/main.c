#include <windows.h>
#include <gl/gl.h>
#include <math.h>
#define pW 40
#define pH 60
#define FRICTION 0.97f
#define RESTITUTION 0.3f
#define B_RESTITUTION 0.83f
#define STRIKE_FORCE 0.62f
#define BALL_RADIUS 1.0f
#define MIN_SPEED 0.001f
#define STB_IMAGE_IMPLEMENTATION
#include "../stb-master/stb-master/stb_image.h"
#include <GL/glu.h>


float camX = 0.0f, camY = 0.0f, camZ = 60.0f;
float camYaw = -90.0f;
float camPitch = 0.0f;

LRESULT CALLBACK WindowProc(HWND, UINT, WPARAM, LPARAM);
void EnableOpenGL(HWND hwnd, HDC*, HGLRC*);
void DisableOpenGL(HWND, HDC, HGLRC);
GLuint ballTextures[16];
GLuint tableTexture;
GLuint woodtex;
int windowWidth = 1200;
int windowHeight = 950;
GLfloat lightPos[4] = {0.0f, 0.0f, 2.0f, 0.0f};



GLuint LoadTexture(const char* filename) {
    int width, height, channels;
    unsigned char* image = stbi_load(filename, &width, &height, &channels, STBI_rgb_alpha);

    if (!image) {
        printf("Failed to load texture: %s\n", filename);
        return 0;
    }

    GLuint textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0,
                 GL_RGBA, GL_UNSIGNED_BYTE, image);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    stbi_image_free(image);
    return textureID;
}

void Parall(float x, float y, float z, float dx, float dy, float dz)
{
    float vertices[] = {
        x,     y,     z,
        x + dx, y,     z,
        x,     y + dy, z,
        x + dx, y + dy, z,
        x,     y,     z + dz,
        x + dx, y,     z + dz,
        x,     y + dy, z + dz,
        x + dx, y + dy, z + dz
    };

    float normals[] = {
        -1,-1,-1, 1,-1,-1, -1,1,-1, 1,1,-1,
        -1,-1,1, 1,-1,1, -1,1,1, 1,1,1
    };

    float repeatX = dx / 16.0f;
    float repeatY = dy / 16.0f;

    float texCoords[] = {
        0.0, 0.0,
        repeatX, 0.0,
        0.0, repeatY,
        repeatX, repeatY,
        0.0, 0.0,
        repeatX, 0.0,
        0.0, repeatY,
        repeatX, repeatY
    };

    GLubyte strips[] = {
        4, 5, 6, 7,
        1, 0, 3, 2,
        5, 1, 7, 3,
        0, 4, 2, 6,
        6, 7, 2, 3,
        0, 1, 4, 5
    };

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    glVertexPointer(3, GL_FLOAT, 0, vertices);
    glNormalPointer(GL_FLOAT, 0, normals);
    glTexCoordPointer(2, GL_FLOAT, 0, texCoords);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    for(int i = 0; i < 6; i++) {
        glDrawElements(GL_TRIANGLE_STRIP, 4, GL_UNSIGNED_BYTE, &strips[i*4]);
    }

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
}

void Parall22(float x, float y, float z, float dx, float dy, float dz)
{
    float vertices[] = {
        x - dy,     y,     z,
        x + dx + dy, y,     z,
        x,     y + dy, z,
        x + dx, y + dy, z,
        x - dy,     y,     z + dz,
        x + dx + dy, y,     z + dz,
        x,     y + dy, z + dz,
        x + dx, y + dy, z + dz
    };

    float normals[] = {
        -1,-1,-1, 1,-1,-1, -1,1,-1, 1,1,-1,
        -1,-1,1, 1,-1,1, -1,1,1, 1,1,1
    };

    float repeatX = dx / 2.0f;
    float repeatY = dy / 2.0f;

    float texCoords[] = {
        0.0, 0.0,
        repeatX, 0.0,
        0.0, repeatY,
        repeatX, repeatY,
        0.0, 0.0,
        repeatX, 0.0,
        0.0, repeatY,
        repeatX, repeatY
    };

    GLubyte strips[] = {
        4, 5, 6, 7,
        1, 0, 3, 2,
        5, 1, 7, 3,
        0, 4, 2, 6,
        6, 7, 2, 3,
        0, 1, 4, 5
    };

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    glVertexPointer(3, GL_FLOAT, 0, vertices);
    glNormalPointer(GL_FLOAT, 0, normals);
    glTexCoordPointer(2, GL_FLOAT, 0, texCoords);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    for(int i = 0; i < 6; i++) {
        glDrawElements(GL_TRIANGLE_STRIP, 4, GL_UNSIGNED_BYTE, &strips[i*4]);
    }

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
}

void ParS2(float x, float y, float z, float dx, float dy, float dz)
{
    float vertices[] = {
        x - 1,     y,     z,
        x + dx + 1, y,     z,
        x,     y + dy, z,
        x + dx, y + dy, z,
        x - 1,     y,     z + dz,
        x + dx + 1, y,     z + dz,
        x,     y + dy, z + dz,
        x + dx, y + dy, z + dz
    };

    float normals[] = {
        -1,-1,-1, 1,-1,-1, -1,1,-1, 1,1,-1,
        -1,-1,1, 1,-1,1, -1,1,1, 1,1,1
    };

    float repeatX = dx / 2.0f;
    float repeatY = dy / 2.0f;

    float texCoords[] = {
        0.0, 0.0,
        repeatX, 0.0,
        0.0, repeatY,
        repeatX, repeatY,
        0.0, 0.0,
        repeatX, 0.0,
        0.0, repeatY,
        repeatX, repeatY
    };

    GLubyte strips[] = {
        4, 5, 6, 7,
        1, 0, 3, 2,
        5, 1, 7, 3,
        0, 4, 2, 6,
        6, 7, 2, 3,
        0, 1, 4, 5
    };

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    glVertexPointer(3, GL_FLOAT, 0, vertices);
    glNormalPointer(GL_FLOAT, 0, normals);
    glTexCoordPointer(2, GL_FLOAT, 0, texCoords);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    for(int i = 0; i < 6; i++) {
        glDrawElements(GL_TRIANGLE_STRIP, 4, GL_UNSIGNED_BYTE, &strips[i*4]);
    }

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
}


void Parall33(float x, float y, float z, float dx, float dy, float dz)
{
    float vertices[] = {
        x,     y-dx,     z,
        x + dx, y,     z,
        x,     y + dy + dx, z,
        x + dx, y + dy, z,
        x,     y - dx,     z + dz,
        x + dx, y,     z + dz,
        x,     y + dy + dx, z + dz,
        x + dx, y + dy, z + dz
    };

    float normals[] = {
        -1,-1,-1, 1,-1,-1, -1,1,-1, 1,1,-1,
        -1,-1,1, 1,-1,1, -1,1,1, 1,1,1
    };

    float repeatX = dx / 2.0f;
    float repeatY = dy / 2.0f;

    float texCoords[] = {
        0.0, 0.0,
        repeatX, 0.0,
        0.0, repeatY,
        repeatX, repeatY,
        0.0, 0.0,
        repeatX, 0.0,
        0.0, repeatY,
        repeatX, repeatY
    };

    GLubyte strips[] = {
        4, 5, 6, 7,
        1, 0, 3, 2,
        5, 1, 7, 3,
        0, 4, 2, 6,
        6, 7, 2, 3,
        0, 1, 4, 5
    };

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    glVertexPointer(3, GL_FLOAT, 0, vertices);
    glNormalPointer(GL_FLOAT, 0, normals);
    glTexCoordPointer(2, GL_FLOAT, 0, texCoords);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    for(int i = 0; i < 6; i++) {
        glDrawElements(GL_TRIANGLE_STRIP, 4, GL_UNSIGNED_BYTE, &strips[i*4]);
    }

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
}

void ParS(float x, float y, float z, float dx, float dy, float dz)
{
    float vertices[] = {
        x,     y-1,     z,
        x + dx, y,     z,
        x,     y + dy + dx, z,
        x + dx, y + dy, z,
        x,     y - dx,     z + dz,
        x + dx, y,     z + dz,
        x,     y + dy + 1, z + dz,
        x + dx, y + dy, z + dz
    };

    float normals[] = {
        -1,-1,-1, 1,-1,-1, -1,1,-1, 1,1,-1,
        -1,-1,1, 1,-1,1, -1,1,1, 1,1,1
    };

    float repeatX = dx / 2.0f;
    float repeatY = dy / 2.0f;

    float texCoords[] = {
        0.0, 0.0,
        repeatX, 0.0,
        0.0, repeatY,
        repeatX, repeatY,
        0.0, 0.0,
        repeatX, 0.0,
        0.0, repeatY,
        repeatX, repeatY
    };

    GLubyte strips[] = {
        4, 5, 6, 7,
        1, 0, 3, 2,
        5, 1, 7, 3,
        0, 4, 2, 6,
        6, 7, 2, 3,
        0, 1, 4, 5
    };

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    glVertexPointer(3, GL_FLOAT, 0, vertices);
    glNormalPointer(GL_FLOAT, 0, normals);
    glTexCoordPointer(2, GL_FLOAT, 0, texCoords);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    for(int i = 0; i < 6; i++) {
        glDrawElements(GL_TRIANGLE_STRIP, 4, GL_UNSIGNED_BYTE, &strips[i*4]);
    }

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
}

void DrawCircle(int count, float r, float x0, float y0, float z0)
{
    float x, y;
    float da = M_PI * 2.0 / count;
    glBegin(GL_TRIANGLE_FAN);
        glColor3f(0,0,0);
        glVertex3f(x0,y0,z0);
        for(int i = 0; i <= count; i++)
        {
            x = x0 + r*sin(da * i);
            y = y0 + r*cos(da * i);
            glVertex3f(x, y, z0);
        }
    glEnd();

}

void HalfwCylin(int segments, float radius, float thickness, float x0, float y0, float z0, float height, float alfa, float beta)
{
    float da = M_PI * beta * radius / segments;
    float dz = height / segments;

    glBegin(GL_QUADS);

    for (int j = 0; j < segments; j++)
    {
        float z1 = z0 + dz * j;
        float z2 = z0 + dz * (j + 1);

        for (int i = 0; i < segments; i++)
        {
            float angle1 = da * i;
            float angle2 = da * (i + 1);

            float rot_angle1 = angle1 + alfa;
            float rot_angle2 = angle2 + alfa;

            float sin_a1 = sin(rot_angle1);
            float cos_a1 = cos(rot_angle1);
            float sin_a2 = sin(rot_angle2);
            float cos_a2 = cos(rot_angle2);

            glNormal3f(sin_a1, cos_a1, 0.0f);
            glVertex3f(x0 + radius * sin_a1, y0 + radius * cos_a1, z1);
            glVertex3f(x0 + radius * sin_a1, y0 + radius * cos_a1, z2);
            glVertex3f(x0 + radius * sin_a2, y0 + radius * cos_a2, z2);
            glVertex3f(x0 + radius * sin_a2, y0 + radius * cos_a2, z1);

            glNormal3f(-sin_a1, -cos_a1, 0.0f);
            glVertex3f(x0 + (radius - thickness) * sin_a1, y0 + (radius - thickness) * cos_a1, z1);
            glVertex3f(x0 + (radius - thickness) * sin_a1, y0 + (radius - thickness) * cos_a1, z2);
            glVertex3f(x0 + (radius - thickness) * sin_a2, y0 + (radius - thickness) * cos_a2, z2);
            glVertex3f(x0 + (radius - thickness) * sin_a2, y0 + (radius - thickness) * cos_a2, z1);

            if (i == 0 || i == segments - 1)
            {
                float side_nx = (i == 0) ? -cos_a1 : cos_a1;
                float side_ny = (i == 0) ? sin_a1 : -sin_a1;

                glNormal3f(side_nx, side_ny, 0.0f);
                glVertex3f(x0 + radius * sin_a1, y0 + radius * cos_a1, z1);
                glVertex3f(x0 + radius * sin_a1, y0 + radius * cos_a1, z2);
                glVertex3f(x0 + (radius - thickness) * sin_a1, y0 + (radius - thickness) * cos_a1, z2);
                glVertex3f(x0 + (radius - thickness) * sin_a1, y0 + (radius - thickness) * cos_a1, z1);
            }
        }
    }

    for (int i = 0; i < segments; i++)
    {
        float angle1 = da * i;
        float angle2 = da * (i + 1);

        float rot_angle1 = angle1 + alfa;
        float rot_angle2 = angle2 + alfa;

        float sin_a1 = sin(rot_angle1);
        float cos_a1 = cos(rot_angle1);
        float sin_a2 = sin(rot_angle2);
        float cos_a2 = cos(rot_angle2);

        glNormal3f(0.0f, 0.0f, -1.0f);
        glVertex3f(x0 + radius * sin_a1, y0 + radius * cos_a1, z0);
        glVertex3f(x0 + radius * sin_a2, y0 + radius * cos_a2, z0);
        glVertex3f(x0 + (radius - thickness) * sin_a2, y0 + (radius - thickness) * cos_a2, z0);
        glVertex3f(x0 + (radius - thickness) * sin_a1, y0 + (radius - thickness) * cos_a1, z0);

        glNormal3f(0.0f, 0.0f, 1.0f);
        glVertex3f(x0 + radius * sin_a1, y0 + radius * cos_a1, z0 + height);
        glVertex3f(x0 + (radius - thickness) * sin_a1, y0 + (radius - thickness) * cos_a1, z0 + height);
        glVertex3f(x0 + (radius - thickness) * sin_a2, y0 + (radius - thickness) * cos_a2, z0 + height);
        glVertex3f(x0 + radius * sin_a2, y0 + radius * cos_a2, z0 + height);
    }

    glEnd();
}


void Quad(float x, float y, float z, float dx, float dy)
{
    float normal[] = {-1,-1,3, 1,-1,3, 1,1,3, -1,1,3};
    float vertex[] = {x, y, z, x + dx, y, z, x+dx, y+dy, z, x, y+dy, z};
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, vertex);
        glNormalPointer(GL_FLOAT, 0, normal);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
}

void DrawSphere(float radius, int sectors, int stacks, float x0, float y0, float z0)
{
    const float PI = 3.14159265f;
    float sectorStep = 2 * PI / sectors;
    float stackStep = PI / stacks;

    int vertexCount = (stacks + 1) * (sectors + 1);
    float vertices[vertexCount * 3];
    float normals[vertexCount * 3];
    float texCoords[vertexCount * 2];

    int index = 0, normalIndex = 0, texIndex = 0;
    for(int i = 0; i <= stacks; ++i)
    {
        float stackAngle = PI / 2 - i * stackStep;
        float xy = radius * cosf(stackAngle);
        float z = radius * sinf(stackAngle);

        for(int j = 0; j <= sectors; ++j)
        {
            float sectorAngle = j * sectorStep;
            float cosSector = cosf(sectorAngle);
            float sinSector = sinf(sectorAngle);

            vertices[index++] = x0 + xy * cosSector;
            vertices[index++] = y0 + xy * sinSector;
            vertices[index++] = z0 + z;

            normals[normalIndex++] = cosSector * cosf(stackAngle);
            normals[normalIndex++] = sinSector * cosf(stackAngle);
            normals[normalIndex++] = sinf(stackAngle);

            texCoords[texIndex++] = (float)j / sectors;
            texCoords[texIndex++] = (float)i / stacks;
        }
    }

    int indexCount = stacks * (sectors + 1) * 2;
    unsigned int indices[indexCount];

    index = 0;
    for(int i = 0; i < stacks; ++i)
    {
        int k1 = i * (sectors + 1);
        int k2 = k1 + sectors + 1;

        for(int j = 0; j <= sectors; ++j)
        {
            indices[index++] = k1 + j;
            indices[index++] = k2 + j;
        }
    }

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);

    glVertexPointer(3, GL_FLOAT, 0, vertices);
    glNormalPointer(GL_FLOAT, 0, normals);
    glTexCoordPointer(2, GL_FLOAT, 0, texCoords);

    for(int i = 0; i < stacks; ++i)
    {
        glDrawElements(GL_TRIANGLE_STRIP, (sectors + 1) * 2,
                      GL_UNSIGNED_INT, &indices[(sectors + 1) * 2 * i]);
    }

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
}

typedef struct {
    float x, y, z;
    float dx, dy;
    float r;
    float vx, vy;
    float rotationAngleX;
    float rotationAngleZ;
    int isMoving;
    int isPocketed;
} TBall;


TBall bitok;
TBall balls[15];



void TBall_init(TBall *obj, float x1, float y1, float z1, float dx1, float dy1, float r1) {
    obj->x = x1;
    obj->y = y1;
    obj->z = z1;
    obj->dx = dx1;
    obj->dy = dy1;
    obj->r = r1;
    obj->vx = 0;
    obj->vy = 0;
    obj->isMoving = 0;
    obj->rotationAngleX = 0;
    obj->rotationAngleZ = 0;
    obj->isPocketed = 0;
}



int CheckBallCollision(TBall *ball1, TBall *ball2) {
    float dx = ball1->x - ball2->x;
    float dy = ball1->y - ball2->y;
    float distanceSquared = dx*dx + dy*dy;
    float minDistance = 2 * BALL_RADIUS;
    return distanceSquared <= minDistance * minDistance;
}


void HandleBallCollision(TBall *ball1, TBall *ball2) {

    float dx = ball2->x - ball1->x;
    float dy = ball2->y - ball1->y;
    float distance = sqrt(dx*dx + dy*dy);


    float nx = dx / distance;
    float ny = dy / distance;


    float dvx = ball2->vx - ball1->vx;
    float dvy = ball2->vy - ball1->vy;


    float velocity_along_normal = dvx * nx + dvy * ny;


    if (velocity_along_normal > 0) return;
    float current_restitution = B_RESTITUTION;
    if (ball1 == &bitok || ball2 == &bitok) {
        current_restitution *= 1.1f;
    }

    float impulse = -(1 + current_restitution) * velocity_along_normal / 2.0f;

    ball1->vx -= impulse * nx;
    ball1->vy -= impulse * ny;
    ball2->vx += impulse * nx;
    ball2->vy += impulse * ny;


    ball1->vx *= 0.96f;
    ball1->vy *= 0.96f;
    ball2->vx *= 0.96f;
    ball2->vy *= 0.96f;


    float overlap = 2 * BALL_RADIUS - distance;
    if (overlap > 0) {
        float moveX = overlap * nx * 0.5f;
        float moveY = overlap * ny * 0.5f;

        ball1->x -= moveX;
        ball1->y -= moveY;
        ball2->x += moveX;
        ball2->y += moveY;
    }

    ball1->isMoving = 1;
    ball2->isMoving = 1;
}


void HandleWallCollision(TBall *ball)
{

    float left = -15.55f + BALL_RADIUS;
    float right = 15.55f - BALL_RADIUS;
    float bottom = -30.1f + BALL_RADIUS;
    float top = 30.1f - BALL_RADIUS;


    if (ball->x < left || ball->x > right)
    {
        ball->vx = -ball->vx * RESTITUTION * 0.8;
        ball->x = (ball->x < left) ? left : right;
        ball->isMoving = 1;
    }


    if (ball->y < bottom || ball->y > top)
    {
        ball->vy = -ball->vy * RESTITUTION * 0.8;
        ball->y = (ball->y < bottom) ? bottom : top;
        ball->isMoving = 1;
    }
}

void CheckPocketCollision(TBall *ball) {
    struct Pocket {
        float x, y, radius;
    } pockets[] = {
        {-15.39f, -29.82f, 1.2f},
        {15.39f, -29.82f, 1.2f},
        {-15.39f, 29.82f, 1.2f},
        {15.39f, 29.82f, 1.2f},
        {-15.65f, 0.0f, 1.25f},
        {15.65f, 0.0f, 1.25f}
    };

    for (int i = 0; i < 6; i++) {
        float dx = ball->x - pockets[i].x;
        float dy = ball->y - pockets[i].y;
        float distance = sqrt(dx*dx + dy*dy);
        if (distance < pockets[i].radius) {
            if (ball == &bitok) {

                ball->x = 0.0f;
                ball->y = -14.55f;
                ball->vx = ball->vy = 0.0f;
                ball->isMoving = 0;
            }
            else if (!ball->isPocketed) {
                ball->isPocketed = 1;
                ball->vx = ball->vy = 0;
                ball->isMoving = 0;
            }
            break;
        }

    }
}

void UpdateBallPosition(TBall *ball) {
    if (ball->isPocketed) return;


    ball->vx *= FRICTION;
    ball->vy *= FRICTION;

    float speed = sqrt(ball->vx * ball->vx + ball->vy * ball->vy);
    if (speed < 0.1f) {
        ball->vx *= 0.94f;
        ball->vy *= 0.94f;
    }

    float oldX = ball->x;
    float oldY = ball->y;

    ball->x += ball->vx;
    ball->y += ball->vy;

    HandleWallCollision(ball);

    float distanceMoved = sqrt((ball->x - oldX)*(ball->x - oldX) +
                     (ball->y - oldY)*(ball->y - oldY));
    if (distanceMoved > 0) {
        float circumference = 2.0f * M_PI * BALL_RADIUS;
        float angleChange = (distanceMoved / circumference) * 360.0f;

        if (ball->vx != 0 || ball->vy != 0) {
            ball->rotationAngleZ = -atan2(ball->vx, ball->vy) * 180.0f / M_PI;
            ball->rotationAngleX -= angleChange * 0.95;

            if (ball->rotationAngleX > 360.0f) ball->rotationAngleX -= 360.0f;
            if (ball->rotationAngleX < 0.0f) ball->rotationAngleX += 360.0f;
        }
    }


    if (speed < MIN_SPEED) {
        ball->vx = ball->vy = 0;
        ball->isMoving = 0;

    }

    CheckPocketCollision(ball);
}

void StrikeCueBall(TBall *cueBall, float mouseX, float mouseY)
{
    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];
    GLdouble winX, winY, winZ;
    GLdouble worldX, worldY, worldZ;

    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);

    winX = mouseX;
    winY = viewport[3] - mouseY;

    glReadPixels(winX, winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);

    GLdouble worldStartX, worldStartY, worldStartZ;
    gluUnProject(winX, winY, 0.0, modelview, projection, viewport, &worldStartX, &worldStartY, &worldStartZ);

    GLdouble worldEndX, worldEndY, worldEndZ;
    gluUnProject(winX, winY, 1.0, modelview, projection, viewport, &worldEndX, &worldEndY, &worldEndZ);

    float rayDirX = worldEndX - worldStartX;
    float rayDirY = worldEndY - worldStartY;

    float t = (1.0 - worldStartZ) / (worldEndZ - worldStartZ);
    float targetX = worldStartX + t * rayDirX;
    float targetY = worldStartY + t * rayDirY;

    float dx = targetX - cueBall->x;
    float dy = targetY - cueBall->y;
    float distance = sqrt(dx*dx + dy*dy);

    if (distance > 0) {
        dx /= distance;
        dy /= distance;

        const float MAX_DISTANCE = 30.0f;
        const float MIN_FORCE = 0.7f;
        const float MAX_FORCE = 2.5f;

        if (distance > MAX_DISTANCE) distance = MAX_DISTANCE;

        float forceFactor = MIN_FORCE + pow(distance / MAX_DISTANCE, 0.7f) * (MAX_FORCE - MIN_FORCE);

        cueBall->vx = -dx * STRIKE_FORCE * forceFactor;
        cueBall->vy = -dy * STRIKE_FORCE * forceFactor;
        cueBall->isMoving = 1;
    }
}



void UpdateAllBalls() {

    for (int i = 0; i < 16; i++) {
        TBall *ball = (i == 15) ? &bitok : &balls[i];
        if (!ball->isPocketed) {
            UpdateBallPosition(ball);
        }
    }

    for (int i = 0; i < 16; i++) {
        TBall *ball1 = (i == 15) ? &bitok : &balls[i];
        if (ball1->isPocketed) continue;

        for (int j = i + 1; j < 16; j++) {
            TBall *ball2 = (j == 15) ? &bitok : &balls[j];
            if (ball2->isPocketed) continue;

            if (CheckBallCollision(ball1, ball2)) {
                HandleBallCollision(ball1, ball2);
            }
        }
    }
}

void Game_init()
{
    ballTextures[0] = LoadTexture("Шар1.jpg");
    ballTextures[1] = LoadTexture("Шар2.jpg");
    ballTextures[2] = LoadTexture("Шар3.jpg");
    ballTextures[3] = LoadTexture("Шар4.jpg");
    ballTextures[4] = LoadTexture("Шар5.jpg");
    ballTextures[5] = LoadTexture("Шар6.jpg");
    ballTextures[6] = LoadTexture("Шар7.jpg");
    ballTextures[7] = LoadTexture("Шар8.jpg");
    ballTextures[8] = LoadTexture("Шар9.jpg");
    ballTextures[9] = LoadTexture("Шар10.jpg");
    ballTextures[10] = LoadTexture("Шар11.jpg");
    ballTextures[11] = LoadTexture("Шар12.jpg");
    ballTextures[12] = LoadTexture("Шар13.jpg");
    ballTextures[13] = LoadTexture("Шар14.jpg");
    ballTextures[14] = LoadTexture("Шар15.jpg");
    ballTextures[15] = LoadTexture("Шар0.jpg");
    TBall_init(&bitok, 0.0, -14.55, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+3, -4, 21.63, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+11, -2, 21.63, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+6, 0, 21.63, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+13, 2, 21.63, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+9, 4, 21.63, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+8, -3, 19.86, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+10, -1, 19.86, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+4, 1, 19.86, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+0, 3, 19.86, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+5, -2, 18.09, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+7, 0, 18.09, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+2, 2, 18.09, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+12, -1, 16.32, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+14, 1, 16.32, 1.0, 0.0, 0.0, 1.0);
    TBall_init(balls+1, 0, 14.55, 1.0, 0.0, 0.0, 1.0);
    tableTexture = LoadTexture("Сукно3.jpg");
    woodtex = LoadTexture("Дерево2.jpg");
}


void TBall_Show(TBall *obj, int textureIndex)
{
    if (obj->isPocketed) return;
    glPushMatrix();
    glTranslatef(obj->x, obj->y, obj->z);

    glRotatef(obj->rotationAngleZ, 0.0f, 0.0f, 1.0f);
    glRotatef(obj->rotationAngleX, 1.0f, 0.0f, 0.0f);

    glScalef(obj->r, obj->r, obj->r);


        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, ballTextures[textureIndex]);

        float matAmbient[] = {0.7f, 0.7f, 0.7f, 1.0f};
        float matDiffuse[] = {1.0f, 1.0f, 1.0f, 1.0f};
        float matSpecular[] = {1.0f, 1.0f, 1.0f, 1.0f};
        float matShininess = 100.0f;

        glMaterialfv(GL_FRONT, GL_AMBIENT, matAmbient);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, matDiffuse);
        glMaterialfv(GL_FRONT, GL_SPECULAR, matSpecular);
        glMaterialf(GL_FRONT, GL_SHININESS, matShininess);

        DrawSphere(1.0, 40, 40, 0, 0, 0);

        glDisable(GL_TEXTURE_2D);
    glPopMatrix();
}

void DrawBorderShadow() {

    glPushMatrix();


    glDisable(GL_LIGHTING);
    glColor4f(0.15f, 0.15f, 0.15f, 0.6f);

    ParS(-15.48, -27.7, 0, 1.08, 25.38, 0.05f);
    ParS(-15.48, 2.31, 0, 1.08, 25.45, 0.05f);
    ParS(14.42, -1.25, 0, 1.08, -27.55, 0.05f);
    ParS(14.42, 28.8, 0, 1.08, -27.55, 0.05f);
    ParS2(-13.27, 30.1, 0, 26.54, -1.08, 0.05f);
    ParS2(-13.27, -30.1, 0, 26.54, 1.08, 0.05f);

    glEnable(GL_LIGHTING);
    glPopMatrix();
}

void ShowWorld()
{
    glNormal3f(0,0,1);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, tableTexture);

    float tableAmbient[] = {0.1f, 0.15f, 0.1f, 1.0f};
    float tableDiffuse[] = {0.0f, 0.3f, 0.15f, 1.0f};
    float tableSpecular[] = {0.05f, 0.05f, 0.05f, 1.0f};
    float tableShininess = 5.0f;

    glMaterialfv(GL_FRONT, GL_AMBIENT, tableAmbient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, tableDiffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, tableSpecular);
    glMaterialf(GL_FRONT, GL_SHININESS, tableShininess);

    glEnable(GL_TEXTURE_COORD_ARRAY);

    glColor3f(1.0, 1.0, 1.0);
    Parall(-16.5, -30.5, -1, 33., 61., 1);

    glDisable(GL_TEXTURE_2D);
    glDisable(GL_TEXTURE_COORD_ARRAY);

    glMaterialfv(GL_FRONT, GL_AMBIENT, tableAmbient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, tableDiffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, tableSpecular);
    glMaterialf(GL_FRONT, GL_SHININESS, tableShininess);


    glDepthMask(GL_FALSE);
    glDisable(GL_DEPTH_TEST);
    DrawBorderShadow();
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);


    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, tableTexture);

    glMaterialfv(GL_FRONT, GL_AMBIENT, tableAmbient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, tableDiffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, tableSpecular);
    glMaterialf(GL_FRONT, GL_SHININESS, tableShininess);

    glEnable(GL_TEXTURE_COORD_ARRAY);

    glColor3f(1.0, 1.0, 1.0);
    Parall22(-13.27, -30.1, 0, 26.54, 1, 1.236);
    Parall33(-15.5, -27.8, 0, 1, 25.55, 1.236);
    Parall33(-15.5, 2.25, 0, 1, 25.55, 1.236);
    Parall33(14.5, -1.25, 0, 1, -27.55, 1.236);
    Parall33(14.5, 28.8, 0, 1, -27.55, 1.236);
    Parall22(13.27, 30.1, 0, -26.54, -1, 1.236);

    glDisable(GL_TEXTURE_2D);
    glDisable(GL_TEXTURE_COORD_ARRAY);

    glMaterialfv(GL_FRONT, GL_AMBIENT, tableAmbient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, tableDiffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, tableSpecular);
    glMaterialf(GL_FRONT, GL_SHININESS, tableShininess);
    DrawCircle(40, 1.2, -15.45, -30.0, 0.2);
    DrawCircle(40, 1.2, 15.45, -30.0, 0.2);
    DrawCircle(40, 1.2, -15.45, 30.0, 0.2);
    DrawCircle(40, 1.2, 15.45, 30.0, 0.2);
    DrawCircle(40, 1.25, -16.475, 0.0, 0.2);
    DrawCircle(40, 1.25, 16.475, 0.0, 0.2);


    glColor3f(0.5,0.5,0.5);
    HalfwCylin(30, 1.2, -0.2, -15.45, -30.0, 0.0, 1.236, 20.52, 1.2);
    HalfwCylin(30, 1.2, -0.2, 15.45, -30.0, 0.0, 1.236, 0.04, 1.2);
    HalfwCylin(30, 1.2, -0.2, -15.45, 30.00, 0.0, 1.236, 15.8, 1.2);
    HalfwCylin(30, 1.2, -0.2, 15.45, 30.00, 0.0, 1.236, 17.37, 1.2);
    HalfwCylin(30, 1.25, -0.2, -16.475, 0.0, 0.0, 1.236, -15.69, 0.8);
    HalfwCylin(30, 1.25, -0.2, 16.475, 0.0, 0.0, 1.236, -0.0, 0.8);
    glColor3f(0.4,0.15,0.07);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, woodtex);

    float woodAmbient[] = {0.2, 0.15, 0.1, 1.0};
    float woodDiffuse[] = {0.6, 0.4, 0.2, 1.0};
    float woodSpecular[] = {0.1, 0.1, 0.1, 1.0};
    float woodShininess = 10.0;

    glMaterialfv(GL_FRONT, GL_AMBIENT, woodAmbient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, woodDiffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, woodSpecular);
    glMaterialf(GL_FRONT, GL_SHININESS, woodShininess);

    glEnable(GL_TEXTURE_COORD_ARRAY);

    Parall(-14.27, -30.1, 0, 28.54, -1, 1.236);
    Parall(-15.5, -28.8, 0, -1, 27.55, 1.236);
    Parall(-15.5, 1.25, 0, -1, 27.55, 1.236);
    Parall(15.5, -1.25, 0, 1, -27.55, 1.236);
    Parall(15.5, 28.8, 0, 1, -27.55, 1.236);
    Parall(14.27, 31.1, 0, -28.54, -1, 1.236);

    glDisable(GL_TEXTURE_2D);
    glDisable(GL_TEXTURE_COORD_ARRAY);

    glMaterialfv(GL_FRONT, GL_AMBIENT, woodAmbient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, woodDiffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, woodSpecular);
    glMaterialf(GL_FRONT, GL_SHININESS, woodShininess);
    glColor3f(1,1,1);
    TBall_Show(&bitok, 15);
    TBall_Show(&balls[0], 0);
    TBall_Show(&balls[1], 1);
    TBall_Show(&balls[2], 2);
    TBall_Show(&balls[3], 3);
    TBall_Show(&balls[4], 4);
    TBall_Show(&balls[5], 5);
    TBall_Show(&balls[6], 6);
    TBall_Show(&balls[7], 7);
    TBall_Show(&balls[8], 8);
    TBall_Show(&balls[9], 9);
    TBall_Show(&balls[10], 10);
    TBall_Show(&balls[11], 11);
    TBall_Show(&balls[12], 12);
    TBall_Show(&balls[13], 13);
    TBall_Show(&balls[14], 14);
}




void MoveCamera()
{
    float moveSpeed = 0.3f;
    float rotateSpeed = 1.4f;

    if (GetKeyState('W') < 0) {
        camX += sin(camYaw * M_PI/180) * moveSpeed;
        camY += cos(camYaw * M_PI/180) * moveSpeed;
    }
    if (GetKeyState('S') < 0) {
        camX -= sin(camYaw * M_PI/180) * moveSpeed;
        camY -= cos(camYaw * M_PI/180) * moveSpeed;
    }
    if (GetKeyState('A') < 0) {
        camX -= cos(camYaw * M_PI/180) * moveSpeed;
        camY += sin(camYaw * M_PI/180) * moveSpeed;
    }
    if (GetKeyState('D') < 0) {
        camX += cos(camYaw * M_PI/180) * moveSpeed;
        camY -= sin(camYaw * M_PI/180) * moveSpeed;
    }
        if (GetKeyState('X') < 0) {

        camZ += 2 * moveSpeed;
    }
    if (GetKeyState('C') < 0) {
        camZ -= 2 * moveSpeed;
    }

    if (GetKeyState(VK_LEFT) < 0) camYaw -= rotateSpeed;
    if (GetKeyState(VK_RIGHT) < 0) camYaw += rotateSpeed;
    if (GetKeyState(VK_UP) < 0) camPitch += rotateSpeed;
    if (GetKeyState(VK_DOWN) < 0) camPitch -= rotateSpeed;

    if (camPitch > 85.0f) camPitch = 85.0f;
    if (camPitch < -85.0f) camPitch = -85.0f;

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glRotatef(camPitch, 1.0, 0.0, 0.0);
    glRotatef(camYaw, 0.0, 0.0, 1.0);
    glTranslatef(-camX, -camY, -camZ);
}

int WINAPI WinMain(HINSTANCE hInstance,
                   HINSTANCE hPrevInstance,
                   LPSTR lpCmdLine,
                   int nCmdShow)
{
    WNDCLASSEX wcex;
    HWND hwnd;
    HDC hDC;
    HGLRC hRC;
    MSG msg;
    BOOL bQuit = FALSE;
    float theta = 0.0f;

    /* register window class */
    wcex.cbSize = sizeof(WNDCLASSEX);
    wcex.style = CS_OWNDC;
    wcex.lpfnWndProc = WindowProc;
    wcex.cbClsExtra = 0;
    wcex.cbWndExtra = 0;
    wcex.hInstance = hInstance;
    wcex.hIcon = LoadIcon(NULL, IDI_APPLICATION);
    wcex.hCursor = LoadCursor(NULL, IDC_ARROW);
    wcex.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
    wcex.lpszMenuName = NULL;
    wcex.lpszClassName = "GLSample";
    wcex.hIconSm = LoadIcon(NULL, IDI_APPLICATION);;


    if (!RegisterClassEx(&wcex))
        return 0;

    /* create main window */
    hwnd = CreateWindowEx(0,
                          "GLSample",
                          "OpenGL Sample",
                          WS_OVERLAPPEDWINDOW,
                          CW_USEDEFAULT,
                          CW_USEDEFAULT,
                          windowWidth,
                          windowHeight,
                          NULL,
                          NULL,
                          hInstance,
                          NULL);

    ShowWindow(hwnd, nCmdShow);

    /* enable OpenGL for the window */
EnableOpenGL(hwnd, &hDC, &hRC);

glEnable(GL_BLEND);
glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
glEnable(GL_DEPTH_TEST);
glDepthFunc(GL_LEQUAL);
glViewport(0, 0, windowWidth, windowHeight);

glMatrixMode(GL_PROJECTION);
glLoadIdentity();
gluPerspective(45.0f, (float)windowWidth/windowHeight, 0.1f, 100.0f);
glMatrixMode(GL_MODELVIEW);


float aspectRatio = (float)windowWidth / (float)windowHeight;
if (windowWidth <= windowHeight) {
    glFrustum(-1.0, 1.0, -1.0/aspectRatio, 1.0/aspectRatio, 1.0, 100.0);
} else {
    glFrustum(-1.0*aspectRatio, 1.0*aspectRatio, -1.0, 1.0, 1.0, 100.0);
}

glTranslatef(0, 0, -18);
glMatrixMode(GL_MODELVIEW);
glLoadIdentity();
glEnable(GL_DEPTH_TEST);
glEnable(GL_BLEND);
glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);


    float lightDiffuse[] = {0.97f, 0.9f, 0.8f, 1.0f};
    float lightSpecular[] = {0.97f, 0.9f, 0.8f, 1.0f};
    float lightPosition[] = {0.0f, 0.0f, 2.0f, 0.0f};

    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);


    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
    /* program main loop */
    Game_init();
    while (!bQuit)
    {
        /* check for messages */
        if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
        {
            /* handle or dispatch messages */
            if (msg.message == WM_QUIT)
            {
                bQuit = TRUE;
            }
            else
            {
                TranslateMessage(&msg);
                DispatchMessage(&msg);
            }
        }
        else
        {
            /* OpenGL animation code goes here */

            glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);



            glPushMatrix();
                //glRotatef(theta, 0, 1, 0);
                float position[] = {0,0,2,0};
                glLightfv(GL_LIGHT0, GL_POSITION, position);
                glScalef(1.0, 1.0, 1.0);
                glColor3f(1,1,1);
                Quad(-53.1,-36.55,90,1,1);


            glPopMatrix();
            MoveCamera();
            ShowWorld();
            //theta += 0.1f;

            SwapBuffers(hDC);

            Sleep(16);  // ~60 FPS
            UpdateAllBalls();
        }
    }

    /* shutdown OpenGL */
    DisableOpenGL(hwnd, hDC, hRC);

    /* destroy the window explicitly */
    DestroyWindow(hwnd);

    return msg.wParam;
}

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
     switch (uMsg)
    {
case WM_SIZE: {
    windowWidth = LOWORD(lParam);
    windowHeight = HIWORD(lParam);
    glViewport(0, 0, windowWidth, windowHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, (float)windowWidth/windowHeight, 0.1f, 100.0f);
    glMatrixMode(GL_MODELVIEW);
    break;
}

        case WM_CLOSE:
            PostQuitMessage(0);
            break;

        case WM_DESTROY:
            return 0;

        case WM_KEYDOWN:
{
    switch (wParam)
    {
        case VK_ESCAPE:
            PostQuitMessage(0);
            break;
        // Остальные клавиши обрабатываются в MoveCamera()
    }
    break;
}
        case WM_LBUTTONDOWN: {
        int mouseX = LOWORD(lParam);
        int mouseY = HIWORD(lParam);
        StrikeCueBall(&bitok, mouseX, mouseY);

    break;
        }
        break;


        default:
            return DefWindowProc(hwnd, uMsg, wParam, lParam);
    }

    return 0;
}

void EnableOpenGL(HWND hwnd, HDC* hDC, HGLRC* hRC)
{
    PIXELFORMATDESCRIPTOR pfd;

    int iFormat;

    /* get the device context (DC) */
    *hDC = GetDC(hwnd);

    /* set the pixel format for the DC */
    ZeroMemory(&pfd, sizeof(pfd));

    pfd.nSize = sizeof(pfd);
    pfd.nVersion = 1;
    pfd.dwFlags = PFD_DRAW_TO_WINDOW |
                  PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
    pfd.iPixelType = PFD_TYPE_RGBA;
    pfd.cColorBits = 24;
    pfd.cDepthBits = 16;
    pfd.iLayerType = PFD_MAIN_PLANE;

    iFormat = ChoosePixelFormat(*hDC, &pfd);

    SetPixelFormat(*hDC, iFormat, &pfd);

    /* create and enable the render context (RC) */
    *hRC = wglCreateContext(*hDC);

    wglMakeCurrent(*hDC, *hRC);
}

void DisableOpenGL (HWND hwnd, HDC hDC, HGLRC hRC)
{
    wglMakeCurrent(NULL, NULL);
    wglDeleteContext(hRC);
    ReleaseDC(hwnd, hDC);
}

