#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

//Global Variables
const int WIDTH = 80;
const int HEIGHT = 23;
const int pixelsLength = 80*23;
float pixels[80*23];
const char * cols = " .,:;$#";

struct mat {
    float m[4][4];
};

struct vector {
    int x;
    int y;
    int z;
    int w;
};

struct tri {
    struct vector p1;
    struct vector p2;
    struct vector p3;
};

struct camera {
    float pNear;
    float pFar;
    float fov;
    float fovRad;
};

float sign(float in) {
    if (in > 0) {
        return 1;
    } else {
        return -1;
    }
}

int ffloor(float in) {
    if (in > 0) {
        return (int) in;
    } else {
        return (int) in - 1;
    }
}

int fceil(float in) {
    return ffloor(in) + 1;
}

int fround(float in) {
    if (in - ffloor(in) < 0.5) {
        return ffloor(in);
    } else {
        return fceil(in);
    }
}

void initOne(struct mat * in) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i == j) {
                in->m[i][j] = 1;
            } else {
                in->m[i][j] = 0;
            }
        }
    }
}

void initZero(struct mat * in) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            in->m[i][j] = 0;
        }
    }
}

void initCamera(struct camera *cam) {
    cam->fov = 3.141592/1.5;
    cam->pFar = 1000;
    cam->pNear = 0.1;
    cam->fovRad = 1/tanf(cam->fov/2);
}

void draw() {
    for (int j = 0; j<HEIGHT; j++) {
        for (int i = 0; i<WIDTH; i++) {
            float tone = pixels[i + j*WIDTH]/255;
            if  (tone < 0.14) {
                printf(" ");
            } else if (tone < 0.28) {
                printf(".");
            } else if (tone < 0.42) {
                printf(",");
            } else if (tone < 0.56) {
                printf(":");
            } else if (tone < 0.7) {
                printf(";");
            } else if (tone < 0.84) {
                printf("$");
            } else {
                printf("#");
            }
        }
        printf("\n");
    }
}

void set(float x, float y, float col) {
    
    if ((x >= 0) && (y >= 0) && (x <= WIDTH) && (y < HEIGHT)) {
        pixels[(int) (x + (((int) y) * WIDTH))] = col;
    }
}

void line(float x1, float y1, float x2, float y2) {
    float xmove = (x2 - x1);
    float ymove = (y2 - y1);

    int dist;

    if (abs(xmove) > abs(ymove)) {
        ymove = ymove / abs(xmove);
        dist = (int) xmove;
        xmove = sign(xmove);
    } else {
        xmove = xmove / abs(ymove);
        dist = (int) ymove;
        ymove = sign(ymove);
    }

    float tx = x1;
    float ty = y1;

    for (int i = 0; i<abs(dist); i++) {
        set(tx, ty, 255);
        tx += xmove;
        ty += ymove;
    }
}

struct mat multimat(struct mat one, struct mat two) {
    struct mat out;
    for (int i = 0; i<4; i++) {
        for (int j = 0; j<4; j++) {
            out.m[i][j] = 0;
            for (int k = 0; k<4; k++) {
                out.m[i][j] +=  one.m[k][j] * two.m[i][k];
            }
        }
    }
    return out;
}

struct vector vectormultimat(struct vector vec, struct mat ma) {
    struct vector out = {
        .x = vec.x * ma.m[0][0] + vec.y * ma.m[1][0] + vec.z * ma.m[2][0] + ma.m[3][0],
        .y = vec.y * ma.m[0][1] + vec.y * ma.m[1][1] + vec.z * ma.m[2][1] + ma.m[3][1],
        .z = vec.x * ma.m[0][2] + vec.y * ma.m[1][2] + vec.z * ma.m[2][2] + ma.m[3][2],
        .w = vec.y * ma.m[0][3] + vec.y * ma.m[1][3] + vec.z * ma.m[2][3] + ma.m[3][3],
    };
    return out;
}

struct vector pointmultimat(struct vector vec, struct mat ma) {
    struct vector out = vectormultimat(vec, ma);
    out.x /= out.w;
    out.y /= out.w;
    out.z /= out.w;
    return out;
}

struct tri trianglemultmat(struct tri tr, struct mat ma) {
    struct tri out = {
        .p1 = pointmultimat(tr.p1, ma),
        .p2 = pointmultimat(tr.p2, ma),
        .p3 = pointmultimat(tr.p3, ma),
    };
    return out;
}

struct tri * meshmultmat(struct tri in [], struct tri out [], int len, struct mat ma) {
    for (int i = 0; i<len; i++) {
        out[i] = trianglemultmat(in[i], ma);
    }
    return out;
}

void rotateMesh(struct tri * mesh, int len, float xrot, float yrot, float zrot) {
    struct mat m;
    initOne(&m);
    m.m[1][1] = cosf(xrot);
    m.m[2][2] = cosf(xrot);
    m.m[2][1] = -sinf(xrot);
    m.m[1][2] = sinf(xrot);
    meshmultmat(mesh, mesh, len, m);
    initOne(&m);
    m.m[0][0] = cosf(yrot);
    m.m[2][2] = cosf(yrot);
    m.m[2][0] = sinf(yrot);
    m.m[0][2] = -sinf(yrot);
    meshmultmat(mesh, mesh, len, m);
    initOne(&m);
    m.m[0][0] = cosf(zrot);
    m.m[1][0] = -sinf(zrot);
    m.m[1][1] = cosf(zrot);
    m.m[0][1] = sinf(zrot);
    meshmultmat(mesh, mesh, len, m);
}

void getProjectionMat(struct camera cam, struct mat * out) {
    initZero(out);
    out->m[0][0] = cam.fovRad/(((float) WIDTH) / ((float) HEIGHT));
    out->m[1][1] = cam.fovRad;
    out->m[2][2] = cam.pFar/(cam.pFar - cam.pNear);
    out->m[3][2] = -(cam.pFar * cam.pNear)/(cam.pFar - cam.pNear);
    out->m[2][3] = 1;
}

void triangle(float x1, float y1, float x2, float y2, float x3, float y3) {
    line(x1, y1, x2, y2);
    line(x2, y2, x3, y3);
    line(x3, y3, x1, y1);
}

void wireMesh(struct tri in [], int len, struct camera cam) {
    struct tri mapped [len];
    struct mat proj;
    getProjectionMat(cam, &proj);
    meshmultmat(in, mapped, len, proj);
    for (int i = 0; i<len; i++) {
        triangle(mapped[i].p1.x, mapped[i].p1.y, mapped[i].p2.x, mapped[i].p2.y, mapped[i].p3.x, mapped[i].p3.y);
    }
}

void ftriangle(float x1, float y1, float x2, float y2, float x3, float y3) {

}

void background(float col) {
    for (int i = 0; i<pixelsLength; i++) {
        pixels[i] = col;
    }
}

void clear() {
    printf("\n\nClear\n\n");    
}

void wait(int t) {
    for (int i = 0; i<t*100000; i++) {
        int temp = i*i;
    }
}

void loop(int frameCount) {
    background(0);
    struct camera cam;
    initCamera(&cam);
    struct tri mesh [2];

    struct tri one = {
        .p1 = {
            .x = 0,
            .y = 0,
            .z = 0,
        },
        .p2 = {
            .x = 0,
            .y = 20,
            .z = 0,
        },
        .p3 = {
            .x = 50,
            .y = 0,
            .z = 0,
        },
    };

    mesh[0].p1.x = 4;
    mesh[0].p1.y = 4;
    mesh[0].p1.z = 4;

    mesh[0].p2.x = 50;
    mesh[0].p2.y = 4;
    mesh[0].p2.z = 4;

    mesh[0].p3.x = 4;
    mesh[0].p3.y = 20;
    mesh[0].p3.z = 4;

    float rot = ((float) frameCount) / 10;
    rotateMesh(mesh, 1, 0.0f, 0.0f, rot);
    wireMesh(mesh, 1, cam);
    
    /*line(15, 7, 25, 10);
    triangle(10, 2, 20, 18, 30, 10);*/
}

int main()  { 
    struct mat m;
    initOne(&m);
    int frameCount = 0;
    do {
        clear();
        loop(frameCount);
        draw();
        wait(1000);
        frameCount++;
    } while (1);
    clear();
    return 0;  
} 
