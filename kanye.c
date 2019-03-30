#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

//Global Variables
const int WIDTH = 50;
const int HEIGHT = 50;
const int pixelsLength = 50*50;
float pixels[50*50];
const char * cols = " .,:;$#";

struct mat {
    float m[4][4];
};

struct vector {
    float x;
    float y;
    float z;
    float w;
};

struct tri {
    struct vector p1;
    struct vector p2;
    struct vector p3;
};

struct camera {
    struct vector pos;
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
    cam->pos.x = 0;
    cam->pos.y = 0;
    cam->pos.z = 0;
    cam->fov = 3.141592/4;
    cam->pFar = 1000;
    cam->pNear = 0.1f;
    cam->fovRad = 1/tanf((cam->fov)/2);
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
        .p2 = pointmultimat	(tr.p2, ma),
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

	//fopen filename

}

void getProjectionMat(struct camera cam, struct mat * out) {
    initZero(out);
    float aspectRatio = ((float) WIDTH) / ((float) HEIGHT);
    out->m[0][0] = cam.fovRad/aspectRatio;
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

void triTri(struct tri t) {
    triangle(t.p1.x * WIDTH/2 + WIDTH/2, t.p1.y * HEIGHT/2 + HEIGHT/2, t.p2.x * WIDTH/2 + WIDTH/2, t.p2.y * HEIGHT/2 + HEIGHT/2, t.p3.x * WIDTH/2 + WIDTH/2, t.p3.y * HEIGHT/2 + HEIGHT/2);
}

struct vector crossProduct(struct vector v1, struct vector v2) {
    struct vector out;

    out.x = v1.y*v2.z - v1.z*v2.y;
    out.y = v1.z*v2.x - v1.x*v2.z;
    out.z = v1.x*v2.y - v1.y*v2.x;

    return out;
}

float dotProduct(struct vector v1, struct vector v2) {
    float out =	v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    return out;
}

struct vector vecSubtract(struct vector v1, struct vector v2) {
    struct vector out = {
	.x = v1.x - v2.x,
	.y = v1.y - v2.y,
	.z = v1.z - v2.z,
    };
    return out;
}

struct vector getNormal(struct tri in) {
    struct vector v1 = vecSubtract(in.p1, in.p2);
    struct vector v2 = vecSubtract(in.p1, in.p3);
    return crossProduct(v1, v2);
}

void wireMesh(struct tri in [], int len, struct camera cam) {
    struct tri mapped [len];
    struct mat proj;
    getProjectionMat(cam, &proj);
    meshmultmat(in, mapped, len, proj);
    for (int i = 0; i<len; i++) {
        if (dotProduct(vecSubtract(in[i].p1,cam.pos), getNormal(in[i])) < 0) triTri(mapped[i]);
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

struct tri mesh [12];
mesh[0].p1.x = 1.0;
mesh[0].p1.y = -1.0;
mesh[0].p1.z = 1.0;
mesh[0].p2.x = -1.0;
mesh[0].p2.y = -1.0;
mesh[0].p2.z = -1.0;
mesh[0].p3.x = 1.0;
mesh[0].p3.y = -1.0;
mesh[0].p3.z = -1.0;
mesh[1].p1.x = -1.0;
mesh[1].p1.y = 1.0;
mesh[1].p1.z = -1.0;
mesh[1].p2.x = 0.999999;
mesh[1].p2.y = 1.0;
mesh[1].p2.z = 1.000001;
mesh[1].p3.x = 1.0;
mesh[1].p3.y = 1.0;
mesh[1].p3.z = -0.999999;
mesh[2].p1.x = 1.0;
mesh[2].p1.y = 1.0;
mesh[2].p1.z = -0.999999;
mesh[2].p2.x = 1.0;
mesh[2].p2.y = -1.0;
mesh[2].p2.z = 1.0;
mesh[2].p3.x = 1.0;
mesh[2].p3.y = -1.0;
mesh[2].p3.z = -1.0;
mesh[3].p1.x = 0.999999;
mesh[3].p1.y = 1.0;
mesh[3].p1.z = 1.000001;
mesh[3].p2.x = -1.0;
mesh[3].p2.y = -1.0;
mesh[3].p2.z = 1.0;
mesh[3].p3.x = 1.0;
mesh[3].p3.y = -1.0;
mesh[3].p3.z = 1.0;
mesh[4].p1.x = -1.0;
mesh[4].p1.y = -1.0;
mesh[4].p1.z = 1.0;
mesh[4].p2.x = -1.0;
mesh[4].p2.y = 1.0;
mesh[4].p2.z = -1.0;
mesh[4].p3.x = -1.0;
mesh[4].p3.y = -1.0;
mesh[4].p3.z = -1.0;
mesh[5].p1.x = 1.0;
mesh[5].p1.y = -1.0;
mesh[5].p1.z = -1.0;
mesh[5].p2.x = -1.0;
mesh[5].p2.y = 1.0;
mesh[5].p2.z = -1.0;
mesh[5].p3.x = 1.0;
mesh[5].p3.y = 1.0;
mesh[5].p3.z = -0.999999;
mesh[6].p1.x = 1.0;
mesh[6].p1.y = -1.0;
mesh[6].p1.z = 1.0;
mesh[6].p2.x = -1.0;
mesh[6].p2.y = -1.0;
mesh[6].p2.z = 1.0;
mesh[6].p3.x = -1.0;
mesh[6].p3.y = -1.0;
mesh[6].p3.z = -1.0;
mesh[7].p1.x = -1.0;
mesh[7].p1.y = 1.0;
mesh[7].p1.z = -1.0;
mesh[7].p2.x = -1.0;
mesh[7].p2.y = 1.0;
mesh[7].p2.z = 1.0;
mesh[7].p3.x = 0.999999;
mesh[7].p3.y = 1.0;
mesh[7].p3.z = 1.000001;
mesh[8].p1.x = 1.0;
mesh[8].p1.y = 1.0;
mesh[8].p1.z = -0.999999;
mesh[8].p2.x = 0.999999;
mesh[8].p2.y = 1.0;
mesh[8].p2.z = 1.000001;
mesh[8].p3.x = 1.0;
mesh[8].p3.y = -1.0;
mesh[8].p3.z = 1.0;
mesh[9].p1.x = 0.999999;
mesh[9].p1.y = 1.0;
mesh[9].p1.z = 1.000001;
mesh[9].p2.x = -1.0;
mesh[9].p2.y = 1.0;
mesh[9].p2.z = 1.0;
mesh[9].p3.x = -1.0;
mesh[9].p3.y = -1.0;
mesh[9].p3.z = 1.0;
mesh[10].p1.x = -1.0;
mesh[10].p1.y = -1.0;
mesh[10].p1.z = 1.0;
mesh[10].p2.x = -1.0;
mesh[10].p2.y = 1.0;
mesh[10].p2.z = 1.0;
mesh[10].p3.x = -1.0;
mesh[10].p3.y = 1.0;
mesh[10].p3.z = -1.0;
mesh[11].p1.x = 1.0;
mesh[11].p1.y = -1.0;
mesh[11].p1.z = -1.0;
mesh[11].p2.x = -1.0;
mesh[11].p2.y = -1.0;
mesh[11].p2.z = -1.0;
mesh[11].p3.x = -1.0;
mesh[11].p3.y = 1.0;
mesh[11].p3.z = -1.0;
    

    float rot = (((float) frameCount) / 10);
	printf("%f\n", rot);

    rotateMesh(mesh, 12, rot, rot, 0.0f);

    struct tri transd [12];

    struct mat translate;

    initOne(&translate);

    translate.m[3][0] = 0.0f;
    translate.m[3][1] = 0.0f;
    translate.m[3][2] = 5.0f + 2*(sinf(rot) + 1);

    meshmultmat(mesh, transd, 12, translate);

    wireMesh(transd, 12, cam);

    

    /*line(15, 7, 25, 10);
    triangle(10, 2, 20, 18, 30, 10);*/
}

int main()  { 
    struct mat m;
    initOne(&m);
    int frameCount = 0;
    do {
        clear();
        loop(frameCount % 62);
        draw();
        wait(100);
        frameCount++;
    } while (1);
    clear();
    return 0;  
} 
