#include "surfaces.h"
#include "vecmat.h"

//Structs
struct _camera {
    vector pos;
    float hrot;
    float vrot;
    float pNear;
    float pFar;
    float fov;
    float fovRad;

    float charRatio;
    float drawDist;
    int width;
    int height;

    surface image;
    surface dists;
} camera;

//functions
void initCamera(camera *cam, int w, int h) {
    cam->hrot = 0;
    cam->vrot = 0;
    cam->pos.x = 0;
    cam->pos.y = 0;
    cam->pos.z = 0;
    cam->fov = 3.141592/4;
    cam->pFar = 1000;
    cam->pNear = 0.1f;
    cam->fovRad = 1/tanf((cam->fov)/2);
    cam->drawDist = 1000;

    cam->width = w;
    cam->height = h;
    cam->charRatio = 2/1;

    initSurf(&(cam->image), w, h);
    initSurf(&(cam->dists), w, h);
}

void clearzs(camera *cam) {
    background(&(cam->dists), cam->drawDist);
}

void getProjectionMat(camera cam, struct mat * out) {
    initMat(out, 0);
    float aspectRatio = ((float) cam.width) / ((float) cam.height);
    out->m[0][0] = cam.fovRad/aspectRatio;
    out->m[1][1] = cam.fovRad;
    out->m[2][2] = cam.pFar/(cam.pFar - cam.pNear);
    out->m[3][2] = -(cam.pFar * cam.pNear)/(cam.pFar - cam.pNear);
    out->m[2][3] = 1;
}

void set3d(camera *cam, float x, float y, float z, colour col) {
    int w = cam->width;
    int h = cam->height;
    if (((x >= 0) && (y >= 0)) && ((x < w) && (y < h))) {
        if (cam->dists.pixels[((int) x) + ((int) y)*w].r > z) {
            cam->dists.pixels[[((int) x) +((int) y)*h].r = z;
            cam->image.pixels[((int) x) + ((int) y)*w] = col;
        }
    }
}

void line3d(camera *cam, float x1, float y1, float z1, float x2, float y2, float z2) {
    float xmove = (x2 - x1);
    float ymove = (y2 - y1);
    float zmove = (z2 - z1);

    int dist;

    if (abs(xmove) > abs(ymove)) {
        ymove = ymove / abs(xmove);
        zmove = zmove / abs(xmove);
        dist = (int) xmove;
        xmove = sign(xmove);
    } else {
        xmove = xmove / abs(ymove);
        zmove = zmove / abs(ymove);
        dist = (int) ymove;
        ymove = sign(ymove);
    }

    float tx = x1;
    float ty = y1;
    float tz = z1;

    for (int i = 0; i<abs(dist); i++) {
        set3d(cam, tx, ty, tz, mColour(255));
        tx += xmove;
        ty += ymove;
        tz += zmove;
    }
}

void triTri(camera *cam, tri t) {
    float charRatio = cam->charRatio;
    float w = cam->width;
    float h = cam->height;
    float x1 = (t.p1.x * w/2) * charRatio + w/2;
    float y1 = t.p1.y * h/2 + h/2;
    float x2 = (t.p2.x * w/2 ) * charRatio + w/2;
    float y2 = t.p2.y * h/2 + h/2;
    float x3 = (t.p3.x * w/2) * charRatio + w/2;
    float y3 = t.p3.y * h/2 + h/2;
    line3d(cam, x1, y1, t.p1.z, x2, y2, t.p2.z);
    line3d(cam, x2, y2, t.p2.z, x3, y3, t.p3.z);
    line3d(cam, x3, y3, t.p3.z, x1, y1, t.p1.z);
}