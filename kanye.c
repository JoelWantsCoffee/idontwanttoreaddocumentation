#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

//Global Variables
const int WIDTH = 100;
const int HEIGHT = 50;
const int pixelsLength = 100*50;
float pixels[100*50];
float zs[100*50];
const char * cols = " .,:;$#";
const float charRatio = 2/1;
float DRAWDIST = 1000;


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

void clearzs() {
    for (int i = 0; i<WIDTH*HEIGHT; i++) {
	zs[i] = DRAWDIST;
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
            float tone = ((float) pixels[i + j*WIDTH])/255;	
	    if  (tone < 0.14) {
                printf("%c", cols[0]);
            } else if (tone < 0.28) {
                printf("%c", cols[1]);
            } else if (tone < 0.42) {
                printf("%c", cols[2]);
            } else if (tone < 0.56) {
                printf("%c", cols[3]);
            } else if (tone < 0.7) {
                printf("%c", cols[4]);
            } else if (tone < 0.84) {
                printf("%c", cols[5]);
            } else {
                printf("%c", cols[6]);
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

void set3d(float x, float y, float z, float col) {
    if (((x >= 0) && (y >= 0)) && ((x <= WIDTH) && (y < HEIGHT))) {
	if (zs[((int) x) + ((int) y)*WIDTH] > z) {
	    zs[((int) x) +((int) y)*WIDTH] = z;
	    pixels[((int) x) + ((int) y)*WIDTH] = col;
	}
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
    triangle((t.p1.x * WIDTH/2) * charRatio + WIDTH/2, t.p1.y * HEIGHT/2 + HEIGHT/2, (t.p2.x * WIDTH/2 ) * charRatio + WIDTH/2, t.p2.y * HEIGHT/2 + HEIGHT/2, (t.p3.x * WIDTH/2) * charRatio + WIDTH/2, t.p3.y * HEIGHT/2 + HEIGHT/2);
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

void fillTriangle3d(float x1, float y1, float x2, float y2, float x3, float y3, float z, float col) {
  struct vector vecs [3];

  vecs[0].x = x1;
  vecs[0].y = y1;
  vecs[1].x = x2;
  vecs[1].y = y2;
  vecs[2].x = x3;
  vecs[2].y = y3;

  int top = 0;
  int bottom = 0;
  for (int i = 0; i<3;i++) {
    if (vecs[i].y > vecs[top].y) top = i;
    if (vecs[i].y < vecs[bottom].y) bottom = i;
  }
  
  int middle = 0;
  
  for (int i = 0; i<3; i++) {
    if ((bottom != i) && (top != i)) middle = i;
  }
  
  struct vector v1 = vecSubtract(vecs[bottom], vecs[middle]);
  struct vector v2 = vecSubtract(vecs[bottom], vecs[top]);
  
  float v1a = v1.x/v1.y;
  float v2a = v2.x/v2.y;
  
  struct vector np = {
	  .x = (vecs[middle].y-vecs[bottom].y)*v2a + vecs[bottom].x,
	  .y =  vecs[middle].y,
  };
  
  for (int i = 0; i<=vecs[middle].y-vecs[bottom].y; i++) {
    float v1x = v1a*i;
    float v2x = v2a*i;
    int dir;
    
    
    //point(floor(vecs[bottom].x + v1x) , floor(i + vecs[bottom].y));
    //point(floor(vecs[bottom].x + v2x) , floor(i + vecs[bottom].y));
    
    int x1 = round(vecs[bottom].x + v1x);
    int x2 = round(vecs[bottom].x + v2x);
    
    int y = floor(i + vecs[bottom].y);
    
    
    if (x1 > x2) {
      dir = -1;
    } else {
      dir = 1;
    }
    
    
    
    for (int j = x1; j != x2+dir; j+=dir) {
      set3d(j, y, z, col);
    }
  }
  
  struct vector vv1 = vecSubtract(vecs[middle], vecs[top]);
  struct vector vv2 = vecSubtract(np, vecs[top]);
  
  float vv1a = vv1.x/vv1.y;  
  float vv2a = vv2.x/vv2.y;
  
  for (int i = 0; i<(vecs[top].y - vecs[middle].y); i++) {
    float v1x = vv1a*i;
    float v2x = vv2a*i;
    int dir;
    
    int x1 = round(vecs[middle].x + v1x);
    int x2 = round(np.x + v2x);
    
    int y = floor(i + vecs[middle].y);
    
    if (x1 > x2) {
      dir = -1;
    } else {
      dir = 1;
    }
    
    
    for (int j = x1; j != x2+dir; j+=dir) {
      set3d(j, y, z, col);
    }
  }
}

void fillTri(struct tri t, float col) {
    float z = (t.p1.z + t.p2.z + t.p3.z)/3;
    fillTriangle3d((t.p1.x * WIDTH/2) * charRatio + WIDTH/2, t.p1.y * HEIGHT/2 + HEIGHT/2, (t.p2.x * WIDTH/2 ) * charRatio + WIDTH/2, t.p2.y * HEIGHT/2 + HEIGHT/2, (t.p3.x * WIDTH/2) * charRatio + WIDTH/2, t.p3.y * HEIGHT/2 + HEIGHT/2, z, col);
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

void fillMesh(struct tri in [], int len, struct camera cam, struct vector light) {
    struct tri mapped [len];
    struct mat proj;
    getProjectionMat(cam, &proj);
    meshmultmat(in, mapped, len, proj);
    for (int i = 0; i<len; i++) {
	float col = (dotProduct(getNormal(in[i]), light) + 1)*0.5*255;
        if (dotProduct(vecSubtract(in[i].p1,cam.pos), getNormal(in[i])) < 0) fillTri(mapped[i], col);
    }
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
    clearzs();
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
    

    struct vector light = {
	.x = 0.6f,
	.y = 0,
	.z = -0.4,
    };

    float rot = (((float) frameCount) / 10);

    rotateMesh(mesh, 12, rot, rot, 0.0f);

    struct tri transd [12];

    struct mat translate;

    initOne(&translate);

    translate.m[3][0] = 0.0f;
    translate.m[3][1] = 0.0f;
    translate.m[3][2] = 5.0f + 2*(sinf(rot) + 1);

    meshmultmat(mesh, transd, 12, translate);

    fillMesh(transd, 12, cam, light);

    wireMesh(transd, 12, cam);

    

    /*line(15, 7, 25, 10);
    triangle(10, 2, 20, 18, 30, 10);*/
}

void run(void (*f)(int), int frameCount, float frameRate) {
    clock_t mt;
    mt = clock();
    clock_t t;
    t = clock();
    f(frameCount);
    while ( ((double)t - (double)mt)/CLOCKS_PER_SEC < 1/frameRate ) {
	t = clock();
    }
}

void engine(int frameCount) {
    //clear();
    loop(frameCount);
    draw();
}

int main()  {
    
    int frameCount = 0;
    do {
	run(engine, frameCount, 30);
        frameCount++;
    } while (1);
    clear();
    return 0;  
} 
