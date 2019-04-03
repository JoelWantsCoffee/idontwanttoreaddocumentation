#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <locale.h>

//Global Variables
const int WIDTH = 299;
const int HEIGHT = 90;
const int pixelsLength = 299*90;
float pixels[299*90];
char pcols[299*90];
float zs[299*90];
const char * cols = " .,>x#@";//" ▏▎▍▌▋▊▉█"; // " ░▒▓█";  
const int colsLen = 7;
const float charRatio = 2/1;
float DRAWDIST = 1000;
/*
float ROTX;
float ROTY;
float ROTZ;
float TRANX;
float TRANY;
float TRANZ;
float LX;
float LY;
float LZ;
*/

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

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

float sign(float in) {
    if (in >= 0) {
        return 1;
    } else {
        return -1;
    }
}

int ffloor(float in) {
    if (in >= 0) {
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
            float tone = ((float) pixels[i + j*WIDTH])/256;	
	        printf("%c",cols[ffloor(tone*((float) colsLen))]);
	    }
        printf("\n");
    }
    fflush(stdout);
}

void set(float x, float y, float col) {
    if ((x >= 0) && (y >= 0) && (x < WIDTH) && (y < HEIGHT)) {
        pixels[(int) (x + (((int) y) * WIDTH))] = col;
    }
}

void set3d(float x, float y, float z, float col) {
    if (((x >= 0) && (y >= 0)) && ((x < WIDTH) && (y < HEIGHT))) {
	if (zs[((int) x) + ((int) y)*WIDTH] > z) {
	    zs[((int) x) +((int) y)*WIDTH] = z;
	    pixels[((int) x) + ((int) y)*WIDTH] = col;
	}
    }
}

void line3d(float x1, float y1, float z1, float x2, float y2, float z2) {
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
        set3d(tx, ty, tz, 255);
        tx += xmove;
        ty += ymove;
        tz += zmove;
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
    float x1 = (t.p1.x * WIDTH/2) * charRatio + WIDTH/2;
    float y1 = t.p1.y * HEIGHT/2 + HEIGHT/2;
    float x2 = (t.p2.x * WIDTH/2 ) * charRatio + WIDTH/2;
    float y2 = t.p2.y * HEIGHT/2 + HEIGHT/2;
    float x3 = (t.p3.x * WIDTH/2) * charRatio + WIDTH/2;
    float y3 = t.p3.y * HEIGHT/2 + HEIGHT/2;
    line3d(x1, y1, t.p1.z, x2, y2, t.p2.z);
    line3d(x2, y2, t.p2.z, x3, y3, t.p3.z);
    line3d(x3, y3, t.p3.z, x1, y1, t.p1.z);
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

struct vector vecAdd(struct vector v1, struct vector v2) {
    struct vector out = {
	.x = v1.x + v2.x,
	.y = v1.y + v2.y,
	.z = v1.z + v2.z,
    };
    return out;
}

struct vector scale(struct vector v, float s) {
    struct vector out = {
	.x = v.x * s,
	.y = v.y * s,
	.z = v.z * s,	
    };
    return out;
}

void vecNormalise(struct vector * v) {
    float largestValue = fmax(abs(v->x), abs(v->y));
    largestValue = fmax(largestValue, abs(v->z));
    v->x /= largestValue;
    v->y /= largestValue;
    v->z /= largestValue;
}

struct vector vecIntersectPlane(struct vector p1, struct vector p2, struct vector *l1, struct vector *l2) {
	vecNormalise(&p2);
	float pd = -dotProduct(p2, p1);
	float ad = dotProduct(*l1, p2);
	float bd = dotProduct(*l2, p2);
	float t = (-pd - ad) / (bd - ad);
	struct vector wl = vecSubtract(*l2, *l1);
	struct vector pl = scale(wl, t);
	return vecAdd(*l1, pl);
}

float distance(struct vector p, struct vector n, struct vector *v) {
    return (n.x)*(v->x) + (n.y)*(v->y) + (n.z)*(v->z) - dotProduct(n, p);
}

int triClip(struct vector pl, struct vector pn, struct tri in, struct tri *out1, struct tri *out2) {
    vecNormalise(&pn);
    
    struct vector *ipts [3];
    struct vector *opts [3];
    int ic = 0;
    int oc = 0;

    float d0 = distance(pl, pn, &in.p1);
    float d1 = distance(pl, pn, &in.p2);
    float d2 = distance(pl, pn, &in.p3);

    if (d0 >= 0) {
	ipts[ic++] = &in.p1;
    } else {
	opts[oc++] = &in.p1;
    }

    if (d1 >= 0) {
	ipts[ic++] = &in.p2;
    } else {
	opts[oc++] = &in.p2;
    }
    
    if (d2 >= 0) {
	ipts[ic++] = &in.p3;
    } else {
	opts[oc++] = &in.p3;
    }


    if (ic == 0) {
	return 0;
    } else if (ic == 3) {
	out1->p1 = in.p1;
	out1->p2 = in.p2;
	out1->p3 = in.p3;
	return 1;
    } else if (ic == 2) {
	out1->p1 = *ipts[0];
	out1->p2 = vecIntersectPlane(pl,pn,ipts[0],opts[0]);
	out1->p3 = *ipts[1];

	out2->p1 = *ipts[1];
	out2->p2 = out1->p2;
	out2->p3 = vecIntersectPlane(pl,pn,ipts[1],opts[0]);
	return 2;
    } else if (ic == 1) {
	out1->p1 = *ipts[0];
	out1->p2 = vecIntersectPlane(pl,pn,ipts[0],opts[0]);
	out1->p3 = vecIntersectPlane(pl,pn,ipts[0],opts[1]);
	return 1;
    }
}

void fillTriangle3d(float x1, float y1, float z1, float x2, float y2, float z2,float x3, float y3, float z3, float col) {
  struct vector vecs [3];
  int con = 1;
  /*if (!ffloor(x1 - x2) && !ffloor(y1 - y2)) con = 0;
  if (!ffloor(x2 - x3) && !ffloor(y2 - y3)) con = 0;
  if (!ffloor(x3 - x1) && !floor(y3 - y1)) con = 0;*/
  if (con) {
  vecs[0].x = x1;
  vecs[0].y = y1;
  vecs[0].z = z1;
  vecs[1].x = x2;
  vecs[1].y = y2;
  vecs[1].z = z2;
  vecs[2].x = x3;
  vecs[2].y = y3;
  vecs[2].z = z3;
/*
  struct vector *minni;
  struct vector *middi;
  struct vector *maxxi;

  if (vecs[0].y < vecs[1].y) {
	if (vecs[0].y < vecs[2].y) {
	    minni = &vecs[0];
	    if (vecs[2].y < vecs[1].y) {
		middi = &vecs[2];
		maxxi = &vecs[1];
	    } else {
		maxxi = &vecs[2];
		middi = &vecs[1];
	    }
	} else {
	    minni = &vecs[2];
	    middi = &vecs[0];
	    maxxi = &vecs[1];
	}
  } else {
	if (vecs[1].y < vecs[2].y) {
	    minni = &vecs[1];
	    if (vecs[0].y < vecs[2].y) {
		middi = &vecs[0];
		maxxi = &vecs[2];
	    } else {
		maxxi = &vecs[0];
		middi = &vecs[2];
	    }
	} else {
	    minni = &vecs[2];
	    middi = &vecs[1];
	    maxxi = &vecs[0];
	}
  }
  for (int i = ceil(minni->y); i<ceil(maxxi->y); i++) {
	if (i < middi->y) {
		struct vector tmax = scale(vecSubtract(*maxxi,*minni),  (i - minni->y)/(maxxi->y - minni->y));
		struct vector tmid = scale(vecSubtract(*middi,*minni),  (i - minni->y)/(middi->y - minni->y));
		int xbase = tmid.x;
		int xend = tmax.x;
		for (int j = xbase; j < xend; j++) {
			struct vector zvec = vecSubtract(tmax, tmid);
			float z = tmid.z + zvec.z * (j-xbase)/(xend-xbase);
			set3d(j+minni->x, i, z, col);
		}
	} else {
		struct vector tmin = scale(vecSubtract(*minni, *maxxi), (maxxi->y - minni->y)/(i - minni->y));
	        struct vector tmid = scale(vecSubtract(*middi, *maxxi), (maxxi->y - middi->y)/(i - middi->y));	
		int xbase = tmid.x;
		int xend = tmin.x;
		for (int j = xbase; j <= xend; j++) {
			struct vector zvec = vecSubtract(tmin, tmid);
			float z = tmid.z + zvec.z * (j-xbase)/(xend-xbase);
			set3d(j+minni->x, i, z, col);
		}
	}
  }*/
  }
  if (con) {
  int top = 0;
  int bottom = 0;
  for (int i = 0; i<3;i++) {
    if (vecs[i].y > vecs[top].y) top = i;
    if (vecs[i].y < vecs[bottom].y) bottom = i;
  }
  
  float z = MAX(MAX(z1, z2), z3);

  int middle = 0;
  
  for (int i = 0; i<3; i++) {
    if ((bottom != i) && (top != i)) middle = i;
  }
  
  struct vector v1 = vecSubtract(vecs[bottom], vecs[middle]);
  struct vector v2 = vecSubtract(vecs[bottom], vecs[top]);
  
  if (!floor(v1.y)) v1.y = 1;
  if (!floor(v2.y)) v2.y = 1;

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
    for (int j = MIN(x1, x2); j<=MAX(x1, x2); j++) {
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
     
    for (int j = MIN(x1, x2); j<=MAX(x1, x2); j++) {
	set3d(j, y, z, col);
    }
  }
  }
}

void fillTri(struct tri t, float col) {
    if (!((t.p1.z >= DRAWDIST) || (t.p2.z >= DRAWDIST) || (t.p3.z >= DRAWDIST))) {
        fillTriangle3d((t.p1.x * WIDTH/2) * charRatio + WIDTH/2, t.p1.y * HEIGHT/2 + HEIGHT/2, t.p1.z, (t.p2.x * WIDTH/2 ) * charRatio + WIDTH/2, t.p2.y * HEIGHT/2 + HEIGHT/2, t.p2.z, (t.p3.x * WIDTH/2) * charRatio + WIDTH/2, t.p3.y * HEIGHT/2 + HEIGHT/2, t.p3.z, col);
    }
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
        if (dotProduct(vecSubtract(in[i].p1,cam.pos), getNormal(in[i])) > 0) triTri(mapped[i]);
    }
}

float magnitude(struct vector in) {
    return sqrt(in.x*in.x + in.y*in.y + in.z*in.z);
}

float sudoTheta(struct vector v1, struct vector v2) {
    return (dotProduct(v1, v2)/(magnitude(v1), magnitude(v2)));
}

void fillMesh(struct tri in [], int len, struct camera cam, struct vector light) {
    struct tri mapped [len];
    struct mat proj;
    getProjectionMat(cam, &proj);
    meshmultmat(in, mapped, len, proj);
    struct vector pl = {
    	.x = 0,
	.y = 0,
	.z = cam.pNear,
    };
    struct vector pn = {
    	.x = 0,
	.y = 0,
	.z = 1,
    };
    for (int i = 0; i<len; i++) {
	float col = (sudoTheta(getNormal(in[i]), light) + 1)*0.5f*255;
        if (dotProduct(vecSubtract(in[i].p1,cam.pos), getNormal(in[i])) < 0) {
	    struct tri out [2];    
	    int ntris = triClip(pl, pn, in[i], &out[0], &out[1]);
	    for (int j = 0; j<ntris; j++) {
	    	struct tri m = trianglemultmat(out[j], proj);
		for (int k = 0; k < ntris; k++) {
		    struct vector pl = {.x = -1, .y = 0, .z = 0,};
		    struct vector pn = {.x = 1,	.y = 0,	.z = 0,};
		    struct tri nout [2];
		    int nntris = triClip(pl, pn, m, &nout[0], &nout[1]);
		    for (int l = 0; l<nntris; l++) {
			struct vector pl = {.x = 1, .y = 0, .z = 0,};
			struct vector pn = {.x = -1, .y = 0, .z = 0,};
			struct tri nnout [2];
			int nnntris = triClip(pl, pn, nout[l], &nnout[0], &nnout[1]);
			for (int m = 0; m<nnntris; m++) {
			    struct vector pl = {.x = 0, .y = -1, .z = 0,};
			    struct vector pn = {.x = 0, .y = 1, .z = 0,};
			    struct tri nout [2];
			    int nntris = triClip(pl, pn, nnout[l], &nout[0], &nout[1]);
			    for (int n = 0; n<nntris; n++) {
				struct vector pl = {.x = 0, .y = 1, .z = 0,};
				struct vector pn = {.x = 0, .y = -1, .z = 0,};
				struct tri nnout [2];
				int nnntris = triClip(pl, pn, nout[l], &nnout[0], &nnout[1]);
				for (int o = 0; o<nnntris; o++) {
				    fillTri(nnout[o], col);
				}
			    }
			}
		    }
		}
	    }
	}
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

    int objlen = 311;

    struct vector light = {
	    .x = 0,
	    .y = -1,
	    .z = -0.1,
    };

    float rot = (((float) frameCount) / 10);

    rotateMesh(mesh, objlen, 3.1415f, rot, 0);

    struct tri transd [objlen];

    struct mat translate;

    initOne(&translate);

    translate.m[3][0] = 0;
    translate.m[3][1] = 2;
    translate.m[3][2] = 10 + 4*(sinf(rot) + 0);


    meshmultmat(mesh, transd, objlen, translate);

    fillMesh(transd, objlen, cam, light);

    wireMesh(transd, objlen, cam);

    /*set(3, 2, 255);
    line(15, 7, 25, 10);
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
    //printf("\033[47:47m");    
    int frameCount = 0;
    do {
	run(engine, frameCount % 62, 30);
        frameCount++;
    } while (1);
    clear();
    return 0;  
} 

