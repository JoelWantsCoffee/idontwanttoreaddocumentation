#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <locale.h>
#include <pthread.h>
#include <termios.h>
#include <X11/Xlib.h>

//Global Variables
const int WIDTH = 1234;
const int HEIGHT = 160;
const int pixelsLength = 1234*160;
float pixels[1234*160];
char pcols[1234*160];
float zs[1234*160];
const char * cols = " .,>x#@";//" ▏▎▍▌▋▊▉█"; // " ░▒▓█";  
const int colsLen = 7;
const float charRatio = 2/1;
float DRAWDIST = 1000;

const char * keys = "wasd";
int keyDown [8];
const int keyDownLen = 8;


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
    float hrot;
    float vrot;
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
    cam->hrot = 0;
    cam->vrot = 0;
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
	        putchar(cols[ffloor(tone*((float) colsLen))]);
	    }
        putchar('\n');
    }
    fflush(stdout);
}
/*
void drawXLib(Display disp) {
    for (int j = 0; j<HEIGHT; j++) {
        for (int i = 0; i<WIDTH; i++) {
            XColor xcolour;
            xcolour.red = pixels[i + j*WIDTH]);	xcolour.green = pixels[i + j*WIDTH]); xcolour.green = pixels[i + j*WIDTH]);	
            xcolour.flags = DoRed | DoGreen | DoBlue;
            XDrawPoint(disp, NULL, NULL, i,j);
	        
	    }
        XFlush(disp);
    }
    XDrawPoint(disp, )
}*/

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
    struct mat translate;
    initOne(&translate);
    translate.m[3][0] = cam.pos.x;
    translate.m[3][1] = cam.pos.y;
    translate.m[3][2] = -cam.pos.z;
    meshmultmat(in, mapped, len, translate);
    meshmultmat(mapped, mapped, len, proj);
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
    struct mat proj;
    getProjectionMat(cam, &proj);
    struct mat translate;
    initOne(&translate);
    translate.m[3][0] = cam.pos.x;
    translate.m[3][1] = cam.pos.y;
    translate.m[3][2] = -cam.pos.z;
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
        if (dotProduct(vecSubtract(cam.pos, in[i].p1), getNormal(in[i])) > 0) {
	    struct tri out [2];    
	    int ntris = triClip(pl, pn, trianglemultmat(in[i], translate), &out[0], &out[1]);
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

void loop(int frameCount, struct camera * cam) {
    background(0);
    clearzs();
    struct tri mesh [370];
    {
        mesh[0].p1.x = 0.007242;
        mesh[0].p1.y = -0.264548;
        mesh[0].p1.z = 0.31558;
        mesh[0].p2.x = 0.023642;
        mesh[0].p2.y = -0.252005;
        mesh[0].p2.z = 0.316736;
        mesh[0].p3.x = 0.007313;
        mesh[0].p3.y = -0.265396;
        mesh[0].p3.z = 0.327658;
        mesh[1].p1.x = 0.002083;
        mesh[1].p1.y = -0.267639;
        mesh[1].p1.z = 0.312103;
        mesh[1].p2.x = 0.024291;
        mesh[1].p2.y = -0.263459;
        mesh[1].p2.z = 0.302255;
        mesh[1].p3.x = 0.007242;
        mesh[1].p3.y = -0.264548;
        mesh[1].p3.z = 0.31558;
        mesh[2].p1.x = 0.011163;
        mesh[2].p1.y = -0.226997;
        mesh[2].p1.z = 0.334062;
        mesh[2].p2.x = 0.029393;
        mesh[2].p2.y = -0.219512;
        mesh[2].p2.z = 0.325764;
        mesh[2].p3.x = 0.01717;
        mesh[2].p3.y = -0.220033;
        mesh[2].p3.z = 0.328763;
        mesh[3].p1.x = 0.023642;
        mesh[3].p1.y = -0.252005;
        mesh[3].p1.z = 0.316736;
        mesh[3].p2.x = 0.029393;
        mesh[3].p2.y = -0.219512;
        mesh[3].p2.z = 0.325764;
        mesh[3].p3.x = 0.011163;
        mesh[3].p3.y = -0.226997;
        mesh[3].p3.z = 0.334062;
        mesh[4].p1.x = 0.031688;
        mesh[4].p1.y = -0.226472;
        mesh[4].p1.z = 0.28078;
        mesh[4].p2.x = 0.037105;
        mesh[4].p2.y = -0.210377;
        mesh[4].p2.z = 0.330419;
        mesh[4].p3.x = 0.029393;
        mesh[4].p3.y = -0.219512;
        mesh[4].p3.z = 0.325764;
        mesh[5].p1.x = 0.031688;
        mesh[5].p1.y = -0.226472;
        mesh[5].p1.z = 0.28078;
        mesh[5].p2.x = 0.029393;
        mesh[5].p2.y = -0.219512;
        mesh[5].p2.z = 0.325764;
        mesh[5].p3.x = 0.023642;
        mesh[5].p3.y = -0.252005;
        mesh[5].p3.z = 0.316736;
        mesh[6].p1.x = 0.007313;
        mesh[6].p1.y = -0.265396;
        mesh[6].p1.z = 0.327658;
        mesh[6].p2.x = 0.023642;
        mesh[6].p2.y = -0.252005;
        mesh[6].p2.z = 0.316736;
        mesh[6].p3.x = 0.011163;
        mesh[6].p3.y = -0.226997;
        mesh[6].p3.z = 0.334062;
        mesh[7].p1.x = 0.048469;
        mesh[7].p1.y = -0.158545;
        mesh[7].p1.z = 0.319398;
        mesh[7].p2.x = 0.037105;
        mesh[7].p2.y = -0.210377;
        mesh[7].p2.z = 0.330419;
        mesh[7].p3.x = 0.038391;
        mesh[7].p3.y = -0.156011;
        mesh[7].p3.z = 0.280595;
        mesh[8].p1.x = 0.031688;
        mesh[8].p1.y = -0.226472;
        mesh[8].p1.z = 0.28078;
        mesh[8].p2.x = 0.038391;
        mesh[8].p2.y = -0.156011;
        mesh[8].p2.z = 0.280595;
        mesh[8].p3.x = 0.037105;
        mesh[8].p3.y = -0.210377;
        mesh[8].p3.z = 0.330419;
        mesh[9].p1.x = 0.035308;
        mesh[9].p1.y = -0.137464;
        mesh[9].p1.z = 0.306814;
        mesh[9].p2.x = 0.048469;
        mesh[9].p2.y = -0.158545;
        mesh[9].p2.z = 0.319398;
        mesh[9].p3.x = 0.038391;
        mesh[9].p3.y = -0.156011;
        mesh[9].p3.z = 0.280595;
        mesh[10].p1.x = 0.048469;
        mesh[10].p1.y = -0.158545;
        mesh[10].p1.z = 0.319398;
        mesh[10].p2.x = 0.035308;
        mesh[10].p2.y = -0.137464;
        mesh[10].p2.z = 0.306814;
        mesh[10].p3.x = 0.053788;
        mesh[10].p3.y = -0.129921;
        mesh[10].p3.z = 0.303985;
        mesh[11].p1.x = 0.049253;
        mesh[11].p1.y = -0.135547;
        mesh[11].p1.z = 0.325422;
        mesh[11].p2.x = 0.039196;
        mesh[11].p2.y = -0.156989;
        mesh[11].p2.z = 0.339227;
        mesh[11].p3.x = 0.066431;
        mesh[11].p3.y = -0.16482;
        mesh[11].p3.z = 0.292051;
        mesh[12].p1.x = 0.039196;
        mesh[12].p1.y = -0.156989;
        mesh[12].p1.z = 0.339227;
        mesh[12].p2.x = -7.79E-4;
        mesh[12].p2.y = -0.146982;
        mesh[12].p2.z = 0.344462;
        mesh[12].p3.x = 0.0;
        mesh[12].p3.y = -0.195958;
        mesh[12].p3.z = 0.354306;
        mesh[13].p1.x = 0.043042;
        mesh[13].p1.y = -0.1738;
        mesh[13].p1.z = 0.335948;
        mesh[13].p2.x = 0.039196;
        mesh[13].p2.y = -0.156989;
        mesh[13].p2.z = 0.339227;
        mesh[13].p3.x = 0.0;
        mesh[13].p3.y = -0.195958;
        mesh[13].p3.z = 0.354306;
        mesh[14].p1.x = 0.037105;
        mesh[14].p1.y = -0.210377;
        mesh[14].p1.z = 0.330419;
        mesh[14].p2.x = 0.024632;
        mesh[14].p2.y = -0.221933;
        mesh[14].p2.z = 0.335448;
        mesh[14].p3.x = 0.029393;
        mesh[14].p3.y = -0.219512;
        mesh[14].p3.z = 0.325764;
        mesh[15].p1.x = 0.037105;
        mesh[15].p1.y = -0.210377;
        mesh[15].p1.z = 0.330419;
        mesh[15].p2.x = 0.011163;
        mesh[15].p2.y = -0.226997;
        mesh[15].p2.z = 0.334062;
        mesh[15].p3.x = 0.024632;
        mesh[15].p3.y = -0.221933;
        mesh[15].p3.z = 0.335448;
        mesh[16].p1.x = 0.024632;
        mesh[16].p1.y = -0.221933;
        mesh[16].p1.z = 0.335448;
        mesh[16].p2.x = 0.011163;
        mesh[16].p2.y = -0.226997;
        mesh[16].p2.z = 0.334062;
        mesh[16].p3.x = 0.01717;
        mesh[16].p3.y = -0.220033;
        mesh[16].p3.z = 0.328763;
        mesh[17].p1.x = 0.031688;
        mesh[17].p1.y = -0.226472;
        mesh[17].p1.z = 0.28078;
        mesh[17].p2.x = 0.01869;
        mesh[17].p2.y = -0.214564;
        mesh[17].p2.z = 0.260349;
        mesh[17].p3.x = 0.033794;
        mesh[17].p3.y = -0.169259;
        mesh[17].p3.z = 0.248461;
        mesh[18].p1.x = 0.01869;
        mesh[18].p1.y = -0.214564;
        mesh[18].p1.z = 0.260349;
        mesh[18].p2.x = 0.031688;
        mesh[18].p2.y = -0.226472;
        mesh[18].p2.z = 0.28078;
        mesh[18].p3.x = 0.025034;
        mesh[18].p3.y = -0.224605;
        mesh[18].p3.z = 0.276978;
        mesh[19].p1.x = 0.025034;
        mesh[19].p1.y = -0.224605;
        mesh[19].p1.z = 0.276978;
        mesh[19].p2.x = 0.020072;
        mesh[19].p2.y = -0.245234;
        mesh[19].p2.z = 0.286429;
        mesh[19].p3.x = 0.01869;
        mesh[19].p3.y = -0.214564;
        mesh[19].p3.z = 0.260349;
        mesh[20].p1.x = -0.022417;
        mesh[20].p1.y = -0.136843;
        mesh[20].p1.z = 0.326557;
        mesh[20].p2.x = -7.79E-4;
        mesh[20].p2.y = -0.146982;
        mesh[20].p2.z = 0.344462;
        mesh[20].p3.x = 0.023841;
        mesh[20].p3.y = -0.138332;
        mesh[20].p3.z = 0.326807;
        mesh[21].p1.x = -0.00563;
        mesh[21].p1.y = -0.191954;
        mesh[21].p1.z = 0.225786;
        mesh[21].p2.x = 0.016642;
        mesh[21].p2.y = -0.186885;
        mesh[21].p2.z = 0.224214;
        mesh[21].p3.x = 0.01869;
        mesh[21].p3.y = -0.214564;
        mesh[21].p3.z = 0.260349;
        mesh[22].p1.x = 0.020072;
        mesh[22].p1.y = -0.245234;
        mesh[22].p1.z = 0.286429;
        mesh[22].p2.x = 0.023713;
        mesh[22].p2.y = -0.253701;
        mesh[22].p2.z = 0.282565;
        mesh[22].p3.x = 0.00371;
        mesh[22].p3.y = -0.265544;
        mesh[22].p3.z = 0.295988;
        mesh[23].p1.x = 0.040152;
        mesh[23].p1.y = -0.135693;
        mesh[23].p1.z = 0.020709;
        mesh[23].p2.x = 0.046714;
        mesh[23].p2.y = -0.146287;
        mesh[23].p2.z = 0.006232;
        mesh[23].p3.x = 0.058952;
        mesh[23].p3.y = -0.120225;
        mesh[23].p3.z = 0.013241;
        mesh[24].p1.x = 0.055855;
        mesh[24].p1.y = -0.084068;
        mesh[24].p1.z = 0.206511;
        mesh[24].p2.x = 0.060234;
        mesh[24].p2.y = -0.111718;
        mesh[24].p2.z = 0.141464;
        mesh[24].p3.x = 0.052896;
        mesh[24].p3.y = -0.077728;
        mesh[24].p3.z = 0.169219;
        mesh[25].p1.x = 0.047802;
        mesh[25].p1.y = -0.086047;
        mesh[25].p1.z = 0.120183;
        mesh[25].p2.x = 0.044295;
        mesh[25].p2.y = -0.086567;
        mesh[25].p2.z = 0.047665;
        mesh[25].p3.x = 0.032255;
        mesh[25].p3.y = -0.08755;
        mesh[25].p3.z = 0.048576;
        mesh[26].p1.x = 0.047802;
        mesh[26].p1.y = -0.086047;
        mesh[26].p1.z = 0.120183;
        mesh[26].p2.x = 0.060234;
        mesh[26].p2.y = -0.111718;
        mesh[26].p2.z = 0.141464;
        mesh[26].p3.x = 0.044295;
        mesh[26].p3.y = -0.086567;
        mesh[26].p3.z = 0.047665;
        mesh[27].p1.x = 0.024929;
        mesh[27].p1.y = -0.10802;
        mesh[27].p1.z = 0.027315;
        mesh[27].p2.x = 0.01653;
        mesh[27].p2.y = -0.109257;
        mesh[27].p2.z = 0.120113;
        mesh[27].p3.x = 0.032255;
        mesh[27].p3.y = -0.08755;
        mesh[27].p3.z = 0.048576;
        mesh[28].p1.x = 0.037948;
        mesh[28].p1.y = -0.084855;
        mesh[28].p1.z = 0.13526;
        mesh[28].p2.x = 0.01653;
        mesh[28].p2.y = -0.109257;
        mesh[28].p2.z = 0.120113;
        mesh[28].p3.x = 0.020677;
        mesh[28].p3.y = 0.005463;
        mesh[28].p3.z = 0.12829;
        mesh[29].p1.x = 0.040651;
        mesh[29].p1.y = -0.136387;
        mesh[29].p1.z = 0.244256;
        mesh[29].p2.x = 0.059678;
        mesh[29].p2.y = -0.129737;
        mesh[29].p2.z = 0.199429;
        mesh[29].p3.x = 0.055855;
        mesh[29].p3.y = -0.084068;
        mesh[29].p3.z = 0.206511;
        mesh[30].p1.x = 0.037183;
        mesh[30].p1.y = -0.008779;
        mesh[30].p1.z = 0.249569;
        mesh[30].p2.x = 0.035825;
        mesh[30].p2.y = 0.085919;
        mesh[30].p2.z = 0.251555;
        mesh[30].p3.x = -0.010543;
        mesh[30].p3.y = -0.015599;
        mesh[30].p3.z = 0.262292;
        mesh[31].p1.x = -0.014487;
        mesh[31].p1.y = 0.080367;
        mesh[31].p1.z = 0.261155;
        mesh[31].p2.x = -0.010543;
        mesh[31].p2.y = -0.015599;
        mesh[31].p2.z = 0.262292;
        mesh[31].p3.x = 0.035825;
        mesh[31].p3.y = 0.085919;
        mesh[31].p3.z = 0.251555;
        mesh[32].p1.x = 0.025447;
        mesh[32].p1.y = -0.16317;
        mesh[32].p1.z = 0.14662;
        mesh[32].p2.x = -0.00563;
        mesh[32].p2.y = -0.191954;
        mesh[32].p2.z = 0.225786;
        mesh[32].p3.x = -0.016354;
        mesh[32].p3.y = -0.170451;
        mesh[32].p3.z = 0.160189;
        mesh[33].p1.x = 0.003674;
        mesh[33].p1.y = -0.085117;
        mesh[33].p1.z = 0.272056;
        mesh[33].p2.x = -0.006761;
        mesh[33].p2.y = -0.10839;
        mesh[33].p2.z = 0.287993;
        mesh[33].p3.x = 0.021506;
        mesh[33].p3.y = -0.10025;
        mesh[33].p3.z = 0.266035;
        mesh[34].p1.x = 0.044295;
        mesh[34].p1.y = -0.086567;
        mesh[34].p1.z = 0.047665;
        mesh[34].p2.x = 0.024485;
        mesh[34].p2.y = -0.097618;
        mesh[34].p2.z = 0.011241;
        mesh[34].p3.x = 0.032255;
        mesh[34].p3.y = -0.08755;
        mesh[34].p3.z = 0.048576;
        mesh[35].p1.x = 0.031585;
        mesh[35].p1.y = -0.136678;
        mesh[35].p1.z = 0.002046;
        mesh[35].p2.x = 0.046714;
        mesh[35].p2.y = -0.146287;
        mesh[35].p2.z = 0.006232;
        mesh[35].p3.x = 0.040152;
        mesh[35].p3.y = -0.135693;
        mesh[35].p3.z = 0.020709;
        mesh[36].p1.x = 0.024929;
        mesh[36].p1.y = -0.10802;
        mesh[36].p1.z = 0.027315;
        mesh[36].p2.x = 0.024485;
        mesh[36].p2.y = -0.097618;
        mesh[36].p2.z = 0.011241;
        mesh[36].p3.x = 0.024813;
        mesh[36].p3.y = -0.137638;
        mesh[36].p3.z = 0.006129;
        mesh[37].p1.x = 0.041861;
        mesh[37].p1.y = -0.117404;
        mesh[37].p1.z = 0.027621;
        mesh[37].p2.x = 0.058952;
        mesh[37].p2.y = -0.120225;
        mesh[37].p2.z = 0.013241;
        mesh[37].p3.x = 0.051727;
        mesh[37].p3.y = -0.106852;
        mesh[37].p3.z = 0.0213;
        mesh[38].p1.x = 0.058952;
        mesh[38].p1.y = -0.120225;
        mesh[38].p1.z = 0.013241;
        mesh[38].p2.x = 0.041861;
        mesh[38].p2.y = -0.117404;
        mesh[38].p2.z = 0.027621;
        mesh[38].p3.x = 0.040152;
        mesh[38].p3.y = -0.135693;
        mesh[38].p3.z = 0.020709;
        mesh[39].p1.x = 0.046714;
        mesh[39].p1.y = -0.146287;
        mesh[39].p1.z = 0.006232;
        mesh[39].p2.x = 0.056245;
        mesh[39].p2.y = -0.133422;
        mesh[39].p2.z = 0.00208;
        mesh[39].p3.x = 0.058952;
        mesh[39].p3.y = -0.120225;
        mesh[39].p3.z = 0.013241;
        mesh[40].p1.x = 0.028943;
        mesh[40].p1.y = 0.117552;
        mesh[40].p1.z = 0.068098;
        mesh[40].p2.x = 0.029122;
        mesh[40].p2.y = 0.138151;
        mesh[40].p2.z = 0.027991;
        mesh[40].p3.x = 0.026391;
        mesh[40].p3.y = 0.117277;
        mesh[40].p3.z = 0.021774;
        mesh[41].p1.x = 0.036214;
        mesh[41].p1.y = 0.086871;
        mesh[41].p1.z = 0.00329;
        mesh[41].p2.x = 0.039372;
        mesh[41].p2.y = 0.091188;
        mesh[41].p2.z = 0.014932;
        mesh[41].p3.x = 0.035054;
        mesh[41].p3.y = 0.113693;
        mesh[41].p3.z = 0.02676;
        mesh[42].p1.x = 0.026391;
        mesh[42].p1.y = 0.117277;
        mesh[42].p1.z = 0.021774;
        mesh[42].p2.x = 0.025813;
        mesh[42].p2.y = 0.096198;
        mesh[42].p2.z = -1.51E-4;
        mesh[42].p3.x = 0.035054;
        mesh[42].p3.y = 0.113693;
        mesh[42].p3.z = 0.02676;
        mesh[43].p1.x = 0.037267;
        mesh[43].p1.y = 0.142487;
        mesh[43].p1.z = 0.021146;
        mesh[43].p2.x = 0.055369;
        mesh[43].p2.y = 0.132142;
        mesh[43].p2.z = 0.02032;
        mesh[43].p3.x = 0.053684;
        mesh[43].p3.y = 0.125516;
        mesh[43].p3.z = 0.002388;
        mesh[44].p1.x = 0.030553;
        mesh[44].p1.y = 0.149125;
        mesh[44].p1.z = 0.079185;
        mesh[44].p2.x = 0.028943;
        mesh[44].p2.y = 0.117552;
        mesh[44].p2.z = 0.068098;
        mesh[44].p3.x = 0.011832;
        mesh[44].p3.y = 0.073182;
        mesh[44].p3.z = 0.136522;
        mesh[45].p1.x = 0.035775;
        mesh[45].p1.y = 0.072902;
        mesh[45].p1.z = 0.121575;
        mesh[45].p2.x = 0.011832;
        mesh[45].p2.y = 0.073182;
        mesh[45].p2.z = 0.136522;
        mesh[45].p3.x = 0.028943;
        mesh[45].p3.y = 0.117552;
        mesh[45].p3.z = 0.068098;
        mesh[46].p1.x = 0.047319;
        mesh[46].p1.y = 0.12642;
        mesh[46].p1.z = 0.207304;
        mesh[46].p2.x = 0.058191;
        mesh[46].p2.y = 0.09789;
        mesh[46].p2.z = 0.14702;
        mesh[46].p3.x = 0.051918;
        mesh[46].p3.y = 0.123289;
        mesh[46].p3.z = 0.067837;
        mesh[47].p1.x = 0.055679;
        mesh[47].p1.y = -0.148737;
        mesh[47].p1.z = 0.164846;
        mesh[47].p2.x = 0.048408;
        mesh[47].p2.y = -0.136161;
        mesh[47].p2.z = 0.125389;
        mesh[47].p3.x = 0.060234;
        mesh[47].p3.y = -0.111718;
        mesh[47].p3.z = 0.141464;
        mesh[48].p1.x = 0.024219;
        mesh[48].p1.y = -0.128092;
        mesh[48].p1.z = 0.121928;
        mesh[48].p2.x = 0.01653;
        mesh[48].p2.y = -0.109257;
        mesh[48].p2.z = 0.120113;
        mesh[48].p3.x = 0.024929;
        mesh[48].p3.y = -0.10802;
        mesh[48].p3.z = 0.027315;
        mesh[49].p1.x = 0.032255;
        mesh[49].p1.y = -0.08755;
        mesh[49].p1.z = 0.048576;
        mesh[49].p2.x = 0.01653;
        mesh[49].p2.y = -0.109257;
        mesh[49].p2.z = 0.120113;
        mesh[49].p3.x = 0.047802;
        mesh[49].p3.y = -0.086047;
        mesh[49].p3.z = 0.120183;
        mesh[50].p1.x = 0.024929;
        mesh[50].p1.y = -0.10802;
        mesh[50].p1.z = 0.027315;
        mesh[50].p2.x = 0.048408;
        mesh[50].p2.y = -0.136161;
        mesh[50].p2.z = 0.125389;
        mesh[50].p3.x = 0.024219;
        mesh[50].p3.y = -0.128092;
        mesh[50].p3.z = 0.121928;
        mesh[51].p1.x = 0.041861;
        mesh[51].p1.y = -0.117404;
        mesh[51].p1.z = 0.027621;
        mesh[51].p2.x = 0.024929;
        mesh[51].p2.y = -0.10802;
        mesh[51].p2.z = 0.027315;
        mesh[51].p3.x = 0.040152;
        mesh[51].p3.y = -0.135693;
        mesh[51].p3.z = 0.020709;
        mesh[52].p1.x = 0.046396;
        mesh[52].p1.y = 0.140708;
        mesh[52].p1.z = 0.081159;
        mesh[52].p2.x = 0.032381;
        mesh[52].p2.y = 0.136286;
        mesh[52].p2.z = 0.15941;
        mesh[52].p3.x = 0.051918;
        mesh[52].p3.y = 0.123289;
        mesh[52].p3.z = 0.067837;
        mesh[53].p1.x = 0.053684;
        mesh[53].p1.y = 0.125516;
        mesh[53].p1.z = 0.002388;
        mesh[53].p2.x = 0.055369;
        mesh[53].p2.y = 0.132142;
        mesh[53].p2.z = 0.02032;
        mesh[53].p3.x = 0.057539;
        mesh[53].p3.y = 0.09589;
        mesh[53].p3.z = 3.54E-4;
        mesh[54].p1.x = 0.023355;
        mesh[54].p1.y = 0.127603;
        mesh[54].p1.z = 0.002924;
        mesh[54].p2.x = 0.029122;
        mesh[54].p2.y = 0.138151;
        mesh[54].p2.z = 0.027991;
        mesh[54].p3.x = 0.037267;
        mesh[54].p3.y = 0.142487;
        mesh[54].p3.z = 0.021146;
        mesh[55].p1.x = 0.051918;
        mesh[55].p1.y = 0.123289;
        mesh[55].p1.z = 0.067837;
        mesh[55].p2.x = 0.035054;
        mesh[55].p2.y = 0.113693;
        mesh[55].p2.z = 0.02676;
        mesh[55].p3.x = 0.051959;
        mesh[55].p3.y = 0.116839;
        mesh[55].p3.z = 0.027809;
        mesh[56].p1.x = 0.029122;
        mesh[56].p1.y = 0.138151;
        mesh[56].p1.z = 0.027991;
        mesh[56].p2.x = 0.028943;
        mesh[56].p2.y = 0.117552;
        mesh[56].p2.z = 0.068098;
        mesh[56].p3.x = 0.030553;
        mesh[56].p3.y = 0.149125;
        mesh[56].p3.z = 0.079185;
        mesh[57].p1.x = 0.037183;
        mesh[57].p1.y = -0.008779;
        mesh[57].p1.z = 0.249569;
        mesh[57].p2.x = 0.051207;
        mesh[57].p2.y = 0.016414;
        mesh[57].p2.z = 0.221731;
        mesh[57].p3.x = 0.035825;
        mesh[57].p3.y = 0.085919;
        mesh[57].p3.z = 0.251555;
        mesh[58].p1.x = 0.011832;
        mesh[58].p1.y = 0.073182;
        mesh[58].p1.z = 0.136522;
        mesh[58].p2.x = 0.009345;
        mesh[58].p2.y = 0.108126;
        mesh[58].p2.z = 0.144102;
        mesh[58].p3.x = 0.021511;
        mesh[58].p3.y = 0.123621;
        mesh[58].p3.z = 0.140422;
        mesh[59].p1.x = 0.046714;
        mesh[59].p1.y = -0.146287;
        mesh[59].p1.z = 0.006232;
        mesh[59].p2.x = 0.031585;
        mesh[59].p2.y = -0.136678;
        mesh[59].p2.z = 0.002046;
        mesh[59].p3.x = 0.024485;
        mesh[59].p3.y = -0.097618;
        mesh[59].p3.z = 0.011241;
        mesh[60].p1.x = 0.051727;
        mesh[60].p1.y = -0.106852;
        mesh[60].p1.z = 0.0213;
        mesh[60].p2.x = 0.058952;
        mesh[60].p2.y = -0.120225;
        mesh[60].p2.z = 0.013241;
        mesh[60].p3.x = 0.049755;
        mesh[60].p3.y = -0.105056;
        mesh[60].p3.z = 0.003318;
        mesh[61].p1.x = 0.023355;
        mesh[61].p1.y = 0.127603;
        mesh[61].p1.z = 0.002924;
        mesh[61].p2.x = 0.053684;
        mesh[61].p2.y = 0.125516;
        mesh[61].p2.z = 0.002388;
        mesh[61].p3.x = 0.038873;
        mesh[61].p3.y = 0.105471;
        mesh[61].p3.z = 5.0E-5;
        mesh[62].p1.x = 0.058191;
        mesh[62].p1.y = 0.09789;
        mesh[62].p1.z = 0.14702;
        mesh[62].p2.x = 0.035775;
        mesh[62].p2.y = 0.072902;
        mesh[62].p2.z = 0.121575;
        mesh[62].p3.x = 0.051918;
        mesh[62].p3.y = 0.123289;
        mesh[62].p3.z = 0.067837;
        mesh[63].p1.x = 0.021506;
        mesh[63].p1.y = -0.10025;
        mesh[63].p1.z = 0.266035;
        mesh[63].p2.x = 0.035308;
        mesh[63].p2.y = -0.137464;
        mesh[63].p2.z = 0.306814;
        mesh[63].p3.x = 0.040651;
        mesh[63].p3.y = -0.136387;
        mesh[63].p3.z = 0.244256;
        mesh[64].p1.x = -0.020833;
        mesh[64].p1.y = -0.125757;
        mesh[64].p1.z = 0.124804;
        mesh[64].p2.x = 0.025447;
        mesh[64].p2.y = -0.16317;
        mesh[64].p2.z = 0.14662;
        mesh[64].p3.x = -0.016354;
        mesh[64].p3.y = -0.170451;
        mesh[64].p3.z = 0.160189;
        mesh[65].p1.x = 0.055855;
        mesh[65].p1.y = -0.084068;
        mesh[65].p1.z = 0.206511;
        mesh[65].p2.x = 0.037183;
        mesh[65].p2.y = -0.008779;
        mesh[65].p2.z = 0.249569;
        mesh[65].p3.x = 0.021506;
        mesh[65].p3.y = -0.10025;
        mesh[65].p3.z = 0.266035;
        mesh[66].p1.x = 0.002484;
        mesh[66].p1.y = 0.265509;
        mesh[66].p1.z = 0.210879;
        mesh[66].p2.x = 0.004478;
        mesh[66].p2.y = 0.178564;
        mesh[66].p2.z = 0.232537;
        mesh[66].p3.x = 0.017448;
        mesh[66].p3.y = 0.138736;
        mesh[66].p3.z = 0.227783;
        mesh[67].p1.x = -1.45E-4;
        mesh[67].p1.y = 0.265834;
        mesh[67].p1.z = 0.193142;
        mesh[67].p2.x = 0.017448;
        mesh[67].p2.y = 0.138736;
        mesh[67].p2.z = 0.227783;
        mesh[67].p3.x = -0.003326;
        mesh[67].p3.y = 0.138226;
        mesh[67].p3.z = 0.212925;
        mesh[68].p1.x = -3.26E-4;
        mesh[68].p1.y = 0.141988;
        mesh[68].p1.z = 0.196277;
        mesh[68].p2.x = 0.017448;
        mesh[68].p2.y = 0.138736;
        mesh[68].p2.z = 0.227783;
        mesh[68].p3.x = 0.047319;
        mesh[68].p3.y = 0.12642;
        mesh[68].p3.z = 0.207304;
        mesh[69].p1.x = 0.021506;
        mesh[69].p1.y = -0.10025;
        mesh[69].p1.z = 0.266035;
        mesh[69].p2.x = -0.006761;
        mesh[69].p2.y = -0.10839;
        mesh[69].p2.z = 0.287993;
        mesh[69].p3.x = 0.023841;
        mesh[69].p3.y = -0.138332;
        mesh[69].p3.z = 0.326807;
        mesh[70].p1.x = -0.027136;
        mesh[70].p1.y = 0.041664;
        mesh[70].p1.z = 0.132141;
        mesh[70].p2.x = 0.005073;
        mesh[70].p2.y = 0.031262;
        mesh[70].p2.z = 0.135312;
        mesh[70].p3.x = 0.020677;
        mesh[70].p3.y = 0.005463;
        mesh[70].p3.z = 0.12829;
        mesh[71].p1.x = 0.009345;
        mesh[71].p1.y = 0.108126;
        mesh[71].p1.z = 0.144102;
        mesh[71].p2.x = -3.26E-4;
        mesh[71].p2.y = 0.141988;
        mesh[71].p2.z = 0.196277;
        mesh[71].p3.x = 0.032381;
        mesh[71].p3.y = 0.136286;
        mesh[71].p3.z = 0.15941;
        mesh[72].p1.x = 0.052896;
        mesh[72].p1.y = -0.077728;
        mesh[72].p1.z = 0.169219;
        mesh[72].p2.x = 0.0337;
        mesh[72].p2.y = 0.048755;
        mesh[72].p2.z = 0.136203;
        mesh[72].p3.x = 0.051207;
        mesh[72].p3.y = 0.016414;
        mesh[72].p3.z = 0.221731;
        mesh[73].p1.x = 0.051918;
        mesh[73].p1.y = 0.123289;
        mesh[73].p1.z = 0.067837;
        mesh[73].p2.x = 0.032381;
        mesh[73].p2.y = 0.136286;
        mesh[73].p2.z = 0.15941;
        mesh[73].p3.x = 0.047319;
        mesh[73].p3.y = 0.12642;
        mesh[73].p3.z = 0.207304;
        mesh[74].p1.x = 0.0337;
        mesh[74].p1.y = 0.048755;
        mesh[74].p1.z = 0.136203;
        mesh[74].p2.x = 0.052896;
        mesh[74].p2.y = -0.077728;
        mesh[74].p2.z = 0.169219;
        mesh[74].p3.x = 0.020677;
        mesh[74].p3.y = 0.005463;
        mesh[74].p3.z = 0.12829;
        mesh[75].p1.x = 0.025447;
        mesh[75].p1.y = -0.16317;
        mesh[75].p1.z = 0.14662;
        mesh[75].p2.x = -0.020833;
        mesh[75].p2.y = -0.125757;
        mesh[75].p2.z = 0.124804;
        mesh[75].p3.x = 0.01653;
        mesh[75].p3.y = -0.109257;
        mesh[75].p3.z = 0.120113;
        mesh[76].p1.x = 0.055679;
        mesh[76].p1.y = -0.148737;
        mesh[76].p1.z = 0.164846;
        mesh[76].p2.x = 0.059678;
        mesh[76].p2.y = -0.129737;
        mesh[76].p2.z = 0.199429;
        mesh[76].p3.x = 0.033794;
        mesh[76].p3.y = -0.169259;
        mesh[76].p3.z = 0.248461;
        mesh[77].p1.x = 0.055679;
        mesh[77].p1.y = -0.148737;
        mesh[77].p1.z = 0.164846;
        mesh[77].p2.x = 0.016642;
        mesh[77].p2.y = -0.186885;
        mesh[77].p2.z = 0.224214;
        mesh[77].p3.x = 0.025447;
        mesh[77].p3.y = -0.16317;
        mesh[77].p3.z = 0.14662;
        mesh[78].p1.x = 0.033794;
        mesh[78].p1.y = -0.169259;
        mesh[78].p1.z = 0.248461;
        mesh[78].p2.x = 0.016642;
        mesh[78].p2.y = -0.186885;
        mesh[78].p2.z = 0.224214;
        mesh[78].p3.x = 0.055679;
        mesh[78].p3.y = -0.148737;
        mesh[78].p3.z = 0.164846;
        mesh[79].p1.x = 0.060234;
        mesh[79].p1.y = -0.111718;
        mesh[79].p1.z = 0.141464;
        mesh[79].p2.x = 0.055855;
        mesh[79].p2.y = -0.084068;
        mesh[79].p2.z = 0.206511;
        mesh[79].p3.x = 0.059678;
        mesh[79].p3.y = -0.129737;
        mesh[79].p3.z = 0.199429;
        mesh[80].p1.x = 0.043042;
        mesh[80].p1.y = -0.1738;
        mesh[80].p1.z = 0.335948;
        mesh[80].p2.x = 0.066431;
        mesh[80].p2.y = -0.16482;
        mesh[80].p2.z = 0.292051;
        mesh[80].p3.x = 0.039196;
        mesh[80].p3.y = -0.156989;
        mesh[80].p3.z = 0.339227;
        mesh[81].p1.x = 0.053788;
        mesh[81].p1.y = -0.129921;
        mesh[81].p1.z = 0.303985;
        mesh[81].p2.x = 0.066431;
        mesh[81].p2.y = -0.16482;
        mesh[81].p2.z = 0.292051;
        mesh[81].p3.x = 0.063421;
        mesh[81].p3.y = -0.164003;
        mesh[81].p3.z = 0.271249;
        mesh[82].p1.x = 0.051207;
        mesh[82].p1.y = 0.016414;
        mesh[82].p1.z = 0.221731;
        mesh[82].p2.x = 0.035775;
        mesh[82].p2.y = 0.072902;
        mesh[82].p2.z = 0.121575;
        mesh[82].p3.x = 0.058191;
        mesh[82].p3.y = 0.09789;
        mesh[82].p3.z = 0.14702;
        mesh[83].p1.x = 0.046396;
        mesh[83].p1.y = 0.140708;
        mesh[83].p1.z = 0.081159;
        mesh[83].p2.x = 0.051918;
        mesh[83].p2.y = 0.123289;
        mesh[83].p2.z = 0.067837;
        mesh[83].p3.x = 0.055369;
        mesh[83].p3.y = 0.132142;
        mesh[83].p3.z = 0.02032;
        mesh[84].p1.x = 0.017448;
        mesh[84].p1.y = 0.138736;
        mesh[84].p1.z = 0.227783;
        mesh[84].p2.x = -0.014487;
        mesh[84].p2.y = 0.080367;
        mesh[84].p2.z = 0.261155;
        mesh[84].p3.x = 0.035825;
        mesh[84].p3.y = 0.085919;
        mesh[84].p3.z = 0.251555;
        mesh[85].p1.x = 0.017448;
        mesh[85].p1.y = 0.138736;
        mesh[85].p1.z = 0.227783;
        mesh[85].p2.x = 0.035825;
        mesh[85].p2.y = 0.085919;
        mesh[85].p2.z = 0.251555;
        mesh[85].p3.x = 0.047319;
        mesh[85].p3.y = 0.12642;
        mesh[85].p3.z = 0.207304;
        mesh[86].p1.x = 0.023841;
        mesh[86].p1.y = -0.138332;
        mesh[86].p1.z = 0.326807;
        mesh[86].p2.x = 0.039196;
        mesh[86].p2.y = -0.156989;
        mesh[86].p2.z = 0.339227;
        mesh[86].p3.x = 0.049253;
        mesh[86].p3.y = -0.135547;
        mesh[86].p3.z = 0.325422;
        mesh[87].p1.x = 0.032381;
        mesh[87].p1.y = 0.136286;
        mesh[87].p1.z = 0.15941;
        mesh[87].p2.x = 0.046396;
        mesh[87].p2.y = 0.140708;
        mesh[87].p2.z = 0.081159;
        mesh[87].p3.x = 0.021511;
        mesh[87].p3.y = 0.123621;
        mesh[87].p3.z = 0.140422;
        mesh[88].p1.x = 0.030553;
        mesh[88].p1.y = 0.149125;
        mesh[88].p1.z = 0.079185;
        mesh[88].p2.x = 0.021511;
        mesh[88].p2.y = 0.123621;
        mesh[88].p2.z = 0.140422;
        mesh[88].p3.x = 0.046396;
        mesh[88].p3.y = 0.140708;
        mesh[88].p3.z = 0.081159;
        mesh[89].p1.x = 0.046396;
        mesh[89].p1.y = 0.140708;
        mesh[89].p1.z = 0.081159;
        mesh[89].p2.x = 0.037267;
        mesh[89].p2.y = 0.142487;
        mesh[89].p2.z = 0.021146;
        mesh[89].p3.x = 0.030553;
        mesh[89].p3.y = 0.149125;
        mesh[89].p3.z = 0.079185;
        mesh[90].p1.x = 0.011832;
        mesh[90].p1.y = 0.073182;
        mesh[90].p1.z = 0.136522;
        mesh[90].p2.x = 0.035775;
        mesh[90].p2.y = 0.072902;
        mesh[90].p2.z = 0.121575;
        mesh[90].p3.x = 0.005073;
        mesh[90].p3.y = 0.031262;
        mesh[90].p3.z = 0.135312;
        mesh[91].p1.x = 0.024485;
        mesh[91].p1.y = -0.097618;
        mesh[91].p1.z = 0.011241;
        mesh[91].p2.x = 0.049755;
        mesh[91].p2.y = -0.105056;
        mesh[91].p2.z = 0.003318;
        mesh[91].p3.x = 0.046714;
        mesh[91].p3.y = -0.146287;
        mesh[91].p3.z = 0.006232;
        mesh[92].p1.x = 0.025034;
        mesh[92].p1.y = -0.224605;
        mesh[92].p1.z = 0.276978;
        mesh[92].p2.x = 0.031688;
        mesh[92].p2.y = -0.226472;
        mesh[92].p2.z = 0.28078;
        mesh[92].p3.x = 0.020072;
        mesh[92].p3.y = -0.245234;
        mesh[92].p3.z = 0.286429;
        mesh[93].p1.x = -0.0117;
        mesh[93].p1.y = -0.264624;
        mesh[93].p1.z = 0.32503;
        mesh[93].p2.x = -0.005389;
        mesh[93].p2.y = -0.264282;
        mesh[93].p2.z = 0.316817;
        mesh[93].p3.x = 0.002083;
        mesh[93].p3.y = -0.267639;
        mesh[93].p3.z = 0.312103;
        mesh[94].p1.x = -0.0117;
        mesh[94].p1.y = -0.264624;
        mesh[94].p1.z = 0.32503;
        mesh[94].p2.x = -0.009455;
        mesh[94].p2.y = -0.266105;
        mesh[94].p2.z = 0.31138;
        mesh[94].p3.x = -0.005389;
        mesh[94].p3.y = -0.264282;
        mesh[94].p3.z = 0.316817;
        mesh[95].p1.x = 0.002083;
        mesh[95].p1.y = -0.267639;
        mesh[95].p1.z = 0.312103;
        mesh[95].p2.x = 0.007313;
        mesh[95].p2.y = -0.265396;
        mesh[95].p2.z = 0.327658;
        mesh[95].p3.x = -0.0117;
        mesh[95].p3.y = -0.264624;
        mesh[95].p3.z = 0.32503;
        mesh[96].p1.x = -0.026484;
        mesh[96].p1.y = -0.259787;
        mesh[96].p1.z = 0.297418;
        mesh[96].p2.x = -0.01772;
        mesh[96].p2.y = -0.262421;
        mesh[96].p2.z = 0.287041;
        mesh[96].p3.x = -0.009455;
        mesh[96].p3.y = -0.266105;
        mesh[96].p3.z = 0.31138;
        mesh[97].p1.x = -0.017193;
        mesh[97].p1.y = -0.210889;
        mesh[97].p1.z = 0.257864;
        mesh[97].p2.x = -0.024261;
        mesh[97].p2.y = -0.22799;
        mesh[97].p2.z = 0.280917;
        mesh[97].p3.x = -0.030009;
        mesh[97].p3.y = -0.229325;
        mesh[97].p3.z = 0.281517;
        mesh[98].p1.x = -0.02214;
        mesh[98].p1.y = -0.257645;
        mesh[98].p1.z = 0.316165;
        mesh[98].p2.x = -0.009896;
        mesh[98].p2.y = -0.227817;
        mesh[98].p2.z = 0.334942;
        mesh[98].p3.x = -0.029006;
        mesh[98].p3.y = -0.229525;
        mesh[98].p3.z = 0.311337;
        mesh[99].p1.x = -0.032359;
        mesh[99].p1.y = -0.171143;
        mesh[99].p1.z = 0.240091;
        mesh[99].p2.x = -0.030009;
        mesh[99].p2.y = -0.229325;
        mesh[99].p2.z = 0.281517;
        mesh[99].p3.x = -0.039008;
        mesh[99].p3.y = -0.158249;
        mesh[99].p3.z = 0.27949;
        mesh[100].p1.x = -0.030009;
        mesh[100].p1.y = -0.229325;
        mesh[100].p1.z = 0.281517;
        mesh[100].p2.x = -0.039859;
        mesh[100].p2.y = -0.204217;
        mesh[100].p2.z = 0.325851;
        mesh[100].p3.x = -0.039008;
        mesh[100].p3.y = -0.158249;
        mesh[100].p3.z = 0.27949;
        mesh[101].p1.x = -0.039859;
        mesh[101].p1.y = -0.204217;
        mesh[101].p1.z = 0.325851;
        mesh[101].p2.x = -0.048469;
        mesh[101].p2.y = -0.158545;
        mesh[101].p2.z = 0.319398;
        mesh[101].p3.x = -0.039008;
        mesh[101].p3.y = -0.158249;
        mesh[101].p3.z = 0.27949;
        mesh[102].p1.x = -0.037814;
        mesh[102].p1.y = -0.136474;
        mesh[102].p1.z = 0.309237;
        mesh[102].p2.x = -0.048469;
        mesh[102].p2.y = -0.158545;
        mesh[102].p2.z = 0.319398;
        mesh[102].p3.x = -0.053788;
        mesh[102].p3.y = -0.129921;
        mesh[102].p3.z = 0.303985;
        mesh[103].p1.x = -0.039197;
        mesh[103].p1.y = -0.156989;
        mesh[103].p1.z = 0.339227;
        mesh[103].p2.x = -7.79E-4;
        mesh[103].p2.y = -0.146982;
        mesh[103].p2.z = 0.344462;
        mesh[103].p3.x = -0.022417;
        mesh[103].p3.y = -0.136843;
        mesh[103].p3.z = 0.326557;
        mesh[104].p1.x = -0.028845;
        mesh[104].p1.y = -0.219951;
        mesh[104].p1.z = 0.330842;
        mesh[104].p2.x = -0.021594;
        mesh[104].p2.y = -0.217944;
        mesh[104].p2.z = 0.346372;
        mesh[104].p3.x = -0.039859;
        mesh[104].p3.y = -0.204217;
        mesh[104].p3.z = 0.325851;
        mesh[105].p1.x = -0.009896;
        mesh[105].p1.y = -0.227817;
        mesh[105].p1.z = 0.334942;
        mesh[105].p2.x = -0.021594;
        mesh[105].p2.y = -0.217944;
        mesh[105].p2.z = 0.346372;
        mesh[105].p3.x = -0.028845;
        mesh[105].p3.y = -0.219951;
        mesh[105].p3.z = 0.330842;
        mesh[106].p1.x = -0.009896;
        mesh[106].p1.y = -0.227817;
        mesh[106].p1.z = 0.334942;
        mesh[106].p2.x = -0.028845;
        mesh[106].p2.y = -0.219951;
        mesh[106].p2.z = 0.330842;
        mesh[106].p3.x = -0.018638;
        mesh[106].p3.y = -0.221672;
        mesh[106].p3.z = 0.328157;
        mesh[107].p1.x = -0.017193;
        mesh[107].p1.y = -0.210889;
        mesh[107].p1.z = 0.257864;
        mesh[107].p2.x = -0.013919;
        mesh[107].p2.y = -0.25944;
        mesh[107].p2.z = 0.288499;
        mesh[107].p3.x = -0.024261;
        mesh[107].p3.y = -0.22799;
        mesh[107].p3.z = 0.280917;
        mesh[108].p1.x = -0.043115;
        mesh[108].p1.y = -0.141878;
        mesh[108].p1.z = 0.237418;
        mesh[108].p2.x = -0.032359;
        mesh[108].p2.y = -0.171143;
        mesh[108].p2.z = 0.240091;
        mesh[108].p3.x = -0.039008;
        mesh[108].p3.y = -0.158249;
        mesh[108].p3.z = 0.27949;
        mesh[109].p1.x = -0.030009;
        mesh[109].p1.y = -0.229325;
        mesh[109].p1.z = 0.281517;
        mesh[109].p2.x = -0.032359;
        mesh[109].p2.y = -0.171143;
        mesh[109].p2.z = 0.240091;
        mesh[109].p3.x = -0.00563;
        mesh[109].p3.y = -0.191954;
        mesh[109].p3.z = 0.225786;
        mesh[110].p1.x = 0.00371;
        mesh[110].p1.y = -0.265544;
        mesh[110].p1.z = 0.295988;
        mesh[110].p2.x = -0.013919;
        mesh[110].p2.y = -0.25944;
        mesh[110].p2.z = 0.288499;
        mesh[110].p3.x = 5.65E-4;
        mesh[110].p3.y = -0.260455;
        mesh[110].p3.z = 0.278232;
        mesh[111].p1.x = -0.01772;
        mesh[111].p1.y = -0.262421;
        mesh[111].p1.z = 0.287041;
        mesh[111].p2.x = -0.024261;
        mesh[111].p2.y = -0.22799;
        mesh[111].p2.z = 0.280917;
        mesh[111].p3.x = -0.013919;
        mesh[111].p3.y = -0.25944;
        mesh[111].p3.z = 0.288499;
        mesh[112].p1.x = -0.030009;
        mesh[112].p1.y = -0.229325;
        mesh[112].p1.z = 0.281517;
        mesh[112].p2.x = -0.024261;
        mesh[112].p2.y = -0.22799;
        mesh[112].p2.z = 0.280917;
        mesh[112].p3.x = -0.01772;
        mesh[112].p3.y = -0.262421;
        mesh[112].p3.z = 0.287041;
        mesh[113].p1.x = 5.65E-4;
        mesh[113].p1.y = -0.260455;
        mesh[113].p1.z = 0.278232;
        mesh[113].p2.x = -0.013919;
        mesh[113].p2.y = -0.25944;
        mesh[113].p2.z = 0.288499;
        mesh[113].p3.x = -0.017193;
        mesh[113].p3.y = -0.210889;
        mesh[113].p3.z = 0.257864;
        mesh[114].p1.x = -0.039606;
        mesh[114].p1.y = 0.087492;
        mesh[114].p1.z = 0.004376;
        mesh[114].p2.x = -0.047607;
        mesh[114].p2.y = 0.104629;
        mesh[114].p2.z = 0.020424;
        mesh[114].p3.x = -0.04562;
        mesh[114].p3.y = 0.088864;
        mesh[114].p3.z = -0.001593;
        mesh[115].p1.x = -0.059617;
        mesh[115].p1.y = -0.1185;
        mesh[115].p1.z = 0.120599;
        mesh[115].p2.x = -0.056427;
        mesh[115].p2.y = -0.085244;
        mesh[115].p2.z = 0.187472;
        mesh[115].p3.x = -0.038493;
        mesh[115].p3.y = -0.084398;
        mesh[115].p3.z = 0.135693;
        mesh[116].p1.x = -0.051674;
        mesh[116].p1.y = -0.109716;
        mesh[116].p1.z = 0.023234;
        mesh[116].p2.x = -0.047173;
        mesh[116].p2.y = -0.08777;
        mesh[116].p2.z = 0.116024;
        mesh[116].p3.x = -0.040724;
        mesh[116].p3.y = -0.084142;
        mesh[116].p3.z = 0.048356;
        mesh[117].p1.x = -0.047107;
        mesh[117].p1.y = -0.005026;
        mesh[117].p1.z = 0.233809;
        mesh[117].p2.x = -0.056427;
        mesh[117].p2.y = -0.085244;
        mesh[117].p2.z = 0.187472;
        mesh[117].p3.x = -0.043115;
        mesh[117].p3.y = -0.141878;
        mesh[117].p3.z = 0.237418;
        mesh[118].p1.x = 0.01653;
        mesh[118].p1.y = -0.109257;
        mesh[118].p1.z = 0.120113;
        mesh[118].p2.x = -0.020833;
        mesh[118].p2.y = -0.125757;
        mesh[118].p2.z = 0.124804;
        mesh[118].p3.x = -0.017058;
        mesh[118].p3.y = -0.107845;
        mesh[118].p3.z = 0.120886;
        mesh[119].p1.x = -0.031352;
        mesh[119].p1.y = -0.118191;
        mesh[119].p1.z = 0.050941;
        mesh[119].p2.x = -0.029086;
        mesh[119].p2.y = -0.121207;
        mesh[119].p2.z = 0.020285;
        mesh[119].p3.x = -0.021311;
        mesh[119].p3.y = -0.109992;
        mesh[119].p3.z = 0.007536;
        mesh[120].p1.x = -0.04676;
        mesh[120].p1.y = -0.146333;
        mesh[120].p1.z = 0.006207;
        mesh[120].p2.x = -0.02965;
        mesh[120].p2.y = -0.136573;
        mesh[120].p2.z = 0.00164;
        mesh[120].p3.x = -0.029086;
        mesh[120].p3.y = -0.121207;
        mesh[120].p3.z = 0.020285;
        mesh[121].p1.x = -0.061288;
        mesh[121].p1.y = -0.130257;
        mesh[121].p1.z = 0.005966;
        mesh[121].p2.x = -0.043082;
        mesh[121].p2.y = -0.121347;
        mesh[121].p2.z = 0.0261;
        mesh[121].p3.x = -0.051674;
        mesh[121].p3.y = -0.109716;
        mesh[121].p3.z = 0.023234;
        mesh[122].p1.x = -0.039606;
        mesh[122].p1.y = 0.087492;
        mesh[122].p1.z = 0.004376;
        mesh[122].p2.x = -0.033126;
        mesh[122].p2.y = 0.084624;
        mesh[122].p2.z = 0.003634;
        mesh[122].p3.x = -0.047607;
        mesh[122].p3.y = 0.104629;
        mesh[122].p3.z = 0.020424;
        mesh[123].p1.x = -0.023356;
        mesh[123].p1.y = 0.126501;
        mesh[123].p1.z = 0.002503;
        mesh[123].p2.x = -0.027847;
        mesh[123].p2.y = 0.120028;
        mesh[123].p2.z = 0.02675;
        mesh[123].p3.x = -0.022864;
        mesh[123].p3.y = 0.094117;
        mesh[123].p3.z = -3.86E-4;
        mesh[124].p1.x = -0.052415;
        mesh[124].p1.y = 0.095648;
        mesh[124].p1.z = -3.5E-5;
        mesh[124].p2.x = -0.047607;
        mesh[124].p2.y = 0.104629;
        mesh[124].p2.z = 0.020424;
        mesh[124].p3.x = -0.057904;
        mesh[124].p3.y = 0.124382;
        mesh[124].p3.z = 0.020753;
        mesh[125].p1.x = -0.052415;
        mesh[125].p1.y = 0.095648;
        mesh[125].p1.z = -3.5E-5;
        mesh[125].p2.x = -0.04562;
        mesh[125].p2.y = 0.088864;
        mesh[125].p2.z = -0.001593;
        mesh[125].p3.x = -0.047607;
        mesh[125].p3.y = 0.104629;
        mesh[125].p3.z = 0.020424;
        mesh[126].p1.x = -0.012166;
        mesh[126].p1.y = 0.073247;
        mesh[126].p1.z = 0.136132;
        mesh[126].p2.x = -0.045197;
        mesh[126].p2.y = 0.121199;
        mesh[126].p2.z = 0.055852;
        mesh[126].p3.x = -0.023133;
        mesh[126].p3.y = 0.120981;
        mesh[126].p3.z = 0.078476;
        mesh[127].p1.x = -0.052194;
        mesh[127].p1.y = 0.037269;
        mesh[127].p1.z = 0.199662;
        mesh[127].p2.x = -0.051287;
        mesh[127].p2.y = 0.123923;
        mesh[127].p2.z = 0.200622;
        mesh[127].p3.x = -0.05208;
        mesh[127].p3.y = 0.124926;
        mesh[127].p3.z = 0.07032;
        mesh[128].p1.x = -0.024667;
        mesh[128].p1.y = -0.101688;
        mesh[128].p1.z = 0.262395;
        mesh[128].p2.x = -0.047107;
        mesh[128].p2.y = -0.005026;
        mesh[128].p2.z = 0.233809;
        mesh[128].p3.x = -0.043115;
        mesh[128].p3.y = -0.141878;
        mesh[128].p3.z = 0.237418;
        mesh[129].p1.x = -0.04672;
        mesh[129].p1.y = 0.093756;
        mesh[129].p1.z = 0.235686;
        mesh[129].p2.x = -0.052194;
        mesh[129].p2.y = 0.037269;
        mesh[129].p2.z = 0.199662;
        mesh[129].p3.x = -0.047107;
        mesh[129].p3.y = -0.005026;
        mesh[129].p3.z = 0.233809;
        mesh[130].p1.x = -0.031352;
        mesh[130].p1.y = -0.118191;
        mesh[130].p1.z = 0.050941;
        mesh[130].p2.x = -0.043082;
        mesh[130].p2.y = -0.121347;
        mesh[130].p2.z = 0.0261;
        mesh[130].p3.x = -0.029086;
        mesh[130].p3.y = -0.121207;
        mesh[130].p3.z = 0.020285;
        mesh[131].p1.x = -0.031246;
        mesh[131].p1.y = 0.135959;
        mesh[131].p1.z = 0.153033;
        mesh[131].p2.x = -0.037203;
        mesh[131].p2.y = 0.138321;
        mesh[131].p2.z = 0.100764;
        mesh[131].p3.x = -0.05208;
        mesh[131].p3.y = 0.124926;
        mesh[131].p3.z = 0.07032;
        mesh[132].p1.x = -0.031653;
        mesh[132].p1.y = 0.105782;
        mesh[132].p1.z = 0.018291;
        mesh[132].p2.x = -0.047607;
        mesh[132].p2.y = 0.104629;
        mesh[132].p2.z = 0.020424;
        mesh[132].p3.x = -0.033126;
        mesh[132].p3.y = 0.084624;
        mesh[132].p3.z = 0.003634;
        mesh[133].p1.x = -0.030328;
        mesh[133].p1.y = 0.140036;
        mesh[133].p1.z = 0.018044;
        mesh[133].p2.x = -0.027847;
        mesh[133].p2.y = 0.120028;
        mesh[133].p2.z = 0.02675;
        mesh[133].p3.x = -0.023356;
        mesh[133].p3.y = 0.126501;
        mesh[133].p3.z = 0.002503;
        mesh[134].p1.x = -0.045197;
        mesh[134].p1.y = 0.121199;
        mesh[134].p1.z = 0.055852;
        mesh[134].p2.x = -0.05208;
        mesh[134].p2.y = 0.124926;
        mesh[134].p2.z = 0.07032;
        mesh[134].p3.x = -0.057904;
        mesh[134].p3.y = 0.124382;
        mesh[134].p3.z = 0.020753;
        mesh[135].p1.x = -0.045197;
        mesh[135].p1.y = 0.121199;
        mesh[135].p1.z = 0.055852;
        mesh[135].p2.x = -0.027847;
        mesh[135].p2.y = 0.120028;
        mesh[135].p2.z = 0.02675;
        mesh[135].p3.x = -0.023133;
        mesh[135].p3.y = 0.120981;
        mesh[135].p3.z = 0.078476;
        mesh[136].p1.x = -0.05208;
        mesh[136].p1.y = 0.124926;
        mesh[136].p1.z = 0.07032;
        mesh[136].p2.x = -0.037603;
        mesh[136].p2.y = 0.15207;
        mesh[136].p2.z = 0.084055;
        mesh[136].p3.x = -0.057904;
        mesh[136].p3.y = 0.124382;
        mesh[136].p3.z = 0.020753;
        mesh[137].p1.x = -0.044942;
        mesh[137].p1.y = 0.141872;
        mesh[137].p1.z = 0.019356;
        mesh[137].p2.x = -0.057904;
        mesh[137].p2.y = 0.124382;
        mesh[137].p2.z = 0.020753;
        mesh[137].p3.x = -0.037603;
        mesh[137].p3.y = 0.15207;
        mesh[137].p3.z = 0.084055;
        mesh[138].p1.x = -0.047107;
        mesh[138].p1.y = -0.005026;
        mesh[138].p1.z = 0.233809;
        mesh[138].p2.x = -0.010543;
        mesh[138].p2.y = -0.015599;
        mesh[138].p2.z = 0.262292;
        mesh[138].p3.x = -0.014487;
        mesh[138].p3.y = 0.080367;
        mesh[138].p3.z = 0.261155;
        mesh[139].p1.x = 0.011832;
        mesh[139].p1.y = 0.073182;
        mesh[139].p1.z = 0.136522;
        mesh[139].p2.x = -0.012166;
        mesh[139].p2.y = 0.073247;
        mesh[139].p2.z = 0.136132;
        mesh[139].p3.x = -0.009247;
        mesh[139].p3.y = 0.106663;
        mesh[139].p3.z = 0.145238;
        mesh[140].p1.x = -0.030099;
        mesh[140].p1.y = 0.074776;
        mesh[140].p1.z = 0.120707;
        mesh[140].p2.x = -0.027136;
        mesh[140].p2.y = 0.041664;
        mesh[140].p2.z = 0.132141;
        mesh[140].p3.x = -0.047097;
        mesh[140].p3.y = 0.052817;
        mesh[140].p3.z = 0.162408;
        mesh[141].p1.x = -0.052194;
        mesh[141].p1.y = 0.037269;
        mesh[141].p1.z = 0.199662;
        mesh[141].p2.x = -0.045197;
        mesh[141].p2.y = 0.121199;
        mesh[141].p2.z = 0.055852;
        mesh[141].p3.x = -0.047097;
        mesh[141].p3.y = 0.052817;
        mesh[141].p3.z = 0.162408;
        mesh[142].p1.x = -0.039008;
        mesh[142].p1.y = -0.158249;
        mesh[142].p1.z = 0.27949;
        mesh[142].p2.x = -0.024667;
        mesh[142].p2.y = -0.101688;
        mesh[142].p2.z = 0.262395;
        mesh[142].p3.x = -0.043115;
        mesh[142].p3.y = -0.141878;
        mesh[142].p3.z = 0.237418;
        mesh[143].p1.x = -0.017841;
        mesh[143].p1.y = 0.13652;
        mesh[143].p1.z = 0.230347;
        mesh[143].p2.x = -0.04672;
        mesh[143].p2.y = 0.093756;
        mesh[143].p2.z = 0.235686;
        mesh[143].p3.x = -0.014487;
        mesh[143].p3.y = 0.080367;
        mesh[143].p3.z = 0.261155;
        mesh[144].p1.x = -0.042681;
        mesh[144].p1.y = -0.137565;
        mesh[144].p1.z = 0.125192;
        mesh[144].p2.x = -0.020833;
        mesh[144].p2.y = -0.125757;
        mesh[144].p2.z = 0.124804;
        mesh[144].p3.x = -0.016354;
        mesh[144].p3.y = -0.170451;
        mesh[144].p3.z = 0.160189;
        mesh[145].p1.x = -0.024667;
        mesh[145].p1.y = -0.101688;
        mesh[145].p1.z = 0.262395;
        mesh[145].p2.x = -0.010543;
        mesh[145].p2.y = -0.015599;
        mesh[145].p2.z = 0.262292;
        mesh[145].p3.x = -0.047107;
        mesh[145].p3.y = -0.005026;
        mesh[145].p3.z = 0.233809;
        mesh[146].p1.x = -0.010941;
        mesh[146].p1.y = 0.235894;
        mesh[146].p1.z = 0.207299;
        mesh[146].p2.x = -0.017841;
        mesh[146].p2.y = 0.13652;
        mesh[146].p2.z = 0.230347;
        mesh[146].p3.x = 0.002484;
        mesh[146].p3.y = 0.265509;
        mesh[146].p3.z = 0.210879;
        mesh[147].p1.x = -0.017841;
        mesh[147].p1.y = 0.13652;
        mesh[147].p1.z = 0.230347;
        mesh[147].p2.x = 0.004478;
        mesh[147].p2.y = 0.178564;
        mesh[147].p2.z = 0.232537;
        mesh[147].p3.x = 0.002484;
        mesh[147].p3.y = 0.265509;
        mesh[147].p3.z = 0.210879;
        mesh[148].p1.x = -0.010941;
        mesh[148].p1.y = 0.235894;
        mesh[148].p1.z = 0.207299;
        mesh[148].p2.x = 0.002484;
        mesh[148].p2.y = 0.265509;
        mesh[148].p2.z = 0.210879;
        mesh[148].p3.x = -1.45E-4;
        mesh[148].p3.y = 0.265834;
        mesh[148].p3.z = 0.193142;
        mesh[149].p1.x = -0.017841;
        mesh[149].p1.y = 0.13652;
        mesh[149].p1.z = 0.230347;
        mesh[149].p2.x = -3.26E-4;
        mesh[149].p2.y = 0.141988;
        mesh[149].p2.z = 0.196277;
        mesh[149].p3.x = -0.051287;
        mesh[149].p3.y = 0.123923;
        mesh[149].p3.z = 0.200622;
        mesh[150].p1.x = -0.027136;
        mesh[150].p1.y = 0.041664;
        mesh[150].p1.z = 0.132141;
        mesh[150].p2.x = -0.012945;
        mesh[150].p2.y = 0.045589;
        mesh[150].p2.z = 0.13692;
        mesh[150].p3.x = 0.005073;
        mesh[150].p3.y = 0.031262;
        mesh[150].p3.z = 0.135312;
        mesh[151].p1.x = -0.031246;
        mesh[151].p1.y = 0.135959;
        mesh[151].p1.z = 0.153033;
        mesh[151].p2.x = -3.26E-4;
        mesh[151].p2.y = 0.141988;
        mesh[151].p2.z = 0.196277;
        mesh[151].p3.x = -0.009247;
        mesh[151].p3.y = 0.106663;
        mesh[151].p3.z = 0.145238;
        mesh[152].p1.x = -0.031246;
        mesh[152].p1.y = 0.135959;
        mesh[152].p1.z = 0.153033;
        mesh[152].p2.x = -0.05208;
        mesh[152].p2.y = 0.124926;
        mesh[152].p2.z = 0.07032;
        mesh[152].p3.x = -0.051287;
        mesh[152].p3.y = 0.123923;
        mesh[152].p3.z = 0.200622;
        mesh[153].p1.x = -0.038493;
        mesh[153].p1.y = -0.084398;
        mesh[153].p1.z = 0.135693;
        mesh[153].p2.x = -0.047097;
        mesh[153].p2.y = 0.052817;
        mesh[153].p2.z = 0.162408;
        mesh[153].p3.x = -0.027136;
        mesh[153].p3.y = 0.041664;
        mesh[153].p3.z = 0.132141;
        mesh[154].p1.x = -0.038493;
        mesh[154].p1.y = -0.084398;
        mesh[154].p1.z = 0.135693;
        mesh[154].p2.x = -0.056427;
        mesh[154].p2.y = -0.085244;
        mesh[154].p2.z = 0.187472;
        mesh[154].p3.x = -0.047097;
        mesh[154].p3.y = 0.052817;
        mesh[154].p3.z = 0.162408;
        mesh[155].p1.x = -0.042681;
        mesh[155].p1.y = -0.137565;
        mesh[155].p1.z = 0.125192;
        mesh[155].p2.x = -0.016354;
        mesh[155].p2.y = -0.170451;
        mesh[155].p2.z = 0.160189;
        mesh[155].p3.x = -0.057323;
        mesh[155].p3.y = -0.147545;
        mesh[155].p3.z = 0.177989;
        mesh[156].p1.x = -0.00563;
        mesh[156].p1.y = -0.191954;
        mesh[156].p1.z = 0.225786;
        mesh[156].p2.x = -0.057323;
        mesh[156].p2.y = -0.147545;
        mesh[156].p2.z = 0.177989;
        mesh[156].p3.x = -0.016354;
        mesh[156].p3.y = -0.170451;
        mesh[156].p3.z = 0.160189;
        mesh[157].p1.x = -0.043115;
        mesh[157].p1.y = -0.141878;
        mesh[157].p1.z = 0.237418;
        mesh[157].p2.x = -0.056427;
        mesh[157].p2.y = -0.085244;
        mesh[157].p2.z = 0.187472;
        mesh[157].p3.x = -0.061648;
        mesh[157].p3.y = -0.116842;
        mesh[157].p3.z = 0.194265;
        mesh[158].p1.x = -0.061648;
        mesh[158].p1.y = -0.116842;
        mesh[158].p1.z = 0.194265;
        mesh[158].p2.x = -0.059617;
        mesh[158].p2.y = -0.1185;
        mesh[158].p2.z = 0.120599;
        mesh[158].p3.x = -0.057323;
        mesh[158].p3.y = -0.147545;
        mesh[158].p3.z = 0.177989;
        mesh[159].p1.x = -0.021311;
        mesh[159].p1.y = -0.109992;
        mesh[159].p1.z = 0.007536;
        mesh[159].p2.x = -0.025373;
        mesh[159].p2.y = -0.09533;
        mesh[159].p2.z = 0.018332;
        mesh[159].p3.x = -0.031352;
        mesh[159].p3.y = -0.118191;
        mesh[159].p3.z = 0.050941;
        mesh[160].p1.x = -0.044853;
        mesh[160].p1.y = -0.189594;
        mesh[160].p1.z = 0.271944;
        mesh[160].p2.x = -0.066432;
        mesh[160].p2.y = -0.16482;
        mesh[160].p2.z = 0.29205;
        mesh[160].p3.x = -0.063422;
        mesh[160].p3.y = -0.164003;
        mesh[160].p3.z = 0.271249;
        mesh[161].p1.x = -0.043043;
        mesh[161].p1.y = -0.1738;
        mesh[161].p1.z = 0.335948;
        mesh[161].p2.x = -0.066432;
        mesh[161].p2.y = -0.16482;
        mesh[161].p2.z = 0.29205;
        mesh[161].p3.x = -0.044853;
        mesh[161].p3.y = -0.189594;
        mesh[161].p3.z = 0.271944;
        mesh[162].p1.x = -0.052194;
        mesh[162].p1.y = 0.037269;
        mesh[162].p1.z = 0.199662;
        mesh[162].p2.x = -0.04672;
        mesh[162].p2.y = 0.093756;
        mesh[162].p2.z = 0.235686;
        mesh[162].p3.x = -0.051287;
        mesh[162].p3.y = 0.123923;
        mesh[162].p3.z = 0.200622;
        mesh[163].p1.x = -0.023133;
        mesh[163].p1.y = 0.120981;
        mesh[163].p1.z = 0.078476;
        mesh[163].p2.x = -0.031246;
        mesh[163].p2.y = 0.135959;
        mesh[163].p2.z = 0.153033;
        mesh[163].p3.x = -0.009247;
        mesh[163].p3.y = 0.106663;
        mesh[163].p3.z = 0.145238;
        mesh[164].p1.x = -0.037203;
        mesh[164].p1.y = 0.138321;
        mesh[164].p1.z = 0.100764;
        mesh[164].p2.x = -0.023133;
        mesh[164].p2.y = 0.120981;
        mesh[164].p2.z = 0.078476;
        mesh[164].p3.x = -0.037603;
        mesh[164].p3.y = 0.15207;
        mesh[164].p3.z = 0.084055;
        mesh[165].p1.x = -0.045197;
        mesh[165].p1.y = 0.121199;
        mesh[165].p1.z = 0.055852;
        mesh[165].p2.x = -0.030099;
        mesh[165].p2.y = 0.074776;
        mesh[165].p2.z = 0.120707;
        mesh[165].p3.x = -0.047097;
        mesh[165].p3.y = 0.052817;
        mesh[165].p3.z = 0.162408;
        mesh[166].p1.x = -0.05208;
        mesh[166].p1.y = 0.124926;
        mesh[166].p1.z = 0.07032;
        mesh[166].p2.x = -0.037203;
        mesh[166].p2.y = 0.138321;
        mesh[166].p2.z = 0.100764;
        mesh[166].p3.x = -0.037603;
        mesh[166].p3.y = 0.15207;
        mesh[166].p3.z = 0.084055;
        mesh[167].p1.x = -0.012945;
        mesh[167].p1.y = 0.045589;
        mesh[167].p1.z = 0.13692;
        mesh[167].p2.x = -0.027136;
        mesh[167].p2.y = 0.041664;
        mesh[167].p2.z = 0.132141;
        mesh[167].p3.x = -0.030099;
        mesh[167].p3.y = 0.074776;
        mesh[167].p3.z = 0.120707;
        mesh[168].p1.x = 0.009345;
        mesh[168].p1.y = 0.108126;
        mesh[168].p1.z = 0.144102;
        mesh[168].p2.x = -0.009247;
        mesh[168].p2.y = 0.106663;
        mesh[168].p2.z = 0.145238;
        mesh[168].p3.x = -3.26E-4;
        mesh[168].p3.y = 0.141988;
        mesh[168].p3.z = 0.196277;
        mesh[169].p1.x = -0.050702;
        mesh[169].p1.y = -0.129375;
        mesh[169].p1.z = 6.69E-4;
        mesh[169].p2.x = -0.021311;
        mesh[169].p2.y = -0.109992;
        mesh[169].p2.z = 0.007536;
        mesh[169].p3.x = -0.04676;
        mesh[169].p3.y = -0.146333;
        mesh[169].p3.z = 0.006207;
        mesh[170].p1.x = -0.050702;
        mesh[170].p1.y = -0.129375;
        mesh[170].p1.z = 6.69E-4;
        mesh[170].p2.x = -0.061288;
        mesh[170].p2.y = -0.130257;
        mesh[170].p2.z = 0.005966;
        mesh[170].p3.x = -0.048878;
        mesh[170].p3.y = -0.099334;
        mesh[170].p3.z = 0.010049;
        mesh[171].p1.x = -0.043016;
        mesh[171].p1.y = 0.11447;
        mesh[171].p1.z = -4.97E-4;
        mesh[171].p2.x = -0.023356;
        mesh[171].p2.y = 0.126501;
        mesh[171].p2.z = 0.002503;
        mesh[171].p3.x = -0.03108;
        mesh[171].p3.y = 0.104637;
        mesh[171].p3.z = 7.6E-5;
        mesh[172].p1.x = -0.052415;
        mesh[172].p1.y = 0.095648;
        mesh[172].p1.z = -3.5E-5;
        mesh[172].p2.x = -0.055813;
        mesh[172].p2.y = 0.124066;
        mesh[172].p2.z = 7.64E-4;
        mesh[172].p3.x = -0.043016;
        mesh[172].p3.y = 0.11447;
        mesh[172].p3.z = -4.97E-4;
        mesh[173].p1.x = -0.03108;
        mesh[173].p1.y = 0.104637;
        mesh[173].p1.z = 7.6E-5;
        mesh[173].p2.x = -0.022864;
        mesh[173].p2.y = 0.094117;
        mesh[173].p2.z = -3.86E-4;
        mesh[173].p3.x = -0.031653;
        mesh[173].p3.y = 0.105782;
        mesh[173].p3.z = 0.018291;
        mesh[174].p1.x = -0.039606;
        mesh[174].p1.y = 0.087492;
        mesh[174].p1.z = 0.004376;
        mesh[174].p2.x = -0.03108;
        mesh[174].p2.y = 0.104637;
        mesh[174].p2.z = 7.6E-5;
        mesh[174].p3.x = -0.033126;
        mesh[174].p3.y = 0.084624;
        mesh[174].p3.z = 0.003634;
        mesh[175].p1.x = 0.007313;
        mesh[175].p1.y = -0.265396;
        mesh[175].p1.z = 0.327658;
        mesh[175].p2.x = 0.002083;
        mesh[175].p2.y = -0.267639;
        mesh[175].p2.z = 0.312103;
        mesh[175].p3.x = 0.007242;
        mesh[175].p3.y = -0.264548;
        mesh[175].p3.z = 0.31558;
        mesh[176].p1.x = 0.00371;
        mesh[176].p1.y = -0.265544;
        mesh[176].p1.z = 0.295988;
        mesh[176].p2.x = 0.023713;
        mesh[176].p2.y = -0.253701;
        mesh[176].p2.z = 0.282565;
        mesh[176].p3.x = 0.024291;
        mesh[176].p3.y = -0.263459;
        mesh[176].p3.z = 0.302255;
        mesh[177].p1.x = 0.002083;
        mesh[177].p1.y = -0.267639;
        mesh[177].p1.z = 0.312103;
        mesh[177].p2.x = 0.00371;
        mesh[177].p2.y = -0.265544;
        mesh[177].p2.z = 0.295988;
        mesh[177].p3.x = 0.024291;
        mesh[177].p3.y = -0.263459;
        mesh[177].p3.z = 0.302255;
        mesh[178].p1.x = 0.023642;
        mesh[178].p1.y = -0.252005;
        mesh[178].p1.z = 0.316736;
        mesh[178].p2.x = 0.007242;
        mesh[178].p2.y = -0.264548;
        mesh[178].p2.z = 0.31558;
        mesh[178].p3.x = 0.024291;
        mesh[178].p3.y = -0.263459;
        mesh[178].p3.z = 0.302255;
        mesh[179].p1.x = -0.009896;
        mesh[179].p1.y = -0.227817;
        mesh[179].p1.z = 0.334942;
        mesh[179].p2.x = 0.011163;
        mesh[179].p2.y = -0.226997;
        mesh[179].p2.z = 0.334062;
        mesh[179].p3.x = 0.020438;
        mesh[179].p3.y = -0.217917;
        mesh[179].p3.z = 0.346311;
        mesh[180].p1.x = 0.023713;
        mesh[180].p1.y = -0.253701;
        mesh[180].p1.z = 0.282565;
        mesh[180].p2.x = 0.023642;
        mesh[180].p2.y = -0.252005;
        mesh[180].p2.z = 0.316736;
        mesh[180].p3.x = 0.024291;
        mesh[180].p3.y = -0.263459;
        mesh[180].p3.z = 0.302255;
        mesh[181].p1.x = 0.023642;
        mesh[181].p1.y = -0.252005;
        mesh[181].p1.z = 0.316736;
        mesh[181].p2.x = 0.023713;
        mesh[181].p2.y = -0.253701;
        mesh[181].p2.z = 0.282565;
        mesh[181].p3.x = 0.031688;
        mesh[181].p3.y = -0.226472;
        mesh[181].p3.z = 0.28078;
        mesh[182].p1.x = -0.017193;
        mesh[182].p1.y = -0.210889;
        mesh[182].p1.z = 0.257864;
        mesh[182].p2.x = -0.00563;
        mesh[182].p2.y = -0.191954;
        mesh[182].p2.z = 0.225786;
        mesh[182].p3.x = 0.01869;
        mesh[182].p3.y = -0.214564;
        mesh[182].p3.z = 0.260349;
        mesh[183].p1.x = 0.043042;
        mesh[183].p1.y = -0.1738;
        mesh[183].p1.z = 0.335948;
        mesh[183].p2.x = 0.020438;
        mesh[183].p2.y = -0.217917;
        mesh[183].p2.z = 0.346311;
        mesh[183].p3.x = 0.037105;
        mesh[183].p3.y = -0.210377;
        mesh[183].p3.z = 0.330419;
        mesh[184].p1.x = 0.038391;
        mesh[184].p1.y = -0.156011;
        mesh[184].p1.z = 0.280595;
        mesh[184].p2.x = 0.031688;
        mesh[184].p2.y = -0.226472;
        mesh[184].p2.z = 0.28078;
        mesh[184].p3.x = 0.033794;
        mesh[184].p3.y = -0.169259;
        mesh[184].p3.z = 0.248461;
        mesh[185].p1.x = 0.043042;
        mesh[185].p1.y = -0.1738;
        mesh[185].p1.z = 0.335948;
        mesh[185].p2.x = 0.037105;
        mesh[185].p2.y = -0.210377;
        mesh[185].p2.z = 0.330419;
        mesh[185].p3.x = 0.048469;
        mesh[185].p3.y = -0.158545;
        mesh[185].p3.z = 0.319398;
        mesh[186].p1.x = 0.035308;
        mesh[186].p1.y = -0.137464;
        mesh[186].p1.z = 0.306814;
        mesh[186].p2.x = 0.023841;
        mesh[186].p2.y = -0.138332;
        mesh[186].p2.z = 0.326807;
        mesh[186].p3.x = 0.053788;
        mesh[186].p3.y = -0.129921;
        mesh[186].p3.z = 0.303985;
        mesh[187].p1.x = 0.053788;
        mesh[187].p1.y = -0.129921;
        mesh[187].p1.z = 0.303985;
        mesh[187].p2.x = 0.023841;
        mesh[187].p2.y = -0.138332;
        mesh[187].p2.z = 0.326807;
        mesh[187].p3.x = 0.049253;
        mesh[187].p3.y = -0.135547;
        mesh[187].p3.z = 0.325422;
        mesh[188].p1.x = 0.049253;
        mesh[188].p1.y = -0.135547;
        mesh[188].p1.z = 0.325422;
        mesh[188].p2.x = 0.066431;
        mesh[188].p2.y = -0.16482;
        mesh[188].p2.z = 0.292051;
        mesh[188].p3.x = 0.053788;
        mesh[188].p3.y = -0.129921;
        mesh[188].p3.z = 0.303985;
        mesh[189].p1.x = 0.043042;
        mesh[189].p1.y = -0.1738;
        mesh[189].p1.z = 0.335948;
        mesh[189].p2.x = 0.048469;
        mesh[189].p2.y = -0.158545;
        mesh[189].p2.z = 0.319398;
        mesh[189].p3.x = 0.063421;
        mesh[189].p3.y = -0.164003;
        mesh[189].p3.z = 0.271249;
        mesh[190].p1.x = 0.063421;
        mesh[190].p1.y = -0.164003;
        mesh[190].p1.z = 0.271249;
        mesh[190].p2.x = 0.044853;
        mesh[190].p2.y = -0.189593;
        mesh[190].p2.z = 0.271944;
        mesh[190].p3.x = 0.043042;
        mesh[190].p3.y = -0.1738;
        mesh[190].p3.z = 0.335948;
        mesh[191].p1.x = 0.039196;
        mesh[191].p1.y = -0.156989;
        mesh[191].p1.z = 0.339227;
        mesh[191].p2.x = 0.023841;
        mesh[191].p2.y = -0.138332;
        mesh[191].p2.z = 0.326807;
        mesh[191].p3.x = -7.79E-4;
        mesh[191].p3.y = -0.146982;
        mesh[191].p3.z = 0.344462;
        mesh[192].p1.x = 0.020438;
        mesh[192].p1.y = -0.217917;
        mesh[192].p1.z = 0.346311;
        mesh[192].p2.x = 0.043042;
        mesh[192].p2.y = -0.1738;
        mesh[192].p2.z = 0.335948;
        mesh[192].p3.x = 0.0;
        mesh[192].p3.y = -0.195958;
        mesh[192].p3.z = 0.354306;
        mesh[193].p1.x = 0.011163;
        mesh[193].p1.y = -0.226997;
        mesh[193].p1.z = 0.334062;
        mesh[193].p2.x = 0.037105;
        mesh[193].p2.y = -0.210377;
        mesh[193].p2.z = 0.330419;
        mesh[193].p3.x = 0.020438;
        mesh[193].p3.y = -0.217917;
        mesh[193].p3.z = 0.346311;
        mesh[194].p1.x = 0.01869;
        mesh[194].p1.y = -0.214564;
        mesh[194].p1.z = 0.260349;
        mesh[194].p2.x = 0.020072;
        mesh[194].p2.y = -0.245234;
        mesh[194].p2.z = 0.286429;
        mesh[194].p3.x = 5.65E-4;
        mesh[194].p3.y = -0.260455;
        mesh[194].p3.z = 0.278232;
        mesh[195].p1.x = 0.020072;
        mesh[195].p1.y = -0.245234;
        mesh[195].p1.z = 0.286429;
        mesh[195].p2.x = 0.00371;
        mesh[195].p2.y = -0.265544;
        mesh[195].p2.z = 0.295988;
        mesh[195].p3.x = 5.65E-4;
        mesh[195].p3.y = -0.260455;
        mesh[195].p3.z = 0.278232;
        mesh[196].p1.x = 0.023841;
        mesh[196].p1.y = -0.138332;
        mesh[196].p1.z = 0.326807;
        mesh[196].p2.x = -0.006761;
        mesh[196].p2.y = -0.10839;
        mesh[196].p2.z = 0.287993;
        mesh[196].p3.x = -0.022417;
        mesh[196].p3.y = -0.136843;
        mesh[196].p3.z = 0.326557;
        mesh[197].p1.x = 0.01869;
        mesh[197].p1.y = -0.214564;
        mesh[197].p1.z = 0.260349;
        mesh[197].p2.x = 0.016642;
        mesh[197].p2.y = -0.186885;
        mesh[197].p2.z = 0.224214;
        mesh[197].p3.x = 0.033794;
        mesh[197].p3.y = -0.169259;
        mesh[197].p3.z = 0.248461;
        mesh[198].p1.x = 0.020072;
        mesh[198].p1.y = -0.245234;
        mesh[198].p1.z = 0.286429;
        mesh[198].p2.x = 0.031688;
        mesh[198].p2.y = -0.226472;
        mesh[198].p2.z = 0.28078;
        mesh[198].p3.x = 0.023713;
        mesh[198].p3.y = -0.253701;
        mesh[198].p3.z = 0.282565;
        mesh[199].p1.x = 0.039372;
        mesh[199].p1.y = 0.091188;
        mesh[199].p1.z = 0.014932;
        mesh[199].p2.x = 0.036214;
        mesh[199].p2.y = 0.086871;
        mesh[199].p2.z = 0.00329;
        mesh[199].p3.x = 0.049285;
        mesh[199].p3.y = 0.093207;
        mesh[199].p3.z = 1.43E-4;
        mesh[200].p1.x = 0.039372;
        mesh[200].p1.y = 0.091188;
        mesh[200].p1.z = 0.014932;
        mesh[200].p2.x = 0.051959;
        mesh[200].p2.y = 0.116839;
        mesh[200].p2.z = 0.027809;
        mesh[200].p3.x = 0.035054;
        mesh[200].p3.y = 0.113693;
        mesh[200].p3.z = 0.02676;
        mesh[201].p1.x = 0.044295;
        mesh[201].p1.y = -0.086567;
        mesh[201].p1.z = 0.047665;
        mesh[201].p2.x = 0.049755;
        mesh[201].p2.y = -0.105056;
        mesh[201].p2.z = 0.003318;
        mesh[201].p3.x = 0.024485;
        mesh[201].p3.y = -0.097618;
        mesh[201].p3.z = 0.011241;
        mesh[202].p1.x = 0.060234;
        mesh[202].p1.y = -0.111718;
        mesh[202].p1.z = 0.141464;
        mesh[202].p2.x = 0.037948;
        mesh[202].p2.y = -0.084855;
        mesh[202].p2.z = 0.13526;
        mesh[202].p3.x = 0.052896;
        mesh[202].p3.y = -0.077728;
        mesh[202].p3.z = 0.169219;
        mesh[203].p1.x = 0.060234;
        mesh[203].p1.y = -0.111718;
        mesh[203].p1.z = 0.141464;
        mesh[203].p2.x = 0.048408;
        mesh[203].p2.y = -0.136161;
        mesh[203].p2.z = 0.125389;
        mesh[203].p3.x = 0.051727;
        mesh[203].p3.y = -0.106852;
        mesh[203].p3.z = 0.0213;
        mesh[204].p1.x = 0.051727;
        mesh[204].p1.y = -0.106852;
        mesh[204].p1.z = 0.0213;
        mesh[204].p2.x = 0.048408;
        mesh[204].p2.y = -0.136161;
        mesh[204].p2.z = 0.125389;
        mesh[204].p3.x = 0.041861;
        mesh[204].p3.y = -0.117404;
        mesh[204].p3.z = 0.027621;
        mesh[205].p1.x = 0.044295;
        mesh[205].p1.y = -0.086567;
        mesh[205].p1.z = 0.047665;
        mesh[205].p2.x = 0.060234;
        mesh[205].p2.y = -0.111718;
        mesh[205].p2.z = 0.141464;
        mesh[205].p3.x = 0.051727;
        mesh[205].p3.y = -0.106852;
        mesh[205].p3.z = 0.0213;
        mesh[206].p1.x = 0.025447;
        mesh[206].p1.y = -0.16317;
        mesh[206].p1.z = 0.14662;
        mesh[206].p2.x = 0.048408;
        mesh[206].p2.y = -0.136161;
        mesh[206].p2.z = 0.125389;
        mesh[206].p3.x = 0.055679;
        mesh[206].p3.y = -0.148737;
        mesh[206].p3.z = 0.164846;
        mesh[207].p1.x = 0.052896;
        mesh[207].p1.y = -0.077728;
        mesh[207].p1.z = 0.169219;
        mesh[207].p2.x = 0.037948;
        mesh[207].p2.y = -0.084855;
        mesh[207].p2.z = 0.13526;
        mesh[207].p3.x = 0.020677;
        mesh[207].p3.y = 0.005463;
        mesh[207].p3.z = 0.12829;
        mesh[208].p1.x = 0.01653;
        mesh[208].p1.y = -0.109257;
        mesh[208].p1.z = 0.120113;
        mesh[208].p2.x = -0.027136;
        mesh[208].p2.y = 0.041664;
        mesh[208].p2.z = 0.132141;
        mesh[208].p3.x = 0.020677;
        mesh[208].p3.y = 0.005463;
        mesh[208].p3.z = 0.12829;
        mesh[209].p1.x = 0.032255;
        mesh[209].p1.y = -0.08755;
        mesh[209].p1.z = 0.048576;
        mesh[209].p2.x = 0.024485;
        mesh[209].p2.y = -0.097618;
        mesh[209].p2.z = 0.011241;
        mesh[209].p3.x = 0.024929;
        mesh[209].p3.y = -0.10802;
        mesh[209].p3.z = 0.027315;
        mesh[210].p1.x = 0.049755;
        mesh[210].p1.y = -0.105056;
        mesh[210].p1.z = 0.003318;
        mesh[210].p2.x = 0.044295;
        mesh[210].p2.y = -0.086567;
        mesh[210].p2.z = 0.047665;
        mesh[210].p3.x = 0.051727;
        mesh[210].p3.y = -0.106852;
        mesh[210].p3.z = 0.0213;
        mesh[211].p1.x = 0.040152;
        mesh[211].p1.y = -0.135693;
        mesh[211].p1.z = 0.020709;
        mesh[211].p2.x = 0.024813;
        mesh[211].p2.y = -0.137638;
        mesh[211].p2.z = 0.006129;
        mesh[211].p3.x = 0.031585;
        mesh[211].p3.y = -0.136678;
        mesh[211].p3.z = 0.002046;
        mesh[212].p1.x = 0.056245;
        mesh[212].p1.y = -0.133422;
        mesh[212].p1.z = 0.00208;
        mesh[212].p2.x = 0.049755;
        mesh[212].p2.y = -0.105056;
        mesh[212].p2.z = 0.003318;
        mesh[212].p3.x = 0.058952;
        mesh[212].p3.y = -0.120225;
        mesh[212].p3.z = 0.013241;
        mesh[213].p1.x = 0.017448;
        mesh[213].p1.y = 0.138736;
        mesh[213].p1.z = 0.227783;
        mesh[213].p2.x = 0.004478;
        mesh[213].p2.y = 0.178564;
        mesh[213].p2.z = 0.232537;
        mesh[213].p3.x = -0.014487;
        mesh[213].p3.y = 0.080367;
        mesh[213].p3.z = 0.261155;
        mesh[214].p1.x = 0.029122;
        mesh[214].p1.y = 0.138151;
        mesh[214].p1.z = 0.027991;
        mesh[214].p2.x = 0.023355;
        mesh[214].p2.y = 0.127603;
        mesh[214].p2.z = 0.002924;
        mesh[214].p3.x = 0.026391;
        mesh[214].p3.y = 0.117277;
        mesh[214].p3.z = 0.021774;
        mesh[215].p1.x = 0.046396;
        mesh[215].p1.y = 0.140708;
        mesh[215].p1.z = 0.081159;
        mesh[215].p2.x = 0.055369;
        mesh[215].p2.y = 0.132142;
        mesh[215].p2.z = 0.02032;
        mesh[215].p3.x = 0.037267;
        mesh[215].p3.y = 0.142487;
        mesh[215].p3.z = 0.021146;
        mesh[216].p1.x = 0.025813;
        mesh[216].p1.y = 0.096198;
        mesh[216].p1.z = -1.51E-4;
        mesh[216].p2.x = 0.036214;
        mesh[216].p2.y = 0.086871;
        mesh[216].p2.z = 0.00329;
        mesh[216].p3.x = 0.035054;
        mesh[216].p3.y = 0.113693;
        mesh[216].p3.z = 0.02676;
        mesh[217].p1.x = 0.037267;
        mesh[217].p1.y = 0.142487;
        mesh[217].p1.z = 0.021146;
        mesh[217].p2.x = 0.053684;
        mesh[217].p2.y = 0.125516;
        mesh[217].p2.z = 0.002388;
        mesh[217].p3.x = 0.023355;
        mesh[217].p3.y = 0.127603;
        mesh[217].p3.z = 0.002924;
        mesh[218].p1.x = 0.057539;
        mesh[218].p1.y = 0.09589;
        mesh[218].p1.z = 3.54E-4;
        mesh[218].p2.x = 0.055369;
        mesh[218].p2.y = 0.132142;
        mesh[218].p2.z = 0.02032;
        mesh[218].p3.x = 0.051959;
        mesh[218].p3.y = 0.116839;
        mesh[218].p3.z = 0.027809;
        mesh[219].p1.x = 0.057539;
        mesh[219].p1.y = 0.09589;
        mesh[219].p1.z = 3.54E-4;
        mesh[219].p2.x = 0.051959;
        mesh[219].p2.y = 0.116839;
        mesh[219].p2.z = 0.027809;
        mesh[219].p3.x = 0.049285;
        mesh[219].p3.y = 0.093207;
        mesh[219].p3.z = 1.43E-4;
        mesh[220].p1.x = 0.049285;
        mesh[220].p1.y = 0.093207;
        mesh[220].p1.z = 1.43E-4;
        mesh[220].p2.x = 0.051959;
        mesh[220].p2.y = 0.116839;
        mesh[220].p2.z = 0.027809;
        mesh[220].p3.x = 0.039372;
        mesh[220].p3.y = 0.091188;
        mesh[220].p3.z = 0.014932;
        mesh[221].p1.x = 0.028943;
        mesh[221].p1.y = 0.117552;
        mesh[221].p1.z = 0.068098;
        mesh[221].p2.x = 0.051918;
        mesh[221].p2.y = 0.123289;
        mesh[221].p2.z = 0.067837;
        mesh[221].p3.x = 0.035775;
        mesh[221].p3.y = 0.072902;
        mesh[221].p3.z = 0.121575;
        mesh[222].p1.x = 0.024929;
        mesh[222].p1.y = -0.10802;
        mesh[222].p1.z = 0.027315;
        mesh[222].p2.x = 0.041861;
        mesh[222].p2.y = -0.117404;
        mesh[222].p2.z = 0.027621;
        mesh[222].p3.x = 0.048408;
        mesh[222].p3.y = -0.136161;
        mesh[222].p3.z = 0.125389;
        mesh[223].p1.x = 0.024813;
        mesh[223].p1.y = -0.137638;
        mesh[223].p1.z = 0.006129;
        mesh[223].p2.x = 0.040152;
        mesh[223].p2.y = -0.135693;
        mesh[223].p2.z = 0.020709;
        mesh[223].p3.x = 0.024929;
        mesh[223].p3.y = -0.10802;
        mesh[223].p3.z = 0.027315;
        mesh[224].p1.x = 0.028943;
        mesh[224].p1.y = 0.117552;
        mesh[224].p1.z = 0.068098;
        mesh[224].p2.x = 0.035054;
        mesh[224].p2.y = 0.113693;
        mesh[224].p2.z = 0.02676;
        mesh[224].p3.x = 0.051918;
        mesh[224].p3.y = 0.123289;
        mesh[224].p3.z = 0.067837;
        mesh[225].p1.x = 0.026391;
        mesh[225].p1.y = 0.117277;
        mesh[225].p1.z = 0.021774;
        mesh[225].p2.x = 0.035054;
        mesh[225].p2.y = 0.113693;
        mesh[225].p2.z = 0.02676;
        mesh[225].p3.x = 0.028943;
        mesh[225].p3.y = 0.117552;
        mesh[225].p3.z = 0.068098;
        mesh[226].p1.x = 0.037183;
        mesh[226].p1.y = -0.008779;
        mesh[226].p1.z = 0.249569;
        mesh[226].p2.x = 0.003674;
        mesh[226].p2.y = -0.085117;
        mesh[226].p2.z = 0.272056;
        mesh[226].p3.x = 0.021506;
        mesh[226].p3.y = -0.10025;
        mesh[226].p3.z = 0.266035;
        mesh[227].p1.x = -0.009247;
        mesh[227].p1.y = 0.106663;
        mesh[227].p1.z = 0.145238;
        mesh[227].p2.x = 0.009345;
        mesh[227].p2.y = 0.108126;
        mesh[227].p2.z = 0.144102;
        mesh[227].p3.x = 0.011832;
        mesh[227].p3.y = 0.073182;
        mesh[227].p3.z = 0.136522;
        mesh[228].p1.x = 0.005073;
        mesh[228].p1.y = 0.031262;
        mesh[228].p1.z = 0.135312;
        mesh[228].p2.x = -3.31E-4;
        mesh[228].p2.y = 0.036278;
        mesh[228].p2.z = 0.123064;
        mesh[228].p3.x = 0.011832;
        mesh[228].p3.y = 0.073182;
        mesh[228].p3.z = 0.136522;
        mesh[229].p1.x = 0.002484;
        mesh[229].p1.y = 0.265509;
        mesh[229].p1.z = 0.210879;
        mesh[229].p2.x = 0.017448;
        mesh[229].p2.y = 0.138736;
        mesh[229].p2.z = 0.227783;
        mesh[229].p3.x = -1.45E-4;
        mesh[229].p3.y = 0.265834;
        mesh[229].p3.z = 0.193142;
        mesh[230].p1.x = 0.032381;
        mesh[230].p1.y = 0.136286;
        mesh[230].p1.z = 0.15941;
        mesh[230].p2.x = -3.26E-4;
        mesh[230].p2.y = 0.141988;
        mesh[230].p2.z = 0.196277;
        mesh[230].p3.x = 0.047319;
        mesh[230].p3.y = 0.12642;
        mesh[230].p3.z = 0.207304;
        mesh[231].p1.x = 0.033794;
        mesh[231].p1.y = -0.169259;
        mesh[231].p1.z = 0.248461;
        mesh[231].p2.x = 0.040651;
        mesh[231].p2.y = -0.136387;
        mesh[231].p2.z = 0.244256;
        mesh[231].p3.x = 0.038391;
        mesh[231].p3.y = -0.156011;
        mesh[231].p3.z = 0.280595;
        mesh[232].p1.x = 0.024813;
        mesh[232].p1.y = -0.137638;
        mesh[232].p1.z = 0.006129;
        mesh[232].p2.x = 0.024485;
        mesh[232].p2.y = -0.097618;
        mesh[232].p2.z = 0.011241;
        mesh[232].p3.x = 0.031585;
        mesh[232].p3.y = -0.136678;
        mesh[232].p3.z = 0.002046;
        mesh[233].p1.x = 0.036214;
        mesh[233].p1.y = 0.086871;
        mesh[233].p1.z = 0.00329;
        mesh[233].p2.x = 0.025813;
        mesh[233].p2.y = 0.096198;
        mesh[233].p2.z = -1.51E-4;
        mesh[233].p3.x = 0.038873;
        mesh[233].p3.y = 0.105471;
        mesh[233].p3.z = 5.0E-5;
        mesh[234].p1.x = 0.030553;
        mesh[234].p1.y = 0.149125;
        mesh[234].p1.z = 0.079185;
        mesh[234].p2.x = 0.011832;
        mesh[234].p2.y = 0.073182;
        mesh[234].p2.z = 0.136522;
        mesh[234].p3.x = 0.021511;
        mesh[234].p3.y = 0.123621;
        mesh[234].p3.z = 0.140422;
        mesh[235].p1.x = 0.037948;
        mesh[235].p1.y = -0.084855;
        mesh[235].p1.z = 0.13526;
        mesh[235].p2.x = 0.047802;
        mesh[235].p2.y = -0.086047;
        mesh[235].p2.z = 0.120183;
        mesh[235].p3.x = 0.01653;
        mesh[235].p3.y = -0.109257;
        mesh[235].p3.z = 0.120113;
        mesh[236].p1.x = 0.037948;
        mesh[236].p1.y = -0.084855;
        mesh[236].p1.z = 0.13526;
        mesh[236].p2.x = 0.060234;
        mesh[236].p2.y = -0.111718;
        mesh[236].p2.z = 0.141464;
        mesh[236].p3.x = 0.047802;
        mesh[236].p3.y = -0.086047;
        mesh[236].p3.z = 0.120183;
        mesh[237].p1.x = 0.038391;
        mesh[237].p1.y = -0.156011;
        mesh[237].p1.z = 0.280595;
        mesh[237].p2.x = 0.040651;
        mesh[237].p2.y = -0.136387;
        mesh[237].p2.z = 0.244256;
        mesh[237].p3.x = 0.035308;
        mesh[237].p3.y = -0.137464;
        mesh[237].p3.z = 0.306814;
        mesh[238].p1.x = -3.26E-4;
        mesh[238].p1.y = 0.141988;
        mesh[238].p1.z = 0.196277;
        mesh[238].p2.x = -0.003326;
        mesh[238].p2.y = 0.138226;
        mesh[238].p2.z = 0.212925;
        mesh[238].p3.x = 0.017448;
        mesh[238].p3.y = 0.138736;
        mesh[238].p3.z = 0.227783;
        mesh[239].p1.x = 0.003674;
        mesh[239].p1.y = -0.085117;
        mesh[239].p1.z = 0.272056;
        mesh[239].p2.x = 0.037183;
        mesh[239].p2.y = -0.008779;
        mesh[239].p2.z = 0.249569;
        mesh[239].p3.x = -0.010543;
        mesh[239].p3.y = -0.015599;
        mesh[239].p3.z = 0.262292;
        mesh[240].p1.x = 0.023841;
        mesh[240].p1.y = -0.138332;
        mesh[240].p1.z = 0.326807;
        mesh[240].p2.x = 0.035308;
        mesh[240].p2.y = -0.137464;
        mesh[240].p2.z = 0.306814;
        mesh[240].p3.x = 0.021506;
        mesh[240].p3.y = -0.10025;
        mesh[240].p3.z = 0.266035;
        mesh[241].p1.x = 0.021511;
        mesh[241].p1.y = 0.123621;
        mesh[241].p1.z = 0.140422;
        mesh[241].p2.x = 0.009345;
        mesh[241].p2.y = 0.108126;
        mesh[241].p2.z = 0.144102;
        mesh[241].p3.x = 0.032381;
        mesh[241].p3.y = 0.136286;
        mesh[241].p3.z = 0.15941;
        mesh[242].p1.x = 0.055855;
        mesh[242].p1.y = -0.084068;
        mesh[242].p1.z = 0.206511;
        mesh[242].p2.x = 0.052896;
        mesh[242].p2.y = -0.077728;
        mesh[242].p2.z = 0.169219;
        mesh[242].p3.x = 0.051207;
        mesh[242].p3.y = 0.016414;
        mesh[242].p3.z = 0.221731;
        mesh[243].p1.x = 0.051207;
        mesh[243].p1.y = 0.016414;
        mesh[243].p1.z = 0.221731;
        mesh[243].p2.x = 0.037183;
        mesh[243].p2.y = -0.008779;
        mesh[243].p2.z = 0.249569;
        mesh[243].p3.x = 0.055855;
        mesh[243].p3.y = -0.084068;
        mesh[243].p3.z = 0.206511;
        mesh[244].p1.x = 0.021506;
        mesh[244].p1.y = -0.10025;
        mesh[244].p1.z = 0.266035;
        mesh[244].p2.x = 0.040651;
        mesh[244].p2.y = -0.136387;
        mesh[244].p2.z = 0.244256;
        mesh[244].p3.x = 0.055855;
        mesh[244].p3.y = -0.084068;
        mesh[244].p3.z = 0.206511;
        mesh[245].p1.x = 0.033794;
        mesh[245].p1.y = -0.169259;
        mesh[245].p1.z = 0.248461;
        mesh[245].p2.x = 0.059678;
        mesh[245].p2.y = -0.129737;
        mesh[245].p2.z = 0.199429;
        mesh[245].p3.x = 0.040651;
        mesh[245].p3.y = -0.136387;
        mesh[245].p3.z = 0.244256;
        mesh[246].p1.x = 0.025447;
        mesh[246].p1.y = -0.16317;
        mesh[246].p1.z = 0.14662;
        mesh[246].p2.x = 0.016642;
        mesh[246].p2.y = -0.186885;
        mesh[246].p2.z = 0.224214;
        mesh[246].p3.x = -0.00563;
        mesh[246].p3.y = -0.191954;
        mesh[246].p3.z = 0.225786;
        mesh[247].p1.x = 0.060234;
        mesh[247].p1.y = -0.111718;
        mesh[247].p1.z = 0.141464;
        mesh[247].p2.x = 0.059678;
        mesh[247].p2.y = -0.129737;
        mesh[247].p2.z = 0.199429;
        mesh[247].p3.x = 0.055679;
        mesh[247].p3.y = -0.148737;
        mesh[247].p3.z = 0.164846;
        mesh[248].p1.x = 0.024219;
        mesh[248].p1.y = -0.128092;
        mesh[248].p1.z = 0.121928;
        mesh[248].p2.x = 0.048408;
        mesh[248].p2.y = -0.136161;
        mesh[248].p2.z = 0.125389;
        mesh[248].p3.x = 0.01653;
        mesh[248].p3.y = -0.109257;
        mesh[248].p3.z = 0.120113;
        mesh[249].p1.x = 0.048408;
        mesh[249].p1.y = -0.136161;
        mesh[249].p1.z = 0.125389;
        mesh[249].p2.x = 0.025447;
        mesh[249].p2.y = -0.16317;
        mesh[249].p2.z = 0.14662;
        mesh[249].p3.x = 0.01653;
        mesh[249].p3.y = -0.109257;
        mesh[249].p3.z = 0.120113;
        mesh[250].p1.x = 0.044853;
        mesh[250].p1.y = -0.189593;
        mesh[250].p1.z = 0.271944;
        mesh[250].p2.x = 0.063421;
        mesh[250].p2.y = -0.164003;
        mesh[250].p2.z = 0.271249;
        mesh[250].p3.x = 0.066431;
        mesh[250].p3.y = -0.16482;
        mesh[250].p3.z = 0.292051;
        mesh[251].p1.x = 0.048469;
        mesh[251].p1.y = -0.158545;
        mesh[251].p1.z = 0.319398;
        mesh[251].p2.x = 0.053788;
        mesh[251].p2.y = -0.129921;
        mesh[251].p2.z = 0.303985;
        mesh[251].p3.x = 0.063421;
        mesh[251].p3.y = -0.164003;
        mesh[251].p3.z = 0.271249;
        mesh[252].p1.x = 0.043042;
        mesh[252].p1.y = -0.1738;
        mesh[252].p1.z = 0.335948;
        mesh[252].p2.x = 0.044853;
        mesh[252].p2.y = -0.189593;
        mesh[252].p2.z = 0.271944;
        mesh[252].p3.x = 0.066431;
        mesh[252].p3.y = -0.16482;
        mesh[252].p3.z = 0.292051;
        mesh[253].p1.x = 0.051207;
        mesh[253].p1.y = 0.016414;
        mesh[253].p1.z = 0.221731;
        mesh[253].p2.x = 0.0337;
        mesh[253].p2.y = 0.048755;
        mesh[253].p2.z = 0.136203;
        mesh[253].p3.x = 0.035775;
        mesh[253].p3.y = 0.072902;
        mesh[253].p3.z = 0.121575;
        mesh[254].p1.x = 0.047319;
        mesh[254].p1.y = 0.12642;
        mesh[254].p1.z = 0.207304;
        mesh[254].p2.x = 0.051207;
        mesh[254].p2.y = 0.016414;
        mesh[254].p2.z = 0.221731;
        mesh[254].p3.x = 0.058191;
        mesh[254].p3.y = 0.09789;
        mesh[254].p3.z = 0.14702;
        mesh[255].p1.x = 0.051918;
        mesh[255].p1.y = 0.123289;
        mesh[255].p1.z = 0.067837;
        mesh[255].p2.x = 0.051959;
        mesh[255].p2.y = 0.116839;
        mesh[255].p2.z = 0.027809;
        mesh[255].p3.x = 0.055369;
        mesh[255].p3.y = 0.132142;
        mesh[255].p3.z = 0.02032;
        mesh[256].p1.x = 0.051207;
        mesh[256].p1.y = 0.016414;
        mesh[256].p1.z = 0.221731;
        mesh[256].p2.x = 0.047319;
        mesh[256].p2.y = 0.12642;
        mesh[256].p2.z = 0.207304;
        mesh[256].p3.x = 0.035825;
        mesh[256].p3.y = 0.085919;
        mesh[256].p3.z = 0.251555;
        mesh[257].p1.x = 0.030553;
        mesh[257].p1.y = 0.149125;
        mesh[257].p1.z = 0.079185;
        mesh[257].p2.x = 0.037267;
        mesh[257].p2.y = 0.142487;
        mesh[257].p2.z = 0.021146;
        mesh[257].p3.x = 0.029122;
        mesh[257].p3.y = 0.138151;
        mesh[257].p3.z = 0.027991;
        mesh[258].p1.x = 0.0337;
        mesh[258].p1.y = 0.048755;
        mesh[258].p1.z = 0.136203;
        mesh[258].p2.x = 0.020677;
        mesh[258].p2.y = 0.005463;
        mesh[258].p2.z = 0.12829;
        mesh[258].p3.x = 0.005073;
        mesh[258].p3.y = 0.031262;
        mesh[258].p3.z = 0.135312;
        mesh[259].p1.x = 0.0337;
        mesh[259].p1.y = 0.048755;
        mesh[259].p1.z = 0.136203;
        mesh[259].p2.x = 0.005073;
        mesh[259].p2.y = 0.031262;
        mesh[259].p2.z = 0.135312;
        mesh[259].p3.x = 0.035775;
        mesh[259].p3.y = 0.072902;
        mesh[259].p3.z = 0.121575;
        mesh[260].p1.x = 0.056245;
        mesh[260].p1.y = -0.133422;
        mesh[260].p1.z = 0.00208;
        mesh[260].p2.x = 0.046714;
        mesh[260].p2.y = -0.146287;
        mesh[260].p2.z = 0.006232;
        mesh[260].p3.x = 0.049755;
        mesh[260].p3.y = -0.105056;
        mesh[260].p3.z = 0.003318;
        mesh[261].p1.x = 0.053684;
        mesh[261].p1.y = 0.125516;
        mesh[261].p1.z = 0.002388;
        mesh[261].p2.x = 0.049285;
        mesh[261].p2.y = 0.093207;
        mesh[261].p2.z = 1.43E-4;
        mesh[261].p3.x = 0.038873;
        mesh[261].p3.y = 0.105471;
        mesh[261].p3.z = 5.0E-5;
        mesh[262].p1.x = 0.025813;
        mesh[262].p1.y = 0.096198;
        mesh[262].p1.z = -1.51E-4;
        mesh[262].p2.x = 0.023355;
        mesh[262].p2.y = 0.127603;
        mesh[262].p2.z = 0.002924;
        mesh[262].p3.x = 0.038873;
        mesh[262].p3.y = 0.105471;
        mesh[262].p3.z = 5.0E-5;
        mesh[263].p1.x = 0.053684;
        mesh[263].p1.y = 0.125516;
        mesh[263].p1.z = 0.002388;
        mesh[263].p2.x = 0.057539;
        mesh[263].p2.y = 0.09589;
        mesh[263].p2.z = 3.54E-4;
        mesh[263].p3.x = 0.049285;
        mesh[263].p3.y = 0.093207;
        mesh[263].p3.z = 1.43E-4;
        mesh[264].p1.x = 0.026391;
        mesh[264].p1.y = 0.117277;
        mesh[264].p1.z = 0.021774;
        mesh[264].p2.x = 0.023355;
        mesh[264].p2.y = 0.127603;
        mesh[264].p2.z = 0.002924;
        mesh[264].p3.x = 0.025813;
        mesh[264].p3.y = 0.096198;
        mesh[264].p3.z = -1.51E-4;
        mesh[265].p1.x = 0.038873;
        mesh[265].p1.y = 0.105471;
        mesh[265].p1.z = 5.0E-5;
        mesh[265].p2.x = 0.049285;
        mesh[265].p2.y = 0.093207;
        mesh[265].p2.z = 1.43E-4;
        mesh[265].p3.x = 0.036214;
        mesh[265].p3.y = 0.086871;
        mesh[265].p3.z = 0.00329;
        mesh[266].p1.x = 0.01717;
        mesh[266].p1.y = -0.220033;
        mesh[266].p1.z = 0.328763;
        mesh[266].p2.x = 0.029393;
        mesh[266].p2.y = -0.219512;
        mesh[266].p2.z = 0.325764;
        mesh[266].p3.x = 0.024632;
        mesh[266].p3.y = -0.221933;
        mesh[266].p3.z = 0.335448;
        mesh[267].p1.x = -0.02214;
        mesh[267].p1.y = -0.257645;
        mesh[267].p1.z = 0.316165;
        mesh[267].p2.x = -0.009455;
        mesh[267].p2.y = -0.266105;
        mesh[267].p2.z = 0.31138;
        mesh[267].p3.x = -0.0117;
        mesh[267].p3.y = -0.264624;
        mesh[267].p3.z = 0.32503;
        mesh[268].p1.x = 0.002083;
        mesh[268].p1.y = -0.267639;
        mesh[268].p1.z = 0.312103;
        mesh[268].p2.x = -0.005389;
        mesh[268].p2.y = -0.264282;
        mesh[268].p2.z = 0.316817;
        mesh[268].p3.x = 0.00371;
        mesh[268].p3.y = -0.265544;
        mesh[268].p3.z = 0.295988;
        mesh[269].p1.x = -0.009455;
        mesh[269].p1.y = -0.266105;
        mesh[269].p1.z = 0.31138;
        mesh[269].p2.x = 0.00371;
        mesh[269].p2.y = -0.265544;
        mesh[269].p2.z = 0.295988;
        mesh[269].p3.x = -0.005389;
        mesh[269].p3.y = -0.264282;
        mesh[269].p3.z = 0.316817;
        mesh[270].p1.x = -0.009455;
        mesh[270].p1.y = -0.266105;
        mesh[270].p1.z = 0.31138;
        mesh[270].p2.x = -0.01772;
        mesh[270].p2.y = -0.262421;
        mesh[270].p2.z = 0.287041;
        mesh[270].p3.x = 0.00371;
        mesh[270].p3.y = -0.265544;
        mesh[270].p3.z = 0.295988;
        mesh[271].p1.x = -0.028845;
        mesh[271].p1.y = -0.219951;
        mesh[271].p1.z = 0.330842;
        mesh[271].p2.x = -0.039859;
        mesh[271].p2.y = -0.204217;
        mesh[271].p2.z = 0.325851;
        mesh[271].p3.x = -0.018638;
        mesh[271].p3.y = -0.221672;
        mesh[271].p3.z = 0.328157;
        mesh[272].p1.x = 0.0;
        mesh[272].p1.y = -0.195958;
        mesh[272].p1.z = 0.354306;
        mesh[272].p2.x = -0.021594;
        mesh[272].p2.y = -0.217944;
        mesh[272].p2.z = 0.346372;
        mesh[272].p3.x = 0.020438;
        mesh[272].p3.y = -0.217917;
        mesh[272].p3.z = 0.346311;
        mesh[273].p1.x = -0.021594;
        mesh[273].p1.y = -0.217944;
        mesh[273].p1.z = 0.346372;
        mesh[273].p2.x = -0.009896;
        mesh[273].p2.y = -0.227817;
        mesh[273].p2.z = 0.334942;
        mesh[273].p3.x = 0.020438;
        mesh[273].p3.y = -0.217917;
        mesh[273].p3.z = 0.346311;
        mesh[274].p1.x = -0.030009;
        mesh[274].p1.y = -0.229325;
        mesh[274].p1.z = 0.281517;
        mesh[274].p2.x = -0.01772;
        mesh[274].p2.y = -0.262421;
        mesh[274].p2.z = 0.287041;
        mesh[274].p3.x = -0.026484;
        mesh[274].p3.y = -0.259787;
        mesh[274].p3.z = 0.297418;
        mesh[275].p1.x = -0.018638;
        mesh[275].p1.y = -0.221672;
        mesh[275].p1.z = 0.328157;
        mesh[275].p2.x = -0.029006;
        mesh[275].p2.y = -0.229525;
        mesh[275].p2.z = 0.311337;
        mesh[275].p3.x = -0.009896;
        mesh[275].p3.y = -0.227817;
        mesh[275].p3.z = 0.334942;
        mesh[276].p1.x = -0.029006;
        mesh[276].p1.y = -0.229525;
        mesh[276].p1.z = 0.311337;
        mesh[276].p2.x = -0.026484;
        mesh[276].p2.y = -0.259787;
        mesh[276].p2.z = 0.297418;
        mesh[276].p3.x = -0.02214;
        mesh[276].p3.y = -0.257645;
        mesh[276].p3.z = 0.316165;
        mesh[277].p1.x = -0.02214;
        mesh[277].p1.y = -0.257645;
        mesh[277].p1.z = 0.316165;
        mesh[277].p2.x = -0.026484;
        mesh[277].p2.y = -0.259787;
        mesh[277].p2.z = 0.297418;
        mesh[277].p3.x = -0.009455;
        mesh[277].p3.y = -0.266105;
        mesh[277].p3.z = 0.31138;
        mesh[278].p1.x = -0.030009;
        mesh[278].p1.y = -0.229325;
        mesh[278].p1.z = 0.281517;
        mesh[278].p2.x = -0.029006;
        mesh[278].p2.y = -0.229525;
        mesh[278].p2.z = 0.311337;
        mesh[278].p3.x = -0.039859;
        mesh[278].p3.y = -0.204217;
        mesh[278].p3.z = 0.325851;
        mesh[279].p1.x = -0.018638;
        mesh[279].p1.y = -0.221672;
        mesh[279].p1.z = 0.328157;
        mesh[279].p2.x = -0.039859;
        mesh[279].p2.y = -0.204217;
        mesh[279].p2.z = 0.325851;
        mesh[279].p3.x = -0.029006;
        mesh[279].p3.y = -0.229525;
        mesh[279].p3.z = 0.311337;
        mesh[280].p1.x = -0.0117;
        mesh[280].p1.y = -0.264624;
        mesh[280].p1.z = 0.32503;
        mesh[280].p2.x = 0.007313;
        mesh[280].p2.y = -0.265396;
        mesh[280].p2.z = 0.327658;
        mesh[280].p3.x = -0.009896;
        mesh[280].p3.y = -0.227817;
        mesh[280].p3.z = 0.334942;
        mesh[281].p1.x = 0.011163;
        mesh[281].p1.y = -0.226997;
        mesh[281].p1.z = 0.334062;
        mesh[281].p2.x = -0.009896;
        mesh[281].p2.y = -0.227817;
        mesh[281].p2.z = 0.334942;
        mesh[281].p3.x = 0.007313;
        mesh[281].p3.y = -0.265396;
        mesh[281].p3.z = 0.327658;
        mesh[282].p1.x = -0.0117;
        mesh[282].p1.y = -0.264624;
        mesh[282].p1.z = 0.32503;
        mesh[282].p2.x = -0.009896;
        mesh[282].p2.y = -0.227817;
        mesh[282].p2.z = 0.334942;
        mesh[282].p3.x = -0.02214;
        mesh[282].p3.y = -0.257645;
        mesh[282].p3.z = 0.316165;
        mesh[283].p1.x = -0.030009;
        mesh[283].p1.y = -0.229325;
        mesh[283].p1.z = 0.281517;
        mesh[283].p2.x = -0.026484;
        mesh[283].p2.y = -0.259787;
        mesh[283].p2.z = 0.297418;
        mesh[283].p3.x = -0.029006;
        mesh[283].p3.y = -0.229525;
        mesh[283].p3.z = 0.311337;
        mesh[284].p1.x = -0.022417;
        mesh[284].p1.y = -0.136843;
        mesh[284].p1.z = 0.326557;
        mesh[284].p2.x = -0.039008;
        mesh[284].p2.y = -0.158249;
        mesh[284].p2.z = 0.27949;
        mesh[284].p3.x = -0.037814;
        mesh[284].p3.y = -0.136474;
        mesh[284].p3.z = 0.309237;
        mesh[285].p1.x = -0.039859;
        mesh[285].p1.y = -0.204217;
        mesh[285].p1.z = 0.325851;
        mesh[285].p2.x = -0.021594;
        mesh[285].p2.y = -0.217944;
        mesh[285].p2.z = 0.346372;
        mesh[285].p3.x = -0.043043;
        mesh[285].p3.y = -0.1738;
        mesh[285].p3.z = 0.335948;
        mesh[286].p1.x = -0.048469;
        mesh[286].p1.y = -0.158545;
        mesh[286].p1.z = 0.319398;
        mesh[286].p2.x = -0.039859;
        mesh[286].p2.y = -0.204217;
        mesh[286].p2.z = 0.325851;
        mesh[286].p3.x = -0.043043;
        mesh[286].p3.y = -0.1738;
        mesh[286].p3.z = 0.335948;
        mesh[287].p1.x = -0.048469;
        mesh[287].p1.y = -0.158545;
        mesh[287].p1.z = 0.319398;
        mesh[287].p2.x = -0.037814;
        mesh[287].p2.y = -0.136474;
        mesh[287].p2.z = 0.309237;
        mesh[287].p3.x = -0.039008;
        mesh[287].p3.y = -0.158249;
        mesh[287].p3.z = 0.27949;
        mesh[288].p1.x = -0.053788;
        mesh[288].p1.y = -0.129921;
        mesh[288].p1.z = 0.303985;
        mesh[288].p2.x = -0.022417;
        mesh[288].p2.y = -0.136843;
        mesh[288].p2.z = 0.326557;
        mesh[288].p3.x = -0.037814;
        mesh[288].p3.y = -0.136474;
        mesh[288].p3.z = 0.309237;
        mesh[289].p1.x = -0.049253;
        mesh[289].p1.y = -0.135547;
        mesh[289].p1.z = 0.325422;
        mesh[289].p2.x = -0.022417;
        mesh[289].p2.y = -0.136843;
        mesh[289].p2.z = 0.326557;
        mesh[289].p3.x = -0.053788;
        mesh[289].p3.y = -0.129921;
        mesh[289].p3.z = 0.303985;
        mesh[290].p1.x = -0.049253;
        mesh[290].p1.y = -0.135547;
        mesh[290].p1.z = 0.325422;
        mesh[290].p2.x = -0.066432;
        mesh[290].p2.y = -0.16482;
        mesh[290].p2.z = 0.29205;
        mesh[290].p3.x = -0.039197;
        mesh[290].p3.y = -0.156989;
        mesh[290].p3.z = 0.339227;
        mesh[291].p1.x = -0.053788;
        mesh[291].p1.y = -0.129921;
        mesh[291].p1.z = 0.303985;
        mesh[291].p2.x = -0.066432;
        mesh[291].p2.y = -0.16482;
        mesh[291].p2.z = 0.29205;
        mesh[291].p3.x = -0.049253;
        mesh[291].p3.y = -0.135547;
        mesh[291].p3.z = 0.325422;
        mesh[292].p1.x = -0.043043;
        mesh[292].p1.y = -0.1738;
        mesh[292].p1.z = 0.335948;
        mesh[292].p2.x = -0.044853;
        mesh[292].p2.y = -0.189594;
        mesh[292].p2.z = 0.271944;
        mesh[292].p3.x = -0.063422;
        mesh[292].p3.y = -0.164003;
        mesh[292].p3.z = 0.271249;
        mesh[293].p1.x = -0.063422;
        mesh[293].p1.y = -0.164003;
        mesh[293].p1.z = 0.271249;
        mesh[293].p2.x = -0.048469;
        mesh[293].p2.y = -0.158545;
        mesh[293].p2.z = 0.319398;
        mesh[293].p3.x = -0.043043;
        mesh[293].p3.y = -0.1738;
        mesh[293].p3.z = 0.335948;
        mesh[294].p1.x = -0.039197;
        mesh[294].p1.y = -0.156989;
        mesh[294].p1.z = 0.339227;
        mesh[294].p2.x = 0.0;
        mesh[294].p2.y = -0.195958;
        mesh[294].p2.z = 0.354306;
        mesh[294].p3.x = -7.79E-4;
        mesh[294].p3.y = -0.146982;
        mesh[294].p3.z = 0.344462;
        mesh[295].p1.x = -0.043043;
        mesh[295].p1.y = -0.1738;
        mesh[295].p1.z = 0.335948;
        mesh[295].p2.x = 0.0;
        mesh[295].p2.y = -0.195958;
        mesh[295].p2.z = 0.354306;
        mesh[295].p3.x = -0.039197;
        mesh[295].p3.y = -0.156989;
        mesh[295].p3.z = 0.339227;
        mesh[296].p1.x = 0.0;
        mesh[296].p1.y = -0.195958;
        mesh[296].p1.z = 0.354306;
        mesh[296].p2.x = -0.043043;
        mesh[296].p2.y = -0.1738;
        mesh[296].p2.z = 0.335948;
        mesh[296].p3.x = -0.021594;
        mesh[296].p3.y = -0.217944;
        mesh[296].p3.z = 0.346372;
        mesh[297].p1.x = -0.022417;
        mesh[297].p1.y = -0.136843;
        mesh[297].p1.z = 0.326557;
        mesh[297].p2.x = -0.024667;
        mesh[297].p2.y = -0.101688;
        mesh[297].p2.z = 0.262395;
        mesh[297].p3.x = -0.039008;
        mesh[297].p3.y = -0.158249;
        mesh[297].p3.z = 0.27949;
        mesh[298].p1.x = -0.017193;
        mesh[298].p1.y = -0.210889;
        mesh[298].p1.z = 0.257864;
        mesh[298].p2.x = -0.030009;
        mesh[298].p2.y = -0.229325;
        mesh[298].p2.z = 0.281517;
        mesh[298].p3.x = -0.00563;
        mesh[298].p3.y = -0.191954;
        mesh[298].p3.z = 0.225786;
        mesh[299].p1.x = -0.01772;
        mesh[299].p1.y = -0.262421;
        mesh[299].p1.z = 0.287041;
        mesh[299].p2.x = -0.013919;
        mesh[299].p2.y = -0.25944;
        mesh[299].p2.z = 0.288499;
        mesh[299].p3.x = 0.00371;
        mesh[299].p3.y = -0.265544;
        mesh[299].p3.z = 0.295988;
        mesh[300].p1.x = -0.017193;
        mesh[300].p1.y = -0.210889;
        mesh[300].p1.z = 0.257864;
        mesh[300].p2.x = 0.01869;
        mesh[300].p2.y = -0.214564;
        mesh[300].p2.z = 0.260349;
        mesh[300].p3.x = 5.65E-4;
        mesh[300].p3.y = -0.260455;
        mesh[300].p3.z = 0.278232;
        mesh[301].p1.x = -0.043082;
        mesh[301].p1.y = -0.121347;
        mesh[301].p1.z = 0.0261;
        mesh[301].p2.x = -0.061288;
        mesh[301].p2.y = -0.130257;
        mesh[301].p2.z = 0.005966;
        mesh[301].p3.x = -0.04676;
        mesh[301].p3.y = -0.146333;
        mesh[301].p3.z = 0.006207;
        mesh[302].p1.x = -0.059617;
        mesh[302].p1.y = -0.1185;
        mesh[302].p1.z = 0.120599;
        mesh[302].p2.x = -0.061648;
        mesh[302].p2.y = -0.116842;
        mesh[302].p2.z = 0.194265;
        mesh[302].p3.x = -0.056427;
        mesh[302].p3.y = -0.085244;
        mesh[302].p3.z = 0.187472;
        mesh[303].p1.x = -0.051674;
        mesh[303].p1.y = -0.109716;
        mesh[303].p1.z = 0.023234;
        mesh[303].p2.x = -0.059617;
        mesh[303].p2.y = -0.1185;
        mesh[303].p2.z = 0.120599;
        mesh[303].p3.x = -0.047173;
        mesh[303].p3.y = -0.08777;
        mesh[303].p3.z = 0.116024;
        mesh[304].p1.x = -0.040724;
        mesh[304].p1.y = -0.084142;
        mesh[304].p1.z = 0.048356;
        mesh[304].p2.x = -0.017058;
        mesh[304].p2.y = -0.107845;
        mesh[304].p2.z = 0.120886;
        mesh[304].p3.x = -0.025373;
        mesh[304].p3.y = -0.09533;
        mesh[304].p3.z = 0.018332;
        mesh[305].p1.x = -0.057323;
        mesh[305].p1.y = -0.147545;
        mesh[305].p1.z = 0.177989;
        mesh[305].p2.x = -0.059617;
        mesh[305].p2.y = -0.1185;
        mesh[305].p2.z = 0.120599;
        mesh[305].p3.x = -0.042681;
        mesh[305].p3.y = -0.137565;
        mesh[305].p3.z = 0.125192;
        mesh[306].p1.x = -0.047097;
        mesh[306].p1.y = 0.052817;
        mesh[306].p1.z = 0.162408;
        mesh[306].p2.x = -0.056427;
        mesh[306].p2.y = -0.085244;
        mesh[306].p2.z = 0.187472;
        mesh[306].p3.x = -0.052194;
        mesh[306].p3.y = 0.037269;
        mesh[306].p3.z = 0.199662;
        mesh[307].p1.x = -0.017058;
        mesh[307].p1.y = -0.107845;
        mesh[307].p1.z = 0.120886;
        mesh[307].p2.x = -0.038493;
        mesh[307].p2.y = -0.084398;
        mesh[307].p2.z = 0.135693;
        mesh[307].p3.x = -0.027136;
        mesh[307].p3.y = 0.041664;
        mesh[307].p3.z = 0.132141;
        mesh[308].p1.x = 0.01653;
        mesh[308].p1.y = -0.109257;
        mesh[308].p1.z = 0.120113;
        mesh[308].p2.x = -0.017058;
        mesh[308].p2.y = -0.107845;
        mesh[308].p2.z = 0.120886;
        mesh[308].p3.x = -0.027136;
        mesh[308].p3.y = 0.041664;
        mesh[308].p3.z = 0.132141;
        mesh[309].p1.x = -0.006761;
        mesh[309].p1.y = -0.10839;
        mesh[309].p1.z = 0.287993;
        mesh[309].p2.x = 0.003674;
        mesh[309].p2.y = -0.085117;
        mesh[309].p2.z = 0.272056;
        mesh[309].p3.x = -0.024667;
        mesh[309].p3.y = -0.101688;
        mesh[309].p3.z = 0.262395;
        mesh[310].p1.x = -0.048878;
        mesh[310].p1.y = -0.099334;
        mesh[310].p1.z = 0.010049;
        mesh[310].p2.x = -0.040724;
        mesh[310].p2.y = -0.084142;
        mesh[310].p2.z = 0.048356;
        mesh[310].p3.x = -0.025373;
        mesh[310].p3.y = -0.09533;
        mesh[310].p3.z = 0.018332;
        mesh[311].p1.x = -0.040724;
        mesh[311].p1.y = -0.084142;
        mesh[311].p1.z = 0.048356;
        mesh[311].p2.x = -0.048878;
        mesh[311].p2.y = -0.099334;
        mesh[311].p2.z = 0.010049;
        mesh[311].p3.x = -0.051674;
        mesh[311].p3.y = -0.109716;
        mesh[311].p3.z = 0.023234;
        mesh[312].p1.x = -0.031352;
        mesh[312].p1.y = -0.118191;
        mesh[312].p1.z = 0.050941;
        mesh[312].p2.x = -0.051674;
        mesh[312].p2.y = -0.109716;
        mesh[312].p2.z = 0.023234;
        mesh[312].p3.x = -0.043082;
        mesh[312].p3.y = -0.121347;
        mesh[312].p3.z = 0.0261;
        mesh[313].p1.x = -0.043082;
        mesh[313].p1.y = -0.121347;
        mesh[313].p1.z = 0.0261;
        mesh[313].p2.x = -0.04676;
        mesh[313].p2.y = -0.146333;
        mesh[313].p2.z = 0.006207;
        mesh[313].p3.x = -0.029086;
        mesh[313].p3.y = -0.121207;
        mesh[313].p3.z = 0.020285;
        mesh[314].p1.x = -0.029086;
        mesh[314].p1.y = -0.121207;
        mesh[314].p1.z = 0.020285;
        mesh[314].p2.x = -0.02965;
        mesh[314].p2.y = -0.136573;
        mesh[314].p2.z = 0.00164;
        mesh[314].p3.x = -0.021311;
        mesh[314].p3.y = -0.109992;
        mesh[314].p3.z = 0.007536;
        mesh[315].p1.x = -0.048878;
        mesh[315].p1.y = -0.099334;
        mesh[315].p1.z = 0.010049;
        mesh[315].p2.x = -0.021311;
        mesh[315].p2.y = -0.109992;
        mesh[315].p2.z = 0.007536;
        mesh[315].p3.x = -0.050702;
        mesh[315].p3.y = -0.129375;
        mesh[315].p3.z = 6.69E-4;
        mesh[316].p1.x = -0.04676;
        mesh[316].p1.y = -0.146333;
        mesh[316].p1.z = 0.006207;
        mesh[316].p2.x = -0.061288;
        mesh[316].p2.y = -0.130257;
        mesh[316].p2.z = 0.005966;
        mesh[316].p3.x = -0.050702;
        mesh[316].p3.y = -0.129375;
        mesh[316].p3.z = 6.69E-4;
        mesh[317].p1.x = -0.014487;
        mesh[317].p1.y = 0.080367;
        mesh[317].p1.z = 0.261155;
        mesh[317].p2.x = 0.004478;
        mesh[317].p2.y = 0.178564;
        mesh[317].p2.z = 0.232537;
        mesh[317].p3.x = -0.017841;
        mesh[317].p3.y = 0.13652;
        mesh[317].p3.z = 0.230347;
        mesh[318].p1.x = -0.030328;
        mesh[318].p1.y = 0.140036;
        mesh[318].p1.z = 0.018044;
        mesh[318].p2.x = -0.023133;
        mesh[318].p2.y = 0.120981;
        mesh[318].p2.z = 0.078476;
        mesh[318].p3.x = -0.027847;
        mesh[318].p3.y = 0.120028;
        mesh[318].p3.z = 0.02675;
        mesh[319].p1.x = -0.057904;
        mesh[319].p1.y = 0.124382;
        mesh[319].p1.z = 0.020753;
        mesh[319].p2.x = -0.047607;
        mesh[319].p2.y = 0.104629;
        mesh[319].p2.z = 0.020424;
        mesh[319].p3.x = -0.045197;
        mesh[319].p3.y = 0.121199;
        mesh[319].p3.z = 0.055852;
        mesh[320].p1.x = -0.03108;
        mesh[320].p1.y = 0.104637;
        mesh[320].p1.z = 7.6E-5;
        mesh[320].p2.x = -0.031653;
        mesh[320].p2.y = 0.105782;
        mesh[320].p2.z = 0.018291;
        mesh[320].p3.x = -0.033126;
        mesh[320].p3.y = 0.084624;
        mesh[320].p3.z = 0.003634;
        mesh[321].p1.x = -0.022864;
        mesh[321].p1.y = 0.094117;
        mesh[321].p1.z = -3.86E-4;
        mesh[321].p2.x = -0.027847;
        mesh[321].p2.y = 0.120028;
        mesh[321].p2.z = 0.02675;
        mesh[321].p3.x = -0.031653;
        mesh[321].p3.y = 0.105782;
        mesh[321].p3.z = 0.018291;
        mesh[322].p1.x = -0.043016;
        mesh[322].p1.y = 0.11447;
        mesh[322].p1.z = -4.97E-4;
        mesh[322].p2.x = -0.030328;
        mesh[322].p2.y = 0.140036;
        mesh[322].p2.z = 0.018044;
        mesh[322].p3.x = -0.023356;
        mesh[322].p3.y = 0.126501;
        mesh[322].p3.z = 0.002503;
        mesh[323].p1.x = -0.057904;
        mesh[323].p1.y = 0.124382;
        mesh[323].p1.z = 0.020753;
        mesh[323].p2.x = -0.044942;
        mesh[323].p2.y = 0.141872;
        mesh[323].p2.z = 0.019356;
        mesh[323].p3.x = -0.055813;
        mesh[323].p3.y = 0.124066;
        mesh[323].p3.z = 7.64E-4;
        mesh[324].p1.x = -0.055813;
        mesh[324].p1.y = 0.124066;
        mesh[324].p1.z = 7.64E-4;
        mesh[324].p2.x = -0.052415;
        mesh[324].p2.y = 0.095648;
        mesh[324].p2.z = -3.5E-5;
        mesh[324].p3.x = -0.057904;
        mesh[324].p3.y = 0.124382;
        mesh[324].p3.z = 0.020753;
        mesh[325].p1.x = -0.012166;
        mesh[325].p1.y = 0.073247;
        mesh[325].p1.z = 0.136132;
        mesh[325].p2.x = -0.030099;
        mesh[325].p2.y = 0.074776;
        mesh[325].p2.z = 0.120707;
        mesh[325].p3.x = -0.045197;
        mesh[325].p3.y = 0.121199;
        mesh[325].p3.z = 0.055852;
        mesh[326].p1.x = -0.059617;
        mesh[326].p1.y = -0.1185;
        mesh[326].p1.z = 0.120599;
        mesh[326].p2.x = -0.051674;
        mesh[326].p2.y = -0.109716;
        mesh[326].p2.z = 0.023234;
        mesh[326].p3.x = -0.042681;
        mesh[326].p3.y = -0.137565;
        mesh[326].p3.z = 0.125192;
        mesh[327].p1.x = -0.051674;
        mesh[327].p1.y = -0.109716;
        mesh[327].p1.z = 0.023234;
        mesh[327].p2.x = -0.031352;
        mesh[327].p2.y = -0.118191;
        mesh[327].p2.z = 0.050941;
        mesh[327].p3.x = -0.042681;
        mesh[327].p3.y = -0.137565;
        mesh[327].p3.z = 0.125192;
        mesh[328].p1.x = -0.020833;
        mesh[328].p1.y = -0.125757;
        mesh[328].p1.z = 0.124804;
        mesh[328].p2.x = -0.025373;
        mesh[328].p2.y = -0.09533;
        mesh[328].p2.z = 0.018332;
        mesh[328].p3.x = -0.017058;
        mesh[328].p3.y = -0.107845;
        mesh[328].p3.z = 0.120886;
        mesh[329].p1.x = -0.040724;
        mesh[329].p1.y = -0.084142;
        mesh[329].p1.z = 0.048356;
        mesh[329].p2.x = -0.047173;
        mesh[329].p2.y = -0.08777;
        mesh[329].p2.z = 0.116024;
        mesh[329].p3.x = -0.017058;
        mesh[329].p3.y = -0.107845;
        mesh[329].p3.z = 0.120886;
        mesh[330].p1.x = -0.031352;
        mesh[330].p1.y = -0.118191;
        mesh[330].p1.z = 0.050941;
        mesh[330].p2.x = -0.020833;
        mesh[330].p2.y = -0.125757;
        mesh[330].p2.z = 0.124804;
        mesh[330].p3.x = -0.042681;
        mesh[330].p3.y = -0.137565;
        mesh[330].p3.z = 0.125192;
        mesh[331].p1.x = -0.048878;
        mesh[331].p1.y = -0.099334;
        mesh[331].p1.z = 0.010049;
        mesh[331].p2.x = -0.061288;
        mesh[331].p2.y = -0.130257;
        mesh[331].p2.z = 0.005966;
        mesh[331].p3.x = -0.051674;
        mesh[331].p3.y = -0.109716;
        mesh[331].p3.z = 0.023234;
        mesh[332].p1.x = -0.045197;
        mesh[332].p1.y = 0.121199;
        mesh[332].p1.z = 0.055852;
        mesh[332].p2.x = -0.047607;
        mesh[332].p2.y = 0.104629;
        mesh[332].p2.z = 0.020424;
        mesh[332].p3.x = -0.031653;
        mesh[332].p3.y = 0.105782;
        mesh[332].p3.z = 0.018291;
        mesh[333].p1.x = -0.031653;
        mesh[333].p1.y = 0.105782;
        mesh[333].p1.z = 0.018291;
        mesh[333].p2.x = -0.027847;
        mesh[333].p2.y = 0.120028;
        mesh[333].p2.z = 0.02675;
        mesh[333].p3.x = -0.045197;
        mesh[333].p3.y = 0.121199;
        mesh[333].p3.z = 0.055852;
        mesh[334].p1.x = -0.030328;
        mesh[334].p1.y = 0.140036;
        mesh[334].p1.z = 0.018044;
        mesh[334].p2.x = -0.037603;
        mesh[334].p2.y = 0.15207;
        mesh[334].p2.z = 0.084055;
        mesh[334].p3.x = -0.023133;
        mesh[334].p3.y = 0.120981;
        mesh[334].p3.z = 0.078476;
        mesh[335].p1.x = -0.04672;
        mesh[335].p1.y = 0.093756;
        mesh[335].p1.z = 0.235686;
        mesh[335].p2.x = -0.047107;
        mesh[335].p2.y = -0.005026;
        mesh[335].p2.z = 0.233809;
        mesh[335].p3.x = -0.014487;
        mesh[335].p3.y = 0.080367;
        mesh[335].p3.z = 0.261155;
        mesh[336].p1.x = -0.009247;
        mesh[336].p1.y = 0.106663;
        mesh[336].p1.z = 0.145238;
        mesh[336].p2.x = -0.012166;
        mesh[336].p2.y = 0.073247;
        mesh[336].p2.z = 0.136132;
        mesh[336].p3.x = -0.023133;
        mesh[336].p3.y = 0.120981;
        mesh[336].p3.z = 0.078476;
        mesh[337].p1.x = -0.012166;
        mesh[337].p1.y = 0.073247;
        mesh[337].p1.z = 0.136132;
        mesh[337].p2.x = -3.31E-4;
        mesh[337].p2.y = 0.036278;
        mesh[337].p2.z = 0.123064;
        mesh[337].p3.x = -0.012945;
        mesh[337].p3.y = 0.045589;
        mesh[337].p3.z = 0.13692;
        mesh[338].p1.x = -0.051287;
        mesh[338].p1.y = 0.123923;
        mesh[338].p1.z = 0.200622;
        mesh[338].p2.x = -3.26E-4;
        mesh[338].p2.y = 0.141988;
        mesh[338].p2.z = 0.196277;
        mesh[338].p3.x = -0.031246;
        mesh[338].p3.y = 0.135959;
        mesh[338].p3.z = 0.153033;
        mesh[339].p1.x = -0.04676;
        mesh[339].p1.y = -0.146333;
        mesh[339].p1.z = 0.006207;
        mesh[339].p2.x = -0.021311;
        mesh[339].p2.y = -0.109992;
        mesh[339].p2.z = 0.007536;
        mesh[339].p3.x = -0.02965;
        mesh[339].p3.y = -0.136573;
        mesh[339].p3.z = 0.00164;
        mesh[340].p1.x = -0.021311;
        mesh[340].p1.y = -0.109992;
        mesh[340].p1.z = 0.007536;
        mesh[340].p2.x = -0.048878;
        mesh[340].p2.y = -0.099334;
        mesh[340].p2.z = 0.010049;
        mesh[340].p3.x = -0.025373;
        mesh[340].p3.y = -0.09533;
        mesh[340].p3.z = 0.018332;
        mesh[341].p1.x = -0.052194;
        mesh[341].p1.y = 0.037269;
        mesh[341].p1.z = 0.199662;
        mesh[341].p2.x = -0.05208;
        mesh[341].p2.y = 0.124926;
        mesh[341].p2.z = 0.07032;
        mesh[341].p3.x = -0.045197;
        mesh[341].p3.y = 0.121199;
        mesh[341].p3.z = 0.055852;
        mesh[342].p1.x = -0.017058;
        mesh[342].p1.y = -0.107845;
        mesh[342].p1.z = 0.120886;
        mesh[342].p2.x = -0.047173;
        mesh[342].p2.y = -0.08777;
        mesh[342].p2.z = 0.116024;
        mesh[342].p3.x = -0.038493;
        mesh[342].p3.y = -0.084398;
        mesh[342].p3.z = 0.135693;
        mesh[343].p1.x = -0.047173;
        mesh[343].p1.y = -0.08777;
        mesh[343].p1.z = 0.116024;
        mesh[343].p2.x = -0.059617;
        mesh[343].p2.y = -0.1185;
        mesh[343].p2.z = 0.120599;
        mesh[343].p3.x = -0.038493;
        mesh[343].p3.y = -0.084398;
        mesh[343].p3.z = 0.135693;
        mesh[344].p1.x = -0.024667;
        mesh[344].p1.y = -0.101688;
        mesh[344].p1.z = 0.262395;
        mesh[344].p2.x = 0.003674;
        mesh[344].p2.y = -0.085117;
        mesh[344].p2.z = 0.272056;
        mesh[344].p3.x = -0.010543;
        mesh[344].p3.y = -0.015599;
        mesh[344].p3.z = 0.262292;
        mesh[345].p1.x = -0.010941;
        mesh[345].p1.y = 0.235894;
        mesh[345].p1.z = 0.207299;
        mesh[345].p2.x = -0.003326;
        mesh[345].p2.y = 0.138226;
        mesh[345].p2.z = 0.212925;
        mesh[345].p3.x = -0.017841;
        mesh[345].p3.y = 0.13652;
        mesh[345].p3.z = 0.230347;
        mesh[346].p1.x = -1.45E-4;
        mesh[346].p1.y = 0.265834;
        mesh[346].p1.z = 0.193142;
        mesh[346].p2.x = -0.003326;
        mesh[346].p2.y = 0.138226;
        mesh[346].p2.z = 0.212925;
        mesh[346].p3.x = -0.010941;
        mesh[346].p3.y = 0.235894;
        mesh[346].p3.z = 0.207299;
        mesh[347].p1.x = 0.011832;
        mesh[347].p1.y = 0.073182;
        mesh[347].p1.z = 0.136522;
        mesh[347].p2.x = -3.31E-4;
        mesh[347].p2.y = 0.036278;
        mesh[347].p2.z = 0.123064;
        mesh[347].p3.x = -0.012166;
        mesh[347].p3.y = 0.073247;
        mesh[347].p3.z = 0.136132;
        mesh[348].p1.x = -0.012945;
        mesh[348].p1.y = 0.045589;
        mesh[348].p1.z = 0.13692;
        mesh[348].p2.x = -3.31E-4;
        mesh[348].p2.y = 0.036278;
        mesh[348].p2.z = 0.123064;
        mesh[348].p3.x = 0.005073;
        mesh[348].p3.y = 0.031262;
        mesh[348].p3.z = 0.135312;
        mesh[349].p1.x = -0.024667;
        mesh[349].p1.y = -0.101688;
        mesh[349].p1.z = 0.262395;
        mesh[349].p2.x = -0.022417;
        mesh[349].p2.y = -0.136843;
        mesh[349].p2.z = 0.326557;
        mesh[349].p3.x = -0.006761;
        mesh[349].p3.y = -0.10839;
        mesh[349].p3.z = 0.287993;
        mesh[350].p1.x = -0.047107;
        mesh[350].p1.y = -0.005026;
        mesh[350].p1.z = 0.233809;
        mesh[350].p2.x = -0.052194;
        mesh[350].p2.y = 0.037269;
        mesh[350].p2.z = 0.199662;
        mesh[350].p3.x = -0.056427;
        mesh[350].p3.y = -0.085244;
        mesh[350].p3.z = 0.187472;
        mesh[351].p1.x = -0.043115;
        mesh[351].p1.y = -0.141878;
        mesh[351].p1.z = 0.237418;
        mesh[351].p2.x = -0.057323;
        mesh[351].p2.y = -0.147545;
        mesh[351].p2.z = 0.177989;
        mesh[351].p3.x = -0.032359;
        mesh[351].p3.y = -0.171143;
        mesh[351].p3.z = 0.240091;
        mesh[352].p1.x = -0.00563;
        mesh[352].p1.y = -0.191954;
        mesh[352].p1.z = 0.225786;
        mesh[352].p2.x = -0.032359;
        mesh[352].p2.y = -0.171143;
        mesh[352].p2.z = 0.240091;
        mesh[352].p3.x = -0.057323;
        mesh[352].p3.y = -0.147545;
        mesh[352].p3.z = 0.177989;
        mesh[353].p1.x = -0.057323;
        mesh[353].p1.y = -0.147545;
        mesh[353].p1.z = 0.177989;
        mesh[353].p2.x = -0.043115;
        mesh[353].p2.y = -0.141878;
        mesh[353].p2.z = 0.237418;
        mesh[353].p3.x = -0.061648;
        mesh[353].p3.y = -0.116842;
        mesh[353].p3.z = 0.194265;
        mesh[354].p1.x = -0.025373;
        mesh[354].p1.y = -0.09533;
        mesh[354].p1.z = 0.018332;
        mesh[354].p2.x = -0.020833;
        mesh[354].p2.y = -0.125757;
        mesh[354].p2.z = 0.124804;
        mesh[354].p3.x = -0.031352;
        mesh[354].p3.y = -0.118191;
        mesh[354].p3.z = 0.050941;
        mesh[355].p1.x = -0.043043;
        mesh[355].p1.y = -0.1738;
        mesh[355].p1.z = 0.335948;
        mesh[355].p2.x = -0.039197;
        mesh[355].p2.y = -0.156989;
        mesh[355].p2.z = 0.339227;
        mesh[355].p3.x = -0.066432;
        mesh[355].p3.y = -0.16482;
        mesh[355].p3.z = 0.29205;
        mesh[356].p1.x = -0.053788;
        mesh[356].p1.y = -0.129921;
        mesh[356].p1.z = 0.303985;
        mesh[356].p2.x = -0.063422;
        mesh[356].p2.y = -0.164003;
        mesh[356].p2.z = 0.271249;
        mesh[356].p3.x = -0.066432;
        mesh[356].p3.y = -0.16482;
        mesh[356].p3.z = 0.29205;
        mesh[357].p1.x = -0.063422;
        mesh[357].p1.y = -0.164003;
        mesh[357].p1.z = 0.271249;
        mesh[357].p2.x = -0.053788;
        mesh[357].p2.y = -0.129921;
        mesh[357].p2.z = 0.303985;
        mesh[357].p3.x = -0.048469;
        mesh[357].p3.y = -0.158545;
        mesh[357].p3.z = 0.319398;
        mesh[358].p1.x = -0.023133;
        mesh[358].p1.y = 0.120981;
        mesh[358].p1.z = 0.078476;
        mesh[358].p2.x = -0.037203;
        mesh[358].p2.y = 0.138321;
        mesh[358].p2.z = 0.100764;
        mesh[358].p3.x = -0.031246;
        mesh[358].p3.y = 0.135959;
        mesh[358].p3.z = 0.153033;
        mesh[359].p1.x = -0.017841;
        mesh[359].p1.y = 0.13652;
        mesh[359].p1.z = 0.230347;
        mesh[359].p2.x = -0.051287;
        mesh[359].p2.y = 0.123923;
        mesh[359].p2.z = 0.200622;
        mesh[359].p3.x = -0.04672;
        mesh[359].p3.y = 0.093756;
        mesh[359].p3.z = 0.235686;
        mesh[360].p1.x = -0.017841;
        mesh[360].p1.y = 0.13652;
        mesh[360].p1.z = 0.230347;
        mesh[360].p2.x = -0.003326;
        mesh[360].p2.y = 0.138226;
        mesh[360].p2.z = 0.212925;
        mesh[360].p3.x = -3.26E-4;
        mesh[360].p3.y = 0.141988;
        mesh[360].p3.z = 0.196277;
        mesh[361].p1.x = -0.022417;
        mesh[361].p1.y = -0.136843;
        mesh[361].p1.z = 0.326557;
        mesh[361].p2.x = -0.049253;
        mesh[361].p2.y = -0.135547;
        mesh[361].p2.z = 0.325422;
        mesh[361].p3.x = -0.039197;
        mesh[361].p3.y = -0.156989;
        mesh[361].p3.z = 0.339227;
        mesh[362].p1.x = -0.030328;
        mesh[362].p1.y = 0.140036;
        mesh[362].p1.z = 0.018044;
        mesh[362].p2.x = -0.055813;
        mesh[362].p2.y = 0.124066;
        mesh[362].p2.z = 7.64E-4;
        mesh[362].p3.x = -0.044942;
        mesh[362].p3.y = 0.141872;
        mesh[362].p3.z = 0.019356;
        mesh[363].p1.x = -0.055813;
        mesh[363].p1.y = 0.124066;
        mesh[363].p1.z = 7.64E-4;
        mesh[363].p2.x = -0.030328;
        mesh[363].p2.y = 0.140036;
        mesh[363].p2.z = 0.018044;
        mesh[363].p3.x = -0.043016;
        mesh[363].p3.y = 0.11447;
        mesh[363].p3.z = -4.97E-4;
        mesh[364].p1.x = -0.030328;
        mesh[364].p1.y = 0.140036;
        mesh[364].p1.z = 0.018044;
        mesh[364].p2.x = -0.044942;
        mesh[364].p2.y = 0.141872;
        mesh[364].p2.z = 0.019356;
        mesh[364].p3.x = -0.037603;
        mesh[364].p3.y = 0.15207;
        mesh[364].p3.z = 0.084055;
        mesh[365].p1.x = -0.012945;
        mesh[365].p1.y = 0.045589;
        mesh[365].p1.z = 0.13692;
        mesh[365].p2.x = -0.030099;
        mesh[365].p2.y = 0.074776;
        mesh[365].p2.z = 0.120707;
        mesh[365].p3.x = -0.012166;
        mesh[365].p3.y = 0.073247;
        mesh[365].p3.z = 0.136132;
        mesh[366].p1.x = -0.052415;
        mesh[366].p1.y = 0.095648;
        mesh[366].p1.z = -3.5E-5;
        mesh[366].p2.x = -0.043016;
        mesh[366].p2.y = 0.11447;
        mesh[366].p2.z = -4.97E-4;
        mesh[366].p3.x = -0.03108;
        mesh[366].p3.y = 0.104637;
        mesh[366].p3.z = 7.6E-5;
        mesh[367].p1.x = -0.03108;
        mesh[367].p1.y = 0.104637;
        mesh[367].p1.z = 7.6E-5;
        mesh[367].p2.x = -0.023356;
        mesh[367].p2.y = 0.126501;
        mesh[367].p2.z = 0.002503;
        mesh[367].p3.x = -0.022864;
        mesh[367].p3.y = 0.094117;
        mesh[367].p3.z = -3.86E-4;
        mesh[368].p1.x = -0.039606;
        mesh[368].p1.y = 0.087492;
        mesh[368].p1.z = 0.004376;
        mesh[368].p2.x = -0.04562;
        mesh[368].p2.y = 0.088864;
        mesh[368].p2.z = -0.001593;
        mesh[368].p3.x = -0.03108;
        mesh[368].p3.y = 0.104637;
        mesh[368].p3.z = 7.6E-5;
        mesh[369].p1.x = -0.04562;
        mesh[369].p1.y = 0.088864;
        mesh[369].p1.z = -0.001593;
        mesh[369].p2.x = -0.052415;
        mesh[369].p2.y = 0.095648;
        mesh[369].p2.z = -3.5E-5;
        mesh[369].p3.x = -0.03108;
        mesh[369].p3.y = 0.104637;
        mesh[369].p3.z = 7.6E-5; 
    }

    int objlen = 370;

    struct vector light = {
	    .x = 0,
	    .y = -1,
	    .z = -0.1,
    };

    float rot = (((float) frameCount) / 10);

    rotateMesh(mesh, objlen, 3.1415f/2, rot +  3.1415f/2, 0);

    struct tri transd [objlen];

    struct mat translate;

    initOne(&translate);

    translate.m[0][0] = 10;
    translate.m[1][1] = 10;
    translate.m[2][2] = 10;

    translate.m[3][0] = 0;
    translate.m[3][1] = 2.5;
    translate.m[3][2] = 9;// + 4*(sinf(rot) + 0);


    if (keyDown[0]) cam->pos.z++;
    if (keyDown[2]) cam->pos.z--;
    if (keyDown[1]) cam->pos.x++;
    if (keyDown[3]) cam->pos.x--;

    if (keyDown[4]) cam->pos.y++;
    if (keyDown[6]) cam->pos.y--;

    //for (int i = 0; i<sizeof(keyDown); i++) keyDown[i] = 0;

    meshmultmat(mesh, transd, objlen, translate);

    fillMesh(transd, objlen, *cam, light);

    wireMesh(transd, objlen, *cam);
}

void run(void (*f)(int, struct camera *), int frameCount, float frameRate, struct camera * cam) {
    clock_t mt;
    mt = clock();
    clock_t t;
    t = clock();
    f(frameCount, cam);
    while ( ((double)t - (double)mt)/CLOCKS_PER_SEC < 1/frameRate ) {
	    t = clock();
    }
}

void engine(int frameCount, struct camera * cam) {
    //clear();
    loop(frameCount, cam);
    draw();
}

void* getIns() {
    char b[5];
    b[0] = ' ';
    int size = sizeof(keyDown);
    //system("xset r rate 10 100");
    read(STDIN_FILENO, b, 5);
    for (int i = 0; i<size; i++) {
        if (b[0] == keys[i]) {
            keyDown[i] = 1;
        } else {
            keyDown[i] = 0;
        }
    }
    //system("xset r rate 220 20");
    fflush(stdout);
}

void getIn(XEvent *e) {
    XKeyEvent *ke = (XKeyEvent *) e;
    int on = (e->type == KeyPress);
    if (ke->keycode == 25) keyDown[0] = on;
    if (ke->keycode == 38) keyDown[1] = on;
    if (ke->keycode == 39) keyDown[2] = on;
    if (ke->keycode == 40) keyDown[3] = on;

    if (ke->keycode == 31) keyDown[4] = on;
    if (ke->keycode == 44) keyDown[5] = on;
    if (ke->keycode == 45) keyDown[6] = on;
    if (ke->keycode == 46) keyDown[7] = on;
	//printf("coolio %c", ke->keycode);
}

int main()  {
    //printf("\033[47:47m");
    Display *disp;
    Window win;
    XEvent event;
    
    disp = XOpenDisplay(NULL);

    int screen = DefaultScreen(disp);

    win = XCreateSimpleWindow(disp, RootWindow(disp, screen), WIDTH, HEIGHT, WIDTH, HEIGHT, 1, BlackPixel(disp, screen), WhitePixel(disp, screen));

    static struct termios told, tnew;
    tcgetattr( STDIN_FILENO, &told);
    tnew = told;
    tnew.c_lflag &= ~(ICANON | ECHO);
    tcsetattr( STDIN_FILENO, TCSANOW, &tnew);
    pthread_t getInputs;
    int frameCount = 0;

    struct camera cam;
    initCamera(&cam);

    XSelectInput(disp, win, ExposureMask | KeyPressMask | KeyReleaseMask);
    XMapWindow(disp, win);
    //pthread_create(&getInputs, NULL, getIns, (void*)&args);
    do {
	    //pthread_create(&getInputs, NULL, getIns, (void*)&args);
	    //if (fgetc(stdin) == 'w') printf("cool");
	    run(engine, frameCount % 62, 30, &cam);
        //XNextEvent(disp, &event);
        while (XCheckWindowEvent(disp, win, ExposureMask | KeyPressMask | KeyReleaseMask, &event)) {
            printf("\n");
            switch (event.type) {
                case KeyRelease:
                    getIn(&event);
	                break;
                case KeyPress:
                    getIn(&event);
	                break;
            }
            break;
        }
	    frameCount++;
        //printf(":%i %i:", keyDown[0], keyDown[2]);
        //printf("cam pos z: %f", cam.pos.z);
    } while (1);

    tcsetattr( STDIN_FILENO, TCSANOW, &told);

    clear();
    return 0;  
} 
