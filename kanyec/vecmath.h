#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//Structs
typedef struct _mat {
    float m[4][4];
} mat;

typedef struct _vector {
    float x;
    float y;
    float z;
    float w;
} vector;

typedef struct _tri {
    vector p1;
    vector p2;
    vector p3;
} tri;

typedef struct _face {
    vector *p1;
    vector *p2;
    vector *p3;
} face;

typedef struct _mesh {
    vector * pts;
    size_t ptCount;
    face * faces;
    size_t faceCount;
} mesh;


//Functions
void initMat(mat * in, float scale) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i == j) {
                in->m[i][j] = scale;
            } else {
                in->m[i][j] = 0;
            }
        }
    }
}

void initMesh(mesh * in, int ptLen, int faLen) {
    in->ptCount = ptLen;
    in->faceCount = faLen;
    in->pts = malloc(sizeof(vector) * ptLen);
    in->faces = malloc(sizeof(face) * faLen);
}

void freeMesh(mesh * me) {
    free(me->pts);
    free(me->faces);
}

mat matmultmat(mat one, mat two) {
    mat out;
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

vector vectormultmat(vector vec, mat ma) {
    vector out = {
        .x = vec.x * ma.m[0][0] + vec.y * ma.m[1][0] + vec.z * ma.m[2][0] + ma.m[3][0],
        .y = vec.x * ma.m[0][1] + vec.y * ma.m[1][1] + vec.z * ma.m[2][1] + ma.m[3][1],
        .z = vec.x * ma.m[0][2] + vec.y * ma.m[1][2] + vec.z * ma.m[2][2] + ma.m[3][2],
        .w = vec.x * ma.m[0][3] + vec.y * ma.m[1][3] + vec.z * ma.m[2][3] + ma.m[3][3],
    };
    return out;
}

vector crossProduct(vector v1, vector v2) {
    vector out = {
        .x = v1.y*v2.z - v1.z*v2.y,
        .y = v1.z*v2.x - v1.x*v2.z,
        .z = v1.x*v2.y - v1.y*v2.x,
    };
    return out;
}

float dotProduct(vector v1, vector v2) {
    return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

vector vecSubtract(vector v1, vector v2) {
    vector out = {
        .x = v1.x - v2.x,
        .y = v1.y - v2.y,
        .z = v1.z - v2.z,
    };
    return out;
}

vector vecAdd(vector v1, vector v2) {
    vector out = {
        .x = v1.x + v2.x,
        .y = v1.y + v2.y,
        .z = v1.z + v2.z,
    };
    return out;
}

vector scale(vector v, float s) {
    vector out = {
        .x = v.x * s,
        .y = v.y * s,
        .z = v.z * s,	
    };
    return out;
}

void vecNormalise(vector * v) {
    float largestValue = MAX( abs(v->x), abs(v->y));
    largestValue = MAX(largestValue, abs(v->z));
    v->x /= largestValue;
    v->y /= largestValue;
    v->z /= largestValue;
}

tri ftot(face in) {
    tri out = {
        .p1 = *in.p1,
        .p2 = *in.p2,
        .p3 = *in.p3,
    };
    return out;
}

//second degree functions
vector vectormultmatScaled(vector vec, mat ma) {
    vector out = vectormultmat(vec, ma);
    out = scale(out, 1/(out.w));
    return out;
}

vector getNormal(tri in) {
    vector v1 = vecSubtract(in.p1, in.p2);
    vector v2 = vecSubtract(in.p1, in.p3);
    return crossProduct(v1, v2);
}

tri trianglemultmat(tri tr, mat ma) {
    tri out = {
        .p1 = vectormultmatScaled(tr.p1, ma),
        .p2 = vectormultmatScaled(tr.p2, ma),
        .p3 = vectormultmatScaled(tr.p3, ma),
    };
    return out;
}

void meshmultmat(mesh *me, mat ma) {
    for (int i = 0; i<me->ptCount; i++) me->pts[i] = vectormultmatScaled(me->pts[i], ma);
}

float magnitude(vector in) {
    return sqrt(in.x*in.x + in.y*in.y + in.z*in.z);
}

float sudoTheta(vector v1, vector v2) {
    return (dotProduct(v1, v2)/(magnitude(v1), magnitude(v2)));
}

float distance(vector p, vector n, vector *v) {
    return (n.x)*(v->x) + (n.y)*(v->y) + (n.z)*(v->z) - dotProduct(n, p);
}

vector vecIntersectPlane(vector p1, vector p2, vector *l1, vector *l2) {
	vecNormalise(&p2);
	float pd = -dotProduct(p2, p1);
	float ad = dotProduct(*l1, p2);
	float bd = dotProduct(*l2, p2);
	float t = (-pd - ad) / (bd - ad);
	vector wl = vecSubtract(*l2, *l1);
	vector pl = scale(wl, t);
	return vecAdd(*l1, pl);
}

void translateMesh(mesh *me, float x, float y, float z) {
    mat m;
    initMat(&m, 1);
    m.m[3][0] = x;
    m.m[3][1] = y;
    m.m[3][2] = z;
    meshmultmat(me, m);
}

void rotateMesh(mesh *me, float xrot, float yrot, float zrot) {
    mat m;
    initMat(&m, 1);
    m.m[1][1] = cosf(xrot);
    m.m[2][2] = cosf(xrot);
    m.m[2][1] = -sinf(xrot);
    m.m[1][2] = sinf(xrot);
    meshmultmat(me, m);
    initMat(&m, 1);
    m.m[0][0] = cosf(yrot);
    m.m[2][2] = cosf(yrot);
    m.m[2][0] = sinf(yrot);
    m.m[0][2] = -sinf(yrot);
    meshmultmat(me, m);
    initMat(&m, 1);
    m.m[0][0] = cosf(zrot);
    m.m[1][0] = -sinf(zrot);
    m.m[1][1] = cosf(zrot);
    m.m[0][1] = sinf(zrot);
    meshmultmat(me, m);
}

int clipFace(vector pl, vector pn, face in, face *out [2], vector *op [2], int *nPts) {
    vecNormalise(&pn);

    vector *ipts [3];
    vector *opts [3];
    int ic = 0;
    int oc = 0;

    float d0 = distance(pl, pn, in.p1);
    float d1 = distance(pl, pn, in.p2);
    float d2 = distance(pl, pn, in.p3);

    if (d0 >= 0) {
	    ipts[ic++] = in.p1;
    } else {
	    opts[oc++] = in.p1;
    }

    if (d1 >= 0) {
    	ipts[ic++] = in.p2;
    } else {
	    opts[oc++] = in.p2;
    }
    
    if (d2 >= 0) {
    	ipts[ic++] = in.p3;
    } else {
	    opts[oc++] = in.p3;
    }


    if (ic == 0) {
        *nPts = 0;
	    return 0;  
    } else if (ic == 3) {
        *out[0] = in;
        *nPts = 0;

        return 1;
    } else if (ic == 2) {
        *op[0] = vecIntersectPlane(pl,pn,ipts[0],opts[0]);
        *op[1] = vecIntersectPlane(pl,pn,ipts[1],opts[0]);


        out[0]->p1 = ipts[0];
        out[0]->p2 = op[0];
        out[0]->p3 = ipts[1];

        out[1]->p1 = ipts[1];
        out[1]->p2 = out[0]->p2;
        out[1]->p3 = op[1];

        *nPts = 2;

        return 2;
    } else if (ic == 1) {
        *op[0] = vecIntersectPlane(pl,pn,ipts[0],opts[0]);
        *op[1] = vecIntersectPlane(pl,pn,ipts[0],opts[1]);

        out[0]->p1 = ipts[0];
        out[0]->p2 = op[0];
        out[0]->p3 = op[1];

        *nPts = 2;

        return 1;
    }
}

void clipMesh(mesh *me, vector pl, vector pn) {
    mesh out;
    int ptCount = 0;
    int faCount = 0;
    face o [2];
    face *po [2];
    po[0] = o;
    po[1] = o+1;
    vector p [2];
    vector *pp [2];
    pp[0] = p;
    pp[1] = p+1;
    int t = 0;
    for (int i = 0; i<me->faceCount; i++) {
        faCount += clipFace(pl, pn, me->faces[i], po, pp, &t);
        ptCount += t;
    }

    initMesh(&out, me->ptCount + ptCount, me->faceCount + faCount);

    ptCount = 0;
    faCount = 0;

    for (int i = 0; i<me->ptCount; i++) out.pts[i] = me->pts[i];

    for (int i = 0; i<me->faceCount; i++) {
        pp[0] = &out.pts[me->ptCount + ptCount];
        pp[1] = &out.pts[me->ptCount + ptCount + 1];

        po[0] = &out.faces[faCount];
        po[1] = &out.faces[faCount + 1];

        int faCountInc = clipFace(pl, pn, me->faces[i], po, pp, &t);

        faCount += faCountInc;
        ptCount += t;
    }

    *me = out;
    freeMesh(&out);
}