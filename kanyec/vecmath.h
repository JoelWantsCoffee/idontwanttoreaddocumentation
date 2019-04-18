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
void matInit(mat * in, float scale) {
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
        .y = vec.y * ma.m[0][1] + vec.y * ma.m[1][1] + vec.z * ma.m[2][1] + ma.m[3][1],
        .z = vec.x * ma.m[0][2] + vec.y * ma.m[1][2] + vec.z * ma.m[2][2] + ma.m[3][2],
        .w = vec.y * ma.m[0][3] + vec.y * ma.m[1][3] + vec.z * ma.m[2][3] + ma.m[3][3],
    };
    return out;
}

vector crossProduct(vector v1, vector v2) {
    vector out;

    out.x = v1.y*v2.z - v1.z*v2.y;
    out.y = v1.z*v2.x - v1.x*v2.z;
    out.z = v1.x*v2.y - v1.y*v2.x;

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
    for (int i = 0; i<me->ptCount; i++) me->pts[i] = vectormultmat(me->pts[i], ma);
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

void rotateMesh(mesh *me, float xrot, float yrot, float zrot) {
    mat m;
    matInit(&m, 1);
    m.m[1][1] = cosf(xrot);
    m.m[2][2] = cosf(xrot);
    m.m[2][1] = -sinf(xrot);
    m.m[1][2] = sinf(xrot);
    meshmultmat(me, m);
    matInit(&m, 1);
    m.m[0][0] = cosf(yrot);
    m.m[2][2] = cosf(yrot);
    m.m[2][0] = sinf(yrot);
    m.m[0][2] = -sinf(yrot);
    meshmultmat(me, m);
    matInit(&m, 1);
    m.m[0][0] = cosf(zrot);
    m.m[1][0] = -sinf(zrot);
    m.m[1][1] = cosf(zrot);
    m.m[0][1] = sinf(zrot);
    meshmultmat(me, m);
}