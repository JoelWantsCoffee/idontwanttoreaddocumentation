//#include "surfaces.h"
//#include "vecmat.h"
#define sign(a) (((a)<(0))?(-1):(1))
//Structs
typedef struct _camera {
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
void wipe(camera *cam) {
    background(&(cam->dists), mColour(cam->drawDist));
    background(&(cam->image), mColour(0));
}

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

    wipe(cam);
}

void freeCamera(camera *cam) {
    freeSurf(&cam->image);
    freeSurf(&cam->dists);
}

void getProjectionMat(camera cam, mat * out) {
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
        if ((cam->dists.pixels[((int) x) + ((int) y)*w].r) > z) {
            cam->dists.pixels[((int) x) + ((int) y)*h].r = z;
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

void fillTri3d(camera * cam, float x1, float y1, float z1, float x2, float y2, float z2,float x3, float y3, float z3, colour col) {
    vector vecs [3];
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

        int top = 0;
    int bottom = 0;
    for (int i = 0; i<3;i++) {
        if (vecs[i].y > vecs[top].y) top = i;
        if (vecs[i].y < vecs[bottom].y) bottom = i;
    }
    
    float z = (z1 + z2 + z3)/3;

    int middle = 0;
    
    for (int i = 0; i<3; i++) {
        if ((bottom != i) && (top != i)) middle = i;
    }
    
    vector v1 = vecSubtract(vecs[bottom], vecs[middle]);
    vector v2 = vecSubtract(vecs[bottom], vecs[top]);
    
    if (!floor(v1.y)) v1.y = 1;
    if (!floor(v2.y)) v2.y = 1;

    float v1a = v1.x/v1.y;
    float v2a = v2.x/v2.y;
    
    vector np = {
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
            set3d(cam, j, y, z, col);
        }
    }
    
    vector vv1 = vecSubtract(vecs[middle], vecs[top]);
    vector vv2 = vecSubtract(np, vecs[top]);
    
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
        set3d(cam, j, y, z, col);
        }
    }
    }
}

void tri3d(camera *cam, tri t) {
    float charRatio = cam->charRatio;
    float w = (float) cam->width;
    float h = (float) cam->height;
    float x1 = (t.p1.x * w/2) * charRatio + w/2;
    float y1 = t.p1.y * h/2 + h/2;
    float x2 = (t.p2.x * w/2 ) * charRatio + w/2;
    float y2 = t.p2.y * h/2 + h/2;
    float x3 = (t.p3.x * w/2) * charRatio + w/2;
    float y3 = t.p3.y * h/2 + h/2;

    //fillTri3d(cam, x1, y1, t.p1.z, x2, y2, t.p2.z, x3, y3, t.p3.z, mColour(t.var));

    
    line3d(cam, x1, y1, t.p1.z, x2, y2, t.p2.z);
    line3d(cam, x2, y2, t.p2.z, x3, y3, t.p3.z);
    line3d(cam, x3, y3, t.p3.z, x1, y1, t.p1.z);
}

void clipMeshToCam(mesh *me) {
    vector pl1 = {.x = -0.5, .y = 0, .z = 0};
    vector pn1 = {.x = 1, .y = 0, .z = 0};
    vector pl2 = {.x = 0.5, .y = 0, .z = 0};
    vector pn2 = {.x = -1, .y = 0, .z = 0};
    vector pl3 = {.x = 0, .y = 1, .z = 0};
    vector pn3 = {.x = 0, .y = -1, .z = 0};
    vector pl4 = {.x = 0, .y = -1, .z = 0};
    vector pn4 = {.x = 0, .y = 1, .z = 0};
    clipMesh(me, pl1, pn1);
    clipMesh(me, pl2, pn2);
    clipMesh(me, pl3, pn3);
    clipMesh(me, pl4, pn4);
}

void mesh3d(camera *cam, mesh me) {
    vector l = {.x = 1, .y = 1, .z = -1};
    vecNormalise(&l);

    for (int i = 0; i<me.faceCount; i++) me.faces[i].var = (sudoTheta(getNormal(ftot(me.faces[i])), l) + 1)*0.5*255; 

    translateMesh(&me, cam->pos.x, cam->pos.y, -cam->pos.z);
    rotateMesh(&me, cam->vrot, cam->hrot, 0);
    mat proj;
    vector pl = {.x = 0, .y = 0, .z = 0.1};
    vector pn = {.x = 0, .y = 0, .z = 1};

    

    for (int i = 0; i<me.faceCount; i++)  {
        if (sudoTheta(getNormal(ftot(me.faces[i])), *me.faces[i].p1) > 0) me.faces[i].var = -1;
    }

    getProjectionMat(*cam, &proj);
    //clipMesh(&me, pl, pn);
    meshmultmat(&me, proj);
    clipMeshToCam(&me);
    vector zero = {.x = 0, .y = 0, .z = 0};
    for (int i = 0; i<me.faceCount; i++) {
        if (me.faces[i].var >= 0) tri3d(cam, ftot(me.faces[i]));
    }
}

int pullSubString(char * in, int inLen, char * out, char b, int index) {
    /*printf("\n---------------stat--------------\n");
    printf("Char: (%c)\n", b);
    printf("In len: (%d)\n", inLen);
    printf("Index: (%d)\n", index);
    printf("in: %s", in);*/
    int c = 0;
    for (int i = 0; i<inLen; i++) {
        //printf("(%c vs %c)\n", in[i], b);
        if (c == index) {
            int j;
            for (j = 0; (!(in[i + j] == b) && (j+i+1 < inLen)); j++) {
                out[j] = in[i+j];
            }
            out[j] = '\0';
            return j+1;
        }

        if (in[i] == b) c++;
    }
    return 0;
}

void inportObj(mesh *me, char * fileName) {
    FILE *fp;

    fp = fopen(fileName, "r");

    int len = 0;

    char ch;

    while((ch = fgetc(fp)) != EOF) {
        len++;
    }

    char t [len + 1];
    int c = 0;

    rewind(fp);

    while((ch = fgetc(fp)) != EOF) {
        t[c] = ch;
        c++;
    }

    fclose(fp);

    if (t[c-1] == '\n') c--;

    t[c] = '\0';
    int ptCount = 0;
    int faCount = 0;

    int lCount = 0;

    for (int c = 0; c<len+1; c++) {
        if (t[c] == '\n') lCount++;
    }

    for (int i = 0; i<lCount+1; i++) {
        char line [100];
        int l = pullSubString(t, len, line, '\n', i);
        if (l > 2) {
            if ((line[0] == 'v') && (line[1] = ' ')) ptCount++;
            if ((line[0] == 'f') && (line[1] = ' ')) faCount++;
        }
    }

    //printf("pt count: %d, f count %d", ptCount, faCount);

    initMesh(me, ptCount, faCount);

    //printf("whot %d", lCount);
    //fflush(stdout);

    int ptTally = 0;
    int faTally = 0;

    for (int i = 0; i<lCount+1; i++) {


        char line [100];
        int l = pullSubString(t, len, line, '\n', i);
        //printf("str: (%s)\n", line);
        //fflush(stdout);

        if (line[0] == 'v') {
            //printf("tally: %d", ptTally);
            char num [100];
            pullSubString(line, l, num, ' ', 1);
            me->pts[ptTally].x = (float) atof(num);
            pullSubString(line, l, num, ' ', 2);
            me->pts[ptTally].y = (float) atof(num);
            pullSubString(line, l, num, ' ', 3);
            me->pts[ptTally].z = (float) atof(num);
            me->pts[ptTally].w = 0;
            ptTally++;
        }
        if (line[0] == 'f') {
            char num [100];
            pullSubString(line, l, num, ' ', 1);
            me->faces[faTally].p1 = &(me->pts[atoi(num) - 1]);
            pullSubString(line, l, num, ' ', 2);
            me->faces[faTally].p2 = &(me->pts[atoi(num) - 1]);
            pullSubString(line, l, num, ' ', 3);
            me->faces[faTally].p3 = &(me->pts[atoi(num) - 1]);
            faTally++;
        }
    }
}

