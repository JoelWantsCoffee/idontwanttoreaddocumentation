#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <termios.h>
#include <time.h>
#include <unistd.h>

#include "vecmath.h"
#include "surfaces.h"
#include "cam3d.h"
#include "shells.h"

void loop(input *ins, int fc, void ** arg) {
    camera c;
    initCamera(&c, 150, 50);
    mesh me;
    //inportObj(&me, "scene.obj");

    dupeMesh((mesh *) arg[1], &me);

    //rotateMesh(&me, 0, 3.1415, 3.1415);

    vector mv = {
        .x = 0,
        .y = 0,
        .z = 0,
        .w = 0,
    };

    
    if (ins->on[4]) ((vector*) arg)->y += 0.1;
    if (ins->on[6]) ((vector*) arg)->y -= 0.1;

    

    if (ins->on[5]) ((vector*) arg)->w += 0.1;
    if (ins->on[7]) ((vector*) arg)->w -= 0.1;

    float hrot = ((vector*) arg)->w;

    if (ins->on[1] + ins->on[3]) hrot += 3.1415926/2;


    mv.z = cos(hrot)*0.1;
    mv.x = sin(hrot)*0.1;

    if (ins->on[1] + ins->on[3]) hrot -= 3.1415926/2;

    
    if (ins->on[2] + ins->on[3]) {
        mv = scale(mv, -1);
    }

    if (!(ins->on[0] + ins->on[1] + ins->on[2] + ins->on[3])) mv = scale(mv, 0);

    *((vector*) arg) = vecAdd(*((vector*) arg), mv);

    c.pos.z = ((vector*) arg)->z;
    c.pos.x = ((vector*) arg)->x;
    c.pos.y = ((vector*) arg)->y;
    c.hrot = ((vector*) arg)->w;

    mesh3d(&c, me);

    drawStipple(c.image, 20);

    freeCamera(&c);
    freeMesh(&me);

    clearIns(ins);
}

int main() {
    mesh me;

    vector pos = {
        .x = 0,
        .y = 0,
        .z = -5,
        .w = 0,
    };

    void * arg [2];

    inportObj(&me, "sphere.obj");

    arg[0] = &pos;
    arg[1] = &me;

    

    liveInputs();
    shell sh;
    initShell(&sh, &loop);

    while (1) {
        getIns(&sh.ins);
        run(&sh, 0, arg, 30);
    }
    return 0;
}