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

void loop(input *ins, int fc, void * arg) {
    camera c;
    initCamera(&c, 400, 100);
    mesh me;
    inportObj(&me, "cube.obj");

    if (ins->on[0]) ((vector*) arg)->z += 0.1;
    if (ins->on[1]) ((vector*) arg)->x += 0.1;
    if (ins->on[2]) ((vector*) arg)->z -= 0.1;
    if (ins->on[3]) ((vector*) arg)->x -= 0.1;

    if (ins->on[4]) ((vector*) arg)->y += 0.1;
    if (ins->on[6]) ((vector*) arg)->y -= 0.1;

    c.pos.z = ((vector*) arg)->z;
    c.pos.x = ((vector*) arg)->x;
    c.pos.y = ((vector*) arg)->y;

    mesh3d(&c, me);

    draw(c.image);

    freeCamera(&c);
    freeMesh(&me);

    clearIns(ins);
}

int main() {
    mesh me;

    vector pos = {
        .x = 0,
        .y = 0,
        .z = 0,
    };

    liveInputs();
    shell sh;
    initShell(&sh, &loop);

    while (1) {
        getIns(&sh.ins);
        run(&sh, 0, &pos, 30);
    }
    return 0;
}