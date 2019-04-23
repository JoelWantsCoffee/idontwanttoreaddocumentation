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

void loop(input *ins, int fc, int arg) {
    camera c;
    initCamera(&c, 80, 40);

    mesh me;
    
    initMesh(&me, 3, 1);

    {
        me.pts[0].x = -0.2;
        me.pts[0].y = -0.3;
        me.pts[0].z = 1;
        me.pts[1].x = 0.4;
        me.pts[1].y = -0.4;
        me.pts[1].z = 1;
        me.pts[2].x = 0.2;
        me.pts[2].y = 0.4;
        me.pts[2].z = 1;

        me.faces[0].p1 = &me.pts[0]; 
        me.faces[0].p2 = &me.pts[1]; 
        me.faces[0].p3 = &me.pts[2]; 
    }
    
    if (ins->on[0]) {
        c.pos.z = -0.1;
    }
    if (ins->on[2]) {
       c.pos.z = 0.1;
    }
    
    mesh3d(&c, me);

    draw(c.image);

    clearIns(ins);
}

int main() {
    liveInputs();
    shell sh;
    initShell(&sh, &loop);
    while (1) {
        getIns(&sh.ins);
        run(&sh, 0, 0, 30);
    }
    return 0;
}