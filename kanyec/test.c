#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include "vecmath.h"
#include "surfaces.h"
#include "cam3d.h"
#include "running.h"

int main() {
    camera c;
    initCamera(&c, 80, 40);

    mesh me;
    
    initMesh(&me, 3, 1);

    {
        me.pts[0].x = -0.4;
        me.pts[0].y = -0.4;
        me.pts[0].z = 10;
        me.pts[1].x = 0.4;
        me.pts[1].y = -0.4;
        me.pts[1].z = 10;
        me.pts[2].x = 0.4;
        me.pts[2].y = 0.4;
        me.pts[2].z = 10;

        me.faces[0].p1 = &me.pts[0]; 
        me.faces[0].p2 = &me.pts[1]; 
        me.faces[0].p3 = &me.pts[2]; 
    }

    mesh3d(&c, me);
    
    draw(c.image);
    return 0;
}