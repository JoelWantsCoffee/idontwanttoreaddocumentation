#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include "vecmath.h"
#include "surfaces.h"

int main() {
    surface surf;
    initSurf(&surf, 100, 50);
    background(&surf, mColour(255));
    draw(surf);
    return 0;
}