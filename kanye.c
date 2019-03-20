#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

//Global Variables
const int WIDTH = 80;
const int HEIGHT = 23;
const int pixelsLength = 80*23;
float pixels[80*23];
const char * cols = " .,:;$#";

struct mat {
    float m[4][4];
};

float sign(float in) {
    if (in > 0) {
        return 1;
    } else {
        return 0;
    }
}

void initOne(struct mat in) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i == j) {
                in.m[i][j] = 1;
            } else {
                in.m[i][j] = 0;
            }
        }
    }
}

void draw() {
    for (int j = 0; j<HEIGHT; j++) {
        for (int i = 0; i<WIDTH; i++) {
            float tone = pixels[i + j*WIDTH]/255;
            if  (tone < 0.14) {
                printf(" ");
            } else if (tone < 0.28) {
                printf(".");
            } else if (tone < 0.42) {
                printf(",");
            } else if (tone < 0.56) {
                printf(":");
            } else if (tone < 0.7) {
                printf(";");
            } else if (tone < 0.84) {
                printf("$");
            } else {
                printf("#");
            }
        }
        printf("\n");
    }
}

void set(float x, float y, float col) {
    if ((x < WIDTH) && (y < HEIGHT)) {
        pixels[(int) (x + (((int) y) * WIDTH))] = col;
    }
}

void line(float x1, float y1, float x2, float y2) {
    float xmove = (x2 - x1);
    float ymove = (y2 - y1);

    float ax = abs(xmove);
    float ay = abs(ymove);

    if (ax > ay) {
        ymove = ymove / abs(xmove);
        xmove = 1;
    } else {
        xmove = xmove / abs(ymove);
        ymove = 1;
    }

    float tx = x1;
    float ty = y1;

    do {
        set(tx, ty, 255);
        tx += xmove;
        ty += ymove;
    } while ((((int) tx) != ((int)x2)) || ((int)ty) != ((int) y2));
}

void triangle(float x1, float y1, float x2, float y2, float x3, float y3) {
    line(x1, y1, x2, y2);
    line(x2, y2, x3, y3);
    line(x3, y3, x1, y1);
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

void loop() {
    background(0);
    triangle(10, 2, 20, 18, 30, 10);
}

int main()  { 
    struct mat m;
    initOne(m);
    do {
        clear();
        loop();
        draw();
        wait(1000);
    } while (1);
    clear();
    return 0;  
} 
