#define sign(a) (((a)<(0))?(-1):(1))

//structs
typedef struct _colour {
    int r;
    int g;
    int b;
} colour;

typedef struct _surface {
    colour lastCol;
    int width;
    int height;
    size_t pixelsLength;
    char cols [9];
    colour * pixels;
} surface;

//functions
void initSurf(surface *surf, int w, int h) {
    strcpy(surf->cols, " .,>xX#@");
    surf->pixelsLength = w*h;
    surf->width = w;
    surf->height = h;
    surf->pixels = malloc(w * h * sizeof(colour));
}

void freeSurf(surface *surf) {
    free(surf->pixels);
}

void draw(surface surf) { 
    for (int j = 0; j<surf.height; j++) {
        for (int i = 0; i<surf.width; i++) {
            colour col = surf.pixels[i + j*surf.width];
            float tone = ( ((float)col.r + col.g + col.b)/(3*256) );	
	        putchar(surf.cols[(int) floor(tone*(sizeof(surf.cols) - 1)) ]);
        }
        putchar('\n');
    }
    fflush(stdout);
}

void drawStipple(surface surf, int res) { 
    for (int j = 0; j<surf.height; j++) {
        for (int i = 0; i<surf.width; i++) {
            colour col = surf.pixels[i + j*surf.width];
            int tone = (int) floor( (((float)col.r + col.g + col.b)/(3*256) - 0.5)*res);
            int num = (i + j) % MAX(abs(tone), 2);
            num = (num>0); 
            if (num - (tone < 0)) {
                putchar(surf.cols[(int) floor(( ((float)col.r + col.g + col.b)/(3*256) )*(sizeof(surf.cols) - 1)) ]);
            } else {
                putchar(surf.cols[0]);
            }
        }
        putchar('\n');
    }
    fflush(stdout);
}

void setCol(surface *surf, colour col) {
    surf->lastCol = col;
}

void set(surface *surf, int x, int y) {
    if ((x >= 0) && (y >= 0) && (x < surf->width) && (y < surf->height)) {
        surf->pixels[(int) (x + (((int) y) * surf->width))] = surf->lastCol;
    }
}

void line(surface *surf, float x1, float y1, float x2, float y2) {
    float xmove = (x2 - x1);
    float ymove = (y2 - y1);

    int dist;

    if ( abs(xmove) > abs(ymove)) {
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
        set(surf, tx, ty);
        tx += xmove;
        ty += ymove;
    }
}

void triangle(surface *surf, float x1, float y1, float x2, float y2, float x3, float y3) {
    line(surf, x1, y1, x2, y2);
    line(surf, x2, y2, x3, y3);
    line(surf, x3, y3, x1, y1);
}

void fillTriangle(surface *surf, float x1, float y1, float x2, float y2, float x3, float y3) {
    line(surf, x1, y1, x2, y2);
    line(surf, x2, y2, x3, y3);
    line(surf, x3, y3, x1, y1);
}

void background(surface *surf, colour col) {
    for (int i = 0; i<surf->pixelsLength; i++) {
        surf->pixels[i] = col;
    }
}

colour fColour(int r, int g, int b) {
    colour col = {
        .r = r,
        .g = g,
        .b = b,
    };
    return col;
}

colour mColour(int v) {
    colour col = {
        .r = v,
        .g = v,
        .b = v,
    };
    return col;
}