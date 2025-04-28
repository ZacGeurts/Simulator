#include "renderer.h"
#include "equations.h"
#include <GL/glew.h>
#include <SDL2/SDL_opengl.h>
#include <GL/glu.h>
#include <SDL2/SDL_ttf.h>
#include <cmath>
#include <string>
#include <algorithm>
#include <vector>
#include <omp.h>

// Marching cubes edge table (unchanged)
static const int edgeTable[256] = {
    0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa, 0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66, 0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55, 0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0xcc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x55, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0xff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x66, 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa, 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33, 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99, 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
};

// Marching cubes triangle table (unchanged)
static const int triTable[256][16] = {
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 4, 5, 2, 4, 2, 11, 2, 3, 11, -1, -1, -1, -1},
    {1, 9, 4, 1, 4, 5, 2, 1, 5, 2, 5, 8, 2, 8, 11, -1},
    {3, 11, 10, 3, 10, 1, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 4, 9, 1, 0, 8, 1, 8, 10, 10, 8, 11, -1, -1, -1, -1},
    {4, 0, 3, 4, 3, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1},
    {11, 4, 9, 11, 9, 10, 8, 11, 10, 8, 10, 5, 8, 5, 4, -1},
    {5, 7, 8, 5, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 7, 4, 4, 7, 3, 4, 3, 0, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 5, 7, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 4, 1, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 5, 7, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 7, 4, 5, 3, 0, 1, 3, 1, 2, 3, 2, 10, -1},
    {0, 2, 4, 4, 2, 10, 4, 10, 5, 5, 10, 9, -1, -1, -1, -1},
    {7, 4, 8, 3, 2, 10, 3, 10, 5, 3, 5, 9, -1, -1, -1, -1},
    {2, 3, 11, 7, 8, 5, 7, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    {7, 4, 5, 3, 0, 11, 0, 2, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 9, 4, 5, 7, 4, 7, 8, -1, -1, -1, -1},
    {9, 4, 5, 7, 3, 1, 7, 1, 2, 7, 2, 11, -1, -1, -1, -1},
    {5, 7, 4, 3, 11, 10, 3, 10, 1, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 0, 10, 0, 11, 7, 4, 5, 7, 5, 8, -1, -1, -1, -1},
    {3, 11, 0, 0, 11, 10, 9, 0, 10, 4, 5, 7, 4, 7, 8, -1},
    {7, 4, 5, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 6, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 9, 10, 5, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 5, 2, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 2, 5, 0, 5, 6, 0, 6, 9, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 2, 8, 2, 6, 8, 6, 9, 6, 5, 9, -1, -1, -1, -1},
    {3, 11, 2, 10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 6, 10, 0, 11, 2, 0, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 9, 10, 5, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 5, 6, 1, 9, 8, 1, 8, 11, 1, 11, 2, -1, -1, -1, -1},
    {3, 6, 5, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {3, 6, 5, 0, 3, 5, 0, 5, 9, 0, 9, 11, -1, -1, -1, -1},
    {9, 8, 11, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 6, 4, 4, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 5, 6, 4, 6, 7, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 4, 5, 6, 4, 6, 7, -1, -1, -1, -1, -1, -1, -1},
    {4, 5, 6, 4, 6, 7, 8, 3, 1, 8, 1, 9, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, 4, 5, 6, 4, 6, 7, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 5, 2, 6, 5, 4, 5, 6, 4, 6, 7, -1},
    {4, 5, 6, 4, 6, 7, 0, 2, 9, 2, 5, 9, -1, -1, -1, -1},
    {4, 5, 6, 4, 6, 7, 3, 2, 9, 3, 9, 8, -1, -1, -1, -1},
    {3, 11, 2, 4, 5, 6, 4, 6, 7, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 4, 5, 6, 4, 6, 7, 2, 0, 11, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 11, 4, 5, 6, 4, 6, 7, -1, -1, -1, -1},
    {4, 5, 6, 4, 6, 7, 2, 1, 9, 2, 9, 11, 11, 9, 8, -1},
    {3, 5, 1, 3, 6, 5, 4, 5, 6, 4, 6, 7, -1, -1, -1, -1},
    {4, 5, 6, 4, 6, 7, 0, 8, 11, 1, 0, 11, 1, 11, 5, -1},
    {0, 9, 11, 0, 3, 11, 4, 5, 6, 4, 6, 7, -1, -1, -1, -1},
    {4, 5, 6, 4, 6, 7, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 7, 10, 7, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 10, 6, 7, 10, 7, 11, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 11, 10, 6, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 9, 11, 10, 6, 11, 6, 7, -1, -1, -1, -1},
    {1, 2, 10, 6, 7, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 10, 6, 7, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 9, 0, 2, 10, 9, 6, 7, 11, -1, -1, -1, -1, -1, -1, -1},
    {6, 7, 11, 3, 2, 10, 3, 10, 9, 3, 9, 8, -1, -1, -1, -1},
    {2, 3, 7, 2, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 2, 6, 7, 2, 7, 11, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 7, 2, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 8, 1, 8, 11, 2, 1, 11, 2, 11, 6, 2, 6, 7, -1},
    {1, 3, 6, 1, 6, 10, 6, 7, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 1, 0, 11, 1, 11, 10, 6, 7, 11, 6, 11, 1, -1},
    {0, 9, 3, 9, 7, 3, 7, 6, 3, 10, 6, 7, -1, -1, -1, -1},
    {6, 7, 11, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 10, 7, 10, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 10, 0, 10, 6, 0, 6, 7, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 10, 6, 8, 6, 7, -1, -1, -1, -1, -1, -1, -1},
    {3, 1, 9, 3, 9, 10, 3, 10, 6, 3, 6, 7, -1, -1, -1, -1},
    {1, 2, 10, 8, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 1, 1, 3, 2, 2, 3, 10, 6, 7, 10, -1, -1, -1, -1},
    {0, 2, 9, 9, 2, 10, 8, 6, 7, 8, 7, 10, -1, -1, -1, -1},
    {3, 2, 10, 3, 10, 6, 3, 6, 7, 9, 3, 7, -1, -1, -1, -1},
    {2, 3, 7, 2, 7, 6, 8, 6, 7, -1, -1, -1, -1, -1, -1, -1},
    {0, 2, 7, 0, 7, 6, 0, 6, 8, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 7, 2, 7, 6, 8, 6, 7, -1, -1, -1, -1},
    {1, 9, 7, 1, 7, 6, 1, 6, 2, 2, 6, 3, -1, -1, -1, -1},
    {1, 3, 6, 1, 6, 10, 7, 8, 6, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 2, 0, 2, 10, 0, 10, 6, 0, 6, 7, -1, -1, -1, -1},
    {0, 9, 3, 3, 9, 7, 6, 3, 7, 10, 6, 7, 8, 6, 7, -1},
    {9, 7, 10, 7, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 0, 8, 3, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {0, 2, 9, 2, 10, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, 3, 2, 10, 3, 10, 9, 3, 9, 8, -1, -1, -1, -1},
    {2, 3, 7, 2, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 2, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 7, 2, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 8, 1, 8, 3, 2, 1, 3, 2, 3, 7, 2, 7, 6, -1},
    {1, 3, 6, 1, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 10, 6, 1, 0, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 9, 3, 9, 6, 3, 6, 10, 3, -1, -1, -1, -1, -1, -1, -1},
    {6, 10, 9, 6, 9, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 10, 7, 10, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 10, 0, 10, 6, 0, 6, 7, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 10, 6, 8, 6, 7, -1, -1, -1, -1, -1, -1, -1},
    {3, 1, 9, 3, 9, 10, 3, 10, 6, 3, 6, 7, -1, -1, -1, -1},
    {1, 2, 10, 8, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 1, 1, 3, 2, 2, 3, 10, 6, 7, 10, -1, -1, -1, -1},
    {0, 2, 9, 9, 2, 10, 8, 6, 7, 8, 7, 10, -1, -1, -1, -1},
    {3, 2, 10, 3, 10, 6, 3, 6, 7, 9, 3, 7, -1, -1, -1, -1},
    {2, 3, 7, 2, 7, 6, 8, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {0, 2, 7, 0, 7, 6, 0, 6, 8, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 7, 2, 7, 6, 8, 2, 7, -1, -1, -1, -1},
    {1, 9, 3, 1, 3, 7, 1, 7, 6, 1, 6, 2, -1, -1, -1, -1},
    {1, 3, 6, 1, 6, 10, 8, 6, 7, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 2, 0, 2, 10, 0, 10, 6, 0, 6, 7, -1, -1, -1, -1},
    {0, 9, 3, 3, 9, 7, 6, 3, 7, 10, 6, 7, 8, 6, 7, -1},
    {9, 7, 10, 7, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};

// Font and TTF management
static TTF_Font* font = nullptr;
static bool ttf_initialized = false;

// Initialize SDL_ttf
bool init_ttf() {
    if (TTF_Init() == -1) {
        log_message("SDL_ttf initialization failed: " + std::string(TTF_GetError()));
        return false;
    }
    ttf_initialized = true;

    // Try loading font from multiple possible paths
    const char* font_paths[] = {
        "../DejaVuSans.ttf",
        "DejaVuSans.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "C:/Windows/Fonts/DejaVuSans.ttf"
    };
    for (const char* path : font_paths) {
        font = TTF_OpenFont(path, 24);
        if (font) {
            return true;
        }
    }

    log_message("Failed to load DejaVuSans.ttf from any path: " + std::string(TTF_GetError()));
    TTF_Quit();
    ttf_initialized = false;
    return false;
}

// Clean up SDL_ttf
void cleanup_ttf() {
    if (font) {
        TTF_CloseFont(font);
        font = nullptr;
    }
    if (ttf_initialized) {
        TTF_Quit();
        ttf_initialized = false;
    }
}

static const int edge_to_vertices[12][2] = {
    {0,1}, {1,2}, {2,3}, {3,0}, {4,5}, {5,6}, {6,7}, {7,4}, {0,4}, {1,5}, {2,6}, {3,7}
};

static const float vertex_positions[8][3] = {
    {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}
};

// Render text to texture
GLuint render_text_to_texture(const std::string& text, SDL_Color color, int& w, int& h) {
    if (!font || text.empty()) {
        w = 0;
        h = 0;
        return 0;
    }

    // Render text to surface
    SDL_Surface* surface = TTF_RenderText_Blended(font, text.c_str(), color);
    if (!surface) {
        log_message("Failed to render text '" + text + "': " + std::string(TTF_GetError()));
        w = 0;
        h = 0;
        return 0;
    }

    // Ensure surface format is RGBA
    SDL_Surface* converted = SDL_ConvertSurfaceFormat(surface, SDL_PIXELFORMAT_RGBA32, 0);
    SDL_FreeSurface(surface);
    if (!converted) {
        log_message("Failed to convert surface to RGBA: " + std::string(SDL_GetError()));
        w = 0;
        h = 0;
        return 0;
    }

    w = converted->w;
    h = converted->h;

    // Create OpenGL texture
    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, converted->w, converted->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, converted->pixels);

    SDL_FreeSurface(converted);
    glBindTexture(GL_TEXTURE_2D, 0);
    return texture;
}

void render(const Simulation& sim) {
    static bool ttf_init_done = false;
    if (!ttf_init_done) {
        if (!init_ttf()) {
            log_message("Failed to initialize TTF. Text rendering disabled.");
        }
        ttf_init_done = true;
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, 800.0 / 600.0, 0.01, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Compute camera position with tilt
    float cam_dist = 5.0f / camera_zoom;
    // Initial position in x-y plane (rotation around z-axis)
    float cam_x = cam_dist * std::cos(camera_angle);
    float cam_y = cam_dist * std::sin(camera_angle);
    float cam_z = cam_dist * 0.5f;
    // Apply tilt (rotation around x-axis)
    float tilted_y = cam_y * std::cos(camera_tilt) - cam_z * std::sin(camera_tilt);
    float tilted_z = cam_y * std::sin(camera_tilt) + cam_z * std::cos(camera_tilt);
    gluLookAt(cam_x, tilted_y, tilted_z, 0, 0, 0, 0, 0, 1);

    // Compute scalar ranges
    double max_scalar = 1e-10, min_scalar = 1e10;
    double max_4d = 1e-10, min_4d = 1e10;
    for (int i = 0; i < sim.size; ++i) {
        for (int j = 0; j < sim.size; ++j) {
            for (int k = 0; k < sim.size; ++k) {
                max_scalar = std::max(max_scalar, sim.computed_scalar[i][j][k]);
                min_scalar = std::min(min_scalar, sim.computed_scalar[i][j][k]);
                for (int w = 0; w < sim.size; ++w) {
                    max_4d = std::max(max_4d, sim.scalar_4d[i][j][k][w]);
                    min_4d = std::min(min_4d, sim.scalar_4d[i][j][k][w]);
                }
            }
        }
    }
    double range = max_scalar - min_scalar;
    if (range == 0) range = 1.0;
    double range_4d = max_4d - min_4d;
    if (range_4d == 0) range_4d = 1.0;
    double iso_level = min_scalar + range * (0.3 + 0.2 * std::sin(sim.t));

    int w_slice = static_cast<int>(sim.t * 2) % sim.size;

    // Get visualization parameters from current equation
    float vis_scale = 1.0f;
    float vis_color_intensity = 1.0f;
    if (!sim.equations.empty()) {
        if (auto* eq = dynamic_cast<CustomEquation*>(sim.equations[sim.current_equation])) {
            vis_scale = static_cast<float>(eq->getVisScale());
            vis_color_intensity = static_cast<float>(eq->getVisColorIntensity());
        }
    }

    // Render particles
    if (current_display_mode == PARTICLES || current_display_mode == HYBRID) {
        glPointSize(8.0);
        glBegin(GL_POINTS);
        for (const auto& p : sim.particles) {
            float speed = std::sqrt(p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
            float r = std::min(speed / 5.0f, 1.0f) * vis_color_intensity;
            float g = 0.5f * vis_color_intensity;
            float b = (1.0f - r) * vis_color_intensity;
            glColor3f(r, g, b);
            glVertex3f(p.x * vis_scale, p.y * vis_scale, p.z * vis_scale);
        }
        glEnd();
        glBegin(GL_LINES);
        glColor3f(0.5f, 0.5f, 1.0f);
        for (const auto& p : sim.particles) {
            glVertex3f(p.x * vis_scale, p.y * vis_scale, p.z * vis_scale);
            glVertex3f((p.x - p.vx * sim.dt * 10) * vis_scale, (p.y - p.vy * sim.dt * 10) * vis_scale, (p.z - p.vz * sim.dt * 10) * vis_scale);
        }
        glEnd();
    }

    // Render grid-based data
    if (current_display_mode != PARTICLES && sim.size > 0) {
        if (current_display_mode == POINTS || current_display_mode == HYBRID) {
            glPointSize(6.0);
            glBegin(GL_POINTS);
            for (int i = 0; i < sim.size; ++i) {
                for (int j = 0; j < sim.size; ++j) {
                    for (int k = 0; k < sim.size; ++k) {
                        float x = static_cast<float>((i - sim.size / 2) * sim.dx);
                        float y = static_cast<float>((j - sim.size / 2) * sim.dx);
                        float z = static_cast<float>((k - sim.size / 2) * sim.dx);
                        float value = (sim.computed_scalar[i][j][k] - min_scalar) / range;
                        float r = value * vis_color_intensity;
                        float g = (1.0f - value) * vis_color_intensity;
                        float b = 0.5f * vis_color_intensity;
                        glColor3f(r, g, b);
                        glVertex3f(x * vis_scale, y * vis_scale, z * vis_scale);
                    }
                }
            }
            glEnd();
        }

        if (current_display_mode == ISOSURFACE) {
            // Initialize GLEW if not already done
            static bool glew_initialized = false;
            if (!glew_initialized) {
                GLenum err = glewInit();
                if (err != GLEW_OK) {
                    log_message("GLEW initialization failed: " + std::string(reinterpret_cast<const char*>(glewGetErrorString(err))));
                    return;
                }
                glew_initialized = true;
            }

            // Validate simulation size and w_slice
            if (sim.size <= 0 || w_slice < 0 || w_slice >= sim.size) {
                log_message("Invalid simulation size or w_slice: size=" + std::to_string(sim.size) + ", w_slice=" + std::to_string(w_slice));
                return;
            }

            // Precompute vertex data
            std::vector<std::vector<std::vector<float>>> colors(sim.size, std::vector<std::vector<float>>(sim.size, std::vector<float>(sim.size)));
            std::vector<std::vector<std::vector<std::array<float, 3>>>> normals(sim.size, std::vector<std::vector<std::array<float, 3>>>(sim.size, std::vector<std::array<float, 3>>(sim.size)));
            #pragma omp parallel for
            for (int i = 0; i < sim.size; ++i) {
                for (int j = 0; j < sim.size; ++j) {
                    for (int k = 0; k < sim.size; ++k) {
                        colors[i][j][k] = (sim.scalar_4d[i][j][k][w_slice] - min_4d) / range_4d;
                        float nx = 0, ny = 0, nz = 0;
                        if (i > 0 && i < sim.size - 1)
                            nx = (sim.computed_scalar[i + 1][j][k] - sim.computed_scalar[i - 1][j][k]) / (2 * sim.dx);
                        if (j > 0 && j < sim.size - 1)
                            ny = (sim.computed_scalar[i][j + 1][k] - sim.computed_scalar[i][j - 1][k]) / (2 * sim.dx);
                        if (k > 0 && k < sim.size - 1)
                            nz = (sim.computed_scalar[i][j][k + 1] - sim.computed_scalar[i][j][k - 1]) / (2 * sim.dx);
                        float norm = std::sqrt(nx * nx + ny * ny + nz * nz);
                        normals[i][j][k] = {norm > 0 ? nx / norm : 0, norm > 0 ? ny / norm : 0, norm > 0 ? nz / norm : 0};
                    }
                }
            }

            // Vertex data for VBO, collected per thread
            std::vector<std::vector<float>> thread_vertices(omp_get_max_threads());
            const int edge_to_vertices[12][2] = {
                {0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}
            };
            static int cube_error_count = 0;
            static int edge_error_count = 0;
            const int max_errors_per_frame = 5;

            #pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < sim.size - 1; ++i) {
                int thread_id = omp_get_thread_num();
                thread_vertices[thread_id].reserve(1000 * 9);
                for (int j = 0; j < sim.size - 1; ++j) {
                    for (int k = 0; k < sim.size - 1; ++k) {
                        float values[8] = {
                            static_cast<float>(sim.computed_scalar[i][j][k]),
                            static_cast<float>(sim.computed_scalar[i + 1][j][k]),
                            static_cast<float>(sim.computed_scalar[i + 1][j + 1][k]),
                            static_cast<float>(sim.computed_scalar[i][j + 1][k]),
                            static_cast<float>(sim.computed_scalar[i][j][k + 1]),
                            static_cast<float>(sim.computed_scalar[i + 1][j][k + 1]),
                            static_cast<float>(sim.computed_scalar[i + 1][j + 1][k + 1]),
                            static_cast<float>(sim.computed_scalar[i][j + 1][k + 1])
                        };
                        float vertices[8][3] = {
                            {static_cast<float>((i - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k - sim.size / 2) * sim.dx * vis_scale)},
                            {static_cast<float>((i + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k - sim.size / 2) * sim.dx * vis_scale)},
                            {static_cast<float>((i + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k - sim.size / 2) * sim.dx * vis_scale)},
                            {static_cast<float>((i - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k - sim.size / 2) * sim.dx * vis_scale)},
                            {static_cast<float>((i - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k + 1 - sim.size / 2) * sim.dx * vis_scale)},
                            {static_cast<float>((i + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k + 1 - sim.size / 2) * sim.dx * vis_scale)},
                            {static_cast<float>((i + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k + 1 - sim.size / 2) * sim.dx * vis_scale)},
                            {static_cast<float>((i - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((j + 1 - sim.size / 2) * sim.dx * vis_scale), static_cast<float>((k + 1 - sim.size / 2) * sim.dx * vis_scale)}
                        };
                        int vi[8] = {i, i + 1, i + 1, i, i, i + 1, i + 1, i};
                        int vj[8] = {j, j, j + 1, j + 1, j, j, j + 1, j + 1};
                        int vk[8] = {k, k, k, k, k + 1, k + 1, k + 1, k + 1};

                        // Compute cube index
                        int cube_index = 0;
                        for (int v = 0; v < 8; ++v) {
                            if (values[v] >= iso_level) cube_index |= (1 << v);
                        }
                        if (cube_index < 0 || cube_index > 255) {
                            #pragma omp critical
                            {
                                if (cube_error_count < max_errors_per_frame) {
                                    log_message("Invalid cube_index: " + std::to_string(cube_index) + " at i=" + std::to_string(i) + ", j=" + std::to_string(j) + ", k=" + std::to_string(k));
                                    cube_error_count++;
                                }
                            }
                            continue;
                        }

                        // Check if cube contributes triangles using edgeTable
                        if (edgeTable[cube_index] == 0) {
                            continue;
                        }

                        // Process triangles
                        for (int t = 0; t < 16 && triTable[cube_index][t] != -1; t += 3) {
                            for (int v = 0; v < 3; ++v) {
                                int edge = triTable[cube_index][t + v];
                                if (edge < 0 || edge >= 12) {
                                    #pragma omp critical
                                    {
                                        if (edge_error_count < max_errors_per_frame) {
                                            log_message("Invalid edge index: " + std::to_string(edge) + " at cube_index=" + std::to_string(cube_index) + ", t=" + std::to_string(t) + ", i=" + std::to_string(i) + ", j=" + std::to_string(j) + ", k=" + std::to_string(k));
                                            edge_error_count++;
                                        }
                                    }
                                    continue;
                                }
                                int v0 = edge_to_vertices[edge][0];
                                int v1 = edge_to_vertices[edge][1];
                                float t0 = 0.5f;
                                if (std::abs(values[v1] - values[v0]) > 1e-6) {
                                    t0 = (iso_level - values[v0]) / (values[v1] - values[v0]);
                                    t0 = std::max(0.0f, std::min(1.0f, t0));
                                }
                                float vert[3], norm[3];
                                float color = colors[vi[v0]][vj[v0]][vk[v0]] + t0 * (colors[vi[v1]][vj[v1]][vk[v1]] - colors[vi[v0]][vj[v0]][vk[v0]]);
                                for (int d = 0; d < 3; ++d) {
                                    vert[d] = vertices[v0][d] + t0 * (vertices[v1][d] - vertices[v0][d]);
                                    norm[d] = normals[vi[v0]][vj[v0]][vk[v0]][d] + t0 * (normals[vi[v1]][vj[v1]][vk[v1]][d] - normals[vi[v0]][vj[v0]][vk[v0]][d]);
                                }
                                float norm_len = std::sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
                                if (norm_len > 0) {
                                    norm[0] /= norm_len;
                                    norm[1] /= norm_len;
                                    norm[2] /= norm_len;
                                }
                                // Append to thread-local vertex data
                                thread_vertices[thread_id].insert(thread_vertices[thread_id].end(), {
                                    vert[0], vert[1], vert[2],
                                    norm[0], norm[1], norm[2],
                                    std::max(0.2f, color * vis_color_intensity), 0.5f, std::max(0.2f, (1.0f - color) * vis_color_intensity)
                                });
                            }
                        }
                    }
                }
            }

            // Reset error counts for next frame
            cube_error_count = 0;
            edge_error_count = 0;

            // Merge thread-local vertex data
            std::vector<float> vertex_data;
            size_t total_size = 0;
            for (const auto& tv : thread_vertices) {
                total_size += tv.size();
            }
            vertex_data.reserve(total_size);
            for (const auto& tv : thread_vertices) {
                vertex_data.insert(vertex_data.end(), tv.begin(), tv.end());
            }

            // Setup VBO
            static GLuint vbo = 0;
            if (vbo == 0) {
                glGenBuffers(1, &vbo);
            }
            glBindBuffer(GL_ARRAY_BUFFER, vbo);
            glBufferData(GL_ARRAY_BUFFER, vertex_data.size() * sizeof(float), vertex_data.data(), GL_STATIC_DRAW);

            // Setup lighting
            glEnable(GL_LIGHTING);
            glEnable(GL_LIGHT0);
            GLfloat light_pos[] = {cam_x, cam_y, cam_z, 1.0f};
            GLfloat diffuse[] = {0.8f, 0.8f, 0.8f, 1.0f};
            GLfloat specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
            glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
            glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
            glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

            // Render VBO
            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_NORMAL_ARRAY);
            glEnableClientState(GL_COLOR_ARRAY);
            glVertexPointer(3, GL_FLOAT, 9 * sizeof(float), (void*)0);
            glNormalPointer(GL_FLOAT, 9 * sizeof(float), (void*)(3 * sizeof(float)));
            glColorPointer(3, GL_FLOAT, 9 * sizeof(float), (void*)(6 * sizeof(float)));
            glDrawArrays(GL_TRIANGLES, 0, vertex_data.size() / 9);
            glDisableClientState(GL_VERTEX_ARRAY);
            glDisableClientState(GL_NORMAL_ARRAY);
            glDisableClientState(GL_COLOR_ARRAY);
            glBindBuffer(GL_ARRAY_BUFFER, 0);

            glDisable(GL_LIGHTING);
            glDisable(GL_COLOR_MATERIAL);
		}
		
		if (current_display_mode == SURFACE) {
            glEnable(GL_LIGHTING);
            glEnable(GL_LIGHT0);
            GLfloat light_pos[] = {cam_x, cam_y, cam_z, 1.0f};
            glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

            // Render 2D surface (x-y plane, z = scalar value)
            int k = sim.size / 2;
            glBegin(GL_QUADS);
            for (int i = 0; i < sim.size - 1; ++i) {
                for (int j = 0; j < sim.size - 1; ++j) {
                    float x0 = static_cast<float>((i - sim.size / 2) * sim.dx);
                    float y0 = static_cast<float>((j - sim.size / 2) * sim.dx);
                    float x1 = static_cast<float>((i + 1 - sim.size / 2) * sim.dx);
                    float y1 = static_cast<float>((j + 1 - sim.size / 2) * sim.dx);

                    float values[4] = {
                        static_cast<float>((sim.computed_scalar[i][j][k] - min_scalar) / range),
                        static_cast<float>((sim.computed_scalar[i + 1][j][k] - min_scalar) / range),
                        static_cast<float>((sim.computed_scalar[i + 1][j + 1][k] - min_scalar) / range),
                        static_cast<float>((sim.computed_scalar[i][j + 1][k] - min_scalar) / range)
                    };
                    float z[4] = {
                        static_cast<float>(values[0] * sim.dx * vis_scale),
                        static_cast<float>(values[1] * sim.dx * vis_scale),
                        static_cast<float>(values[2] * sim.dx * vis_scale),
                        static_cast<float>(values[3] * sim.dx * vis_scale)
                    };

                    // Compute normal
                    float vec1[3] = {x1 - x0, y0 - y0, z[1] - z[0]};
                    float vec2[3] = {x0 - x0, y1 - y0, z[3] - z[0]};
                    float normal[3] = {
                        vec1[1] * vec2[2] - vec1[2] * vec2[1],
                        vec1[2] * vec2[0] - vec1[0] * vec2[2],
                        vec1[0] * vec2[1] - vec1[1] * vec2[0]
                    };
                    float norm_len = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
                    if (norm_len > 0) {
                        normal[0] /= norm_len;
                        normal[1] /= norm_len;
                        normal[2] /= norm_len;
                    }

                    glNormal3fv(normal);

                    // Vertex 0
                    glColor3f(values[0] * vis_color_intensity, 0.5f, (1.0f - values[0]) * vis_color_intensity);
                    glVertex3f(x0 * vis_scale, y0 * vis_scale, z[0]);

                    // Vertex 1
                    glColor3f(values[1] * vis_color_intensity, 0.5f, (1.0f - values[1]) * vis_color_intensity);
                    glVertex3f(x1 * vis_scale, y0 * vis_scale, z[1]);

                    // Vertex 2
                    glColor3f(values[2] * vis_color_intensity, 0.5f, (1.0f - values[2]) * vis_color_intensity);
                    glVertex3f(x1 * vis_scale, y1 * vis_scale, z[2]);

                    // Vertex 3
                    glColor3f(values[3] * vis_color_intensity, 0.5f, (1.0f - values[3]) * vis_color_intensity);
                    glVertex3f(x0 * vis_scale, y1 * vis_scale, z[3]);
                }
            }
            glEnd();

            glDisable(GL_LIGHTING);
            glDisable(GL_COLOR_MATERIAL);
        }

        if (current_display_mode == WIREFRAME) {
            glBegin(GL_LINES);
            glColor3f(0.5f, 0.5f, 1.0f);
            for (int i = 0; i < sim.size; ++i) {
                for (int j = 0; j < sim.size; ++j) {
                    for (int k = 0; k < sim.size; ++k) {
                        float x = static_cast<float>((i - sim.size / 2) * sim.dx);
                        float y = static_cast<float>((j - sim.size / 2) * sim.dx);
                        float z = static_cast<float>((k - sim.size / 2) * sim.dx);
                        if (i < sim.size - 1) {
                            glVertex3f(x * vis_scale, y * vis_scale, z * vis_scale);
                            glVertex3f((x + sim.dx) * vis_scale, y * vis_scale, z * vis_scale);
                        }
                        if (j < sim.size - 1) {
                            glVertex3f(x * vis_scale, y * vis_scale, z * vis_scale);
                            glVertex3f(x * vis_scale, (y + sim.dx) * vis_scale, z * vis_scale);
                        }
                        if (k < sim.size - 1) {
                            glVertex3f(x * vis_scale, y * vis_scale, z * vis_scale);
                            glVertex3f(x * vis_scale, y * vis_scale, (z + sim.dx) * vis_scale);
                        }
                    }
                }
            }
            glEnd();
        }
    }
	
	// Render XYZ arrows
    glBegin(GL_LINES);
    // X-axis (red)
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(-1.0f * vis_scale, -1.0f * vis_scale, -1.0f * vis_scale);
    glVertex3f(1.0f * vis_scale, -1.0f * vis_scale, -1.0f * vis_scale);
    // Y-axis (green)
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(-1.0f * vis_scale, -1.0f * vis_scale, -1.0f * vis_scale);
    glVertex3f(-1.0f * vis_scale, 1.0f * vis_scale, -1.0f * vis_scale);
    // Z-axis (blue)
    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(-1.0f * vis_scale, -1.0f * vis_scale, -1.0f * vis_scale);
    glVertex3f(-1.0f * vis_scale, -1.0f * vis_scale, 1.0f * vis_scale);
    glEnd();

    // Render text (equation name and simulation time)
    if (font && !sim.equations.empty()) {
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrtho(0, 800, 0, 600, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        glEnable(GL_TEXTURE_2D);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        std::string text = sim.equations[sim.current_equation]->name() + " | t = " + std::to_string(sim.t);
        SDL_Color color = {255, 255, 255, 255};
        int w, h;
        GLuint texture = render_text_to_texture(text, color, w, h);
        if (texture) {
            glBindTexture(GL_TEXTURE_2D, texture);
            glBegin(GL_QUADS);
            glTexCoord2f(0, 1); glVertex2f(10, 590 - h);
            glTexCoord2f(1, 1); glVertex2f(10 + w, 590 - h);
            glTexCoord2f(1, 0); glVertex2f(10 + w, 590);
            glTexCoord2f(0, 0); glVertex2f(10, 590);
            glEnd();
            glDeleteTextures(1, &texture);
        }

        glDisable(GL_TEXTURE_2D);
        glDisable(GL_BLEND);
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
    }
}