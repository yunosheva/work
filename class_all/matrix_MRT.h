#pragma once
#include <iostream> 
#include <vector> 
#include <fstream>
#include <cmath>
#include <sstream>
#include <omp.h>

using namespace std;

const int kMax = 19;

const double M[kMax][kMax] = {
   {  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  1.0,  1.0, 1.0, 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0}, // rho
   {-30.0, -11.0, -11.0, -11.0, -11.0, -11.0, -11.0,  8.0,  8.0, 8.0, 8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0}, // e
   { 12.0,  -4.0,  -4.0,  -4.0,  -4.0,  -4.0,  -4.0,  1.0,  1.0, 1.0, 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0}, // eps
   {  0.0,   1.0,  -1.0,   0.0,   0.0,   0.0,   0.0,  1.0, -1.0, 1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0}, // jx
   {  0.0,  -4.0,   4.0,   0.0,   0.0,   0.0,   0.0,  1.0, -1.0, 1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0}, // qx
   {  0.0,   0.0,   0.0,   1.0,  -1.0,   0.0,   0.0,  1.0,  1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0}, // jy
   {  0.0,   0.0,   0.0,  -4.0,   4.0,   0.0,   0.0,  1.0,  1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0}, // qy
   {  0.0,   0.0,   0.0,   0.0,   0.0,   1.0,  -1.0,  0.0,  0.0, 0.0, 0.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0, -1.0}, // jz
   {  0.0,   0.0,   0.0,   0.0,   0.0,  -4.0,   4.0,  0.0,  0.0, 0.0, 0.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0, -1.0}, // qz
   {  0.0,   2.0,   2.0,  -1.0,  -1.0,  -1.0,  -1.0,  1.0,  1.0, 1.0, 1.0,  1.0,  1.0,  1.0,  1.0, -2.0, -2.0, -2.0, -2.0}, // pxx
   {  0.0,  -4.0,  -4.0,   2.0,   2.0,   2.0,   2.0,  1.0,  1.0, 1.0, 1.0,  1.0,  1.0,  1.0,  1.0, -2.0, -2.0, -2.0, -2.0}, // pi_xx
   {  0.0,   0.0,   0.0,   1.0,   1.0,  -1.0,  -1.0,  1.0,  1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0}, // p_ww
   {  0.0,   0.0,   0.0,  -2.0,  -2.0,   2.0,   2.0,  1.0,  1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0}, // pi_ww
   {  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  1.0, -1.0, -1.0, 1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // p_xy
   {  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0}, // p_yz
   {  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0, 0.0,  1.0, -1.0, -1.0,  1.0,  0.0,  0.0,  0.0,  0.0}, // p_xz
   {  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  1.0, -1.0, 1.0, -1.0, -1.0,  1.0, -1.0,  1.0,  0.0,  0.0,  0.0,  0.0}, // mx
   {  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, -1.0, -1.0, 1.0, 1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0}, // my
   {  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0, 0.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0,  1.0,  1.0}, // mz
};

const double Mi[kMax][kMax] = {
   {1.0 / 19.0,  -5.0 / 399.0,   1.0 / 21.0,   0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
   {1.0 / 19.0, -11.0 / 2394.0, -1.0 / 63.0,   1.0 / 10.0, -1.0 / 10.0, 0.0, 0.0,  0.0, 0.0, 1.0 / 18.0, -1.0 / 18.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0},
   {1.0 / 19.0, -11.0 / 2394.0, -1.0 / 63.0,  -1.0 / 10.0,  1.0 / 10.0, 0.0, 0.0, 0.0, 0.0, 1.0 / 18.0, -1.0 / 18.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0},
   {1.0 / 19.0, -11.0 / 2394.0, -1.0 / 63.0,   0.0, 0.0, 1.0 / 10.0, -1.0 / 10.0,  0.0, 0.0, -1.0 / 36.0,  1.0 / 36.0, 1.0 / 12.0, -1.0 / 12.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0},
   {1.0 / 19.0, -11.0 / 2394.0, -1.0 / 63.0,   0.0, 0.0, -1.0 / 10.0, 1.0 / 10.0, 0.0, 0.0, -1.0 / 36.0,  1.0 / 36.0, 1.0 / 12.0, -1.0 / 12.0, 0.0, 0.0, 0.0,  0.0,  0.0, 0.0},
   {1.0 / 19.0, -11.0 / 2394.0, -1.0 / 63.0,   0.0, 0.0, 0.0, 0.0,  1.0 / 10.0, -1.0 / 10.0, -1.0 / 36.0,  1.0 / 36.0, -1.0 / 12.0, 1.0 / 12.0, 0.0,  0.0, 0.0,  0.0, 0.0, 0.0},
   {1.0 / 19.0, -11.0 / 2394.0, -1.0 / 63.0,   0.0, 0.0, 0.0, 0.0, -1.0 / 10.0,  1.0 / 10.0, -1.0 / 36.0,  1.0 / 36.0, -1.0 / 12.0, 1.0 / 12.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0},
   {1.0 / 19.0,   4.0 / 1197.0,  1.0 / 252.0,  1.0 / 10.0,  1.0 / 40.0, 1.0 / 10.0, 1.0 / 40.0,  0.0, 0.0,  1.0 / 36.0,  1.0 / 72.0, 1.0 / 12.0,  1.0 / 24.0, 1.0 / 4.0, 0.0,  0.0, 1.0 / 8.0, -1.0 / 8.0, 0.0},
   {1.0 / 19.0,   4.0 / 1197.0,  1.0 / 252.0, -1.0 / 10.0, -1.0 / 40.0, 1.0 / 10.0, 1.0 / 40.0,  0.0, 0.0,  1.0 / 36.0,  1.0 / 72.0, 1.0 / 12.0,  1.0 / 24.0, -1.0 / 4.0, 0.0,  0.0, -1.0 / 8.0, -1.0 / 8.0, 0.0},
   {1.0 / 19.0,   4.0 / 1197.0,  1.0 / 252.0,  1.0 / 10.0,  1.0 / 40.0, -1.0 / 10.0, -1.0 / 40.0,  0.0, 0.0,  1.0 / 36.0,  1.0 / 72.0, 1.0 / 12.0,  1.0 / 24.0, -1.0 / 4.0, 0.0,  0.0, 1.0 / 8.0, 1.0 / 8.0, 0.0},
   {1.0 / 19.0,   4.0 / 1197.0,  1.0 / 252.0, -1.0 / 10.0, -1.0 / 40.0, -1.0 / 10.0, -1.0 / 40.0,  0.0, 0.0,  1.0 / 36.0,  1.0 / 72.0, 1.0 / 12.0,  1.0 / 24.0, 1.0 / 4.0, 0.0,  0.0,-1.0 / 8.0, 1.0 / 8.0, 0.0},
   {1.0 / 19.0,   4.0 / 1197.0,  1.0 / 252.0,  1.0 / 10.0,  1.0 / 40.0,  0.0, 0.0,  1.0 / 10.0,  1.0 / 40.0,  1.0 / 36.0,  1.0 / 72.0, -1.0 / 12.0, -1.0 / 24.0, 0.0, 0.0,  1.0 / 4.0, -1.0 / 8.0, 0.0, 1.0 / 8.0},
   {1.0 / 19.0,   4.0 / 1197.0,  1.0 / 252.0, -1.0 / 10.0, -1.0 / 40.0,  0.0, 0.0,  1.0 / 10.0,  1.0 / 40.0,  1.0 / 36.0,  1.0 / 72.0, -1.0 / 12.0, -1.0 / 24.0, 0.0, 0.0, -1.0 / 4.0,  1.0 / 8.0, 0.0, 1.0 / 8.0},
   {1.0 / 19.0,   4.0 / 1197.0,  1.0 / 252.0,  1.0 / 10.0,  1.0 / 40.0,  0.0, 0.0, -1.0 / 10.0, -1.0 / 40.0,  1.0 / 36.0,  1.0 / 72.0, -1.0 / 12.0, -1.0 / 24.0, 0.0, 0.0, -1.0 / 4.0, -1.0 / 8.0, 0.0, -1.0 / 8.0},
   {1.0 / 19.0,   4.0 / 1197.0,  1.0 / 252.0, -1.0 / 10.0, -1.0 / 40.0,  0.0, 0.0, -1.0 / 10.0, -1.0 / 40.0,  1.0 / 36.0,  1.0 / 72.0, -1.0 / 12.0, -1.0 / 24.0, 0.0, 0.0,  1.0 / 4.0,  1.0 / 8.0, 0.0, -1.0 / 8.0},
   {1.0 / 19.0,   4.0 / 1197.0,  1.0 / 252.0,  0.0, 0.0, 1.0 / 10.0, 1.0 / 40.0,  1.0 / 10.0,  1.0 / 40.0, -1.0 / 18.0, -1.0 / 36.0, 0.0, 0.0, 0.0,  1.0 / 4.0, 0.0,  0.0, 1.0 / 8.0, -1.0 / 8.0},
   {1.0 / 19.0,   4.0 / 1197.0,  1.0 / 252.0,  0.0, 0.0, -1.0 / 10.0, -1.0 / 40.0,  1.0 / 10.0,  1.0 / 40.0, -1.0 / 18.0, -1.0 / 36.0, 0.0, 0.0, 0.0, -1.0 / 4.0, 0.0,  0.0, -1.0 / 8.0, -1.0 / 8.0},
   {1.0 / 19.0,   4.0 / 1197.0,  1.0 / 252.0,  0.0, 0.0, 1.0 / 10.0, 1.0 / 40.0, -1.0 / 10.0, -1.0 / 40.0, -1.0 / 18.0, -1.0 / 36.0, 0.0, 0.0, 0.0, -1.0 / 4.0, 0.0,  0.0, 1.0 / 8.0, 1.0 / 8.0},
   {1.0 / 19.0,   4.0 / 1197.0,  1.0 / 252.0,  0.0, 0.0, -1.0 / 10.0, -1.0 / 40.0, -1.0 / 10.0, -1.0 / 40.0, -1.0 / 18.0, -1.0 / 36.0, 0.0, 0.0,0.0,  1.0 / 4.0, 0.0,  0.0, -1.0 / 8.0, 1.0 / 8.0} };