#pragma once

#include <cuda_runtime.h>
#include<iostream>
#include "Matrix3x3.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//Render options
__device__ const int RESOLUTION[2] = {1024,786};
__device__ const int iterations = 250;
__device__ const double stepsize = 0.16;

//Geometry
//Position of the camera in units of r_S
//x,z are the disc's plane
// Vector3D LOOKAT(0.0, 0.0, 0.0);
// Vector3D CAMERA_POS(0.0, 1.0, -20.0);
// Vector3D UPVEC(0.2, 1.0, 0.0);
// __device__ const double FOV = 2;

//Disk radius
__device__ const double disk_inner = 2.6, disk_outer = 14.0;
__device__ const double disk_inner_sq = 6.76;
__device__ const double disk_outer_sq = 196.0;

//Material
__device__ const bool sRGB_in = 1, sRGB_out = 0;

__device__ const int sky_width = 4096;
__device__ const int sky_height = 2048;
__device__ double bg_image[sky_width][sky_height][3];

const Vector3D white(1.0, 1.0, 1.0);
const Vector3D black(0.0, 0.0, 0.0);

__device__ const int disk_width = 2512;
__device__ const int disk_height = 400;
__device__ double ad_image[disk_width][disk_height][3];


