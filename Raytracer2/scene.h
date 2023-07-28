#include<iostream>
#include "Vector3D.h"
#include "Matrix3x3.h"
#include <cmath>

//Render options
const int RESOLUTION[2] = {1024,786};
const int iterations = 250;
const double stepsize = 0.16;

//Geometry
//Position of the camera in units of r_S
//x,z are the disc's plane
Vector3D LOOKAT(0.0, 0.0, 0.0);
Vector3D CAMERA_POS(0.0, 1.0, -20.0);
Vector3D UPVEC(0.2, 1.0, 0.0);
const double FOV = 2;

//Disk radius
const double disk_inner = 2.6, disk_outer = 14.0;
const double disk_inner_sq = disk_inner*disk_inner;
const double disk_outer_sq = disk_outer*disk_outer;

//Material
const bool sRGB_in = 1, sRGB_out = 1;

const int sky_width = 4096;
const int sky_height = 2048;
double bg_image[sky_width][sky_height][3];

const Vector3D white(1.0, 1.0, 1.0);
const Vector3D black(0.0, 0.0, 0.0);

const int disk_width = 2512;
const int disk_height = 400;
double ad_image[disk_width][disk_height][3];
