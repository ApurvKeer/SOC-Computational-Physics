// #include "Vector3D.h"
// #include "Matrix3x3.h"
#include <iostream>
#include <fstream>
#include "utils.h"

using namespace std;

__host__ __device__ void rgbtosrgb(double* arr) {
    const double threshold = 0.0031308f;
    for (int i = 0; i < 3; ++i) {
        if (arr[i] > threshold) {
            arr[i] = std::pow(arr[i], 1.0f / 2.4f) * 1.055f - 0.055f;
        } else {
            arr[i] *= 12.92f;
        }
    }
}

__host__ __device__ void srgbtorgb(double* arr) {
    const double threshold = 0.04045f;
    for (int i = 0; i < 3; ++i) {
        if (arr[i] > threshold) {
            arr[i] = std::pow((arr[i] + 0.055f) / 1.055f, 2.4f);
        } else {
            arr[i] /= 12.92f;
        }
    }
}

__host__ __device__ double* bgimage_initialize(bool sRGB_in) {
    double* bg_image = new double[sky_width * sky_height * 3];
    
    for (int x = 0; x < sky_width; ++x) {
        for (int y = 0; y < sky_height; ++y) {
            double rgb[3];
            int pixel = 3 * (x + y * sky_width);
            rgb[0] = 0.0;
            rgb[1] = 0.0;
            rgb[2] = 0.0;
            if (sRGB_in) {
                srgbtorgb(rgb);
            }
            bg_image[pixel] = rgb[0];
            bg_image[pixel + 1] = rgb[1];
            bg_image[pixel + 2] = rgb[2];
        }
    }
    return bg_image;
}

__host__ __device__ double* adimage_initialize(bool sRGB_in, unsigned char *img){
    double* ad_image = new double[disk_width * disk_height * 3];
    
    for (int x = 0; x < disk_width; ++x) {
        for (int y = 0; y < disk_height; ++y) {
            double rgb[3];
            int pixel = 3 * (x * disk_height + y);
            rgb[0] = 1.0;
            rgb[1] = 1.0;
            rgb[2] = 1.0;
            if (sRGB_in) {
                srgbtorgb(rgb);
            }
            ad_image[pixel] = rgb[0];
            ad_image[pixel + 1] = rgb[1];
            ad_image[pixel + 2] = rgb[2];
        }
    }
    return ad_image;
}




// int main(){
//     std::ofstream MyFile("2.ppm");

//     MyFile << "P3\n" << RESOLUTION[0] << ' ' << RESOLUTION[1] << "\n255\n";

//     bgimage_initialize();

//     for(int x = 0; x < RESOLUTION[0]; ++x){
//         for(int y = 0; y < RESOLUTION[1]; ++y){
//             int ir = static_cast<int>(255.999 * bg_image[x][y][0]);
//             int ig = static_cast<int>(255.999 * bg_image[x][y][1]);
//             int ib = static_cast<int>(255.999 * bg_image[x][y][2]);

//             MyFile << ir << ' ' << ig << ' ' << ib << '\n'; 
//         }
//     }
//     MyFile.close();
// }


