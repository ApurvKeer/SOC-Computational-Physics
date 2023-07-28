// #include "Vector3D.h"
// #include "Matrix3x3.h"
#include <iostream>
#include <fstream>
#include "scene.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image/stb_image_write.h"

using namespace std;

void rgbtosrgb(double* arr) {
    const double threshold = 0.0031308;
    for (int i = 0; i < 3; ++i) {
        if (arr[i] > threshold) {
            arr[i] = std::pow(arr[i], 1.0 / 2.4) * 1.055 - 0.055;
        } else {
            arr[i] *= 12.92;
        }
    }
}

void srgbtorgb(double* arr) {
    const double threshold = 0.04045;
    for (int i = 0; i < 3; ++i) {
        if (arr[i] > threshold) {
            arr[i] = std::pow((arr[i] + 0.055) / 1.055, 2.4);
        } else {
            arr[i] /= 12.92;
        }
    }
}

void bgimage_initialize() {
    int w, h, c;
    unsigned char *img = stbi_load("bgedit.jpg", &w, &h, &c, 0);
    for (int x = 0; x < sky_width; ++x) {
        for (int y = 0; y < sky_height; ++y) {
            double rgb[3];
            int pixel = 3*(x + y*sky_width);
            rgb[0] = double(img[pixel])/255.0;
            rgb[1] = double(img[pixel + 1])/255.0;
            rgb[2] = double(img[pixel + 2])/255.0;
            if (sRGB_in) {
                srgbtorgb(rgb);
            }
            bg_image[x][y][0] = rgb[0];
            bg_image[x][y][1] = rgb[1];
            bg_image[x][y][2] = rgb[2];
        }
    }
}

void adimage_initialize(){
    int w, h, c;
    unsigned char *img = stbi_load("adisk.jpg", &w, &h, &c, 0);
    for (int x = 0; x < w; ++x) {
        for (int y = 0; y < h; ++y) {
            double rgb[3];
            int pixel = 3*(x*h + y);
            rgb[0] = double(img[pixel])/255.0;
            rgb[1] = double(img[pixel + 1])/255.0;
            rgb[2] = double(img[pixel + 2])/255.0;
            if (sRGB_in) {
                srgbtorgb(rgb);
            }
            ad_image[x][y][0] = rgb[0];
            ad_image[x][y][1] = rgb[1];
            ad_image[x][y][2] = rgb[2];
        }
    }
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


