#include <iostream>
#include <cmath>
#include "texture.h"
#include "Vector3D.h"
#include "Matrix3x3.h"
#include <random>

Vector3D transform(Vector3D v, Matrix3x3 A){
    double w0 = v[0]*A(0,0) + v[1]*A(1,0) + v[2]*A(2,0);
    double w1 = v[0]*A(0,1) + v[1]*A(1,1) + v[2]*A(2,1);
    double w2 = v[0]*A(0,2) + v[1]*A(1,2) + v[2]*A(2,2);
    Vector3D w(w0, w1, w2);
    return w;
}

double SetInitialConditions(Vector3D& point, Vector3D& velocity) {
    Vector3D c = cross(point, velocity);
    return c.norm2();
}

Vector3D dxdt(Vector3D point, Vector3D velocity){
    return velocity;
}

Vector3D dvdt(Vector3D point, Vector3D velocity, double h2){
    Vector3D accl = -1.5 * h2 * point / pow(point.norm2(), 2.5);
    // Vector3D accl(0.0, 0.0, 0.0);
    return accl;
}

void RK4(Vector3D& point, Vector3D& velocity, double h2, double step){
    Vector3D k1 = dxdt(point, velocity);
    Vector3D l1 = dvdt(point, velocity, h2);
    Vector3D k2 = dxdt(point + 0.5*step*k1, velocity + 0.5*step*l1);
    Vector3D l2 = dvdt(point + 0.5*step*k1, velocity + 0.5*step*l1, h2);
    Vector3D k3 = dxdt(point + 0.5*step*k2, velocity + 0.5*step*l2);
    Vector3D l3 = dvdt(point + 0.5*step*k2, velocity + 0.5*step*l2, h2);
    Vector3D k4 = dxdt(point + step*k3, velocity + step*l3);
    Vector3D l4 = dvdt(point + step*k3, velocity + step*l3, h2);

    point = point + (1.0/6.0)*step*(k1 + 2*k2 + 2*k3 + k4);
    velocity = velocity + (1.0/6.0)*step*(l1 + 2*l2 + 2*l3 + l4);
}

Vector3D blendcolors(Vector3D cb, double balpha, Vector3D ca, double aalpha){
    Vector3D c = ca + cb * (balpha*(1.0 - aalpha));
    return c;
}

double blendalpha(double balpha, double aalpha){
    return aalpha + balpha*(1.0 - aalpha);
}

double mod(double a, double b){
    if(a < b){
        return a;
    }
    else{
        int n = static_cast<int>(a/b);
        return a - n*b;
    }
}

// Vector3D random_color(int x, int y){
//     Vector3D r(double(x)/RESOLUTION[0], double(y)/RESOLUTION[1], 0.25);
//     return r;
// }

// void bg_image_initial(){
//     for(int x=0; x < RESOLUTION[0]; ++x){
//         for(int y=0; y < RESOLUTION[1]; ++y){
//             Vector3D r;
//             if(rand()%100 > 95){
//                 r = white;
//             }
//             else{
//                 r = black;
//             }
//             bg_image[x][y][0] = r[0];
//             bg_image[x][y][1] = r[1];
//             bg_image[x][y][2] = r[2];
//         }
//     }
// }

double clip(double x, double min, double max){
    if(x < min){
        return min;
    }
    else if (x > max){
        return max;
    }
    else{
        return x;
    }
    
}

Vector3D texture_lookup_bg(double* uvarrin){
    double uvarr[2] = {clip(uvarrin[0], 0.0, 0.999), clip(uvarrin[1], 0.0, 0.999)};

    uvarr[0] *= double(sky_height);
    uvarr[1] *= double(sky_width);

    int x = int(uvarr[1]);
    int y = int(uvarr[0]);

    Vector3D image(bg_image[x][y][0], bg_image[x][y][1], bg_image[x][y][2]);
    return image;
}

Vector3D texture_lookup_ad(double* uvarrin){
    double uvarr[2] = {clip(uvarrin[0], 0.0, 0.999), clip(uvarrin[1], 0.0, 0.999)};

    uvarr[0] *= double(disk_height);
    uvarr[1] *= double(disk_width);

    int x = int(uvarr[1]);
    int y = int(uvarr[0]);

    Vector3D image(ad_image[x][y][0], ad_image[x][y][1], ad_image[x][y][2]);
    return image;
}
