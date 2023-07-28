#include <iostream>
#include <fstream>
#include "Vector3D.h"
#include "Matrix3x3.h"
#include "utils.h"
// #include "texture.h"

int main(){
    Vector3D FRONTVEC(LOOKAT - CAMERA_POS);
    FRONTVEC.normalize();
    Vector3D LEFTVEC(cross(UPVEC, FRONTVEC));
    LEFTVEC.normalize();

    Vector3D NUPVEC(cross(FRONTVEC, LEFTVEC));

    double data[9] = {LEFTVEC[0], NUPVEC[0], FRONTVEC[0], 
                        LEFTVEC[1], NUPVEC[1], FRONTVEC[1],
                        LEFTVEC[2], NUPVEC[2], FRONTVEC[2]};

    Matrix3x3 viewmatrix(data);

    viewmatrix.T();

    // int numPixels = RESOLUTION[0]*RESOLUTION[1];
    // int pixelindices[numPixels];
    // for(int i = 0; i < numPixels; i++){
    //     pixelindices[i] = i;
    // }

    std::ofstream MyFile("3.ppm");

    MyFile << "P3\n" << RESOLUTION[0] << ' ' << RESOLUTION[1] << "\n255\n";

    Vector3D object_color(0.0, 0.0, 0.0);
    double object_alpha = 0.0;
    bgimage_initialize();
    adimage_initialize();

    for(int y = 0; y < RESOLUTION[1] ; y++){
        for(int x = 0; x < RESOLUTION[0]; x++){
            Vector3D view = Vector3D(((float)x / RESOLUTION[0] - 0.5f) * FOV,
            ((-(float)y / RESOLUTION[1] + 0.5f) * RESOLUTION[1] / RESOLUTION[0]) * FOV, 1.0);

            view = transform(view, viewmatrix);
            Vector3D normView = view / view.norm();

            Vector3D velocity = normView;
            Vector3D point = CAMERA_POS;
            double pointsq = point.norm2();

            double h2 = SetInitialConditions(point, velocity);

            for(int iter=0; iter < iterations; ++iter){
                Vector3D prevpoint = point;
                double prevsqnorm = pointsq;

                RK4(point, velocity, h2, stepsize);

                pointsq = point.norm2();

                //accretion disk
                int diskmask = 0;

                Vector3D temp(0.0, 0.0, 0.0);
                object_color = temp;
                object_alpha = 0.0;

                if(point[1]*prevpoint[1] < 0){
                    if(pointsq > disk_inner_sq && pointsq < disk_outer_sq){
                        diskmask = 1;
                    }
                }

                if(diskmask == 1){
                    double lambda = - point[1]/velocity[1];
                    Vector3D col_point = point + lambda*velocity;
                    double colpointnorm = col_point.norm();

                    double phi = atan(col_point[0]/point[2]);
                    double uv[2] = {mod(phi+2.0*M_PI, 2.0*M_PI)/(2.0*M_PI), (colpointnorm - disk_inner)/(disk_outer - disk_inner)};
                    clip(uv[0], 0.0, 1.0);
                    clip(uv[1], 0.0, 1.0);

                    // Vector3D disk_color(1.0, 1.0, 0.98);
                    // double disk_alpha = 1.0;

                    Vector3D disk_color;
                    double disk_alpha;
                    disk_color = texture_lookup_ad(uv);
                    disk_alpha = diskmask * clip(disk_color.norm2()/3.0, 0.0, 1.0);

                    object_color = blendcolors(disk_color, disk_alpha, object_color, object_alpha);
                    object_alpha = blendalpha(disk_alpha, object_alpha);

                    //std::cout << "Hit disk \n";
                    break;
                }

                //event horizon
                int maskhorizon = 0;
                if(pointsq < 1 && prevsqnorm > 1){
                    maskhorizon = 1;
                }
                if(maskhorizon == 1){
                    Vector3D horizon_color(0.0, 0.0, 0.0);
                    double horizon_alpha = 1.0;

                    object_color = blendcolors(horizon_color, horizon_alpha, object_color, object_alpha);
                    object_alpha = blendalpha(horizon_alpha, object_alpha);
                    //std::cout << "Hit horizon \n";
                    break;
                }
            }

            double vphi = atan(velocity[0]/velocity[2]);
            double vtheta = atan(velocity[1]/velocity.norm());
            // double vuv[2] = {mod((vphi+4.5),(2.0*M_PI))/(2.0*M_PI), (vtheta+M_PI/2.0)/M_PI};

            double vuv[2] = {mod(0.5 + vphi/(2.0*M_PI), 1.0), mod(0.5 - (vtheta/M_PI), 1.0)};
            if(vuv[0] < 0){vuv[0] = 1.0 + vuv[0];}
            if(vuv[1] < 0){vuv[1] = 1.0 + vuv[1];}

            Vector3D bg_color = texture_lookup_bg(vuv);
            // Vector3D bg_color(bg_image[x][y][0], bg_image[x][y][1], bg_image[x][y][2]);
            double bg_alpha = 1.0;

            object_color = blendcolors(bg_color, bg_alpha, object_color, object_alpha);
            object_alpha = blendalpha(bg_alpha, object_alpha);

            double color[3] = {object_color[0], object_color[1], object_color[2]};

            if(sRGB_out){
                rgbtosrgb(color);
            }

            int ir = static_cast<int>(255.999 * color[0]);
            int ig = static_cast<int>(255.999 * color[1]);
            int ib = static_cast<int>(255.999 * color[2]);

            MyFile << ir << ' ' << ig << ' ' << ib << '\n';   

            //std::cout << x << " " << y << "\n"; 
        }
    }
    MyFile.close();
}


