#include <iostream>
#include <time.h>
#include "texture.h"

// limited version of checkCudaErrors from helper_cuda.h in CUDA examples
#define checkCudaErrors(val) check_cuda( (val), #val, __FILE__, __LINE__ )

void check_cuda(cudaError_t result, char const *const func, const char *const file, int const line) {
    if (result) {
        std::cerr << "CUDA error = " << static_cast<unsigned int>(result) << " at " <<
            file << ":" << line << " '" << func << "' \n";
        // Make sure we call CUDA Device Reset before exiting
        cudaDeviceReset();
        exit(99);
    }
}

__global__ void render(Vector3D *fb, int max_x, int max_y) {
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;
    if((x >= max_x) || (y >= max_y)) return;
    int pixel_index = y*max_x + x;

    Vector3D LOOKAT(0.0f, 0.0f, 0.0f);
    Vector3D CAMERA_POS(0.0f, 1.0f, -20.0f);
    Vector3D UPVEC(0.2f, 1.0f, 0.0f);
    double FOV = 2;

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

    Vector3D object_color(0.0f, 0.0f, 0.0f);
    double object_alpha = 0.0f;
    
    Vector3D view = Vector3D(((float)x / RESOLUTION[0] - 0.5f) * FOV,
    ((-(float)y / RESOLUTION[1] + 0.5f) * RESOLUTION[1] / RESOLUTION[0]) * FOV, 1.0f);

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

        Vector3D temp(0.0f, 0.0f, 0.0f);
        object_color = temp;
        object_alpha = 0.0f;

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
            double uv[2] = {mod(phi+2.0f*M_PI, 2.0f*M_PI)/(2.0f*M_PI), (colpointnorm - disk_inner)/(disk_outer - disk_inner)};
            clip(uv[0], 0.0f, 1.0f);
            clip(uv[1], 0.0f, 1.0f);

            // Vector3D disk_color(1.0, 1.0, 0.98);
            // double disk_alpha = 1.0;

            // Vector3D disk_color;
            double disk_alpha;
            // disk_color = texture_lookup_ad(uv);
            Vector3D disk_color(1.0f, 1.0f, 1.0f);
            disk_alpha = diskmask * clip(disk_color.norm2()/3.0f, 0.0f, 1.0f);

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
            Vector3D horizon_color(0.0f, 0.0f, 0.0f);
            double horizon_alpha = 1.0f;

            object_color = blendcolors(horizon_color, horizon_alpha, object_color, object_alpha);
            object_alpha = blendalpha(horizon_alpha, object_alpha);
            //std::cout << "Hit horizon \n";
            break;
        }
    }

    double vphi = atan(velocity[0]/velocity[2]);
    double vtheta = atan(velocity[1]/velocity.norm());
    // double vuv[2] = {mod((vphi+4.5),(2.0*M_PI))/(2.0*M_PI), (vtheta+M_PI/2.0)/M_PI};

    double vuv[2] = {mod(0.5f + vphi/(2.0f*M_PI), 1.0f), mod(0.5f - (vtheta/M_PI), 1.0f)};
    if(vuv[0] < 0){vuv[0] = 1.0f + vuv[0];}
    if(vuv[1] < 0){vuv[1] = 1.0f + vuv[1];}

    // Vector3D bg_color = texture_lookup_bg(vuv, bg_image);
    Vector3D bg_color(0.0f, 0.0f, 0.0f);
    double bg_alpha = 1.0f;

    object_color = blendcolors(bg_color, bg_alpha, object_color, object_alpha);
    object_alpha = blendalpha(bg_alpha, object_alpha);

    double color[3] = {object_color[0], object_color[1], object_color[2]};

    if(sRGB_out){
        rgbtosrgb(color);
    }

    int ir = static_cast<int>(255.999f * color[0]);
    int ig = static_cast<int>(255.999f * color[1]);
    int ib = static_cast<int>(255.999f * color[2]);

    fb[pixel_index] = Vector3D(ir, ig, ib);
}

int main() {
    int nx = 1024;
    int ny = 786;
    int tx = 8;
    int ty = 8;

    std::cerr << "Rendering a " << nx << "x" << ny << " image ";
    std::cerr << "in " << tx << "x" << ty << " blocks.\n";

    int num_pixels = nx*ny;
    size_t fb_size = num_pixels*sizeof(Vector3D);

    // allocate FB
    Vector3D *fb;
    checkCudaErrors(cudaMallocManaged((void **)&fb, fb_size));

    clock_t start, stop;
    start = clock();
    // Render our buffer
    dim3 blocks(nx/tx+1,ny/ty+1);
    dim3 threads(tx,ty);
    render<<<blocks, threads>>>(fb, nx, ny);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());
    stop = clock();
    double timer_seconds = ((double)(stop - start)) / CLOCKS_PER_SEC;
    std::cerr << "took " << timer_seconds << " seconds.\n";

    // Output FB as Image
    std::cout << "P3\n" << nx << " " << ny << "\n255\n";
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            size_t pixel_index = j*nx + i;
            int ir = int(255.99f*fb[pixel_index][0]);
            int ig = int(255.99f*fb[pixel_index][1]);
            int ib = int(255.99f*fb[pixel_index][2]);
            std::cout << ir << " " << ig << " " << ib << "\n";
        }
    }

    checkCudaErrors(cudaFree(fb));
}
