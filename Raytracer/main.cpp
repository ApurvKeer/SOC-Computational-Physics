// Reference : https://locklessinc.com/articles/raytracing/

//Bunch of Libraries
#include <math.h>
#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <iostream>
#include <random>
#include <cstdlib>

#define N	6 //Number of independent variables to be tracked

//Bunch of constants
const double pi = M_PI;
const double a = 0.5;
const int size = 1000;
double inclination = 85;

double r0, theta0;

double a2;
double r_horizon, r_stable_orbit, r_disk;

double L, kappa;

// Coupled differential equations describing motion of photon, null geodesics around a Kerr Black hole
void geodesic(double *y, double *dydx){

	double r, theta, pr, ptheta;

	r = y[0];
	theta = y[1];
	pr = y[4];
	ptheta = y[5];

	double r2 = r*r;
	double twor = 2.0*r;

	double sintheta, costheta;
	sintheta = sin(theta);
	costheta = cos(theta);
	double cos2 = costheta*costheta;
	double sin2 = sintheta*sintheta;

	double sigma = r2+a2*cos2;
	double delta = r2-twor+a2;
	double sd = sigma*delta;
	double siginv = 1.0/sigma;
	double bot = 1.0/sd;

	// Prevent problems along Z azis
	if (sintheta < 1e-8){
		sintheta = 1e-8;
		sin2 = 1e-16;
	}

	dydx[0] = -pr*delta*siginv;
	dydx[1] = -ptheta*siginv;
	dydx[2] = -(twor*a+(sigma-twor)*L/sin2)*bot;
	dydx[3] = -(1.0+(twor*(r2+a2)-twor*a*L)*bot);
	dydx[4] = -(((r-1.0)*(-kappa)+twor*(r2+a2)-2.0*a*L)*bot-2.0*pr*pr*(r-1.0)*siginv);
	dydx[5] = -sintheta*costheta*(L*L/(sin2*sin2)-a2)*siginv;
}

//Setting up initial conditions based on camera view
void initial(double *y0, double *ydot0, double x, double y){

	y0[0] = r0;
	y0[1] = theta0;
	y0[2] = 0;
	y0[3] = 0;
	y0[4] = cos(y)*cos(x);
	y0[5] = sin(y)/r0;

	double sintheta, costheta;
	sintheta = sin(theta0);
	costheta = cos(theta0);
	double cos2 = costheta*costheta;
	double sin2 = sintheta*sintheta;

	double rdot0 = y0[4];
	double thetadot0 = y0[5];

	double r2 = r0 * r0;
	double sigma = r2 + a2*cos2;
	double delta = r2 - 2.0 * r0 + a2;
	double s1 = sigma - 2.0 * r0;

	y0[4]= rdot0*sigma/delta;
	y0[5]= thetadot0*sigma;

	ydot0[0] = rdot0;
	ydot0[1] = thetadot0;
	ydot0[2] = cos(y)*sin(x)/(r0*sin(theta0));

	double phidot0 = ydot0[2];
	double energy2 = s1*(rdot0*rdot0/delta+thetadot0*thetadot0)
					+ delta*sin2*phidot0*phidot0;

	double energy = sqrt(energy2);

	/* Rescale */
	y0[4] = y0[4]/energy;
	y0[5] = y0[5]/energy;

	/* Angular Momentum with E = 1 */
	L = ((sigma*delta*phidot0-2.0*a*r0*energy)*sin2/s1)/energy;

	kappa = y0[5]*y0[5]+a2*sin2+L*L/sin2;

	/* Hack - make sure everything is normalized correctly by a call to geodesic */
	geodesic(y0, ydot0);
}

//One step of the Runge Kutta Fehlberg method
void rkstep(double *y, double *dydx, double h, double *yout, double *yerr){

	int i;

	double ak[N];

	double ytemp1[N], ytemp2[N], ytemp3[N], ytemp4[N], ytemp5[N];

	for (i = 0; i < N; i++){
		double hdydx = h * dydx[i];
		double yi = y[i];
		ytemp1[i] = yi + 0.2 * hdydx;
		ytemp2[i] = yi + (3.0/40.0) * hdydx;
		ytemp3[i] = yi + 0.3 * hdydx;
		ytemp4[i] = yi -(11.0/54.0) * hdydx;
		ytemp5[i] = yi + (1631.0/55296.0) * hdydx;
		yout[i] = yi + (37.0/378.0) * hdydx;
		yerr[i] = ((37.0/378.0)-(2825.0/27648.0)) * hdydx;
	}

	geodesic(ytemp1, ak);

	for (i = 0; i < N; i++){
		double yt = h * ak[i];
		ytemp2[i] += (9.0/40.0) * yt;
		ytemp3[i] -= 0.9 * yt;
		ytemp4[i] += 2.5 * yt;
		ytemp5[i] += (175.0/512.0) * yt;
	}

	geodesic(ytemp2, ak);

	for (i = 0; i < N; i++){
		double yt = h * ak[i];
		ytemp3[i] += 1.2 * yt;
		ytemp4[i] -= (70.0/27.0) * yt;
		ytemp5[i] += (575.0/13824.0) * yt;
		yout[i] += (250.0/621.0) * yt;
		yerr[i] += ((250.0/621.0)-(18575.0/48384.0)) * yt;
	}

	geodesic(ytemp3, ak);

	for (i = 0; i < N; i++){
		double yt = h * ak[i];
		ytemp4[i] += (35.0/27.0) * yt;
		ytemp5[i] += (44275.0/110592.0) * yt;
		yout[i] += (125.0/594.0) * yt;
		yerr[i] += ((125.0/594.0)-(13525.0/55296.0)) * yt;
	}

	geodesic(ytemp4, ak);

	for (i = 0; i < N; i++){
		double yt = h * ak[i];
		ytemp5[i] += (253.0/4096.0) * yt;
		yerr[i] -= (277.0/14336.0) * yt;
	}

	geodesic(ytemp5, ak);

	for (i = 0; i < N; i++){
		double yt = h * ak[i];
		yout[i] += (512.0/1771.0) * yt;
		yerr[i] += ((512.0/1771.0)-0.25) * yt;
	}
}

//The R-K algorithm will dynamically tune the sizes of the steps it takes in order to maintain error within some given tolerance.
double rkqs(double *y, double *dydx, double htry, double escal, double *yscal, double *hdid){

	int i;

	double hnext;

	double errmax, h = htry, htemp;
	double yerr[N], ytemp[N];

	while (true){

		rkstep(y, dydx, h, ytemp, yerr);

		errmax = 0.0;
		for (i = 0; i < N; i++){
			double temp = fabs(yerr[i]/yscal[i]);
			if (temp > errmax) errmax = temp;
		}

		errmax *= escal;
		if (errmax <= 1.0) break;

		htemp = 0.9 * h / sqrt(sqrt(errmax));

		h *= 0.1;

		if (h >= 0.0){
			if (htemp > h) h = htemp;
		}
		else{
			if (htemp < h) h = htemp;
		}
	}

	if (errmax > 1.89e-4){
		hnext = 0.9 * h * pow(errmax, -0.2);
	}
	else{
		hnext = 5.0 * h;
	}

	*hdid = h;

	memcpy(y, ytemp, N * sizeof(double));

	return hnext;
}

//We use a binary search technique by noticing which hemisphere the photon happens to lie in
void binarysearch(double *y, double *dydx, double hbig){

	double hsmall = 0.0;

	int side;
	if (y[1] > pi/2.0){
		side = 1;
	}
	else if (y[1] < pi/2.0){
		side = -1;
	}
	else{
		/* Already at the equator */
		return;
	}

	geodesic(y,dydx);

	while ((y[0] > r_horizon) && (y[0] < r0) && (side != 0)){
		double yout[N], yerr[N];

		double hdiff = hbig - hsmall;

		if (hdiff < 1e-7){
			rkstep(y, dydx, hbig, yout, yerr);

			memcpy(y, yout, N * sizeof(double));

			return;
		}

		double hmid = (hbig + hsmall) / 2;

		rkstep(y, dydx, hmid, yout, yerr);

		if (side * (yout[1] - pi/2.0) > 0){
			hsmall = hmid;
		}
		else{
			hbig = hmid;
		}
	}
}

//color gradient of the accretion disk
void color_ramp(float *rgb, double r){
	int red[3] = {12, 15, 149};
	int orange[3] = {8, 94, 232};
	int black[3] = {0, 0, 0};

	float factor = (r - r_disk)/(r_stable_orbit - r_disk);
	for(int i = 0; i < 3; i++){
		rgb[i] = orange[i]*factor + red[i]*(1 - factor);
	}
	for(int i = 0; i < 3; i++){
		rgb[i] = rgb[i]*factor + black[i]*(1 - factor);
	}
}

//Firing a ray from the observer and checker whether it hit the accretion disk or not
void fire_ray(float *rgb, int x1, int y1){

	double htry = 0.5, escal = 1e11, hdid = 0.0, hnext = 0.0;

	double range = 0.0025 * r_disk / (size - 1.0);

	double y[N], dydx[N], yscal[N], ylaststep[N];

	int side;
	int i;

	initial(y, dydx, (x1 - (size + 1.0) / 2) * range, (y1 - (size + 1.0) / 2) * range);

	while (1){
		memcpy(ylaststep, y, N * sizeof(double));

		geodesic(y, dydx);

		for (i = 0; i < N; i++){
			yscal[i] = fabs(y[i]) + fabs(dydx[i] * htry) + 1.0e-3;
		}

		if (y[1] > pi/2) side = 1;
		else if (y[1] < pi/2) side = -1;
		else side = 0;

		hnext = rkqs(y, dydx, htry, escal, yscal, &hdid);

		if ((y[1]-pi/2)*side < 0){
			memcpy(y, ylaststep, N * sizeof(double));

			binarysearch(y, dydx, hdid);

			/* Did we hit the disk? */
			if ((y[0] <= r_disk) && (y[0] >= r_stable_orbit)){	
				color_ramp(rgb, y[0]);
				return;
			}
		}

		/* Inside the hole, or escaped to infinity */
		if ((y[0] < r_horizon)){
			rgb[0] = 0;
			rgb[1] = 0;
			rgb[2] = 0;
			return;
		}

		if ((y[0] > r0)){ //starry background texture
			if(rand()%100 < 99){
				rgb[0] = 0;
				rgb[1] = 0;
				rgb[2] = 0;
			}
			else{
				rgb[0] = 255;
				rgb[1] = 255;
				rgb[2] = 255;
			}
			return;
		}

		htry = hnext;
	}
}

//////////////////////////////////////

double inner_orbit(void){

	double z1 = 1+cbrt(1-a2)*(cbrt(1+a)+cbrt(1-a));
	double z2 = sqrt(3*a2+z1*z1);
	return 3+z2-sqrt((3-z1)*(3+z1+2*z2));
}

int main(){

	int i, j;
	float color[3];

	r0 = 1000.0;
	theta0 = (pi/180.0) * inclination;

	a2 = a*a;

	r_horizon = 1.0 + sqrt(1.0-a2) + 1e-5;
	r_disk = 40.0; //times the Schwarzschild radius
	r_stable_orbit = inner_orbit();

	//.ppm file header
	std::cout << "P3\n" << size << " " << size << "\n255\n";

	for (j = 0; j < size; j++){
		for (i = 0; i < size; i++){
			fire_ray(color, i, size - j);
			unsigned int b = color[0];
			unsigned int g = color[1];
			unsigned int r = color[2];
			std::cout << r << " " << g << " " << b << "\n";
		}
	}
}

/*To run this file use
$ g++ main.cpp
$ ./a.out>image.ppm
This will store the image in a ppm format. You can open it in VSCode. It has an extension for it.
*/