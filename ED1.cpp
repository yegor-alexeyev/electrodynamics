/* 
 * File:   ED1.cpp
 * Author: Thomas Bronzwaer
 * Year:   2014
 * 
 * This C++ program creates a numerical plot of the Lienard-Wiechert potential
 * due to a moving charge (x(t) and v(t) specified by user).
 * 
 * The same field is also computed using a Fourier method (long-wavelength
 * approximation).
 * 
 * Units are cgs-Gauss. 
 * 
 * Copyright (c) 2014 Thomas Bronzwaer
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy 
 * of this software and associated documentation files (the "Software"), to deal 
 * in the Software without restriction, including without limitation the rights 
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
 * copies of the Software, and to permit persons to whom the Software is 
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in 
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
 * SOFTWARE.
 */

#include <cstdlib>
#include <complex>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>

#include "parameters_ED1.h"
#include "utilities.h"
#include "vector3.h"

using namespace std;

// DEFINITIONS
//////////////

#define _i_ complex<double>(0., 1.) // Imaginary unit
#define c   2.99792458e8           // Speed of light, m/s

// FUNCTIONS
////////////
vector3 position(double t){
    vector3 X;
    X.x = WIDTH/2 + AMPLITUDE * sin(FREQUENCY * t); // Horiz Osc
    //X.x = 0.5 * ACCELERATION * t * t; // Linear
    //X.y = SCALE * AMPLITUDE * sin(FREQUENCY * t); // Vert Osc
    //X.x = 0.;
    /* X.x = 400.; */
    /* X.y = VELOCITY * t - VELOCITY * T_INIT - 300.; */
    //X.y = 0.;
    /* X.y = HEIGHT/2 + AMPLITUDE * cos(FREQUENCY * t); */
    /* X.x = WIDTH/2.; */
    X.y = HEIGHT/2.;
    X.z = 0.;
    return X;
}

vector3 velocity(double t){
    
    // If tau - t <= 0, we're asking for velocity before the simulation begins!
    // And we assume the particle is STATIONARY before t = 0
    // NOTE: we allow a small tolerance due to the binary search algorithm's
    // finite precision.
    /* if (tau <= 1.e-12) */
    /*     return vector3(0., 0., 0.); */
    
    vector3 V;
    /* V.x = -FREQUENCY * AMPLITUDE * sin(FREQUENCY * tau); // Horiz Osc */
    //V.x = ACCELERATION * t; // Linear
    //V.y = SCALE * FREQUENCY * AMPLITUDE * cos(FREQUENCY * t); // Vert Osc
    V.x = AMPLITUDE * FREQUENCY*cos(FREQUENCY * t); // Horiz Osc
    /* V.x = 0.; */
    /* V.y = 0.; */
    /* V.y = HEIGHT/2 - AMPLITUDE * FREQUENCY*sin(FREQUENCY * t); */
    /* V.y = VELOCITY; */
    /* V.x = 0.; */
    V.y = 0.;
    V.z = 0.;
    return V;
}

// We need this to compute the E-field
vector3 acceleration(double t){
    
    // If tau - t <= 0, we're asking for accel. before the simulation begins!
    // And we assume the particle is STATIONARY before t = 0
    /* if (tau <= 1.e-12) */
    /*     return vector3(0., 0., 0.); */
    
    vector3 a;
    /* a.x = -FREQUENCY * FREQUENCY * AMPLITUDE * cos(FREQUENCY * tau); */
    a.x = -AMPLITUDE * FREQUENCY*FREQUENCY*sin(FREQUENCY * t); // Horiz Osc
    /* a.x = 0.; */
    a.y = 0.;
    a.z = 0.;
    return a;
}

typedef double* matrix;

void process( double x0, double y0, matrix result, double t, double rendered_time, int dx, int dy)
{
    /* result[3 * (WIDTH * y + x) + 0] = (int)(t/8.5e-12); */


    const int nearby_x = (int)(x0 + 0.5) + dx;
    const int nearby_y = (int)(y0 + 0.5) + dy;
    if (nearby_x >= 0 && nearby_y >= 0 && nearby_x < WIDTH && nearby_y < HEIGHT) {
        result[WIDTH * nearby_y + nearby_x] = t;
    }
}

void processOctet(double x0, double y0, int x, int y, matrix result, double t, double rendered_time) {
		process(x0, y0, result, t, rendered_time, x, y);
		process(x0, y0, result, t, rendered_time, y, x);
		process(x0, y0, result, t, rendered_time, -x, y);
		process(x0, y0, result, t, rendered_time, -y, x);
		process(x0, y0, result, t, rendered_time, -x, -y);
		process(x0, y0, result, t, rendered_time, -y, -x);
		process(x0, y0, result, t, rendered_time, x, -y);
		process(x0, y0, result, t, rendered_time, y, -x);
}


void process_wave_front(double retarded_time, int radius, matrix result, double rendered_time) {
    int x = 0.0;
    int y = radius;
    const vector3 pos = position(retarded_time);
    double x0 = pos.x;
    double y0 = pos.y;

    while (x <= y) {
        processOctet(x0,y0,x,y, result, retarded_time, rendered_time);
        if ( (2*x+1)*(2*x+1) + (2*y-1)*(2*y-1) <= 4*radius*radius ) {
            x++;
        } else {
            y--;
        }
    }
}


//given a particle trajectory and a timestamp calculate a map  of retarded times for each discrete point in [0..WIDTH,0..HEIGHT]
void calculate_retarded_times(double rendered_time, matrix result)
{
    for (int wave_front_radius = 0; wave_front_radius <= (WIDTH+HEIGHT); wave_front_radius++)
    {
        double retarded_time = rendered_time - wave_front_radius/c;
        process_wave_front(retarded_time, wave_front_radius, result, rendered_time);
    }

}



// Returns the retarded time at the current position/time. 
// Uses binary search between 0 and t.
double ret_time(double t, vector3 r){
    double tau = t / 2.;   
    double divider = 4.;
    
    // Perform binary search NSTEPS times; accuracy ~t/2^NSTEPS
    for (int i = 2; i < NSTEPS; i++){
        if (tau < (t - separation(r, position(tau))/c)){
            tau += t / divider;
        }
        else{
            tau -= t / divider;
        }
        divider *= 2.;
    }
    
    return tau;
}

// Returns the advanced time at the current position/time. 
// Uses binary search between t and 2t ??????
double adv_time(double t, vector3 r){
    
    
    // DOESN'T WORK YET
    
    
    double tau = t * 1.5;   
    double divider = 4.;
    
    // Perform binary search NSTEPS times; accuracy ~t/2^NSTEPS
    for (int i = 2; i < NSTEPS; i++){
        if (tau < (t + separation(r, position(tau))/c)){
            tau += t / divider;
        }
        else{
            tau -= t / divider;
        }
        divider *= 2.;
    }
    
    return tau;
}

// MAIN FUNCTION
////////////////
void display_image(unsigned char* data, int imagecounter)
{

    cv::Mat image(HEIGHT, WIDTH, CV_8UC3, data);
    imshow( "rendering", image );
}

int main() {
    
    // INITIALIZE VARIABLES
    ///////////////////////
    
    // RGB array for creating images
    unsigned char* bitmap = new unsigned char[WIDTH * HEIGHT*3];

    matrix data = new double[WIDTH * HEIGHT];   
    
    /* double k = FREQUENCY / c;           // Wave number */
    /* vector3 xhat = vector3(1., 0., 0.); */
    /* double X, Y, tau, gamma; */
    /* vector3 r, R, R_hat, V, a; */
    /* double PHI, PHI_fourier; */
    /* vector3 A, A_fourier, E, B;         // Quantities of interest */
    int q = 0;                          // Image counter
    /* double RED,GRE,BLU; */
    /* vector3 *A_field = new vector3[WIDTH * HEIGHT]; */
    /* vector3 *E_field = new vector3[WIDTH * HEIGHT]; */
    /* double *PHI_field = new double[WIDTH * HEIGHT]; */
    
    // Set initial time
    double time = T_INIT;
    
    // MAIN COMPUTATION
    ///////////////////
    
    cv::namedWindow( "rendering");

    // for all time steps...

    for (int t = 0; t < TIMESTEPS; t++){        
        memset(data, 0, WIDTH * HEIGHT * sizeof(double));
        calculate_retarded_times(time, data);

        for (int j = 0; j < HEIGHT; j++){
            for (int i = 0; i < WIDTH; i++){               
                double retarded_time = data[WIDTH * j + i];
                if (retarded_time == 0.0) {
                    data[WIDTH * j + i] = 0;
                    continue;
                }
                vector3 r;
                r.x = i;
                r.y = j;
                r.z = 0;
                vector3 R = r - position(retarded_time);
                vector3 V = velocity(retarded_time);
                vector3 a = acceleration(retarded_time);
                const double PHI = 1.0 / (norm(R) - dot(R, V)/c);

                //reusing the matrix
                /* data[WIDTH * j + i] = PHI; */

                vector3 R_hat = normalize(R);
                double mu = norm(V)/c;
                double gamma = 1. / sqrt(1. - mu * mu);
                vector3 E = 
                    (R_hat - (V/c)) * CHARGE / (gamma * gamma *
                    pow((1. - dot(R_hat, (V/c))),3.) * norm(R) * norm(R)) 
                    + 
                    (cross(R_hat, cross((R_hat - V/c), a)) / 
                    (pow((1. - dot(R_hat, (V/c))),3.) * norm(R))) * CHARGE / (c*c)
                    ;
                /* data[WIDTH * j + i] = sqrt(10*norm(E)); */
                data[WIDTH * j + i] = E.x;

            }
        }

        // First pass over all pixels to compute the vector potential...
        /* for (int j = 0; j < HEIGHT; j++){ */
        /*     for (int i = 0; i < WIDTH; i++){ */               
        /*         // COMPUTE POSITION & RETARDED TIME */
        /*         /////////////////////////////////// */

        /*         // Position coordinates at each pixel center */
        /*         r.x = SCALE * ((double) i + .5 - WIDTH / 2.); */
        /*         //r.y = SCALE * ((double) j + .5 + 3413.15 * HEIGHT / 1.); */
        /*         r.y = SCALE * ((double) j + .5 - HEIGHT / 2.); */
        /*         r.z = 0.; */                
                
        /*         tau = ret_time(time, r); */               
        /*         R = r - position(tau); */  
        /*         V = velocity(tau); */
        /*         a = acceleration(tau); */
                
        /*         // COMPUTE LIÉNARD-WIECHERT SCALAR POTENTIAL (PHI) */
        /*         ////////////////////////////////////////////////// */

        /*         PHI = CHARGE / (norm(R) - dot(R, V)/c); */
        /*         PHI_field[WIDTH * j + i] = PHI; */

        /*         // COMPUTE LIÉNARD-WIECHERT VECTOR POTENTIAL (A) */
        /*         //////////////////////////////////////////////// */
                
        /*         A = V / c * PHI; */
        /*         A_field[WIDTH * j + i] = A; */
                
        /*         // FOURIER VECTOR POTENTIAL (A_fourier) */
        /*         /////////////////////////////////////// */
                
        /*         complex<double> complexterm = (-_i_ * exp(_i_ * (k * norm(r) - */ 
        /*                                       FREQUENCY * time))); */
        /*         double realpart = complexterm.real(); */
        /*         A_fourier = xhat * k/norm(r) * CHARGE * AMPLITUDE * realpart; */  
                
        /*         // COMPUTE ELECTRIC FIELD */
        /*         ///////////////////////// */
                
        /*         R_hat = normalize(R); */
        /*         gamma = 1. / sqrt(1. - norm(V) * norm(V) / (c * c)); */
        /*         E = (R_hat - (V/c)) * CHARGE / (gamma * gamma * */
        /*             pow((1. - dot(R_hat, (V/c))),3.) * norm(R) * norm(R)) + */ 
        /*             (cross(R_hat, cross((R_hat - V/c), a)) / */ 
        /*             (pow((1. - dot(R_hat, (V/c))),3.) * norm(R))) * CHARGE / (c*c); */
        /*         E_field[WIDTH * j + i] = E; */

                
        /*     } */      
        /* } */   
        
        
        
        // NOTE - we only need one pass - just compute B from E
        
        
        
        
        // Second pass over all pixels - first we fill A, then we compute curl A
        /* for (int j = 0; j < HEIGHT; j++){ */
        /*     for (int i = 0; i < WIDTH; i++){ */
        /*         // COMPUTE B = curl A */
        /*         ///////////////////// */
                
        /*         // Reinitialize B-vector */
        /*         B.x = 0.; */
        /*         B.y = 0.; */
        /*         B.z = 0.; */
                
        /*         // We ignore the outermost border (b/c of spatial derivative) */
        /*         if (i > 1 && i < (WIDTH-1) && j > 1 && j < (HEIGHT-1)){ */
        /*             // Assume d/dz = 0 (we are in the XY-plane and the fields */  
        /*             // are symmetric about this plane) */
        /*             // We need the following derivatives; */
        /*             double dAzdx = (A_field[WIDTH * j + (i+1)].z - */ 
        /*                             A_field[WIDTH * j + (i-1)].z) / (2.*SCALE); */
        /*             double dAzdy = (A_field[WIDTH * (j+1) + i].z - */ 
        /*                             A_field[WIDTH * (j-1) + i].z) / (2.*SCALE); */
        /*             double dAydx = (A_field[WIDTH * j + (i+1)].y - */ 
        /*                             A_field[WIDTH * j + (i-1)].y) / (2.*SCALE); */
        /*             double dAxdy = (A_field[WIDTH * (j+1) + i].x - */ 
        /*                             A_field[WIDTH * (j-1) + i].x) / (2.*SCALE); */
                    
        /*             // Update B-vector */
        /*             B.x = dAzdy - 0.; */
        /*             B.y = 0. - dAzdx; */
        /*             B.z = dAydx - dAxdy; */             
        /*         } */
                
        /*         // COMPUTE PIXEL COLOR (PLOT) */
        /*         ///////////////////////////// */
                
        /*         double factor = 1.25e14;//6.e10; */
        /*         double color = norm(A_field[WIDTH * j + i]); */
                
        /*         // Draw phi */
        /*         //RED = factor * PHI_field[WIDTH * j + i]; */
        /*         //GRE = factor * PHI_field[WIDTH * j + i]; */
        /*         //BLU = factor * PHI_field[WIDTH * j + i]; */
                
        /*         // Draw the norm of A */
        /*         //RED = factor * color; */
        /*         //GRE = factor * color; */
        /*         //BLU = factor * color; */
                
        /*         // Draw A itself (R = x, G = y, B = z) */
        /*         //RED = factor * abs(A_field[WIDTH * j + i].x); */
        /*         //GRE = factor * abs(A_field[WIDTH * j + i].y); */
        /*         //BLU = factor * abs(A_field[WIDTH * j + i].z); */
                
        /*         // Draw E */ 
        /*         //RED = 0.5 + factor * E_field[WIDTH * j + i].x; */
        /*         //GRE = 0.5 + factor * E_field[WIDTH * j + i].y; */
        /*         //BLU = 0.5 + factor * E_field[WIDTH * j + i].z; */
        /*         // */
                
        /*         // Draw B */
        /*         RED = 0.5 + factor * B.x; */
        /*         GRE = 0.5 + factor * B.y; */
        /*         BLU = 0.5 + factor * B.z; */
        /*         // */
                
        /*         // Draw E + B */
        /*         /1* RED = 0.5 * (0.5 + factor * E_field[WIDTH * j + i].x + *1/ */ 
        /*         /1*              0.5 + factor * B.x); *1/ */
        /*         /1* GRE = 0.5 * (0.5 + factor * E_field[WIDTH * j + i].y + *1/ */ 
        /*         /1*              0.5 + factor * B.y); *1/ */
        /*         /1* BLU = 0.5 * (0.5 + factor * E_field[WIDTH * j + i].z + *1/ */ 
        /*         /1*              0.5 + factor * B.z); *1/ */
                
        /*         // Draw norm(B) */
        /*         //RED = factor * norm(B); */
        /*         //GRE = factor * norm(B); */
        /*         //BLU = factor * norm(B); */

        /*         // Gamma correction & color clamping */
        /*         RED = 255. * pow(RED, GAMMA); */
        /*         GRE = 255. * pow(GRE, GAMMA); */
        /*         BLU = 255. * pow(BLU, GAMMA); */

        /*         data[3 * (WIDTH * j + i) + 0] = max(0., min(RED,255.)); */
        /*         data[3 * (WIDTH * j + i) + 1] = max(0., min(GRE,255.)); */
        /*         data[3 * (WIDTH * j + i) + 2] = max(0., min(BLU,255.)); */ 
                
        /*     } */
        /* } */
        
        for (int j = 0; j < HEIGHT; j++){
            for (int i = 0; i < WIDTH; i++){
                const double value = data[WIDTH * j + i];
                const double scaled_value = value*1000000;
                if (scaled_value < 0) {
                    bitmap[3 * (WIDTH * j + i) + 2] = MIN(255, -scaled_value);
                } else {
                    bitmap[3 * (WIDTH * j + i) + 1] = MIN(255, scaled_value);
                }
            }
        }
        display_image(bitmap, q);
        cv::waitKey(1);                                          
        
        cout << "\nOn time step " << (int) q << " of " << (int) (TIMESTEPS-1);
        
        q++;       
        time += T_INCREMENT;
    }

    // Cleanup & Exit
    /////////////////
    
    delete[] data;
    delete[] bitmap;
    
    return 0;
}
