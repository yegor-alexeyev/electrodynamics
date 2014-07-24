/* 
 * File:   ED1.cpp
 * Author: Thomas Bronzwaer
 * Year:   2014
 * 
 * This C++ program creates a numerical plot of the Lienard-Wiechert potential
 * due to a moving charge (x(t) and v(t) specified by user).
 * 
 * Units are cgs-Gauss. 
 */

#include <cstdlib>
#include <complex>
#include "parameters_ED1.h"
#include "utilities.h"
#include "vector3.h"

using namespace std;

// DEFINITIONS
//////////////

#define _i_ complex<double>(0., 1.)

// FUNCTIONS
////////////

vector3 position(double t){
    vector3 X;
    X.x = AMPLITUDE * cos(FREQUENCY * t); // Horiz Osc
    //X.x = 0.5 * ACCELERATION * t * t; // Linear
    //X.y = SCALE * AMPLITUDE * sin(FREQUENCY * t); // Vert Osc
    //X.x = 0.;
    //X.y = VELOCITY * t;
    X.y = 0.;
    X.z = 0.;
    return X;
}

vector3 velocity(double t){
    vector3 V;
    V.x = -FREQUENCY * AMPLITUDE * sin(FREQUENCY * t); // Horiz Osc
    //V.x = ACCELERATION * t; // Linear
    //V.y = SCALE * FREQUENCY * AMPLITUDE * cos(FREQUENCY * t); // Vert Osc
    //V.x = 0.;
    //V.y = VELOCITY;
    V.y = 0.;
    V.z = 0.;
    return V;
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

// MAIN FUNCTION
////////////////

int main() {
    
    // RGB array for creating images
    unsigned char* data = new unsigned char[WIDTH * HEIGHT * 3];
  
    // MAIN COMPUTATION
    ///////////////////
    
    double X, Y, tau;
    double RED,GRE,BLU;
    vector3 r, R;
    //double time = 1.5e-7;
    double time = 0.0001;
    double PHI = 0.;
    vector3 A, A_fourier;
    int q = 0; //img counter
    
    // for all times
    for (int t = 0; t < TIMESTEPS; t++){
        // For all pixels...
        for (int j = 0; j < HEIGHT; j++){
            for (int i = 0; i < WIDTH; i++){
                
                // COMPUTE POSITION & RETARDED TIME

                // Position coordinates at each pixel center
                r.x = SCALE * ((double) i + .5 - WIDTH / 2.);
                //r.y = SCALE * ((double) j + .5 + 6. * HEIGHT / 1.);
                r.y = SCALE * ((double) j + .5 - HEIGHT / 2.);
                r.z = 0.;                
                
                tau = ret_time(time, r);               
                R = r - position(tau);      
                
                // SCALAR POTENTIAL, PHI
                ////////////////////////
                
                //PHI = 0.;               
                //if (sqrt(dot(R,R)) < c * time){
                if (1){
                    PHI = CHARGE / (norm(R) - dot(R, velocity(tau))/c);
                }
                
                // VECTOR POTENTIAL, A
                //////////////////////
                
                A = velocity(tau) / c * PHI;
                
                double k = 0.62875350658 / (2. * M_PI);
                vector3 xhat = vector3(1., 0., 0.);
                complex<double> complexterm = (-_i_ * exp(_i_ * (k * norm(r) - FREQUENCY * time)));
                double cterm = complexterm.real();
                A_fourier = xhat * k/norm(r) * CHARGE * AMPLITUDE * cos(FREQUENCY * time)
                        * cterm;
                
                // Set RGB at this location
                
                double factor = 56.;
                double color = norm(A_fourier);
                
                RED = factor * color;
                GRE = factor * color;
                BLU = factor * color;

                RED = 255. * pow(RED, GAMMA);
                GRE = 255. * pow(GRE, GAMMA);
                BLU = 255. * pow(BLU, GAMMA);

                data[3 * (WIDTH * j + i) + 0] = min(RED,255.);
                data[3 * (WIDTH * j + i) + 1] = min(GRE,255.);
                data[3 * (WIDTH * j + i) + 2] = min(BLU,255.); 
            }
        }
        
        write_image(data, q);
        q++;       
        time += T_INCREMENT;
    }

    // Cleanup & Exit
    /////////////////
    
    delete[] data;
    
    return 0;
}

