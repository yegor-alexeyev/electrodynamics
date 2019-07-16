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
#define c   C_SPEED
#define pi 3.14159

// FUNCTIONS
////////////
vector3 position(double t){
    vector3 X;
/* X.x=WIDTH/2; */
    X.x = WIDTH/2 + AMPLITUDE * sin(2*pi*FREQUENCY * t); // Horiz Osc
    /* X.y = HEIGHT/2 + AMPLITUDE * cos(2*pi*FREQUENCY * t); */
    /* X.x = WIDTH/2.; */
    X.y = HEIGHT/2.;
    X.z = 0.;
    return X;
}

vector3 velocity(double t){
    
    
    vector3 V;
/* V.x = 0; */
    V.x = AMPLITUDE *2*pi* FREQUENCY*cos(2*pi*FREQUENCY * t); // Horiz Osc
    /* V.y =  - AMPLITUDE *2*pi* FREQUENCY*sin(2*pi*FREQUENCY * t); */
    /* V.x = 0.; */
    V.y = 0.;
    V.z = 0.;
    return V;
}

// We need this to compute the E-field
vector3 acceleration(double t){
    
    
    vector3 a;
    /* a.x = 0; */
    a.x = -AMPLITUDE * 4*pi*pi*FREQUENCY*FREQUENCY*sin(2*pi*FREQUENCY * t); // Horiz Osc
    /* a.y = -AMPLITUDE * 4*pi*pi*FREQUENCY*FREQUENCY*cos(2*pi*FREQUENCY * t); // Horiz Osc */
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
                    data[WIDTH * j + i] =0;
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
                vector3 E = 
                    (R_hat*c - V) * CHARGE * (c*c - norm(V)*norm(V))/ (
                    pow((c - dot(R_hat, V)),3.) * norm(R) * norm(R)) 
                    + 
                    (cross(R_hat, cross(R_hat*c - V, a)) / 
                    (pow(c - dot(R_hat, V),3.) * norm(R))) * CHARGE
                    ;
                /* data[WIDTH * j + i] = E.x; */
                data[WIDTH * j + i] = norm(E);
if (i == (WIDTH/2 - 1) && j == (HEIGHT/2))
{
    cout << "R = " << norm(R) << "; E = " << norm(E)  << endl;
}
if (i == (WIDTH/2 - 2) && j == (HEIGHT/2))
{
    cout << "R = " << norm(R) << "; E = " << norm(E)  << endl;
}
if (i == (WIDTH/2 - 3) && j == (HEIGHT/2))
{
    cout << "R = " << norm(R) << "; E = " << norm(E)  << endl;
}

            }
        }

        
        for (int j = 0; j < HEIGHT; j++){
            for (int i = 0; i < WIDTH; i++){
                const double value = data[WIDTH * j + i];
                /* const double scaled_value = value; */
//TODO remove sqrt
                const double scaled_value = sqrt(value);
                if (scaled_value < 0) {
                    bitmap[3 * (WIDTH * j + i) + 1] = 0;
                    bitmap[3 * (WIDTH * j + i) + 2] = MIN(255, -scaled_value);
                } else {
                    bitmap[3 * (WIDTH * j + i) + 1] = MIN(255, scaled_value);
                    bitmap[3 * (WIDTH * j + i) + 2] = 0;
                }
            }
        }
        display_image(bitmap, q);
        cv::waitKey(1);                                          
        
        cout << "\nOn time step " << (int) q << " of " << (int) (TIMESTEPS-1) << endl;
        
        q++;       
        time += T_INCREMENT;
    }

    // Cleanup & Exit
    /////////////////
    
    delete[] data;
    delete[] bitmap;
    
    return 0;
}
