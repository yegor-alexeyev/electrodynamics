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


double sine(double arc_amount) {
    return sin(arc_amount*2.0*pi);
}

double cosine(double arc_amount) {
    return cos(arc_amount*2.0*pi);
}


struct Particle
{
double phase;
double charge;
double start_position;
    explicit Particle(double charge_, double phase_, double start_position_) {
        phase = phase_;
        charge = charge_;
        start_position = start_position_;
    }

    vector3 position(double t) const {
        vector3 X;
        X.x = start_position - speed_amplitude * cosine(frequency * t + phase)/frequency/2.0/pi; // Horiz Osc
        X.y = HEIGHT/2.;
        X.z = 0.;
        return X;
    }

    vector3 velocity(double t) const {
        
        
        vector3 V;
        V.x = speed_amplitude * sine(frequency * t + phase); // Horiz Osc
        V.y = 0.;
        V.z = 0.;
        return V;
    }

    // We need this to compute the E-field
    vector3 acceleration(double t) const {
        
        
        vector3 a;
        a.x = speed_amplitude *frequency*2.0*pi*cosine(frequency * t + phase); // Horiz Osc
        a.y = 0.;
        a.z = 0.;
        return a;
    }
};

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


void process_wave_front(double retarded_time, const vector3 retarded_position, int radius, matrix result, double rendered_time) {
    int x = 0.0;
    int y = radius;
    double x0 = retarded_position.x;
    double y0 = retarded_position.y;

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


// MAIN FUNCTION
////////////////
void display_image(unsigned char* data, int imagecounter)
{

    cv::Mat image(HEIGHT, WIDTH, CV_8UC3, data);
    imshow( "rendering", image );
}


void superposition_e_field(const Particle& particle, vector<vector3>& data, const double rendered_time)
{
    matrix retarded_time_data = new double[WIDTH * HEIGHT];   
    memset(retarded_time_data, 0, WIDTH * HEIGHT * sizeof(double));

    for (double wave_front_radius = 0; wave_front_radius <= (WIDTH+HEIGHT); wave_front_radius++)
    {
        double retarded_time = rendered_time - wave_front_radius/1.0; //time is scaled in lightspeed units
        const vector3 retarded_position = particle.position(retarded_time);
        process_wave_front(retarded_time, retarded_position, wave_front_radius, retarded_time_data, rendered_time);
    }

    for (int j = 0; j < HEIGHT; j++) {
        for (int i = 0; i < WIDTH; i++) {               
            double retarded_time = retarded_time_data[WIDTH * j + i];
            if (data[WIDTH * j + i] == vector3(-1,0,0)) {
                continue;
            }
            if (retarded_time == 0.0) {
                data[WIDTH * j + i] =vector3(-1,0,0);
                continue;
            }
            vector3 r;
            r.x = i;
            r.y = j;
            r.z = 0;
            vector3 R = r - particle.position(retarded_time);
            vector3 V = particle.velocity(retarded_time);
            vector3 a = particle.acceleration(retarded_time);


            vector3 R_hat = normalize(R);
            vector3 E = ( 
                (R_hat - V) * ( 1.0 - dot(V,V)) / dot(R,R) 
                + 
                cross(R_hat, cross(R_hat - V, a)) / norm(R)
            )* particle.charge / pow(1.0 - dot(R_hat, V),3.) 
                ;
            data[WIDTH * j + i] = data[WIDTH * j + i] + E;

            /* if (j == HEIGHT/2 && i > WIDTH/2 + 160 && i < WIDTH/2 + 170 ) { */
                /* cout << i << " " << j <<" line: " << E.x << " " << E.y << " " << R.x << " time: " << retarded_time << endl; */
                /* cout << i << " " << j <<" line: " << norm(data[WIDTH*j+i]) << endl; */
                /* cout << i << " " << j <<" line: " << E.x  << " " << E.y<< endl; */
            /* } */

        }
    }
/* cout << "done" << endl; */
}

int main() {
    
    // INITIALIZE VARIABLES
    ///////////////////////
    
    // RGB array for creating images
    unsigned char* bitmap = new unsigned char[WIDTH * HEIGHT*3];

    
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
std::vector<Particle> particles;
for (size_t i = 0; i < 1; i++) {
    particles.push_back(Particle(-CHARGE, 0.0, WIDTH/2));
}
for (size_t i = 0; i < 1; i++) {
    particles.push_back(Particle(CHARGE, 0.5, WIDTH/2 ));
}
    
    // MAIN COMPUTATION
    ///////////////////
    
    cv::namedWindow( "rendering");

    // for all time steps...
    vector<vector3> data(WIDTH * HEIGHT);   

    for (int t = 0; t < TIMESTEPS; t++) {
        for (int i = 0; i < WIDTH*HEIGHT; i++) {
            data[i] = vector3(0,0,0);
        }


for (Particle& particle: particles) {
        superposition_e_field(particle, data, time);
}

        vector<double> euclideanData;
        for (vector3& datum: data) {
            if ( datum == vector3(-1,0,0)){
                euclideanData.push_back(0);
            } else {
                euclideanData.push_back(norm(datum));
            }
        }

        /* cv::normalize(euclideanData, euclideanData , 255.0, 0.0, cv::NORM_INF); */

        for (int i = WIDTH/2 - AMPLITUDE; i < WIDTH/2 + AMPLITUDE; i++) {
            bitmap[3 * (WIDTH * HEIGHT/2 + i)] = 255;
        }
        for (int j = 0; j < HEIGHT; j++) {
            for (int i = 0; i < WIDTH; i++) {
                const double value = euclideanData[WIDTH * j + i];
                const double scaled_value = value/5;
                if (scaled_value < 0) {
                    bitmap[3 * (WIDTH * j + i) + 1] = 0;
                    bitmap[3 * (WIDTH * j + i) + 2] = MIN(255, -scaled_value);
                } else {
                    bitmap[3 * (WIDTH * j + i) + 1] = MIN(255, scaled_value);
                    bitmap[3 * (WIDTH * j + i) + 2] = 0;
                }
for (Particle& particle: particles) {
                if ( int(particle.position(time).x) == i && int(particle.position(time).y) == j) {
                    bitmap[3 * (WIDTH * j + i) + 1] = 255;
                    bitmap[3 * (WIDTH * j + i) + 2] = 255;
                }
}
            }
        }
        display_image(bitmap, q);
        int a = cv::waitKey(1);

if (a == '9') {
speed_amplitude -= C_SPEED*0.01;

}         
if (a == '0') {
speed_amplitude += C_SPEED*0.01;
}
        cout << "\nOn time step " << (int) q << " of " << (int) (TIMESTEPS-1) << endl;
        
        q++;       
        time += T_INCREMENT;
    }

    // Cleanup & Exit
    /////////////////
    
    delete[] bitmap;
    
    return 0;
}
