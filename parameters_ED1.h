/* 
 * File:   parameters_ED1.h
 * Author: Thomas Bronzwaer
 *
 * Created on July 7, 2014, 4:03 PM
 * 
 * NOTES ON cgs-Gauss UNITS:
 * 
 * Units are centimeter, gram, second
 */

#ifndef PARAMETERS_ED1_H
#define	PARAMETERS_ED1_H

#include <cstdlib>
#include <complex>
#include <stdio.h>
#include <iostream>

// PARAMETERS
/////////////

#define pi 3.14159
const double C_SPEED =  290.792458;           // Speed of light, m/us

const int    WIDTH           = 800;
const int    HEIGHT          = 600;
//2*frequency*amplitude < c
double speed_amplitude       = 0.1; //fraction of  lightspeed
double frequency = 2/C_SPEED;
const double AMPLITUDE       = speed_amplitude/frequency/2./pi;
const double period_distance       = AMPLITUDE*4;
/* const double CHARGE          = 4.80320425e-10; // Elementary charge, statcoulomb */
const double CHARGE          = 10000000.480320425; 
const int    TIMESTEPS       = 42000;
const double T_INCREMENT     = 0.01/frequency;

const double T_INIT          = 0.0;

const std::string OUTPUT_DIR = "/tmp/electrodynamics";

#endif
