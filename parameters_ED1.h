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

const int    WIDTH           = 800;
const int    HEIGHT          = 600;
const double SCALE           = 2.0;            // cm; scale factor
const double VELOCITY        = 2.448e10;         // cm/s Linear velocity, 0.8333 c
const double ACCELERATION    = 4.e18;          // Used for linear acceleration
//2*pi*frequency*amplitude < c
const double FREQUENCY       = 4.29e6;
const double AMPLITUDE       = 10;
/* const double CHARGE          = 4.80320425e-10; // Elementary charge, statcoulomb */
const double CHARGE          = 4.80320425e-1; // Elementary charge, statcoulomb
const int    NSTEPS          = 50;             // For binary search (ret. time)
const int    TIMESTEPS       = 42000;
const double T_INCREMENT     = 2.5e-8;
const double GAMMA           = 1. / 2.2;

const double T_INIT          = 0.0;

const std::string OUTPUT_DIR = "/tmp/electrodynamics";

#endif
