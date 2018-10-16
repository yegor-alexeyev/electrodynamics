/* 
 * File:   utilities.h
 * Author: Thomas
 *
 * Created on March 21, 2014, 12:43 AM
 * 
 * This header file contains the following utility functions for 3DPauli:
 * 
 * Complex to RGB conversion - for plotting complex functions
 */

#ifndef UTILITIES_H
#define	UTILITIES_H

#include <sstream>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "parameters_ED1.h"
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

using namespace std;

// Returns RGB values for a given complex number
void complex_to_rgb(const complex<double> C, double RGB[]){
    double magnitude = abs(C);
    double phase = arg(C);

    // Normalize the phase.
    phase += M_PI;
    phase /= (2. * M_PI);

    // Normalize the magnitude.
    magnitude /= 1.0;

    // HSV to RGB method from Wikipedia 
    ///////////////////////////////////
    
    double Hprime = phase * 6.;	
    double Ch = magnitude * 1.;	
    double X = Ch * (1. - abs(fmod(Hprime, 2.) - 1.));

    X = max(0., min(255., X));
    Ch = max(0., min(255., Ch));

    if (0.0 <= Hprime && Hprime <= 1.0)
        { RGB[0] = Ch; RGB[1] = X; RGB[2] = 0.; }
    else if (1.0 <= Hprime && Hprime <= 2.0)
        { RGB[0] = X; RGB[1] = Ch; RGB[2] = 0.; }
    else if (2.0 <= Hprime && Hprime <= 3.0)
        { RGB[0] = 0.; RGB[1] = Ch; RGB[2] = X; }
    else if (3.0 <= Hprime && Hprime <= 4.0)
        { RGB[0] = 0.; RGB[1] = X; RGB[2] = Ch; }
    else if (4.0 <= Hprime && Hprime <= 5.0)
        { RGB[0] = X; RGB[1] = 0.; RGB[2] = Ch; }
    else if (5.0 <= Hprime && Hprime <= 6.0)
        { RGB[0] = Ch; RGB[1] = 0.; RGB[2] = X; }
    else 
        { RGB[0] = 0.; RGB[1] = 0.; RGB[2] = 0.; }
    
    // END OF MATERIAL FROM WIKIPEDIA
    /////////////////////////////////
}

static unsigned char * pixel_at (unsigned char * data, int x, int y)
{
    return data + (WIDTH * y + x)*3;
}

/* Write "bitmap" to a PNG file specified by "path"; returns 0 on
   success, non-zero on error. */

static int save_png_to_file (unsigned char *data, const char *path)
{
    FILE * fp;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    size_t x, y;
    png_byte ** row_pointers = NULL;
    /* "status" contains the return value of this function. At first
       it is set to a value which means 'failure'. When the routine
       has finished its work, it is set to a value which means
       'success'. */
    int status = -1;
    /* The following number is set by trial and error only. I cannot
       see where it it is documented in the libpng manual.
    */
    int pixel_size = 3;
    int depth = 8;
    
    fp = fopen (path, "wb");
    if (! fp) {
        goto fopen_failed;
    }

    png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (png_ptr == NULL) {
        goto png_create_write_struct_failed;
    }
    
    info_ptr = png_create_info_struct (png_ptr);
    if (info_ptr == NULL) {
        goto png_create_info_struct_failed;
    }
    
    /* Set up error handling. */

    if (setjmp (png_jmpbuf (png_ptr))) {
        goto png_failure;
    }
    
    /* Set image attributes. */

    png_set_IHDR (png_ptr,
                  info_ptr,
                  WIDTH,
                  HEIGHT,
                  depth,
                  PNG_COLOR_TYPE_RGB,
                  PNG_INTERLACE_NONE,
                  PNG_COMPRESSION_TYPE_DEFAULT,
                  PNG_FILTER_TYPE_DEFAULT);
    
    /* Initialize rows of PNG. */

    row_pointers = (png_byte**)png_malloc (png_ptr, HEIGHT * sizeof (png_byte *));
    for (y = 0; y < HEIGHT; y++) {
        png_byte *row = 
            (png_byte*)png_malloc (png_ptr, sizeof (uint8_t) * WIDTH * pixel_size);
        row_pointers[y] = row;
        for (x = 0; x < WIDTH; x++) {
            unsigned char * pixel = pixel_at (data, x, y);
            *row++ = *pixel;
            *row++ = *(pixel + 1);
            *row++ = *(pixel + 2);
        }
    }
    
    /* Write the image data to "fp". */

    png_init_io (png_ptr, fp);
    png_set_rows (png_ptr, info_ptr, row_pointers);
    png_write_png (png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    /* The routine has successfully written the file, so we set
       "status" to a value which indicates success. */

    status = 0;
    
    for (y = 0; y < HEIGHT; y++) {
        png_free (png_ptr, row_pointers[y]);
    }
    png_free (png_ptr, row_pointers);
    
 png_failure:
 png_create_info_struct_failed:
    png_destroy_write_struct (&png_ptr, &info_ptr);
 png_create_write_struct_failed:
    fclose (fp);
 fopen_failed:
    return status;
}

int write_image(unsigned char *data, int imagecounter){
    stringstream filename;
    filename.str("");
    if (imagecounter<10)
        filename << OUTPUT_DIR << "/frame0000" << (int) imagecounter << ".png";
    else if (imagecounter<100)
        filename << OUTPUT_DIR << "/frame000" << (int) imagecounter << ".png";
    else if (imagecounter<1000)
        filename <<  OUTPUT_DIR << "/frame00" << (int) imagecounter << ".png";
    else if (imagecounter<10000)
        filename <<  OUTPUT_DIR << "/frame0" << (int) imagecounter << ".png";
    else
        filename << OUTPUT_DIR << "/frame" << (int) imagecounter << ".png";

    save_png_to_file(data, filename.str().c_str());


    return 0;
}

#endif
