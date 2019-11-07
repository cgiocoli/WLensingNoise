#ifndef NOISEKAPPA_H_
#define NOISEKAPPA_H_
#include <vector>
#include <gsl/gsl_math.h>
#include <fftw3.h> 
#include "../Moka/distributions.h"
#include "../WeakLMoka/power2D.h"
#include <gsl/gsl_randist.h>

/** 
 * created by:  Margarita Petkova, University of Bologna 2012
 * modified and adapted to MOKA by: Carlo Giocoli, University of Bologna 2012 - (carlo.giocoli@unibo.it)
 */

double* getkappanoise3(std:: vector<double> Pklin, std:: vector<double> llin,
		       double boxlrad,int npix, std:: valarray<float> map,
		       std:: vector<double> Pk, std:: vector<double> l,
		       double filter,long seedi=-1);

double* getkappanoise(std:: vector<double> Pklin, std:: vector<double> llin,
		      double boxlrad,int npix,long seedi=-1);

#endif
