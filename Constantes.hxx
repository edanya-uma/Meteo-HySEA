#ifndef _CONSTANTES_H_
#define _CONSTANTES_H_

#include <iostream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define EPSILON        DBL_EPSILON
#define EARTH_RADIUS   6378136.6
#define MAX_FAULTS     5000   // Maximum number of states in the dynamic deformation file
#define SPONGE_SIZE    4      // Size of the sponge layer (>= 0)
#define SEA_LEVEL      0.0    // Sea level in meters. Used in sponge layer
#define DEFLATE_LEVEL  5      // Level of compression of the NetCDF files (0-9)
//#define EARTH_ANGULAR_VELOCITY 0.0
#define EARTH_ANGULAR_VELOCITY	7.2921159e-5	//NUEVO
#define WATER_DENSITY	1027.0	//kg/(m^3)	NUEVO
#define FRMAX          10.0

#define SEA_SURFACE_FROM_FILE 0
#define OKADA_STANDARD        1
#define DYNAMIC_DEFORMATION   6

using namespace std;

#define NUM_HEBRASX_ARI 8
#define NUM_HEBRASY_ARI 8
#define NUM_HEBRASX_EST 8
#define NUM_HEBRASY_EST 8
#define NUM_HEBRAS_PUNTOS 128

#define iDivUp(a,b)  (((a)%(b) != 0) ? ((a)/(b) + 1) : ((a)/(b)))

typedef struct {
	double2 *areaYCosPhi;
	double *anchoVolumenes;
	double *altoVolumenes;
	double *longitud;
	double *latitud;
} tipoDatosSubmalla;

#endif

