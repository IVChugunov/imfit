/* FILE: func_SersicBPS.cpp ----------------------------------------------- */
/*
 *
 *   Function object class for a photometric model of
 *   B/PS bulge, commonly found in edge-on galaxies.
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x and y.
 *      GetValue() then completes the calculation, using the actual value
 *      of x and y, and returns the result.
 *      So for an image, we expect the user to call Setup() once at
 *      the start, then loop through the pixels of the image, calling
 *      GetValue() to compute the function results for each pixel coordinate
 *      (x,y).
 */

//toreturn open gaussexpx_tmpsave.cpp

/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <gsl/gsl_sf.h>

#include "func_bps.h"
#include "helper_funcs.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int  N_PARAMS = 7;
const char  PARAM_LABELS[][20] = {"PA", "eps", "n", "I_e", "r_e", "k", "re_hx"};
const char  FUNCTION_NAME[] = "Sersic BPS function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char SersicBPS::className[] = "SersicBPS";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

SersicBPS::SersicBPS( )
{
  string  paramName;
  
  nParams = N_PARAMS;
  functionName = FUNCTION_NAME;
  shortFunctionName = className;
  // Set up the vector of parameter labels
  for (int i = 0; i < nParams; i++) {
    paramName = PARAM_LABELS[i];
    parameterLabels.push_back(paramName);
  }
  
  doSubsampling = true;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void SersicBPS::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = (params[0 + offsetIndex] + 90.0) * DEG2RAD;
  eps =  params[1 + offsetIndex];
  n =  params[2 + offsetIndex ];
  I_e = params[3 + offsetIndex ];
  r_e = params[4 + offsetIndex];
  k = params[5 + offsetIndex ];
  re_hx = params[6 + offsetIndex ];

  cosPA = cos(PA);
  sinPA = sin(PA);
  bn = Calculate_bn(n);
  invn = 1.0 / n;

  if (eps > 0) {
    invq = 1.0 / (1.0 - eps);
  } else {
    invq = (1.0 + eps);
  }

  hx = r_e / re_hx; // large re_hx = rapid decline in lobes
  
}

/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// This function calculates and returns the intensity value for a pixel with
// coordinates (x,y), including pixel subsampling if necessary (and if subsampling
// is turned on). The CalculateIntensity() function is called for the actual
// intensity calculation.

double SersicBPS::GetValue( double x, double y )
{

  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp, yp_scaled, x0, r, lobe_factor, I_ser, totalIntensity;

  xp = x_diff * cosPA + y_diff * sinPA;
  yp = -x_diff * sinPA + y_diff * cosPA;
  yp_scaled = yp * invq;

  x0 = std::copysign(yp, xp) / k;

  if (abs(xp) > abs(x0)) { // normal Sersic
    r = sqrt(xp * xp + yp_scaled * yp_scaled);
    lobe_factor = 1;
  } else { // Sersic intensity at ray multiplied by exp decay
    r = sqrt(x0 * x0 + yp_scaled * yp_scaled);
    lobe_factor = exp(-abs((xp - x0) / hx));
  }

  I_ser = I_e * exp(-bn * (pow((r / r_e), invn) - 1.0));
  totalIntensity = I_ser * lobe_factor;

  return totalIntensity;
}

/* END OF FILE: func_bps.cpp ---------------------------------------- */
