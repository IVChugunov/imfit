/* FILE: func_edge-on-broken-disk.cpp ----------------------------------------- */
/* 

// Copyright 2010--2022 by Peter Erwin.
// 
// This file is part of Imfit.
// 
// Imfit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your
// option) any later version.
// 
// Imfit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License along
// with Imfit.  If not, see <http://www.gnu.org/licenses/>.



/* ------------------------ Include Files ------------------------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>

#include "func_edge-on-broken-disk.h"
#include "gsl/gsl_sf_bessel.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

using namespace std;


/* ---------------- Definitions ---------------------------------------- */

const int   N_PARAMS = 7;
const char  PARAM_LABELS[][20] = {"PA", "L_0", "h1", "h2", "r_break", "n", "z_0"};
const char  PARAM_UNITS[][30]  = {"deg (CCW from +y axis)", "counts/voxel",
                                  "pixels", "pixels", "pixels", "", "pixels"};

const char  FUNCTION_NAME[] = "Edge-on Broken Disk function";

const double DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;
const double COSH_LIMIT = 100.0;

const double TABLE_DR = 1.0;        // radial sampling (pixels)

const char EdgeOnBrokenDisk::className[] = "EdgeOnBrokenDisk";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

EdgeOnBrokenDisk::EdgeOnBrokenDisk()
{
  nParams = N_PARAMS;
  functionName = FUNCTION_NAME;
  shortFunctionName = className;

  for (int i = 0; i < nParams; i++) {
    parameterLabels.push_back(PARAM_LABELS[i]);
    parameterUnits.push_back(PARAM_UNITS[i]);
  }

  parameterUnitsExist = true;
  doSubsampling = true;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void EdgeOnBrokenDisk::Setup(double params[], int offsetIndex, double xc, double yc)
{
  x0 = xc;
  y0 = yc;

  PA      = params[0 + offsetIndex];
  L_0     = params[1 + offsetIndex];
  h1      = params[2 + offsetIndex];
  h2      = params[3 + offsetIndex];
  r_break = params[4 + offsetIndex];
  n       = params[5 + offsetIndex];
  z_0     = params[6 + offsetIndex];

  PA_rad = (PA + 90.0) * DEG2RAD;

  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);

  alpha = 2.0 / n;
  scaledZ0 = alpha * z_0;
  two_to_alpha = pow(2.0, alpha);
  Sigma_00_2 = 2.0 * L_0 * h2 * exp(-r_break / h1 + r_break / h2);

  TABLE_SIZE = (int)(r_break / TABLE_DR) + 2; // radial table size
  radialTable.resize(TABLE_SIZE);

  BuildRadialTable();
}


/* ---------------- Luminosity density -------------------------------- */

double EdgeOnBrokenDisk::Luminosity(double R)
{
  if (R < r_break)
    return L_0 * exp(-R / h1);

  return L_0 * exp(-r_break / h1) * exp(-(R - r_break) / h2);
}


/* ---------------- LOS integral using R = r cosh(u) ------------------ */

struct LOSParams
{
  double r;
  double L_0;
  double h1;
  double h2;
  double r_break;
};

/* ---------------- Numerically stable log(cosh(u)) ------------------- */

static inline double LogCosh(double u)
{
  if (u < 100.0)
    return log(cosh(u));
  return u - log(2.0);
}


/* ---------------- LOS integrands ------------------------------- */

double LOSIntegrandInner(double u, void *params)
{
  LOSParams *p = static_cast<LOSParams *>(params);

  const double log_cosh = LogCosh(u);

  // log(R / h1) = log(r / h1) + log(cosh(u))
  const double log_R_over_h = log(p->r / p->h1) + log_cosh;

  const double R_over_h = exp(log_R_over_h);
  const double log_f = log(p->L_0) + log_cosh - R_over_h;

  if (log_f < -100.0)
    return 0.0;

  return exp(log_f);
}


double LOSIntegrandOuter(double u, void *params)
{
  LOSParams *p = static_cast<LOSParams *>(params);

  const double log_cosh = LogCosh(u);

  // R/h2 = (r/h2) cosh(u)
  const double log_R_over_h = log(p->r / p->h2) + log_cosh;

  const double R_over_h = exp(log_R_over_h);

  // L(R) = L_0 exp(-rb/h1) exp(-(R-rb)/h2) = L_0 exp(-rb/h1 + rb/h2) exp(-R/h2)
  const double log_f = log(p->L_0) - p->r_break / p->h1 + p->r_break / p->h2 + log_cosh - R_over_h;

  if (log_f < -100.0)
    return 0.0;

  return exp(log_f);
}

double EdgeOnBrokenDisk::LOSIntegral(double r)
{
  if (r <= 0.0)
    return 2.0 * L_0 * (h1 * (1.0 - exp(-r_break / h1)) + h2 * exp(-r_break / h1));

  LOSParams params;

  params.L_0       = L_0;
  params.h1       = h1;
  params.h2       = h2;
  params.r_break  = r_break;
  params.r        = r;

  gsl_function F_inner;
  gsl_function F_outer;

  F_inner.function = &LOSIntegrandInner;
  F_inner.params   = &params;

  F_outer.function = &LOSIntegrandOuter;
  F_outer.params   = &params;

  gsl_integration_workspace *workspace_inner = gsl_integration_workspace_alloc(1000);
  gsl_integration_workspace *workspace_outer = gsl_integration_workspace_alloc(1000);


  double inner_result = 0.0;
  double inner_error  = 0.0;

  double outer_result = 0.0;
  double outer_error  = 0.0;

  int status_inner = GSL_SUCCESS;
  int status_outer = GSL_SUCCESS;

  if (r < r_break)
  {
    double u_break = acosh(r_break / r);
    status_inner = gsl_integration_qag(&F_inner, 0.0, u_break, 0.0, 1e-8, 1000, GSL_INTEG_GAUSS61, workspace_inner, &inner_result, &inner_error);
    status_outer = gsl_integration_qagiu(&F_outer, u_break, 0.0, 1e-8, 1000, workspace_outer, &outer_result, &outer_error);
  }

  // This is implemented because table is calculated with some reserve.
  else
    status_outer = gsl_integration_qagiu(&F_outer, 0.0, 0.0, 1e-8, 1000, workspace_outer, &outer_result, &outer_error);

  gsl_integration_workspace_free(workspace_inner);
  gsl_integration_workspace_free(workspace_outer);

  if (status_inner != GSL_SUCCESS || status_outer != GSL_SUCCESS)
    fprintf(stderr, "Warning: GSL LOS integration failed at r = %.6g (inner status = %d, outer status = %d)\n", r, status_inner, status_outer);

  return 2.0 * r * (inner_result + outer_result);
}


/* ---------------- Build interpolation table ------------------------- */

void EdgeOnBrokenDisk::BuildRadialTable()
{
  for (int i = 0; i < TABLE_SIZE; i++)
  {
    double r = i * TABLE_DR;
    radialTable[i] = LOSIntegral(r);
  }
}


/* ---------------- Interpolate radial integral ----------------------- */

double EdgeOnBrokenDisk::RadialValue(double r)
{
  double t = r / TABLE_DR;
  int i = (int)t;

  if (i >= TABLE_SIZE - 1)
    return radialTable[TABLE_SIZE - 1];

  double f = t - i;

  return radialTable[i] * (1.0 - f) + radialTable[i + 1] * f;
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ---------------- */

double EdgeOnBrokenDisk::CalculateIntensity(double r, double z)
{
  double I_radial;

  if (r == 0) {
    I_radial = 2.0 * L_0 * (h1 + (h2 - h1) * exp(-r_break / h1));
  }
  else if (r > r_break) {
    double scaledR = r / h2;
    I_radial = Sigma_00_2 * scaledR * gsl_sf_bessel_K1(scaledR);
  }
  else {
    I_radial = RadialValue(r);
  }

  double verticalScaling;

  if ((z / scaledZ0) > COSH_LIMIT)
    verticalScaling = two_to_alpha * exp(-z / z_0);
  else
  {
    double sech = 1.0 / cosh(z / scaledZ0);
    verticalScaling = pow(sech, alpha);
  }

  return I_radial * verticalScaling;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// NOTE: x and y are orthognal image coordinates; R and z are orthogonal
// coordinates *in the component reference frame* (corresponding to xp and
// yp in other, non-edge-on function objects).
// Note that both CalculateIntensity() and CalculateSubsamples() assume that
// R and z are *non-negative*!
double EdgeOnBrokenDisk::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  R, z, totalIntensity;
  int  nSubsamples;
  
  // Calculate R,z (= x,y in component reference frame)
  R = fabs(x_diff*cosPA + y_diff*sinPA);    // "R" is x in the component reference frame
  z = fabs(-x_diff*sinPA + y_diff*cosPA);   // "z" is y in the component reference frame

  nSubsamples = CalculateSubsamples(R, z);
  if (nSubsamples > 1) {
    // Do subsampling
    // start in center of leftmost/bottommost sub-pixel
    double deltaSubpix = 1.0 / nSubsamples;
    double x_sub_start = x - 0.5 + 0.5*deltaSubpix;
    double y_sub_start = y - 0.5 + 0.5*deltaSubpix;
    double theSum = 0.0;
    for (int ii = 0; ii < nSubsamples; ii++) {
      double x_ii = x_sub_start + ii*deltaSubpix;
      for (int jj = 0; jj < nSubsamples; jj++) {
        double y_ii = y_sub_start + jj*deltaSubpix;
        x_diff = x_ii - x0;
        y_diff = y_ii - y0;
        R = fabs(x_diff*cosPA + y_diff*sinPA);
        z = fabs(-x_diff*sinPA + y_diff*cosPA);
        theSum += CalculateIntensity(R, z);
      }
    }
    totalIntensity = theSum / (nSubsamples*nSubsamples);
  }
  else
    totalIntensity = CalculateIntensity(R, z);

  return totalIntensity;
}


/* ---------------- PROTECTED METHOD: CalculateSubsamples ------------------------- */
// Function which determines the number of pixel subdivisions for sub-pixel integration,
// given that the current pixel is a distance of r away from the center of the
// r=0 line and/or z away from the disk plane.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int EdgeOnBrokenDisk::CalculateSubsamples( double r, double z )
{
  int  nSamples = 1;
  int  nSr, nSz;
  double  R_abs = fabs(r);
  double  z_abs = fabs(z);
  
  // based on standard exponential-function subsampling
  if ( (doSubsampling) && ((R_abs < 10.0) || (z_abs < 10.0)) ) {
    if ( ((h1 <= 1.0) && (R_abs <= 1.0)) || ((z_0 <= 1.0) && (z_abs <= 1.0)) ) {
      nSr = min(100, (int)(2 * SUBSAMPLE_R / h1));
      nSz = min(100, (int)(2 * SUBSAMPLE_R / z_0));
      nSamples = max(nSr, nSz);
    }
    else {
      if ((R_abs <= 3.0) || (z_abs <= 3.0))
        nSamples = 2 * SUBSAMPLE_R;
      else {
        nSr = min(100, (int)(2 * SUBSAMPLE_R / R_abs));
        nSz = min(100, (int)(2 * SUBSAMPLE_R / z_abs));
        nSamples = max(nSr, nSz);
      }
    }
  }
  return nSamples;
}


/* END OF FILE: func_edge-on-broken-disk.cpp ---------------------------------- */
