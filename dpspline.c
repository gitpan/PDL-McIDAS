//
//  Copyright (c) 1995-1997 University Corporation for Atmospheric Research
// All rights reserved
//
/* Spline interpolation routines from Numerical Recipes in C */
/**----------------------------------------------------------------------    
* @name       dplspline.c
* 
* Spline interpolation routine.  Parts adapted from Numerical Recipes
*  
* @author     Doug Hunt
* @version    $Revision: 1.2 $
* @cvsinfo    $Id: dpspline.c,v 1.2 2000/01/05 22:26:37 huntd Exp $
* -----------------------------------------------------------------------*/

#include <stdio.h>
#include "dpextC.h"

#define ERR -999

/**----------------------------------------------------------------------    
* @name       dpspline
*
* Using the function described as y = f(x) for each
* y and x in the input vectors, interpolate at points 'in' (x coordinates)
* resulting in the vector 'out' of y coordinates.  M is the
* size of the interpolation window for each
*
* If ever one of the 'in' coordinates is outside of the 'x' array, 
* return for the corresponding 'y' the error value, -999.
*
* @parameter  
* @ input:    x  -- independent variable
* @           y  -- dependent variable
* @           in -- x coordinates to interpolate to
* @           nx -- size of the x and y arrays
* @           nin -- size of 'in' array
* @           y2limit -- limit for second derivatives
* @           m   -- window size
* @ output:   out -- interpolated value
*-----------------------------------------------------------------------*/
int dpspline (double x[], double y[], double in[], int nx, int nin, double y2limit, int m, double out[]) {

  int i, k, jlo;
  double *vector();

  double * y2;  /* The vector of second derivatives computed by 'spline' */
  double * xp  = ((double *)x)-1;  /* One-based pointer to x array */
  double * yp  = ((double *)y)-1;  /* One-based pointer to y array */

  //for (i=0;i<nx;i++) {
  //  printf ("i=%d\tx[i]=%f\ty[i]=%f\n", i, x[i], y[i]);
  //}

  /* 
   * First call the 'spline' routine to compute the 2nd derivatives 
   */
  y2 = vector(1,nx);
  spline (xp, yp, (int)nx, 1e31, 1e31, y2);

  /* loop through the interpolation 'x' locations in 'in' */
  for (i=0;i<nin;i++) {

    /* 
     * Locate the index 'jlo' into the 'x' array such that
     * x[jlo] < in < x[jlo+1].  This is used to find the
     * appropriate length-four array slices of 'x' and 'y' to 
     * pass to the interpolation routine for the current interpolation
     * location, in
     */
    hunt (xp, (unsigned long)nx, in[i], &jlo);
  
    /* 
     * compute starting point in arrays 'x' and 'y' for four point
     * interpolation window.
     */
    k = IMIN(IMAX(jlo-(m-1)/2,1),nx+1-m);
  
    /*
     * If the input value lies on either end point of the
     * x array, the result is equal to that end point.
     */
    if (xp[1] == in[i]) {
    
      out[i] = yp[1];
    
    } else if (xp[nx] == in[i]) {
    
      out[i] = yp[nx];
    
      /* 
       * If interpolation point not contained in 'x' array, then
       * assign the output value to the error value.
       */
    } else if ((jlo == 0) || (jlo == nx)) { 
      
      out[i] = ERR; 
    
      /*
       * If the second derivative for this point is too large, set
       * to the error value.
       */
    } else if (abs(y2[k]) > y2limit) {
      
      out[i] = ERR;
    
      /*
       * If the y value at this point is the error value, set this interpolated
       * point to the error value.
       */
    } else if ((yp[k] == ERR) || (yp[k-1+m] == ERR)) {
    
      out[i] = ERR;
    
      /* 
       * All looks good, do the interpolation...
       */
    } else {
      
      /*
       * Do the interpolation
       */
      splint(&xp[k], &yp[k], &y2[k], m, in[i], &out[i]);
      
    /*
      printf ("i = %d, inp[i] = %f, outp[i] = %f\n", i, inp[i], outp[i]);
      printf ("k-1 = %d, xp[k] = %f, yp[k] = %f, y2[k] = %f\n", k-1, xp[k], yp[k], y2[k]);
    */
      
    }  /* jlo out of range or 2nd derivative too large */

  }
  
  /*
   * Free the y2 array
   */
  free_vector (y2, 1, nx); 
  
  return 0;  /* success */
  
}
