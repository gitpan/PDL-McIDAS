/**----------------------------------------------------------------------    
* @name       dplog1.c
* 
* Log interpolation routine.  Parts adapted from Numerical Recipes
*  
* @author     Doug Hunt
* @version    $Revision: 1.2 $
* @cvsinfo    $Id: dplog1.c,v 1.2 2000/01/05 22:26:37 huntd Exp $
* -----------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>

#define ERR -999

/**----------------------------------------------------------------------    
* @name       dplog1
*
* Using the function described as y = f(x) for each
* y and x in the input vectors, interpolate at points 'in' (x coordinates)
* resulting in the vector 'out' of y coordinates.
*
* If ever one of the 'in' coordinates is outside of the 'x' array, 
* return for the corresponding 'y' the error value, -999.
* @parameter  
* @ input:    x  -- independent variable
* @           y  -- dependent variable
* @           in -- x coordinate to interpolate to
* @           nx -- size of the x and y arrays
* @ output:   out -- interpolated value
*-----------------------------------------------------------------------*/
void dplog1 (double x[], double y[], double in, int nx, double *out) {

  int i, k, jlo;
  int m = 2; /* degree of interpolation */
  double frac;
  
  double * xp  = ((double *)x)-1;  /* One-based pointer to x Vector */
  double * yp  = ((double *)y)-1;  /* One-based pointer to y Vector */

  /* 
   * Locate the index 'jlo' into the 'x' array such that
   * x[jlo] < in < x[jlo+1].  This is used to find the
   * appropriate length-four array slices of 'x' and 'y' to 
   * pass to the interpolation routine for the current interpolation
   * location, in
   */
  hunt (xp, (unsigned long)nx, in, &jlo);

  /*
   * If the input value lies on either end point of the
   * x array, the result is equal to that end point.
   */
  if (xp[1] == in) {
    
    *out = yp[1];
    
  } else if (xp[nx] == in) {
    
    *out = yp[nx];
    
    /* 
     * If interpolation point not contained in 'x' array, then
     * assign the output value to the error value.
     */
  } else if ((jlo == 0) || (jlo == nx)) { 
    
    *out = ERR; 
    
    /* 
     * If the y value at the interpolation point is the error value, 
     * assign the output value to the error value.
     */
  } else if ((yp[jlo] == ERR) || (yp[jlo+1] == ERR)) {
    
    *out = ERR; 
    
    
    /*
     * If any of the values would cause a log to blow up, return error...
     */
  } else if ((xp[jlo] <= 0) || (xp[jlo+1] <= 0) || (in <= 0)) {
    
    *out = ERR;
    
    
    /* 
     * All looks good, do the interpolation...
     */
  } else {
    
    /* 
     * Compute how far along between two 'x' points the current
     * 'in' value is.
     */
    frac = (log(in) - log(xp[jlo])) / (log(xp[jlo+1]) - log(xp[jlo]));
    
    /*
     * Do the interpolation
     */
    *out = yp[jlo] + frac * (yp[jlo+1] - yp[jlo]);
    
  }  /* jlo out of range*/
  
  return;
  
}
