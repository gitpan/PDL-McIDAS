//
//  Copyright (c) 1995-1997 University Corporation for Atmospheric Research
// All rights reserved
//
/* header definitions from numerical recipes */

static int imaxarg1, imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
(imaxarg1) : (imaxarg2))

static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
(iminarg1) : (iminarg2))

void spline(double[], double[], int, double, double, double[]);
void splint(double[], double[], double[], int, double, double *);
void hunt  (double[], int, double, int *);


/*
 * Largest permissible 2nd derivative for spline
 * interpolation.
 */
#define YPLIMIT  100   

