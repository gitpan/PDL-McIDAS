!/**----------------------------------------------------------------------    
! @name       util.f
!
! Various FORTRAN routines for calling by perl
!
! @author     Doug Hunt
! @version    $Revision: 1.2 $
! @since      11/26/97
! @cvsinfo    $Id: util.f,v 1.2 2000/01/05 22:26:41 huntd Exp $
* -----------------------------------------------------------------------*/


!/**----------------------------------------------------------------------    
! @name       trapezoid
!
! Simple integration using trapezoid method.
! y = f(x) is given at n points
! x must be monotonically increasing
!
! @parameter  
! @ input:    n -- number of values in x and y arrays
! @           x, y -- Function values
! @ output:   integral -- value of integral under y=f(x)
! -----------------------------------------------------------------------*/
      subroutine trapezoid (n, x, y, integral)
      implicit none
      integer n, i
      real*8 x(n), y(n), integral, dx

      integral = 0.d0
      do i = 1, n-1
         dx = x(i+1) - x(i)
         integral = integral + (dx * 0.5 * (y(i+1) + y(i)))
      end do

      return
      end


!/**----------------------------------------------------------------------    
! @name       lreg
!
! Linear regression
!
! @parameter  
! @ input:    n -- input vector size
! @           x -- independent variable
! @           y -- dependent variable
! @ output:   z -- linear fit of y=f(x):  z=linear(x)
! @           c0, c1 -- Linear fit parameters:  z = c0 + c1 * x
! -----------------------------------------------------------------------*/
      subroutine lreg(n,x,y,z,c0,c1)
      implicit none
      integer n
      real*8 x(n),y(n),z(n)
      real*8 c0, c1
      
      real*8 a0, a1, a2, ab, b1
      integer i

      a0=dble(n)
      a1=0.d0
      a2=0.d0
      ab=0.d0
      b1=0.d0

      do i=1,n
         a1=a1+x(i)
         b1=b1+y(i)
         a2=a2+x(i)*x(i)
         ab=ab+x(i)*y(i)
      end do

      c0=(a2*b1-a1*ab)/(a2*a0-a1**2)
      c1=(ab*a0-a1*b1)/(a2*a0-a1**2)

      do i=1,n
         z(i)=c0+c1*x(i)
      end do

      return
      end


!/**----------------------------------------------------------------------    
! @name       lregslope
!
! Linear regression, with the y offset contrained to zero
! (the output line is of the form z = c1 * x, not z = c0 + c1 * x)
!
! @parameter  
! @ input:    n -- input vector size
! @           x -- independent variable
! @           y -- dependent variable
! @ output:   z -- linear fit of y=f(x):  z=linear(x)
! @           c1 -- Linear fit parameters:  z = c1 * x
! -----------------------------------------------------------------------*/
      subroutine lregslope(n,x,y,z,c1)
      implicit none
      integer n
      real*8 x(n),y(n),z(n)
      real*8 c1
      
      real*8 a0, a2, ab
      integer i

      a0=dble(n)
      a2=0.d0
      ab=0.d0

      do i=1,n
         a2=a2+x(i)*x(i)
         ab=ab+x(i)*y(i)
      end do

      c1=ab/a2
      do i=1,n
         z(i)=c1*x(i)
      end do

      return
      end



!/**----------------------------------------------------------------------    
! @name       find_geop
!
! Compute a set of geopotential heights for a data set.
!
! @parameter  
! @ input:    alt -- Base geopotential height (km)
! @           t   -- Array of temperature values (C)
! @           p   -- Array of pressure values (mb)
! @           e   -- Array of water Vp values (mb, -999 if undef)
! @ output:   geo -- list of geopotential heights, in km.
! -----------------------------------------------------------------------*/
      subroutine find_geop(n, alt, t, p, e, geo)
      implicit none
      integer n, i
      real*8 alt, t(n), p(n), e(n), geo(n)
      real*8 R, M, G, X, t1, t2, tv1, tv2

! t[0], p[0], e[0] are assumed to correspond to $alt.
! Reference:  Paper by Chris Rocken: "Range Correction for Atmospheric
!             Propagation of E-M Signals from Space" pg 11.

      R = 8314.36               ! joule/deg K-kg-mol, universal gas constant
      M = 28.966                ! molecular weight of air
      G = 9.80665               ! WMO acceleration of gravity, m / sec**2

      geo(1) = alt
      do i=2,n

         if (t(i).eq.-999.) return
         if (p(i).eq.-999.) return

                                ! Convert to Kelvin
         t1 = t(i-1) + 273.16
         t2 = t(i)   + 273.16

                                ! Compute virtual temperatures, auxillary variable X
         tv1 = t1/(1. - 0.379 * (e(i-1)/p(i-1)))
         tv2 = t2/(1. - 0.379 * (e(i)/p(i)))
         X   = (tv2 - tv1)/tv1
         
         if (X .eq. 0.) then
            geo(i) = geo(i-1) 
         else 
            geo(i) = geo(i-1) + abs(
     +       ((R * tv1) / (G * M)) *
     +       ((X * log(p(i)/p(i-1))) / log (1 + X)) * 0.001) 
         endif
      end do
    
      return
      end


!/**----------------------------------------------------------------------    
! @name       geo2msl
!
! Convert a list of geopotential altitudes to mean sea level altitude.
!
! @parameter  
! @ input:    h   -- ref to list of geopotential altitudes (in km)
! @           lat -- latitude  of profile in degrees.
! @           lon -- longitude of profile in degrees.
! @ output:   z   -- list of MSL altitudes, in km.
! -----------------------------------------------------------------------*/
      subroutine geo2msl (h, lat, lon, z)
      implicit none
      
      real*8 h, lat, lon, z  ! note lon not currently used

      real*8 b, a, pi, G, g0, r0
      
      b = 6356.7516d0             ! min earth radius, km
      a = 6378.1363d0             ! max earth radius, km

      pi = 3.14159265358979d0
      lat = lat * (pi/180.d0)           ! in radians
      lon = lon * (pi/180.d0)           ! in radians

  ! These are the equations for g0 and r0 in Chris Rocken's paper.
  ! I am not using them because they imply a standard ellipsoid which
  ! is quite different from our standard ellipsoid. D. Hunt 10/28/96
  !G = 0.0098
  !g0 = 0.001 * 9.780356 * (1+0.0052885 * (sin(lat))**2 - 5.9e-6 * (sin(2*lat))**2)
  !r0 = (2*g0)/(3.085462e-6 + 2.27e-9 * cos(2*lat) - 2e-12*cos(4*lat))

      G = 0.00980665d0          ! WMO reference g value, km/s**2

      g0 = 0.d0
      call gravity (lat, 0.d0, g0)
      g0 = g0 * 0.00001d0             ! convert to km/s**2
 
  ! compute local earth's radius using ellipse equation
      r0 = sqrt ( a**2 * cos(lat)**2 + b**2 * sin(lat)**2)

      if (h.eq.-999.d0) then
         z = -999.d0
      else 
      ! Compute altitude above sea level
         z = (r0 * h) / (((g0*r0)/G) - h)
      endif
      
      return
      end


!/**----------------------------------------------------------------------    
! @name       msl2geo
!
! Convert a list of geometric mean sea level altitudes to geopotential height
!
! @parameter  
! @ input:    z   -- list of geometric altitudes (in km)
! @           lat -- latitude  of profile in degrees
! @           lon -- longitude of profile in degrees.
! @ output:   h   -- list of geopotential altitudes, in km
! -----------------------------------------------------------------------*/
      subroutine msl2geo (z, lat, lon, h)
      implicit none
      real*8 h, lat, lon, z     ! Note lon is currently not used
      real*8 b, a, pi, g0, r0, G

      b = 6356.7516             ! min earth radius, km
      a = 6378.1363             ! max earth radius, km

      pi = 3.14159265358979d0
      lat = lat * (pi/180.d0)        ! in radians
      lon = lon * (pi/180.d0)        ! in radians

      G = 0.00980665d0          ! WMO reference g value, km/s**2

      g0 = 0.d0
      call gravity (lat, 0.d0, g0)
      g0 = g0 * 0.00001d0       ! convert to km/s**2
 
      ! compute local earth's radius using ellipse equation
      r0 = sqrt ( a**2 * cos(lat)**2 + b**2 * sin(lat)**2)

      if (z.eq.-999.d0) then
         h = -999.d0
      else 
         ! Compute geopotential height
         h = (g0/G) * ((r0 * z)/(r0 + z))
      endif
  
      return
      end


!/**----------------------------------------------------------------------    
! @name       oned
!
! Polynomial interpolation given x = [0,1] and line = values of 
! four points on line around x:
!
!    line[1]  line[2]   x   line[3]  line[4] 
!
! x gives fraction of way between line[2] and line[3] for interpolation.
!
! @parameter  
! @ input:    n   -- size of 'line'
! @           x   -- fraction of way between line[2] and line[3] for interpolation.
! @           line-- Four values to interpolate between
! @ output:   out -- Interpolated value
! -----------------------------------------------------------------------*/
      subroutine oned(n, x, line, out)
      implicit none

      integer n
      real*8 x, line(n)
      real*8 a, b, c, d, out
      real*8 err

      err = -999.d0

      if (n.ne.4) stop 'oned:  input vector not equal to 4'

      a = line(1)
      b = line(2)
      c = line(3)
      d = line(4)
      out = err

      if(x.eq.0.)out = b
      if(x.eq.1.)out = c
      if((b.eq.err).or.(c.eq.err))return
      if((a.eq.err).or.(d.eq.err)) go to 2
      out = (1.0-x)*(b+x*(0.5*(c-a)+x*(0.5*(c+a)-b)))+x*(c+(1.0-x)*(0.5
     .     *(b-d)+(1.0-x)*(0.5*(b+d)-c)))
      return
 2    out = b*(1.0-x)+c*x
      if(a.ne.err)out = b+x*(0.5*(c-a)+x*(0.5*(c+a)-b))
      if(d.ne.err)out = c+(1.0-x)*(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c))
      return
      end




