!/**----------------------------------------------------------------------    
! @name       coords.f
!
! Determine sun-fixed latitude and longitude
!
! Input is the time and location of a point on the Earth.
!
! @author     Doug Hunt
! @since      10/24/97
! @version    $Revision: 1.2 $
! @cvsinfo    $Id: coords.f,v 1.2 2000/01/05 22:26:36 huntd Exp $
* -----------------------------------------------------------------------*/


      program foo  ! needed to allow g77 linkage (supplies __MAIN symbol)
      end


      subroutine coords (date, geolat, geolon, sflat, sflon)
!/**----------------------------------------------------------------------    
! @name       coords
!
! Determine sun-fixed latitude and longitude
!
! Input is the time and location of a point on the Earth.
!
! Processing:      - Compute solar orbital elements 
!                  - Compute sun position and velocity vectors
!                  - Convert input lat/lon into ECI unit vector
!                  - Rotate this vector into sun-fixed frame
!                  - Convert the rotated vector back into (sun fixed) lat/lon
!
! @parameter  date/time [empress 14 digit format--double]
! @           geodetic_latlon [degrees--double(2)],
! @return     sunfixed_latlon [degrees--double(2)]
! -----------------------------------------------------------------------*/
**************************************************************************
c
c     uses j. kwok's subroutines to compute the lunar and solar
c     ephemerides for a given julian date.
c     
c     t. kelecy    4-13-87
c     d. hunt      6-13-95
c
      implicit none
      real*8 sflat, sflon, geolat, geolon, tout
      real*8 tin(2), es(7), ys(6), sunpos(3), date, sunvel(3), unit(3)
c
      data es /7*0.d0/ 
      data ys /6*0.0d0/
c
c     Convert input YYYYMMDDHHMMSS empress date/time to tin format
c     (YYYYMMDD. , HHMMSS.SSSS)
c
      tin(1) = int (date/1e6)
      tin(2) = date - (tin(1) * 1e6)

      ! /**
      !  * @call coords.f:julian
      !  compute julian date
      !  */
      call julian(tin,tout)

      ! /**
      !  * @call  coords.f:setsm
      ! compute the solar ephemerides
      !  */
      call setsm(tout,es)
      call ephem(tout,tout,es,ys)

      ! /**
      !  * @call  coords.f:xyz
      !    compute vector from center of earth to sun (sun position vector)
      !    and sun velocity vector.
      !  */
      call xyz(ys, sunpos, sunvel)

      ! /**
      ! @call coords.f:latlon2unit
      !  Compute the unit ECI position vector of the occultation
      !  given the julian time and the lat/lon
      ! */
      call latlon2unit (tout, geolat, geolon, unit)

      ! /**
      ! @call coords.f:sunfixed
      ! determine if the lat, lon are in sunlight
      ! */
      call sunfixed(sunpos, sunvel, unit, sflat, sflon)

      return
      end

c---------------------------------------------------------------------------

      subroutine julian(tin,tout)
!/**----------------------------------------------------------------------    
! @name       julian
!
! Computes julian date when given calander date and time.
!
! @author     John Kwok
! @reference  jpl em 312/86-151, 30 june 1986
! @parameter  tin   array of 2 values containing yyyymmdd. and hhmmss.ssss
! @return     tout  output julian date
! -----------------------------------------------------------------------*/
c
      implicit double precision (a-h,o-z)
      dimension tin(2)
      data d0,d1,d2,d3,d4,d5/2415020.5d0,1.d-4,1.d-2,24.d0,1.44d3
     1,8.64d4/,eps/1.0d-14/
c
c  this is the julian date for 19000101
c
      tout=d0
c
      iy=int(tin(1)*d1-1.9d3)
      im=int(tin(1)*d2-1.9d5)-iy*100
      id=int(tin(1)-1.9d7)-iy*10000-im*100
      ihour=int(tin(2)*d1+eps)
      imin=int(tin(2)*d2+eps)-ihour*100
      sec=tin(2)-ihour*10000-imin*100
      jd=iy*365+(iy-1)/4
      im1=im-1
      if (im1.eq.0) go to 12
      go to (1,2,3,4,5,6,7,8,9,10,11),im1
    1 jd=jd+31
      go to 12
    2 jd=jd+59
      go to 12
    3 jd=jd+90
      go to 12
    4 jd=jd+120
      go to 12
    5 jd=jd+151
      go to 12
    6 jd=jd+181
      go to 12
    7 jd=jd+212
      go to 12
    8 jd=jd+243
      go to 12
    9 jd=jd+273
      go to 12
   10 jd=jd+304
      go to 12
   11 jd=jd+334
   12 continue
      if (iy/4*4-iy.eq.0.and.im.gt.2) jd=jd+1
      jd=jd+id-1
      tout=tout+jd+ihour/d3+imin/d4+sec/d5
      return
      end

c---------------------------------------------------------------------------

      subroutine setsm(t,es)
!/**----------------------------------------------------------------------    
! @name       setsm
! 
! Sets up some orbital elements of the built-in luni-solar
! ephemerides
!
! @author     John Kwok
! @reference  jpl em 312/86-151, 30 june 1986
! @           explanatory supplement to the astronomical ephemeris and the
! @           american ephemeris and nautical almanac
! @parameter  t   current julian date
! @return     es  mean orbital elements of the sun in earth mean equator
! @               and equinox of date
! @           (1)  = a, semi-major axis (km)
! @           (2)  = e, eccentricity
! @           (3)  = i, inclination (rad)
! @           (4)  = capw, longitude of ascending node, =0 by definition
! @           (5)  = w, argument of periapsis (rad)
! -----------------------------------------------------------------------*/
      implicit double precision (a-h,o-z)
      dimension es(7)
      data rjd/2.41502d6/
      dt=t-rjd
      dt1=dt*1.d-4
      dt2=dt1*dt1
      dt3=dt1*dt2
      es(1)=1.496d8
      es(2)=1.675104d-2-1.1444d-5*dt1-9.4d-9*dt2
      es(3)=.4093197474d0-6.217910d-5*dt1-2.1468d-9*dt2+1.7977d-10*dt3
      es(4)=0.d0
      es(5)=4.908229653d0+8.2149855d-7*dt+5.9167d-7*dt2+1.22d-9*dt3
      return
      end

c---------------------------------------------------------------------------

      subroutine ephem(t,tr,es,ys)
!/**----------------------------------------------------------------------    
! @name       ephem
! Computes the orbital elements of a sun and/or a moon
!
! For mercury, venus, and mars, there is no luni-perturbation, user
! should input the sun's orbital elements relative to planet equator
! of epoch.  the inertial x-axis could be planet equinox or the
! direction of the prime meridian at some epoch such as defined by
! the iau.  for earth, the orbital elements of the moon varies too
! much to use two-body ephemeris.  user may use the mean lunar
! ephemeris built into the program.  the short period variations of
! the moon may introduce an error up to two degrees.  as for the sun,
! the mean orbital elements of the sun relative to earth equator and
! equinox are built in also.  but since they don't change as rapidly
! as those of the moon, they are initialize elsewhere and then held
! fixed for the duration of the trajectory propagation.
!
! @author     John Kwok
! @reference  jpl em 312/86-151, 30 june 1986
! @           explanatory supplement to the astronomical ephemeris and the
! @           american ephemeris and nautical almanac
! @parameter  isun   = see subroutine lop    
! @           isrp   = see subroutine lop    
! @           iephem = see subroutine lop    
! @           t      = current julian date   
! @           tr     = reference julian date 
! @           es     = see subroutine lop    
! 
! @return     ys     = six output orbital elements a, e, i, node, w, and m,
! @                    of sun at t (km, rad)
! -----------------------------------------------------------------------*/
      implicit double precision (a-h,o-z)
      dimension es(7),ys(6)
      data rjd/2.41502d6/
      data pi,tpi/3.141592653589793d0,6.283185307179586d0/

      dt=t-rjd
      dt1=dt*1.d-4
      dt2=dt1*dt1
      dt3=dt1*dt2
c
c  sun relative to earth in mean equator of date
c
      do 110 i=1,5
  110 ys(i)=es(i)
      x = 6.256583575d0+1.720196977d-2*dt - 1.9548d-7*dt2-1.22d-9*dt3
      ys(6)= x - (int (x/tpi) * tpi)  ! replacement for dmod function

      return
      end

c---------------------------------------------------------------------------

      subroutine xyz(ys, pos, vel)
!/**----------------------------------------------------------------------    
! @name       xyz
! 
! Computes the xyz position of the sun in ECI coordinates at input time
!
! @author     Doug Hunt
! @reference  GPS Theory and Practice, Hoffmann-Wellenhof et al. pp. 45-4
! @parameter  ys   six output orbital elements a, e, i, node, w, and m,
! @                of sun at t (km, rad)
! @return     xyzpos   the X, Y, and Z position components of the sun at t (km)
! -----------------------------------------------------------------------*/
      implicit none
      real*8 vel(3), pos(3), ys(6), mu, epsi, a, e, i, w, omg, ma, ea1
      real*8 ea, pi, f2, f, r, wf, x, y, z, xv1, xv2, xv, yv1, yv2
      real*8 yv, zv, mag, l

      pi = dacos(-1.d0)

      mu = 1.3271544d11 ! Heliocentric mu value, km**3/sec**2
                        ! Ref:  Bate, Mueller & White, Fundamentals of Astrodynamics
                        ! pg 429.
      epsi = 1e-8

      a = ys(1) ! Semi-major axis of orbit (km)
      e = ys(2) ! eccentricity (unitless)
      i = ys(3) ! inclination  (rad)
      w = ys(5) ! Argument of perigee (rad)
      omg = ys(4) ! Omega, right ascension of ascending node = 0 for sun (rad)
      ma  = ys(6) ! mean anomaly (rad)

c
      ea1 = ma
c
c PERFORM NEWTON-RAPHSON ITERATION TO DETERMINE ECCENTRIC ANOMALY
c

 10   ea = ea1
      if (ea .gt. 2*pi) then ! MODULATE ECCENTRIC ANOMALY IF GREATER THAN 2*PI
         ea = ea - 2*pi
      endif
      ea1 = ea - (ea - e*sin(ea) - ma)/(1 - e*cos(ea))
      if (abs(ea1-ea) .lt. epsi) goto 20
      goto 10
 20   continue

c
c COMPUTE TRUE ANOMALY
c
      f2 = atan(sqrt((1+e)/(1-e))*tan(ea/2))
      f  = 2*f2
c     
c COMPUTE RADIUS
c
      r = a*(1 - e*cos(ea))
c
c COMPUTE  SPECIFIC  ANGULAR MOMENTUM
c
      l = sqrt(mu*a*(1 - e**2))
c
c COMPUTE ARGUMENT OF LATITUDE
c
      wf = w + f
c
c COMPUTE POSITIONS (in km)
c
      x = r*(cos(omg)*cos(wf) - sin(omg)*sin(wf)*cos(i))
      y = r*(sin(omg)*cos(wf) + cos(omg)*sin(wf)*cos(i))
      z = r*(sin(i)*sin(wf))

c
c COMPUTE VELOCITIES (in km/sec)
c
      xv1 = x*l*e*sin(f)/(r*a*(1 - e**2))
      xv2 = (l/r)*(cos(omg)*sin(wf) + sin(omg)*cos(wf)*cos(i))
      xv = xv1 - xv2
      yv1 = y*l*e*sin(f)/(r*a*(1 - e**2))
      yv2 = (l/r)*(sin(omg)*sin(wf) - cos(omg)*cos(wf)*cos(i))
      yv = yv1 - yv2
      zv = z*l*e*sin(f)/(r*a*(1-e**2))+(l/r)*sin(i)*cos(wf)

c
c     Normalize the position vector
c

      mag = sqrt(x**2 + y**2 + z**2)
      pos(1) = x / mag
      pos(2) = y / mag
      pos(3) = z / mag

c
c     Normalize the velocity vector
c

      mag = sqrt(xv**2 + yv**2 + zv**2)
      vel(1) = xv / mag
      vel(2) = yv / mag
      vel(3) = zv / mag

      return
      end

c---------------------------------------------------------------------------

      subroutine latlon2unit (tjul, lat, lon, unit)
!/**----------------------------------------------------------------------    
! @name       latlon2unit
! 
! Compute a unit position vector in ECI coordinates given
! a julian time value and a geodetic lat/lon point.
!
! @author     Doug Hunt
! @parameter  
! @ input:    tjul:  Julian time value
! @           lat, lon:  geocentric latitude and longitude of point on earth in 
! @                      question (deg)
! @ output:   unit:  ECI unit position vector of input lat/lon at tjul.
! -----------------------------------------------------------------------*/

      implicit none
      real*8 lat, lon, pi, unit(3), latr, lonr
      real*8 julfrac, tjul, gangle, jsecs

      pi = dacos(-1.d0)

c
c     Convert lat, lon to a unit position vector in ECI coordinates
c

      julfrac = tjul - dint(tjul) ! fractional part of tjul (measure 
                                  ! of hours since noon GMT)
      jsecs   = (julfrac - 0.5) * 86400   ! Seconds into day
      
      
      call gast (tjul, jsecs, gangle)    ! Compute Greenwich Apparent Siderial
                                         ! Time angle, gangle
      if (lon.lt.0.d0) then
         lonr = lon + 360.d0
      else
         lonr = lon
      endif
      lonr = lonr * (pi / 180)     ! convert lon to rad, 0 to 2pi
      latr = lat * (pi / 180)      ! convert lat to radians
      
      lonr = lonr + gangle         ! convert to ECI by adding Greenwich angle

c
c     Compute an ECI unit vector pointing to the input lat/lon
c
      unit(1) = cos(latr) * cos(lonr)
      unit(2) = cos(latr) * sin(lonr)
      unit(3) = sin(latr)

      return
      end

c---------------------------------------------------------------------------


      subroutine sunfixed(xsun, vsun, unit, sflat, sflon)
!/**----------------------------------------------------------------------    
! @name       sunfixed
! 
! Compute sun fixed lat and lon (degrees) from sun XYZ pos and vel unit
! vectors.
!
! @author     Doug Hunt
! @parameter  
! @ input:    xsun:  XYZ unit position vector of sun (ECI coordinates) 
! @           vsun:  XYZ unit velocity vector of sun (ECI coordinates)
! @           unit:  unit ECI position vector of the occultation
! @ output:   sflat, sflon:  Sun fixed latitude and longitude (deg)
! -----------------------------------------------------------------------*/

      implicit none
      real*8 xsun(3), vsun(3), pi, unit(3)
      real*8 sflat, sflon, r(3,3), xs(3), zsun(3), ysun(3)

      pi = dacos(-1.d0)

c
c     Rotate this unit vector into the Sun-fixed frame of reference
c

c     compute cross product, xsun X vsun, to get zsun
      zsun(1) = xsun(2)*vsun(3) - xsun(3)*vsun(2)
      zsun(2) = xsun(3)*vsun(1) - xsun(1)*vsun(3)
      zsun(3) = xsun(1)*vsun(2) - xsun(2)*vsun(1)

c     compute cross product, zsun X xsun, to get ysun
      ysun(1) = zsun(2)*xsun(3) - zsun(3)*xsun(2)
      ysun(2) = zsun(3)*xsun(1) - zsun(1)*xsun(3)
      ysun(3) = zsun(1)*xsun(2) - zsun(2)*xsun(1)

c     rotation matrix
      r(1,1) = xsun(1)
      r(1,2) = xsun(2)
      r(1,3) = xsun(3)
      r(2,1) = ysun(1)
      r(2,2) = ysun(2)
      r(2,3) = ysun(3)
      r(3,1) = zsun(1)
      r(3,2) = zsun(2)
      r(3,3) = zsun(3)

c     x(3), occultation Cartesian vector in eci coords
c     xs(3), occultation Cartesian vector in sun coords

      xs(1) = r(1,1)*unit(1) + r(1,2)*unit(2) + r(1,3)*unit(3)
      xs(2) = r(2,1)*unit(1) + r(2,2)*unit(2) + r(2,3)*unit(3)
      xs(3) = r(3,1)*unit(1) + r(3,2)*unit(2) + r(3,3)*unit(3)
c
c     now rotate xyz sun vector to lat,lon of sun vector
c
      sflat = asin(xs(3)) * (180/pi)
      sflon = atan2 (xs(2), xs(1)) * (180/pi)

      return
      end

c---------------------------------------------------------------------------

      double precision function get_true_anomaly (mk, e)
!/**----------------------------------------------------------------------    
! @name       get_true_anomaly
! Compute the true anomaly from the mean anomaly and the eccentricity.
! use newton's method to solve the non-linear equation Mk = Ek - e*sin(Ek) for Ek.
! Use eccentric anomaly thus computed to find true anomaly.
!
! @author     Doug Hunt
! @parameter  mk   Mean Anomaly
! @           e    Eccentricity
! @return     true anomaly
! -----------------------------------------------------------------------*/
      implicit none

      integer max, j
      real*8  f, df, dx, tol, mk, e, v, ek

c Compute Ek

      max = 20      ! max iterations
      tol = 1.0d-10 ! solution tolerence

      ek = mk       ! first guess at solution
      
      do j=1,max
         f = mk - (ek - e * sin(ek))
         df= -1 + (e * cos(ek))
         dx = f/df
         ek = ek - dx
         
         if (abs(dx).lt.tol) then
            goto 100
         endif
      end do

 100  continue

c Compute v (true anomaly) from Ek and e

      v = 2.0 * atan ( sqrt( (1.0+e)/(1.0-e) ) * tan(ek/2.0) )
      get_true_anomaly = v
      return
      end
      

c---------------------------------------------------------------------------

      subroutine gast(djd, utco, theta)
!/**----------------------------------------------------------------------    
! @name       gast
!
! This subroutine computes the Greenwich Apparent Siderial
! Time angle given a UTC date and time.
!
! @author     Bill Schreiner, Doug Hunt
! @parameter  djd    Julian day/time (whole day, eg. 10/04/1995 = 2449994.5)
! @           utco   Seconds since start of day
! @return     theta, GAST angle in radians
! -----------------------------------------------------------------------*/
       implicit none

      real*8 djd, utco, theta
      real*8 tu, gmst, utc, pi

      pi = dacos(-1.0d0)      !pi

c
C Coordinate transform from the celestial inertial reference frame to the geo-
C centered Greenwich reference frame.

      tu = (djd - 2451545.0d0)/36525.0d0
      gmst = 24110.54841 + 8640184.812866*tu +
     1       0.093104*tu**2 - 6.2d-6*tu**3       !gmst=Greenwich mean...

      utc = utco * 1.0027379093
      gmst = gmst + utc         !in seconds, without eoe correction.

      ! gmst may be positive or negative.  Subtract or add until
      ! a value between 0 and 86400 is obtained.
      do while (gmst .lt. 0.0d0) 
        gmst = gmst + 86400.0d0
      end do
      do while (gmst .gt. 86400.0d0) 
        gmst = gmst + 86400.0d0
      end do

C gmst = the Greenwich mean sidereal time.
C This gmst is without the corrections from the equation of equinoxes.  For
C GPS/MET applications, the corrections from equation of equinoxes is not 
C necessary because of the accurary needed.

      theta = gmst*2.0d0*pi/86400.0d0 !*** This is the THETA in radians.

      return
      end
