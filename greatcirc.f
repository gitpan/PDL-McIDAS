!/**----------------------------------------------------------------------    
! @name       greatcirc.f
!
! Compute great circle distance between two points on Earth
!
! @author     Doug Hunt, Bill Schreiner
! @version    $Revision: 1.2 $
! @cvsinfo    $Id: greatcirc.f,v 1.2 2000/01/05 22:26:39 huntd Exp $
* -----------------------------------------------------------------------*/

      subroutine greatcirc (lat1, lon1, lat2, lon2, dist)
!/**----------------------------------------------------------------------    
! @name       gravity
!
! Compute great circle distance between two points on Earth.
!
! Processing:      - Compute ECF positions of each lat/lon pair
!                  - Take the dot product of these two vectors to get the cosine of
!                    the angle between them
!                  - Multiply the resulting angle by the average earth radius at
!                    the two points.
!
! @parameter  
! @ input:    lat1, lon1, lat2, lon2 (degrees)
! @ output:   dist, Great circle distance between them, km.
! -----------------------------------------------------------------------*/
      implicit none
      real*8 lat1, lon1, lat2, lon2, dist, ecf1(3), ecf2(3), r1, r2
      real*8 unit1(3), unit2(3), costheta, ht

      integer itype

      ! Convert lat/lon pairs to ECF
      ht = 0.d0
      itype = 0 ! convert from lat/lon/ht to ecf
      call xyzsub (itype, ecf1(1), ecf1(2), ecf1(3), lat1, lon1, ht) 
      call xyzsub (itype, ecf2(1), ecf2(2), ecf2(3), lat2, lon2, ht) 

      ! Compute unit vectors
      r1 = sqrt (ecf1(1)**2 + ecf1(2)**2 + ecf1(3)**2)
      unit1(1) = ecf1(1)/r1
      unit1(2) = ecf1(2)/r1
      unit1(3) = ecf1(3)/r1

      r2 = sqrt (ecf2(1)**2 + ecf2(2)**2 + ecf2(3)**2)
      unit2(1) = ecf2(1)/r2
      unit2(2) = ecf2(2)/r2
      unit2(3) = ecf2(3)/r2


      ! Find the angle between these vectors
      costheta = unit1(1) * unit2(1) +
     +           unit1(2) * unit2(2) +
     +           unit1(3) * unit2(3)  

      ! compute distance using average earth radius of two positions
      dist = dacos(costheta) * ((r1 + r2)/2.d0)

      return
      end

      subroutine xyzsub(i,x,y,z,lat,lon,ht)
c This program converts (x,y,z) to (lat,lon,ht) and vice versa
c Equations obtained from ASEN5004 at CU.
c XYZ in km, lat/lon in degrees
c
      implicit real*8(a-h,o-z)
      double precision lat,lon,latr,lonr
c
c WGS-84 constants needed for conversion
      a = 6378.137d0
      fi = 298.257223563
      f = 1.d0/fi
      e2 = 2.d0*f - f*f
      pi = 4.d0*atan(1.d0)
      rad = pi/180.d0
c
      if (i.eq.1) goto 100
c
c......................................
c Convert from lat,lon,ht to x,y,z
c......................................
      latr = lat*rad
      lonr = lon*rad
c
      b = dsqrt(1.d0 - e2*(sin(latr))**2)
      x = (a/b + ht)*cos(latr)*cos(lonr)
      y = (a/b + ht)*cos(latr)*sin(lonr)
      z = (a*(1.d0 - e2)/b + ht)*sin(latr)
c
      goto 1000
c
c......................................
c Convert from  x,y,z  to lat,lon,ht
c......................................
c
c Iterative method described in ASEN5004
 100  lon = atan2(y,x)
c
c initial value of phi for latitude iteration
      n = 0
      dxy = dsqrt(x*x + y*y)
      phi0 = atan(z/(1.d0 - e2)/dxy)
      ht0 = 0.0d0
c
c iterate on the latitude
 10   continue
      n = n + 1
      if (n.gt.10) then
       print*,'more than 10 iterations'
       goto 1000
      endif
c
      xn = a/dsqrt(1.d0 - e2*(sin(phi0)**2))
      tphi = z/dxy*(1.d0 + e2*xn*sin(phi0)/z)
      phi = atan(tphi)
c
c  solve for the height
      xn = a/dsqrt(1.d0 - e2*(sin(phi)**2))
      ht = dxy/cos(phi) - xn
c
c  check latitude and height tolerence
      if (dabs(phi - phi0).gt.1.d-10.or.dabs(ht - ht0).gt.1.d-4) then
       phi0 = phi
       ht0 = ht
       goto 10
      endif
c
      lat = phi
c
c
 1000 continue
      return
      end
