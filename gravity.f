!/**----------------------------------------------------------------------    
! @name       gravity.f
!
! Compute the acceleration due to the Earth's gravity at any latitude/altitude
!
! @author     Bill Schreiner
! @since      5/95
! @version    $Revision: 1.2 $
! @cvsinfo    $Id: gravity.f,v 1.2 2000/01/05 22:26:38 huntd Exp $
* -----------------------------------------------------------------------*/

      subroutine gravity(xlat,alt,galt)
!/**----------------------------------------------------------------------    
! @name       gravity
!
! This subroutine computes the Earth's gravity at any altitude
! and latitude.  The model assumes the Earth is an oblate 
! spheriod rotating at a the Earth's spin rate.  The model
! was taken from "Geophysical Geodesy, Kurt Lambeck, 1988".
!
! @parameter  
! @ input:    xlat, latitude in radians
! @           alt,  altitude above the reference ellipsiod, km
! @ output:   galt, gravity at the given lat and alt, cm/sec
! -----------------------------------------------------------------------*/

      implicit double precision (a-h,o-z)
c
      xmu = 398600.4415d0	! km^3/s^2
      ae = 6378.1363d0		! km
      f = 1.0d0/298.2564d0	!
      w = 7.292115d-05		! rad/s
      xm = 0.003468d0		!
c     f2 = -f + 5.0/2.0*xm - 17.0/14.0*f*xm + 15.0/4.0*xm**2
c     f4 = -f**2/2.0 + 5.0/2.0*f*xm
      f2 = 5.3481622134089D-03    
      f4 = 2.3448248012911D-05
c
c compute gravity at the equator, km/s2
      ge = xmu/ae**2/(1.0d0 - f + 1.5d0*xm - 15.0/14.0*xm*f)
c
c compute gravity at any latitude, km/s2
      g = ge*(1.0d0 + f2*(sin(xlat))**2 - 1.0/4.0*f4*(sin(2.0*xlat))**2)
c
c compute gravity at any latitude and at any height, km/s2
      galt = g - 2.0*ge*alt/ae*(1.0 + f + xm + (-3.0*f + 5.0/2.0*xm)*
     &                       (sin(xlat))**2) + 3.0*ge*alt**2/ae**2
c
      galt = galt*1.0d5		! convert from km/s2 to cm/s2
c
      return
      end
