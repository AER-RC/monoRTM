!     path:      $HeadURL: https://svn.aer.com/svn/aer/project/RD/LBLRTM/trunk/src/phys_consts.f90 $
!     author:    $Author: malvarad $
!     revision:  $Revision: 16421 $
!     created:   $Date: 2012-10-19 09:20:37 -0400 (Fri, 19 Oct 2012) $
!
!  --------------------------------------------------------------------------
! |  Copyright �, Atmospheric and Environmental Research, Inc., 2012         |
! |                                                                          |
! |  All rights reserved. This source code is part of the LBLRTM software    |
! |  and is designed for scientific and research purposes. Atmospheric and   |
! |  Environmental Research, Inc. (AER) grants USER the right to download,   |
! |  install, use and copy this software for scientific and research         |
! |  purposes only. This software may be redistributed as long as this       |
! |  copyright notice is reproduced on any copy made and appropriate         |
! |  acknowledgment is given to AER. This software or any modified version   |
! |  of this software may not be incorporated into proprietary software or   |
! |  commercial software offered for sale.                                   |
! |                                                                          |
! |  This software is provided as is without any express or implied          |
! |  warranties.                                                             |
! |                       (http://www.rtweb.aer.com/)                        |
!  --------------------------------------------------------------------------
!
MODULE phys_consts   ! Physical constants

  implicit none
!                                                                       
!    Constants from NIST May 2010
!                                                                       
!     units are generally cgs                                           
!                                                                       
      real, parameter  :: PI = 3.1415926535898
      real, parameter  :: PLANCK = 6.62606876E-27    !Planck constant          | [g.cm^2/s]
      real, parameter  :: BOLTZ = 1.3806503E-16      !Boltzman Constant        | [g.cm^2/(s^2.K)]
      real, parameter  :: CLIGHT = 2.99792458E+10    !Speed of Light           | [cm/s]
      real, parameter  :: AVOGAD = 6.02214199E+23    !Avogadro's Number        | [molec./mole]
      real, parameter  :: ALOSMT = 2.6867775E+19 
      real, parameter  :: GASCON = 8.314472E+07      !Gas constant             | [g.cm^2/(s^2.K.mole)]
      real, parameter  :: RADCN1 = 1.191042722E-12   !1st radiation constant   |
      real, parameter  :: RADCN2 = 1.4387752         !2nd radiation constant   |
!                                                                       
!     The first and second radiation constants are taken from NIST.     
!     They were previously obtained from the relations:                 
!                            RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07    
!                            RADCN2 = PLANCK*CLIGHT/BOLTZ               

END MODULE phys_consts
