!***********************************************************************************************
!
! Filename:     constants.F90
!
! Module name:
!
! Description:  Catalog of physical constants
!
! Module procedures:
!
! References:
!
!***************************************** Change log ******************************************
!
! Creator:              Hari Nair
! Creation date:        Dec 05, 2005
! Modification
!
!    Date:  02/11/08 	Developer: mcduffie
!    Description:       Add constant for gas absorption specific humidity correction (bugzilla 174)
!
!    Date:  mm/dd/yy 	Developer: username
!    Description:
!
!***********************************************************************************************
!
!              Copyright 2005, by the California Institute of Technology
!        ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
!      Any commercial use must be negotiated with the Office of Technology Transfer
!                       at the California Institute of Technology.
!
!        This software may be subject to U.S. export control laws and regulations.
!        By accepting this document, the user agrees to comply with all applicable
!                          U.S. export laws and regulations.
!    User has the responsibility to obtain export licenses, or other export authority
! as may be required before exporting such information to foreign countries or providing access
!                                to foreign persons.
!!END********************************************************************************************

module constants

  implicit none

  ! Astronomical unit, in km (http://neo.jpl.nasa.gov/glossary/au.html)
  double precision, parameter :: AU = 1.49597870691D8     
  ! Boltzmann's constant, in J K^{-1}      
  double precision, parameter :: BOLTZMANN = 1.3806505D-23
  ! Avogadro's number, molecules/mole
  double precision, parameter :: NA = 6.02214179D23
  ! PI?  Never heard of it.
  double precision, parameter :: PI = 3.14159265358979323846D0
  ! Planck's constant, in J s
  double precision, parameter :: PLANCK = 6.6260693D-34
  ! Earth equatorial radius, in km
  double precision, parameter :: RADIUS = 6378.137D0
  ! Gas constant, in J g^{-1} K^{-1}
  double precision, parameter :: RATM = 8314.D0/28.96D0
  ! Gravitational acceleration at the earth's surface, m/s^2
  double precision, parameter :: SGRAV = 9.80665D0
  ! Speed of light, m s^{-1}
  double precision, parameter :: SPEED_OF_LIGHT = 2.99792458D8
  ! In radians at 1 AU (Chollet and Sinceac, Astronomy and Astrophysics Supplement, v.139, p.219-229
  double precision, parameter :: SOLAR_ANGULAR_RADIUS = (959.44D0 / 3600) * (PI / 180)
  ! molar weight of dry air, in g
  double precision, parameter :: WGT = 28.9600D0
  ! specific humidity factor
  double precision, parameter :: WH2O = 18.01528D0
  
  ! Not physical constants, but used by the code in many places
  integer, parameter :: MAX_LINE_LENGTH = 2048
  integer, parameter :: MAX_LARGE_MATRIX_LOG_LIMIT = 100
  
end module constants
