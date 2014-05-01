     *****************************************************
     **  prepared April 5, 2002 by Robert Koelemeijer   **
     *****************************************************

This database is provided on an "as is" basis. If you plan to use this 
database, please ask Dr. R. Koelemeijer for permission. Distribution
of this database without contacting R. Koelemeijer is not permitted.


** Contact ** 
Dr. Robert Koelemeijer
Space Research Organization Netherlands
Sorbonnelaan 2
3584 CA Utrecht
The Netherlands
tel: +31 30 253 8573
fax: +31 30 254 0860
email: R.B.A.Koelemeijer@sron.nl


** outline of method **
The database contains minimum Lambert-equivalent reflectivities (MLER)
values at eleven 1-nm wide wavelength bins centered at 335.0, 380.0, 416.0, 
440.0, 463.0, 494.5, 555.0, 610.0, 670.0, 758.0, 772.0 nm.
Lambert-equivalent reflectivities (LERs) are derived from the GOME data 
using the Doubling-Adding KNMI (DAK) polarized radiative transfer code.
The model consists of an atmosphere for which Rayleigh scattering and 
ozone absorption is taken into account, and which is bounded below by a 
Lambertian surface. 
For each GOME pixel and wavelength, the LER is found as the value of the 
surface albedo in the radiative transfer model needed to match the observed 
top-of-atmosphere reflectivity. The GOME data considered were acquired in 
the period June 27, 1995  -- December 31, 2000.

The LER spectra are then binned per month and in grid-cells of 1*1 degree. 
The MLER is then determined as the minimum LER in each grid-cell and 
each month, over a period of 5.5 years, and preforming corrections for 
missing data, cloud contamination, etc (see FLAG-files).
The minimum is searched at 670 nm, and the corresponding values pertaining
to this GOME pixel are stored in the database, to ensure that the MLER values 
of different wavelengths are derived from the same GOME spectrum. 

In addition to the database with monthly minimum values, a database
was made with annual minimum values, containing the minimum MLER value that
occurred in each of the twelve months for each grid cell and wavelength.
In this database of annual minimum LER values it is no longer the case 
that the values are derived from the same GOME spectrum, as
minima at different wavelength generally occur in different months.
Both the monthly minimum and the annual minimum databases are accompanied 
by a database of flags, which keeps track of corrections that have been 
applied to each grid cell.


** Files **
sacspecTOTL01_335.dat -- monthly minimum LER values (MLER) for month 01
                         (January) and wavelength 335 nm
sacspecALLM335.dat    -- data for 335 nm, annual minimum values
sacspecFLAG01.dat     -- file with flags corresponding to the TOTL-files
sacspecFLAG335.dat    -- file with flags corresponding to the ALLM-files

The MLER values are three digit integers, and should be divided by 1000 to 
obtain the MLER in the range (0,1). The format resembles the well-known 
TOMS format.


** Fortran subroutine to read the files **
      subroutine rd_cscdb(filename,header,iler)
      implicit none
      integer i,j
      character*100 filename
      character*75 header(3)
      character*14 string
      integer iler(180,360)
      open(10,file=filename)
      do i=1,3
         read(10,'(a)') header(i)
      enddo
      do i=1,180
         read(10,91) (iler(i,j),j=1,360),string
      enddo
      close(10)
  91  format(14(25i3/),10i3,a14)
      return
      end


** Meaning of values in flag-files **
-----------------
FLAG	EXPLANATION
-----------------
0     OK, no corrections applied
1     Residual clouds above ocean detected => replaced by 
      weighted average of surrounding pixels in 5*5 degr area
2     Monthly variation threshold exceeded => replaced by
	maximum of adjacent months
3     Missing data, filled in by nearest month with data (in case of
	polar regions which are observed only part of the year), or with 
	average of adjacent pixels at same latitude (in case of ERS-2 tape 
	recorder overflow gap)
4     Missing data, filled in by surrounding pixels in 3*3 degr area
5     Missing data throughout the year. Data copied from other location
      with similar surface properties (polar areas)

VALUE+10 Residual cloud contamination likely in final data

REFERENCE 
Koelemeijer, R. B. A.,  J. F. de Haan, and P. Stammes, 
A database of spectral surface reflectivity in the range 335--772 nm
derived from 5.5 years of GOME observations, J. Geophys. Res., submitted.


