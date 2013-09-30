c---------------------------------------------------
c parameters for mesh
c---------------------------------------------------
	implicit none
	integer ngrid, natmx
	parameter (ngrid = 65)
	parameter (natmx = 40000)
c---------------------------------------------------
 	integer atmp(ngrid,ngrid,ngrid)
	integer nat
c
	real*4 oldmid(3),scale
	real*4 atmcrd(3,natmx),atmrad(natmx)
	real*4 epsmp(ngrid,ngrid,ngrid)
	integer neps(ngrid,ngrid,ngrid),nneigh(0:6)
	integer trilist(4,3,ngrid**3),trilink(ngrid**3),ntri,nlink
	integer trinext(3,ngrid**3),cuttri(ngrid**3)
	integer ecutlink(3,ngrid,ngrid,ngrid)
c---------------------------------------------------
      character*5 atm
c---------------------------------------------------
	common
     &	/maps/  epsmp,neps,trilist,trilink,ntri,ecutlink,trinext,nlink,atmp
     &	/scale/ scale,oldmid
     &      /atoms/ atmcrd,atmrad,nat
c---------------------------------------------------
