	program trisrf
c
c Produces a mesh representation of the solvent accessible
c surface of a molecule, 
c outputs *.usr format files for display in INSIGHT
c takes a pdb file, assigns radii, 
c assigns an atomic density function, based on overlapping Gaussians
c then contours it
c to use less/memory/time, reduce parameter ngrid from 65
c to some smaller ODD number in trisrf.h
c kas aug 2003- investigate neighbors for surface mapping
c
	include "trisrf.h"
c------------------------------------------------------------------
	character*5 probe
	character*6 head,radstr
	character*40 filnam
	character*80 line
	character*24 crdstr
      character*60 ntitle
      character*10 phead


c------------------------------------------------------------------
	real*4 cmin(3),cmax(3),cran(3),xo(3),perfil
	real*4 radprb,rad,temp,bran,rmidg,gran
	real*4 clevel
c------------------------------------------------------------------
	integer iargc,nrec
	integer midg,nsurf,ineigh
	integer i,j,k
c------------------------------------------------------------------
c
c parse command line: 
c must be a pdb file name- 
c
	radprb = 1.8
	clevel = 5.
	if(iargc().lt.1)then
	  print * ,'Usage: trisrf pdb_file [coutour_level]'
	else if(iargc().eq.2)then
	  call getarg(2,probe)
	  read(probe,'(f5.2)',err=904)clevel
	endif
	call getarg(1,filnam)
	open(unit=10,file=filnam,err=901,status='old')
c	perfil = 25.0
c	perfil = 45.0
	perfil = 75.0
	print *,'% fill? '
	read(5,*)perfil
	write(6,*)filnam
	write(6,*)'grid size                  :',ngrid
	write(6,*)'percent of box to be filled:',perfil
	write(6,*)'probe radius (A)           :',radprb
	write(6,*)'contour level              :',clevel
c
c read atom coordinates through, assign radii and compute scale
c
	write(6,*)'   '
	write(6,*)'atomic coordinates being read from file'
	write(6,*)filnam
	write(6,*)'...'
	nat = 0
	nrec = 0
	do k = 1,3
	  cmin(k) = 1.e6
	  cmax(k) = -1.e6
	end do
	open(unit=19,file='trisrf.pdb',status='unknown')
	open(17,file='trisrf.tri')
	write(17,'(A)')'PDB_RECORD'
	open(20,file='trisrf.rec')
	write(20,'(A)')'VERTEX_PDB_RECORD'
102   continue
        read(10,204,end=302)line
204     format(A80)
        head = line(1:6)
	  call up(head,6)
c
c skip header lines
c
        if((head.ne.'ATOM  ').and.(head.ne.'HETATM')) go to 102
	  nrec = nrec + 1
	  atm = line(12:16)
	  call elb(atm,5)
	  call up(atm,5)
	  if(atm(1:1).eq.'C')then
	    rad = 1.90
	  else if(atm(1:1).eq.'O')then
	    rad = 1.60
	  else if(atm(1:1).eq.'N')then
	    rad = 1.65
	  else if(atm(1:1).eq.'H')then
	    rad = 0.0 
	  else if(atm(1:1).eq.'1')then
	    rad = 0.0
	  else if(atm(1:1).eq.'2')then
	    rad = 0.0
	  else if(atm(1:1).eq.'3')then
	    rad = 0.0
	  else if(atm(1:1).eq.'4')then
	    rad = 0.0
	  else if(atm(1:1).eq.'P')then
	    rad = 1.90
	  else if(atm(1:1).eq.'S')then
	    rad = 1.90
	  else
	    print *,'unknown atom radius -> 0: ',atm
	    rad = 0.
	  end if
	if(rad.gt.0.) then
	  nat = nat + 1
	  if(nat.gt.natmx) then
	    write(6,*)' maximum # of atom records exceeded'
	    write(6,*)nat
	    write(6,*)' - increase natmx'
	    stop 
	  end if
	  atmrad(nat)=rad
	  crdstr = line(31:54)
        read(crdstr,205,err=905)xo
	  do k = 1,3
	    atmcrd(k,nat) = xo(k)
	  end do
c	  print *,xo
205     format(3F8.3)
c
c find largest extent of molecule including radii of atoms
c
	  do i = 1,3
	    temp = xo(i) - rad
	    cmin(i) = amin1(cmin(i),temp)
	    temp = xo(i) + rad
	    cmax(i) = amax1(cmax(i),temp)
	  end do
        write(radstr,'(f6.2)')rad
	  line(55:60) = radstr
	  write(19,204)line
	  write(17,204)line
	end if
      go to 102
302   continue		
c	end of file
	close(10)
	close(19)
	write(17,'(A)')'END PDB_RECORD'
	write(17,'(A)')'TRIANGLE_XYZ'
	if(nat.gt.0)then
	  print *,nrec,' Coordinate records read'
	  print *,nat,' Atoms with non-zero radius found'
	else
	  print *,'Error: no atom records read'
	  stop
	end if
c
c calculate scale to fill required percentage of box- allow for probe radius
c added to vdw radius 
c
	do i = 1,3
	  cran(i) = cmax(i) - cmin(i) + 2.*abs(radprb)
	  oldmid(i) = (cmax(i) + cmin(i))/2.
	end do
	bran = 0.0
	do i = 1,3
	  bran = amax1(bran,cran(i))
	end do
	midg = (ngrid+1)/2
	rmidg = midg
	gran = ngrid - 1
	scale = gran*perfil/(100.*bran)
	write(6,*)'scale   (grids/A): ',scale
c	print *,'enter new scale'
c	read(5,*)scale

	write(6,*)'  '
	write(6,*)'xmin,xmax     (A):',cmin(1),cmax(1)
	write(6,*)'ymin,ymax     (A):',cmin(2),cmax(2)
	write(6,*)'zmin,zmax     (A):',cmin(3),cmax(3)
	write(6,*)'x,y,z range   (A): ',cran
	write(6,*)'scale   (grids/A): ',scale
	write(6,*)'object centre (A): ',oldmid
	write(6,*)'  '
	call gauden(abs(radprb),clevel)
      open(21,file='trisrf.phi',form='unformatted')
      phead = 'gauden    '
      ntitle = 'from trisrf '
      write(21)' now starting phimap'
      write(21)phead,ntitle
      write(21)epsmp
      write(21) ' end of phi map '
      write(21)scale,oldmid
      close(21)

c
c analyse neighbors
c
	do i = 0,6
	  nneigh(i) = 0
	end do
	nsurf = 0
      do k = 2,ngrid-1
      do j = 2,ngrid-1
      do i = 2,ngrid-1
	  neps(i,j,k) = 0
	  ineigh = 0
	  if(epsmp(i+1,j,k).gt.clevel)ineigh = ineigh + 1
	  if(epsmp(i-1,j,k).gt.clevel)ineigh = ineigh + 1
	  if(epsmp(i,j+1,k).gt.clevel)ineigh = ineigh + 1
	  if(epsmp(i,j-1,k).gt.clevel)ineigh = ineigh + 1
	  if(epsmp(i,j,k+1).gt.clevel)ineigh = ineigh + 1
	  if(epsmp(i,j,k-1).gt.clevel)ineigh = ineigh + 1
	  if(epsmp(i,j,k).gt.clevel)then
	    if((ineigh.gt.0).and.(ineigh.lt.6))then
		neps(i,j,k) = 1
	      nneigh(ineigh) = nneigh(ineigh) + 1
	      nsurf = nsurf + 1
	    end if
	  end if
      end do
      end do
      end do
	print *,'# of surface points: ',nsurf
	print *,'with # of neighbors inside,freq: '
	do i = 0,6
	  print *,i,nneigh(i)
	end do
	print *,' '
	do i = 0,6
	  nneigh(i) = 0
	end do
	nsurf = 0
      do k = 2,ngrid-1
      do j = 2,ngrid-1
      do i = 2,ngrid-1
	  if(neps(i,j,k).eq.1)then
	    nsurf = nsurf + 1
          ineigh = neps(i+1,j,k)
     &        + neps(i,j+1,k)
     &        + neps(i,j,k+1)
     &        + neps(i-1,j,k)
     &        + neps(i,j-1,k)
     &        + neps(i,j,k-1)
	    nneigh(ineigh) = nneigh(ineigh) + 1
	  end if
      end do
      end do
      end do
	print *,'# of surface points: ',nsurf
	print *,'with # of surface neighbors,freq: '
	do i = 0,6
	  print *,i,nneigh(i)
	end do
	call putsrf()
	call boxes(clevel)
	write(17,'(A)')'END TRIANGLE_XYZ'
	close(17)
	write(20,'(A)')'END VERTEX_PDB_RECORD'
	close(20)

	stop
901   print *,'can`t find coordinate file ',filnam
	stop
904	print *,'Error reading probe radius of ',probe
	stop
905   print *,'Error reading the atom file'
	stop

	end
c------------------------------------------------------------------
      subroutine elb(txt,len)
c
c eliminate leading blanks from a character string
c
      character*(*) txt
      character*80 save

      do 9000 i=1,len
        if(txt(i:i).ne.' ') then
          nonb = i
          go to 100
        end if
9000  continue
      return
100   continue
      save = txt(nonb:len)
      txt = save
      return
      end
c------------------------------------------------------------
      subroutine up(txt,len)
c
c convert character string to upper case
c
      character*(*) txt
      character*80 save
      character*26 ualpha,lalpha
      data ualpha /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data lalpha /'abcdefghijklmnopqrstuvwxyz'/

      do 9000 i=1,len
        if((txt(i:i).ge.'a').and.(txt(i:i).le.'z')) then
          match = index(lalpha,txt(i:i))
          save(i:i) = ualpha(match:match)
        else
          save(i:i) = txt(i:i)
        end if
c      end do
9000	continue

      txt = save
      return
      end
c----------------------------------------------------------
	subroutine gauden(radprb,clevel)
c
c add gaussian density function to each grid point from each atom
c
c
c----------------------------------------------------------
	include "trisrf.h"
c----------------------------------------------------------
	real*4 radprb,rmidg,rad,rad2,d2,vin,surf
	real*4 xn(3),ismin(3),ismax(3)
	real*4 atdist(ngrid,ngrid,ngrid)
c----------------------------------------------------------
	integer i,j,k,ix,iy,iz
	integer nvin,nsurf,ineigh
	real*4 gaucnst,gauwidth,clevel
c----------------------------------------------------------
	data gaucnst / 10. /
c----------------------------------------------------------
	write(6,*)'   '
	write(6,*)'initializing SAV map...'
	write(6,*)'   '
	do k = 1,ngrid
	    do j = 1,ngrid
		do i = 1,ngrid
		  epsmp(i,j,k) = 0.
		  atmp(i,j,k) = 0.
		  atdist(i,j,k) = 1.e6
	end do
	end do
	end do
c
	rmidg = (ngrid+1)/2.
	do i = 1,nat
c
c scale atoms to grid space
c
	  do k = 1,3
	    xn(k) = (atmcrd(k,i)-oldmid(k))*scale + rmidg
	  end do
c
c scale radius to grid
c
	  gauwidth = (scale*atmrad(i))**2
	  rad = (atmrad(i)+radprb)*scale
	  rad2 = rad*rad
c
c find upper and lower grid index limits
c that fully enclose sphere of atom
c ensure they lie within grid
c
        do k = 1,3
          ismin(k) = (xn(k) - rad - 1.)
	    ismin(k) = min(ismin(k),ngrid)
	    ismin(k) = max(ismin(k),1)
          ismax(k) = (xn(k) + rad + 1.)
	    ismax(k) = min(ismax(k),ngrid)
	    ismax(k) = max(ismax(k),1)
	  end do
c
c calculate grid's distance from
c atom centre. if  less than atom rad + probe rad, add density in
c
        do iz =  ismin(3),ismax(3)
          do iy =  ismin(2),ismax(2)
            do ix =  ismin(1),ismax(1)
              d2 = (ix-xn(1))**2 + (iy-xn(2))**2 + (iz-xn(3))**2
              if(d2.le.rad2) then
                epsmp(ix,iy,iz) = epsmp(ix,iy,iz) + gaucnst*exp(-d2/gauwidth)
                epsmp(ix,iy,iz) = min(gaucnst,epsmp(ix,iy,iz))
		    if(d2.lt.atdist(ix,iy,iz))then
		      atdist(ix,iy,iz) = d2
			atmp(ix,iy,iz) = i
		    end if
              end if                  
	  end do
	  end do
	  end do
	end do
	call slice(clevel)
	nvin = 0
	nsurf = 0
	do k = 2,ngrid-1
	do j = 2,ngrid-1
	do i = 2,ngrid-1
	  if(epsmp(i,j,k).gt.clevel)nvin = nvin + 1
	  ineigh = 0
	  if(epsmp(i+1,j,k).gt.clevel)ineigh = ineigh + 1
	  if(epsmp(i-1,j,k).gt.clevel)ineigh = ineigh + 1
	  if(epsmp(i,j+1,k).gt.clevel)ineigh = ineigh + 1
	  if(epsmp(i,j-1,k).gt.clevel)ineigh = ineigh + 1
	  if(epsmp(i,j,k+1).gt.clevel)ineigh = ineigh + 1
	  if(epsmp(i,j,k-1).gt.clevel)ineigh = ineigh + 1
	  if((ineigh.gt.0).and.(ineigh.lt.6)) nsurf = nsurf + 1
	end do
	end do
	end do
	vin = nvin/(scale**3)
c sqrt 3. since mean orientation of surface elements of sphere 
c is at angle to grid
	surf = nsurf/(sqrt(3.)*scale**2)  
	print *,'# of points inside gaussian volume:',nvin
	print *,' giving volume (A^3) of: ',vin
	print *,'# of points on surface of gaussian volume:',nsurf
	print *,' giving surface area (A^2) of: ',surf
	print *,'using contour level of :',clevel
	return
	end


	subroutine slice(clevel)
	include "trisrf.h"
c-------------------------------------------------------------------------
c-------------------------------------------------------------------------
	integer irow(ngrid),i,j,nrow,k
	real*4 clevel
c-------------------------------------------------------------------------
	print *,'z-slice # (1-65)>> '
	read(5,*)nrow
	do i = 1,ngrid
	  do j = 1,ngrid
	    if(epsmp(i,j,nrow).gt.clevel)then
		irow(j) = 1
	    else
		irow(j) = 0
	    end if
	  end do
	  write(6,'(65I1)')(irow(j),j=1,ngrid)
	end do
	return
	end

	subroutine putsrf()
c put surface points out in insight *.usr format
c-------------------------------------------------------
	include "trisrf.h"
c-------------------------------------------------------
	real*4 xn(3),rmidg
c-------------------------------------------------------
	integer nsrf,icol,ix1,iy1,iz1
c-------------------------------------------------------
	write(6,*)' putting surface points out...'
	open(15,file='trisrf.usr')
	rmidg = (ngrid+1)/2.
	nsrf = 0
	icol = 120
	write(15,'(4A)')'DOTS'
	do iz1 = 2,ngrid-1
	  xn(3) = (iz1-rmidg)/scale + oldmid(3)
	  do iy1 = 2,ngrid-1
	    xn(2) = (iy1-rmidg)/scale + oldmid(2)
	    do ix1 = 2,ngrid-1
	      xn(1) = (ix1-rmidg)/scale + oldmid(1)
c		nsrf = nsrf + neps(ix1,iy1,iz1)
		if(neps(ix1,iy1,iz1).eq.1)then
c
c a point on boundary-write it
c
		  nsrf = nsrf + 1
	        write(15,'(1x,3f9.3,i5)')xn,icol
		end if	
	    end do
	  end do
	end do
	close(15)
	print *,nsrf,' surface points written to trisrf.srf'
	return
	end

	subroutine boxes(clevel)
c analyse and contour boxes
c-------------------------------------------------------
	include "trisrf.h"
c-------------------------------------------------------
	real*4 clevel
	integer ecuts(3,ngrid,ngrid,ngrid),ecutbx(3,2,2,2)
	integer ibase(3)
	integer i,j,k,itype(0:8),ihigh,l
	integer icut(0:12)
c-------------------------------------------------------
	integer nsrf,icol,ix1,iy1,iz1,ncut
	integer i1,i2,i3,i4,i5,i6,i7,i8
	integer idist(3),j1,k1,l1,itemp,iwant
	integer jloop(4),nbox
c-------------------------------------------------------
	open(16,file='triline.usr')
	open(18,file='trilinel.dat')
	write(16,'(a4)')'LINE'
c
c generate cut points- 0 no cut 1 cut with high point toward higher index
c -1 or 2 high point to lower index
c
cc debug test all box types
c	nbox = 0
c	do i1 = 0,1
c	 epsmp(2,2,2) = i1
c	do i2 = 0,1
c	 epsmp(3,2,2) = i2
c	do i3 = 0,1
c	 epsmp(2,3,2) = i3
c	do i4 = 0,1
c	 epsmp(3,3,2) = i4
c	do i5 = 0,1
c	 epsmp(2,2,3) = i5
c	do i6 = 0,1
c	 epsmp(3,2,3) = i6
c	do i7 = 0,1
c	 epsmp(2,3,3) = i7
c	do i8 = 0,1
c	 epsmp(3,3,3) = i8
c
c	nbox = nbox + 1
	do i = 0,8
	  itype(i) = 0
	end do
c
c initialize triangle/point lists
c
	ncut = 0
	ntri = 0
	nlink = 0
	do i = 1,ngrid**3
	  trilink(i) = 0
	  do j = 1,3
	    trinext(j,i) = 0
	    do k = 1,4
	      trilist(k,j,i) = 0
	    end do
	  end do
	end do
	do i = 1,ngrid
	do j = 1,ngrid
	do k = 1,ngrid
	  do l = 1,3
	    ecuts(l,i,j,k) = 0
	    ecutlink(l,i,j,k) = 0
	  end do
	end do
	end do
	end do
cc debug test
c	do i = 2,3
c	do j = 2,3
c	do k = 2,3
	do i = 1,ngrid-1
	do j = 1,ngrid-1
	do k = 1,ngrid-1
	  ihigh = 0

	  if(epsmp(i,j,k).gt.clevel)ihigh = ihigh + 1 
	  if(epsmp(i+1,j,k).gt.clevel)ihigh = ihigh + 1 
	  if(epsmp(i,j+1,k).gt.clevel)ihigh = ihigh + 1 
	  if(epsmp(i+1,j+1,k).gt.clevel)ihigh = ihigh + 1
	  if(epsmp(i,j,k+1).gt.clevel)ihigh = ihigh + 1 
	  if(epsmp(i+1,j,k+1).gt.clevel)ihigh = ihigh + 1 
	  if(epsmp(i,j+1,k+1).gt.clevel)ihigh = ihigh + 1 
	  if(epsmp(i+1,j+1,k+1).gt.clevel)ihigh = ihigh + 1
	  itype(ihigh) = itype(ihigh) + 1
	  if((epsmp(i,j,k).lt.clevel).and.(epsmp(i+1,j,k).gt.clevel))then
	    ecuts(1,i,j,k) = 1
	    ncut = ncut + 1
	  else if((epsmp(i,j,k).gt.clevel).and.(epsmp(i+1,j,k).lt.clevel))then
	    ecuts(1,i,j,k) = -1
	    ncut = ncut + 1
	  else
	    ecuts(1,i,j,k) = 0
	  end if
	  if((epsmp(i,j,k).lt.clevel).and.(epsmp(i,j+1,k).gt.clevel))then
	    ecuts(2,i,j,k) = 1
	    ncut = ncut + 1
	  else if((epsmp(i,j,k).gt.clevel).and.(epsmp(i,j+1,k).lt.clevel))then
	    ecuts(2,i,j,k) = -1
	    ncut = ncut + 1
	  else
	  end if
	  if((epsmp(i,j,k).lt.clevel).and.(epsmp(i,j,k+1).gt.clevel))then
	    ecuts(3,i,j,k) = 1
	    ncut = ncut + 1
	  else if((epsmp(i,j,k).gt.clevel).and.(epsmp(i,j,k+1).lt.clevel))then
	    ecuts(3,i,j,k) = -1
	    ncut = ncut + 1
	  else
	  end if
	end do
	end do
	end do
	print *,'# of edge cuts found: ',ncut
	print *,'# of boxes with high points '
	nsrf = 0
	do i = 0,8
	  print *,i,itype(i)
	  nsrf = nsrf + itype(i)
	end do
	nsrf = nsrf - itype(0) - itype(8)
	print *,'number of surface boxes: ',nsrf
	do i = 0,12
	  icut(i) = 0
	end do
cc debug
c	do i = 2,2
c	do j = 2,2
c	do k = 2,2
	do i = 1,ngrid-2
	do j = 1,ngrid-2
	do k = 1,ngrid-2
c debug
c	  print *,'ecuts: ',(ecuts(1,i,j,k)) ,(ecuts(2,i,j,k)),(ecuts(3,i,j,k))
c        print *,          (ecuts(2,i+1,j,k)) , (ecuts(3,i+1,j,k)) 
c        print *,          (ecuts(1,i,j+1,k)) , (ecuts(3,i,j+1,k)) 
c        print *,          (ecuts(1,i,j,k+1)) , (ecuts(2,i,j,k+1)) 
c        print *,      (ecuts(3,i+1,j+1,k))  , (ecuts(2,i+1,j,k+1)) , (ecuts(1,i,j+1,k+1)) 
	  ncut = abs(ecuts(1,i,j,k)) + abs(ecuts(2,i,j,k)) + abs(ecuts(3,i,j,k)) +
     &          abs(ecuts(2,i+1,j,k)) + abs(ecuts(3,i+1,j,k)) +
     &          abs(ecuts(1,i,j+1,k)) + abs(ecuts(3,i,j+1,k)) +
     &          abs(ecuts(1,i,j,k+1)) + abs(ecuts(2,i,j,k+1)) +
     &          abs(ecuts(3,i+1,j+1,k))  + abs(ecuts(2,i+1,j,k+1)) + abs(ecuts(1,i,j+1,k+1)) 
	  if(ncut.ne.0)write(10,*)i,j,k,ncut
	  icut(ncut) = icut(ncut) + 1
	  if((ncut.gt.0).and.(ncut.le.12))then
	    ecutbx(1,1,1,1) = ecuts(1,i,j,k) 
	    ecutbx(2,1,1,1) = ecuts(2,i,j,k) 
	    ecutbx(3,1,1,1) = ecuts(3,i,j,k) 
	    ecutbx(2,2,1,1) = ecuts(2,i+1,j,k) 
	    ecutbx(3,2,1,1) = ecuts(3,i+1,j,k) 
	    ecutbx(1,1,2,1) = ecuts(1,i,j+1,k) 
	    ecutbx(3,1,2,1) = ecuts(3,i,j+1,k) 
	    ecutbx(1,1,1,2) = ecuts(1,i,j,k+1) 
	    ecutbx(2,1,1,2) = ecuts(2,i,j,k+1) 
	    ecutbx(3,2,2,1) = ecuts(3,i+1,j+1,k)  
	    ecutbx(2,2,1,2) = ecuts(2,i+1,j,k+1) 
	    ecutbx(1,1,2,2) = ecuts(1,i,j+1,k+1) 
c	print *,'ecutbx in boxes: ',ecutbx
	    ibase(1) = i
	    ibase(2) = j
	    ibase(3) = k
	    call wrapbox(ecutbx,ibase,jloop,clevel)
	  end if
	end do
	end do
	end do
	nsrf = 0
	print *,'# of boxes with cuts'
      do i = 0,12
        print *,i,icut(i)
        nsrf = nsrf + icut(i)
      end do
      nsrf = nsrf - icut(0) 
      print *,'number of surface boxes: ',nsrf
	print *,'number of triangles: ',ntri
	print *,'neighbors: '
	open(19,file='trinext.dat')
	do i = 1,ntri
	  write(19,'(4i8)')i,(trinext(k,i),k=1,3)
	end do
	close(16)
	close(18)
	close(19)
	return
	end 

	subroutine wrapbox(ecutbx,ibase,jloop,clevel)
	implicit none
	real*4 clevel
	integer ecutbx(3,2,2,2),ibase(3),jloop(4),kloop
	integer i,j,k,i1,j1,k1
c-------------------------------------------------------
	integer number(3,2,2,2) ! maps index #'s of edge cuts to order number
	integer idx(4,12) ! maps edge order number to index #'s of edge
	integer next(3,2,12) ! give edge number, direction
	! from inside box with high point to right, gives edge number of
	! next edge to left, ahead, to right
	integer ecutloc(12),istart,ifound
	integer indx(4),idir,inext,icurr,nloop,nleft
	integer itriangle(3),it
c-------------------------------------------------------
	data number / 
     & 1,2,3,   ! 000
     & 0,4,7,   ! 100
     & 8,0,5,   ! 010
     & 0,0,10,  ! 110
     & 6,9,0,   ! 001
     & 0,12,0,  ! 101
     & 11,0,0,  ! 011
     & 0,0,0    ! 111
     &            /
c-------------------------------------------------------
	data idx / 
     & 1,1,1,1,
     & 2,1,1,1,
     & 3,1,1,1,
     & 2,2,1,1,
     & 3,1,2,1,
     & 1,1,1,2,
     & 3,2,1,1,
     & 1,1,2,1,
     & 2,1,1,2,
     & 3,2,2,1,
     & 1,1,2,2,
     & 2,2,1,2
     &            /
c-------------------------------------------------------
	data next / 
     & 2,8,4,7,6,3,
     & 3,9,5,8,4,1,
     & 1,7,6,9,5,2,
     & 1,2,8,10,12,7,
     & 2,3,9,11,10,8,
     & 3,1,7,12,11,9,
     & 4,10,12,6,3,1,
     & 5,11,10,4,1,2,
     & 6,12,11,5,2,3,
     & 8,5,11,12,7,4,
     & 9,6,12,10,8,5,
     & 7,4,10,11,9,6
     &            /
c-------------------------------------------------------
c	print *,'in wrap box: '
c	print *,'ibase: ',ibase
c	print *,'ecutbx: ',ecutbx
	nleft = 0
	do i = 1,4
	  jloop(i) = 0
	end do
	kloop = 0
	do i = 1,12
	  ecutloc(i) = ecutbx(idx(1,i),idx(2,i),idx(3,i),idx(4,i))
	  if(ecutloc(i).ne.0)nleft = nleft + 1
	end do
c	call boxpic(ecutloc)
c
c find 1st cut point
c
300	continue
	istart = 0
	ifound = 0
	do while((ifound.eq.0).and.(istart.lt.12))
	  istart = istart + 1
	  ifound = ecutloc(istart)
	end do
c	print *,'first cut at point',istart
c	print *,'starting chase...'
	if(ifound.eq.-1)then
	  idir = 2
	else if(ifound.eq.1)then
	  idir = 1
	else
	  print *,'wrong cut type: ',istart,ifound
	  return
	end if
	icurr = istart
c	ecutloc(icurr) = 0
	nloop = 0
	it = 1
	itriangle(it) = icurr
c	print *,'istart/current: ',icurr
100	continue
c look at next 3 edges left, ahead, right
c always turn towards high point
	  do k1 = 3,1,-1
	    inext = next(k1,idir,icurr)
c	    print *,'k1,idir,inext: ',k1,idir,inext
c	    print *,'idx: ',idx(1,inext),idx(2,inext),idx(3,inext),idx(4,inext)
c	    print *, ecutloc(inext)
	    ifound = ecutloc(inext)
	    if(ifound.ne.0)then
	      icurr = inext
		ecutloc(icurr) = 0
		nleft = nleft - 1
c	      call boxpic(ecutloc)
		nloop = nloop + 1
	      if(icurr.eq.istart)then
c		  print *,'closing loop'
		  goto 200
		else
		  it = it + 1
		  itriangle(it) = icurr
		  if(it.eq.3)then
c		    print *,'new triangle: ',itriangle
		    call puttri(itriangle,ibase,idx,clevel)
		    itriangle(2) = itriangle(3)
		    it = 2
		  end if
		  idir = 1
	        if(ifound.eq.-1)idir = 2
c	        print *,'current edge,direction: ',icurr,idir
		  goto 100
		endif
	    endif
	  end do
200	continue
	kloop = kloop + 1
	jloop(kloop) = nloop
c	print *,'done loop of ',nloop
c	print *,'cut points left: ',nleft
	if(nleft.gt.0)goto 300
	return
	end

	subroutine boxpic(ecutloc)
	implicit none
	integer ecutloc(12)
c picture of box, top
	write(6,'(''    '',''*.....'',i2,''.....*'')')ecutloc(11)
	write(6,'(''   /|           /|'')')
	write(6,'('' '',I2,'' |         '',i2,'' |'')')ecutloc(9),ecutloc(12)
	write(6,'('' /  |         /  |  '')')
	write(6,'(''*.....'',i2,''.....*   |'')')ecutloc(6)
	write(6,'(''|   |        |   |  '')')
	write(6,'(''|  '',I2,''        |  '',i2)')ecutloc(5),ecutloc(10)
	write(6,'(''|   |        |   |  '')')
c	write(6,'(''|   |        |   |  '')')
	write(6,'(I2,''  |        '',i2,''  |'')')ecutloc(3),ecutloc(7)
	write(6,'(''|   |        |   |  '')')
c	write(6,'(''|   |        |   |  '')')
c bottom
	write(6,'(''|   '',''*.....'',i2,''.|...*'')')ecutloc(8)
	write(6,'(''|  /         |  /'')')
	write(6,'(''|'',I2,''          |'',i2)')ecutloc(2),ecutloc(4)
	write(6,'(''|/           |/  '')')
	write(6,'(''*.....'',i2,''.....*'')')ecutloc(1)
	return
	end

	subroutine puttri(itriangle,ibase,idx,clevel)
	include "trisrf.h"
c-----------------------------------
	integer ibase(3)
	integer idx(4,12) ! maps edge order number to index #'s of edge
	integer itriangle(3),it
c-----------------------------------
	character*1 point,line
c-----------------------------------
	integer i,j,k,l
      integer ix1,iy1,iz1,ilink,ilink1
	integer jlink,imatch,itri,jtri
c-----------------------------------
	real*4 xn(3,3),rmidg
	real*4 clevel,eps1,eps2,dx,ratm
	integer iatm(3)
c-----------------------------------
	data point / 'P' /
	data line / 'L' /
c
c generate coords of triangle from indices
c
c-------------------------------------------------------
	ntri = ntri + 1
      rmidg = (ngrid+1)/2.
c	print *,'itriangle in pout: ',itriangle
cc debug
c	scale = 1.
c
c gather real indices of cut point
c
	do it = 1,3
	  l = idx(1,itriangle(it))
	  i = idx(2,itriangle(it)) - 1 + ibase(1)
	  j = idx(3,itriangle(it)) - 1 + ibase(2)
	  k = idx(4,itriangle(it)) - 1 + ibase(3)
c	print *,'indices l,i,j,k: ',l,i,j,k
c
c generate real space coords of triangle
c
        xn(1,it) = (i-rmidg)/scale + oldmid(1)
        xn(2,it) = (j-rmidg)/scale + oldmid(2)
        xn(3,it) = (k-rmidg)/scale + oldmid(3)
	  eps1 = epsmp(i,j,k)
	  if(l.eq.1)then
	    eps2 = epsmp(i+1,j,k)
	    dx = (clevel - eps1)/(eps2-eps1)
	    dx = min(dx,0.95)
	    dx = max(dx,0.05)
	    if(eps2.gt.clevel)then
	      iatm(it) = atmp(i+1,j,k)
	    else
	      iatm(it) = atmp(i,j,k)
	    end if
	    xn(1,it) = xn(1,it) + dx/scale
	  else if(l.eq.2)then
	    eps2 = epsmp(i,j+1,k)
	    dx = (clevel - eps1)/(eps2-eps1)
	    dx = min(dx,0.95)
	    dx = max(dx,0.05)
	    if(eps2.gt.clevel)then
	      iatm(it) = atmp(i,j+1,k)
	    else
	      iatm(it) = atmp(i,j,k)
	    end if
	    xn(2,it) = xn(2,it) + dx/scale
	  else if(l.eq.3)then
	    eps2 = epsmp(i,j,k+1)
	    dx = (clevel - eps1)/(eps2-eps1)
	    dx = min(dx,0.95)
	    dx = max(dx,0.05)
	    if(eps2.gt.clevel)then
	      iatm(it) = atmp(i,j,k+1)
	    else
	      iatm(it) = atmp(i,j,k)
	    end if
	    xn(3,it) = xn(3,it) + dx/scale
	  end if
c debug
c	  if(ntri.eq.3632)then
c	print *,'itriangle in pout: ',itriangle
c	print *,'indices l,i,j,k: ',l,i,j,k
c	print *,'eps1,eps2: ',eps1,eps2
c	  end if
	  if((dx.lt.0.01).or.(dx.gt.0.99))then
	    print *,'point ',it,' of triangle ',ntri,' is awfully close to a grid point ',dx
	    
	  end if
c
c enter triangle in triangle lists and try to find neighbors
c
	  trilist(1,it,ntri) = l
	  trilist(2,it,ntri) = i
	  trilist(3,it,ntri) = j
	  trilist(4,it,ntri) = k
c enter triangle in this cut point's linked list
	  nlink = nlink + 1
	  if(ecutlink(l,i,j,k).eq.0)then
	    ecutlink(l,i,j,k) = nlink
	    cuttri(nlink) = ntri
	    trilink(nlink) = 0
	  else
	    ilink = ecutlink(l,i,j,k)
	    do while(ilink.ne.0)
		ilink1 = ilink
	      ilink = trilink(ilink1)
	    end do
	    cuttri(nlink) = ntri
	    trilink(ilink1) = nlink
	  end if
	end do
c
c to find neighbors of this triangle, examine each pair of cutpoints
c list of triangles to find a common triangle, note that list of triangles in ascending order
c
c	do it = 1,3
c	  do i = 1,4
c	    print *,'trilist: ',trilist(i,it,ntri)
c	  end do
c	end do
	do i = 1,3
	  ix1 = i
	  iy1 = i + 1
	  if(iy1.eq.4)iy1 = 1

	ilink = ecutlink(trilist(1,ix1,ntri),trilist(2,ix1,ntri),trilist(3,ix1,ntri),trilist(4,ix1,ntri))
	jlink = ecutlink(trilist(1,iy1,ntri),trilist(2,iy1,ntri),trilist(3,iy1,ntri),trilist(4,iy1,ntri))
c	print *,'ilink, jlink: ',ilink,jlink
	itri = 0
	jtri = 0
	imatch = 0
	do while((ilink.ne.0).and.(itri.lt.ntri))
	  itri = cuttri(ilink)
	  ilink = trilink(ilink)
	  do while((jlink.ne.0).and.(jtri.lt.itri))
	    jtri = cuttri(jlink)
	    jlink = trilink(jlink)
	  end do
	  if((itri.ne.ntri).and.(itri.eq.jtri))then
	    imatch = jtri
	    itri = ntri
	  end if
	end do
	if(imatch.ne.0)then
c	  print *,'have a match for points 1, 2 of ',ntri,imatch
	  if(trinext(1,ntri).eq.0)then
	    trinext(1,ntri) = jtri
	  else if(trinext(2,ntri).eq.0)then
	    trinext(2,ntri) = jtri
	  else if(trinext(3,ntri).eq.0)then
	    trinext(3,ntri) = jtri
	  else
	    print * ,'more than 3 neighbors! for ',ntri
	  end if
	  if(trinext(1,jtri).eq.0)then
	    trinext(1,jtri) = ntri
	  else if(trinext(2,jtri).eq.0)then
	    trinext(2,jtri) = ntri
	  else if(trinext(3,jtri).eq.0)then
	    trinext(3,jtri) = ntri
	  else
	    print * ,'more than 3 neighbors! for ',jtri
	  end if
	end if
	end do
c
c write to insight type plot file
c
	write(16,'(1x,3f9.3,1x,A1)')(xn(k,3),k=1,3),point
	do it = 1,3
	  write(16,'(1x,3f9.3,1x,A1)')(xn(k,it),k=1,3),line
c	  write(17,'(1x,3f9.3,1x,A1)')(xn(k,it),k=1,3)
	  ratm = iatm(it)
	  write(17,'(1x,3f9.3,f15.6)')(xn(k,it),k=1,3),ratm
	  write(20,'(2i8)')ntri,iatm(it)
	end do
	write(18,'(1x,3f9.3,1x,A1)')(xn(k,1),k=1,3),point
	write(18,'(1x,3f9.3,1x,A1)')(xn(k,2),k=1,3),point
	write(18,'(1x,3f9.3,1x,A1)')(xn(k,2),k=1,3),point
	write(18,'(1x,3f9.3,1x,A1)')(xn(k,3),k=1,3),point
	write(18,'(1x,3f9.3,1x,A1)')(xn(k,3),k=1,3),point
	write(18,'(1x,3f9.3,1x,A1)')(xn(k,1),k=1,3),point

	return
	end
