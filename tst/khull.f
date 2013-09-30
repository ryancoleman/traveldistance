        program khull
c
c generate points/triangles/edges of convex hull read in a *.tst file output from meshsrfA/B,
c output another *.tst file or tri   file
c
c-------------------------------------------------------
        implicit none
        integer ntrimax
	  integer nseemx
        parameter (ntrimax = 90000 )
        parameter (nseemx = 1000 )
c-------------------------------------------------------
        integer i,j,k,ntri,ntet,k1,i1,j1
        integer i2,i3
        character*80 infile,outfile
	  integer stdin
        character*132 line
c--------------------------------------------
        integer lunit,narg
        integer npoint,np
c
c unique list of points
c
        real*8  xyzpoint(3,ntrimax)
c
c list of points  in each triangle
c
        integer  tripoint(3,ntrimax)
c
c list of each tris neighboring tetrahedra
c
        integer  trineigh(2,ntrimax)
        integer  tetpoint(4,ntrimax) ! points of each tetra
c
c list of each tri neighboring triangles
c
        integer tritri(3,ntrimax)
c
c list of edges
c
	  integer ledge(3,ntrimax),nedge
c
c list of 'seen triangles, with tetra volumes
c
	  integer nsee,isee(nseemx)
	  real*8  seevol(nseemx)
c
c convex hull arrays
c
        integer inchull(ntrimax),nptchull,ifnd,jfnd,kfnd
	  integer ndo
        real*8 tpoint(3),tpoint1(3),color(3)
c
	  real*8 xmid(3),xmin(3),xmax(3),dxyz(3),rdot
	  integer it1,it2,it3,it4
        real*8 vnorm(3),v1(3),v2(3),v3(3),v4(3),rnorm,off
	  real*8 v1max,v2max,v3max,v4max
	  real*8 tvol,tvolmx,thi,thimx,eps,tvoltot,Atot
	  real*8 tvolmn,pyrvol
	  integer ncapt,naddt,nnew
	  real*4 rtod
	  real*4 beg,fin,secnds
	  integer fndtri
	  integer fndedge
	  integer ntrisurf
	  integer nerr
        data color / 1., 0., 0. /
	  data eps / 1.e-3 /
c--------------------------------------------
	  beg = secnds(0.0)
        rtod = 180./acos(-1.)
	  print *,'----------'
	  print *,'KAS implementation of convex hull algorithm,' 
	  print *,'output in qhull format so can be used by travel depth'
	  print *,'input is either from tst file, or if cant find file, '
	  print *,'read input from stdin using qhull format'
	  print *,'----------'
	  narg = iargc()
        lunit = 10
        if(narg.lt.2)then
          print *,'USAGE: khull (input_tst_file) TO outfile'
          stop
	  else if(narg.eq.2)then
          call getarg(2,outfile)
	    stdin = 1 
	  else if(narg.ge.3)then
          call getarg(1,infile)
	    print *,'input file: ',infile
          call getarg(narg,outfile)
	    stdin = 0
        end if
	  if(stdin.eq.0)then
          open(lunit,file=infile,status='old',err=904) ! check if file exists, if not swith to stdin
	    close(lunit)
	    goto 905
904	    stdin = 1
905       continue
	  end if
	  print *,'output file: ',outfile
      npoint = 0
	do k = 1,3
	  xmid(k) = 0.
	  xmin(k) = 1.e6
	  xmax(k) = -1.e6
	end do
	if(stdin.eq.0)then
c--------------------------input from tst file
	  print *,'reading from tst file: ',infile
        open(lunit,file=infile,status='old',err=900)
c
c get list of points from triangle structure file (*.tst)
c
100     read(lunit,'(A)',end=200)line
        if(line(1:14).ne.'POINT_XYZ LIST')goto 100
        print *,line
        read(lunit,'(A)',end=900)line
        do while(line(1:3).ne.'END')
          npoint = npoint + 1
          read(line,'(i8,3f9.4,f15.6)',err=200)i,(xyzpoint(k,npoint),k=1,3)
	    do k = 1,3
	      xmid(k)  = xmid(k) + xyzpoint(k,npoint)
	      xmin(k)  = min(xmin(k),xyzpoint(k,npoint))
	      xmax(k)  = max(xmax(k),xyzpoint(k,npoint))
	    end do
	    inchull(npoint) = 0
          read(lunit,'(A)',end=200)line
        end do
200     close(lunit)
        print *,'POINT_XYZ LIST records: ',npoint
c--------------------------end input from tst file
	else
c--------------------------
	  print *,'reading input from stdin in qhull format'
	  read(5,*,err=901)i1
	  print *,'data dimension: ',i1
	  read(5,*,err=901)npoint
	  do i = 1,npoint
          read(5,*,err=901)(xyzpoint(k,i),k=1,3)
	    do k = 1,3
	      xmid(k)  = xmid(k) + xyzpoint(k,i)
	      xmin(k)  = min(xmin(k),xyzpoint(k,i))
	      xmax(k)  = max(xmax(k),xyzpoint(k,i))
	    end do
	    inchull(npoint) = 0
	  end do
        print *,'# qhull style input records: ',npoint
c--------------------------end input from stdin in qhull format
	end if
	do k = 1,3
	  xmid(k) = xmid(k)/npoint
	end do
	print *,'xmin: ',xmin
	print *,'xmid: ',xmid
	print *,'xmax: ',xmax
c
c check for duplicate points, collinear points
c
	print *,'checking for duplicate points'
	call checkdup(xyzpoint,npoint,ntrimax,inchull)
c201	call checklin(xyzpoint,npoint,ntrimax,inchull,nlinear)
c	if(nlinear.gt.0)goto 201
c
c put out POINTS in form pymol can display
c
	open(23,file='khullp.py')
	write(23,'(a)')'from pymol.cgo import *'
	write(23,'(a)')'from pymol import cmd'
	write(23,'(a)')'import math'
	write(23,'(a)')'obj = ['
	write(23,'(a)')'    LINEWIDTH, 5.0,'
	write(23,'(a)')'   BEGIN, POINTS,'
	do i = 1,npoint
        if(inchull(i).lt.0)then
	    write(23,'(a)')'   COLOR, 0., 0., 1.,'
        else if(inchull(i).eq.0)then
	    write(23,'(a)')'   COLOR, 0., 1., 0.,'
        else 
	    write(23,'(a)')'   COLOR, 1., 1., 1.,'
        end if 
	  do k = 1,3 
	    tpoint(k) = xyzpoint(k,i)
	  end do
	  write(23,'(a,3(f8.3,'',''))')'   VERTEX,',tpoint
	end do
	write(23,'(a)')'   END'
	write(23,'(a)')'   ]'
	write(23,'(a)')"cmd.load_cgo(obj,'khullp')"
	close(23)
c
c find points of initial tetra
c and
c initialize triangle arrays 
c
        ntri = 0
	  do i = 1,ntrimax
          do j = 1,4
            tetpoint(j,i) = 0
          end do
          do j = 1,3
            tripoint(j,i) = 0
            ledge(j,i) = 0
          end do
          do j = 1,2
	      trineigh(j,i) = 0
          end do
	  end do
        nptchull = 0
	  nedge = 0
c construct trial vectors pointing in roughly 4 equi-spaced directions
c and find 4 ponts that most align along them
	  v1max  = 0.
	  v2max = 0.
	  v3max = 0.
	  v4max = 0.
	  v1(1) = 1.1
	  v1(2) = 1.2
	  v1(3) = 0.
c
	  v2(1) = 1.
	  v2(2) = -1.3
	  v2(3) = 0.
c
	  v3(1) = -1.25
	  v3(2) = 0.
	  v3(3) = 1.05
c
	  v4(1) = -1.35
	  v4(2) = 0.
	  v4(3) = -1.15
c
        do i = 1,npoint
	    if(inchull(i).ge.0)then 
c	    if(inchull(i).lt.0)then  ! debug
	      do k = 1,3
	        dxyz(k) = xyzpoint(k,i) - xmid(k)
	      end do
		call dot(dxyz,v1,rdot)
		if(rdot.gt.v1max)then
		  it1 = i
		  v1max = rdot
		end  if
		call dot(dxyz,v2,rdot)
		if(rdot.gt.v2max)then
		  it2 = i
		  v2max = rdot
		end  if
		call dot(dxyz,v3,rdot)
		if(rdot.gt.v3max)then
		  it3 = i
		  v3max = rdot
		end  if
		call dot(dxyz,v4,rdot)
		if(rdot.gt.v4max)then
		  it4 = i
		  v4max = rdot
		end  if
	    end if
        end do
c	it1 = 1
c	it2 = 2
c	it3 = 3
c	it4 = 4
	print *,'initial tetra: ',it1,it2,it3,it4
	write(6,'(3f9.4)')(xyzpoint(k,it1),k=1,3)
	write(6,'(3f9.4)')(xyzpoint(k,it2),k=1,3)
	write(6,'(3f9.4)')(xyzpoint(k,it3),k=1,3)
	write(6,'(3f9.4)')(xyzpoint(k,it4),k=1,3)
	ntet = 1
c tetra points stored in right hand screw order
	tetpoint(1,ntet)  = it1
	tetpoint(2,ntet)  = it2
	tetpoint(3,ntet)  = it3
	tetpoint(4,ntet)  = it4
	nptchull = 4
	inchull(it1) = ntet
	inchull(it2) = ntet
	inchull(it3) = ntet
	inchull(it4) = ntet
	call tetvol(xyzpoint,ntrimax,it1,it2,it3,it4,tvoltot,thi)
	print *,'initial tvol, thi: ',tvoltot,thi
	ntri = 4
	do k = 1,4
	  trineigh(1,k) = ntet
	end do
c
c 4 tri faces of initial tetra, in anticlockwise order from outside view
c
	tripoint(1,1) = it2
	tripoint(2,1) = it3
	tripoint(3,1) = it4
c
	tripoint(1,2) = it1
	tripoint(2,2) = it4
	tripoint(3,2) = it3
c
	tripoint(1,3) = it4
	tripoint(2,3) = it1
	tripoint(3,3) = it2
c
	tripoint(1,4) = it3
	tripoint(2,4) = it2
	tripoint(3,4) = it1
c
c update triangle neighbors- anticlock order from outside
c
	tritri(1,1) = 2
	tritri(2,1) = 3
	tritri(3,1) = 4
	tritri(1,2) = 1
	tritri(2,2) = 4
	tritri(3,2) = 3
	tritri(1,3) = 4
	tritri(2,3) = 1
	tritri(3,3) = 2
	tritri(1,4) = 3
	tritri(2,4) = 2
	tritri(3,4) = 1
c
c update edge list
c
	nedge = 6
	ledge(1,1) = it1
	ledge(2,1) = it2
	ledge(3,1) = 2  ! number of exposed triangles attached to this edge
	ledge(1,2) = it1
	ledge(2,2) = it3
	ledge(3,2) = 2  ! number of exposed triangles attached to this edge
	ledge(1,3) = it1
	ledge(2,3) = it4
	ledge(3,3) = 2  ! number of exposed triangles attached to this edge
	ledge(1,4) = it2
	ledge(2,4) = it3
	ledge(3,4) = 2  ! number of exposed triangles attached to this edge
	ledge(1,5) = it2
	ledge(2,5) = it4
	ledge(3,5) = 2  ! number of exposed triangles attached to this edge
	ledge(1,6) = it3
	ledge(2,6) = it4
	ledge(3,6) = 2  ! number of exposed triangles attached to this edge
c
c check face norms
c
	do k  = 1,3
	  xmid(k) = (xyzpoint(k,it1)+xyzpoint(k,it2)+xyzpoint(k,it3)+xyzpoint(k,it4))/4.
	end do
	do i = 1,ntri
	  do k = 1,3
	    v1(k) = xyzpoint(k,tripoint(2,i)) - xyzpoint(k,tripoint(1,i)) 
	    v2(k) = xyzpoint(k,tripoint(3,i)) - xyzpoint(k,tripoint(1,i)) 
	    v4(k) = xmid(k) - xyzpoint(k,tripoint(1,i)) 
	  end do
	  call norm(v1,rnorm)
	  call norm(v2,rnorm)
	  call cross(v1,v2,v3)
	  call norm(v3,rnorm)
	  call dot(v3,v4,rdot)
	  write(6,'(i6,a,4f9.4)')i,'tri norm: ',(v3(k),k=1,3),rdot
	  if(rdot.gt.0.)then
	    print *,'norm pointing in!!, swapping order'
	    i1 = tripoint(1,i)
	    tripoint(1,i) = tripoint(3,i)
	    tripoint(3,i) = i1
	  end if
	end do
	print *,'checking edges of initial tetrahedron...'
	call chkedge(ledge,nedge,tripoint,trineigh,ntri,xyzpoint,ntrimax,nerr)
	print *,'# errors: ',nerr
c	goto 800 ! debug
c
c check each point against each surface tri (i.e those with only 1 tetra neighbor)
c if new tetra has + height/volume, add tetrahedron and new faces
c
	np = npoint/20
	print *,'0%              100%'
	do i = 1,npoint
	  if(mod(i,np).eq.0)then
c	    print *,'doing point ',i
	    write(6,'(a,$)')'.'
	  end if
	  if(inchull(i).eq.0)then
	    ncapt = 0
	    naddt = ntri
	    tvolmn = 1.e6
	    tvolmx  = -1.e6
	    thimx  = -1.e6
	    nsee = 0
	    pyrvol = 0.
	    do ndo = 1,ntri
	      if(trineigh(2,ndo).eq.0)then
	        it1 = tripoint(1,ndo)
	        it2 = tripoint(2,ndo)
	        it3 = tripoint(3,ndo)
		  call tetvol(xyzpoint,ntrimax,it1,it2,it3,i,tvol,thi)
c		  if(thi.gt.eps)then
c		  if(thi.ge.0.)then
c		    print *,'add tet: ',i,ndo,tvol,thi,it1,it2,it3
c		  if(tvol.gt.eps)then
		  if(tvol.gt.1.e-6)then
c store list of seen triangles
		    nsee = nsee + 1
		    if(nsee.gt.nseemx)then
		      print *,'nsee > nseemx'
			stop
		    endif
		    isee(nsee) = ndo
		    seevol(nsee) = tvol
		    pyrvol = pyrvol + tvol
		  else
c		    print *,'No tet: ',i,ndo,tvol,thi,it1,it2,it3
		  end if ! if this triangle can be 'seen' by current point
		end if ! if surface triangle
	    end do ! list of triangles
c	    if(nsee.gt.0)then
c	      print *,' # of visible triangles, pyramid volume: ',nsee,pyrvol
c	    end if
c
c process list of seen triangles
c
	    do j = 1,nsee
	      ndo = isee(j)
		tvol = seevol(j)
	      it1 = tripoint(1,ndo)
	      it2 = tripoint(2,ndo)
	      it3 = tripoint(3,ndo)
c add new tetrahedron
		    tvoltot = tvoltot + tvol
		    tvolmn = min(tvolmn,tvol)
		    tvolmx  = max(tvolmx,tvol)
		    thimx  = max(thimx,thi)
		    ntet  = ntet + 1
		    tetpoint(1,ntet) = it1
		    tetpoint(2,ntet) = it2
		    tetpoint(3,ntet) = it3
		    tetpoint(4,ntet) = i
c assign point to convex hull
		    if(inchull(i).eq.0)then
		      inchull(i) = ntet
		      nptchull  = nptchull + 1
		    end if
c this surface tri gets covered by new tet
		    trineigh(2,ndo) = ntet
c		    print *,'cap base tri: ',it1,it2,it3
c check remaining 3 faces to see if they are new
	          ifnd = fndtri(tripoint,ntrimax,naddt,it2,it3,i)
		    if(ifnd.eq.0)then ! new tri
		      naddt = naddt + 1
			tripoint(1,naddt) = it2
			tripoint(2,naddt) = it3
			tripoint(3,naddt) = i
			trineigh(1,naddt) = ntet
c			call addneigh(tritri,tripoint,ntrimax,naddt,ndo,it2,it3)
c			print *,'add tri: ',it2,it3,i
		    else
			trineigh(2,ifnd) = ntet
c			print *,'cap side tri: ',it2,it3,i
			ncapt = ncapt + 1
		    end  if
	          jfnd = fndtri(tripoint,ntrimax,naddt,it1,i,it3)
		    if(jfnd.eq.0)then ! new tri
		      naddt = naddt + 1
			tripoint(1,naddt) = it1
			tripoint(2,naddt) = i
			tripoint(3,naddt) = it3
			trineigh(1,naddt) = ntet
c			call addneigh(tritri,tripoint,ntrimax,naddt,ndo,it1,it3)
c			print *,'add tri: ',it1,i,it3
		    else
			trineigh(2,jfnd) = ntet
c			print *,'cap side tri: ',it1,i,it3
			ncapt = ncapt + 1
		    end  if
	          kfnd = fndtri(tripoint,ntrimax,naddt,i,it1,it2)
		    if(kfnd.eq.0)then ! new tri
		      naddt = naddt + 1
			tripoint(1,naddt) = i
			tripoint(2,naddt) = it1
			tripoint(3,naddt) = it2
			trineigh(1,naddt) = ntet
c			call addneigh(tritri,tripoint,ntrimax,naddt,ndo,it1,it2)
c			print *,'add tri: ',i,it1,it2
		    else
			trineigh(2,kfnd) = ntet
c			print *,'cap side tri: ',i,it1,it2
			ncapt = ncapt + 1
		    end  if
c
c add new edges
c
		    ifnd = fndedge(ledge,ntrimax,nedge,i,it1)
		    if(ifnd.eq.0)then
		      nedge = nedge + 1
			ledge(1,nedge) = i
			ledge(2,nedge) = it1
			ledge(3,nedge) = 2
		    end if
		    jfnd = fndedge(ledge,ntrimax,nedge,i,it2)
		    if(jfnd.eq.0)then
		      nedge = nedge + 1
			ledge(1,nedge) = i
			ledge(2,nedge) = it2
			ledge(3,nedge) = 2
		    end if
		    kfnd = fndedge(ledge,ntrimax,nedge,i,it3)
		    if(kfnd.eq.0)then
		      nedge = nedge + 1
			ledge(1,nedge) = i
			ledge(2,nedge) = it3
			ledge(3,nedge) = 2
		    end if
	    end do
	    nnew = (naddt - ntri) - ncapt
	    ntri   = naddt
c	    if(nnew.gt.0)print *,nnew,' sided pyramid added for point',i,tvolmn,ntet
c	    print *,i,'max vol, hi: ',tvolmx,thimx
c
c check edges
c
c	    if(nnew.gt.0)then
c	      call chkedge(ledge,nedge,tripoint,trineigh,ntri,xyzpoint,ntrimax,nerr)
c	      if(nerr.gt.0)print *,'# of errors: ',nerr
c	    end if
c	    if(nerr.gt.0)goto 800
cc debug	
c	      if(ntet.ge.238)then
c		  goto 800
c		end if
cc debug	
	  end if
	end do ! list of points
	call chkedge(ledge,nedge,tripoint,trineigh,ntri,xyzpoint,ntrimax,nerr)
	print *,'# of errors: ',nerr
800	print *,'# of points examined: ',nptchull
	ntrisurf = 0
	do i  = 1,npoint
	  inchull(i) = 0
	end  do
	nptchull = 0
	Atot = 0.
	do i = 1,ntri
	  if(trineigh(2,i).eq.0)then
	    ntrisurf = ntrisurf + 1
	    do k = 1,3
	      v1(k) = xyzpoint(k,tripoint(2,i)) - xyzpoint(k,tripoint(1,i))
	      v2(k) = xyzpoint(k,tripoint(3,i)) - xyzpoint(k,tripoint(1,i))
	    end do
	    call cross(v1,v2,v3)
	    call norm(v3,rnorm)
	    Atot = Atot + 0.5*rnorm
	    do j  = 1,3
            i1 = tripoint(j,i)
	      if(inchull(i1).eq.0)then
	        inchull(i1)   = 1
	        nptchull = nptchull +1
	      end if
	    end do
	  end if
	end do
	print *,'total tri, tetra: ',ntri,ntet
	print *,'# convex hull pts, tri: ',nptchull,ntrisurf
	print *,'convex hull volume, Area: ',tvoltot,Atot

	fin = secnds(beg)
	print *,'cpu time: ',fin
c      print *,'--------'
c      print *,'points'
c      do i=1,npoint
c        if(inchull(i).gt.0)print *,i,inchull(i)
c      end do
      print *,'--------'
      open(17,file='khull.tri')
	write(17,'(a)')'PDB_RECORD'
	write(17,'(a)')'END PDB_RECORD'
      write(17,'(A)')'TRIANGLE_XYZ'
      do i = 1,ntri
	  if(trineigh(2,i).eq.0)then
        do j = 1,3
          i1 = tripoint(j,i)
          write(17,'(1x,3f9.3,1x,A1)')(xyzpoint(k,i1),k=1,3)
        end do
	  end if
      end do
      write(17,'(A)')'END TRIANGLE_XYZ'
      close(17)
c
c put out lines in form pymol can display
c
	color(3) = 1.
	color(2) = 0.
	color(1) = 0.
	open(23,file='khull.py')
	write(23,'(a)')'from pymol.cgo import *'
	write(23,'(a)')'from pymol import cmd'
	write(23,'(a)')'import math'
	write(23,'(a)')'obj = ['
	write(23,'(a)')'   BEGIN, LINES,'
	do  i = 1,ntri
	  if(trineigh(2,i).eq.0)then
	  do j = 1,3
	    if(j.eq.3)then
	      j1 = 1
	    else
	      j1 = j + 1
	    end if
	    i1 = tripoint(j,i)
	    i2  = tripoint(j1,i)
	    do k = 1,3
	      tpoint(k) = xyzpoint(k,i1)
	      tpoint1(k) = xyzpoint(k,i2)
	    end do
          write(23,'(a,3(f8.3,'',''))')'   COLOR,',(color(k),k=1,3)
          write(23,'(a,3(f8.3,'',''))')'   VERTEX,',tpoint
          write(23,'(a,3(f8.3,'',''))')'   VERTEX,',tpoint1
	  end do
	  end if
	end do
      write(23,'(a)')'   END'
      write(23,'(a)')'   ]'
      write(23,'(a)')"cmd.load_cgo(obj,'khull')"
      close(23)
c
c output in qhull format
c note qhull indexes points from 0, so need to subtract 1
c and also  clockwise order looking from outside
c
      open(17,file=outfile)
	off = 0.
	write(17,'(i8)')ntrisurf
	do i = 1,ntri
	  if(trineigh(2,i).eq.0)then
          write(17,'(3i8)')(tripoint(k,i)-1,k=3,1,-1)
	  end if
      end do
	write(17,'(a)')'4'
	write(17,'(i8)')ntrisurf
	do i = 1,ntri
	  if(trineigh(2,i).eq.0)then
	    do k = 1,3
	      v1(k) = xyzpoint(k,tripoint(2,i)) - xyzpoint(k,tripoint(1,i))
	      v2(k) = xyzpoint(k,tripoint(3,i)) - xyzpoint(k,tripoint(1,i))
	    end do
	    call cross(v1,v2,v3)
	    call norm(v3,rnorm)
          write(17,'(4f10.5)')(v3(k),k=1,3),off
	  end if
      end do
      close(17)
	
      stop
c
c----------------------------------------------------
900     print *,'cant find triangle structure file: ',infile
        stop
901     print *,'error reading qhull style input from stdin: ',i
        stop
        end


c------------------------------------------------------
      SUBROUTINE CROSS(U,V,W)
	implicit none
      real*8 u(3),v(3),w(3)
      w(1) = u(2)*v(3) - u(3)*v(2)
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)
      return
      end
c------------------------------------------------------
      SUBROUTINE DOT(U,V,RDOT)
	implicit none
      real*8 u(3),v(3),rdot
      rdot = u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
      return
      end
c------------------------------------------------------
      SUBROUTINE NORM(U,SIZE)
	implicit none
      real*8 u(3),size
      call dot(u,u,size)
      if(size.ne.0.) then
        size = sqrt(size)
        u(1) = u(1)/size
        u(2) = u(2)/size
        u(3) = u(3)/size
      else
c        type *,'zero vector'
      end if
      return
      end
c-------------------------------------------------------

      function ran1(iseed)
      iseed = iseed*69069 + 1
        ran1 = real(ishft(iseed,-8))*0.5**24
        return
      end



	function fndtri(tripoint,ntrimax,ntri,ic1,ic2,ic3)
c
c search triangle point list for ic1,ic2 & ic3
c return triangle number
c	
	implicit none
	integer fndtri
	integer ntrimax,ntri,ic1,ic2,ic3
      integer  tripoint(3,ntrimax)
	integer i,j,k,kfnd
c-----------------
	fndtri = 0
      do i = 1,ntri
        kfnd = 0
        if((ic1.eq.tripoint(1,i)).or.(ic1.eq.tripoint(2,i)).or.(ic1.eq.tripoint(3,i)))kfnd = kfnd + 1
        if((ic2.eq.tripoint(1,i)).or.(ic2.eq.tripoint(2,i)).or.(ic2.eq.tripoint(3,i)))kfnd = kfnd + 1
        if((ic3.eq.tripoint(1,i)).or.(ic3.eq.tripoint(2,i)).or.(ic3.eq.tripoint(3,i)))kfnd = kfnd + 1
        if(kfnd.eq.3)then
	    fndtri = i
	    return
	  end if
      end do
	return
	end


	subroutine checkdup(xyzpoint,npoint,ntrimax,inchull)
	implicit none
	integer ntrimax
	integer npoint, inchull(ntrimax)
	real*8 xyzpoint(3,ntrimax)
	integer i,j,k,ndup,kfnd
c-----------------
	ndup = 0
	do i = 1,npoint
	  if(inchull(i).lt.0)goto 401
	  do j = i+1,npoint
	    if(inchull(j).lt.0)goto 402
	    kfnd = 0
	    do k = 1,3
	      if(xyzpoint(k,i).eq.xyzpoint(k,j))kfnd = kfnd + 1
	    end do
	    if(kfnd.eq.3)then
		ndup = ndup + 1
		inchull(j) = -1
	      print *,'duplicate points: ',i,j
	    end if
402	    continue
	  end do
401	  continue
	end do
	print *,'# of duplicate points: ',ndup
	return
	end

	subroutine checklin(xyzpoint,npoint,ntrimax,inchull,nlinear)
	implicit none
	integer ntrimax
	integer npoint, inchull(ntrimax)
	real*8 xyzpoint(3,ntrimax)
      real*8 vnorm(3),v1(3),v2(3),v3(3),rnorm
	real*8 rdot1,rdot2,rdot3
	integer i,j,k,nlinear,k1
	real*8 eps,tol
	data eps / 1.e-4 /
c-----------------
	nlinear = 0
	tol = 1.-eps
	do i = 1,npoint
	  if(inchull(i).lt.0)goto 401
	  do j = i+1,npoint
	    if(inchull(j).lt.0)goto 402
	    do k1 = 1,3
	      v1(k1)   = xyzpoint(k1,j) - xyzpoint(k1,i)
	    end do
	    call norm(v1,rnorm)
	    do k = j+1,npoint
	      if(inchull(k).lt.0)goto 403
	      do k1 = 1,3
	        v2(k1)   = xyzpoint(k1,k) - xyzpoint(k1,j)
	        v3(k1)   = xyzpoint(k1,i) - xyzpoint(k1,k)
	      end do
	      call norm(v2,rnorm)
	      call norm(v3,rnorm)
		call dot(v1,v2,rdot1)
		call dot(v2,v3,rdot2)
		call dot(v3,v1,rdot3)
		if(abs(rdot1).gt.tol)then
		  nlinear = nlinear + 1
		  print *,'lin:  ',i,j,k
		  if(rdot1.gt.0.)then ! j is middle point
		    inchull(j) = -1
		    goto 402
		  end if
		  if(rdot2.gt.0.)then ! k is middle point
		    inchull(k) = -1
		    goto 403
		  end if
		  if(rdot3.gt.0.)then ! i is middle point
		    inchull(i) = -1
		    goto 401
		  end if
		end if
403	      continue
	    end do
402	    continue
	  end do
401	  continue
	end do
	print *,'# collinear points: ',nlinear
	return
	end

	subroutine tetvol(xyzpoint,ntrimax,it1,it2,it3,it4,tvol,thi)
c
c volume, height of tetra formed by it1,it2,it3,i
c
	implicit none
	integer ntrimax
	integer it1,it2,it3,it4,i,j,k
	real*8 xyzpoint(3,ntrimax)
      real*8 vnorm(3),v1(3),v2(3),v3(3),rnorm
	real*8 rdot1,rdot2,rdot3
	real*8 tvol,thi
c-----------------
	do k = 1,3
	  v1(k)  = xyzpoint(k,it2) - xyzpoint(k,it1)
	  v2(k)  = xyzpoint(k,it3) - xyzpoint(k,it1)
	  v3(k)  = xyzpoint(k,it4) - xyzpoint(k,it1)
	end  do
	call cross(v1,v2,vnorm)
	call dot(vnorm,v3,rdot1)
	tvol = rdot1/6.
	call norm(vnorm,rnorm)
	call dot(vnorm,v3,rdot2)
	thi = rdot2
	return
	end

	subroutine addneigh(tritri,tripoint,ntrimax,newt,oldt,pold1,pold2)
c
c replace triangle oldt with newt in oldt's neighbors' lists
c
	implicit none
	integer ntrimax
	integer tritri(3,ntrimax)
	integer tripoint(3,ntrimax)
	integer newt,oldt,pold1,pold2
	integer i,j,k,iside,jside
	integer i1,i2
	integer neigh
c
c check all neighbors for one which has edge pold1-pold2
	neigh = 0
	do i = 1,3
	  iside = tritri(i,oldt)
	  i1 = tripoint(1,iside)
	  i2 = tripoint(2,iside)
	  if(((i1.eq.pold1).and.(i2.eq.pold2)).or.((i2.eq.pold1).and.(i1.eq.pold2)))neigh = iside
	  i1 = tripoint(2,iside)
	  i2 = tripoint(3,iside)
	  if(((i1.eq.pold1).and.(i2.eq.pold2)).or.((i2.eq.pold1).and.(i1.eq.pold2)))neigh = iside
	  i1 = tripoint(3,iside)
	  i2 = tripoint(1,iside)
	  if(((i1.eq.pold1).and.(i2.eq.pold2)).or.((i2.eq.pold1).and.(i1.eq.pold2)))neigh = iside
	end do
	if(neigh.eq.0)then
	  print *,'ERROR: cannot find edge ',pold1,pold2,' in ',oldt,' neighbors'
	  stop
	else
	  jside = 0
	  do i = 1,3
	    if(tritri(i,neigh).eq.oldt)jside = i
	  end do
	  if(jside.eq.0)then
	    print *,'ERROR: cannot find ',oldt,' as neighbor of ',neigh
	    stop
	  else
	    print *,'replacing neighbor ',oldt,' for ',neigh,' with ',newt
	    tritri(jside,neigh) = newt
	    tritri(jside,neigh) = newt
	  end if
	end if
	return
	end

	function fndedge(ledge,ntrimax,nedge,ic1,ic2)
c
c search edge point list for ic1,ic2 
c return edge number
c	
	implicit none
	integer fndedge
	integer ntrimax,nedge,ic1,ic2
      integer  ledge(3,ntrimax)
	integer i,j,k,kfnd
c-----------------
	fndedge = 0
      do i = 1,nedge
        kfnd = 0
        if((ic1.eq.ledge(1,i)).or.(ic1.eq.ledge(2,i)))kfnd = kfnd + 1
        if((ic2.eq.ledge(1,i)).or.(ic2.eq.ledge(2,i)))kfnd = kfnd + 1
        if(kfnd.eq.2)then
	    fndedge = i
	    return
	  end if
      end do
	return
	end

	subroutine chkedge(ledge,nedge,tripoint,trineigh,ntri,xyzpoint,ntrimax,nerr)
c check edges
	implicit none
	integer ntrimax,nedge,ntri
	integer nerr
	integer nbury,nsurf
      integer  ledge(3,ntrimax)
      integer  trineigh(2,ntrimax),tripoint(3,ntrimax)
	real*8 xyzpoint(3,ntrimax)
	real*8 tvol,thi
	integer ie1,ie2,it1,it2,it3,it4,nexp,iexp(2),ifnd 
	integer i,j,k
	real*8 eps
	data eps / 1.e-4 /
c
c for each edge, if buried (ledge(3,i) = 0) skip
c else find all triangles bounding this edge , that are exposed
c if none, then bury this edge, if 2 then check angle, else error
c
	nerr = 0
	nbury = 0
	nsurf = 0
	do i = 1,nedge
	  if(ledge(3,i).ne.0)then
	    ie1 = ledge(1,i)
	    ie2 = ledge(2,i)
	    nexp = 0
	    do j = 1,ntri
	      if(trineigh(2,j).eq.0)then
		  it1 = tripoint(1,j)
		  it2 = tripoint(2,j)
		  it3 = tripoint(3,j)
	        ifnd = 0
		  if((it1.eq.ie1).or.(it1.eq.ie2))ifnd = ifnd + 1
		  if((it2.eq.ie1).or.(it2.eq.ie2))ifnd = ifnd + 1
		  if((it3.eq.ie1).or.(it3.eq.ie2))ifnd = ifnd + 1
		  if(ifnd.gt.2)then
		    print *,'error edge : ',ie1,ie2
		    print *,'overmatches triangle : ',it1,it2,it3
		    nerr = nerr + 1
		  else 
		    if(ifnd.eq.2) then
		      nexp = nexp + 1
			if(nexp.gt.2)then
			  print *,'error: edge has > 2 exposed triangles'
		        nerr = nerr + 1
c			  stop
			else
			  iexp(nexp) = j
			end if
		    end if
		  end if
		end if ! exposed triangles
	    end do ! triangle list
	    if(nexp.eq.0)then
c		print *,'burying edge: ',ie1,ie2
	      ledge(3,i) = 0
		nbury = nbury + 1
	    else if(nexp.eq.1)then
	      print *,'error: edge has only 1 exposed triangle: ',ie1,ie2,iexp(1)
		nerr = nerr + 1
c	      stop
	    else
c reached here, edge has 2 exposed triangles- check that angle is convex
		nsurf = nsurf + 1
		it1 = tripoint(1,iexp(1))
		it2 = tripoint(2,iexp(1))
		it3 = tripoint(3,iexp(1))
		it4 = 0
		do k = 1,3
		  if((tripoint(k,iexp(2)).ne.ie1).and.(tripoint(k,iexp(2)).ne.ie2))it4 = tripoint(k,iexp(2))
		end do
		if(it4.eq.0)then
		  print *,'error: failed to find non edge point in 2nd triangle: ',iexp(2)
		  stop
		end if
		call tetvol(xyzpoint,ntrimax,it1,it2,it3,it4,tvol,thi)
c		if(thi.gt.eps)then
		if(thi.gt.1.e-6)then
		  print *,'WARNING: edge is concave: ',ie1,ie2,iexp,thi
		nerr = nerr + 1
		nerr = nerr + 1
		else
c		  print *,'OK: edge is convex: ',ie1,ie2,iexp,thi
		end if
	    end if
	  end if ! exposed edges
	end do ! edge list
	print *,'nedge,nbury,nsurf: ',nedge,nbury,nsurf

	return
	end
