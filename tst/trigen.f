	program trigen
c
c read in output from
c mesh srf and check triangulation
c generate triangluated surface structure, normals, 
c point and atom 'labels'
c
c input is triples of points forming triangles
c in format 1X, 3f9.3
c
	implicit none
	integer ntrimax,natmax
	integer nnmax,nlinemx
	parameter (ntrimax = 90000 )
	parameter (natmax = 40000 )
	parameter (nnmax = 12 )
	parameter ( nlinemx = ntrimax/10 )
c--------------------------------------------
	character*80 trifile,line,tstfile,recfile
c--------------------------------------------
	real*4 triangle(4,3,ntrimax),tpoint(3)
	real*4 tnorm(3),tnorm1(3),tnorm2(3)
	real*4 tpoint1(3)
	real*4 tpoint2(3)
	real*4 tri1(3,3)
	real*4 v1(3),v2(3),v3(3)
	real*4 orderav
	character*80 pdbrec(natmax)
c--------------------------------------------
	integer i,j,k,ntri,i1,j1,k1,i2,j2,k2
	integer itri
	integer lunit,lunit1
	integer nrec
	integer pair(2,nnmax)
c
c triangle neighbor list: neighbor on 1-2 side in 1st, 1-3 side 2nd, and 2-3 side 3rd
c
	integer trineigh(3,ntrimax)
c
c stack and list for running through each triangle once
c
c
c unique list of points
c
	real*4  xyzpoint(3,ntrimax),xyznorm(3,ntrimax),xyzav(3),rnorm
	real*4 curve(ntrimax),dpoint(3),dotdn,dotdd,crvmin,crvmax,crvav
	real*4 curve1(ntrimax)
	real*4 pproperty(ntrimax)
	integer pntcurv(nlinemx),npntcurv
	integer npoint,neuler
c
c general size etc
c
	real*4  xyzmax(3),xyzmin(3),xyzran(3),xyzmid(3),rog
c
c list of points  in each triangle
c	
	integer tripoint(3,ntrimax)
c	
c list of each point's neighboring points
c	
	integer pointneigh(nnmax,ntrimax)
	integer pointorder(ntrimax)
c	
c list of each points neighboring triangles
c
	integer pointtri(nnmax,ntrimax)
	integer orderfrq(nnmax)
c
c subset each point belongs to
c
	integer isubset(ntrimax),nsubset,nhandle,nedge,npoint_good,nedge_good
c
	real*8 Atot,Vtot
	real*4 Atri,Vtri
c
c	pdb atom record associated with each triangle, vertex
c
	real*4 pdbxyz(3),pdbrad,trimid(3),dist(3),distmid
	real*4 distmin(3),distmidmin
	integer pointrec(ntrimax),trirec(ntrimax),vertex_pdb(3,ntrimax)
	include 'tricol.h'
	real*4 vmid,vup,vlow,val

c--------------------------------------------
c	read(5,*)trifile
	if(iargc().lt.2)then
	  print *,'generate triangle structure from meshsrf output list of triangle points'
	  print *,'USAGE: trigen meshsrf_input_file (*.tri) output triangle structure file (*.tst)'
	  stop
	end if
	call getarg(1,trifile)
	call getarg(2,tstfile)
	lunit = 10
	open(lunit,file=trifile,status='old',err=900)
100	read(lunit,'(A)',end=200)line
	if(line(1:).eq.'PDB_RECORD')then
	  call getpdbrec(lunit,pdbrec,natmax,nrec)
	else if(line(1:).eq.'TRIANGLE_XYZ')then
	  call gettrixyz(lunit,triangle,ntrimax,ntri)
	  open(9,file='trisrf.rec',status='old',err=902)
	  print *,'reading vertex to pdb atom mapping...'
	  read(9,'(A)',end=200)line
	  print *,line
	  do i = 1,ntri
	    do k = 1,3
	      read(9,'(2i8)',err=902)i1,vertex_pdb(k,i)
	    end do
	  end do
	  close(9)
	  goto 100
902	  print *,'cant find vertex to pdb atom mapping file trisrf.rec'
	  print *,' or records missing after: ',i
	  do i = 1,ntri
	    do k = 1,3
	      vertex_pdb(k,i) = 0
	    end do
	  end do
	end if
	goto 100
200	close(lunit)
	print *,'generating vertex list...'
	npoint = 0
	do i = 1,ntrimax
	  pointorder(i) = 0
	end do
	do i = 1,ntri
	  do j = 1,3 ! triangle vertexes
	    do i1 = 1,npoint
	      if(triangle(1,j,i).ne.xyzpoint(1,i1))goto 300
	      if(triangle(2,j,i).ne.xyzpoint(2,i1))goto 300
	      if(triangle(3,j,i).ne.xyzpoint(3,i1))goto 300
c if we reach here, found previous point
	      tripoint(j,i) = i1
	      pointorder(i1) = pointorder(i1) + 1
	      pointtri(pointorder(i1),i1) = i
	      goto 400 ! next triangle point
300	      continue ! not a match- go to next vertex of current triangle
	    end do
c if we reach here, it is a new point
	    npoint = npoint + 1
	    do k = 1,3
	      xyzpoint(k,npoint) = triangle(k,j,i)
	      pproperty(npoint) = triangle(4,j,i)
	    end do
	    tripoint(j,i) = npoint
	    pointorder(npoint) = pointorder(npoint) + 1
	    pointtri(pointorder(npoint),npoint) = i
	    pointrec(npoint) = vertex_pdb(j,i)
400	    continue ! with next triangle vertex
	  end do
c map triangle to  atom record - if 2/3 points map to same point, use that, otherwise use first point'
	  if(pointrec(tripoint(2,i)).eq.pointrec(tripoint(3,i)))then
	    trirec(i) = pointrec(tripoint(2,i))
	  else 
	    trirec(i) = pointrec(tripoint(1,i))
	  end if
	end do
c	open(21,file='trigen.dat')
c	do i = 1,ntri
c	  do j = 1,3
c	    tpoint(j) = triangle(j,1,i)
c	    tpoint1(j) = triangle(j,2,i)
c	    tpoint2(j) = triangle(j,3,i)
c	  end do
c	  write(21,'(a,3(3f8.3,a,3f8.3,a))')'draw triangle { ',tpoint,' } { ',tpoint1,' } { ',tpoint2,' }'
c	end do
c	close(21)
	print *,'# of unique points: ',npoint
c	print *,'xyzpoint list:'
c	do i = 1,npoint
c	  write(6,'(i6,3f9.4,f15.6)')i,(xyzpoint(k,i),k=1,3)
c	end do
c	print *,' '
c	print *,'triangle point list (clockwise really, this was checked):'
c	do i = 1,ntri
c	  write(6,'(i6,3i5)')i,(tripoint(k,i),k=1,3)
c	end do
c	print *,' '
c	print *,'point triangle list (clockwise):'
	do i = 1,nnmax
	  orderfrq(i) = 0
	end do
	nedge = 0
	do i = 1,npoint
	  nedge = nedge + pointorder(i)
	  if(pointorder(i).gt.nnmax)then
	    print *,'WARNING: max point order exceeded: ',pointorder(i)
	  else
	    orderfrq(pointorder(i))  = orderfrq(pointorder(i)) + 1
	  end if
c	  write(6,'(i6,11i5)')i,(pointtri(k,i),k=1,pointorder(i)),pointorder(i)
	end do
	print *,' '
	orderav = 3.*ntri/npoint
	print *,'mean point order by F/V = ',orderav
	print *,'Point order frequency: '
	order av = 0
	do i = 1,nnmax
	  orderav = orderav + i*orderfrq(i)
	  write(6,'(2I6)')i,orderfrq(i)
	end do
	orderav = orderav/npoint
	nedge = nedge / 2
	print *,'mean point order = ',orderav
	print *,' '
c 
c generate point neighbor list
c
	print *,'generating vertex conectivity...'
	do i = 1,npoint
	  do i1 = 1,pointorder(i)
c for each point, get the triangles it belongs to and pull out other pair of points
c in cyclic order, because triangle vertices are ordered clockwise
	    itri = pointtri(i1,i)
	    do j = 1,3 ! all triangle points
	      j1 = tripoint(j,itri)
	      if(j1.eq.i)then
	        if(j.eq.1)then
		  pair(1,i1) = tripoint(2,itri)
		  pair(2,i1) = tripoint(3,itri)
	        else if(j.eq.2)then
		  pair(1,i1) = tripoint(3,itri)
		  pair(2,i1) = tripoint(1,itri)
		else ! j = 3
		  pair(1,i1) = tripoint(1,itri)
		  pair(2,i1) = tripoint(2,itri)
		end if
	      end if
	    end do
	  end do
c
c	  print *,'point neighbors: ',i
c	  do i1 = 1,pointorder(i)
c	    print *,pair(1,i1),pair(2,i1)
c	  end do
	  i1 = 1
	  j1 = 1
	  do while (i1.le.pointorder(i))
	    pointneigh(i1,i) = pair(1,j1)
	    k1 = pair(2,j1)
c	    print *,'next in neighborhood: ',pair(1,j1),k1
	    j1 = 1
	    do while(pair(1,j1).ne.k1)
	      j1 = j1 + 1
c		print *,'j1: ',j1
	    end do
	    i1 = i1 + 1
	  end do
	end do
	call dostack(ntrimax,nnmax,pointneigh,pointorder,npoint,nsubset,isubset)
	nedge_good = 3*ntri/2
	print *,'Faces, Vertices, Edges: ',ntri,npoint,nedge
	print *,'# of surfaces (S): ',nsubset
	npoint_good = ntri/2 + 2*nsubset
	print *,'expected, actual points: ',npoint_good,npoint
	print *,'expected, actual edges: ',nedge_good,nedge
	if(npoint.ne.npoint_good)then
	print *,'WARNING: atucal number of points defferent from expected'
	end if
	neuler = ntri + npoint - 3*ntri/2 - 2*nsubset
	nhandle = -neuler/2
	print *,'F + V - E - 2S = ',neuler
	print *,'# of handles: ',nhandle


c	print *,'point neighbor list (clockwise):'
	print *,'generating outward normals...'
c note point neighbors are ordered clockwise
c calculate normal of each triangle surrounding point i, and average normal vector
	do i = 1,npoint
	  do k = 1,3
	    xyzav(k) = 0
	  end do
	  do j = 1,pointorder(i)-1
	    do k = 1,3
	      v1(k) = xyzpoint(k,pointneigh(j,i)) - xyzpoint(k,i)
	      v2(k) = xyzpoint(k,pointneigh(j+1,i)) - xyzpoint(k,i)
	    end do
	    call cross(v2,v1,v3)
	    call norm(v3,rnorm)
	    do k = 1,3
	       xyzav(k) = xyzav(k) + v3(k)
	    end do
	  end do
	  do k = 1,3
	      v1(k) = xyzpoint(k,pointneigh(pointorder(i),i)) - xyzpoint(k,i)
	      v2(k) = xyzpoint(k,pointneigh(1,i)) - xyzpoint(k,i)
	  end do
	  call cross(v2,v1,v3)
	  do k = 1,3
	    xyzav(k) = xyzav(k) + v3(k)
	  end do
	  do k = 1,3
	    xyzav(k) = xyzav(k)/pointorder(i)
	  end do
	  call norm(xyzav,rnorm)
cc 1st neighbor code
c	  do k = 1,3
c	    xyznorm(k,i) = xyzav(k)
c	  end do
cc	  write(6,'(i6,11i5)')i,(pointneigh(k,i),k=1,pointorder(i)),pointorder(i)
c	end do
c 2nd neighbor neighbor code
	  do k = 1,3
	    triangle(k,1,i) = xyzav(k) ! triangle( is reused as temporary storage
	  end do
	end do
c now replace each normal by average of itself plus neighboring normals
	print *,'normals are average of vertex and 1st neighbors'
	do i = 1,npoint
	  do k = 1,3
	    xyzav(k) = 0
	  end do
	  do j = 1,pointorder(i)
	    j1 = pointneigh(j,i)
	    do k = 1,3
	       xyzav(k) = xyzav(k) + triangle(k,1,j1)
	    end do
	  end do
	  do k = 1,3
	    xyzav(k) = 0.1*triangle(k,1,i) + 0.9*xyzav(k)/pointorder(i)
	  end do
	  call norm(xyzav,rnorm)
	  do k = 1,3
	    xyznorm(k,i) = xyzav(k)
	  end do
	end do
c	print *,'Normals '
c	do i = 1,npoint
c	  write(6,'(i6,3f9.4)')i,(xyznorm(k,i),k=1,3)
c	end do
c
c compute radius of curvature (actually 1/r so flat = 0, not infinity!)
c
c treat each line joining point i to its neighbor points as a segment of polygon
c inscribed in a circle, and compute 1/radius of that circle- this defintiion relys on normal
c vector
c now expand to second neighbors, as 1st neighbors too close to get good curv. estimate
	crvmin = 1.e6
	crvmax = -1.e6
	crvav = 0.
c 2nd neighbor code
	do i = 1,npoint
	  npntcurv = 0
	  curve(i) = 0.
	  do j = 1,pointorder(i)
	    i1 = pointneigh(j,i)
c collect list of 2nd neighbors
	    do j1  = 1,pointorder(i1)
	      i2 = pointneigh(j1,i1)
		if(i2.ne.i)then
		  npntcurv = npntcurv + 1
		  pntcurv(npntcurv) = i2
		end if
	    end do
	  end do
	  do j = 1,npntcurv
	    do k = 1,3
	      dpoint(k) = xyzpoint(k,pntcurv(j)) - xyzpoint(k,i)
	      v1(k) = xyznorm(k,i)
	    end do
	    call dot(dpoint,v1,dotdn)
	    call dot(dpoint,dpoint,dotdd)
c	    print *,'i,dotdn,dotdd: ',i,j,dotdn,dotdd
	    curve(i) = curve(i) - 2*dotdn/dotdd
	  end do
	  curve(i) = curve(i)/npntcurv
	  crvmin = min(crvmin,curve(i))
	  crvmax = max(crvmax,curve(i))
	  crvav = crvav + curve(i)
	end do
cc 1st neighbor code
c	do i = 1,npoint
c	  curve(i) = 0.
c	  do j = 1,pointorder(i)
c	    do k = 1,3
c	      dpoint(k) = xyzpoint(k,pointneigh(j,i)) - xyzpoint(k,i)
c	      v1(k) = xyznorm(k,i)
c	    end do
c	    call dot(dpoint,v1,dotdn)
c	    call dot(dpoint,dpoint,dotdd)
cc	    print *,'i,dotdn,dotdd: ',i,j,dotdn,dotdd
c	    curve(i) = curve(i) - 2*dotdn/dotdd
c	  end do
c	  curve(i) = curve(i)/pointorder(i)
c	  crvmin = min(crvmin,curve(i))
c	  crvmax = max(crvmax,curve(i))
c	  crvav = crvav + curve(i)
c	end do
	crvav = crvav/npoint
	print *,'max,min, av curvature: ',crvmin,crvmax,crvav
c	print *,'curvature '
c	do i = 1,npoint
c	  write(6,'(i6,3f9.4)')i,curve(i)
c	end do
c
c smooth curvature by averaging
c
	print *,'smoothing curvature...'
	do i = 1,npoint
	  curve1(i) = 0.0
	  do j = 1,pointorder(i)
	    curve1(i) = curve1(i) + curve(j)
	  end do
	  curve1(i) = 0.8*curve(i) + 0.2*curve1(i)/pointorder(i)
	end do
	do i = 1,npoint
	  curve(i) = 0.0
	  do j = 1,pointorder(i)
	    curve(i) = curve(i) + curve1(j)
	  end do
	  curve(i) = 0.8*curve1(i) + 0.2*curve(i)/pointorder(i)
	end do
c
c find closest atom for each triangle, vertex
c
	goto 501 ! debug
	print *,'associating triangles/vertices with atoms...' 
	do i = 1,ntri
	  if(mod(i,100).eq.0)then
	    print *,'triangle ',i
	  endif
	  do k = 1,3
	    tri1(k,1) = xyzpoint(k,tripoint(1,i))
	    tri1(k,2) = xyzpoint(k,tripoint(2,i))
	    tri1(k,3) = xyzpoint(k,tripoint(3,i))
	    trimid(k) = (tri1(k,1) +  tri1(k,2) +  tri1(k,3))/3.
	    distmin(k) = 1.e8
	  end do
	  distmidmin = 1.e8
	  do i1 = 1,nrec
	    line = pdbrec(i1)
	    read(line,'(30x,3f8.3,f6.2)')pdbxyz,pdbrad
	    do j1 = 1,3
	      dist(j1) = (pdbxyz(1)-tri1(1,j1))**2 +  (pdbxyz(2)-tri1(2,j1))**2 + (pdbxyz(3)-tri1(3,j1))**2
	      if(dist(j1).lt.distmin(j1))then
	        distmin(j1) = dist(j1)
		pointrec(tripoint(j1,i)) = i1
	      end if
	    end do
	    distmid = (pdbxyz(1)-trimid(1))**2 +  (pdbxyz(2)-trimid(2))**2 + (pdbxyz(3)-trimid(3))**2
	    if(distmid.lt.distmidmin)then
	      distmidmin = distmid
	      trirec(i) = i1
	    end if
	  end do
	end do
501	continue
c
c radius of gyration- just a check for consistency with meshchk
c original area, volume
c
c
c find range of points 
c
        do k = 1,3
          xyzmin(k) = 1.e6
          xyzmax(k) = -1.e6
          xyzmid(k) = 0.
        end do
        do i = 1,npoint
          do k = 1,3
            xyzmid(k) = xyzmid(k) + xyzpoint(k,i)
            xyzmin(k) = min(xyzmin(k),xyzpoint(k,i))
            xyzmax(k) = max(xyzmax(k),xyzpoint(k,i))
          end do
        end do
        print *,' '
        print *,'min,max,range,mid of points'
        do k = 1,3
          xyzmid(k) = xyzmid(k)/npoint
          xyzran(k) = xyzmax(k) - xyzmin(k)
          write(6,'(4f9.4)')xyzmin(k),xyzmax(k),xyzran(k),xyzmid(k)
        end do
	rog = 0
	do i = 1,npoint
	  rog = rog + (xyzpoint(1,i)-xyzmid(1))**2 +  (xyzpoint(2,i)-xyzmid(2))**2 + (xyzpoint(3,i)-xyzmid(3))**2 
	end do
	rog = sqrt(rog/npoint)
	print *,'radius of gyration: ',rog
	Atot = 0.
	Vtot = 0.
	do i = 1,ntri
	  do k = 1,3
	    do j = 1,3
	      tri1(j,k) = xyzpoint(j,tripoint(k,i))
	    end do
	  end do
	  call triVA(tri1,xyzmid,Atri,Vtri)
	  Atot = Atot + Atri
	  Vtot = Vtot + Vtri
	end do
	print *,'Original Area, Volume: ',Atot,Vtot
	lunit1 = 11
	open(lunit1,file=tstfile)
	call putstruct(lunit1,ntrimax,nnmax,xyzpoint,npoint,trineigh,
     &   tripoint,pointneigh,pointtri,pointorder,ntri,xyznorm)
	call putrec(lunit1,ntrimax,trirec,ntri,pointrec,npoint)
	call putpdbrec(lunit1,pdbrec,natmax,nrec)
	write(lunit1,'(A)')'CURVATURE_XYZ'
	do i = 1,npoint
	  write(lunit1,'(i8,f15.6)')i,curve(i)
	end do
	write(lunit1,'(A)')'END CURVATURE_XYZ'
	write(lunit1,'(A)')'PROPERTY_XYZ'
	do i = 1,npoint
	  write(lunit1,'(i8,f15.6)')i,pproperty(i)
	end do
	write(lunit1,'(A)')'END PROPERTY_XYZ'
	close(lunit1)
c
c put out triangles in form pymol can display them
c
c	vup = 0.4
c	vlow = -0.4
c	vmid = 0.0
	vup = 1000.
	vlow = 0.0
	vmid = 500.0
	open(23,file='trigen.py')
	write(23,'(a)')'from pymol.cgo import *'
	write(23,'(a)')'from pymol import cmd'
	write(23,'(a)')'import math'
	write(23,'(a)')'obj = ['
	write(23,'(a)')'   BEGIN, TRIANGLES,'
	do i = 1,ntri
	  do k = 1,3 ! xyz and vertex index
c use original ms for display
	    tpoint(k) = xyzpoint(k,tripoint(1,i))
	    tpoint1(k) = xyzpoint(k,tripoint(2,i))
	    tpoint2(k) = xyzpoint(k,tripoint(3,i))
	    tnorm(k) = xyznorm(k,tripoint(1,i))
	    tnorm1(k) = xyznorm(k,tripoint(2,i))
	    tnorm2(k) = xyznorm(k,tripoint(3,i))
c	    val = curve(tripoint(k,i))
	    call color_blend(green,grey,white,vup,vlow,vmid,val,color(1,k))
	    val = pproperty(tripoint(k,i))
c	    call color_blend(orange,cyan,yellow,vup,vlow,vmid,val,color(1,k))
	    call color_blend(orange,violet,green,vup,vlow,vmid,val,color(1,k))
c            if(curve(tripoint(k,i)).lt.-0.2)then
c                do k1 = 1,3
c                  color(k1,k) = violet(k1)
c                end do
c            else if(curve(tripoint(j,i)).gt.0.2)then
c                do k1 = 1,3
c                  color(k1,k) = green(k1)
c                end do
c            else
c                do k1 = 1,3
c                  color(k1,k) = white(k1)
c                end do
c            end if
	  end do
          write(23,'(a,3(f8.3,'',''))')'   COLOR,',(color(k,1),k=1,3)
	  write(23,'(a,3(f8.3,'',''))')'   NORMAL,',tnorm
	  write(23,'(a,3(f8.3,'',''))')'   VERTEX,',tpoint
          write(23,'(a,3(f8.3,'',''))')'   COLOR,',(color(k,2),k=1,3)
	  write(23,'(a,3(f8.3,'',''))')'   NORMAL,',tnorm1
	  write(23,'(a,3(f8.3,'',''))')'   VERTEX,',tpoint1
          write(23,'(a,3(f8.3,'',''))')'   COLOR,',(color(k,3),k=1,3)
	  write(23,'(a,3(f8.3,'',''))')'   NORMAL,',tnorm2
	  write(23,'(a,3(f8.3,'',''))')'   VERTEX,',tpoint2
	end do
	write(23,'(a)')'   END'
	write(23,'(a)')'   ]'
	write(23,'(a)')"cmd.load_cgo(obj,'trigen')"
	close(23)
c
c put out lines in form pymol can display
c
	open(23,file='triline.py')
	write(23,'(a)')'from pymol.cgo import *'
	write(23,'(a)')'from pymol import cmd'
	write(23,'(a)')'import math'
	write(23,'(a)')'obj = ['
	write(23,'(a)')'   BEGIN, LINES,'
	do i = 1,npoint
	  do j = 1,pointorder(i)
	    i1 = pointneigh(j,i)
	    if(i1.gt.1)then
	      do k = 1,3 
	        tpoint(k) = xyzpoint(k,i)
	        tpoint1(k) = xyzpoint(k,i1)
	      end do
	      val = curve(i)
	      call color_blend(green,grey,white,vup,vlow,vmid,val,color(1,1))
              write(23,'(a,3(f8.3,'',''))')'   COLOR,',(color(k,1),k=1,3)
	      write(23,'(a,3(f8.3,'',''))')'   VERTEX,',tpoint
	      write(23,'(a,3(f8.3,'',''))')'   VERTEX,',tpoint1
	    end if
	  end do
	end do
	write(23,'(a)')'   END'
	write(23,'(a)')'   ]'
	write(23,'(a)')"cmd.load_cgo(obj,'triline')"
	close(23)
	stop
c----------------------------------------------------
900	print *,'cant find triangle file: ',trifile
	stop
901	print *,'error reading triangle file at line: ',ntri
	stop
	end


	subroutine norm(xyz,rnorm)
	implicit none
	real*4 xyz(3),rnorm
	rnorm = xyz(1)**2 + xyz(2)**2 + xyz(3)**2
	if(rnorm.gt.1.e-6)then
	rnorm = sqrt(rnorm)
	xyz(1) = xyz(1) / rnorm
	xyz(2) = xyz(2) / rnorm
	xyz(3) = xyz(3) / rnorm
	end if
	return
	end 


	subroutine triVA(tri,xyzmid,Atri,Vtri)
	real*4 tri(3,3),xyzmid(3)
	real*4 Atri,Vtri
	real*4 v1(3),v2(3),v3(3),rcen(3)
	do k = 1,3
	  v1(k) = tri(k,3) - tri(k,1)
	  v2(k) = tri(k,2) - tri(k,1)
	  rcen(k) = tri(k,1) - xyzmid(1)
	end do
	call cross(v1,v2,v3)
	call dot(v3,rcen,Vtri)
	call norm(v3,Atri)
	Vtri = Vtri/6.
	Atri = Atri/2.
	return
	end
c------------------------------------------------------
      SUBROUTINE CROSS(U,V,W)
      dimension u(3),v(3),w(3)
      w(1) = u(2)*v(3) - u(3)*v(2)
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)
      return
      end
c------------------------------------------------------
      SUBROUTINE DOT(U,V,RDOT)
      dimension u(3),v(3)
      rdot = 0.0
        rdot = u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
      return
      end
c------------------------------------------------------
