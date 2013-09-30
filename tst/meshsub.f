c
c subroutines for various mesh programs
c
	subroutine getpdbrec(lunit,pdbrec,natmax,nrec)
c read pdb records from a meshsrf triangle file until end encoutered
	implicit none
	integer nrec,natmax,lunit
	character*80 pdbrec(natmax),line
c------------------------------------------------------
	nrec = 0
100	read(lunit,'(A)',end=900)line
	if(line(1:3).eq.'END')then
	  print *,'pdb records read: ',nrec
	  return
	else
	  nrec = nrec + 1
	  if(nrec.gt.natmax)then
	    print *,'exceeded atom record max- increase natmax: ',natmax
	    stop
	  end if
	  pdbrec(nrec) = line
	end if
	goto 100
900	print *,'unexpected end of PDB records at # ',nrec
	close(lunit)
	return
	end

	subroutine putpdbrec(lunit,pdbrec,natmax,nrec)
c put pdb records intoa triangle structure file
	implicit none
	integer nrec,natmax,lunit
	character*80 pdbrec(natmax),line
	integer i
c------------------------------------------------------
	write(lunit,'(A)')'PDB_RECORD'
	do i = 1,nrec
	  write(lunit,'(A)')pdbrec(i)
	end do
	write(lunit,'(A)')'END PDB_RECORD'
	return
	end

	subroutine gettrixyz(lunit,triangle,ntrimax,ntri)
c read triangle coords from a meshsrf triangle file until end encoutered
	implicit none
	integer ntrimax,ntri,lunit,nline,j,k
	character*80 line
	real*4 triangle(4,3,ntrimax)
	real*4 tpoint(3),property
c--------------------------------------------
	ntri = 0 
	nline = 0
	j = 3
100	read(lunit,'(A)',end=900)line
	if(line(1:3).eq.'END')then
	  print *,'triangle records read: ',nline
          print *,'# of triangles: ',ntri
          if(j.ne.3)then
            print *,'warning not 3N lines in file!'
          end if
	  return
	else
	  nline = nline + 1
c try reading 4th column or not
	  read(line,'(1x,3f9.3,f15.6)',err=901)tpoint,property
	  goto 902
901	  read(line,'(1x,3f9.3)',err=900)tpoint
	  property = 0.
902	  continue
          j = j + 1
          if(j.eq.4)then
            j = 1
            ntri = ntri + 1
	    if(ntri.gt.ntrimax)then
	      print *,'exceeded triangle record max- increase ntrimax: ',ntrimax
	      stop
	    end if
	  end if
          do k = 1,3
            triangle(k,j,ntri) = tpoint(k)
          end do
          triangle(4,j,ntri) = property
	end if
        goto 100
900	print *,'unexpected end of triangle xyz records at # ',nline
	close(lunit)
	return
	end

	subroutine putstruct(lunit,ntrimax,nnmax,xyzpoint,npoint,trineigh,
     &   tripoint,pointneigh,pointtri,pointorder,ntri,xyznorm)
c
c put all neighbor lists into file
c
	implicit none
	integer ntrimax,nnmax
	integer ntri,npoint
	integer trineigh(3,ntrimax)
	integer tripoint(3,ntrimax)
	integer pointneigh(nnmax,ntrimax)
	integer pointtri(nnmax,ntrimax)
	integer pointorder(ntrimax)
	real*4  xyzpoint(3,ntrimax)
	real*4  xyznorm(3,ntrimax)
	integer i,k,lunit
c---------------------------------------------------
	write(lunit,'(A)')'TRIANGLE_NEIGHBOR LIST (CLOCKWISE)'
	do i = 1,ntri
	  write(lunit,'(4i8)')i,(trineigh(k,i),k=1,3)
	end do
	write(lunit,'(A)')'END TRIANGLE_NEIGHBOR LIST'
	write(lunit,'(A)')'POINT_XYZ LIST'
	do i = 1,npoint
	  write(lunit,'(i8,3f9.4,f15.6)')i,(xyzpoint(k,i),k=1,3)
	end do
	write(lunit,'(A)')'END POINT_XYZ LIST'
	write(lunit,'(A)')'TRIANGLE_POINT LIST (CLOCKWISE really, this was checked)'
	do i = 1,ntri
	  write(lunit,'(i8,3i8)')i,(tripoint(k,i),k=1,3)
	end do
	write(lunit,'(A)')'END TRIANGLE_POINT LIST'
	write(lunit,'(A)')'POINT_TRIANGLE LIST (CLOCKWISE)'
	do i = 1,npoint
	  write(lunit,'(i8,13i8)')i,pointorder(i),(pointtri(k,i),k=1,pointorder(i))
	end do
	write(lunit,'(A)')'END POINT_TRIANGLE LIST'
	write(lunit,'(A)')'POINT_NEIGHBOR LIST (CLOCKWISE)'
	do i = 1,npoint
	  write(lunit,'(i8,13i8)')i,pointorder(i),(pointneigh(k,i),k=1,pointorder(i))
	end do
	write(lunit,'(A)')'END POINT_NEIGHBOR LIST (CLOCKWISE)'
	write(lunit,'(A)')'NORM_XYZ LIST'
	do i = 1,npoint
	  write(lunit,'(i8,3f9.4,f15.6)')i,(xyznorm(k,i),k=1,3)
	end do
	write(lunit,'(A)')'END NORM_XYZ LIST'
	return
	end

	subroutine getstruct(lunit,ntrimax,nnmax,xyzpoint,npoint,trineigh,
     &   tripoint,pointneigh,pointtri,pointorder,ntri,xyznorm)
c
c gets neighbor lists from file- note this starts after
c findind a record: 
c	TRIANGLE_NEIGHBOR LIST (CLOCKWISE)
c
	implicit none
	integer ntrimax,nnmax
	integer ntri,npoint
	integer trineigh(3,ntrimax)
	integer tripoint(3,ntrimax)
	integer pointneigh(nnmax,ntrimax)
	integer pointtri(nnmax,ntrimax)
	integer pointorder(ntrimax)
	real*4  xyzpoint(3,ntrimax)
	real*4  xyznorm(3,ntrimax)
	integer i,j,k,lunit
	character*132 line
	integer ntri1,npoint1
c---------------------------------------------------
	read(lunit,'(A)',end=900)line
	ntri = 0
	do while(line(1:3).ne.'END')
	  ntri = ntri + 1
	  read(line,'(4i8)',err=900)i,(trineigh(k,ntri),k=1,3)
	  read(lunit,'(A)',end=900)line
	end do
	ntri1 = ntri
	print *,'TRIANGLE_NEIGHBOR LIST records: ',ntri
c
	read(lunit,'(A)',end=900)line
	print *,line
	if(line(1:14).ne.'POINT_XYZ LIST')goto 900
	read(lunit,'(A)',end=900)line
	npoint = 0
	do while(line(1:3).ne.'END')
	  npoint = npoint + 1
	  read(line,'(i8,3f9.4,f15.6)',err=900)i,(xyzpoint(k,npoint),k=1,3)
	  read(lunit,'(A)',end=900)line
	end do
	npoint1 = npoint
	print *,'POINT_XYZ LIST records: ',npoint
c
	read(lunit,'(A)',end=900)line
	print *,line
	if(line(1:19).ne.'TRIANGLE_POINT LIST')goto 900
	read(lunit,'(A)',end=900)line
	ntri = 0
	do while(line(1:3).ne.'END')
	  ntri = ntri + 1
	  read(line,'(4i8)',err=900)i,(tripoint(k,ntri),k=1,3)
	  read(lunit,'(A)',end=900)line
	end do
	print *,'TRIANGLE_POINT LIST records: ',ntri
	if(ntri.ne.ntri1)then
	  print *,'warning, different triangle #s:',ntri1,ntri
	end if
c
	read(lunit,'(A)',end=900)line
	print *,line
	if(line(1:14).ne.'POINT_TRIANGLE')goto 900
	read(lunit,'(A)',end=900)line
	i = 0
	do while(line(1:3).ne.'END')
	  i = i + 1
	  read(line,'(i8,13i8)')j,pointorder(i),(pointtri(k,i),k=1,pointorder(i))
cc debug
c	  if(pointorder(i).ge.9)then
c	    print *,line
c	    write(6,'(2i4,13i8)')i,pointorder(i),(pointtri(k,i),k=1,pointorder(i))
c	  end if
c debug
	  read(lunit,'(A)',end=900)line
	end do
	print *,'POINT_TRIANGLE LIST records: ',i
c
	read(lunit,'(A)',end=900)line
	print *,line
	if(line(1:11).ne.'POINT_NEIGH')goto 900
	read(lunit,'(A)',end=900)line
	i = 0
	do while(line(1:3).ne.'END')
	  i = i + 1
	  read(line,'(i8,13i8)')j,pointorder(i),(pointneigh(k,i),k=1,pointorder(i))
cc debug
c	  if(pointorder(i).ge.9)then
c	    print *,line
c	    write(6,'(2i4,13i8)')i,pointorder(i),(pointneigh(k,i),k=1,pointorder(i))
c	  end if
cc debug
	  read(lunit,'(A)',end=900)line
	end do
	print *,'POINT_NEIGH LIST records: ',i
	read(lunit,'(A)',end=900)line
	print *,line
	if(line(1:14).ne.'NORM_XYZ LIST')goto 900
	read(lunit,'(A)',end=900)line
	npoint = 0
	do while(line(1:3).ne.'END')
	  npoint = npoint + 1
	  read(line,'(i8,3f9.4,f15.6)',err=900)i,(xyznorm(k,npoint),k=1,3)
	  read(lunit,'(A)',end=900)line
	end do
	print *,'NORM_XYZ LIST records: ',npoint

	return
900	print *,'premature error or end of triangle structure records: ',line
	end

	subroutine putrec(lunit,ntrimax,trirec,ntri,pointrec,npoint)
c
c put links to atoms for triangles and points
c
	implicit none
	integer ntrimax
	integer ntri,npoint
	integer trirec(ntrimax),pointrec(ntrimax)
	integer i,k,lunit
c---------------------------------------------------
	write(lunit,'(A)')'POINT_PDB_RECORD'
	do i = 1,npoint
	  write(lunit,'(2i8)')i,pointrec(i)
	end do
	write(lunit,'(A)')'END POINT_PDB_RECORD'
	write(lunit,'(A)')'TRIANGLE_PDB_RECORD'
	do i = 1,ntri
	  write(lunit,'(2i8)')i,trirec(i)
	end do
	write(lunit,'(A)')'END TRIANGLE_PDB_RECORD'
	return
	end

	subroutine getrec(lunit,ntrimax,trirec,ntri,pointrec,npoint)
c
c get links to atoms for triangles and points
c
	implicit none
	integer ntrimax
	integer ntri,npoint
	integer trirec(ntrimax),pointrec(ntrimax)
	integer i,k,lunit
	character*132 line
c---------------------------------------------------
	npoint = 0
	print *,'reading point/triangle to atom mapping...'
	read(lunit,'(A)',end=900)line
	do while(line(1:3).ne.'END')
	  npoint = npoint + 1
	  read(line,'(2i8)',err=900)i,pointrec(npoint)
	  read(lunit,'(A)',end=900)line
	end do
	write(6,'(A)')line
	read(lunit,'(A)',end=900)line
	write(6,'(A)')line
	ntri = 0
	read(lunit,'(A)',end=900)line
	do while(line(1:3).ne.'END')
	  ntri = ntri + 1
	  read(line,'(2i8)',err=900)i,trirec(ntri)
	  read(lunit,'(A)',end=900)line
	end do
	write(6,'(A)')line
	print *,'point,triangle mappings read: ',npoint,ntri
	return
900	print *,'error reading point/triangle mappings or early end of file: ',npoint,ntri
	return
	end


	subroutine dostack(ntrimax,nnmax,pointneigh,pointorder,npoint,nsubset,isubset)
	implicit none
        integer ntrimax
        integer nnmax
c      
c list of each point's neighboring points
c      
	integer npoint
        integer pointneigh(nnmax,ntrimax)
        integer pointorder(ntrimax)
	integer isubset(ntrimax)
c
c stack and list for running through each triangle once
c
        integer idone(ntrimax),istack(ntrimax),nstack
        integer icurr,inext,ndone,i,k,nsubset
c----------------------------------------------------------------
	ndone = 0
	nsubset = 1
	do i = 1,npoint
	  idone(i) = 0
	end do
100	i = 1
	do while((idone(i).eq.1).and.(i.le.npoint))
	  i = i + 1
	end do
	nstack = 1
	istack(nstack) = i
c
c run thru points
c
	do while(nstack.gt.0)
c
c pop next point off stack
c
	  icurr = istack(nstack)
	  nstack = nstack - 1
	  if(idone(icurr).eq.0)then
c
c if new point
c
	    idone(icurr) = 1
	    isubset(icurr) = nsubset
	    ndone = ndone + 1
c	    print *,'doing point: ',icurr
c
c push its neighbors onto stack if they haven't been done
c
	    do k = 1,pointorder(icurr)
	      inext = pointneigh(k,icurr)
	      if(idone(inext).eq.0)then
	        nstack = nstack + 1
	        istack(nstack) = inext
	      end if
	    end do
	  end if
	end do
	print *,'# done, # points; ',ndone,npoint
	if(ndone.lt.npoint)then
	  nsubset = nsubset + 1
	  goto 100
	end if
	print *,' # of point subsets: ',nsubset,' > one means cavities'
	return
	end

        subroutine getproperty(lunit,property,ntrimax,npoint)
	implicit none
	integer ntrimax
	integer npoint
	real*4 property(ntrimax)
	integer i,k,lunit
	character*132 line
c-------------------------------------------------------
	npoint = 0
	print *,'reading point property'
	read(lunit,'(A)',end=900)line
	do while(line(1:3).ne.'END')
	  npoint = npoint + 1
	  read(line,'(i8,f15.6)',err=900)i,property(npoint)
	  read(lunit,'(A)',end=900)line
	end do
	write(6,'(A)')line
	print *,'point properties read: ',npoint
	return
900	print *,'error reading point properties or early end of file: ',npoint
	return
	end

        subroutine maptotri(tripoint,ntrimax,property,npoint,propertyt,ntri)
	implicit none
	integer ntrimax
	integer npoint,ntri
	integer tripoint(3,ntrimax)
	real*4 property(ntrimax),propertyt(ntrimax)
	integer i,j,k
c----------------------------------------
	do i = 1,ntri
	  propertyt(i) = 0.
	  do j = 1,3
	    propertyt(i) = propertyt(i) + property(tripoint(j,i))/3.
	  end do
	end do
	return
	end

	subroutine color_blend(colup,collow,colmid,vup,vlow,vmid,val,col)

	real*4 colup(3),collow(3),colmid(3),col(3),vup,vlow,vmid,val
	real*4 dval
	integer i,j,k
c---------------------------------------------------------------
	if(val.le.vlow)then
	  do k = 1,3
	    col(k) = collow(k)
	  end do
	else if(val.ge.vup)then
	  do k = 1,3
	    col(k) = colup(k)
	  end do
	else
c	  do k = 1,3
c	    col(k) = colmid(k)
c	  end do
	  if(val.le.vmid)then
	    dval  = (vmid - val)/(vmid - vlow)
	    do k = 1,3
	      col(k) = dval*collow(k) + (1.-dval)*colmid(k)
	    end do
	  else
	    dval  = (val - vmid)/(vup - vmid)
	    do k = 1,3
	      col(k) = dval*colup(k) + (1.-dval)*colmid(k)
	    end do
	  end if
	end if
c	print *,'in color blend colup : ',colup
c	print *,'in color blend collow : ',collow
c	print *,'in color blend colmid : ',colmid
c	print *,'in color blend: ',vup,vlow,vmid,val,dval,col



	return
	end
