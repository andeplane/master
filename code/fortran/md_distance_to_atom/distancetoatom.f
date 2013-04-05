      program distancetoatom
c     based on the smarter algorithm to mark all parts of an array 
c     within a given distance of an Atom - this is also more flexible
c     in allowing different distances for different atoms

c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)

c     Data for storing the set
      integer maxpart
      parameter(maxpart=1000000)
      integer npart ! actual nr of particles
      real*8 coord(maxpart,3) ! point positions

c     Data for matrix
      integer nnmax
      parameter(nnmax=310) ! max set is nnmax cuber
      real*8 atomdist(nnmax,nnmax,nnmax)
      real*8 matboxsize ! size of matrix boxes

c     Data for analysis
      real*8 partdist ! distance to point, also used for boxlist

      character*50 infile
      character*5 tmpstr
      real*4 rout1

c-----------------------------------------------------------------------

c     Read in parameters
      write(*,*) 'Filename for data'
      read(*,*) infile
      write(*,*) 'Reading points from ',infile
      write(*,*) 'Matrix pixel length (in Angstroms)'
      read(*,*) matboxsize
      write(*,*) 'Matrix pixel length = ',matboxsize,' (Angstroms)'
      write(*,*) 'Matrix dimensions'
      read(*,*) ndx,ndy,ndz
      if ((ndx.lt.1).or.(ndx.gt.nnmax)) then
         write(*,*) 'Nx out of range ',ndx
         stop
      end if
      if ((ndy.lt.1).or.(ndy.gt.nnmax)) then
         write(*,*) 'Ny out of range ',ndy
         stop
      end if
      if ((ndz.lt.1).or.(ndz.gt.nnmax)) then
         write(*,*) 'Nz out of range ',ndz
         stop
      end if
      write(*,*) 'Notice that points outside matrix will be ignored'
      write(*,*) 'Max distance around atom for analysis'
      read(*,*) partdist
      write(*,*) 'Max distance around atom (in Angstroms) = ',partdist
c     Convert partdist to box size in atomdist
      partdist = partdist/matboxsize

      read(*,*) nperiod
      if (nperiod.le.0) then
         write(*,*) 'Not period boundaries'
      else
         write(*,*) 'Periodic boundaries'
      end if

      xsize = ndx*matboxsize
      xsize2 = xsize/2
      ysize = ndy*matboxsize
      ysize2 = ysize/2
      zsize = ndz*matboxsize
      zsize2 = zsize/2

c     Reset matrix
      do ix = 1,ndx
         do iy = 1,ndy
            do iz = 1,ndz
               atomdist(ix,iy,iz) = xsize2
            end do
         end do
      end do

c     Read points from file
      xmin = 1e10
      xmax = -xmin
      ymin = 1e10
      ymax = -ymin
      zmin = 1e10
      zmax = -zmin
      open(unit=10,file=infile)
      read(10,*) npart
      write(*,*) 'Nr of particles = ',npart
      read(10,*) tmpstr
      do i = 1,npart
         read(10,*) tmpstr,x,y,z,id
         coord(i,1) = x
         coord(i,2) = y
         coord(i,3) = z
         if (x.gt.xmax) xmax = x
         if (x.lt.xmin) xmin = x
         if (y.gt.ymax) ymax = y
         if (y.lt.ymin) ymin = y
         if (z.gt.zmax) zmax = z
         if (z.lt.zmin) zmin = z
      end do
      close(10)
      write(*,*) 'xmin,xmax = ',xmin,xmax
      write(*,*) 'ymin,ymax = ',ymin,ymax
      write(*,*) 'zmin,zmax = ',zmin,zmax

c     Loop through all atoms and mark them in the matrix
      nstep = 100
      rad = partdist/matboxsize
      nrad = rad
      nrad = nrad + 1
      do i = 1,npart
         if (mod(i,nstep).eq.0) then
            write(*,*) 'Atom nr ',i
         end if
         x = coord(i,1)
         y = coord(i,2)
         z = coord(i,3)
         cx = x/matboxsize
         nx = cx
         nx = nx + 1
         cy = y/matboxsize
         ny = cy
         ny = ny + 1
         cz = z/matboxsize
         nz = cz
         nz = nz + 1
         do ix = -nrad,nrad
            do iy = -nrad,nrad
               do iz = -nrad,nrad
                  nnx = nx + ix
                  nny = ny + iy
                  nnz = nz + iz
                  xx = ((nnx-1)+0.5)*matboxsize 
                  yy = ((nny-1)+0.5)*matboxsize 
                  zz = ((nnz-1)+0.5)*matboxsize 
                  dx = xx - x
                  dy = yy - y
                  dz = zz - z
                  if (nperiodic.gt.0) then
                     if (dx.gt.xsize2) dx = dx - xsize
                     if (dx.lt.-xsize2) dx = dx + xsize
                     if (dy.gt.ysize2) dy = dy - ysize
                     if (dy.lt.-ysize2) dy = dy + ysize
                     if (dz.gt.zsize2) dz = dz - zsize
                     if (dz.lt.-zsize2) dz = dz + zsize
                  end if
                  rr = sqrt(dx*dx + dy*dy + dz*dz)
                  if (rr.lt.partdist) then
c     Inside this radius - check that it is not outside the matrix
                     if ((nnx.ge.1).and.(nnx.le.ndx).and.
     &                    (nny.ge.1).and.(nny.le.ndy).and.
     &                    (nnz.ge.1).and.(nnz.le.ndz)) then
                        rrr = atomdist(nnx,nny,nnz)
                        if (rr.lt.rrr) then
                           atomdist(nnx,nny,nnz) = rr
                        end if
                     end if
                  end if
               end do
            end do
         end do
      end do

      rmin = xsize2
      rmax = -rmin
c     Output value of matrix to vtk voxel file
      open(unit=20,file='tmp.vtk')
      write(20,90) '# vtk DataFile Version 2.0'
 90   format(a)
 91   format(a,3i8)
 92   format(a,i8)
      write(20,90) 'structured_point'
      write(20,90) 'ASCII'
      write(20,*)
      write(20,90) 'DATASET STRUCTURED_POINTS'
      write(20,91) 'DIMENSIONS',ndx,ndy,ndz
      write(20,90) 'ORIGIN 0.0 0.0 0.0'
      write(20,90) 'SPACING 1 1 1'
      nn = ndx*ndy*ndz
      write(20,92) 'POINT_DATA ',nn
      write(20,90) 'SCALARS atomdist double'
      write(20,90) 'LOOKUP_TABLE default'
      write(20,*)
      do ix = 1,ndx
         do iy = 1,ndy
            do iz = 1,ndz
               rr = atomdist(ix,iy,iz)
               rout1 = rr
               write(20,*) rout1
               if (rr.lt.rmin) rmin = rr
               if (rr.gt.rmax) rmax = rr
            end do
         end do
      end do
      close(20)
      write(*,*) 'Min dist = ',rmin
      write(*,*) 'Max dist = ',rmax
      
      end

