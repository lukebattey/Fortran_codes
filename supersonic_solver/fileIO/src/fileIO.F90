!-----------------------------------------------------------------------------------------------!
!                     ASCII Format Tecplot File IO Program                                      !
!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
! Filename   : Driver.f90                                                                       !
! Author     : R. Ranjan                                                                        !
!              School of Aerospace Engineering                                                  !
!              Georgia Institute of Technology                                                  !
!              Atlanta, Georgia, US                                                             !
!              3/15/2017                                                                        !
! Description: This is the driver program to define a function in 2D over a rectangular grid    !
!              and save the data in ASCII format, which can be loaded in Tecplot software.      ! 
!-----------------------------------------------------------------------------------------------!
program fileIO
  implicit none

  double precision, allocatable, dimension(:,:) :: x, y                         ! 2D array to store grid
  double precision, allocatable, dimension(:,:) :: q                            ! 2D array to store function
  double precision                              :: xmin, xmax, ymin, ymax       ! Domain extents in x and y directions
  double precision                              :: dx, dy                       ! Grid size in x and y directions
  double precision                              :: time                         ! Time of the solution 
  double precision, parameter                   :: pi = 4.0d0*atan(1.0d0)       ! Value of constant pi
  character (100)                               :: fname                        ! Name of the file to which data uis saved
  integer                                       :: imax, jmax                   ! Number of points along x and y directions
  integer                                       :: i, j                         ! Loop variables
  integer                                       :: i1, i2, j1, j2               ! Internal variables
  integer                                       :: stID                         ! Strand ID for tecplot data to understand it is a time varying data
  integer                                       :: lout                         ! File IO unit

  ! Set domain size
  print *, 'Setting up the problem...'
  xmin = 0.0d0
  xmax = 6.0d0
  ymin = 0.0d0
  ymax = 2.0d0

  ! Set grid size
  imax = 101    
  jmax = 51      

  dx = (xmax - xmin)/(imax - 1.0d0)
  dy = (ymax - ymin)/(jmax - 1.0d0)

  ! Allocate arrays
  print *, 'Allocating dynamic memory...'

  allocate(x(1:imax,1:jmax))
  allocate(y(1:imax,1:jmax))
  allocate(q(1:imax,1:jmax))

  ! Compute the grid
  print *, 'Generating the grid...'
  x(1,:) = xmin
  do i=2,imax
     x(i,:) = x(i-1,:) + dx
  enddo

  y(:,1) = ymin
  do j=2,jmax
     y(:,j) = y(:,j-1) + dy
  enddo

  ! Set the function
  print *, 'Computing the function...'
  do j=1,jmax
     do i=1,imax
        q(i,j) = sin(2*pi*x(i,j)) * cos(2*pi*y(i,j))
     enddo
  enddo

  ! Save ascii Tecplot file
  print *, 'Saving ASCII tecplot file...'
  lout = 151
  time = 0.0d0
  stID = 1
  i1 = 1
  i2 = imax
  j1 = 1
  j2 = jmax

  ! Open file
  write(fname,'("./data/fileIO.dat")')
  open(unit = lout, file=trim(fname),status='unknown', form='formatted')
  
  ! Write the header
  write(lout,'(A30)') 'VARIABLES = "X" "Y" "Q"'
  write(lout,'("ZONE SolutionTime= ",F13.8, ", StrandId=",I5,", I=", I5, ", J= ", I5, &
       ", ZONETYPE=Ordered DATAPACKING=BLOCK")')  time, stID, (i2-i1+1), (j2-j1+1)

  ! Write the variables
  call WriteTecVar(lout, i1, i2, j1, j2, x(i1:i2,j1:j2)) 
  call WriteTecVar(lout, i1, i2, j1, j2, y(i1:i2,j1:j2)) 
  call WriteTecVar(lout, i1, i2, j1, j2, q(i1:i2,j1:j2)) 

  ! Close file
  close(lout)

  print *, 'Finished saving tecplot file...'

  ! Deallocate arrays
  print *, 'Performing cleanup...'
  deallocate(x)
  deallocate(y)
  deallocate(q)

  print *, 'DONE.'
  
  stop
end program fileIO
!-----------------------------------------------------------------------------------------------!
! End of File                                                                                   !
!-----------------------------------------------------------------------------------------------!
