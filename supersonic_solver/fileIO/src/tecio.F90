!-----------------------------------------------------------------------------------------------!
! Filename   : tecio.f90                                                                        !
! Author     : R. Ranjan                                                                        !
!              School of Aerospace Engineering                                                  !
!              Georgia Institute of Technology                                                  !
!              Atlanta, Georgia, US                                                             !
!              3/15/2017                                                                        !
! Description: Subroutine to write data in ASCII block format for Tecplot                       !
!              - WriteTecVar                                                                    !
!-----------------------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------------------!
! WriteTecVar                                                                                   !
! Arguments   : flunit, i1, i2, j1, j2, qvar                                                    !
!               flunit    integer    File id to which the data is written                       !
!               i1        integer    Start index of the data in x direction                     !
!               i2        integer    End index of the data in x direction                       !
!               j1        integer    Start index of the data in y direction                     !
!               j2        integer    End index of the data in x direction                       !
!               qvar      double     2D array having the data that need to be written           !
! Description : The subroutine writes the 2D data in ASCII format to the file ID passed to this !
!               subroutine.                                                                     !
!-----------------------------------------------------------------------------------------------!  
subroutine WriteTecVar(flunit, i1, i2, j1, j2, qvar)
  implicit none

  integer, intent(in) :: flunit
  integer, intent(in) :: i1, i2, j1, j2
  double precision, intent(in), dimension(i1:i2,j1:j2) :: qvar
  integer :: i, j, k

  do j = j1, j2
     do i = i1, i2
        write(flunit, *), qvar(i,j)
     enddo
  enddo
  
  return
end subroutine WriteTecVar
!-----------------------------------------------------------------------------------------------!  
! End of File                                                                                   !
!-----------------------------------------------------------------------------------------------!  
