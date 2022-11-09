program mdft_plugin
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  read_wf = .True.
  touch read_wf 
  call routine_grid_pts_mdft
  call routine_mos_grid_mdft
  call routine_den_grid_mdft
  call routine_dipolar_moment_mdft
  call routine_pot_mm_in_aos
  call routine_pot_mm_in_mos
  call routine_ao_overlap
!  call routine_ao_to_mo
end

subroutine routine_grid_pts_mdft
  BEGIN_DOC
  ! Print the points grid of MDFT in a file.dat 
  END_DOC
  implicit none
  integer :: i,j
  open(10, file="Point_grid_mdft.dat") 
  do i = 1,n_pts_mdft_tot
   write(10,'(100(F10.5,X))')grid_pts_mdft_xyz(:,i)
  enddo
  close(10)
end

subroutine routine_mos_grid_mdft
  BEGIN_DOC
  ! Print the MO grid of MDFT in a file.dat
  END_DOC
  implicit none
  integer :: i,j
  open(11, file="MOs_grid_mdft.dat")
  do i = 1, n_pts_mdft_tot
    write(11,'(100(F10.5,X))')mos_grid_mdft(:,i)
  enddo
  close(11)
end

subroutine routine_den_grid_mdft
  BEGIN_DOC
  ! Print the point grid of MDFT in a file.dat
  END_DOC
  implicit none
  integer :: i,j
  open(12, file="Density_grid_mdft.dat")
  do i = 1, n_pts_mdft_tot
    write(12,'(100(F10.5,X))')density_grid_mdft(:,i)
  enddo
  close(12)
end

subroutine routine_dipolar_moment_mdft
  BEGIN_DOC
  ! Print the dipolar moment of the molecule in the x/y/z direction in a file.dat 
  END_DOC
  implicit none
  integer :: i,j
  open(13, file="Dipolar_moment_mdft.dat")
  write(13,'(100(F10.5,X))') dipolar_moment_xyz(1), dipolar_moment_xyz(2), dipolar_moment_xyz(3)
  close(13)
end

subroutine routine_pot_mm_in_aos
  BEGIN_DOC
  ! Print the potential from MDFT program in ao basis in a file.dat
  END_DOC
  implicit none
  integer :: i,j
  open(14, file="pot_mm_in_aos.dat")
  do i = 1, ao_num
    write(14,'(100(F10.5,X))') pot_mm_in_aos(:,i)
  enddo
  close(14)
end

subroutine routine_ao_overlap
  BEGIN_DOC
  ! Print the potential from MDFT program in ao base in a file.dat
  END_DOC
  implicit none
  integer :: i,j
  open(20, file="ao_overlap.dat")
  do i = 1, ao_num
    write(20,'(100(F10.5,X))') pot_mm_in_aos(:,i)
    write(20,'(100(F10.5,X))') ao_overlap(:,i)
  enddo
  close(20)
end

subroutine routine_pot_mm_in_mos
  BEGIN_DOC
  ! Print the potential from MDFT program in mo base in a file.dat
  END_DOC
  implicit none
  integer :: i,j
  open(15, file="pot_mm_in_mos.dat")
  do i = 1, mo_num
    write(15,'(100(F10.5,X))') pot_mm_in_mos(:,i)
  enddo
  close(15)
end

subroutine routine_ao_to_mo
  BEGIN_DOC
  ! Print the potential from MDFT program in ao base in a file.dat
  END_DOC
  implicit none
  integer :: i,j,k
  open(21, file="ao_to_mo.dat")
!  call ao_to_mo()
!  do i = 1,
!    write(21,'(100(F10.5,X))') 
!  enddo
  close(21)
end

