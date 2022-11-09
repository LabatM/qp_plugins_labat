BEGIN_PROVIDER [ integer, n_pts_mdft_x ]
 implicit none
 BEGIN_DOC
 ! Number of grid points for MDFT in the x direction
 END_DOC
 n_pts_mdft_x = 40
END_PROVIDER

BEGIN_PROVIDER [ integer, n_pts_mdft_y ]
 implicit none
 BEGIN_DOC
 ! Number of grid points for MDFT in the y direction
 END_DOC
 n_pts_mdft_y = 40
END_PROVIDER

BEGIN_PROVIDER [ integer, n_pts_mdft_z ]
 implicit none
 BEGIN_DOC 
 ! Number of grid points for MDFT in the z direction
 END_DOC
 n_pts_mdft_z = 40
END_PROVIDER 

BEGIN_PROVIDER [ integer, n_pts_mdft_xyz, (3) ]
 implicit none
 BEGIN_DOC
 ! Number of grid points for MDFT in the x/y/z directions
 END_DOC
 n_pts_mdft_xyz (1) = n_pts_mdft_x
 n_pts_mdft_xyz (2) = n_pts_mdft_y
 n_pts_mdft_xyz (3) = n_pts_mdft_z
END_PROVIDER

BEGIN_PROVIDER [ integer, n_pts_mdft_tot ]
 implicit none
 BEGIN_DOC
 ! Total number of grid points for MDFT
 END_DOC
 n_pts_mdft_tot = (n_pts_mdft_x+1)*(n_pts_mdft_y+1)*(n_pts_mdft_z+1)
END_PROVIDER

BEGIN_PROVIDER [ double precision, grid_length_mdft_x ]
 implicit none
  BEGIN_DOC
  ! Grid length for MDFT in the x direction
  END_DOC
  grid_length_mdft_x = 40.d0
END_PROVIDER

BEGIN_PROVIDER [ double precision, grid_length_mdft_y ]
 implicit none
  BEGIN_DOC
  ! Grid length for MDFT in the y direction
  END_DOC
  grid_length_mdft_y = 40.d0
END_PROVIDER

BEGIN_PROVIDER [ double precision, grid_length_mdft_z ]
 implicit none
  BEGIN_DOC
  ! Grid length for MDFT in the z direction
  END_DOC
  grid_length_mdft_z = 40.d0
END_PROVIDER

BEGIN_PROVIDER [ double precision, grid_length_mdft_xyz, (3) ]
 implicit none
  BEGIN_DOC
  ! Grid length for MDFT in the x/y/z directions
  END_DOC
  grid_length_mdft_xyz(1) = grid_length_mdft_x
  grid_length_mdft_xyz(2) = grid_length_mdft_y
  grid_length_mdft_xyz(3) = grid_length_mdft_z
END_PROVIDER

 BEGIN_PROVIDER [ double precision, grid_step_mdft_xyz, (3) ]
&BEGIN_PROVIDER [ double precision, grid_volume_mdft ]
 implicit none
  BEGIN_DOC
  ! Grid steps for MDFT in the x/y/z directions and the volume of the grid 
  END_DOC
  integer :: i
  do i=1,3
  if (n_pts_mdft_xyz(i).eq.0) then
   grid_step_mdft_xyz(i) = 0.d0
  else
   grid_step_mdft_xyz(i) = grid_length_mdft_xyz(i)/n_pts_mdft_xyz(i)
  end if
  enddo
  grid_volume_mdft = grid_step_mdft_xyz(1) * grid_step_mdft_xyz(2) * grid_step_mdft_xyz(3)
END_PROVIDER

BEGIN_PROVIDER [ double precision, grid_pts_mdft_x, (n_pts_mdft_x+1) ] 
 implicit none
  BEGIN_DOC
  ! Grid points for MDFT in the x direction
  END_DOC
  integer :: i
  do i = 1, n_pts_mdft_x+1
   grid_pts_mdft_x(i) = -grid_length_mdft_xyz(1)*0.5+(i-1)*grid_step_mdft_xyz(1)
  enddo
END_PROVIDER
BEGIN_PROVIDER [ double precision, grid_pts_mdft_y, (n_pts_mdft_y+1) ]
 implicit none
  BEGIN_DOC
  ! Grid points for MDFT in the y direction
  END_DOC
  integer :: i
  do i = 1, n_pts_mdft_y+1
   grid_pts_mdft_y(i) = -grid_length_mdft_xyz(2)*0.5+(i-1)*grid_step_mdft_xyz(2)
  enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, grid_pts_mdft_z, (n_pts_mdft_z+1) ]
 implicit none
  BEGIN_DOC
  ! Grid points for MDFT in the z direction
  END_DOC
  integer :: i
  do i = 1, n_pts_mdft_z+1
    grid_pts_mdft_z(i) = -grid_length_mdft_xyz(3)*0.5+(i-1)*grid_step_mdft_xyz(3) 
  enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, grid_pts_mdft_xyz, (3,n_pts_mdft_tot) ]
 implicit none
  BEGIN_DOC
  ! Grid containing all the points for MDFT in the x/y/z directions
  END_DOC
  integer :: i,j,k,l
  double precision :: wall0,wall1
  print*,'Providing the point grid...'
  call wall_time(wall0)
  l = 0
  do i = 1, n_pts_mdft_x+1
    do j = 1, n_pts_mdft_y+1
      do k = 1, n_pts_mdft_z+1
        l += 1
        grid_pts_mdft_xyz(1,l) = grid_pts_mdft_x(i)
        grid_pts_mdft_xyz(2,l) = grid_pts_mdft_y(j)
        grid_pts_mdft_xyz(3,l) = grid_pts_mdft_z(k) 
      enddo
    enddo
  enddo 
  call wall_time(wall1)
  print*, 'Time to provide the points grid = ', wall1 - wall0, 's'
END_PROVIDER



