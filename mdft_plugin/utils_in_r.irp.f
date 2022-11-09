BEGIN_PROVIDER [ double precision, aos_grid_mdft, (ao_num,n_pts_mdft_tot) ]
BEGIN_DOC
  ! Grid of atomic orbitals in the x/y/Z directions
  END_DOC
  implicit none
  integer :: i,j
  double precision, allocatable :: aos_array(:)
  double precision :: wall0,wall1
  allocate(aos_array(ao_num))
  print*,'Providing the AOs grid...'
  call wall_time(wall0)
  do i = 1, n_pts_mdft_tot
    call give_all_aos_at_r(grid_pts_mdft_xyz(:,i),aos_array)
    do j = 1, ao_num
      aos_grid_mdft(j,i) = aos_array(j)
    enddo
  enddo
  call wall_time(wall1)
  print*,'Time to provide the AOs grid = ', wall1 - wall0, 's'
END_PROVIDER

BEGIN_PROVIDER [ double precision, mos_grid_mdft, (mo_num,n_pts_mdft_tot) ]
  BEGIN_DOC
  ! Grid of molecular orbitals in the x/y/Z directions 
  END_DOC
  implicit none
  integer :: i,j
  double precision, allocatable :: mos_array(:)
  double precision :: wall0,wall1
  allocate(mos_array(mo_num))
  print*,'Providing the MOs grid...'
  call wall_time(wall0)
  do i = 1, n_pts_mdft_tot
    call give_all_mos_at_r(grid_pts_mdft_xyz(:,i),mos_array)
    do j = 1, mo_num
      mos_grid_mdft(j,i) = mos_array(j)
    enddo
  enddo
  call wall_time(wall1)
  print*,'Time to provide the MOs grid = ', wall1 - wall0, 's'
END_PROVIDER

BEGIN_PROVIDER [ double precision, density_grid_mdft, (N_states, n_pts_mdft_tot) ]
  BEGIN_DOC
  ! Grid containing the densities for all grid points and states 
  END_DOC
  implicit none
  integer :: i,j,k,l
  double precision, allocatable :: dm_a(:), dm_b(:)
  double precision :: wall0,wall1
  allocate ( dm_a(N_states), dm_b(N_states) ) 
  print*,'Providing the density grid...'
  call wall_time(wall0)
  density_grid_mdft = 0.d0
  do i = 1, n_pts_mdft_tot
    call dm_dft_alpha_beta_at_r(grid_pts_mdft_xyz(:,i),dm_a,dm_b)
    do j = 1, N_states
      density_grid_mdft(j,i) = dm_a(j) + dm_b(j)
    enddo
  enddo
  call wall_time(wall1)
  print*,'Time to provide the density grid = ', wall1 - wall0, 's'
END_PROVIDER

BEGIN_PROVIDER [double precision, dipolar_moment_x]
  BEGIN_DOC
  ! Dipolar moment of the molecule in the x direction 
  END_DOC
  implicit none
  integer :: i,j
  dipolar_moment_x = 0.d0
  do i = 1, nucl_num
    dipolar_moment_x += nucl_charge (i) * nucl_coord_transp (1,i)
  enddo
  do i = 1, mo_num
    do j = 1, mo_num
      dipolar_moment_x += - mo_dipole_z (j,i) * dm_tot_e_wft (j,i)
    enddo
  enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, dipolar_moment_y]
  BEGIN_DOC
  ! Dipolar moment of the molecule in the y direction 
  END_DOC
  implicit none
  integer :: i,j
  dipolar_moment_y = 0.d0
  do i = 1, nucl_num
    dipolar_moment_y += nucl_charge (i) * nucl_coord_transp (2,i)
  enddo
  do i = 1, mo_num
    do j = 1, mo_num
      dipolar_moment_y += - mo_dipole_z (j,i) * dm_tot_e_wft (j,i)
    enddo
  enddo
END_PROVIDER


BEGIN_PROVIDER [double precision, dipolar_moment_z]
  BEGIN_DOC
  ! Dipolar moment of the molecule in the z direction  
  END_DOC
  implicit none
  integer :: i,j
  dipolar_moment_z = 0.d0
  do i = 1, nucl_num
    dipolar_moment_z += nucl_charge (i) * nucl_coord_transp (3,i)
  enddo
  do i = 1, mo_num
    do j = 1, mo_num
      dipolar_moment_z += - mo_dipole_z (j,i) * dm_tot_e_wft (j,i)
    enddo
  enddo 
END_PROVIDER  

BEGIN_PROVIDER [double precision, dipolar_moment_xyz, (3)]
  BEGIN_DOC
  ! Dipolar moment of the molecule in the x/y/z directions 
  END_DOC
  implicit none
    dipolar_moment_xyz(1) = dipolar_moment_x
    dipolar_moment_xyz(2) = dipolar_moment_y    
    dipolar_moment_xyz(3) = dipolar_moment_z
END_PROVIDER

!!!!!! Lecture du pot V issu du code MDFT --- A voir si il faut créer un autre fichier provider !!!!!!

BEGIN_PROVIDER [ double precision, grid_weight_mdft_at_r, (n_pts_mdft_tot) ]
  BEGIN_DOC
  ! Grid weight at each point r (For now, it's assumed to be equal to the cubic grid volume)
  END_DOC
  integer :: i
  do i = 1, n_pts_mdft_tot
    grid_weight_mdft_at_r(i) = grid_volume_mdft
  enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, pot_mm_at_r, (n_pts_mdft_tot) ]
  BEGIN_DOC
  ! Potential from mdft program at each point r (For now, it's assumed to be equal to 1)
  END_DOC
  integer :: i
  do i = 1, n_pts_mdft_tot
    pot_mm_at_r(i) = 1.d0
  enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, pot_mm_in_aos, (ao_num, ao_num) ]
  BEGIN_DOC
  ! Potential from mdft program in ao base
  END_DOC
  integer :: i,j,k
  do i = 1, ao_num
   do j = 1, ao_num
     pot_mm_in_aos(j,i) = 0.d0
     do k = 1, n_pts_mdft_tot
       pot_mm_in_aos(j,i) += grid_weight_mdft_at_r(k) * aos_grid_mdft(j,k) * aos_grid_mdft(i,k) * pot_mm_at_r(k)
     enddo
    enddo
  enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, pot_mm_in_mos, (mo_num, mo_num) ]
  BEGIN_DOC
  ! Potential from mdft program in mo basis
  END_DOC
  integer :: i,j,k,l
  do i = 1, mo_num
   do j = 1, mo_num
     pot_mm_in_aos(j,i) = 0.d0
     do k = 1, ao_num
       do l = 1, ao_num
       pot_mm_in_mos(j,i) += mo_coef(k,j) * mo_coef(l,i) * pot_mm_in_aos(j,i)
       enddo
     enddo
    enddo
  enddo
END_PROVIDER
