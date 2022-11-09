BEGIN_PROVIDER [ double precision, dm_tot_e_wft, (mo_num, mo_num)] 
 implicit none
 BEGIN_DOC
 ! total density matrix for electrons in the WFT 
 END_DOC
 integer :: i,j
  do i = 1, mo_num
   do j = 1, mo_num
    dm_tot_e_wft(j,i) = one_e_dm_mo_alpha_average(j,i)+one_e_dm_mo_beta_average(j,i)
   enddo
  enddo
END_PROVIDER

! call give_all_aos_at_r(r,aos_array) 
!      give_all_mos_at_r_bis(r, mos_array) ! give_all_aos_at_r, mo_coef
! call give_all_mos_at_r(r, mos_array)

