!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS) Coupler.
!*
!* FMS Coupler is free software: you can redistribute it and/or modify
!* it under the terms of the GNU Lesser General Public License as
!* published by the Free Software Foundation, either version 3 of the
!* License, or (at your option) any later version.
!*
!* FMS Coupler is distributed in the hope that it will be useful, but
!* WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!* General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS Coupler.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module atm_ice_wave_exchange_mod

  !! FMS
  use FMS
  use atmos_model_mod,     only: atmos_data_type, land_ice_atmos_boundary_type
  use land_model_mod,      only: land_data_type
  use ice_model_mod,       only: ice_data_type
  use ocean_model_mod,     only: ocean_public_type
  use wave_model_mod,      only: wave_data_type, atmos_wave_boundary_type, ice_wave_boundary_type
  implicit none
  private


  !---- exchange grid maps -----

  type(FmsXgridXmap_type), save :: xmap_atm_wav
  integer, save         :: n_xgrid_atm_wav
  type(FmsXgridXmap_type), save :: xmap_ice_wav
  integer, save         :: n_xgrid_ice_wav

  ! Exchange grid indices
  integer :: X2_GRID_ATM, X2_GRID_WAV

  public :: atm_wave_exchange_init, atm_to_wave, ice_wave_exchange_init, ice_to_wave

  integer :: cplClock, fluxLandIceClock
  logical :: do_runoff
  real    :: Dt_cpl
contains

  subroutine atm_wave_exchange_init(Atm, Wav, Atmos_wave_boundary)
    type(atmos_data_type),          intent(in)    :: Atm !< A derived data type to specify atmospheric boundary data
    type(wave_data_type),           intent(inout) :: Wav !< A derived data type to specify wave boundary data
    type(atmos_wave_boundary_type), intent(inout) :: atmos_wave_boundary !< A derived data type to specify properties
                                                                         !! passed from atmos to waves
    !-----------  local variables  --------------
    integer :: is, ie, js, je

    call fms_xgrid_setup_xmap(xmap_atm_wav, (/ 'ATM', 'WAV' /),       &
         (/ atm%Domain, Wav%Domain /),                    &
         "INPUT/grid_spec.nc", Atm%grid )
    ! exchange grid indices
    X2_GRID_ATM = 1; X2_GRID_WAV = 2;
    n_xgrid_atm_wav = max(fms_xgrid_count(xmap_atm_wav),1)
    call fms_mpp_domains_get_compute_domain( Wav%domain, is, ie, js, je )
    !allocate atmos_wave_boundary
    allocate( atmos_wave_boundary%wavgrd_u10_mpp(is:ie,js:je,1) )
    allocate( atmos_wave_boundary%wavgrd_v10_mpp(is:ie,js:je,1) )
    !variables from wave to atmospheric model, added by Biao
    allocate( Wav%hs(is:ie,js:je,1) )
    allocate( Wav%ust_wav(is:ie,js:je,1) )
    allocate( Wav%ustdir_wav(is:ie,js:je,1) )
    allocate( Wav%charn_wav(is:ie,js:je,1) )

    atmos_wave_boundary%wavgrd_u10_mpp(:,:,:) = 0.0
    atmos_wave_boundary%wavgrd_v10_mpp(:,:,:) = 0.0
    !variables from wave to atmospheric model, added by Biao
    Wav%hs(:,:,1) = 0.0
    Wav%ust_wav(:,:,1) = 0.0
    Wav%ustdir_wav(:,:,1) = 0.0
    Wav%charn_wav(:,:,1) = 0.0


  end subroutine atm_wave_exchange_init

  subroutine ice_wave_exchange_init(Ice, Wav, Ice_wave_boundary)
    type(ice_data_type),          intent(in)      :: Ice !< A derived data type to specify ocean/ice boundary data
    type(wave_data_type),           intent(inout) :: Wav !< A derived data type to specify wave boundary data
    type(ice_wave_boundary_type), intent(inout)   :: Ice_wave_boundary !< A derived data type to specify properties
                                                                     !! passed from atmos to waves
    !-----------  local variables  --------------
    integer :: is, ie, js, je

    call fms_xgrid_setup_xmap(xmap_ice_wav, (/ 'WAV', 'OCN' /),       &
         (/ Wav%Domain, Ice%Domain /),                    &
         "INPUT/grid_spec.nc" )
    ! exchange grid indices
    n_xgrid_ice_wav = max(fms_xgrid_count(xmap_ice_wav),1)
    call fms_mpp_domains_get_compute_domain( Wav%domain, is, ie, js, je )
    !allocate land_ice_boundary
    allocate( ice_wave_boundary%wavgrd_ucurr_mpp(is:ie,js:je,1) )
    ice_wave_boundary%wavgrd_ucurr_mpp(:,:,:) = 0.0
    allocate( ice_wave_boundary%wavgrd_vcurr_mpp(is:ie,js:je,1) )
    ice_wave_boundary%wavgrd_vcurr_mpp(:,:,:) = 0.0
    allocate( ice_wave_boundary%wavgrd_Ice_mpp(is:ie,js:je,1) )
    ice_wave_boundary%wavgrd_Ice_mpp(:,:,:) = 0.0

    allocate( wav%ustkb_mpp(is:ie,js:je,wav%num_stk_bands) )
    wav%ustkb_mpp(:,:,:) = 0.0
    allocate( wav%vstkb_mpp(is:ie,js:je,wav%num_stk_bands) )
    wav%vstkb_mpp(:,:,:) = 0.0
    allocate( wav%tauox_wav(is:ie,js:je,1) )
    wav%tauox_wav(:,:,1) =0.0 
    allocate( wav%tauoy_wav(is:ie,js:je,1) )
    wav%tauoy_wav(:,:,1) =0.0

    ! This are a temporary and costly trick to make MPI work
    allocate( wav%glob_loc_X(is:ie,js:je) )
    wav%glob_loc_X(:,:) = 0
    allocate( wav%glob_loc_Y(is:ie,js:je) )
    wav%glob_loc_Y(:,:) = 0

    call fms_mpp_domains_get_compute_domain( Ice%domain, is, ie, js, je )
    allocate( ice_wave_boundary%icegrd_ustkb_mpp(is:ie,js:je,1,wav%num_stk_bands) )
    ice_wave_boundary%icegrd_ustkb_mpp(:,:,:,:) = 0.0
    allocate( ice_wave_boundary%icegrd_vstkb_mpp(is:ie,js:je,1,wav%num_stk_bands) )
    ice_wave_boundary%icegrd_vstkb_mpp(:,:,:,:) = 0.0

    return
  end subroutine ice_wave_exchange_init

  !> Does atmosphere TO wave operations (could do wave TO atmosphere operations too).
  subroutine atm_to_wave( Time, land_ice_atmos_boundary, Wav, Atmos_wave_Boundary )
    use atm_land_ice_flux_exchange_mod, only: id_hs_wav, id_ust_wav, id_ustdir_wav, id_charn_wav, id_un_ref, id_vn_ref 
    type(FmsTime_type),                intent(in) :: Time !< Current time
    type(land_ice_atmos_boundary_type),intent(in) :: land_ice_atmos_boundary
    type(wave_data_type),            intent(in) :: Wav
    type(atmos_wave_boundary_type), intent(inout):: Atmos_wave_Boundary

! ---- local vars, added by Biao ----------------------------------------------------------
    real, dimension(size(Land_Ice_Atmos_Boundary%t,1),size(Land_Ice_Atmos_Boundary%t,2)) :: diag_atm
    real, dimension(n_xgrid_atm_wav) :: &
         ex_Unref_atm,    &
         ex_Vnref_atm,    &
         ex_hs_wav,       &
         ex_ust_wav,      &
         ex_ustdir_wav,   &
         ex_charn_wav


    integer :: remap_method
    logical :: used

    remap_method = 1

    ! initilize variables on Exchange grid
    ex_Unref_atm  = 0.0
    ex_Vnref_atm  = 0.0
    ex_hs_wav     = 0.0
    ex_ust_wav    = 0.0
    ex_ustdir_wav = 0.0
    ex_charn_wav  = 0.0

    !> Put atmospheric variables onto exchange grid, modified by Biao
    call fms_xgrid_put_to_xgrid (land_ice_atmos_boundary%un_ref , 'ATM', ex_Unref_atm , xmap_atm_wav, remap_method=remap_method, complete=.false.)
    call fms_xgrid_put_to_xgrid (land_ice_atmos_boundary%vn_ref , 'ATM', ex_Vnref_atm , xmap_atm_wav, remap_method=remap_method, complete=.true.)

    !> Put wave related variables onto exchange grid, added by Biao
    call fms_xgrid_put_to_xgrid (Wav%hs, 'WAV', ex_hs_wav , xmap_atm_wav)
    call fms_xgrid_put_to_xgrid (Wav%ust_wav, 'WAV', ex_ust_wav , xmap_atm_wav)
    call fms_xgrid_put_to_xgrid (Wav%ustdir_wav, 'WAV', ex_ustdir_wav , xmap_atm_wav) 
    call fms_xgrid_put_to_xgrid (Wav%charn_wav, 'WAV', ex_charn_wav , xmap_atm_wav)

    if (Wav%pe) then
       call fms_xgrid_get_from_xgrid(Atmos_Wave_Boundary%wavgrd_u10_mpp, 'WAV', ex_Unref_atm, xmap_atm_wav)
       call fms_xgrid_get_from_xgrid(Atmos_Wave_Boundary%wavgrd_v10_mpp, 'WAV', ex_Vnref_atm, xmap_atm_wav)
    endif

    !------- output diagnostic variables from wave, added by Biao-----------

    if ( id_hs_wav > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_hs_wav, xmap_atm_wav)
       used = fms_diag_send_data ( id_hs_wav, diag_atm, Time )
    endif

    if ( id_ust_wav > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_ust_wav, xmap_atm_wav)
       used = fms_diag_send_data ( id_ust_wav, diag_atm, Time )
    endif    

    if ( id_ustdir_wav > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_ustdir_wav, xmap_atm_wav)
       used = fms_diag_send_data ( id_ustdir_wav, diag_atm, Time )
    endif

    if ( id_charn_wav > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_charn_wav, xmap_atm_wav)
       used = fms_diag_send_data ( id_charn_wav, diag_atm, Time )
    endif

    if ( id_un_ref > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_Unref_atm, xmap_atm_wav)
       used = fms_diag_send_data ( id_un_ref, diag_atm, Time )
    endif

    if ( id_vn_ref > 0 ) then
       call fms_xgrid_get_from_xgrid (diag_atm, 'ATM', ex_Vnref_atm, xmap_atm_wav)
       used = fms_diag_send_data ( id_vn_ref, diag_atm, Time )
    endif

 
  end subroutine atm_to_wave

  !> Does both ice TO wave and wave TO ice exchange grid operations.
  subroutine ice_to_wave( Time, Ice, Wav, Ice_wave_Boundary )
    type(FmsTime_type),                intent(in) :: Time !< Current time
    type(ice_data_type),            intent(in) :: Ice !< The ice module container
    type(wave_data_type),           intent(in) :: Wav !< The wave module container
    type(ice_wave_boundary_type), intent(inout):: Ice_wave_Boundary !< The ice-wave boundary container

    ! ---- local vars, added by Biao ----------------------------------------------------------
    real, dimension(size(Ice%part_size,1),size(Ice%part_size,2),size(Ice%part_size,3)) :: ice_frac
    real, dimension(size(Ice%part_size,1),size(Ice%part_size,2)) :: ice_frac_2d
    real, dimension(n_xgrid_ice_wav) :: &
         ex_ucurr,  & ! Exchange grid x-current
         ex_vcurr,  & ! Exchange grid y-current
         ex_ustokes,& ! Exchange grid x-Stokes drift
         ex_vstokes,& ! Exchange grid y-Stokes drift
         ex_ice_frac,&
         ex_ust_wav,      &
         ex_charn_wav,    &
         ex_tauox_wav,    &
         ex_tauoy_wav

    integer :: remap_method ! Interpolation method (todo: list options)
    integer :: i_stk
    integer :: i, j, ncat

    remap_method = 1

    ! initilize variables on Exchange grid
    ex_ucurr      = 0.0
    ex_vcurr      = 0.0
    ex_ustokes    = 0.0
    ex_vstokes    = 0.0
    ex_ice_frac   = 0.0
    ex_ust_wav    = 0.0
    ex_charn_wav  = 0.0
    ex_tauox_wav  = 0.0
    ex_tauoy_wav  = 0.0 

    ice_frac      = 0.0
    ice_frac_2d   = 0.0
    ncat = size(Ice%part_size, 3)
    ice_frac_2d = sum(Ice%part_size(:,:,2:ncat), dim=3)
    ice_frac_2d = max(0.0, min(1.0, ice_frac_2d))
    ice_frac(:,:,1) = ice_frac_2d

    ! -> Put Ocean (ice) parameters onto exchange grid
    call fms_xgrid_put_to_xgrid (Ice%u_surf, 'OCN', ex_ucurr , xmap_ice_wav)
    call fms_xgrid_put_to_xgrid (Ice%v_surf, 'OCN', ex_vcurr , xmap_ice_wav)
    call fms_xgrid_put_to_xgrid (ice_frac ,  'OCN', ex_ice_frac , xmap_ice_wav)

    ! -> Put Wave parameters onto exchange grid, added by Biao
    call fms_xgrid_put_to_xgrid (Wav%ust_wav(:,:,1), 'WAV', ex_ust_wav , xmap_ice_wav)
    call fms_xgrid_put_to_xgrid (Wav%charn_wav(:,:,1), 'WAV', ex_charn_wav , xmap_ice_wav)
    call fms_xgrid_put_to_xgrid (Wav%tauox_wav(:,:,1), 'WAV', ex_tauox_wav , xmap_ice_wav)
    call fms_xgrid_put_to_xgrid (Wav%tauoy_wav(:,:,1), 'WAV', ex_tauoy_wav , xmap_ice_wav)    

    ! -> Only on wave-PEs, bring wave information off exchange grid
    if (Wav%pe) then
       call fms_xgrid_get_from_xgrid(Ice_Wave_Boundary%wavgrd_ucurr_mpp(:,:,1), 'WAV', ex_ucurr, xmap_ice_wav)
       call fms_xgrid_get_from_xgrid(Ice_Wave_Boundary%wavgrd_vcurr_mpp(:,:,1), 'WAV', ex_vcurr, xmap_ice_wav)
       call fms_xgrid_get_from_xgrid(Ice_Wave_Boundary%wavgrd_Ice_mpp(:,:,1), 'WAV', ex_ice_frac, xmap_ice_wav)
    endif
    
    if (Ice%pe) then
        call fms_xgrid_get_from_xgrid(Ice%ust_wav, 'OCN', ex_ust_wav, xmap_ice_wav)
        call fms_xgrid_get_from_xgrid(Ice%charn_wav, 'OCN', ex_charn_wav, xmap_ice_wav)
    endif

    do i_stk = 1,wav%num_Stk_bands
      call fms_xgrid_put_to_xgrid (Wav%ustkb_mpp(:,:,i_stk) , 'WAV', ex_ustokes , xmap_ice_wav)
      call fms_xgrid_put_to_xgrid (Wav%vstkb_mpp(:,:,i_stk) , 'WAV', ex_vstokes , xmap_ice_wav)
      ! -> Only on ice-PEs, bring ice information off exchange grid
      if (Ice%pe) then
        call fms_xgrid_get_from_xgrid(Ice_Wave_Boundary%icegrd_ustkb_mpp(:,:,:,i_stk), 'OCN', ex_ustokes, xmap_ice_wav)
        call fms_xgrid_get_from_xgrid(Ice_Wave_Boundary%icegrd_vstkb_mpp(:,:,:,i_stk), 'OCN', ex_vstokes, xmap_ice_wav)
      endif

    enddo

    return
  end subroutine ice_to_wave

end module atm_ice_wave_exchange_mod
