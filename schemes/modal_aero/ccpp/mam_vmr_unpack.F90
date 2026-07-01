! Unpack MAM microphysics VMR tendencies back into CCPP constituent (mass / number
! mixing ratio) tendencies.
!
! CCPP analog of CAM's vmr2mmr at the bottom of the aerosol cluster: the cluster
! schemes return molar-mixing-ratio tendencies (dqdt), which are converted to
! mixing-ratio tendencies and accumulated into ccpp_constituent_tendencies for the
! framework's apply_constituent_tendencies. Inverse of mam_vmr_pack:
!
!     d(q)/dt = (adv_mass / mbar) * d(vmr)/dt
!
! adv_mass = molar_mass [g mol-1], read from the registered constituent props for
! every species (number tracers carry CAM's nominal cnst_mw = 1.0074 g mol-1),
! consistent with mam_vmr_pack. mbar is held fixed over the step, as in CAM.
!
! For the gasaerexch (+ rename) milestone the cluster members are tendency-return, so
! the total interstitial tendency is exactly dqdt and the per-constituent flag is
! dotend. NOTES for the cluster extension:
!   - newnuc/coag mutate vmr in place rather than returning a tendency; when added,
!     this unpack switches to a saved (vmr_final - vmr_initial)/deltat delta.
!   - the cloud-borne tendency (dqqcwdt / dotendqqcw) is added together with rename.
module mam_vmr_unpack

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: mam_vmr_unpack_run

contains

!> \section arg_table_mam_vmr_unpack_run Argument Table
!! \htmlinclude mam_vmr_unpack_run.html
  subroutine mam_vmr_unpack_run(ncol, pver, num_q, const_props, mbar, &
                                dqdt, dotend, const_tend, errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    integer,                           intent(in)    :: ncol
    integer,                           intent(in)    :: pver
    integer,                           intent(in)    :: num_q
    type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)     ! (num_q)
    real(kind_phys),                   intent(in)    :: mbar(:,:)          ! (ncol,pver) [g mol-1]
    real(kind_phys),                   intent(in)    :: dqdt(:,:,:)        ! (ncol,pver,num_q) vmr tendency
    logical,                           intent(in)    :: dotend(:)          ! (num_q)
    real(kind_phys),                   intent(inout) :: const_tend(:,:,:)  ! (ncol,pver,num_q) constituent tendency
    character(len=*),                  intent(out)   :: errmsg
    integer,                           intent(out)   :: errflg

    integer            :: m
    real(kind_phys)    :: molar_mass    ! [kg mol-1]
    real(kind_phys)    :: adv_mass      ! [g mol-1]
    character(len=256) :: cname

    errmsg = ''
    errflg = 0

    do m = 1, num_q
       if (.not. dotend(m)) cycle
       call const_props(m)%molar_mass(molar_mass, errflg, errmsg)
       if (errflg /= 0) return
       ! A registered molar mass is required to convert the tendency back (number
       ! carries CAM's nominal cnst_mw); a missing one is the unset sentinel (huge).
       if (molar_mass > 1.0e30_kind_phys) then
          call const_props(m)%standard_name(cname, errflg, errmsg)
          errflg = 1
          errmsg = 'mam_vmr_unpack_run: constituent '//trim(cname)// &
               ' was registered without a molar_mass; cannot convert vmr<->mmr'
          return
       end if
       adv_mass = molar_mass * 1.0e3_kind_phys
       const_tend(:ncol, :, m) = const_tend(:ncol, :, m) + &
            (adv_mass / mbar(:ncol, :)) * dqdt(:ncol, :, m)
    end do

  end subroutine mam_vmr_unpack_run

end module mam_vmr_unpack
