! Constituent-ordered chemistry VMR workspace shared by the chemistry and
! aerosol clusters: mean wet atmospheric mass (chem_mean_mass), MMR -> VMR
! pack (chem_vmr_pack), cluster-entry checkpoint (chem_vmr_checkpoint_entry),
! and the difference-unpack into constituent tendencies (chem_vmr_unpack).
module chem_vmr_workspace
  implicit none
  private

  public :: chem_mean_mass_run
  public :: chem_vmr_pack_run
  public :: chem_vmr_checkpoint_entry_run
  public :: chem_vmr_unpack_run

contains

!> \section arg_table_chem_mean_mass_run Argument Table
!! \htmlinclude arg_table_chem_mean_mass_run.html
  subroutine chem_mean_mass_run(ncol, pver, mwdry, mbar, errmsg, errflg)
    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    integer,          intent(in)  :: pver
    real(kind_phys),  intent(in)  :: mwdry        ! dry-air molecular weight [g mol-1]
    real(kind_phys),  intent(out) :: mbar(:,:)    ! (ncol,pver) mean wet atmospheric mass [g mol-1]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    ! FHIST easy path: mean wet atmospheric mass = dry-air molecular weight.
    mbar(:ncol, :) = mwdry

  end subroutine chem_mean_mass_run

!> \section arg_table_chem_vmr_pack_run Argument Table
!! \htmlinclude arg_table_chem_vmr_pack_run.html
  subroutine chem_vmr_pack_run(ncol, pver, num_q, const_props, q, mbar, &
                               vmr, loffset, errmsg, errflg)

    use ccpp_kinds, only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use ieee_arithmetic,           only: ieee_value, ieee_signaling_nan
    use chem_vmr_metadata,         only: chem_vmr_slot_kind, &
                                         CHEM_VMR_SLOT_SOLVED, CHEM_VMR_SLOT_INVARIANT

    integer,                           intent(in)    :: ncol
    integer,                           intent(in)    :: pver
    integer,                           intent(in)    :: num_q
    type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)   ! (num_q)
    real(kind_phys),                   intent(in)    :: q(:,:,:)         ! (ncol,pver,num_q) mass/number mixing ratio
    real(kind_phys),                   intent(in)    :: mbar(:,:)        ! (ncol,pver) mean wet atmospheric mass [g mol-1]
    real(kind_phys),                   intent(inout) :: vmr(:,:,:)       ! (ncol,pver,num_q) molar mixing ratio
    integer,                           intent(out)   :: loffset
    character(len=*),                  intent(out)   :: errmsg
    integer,                           intent(out)   :: errflg

    integer         :: m
    real(kind_phys) :: molar_mass    ! [kg mol-1] from constituent props
    real(kind_phys) :: adv_mass      ! [g mol-1] advected molar mass
    real(kind_phys) :: nan_poison    ! signaling NaN for non-workspace constituents

    errmsg = ''
    errflg = 0

    ! The workspace VMR indices are the same as CCPP constituent indices.
    loffset = 0

    nan_poison = ieee_value(1.0_kind_phys, ieee_signaling_nan)

    do m = 1, num_q
       select case (chem_vmr_slot_kind(m))
       case (CHEM_VMR_SLOT_SOLVED)
          call const_props(m)%molar_mass(molar_mass, errflg, errmsg)
          if (errflg /= 0) return

          if (molar_mass > 1.0e30_kind_phys) then
             errflg = 1
             write(errmsg,'(a,i0,a)') 'chem_vmr_pack_run: solved slot ', m, &
                  ' has no registered molar mass'
             return
          end if
          adv_mass = molar_mass * 1.0e3_kind_phys        ! kg mol-1 -> g mol-1
          vmr(:ncol, :, m) = mbar(:ncol, :) * q(:ncol, :, m) / adv_mass
       case (CHEM_VMR_SLOT_INVARIANT)
          ! vmr supplied directly (ic read / provider); leave untouched
       case default
          ! not a chemistry-aerosol constituent: poison to avoid accidental use
          vmr(:ncol, :, m) = nan_poison
       end select
    end do

  end subroutine chem_vmr_pack_run

!> \section arg_table_chem_vmr_checkpoint_entry_run Argument Table
!! \htmlinclude arg_table_chem_vmr_checkpoint_entry_run.html
  subroutine chem_vmr_checkpoint_entry_run(ncol, vmr, vmr_entry, &
                                           errmsg, errflg)

    use ccpp_kinds, only: kind_phys

    integer,          intent(in)  :: ncol
    real(kind_phys),  intent(in)  :: vmr(:,:,:)        ! (ncol,pver,num_q) working molar mixing ratio
    real(kind_phys),  intent(out) :: vmr_entry(:,:,:)  ! (ncol,pver,num_q) cluster-entry checkpoint copy
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    errmsg = ''
    errflg = 0

    vmr_entry(:ncol,:,:) = vmr(:ncol,:,:)

  end subroutine chem_vmr_checkpoint_entry_run

!> \section arg_table_chem_vmr_unpack_run Argument Table
!! \htmlinclude arg_table_chem_vmr_unpack_run.html
  subroutine chem_vmr_unpack_run(ncol, pver, num_q, deltat, const_props, mbar, &
                                 vmr, vmr_entry, const_tend, errmsg, errflg)

    use ccpp_kinds, only: kind_phys
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t
    use chem_vmr_metadata,         only: chem_vmr_slot_kind, CHEM_VMR_SLOT_SOLVED

    integer,                           intent(in)    :: ncol
    integer,                           intent(in)    :: pver
    integer,                           intent(in)    :: num_q
    real(kind_phys),                   intent(in)    :: deltat
    type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)     ! (num_q)
    real(kind_phys),                   intent(in)    :: mbar(:,:)          ! (ncol,pver) [g mol-1]
    real(kind_phys),                   intent(in)    :: vmr(:,:,:)         ! (ncol,pver,num_q) post-cluster vmr
    real(kind_phys),                   intent(in)    :: vmr_entry(:,:,:)   ! (ncol,pver,num_q) cluster-entry vmr
    real(kind_phys),                   intent(inout) :: const_tend(:,:,:)  ! (ncol,pver,num_q) constituent tendency
    character(len=*),                  intent(out)   :: errmsg
    integer,                           intent(out)   :: errflg

    integer            :: m
    real(kind_phys)    :: molar_mass    ! [kg mol-1]
    real(kind_phys)    :: adv_mass      ! [g mol-1]

    errmsg = ''
    errflg = 0

    do m = 1, num_q
       if (chem_vmr_slot_kind(m) /= CHEM_VMR_SLOT_SOLVED) cycle
       call const_props(m)%molar_mass(molar_mass, errflg, errmsg)
       if (errflg /= 0) return
       ! chem_vmr_metadata_init asserted a registered molar mass on every
       ! solved slot; guard against the unset sentinel (huge) regardless
       if (molar_mass > 1.0e30_kind_phys) then
          errflg = 1
          write(errmsg,'(a,i0,a)') 'chem_vmr_unpack_run: solved slot ', m, &
               ' has no registered molar mass'
          return
       end if
       adv_mass = molar_mass * 1.0e3_kind_phys
       const_tend(:ncol, :, m) = const_tend(:ncol, :, m) + &
            (adv_mass / mbar(:ncol, :)) * &
            ((vmr(:ncol, :, m) - vmr_entry(:ncol, :, m)) / deltat)
    end do

  end subroutine chem_vmr_unpack_run

end module chem_vmr_workspace
