! Conservative sigma-to-sigma vertical interpolation for HEMCO data in
! direct-to-physics-grid mode.
!
! In legacy intermediate-grid mode, MESSy NCREGRID handles vertical
! interpolation as part of combined horizontal+vertical regridding. In direct
! mode, ESMF does horizontal regridding per-level and the vertical step is
! split out here as a per-column overlap-weighted operation.
!
! Sigma edges are ordered surface (index 1, sigma~1) to TOA (index nlev+1,
! sigma~0), matching HEMCO's L=1=surface convention.
module hco_vertregrid_mod
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: HCO_VertRegrid_Column
  public :: HCO_VertRegrid_3D

contains

  ! Conservative sigma-to-sigma interpolation for a single column. Conserves
  ! the column-integrated quantity for layer-mean intensive inputs (e.g.
  ! mixing ratios, kg/m2/s per layer) via overlap-weighted averaging.
  ! Source layers extending past the target grid edges contribute their
  ! boundary value (no data is lost).
  subroutine HCO_VertRegrid_Column(nlev_src, sig_src, data_src, &
                                   nlev_tgt, sig_tgt, data_tgt)
    integer,         intent(in)  :: nlev_src
    real(kind_phys), intent(in)  :: sig_src(nlev_src + 1)
    real(kind_phys), intent(in)  :: data_src(nlev_src)
    integer,         intent(in)  :: nlev_tgt
    real(kind_phys), intent(in)  :: sig_tgt(nlev_tgt + 1)
    real(kind_phys), intent(out) :: data_tgt(nlev_tgt)

    integer  :: L_tgt, L_src
    real(kind_phys) :: tgt_bot, tgt_top
    real(kind_phys) :: src_bot, src_top
    real(kind_phys) :: overlap
    real(kind_phys) :: tgt_thickness
    real(kind_phys) :: weighted_sum

    do L_tgt = 1, nlev_tgt
      ! Target layer bounds (sigma decreases with altitude).
      tgt_bot = sig_tgt(L_tgt)
      tgt_top = sig_tgt(L_tgt + 1)
      tgt_thickness = tgt_bot - tgt_top

      weighted_sum = 0.0_kind_phys

      do L_src = 1, nlev_src
        src_bot = sig_src(L_src)
        src_top = sig_src(L_src + 1)

        ! Both intervals run high-sigma (bottom) to low-sigma (top).
        overlap = max(0.0_kind_phys, min(tgt_bot, src_bot) - max(tgt_top, src_top))

        if (overlap > 0.0_kind_phys) then
          weighted_sum = weighted_sum + data_src(L_src)*overlap
        end if
      end do

      if (tgt_thickness > 0.0_kind_phys) then
        data_tgt(L_tgt) = weighted_sum/tgt_thickness
      else
        data_tgt(L_tgt) = 0.0_kind_phys
      end if
    end do

  end subroutine HCO_VertRegrid_Column

  ! Per-column conservative regrid for a 3-D field. Source sigma edges are
  ! either uniform across columns (sig_src_1d) or per-column (sig_src_3d);
  ! one of the two must be present. Target sigma edges are always per-column.
  subroutine HCO_VertRegrid_3D(ncol, nlev_src, nlev_tgt, &
                               data_src, data_tgt, &
                               sig_tgt, &
                               sig_src_1d, sig_src_3d)
    integer,         intent(in)            :: ncol
    integer,         intent(in)            :: nlev_src
    integer,         intent(in)            :: nlev_tgt
    real(kind_phys), intent(in)            :: data_src(ncol, nlev_src)
    real(kind_phys), intent(in)            :: sig_tgt(ncol, nlev_tgt + 1)
    real(kind_phys), intent(out)           :: data_tgt(ncol, nlev_tgt)
    real(kind_phys), intent(in), optional  :: sig_src_1d(nlev_src + 1)
    real(kind_phys), intent(in), optional  :: sig_src_3d(ncol, nlev_src + 1)

    integer :: I

    do I = 1, ncol
      if (present(sig_src_3d)) then
        call HCO_VertRegrid_Column(nlev_src, sig_src_3d(I, :), data_src(I, :), &
                                   nlev_tgt, sig_tgt(I, :),    data_tgt(I, :))
      else if (present(sig_src_1d)) then
        call HCO_VertRegrid_Column(nlev_src, sig_src_1d,       data_src(I, :), &
                                   nlev_tgt, sig_tgt(I, :),    data_tgt(I, :))
      else
        ! Neither source sigma provided - should not happen; zero out.
        data_tgt(I, :) = 0.0_kind_phys
      end if
    end do

  end subroutine HCO_VertRegrid_3D
end module hco_vertregrid_mod
