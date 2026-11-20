! Saves the current stratiform (ast) and convective (concld) cloud fractions
! into their lagged previous-timestep state variables (cldo and concld_prev)
! for droplet activation (where ast - cldo drives activation in expanding cloud)
! and Park macrophysics on the next timestep.
module stratiform_cloud_fraction_save
  implicit none
  private

  public :: stratiform_cloud_fraction_save_run

contains

!> \section arg_table_stratiform_cloud_fraction_save_run Argument Table
!! \htmlinclude stratiform_cloud_fraction_save_run.html
  subroutine stratiform_cloud_fraction_save_run(ncol, pver, top_lev, ast, concld, &
                                                cldo, concld_prev, errmsg, errflg)

    use ccpp_kinds, only: kind_phys

    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: top_lev          ! top vertical level for cloud physics [index]
    real(kind_phys),  intent(in)    :: ast(:,:)         ! stratiform cloud fraction [fraction]
    real(kind_phys),  intent(in)    :: concld(:,:)      ! convective cloud fraction [fraction]
    real(kind_phys),  intent(inout) :: cldo(:,:)        ! previous-timestep stratiform cloud fraction [fraction]
    real(kind_phys),  intent(inout) :: concld_prev(:,:) ! previous-timestep convective cloud fraction [fraction]

    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    errmsg = ''
    errflg = 0

    cldo(:ncol,top_lev:pver) = ast(:ncol,top_lev:pver)
    concld_prev(:ncol,top_lev:pver) = concld(:ncol,top_lev:pver)

  end subroutine stratiform_cloud_fraction_save_run

end module stratiform_cloud_fraction_save
