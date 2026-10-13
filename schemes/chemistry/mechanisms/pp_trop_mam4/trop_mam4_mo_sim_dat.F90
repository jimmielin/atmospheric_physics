
      module trop_mam4_mo_sim_dat

      private
      public :: set_sim_dat

      contains

      subroutine set_sim_dat( errmsg, errflg )

      use trop_mam4_chem_mods,     only : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass, crb_mass
      use trop_mam4_chem_mods,     only : diag_map
      use trop_mam4_chem_mods,     only : phtcnt, rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
      use trop_mam4_chem_mods,     only : pht_alias_lst, pht_alias_mult
      use trop_mam4_chem_mods,     only : extfrc_lst, inv_lst, slvd_lst
      use trop_mam4_chem_mods,     only : enthalpy_cnt, cph_enthalpy, cph_rid, num_rnts, rxntot
      use trop_mam4_chem_mods,     only : solsym
      use trop_mam4_chem_mods,     only : frc_from_dataset
      use ccpp_kinds,    only : r8 => kind_phys

      implicit none

!--------------------------------------------------------------
!      ... local variables
!--------------------------------------------------------------
      integer :: ios
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      errmsg = ' '
      errflg = 0

      clscnt(:) = (/      0,     0,     0,    26,     0 /)

      cls_rxt_cnt(:,4) = (/      1,     9,     0,    26 /)

      solsym(: 26) = (/ 'bc_a1           ','bc_a4           ','DMS             ','dst_a1          ','dst_a2          ', &
                        'dst_a3          ','H2O2            ','H2SO4           ','ncl_a1          ','ncl_a2          ', &
                        'ncl_a3          ','num_a1          ','num_a2          ','num_a3          ','num_a4          ', &
                        'pom_a1          ','pom_a4          ','SO2             ','so4_a1          ','so4_a2          ', &
                        'so4_a3          ','soa_a1          ','soa_a2          ','SOAE            ','SOAG            ', &
                        'H2O             ' /)

      adv_mass(: 26) = (/    12.011000_r8,    12.011000_r8,    62.132400_r8,   135.064039_r8,   135.064039_r8, &
                            135.064039_r8,    34.013600_r8,    98.078400_r8,    58.442468_r8,    58.442468_r8, &
                             58.442468_r8,     1.007400_r8,     1.007400_r8,     1.007400_r8,     1.007400_r8, &
                             12.011000_r8,    12.011000_r8,    64.064800_r8,   115.107340_r8,   115.107340_r8, &
                            115.107340_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                             18.014200_r8 /)

      crb_mass(: 26) = (/    12.011000_r8,    12.011000_r8,    24.022000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                             12.011000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8, &
                              0.000000_r8 /)

      fix_mass(:  7) = (/ 0.00000000_r8, 31.9988000_r8, 17.0068000_r8, 47.9982000_r8, 62.0049400_r8, &
                          33.0062000_r8, 28.0134800_r8 /)

      clsmap(: 26,4) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                            21,  22,  23,  24,  25,  26 /)

      permute(: 26,4) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                             11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                             21,  22,  23,  24,  25,  26 /)

      diag_map(: 26) = (/    1,   2,   3,   5,   6,   7,   8,  10,  11,  12, &
                            13,  14,  15,  16,  17,  18,  19,  21,  22,  23, &
                            24,  25,  26,  27,  29,  30 /)

      extfrc_lst(:  9) = (/ 'SO2             ','so4_a1          ','so4_a2          ','pom_a4          ','bc_a4           ', &
                            'H2O             ','num_a1          ','num_a2          ','num_a4          ' /)

      frc_from_dataset(:  9) = (/ .true., .true., .true., .true., .true., &
                                  .true., .true., .true., .true. /)

      inv_lst(:  7) = (/ 'M               ', 'O2              ', 'OH              ', 'O3              ', 'NO3             ', &
                         'HO2             ', 'N2              ' /)

      if( allocated( rxt_tag_lst ) ) then
         deallocate( rxt_tag_lst )
      end if
      allocate( rxt_tag_lst(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         errflg = 1
         write(errmsg,'(a,i0)') 'set_sim_dat: failed to allocate rxt_tag_lst; error = ', ios
         return
      end if
      if( allocated( rxt_tag_map ) ) then
         deallocate( rxt_tag_map )
      end if
      allocate( rxt_tag_map(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         errflg = 1
         write(errmsg,'(a,i0)') 'set_sim_dat: failed to allocate rxt_tag_map; error = ', ios
         return
      end if
      rxt_tag_lst(     1:    10) = (/ 'jh2o2                           ', 'jsoa_a1                         ', &
                                      'jsoa_a2                         ', 'OH_H2O2                         ', &
                                      'usr_HO2_HO2                     ', 'DMS_NO3                         ', &
                                      'DMS_OHa                         ', 'SO2_OH_M                        ', &
                                      'usr_DMS_OH                      ', 'SOAE_tau                        ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10 /)
      if( allocated( pht_alias_lst ) ) then
         deallocate( pht_alias_lst )
      end if
      allocate( pht_alias_lst(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         errflg = 1
         write(errmsg,'(a,i0)') 'set_sim_dat: failed to allocate pht_alias_lst; error = ', ios
         return
      end if
      if( allocated( pht_alias_mult ) ) then
         deallocate( pht_alias_mult )
      end if
      allocate( pht_alias_mult(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         errflg = 1
         write(errmsg,'(a,i0)') 'set_sim_dat: failed to allocate pht_alias_mult; error = ', ios
         return
      end if
      pht_alias_lst(:,1) = (/ '                ', '                ', '                ' /)
      pht_alias_lst(:,2) = (/ '                ', 'jno2            ', 'jno2            ' /)
      pht_alias_mult(:,1) = (/ 1._r8, 1._r8, 1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8, .0004_r8, .0004_r8 /)
      if( allocated( num_rnts ) ) then
         deallocate( num_rnts )
      end if
      allocate( num_rnts(rxntot-phtcnt),stat=ios )
      if( ios /= 0 ) then
         errflg = 1
         write(errmsg,'(a,i0)') 'set_sim_dat: failed to allocate num_rnts; error = ', ios
         return
      end if
      num_rnts(:) = (/      2,     2,     2,     2,     3,     2,     1 /)

      end subroutine set_sim_dat

      end module trop_mam4_mo_sim_dat
