# hplin/cam5_mam4_full — assembly record (atmospheric_physics)

**EPHEMERAL INTEGRATION BRANCH — never a PR.** Rebuild it from the recipe
below whenever an ingredient moves; do not merge updates into it, and never
cherry-pick from this branch back into the unit branches. Master copy of
this record + the full fix register:
`~/devel/design_docs/mam_project_scratchpad/cam5_mam4_assembly_scoping.md`.

Assembled 2026-07-16.

## Ingredients (exact heads)

| Ingredient | Head | Content |
|-|-|-|
| `hplin/modal_aero_rebased_on_bulk_aero_3` (base) | `03f4fe7` | BAM + MAM + microp_aero + emissions + cloud_water + MVP sulfur chemistry; includes FIX-5 (hetfrz flag std name -> do_heterogeneous_ice_nucleation) |
| `hplin/uwshcu` | `50490d9` | UW shallow convection + suite_cam5 UW block |
| `hplin/cam5_macrop` | `ce806db` | Park macrophysics (+ upstream main b98c7f5, 979c6bc ride along) |
| `pumas_round3` | `25e1c7f` | Jesse's PUMAS CCPP + FIX-1 (dims_post std names), FIX-2 (naai/dust-dim renames), FIX-3 (optics limiter) |

## Merge order and conflict resolutions

```
git checkout -b hplin/cam5_mam4_full 03f4fe7
git merge hplin/uwshcu       # CONFLICT suites/suite_cam5.xml: keep UW block
                             # (theirs) THEN upstream's Park-placeholder +
                             # compute_cloud_fraction (ours)
git merge hplin/cam5_macrop  # clean
git merge pumas_round3       # CONFLICT clamp_number_concentrations.F90:
                             # take pumas_round3 (4 species, no graupel;
                             # 5-species variant is wrong, FIX-6)
```

rerere is enabled in this repo; both resolutions are recorded and replay on
a rebuild.

## Submodules (checked out manually, not by git)

- `schemes/pumas/pumas`: nusbaume/PUMAS @ `c4aec4a` (gitlink corrected to
  match in `b795cc8`; it had been stale at `c04b38ad`, whose two
  pumas_*cloud_liquid* names mismatch our interstitials).
- `schemes/rrtmgp/ext`, `schemes/musica/...`: unchanged from base.

## Post-merge work on this branch

- suite_cam5.xml GAP-E assembly (separate commits; insertion manifest in the
  scoping doc).
- `2dbfe2b` FIX-12: stratiform snow rate std name ->
  lwe_large_scale_snowfall_rate_at_surface (rk_stratiform +
  set_surface_coupling_vars).
- `d8029e4` FIX-15: Beljaars stress dummies -> taux_beljaars/tauy_beljaars
  (capgen local-vs-dummy collision with the vdiff stress pair).
- `b795cc8` FIX-17: backport of pumas_round3-tip b4b post-interstitial
  fixes (pumasr3 9020f31+70cc4db): pumas_post_main recomputes
  dei/pgam/lamc/des/degrau for radiation from updated constituents;
  suite_cam5 moves pumas_post_main + optics limiter BEFORE
  apply_heating_rate. Mechanically reverse-swept from the new
  standard-name dialect using
  `~/devel/design_docs/mam_project_scratchpad/tools/std_name_sweep.py`
  with `tools/waves/renames_2026-07_pumas_round3.tsv` (48 pairs; ledger
  validated complete against pumasr3 tip). Use the same ledger with
  `sweep --map ... [--reverse]` for any further cross-dialect port.
- `c1e6e2b` FIX-19 (user): removed the duplicate ZMDLF add/out pair from
  park_macrophysics_diagnostics (convect_shallow_diagnostics already
  registers it; same total detrainment value, so field content unchanged).
- `b29863f` FIX-18: chem_srf_emissions + chem_extfrc treat the host
  'UNSET' sentinel as end-of-list. The namelist generator initializes
  char namelist vars to 'UNSET', not blanks, so CAM's exit-on-blank parse
  loop ran on into the unused elements and produced an empty species name.
- `625d44c` FIX-28: mam_mode_metadata resolves the renaming pairs per
  modal_accum_coarse_exch — one (aitken->accum) when false, three when
  true — mirroring CAM's two init routines. It always built three, but the
  cam5 suite sets the flag false, so the run reached
  modal_aero_rename_no_acc_crs_sub, which supports one pair only.
- `cf9e211` FIX-27: rrtmgp_lw_calculate_fluxes + rrtmgp_post declare
  flwds/netsw with the `_to_coupler` standard names, so they resolve to
  the registry cam_out variables instead of group-locals that the cap
  allocated, filled and threw away (leaving Faxa_lwdn at 0 -> CLM's
  "Longwave down ... is negative or zero"). The SW siblings already did
  this. Latent upstream: ours is the first active-land CAM-SIMA run.
- `eef79f9` FIX-26: aerosol_optics constructs the volcanic radius
  constituent name per bin from the optics type (volcanic_radius1 ->
  VOLC_RAD_GEOM1), as CAM does. It had looked up one unsuffixed
  VOLC_RAD_GEOM outside the bin loop, which only matches the legacy
  single-field strataero file; our 3-mode modal file registers
  VOLC_RAD_GEOM1/2/3 and needs a radius per mode.

## MUST PORT BACK to the unit branches (do NOT lose these)

These are [ours] fixes made HERE first; this branch is ephemeral, so they
die with it unless ported. Reapply them on the unit branch — do not
cherry-pick FROM the octopus (rebuild-don't-maintain).

- FIX-18 `b29863f` -> `hplin/modal_aero_rebased_on_bulk_aero_3`
  (chem_srf_emissions.F90 + chem_extfrc.F90; rides the emissions unit PR).
- FIX-26 `eef79f9` -> `hplin/modal_aero_rebased_on_bulk_aero_3`
  (aerosol_optics.F90: per-mode VOLC_RAD_GEOM lookup; rides the
  aerosol-optics/strataero unit PR). Not octopus-specific — any modal
  (3- or 5-mode) prescribed_strataero file hits it.
- FIX-28 `625d44c` -> `hplin/modal_aero_rebased_on_bulk_aero_3`
  (mam_mode_metadata.{F90,meta}: resolve renaming pairs per
  modal_accum_coarse_exch). Not octopus-specific — the no_acc_crs path is
  broken for anyone who selects it.
- FIX-19 `c1e6e2b` -> `hplin/cam5_macrop` (park_macrophysics_diagnostics
  ZMDLF removal; the register notes the ZMDLF-ownership question to raise
  with the team when that branch goes upstream).
- FIX-12 `2dbfe2b` + FIX-15 `d8029e4` + FIX-27 `cf9e211` -> small upstream
  atmos_phys PR (rk_stratiform + set_surface_coupling_vars; beljaars_drag;
  rrtmgp_lw_calculate_fluxes + rrtmgp_post). All are upstream-main schemes.

Anything else changed here to make the run work MUST be added to the fix
register in the scoping doc with a durable home (our unit branches /
pumas_round3 PR to Cheryl+Jesse / recorded-here-only).
