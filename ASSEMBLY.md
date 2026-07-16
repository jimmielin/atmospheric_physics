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

- `schemes/pumas/pumas`: nusbaume/PUMAS @ `c4aec4a` (gitlink `c04b38ad`),
  from pumas_round3's .gitmodules (taken as-is).
- `schemes/rrtmgp/ext`, `schemes/musica/...`: unchanged from base.

## Post-merge work on this branch

- suite_cam5.xml GAP-E assembly (separate commits; insertion manifest in the
  scoping doc).

Anything else changed here to make the run work MUST be added to the fix
register in the scoping doc with a durable home (our unit branches /
pumas_round3 PR to Cheryl+Jesse / recorded-here-only).
