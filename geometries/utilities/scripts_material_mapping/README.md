# Material mapping

Tools to build and validate the ACTS surface **material map** for the ALICE 3
tracking geometry (Gen1 or Gen3). The whole chain is driven by
`run_material_mapping.sh`:

```bash
./run_material_mapping.sh --geometry <geom_dir> --geomGeneration <gen1|gen3> -n <nEvents>
```

`<geom_dir>` holds `o2sim_geometry.gdml` (+ `o2sim_geometry.root` for Gen3).

*Example command:*
`./run_material_mapping.sh --geometry ../../geometry_2026_06_10_segmented_withoutRICH_noEndcapTOFs_Gen3/geom  --geomGeneration gen3 --events 1000`

## Pipeline

1. **Record** (`material_recording.py`) — a Geant4 geantino scan over the
   detailed GDML records the real material along each ray →
   `geant4_material_tracks.root`. Geometry-generation agnostic; this is the
   ground truth.

2. **(Gen1 only) receiver config** (`geometry.py`,
   `GeometryVisualisationAndMaterialHandling.py`, `writeMapConfig.py`,
   `modifyConfigMap.py`, `configureMap.py`) — declare which surfaces receive
   material and set their binning → `geometry-map.json`. Gen3 declares receivers
   and binning inside the C++ builder, so these steps are skipped.

3. **Map** (`material_mapping.py --geomGeneration ...`) — assign the recorded
   geant material onto the tracking-geometry surfaces and bin it →
   `material_map.json` (the map), `material_mapped.root` (assigned material per
   surface), `material_unmapped.root` (leftover not attached to any surface).

4. **Validate** (`material_validation.py --geomGeneration ...`) — apply the map
   to the geometry and re-integrate it along tracks → `propagation-material.root`
   (the binned map as seen by the assigner).

5. **Global validation** (`Mat_map.C`, `CheckMaterial.cpp` → `CheckMaterial.exe`)
   — compare the geantino truth against the mapped/validated material (X/X0 vs
   eta, phi, r–z) → `Validation/`.

6. **Per-surface validation** (`material_validation_per_surface.C`) — check the
   binning fidelity per surface: `propagation-material.root` (binned map) vs
   `material_mapped.root` (assigned map), with 1D (z/r, phi), 2D, and ratio
   plots → `Validation/surfaces/`.

The resulting `material_map.json` is what the reconstruction/simulation chain
loads onto the geometry.

## Notes

- `--geomGeneration gen1|gen3` selects the geometry build and, for Gen3, skips
  the `geometry-map.json` receiver-config steps (binning is fixed in the
  builder).
</content>
</invoke>
