#!/bin/bash

# Default value
GEOMETRY="invalid"
NEVENTS=1000
GEOMGENERATION="gen1"

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --geometry)
            GEOMETRY="$2"
            shift 2
            ;;
        -n|--events)
            NEVENTS="$2"
            shift 2
            ;;
        --geomGeneration)
            GEOMGENERATION="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if [[ "$GEOMGENERATION" != "gen1" && "$GEOMGENERATION" != "gen3" ]]; then
    echo "Invalid --geomGeneration: $GEOMGENERATION (expected gen1 or gen3)"
    exit 1
fi

# Use the geometry variable
echo "Running with geometry: $GEOMETRY"
echo "Number of events: $NEVENTS"
echo "Geometry generation: $GEOMGENERATION"


# record material hits via geantino scan (geometry-generation agnostic: Geant4
# over the GDML, never touches the ACTS tracking geometry)
OUTPUT_GEANTINO_HITS="geant4_material_tracks.root"
python material_recording.py --input ${GEOMETRY}/o2sim_geometry.gdml -n $NEVENTS \
    --output ${OUTPUT_GEANTINO_HITS}

MATERIAL_MAP_OUTPUT_BASE="material"

if [[ "$GEOMGENERATION" == "gen1" ]]; then
    # --- Gen1-only preamble --------------------------------------------------
    # Gen1 declares the material receiver surfaces (and their binning) through a
    # geometry-map.json that is fed to material_mapping.py as --matconfig. Gen3
    # declares them inside the C++ builder (ProtoSurfaceMaterial, fixed binning),
    # so none of the config-generation steps below apply to gen3.

    # visualization of ACTS tracking geometry + write geometry-map.json
    python geometry.py --geo-dir $GEOMETRY
    python $ACTS_SOURCE_DIR/Examples/Scripts/MaterialMapping/GeometryVisualisationAndMaterialHandling.py --geometry geometry-map.json

    # tune binning
    python $ACTS_SOURCE_DIR/Examples/Scripts/MaterialMapping/writeMapConfig.py geometry-map.json config-map.json
    python modifyConfigMap.py config-map.json 80 80 False
    python $ACTS_SOURCE_DIR/Examples/Scripts/MaterialMapping/configureMap.py geometry-map.json config-map-binned.json
fi

# map material
python material_mapping.py --geomGeneration ${GEOMGENERATION} --geo-dir ${GEOMETRY} -n $NEVENTS \
    --input ${OUTPUT_GEANTINO_HITS} --output ${MATERIAL_MAP_OUTPUT_BASE}

# validation
VALIDATION_OUTPUT_BASE="propagation-material"
python material_validation.py --geomGeneration ${GEOMGENERATION} --geo-dir ${GEOMETRY} -n $NEVENTS \
    --map "${MATERIAL_MAP_OUTPUT_BASE}_map.json" -o ${VALIDATION_OUTPUT_BASE}
mkdir -p Validation
root -l -b -q "$ACTS_SOURCE_DIR/Examples/Scripts/MaterialMapping/Mat_map.C(\
    \"${VALIDATION_OUTPUT_BASE}.root\",\
    \"${MATERIAL_MAP_OUTPUT_BASE}_mapped.root\",\
    \"Validation\")"

# Check if the code for making 2D maps have been compiled
if [ ! -f CheckMaterial.exe ]; then
    make
fi

# Run 2D Material maps
./CheckMaterial.exe --G4Input geant4_material_tracks.root --TrkInput ${VALIDATION_OUTPUT_BASE}.root

# --- LAST STEP: per-surface material validation (gen1 and gen3) -------------
# Self-contained macro (no geometry-map.json). Compares the binned map against
# the assigned map (binning fidelity), per surface:
#   A = ${VALIDATION_OUTPUT_BASE}.root  (binned map, assigner)
#   B = ${MATERIAL_MAP_OUTPUT_BASE}_mapped.root  (assigned map, mapping output)
# Output in material_mapping_validation/surfaces/.
root -l -b -q "material_validation_per_surface.C(\
    \"${VALIDATION_OUTPUT_BASE}.root\",\
    \"${MATERIAL_MAP_OUTPUT_BASE}_mapped.root\",\
    \"material_mapping_validation\",\
    \"${GEOMGENERATION}\")"

