python sensitive_layer_overlap_check.py  o2sim_geometry.gdml --eta 0.01 \
    --phi-step 0.0001 --png out.png --csv out.csv --no-show \
    --root sensitive_scan.root  \
    --splitIntoRadialLayers  --maxAllowedHoleSize 0.28 \
    --materials TRK_SILICON FT3_SILICON

# python sensitive_layer_overlap.py o2sim_geometry.gdml --eta -0.7 --materials TRK_SILICON FT3_SILICON