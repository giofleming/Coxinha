import material_merge_rz_view as mrz
from pathlib import Path
import sys
import acts
from AliceActsPythonBindings import buildALICE3GeometryGen3

tg = buildALICE3GeometryGen3(
    # "/Users/igor/ALICE3/actso2/geometries/geometry_2026_06_10_segmented_withoutRICH_noEndcapTOFs_Gen3/geom/o2sim_geometry.root",
    # "/Users/igor/ALICE3/actso2/geometries/geometry_2026_07_31_segmented_with_OT_realistic_staves_eosCard_1.78perc_Gen3/geom/o2sim_geometry.root",
    # "/Users/igor/ALICE3/actso2/geometries/geometry_2026_07_31_segmented_with_OT_realistic_staves_eosCard_1.78perc_Gen3/geom/o2sim_geometry.root",
    # "/Users/igor/ALICE3/actso2/geometries/geometry_2026_08_06_realistic_OT_with_eosCards_1.78percBarrel_4percDiscs_Gen3/geom/o2sim_geometry.root",
    "/Users/igor/ALICE3/actso2/geometries/geometry_2026_08_06_realistic_OT_with_eosCards_4percBarrel_4percDiscs_Gen3/geom/o2sim_geometry.root",
    acts.logging.DEBUG,          # DEBUG prints the per-layer discovery
    "alice3_blueprint.dot",
)


gctx = acts.GeometryContext.dangerouslyDefaultConstruct()

# obj file for 3D viewer
# vis = acts.ObjVisualization3D()
# tg.visualize(vis, gctx)
# vis.write("alice3_gen3.obj")


# rz projection of all cylinder volumes - the most useful validation view
sys.path.insert(0, "/Users/igor/alice/ACTS/Examples/Scripts/Python")

volumes, markers = mrz.collect_rz(tg, gctx)
# mrz.render(volumes, markers, Path("alice3_rz.svg"), title="ALICE 3 Gen3 (r-z)")
mrz.render(volumes, markers, Path("alice3_rz.pdf"), title="ALICE 3 Gen3 (r-z)")
