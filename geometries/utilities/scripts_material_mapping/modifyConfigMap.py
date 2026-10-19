# This script only changes binning of approach = 1;
# From terminal run: python changeApproach1.py config-map.json <desired_binning_binPhi,(int)> <desired_binning_binZ,(int)>
# Returns: config-map-binned.json

import sys
import json

if len(sys.argv) < 5:
    conFile = "config-map.json"
    bins = {"binPhi": sys.argv[1], "binZ": sys.argv[2]}
    mapOnBoundaries = True
else:
    conFile = sys.argv[1]
    bins = {"binPhi": sys.argv[2], "binZ": sys.argv[3]}
    mapOnBoundaries = sys.argv[4].lower() in ("true", "1", "yes")

print("Using: conFile=", conFile, "bins=", bins, "mapOnBoundaries (type)=",
      mapOnBoundaries, "(", type(mapOnBoundaries), ")")

with open(conFile, "r") as con_file:
    config = json.load(con_file)
    # bins = {"binPhi":sys.argv[2],"binZ":sys.argv[3]}

surfaces = config["Surfaces"]  # The surfaces dictionary containing keys "1", "2"...

for surface in surfaces:
    for object in surfaces[surface]:  # An object in the list 'surface' contains a dictionary of properties
        if "approach" in object and object["approach"] == 1:
            object["value"]["material"]["mapMaterial"] = True
            # object["value"]["material"]["mappingType"] = "PreMapping" # IA
            for bin_property in object["value"]["material"]["binUtility"]["binningdata"]:
                if bin_property["value"] == "binPhi":
                    bin_property["bins"] = int(bins["binPhi"])
                else:
                    bin_property["bins"] = int(bins["binZ"])

for surface in surfaces:
    for object in surfaces[surface]:  # An object in the list 'surface' contains a dictionary of properties
        if "boundary" in object:
            object["value"]["material"]["mapMaterial"] = mapOnBoundaries
            # object["value"]["material"]["mappingType"] = "PreMapping" # IA
            for bin_property in object["value"]["material"]["binUtility"]["binningdata"]:
                if bin_property["value"] == "binPhi":
                    bin_property["bins"] = int(bins["binPhi"])
                else:
                    bin_property["bins"] = int(bins["binZ"])


config["Surfaces"] = surfaces
con_file = json.dumps(config, indent=4)
with open('config-map-binned.json', 'w') as file:
    file.write(con_file)
