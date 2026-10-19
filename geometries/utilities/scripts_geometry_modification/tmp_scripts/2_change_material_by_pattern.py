import xml.etree.ElementTree as ET
import re
from pathlib import Path

input_gdml = "o2sim_geometry.gdml"


#####################
#####################
# Step 1: Create a copy of TRK_SILICON --> Technical_SILICON (skip if it already exists)
# (we want to distinguish sensitive material from "bulk" -- important for proper production of hits in Geant4)

NEW_NAME = "Technical_SILICON"

p = Path(input_gdml)
text = p.read_text(encoding="utf-8")


exists_pat = re.compile(
    rf'\bname\s*=\s*["\']{re.escape(NEW_NAME)}["\']',
    re.IGNORECASE,
)

already_exists = False
if exists_pat.search(text):
    print(f"{NEW_NAME} already exists. Nothing done.")
    # raise SystemExit(0)
    already_exists = True


# find TRK_SILICON block
pattern = re.compile(
    r'(<material\s+name="TRK_SILICON"[\s\S]*?</material>)',
    re.IGNORECASE,
)
m = pattern.search(text)
if not m:
    print("Material TRK_SILICON not found")
    exit()
orig_block = m.group(1)
# copy it with new name
copy_block = orig_block.replace(
    'name="TRK_SILICON"',
    f'name="{NEW_NAME}"',
    1,
)
# paste after the original block
insert_pos = m.end(1)
new_text = text[:insert_pos] + "\n\n" + copy_block + text[insert_pos:]
p.write_text(new_text, encoding="utf-8")
print(f"Done. Created {NEW_NAME}.")


#####################
#####################
# Step 2: Replace materials of some volumes

# input_gdml  = "o2sim_geometry.gdml"
# output_gdml = "o2sim_geometry_modified.gdml"

TARGET_KEYWORDS = (
    "TRKDeadzone",
    "TRKMetalStack",
    "TRKChip",
    "TRKModule",
    "ITOFChip",
    "OTOFChip",
)
OLD_MATERIAL = "TRK_SILICON"
# NEW_MATERIAL = NEW_NAME

tree = ET.parse(input_gdml)
root = tree.getroot()

n_replaced = 0

for volume in root.iter("volume"):
    name = volume.get("name", "")

    if "TRKSensor" in name:
        continue

    if any(key in name for key in TARGET_KEYWORDS):
        material = volume.find("materialref")
        if material is not None and (
            material.get("ref") == OLD_MATERIAL or material.get("ref") == "TF3_SILICON"
        ):
            material.set("ref", NEW_NAME)
            n_replaced += 1

    # if "ITOFSensor" in name:
    #     volume.set("name", "ITOFTRKSensor")
    #     print(f"found 'ITOFSensor' --> replaced with 'ITOFTRKSensor'.")

# --- временная запись ---
tree.write(input_gdml, encoding="utf-8", xml_declaration=True)

# --- FIX: убрать пробел перед "/>" ---
with open(input_gdml, "r", encoding="utf-8") as f:
    content = f.read()

content = content.replace(" />", "/>")

with open(input_gdml, "w", encoding="utf-8") as f:
    f.write(content)

print(f"Done. Replaced material in {n_replaced} volumes.")
