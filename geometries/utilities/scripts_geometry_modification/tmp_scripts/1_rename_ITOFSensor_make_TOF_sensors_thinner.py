import re
from pathlib import Path

p = Path("o2sim_geometry.gdml")


# Step 1: rename ITOFSensor to see it in OuterPixels/Barrel ACTS volume

text = p.read_text(encoding="utf-8")
text = text.replace("ITOFSensor", "ITOFTRKSensor")
# p.write_text(text, encoding="utf-8")
print("ITOFSensor --> ITOFTRKSensor replaced.")

p.write_text(text, encoding="utf-8")

# Step 2: now make this sensitive layer thinner

p = Path("o2sim_geometry.gdml")
text = p.read_text(encoding="utf-8")

# 1) find solidref in ITOFTRKSensor
m = re.search(
    r'<volume\s+name="ITOFTRKSensor".*?<solidref\s+ref="([^"]+)"',
    text,
    re.DOTALL,
)
if not m:
    print("ITOFTRKSensor volume not found")
    exit()

solid_name = m.group(1)
print("Making ITOFTRKSensor thinner: solid =", solid_name)

# 2) replace rmax
pattern = re.compile(
    rf'(<tube\b[^>]*name="{re.escape(solid_name)}"[^>]*\brmax=")([^"]*)(")',
    re.IGNORECASE,
)
text, n = pattern.subn(r"\g<1>19.01\g<3>", text)
print(f"Done. Updated rmax ITOFTRKSensor, {n} tube(s).")


# 1) find solidref in OTOFSensor
m = re.search(
    r'<volume\s+name="OTOFSensor".*?<solidref\s+ref="([^"]+)"',
    text,
    re.DOTALL,
)
if not m:
    print("OTOFSensor volume not found")
    exit()

solid_name = m.group(1)
print("Making OTOFSensor thinner: solid =", solid_name)

# 2) replace rmax
pattern = re.compile(
    rf'(<tube\b[^>]*name="{re.escape(solid_name)}"[^>]*\brmax=")([^"]*)(")',
    re.IGNORECASE,
)
text, n = pattern.subn(r"\g<1>85.01\g<3>", text)
print(f"Done. Updated rmax OTOFSensor, {n} tube(s).")


###
p.write_text(text, encoding="utf-8")
