import re
from pathlib import Path


p = Path("o2sim_geometry.gdml")

text = p.read_text(encoding="utf-8")

pattern = re.compile(
    r'(<position\b[^>]*name="TRKSensor\d+_1inTRKChip\d+pos"[^>]*\by=")([^"]*)(")',
    re.IGNORECASE,
)

text, n = pattern.subn(r'\g<1>0\g<3>', text)

p.write_text(text, encoding="utf-8")

print(f"Done. Replaced {n} entries for TRKSensor.")


pattern = re.compile(
    r'(<position\b[^>]*name="TRKMetalStack\d+_1inTRKChip\d+pos"[^>]*\by=")([^"]*)(")',
    re.IGNORECASE,
)

text, n = pattern.subn(r'\g<1>0\g<3>', text)

p.write_text(text, encoding="utf-8")

print(f"Done. Replaced {n} entries TRKMetalStack.")
