"""Compact view of the blind dissent packet (display only; no key)."""

import re
import sys

COMPARATIVE = re.compile(r"minus|difference|regioisomer|favou?r|major|selectiv|orientation.*(lower|higher)", re.IGNORECASE)
code = None
item = {}


def flush():
    if not item:
        return
    text = item.get("meaning", "")
    if COMPARATIVE.search(text) or "category" in item.get("head", ""):
        print("  *", item["head"])
        print("    M:", text[:420])
        print("    B:", item.get("basis", "")[:520])
    else:
        print("  .", item["head"][:110], "|", text[:110])


for line in open(sys.argv[1]):
    line = line.rstrip("\n")
    if line.startswith("====="):
        flush()
        item = {}
        print(line)
    elif line.startswith("  - "):
        flush()
        item = {"head": line[4:]}
    elif line.startswith("    meaning: "):
        item["meaning"] = line[13:]
    elif line.startswith("    basis: "):
        item["basis"] = line[11:]
flush()
