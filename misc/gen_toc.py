#!/usr/bin/env python3

import sys

filename = sys.argv[1] if len(sys.argv) > 2 else "README.md"

with open("README.md") as f:
    for line in f:
        line = line.strip()
        indent = 0
        if not line.startswith("#"):
            continue
        if not line.startswith("##"):
            continue
        elif not line.startswith("###"):
            indent = 0
            body = line.split("##")[1]
        elif not line.startswith("####"):
            indent = 1
            body = line.split("###")[1]
        else:
            continue
        body = body.strip()
        if body.lower() == "table of contents":
            continue

        print(
            "{indent}- [{name}](#{anchor})".format(
                indent="    " * indent,
                name=body,
                anchor=body.lower().replace(" ", "-").replace(".", ""),
            )
        )
