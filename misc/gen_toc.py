#!/usr/bin/env python3

# /* TeNeS - Massively parallel tensor network solver /
# / Copyright (C) 2019- The University of Tokyo */
#
# /* This program is free software: you can redistribute it and/or modify /
# / it under the terms of the GNU General Public License as published by /
# / the Free Software Foundation, either version 3 of the License, or /
# / (at your option) any later version. */
#
# /* This program is distributed in the hope that it will be useful, /
# / but WITHOUT ANY WARRANTY; without even the implied warranty of /
# / MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the /
# / GNU General Public License for more details. */
#
# /* You should have received a copy of the GNU General Public License /
# / along with this program. If not, see http://www.gnu.org/licenses/. */

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
