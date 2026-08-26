#!/bin/sh
# Generate the Doxygen API documentation of TeNeS.
#
# Usage: ./generate.sh
#
# Runs doxygen with the Doxyfile in this directory (the working directory
# does not matter; the script changes into its own directory first).
# The result is written to doxygen_out/ next to the Doxyfile:
#   doxygen_out/html/index.html  - the HTML documentation
#   doxygen_out/latex/           - LaTeX sources (run make there for a PDF)

set -eu

cd "$(dirname "$0")"
doxygen Doxyfile

echo
echo "Generated: $(pwd)/doxygen_out/html/index.html"
