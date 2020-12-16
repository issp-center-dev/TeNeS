.. highlight:: none

.. _sec-input-summary:

Short summary for input files of TeNeS
---------------------------------------

Input files of TeNeS is written in `TOML <https://github.com/toml-lang/toml/blob/master/versions/ja/toml-v0.5.0.md>`__ format
and each file has some sections.
``tenes_simple`` and ``tenes_std`` read some sections and generate an input file for ``tenes_std`` and ``tenes``, respectively.
``tenes`` reads some sections and performs simulation.

For example, ``tenes_simple`` reads ``model`` and ``lattice`` sections and generates ``tensor``, ``observable``, and ``hamiltonian`` ones.
Additionary, this copies ``parameter`` and ``correlation`` sections.

The following table summarizes how each tool deal with sections.

.. csv-table::
  :header: "Section", ``tenes_simple``, ``tenes_std``, ``tenes``
  :widths: 15, 15, 15, 10

  ``parameter``,   "copy", "in / copy", "in"
  ``model``,       "in",   "",        ""
  ``lattice``,     "in",   "",        ""
  ``tensor``,      "out",  "in / copy", "in"
  ``observable``,  "out",  "copy",    "in"
  ``correlation``, "copy", "copy",    "in"
  ``hamiltonian``, "out",  "in",      ""
  ``evolution``,   "",     "out",     "in"

- "in"

  - Tool uses this section as input

- "out"

  - Tool generates this section in output (= next input)

- "copy"

  - Tool copies this section into output (= next input)

