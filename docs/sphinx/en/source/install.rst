
Install
-------------------

Download
===================
You can download the source code for``TeNeS``from the `GitHub page <https://github.com/issp-center-dev/TeNeS>`_ . 
If you have git installed on your machine, type the following command to start download:

``$ git clone https://github.com/issp-center-dev/TeNeS``


Prerequisites
======================
The following tools are required for building TeNeS.

1. C++11 compiler
2. CMake (>=2.8.14)
3. MPI and ScaLAPACK

``TeNeS`` depends on the following libraries, but these are downloaded automatically through the build process.

1. `mptensor <https://github.com/smorita/mptensor>`_ 
2. `cpptoml <https://github.com/skystrife/cpptoml>`_
3. `sanitizers-cmake <https://github.com/arsenm/sanitizers-cmake>`_

ScaLAPACK must be installed by yourself. If you use homebrew in Mac, type the following command:

.. code::

   brew install scalapack

For others, see the web page of ScaLAPACK.

For ``tenes_simple`` which generates the input file for ``tenes``, 
the following libraries are needed.

1. Python (>= ver.3 is recommended)
2. numpy
3. toml


Install
======================

1. Build ``TeNeS`` by typing the following commands:

::

  $ mkdir build
  $ cd build
  $ cmake ..
  $ make

The executable file ``tests``  will be generated in  ``build/src`` directory.
  
2. Install ``TeNeS`` by typing the following commands:

::

  $ cmake -DCMAKE_INSTALL_PREFIX=<path to install to> ../
  $ make
  $ make install

The above installs ``tenes`` into the ``<path to install to>/bin`` . The default value of the ``<path to install to>`` is ``/usr/local``. 

.. admonition:: Specify compiler

   ``CMake`` detects your compiler automatically but sometimes this is not what you want. In this case, you can specify the compiler by the following way,

   ::

      $ cmake -DCMAKE_CXX_COMPILER=<path to your compiler> ../


.. admonition:: Use the pre-built mptensor

   ``TeNeS`` is based on the parallerized tensor library ``mptensor``. The build system of ``TeNeS`` installs this automatically, but you can use the extra pre-built mptensor by the following way.
   ::

      $ cmake -DMPTENSOR_ROOT=<path to mptensor> ../
