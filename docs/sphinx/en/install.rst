
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

``TeNeS`` depends on the following libraries, but these are downloaded automatically through the build process.

1. `mptensor <https://github.com/smorita/mptensor>`_ 
2. `cpptoml <https://github.com/skystrife/cpptoml>`_
3. `sanitizers-cmake <https://github.com/arsenm/sanitizers-cmake>`_

TeNeS can use MPI and ScaLAPACK for parallel operations of tensors.
MPI and ScaLAPACK must be installed by yourself. If you use homebrew on macOS, for example, type the following command:

.. code::

   brew install open-mpi scalapack

For others, see the official instruction of some MPI implementation (e.g., OpenMPI) and ScaLAPACK.

For ``tenes_simple`` which generates the input file for ``tenes``, 
Python3 is required.
Additionary, the following python packages are required.

1. numpy
2. toml


Install
======================

1. Build ``TeNeS`` by typing the following commands:

::

  $ mkdir build
  $ cd build
  $ cmake  -DCMAKE_INSTALL_PREFIX=<path to install to> ..
  $ make

The executable file ``tests``  will be generated in  ``build/src`` directory.
The default value of the ``<path to install to>`` is ``/usr/local``. 

2. Install ``TeNeS`` by typing the following commands:

::

  $ make install

In this case, ``tenes`` is installed into the ``<path to install to>/bin`` . 

.. admonition:: Disable MPI/ScaLAPACK parallelization

  If you want to disable MPI/ScaLAPACK parallelization, pass ``-DENABLE_MPI=OFF`` option to ``cmake`` command.

.. admonition:: Specify compiler

   ``CMake`` detects your compiler automatically but sometimes this is not what you want. In this case, you can specify the compiler by the following way,

   ::

      $ cmake -DCMAKE_CXX_COMPILER=<path to your compiler> ../


.. admonition:: Use the pre-built mptensor

   ``TeNeS`` is based on the parallerized tensor library ``mptensor``. The build system of ``TeNeS`` installs this automatically, but if you want to use the specific version of the mptensor, please add the following option in cmake.
   ::

      $ cmake -DMPTENSOR_ROOT=<path to mptensor> ../
