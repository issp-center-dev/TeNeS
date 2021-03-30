
Install
-------------------

Download
===================
You can download the source code for TeNeS from the `GitHub page <https://github.com/issp-center-dev/TeNeS>`_ . 
If you have git installed on your machine, type the following command to start download:

``$ git clone https://github.com/issp-center-dev/TeNeS``


Prerequisites
======================
The following tools are required for building TeNeS.

1. C++11 compiler
2. CMake (>=3.6.0)

TeNeS depends on the following libraries, but these are downloaded automatically through the build process.

1. `mptensor <https://github.com/smorita/mptensor>`_ 
2. `cpptoml <https://github.com/skystrife/cpptoml>`_

TeNeS can use MPI and ScaLAPACK for parallel operations of tensors.
MPI and ScaLAPACK must be installed by yourself.
For example, if you use Debian GNU/Linux (or Debian based system such as Ubuntu) and have root priviledges,
you can easily install them by the following:

.. code::

   sudo apt install openmpi-bin libopenmpi-dev libscalapack-mpi-dev

For others, see the official instruction of some MPI implementation and ScaLAPACK.

Python3 is required for the input file generators, ``tenes_simple`` and ``tenes_std`` .
Additionary, the following python packages are also required.

1. numpy
2. scipy
3. toml


Install
======================

1. Build TeNeS by typing the following commands (Some environment such as CentOS provides CMake3 as ``cmake3``)::

::

  $ mkdir build
  $ cd build
  $ cmake  -DCMAKE_INSTALL_PREFIX=<path to install to> ..
  $ make

The default value of the ``<path to install to>`` is ``/usr/local``. 

.. admonition:: Parallel Build

  The ``make`` command accepts ``-j <num>`` options and then uses ``<num>`` processes for a parallel building.
  This reduces the time to build TeNeS drastically.

The executable file ``tenes``  will be generated in  ``build/src`` directory.
By typing the following command, tests for ``tenes`` can be done.

::
 
  $ make tests

2. Install TeNeS by typing the following commands:

::

  $ make install

In this case, ``tenes``, ``tenes_std`` and ``tenes_simple`` are installed into the ``<path to install to>/bin`` . 

.. admonition:: Disable MPI/ScaLAPACK parallelization

  If you want to disable MPI/ScaLAPACK parallelization, pass ``-DENABLE_MPI=OFF`` option to ``cmake`` command.
  On macOS, some functions of ScaLAPACK are incompatible with the system's BLAS and LAPACK,
  and TeNeS ends in error. It is recommended to disable MPI parallel.

.. admonition:: Specify compiler

   CMake detects your compiler automatically but sometimes this does not work. In this case, you can specify the compiler by the following way,

   ::

      $ cmake -DCMAKE_CXX_COMPILER=<path to your compiler> ../


.. admonition:: Specify ScaLAPACK

  CMake detects your ScaLAPACK library automatically but may fail. In this case, you can specify the ScaLAPACK library (``<path>/lib/libscalapack.so``) by the following way,

  ::

    $ cmake -DSCALAPACK_ROOT=<path> ../

.. admonition:: Use the pre-built mptensor

   TeNeS is based on the parallerized tensor library ``mptensor``. The build system of TeNeS installs this automatically, but if you want to use the specific version of the mptensor (``<path>/lib/libmptensor.a``), please add the following option in cmake.
   ::

      $ cmake -DMPTENSOR_ROOT=<path> ../


.. admonition:: Specify Python interpreter

   TeNeS tools (``tenes_simple`` and ``tenes_std``) use ``python3`` interpreter which is found in ``PATH`` via ``/usr/bin/env python3``.
   Please make sure that ``python3`` command invokes the interpreter which you want to use, for example, by using ``type python3``.

   If you want to fix the interpreter (or ``/usr/bin/env`` does not exist), you can specify the interpreter by the following way,

   ::

      $ cmake -DTENES_PYTHON_EXECUTABLE=<path to your interpreter> ../
