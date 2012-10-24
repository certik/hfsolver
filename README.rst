Atomic and Molecular Solver
===========================

``hfsolver`` is a quantum chemistry solver. Features:

* Roothaan-Hartree-Fock equations for closed shell molecules
* radial Roothaan-Hartree-Fock equations for closed shell atoms
* second-order single particle many body Green's function (total energy,
  ionization potentials and electron affinities)
* many body perturbation theory (MBPT), order 2, 3 and 4
* Debye screening of electron-nucleus and electron-electron interactions

Bases:

* finite element (FE), Slater Type Orbitals (STO) and Gaussian Type Orbitals
  (GTO) basis for atoms

* GTO for molecules. It can use `Libint <http://sourceforge.net/p/libint>`_ for
  two-particle integrals or the built-in code based on
  `PyQuante <http://pyquante.sourceforge.net/>`_.

How to build
------------

Only cmake, Lapack and gfortran is needed::

    cmake .
    make

Additional features can be turned on in ``CMakeCache.txt`` or on the command
line. See the ``CMakeLists.txt`` for available options.

License
-------

All code is BSD licensed, except for files taken from other projects (some use
a BSD license, some do not). See the
`LICENSE <https://github.com/certik/hfsolver/blob/master/LICENSE>`_ file for
more information.
