============
DET_COMPRESS
============

This utility compresses multideterminant expansions using algebraic
properties of determinants to combine them and produce a shorter
expansion.  The resulting expansion contains orbitals which are
linear combinations of orbitals as matrix elements.

The methodology used by DET_COMPRESS was developed by
Gihan L. Weerasinghe, Pablo Lopez Rios, and Richard J. Needs, and
is described in a paper available at arXiv.org at
http://arxiv.org/abs/1311.3697 .

DET_COMPRESS is distributed both alongside CASINO and separately,
under a different license in each case.  For CASINO see

  http://vallico.net/casinoqmc/

and for the independent DET_COMPRESS distribution see

  http://github.com/plopezrios/det_compress

In what follows, "SRC", "LIB", "EXAMPLES" and "INTERFACE" refer to
the following directories in the CASINO distribution:
  SRC       = utils/det_compress/
  LIB       = lib/
  EXAMPLES  = examples/det_compress/
  INTERFACE = utils/det_compress/interface/
and to the following directories in the stand-alone DET_COMPRESS
distribution:
  SRC       = src/
  LIB       = lib/
  EXAMPLES  = examples/
  INTERFACE = interface/

Note on licensing: different parts of this distribution are covered
by different licenses.  See the bottom of this README file for
licensing information.

PLR 11.2013


===========
Compilation
===========

* The version of DET_COMPRESS distributed with CASINO is integrated
  into the CASINO build system and requires no indepedent set-up.
  This section applies mostly to the stand-alone DET_COMPRESS
  distribution.  Note that the interface demo programs in the
  INTERFACE directory are not compiled by default in the CASINO
  version of DET_COMPRESS.

* By default the utility is configured to work with gcc and gfortran.
  If you have gcc, gfortran and GNU make installed you simply need to
  type 'make' under the root directory of the DET_COMPRESS
  distribution.  The utility binary will appear as bin/det_compress.
  This will also compile the interface demo programs under the
  INTERFACE directory and place them in their respective
  INTERFACE/*/bin directories.

* To use other compilers and/or compiler options, override the FC,
  FFLAGS, CC, and CFLAGS variables when running make with, e.g.,

    make FC="ifort" FFLAGS="-O3" CC="icc" CFLAGS="-O3"

  The Fortran compiler must support the Fortran 95 standard.

* Some C compilers represent function names differently from others,
  and this affects how Fortran code should interface with C.  If you
  get compilation errors, you may want to try passing one of these to
  'make':

    # IBM compilers typically need:
    CFLAGS_FC_INTERFACE="-DF90_NO_UNDERSCORE"

    # PathScale compiler typically needs:
    CFLAGS_FC_INTERFACE="-DF90_DOUBLE_UNDERSCORE"

    # Some old, rare compilers have been known to need:
    CFLAGS_FC_INTERFACE="-DF90_NO_UNDERSCORE -DF90_CAPITALS"

* We distribute the source of the LP_SOLVE library along with
  DET_COMPRESS to simplify the compilation process.  However if you
  wish to use a system-wide installation of LP_SOLVE you can do so by
  passing the following variables to 'make':

    # Location of liblpsolve*.a:
    LPSOLVE_LIBDIR="/usr/lib"

    # Name of liblpsolve*.a (without the 'lib' prefix and the '.a'
    # suffix):
    LPSOLVE_LIBNAME="lpsolve55"

    # Location of LPSOLVE's headers:
    LPSOLVE_INCDIR="/usr/include/lpsolve"

    # Linker flags required to link dependencies of LPSOLVE:
    LDLPSOLVE_DEPS="-lcolamd -ldl"

  For example, the following works under Ubuntu 13.10, having
  installed the package liblpsolve55-dev:

    make LPSOLVE_LIBNAME=lpsolve55 LPSOLVE_LIBDIR=/usr/lib
      LPSOLVE_INCDIR=/usr/include/lpsolve
      LDLPSOLVE_DEPS="-lcolamd -ldl"

* Test the compilation by running an example, e.g., change into one of
  the subdirectories of the EXAMPLES directory and in the case of the
  stand-alone destribution run

    ../../bin/det_compress

  or in the CASINO version, simply run

    det_compress

  Then press 'd' followed by [Enter] to compress the expansion.


=====
Usage
=====

* DET_COMPRESS must be run in a directory that contains an mdet.casl
  file describing a multi-determinant expansion, and it will produce a
  cmdet.casl file describing the compressed multi-determinant
  expansion.

* An mdet.casl file can be produced by CASINO by setting 'runtype' to
  'gen_mdet_casl' in the input file.  The cmdet.casl file is
  automatically picked up by CASINO on subsequent runs in the same
  directory.

* The CASL file format is used by DET_COMPRESS for data exchange.  See
  the INTERFACE directory for sample code that uses the casl.f90
  module to write/read the input/output files of DET_COMPRESS, for use
  with codes other than CASINO.

* DET_COMPRESS requires an option to be passed on the standard input
  to determine the operational level of the compression:

  - 'a' for deduplication only

  - 'b' for deduplication and compression using the greedy algorithm
    and the simple iterative method

  - 'c' for deduplication and compression using LP_SOLVE and the
    simple iterative algorithm

  - 'd' for deduplication and compression using LP_SOLVE and the
    unified iteration algorithm

  The different levels yield increasingly shorter compressed
  expansions and take increasingly more CPU time.

  Levels 'a' and 'b' scale quadratically with initial expansion size.
  Levels 'c' and 'd' use LP_SOLVE, which is in principle non-
  polynomial.  Level 'd' uses the unified iteration algorithm, which
  involves more non-polynomial operations.

  In practice neither LP_SOLVE nor the unified iteration algorithm
  seem to make the computation intractable, even on modest CPUs.  We
  recommend that the 'd' mode be used, resorting to cheaper modes only
  if the compression takes too long.

* Multi-determinant expansions are usually arranged in groups of
  determinants (Configuration State Functions, CSFs).  During wave
  function optimization in QMC it is common practice to keep the
  ratios between the determinant coefficients in each CSFs fixed,
  so that there is effectively one optimizable parameter per CSF.

  DET_COMPRESS produces compressed expansions in a form that exposes
  the same set of optimizable parameters as the uncompressed expansion
  does.  For this to be possible it is necessary to prevent
  "accidental" compression operations that require determinant
  coefficients in different CSFs to be fixed at their initial ratio.

  This restriction can in principle worsen the compression ratio.
  The main SRC/det_compress.f90 code contains an internal flag,
  IGNORE_COEFF_LABELS, which allows the code to ignore this
  restriction if set to .true. (the code then effectively treats all
  the determinants as if they were part of the same CSF).

  With this modification the code may produce somewhat shorter
  expansions, at the expense of modifying their variational freedom
  if optimized within QMC.  While it is not advisable to optimize
  the coefficients in such compressed expansions, one might want to
  use them for DMC calculations, for example.

  However in our testing we have not encountered any cases where the
  code with IGNORE_COEFF_LABELS=.false. yields larger compressed
  expansions than with IGNORE_COEFF_LABELS=.true..  It remains
  possible in principle, but the gains are likely to be very small,
  and the safest option is to stick to the default.


====================
DET_COMPRESS license
====================

This version of DET_COMPRESS (except the contents of the LIB/lpsolve
and LIB/colamd directories, see below) is dual-licensed:

* DET_COMPRESSED is licensed under the proprietary CASINO license when
  distributed as part of the CASINO quantum Monte Carlo package.  See
  CASINO/doc/academic_consent_form.pdf for the license terms.

* DET_COMPRESSED is licensed under the GPLv3 when distributed
  directly on its own.  See http://www.gnu.org/licenses/gpl-3.0.html
  for the license terms.

Subsequent versions of DET_COMPRESS may be licensed differently
at the authors' discretion.  This version of DET_COMPRESS will in any
case remain licensed as stated above.


================
LP_SOLVE license
================

We use and distribute the unmodified upstream source code of the
LP_SOLVE library with DET_COMPRESS under LIB/lpsolve and LIB/colamd.
LP_SOLVE as distributed with DET_COMPRESS is licensed under the
LGPLv2.1.  See http://www.gnu.org/licenses/lgpl-2.1.html for the
license terms.

The stand-alone GPLv3 distribution of DET_COMPRESS complies with the
licensing terms for using and distributing LP_SOLVE.

When distributed under the CASINO license, DET_COMPRESS classifies as
a "work that uses the library" in the terms of the LGPLv2.1.  While
DET_COMPRESS itself is covered by the propietary CASINO license terms,
it can be compiled against a modified version of the LP_SOLVE library,
satisfying the requirements of the LGPLv2.1.  We believe in good
faith that our use and distribution of the LP_SOLVE library in CASINO
complies with the terms of the LGPLv2.1.
