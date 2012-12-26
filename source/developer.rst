======================
CHARMM Developer Guide
======================

This document provides a basic guide for understanding CHARMM's
architecture, implementation, and development protocols and tools.
Prospective developers are urged to familiarize themselves with
its contents.

.. _develop_implement:

CHARMM Implementation and Management
------------------------------------

CHARMM is implemented as a single program package, which is
developed for use on a variety of platforms.  The single source
structure makes the program easier to handle and promotes the
program's integrity.

CHARMM was originally written in FLECS, FORTRAN77 and C languages.
Before FORTRAN77, FLECS allowed us to use a variety of
control constructs, e.g., WHEN-ELSE, WHILE, UNLESS, etc.  A FLECS to
FORTRAN translator was used to process FLECS source code to produce
FORTRAN source.  With CHARMM 23, the FLECS source code was converted
to FORTRAN 77. CHARMM 23f2 and later versions are fully in FORTRAN,
except for some machine-specific codes written in C.  All new
code should be written in FORTRAN 95.

Since CHARMM version 22, all files are maintained by utilizing
software engineering tools.  The Subversion utility is used to maintain
the CHARMM source code, documentation and other supporting files.
The Subversion repository resides on charmm.hanyang.ac.kr.  The CHARMM
manager controls both the developmental and released versions of the code.
He schedules contributions from all CHARMM developers.


.. _develop_directories:

CHARMM Directory Structure
--------------------------

CHARMM files are organized in the following directories.  UNIX pathnames
are used throughout the document.  ~/ is the parent directory that
contains the CHARMM main directory, ~/cnnXm.  nn is the version
number, X is the version trunk designator (a for alpha or
developmental, b for beta release and c for gamma or general release)
and m is the revision number.  For example, c24b1 is CHARMM version 24
beta release revision 1.

==================  ===================================================
Directory           Purpose
==================  ===================================================
~/cnnXm             The main directory of the current CHARMM version.
                    The install.com installation script runs in this
                    directory.

~/cnnXm/source      Source files.

~/cnnXm/doc         Documentation

~/cnnXm/test        Testcases

~/cnnXm/toppar      Standard topology and parameter files.

~/cnnXm/support     Holds various support programs and data files for
                    CHARMM. See :doc:`support`.

~/cnnXm/tool        Contains the preprocessor, prefx, and other
                    CHARMM processing/management tools.

~/cnnXm/build       Contains Makefile, module makefiles and the log
                    file of the install make command for each machine
                    in the subdirectory named after the machine type.

~/cnnXm/lib         Contains library files

~/cnnXm/exec        Will hold executables
==================  ===================================================


.. _develop_standards:

Standards (rules) for writing CHARMM code
-----------------------------------------

Because CHARMM is developed by various groups, there are a number of
standards to which all contributors must adhere in order for the
program to remain modifiable, usable, and transferable.  The rules
which have been established towards this end are listed below.

1) CHARMM module and subroutine structure:

   Fortran modules should be used for all code when possible. The
   general form of a CHARMM fortran module is:

   ::

     module SAMPLEMOD
      ! Comments describing the general function of the module
        use chm_kinds
        implicit none

        ! Declarations of public variables
        ! Declarations of private variables

      contains

        subroutine SUBROUTINE1(arg1,arg2,...)
           .
        end subroutine SUBROUTINE1

        subroutine SUBROUTINE2(arg1,arg2,...)
           .
        end subroutine SUBROUTINE2
           .
           .

        !----- other subroutines and functions --------------!

      end module SAMPLEMOD

   The subroutines in a CHARMM module might include, for example,
   setup subroutines, memory allocation/deallocation subroutines,
   output subroutines, and the principal subroutines for the
   method.

   A subroutine in CHARMM should have the following general form:

   ::

      SUBROUTINE DOTHIS(ARG1,ARG2,....
      ! A comment which describes the purpose of this subroutine.
      ! This may include important variables and what their use is
      ! to aid in understanding and modifying the routine.
      ! A description of all passed arrays and arguments if
      ! users need to call this routine.

        use chm_kinds
        use dimens (if dimensioned common blocks are included)
        use number (if commonly used real numbers are used)
        use somemodule, only: var1, var2, func1, sub1

        implicit none

        !---- declare all passed variables here --------
        real(chm_real), intent(in) :: arg1
        integer(chm_int), intent(out) :: arg2

        !---- local ------------
        Declarations of ALL local variables and parameters.
           .
        data statements at end of declarations.

        !----- begin -----------
           .
        Code (liberally documented through comments)
           .
      end subroutine DOTHIS

   Note that data statements, if present, come after all declarations
   and parameter statements, but before the first line of executable code.

   The use of subroutines outside of modules should be avoided because
   the compiler does not check their arguments.


2) All code should be written clearly.  Since the code must be
   largely self-documenting, clarity should not be sacrificed for
   insignificant gains in efficiency.  Variable, function, subroutine
   and module names should be chosen with care so as to help illustrate
   their purpose.  Avoid using single letter variable names except
   for scratch variables in simple loop constructs.  Comments should be
   used where the function of the code is not obvious. Define/explain
   important variables.  Use the appropriate "intent" attribute (in, out,
   inout) in variable declarations wherever possible.

3) Input/Output

   a) The RDCMND routine should be used to read lines from the
      command stream.  XTRANE should be called to be sure that the
      entire command line is parsed.
   b) Short outputs, messages, warnings, and error should be sent to
      unit OUTU (accessed by USE stream) for output.
   c) Any warning and error message should state which subroutine
      generated it.
   d) PSF and parameter unformatted I/O file formats must remain
      upward compatible. Use an ICNTRL array element to indicate
      which version of CHARMM wrote the file. Such upward
      compatibility must be maintained only across release
      versions of CHARMM. In other words, a file format for the
      developmental version may be freely changed until a new release
      version is generated, at which point all future versions must
      be able to read it.
   e) Use as many significant digits as needed but not more.
      Do not use list directed output (use of "*" for format) or print.
      In particular,

      ::

         write(outu,*) var1, var2
         print *, "Value is", e_variable

      should not be used. It makes output unreadable and makes testing on
      different machines difficult.
   g) All output must be performed based on the PRNLEV value.  This
      is used, for example, to restrict I/O to node_0 in parallel
      implementations of the code.  For example,

      ::

        write (OUTU,'(FORMAT)') ITEMs

      should be coded as

      ::

        if (PRNLEV.GT.2) write (OUTU,'(FORMAT)') ITEMs

      where N is an appropriate print level (always >= 2; see also
      :doc:`miscom`).

4) All error conditions must terminate with a CALL WRNDIE(...) statement;
   direct calls to DIE should not be used. The first argument is the
   severity of the error, the second argument must contain the source
   file in angle brackets, and the routine name of the location of
   the error condition. The third argument contains a string
   describing the error condition.

   ::
     call WRNDIE(-3,'<tamd.src> tamd_allocate', &
       'Failed to allocate memory for islct array')

4) All error conditions must terminate with a CALL WRNDIE(...) statement;
    direct calls to DIE should not be used. The first argument is the
    severity of the error, the second argument must contain the source
    file in angle brackets, and the routine name of the location of
    the error condition. The third argument contains a string
    describing the error condition.

      call WRNDIE(-3,'<tamd.src> tamd_allocate', &
           'Failed to allocate memory for islct array')

5) A. Use subroutines chmalloc and chmdealloc from module memory
      to allocate and deallocate memory.

      Example:

      ::

        module mymod

          use chm_kinds
          integer, allocatable, dimension(:) :: myints
          real(chm_real), allocatable, dimension(:) :: myreals

        contains

          subroutine mystart(natom)
            use memory
            integer, intent(in) :: natom
            call chmalloc('mymod.src', 'mystart', 'myints', natom, intg=myints)
            call chmalloc('mymod.src', 'mystart', 'myreals', natom, crl=myreals)
          end subroutine mystart

          subroutine myfinish(natom)
            use memory
            integer, intent(in) :: natom
            call chmdealloc('mymod.src', 'myfinish', 'myreals', natom, crl=myreals)
            call chmdealloc('mymod.src', 'myfinish', 'myints', natom, intg=myints)
          end subroutine myfinish

        end module mymod

   B. If possible, combine memory allocation statements in a
      separate subroutine.  Do the same for memory deallocation.

6) Arrays should be dimensioned with a constant or integer variable
   where possible. Do not use (*) when dimensioning arrays passed as
   arguments; be sure to pass the dimensioning information
   also. Arrays cannot change shape or type when passed.

   ::

     subroutine NEWFORCE(natom,force,fsum)
       use chm_kinds
       implicit none
       integer(chm_int),intent(in) :: natom
       real(chm_real),intent(out) :: fsum
       real(chm_real),intent(inout),dimension(3,natom) :: force
          ---or---
       real(chm_real),intent(inout) :: force(3,natom)

   Check for array overflows.

7) Error checking in general should be as complete as possible.
   Consider checking for overflows (reciprocals of very small
   numbers, exponentials of very large numbers, etc.), square roots
   of negative numbers, arccosine or arcsine of numbers of absolute
   value greater than one, etc. Code should contain checks for error
   conditions where it will not impact performance. Use compiler
   flags for checking, then remove checking flags for production
   compiling and code submission.

7) The code should not use non-standard Fortran 95 features.
   Such features must be restricted to the machine dependent modules,
   or encapsulated in "##IF - ##ELSE - ##ENDIF" preprocessor
   constructs.

8) * Do not use obsolescent or deprecated Fortran constructs, only use
     what is in the current standard.

   * Do not create common blocks, and use the established common blocks
     only when necessary.

   * Do not use entry points.

   * Do not use computed or assigned goto statements.

   * Avoid "goto" when other constructs are available such as "if," "case,"
     "cycle," and "exit."

9) Functions should never be called with a "call" statement.

10) The generic form of an intrinsic function should be used whenever
    possible.  For example, use SQRT(DP) rather than DSQRT(DP).

11) Real or integer constants should be defined as parameters.

    ::

      real(chm_real),parameter :: ONE=1.0D0, THREE=3.0D0, FIVE=5.0D0, &
                                  SEVEN=7.0D0
      integer(chm_int),parameter :: MAXATM=99999

    See ltm/number_ltm.src for frequently used numbers.
    Real numbers may never be placed in a calling sequence.
    All physical constants should be declared as parameters in
    ltm/consta_ltm.src

    Constants and numbers can be used by adding

    ::

      use consta

    to routines or modules.

12) Routines should have no implicitly declared variables.  This means that
    all variables and arrays, whether passed or not, must be explicitly
    declared. Each module must thus contain

    ::

      implicit none

    as the first statement after any "use" statements and before any
    declarations. This obviates the need for having the implicit none
    statement in each of the subroutines contained in the module.
    A routine that is not contained in a module must have the
    implicit none statement as the first statement after any "use"
    statements and before any declarations.

13) There is no upper or lower case rule in CHARMM except that any
    variable or subroutine name should have consistent case throughout
    all of CHARMM code. It is helpful to use all caps for parameter
    variables.

14) No tabs should appear in code or documentation.

15) All strings must be stored in character variables. Use
    character(len=<n>) for character declarations, where <n> is the
    length of the string. a "*" may be used for the length of a passed
    string but this is discouraged.

16) For routine command parsing, the keyword parsing functions INDXA,
    GTRMA, GTRMF, GTRMI, and NEXTA4 should be used.

17) The recommended length for names of frequently used variables is
    4-10 characters.  Avoid single letter variables, except as indices
    for simple loops. Avoid overriding standard Fortan words such as
    SUM or TYPE.

18) All variables must be declared with a kind. Use the kinds available
    in chm_kinds_ltm.src, the CHARMM variable kinds (chm_kinds) module.  If a
    new kind is used, add it to the chm_kinds module.

19) All variables must be initialized before first use. Most
    initialization is done in the setup or initialization routine in a
    module.

20) Other coding conventions make it easier to search through text for
    particular strings using the SEARCH, fpat, or grep commands.
    Poorly placed spaces can make it very difficult to maintain code.
    There cannot be a space within a variable name.  Here are some
    other examples;

    =================            ===================
     Good                         Please Avoid
    =================            ===================
     GOTO                         GO TO
     |CALL DOSOME(...             CALL  DOSOME(...
     |                            CALL DOSOME (...
     ARRAY(5) = 20                |
     ARRAY(5)=20                  |ARRAY (5) = 20
    =================            ===================

21) The ltm directory contains modules that have no dependencies on
    other modules, except for the chm_kinds module. Examples are modules
    containing global parameters and variables.

22) Avoid duplicating code. If a defect is found in some code, and there
    are multiple copies of it, the remedy must be applied to all copies.
    If code is not duplicated, all its users benefit from any improvement
    to it.


.. _develop_tools:

CHARMM Developer Tools
----------------------

CHARMM is available on a variety of computational devices and we
strongly support multiplatform development efforts.  CHARMM tools are
utility programs/procedures for installation, modification,
optimization, etc.  The preprocessor PREFX and utility procedures
for makefile generation are located in ~/cnnXm/tool.  The
FLECS to FORTRAN translator FLEXFORT is no longer needed since CHARMM
c23f2 and was removed from this and later distribution versions.

_develop_makemod:

Module Makefiles and Optimization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The installation script install.com works with a set of makefiles in
~/cnnXm/build/{machine_type}.  These makefiles play the key role in
developing, optimizing and porting CHARMM code on the machine you are
working with.

(1) Porting to Other Machines

    You may begin with the given set of makefiles for a machine close
    in the architecture to the one to which you intend to port CHARMM.
    First you have to decide a name for the machine platform.  For
    example, IBMRS was chosen for the IBM RS/6000 series.

    ::

      cp -r  ~/cnnXm/build/{closely_related_machine_type} \
             ~/cnnXm/build/{your_chosen_machine_type}

   Then delete Makefile in the new build directory and remane
   Makefile_{closely_related_machine_type} to Makefile_{your_machine_type}.
   You may have to modify compile commands and compiler flags in the
   Makefile template.

   Study carefully ~/cnnXm/install.com and modify it if necessary.
   In most cases, you just need to correct echo messages to address your
   machine properly.  Then issue the install.com command.

(2) Optimization

    Once the makefiles are working properly, you can carry out a
    compiler level optimization for the CHARMM version.  FORTRAN compile
    macro's are defined in Makefile_{machine_type}, e.g., $(FC1), $(FC2),
    $(FC3), etc.  Compiler options are bound to these compile macros.  You
    may inspect each module makefiles and set a proper compile command for
    a given FORTRAN source.  For example, the following are the default
    optimization flags for the c24b1 release.  Most of source files are
    compiled by $(FC2) except

    ::

        build/convex/energy.mk
             $(FC0) ehbond.f
             $(FCR) enefst2.f
             $(FCR) enefst2q.f
             $(FC3) enefvect.f

        build/convex/image.mk
             $(FCR) imnbf2p.f
             $(FC3) imnbfp.f
             $(FC0) nbondm.f

        build/convex/manip.mk
             $(FC0) corman.f
             $(FC3) fshake.f
             $(FCR) fshake2.f

        build/convex/nbonds.mk
             $(FCR) enbf2.f
             $(FCR) enbf3.f
             $(FCR) enbf4.f
             $(FCR) enbf5.f
             $(FC3) ewaldf.f
             $(FCR) ewaldf2.f
             $(FCR) nbndf2p.f
             $(FC3) nbndfp.f

        build/convex/quantum.mk
             $(FC0) qmdata.f
             $(FC0) qmene.f
             $(FC0) qmjunc.f
             $(FC0) qmpac.f
             $(FC0) qmset.f


(3) Generating Module Makefiles

    We have included scripts that find all module dependencies.  When you
    want to create the full set of module makefiles, you may use setmk.com
    in ~/cnnXm/tool.

    ::

        setmk.com UNX

    This generates makefiles in ~/cnnXm/build/UNX, which install.com
    copies to ~/cnnXm/build/{machine-type} as needed.

(4) Usage Note on makemod

    When you generate module.mk files from scratch, the FORTRAN
    compile macro $(FC2) is used for all source files.  In order to set
    the compiler option for further optimization, you have to modify
    the module makefiles to set the macro manually.


.. _develop_modify:

The procedure for modifying anything in CHARMM
----------------------------------------------

This procedure describes the steps which should be taken when
modifying a source file in CHARMM.  When you are developing CHARMM
source code, always maintain close contact with the CHARMM manager
and other developers.  Inform them your development plan and which
files you are working on.  See :ref:`develop_checkin` for the procedures
to follow when submitting your developmental code to the CHARMM manager.

1) Get a copy of the current development code.  If you are a CHARMM
   developer and plan to integrate your program into CHARMM in the
   future, make sure that you obtain the most current revision from
   Subversion.  Check with the CHARMM manager.

2) Once you obtain the code, you are branching out from the main
   CHARMM source code control system.  You should record details
   of modification so that you may reproduce them when you check your
   files in with the CHARMM manager.

3) While you make modifications and debug them, follow the guidelines
   in *note Standards::, so that CHARMM code will be consistent.
   If your modification does not involve any changes in the source file
   directory structure and makes no changes in USE statements,
   you may use the module makefiles supplied (with the extension .mk)
   in ~/cnnXm/build/UNX.  If you add/remove any source files,
   reorganize them, modify any USE statements or are porting to
   a machine that is not already supported, you have to build the
   relevant module make files.  See *note makemod:: for more
   information on makemod.

4) In your local ~/cnnXm directory, you may issue the install.com
   command to build the library and the executable.
   See *note Install: (install.doc)Install.
   Your library is built in ~/cnnXm/lib/{machine_type} and the
   executable will be in ~/cnnXm/exec/{machine_type}.  You may find
   the log file {machine_type}.log in ~/cnnXm/build/{machine_type}.

5) If your modification involves a new feature, you should either
   modify an existing test case or make a new test case to demonstrate
   and check its operation.  See *note testing: (testcase.doc), for a
   description of the tests currently available.  If you add a new test
   case, update the ~charmm/doc/testcase.doc file.  Each new test case
   should run in 1 second or less on a modern workstation.

6) If your change involves adding or modifying a command or adding or
   modifying a feature, modify the existing documentation or if none is
   available, create new documentation.  Make sure that the emacs info
   program can read the document and the format of your documentation
   is consistent with other documents. See *note Document::.




.. _develop_document:

How to Document CHARMM Commands and Features
--------------------------------------------

Documentation is an integral part of CHARMM developments.  In order to
document commands and features under development in a consistent
manner, the following documentation format is recommended.  All
documentation should be accessible (readable) through the emacs info
facility.  If you do not know how to insert the info directives, ask the
CHARMM manager for assistance.  If a new functional module is being
introduced, a new .doc file should be created.  For modifications or
extensions of existing CHARMM modules/functions, the preexisting .doc
file should be revised.

Each documentation file, with the extension .doc, should contain

1) One brief paragraph describing the motivation, theory, or
   procedure relevant to the feature being documented.  Here,
   a few references can be given.

2) A table of contents of the documentation (to serve as the info
   menu).

3) The command syntax.

4) A complete description of all the commands, sub-commands, and
   command options.  The syntax, defaults and file names involved
   should be described.  A brief description of what the command
   accomplishes should also be given.  The order in which various
   commands should be invoked should be described.  Relevant commands
   and subcommands can be cross-referenced with a key.

5) One or two examples involving the concepts and commands described
   (No output listing).

The same notation should be followed throughout the documentation.

::

   [...]   optional, can be present only once, if at all.

   {...}   can be repeated any number of times, must be present
           at least once.

   [{...}] or [{...}]  can either be missing or be present any
                       number of times.

   n{...}  must be present exactly n times.

   <A|b>   either  A or B must be present.

Syntax definitions will use literal keywords such as VIBRan,
READ, MINI, VERLet, etc.  These are to be typed as such.

Syntax definitions can also use dummy keywords such as atom_name,
atom_index and atom_type.  The meaning and variable type can be listed
just after the syntax notation.

For literal keywords the documentation and examples will use
uppercase characters immediately followed by zero or more lower case
characters.  Dummy keywords will be written in all lower case.


.. _develop_api:

API Documentation
-----------------

As a complement to the required command and feature documentation,
the free Doxygen tool (http://www.doxygen.org/) can automatically
generate low-level documentation from the source code. Doxygen reads
declarations of modules, subroutines, and variables, and generates HTML
with cross-references showing the relationships between them.

Doxygen also recognizes Fortran comments beginning with !> as text to
include in the HTML output. Some examples of such comments are in
source/pert/lambdadyn.src.

Before you run Doxygen on CHARMM code, download the Doxygen 1.6
source distribution and apply tool/doxygen.patch as follows:

::

     cd doxygen-1.6
     patch -p0 < ../charmm/tool/doxygen.patch

and follow the Doxygen installation instructions. Then the shell command

::

     cd ../charmm
     doxygen tool/Doxyfile

generates documentation in a directory named apidoc. A good starting
point for browsing is apidoc/html/files.html.

If you have Graphviz (http://www.graphviz.org/) installed, Doxygen
can generate call diagrams for each subroutine. To do this, edit
tool/Doxyfile and change CALL_GRAPH and CALLER_GRAPH to YES.




.. _develop_checkin:

Checkin Procedure in CHARMM Management System
---------------------------------------------

CHARMM is maintained as a single-source software package.
The following rules have been established to minimize conflicts
and delays and to allow for error-free integration of CHARMM
developments from many scientists.

1) It is always wise to inform the CHARMM manager about your
   development plan and timetable so that he may better arrange the
   administrative schedule and also prevent you from duplicating the work
   of others.  The list of files you are working on and the nature of
   the modifications should be reported in advance.

2) Developments must be submitted to the CHARMM manager by
   December 30 for inclusion in the February distributions or
   June 30 for the August distributions.

3) Be certain that the submitted developments are based on the most
   recent development version of CHARMM.  The check-in package
   (see below) should compile out-of-the-box when merged with the
   base version.

4) Prior to check-in, run "test.com" after integration of the
   check-in package with the appropriate base version of CHARMM.
   Compare the results in the output directory with those obtained
   when using the base version.  If you are introducing new
   preprocessor (pref.dat) keywords to control the compilation of
   code for new features (which is recommended), check the test.com
   results for executables produced both with and without the new
   keyword(s).

5) The check-in or submission package should include the following:

   *  modified source files (using the CHARMM source directory structure),
   *  updated or new documentation (doc/*.doc) files, and
   *  updated testcase files

   These files should be assembled in a tar archive conforming to the
   CHARMM distribution directory structure.

6) Post a completed project form (http://charmm.hanyang.ac.kr/172) to
   the CHARMM development bulletin board (http://charmm.hanyang.ac.kr/).
   The project form should contain a succinct description of the
   submitted modifications and new features, the name(s) and institutional
   affiliations of the developer(s), the date, the base CHARMM version,
   new preprocessor keywords, and lists of the new and modified files.
   Upload your tar archive as an attachment to the project form.