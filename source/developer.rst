======================
CHARMM Developer Guide
======================

This is to provide a guide to someone who wants to understand how
CHARMM is implemented, and a variety of rules that should be followed
by anyone who wishes to modify it.  Anyone who wishes to modify CHARMM
is advised to read through everything in this document.

.. _develop_implement:

CHARMM Implementation and Management
------------------------------------

CHARMM is implemented as a single program package, which is
developed on a variety of platforms.  As a result, it includes some
machine specific implementations and makes heavy use of the virtual
memory capabilities.  By placing everything together, the task of
modifying the program is made more reliable because errors in
modifying the program are more likely to be noticed.  The single
source package concept helps us to maintain integrity of CHARMM as the
paradigmatic macromolecular research software system running on a
variety of platforms.

CHARMM was originally written in FLECS, FORTRAN77 and C languages.
In the past, before FORTRAN77, FLECS allowed us to use a variety of
control constructs, e.g., WHEN-ELSE, WHILE, UNLESS, etc.  A FLECS to
FORTRAN translator was used to process FLECS source code to produce
FORTRAN source.  FORTRAN77 provide us some structured language
constructs.  We began to program directly in FORTRAN77.  The initial
project in CHARMM23 development was to convert all FLECS source codes
into FORTRAN.  CHARMM 23f2 and later versions are fully in FORTRAN
except some machine specific codes written in the C language.  All new
code should be written in FORTRAN77.

Since CHARMM version 22, all files are maintained by utilizing
software engineering tools.  We use the CVS (Concurrent Versions
System) utility to maintain the CHARMM source code, the documentation
and other supporting files.  The CVS repository resides on
tammy.harvard.edu (Convex C220 running under UNIX).  CVS is a superset
of RCS (Revision Control System); file and code management is carried
out with CVS and RCS commands.  The CHARMM manager will be the sole
owner of all CHARMM files and he will schedule/checkin/merge
contributions from both in-house and remote developers.


.. _develop_directories:

CHARMM Directory Structure
--------------------------

CHARMM files are organized in the following directories.  We use UNIX
pathnames throughout the document.  ~/ is the parent directory that
contains the CHARMM main directory, ~/cnnXm.  nn is the version
number, X is the version trunk designator (a for alpha or
developmental, b for beta release and c for gamma or general release)
and m is the revision number.  For example, c24b1 is CHARMM version 24
beta release revision 1.

==================  ===================================================
Directory           Purpose
==================  ===================================================
~/cnnXm             The main directory of the current CHARMM version.
                    install.com procedure runs in this directory.

~/cnnXm/source      Source and include files.

~/cnnXm/doc         Documentation

~/cnnXm/test        Testcases

~/cnnXm/toppar      Standard topology and parameter files.

~/cnnXm/support     Holds various support programs and data files for
                    CHARMM.  See :doc:`support`.

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

Because CHARMM is implemented by a group, there are a number of
conventions which must be observed in order for the program to remain
modifiable, usable, and transferable.  The rules which have been
established towards this end are listed below.

1) Gross subroutine organization:
   
   All INCLUDE statements are processed by the preprocessor to handle
   machine dependent INCLUDE'ing.  ``##INCLUDE`` is the preprocessor
   keyword.  nn as in charmm_nn denotes the version number.  Each
   subroutine should have the following structure.  Note that any
   data statements come after all declarations and parameter
   statements, but before the first line of executable code. 

   ::
   
            SUBROUTINE DOTHIS(ARG1,ARG2,....
      C
      C     A comment which describes the purpose of this subroutine.
      C     This may include important variables and what their use is
      C     to aid in understanding and modifying the routine.
      C
      C     A description of all passed arrays and arguments if
      C     users need to call this routine. 
      C
      ##INCLUDE '~/charmm_fcm/impnon.fcm' (required)
      ##INCLUDE '~/charmm_fcm/dimens.fcm' (if dimensioned common blocks are included)
      ##INCLUDE '~/charmm_fcm/exfunc.fcm' (if external functions are called)
      ##INCLUDE '~/charmm_fcm/number.fcm' (if commonly used real*8 numbers are used)
      C
            declare all passed variables here
      C
      ##INCLUDE '~/charmm_fcm/what_i_need.fcm'
      ##INCLUDE '~/charmm_fcm/more_i_need.fcm'
      C
      C local
                 Declarations of ALL local variables and parameters.
                  .
                 data statements at end of declarations.
      C
      C begin
                  .
                  .
                 Code (liberally documented through comments)
                  .
                  .
            END

2) All code should be written clearly.  Since the code must be
   largely self-documenting, clarity should not be sacrificed for 
   insignificant gains in efficiency.  Variable names should be
   chosen with care so as to illustrate their purpose.  Avoid using
   one or two letter variable names except for scratch variables.
   Comments should be used where the function of code is not obvious.

3) Input/Output

   a) The RDCMND routine should be used to read lines from the
      command stream.  XTRANE should be called to be sure that the
      entire command line is parsed.
   b) Short outputs, messages, warnings, and error should be sent to
      unit OUTU (accessed by ``##INCLUDE '~/charmm_fcm/stream.fcm'``) for
      output. 
   c) All non-fatal messages should state what subroutine generated it.
   d) PSF and parameter unformatted I/O file formats must remain
      upward compatible. Use an ICNTRL array element to indicate
      which version of CHARMM wrote the file. Such upward
      compatibility must be maintained only across production
      versions of CHARMM. In other words, a file format for the
      developmental version may be freely changed until a new version
      is generated, at which point all future versions must be able
      to read it.
   e) I/O of files should be possible in both card and binary format,
      and routines should exist to interconvert between the two.
   f) use as many significant digits as needed but not more;
      in particular ``WRITE(OUTU,*)`` X should be avoided.
      It makes output unreadable and makes testing on
      different machines difficult.
   g) All output must be performed based on the PRNLEV value.  This
      is to enforce only node_0 to carry out I/O on parallel platforms.
      For example,
      
      ::

        WRITE (OUTU,'(FORMAT)') ITEMs

      should be coded as
      
      ::

       IF (PRNLEV.GT.N) WRITE (OUTU,'(FORMAT)') ITEMs

      where N is an appropriate print level.
   h) ``PRINT``, especially ``PRINT,*`` should NOT be used.

4) All error conditions must terminate with a CALL WRNDIE(...);
   direct calls to DIE should not be used; subroutine DIEWRN should be
   phased out.

5) Large or variable storage requirements must be met on the stack or
   heap. When allocating space for the stack or heap, the appropriate
   space allocation subroutine MUST be called. For example, to
   allocate J integer words off the stack, POINTER=ALLSTK(J) is not
   sufficient.  One must use POINTER=ALLSTK(INTEG4(J)) to ensure
   proper performance across different machines.  This also applies
   when freeing the space. The amount of space required for any
   purpose should NEVER be assumed.  This is essential for the
   portability of CHARMM. 

6) Array overflows should be checked for. Error checking in general
   should be as complete as feasible.  Consider checking for
   overflows (reciprocals of very small numbers, exponentials of very
   large numbers, etc.), square roots of negative numbers, arccosine 
   or arcsine of numbers of absolute value greater than one, etc.

7) The code should use a minimum of non-standard Fortran-77 features.
   Such features MUST be restricted to the machine dependent modules,
   or encapsulated in ``##IF - ##ELSE - ##ENDIF`` preprocessor
   constructs.  The only non-Fortran-77 features we use are the INCLUDE 
   statement and the REAL*8 (and INTEGER*2) designators.

8) All common blocks are to be placed in files and INCLUDE'd into
   the program.  The common blocks should have comments describing
   each variable in the common block so that new users will know
   what's there.  The comments should also give clear relationships
   between the variables, so that redimensioning the common block is
   straightforward.  The device on which the file resides must be
   given as ~/charmm_fcm/ where fcm stands for FORTRAN COMMON file.
   The common block files should be named with lower case and have
   the extension .fcm.  Every variable in every common block must be
   declared within the FCM (INCLUDE'd) file.  No Data statements
   should appear in FCM files, and variables declared in FCM's should
   not then be initialized in any other Data statement within
   subroutines.  Moreover, a variable should not appear in more than
   one FCM if there is a possibility that both FCM's will be used in
   the same subroutine.  The multiple declaration will result in an
   error. 

9) Functions should NEVER be called with a CALL statement. Note that
   ENTRY points of functions are also functions.   Moreover, avoid the 
   use of ENTRY points.

10) The generic use of a function should be used unless there is a
    good reason not to.  For example, use SQRT(DP) rather than DSQRT(DP). 

11) Real constants should be defined in PARAMETER statements.  The
    statement above the PARAMETER statement should declare the
    parameter.  Only parameters defined in the following line should
    be declared in such a position.  Double-precision (REAL*8)
    constants should be PARAMETERized with a D.  For example: 
    
    ::
    
      REAL*8 ONED, THREE, FIVE, SEVEN
      PARAMETER (ONED=1.0D0, THREE=3.0D0, FIVE=5.0D0, SEVEN=7.0D0)
      INTEGER MAXATM
      PARAMETER (MAXATM=99999)
      REAL ONES
      PARAMETER (ONES=1.0)
      
    Declaration and Parameter statements should not use continuation
    cards.  See ``~/charmm_fcm/number.fcm`` for frequently used numbers.
    Real numbers may NEVER be placed in a calling sequence!

12) All routines should be up to the IMPLICIT NONE standard.  This
    means that all variables and arrays, whether passed or not, must 
    be declared. This is accomplished by inserting
    ``##INCLUDE '~/charmm_fcm/impnon.fcm'`` in each routine.  The file
    ~/cnnXm/source/fcm/impnon.fcm may then be modified for testing
    purposes, but should contain only comments for normal usage or for 
    machines without an IMPLICIT NONE statement. (Here ~/charmm_fcm/ is
    logically bound to the directory ~/cnnXm/source/fcm)  All
    elements of common blocks MUST be declared in the appropriate
    common file. 

13) All programming should be done in capital letters.  Only comments
    and character strings may use lower case.  No tabs should appear
    in code or documentation. 

14) All strings must be stored in CHARACTER variables.  Although
    integer and real variables will serve on some machines, this is 
    non-standard and eventually causes problems in transportation.

15) For routine command parsing, the keyword parsing functions INDXA,
    GTRMA, GTRMF, GTRMI, and NEXTA4 should be used.

16) The DIMENSION statement should not be used.  Neither should the
    PRINT statement. 

17) Variable names longer than 10 characters should not be used.  Also,
    1 and 2 letter variables should be avoided in large routines
    (except for loop count variables). To facilitate the use of lischeck,
    subroutine and function names may not exceed 13 characters.

18) Precision variables should be ``REAL*8`` (rather than ``DOUBLE PRECISION``).

19) All variables must be initialized before first use.  This may
    best be done in the routine INIALL, which is called very early in
    every CHARMM run. 

20) Other coding conventions make it easier to search through text for
    particular strings using the SEARCH, fpat, or grep commands.
    Poorly placed spaces can make it very difficult to maintain code.
    Never put a space within a variable name.  Here are some other
    examples; 

        =================            ===================   
         Good                         Please Avoid
        =================            ===================
         GOTO                         GO TO
         CALL DOSOME(...              CALL  DOSOME(...
         ARRAY(5) = 20                CALL DOSOME (...
         ARRAY(5)=20                  ARRAY (5) = 20
        =================            ===================
        

.. _develop_tools:

CHARMM Developer Tools
----------------------

CHARMM is available on a variety of computational devices and we
strongly support multiplatform development efforts.  CHARMM tools are
utility programs/procedures for installation, modification, 
optimization, etc.  In ~/cnnXm/tool, we include the preprocessor
PREFX and utility procedures for module makefile generation.  The
FLECS to FORTRAN translator FLEXFORT is no longer needed since CHARMM
c23f2 and removed from this and later distribution versions.

.. _develop_prefx:

CHARMM Preprocessing
--------------------

There is a CHARMM preprocessor, PREFX (formerly PREFLX), which
reads source files as input and produces fortran files for subsequent
compilation.  The main purpose of this preprocessor is to allow a single
version of the source code to work with all platforms and compile options.
A summary of preflx capabilities:

1.  Allows selective compile of machine specific code
2.  Allows selected features to be not compiled (to reduce memory needs)
3.  Supports a size directive to allows larger (and smaller) versions.
4.  Handles the inclusion of .fcm files in a general manner
5.  Allows alternate include file directory to be specified
6.  Allows code expansion for alternate compiles
    (can move IFs from a DO loop).
7.  Allows comments on source lines following a "!"
8.  Handles the conversion to single precision (CRAY, DEC alpha,...)
9.  Identifies unwanted tabs in the source code
10. Checks for line lengths exceeding 72 for non-comments
11. Allows processing multiple files from a list (Macintosh version).
12. Allows the removal of "IMPLICIT NONE" from source files.

The source files have the extension ".src" and the include files have ".fcm".
These files are processed by the preprocessor (PREFX)

SELECTIVE COMPILATION
^^^^^^^^^^^^^^^^^^^^^

Conditional compilation is controlled by simple directives.  The directives
all start with "##" in the first column.  Global keywords should be in upper
case and have multiple letters (local keywords use single character or lower
case).  ##IF constructs can be nested (up to 40 levels).

::

     ##IF keyword(s)     (match-token) ! process code if any keyword is active.
     ##ELIF keyword(s)   (match-token) ! after and IF or IFN,
                                         alternate processing.
     ##ELSE              (match-token) ! after and IF, ELIF, or IFN,
                                         process the rest.
     ##ERROR 'message'                 ! indicates an error within a
                                         ##IF construct
                                         (usually after a ##ELSE condition).
                                         Note: single quotes are required.
     ##ENDIF             (match-token) ! terminates IF, IFN, ELSE,
                                         or ELIF constructs.
     ##IFN keyword(s)    (match-token) ! process code if no keyword is active.


     keywords:: A set of one or more keyword that may be specified in prefx.dat
                or in an ##EXPAND construct (see below).

     match-token :: unique text string in parentheses; must be the
                    same for each use in an ##IF ... ##ENDIF block

Example (from fcm/dimens.fcm):

::

            INTEGER MAXVEC
      ##IFN VECTOR PARVECT               (maxvec_spec)
            PARAMETER (MAXVEC = 10)     
      ##ELIF LARGE XLARGE                (maxvec_spec)
            PARAMETER (MAXVEC = 4000)
      ##ELIF MEDIUM                      (maxvec_spec)
            PARAMETER (MAXVEC = 2000)
      ##ELIF SMALL                       (maxvec_spec)
            PARAMETER (MAXVEC = 2000)
      ##ELIF XSMALL                      (maxvec_spec)
            PARAMETER (MAXVEC = 1000)
      ##ELSE                             (maxvec_spec)
      ##ERROR 'Unrecognized size directive in DIMENS.FCM.'
      ##ENDIF                            (maxvec_spec)


When multiple keywords are specified, an "OR" condition is implied.
IF an "AND" condition is required, use a nested ##IF construct.
In the example above, MAXVEC will not be 10 if either VECTOR or
PARVECT is specified.
 
The text ".not." may be added before a keyname to test for its inverse.
For example, the following constructs are equivalent:

::

         ##IFN BLOCK                    ##IF .not.BLOCK

but these are not equivalent:

::

         ##IFN BLOCK TSM                ##IF .not.BLOCK  .not.TSM

This is because the the first will select when both are false but the second
will select when either is false.

Selective compilation may also be done using on a single line using
a "!##" construct.  The syntax is:

::

        standard-fortran-line    !## keyword(s)  ! comments

A space is not required between the "!##" and the keyword list.
For example the following constructs are equivalent:

Standard format:      

::

      ##IF LONGLINE
            QLONGL=.TRUE.
      ##ELSE
            QLONGL=.FALSE.
      ##ENDIF

Compact format (with comments):      

::

            QLONGL=.TRUE.     !##LONGLINE       ! specify the QLONGL flag
            QLONGL=.FALSE.    !##.not.LONGLINE  ! based on compilation options

Both "and" and "or" conditions can be used for one line processing:

::

            !##PERT  !##PARALLEL   - An "AND" conditional compile
            !##PERT PARALLEL       - An "OR" conditional compile


.. note::
   To assist the automatic ``##IF`` checking utilities, please do not place
   ``##`` characters in the source code (other than on comment lines) unless
   absolutely necessary.

Keyword listing directives;
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``##KEYWORDS LIST unit``

  Inserts a fortran write lines of all current keywords to the selected
  write unit. "unit" may be a variable name or a number but it is limited
  to 8 characters maximum.

* ``##KEYWORDS FILL count array``

  Fills an integer count variable with the number of current keys and
  also fill a character*12 array in the program with the current keys.
  Count and array variable names are limited to 8 characters maximum.


INCLUDE FILES
^^^^^^^^^^^^^

Common files may be included with the ##INCLUDE directive.  The filename must
follow in single quotes.  A directory may preceed the filename with the UNIX
format.

The keyword PUTFCM causes the contents of the included file to be copied and
processed as well.  This is necessary if ## constructs are present in the
included file.

The keyword FCMDIR may override the specified directory in the include
directive.

The VMS keyword will convert the directory name to VMS format.  There is also
special directory name conversion for the Macintosh version.  An included file
may invoke another include file (up to 20 levels).

Example:

::

   ##INCLUDE '~/charmm_fcm/impnon.fcm'


CODE EXPANSION
^^^^^^^^^^^^^^

For computational intensive routines which are not too large, code expansion
may be used to increase efficiency.  This is achieved by moving constant IF
conditions to the outside of major loops. Code expansion is optional and (if
done properly) the code should function in both expanded and unexpanded forms.
This means that the code should be written and tested in an unexpanded
form and then retested with expansion enabled.

::

   ##EXPAND local-flag(s) .when.  conditional-flag(s)   (identifier)

Expand subcommands control section (immediately following the ##EXPAND):

::

     ##PASS1 flag1 flag2 ... 
     ##PASS2 flag1 flag2 ... 
     ##PASS3 ...   - code sections and conditions for each pass
     ##PASS[n] ... 
     ##EXFIN       - code section for the termination of the expand section
     ##EXEND       - end of expansion specification
     ##ENDEX   (identifier)
     
(the identifier is required and must match the corresponding ##EXPAND).
For each pass, the specified flags are temporarily set (or .not. set)
as requested.  If all of the conditions for the code expansion (flags
specified after the .when. construct) are not set, then all flags from
the ##EXPAND line (before the .when.) are temporarily set and no code
expansion is processed.

Example (from nbonds/enbfs8.src):

::

      ...
      ...
      ...
      C Do block expansion of code
      ##EXPAND  B  forces    .when. BLOCK EXPAND  (expand_block)
      ##PASS1  .not.forces
            IF(QBLOCK .AND. NOFORC) THEN
      ##PASS2  forces
            ELSE IF(QBLOCK) THEN
      ##PASS3 .not.BLOCK  forces
            ELSE
      ##EXFIN
            ENDIF
      ##EXEND
      C
            DO I=1,NATOMX       ! Begin of main loop
      ...
      ...
      ...
               IF (.NOT. NOFORC) THEN     !##B
      ##IF forces
               DX(I)=DX(I)+DTX
               DY(I)=DY(I)+DTY
               DZ(I)=DZ(I)+DTZ
      ##ENDIF
               ENDIF                      !##B
      ...
      ...
      ...
            ENDDO               ! End of main loop
      
      ##ENDEX    (expand_block)
            RETURN
            END

This example will do a multi pass compilation when BOTH the
"EXPAND" and the "BLOCK" keywords are set.  If they are not both
set, then the local flags "B" and "forces" will be set until
the corresponding ##ENDEX is reached.  If the "EXPAND" and "BLOCK"
conditions are met, then the body of the expanded section will be
compiled three times.

::

 PASS1 - additional active flag:             disabled flag: forces
 PASS2 - additional active flag: forces      disabled flag: 
 PASS3 - additional active flag: forces      disabled flag: BLOCK


RESERVED KEYWORDS
^^^^^^^^^^^^^^^^^

The following keywords are reserved:

   ============= =======================================================
   END           The end of keywords in prefx.dat (END is not a keyword)
   SINGLE        Conversion to single precision   (SINGLE is a keyword)
   PUTFCM        Include files are to be copied into fortran files
   VMS           Use VMS directory names (from DEC's DCL)
   REMIMPNON     Remove any "IMPLICIT NONE" lines found in the source
   FCMDIR        Specification of include file directory
   UPPERCASE     Convert all non-text code to uppercase Fortran
   LONGLINE      Allows a longer line output format (>80 characters).
   SAVEFCM       Include all SAVE statements
   EXPAND        Do semi-automatic code expansion
   single-letter reserved for unexpanded compile conditionals
   lower-case    reserved for local compile flags (within a routine)
   ============= =======================================================
   
Other Keyword Rules

- Keyword may not exceed 12 characters in length.
- Global keywords must be all uppercase
- Local  keywords must be all lowercase
- Keywords should otherwise follow fortran standards for naming
- Recommendation: Avoid one and two letter keywords (harder to find)

preflx.dat or pref.dat are the preprocessor instruction data files.
Create a file preflx.dat or pref.dat (with UNIX) that contains 
one or more of the keywords specified below.  On UNIX platforms, install.com
generates the default pref.dat file in build/{machine_type} directory.
"END" keyword stops parsing keywords.  The use of a (Match-Token) can
help to identify the components of ##IF blocks in source files that make
heavy use of ## directives; it should follow any keywords, and must be
appended to all components of a given ##IF block (if it is used).  See
the code for more examples.


LIST OF ALL KEYWORDS IN CHARMM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A complete list of all compile flags and options WITH SUITABLE DESCRIPTIONS
will be found in the documentation file; :doc:`preflx_list`.
It can change much between released CHARMM versions.

The information list here highlights reserved keywords and other basic
information that is unlikely to change between versions.

::

   [1] Include File Directory
       FCMDIR=directory_name   ! point to a particular directory
       FCMDIR=CURRENT          ! use what is specified in the include line.
       FCMDIR=LOCAL            ! use the local directory.

   [2] Machine Type (choose exactly one)
       ALLIANT     = Alliant
       ALPHA       = DEC alpha workstation
       APOLLO      = HP-Apollo, both AEGIS and UNIX
       ARDENT      = Stardent, Titan series
       CONVEX      = Convex Computer
       CRAY        = Cray Research Inc.
       DEC         = DEC ULTRIX
       FUJITSU
       HAL         = Sun computer port - special
       GWS         = Sun Global Works System
       HPUX        = Hewlett-Packard series 700.
       IBM         = IBM-3090 running AIX
       IBMRS       = IBM-RS
       IRIS        = Silicon Graphics
       MACINTOSH   = Apple Macintosh computers (system 7)
       SUN         = Sun Microsystems
       ULTRA       = For modern Sun compilers circa 2000.
       VAX         = Digital Equipment Corp. VAX VMS.

     Other machine descriptors
       IBMMVS      = IBM's MVS platform
       IBMVM       = IBM's VM platform
       GNU         = using GNU Fortran compiler
       CMEM        = A convex option?
       GRAPE       = Use MD-GRAPE-II board to speedup nonbond calculations
       LOBOS       = LoBoS cluster specific code

     Parallel machine types
       ALPHAMP     = DEC Alpha Multi Processor machines
       CM5         = Machine type            = TMC's CM-5 machine
       CSPP        = Convex PA-RISC parallel system (HP chip)
       CSPPMPI     = Convec SPP using proprietary MPI library
       DELTA       = machine type            = Intel delta (Caltech) machine
       IBMSP       = machine type            = IBM's SPn cluster machines
       IBMSP1      = machine type            = IBM's SP1 cluster machines
       INTEL       = machine type            = Intel iPSC Hypercube
       PARAGON     = machine type            = Intel Paragon machine
       SGIMP       = machine type            = SGI Power Challenge
       T3D         = Cray massively parallel (DEC Alpha chip)
       T3E         = Cray massively parallel (DEC Alpha chip)
       TERRA       = multiprocessor DEC Alpha chip system

   [3] Operating system (choose at most one)
       AIX370      = IBM UNIX
       UNIX        = UNIX
       UNICOS      = Cray UNIX
       OS2         = IBM pre-emptive multitasking

   [4] Size directive (must choose exactly one)
       XXLARGE     =360720 atom limit
       XLARGE      =240480 atom limit
       LARGE       = 60120 atom limit
       MEDIUM      = 25140 atom limit
       REDUCE      = 15000 atom limit 
       SMALL       =  6120 atom limit
       XSMALL      =  2040 atom limit

   [5] Machine Architecture (may choose several)
       SCALAR      = machine characteristics = default for scalar machines
       VECTOR      = feature directive *     = Vectorized routines
       PARVECT     = Parallel vector code (multi processor vector machines)
       CRAYVEC     = Fast vector code (standard vector code)
       SINGLE      = specifies single precision version (primarily used for CRAY)
       SGIF90      = Used to compile CHARMM using F90 compiler on SGI machines
       T3ETRAJ     = Used to read t3e trajectories on IEEE machines
                     w/ 32 bit integers

   [6] Parallel CHARMM descriptors  (see :doc:`parallel`)
       (all require the PARALLEL keyword)

   [7] Feature directives  (see preflx_:doc:`used`). See;
          *note preflx (chmdoc/preflx_:doc:`used`).

   [8] Graphics keywords; choose only one (except on Apollo)
       GLDISPLAY   = use the GL display code for the graphics window (*)
       NODISPLAY   = no graphics window; PostScript, other files produced
       NOGRAPHICS  = graphics code not compiled
       XDISPLAY    = use the X11 display code for the graphics window

        (*) the GL code is relatively untested, and may have problems

   [9] Keywords Not for Normal Use  (see preflx_:doc:`used`). See;
          *note preflx (chmdoc/preflx_:doc:`used`).

   [10] Major Blocks that can be Removed, but normally are not. See;
          *note preflx (chmdoc/preflx_:doc:`used`).

   [11] Other Control Directives
       EXPAND      = Do semi-automatic code expansion
       LONGLINE    = Allows a longer line output format (>80 characters).
       SAVEFCM     = Include all SAVE statements in .fcm files
       SINGLE      = Conversion to single precision (SINGLE is a keyword)
       PUTFCM      = Include files are to be copied into fortran files
       VMS         = Use VMS directory names (from DEC's DCL)
       REMIMPNON   = Remove any "IMPLICIT NONE" lines found in the source
       UPPERCASE   = Convert all non-text code to uppercase Fortran

By employing appropriate preprocessor keys, one can generate a
variant of CHARMM for a specific machine with specific features.

.. _develop_makemod:

Module Makefiles and Optimization

The installation script install.com works with a set of makefiles in
~/cnnXm/build/{machine_type}.  These makefiles play the key role in
developing, optimizing and porting CHARMM code on the machine you are
working with.

(1) Porting to Other Machines

    You may begin with the given set of makefiles for a machine close
    in the architecture to the one to which you intend to port CHARMM.
    First you have to decide a name for the machine platform.  For
    example, we have chosen IBMRS for IBM RS/6000 series.
    
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

    Once you make the makefiles working properly, you can carry out a
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

    We have included the makemod script that finds all include file
    dependencies.  The makemod script is used for all source modules
    except main, for which we use mainmake instead.

    ::
    
        makemod [-v] [-n] module_name actual_path sourcetree_path \
                          makefile_name [definition_file_name]

    where -v means verbose and -n means the include files don't have
    includes (This saves time).  The module name is general the module
    directory, e.g., dynamc.  The actual path is generally
    ~/cnnXm/source/{module} and the source tree path is ~/cnnXm/source
    or the like.  We generally use module_name.mk for the makefile name.
    The definition file if specified (we generally don't) will prepend a
    file of definitions to the makefile.  For example, the way to generate
    the makefiles is to:

    ::
    
        cd ~/cnnXm/source/{module}
        makemod -vn {module} `pwd` `cd ..;pwd` {module}.mk

    where {module} is the sub-directory in source.

    When you want to create the full set of module makefiles, you may
    use setmk.com in ~/cnnXm/tool.

    ::
    
        setmk.com your_machine_type


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
source code, always maintain close contacts with the CHARMM manager
and other developers.  Inform them your development plan and which
files you are working on.  See :ref:`_developer_checkin` for checkin procedure
needed when you deposit your code in the CHARMM central library.

1) Get a copy of the current release package.  If you are a CHARMM
   developer and plan to integrate your program into CHARMM in the
   future, make sure that you obtain the most current version.  Check
   with the CHARMM manager.

2) Once, you obtain the package, you are branching out from the main
   CHARMM source code control system.  You should record details
   of modification so that you may REDO them when you check your
   files in the central CHARMM library.

3) While you make modifications and debug them, follow the guidelines
   in :ref:`developer_standards`, so that CHARMM code will be consistent.
   If your modification does not involve any changes in source file
   directory structure and makes no changes in INCLUDE statements,
   you may use the module makefiles supplied (with the extension .mk)
   in ~/cnnXm/build/UNX.  If you add/remove any source files,
   reorganize them, modify any INCLUDE statements or are porting to
   other machine than those already supported, you have to build the
   relevant module make files.  See :ref:`developer_tools` for more
   information on makemod.

4) In your local ~/cnnXm directory, you may issue install.com
   command to build the library and the executable.
   See :doc:`install`.
   Your library is built in ~/cnnXm/lib/{machine_type} and the
   executable will be in ~/cnnXm/exec/{machine_type}.  You may find
   the log file {machine_type}.log in ~/cnnXm/build/{machine_type}.

5) If your modification involves a new feature, you should either
   modify an existing test or make a new test to demonstrate and
   check its operation. See :doc:`testcase`,
   for a description of the tests currently available.  If you add a
   new test, update the ~charmm/doc/testcase.doc file.

6) If your change involves adding or modifying a command or adding or
   modifying a feature, modify existing documentation or if none is
   available, make new documentation.  Make sure that the emacs info
   program can read the document and the format of your documentation
   is consistent with other documents.


.. _develop_document:

How to Document CHARMM Commands and Features
--------------------------------------------

Documentation is an integral part of CHARMM developments.  In order to
document commands and features under development in a consistent
manner, we recommend the following documentation format.  All
documentations should be accessible (readable) through the emacs info
facility.  If you do not know how to insert the info directives, ask the
CHARMM manager for assistance.

Each documentation file, with the extension .doc, should contain

1) One brief paragraph of motivation, theory, procedure or
   whatever is necessary for a particular feature.  Here, some
   references can be given.

2) A table of contents of the documentation (to serve as the info
   menu).

3) The command syntax.

4) Complete description of all the commands and sub-commands.
   Here the syntax, defaults and file names involved would be
   described.  A brief account of what the command accomplishes
   would also be given.  The order in which various commands
   should be invoked would be described.  Relevant commands and
   subcommands can be cross-referenced with a key.

5) One or two examples involving concepts and commands described
   (No output listing).

The same notation should be followed throughout the documentation.

::

   [...]   optional,  can be present only once, if at all.

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


.. _develop_checkin:

Checkin Procedure in CHARMM Management System
---------------------------------------------

We maintain CHARMM as an integrated single source package.  The
following rules have been established to keep CHARMM as such and
to minimize time consuming problems that result from carelessness
and conflicts between developers.

1) It is always wise to inform the CHARMM manager about your
   development plan and time table so that he may arrange CHARMM
   management schedule and prevent you duplicating works done by
   others.  The list of files you are working on and the nature of
   modification should be reported in advance.  You may use the
   template form for such a report, found in support/form/project.form.

2) When you finish your project in code development, make an
   appointment with the CHARMM manager and get the most current
   version of the distribution package.  Then, he will lock the
   central CHARMM library and allow you to work on it.  Normally you
   are given a couple of weeks to incorporate your contributions into
   the main source.  It is very important to plan ahead for the
   appointment.

3) Follow the steps in :ref:`developer_modify` to integrate your modifications
   with the current source.  When you finish incorporating modifications,
   debugging, testing and documenting, prepare to send your final
   version in CHARMM management system.  (a) Modified/added files
   including source, testcase input, documentation and others and
   (b) the change-log file describing the modifications in detail are
   required for your checkin.

4) You and the manager will work together to check the modified
   files in the central CHARMM system.  After the successful
   incorporation, the change-log file will be mailed to all CHARMM 
   developers (charmm-bugs@tammy.harvard.edu) and the new developmental
   version will be established.

