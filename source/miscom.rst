.. py:module:: miscom

######################
Miscellaneous Commands
######################

The commands described in this section are generally more
simple in nature than those of previous sections. Some are perhaps
obsolete, but included for the sake of completeness.

.. index:: misc; syntax
.. _miscom_syntax:

Syntax of miscellaneous commands
--------------------------------

* File handling:

  :: 

   OPEN    UNIT integer NAME filename [WRITe ]   [UNFORMatted]
                                      [READ  ]   [FILE]
                                      [APPEnd]   [FORMatted]
                                                 [CARD]

   LOWEr                                 ! Force the case of output file names
   UPPEr                                 !             "

   CLOSe   UNIT integer  [DISPosition KEEP  ]
                         [DISPosition DELEte]

   REWInd  UNIT integer

   REWFalse ! do not allow rewinds during trajectory I/O 
   REWTrue  ! allow rewinds during trajectory I/O

   INQUire ! get a list of open files and their qualifiers, only from CHARMM
           ! possible

   STREam  [ UNIT integer       ] [ repeat(argument) ]
           [ file_specification ]        ! Call another input file


   OUTUnit  integer                      ! Redirect output to a different unit.

   RETUrn                                ! Return to the previous unit

---------------------------------------------------------------------------

   ::

      DEFIne  keyname  SELE atom_selection END

      SYSTem  "shell-cmd"  ! Execute a shell command (Unix)

      STOP               ! Terminate CHARMM

      USER ....          ! Invoke a user supplied subroutine (USERSB)

      PREF               ! Print out pref.dat keywords that were used during
                         !     compilation of running executable

      [ECHU [integer]]   ! Print remainder of commandline to OUTU or to file open on 
      ECHO ....          ! the unit specified with ECHU integer. Just ECHU resets
                         ! echo-unit to OUTU.
                      
---------------------------------------------------------------------------

* Title manipulation:

  ::

   TITLe  [COPY]           ! Specify the main "write" title

   WRITe TITLe [UNIT int]  ! output a title to the specified unit without
                           ! closing the file and without initial "*"s.
                           ! ***** (only from main command parser)

---------------------------------------------------------------------------

* Control levels:

  ::

   TIME  {integer}      ! Specify the timing level for performance evaluation
         {NOW    }      ! Show current time.
         {DIFF   }      ! Show elapsed time since last timed event

   DATE                 ! Display the current date and time

   WRNLev int [NODE int]! Set the warning print level. Higher values mean
                        ! more warnings printed (default: 5).
                        ! If running in parallel, only node 0 is modified
                        ! unless another node is specified.

   PRNLev int [NODE int]! Set the print level. Higher values mean
                        ! more printout (default: 5) (default: all nodes).

   IOLEv int [NODE int] ! Set the I/O level (not for normal use)

   BOMBlev integer      ! Set the error termination level
   NOBOmb               ! Don't bomb on typing errors (same as BOMBlev -1).

   FASTer  [integer]    ! Specify efficiency level
           [OFF    ]    ! Disable fast routines
           [DEFAult]    ! Use fast routines if possible.
           [ON     ]    ! Use fast routines, otherwise error.
           [SCALar ]    ! Use scalar fast routines.
           [VECTor ]    ! Use vector fast routines.
           [VPAR   ]    ! Use vector/parallel fast routines.
           [CRAYvec]    ! Use vectorized CRAY fast routines.

   LONG                 ! specify long line output (<256 characters)
   SHORT                ! specify short line output (<80 characters)

.. index:: miscom; quick

---------------------------------------------------------------------------

* Quick and simple structure analysis:

  ::

   QUICK { repeat(atom-spec [COMP] ) } ! one atom     - position and projection
   Q                                   ! two atoms    - distance
                                       ! three atoms  - angle
                                       ! four atoms   - dihedral
                                       ! five or more - list positions only

         atom-spec::= { residue-number atom-name  }
                      { segid  resid atom-name    }
                      { BYNUm  atom-number        }
                      { atom-selection [MASS]     }
                      { atom-number ***           }


  If only one atom is specified its position will be printed as well as
  its relationship to the previously defined axis (if any).
  (e.g. :ref:`COOR AXIS <corman_axis>` command).
  
  If the keyword ":chm:`COMP`" immediately follows an atom specification
  (or atom selection), then the comparison coordinate value(s) will be used
  for that atom only.
  
  If atom selections involving multiple atoms are specified, the center
  of geometry or center of mass of each atom selection will be used as the
  coordinate for the analysis. Note that if mass weighting is used, the
  keyword MASS must immediately follow the associated atom selection.
  
  .. note::
     This is the old syntax.  It may be used only if ALL atoms are specified
     in this manner (simple integers) and no :chm:`COMP` feature is allowed.
  
  The QUICK command sets the following substitution paramters
  (for use subsequent commands);
  
  ::
  
      one   atom  specified - @XVAL, @YVAL, @ZVAL
      two   atoms specified - @DIST
      three atoms specified - @THET
      four  atoms specified - @PHI
  
  Some examples:
  
  ::
  
     --- bond distance using atom selections ------------
            quick sele atom aseg 53 HN end sele atom aseg 53 N end
     --- Angle using atom selections ------------
            quick sele atom aseg 53 HN end -
                  sele atom aseg 53 N end  -
                  sele atom aseg 53 CA  end 
     --- Dihedral using atom selections ------------
            quick sele segi buta .and. type C4 end  -
                  sele segi buta .and. type C3 end -
                  sele segi buta .and. type C2 end -
                  sele segi buta .and. type C1 end
     --- Dihedral using segid/resid/atom ------------
            quick buta 1 C4  buta 1 C3 buta 1 C2 buta 1 C1 end
     --- Using simple atom numbers -----
            q 1 2 4   ! bond angle involing atoms 1-2-4
            q 1 2 4 6 ! dihedral involving atoms 1-2-4-6
     --- Using a mixture of formats -----
            q 1 CL 2 N sele ires 4 end MASS B 1 C COMP
               ! dihedral involving atoms:
               ! ires 1 CL -- ires 2 N -- center of mass of ires 4 --
               ! and the comparison coordinate value of atom B 1 C.
     --- distance between center of mass of two segments-----
            quick sele segid A end MASS sele segid B end MASS
     --- distance between an atom and its comparison coordinate value-----
            quick sele atom aseg 53 HN end sele atom aseg 53 HN end COMP
 
---------------------------------------------------------------------------

* RANDom and IRANdom specifications:

  ::

   RANDom  OLDRandom
           CLCG
           UNIForm         [SCALe scale]  [OFFSet offset]  [ASIN]   [ISEEd  iseed]
           GAUSsian sigma                                  [ACOS] 

   IRANdom                 [SERIes int]  [SETUp]  [BEGInt int]  [ENDInt int]  
                           [SEED int] 

---------------------------------------------------------------------------

* Run control:

  ::

   Command line sustitutions:

   SET parameter string                      ! Define a parameter

   CALC parameter arithmetic_expression      ! Evaluate an arithmetic expression

   command  ........ @parameter ........     ! use a parameter in a command
   command  ........ @?parameter ........    ! existance of parameter 

   command  ........ ?energy-term ........   ! use an energy value in a command

   command  ........ ?corman-value ........  ! use a corman value in a command

   SHOW [BUILtins]                           ! list all "?" substitution values.
   SHOW PARAmeters [VERBose]                 ! list contents of parameter table

   IF [parameter] [ EQ ] [ string] [THEN] command ! process a conditional
      [string** ] [ NE ]                     ! (**= single character not allowed
                  [.EQ.]                        unless from @ or ? variables)
                  [.NE.]
                                             !
   IF [parameter] [ GT ] [ value ] [THEN] command ! process a conditional
      [value**  ] [ LT ]                     ! (**= single character not allowed
                  [ GE ]                        unless from @ or ? variables)
                  [ LE ]
                  [ AE ]                     ! AE = almost equal (diff<0.0001)
                  [.GT.]
                  [.LT.]
                  [.GE.]
                  [.LE.]
                  [.AE.]
   The following forms are also allowed, and may be nested
   IF ... THEN
   statements
   ENDIF

   IF ... THEN
   statements
   ELSE
   statements
   ENDIF

   GOTO label                                ! A branching command

   LABEL label                               ! Label (up to 20 characters)
                                               that may be branched to

   INCRement  parameter [ BY value ]         ! Do an addition

   DECRement  parameter [ BY value ]         ! Do a subtraction

   GET        parameter UNIT int             ! read a parameter string

   FORMat  [ (format_spec) ]                 ! Specify a format for encoding.

   TRIM  parameter [ FROM integer ] [ TO integer ]  ! Take a substring

---------------------------------------------------------------------------

  :: 
  
   MMQM [atom-selection] [UNIT integer] [NCHAr integer]
                                        ! Write selected QM atoms together
   GAUSSIAN_HEADER                      ! with the rest of atoms as charges
   <gaussian commands>                  ! as input to GAUSSIAN program
   END
   GAUSSIAN_BASIS
   <optional gaussian general basis set specification or other input>
   END

  NCHAr specifies the number of characters of the atom type that will be output. 
  The default is one (NCHAr=1) such that, for example, for that atom type HG1, 
  only the character H will be printed in the output file.

  If CHARMM is compiled with Q-Chem then MMQM is slightly modified to function
  as a Q-Chem input writer instead of Gaussian. The modified routine should be
  called in the following way.
  
  ::

   MMQM [atom-selection] [UNIT integer]
   $rem section
   $molecule section but do not $end it
   QCHEM_MOLECULE
   $end for $molecule section
   QCHEM_MISC
   add any additional Q-Chem input sections
   END

---------------------------------------------------------------------------

* DEADline commands:

  :: 
  
   DEADline [CPU real] [CLOCk real]                 ! Time limits for job

   [SYNTAX ATLImit]

   ATLIimit alternate_command                       ! Execute if limits reached

---------------------------------------------------------------------------

  ::
  
   For assignment:
   parameter::= string containing alphanumeric or non-alphanumeric characters
   (no white-space (blanks or tabs)
   For substitution:
   parameter::= string-containing- alphanumeric-characters
   parameter::= {string containing lphanumeric or non-alphanumeric characters}
        
   energy-term::= see *note eterm:(chmdoc/energy.doc)Skipe.
   
---------------------------------------------------------------------------

* Convex ONLY:

  ::

   SPECIfy  specify-keywords

     specify-keywords ::=
                         PARAllel [NCPU integer-number-of-cpus] |
                         FLUSh |
                         NOFLush |
                         NBFActor  real-nonbond-memory-factor |
                         FNBL { ON | OFF }



Purpose of the various miscellaneous commands
---------------------------------------------

1) The OPEN command is used to open logical units to specific files specified
   from the input file rather than logical name assignments made prior
   to the run.  This is the recommended procedure to access a file
   within the program.  OPEN can be used to redirect the output that
   appears on unit 6 to different files by opening unit 6 in the middle
   of a run. The APPEnd keyword causes output to be appended to the
   output file; useful if you want to get back to your normal output
   file without sacrificing the first part of it.
   
   The case of filenames opened for WRITE access may be specified with
   the LOWEr or UPPEr commands.

2) The CLOSe command closes a logical unit.  This frees the associated file
   and logical unit so that they can be used for other purposes.  The
   default disposition of the file is KEEP.

3) The REWInd command

   The REWInd command causes the requested logical unit to
   be rewound. When used with the STREam command, a particular sequence can
   be used more than once.

4) The STREam command

   The steam command allows the input of command sequence
   to be shifted to another file. This is useful when parts of an
   input file are to be used many times or used by many different
   calculations. The only input value is the unit number to transfer to.
   In place of a unit number, a file may be specified. Stream files
   must be card format and should begin with a title.
   
   Arguments may be set by the stream command.  Arguments must
   not contain any blanks (or other delimiting characters).  They
   are assigned to the variable IN1, :/RAIN2, IN3, etc..  The command;

   ::
   
      STREam filename  arg1  arg2  arg3  arg4

   is functionally equivalent to;

   ::
   
      SET IN1 arg1
      SET IN2 arg2
      SET IN3 arg3
      SET IN4 arg4
      STREam filename

   This simplifies the use of passed parameters to a stream file.

5) The RETUrn command

   The return command causes the input of command sequence
   to return to the stream that called the current stream. Streams
   may be nested to up to 20 calls. There are no parameters for this command

---------------------------------------------------------------------------

6) The DEFIne command

   This command allows the user to specify selection keywords.
   This command must contain a keyword and an atom selection. The
   keyword may then be used in subsequent atom selections.  The keywords
   may not be abbreviated.

7) The SYSTem command

   Allows shell commands (sh, csh, ksh, etc.) to be executed from
   within CHARMM. NOTE: CHARMM assumes all shell commands are protected
   by double quotes, e.g., system "awk -f file.awk p1=1 p2=3 filename >
   filout" will process the file filename using the awk script file.awk
   and the parameters p1=1, p2=3 placing output in file fileout.

8) The STOP command

   The STOP command causes the program to terminate and to
   ignore all command that follow this command. This is useful for
   making temporary modifications to input files.
   
   .. note::
   
      This command is only available from the main program.

9)  The USER command, see :ref:`Interface <usage_interface>`.

10) The PREF command will prting out the pref.dat keywords that were
    used in the current executable. The purpose is to allow the user
    to probe the executable about whether the feature(s) that are 
    desired were in fact compiled into the executable and whether
    one can expect certain features to work. Currently these are the
    keywords that will be checked for and printed. If keywords other
    than the following were used, they will not be detected or printed.

    ::
    
      ACE ADUMB AIX370 ALLIANT ALPHA ALPHAMP AMBER APOLLO ARDENT
      ASPENER BANBA BLOCK BUFFERED CADPAC CFF CHARMMRATE CM5 CMPI
      COMMEASURE CONCURR CONVEX CRAY CRAYVEC CRAY_1DFFT CSPP
      DEBUG DEBUGGB DELTA DIMB DMCONS DOCK EISPACK ETHER
      FILEINPUT FILEOUTPUT FMA FOURD FSSHK GAMESS GBBLCK GBFIXAT
      GBINLINE GBNOLIST GBSWIT GENBORN GENCOMM GENETIC GLDISPLAY
      GNU GRAPE GWS HMCM HPUX IBM IBMRS IBMSP IBMVM IMCUBES INTEL
      IPRESS IRIS JUNK LARGE LATTICE LDLAN LDM LDMGEN LMC
      LONEPAIR LONGLINE LRST MANYNODES MBOND MC MCSS MMFF MOLVIB
      MPI MTS MULTCAN NEWTIMER NIH NOCORREL NODISPLAY NOGRAPHICS
      NOIMAGES NOLDMUP NOMISC NOPARASWAP NOST2 NOVIBRAN NO_BYCC
      NO_BYCU NO_DQS OLDDYN OS2 OTHERPARSHK PARAFULL PARALLEL
      PARALLELSHK PARASCAL PARVECT PATHINT PBC PBCUBES PBEQ
      PBEWALD PBOUND PERT PM1 PMEPLSMA PNOE POINTER_KEYWORD POLAR
      POSIX PREFMSI PRIMSH PVM PVMC QBLOCK QUANTA QUANTUM REDUCE
      REPDEB REPLICA RGYCONS RISM RXNCOR SAVEFCM SCALAR SCHED
      SGIMP SGMD SHAPES SHMEM SINGLE SOCKET SOFTVDW T3D T3E TERRA
      TIMESTAMP TNPACK TRAVEL TSM UNICOS UNIX UNUSED VAX VECTOR
      XDISPLAY XLARGE XSMALL YAMMP
      
    See also :ref:`subst` for variable substitutions for detecting
    the keywords.

---------------------------------------------------------------------------

11) The TITle command is used to modify TITLEA which is used whenever
    a file is written. This title is normally filled only in the
    CHARMM startup procedure. If the COPY keyword is used, then
    the TITLEB (the title from the most recently read file) is
    copied to TITLEA. Otherwise, a valid title specification should
    follow this command.

---------------------------------------------------------------------------

12) The TIMEr command sets the value of TIMER in COMMON /TIMER/ to the
    specified value.  This variable is used to time different functions
    in the program.

    - 1 will print out the time to evaluate ENERGY.
    - 2 will print out individual component times in ENERGY, and
      the times for various components of the EXEL nbonds update.

13) The WRNLEV command sets the value of the WRNLEV variable in
    COMMON /TIMER/ to the specified value. This is used in WRNDIE
    and elsewhere. Suggested values        for future use:

    ============ ==========================================================
    -5,5         warnings associated with fatal errors (see BOMBlev).
    5            default should print brief warning and error messages
                 for conditions that will affect outcome.
    6            more extensive information on errors and some information
                 on normal partial results and conditions
    7            verbose error messages and more normal processing
                 information for debugging
    8            all information that might be relevent to an error condition
                 plus checking results
    9,10         debugging levels for anything you might concievably want.
    10 or higher for term by term outputs from energy routines, or
                 other tasks where huge amounts of data useful only in
                 debugging might be generated.
    ============ ==========================================================

14) The BOMBlev command sets the level which determines the types of
    errors which will terminate the program. The default is zero.
    A value of -1 is suggested for interactive use. Suggested values are;

    ============ ===========================================================
      -5,-4      Limit exceeded type of errors. Run only as debug.
      -3,-2      Severe errors where results will be incorrect if continued.
      -1         Moderately severe errors, results may be bad.
      0          Parsing type errors. Some important warnings.
      1,2        Serious warnings.
      3,4,5      Assorted minor warnings (see WARNlev for their suppression).
    ============ ===========================================================

15) The FASTer command controls when and whic fast energy routines will be used.

    Certain conditions must be met in order to use the fast routines.  If
    the fast routines are requested and cannot be used, an error message will
    be issued and the slow routines will be substituted.  Also, there is less
    error checking for the fast routines. See :ref:`fast <energy_fast>`.

---------------------------------------------------------------------------

16) The SET command sets up a command line parameter.  The command line
    parameters will be substituted into the command line by the
    command line reader when it encounters the symbol "@".

    A command line parameter token can now be a string rather than just one of the
    single characters 0-9,a-z,A-Z. For substitution a token is indicated by the use
    of the @ character as before.  Arrays can be made by preceding the array
    indices with '@@', e.g. @segid@@j can be used to loop over parameter tokens
    segid1, segid2, ...
    (Note: Pete Steinbach's precursor to PARSUB, called to first substitute
    parameters preceded by '@@'.  Allows parameters to reference array elements.
    Mar 20, 1998)

    The token is end-delimited by any
    non-alphanumeric character. In the case that the token is not found in the
    parameter table, a check is made to see if the first character of the token is
    itself a token in the parameter table. If this single character token is in the
    table, the corresponding value is substituted -- this is the necessary scheme
    to allow backwards compatibility with the old parameter substitution, which
    allowed parameters embedded in strings.  For unambiguous token detection,
    "protect" the token with brackets {} --- this allows for the use of non
    alphanumerics in tokens such as -,_.  To test whether a token is in the
    parameter table, use @?token.  This will substitute 1 if token is in the table,
    0 if not. This is useful (in conjunction with the IF command) for setting
    defaults. (Note that @? takes precedence over any of the built-in parameters
    such as ?ENER etc. --- it is parsed first).
    
    :: 

      SET outfile = myjob 
      OPEN UNIT 1 WRITE CARD NAME @outfile.dat

    In the above example the token is delimited by the "." in the filename
    and the value "myjob" is substituted in place of "@outfile", resulting
    in an unit 1 being attached to the file "myjob.dat".
    To protect a token from surrounding alphanumerics, use brackets, 
    
    :: 
    
       OPEN UNIT 1 WRITE CARD NAME @{outfile}today.dat
       File name becomes "myjobtoday.dat".

    The token is taken to be whatever is delimited by the brackets  --- thus
    the token may in this case may also contain non-alphanumerics.

    ::

       SET max-temp = 500.
       DYNA VERLET FINALT = @{max-temp} ... etc...

    For backwards compatibility, get token, check in table, if not present,
    then drop back to first character of token and check again. 
    Substitute appropriately.

    :: 
    
       SET 1 rdie
       OPEN UNIT 1 WRITE CARD NAME @15.dat
       
    will result in a file named rdie5.dat

    To test the presence of a token in the parameter table use the @? operator.
    If the token is present, the value substituted is 1, if not 0.
    This is useful for setting defaults:
    
    ::
    
       if @?{max-temp} .eq. 0 set max-temp 300. 

    At present the parameter table is dimensioned as follows:
    
    ============================ ===
    Maximum number of parameters 256
    Maximum token length          32
    Maximum value length         128
    ============================ ===
    
    For current sizes use command SHOW PARAmeters VERBose (see below).
    
---------------------------------------------------------------------------

17) The SHOW command prints the available command line substitution 
    parameters.

    SHOW by itself or with BUILtin keyword prints the parameters set internally
    by the program functions, such as ?ENER, ?RMS etc.
    SHOW PARAmeters lists the user defined @ command substitution parameter table.
    The VERBose keyword prints table limits on string sizes for tokens and values.


18) The IF command will optionally execute a command based on the
    value of the parameter used. Example;
    
    ::
    
       IF 1 GT 25.0  PRINT COOR
       
    The "EQ" and "NE" operations only compare strings. Thus the string
    "2.00" would not be equal to "2.0" with these conditions. The options
    requesting a value, do a value comparison.
    The AE option will test if two values are almost equal (difference
    less than 0.0001). This avoids the problem of round off error in
    loop counters (i.e. values like 3.999999).

19) The GOTO command will rewind the current input stream and search for
    the requested label. For the sake of efficiency, frequent use of this
    command (i.e. looping) should not be used with long input files.

20) The LABEL command does nothing except mark the presence of a label
    (up to 20 characters in length) to be used by the GOTO command.

---------------------------------------------------------------------------

21) The INCRement command will modify the selected parameter. If a value
    is not specified, then a value of 1.0 will be used. Example
    
    ::
    
       INCR 1 BY 2.0

22) The DECRement command is identical to the INCRement command except
    that a subtraction is done. The purpose of this command is to allow
    the subtraction of parameters. For example, the sequence;
    
    ::
    
        SET 1 ?ENER
        DECR 1 BY ?HARM
        WRITE TITLE UNIT 30
        * @1
        *
        
    will compute the total energy less the constraint energy and write it
    to a file.

23) The FORMat command allows the user to specify the format for
    ALL subsequent calls to ENCODF. This can be used to format the output
    of titles or other internal strings. Here are some examples;
    
    ================  =========================================================
    FORMat (I5)       - All values will be integers. Good for looping and such.
    FORMat (F12.4)    - Just what it says.
    FORMat            - Reverts to current scheme for ENCODF (1PG14.6) followed
                      by trimmimg
    FORMat (A12)      - Won't work...
    ================  =========================================================

    If an integer format is used, the real value will be rounded to the
    nearest integer. The parenthesis are required around the format specified.
    If several different formats are needed, then the FORMat command should
    precede each different required usage.
    
    .. note::
   
      Not all string manipulation commands call ENCODF.  The SET command
      does not.  The INCRement command does, so the sequence;
      
      ::
      
          FORMat   (f10.5)     ! specify the format
          INCRement  a  by 0.0 ! apply the format to variable "a"
          
      may be used to format a particular variable without modifying its value.
---------------------------------------------------------------------------

24) The TRIM command allows a substring of a parameter to replace the
    same parameter. The FROM value determines the first character to be kept
    (default first nonblank character), and the TO value determines the last
    character to be kept (last nonblank character). If a TO value that is larger
    than the length of the current parameter is used, blanks will be padded at
    the end.
    
    Preceding blanks may be added by;
    
    ::
    
          SET 5           ! set parameter five to the null string
          TRIM 5 to 10    ! convert parameter five to a string with 10 blanks
          SET 6 @5@6      ! add these 10 blanks to parameter six

    This command may be used for general formatting.

---------------------------------------------------------------------------

25) The DEADline command sets CPU and/or clock-time limits. These
    limits are checked in DCNTRL,ECNTRL, and GAUSHS (the parameter-fitting
    routine) at regular intervals. When a deadline has been reached the
    routine exits normally. This is useful when you have to stop computing
    before a given time of day (taking advantage of lower charge during the night
    or some such) or when you want to get some useful results and you are not
    sure that you can actually stay within the CPUlimit in a given batch queue.

    Keyword CPU <real> specifies that  <real> CPUminutes from the
    time the command is given is to be one deadline.

    Keyword CLOCk <real> sets the time HH.MM (in 24-hour format) as
    one deadline. The routine assumes that if the command is issued after
    the specified time, you mean the following day. (If at 6 pm you start
    a job containing the line DEAD CLOC 13.00 CPU 600. your minimization
    will run until 600 CPU-minutes have been used, or until 1 pm the next
    day, whichever comes first.)

26) The ATLImit command can be given at any point in the input file.
    CHARMM checks before reading each command if either of the DEADlines
    (CPU or CLOCk) has been reached. If this is the case the alternate_command
    of the most recent ATLImit command is executed. This would typically be
    a GOTO SHUTdown or some other simple thing, but could be any CHARMM command.
    Currently the alternate_command is limited to 80 characters.

27) Substitutions and punctuation in command input.

    +---------+----------------------------------------------------------+
    |"!"      |Ignore this and all subsequent characters on this line    |
    +---------+----------------------------------------------------------+
    |"-"      |If this is the last character of a line then the following|
    |         |line is a continuation                                    |
    +---------+----------------------------------------------------------+
    |"*"      |As a first character indicates a title line. Alone on a   |
    |         |line indicates a title terminator.                        |
    +---------+----------------------------------------------------------+
    |"$"      |The default delimiter                                     |
    +---------+----------------------------------------------------------+
    |"* % # +"|Atom selection wildcards, alone or in a word              |
    |         |  === ====================================================|
    |         |   \* matches any string of characters (including none),  |
    |         |   %  matches any single character,                       |
    |         |   #  matches any string of digits (including none),      |
    |         |   \+ matches any single digit.                           |
    |         |  === ====================================================|
    +---------+----------------------------------------------------------+
    |"@"      | Command parameter substitution                           |
    +---------+----------------------------------------------------------+
    |"?"      | Energy value substitution                                |
    +---------+----------------------------------------------------------+

---------------------------------------------------------------------------

28) File inquiry.  The inquiry command (from CHARMM) may be used to
    get a list of currently open files.  This is very useful in interactive
    sessions when one has forgotten which FORTRAN units are already assigned.
    The command won't work if the files are assigned outside of CHARMM.

29) Random number generation.  

    1) RANDom command.   The expression  ?RAND  will have a random number
       substituted for it during command line evaluation.  The default is to provide
       a number from a uniform distribution, between 0.0 and 1.0; the RANDom command
       allows modification of the distribution type and specification of other factors.
       The only required keyword is the distribution type, which must be second; for a
       GAUSsian distribution, a value for sigma is required; the default mean is 0.0.

       :: 
       
         RANDom  OLDRandom
         CLCG
         UNIForm         [SCALe scale]  [OFFSet offset]  [ASIN]   [ISEEd  iseed]
         GAUSsian sigma                                  [ACOS] 

       Additional keywords:

       ================  ========================================================
       SCALe  scale      multiply the number by scale
       OFFSet offset     add offset to the number
       ACOS              treat the number as a cosine and return the angle (deg)
       ASIN              treat the number as a sine and return the angle (deg)
       ISEEd iseed       specify a new random seed (integer)
       ================  ========================================================

       Examples:
       
       ::

          RANDOM GAUSS 0.2 SCALE 10.0   !     gaussian  mean of 0.0 with a sigma of 2.
          RANDOM UNIFORM SCALE 360.     !     uniform   0. to 360
          RANDOM UNIFORM ACOS SCALE .5  !     uniform   angles with cosines from 0. to .5
          RAND GAUS 5. OFFS 60.         !     gaussian  mean of 60. with a sigma of 5.
          RAND UNIF ISEED 7734          !     uniform   new random seed

      Subsequent use of ?RAND will substitute a number from the appropriate
      distribution. 

      Note that OLDRandom subcommand sets OLDRNG, which runs "old" random
      number generator instead of "new" CLCG method. CLCG unsets OLDRNG,
      and runs the CLCG random number generator.

   2) IRANdom command.  This command is designed to generate series of random
      integers taken from uniform distributions between user-specified limits.
      Each series or distribution must first be set up with the IRANdom SETUp
      command, in which the lower and upper limits of the distribution, the series
      number, and an integer seed are specified. E.g.

      ::
      
         IRAND SERIES 1 SETUp BEGI 1 ENDI 18  SEED 2346
         IRAND SERIES 2 SETUp BEGI 1 ENDI 402 SEED 4028987

      The random integers for each series are then generated with the commands
      
      ::
      
         IRAND SERIES 1
         IRAND SERIES 2

      The ?iran expression accesses the last random integer generated. 

      The purpose of the multiple series feature is at least two-fold. First, it
      allows users to generate random numbers easily from many different
      distributions during the same CHARMM run (e.g. for use in different parts
      of the same calculation). Second, it may help the user avoid correlations
      between random numbers generated for different parts of a calculation.
      An internal counter, corresponding initially to the seed, is incremented
      by several units with each instance of the IRANdom command; by separating
      the seeds of the various distributions sufficiently, the user can thus avoid
      cross-series correlations.  The use of multiple seeds for a given series
      should be unnecessary and is discouraged.  The IRANdom function has an
      overall period of no less than 10^12 for distribution widths of 10^10 or less.
      IRANdom can also be used to effectively generate random real numbers, through
      a division of the generated integers by a constant, with the use of the
      CALC command.

---------------------------------------------------------------------------

30) The CALC command allows the evaluation of any fortran-admissible
    arithmetic expression.  It supports most of the normal fortran functions
    such as COS, SIN, TAN, EXP, LN, LOG, TANH, etc...  Any number of parenthesis
    nesting is allowed.  The substitution parameters @ is allowed directly.  The
    substitution parameters ? can also be used but the character chain must
    be surrounded by blanks to be properly recognized; e.g., COS( ?pi ) is
    ok but not COS(?pi).  Otherwise, there can be any number of blanks
    between the quantities involved in the arithmetic expression.  See the
    testcase calc.inp for examples.

    .. note::
    
       All transcedental functions work in natural units (not degrees).

---------------------------------------------------------------------------

31) Writing input for GAUSSIAN series of programs. Selected atoms are
    treated as quantum atoms while the rest of the system is put at the
    end of the file in a format ready for CHARGE command within GAUSSIAN
    (must be at least version 92) Gaussian commands are specified after
    GAUSSIAN_HEADER keyword ended by the END keyword, and other input is
    optionally specified after GAUSSIAN_BASIS keyword. If none of the two
    is specified both END keywords must still be present. There is no
    check for the names of atoms not specified according to periodic table
    of elements.  Use RENAme ATOM command to rename CA atoms for
    example. Charges are taken from RTF.

    ::
    
      MMQM [atom-selection] [UNIT integer]
      GAUSSIAN_HEADER
      # 6-31g** charge scf=direct mp2=fulldirect gen
      END
      GAUSSIAN_BASIS
      <optional gaussian general basis set specification or other input>
      END

   
   a) Example using MMQM to write Q-Chem input file:

      ::
      
         MMQM [atom-selection] [UNIT integer]
         $rem
         .... add all rem variables ....
         $end

         $molecule
         0 1 
         QCHEM_MOLECULE
         $end

         QCHEM_MISC
         .... add additional Q-Chem input sections ....
         this is the place to specify an $opt section and add
         constraints for performin entire PES scans
         END

---------------------------------------------------------------------------

32) SPECIfy  specify-keywords             !  Convex ONLY 

    ::
    
       specify-keywords ::=
                            PARAllel [NCPU integer-number-of-cpus] |
                            FLUSh |
                            NOFLush |
                            NBFActor  real-nonbond-memory-factor |
                            FNBL { ON | OFF }

    description:
    
    1. PARAllel - Tells CHARMm to run parallel (where possible). The optional
       NCPU keyword specifies the maximum number of processors to use.  If
       a number is specified that is greater than the maximum allowed for the
       particular machine, a warning message is printed and the number of cpu's
       is set to the maximum.  Note that at startup CHARMm senses the number of
       cpu's and sets NCPU accordingly.

    2. FLUSh -  Specifies the that trajectory; coordinate; dynamics restart
       and other output files should be flushed after each data set is written.
       See below.  This is the default action. The command is provided to reset
    3. NBFActor - When the parallel non-bond list generators allocate
       memory for the temporary arrays used by each thread, the predicted size
       of list array (MXJNB and the like), is divided by the number of cpu's
       and multiplied by NBFACT. The default is 1.5 and has worked well so far.
       If it doesn't the SPECIfy NBFACT <num> command is available to adjust it.
    4. FNBL - FastNonBondListgeneration -  Specifies whether or not to use
       the new non-bond list generation routines. Just included for testing
       and timing purposes.

------------------------------------------------------------------------------

33) ``IOFOrmat [ EXTEned | NOEXtended ]``

    In c30a2, atom numbers can assoume I10 and PSF IDs (SEGID, RESID,
    RES and TYPE) can be character*8.  Atom numbers take I5 in coordinate
    files and I8 in psf files and CHARACTER*4 PSF IDs are used for Normal
    (noextended) I/O operation.  These are expanded to I10 and A8
    respectively.  Noextended format is the default and the expanded
    format is used only when the number of atoms is greater than 100000 or
    any PSF ID is longer than 4 characters.  This command overrides the
    default set:  IOFOrmat EXTEnded enforces the extended format and
    IOFOrmat NOEXtended does the normal (old) format.
