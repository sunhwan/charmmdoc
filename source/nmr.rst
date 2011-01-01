.. py:module::nmr

===================
NMR Analysis Module   
===================

::

   CHARMM Element doc/nmr.doc 1.1
   
   File: NMR, Node: Top, Up: (chmdoc/commands.doc), Next: Syntax


                            NMR Analysis Module   

       The NMR commands may be used to obtain a set of time series for a
   number of NMR properties from a trajectory.  Among the possible
   properties are relaxation rates due to dipole-dipole fluctuations (T1,
   T2, NOE, ROE), chemical shift anisotropy and Deuterium order
   parameters for oriented samples.  The documentation assumes that users
   are already familiar with NMR.  Several textbooks are available for
   users interested in more information.  The NMR command invokes the NMR
   subcommand parser.

       Because several properties are based uppon the position of nuclei
   that may not have been included in the PSF (and the trajectory) the
   module has its own building submodule (see BUILD) to construct atoms.
   For example, the H_alpha on the C_alpha can be constructed without
   invoking HBUILD for T1 and T2 calculations.  

       Everthing is stored on the HEAP and no variables are kept when the
   module is left (there is no nmr.fcm common block).  Everything is
   re-initialized when the module is exited with the END command.


   WARNING: The module has not been used in numerous situations and caution
            should be the rule. In case of doubt it is best to study the
            source code. 

   * Menu:

   * Syntax::      Syntax of the NMR commands
   * Function::    Purpose of each of the commands
   * Examples::    Usage examples of the NMR analysis commands

   
   File: NMR, Node: Syntax, Up: Top, Previous: Top, Next: Function


                                    Syntax

   [SYNTAX NMR functions]

   Syntax:

   NMR   enter the NMR module
   END   exit the NMR module

   Subcommands:

   BUIL name  4 x single-atom-selection [DIST  real] [THETA  real] [DIHE  real]

   CSA        4 x single-atom-selection  -
              [THE1  real] [PHI1  real] [THE2  real] [PHI2  real] -
              [S11  real]  [S22  real] [S33  real]

   DQS        2 x atoms-selection   [CTDI real]     

   SET        [HFIELD  real]  [RTUMBL  real]  [CUT real]
              [GAMMA real [atoms-selection]]

   RTIM       2 x atoms-selection [RTUMBL real]  [HFIELD real] 
              [CSAR DSIGMA real] [CUT real] [STAT]
              {ANIS} [DRAT real] [DTSX real DTSY real DTSZ real]
           

   DYNA       traj-specification -
              [CUT real] [RTUMBL real] [HFIELD real] [TMAX real] -
              [ORIENT { MASS }   atom-selection [NOCOMP] ]
                      { WEIGH }
              [UNIT INTEGER] {quantity} 
              [DSIGma real] [ILISt integer] [MODFree integer] [MFDAta integer]
              [SAVE ]

   DYNA        WRSTAT [ILISt integer] [ISDLIST integer] 

   WRIT       { LIST }   [UNIT  INTEGER]
              { COOR }
              { CSA  }
              { DQS  }
              { RTIM }

   RESE       [KEEP]

   name::= an arbitrary name chosen by the user (four characters)
   atom-selection::= a selection of a group of atoms 
   single-atom-selection::= a selection of a single atom
   quantity::= none or any combination of C(t), R(t), PROP, UVEC, VERBOSE

   traj-specification::= FIRStu int [ NUNIt int ] 
                         [ BEGIn int ]   [ SKIP int ]  STOP int

   NB!  STOP has to be specified so that enough memory may be allocated.
   It is a good idea to specify the other parametters as well (they default to 1).
   The command TRAJ QUERY UNIT int (*note dynamc (chmdoc/dynamc.doc)) may be used
   to obtain necessary data from the trajectory files.

   
   File: NMR, Node: Function, Up: Top, Previous: Syntax, Next: Examples
 

           General discussion regarding the NMR analysis module


   1. RTIM (Relaxation TIMes)
   --------------------------

   This serves as a set-up for the DYNA command but can also be used to calculate
   the relaxation parameters for a rigid body if keyword STAT is specified.

   For the dipole-dipole relaxation rate properties the rotational tumbling
   time are taken in PS, the magnetic field in TESLA (note that 11.74 T yields
   500 MHz for a proton H).  The Relaxation Times are in 1/SEC.
   The gyromagnetic ratio of a nucleus are obtained from the first letter
   in the TYPE() TYPE array and can be modified using the SET commands.
   The default GAMMA constants are taken from table 2.1 in "NMR of Proteins
   and Nucleic Acid" by K. Wuthrich.

        Nucleus   Gamma [RADIAN/(TESLA*SEC)]:
         H        26.75D07
         C        6.73D07
         N        -2.71D07    *NB the sign is important for the spectral densities
         P        10.83D07
 
   All the time series of all particles involved in a RTIM selection are kept
   on the HEAP.  The total number of time series is indicated in the output,
   so is the HEAP required storage.  Very large sets of long trajectory can
   be broken down (one can use the repeated loops of the MISC command to
   do this).

   In the NMR analysis module, the spectral densities are defined as

                     +inf
                     /
          J(W) =     \  COS(W*t) C(t) Dt    =   J(|W|)
                     /
                     0

   following the convention of R.M. Levy et al., JACS 103, 5998 (1981), or
   E.T. Olejniczak et al., JACS 106, 1923 (1983).  Notice that this convention
   differs from other notations such as in "Principles of Nuclear Magnetic
   Resonance in One and Two Dimensions" by R.R. Ernst, G. Bodenhausen and
   A. Wokaun, Oxford 1987, where there is a factor of 2 to account for an
   integral from -inf to +inf (see section 2.3 of that reference). 

   The fast part (decays on time scale of ps) and the slow part (decays from
   the rotational diffusion with RTUMB) are integrated separatly.

                     +TMAX 
                     /                                    
   J_{fast}(W) =     \  COS(W*t) [C(t)-C_{plateau}] Exp[-t/RTUMBL] Dt  
                     /                         
                     0    


                     +inf                     
                     /                         
   J_{slow}(W) =     \  COS(W*t) C_{plateau} Exp[-t/RTUMBL] Dt
                     /                           
                     0                          

   Many papers report different formulas for T1, T2, T1R, NOE and ROE.
   See for instance Levy and M. Karplus p. 445 "Trajectory studies of NMR
   relaxation in flexible molecules", Chap 18, p. 445, American Chemical
   Society 1983. See also, I. Solomon, Phys. Rev. 99, 599 (1955).
   In the NMR module, the expressions used are:

   Spectral densities:  J(0)      J(W1)     J(W2)    J(W1-W2)       J(W1+W2)

   1/T1  = FACT*(3*J(W1)+J(W1-W2)+6*J(W1+W2))
   1/T2  = FACT*(4*J(0)+3*J(W1)+6*J(W2)+J(W1-W2)+6*J(W1+W2))/2
   1/T1R = FACT*(3*J(0)+5*J(W1)+2*J(W1+W2))  (T1 in rotating frame, to be checked)
   NOE   = 1+(GAMMA2/GAMMA1)*(-J(W1-W2)+6*J(W1+W2))/(3*J(W1)+J(W1-W2)+6*J(W1+W2))
   ROE   = FACT*(3*J(W1)+2*J(W1-W2))

   with the prefactor given by:

   FACT  = 1/10 * ((MU0/4*PI)*PLANCK/(TWO*PI)*GAMMA1*GAMMA2)**2 * (PSEC/ANGS**6)

   where PLANCK = 6.62618D-34, ANGS=1.0D-10, PSEC = 1.0D-12, MU0 = 4*PI*1D-07
   the permitivity of vacuum in SI units.  The rates are converted to [1/SEC]
   by the factor FACT. Note that the spectral density contains the distance
   dependent part <r**-3>**2.

   The order parameters are calculated from the average of 
   plateau = 3/4 <Y2/R**3>**2 + 3 <Y1/R**3>**2 + 1/4 <Y0/R*3>**2
   using COMPLEX arithmetics.

   If CSAR DSIG {real} is given the contribution of chemical shift anistropy 
   to the relaxation will also be calculated for bonds less or equal to 1 Angstrom
   length using the value DSIG for the chem. shift anisotropy. For N15 nucleus
   a value of -160 ppm is recommended and the director is approximately
   directed along the N-H bond. A unit vector can be generated using the BUILD
   facility of the NMR module if other axis are desired. The expressions for the
   chemical shift anistropy relaxation were taken from Goldman's book on NMR:

           1/T1 =  (2/15)*(1.0E-06*DSIGMA*W1)**2*J(W1)/<1/R**6>
           1/T2 =  (2/15)*(1.0E-06*DSIGMA*W1)**2*((2/3)*J(0)+(1/2)*J(W1))/<1/R**6>

   where DSIGMA is in ppm and W1 is GAMMA1*HFIELD. The distance dependence in
   J(W) is also removed here.

   The relaxation contribution due to CSA is added to give the total relaxation
   value for the spin pair in the output to ILIST file command when DSIG
   keyword is present in DYNA command.

   Keyword ANIS: Anisotropy is now implemented for an axially symmetric molecule,
   i.e. Dy ~ Dz of the principle axes of interia or diffusion tensor, so that
   Dparallel and Dperpendicular to a long axis can be used. Obtain these values
   from the relaxation data via a program such as ROTDIF (Walker O, Varadan R, 
   Fushman D. 2004. J. Magn. Reson. 168:336-345), or via hydrodynamics
   calculations. DRAT {real} is the ratio Dparalell/Dperp. of the diffusion tensor
   DTSX, DTSY, DTSZ {real} is the diffusion tensor axes, typically using orient
   will align the coordinate set with the longest axes of inertia along x
   (so it would be 1.0 0.0 0.0 but any alignment could be chosen). For alignment
   the coordinates of the comparison set are used.

   Note: Tau1,2,3 are calculated from value of DRAT and RTUMBL With keyword
   STATic {logical} the relaxation parameters will be calculated for global
   tumbling with a correlation time RTUMBL or in case of anisotropy with
   Tau1,2,3 and the structure in the comparison set.

   It should be noted that since global and internal motions are modeled
   separately, that anisotropy has no effect on the correlation functions,
   i.e. S2, but mixes into the calculation of relaxation parameters.
   Equations used follow those in Barbato et al., Biochem. 31, 5269-78 (1992)
   and further description is given in Buck et al., 2005 (submitted to JACS).


   2.  DYNA option
   ---------------

           The DYNAmics command reads in the trajectory from fortran units
   opened with sequential numbers.  *note dynamc (chmdoc/dynamc.doc)

   	FIRSTU is the unit assigned to the first file of the trajectory, 
   and must be specified.  NUNIT gives the number of units to be scanned,
   and defaults to 1.  BEGIN, STOP, and SKIP are used to specify which steps
   in the trajectory are actually used. BEGIN specifies the first step number to
   be used. STOP specifies the last. SKIP is used to select steps
   periodically as follows: only those steps whose step number is evenly
   divisible by STEP are selected. The default value for BEGIN is the first
   step in the trajectory; for STOP, it is the last step in the trajectory;
   and for SKIP, the default is 1.  A similar logic is used in the CORREL
   module (*note (chmdoc/correl.doc) ).

   Keyword CUT can be used to specify a cutoff for the distance between nuclei to
   be included in the calculation.

   	ORIE is used to reorient all coordinate frames of a trajectory with
   respect to the comparison set; if NOCOMP keyword is present, orientation
   will be wrt the first frame of the trajectory piece to be analyzed.
   This is done to obtain the internal dipole-dipole correlation functions
   in the molecular frame assuming internal motions and overall rotation are
   independent.  Overall rotation is assumed to be isotropic and to correspond
   to an exponential correlation function with a characteristic time equal
   to RTUMBL (ps). 

            HFIELD is the magnetic field strength in tesla.  Default = 11.74 Tesla
   which yields a Larmor frequency of 500 MHz for protons.  The value of TMAX
   is the maximum time used to numerically integrate the fast part of the
   internal correlation function.  A simple trapezoidal rule is used.
   The default value of TMAX is 0.0, the correlation function should be examined
   to set a reasonable value for TMAX [for instance, see R. Bruschweiler, 
   B. Roux, M. Blackledge, C. Griesinger, M. Karplus and R. Ernst.
   ``Influence of Rapid Intramolecular Motions on NMR Cross-Relaxation Rates.  
   A Molecular Dynamics Study of Antamanide in Solution'', J. am. Chem. Soc. 
   114, 2289 (1992)].

   If  RTUMB .le. 0.0  then no analytic overall rotation contribution is computed.
      This is to be used with trajectories that retain the overall diffusion.

   Output includes a rough estimate of the effective correlation time for the 
   analyzed (NH) motions, and an entropy estimate using the "diffusion in a cone"
   model (Yang&Kay,JMB263,p369 (1996) "model 3")

   DSIGma adds a CSA contribution to the relaxation rate (see also CSA below)
      The TOTAL rates (and the rates written to the ILISt file) contain this
      CSA contribution, whereas the rates printed immediately after each
      spin-system do not.  

   ILISt  specifies a file for compact writing of relaxation parameters.
     The columns are relaxation rates as defined above (in 1/sec) R1, R2, NOE,
     ROE, R2/R1, <S2>, Sconf, Taue, TMXE, and atom identifiers.
     Here <S2> is the plateau value (generalized order parameter),
     Sconf is an entropy estimate using the  diffusion-in-a-cone model
     (Yang&Kay,JMB263,p369 (1996) "model 3") neglecting alternative Sconf values
     for S2 < 1/64, and using approximation A=-0.11 as suggested by Yang&Kay.
     Taue is the effective correlation time for this motion computed from the
     integral of the correlation function C(t) out to TMXE, the first time when
     C(t) is <= <S2>.  

   MODFree and MFDAta specify files that can be used as input to Art Palmer's
         ModelFree NMR analysis program

   The SAVE keyword adds relaxation parameters for subsequent statistical
      averaging (DYNA WRSTAT)

   The output is written to UNIT.  The output level is controlled by the keywords:
   C(t)    dipole-dipole relaxation correlation functions 
   R(t)    dipole-dipole time series
   PROP    CSA and DQS for solid state NMR properties
   UVEC    unit vectors for CSA and DQS solid state NMR
   VERBOSE all quantites will be written out (including all coordinate frames!)
 
   DYNA WRSTAT is a special form of the command, which simply computes averages
   and standard-deviations of the relaxation parameters that were SAVEd in
   previous DYNA commands, and writes them out to ILIST and ISDLIST, respectively.
   Accumulators are zeroed in preparation for a new round of statistics
   collection. 

   In addition to the correlation functions, relaxation parameters are calculated
   (see above). It should be noted that spin-spin distances and anisotropy
   (specifically the angle of the vector with the long axis) are taken as the
   trajectory average. If a constant distance, e.g. 1.02A for N-H is desired
   you need to alter the source-code. 


   3. other NMR properties supported   
   ---------------------------------

   3.1 CSA (Chemical Shift Anisotropy): 

   Construct the principal axis from a z-matrix
         1          u
          \        /          theta 2-3-u   (theta=0 gives u along 2-3)
           2 --- 3*           phi   1-2-3-u (phi=0 gives a cis)
   "u" is the end of the unit vector indicating a principal axis starting from
   atom 3

   CSA = SUM_{axis_i}  S_ii (Z(i)**2 - 0.5 *(X(i)**2+Y(i)**2) )

   where X(i), Y(i), and Z(i) are the components of the i-th unit vector of the
   chemical shift tensor elements and S_ii is the magnitude of the i-th tensor 
   element.  The chemical shift tensor is a symmetric second rank tensor 
   and is determined by 3 chemical tensor elements and 3 unit vectors.  The
   value of the chemical shift parallel, Z(i)**2, and perpendicular,
   0.5*(X(i)**2+Y(i)**2, are also given independently.

   For example, the N15- chemical shift anisotropy for the peptide backbone
   has been studie by Mai W., Hu W., Wang C., and Cross TA.  
   "Orientational constraints as three-dimensional structural constraints 
   from chemical shift anisotropy: the polypeptide backbone of 
   gramicidin A in a lipid bilayer".  Protein Science (1993) Apr;2(4):532-42.

   CSA  S11 37.0 S22  62.0  S33 202.0   -
        the1   71.0  phi1 180.0 the2  -90.0 phi2 90.0 -
        select resid 2 .and. type C   end  -
        select resid 3 .and. type H   end  -
        select resid 3 .and. type N   end


   3.2 DQS (Deuterium Quadrupol Splitting):

   Construct the unit vector between a pair of atoms and project it onto the 
   reference Z-direction.

   DQS = (3*Z**2-1)/2.0, 

   where Z is the projection along the Z axis of the unit vector of a 
   carbon-deuterium bond.  This particular property could also be easily computed
   from the options of the CORREL module, *note correl: (chmdoc/correl.doc).


   4. BUILD
   --------

        The build command is useful for constructing hydrogen atoms, or
   any other particle, that is involved in the calculation of an NMR propertiy
   but is not present explicitly in the trajectory file.  An example would be
   the NMR relaxation times T1, T2 of the H_alpha, which is not included in
   the extended atom potential function (e.g., in toph19.inp).  The syntax
   is simply a Z-MATRIX input line, where the first three atoms have well-defined
   coordinates.  The name given to the new atom is arbitrary.  By default
   the RESID and RESNAM are the same as that of the first atom-selection and
   the SEGID is called "BUIL".  The atom position is stored starting from
   NATOM+1, at the end of the coordinate list.   The coordinates are
   re-built automatically before computing any NMR property.


   5. WRITE
   --------

        The WRITE command is used to write out most information.   The default
   output is used unless a UNIT number is given (that unit is not closed by
   the NMR module). The keywords LIST (write out all the list of all properties,
   mostly used for debugging), COOR (mostly to have access to the coordinates 
   constructed by the BUILD option), and the NMR properties (CSA, DQS and RTIM).
   The level of printout detail is controlled by PRNLEV (see (chmdoc/misc.doc)).
   This will change in future versions and the printout level will be controlled 
   by direct keywords.  The present levels of printout are:

     PRNLEV         OUTPUT
     0 (default)    normal output for all options and commands
     1              value of DQS, CSA for individual structure 
     2              Value of the spectral densities J(W1)
     3              Larmor Frequencies
     4              Dynamics steps, time and NCOORD
                    Fast and plateau part of the spectral densities
     5              Associated unit-vectors for CSA and DQS
                    COOR ORIENT normal output in DYNAM (angle and axis printed)
     6              Correlation function for relaxation
                    Integrand in calculations of spectral densities
     7              Spin-spin time series used to compute the correlation function
     8              Full spin trajectory


   6. SET
   ------

        The SET command is useful to enter a the value of the gyromagnetic 
   ratio GAMMA for a new type of nucleus (with the atom selection) and add it 
   to the default list of nuclei (the gyromagnetic ratio GAMMA is involved 
   in the relation OMEGA=GAMMA*HFIELD, where OMEGA is the Larmor frequency).  
   The nuclei now supported by the NMR module are: H, C, N, and P. 
   It is also possible to use the SET command to give values for RTUMBL 
   and HFIELD which are kept for the relation calculations.  


   7. RESET
   --------

            Resets all assignements of the NMR module.  Destroys all lists and
   is equivalent to exiting and re-entering the module. 


   8. Miscellaneous command manipulations
   --------------------------------------

   *note misc: (chmdoc/miscom.doc) are supported within the NMR module,
   allowing opening and closing of files, label assignments (e.g., LABEL), 
   and repeated loops (e.g., GOTO), parameter substitutions (e.g., @1, @2, etc...)
   and control (e.g., IF 1 eq 10.0 GOTO LOOP).

   
   File: NMR, Node: Examples, Up: Top, Previous: Function, Next: Top


                                   Examples

           These examples are meant to be a partial guide in setting up
   input files for NMR. The test cases may be examined for a wider
   set of applications.  There is 1 file: nmrtest1.inp which can be submitted 
   through nmrtest.com.


   Example (1)
   -----------

   NMR
   reset  

   ! Relaxation times
   ! H - N pair
   RTIMES  select type N end    select type H end 
   WRITE RTIMS rtumbl 500.0 hfield 11.74 cut 3.5 iwrite 6

   END

   Produces a verbose output of all the N-H dipole-dipole relaxation rates
   within a distance of 3.5 angstroms in the presence of a magnetic field of
   11.74 Tesla and assuming a isotropic tumbling of 500 picoseconds.
   Print out to unit 6.


   Example (2)
   -----------                                     

   NMR 
   reset
   BUILD HA1 select type CA  .and. resid 2 end  dist     1.08 -
             select type C   .and. resid 2 end  theta  109.28 -
             select type N   .and. resid 2 end  dihe  -120.00

   WRITE COOR select segid BUIL .or. resid 2 end

   END

   Build the position of hydrogen  bonded to CA #1 with ZMATRIX syntax and
   print out the coordinates to verify the structure (verification should
   always be done).  The NAME of the atom built is HA1, the RESNAM and the
   RESID are the same as those of the first selected atom, the SEGID is
   called BUIL by default.  The coordinates are added at the end of the
   structure (after NATOM).  The command ZMAT can be called from outside
   the NMR module and supplements the IC table with a "gaussian-like" zmatrix.


   Example (3)
   -----------

   NMR
   reset
   ! Phosphate group chemical shift anisotropy for lipids
   ! from J. Herzfeld et al., Biochem. 17, 2711 (1978).
   CSA  S11 -76.0   S22 -17.0  S33 110.0 -
        the1   180.00  phi1  0.0 the2  90.00  phi2 0.0 -
        select resid 1 .and. type P end -
        select resid 1 .and. type O11 end  -
        select resid 1 .and. type O12 end

   write CSA

   build HA select type C11 .and. resid 1 end  dist     1.08 -
            select type C12 .and. resid 1 end  theta  109.28 -
            select type O12 .and. resid 1 end  dihe   120.00

   build HB select type C11 .and. resid 1 end  dist     1.08 -
            select type C12 .and. resid 1 end  theta  109.28 -
            select type O12 .and. resid 1 end  dihe  -120.00

   DQS select type C* end select type H* end
   write DQS

   END

   Defines the Chemical Shift Anisotropy of a phosphate group in the phospholipids
   DPPC with the experimental principal axis values and print it.
   Construct the coordinates of two hydrogen (deuterium) and calculate
   the order parameters of the static structure.


   Example (4)
   -----------


   open read unformatted unit 50 name nmrtest1.trj
   DYNA nunit      1   firstu    50   begin   100   stop  10000   skip 100 -
        rtumbl 500.0   hfield 11.74   cut     3.5   tmax    3.0   -
        iwrite     6   C(t)           R(t)           -
        orient         select type CA end 


   Calculate the NMR properties from trajectory nmrtest1.trj re-orienting
   all the frames with respect to the carbon CA of the COMP cordinate set.
   For the relaxation correlation function integrals are cut at a TMAX of
   3.0 psec.  Write out the time series and the correlation function.


   Example (5)
   -----------


   ! build the position of chemical shift director with ZMATRIX syntax 
   build X  select type N   .and. resid 2 end  dist     1.00 -
            select type H   .and. resid 2 end  theta    0.00 -
            select type C   .and. resid 2 end  dihe     0.00

   ! build the position of chemical shift director with ZMATRIX syntax 
   build X  select type N   .and. resid 3 end  dist     1.00 -
            select type H   .and. resid 3 end  theta    0.00 -
            select type C   .and. resid 3 end  dihe     0.00


   RTIMES CSAR dsigma 160.0 rtumbl  500.0 hfield 11.74  -
               select type N end  select type X end  

   open read unformatted unit 50 name nmrtest1.trj
   DYNA nunit      1   firstu    50   begin   100   stop  10000   skip 100 -
        rtumbl 500.0   hfield 11.74   cut     3.5   tmax    3.0   -
        iwrite     6   -
        orient         select type CA end 


   Defines fictitious unit vectors with the build facility and calculate
   the chemical shift anisotropy relaxation for N15.  The anisotropy is
   about 160 ppm between the principal axis if a near cylindrical symmetry
   is assumed.


   Example (6)
   -----------
   {see also test/c33test/nmrtest2.inp}

   RTIMES STAT CSAR DSIG 170.0 rtumbl  500.0 hfield 11.74  -
         ANIS DTSX 1.0 DTSY 0.0 DTSZ 0.0 DRAT 1.2 CUT 2.3 -
               select type N end  select type X end  

   Calculates the relaxation parameters for mainchain N-H spin pairs assuming
   a rigid molecule (coordinates in the comparison set) tumbling as a symmetric
   top with the long axis aligned along x (thus DTSX,y,z are 1,0,0) and
   a Dparallel/Deper ratio of 1.2. Dipole-Dipole and CSA contributions are
   calculated


   Example (7)
   -----------

   RTIMES rtumbl  500.0 hfield 11.74  -
         ANIS DTSX 1.0 DTSY 0.0 DTSZ 0.0 DRAT 1.2 CUT 2.3 -
               select type N end  select type X end  

   open read unformatted unit 50 name nmrtest1.trj
   DYNA nunit      1   firstu    50   begin   100   stop  10000   skip 100 -
        rtumbl 500.0   hfield 11.74   cut     2.3   tmax    3.0   -
        iwrite     6  C(t) modf 6 mfda 6 dsig 170.0 -
        orient         select type CA end

   Calculates the relaxation parameters for mainchain N-H spin pairs from the
   trajectory after alignment with the maincain CA in the comparison set.
   Anisotropic tumbling as a symmetric top is modeled with the long axis aligned
   along x (thus DTSX,y,z are 1,0,0) and a Dparallel/Deper ratio of 1.2.
   However, N-H vector angles to the long axis are trajectory averaged. 
   Both correlation functions for the internal motions as well as relaxation
   parameters are calculated.
