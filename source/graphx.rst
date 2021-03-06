.. py:module:: graphx

========
GRAPHICS
========

Graphics is a subparser of charmm, invoked by via the GRAPH command.
All of the miscellaneous commands (:doc:`miscom`), coordinate commands
(:doc:`corman`), and internal coordinate commands (:doc:`intcor`) are
available from the GRAPHX> prompt.  Only the 1st three characters are
used for primary graphics commands, but many of the options require
the 1st four characters.

The graphics facility has been extended to provide general X11
support, and the original Apollo GPR screen display has been dropped; 
a NODISPLAY version can also be built, which will generate all of the
derived files.  The other major enhancement is the production of
PostScript output files, in either color or grayscale; both X11 and
PostScript use the Apollo imaging model.  Additional information on X11
usage tips and compiling for X11 are given at the end of this document.
Finally, a recent addition is the production of input files for POV-Ray,
an excellent freeware ray tracing package for making high quality 
molecular images.  See   http://www.povray.org

Option keywords are indicated by the use of upper case; lower
case terms are variable values, generally real numbers, but decimal
points are not required.  Triplets ( x y z ) are position dependent;
omitted values are assumed to be zero.  Items enclosed in square 
brackets are [optional] but their absence often implies a default 
choice.  Default choices are indicated with an asterisk (*) in syntax 
listings where apropriate.


.. _graphx_summary:

GRAPHX COMMAND SUMMARY
----------------------

Initializing, exiting:

::

   GRAPhx  [XXSZ iwidth] [XYSZ iheight]  [ NVERtex integer ]
           [NOWIndow]

   OFF                 :: disable all graphics and exit (see END, below).
   END                 :: exit from command parser only; preserve window
   RESet               :: revert to default view, scale
   DEFault-colors      :: revert to default colors

   HELp                :: provides a command listing

Available within the graphics subcommand parser:

::

   miscellaneous-command-spec         ! see *note miscom:(chmdoc/miscom.doc).
   IC        ic-subcommand-spec       ! see *note intcor:(chmdoc/intcor.doc).
   COOR      coor-subcommand-spec     ! see *note corman:(chmdoc/corman.doc).
   TRAJ      traj-subcommand-spec     ! see *note dynamc:(chmdoc/dynamc.doc).

These commands affect what is viewed:

::

   DRAw [atom-selection]   :: change displayed atoms, force a redraw

   DISplay [ON]* [MAIN]* [COMP] [BOND]* [VECT] [ATOM] [TEXT] -
           [OFF]
                         [HBONds] [LABEls] [AXES]    :: change displayed objects

   COLor color-name [brightfactor] [COMP] atom-selection
                    [ LEVel izcue-level ]

     color-name =   NONE   RED  BLACk   YELLow GREEn  WHITe
                    BLUE   CYAN MAGEnta GRAY   ORANge BROWn
                    PURPle TURQuoise    CHARtreuse  DKBLue
                       (or an integer from 0 to 15)

     izcue-level = an integer from 0 to 11; 0 == brightest, 11 == dimmest

   HBStyle  [COLOR color-name]  [WIDTH iwidth]  [DASH idash]  :: HBOND style
   LBL  label-type  label-atoms  SIZE label-size COLOr label-color
      label-type  = INIT SEGId RESN RESId TYPE ATNUm CHEM CHARge WEIGht
                    USER user-label
      label-atoms = FIRSt and/or atom-selection
      label-size  = VSMAll SMALl MEDIum LARGe                   default: SMALL
      label-color = color-name ( see COLOR command above ... )  default: YELLOW
      user-label  = up to 8 characters

   LINe iwidth         :: bonds or vectors; pixels
   RADii [DEFaults] [PARam] atom-scale [bond-scale] atom-sel
   AXEs  [XLIM xmin xmax] [YLIM ymin ymax] [ZLIM zmin zmax]
         [DEFAult]  (all limits set to -25, +25 A)
   STEreo [ON/OFF] [dist] [angle]

   TEXt [text-body]    :: set title text, null entry clears
   FONT [ VSMAll | SMALl | MEDIum | LARGe ]  :: text font; default MEDIUM


These commands change the view only:

::

   SCAle  factor [MOL/LAB] [REP int]
   BOXsize size  [MOL/LAB]
   CENter [atom-selection]
   MAXwindow
   POInt  x y z
   ROTate rx ry rz  [MOL/LAB] [REP int]
   TRAnslate x y z  [MOL/LAB] [REP int]
   ZCLip [low] high
   ZCUe [ [low] high ]/[AUTO]

These commands do not affect the display:

::

   AUTo [ON*/OFF]      :: redraw after every command
   ERAse [ON*/OFF]     :: screen clear prior to next drawing
   EXEcute pathname    :: execute a program; no arguments are passed

Output Files: the UNIT must be OPENed first...

::

   PLUto UNIT n                                 :: PLUTO FDAT file
   MAKE  UNIT n                                 :: LIGHT .atm file

   PSC   UNIT n [COLOr] [BWREv] [PORTrait]      :: PostScript file
                [INIT]
           [TERM]

   POV   UNIT n  [ INIT ]                       :: POV-Ray file
                 [ TERM ]
                 [ INCLude n ]  [ UOBJ n ]  
                 [ OBJect  obj-spec ]

(see the description for full syntax details)


.. _graphx_description:

DETAILED COMMAND DESCRIPTION
----------------------------

*  GRAPHX

   ::

      GRAPhx  [XXSZ iwidth]   [XYSZ iheight]     [NVERtex integer]
              [NOWIndow]

   Invoked from the main CHARMM command parser; if already initialized
   (i.e. GRAPHX ... END) the previous graphic states are retained.

   The XXSZ and XYSZ options set the X11 window width and height for
   the duration of the graphics session; a window resize is NOT passed
   back to graphics code at this time.

   The NVERtex option increases the allocated storage for the vertex array
   used by the POV-Ray object code, especially for the MESH and RIBBON
   objects; the default is set to NATOM. but may need to be increased for
   scripts which create a lot of POV objects, especially RIBBON objects.

   The NOWIndow option suppresses the X11 window for batch mode or remote
   usage; ideal for making PostScript files, etc.

*  HELP

   ::

      HELp

   List the available commands and syntax.

*  OFF

   ::

      OFF

   Disable all graphics and exit the graphics subcommand parser; see
   also the END command.

*  END

   ::

      END

   Exit from command parser only; allows other CHARMM commands to be
   performed w/o losing the current display.  Re-invoking graphics
   does not re-initialize the graphics settings.  Useful for trajectory
   movies using MERGE DRAW, e.g.

   ::

      end
      open unit 12 read file name dyn.trj.10
      merge draw firstu 12 nunit 1
      graphics

   Overlays can be created (on screen only) by preceding the END command 
   in the above example with ERAse OFF.

*  RESET

   ::

      RESet

   Restores view settings to the program defaults ( scale, translation
   and rotations).

*  DISPLAY

   ::

      DISplay [ON*]  [MAIN*] [COMP] [VECT] [ATOM] [BOND*] [TEXT] [HBONds] -
              [OFF]
                     [LABEls] [AXES]

   Turns the display of various graphic features on or off; the default
   is DISPLAY ON MAIN BOND TEXT which will show the connectivity of
   the atoms in main coordinate set, using the default CHARMM title.
   The options are:

   ==========  =======================================================
   **ON**/OFF  enable (default) or disable the display of a feature
   **MAIN**    the main coordinate set; displayed by default
   COMP        the comparison coordinate set; both may be displayed
   **BOND**    atom connectivity as atom-colored half-bonds; default
   ATOM        filled circles using current radius (see RADII)
   **TEXT**    title display
   HBONDS      current HBOND list, using double width lines
   LABELS      residue names, atom types, user labels ...
   AXES        lab frame axes; + end solid, labeled; - end dashed
   ==========  =======================================================
   
   Examples:
   
   ::
   
      display atom on          ! enable atom display; MAIN assumed
      display text             ! toggles title display
      display hbonds off       ! disable H-bond display


*  DRAW

   ::
   
      DRAw [atom-selection]

   Forces a redraw when AUTO mode is off; also used to to change which
   set of atoms is currently being displayed.  All display modes and
   output files from PLOT, PSC, PLUTO, XMOLE, and MAKE commands use this
   atom selection.  The initial selection is all atoms and may be restored 
   via:

   ::
   
      draw sele all end

*  COLOR

   ::

      COLor color-name [brightfactor] [COMP] atom-selection
                       [LEVel i]

   Sets the color of individual atoms according to the atom
   selection, using one of the color names below:

   ::
   
      NONE   RED  BLACk   YELLow GREEn  WHITe
      BLUE   CYAN MAGEnta GRAY   ORANge BROWn
      PURPle TURQuoise    CHARtreuse  DKBLue

   The color applies to the main coordinate set unless COMP is
   specified; brightfactor is a relative intensity, 0.0 -- 1.0
   Alternatively, the zcue levels may be set directly to one of the
   indices, an integer in the range 0:11 with 0 being fully bright.

   Setting all atoms to WHITE, and using the brightfactor gives the
   best control over how PostScript grayscale output will look on the
   printed page; note that NONE is background, BLACK is the basic carbon
   color is really gray, and WHITE will be black in PostScript grayscale
   output or in color output with black/white reversal.

   Examples: All carbons for segment s are colored cyan:

   ::
   
      color cyan sele type c* .and. segid s end

   Colored based on weighting array:

   ::
   
      color green lev 11 sele prop wmain .gt. 1.0 end
      color green lev 10 sele prop wmain .gt. 2.0 end
                     :
                    etc.

*  DEFAULT

   ::
   
      DEFault-colors

   Restore the default color assignments, based on element type.

   Example:

   ::
   
      default

*  LINE

   ::

      LINe iwidth  (bonds or vectors; pixels)

   Set the line width for bonds & vectors, in pixels (integer).

   Example:

   ::
   
      line 2

*  RADII

   ::
   
      RADii [DEFaults] [PARam] scale [bond] atom-sel

   Sets the radius for displaying atoms, and for output files produced
   by the PLOT, PSC, XMOLE, and MAKE commands.  The options are:

      ===============  ======================================================
        scale          required if no other options are specified, and
                       assumed to be 1.0 if omitted; performs a relative
                       scale if used by itself, or scales the radii set
                       by the DEFAULT or PARAM options

        DEFAULTS       set radii to a convenient size for display, based
                       on atom type ( C, N, O, ... )

        PARAM          use VDW radii

        bond           value for bond radii for LIGHT program (see MAKE)

        atom-sel       atom selection to apply the radii command to
      ===============  ======================================================
   
   Examples:

   ::
   
      rad 0.8                       ! reduces radii to 80% of current value
      radii param .5                ! set radii to 50% of VDW
      radii param .5  .15           ! bond radii to 0.15 A
      rad 1.5 sele type H* end      ! enlarge all H atoms by 50%

*  LBL

   :: 
   
      LBL  label-type  label-atoms  SIZE label-size  COLOR label-color
      label-type  = INIT SEGID RESN* RESID* TYPE CHEM CHARGE WEIGHT USER user-label
      label-atoms = FIRST* and/or atom-selection
      label-size  = VSMALL SMALL* MEDIUM LARGE
      label-color = color-name ( see COLOR command above; default: YELLOW )
      user-label  = up to 8 characters

   The LBL command identifies which atoms are to be labeled, what
   atom attributes are to be included in the label, and the relative
   size of the labels; the defaults are marked with an asterisk [*].

   The INIT option clears all labels, and any other options are ignored.

   One or more the following attributes may be included in the label,
   by simply including the keyword(s) in the LBL statement:

      ======      ===========================================
      SEGId       segment name (from GENErate; A4)
      RESN        residue name (from the RTF; A4)
      RESId       residue ID, a numeral (A4)
      TYPE        atom type, e.g. N, CA, CB ...  (A4)
      CHEM        atom parameter type code (A4)
      CHARge      atomic charge (G12.4)
      WEIGht      value stored in the weight vector (G12.4)
      USER        arbitrary user-specified text (A8)
      ======      ===========================================
      
   The label length (24 bytes) is such that all attributes may NOT be
   displayed simulaneously; in particular, CHARge and WEIGht may not be
   displayed at the same time for the same atom.

   SIZE is specified by one the keywords VSMALL, SMALL, MEDIUM, LARGE,
   with a default of SMALL.  The COLOR keyword allows setting the label
   color, using the same color names as the COLor command; the default
   label color is yellow.  Each use of the LBL command can create a group
   of labels of a different size and color, for atoms which don't overlap
   with any previous label atom selections.

   The blank delimited word following the keyword USER, up to 8 characters, allows
   the use of any text string as a label for the selected atoms.

   The default is to label the first atom of each residue; an atom
   selection overrides this, unless the FIRSt keyword is present;
   in this case, the first atom of each selected residue is labeled.

   Examples:
   
   ::

      ! the first atom of each residue is labeled by name with normal text
      LBL RESN COLOR CYAN

      ! the first atom of each residue in the segment MAIN is labeled
      ! by name and number with very small text
      LBL RESN RESID FIRST SELE SEGID MAIN END SIZE VSMALL

      ! all oxygen atoms are labeled by charge with small text
      LBL CHARGE SELE TYPE O* END COLOR CHAR

      ! all alpha-carbons are labeled by the weight vector with medium text
      LBL WEIGHT SELE TYPE CA END SIZE MEDIUM

      ! enter a null label; may be used to selectively "blank" labels
      ! in this case, all alpha-carbon labels are set to a string of blanks
      ! for display efficiency, LBL INIT is preferable
      LBL SELE TYPE CA END USER

      ! show the location of formal charges on amino acid side chains
      LBL USER - SELE RESN ASP .AND. TYPE CG END SIZE LARGE COLOR GREEN
      LBL USER - SELE RESN GLU .AND. TYPE CD END SIZE LARGE COLOR GREEN
      LBL USER + SELE RESN ARG .AND. TYPE CZ END SIZE LARGE COLOR GREEN
      LBL USER + SELE RESN LYS .AND. TYPE NZ END SIZE LARGE COLOR GREEN

*  HBSTYLE

   ::
   
      HBStyle  [COLOR color-name]  [WIDTH iwidth]  [DASH idash] 

   Set the style for representing HBONDS; color-name is as for the
   COLOR command, and iwidth and idash are integers in pixel units.
   Specifying HBSTYLE alone resets to the default style, which is
   equivalent to

   ::
   
      HBSTYLE COLOR ORANGE WIDTH 4 DASH 4    (N.B. DASH 0 = solid line)

   If at least one option (COLOR, WIDTH, or DASH) is specified, the
   remaining options are unchanged; thus HBSTYLE COLOR WHITE will
   not reset the WIDTH to 4 pixels, but leave it as it was.

*  AXES

   ::
   
      AXEs  [XLIM xmin xmax] [YLIM ymin ymax] [ZLIM zmin zmax]
            [DEFAult]

   Changes the length of the displayed axes; the default, which is
   restored via the DEFAult keyword, is -25 to +25 A for each dimension.
   Only the axes lines specified are changed, e.g. AXES ZLIM -5 5 only
   changes the endpoints of the Z axes line.

*  STEREO

   ::
   
      STEreo [ON/OFF] [dist] [angle]

   Invoke side-by-side stereo mode for screen display and for output
   files produced by the PLOT and MAKE commands, when the ON keyword
   is used, or when stereo is off.  The  dist  option controls the
   separation between the two images; the  angle  option specifies
   the parallax angle for left and right eye views.  The OFF keyword
   is used to return to mono mode, and is assumed if the command is
   used while in stereo mode.  The default angle is 7 degrees, but the
   dist parameter defaults to an interval derived from the window width;
   wider windows are best for stereo.

   Example:

   ::
   
      stereo on 16.0 7 

*  TEXT

   ::
   
      TEXT [ text-body | "text-body" ]

   Supply the text for the title display; quotes override the conversion
   to all upper case, and a null entry clears the current title.  The
   quotes are removed, and are not displayed as part of the title; the
   current CHARMM title is the default graphics title.

*  FONT

   ::
   
      FONt [ VSMAll | SMALl | MEdium | LArge ]

   Change the title font to one of four sizes; the initial setting
   is MEDIUM, as is the default if no size is specified.

*  SCALE
   
   ::
   
      SCAle  factor [MOL/LAB] [REP int]

   Change the Angstrom/pixel scale factor; initially, 1 A = 32 pixels.
   The repeat factor can provide an ersatz zoom effect.  The default is
   the LAB frame (i.e. scale is done relative to screen center rather
   then the system center).


*  BOXSIZE

   ::

      BOXsize size  [MOL/LAB]

   Set the viewing 'box' to (size)x(size), in Angstroms.

*  CENTER

   ::

      CENter [atom-selection]

   Move the selected atoms to the center of display space.

*  MAXWINDOW

   ::
   
      MAXwindow

   Scales the molecule to fit in the display window.

*  POINT

   :: 
   
      POInt  x y z

   The point specified by  x y z  becomes the center of display space.

*  ROTATE

   ::
   
      ROTate rx ry rz  [MOL/LAB] [REP int]

   Apply a rotation to the viewing transform; does not affect the
   coordinates.  The default is the LAB frame (i.e. a Z rotation
   always spins the screen).  If the MOL (molecule) frame is used,
   then the rotation will be about the origin of the molecular coordinate
   system (i.e. does not depend on the current view matrix).
   
   (NOTE1: the ROTate command uses a left handed system, multiply rotation
   angles by -1 as necessary).
   
   (NOTE2: if more than one rotation angle is given as nonzero, then the
   rotations will occur sequentially, first x, then y, and then z).

   Examples:
   
   ::

      rot 0 90       ! rotate by 90 deg around the y axis
      rotate 180     ! rotate by 180 deg around the x axis
      rot 90 0 90    ! rotate by 90 deg around x, then 90 deg around z


*  TRANSLATE

   ::
   
      TRAnslate x y z  [MOL/LAB] [REP int]

   Apply a translation to the viewing transform; does not affect the
   coordinates.

   Examples:

   ::
   
      tran 2 9.5     ! translate +2 A along x axis, +9.5 along y axis
      tra 0 0 4      ! translate +4 axis along the z axis

*  ZCLIP

   ::
   
      ZCLip [low] high

   Set hither and yon clip limits; atoms outside the limits are not
   displayed, whether selected or not.

   Examples:

   ::
   
      zclip 10       ! atoms outside z = ( -10 .. +10 ) are not displayed
      zclip -5 10


*  ZCUE

   :: 
   
      ZCUe [[low] high ]/[AUTO]

   Controls the z coordinate range over which depth cueing will be
   applied.  AUTO use the displayed atom coordinates to set the limits.

   Examples:

   ::
   
      zcue 10        ! zcue from -10 to +10
      zcue -4 8
      zcue auto

*  AUTO

   :: 

      AUTo [ON/OFF]

   Enables or disables automatic redraw after every command; the
   initial setting is ON, and the command functions as a toggle.
   AUTO OFF is useful for making multiple changes without the
   time required for a redraw.

*  ERASE

   :: 

      ERAse [ON/OFF]

   Enables or disables a screen clear prior to the next drawing; the
   initial setting is ON, and the command functions as a toggle.
   ERASE OFF is useful for overlaying trajectory frames or other
   related collections of structures.  Possibly best used with AUTO OFF,
   and using DRAW for each structure.  Also applies to POV and PSC
   commands; after using the INIT keyword, the file is not closed until 
   one of the following:

   (1) PSC/POV command with TERMinate keyword; only correct method
   (2) the UNIT is CLOSEd; no title or page eject, etc.
   (3) the program terminates; no title or page eject, etc.

   Creates an overlay of 50 consecutive trajectory frames (coord traj 
   open on unit 20; .ps file open on unit 12):

   ::
   
      erase off
      psc unit 12 init
      set n 1
      label frmlp
      coor read file unit 20 ifile -1
      psc unit 12
      let n += 1
      if n .lt. 50 goto frmlp
      psc unit 12 term

*  EXECUTE

   ::
   
      EXEcute pathname

   Execute a program w/o arguments.

   Example:
   
   ::

      exe /bin/ls         ! list the current directory  (unix)

   Equivalent to the SYSTem command (:doc:`miscom`)


.. _graphx_output:

Output Facility
---------------

*  PSC

   ::

      PSC UNIT n [COLOr] [PORTrait] [BWREverse] [INIT | TERM]

   Writes out a PostScript display file using the current viewing transform
   and selected atoms; the default is grayscale in landscape mode, and may
   be changed using the COLOr and/or PORTrait keywords.  Colors are direct
   translations of the RGB color map used for Apollo and X11 displays, or
   arbitrary conversions to grayscale.  Direct control over grayscale may 
   be achieved by setting all atoms to COLOR WHITE, and using the 
   'brightfactor' option of the COLOR command to set the various gray 
   levels; note that full WHITE is black, since the background (paper) is 
   white in grayscale mode.  The default background is black in COLOR mode, but may changed to white
   with BWREverse keyword; atoms or text colored white on the screen will 
   be printed in black.  As always, the UNIT must be opened first 
   (recommended file extension .ps).

   The INIT keyword forces device initialization in ERASE OFF mode, and 
   should only be used in that context; it is not normally required. 
   Likewise, the TERM keyword writes the final part of the PostScript file in ERASE OFF mode,
   and does not draw the structure; only the graphics title and the final 
   few PostScript commands are written to the file.  See the ERAse command 
   for additional details and a usage example.

* PLUTO

   ::

      PLUto UNIT n

   Writes out atom coordinates and connectivity based on the current atom
   selection and view transform, in CSD FDAT format.  

   .. warning::
      999 atom limit!

   Stereo mode settings are ignored, as are radii and color.  Also, atoms
   should be renamed to their element types to get proper radii, etc within
   pluto (e.g. CA is not carbon alpha, it's calcium).  As with PLOT, the
   UNIT must be opened first ( recommended file extension  .fdat ).

   Example:

   ::

      rename atom C1 sele type C  end
      rename atom C2 sele type CA end
      rename atom N1 sele type N  end
      rename atom O1 sele type O  end
      rename atom H1 sele type H1 end
         ! etc., etc., etc.
      graphics
      scale .5
      rot 0 90
      open unit 50 write card name molecule.fdat
      pluto unit 50

*  MAKE

   ::
   
      MAKE  UNIT n

   Writes atom coordinates, radii, and color to a file formatted for
   the LIGHT program (available from the N.I.H.) which produces nice
   ray-traced images using current stereo settings and view transform.

   A 1280x1024 window is the optimum size for best scaling and centering
   of the image produced by LIGHT.

   The UNIT must be opened first; the LIGHT program requires that an
   extension of .atm be use for the file name.

*  POV

   ::
   
      POV  UNIT int  [ INIT ]  [ INCLude int ]
                     [ TERM ]

              Write atom-based objects; bonds, atoms, H-bonds, labels

        UNIT int    :: final POV file; user must open; also open int+1 for stereo
        INCLude int :: read file w. user override of default camera, lighting, ...

        INIT        :: intialize POV file in ERASE OFF mode (overlays)
        TERM        :: close POV file in ERASE OFF mode; other options ignored

   All units must be opened; in STEreo mode, UNIT numbers are used in
   pairs as N, N+1 for the left- and right-eye views.  The atom, bonds,
   H-bonds, and labels currently displayed are written to the file; just 
   as for the X11 display and PostScript output, this is determined by the 
   last atom selection used with a DRAw command, and graphics elements 
   enabled via the DISplay command.  The current viewpoint is used for 
   transforming coordinates for the POV output file, and the graphics 
   SCAle determines the camera distance.  The atom radii are determined 
   from the current graphics RADii, while the bond cylinders are scaled 
   based on the graphics LINe width.

   The POV file produced (recommended extension .pov) is a flat text file, 
   and may be manually edited to further customize the image that will be 
   ultimately created by processing the POV file.  Note that CHARMM only 
   creates an input file for POV-Ray; it is assumed that the user has 
   access to version 3.x of the POV-Ray program.  When using the povray3 
   program on Unix systems, the PPM bitmap output format is recommended.  
   For details on obtaining and using the POV-Ray program, see the web 
   site at  http://www.povray.org

   The INCLude file provides an easy mechanism to add user customization, 
   e.g. changing the light sources or camera position, or adding other 
   POV-Ray objects to be displayed.  Comments in the POV input file 
   created indicate the section which is completely replaced by the 
   userfile; using the CHARMM supplied default setup as a template for 
   INCLude files is strongly encouraged.

   The INIT keyword forces device initialization in ERASE OFF mode, and 
   should only be used in that context; it is not normally required.  
   Likewise, the TERM keyword writes the final part of the POV file in 
   ERASE OFF mode, and does not output a structure; only the graphics 
   title is written to the file.  See the ERAse command for additional 
   details and a PSC (PostScript) usage example.

   Example: write a single POV file 

   ::
   
      open unit 21 write card name protein.pov
      pov unit 21

   Example: write two POV files for a stereo image pair

   ::
   
      open unit 21 write card name protein-l.pov
      open unit 22 write card name protein-r.pov
      pov unit 21


.. _graphx_addendum:

X11 Usage Tips
--------------

On some workstations, the graphics display colors will only be
correct when that window has "focus".

The recommended "focus" policy is pointer-has-focus, i.e. the active
window is whichever one the mouse pointer is in; this permits typing
graphics commands in a mostly obscured terminal window (xterm, winterm,
hpterm, aixterm, etc.).  The alternative is to "lower" the text command
window via the window manager.  Automatic "topping" of the window with
focus is not recommended.

Although the graphics display window can be resized, the size change is
not currently passed back to the drawing routines.

X11 Compiling Problems
----------------------

The default assumptions are:

(1) required include files (Xlib.h, Xutil.h) reside in /usr/include/X11

(2) the library libX11.a resides in /usr/lib/X11

These files *must* exist for the X11 graphics to be compiled, and of
course the keyword XDISPLAY must be in the pref.dat file.  If the files
are NOT in the standard places (such as under HP-UX), the compile option
-I may be needed to point to the directory with the include files, and
the corresponding -L option may be needed to point to the directory 
which contains libX11.a for the final link step.

Other Useful Programs (all "freeware")
--------------------------------------

(1) ghostscript -- besides allowing on screen preview of PostScript files
prior to printing, this useful utility can convert the .ps files output by
the PSC command to other formats (e.g. GIF, PBM) and can also be used to make
an "encapsulated" PostScript (EPS) file.  Available via FTP from almost any
software archive that distributes GNU software (Free Software Foundation).

(2) xfig -- this X11-based technical drawing package can read EPS files, and
provides an ideal way to do "pasteup" of PostScript files from diverse sources
into a composite figure.  Available via FTP from sites carrying X11 software;
hardware specific sites may have versions pre-configured for SGI, HP, etc.

(3) pbmplus, netpbm -- this suite of bitmap manipulation programs offers the
best looking conversion to GIF files, using a PBM file from ghostscript as an
intermediate file type.  The names pbmplus and netpbm indicate minor variants
of the same suite, with netpbm being more recent, according to UseNet lore.
Available at the MIT X11 site, and similar FTP sites.

(4) light -- a ray-tracing program written by BR Brooks at NIH, which uses
the file format ouput by the graphics MAK command; used to produce movies
and slick-looking single images for presentations.  Contact the author.

(5) gnuplot -- 2D and 3D graphics, with X11 and PostScript support (and about
50 other output devices as well); useful for time series plots, etc., which
can be combined with molecular graphics via xfig (see above).  Available at
most GNU sites and from the "home" site at Dartmouth.

(6) ImageMagick -- a powerful suite of image manipulation programs, including
animation capabilities; available at ftp.x.org

(7) MPEG -- mpegencode and mpegplay produce and replay animations; multiple
images in PBM format can be combined into a movie

