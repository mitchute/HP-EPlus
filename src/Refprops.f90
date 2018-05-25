
MODULE FluidProperties_HPSim
        !RS Comment: Renamed FluidProperties_HPSim from FluidProperties as Energy+ already has a different refrigerant file
        !            called FluidProperties with different refrigerant data/modules. Renaming this module makes it easier for integration.

        ! MODULE INFORMATION:
        !       AUTHOR         Mike Turner
        !       DATE WRITTEN   10 December 99
        !       MODIFIED       Rick Strand (April 2000, May 2000)
        !                      Simon Rees  (May, June 2002)
		!                      Kenneth Tang (Nov, 2003)
        !       RE-ENGINEERED  Rick Strand (April 2000, May 2000)
		!		RE-ENGINEERED  Kennneth Tang (Nov,2003)
		!                      Ipseng Iu (Jan 2004)

        ! PURPOSE OF THIS MODULE:
        ! This module contains subroutines which determine and return properties
        ! of materials including enthalpy, quality, specific heat, and density.
        ! The module uses InputProcessor to read the material type and the
        ! associated charts from IN.IDF.  The module is only as powerful as the
        ! amount of data loaded into this file.

        ! METHODOLOGY EMPLOYED:
        ! The module will first check if the current refrigerant has been read
        ! in yet.  If not, it will get the data from IN.IDF and "store" it into
        ! a set of variables.  Any future iterations with that refrigerant will
        ! simply retrieve the data from storage instead of reading from the .IDF
        ! file again.  After the data is made available, the module uses input
        ! temperatures, pressures, and either quality or enthalpy to locate the
        ! state point and choose the proper routine.  Finally, it performs a
        ! double interpolation between temperatures and pressures or qualities
        ! which surround the point on a chart specified by the input conditions.
        ! The program is designed to work on either side of or under the vapor
        ! dome.  This data can be added as needed.
        !
        ! Where properties are invalid at particular pressure/temperature points
        ! in the data input file, zeros have to be inserted. This is necessary
        ! as the data structures are rectangular. The zero values are used to detect
        ! bounds of the data and issue appropriate warnings.
        !
        ! Properties of liquids (e.g. water) can be specified as glycol properties by
        ! supplying the same data for concentrations of 0.0 and 1.0 only.
        !
        ! Temperature data has to be supplied in ascending order only.

        ! REFERENCES:

        ! USE STATEMENTS

USE DataGlobals_HPSimIntegrated, ONLY: MaxNameLength

IMPLICIT NONE                           ! Enforce explicit typing of all variables
PRIVATE
        ! MODULE PARAMETER DEFINITIONS
CHARACTER(len=11) :: Refrig       = "REFRIGERANT"
CHARACTER(len=6)  :: Glycol       = "GLYCOL"
CHARACTER(len=8)  :: Pressure     = "PRESSURE"
CHARACTER(len=8)  :: Enthalpy     = "ENTHALPY"
CHARACTER(len=7)  :: Entropy      = "ENTROPY"
CHARACTER(len=7)  :: Density      = "DENSITY"
CHARACTER(len=12) :: SpecificHeat = "SPECIFICHEAT"
CHARACTER(len=12) :: Conductivity = "CONDUCTIVITY"
CHARACTER(len=9)  :: Viscosity    = "VISCOSITY"
CHARACTER(len=14) :: SurfaceTension  = "SURFACETENSION"
CHARACTER(len=5)  :: Fluid        = "FLUID"
CHARACTER(len=8)  :: GasFluid     = "FLUIDGAS"

        ! DERIVED TYPE DEFINITIONS
TYPE FluidPropsRefrigerantData
  CHARACTER(len=MaxNameLength) :: Name = ' '    ! Name of the refrigerant
  INTEGER :: NumPsPoints       = 0              ! Number of saturation pressure
  REAL, POINTER, DIMENSION(:)  :: PsTemps       ! Temperatures for saturation pressures
  REAL, POINTER, DIMENSION(:)  :: PsfValues     ! Saturation liquid pressures at PsTemps
  REAL, POINTER, DIMENSION(:)  :: PsgValues     ! Saturation vapor pressures at PsTemps
  INTEGER :: NumHPoints        = 0              ! Number of enthalpy points
  REAL, POINTER, DIMENSION(:)  :: HTemps        ! Temperatures for enthalpy points
  REAL, POINTER, DIMENSION(:)  :: HfValues      ! Enthalpy of saturated fluid at HfTemps
  REAL, POINTER, DIMENSION(:)  :: HfgValues     ! Enthalpy of saturated fluid/gas at HfgTemps
  INTEGER :: NumCpPoints       = 0              ! Number of specific heat of fluid points
  REAL, POINTER, DIMENSION(:)  :: CpTemps       ! Temperatures for specific heat points
  REAL, POINTER, DIMENSION(:)  :: CpfValues     ! Specific heat of saturated fluid at CpfTemps
  REAL, POINTER, DIMENSION(:)  :: CpfgValues    ! Specific heat of saturated fluid/gas at CpfgTemps
  INTEGER :: NumRhoPoints      = 0              ! Number of density of fluid points
  REAL, POINTER, DIMENSION(:)  :: RhoTemps      ! Temperatures for density of fluid points
  REAL, POINTER, DIMENSION(:)  :: RhofValues    ! Density of saturated fluid at RhofTemps
  REAL, POINTER, DIMENSION(:)  :: RhofgValues   ! Density of saturated fluid/gas at RhofgTemps
  INTEGER :: NumSPoints      = 0                ! Number of entropy of fluid points
  REAL, POINTER, DIMENSION(:)  :: STemps        ! Temperatures for entropy of fluid points
  REAL, POINTER, DIMENSION(:)  :: SfValues      ! Entropy of saturated fluid at STemps
  REAL, POINTER, DIMENSION(:)  :: SfgValues     ! Entropy of saturated fluid/gas at STemps
  INTEGER :: NumCPoints      = 0                ! Number of entropy of fluid points
  REAL, POINTER, DIMENSION(:)  :: CTemps        ! Temperatures for conductivity of fluid points
  REAL, POINTER, DIMENSION(:)  :: CfValues      ! Conductivity of saturated fluid at CTemps
  REAL, POINTER, DIMENSION(:)  :: CfgValues     ! Conductivity of saturated fluid/gas at CTemps
  INTEGER :: NumDVPoints      = 0               ! Number of viscosity of fluid points
  REAL, POINTER, DIMENSION(:)  :: DVTemps       ! Temperatures for dynamic viscosity of fluid points
  REAL, POINTER, DIMENSION(:)  :: DVfValues     ! Dynamic viscosity of saturated fluid at VTemps
  REAL, POINTER, DIMENSION(:)  :: DVfgValues    ! Dynamic viscosity of saturated fluid/gas at VTemps
  INTEGER :: NumSTPoints      = 0               ! Number of viscosity of fluid points
  REAL, POINTER, DIMENSION(:)  :: STTemps       ! Temperatures for surface tension of fluid points
  REAL, POINTER, DIMENSION(:)  :: STfValues     ! Surface tension of saturated fluid at VTemps
  REAL, POINTER, DIMENSION(:)  :: STgValues     ! Surface tension of saturated vapor at VTemps
  INTEGER :: NumSuperTempPts   = 0              ! Number of temperature points for superheated enthalpy
  INTEGER :: NumSuperPressPts  = 0              ! Number of pressure points for superheated enthalpy
  REAL, POINTER, DIMENSION(:)   :: SHTemps      ! Temperatures for superheated gas
  REAL, POINTER, DIMENSION(:)   :: SHPress      ! Pressures for superheated gas
  REAL, POINTER, DIMENSION(:,:) :: HshValues    ! Enthalpy of superheated gas at HshTemps, HshPress
  REAL, POINTER, DIMENSION(:,:) :: RhoshValues  ! Density of superheated gas at RhoshTemps, RhoshPress
  REAL, POINTER, DIMENSION(:,:) :: SshValues    ! Entropy of superheated gas at SshTemps, SshPress
  REAL, POINTER, DIMENSION(:,:) :: DVshValues   ! Dynamic viscosity of superheated gas at DVshTemps, DVshPress
  REAL, POINTER, DIMENSION(:,:) :: CshValues    ! Conductivity of superheated gas at CshTemps, CshPress
  REAL, POINTER, DIMENSION(:,:) :: CPshValues   ! Specific Heat of superheated gas at CPshTemps, CPshPress
  
  INTEGER :: NumSubcoolTempPts   = 0            ! Number of temperature points for subcooled enthalpy
  INTEGER :: NumSubcoolPressPts  = 0            ! Number of pressure points for subcooled enthalpy
  REAL, POINTER, DIMENSION(:)   :: SCTemps      ! Temperatures for subcooled gas
  REAL, POINTER, DIMENSION(:)   :: SCPress      ! Pressures for subcooled gas
  REAL, POINTER, DIMENSION(:,:) :: HscValues    ! Enthalpy of subcooled gas at HscTemps, HscPress
  REAL, POINTER, DIMENSION(:,:) :: RhoscValues  ! Density of subcooled gas at RhoscTemps, RhoscPress
  REAL, POINTER, DIMENSION(:,:) :: SscValues	! Entropy of subcooled gas at SscTemps, SscPress
  REAL, POINTER, DIMENSION(:,:) :: DVscValues   ! Dynamic viscosity of subcooled gas at DVscTemps, DVscPress
  REAL, POINTER, DIMENSION(:,:) :: CscValues    ! Conductivity of subcooled gas at CscTemps, CscPress
  REAL, POINTER, DIMENSION(:,:) :: CPscValues   ! Specific Heat of subcooled gas at CPscTemps, CPscPress
END TYPE

TYPE (FluidPropsRefrigerantData), ALLOCATABLE, DIMENSION(:) :: RefrigData

        ! INTERFACE BLOCK SPECIFICATIONS
        ! na

        ! MODULE VARIABLE DECLARATIONS
        ! na

LOGICAL :: GetInput = .TRUE.     ! Used to get the input once only
INTEGER :: NumOfRefrigerants = 0 ! Total number of refrigerants input by user
INTEGER :: NumOfGlycols = 0      ! Total number of glycols input by user

        ! ACCESSIBLE SPECIFICATIONS OF MODULE SUBROUTINES OR FUNCTONS:
PRIVATE GetFluidPropertiesData
PUBLIC	Pcrit
PUBLIC  Tcrit
PUBLIC  MW
PUBLIC  PQ
PUBLIC  TQ
PUBLIC  TP
PUBLIC  PH
PUBLIC  PS
PUBLIC  CheckFluidPropertyName
PUBLIC  FindRefrigerant
PUBLIC  FindArrayIndex
PUBLIC  FindArrayIndexSC
PRIVATE GetInterpolatedSatProp

CONTAINS

          ! MODULE SUBROUTINES:

SUBROUTINE GetFluidPropertiesData

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Rick Strand
          !       DATE WRITTEN   April 2000
          !       MODIFIED       May 2002 Simon Rees
          !                         Added saturated pressure data retreaval
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! The purpose of this subroutine is to read in all of the fluid
          ! property data contained in the user input file.

          ! METHODOLOGY EMPLOYED:
          ! Standard EnergyPlus methodology.  Derived type portions are
          ! allocated as necessary as the data is read into the program.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
          ! na

          ! SUBROUTINE PARAMETER DEFINITIONS:
  REAL, PARAMETER :: TempToler = 0.1    ! Some reasonable value for comparisons
  REAL, PARAMETER :: PressToler = 1.0   ! Some reasonable value for comparisons

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
  TYPE FluidTempData
    CHARACTER(len=MaxNameLength) Name       ! Name of the temperature list
    INTEGER NumOfTemps                      ! Number of temperatures in a particular arry
    REAL, POINTER, DIMENSION(:) :: Temps ! Temperature values (degrees C)
  END TYPE

  TYPE(FluidTempData), ALLOCATABLE, DIMENSION(:) :: FluidTemps

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  CHARACTER(len=MaxNameLength),DIMENSION(20) :: Alphas ! Reads string value from input file
  INTEGER :: Loop                    ! DO loop counter (various uses)
  INTEGER :: NumAlphas               ! States which alpha value to read from a "Number" line
  REAL, DIMENSION(500) :: Numbers    ! brings in data from IP
  INTEGER :: NumNumbers              ! States which number value to read from a "Numbers" line
  INTEGER :: Status                  ! Either 1 "object found" or -1 "not found"

  INTEGER :: InData
  INTEGER :: TempLoop
  INTEGER :: NumOfFluidTempArrays
  INTEGER :: NumOfSatFluidPropArrays
  INTEGER :: NumOfSHFluidPropArrays
  INTEGER :: NumOfSCFluidPropArrays
  INTEGER :: NumOfGlyFluidPropArrays
  CHARACTER(len=MaxNameLength) :: TempsName
  LOGICAL :: FirstSHMatch
  INTEGER :: NumOfPressPts
  LOGICAL :: ErrorsFound=.false.
  
  INTEGER, PARAMETER :: r64=KIND(1.0D0)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12) 
  REAL(r64), DIMENSION(500) :: TmpNumbers !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

          ! FLOW:
  NumOfFluidTempArrays    = GetNumObjectsFound('FluidPropertyTemperatures')
  NumOfSatFluidPropArrays = GetNumObjectsFound('FluidPropertySaturated')
  NumOfSHFluidPropArrays  = GetNumObjectsFound('FluidPropertySuperheated')
  NumOfSCFluidPropArrays  = GetNumObjectsFound('FluidPropertySubcooled')
  NumOfGlyFluidPropArrays = GetNumObjectsFound('FluidPropertyConcentration')

  CALL GetObjectItem('FLUIDNAMES',1,Alphas,NumAlphas,TmpNumbers,NumNumbers,Status)
    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
  
          ! Get a count on the number of refrigerants and the number of glycols entered
          ! so that the main derived types can be allocated
  DO Loop = 2, NumAlphas, 2
    IF (SameString(Alphas(Loop),Refrig)) THEN
      NumOfRefrigerants = NumOfRefrigerants + 1
    ELSEIF (SameString(Alphas(Loop),Glycol)) THEN
      NumOfGlycols = NumOfGlycols + 1
    ELSE
      CALL ShowSevereError('GetFluidPropertiesData: Only REFRIGERANT or GLYCOL allowed as fluid types in the FluidNames syntax')
      CALL ShowContinueError('Fluid Type in Error='//TRIM(Alphas(Loop))//', Fluid Name='//TRIM(Alphas(Loop-1)))
      ErrorsFound=.true.
    END IF
  END DO

  IF (NumOfRefrigerants > 0) THEN
      ALLOCATE(RefrigData(NumOfRefrigerants))
  END IF

          ! Take the fluid names and assign them to the appropriate derived type
  NumOfRefrigerants = 0
  DO Loop = 2, NumAlphas, 2
    IF (SameString(Alphas(Loop),Refrig)) THEN
      NumOfRefrigerants = NumOfRefrigerants + 1
      RefrigData(NumOfRefrigerants)%Name = Alphas(Loop-1)
    END IF
  END DO

          ! Read in all of the temperature arrays in the input file
  ALLOCATE(FluidTemps(NumOfFluidTempArrays))

  DO Loop = 1, NumOfFluidTempArrays
    
    CALL GetObjectItem('FluidPropertyTemperatures',Loop,Alphas,NumAlphas, &
                        TmpNumbers,NumNumbers,Status)
    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    
    FluidTemps(Loop)%Name       = Alphas(1)
    FluidTemps(Loop)%NumOfTemps = NumNumbers

    ALLOCATE(FluidTemps(Loop)%Temps(FluidTemps(Loop)%NumOfTemps))
    FluidTemps(Loop)%Temps = Numbers(1:NumNumbers)

  END DO

          ! *************** REFRIGERANTS ***************
          ! Go through each refrigerant found in the fluid names statement and read in the data
          ! Note that every valid fluid must have ALL of the necessary data or a fatal error will
          ! be produced.
  DO Loop = 1, NumOfRefrigerants

          ! For each property, cycle through all the valid input until the proper match is found.

          ! **********    SATURATED DATA SECTION    **********

          ! Get: ***** Saturation Pressure temperatures and data (fluidgas only) *****
          ! This section added by S.J.Rees May 2002.

       TempsName     = ' '
    DO InData = 1, NumOfSatFluidPropArrays

      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    
      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),Pressure) )              .AND. &
           (SameString(Alphas(3),GasFluid ) )           ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
            TempsName = FluidTemps(TempLoop)%Name
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
            RefrigData(Loop)%NumPsPoints = FluidTemps(TempLoop)%NumOfTemps
            ALLOCATE(RefrigData(Loop)%PsTemps(RefrigData(Loop)%NumPsPoints))
			ALLOCATE(RefrigData(Loop)%PsgValues(RefrigData(Loop)%NumPsPoints))

            ! Same number of points so assign the values
            RefrigData(Loop)%PsTemps  = FluidTemps(TempLoop)%Temps
			RefrigData(Loop)%PsgValues = Numbers(1:NumNumbers)
            EXIT ! the TempLoop DO loop
          END IF
        END DO  ! ...end of FluidTemps DO loop
        EXIT ! the InData DO loop

      END IF

          ! If it made it all the way to the last input occurrence and didn't find a match, 
          ! then no sat press data found
      IF (InData == NumOfSatFluidPropArrays) THEN
        CALL ShowSevereError('GetFluidPropertiesData: No fluid saturation pressure found')
        CALL ShowContinueError('Was looking for properties for Refrigerant='//TRIM(RefrigData(Loop)%Name))
        ErrorsFound=.true.
      ENDIF

    END DO  ! ...end of DO loop through all of the input syntax trying to find saturation pressure for this refrigerant

          ! Get: ***** Saturation Pressure temperatures and data (fluid only) *****
          ! This section added by I.S. Iu Nov 2006.

       TempsName     = ' '
    DO InData = 1, NumOfSatFluidPropArrays

      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status) 
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),Pressure) )              .AND. &
           (SameString(Alphas(3),Fluid ) )           ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
            TempsName = FluidTemps(TempLoop)%Name
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
			ALLOCATE(RefrigData(Loop)%PsfValues(RefrigData(Loop)%NumPsPoints))

            ! Same number of points so assign the values
			RefrigData(Loop)%PsfValues = Numbers(1:NumNumbers)
            EXIT ! the TempLoop DO loop
          END IF
        END DO  ! ...end of FluidTemps DO loop
        EXIT ! the InData DO loop

      END IF

          ! If it made it all the way to the last input occurrence and didn't find a match, 
          ! then no sat press data found
      IF (InData == NumOfSatFluidPropArrays) THEN
        CALL ShowSevereError('GetFluidPropertiesData: No fluid saturation pressure found')
        CALL ShowContinueError('Was looking for properties for Refrigerant='//TRIM(RefrigData(Loop)%Name))
        ErrorsFound=.true.
      ENDIF

    END DO  ! ...end of DO loop through all of the input syntax trying to find saturation pressure for this refrigerant

          ! Get: ***** ENTHALPY of SATURATED LIQUID *****
    TempsName     = " "
    DO InData = 1, NumOfSatFluidPropArrays
      
      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),Enthalpy) )              .AND. &
           (SameString(Alphas(3),Fluid ) )                ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
            TempsName = FluidTemps(TempLoop)%Name
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
            RefrigData(Loop)%NumHPoints = FluidTemps(TempLoop)%NumOfTemps
            ALLOCATE(RefrigData(Loop)%HTemps(RefrigData(Loop)%NumHPoints))
            ALLOCATE(RefrigData(Loop)%HfValues(RefrigData(Loop)%NumHPoints))

            ! Same number of points so assign the values
            RefrigData(Loop)%HTemps  = FluidTemps(TempLoop)%Temps
            RefrigData(Loop)%HfValues = Numbers(1:NumNumbers)

            EXIT ! the TempLoop DO loop
          END IF
        END DO  ! ...end of FluidTemps DO loop

        EXIT ! the InData DO loop

      END IF

          ! If it made it all the way to the last input occurrence and didn't find a match, then no sat fluid enthalpy data found
    END DO  ! ...end of DO loop through all of the input syntax trying to find saturated fluid enthalpy for this refrigerant

          ! Get: ***** ENTHALPY of SATURATED LIQUID/VAPOR ***** (difference between Hf and Hg, i.e. Hfg)
    DO InData = 1, NumOfSatFluidPropArrays
      
      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),Enthalpy) )              .AND. &
           (SameString(Alphas(3),GasFluid ) )             ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
           
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
            ALLOCATE(RefrigData(Loop)%HfgValues(RefrigData(Loop)%NumHPoints))

            ! Same number of points so assign the values
            RefrigData(Loop)%HfgValues = Numbers(1:NumNumbers)
            EXIT ! the TempLoop DO loop
          END IF

          ! If it made it all the way to the last temperature array and didn't find a match, then no match was found     
        END DO  ! ...end of FluidTemps DO loop
        EXIT ! the InData DO loop
      END IF
    END DO  ! ...end of DO loop through all of the input syntax trying to find saturated gas/fluid enthalpy for this refrigerant

          ! Get: ***** SPECIFIC HEAT of SATURATED LIQUID *****
    TempsName     = " "
    DO InData = 1, NumOfSatFluidPropArrays

      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),SpecificHeat) )          .AND. &
           (SameString(Alphas(3),Fluid ) )                ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
            TempsName = FluidTemps(TempLoop)%Name
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
            RefrigData(Loop)%NumCpPoints = FluidTemps(TempLoop)%NumOfTemps
            ALLOCATE(RefrigData(Loop)%CpTemps(RefrigData(Loop)%NumCpPoints))
            ALLOCATE(RefrigData(Loop)%CpfValues(RefrigData(Loop)%NumCpPoints))

            ! Same number of points so assign the values
            RefrigData(Loop)%CpTemps  = FluidTemps(TempLoop)%Temps
            RefrigData(Loop)%CpfValues = Numbers(1:NumNumbers)

            EXIT ! the TempLoop DO loop

          END IF
        END DO  ! ...end of FluidTemps DO loop
        EXIT ! the InData DO loop
      END IF
      
    END DO  ! ...end of DO loop through all of the input syntax trying to find saturated fluid Cp for this refrigerant

          ! Get: ***** SPECIFIC HEAT of SATURATED LIQUID/VAPOR ***** (difference between Cpf and Cpg, i.e. Cpfg)
    DO InData = 1, NumOfSatFluidPropArrays

      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),SpecificHeat) )          .AND. &
           (SameString(Alphas(3),GasFluid ) )             ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
            
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
            ALLOCATE(RefrigData(Loop)%CpfgValues(RefrigData(Loop)%NumCpPoints))

            ! Make sure the number of points in the two arrays (temps and values) are the same

            ! Same number of points so assign the values
            RefrigData(Loop)%CpfgValues = Numbers(1:NumNumbers)

            EXIT ! the TempLoop DO loop

          END IF

        END DO  ! ...end of FluidTemps DO loop

        EXIT ! the InData DO loop

      END IF

    END DO  ! ...end of DO loop through all of the input syntax trying to find saturated gas/fluid Cp for this refrigerant

          ! Get: ***** DENSITY of SATURATED LIQUID *****
    TempsName     = " "
    DO InData = 1, NumOfSatFluidPropArrays

      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),Density) )               .AND. &
           (SameString(Alphas(3),Fluid ) )                ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
            TempsName = FluidTemps(TempLoop)%Name
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
            RefrigData(Loop)%NumRhoPoints = FluidTemps(TempLoop)%NumOfTemps
            ALLOCATE(RefrigData(Loop)%RhoTemps(RefrigData(Loop)%NumRhoPoints))
            ALLOCATE(RefrigData(Loop)%RhofValues(RefrigData(Loop)%NumRhoPoints))

            ! Same number of points so assign the values
            RefrigData(Loop)%RhoTemps  = FluidTemps(TempLoop)%Temps
            RefrigData(Loop)%RhofValues = Numbers(1:NumNumbers)

            EXIT ! the TempLoop DO loop

          END IF

        END DO  ! ...end of FluidTemps DO loop

        EXIT ! the InData DO loop

      END IF

          ! If it made it all the way to the last input occurrence and didn't find a match, then no sat fluid density data found

    END DO  ! ...end of DO loop through all of the input syntax trying to find saturated fluid enthalpy for this refrigerant

          ! Get: ***** DENSITY of SATURATED LIQUID/VAPOR ***** (difference between Rhof and Rhog, i.e. Rhofg)
    DO InData = 1, NumOfSatFluidPropArrays

      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),Density) )               .AND. &
           (SameString(Alphas(3),GasFluid ) )             ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
            
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
            ALLOCATE(RefrigData(Loop)%RhofgValues(RefrigData(Loop)%NumRhoPoints))

            ! Make sure the number of points in the two arrays (temps and values) are the same

            ! Same number of points so assign the values
            RefrigData(Loop)%RhofgValues = Numbers(1:NumNumbers)

            EXIT ! the TempLoop DO loop

          END IF

          ! If it made it all the way to the last temperature array and didn't find a match, then no match was found

        END DO  ! ...end of FluidTemps DO loop

        EXIT ! the InData DO loop

      END IF

          ! If it made it all the way to the last input occurrence and didn't find a match, then no sat f/g density data found
      
    END DO  ! ...end of DO loop through all of the input syntax trying to find saturated gas/fluid density for this refrigerant

          ! Get: ***** ENTROPY of SATURATED LIQUID *****
    TempsName     = " "
    DO InData = 1, NumOfSatFluidPropArrays

      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),Entropy) )               .AND. &
           (SameString(Alphas(3),Fluid ) )                ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
            TempsName = FluidTemps(TempLoop)%Name
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
            RefrigData(Loop)%NumSPoints = FluidTemps(TempLoop)%NumOfTemps
            ALLOCATE(RefrigData(Loop)%STemps(RefrigData(Loop)%NumSPoints))
            ALLOCATE(RefrigData(Loop)%SfValues(RefrigData(Loop)%NumSPoints))

            ! Make sure the number of points in the two arrays (temps and values) are the same

            ! Same number of points so assign the values
            RefrigData(Loop)%STemps  = FluidTemps(TempLoop)%Temps
            RefrigData(Loop)%SfValues = Numbers(1:NumNumbers)

            EXIT ! the TempLoop DO loop

          END IF

        END DO  ! ...end of FluidTemps DO loop

        EXIT ! the InData DO loop

      END IF

          ! If it made it all the way to the last input occurrence and didn't find a match, then no sat fluid density data found
      IF (InData == NumOfSatFluidPropArrays) THEN
        CALL ShowSevereError('GetFluidPropertiesData: No saturated fluid entropy found')
        CALL ShowContinueError('Was looking for properties for Refrigerant='//TRIM(RefrigData(Loop)%Name))
        ErrorsFound=.true.
      ENDIF

    END DO  ! ...end of DO loop through all of the input syntax trying to find saturated fluid enthalpy for this refrigerant

          ! Get: ***** ENTROPY of SATURATED LIQUID/VAPOR ***** (difference between Rhof and Rhog, i.e. Rhofg)
    DO InData = 1, NumOfSatFluidPropArrays

      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),Entropy) )               .AND. &
           (SameString(Alphas(3),GasFluid ) )             ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
            
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
            ALLOCATE(RefrigData(Loop)%SfgValues(RefrigData(Loop)%NumSPoints))

            ! Same number of points so assign the values
            RefrigData(Loop)%SfgValues = Numbers(1:NumNumbers)

            EXIT ! the TempLoop DO loop

          END IF

          ! If it made it all the way to the last temperature array and didn't find a match, then no match was found
          
        END DO  ! ...end of FluidTemps DO loop

        EXIT ! the InData DO loop

      END IF

          ! If it made it all the way to the last input occurrence and didn't find a match, then no sat f/g density data found
      IF (InData == NumOfSatFluidPropArrays) THEN
        CALL ShowSevereError('GetFluidPropertiesData: No saturated gas/fluid entropy found')
        CALL ShowContinueError('Was looking for properties for Refrigerant='//TRIM(RefrigData(Loop)%Name))
        ErrorsFound=.true.
      ENDIF

    END DO  ! ...end of DO loop through all of the input syntax trying to find saturated gas/fluid density for this refrigerant

         ! Get: ***** DYNAMIC VISCOSITY of SATURATED LIQUID *****
    TempsName     = " "
    DO InData = 1, NumOfSatFluidPropArrays

      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),Viscosity) )               .AND. &
           (SameString(Alphas(3),Fluid ) )                ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
            TempsName = FluidTemps(TempLoop)%Name
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
            RefrigData(Loop)%NumDVPoints = FluidTemps(TempLoop)%NumOfTemps
            ALLOCATE(RefrigData(Loop)%DVTemps(RefrigData(Loop)%NumDVPoints))
            ALLOCATE(RefrigData(Loop)%DVfValues(RefrigData(Loop)%NumDVPoints))

            ! Same number of points so assign the values
            RefrigData(Loop)%DVTemps  = FluidTemps(TempLoop)%Temps
            RefrigData(Loop)%DVfValues = Numbers(1:NumNumbers)

            EXIT ! the TempLoop DO loop

          END IF

        END DO  ! ...end of FluidTemps DO loop

        EXIT ! the InData DO loop

      END IF

          ! If it made it all the way to the last input occurrence and didn't find a match, then no sat fluid density data found

    END DO  ! ...end of DO loop through all of the input syntax trying to find saturated fluid enthalpy for this refrigerant

          ! Get: ***** VISCOSITY of SATURATED LIQUID/VAPOR ***** (difference between DVf and DVg, i.e. DVfg)
    DO InData = 1, NumOfSatFluidPropArrays

      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status) 
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),Viscosity) )               .AND. &
           (SameString(Alphas(3),GasFluid ) )             ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
           
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
            ALLOCATE(RefrigData(Loop)%DVfgValues(RefrigData(Loop)%NumDVPoints))

            ! Make sure the number of points in the two arrays (temps and values) are the same

            ! Same number of points so assign the values
            RefrigData(Loop)%DVfgValues = Numbers(1:NumNumbers)

            EXIT ! the TempLoop DO loop

          END IF

          ! If it made it all the way to the last temperature array and didn't find a match, then no match was found
       
        END DO  ! ...end of FluidTemps DO loop

        EXIT ! the InData DO loop

      END IF

    END DO  ! ...end of DO loop through all of the input syntax trying to find saturated gas/fluid density for this refrigerant

        ! Get: ***** CONDUCTIVITY of SATURATED LIQUID *****
    TempsName     = " "
    DO InData = 1, NumOfSatFluidPropArrays

      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),Conductivity) )               .AND. &
           (SameString(Alphas(3),Fluid ) )                ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
            TempsName = FluidTemps(TempLoop)%Name
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
            RefrigData(Loop)%NumCPoints = FluidTemps(TempLoop)%NumOfTemps
            ALLOCATE(RefrigData(Loop)%CTemps(RefrigData(Loop)%NumCPoints))
            ALLOCATE(RefrigData(Loop)%CfValues(RefrigData(Loop)%NumCPoints))

            ! Same number of points so assign the values
            RefrigData(Loop)%CTemps  = FluidTemps(TempLoop)%Temps
            RefrigData(Loop)%CfValues = Numbers(1:NumNumbers)

            EXIT ! the TempLoop DO loop

          END IF

          ! If it made it all the way to the last temperature array and didn't find a match, then no match was found
         
        END DO  ! ...end of FluidTemps DO loop

        EXIT ! the InData DO loop

      END IF

          ! If it made it all the way to the last input occurrence and didn't find a match, then no sat fluid density data found
     
    END DO  ! ...end of DO loop through all of the input syntax trying to find saturated fluid enthalpy for this refrigerant

          ! Get: ***** CONDUCTIVITY of SATURATED LIQUID/VAPOR ***** (difference between Cf and Cg, i.e. Cfg)
    DO InData = 1, NumOfSatFluidPropArrays

      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),Conductivity) )               .AND. &
           (SameString(Alphas(3),GasFluid ) )             ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
          
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
            ALLOCATE(RefrigData(Loop)%CfgValues(RefrigData(Loop)%NumCPoints))

            ! Same number of points so assign the values
            RefrigData(Loop)%CfgValues = Numbers(1:NumNumbers)

            EXIT ! the TempLoop DO loop

          END IF

          ! If it made it all the way to the last temperature array and didn't find a match, then no match was found
        
        END DO  ! ...end of FluidTemps DO loop

        EXIT ! the InData DO loop

      END IF

          ! If it made it all the way to the last input occurrence and didn't find a match, then no sat f/g density data found
    
    END DO  ! ...end of DO loop through all of the input syntax trying to find saturated gas/fluid density for this refrigerant

    
	      ! Get: ***** SURFACE TENSION of SATURATED LIQUID *****
    TempsName     = " "
    DO InData = 1, NumOfSatFluidPropArrays

      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),SurfaceTension) )               .AND. &
           (SameString(Alphas(3),Fluid ) )                ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
            TempsName = FluidTemps(TempLoop)%Name
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.
            
            !ISI - 11/07/06
			RefrigData(Loop)%NumSTPoints = FluidTemps(TempLoop)%NumOfTemps
            ALLOCATE(RefrigData(Loop)%STTemps(RefrigData(Loop)%NumSTPoints))
            ALLOCATE(RefrigData(Loop)%STfValues(RefrigData(Loop)%NumSTPoints))
            
            ! Same number of points so assign the values
			RefrigData(Loop)%STTemps  = FluidTemps(TempLoop)%Temps !ISI - 11/07/06
            RefrigData(Loop)%STfValues = Numbers(1:NumNumbers)

            EXIT ! the TempLoop DO loop

          END IF

          ! If it made it all the way to the last temperature array and didn't find a match, then no match was found
         
        END DO  ! ...end of FluidTemps DO loop

        EXIT ! the InData DO loop

      END IF

          ! If it made it all the way to the last input occurrence and didn't find a match, then no sat fluid density data found
     
    END DO  ! ...end of DO loop through all of the input syntax trying to find saturated fluid enthalpy for this refrigerant


	      ! Get: ***** SURFACE TENSION of SATURATED VAPOR *****
    TempsName     = " "
    DO InData = 1, NumOfSatFluidPropArrays

      CALL GetObjectItem('FluidPropertySaturated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IF ( (SameString(Alphas(1),RefrigData(Loop)%Name) ) .AND. &
           (SameString(Alphas(2),SurfaceTension) )               .AND. &
           (SameString(Alphas(3),GasFluid ) )                ) THEN

        DO TempLoop = 1, NumOfFluidTempArrays

          IF (SameString(Alphas(4),FluidTemps(TempLoop)%Name)) THEN
            TempsName = FluidTemps(TempLoop)%Name
            ! At this point, we have found the correct input line and found a match
            ! for the temperature array.  It's time to load up the local derived type.

            !ISI - 11/07/06
			RefrigData(Loop)%NumSTPoints = FluidTemps(TempLoop)%NumOfTemps
            ALLOCATE(RefrigData(Loop)%STgValues(RefrigData(Loop)%NumSTPoints))
            
            ! Same number of points so assign the values
            RefrigData(Loop)%STgValues = Numbers(1:NumNumbers)

            EXIT ! the TempLoop DO loop

          END IF

          ! If it made it all the way to the last temperature array and didn't find a match, then no match was found
         
        END DO  ! ...end of FluidTemps DO loop

        EXIT ! the InData DO loop

      END IF

          ! If it made it all the way to the last input occurrence and didn't find a match, then no sat fluid density data found
     
    END DO  ! ...end of DO loop through all of the input syntax trying to find saturated fluid enthalpy for this refrigerant

!*************************   SUPERHEATED DATA SECTION   **************************************
          ! Get: ***** ENTHALPY of SUPERHEATED GAS  *****
          ! First find the number of pressure value syntax lines have been entered and
          ! make sure that all of the pressure input is linked to the same temperature list
    TempsName     = " "
    FirstSHMatch  = .TRUE.
    NumOfPressPts = 0
    DO InData = 1, NumOfSHFluidPropArrays

      CALL GetObjectItem('FluidPropertySuperheated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
      
      IF ((SameString(Alphas(1),RefrigData(Loop)%Name)).AND.(SameString(Alphas(2),Enthalpy))) THEN
        NumOfPressPts = NumOfPressPts + 1
        IF (FirstSHMatch) THEN
          TempsName = Alphas(3)
          FirstSHMatch = .FALSE.
        ELSE
         
        END IF
      END IF
    END DO
   
          ! Now allocate the arrays and read the data into the proper place
          ! First, allocate the temperature array and transfer the data from the FluidTemp array
    DO TempLoop = 1, NumOfFluidTempArrays
      IF (SameString(TempsName,FluidTemps(TempLoop)%Name)) THEN
        RefrigData(Loop)%NumSuperTempPts = FluidTemps(TempLoop)%NumOfTemps
        ALLOCATE(RefrigData(Loop)%SHTemps(RefrigData(Loop)%NumSuperTempPts))
        RefrigData(Loop)%SHTemps = FluidTemps(TempLoop)%Temps
        EXIT ! the TempLoop DO loop
      END IF
      
    END DO

          ! Next, allocate the pressure related arrays
    RefrigData(Loop)%NumSuperPressPts = NumOfPressPts
    ALLOCATE(RefrigData(Loop)%SHPress(RefrigData(Loop)%NumSuperPressPts))
    ALLOCATE(RefrigData(Loop)%HshValues(RefrigData(Loop)%NumSuperTempPts,RefrigData(Loop)%NumSuperPressPts))

          ! Finally, get the pressure and enthalpy values from the user input
    NumOfPressPts = 0
    DO InData = 1, NumOfSHFluidPropArrays

      CALL GetObjectItem('FluidPropertySuperheated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
      
      IF ((SameString(Alphas(1),RefrigData(Loop)%Name)).AND.(SameString(Alphas(2),Enthalpy))) THEN
        NumOfPressPts = NumOfPressPts + 1
        RefrigData(Loop)%SHPress(NumOfPressPts) = Numbers(1)
        ! a little error trapping

        IF ((NumNumbers-1) == RefrigData(Loop)%NumSuperTempPts) THEN
          RefrigData(Loop)%HshValues(1:RefrigData(Loop)%NumSuperTempPts,NumOfPressPts) = Numbers(2:NumNumbers)
        ELSE
          CALL ShowSevereError('GetFluidPropertiesData: Number of superheated enthalpy data points '// &
                               'not equal to number of temperature points')
          CALL ShowContinueError('Error occurs in Refrigerant Data Name='//TRIM(RefrigData(Loop)%Name))
          ErrorsFound=.true.
        END IF
      END IF
    END DO

          ! Get: ***** DENSITY of SUPERHEATED GAS  *****
          ! First find the number of pressure value syntax lines have been entered and
          ! make sure that all of the pressure input is linked to the same temperature list
          ! Then allocate the arrays and read the data into the proper place
    ALLOCATE(RefrigData(Loop)%RhoshValues(RefrigData(Loop)%NumSuperTempPts,RefrigData(Loop)%NumSuperPressPts))

    NumOfPressPts = 0
    DO InData = 1, NumOfSHFluidPropArrays

      CALL GetObjectItem('FluidPropertySuperheated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
      
      IF ((SameString(Alphas(1),RefrigData(Loop)%Name)).AND.(SameString(Alphas(2),Density))) THEN
        NumOfPressPts = NumOfPressPts + 1       
        
        IF ((NumNumbers-1) == RefrigData(Loop)%NumSuperTempPts) THEN
          RefrigData(Loop)%RhoshValues(1:RefrigData(Loop)%NumSuperTempPts,NumOfPressPts) = Numbers(2:NumNumbers)
        ELSE
          CALL ShowSevereError('GetFluidPropertiesData: Number of superheated density data points not equal to '// &
                               'number of temperature points')
          CALL ShowContinueError('Error occurs in Refrigerant Data Name='//TRIM(RefrigData(Loop)%Name))
          ErrorsFound=.true.
        END IF
      END IF
    END DO
   
   ! Get: ***** ENTROPY of SUPERHEATED GAS  *****
          ! First find the number of pressure value syntax lines have been entered and
          ! make sure that all of the pressure input is linked to the same temperature list
          ! Then allocate the arrays and read the data into the proper place
    ALLOCATE(RefrigData(Loop)%SshValues(RefrigData(Loop)%NumSuperTempPts,RefrigData(Loop)%NumSuperPressPts))

    NumOfPressPts = 0
    DO InData = 1, NumOfSHFluidPropArrays

      CALL GetObjectItem('FluidPropertySuperheated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
      
      IF ((SameString(Alphas(1),RefrigData(Loop)%Name)).AND.(SameString(Alphas(2),Entropy))) THEN
        NumOfPressPts = NumOfPressPts + 1
       
        IF ((NumNumbers-1) == RefrigData(Loop)%NumSuperTempPts) THEN
          RefrigData(Loop)%SshValues(1:RefrigData(Loop)%NumSuperTempPts,NumOfPressPts) = Numbers(2:NumNumbers)
        ELSE
          CALL ShowSevereError('GetFluidPropertiesData: Number of superheated density data points not equal to '// &
                               'number of temperature points')
          CALL ShowContinueError('Error occurs in Refrigerant Data Name='//TRIM(RefrigData(Loop)%Name))
          ErrorsFound=.true.
        END IF
      END IF
    END DO

	! Get: ***** CONDUCTIVITY of SUPERHEATED GAS  *****
          ! First find the number of pressure value syntax lines have been entered and
          ! make sure that all of the pressure input is linked to the same temperature list
          ! Then allocate the arrays and read the data into the proper place
    ALLOCATE(RefrigData(Loop)%CshValues(RefrigData(Loop)%NumSuperTempPts,RefrigData(Loop)%NumSuperPressPts))

    NumOfPressPts = 0
    DO InData = 1, NumOfSHFluidPropArrays

      CALL GetObjectItem('FluidPropertySuperheated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
      
      IF ((SameString(Alphas(1),RefrigData(Loop)%Name)).AND.(SameString(Alphas(2),Conductivity))) THEN
        NumOfPressPts = NumOfPressPts + 1
       IF ((NumNumbers-1) == RefrigData(Loop)%NumSuperTempPts) THEN
          RefrigData(Loop)%CshValues(1:RefrigData(Loop)%NumSuperTempPts,NumOfPressPts) = Numbers(2:NumNumbers)
        ELSE
          CALL ShowSevereError('GetFluidPropertiesData: Number of superheated dynamic viscosity data points not equal to '// &
                               'number of temperature points')
          CALL ShowContinueError('Error occurs in Refrigerant Data Name='//TRIM(RefrigData(Loop)%Name))
          ErrorsFound=.true.
        END IF  
        
      END IF
    END DO

	! Get: ***** DYNAMIC VISCOSITY of SUPERHEATED GAS  *****
          ! First find the number of pressure value syntax lines have been entered and
          ! make sure that all of the pressure input is linked to the same temperature list
          ! Then allocate the arrays and read the data into the proper place
    ALLOCATE(RefrigData(Loop)%DVshValues(RefrigData(Loop)%NumSuperTempPts,RefrigData(Loop)%NumSuperPressPts))

    NumOfPressPts = 0
    DO InData = 1, NumOfSHFluidPropArrays

      CALL GetObjectItem('FluidPropertySuperheated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
      
      IF ((SameString(Alphas(1),RefrigData(Loop)%Name)).AND.(SameString(Alphas(2),Viscosity))) THEN
        NumOfPressPts = NumOfPressPts + 1
        
        IF ((NumNumbers-1) == RefrigData(Loop)%NumSuperTempPts) THEN
          RefrigData(Loop)%DVshValues(1:RefrigData(Loop)%NumSuperTempPts,NumOfPressPts) = Numbers(2:NumNumbers)
        ELSE
          CALL ShowSevereError('GetFluidPropertiesData: Number of superheated dynamic viscosity data points not equal to '// &
                               'number of temperature points')
          CALL ShowContinueError('Error occurs in Refrigerant Data Name='//TRIM(RefrigData(Loop)%Name))
          ErrorsFound=.true.
        END IF
      END IF
    END DO

	! Get: ***** SPECIFIC HEAT of SUPERHEATED GAS  *****
          ! First find the number of pressure value syntax lines have been entered and
          ! make sure that all of the pressure input is linked to the same temperature list
          ! Then allocate the arrays and read the data into the proper place
    ALLOCATE(RefrigData(Loop)%CPshValues(RefrigData(Loop)%NumSuperTempPts,RefrigData(Loop)%NumSuperPressPts))

    NumOfPressPts = 0
    DO InData = 1, NumOfSHFluidPropArrays

      CALL GetObjectItem('FluidPropertySuperheated',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
      
      IF ((SameString(Alphas(1),RefrigData(Loop)%Name)).AND.(SameString(Alphas(2),SpecificHeat))) THEN
        NumOfPressPts = NumOfPressPts + 1
   
        IF ((NumNumbers-1) == RefrigData(Loop)%NumSuperTempPts) THEN
          RefrigData(Loop)%CPshValues(1:RefrigData(Loop)%NumSuperTempPts,NumOfPressPts) = Numbers(2:NumNumbers)
        ELSE
          CALL ShowSevereError('GetFluidPropertiesData: Number of superheated specific heat data points not equal to '// &
                               'number of temperature points')
          CALL ShowContinueError('Error occurs in Refrigerant Data Name='//TRIM(RefrigData(Loop)%Name))
          ErrorsFound=.true.
        END IF
      END IF
    END DO

  !*****************************   SUBCOOLED DATA SECTION   ***********************************
          ! Get: ***** ENTHALPY of SUBCOOLED GAS  *****
          ! First find the number of pressure value syntax lines have been entered and
          ! make sure that all of the pressure input is linked to the same temperature list
    TempsName     = " "
    FirstSHMatch  = .TRUE.
    NumOfPressPts = 0
    DO InData = 1, NumOfSCFluidPropArrays

      CALL GetObjectItem('FluidPropertySubcooled',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
      
      IF ((SameString(Alphas(1),RefrigData(Loop)%Name)).AND.(SameString(Alphas(2),Enthalpy))) THEN
        NumOfPressPts = NumOfPressPts + 1
        IF (FirstSHMatch) THEN
          TempsName = Alphas(3)
          FirstSHMatch = .FALSE.
        ELSE
          IF (.NOT.SameString(TempsName,Alphas(3))) THEN
            CALL ShowSevereError('GetFluidPropertiesData: All subcooled data for the same property must use '// &
                                 'the same temperature list')
            CALL ShowContinueError('Error occurs in Refrigerant Data Name='//TRIM(RefrigData(Loop)%Name))
            CALL ShowContinueError('Expected Temperature name='//TRIM(TempsName)//', Input had name='//TRIM(Alphas(3)))
            ErrorsFound=.true.
          ENDIF
        END IF
      END IF
    END DO

          ! Now allocate the arrays and read the data into the proper place
          ! First, allocate the temperature array and transfer the data from the FluidTemp array
    DO TempLoop = 1, NumOfFluidTempArrays
      IF (SameString(TempsName,FluidTemps(TempLoop)%Name)) THEN
        RefrigData(Loop)%NumSubcoolTempPts = FluidTemps(TempLoop)%NumOfTemps
        ALLOCATE(RefrigData(Loop)%SCTemps(RefrigData(Loop)%NumSubcoolTempPts))
        RefrigData(Loop)%SCTemps = FluidTemps(TempLoop)%Temps
        EXIT ! the TempLoop DO loop
      END IF    
    END DO

          ! Next, allocate the pressure related arrays
    RefrigData(Loop)%NumSubcoolPressPts = NumOfPressPts
    ALLOCATE(RefrigData(Loop)%SCPress(RefrigData(Loop)%NumSubcoolPressPts))
    ALLOCATE(RefrigData(Loop)%HscValues(RefrigData(Loop)%NumSubcoolTempPts,RefrigData(Loop)%NumSubcoolPressPts))

          ! Finally, get the pressure and enthalpy values from the user input
    NumOfPressPts = 0
    DO InData = 1, NumOfSCFluidPropArrays

      CALL GetObjectItem('FluidPropertySubcooled',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
      
      IF ((SameString(Alphas(1),RefrigData(Loop)%Name)).AND.(SameString(Alphas(2),Enthalpy))) THEN
        NumOfPressPts = NumOfPressPts + 1
        RefrigData(Loop)%SCPress(NumOfPressPts) = Numbers(1)
        ! a little error trapping
        IF (NumOfPressPts == 1) THEN
          IF (RefrigData(Loop)%SCPress(NumOfPressPts) <= 0.0) THEN
            CALL ShowSevereError('GetFluidPropertiesData: Negative pressures not allowed in fluid property input data')
            CALL ShowContinueError('Error occurs in Refrigerant Data Name='//TRIM(RefrigData(Loop)%Name))
            ErrorsFound=.true.
          ENDIF
        ELSE
        END IF

        IF ((NumNumbers-1) == RefrigData(Loop)%NumSubcoolTempPts) THEN
          RefrigData(Loop)%HscValues(1:RefrigData(Loop)%NumSubcoolTempPts,NumOfPressPts) = Numbers(2:NumNumbers)
        ELSE
          CALL ShowSevereError('GetFluidPropertiesData: Number of subcooled enthalpy data points '// &
                               'not equal to number of temperature points')
          CALL ShowContinueError('Error occurs in Refrigerant Data Name='//TRIM(RefrigData(Loop)%Name))
          ErrorsFound=.true.
        END IF
      END IF
    END DO

          ! Get: ***** DENSITY of SUBCOOLED GAS  *****
          ! First find the number of pressure value syntax lines have been entered and
          ! make sure that all of the pressure input is linked to the same temperature list
          ! Then allocate the arrays and read the data into the proper place
    ALLOCATE(RefrigData(Loop)%RhoscValues(RefrigData(Loop)%NumSubcoolTempPts,RefrigData(Loop)%NumSubcoolPressPts))

    NumOfPressPts = 0
    DO InData = 1, NumOfSCFluidPropArrays

      CALL GetObjectItem('FluidPropertySubcooled',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
      
      IF ((SameString(Alphas(1),RefrigData(Loop)%Name)).AND.(SameString(Alphas(2),Density))) THEN
        NumOfPressPts = NumOfPressPts + 1
        
        IF ((NumNumbers-1) == RefrigData(Loop)%NumSubcoolTempPts) THEN
          RefrigData(Loop)%RhoscValues(1:RefrigData(Loop)%NumSubcoolTempPts,NumOfPressPts) = Numbers(2:NumNumbers)
        ELSE
          CALL ShowSevereError('GetFluidPropertiesData: Number of subcooled density data points not equal to '// &
                               'number of temperature points')
          CALL ShowContinueError('Error occurs in Refrigerant Data Name='//TRIM(RefrigData(Loop)%Name))
          ErrorsFound=.true.
        END IF
      END IF
    END DO

   ! Get: ***** ENTROPY of SUBCOOLED GAS  *****
          ! First find the number of pressure value syntax lines have been entered and
          ! make sure that all of the pressure input is linked to the same temperature list
          ! Then allocate the arrays and read the data into the proper place
    ALLOCATE(RefrigData(Loop)%SscValues(RefrigData(Loop)%NumSubcoolTempPts,RefrigData(Loop)%NumSubcoolPressPts))

    NumOfPressPts = 0
    DO InData = 1, NumOfSCFluidPropArrays

      CALL GetObjectItem('FluidPropertySubcooled',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
      
      IF ((SameString(Alphas(1),RefrigData(Loop)%Name)).AND.(SameString(Alphas(2),Entropy))) THEN
        NumOfPressPts = NumOfPressPts + 1
        
        IF ((NumNumbers-1) == RefrigData(Loop)%NumSubcoolTempPts) THEN
          RefrigData(Loop)%SscValues(1:RefrigData(Loop)%NumSubcoolTempPts,NumOfPressPts) = Numbers(2:NumNumbers)
        ELSE
          CALL ShowSevereError('GetFluidPropertiesData: Number of subcooled density data points not equal to '// &
                               'number of temperature points')
          CALL ShowContinueError('Error occurs in Refrigerant Data Name='//TRIM(RefrigData(Loop)%Name))
          ErrorsFound=.true.
        END IF
      END IF
    END DO

	! Get: ***** CONDUCTIVITY of SUBCOOLED GAS  *****
          ! First find the number of pressure value syntax lines have been entered and
          ! make sure that all of the pressure input is linked to the same temperature list
          ! Then allocate the arrays and read the data into the proper place
    ALLOCATE(RefrigData(Loop)%CscValues(RefrigData(Loop)%NumSubcoolTempPts,RefrigData(Loop)%NumSubcoolPressPts))

    NumOfPressPts = 0
    DO InData = 1, NumOfSCFluidPropArrays

      CALL GetObjectItem('FluidPropertySubcooled',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
      
      IF ((SameString(Alphas(1),RefrigData(Loop)%Name)).AND.(SameString(Alphas(2),Conductivity))) THEN
        NumOfPressPts = NumOfPressPts + 1
        
        IF ((NumNumbers-1) == RefrigData(Loop)%NumSubcoolTempPts) THEN
          RefrigData(Loop)%CscValues(1:RefrigData(Loop)%NumSubcoolTempPts,NumOfPressPts) = Numbers(2:NumNumbers)
        ELSE
          CALL ShowSevereError('GetFluidPropertiesData: Number of subcooled conductivity data points not equal to '// &
                               'number of temperature points')
          CALL ShowContinueError('Error occurs in Refrigerant Data Name='//TRIM(RefrigData(Loop)%Name))
          ErrorsFound=.true.
        END IF
      END IF
    END DO

	! Get: ***** DYNAMIC VISCOSITY of SUBCOOLED GAS  *****
          ! First find the number of pressure value syntax lines have been entered and
          ! make sure that all of the pressure input is linked to the same temperature list
          ! Then allocate the arrays and read the data into the proper place
    ALLOCATE(RefrigData(Loop)%DVscValues(RefrigData(Loop)%NumSubcoolTempPts,RefrigData(Loop)%NumSubcoolPressPts))

    NumOfPressPts = 0
    DO InData = 1, NumOfSCFluidPropArrays

      CALL GetObjectItem('FluidPropertySubcooled',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
      
      IF ((SameString(Alphas(1),RefrigData(Loop)%Name)).AND.(SameString(Alphas(2),Viscosity))) THEN
        NumOfPressPts = NumOfPressPts + 1
        
        IF ((NumNumbers-1) == RefrigData(Loop)%NumSubcoolTempPts) THEN
          RefrigData(Loop)%DVscValues(1:RefrigData(Loop)%NumSubcoolTempPts,NumOfPressPts) = Numbers(2:NumNumbers)
        ELSE
          CALL ShowSevereError('GetFluidPropertiesData: Number of subcooled dynamic viscosity data points not equal to '// &
                               'number of temperature points')
          CALL ShowContinueError('Error occurs in Refrigerant Data Name='//TRIM(RefrigData(Loop)%Name))
          ErrorsFound=.true.
        END IF
      END IF
    END DO
    
	! Get: ***** SPECIFIC HEAT of SUBCOOLED GAS  *****
          ! First find the number of pressure value syntax lines have been entered and
          ! make sure that all of the pressure input is linked to the same temperature list
          ! Then allocate the arrays and read the data into the proper place
    ALLOCATE(RefrigData(Loop)%CPscValues(RefrigData(Loop)%NumSubcoolTempPts,RefrigData(Loop)%NumSubcoolPressPts))

    NumOfPressPts = 0
    DO InData = 1, NumOfSCFluidPropArrays

      CALL GetObjectItem('FluidPropertySubcooled',InData,Alphas,NumAlphas, &
                          TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
      
      IF ((SameString(Alphas(1),RefrigData(Loop)%Name)).AND.(SameString(Alphas(2),SpecificHeat))) THEN
        NumOfPressPts = NumOfPressPts + 1
       
        IF ((NumNumbers-1) == RefrigData(Loop)%NumSubcoolTempPts) THEN
          RefrigData(Loop)%CPscValues(1:RefrigData(Loop)%NumSubcoolTempPts,NumOfPressPts) = Numbers(2:NumNumbers)
        ELSE
          CALL ShowSevereError('GetFluidPropertiesData: Number of subcooled specific heat data points not equal to '// &
                               'number of temperature points')
          CALL ShowContinueError('Error occurs in Refrigerant Data Name='//TRIM(RefrigData(Loop)%Name))
          ErrorsFound=.true.
        END IF
      END IF
    END DO
   
  END DO    ! ...end of DO loop through all of the refrigerants

  DEALLOCATE(FluidTemps)

  IF (ErrorsFound) THEN
    CALL ShowFatalError('GetFluidPropertiesData: Previous errors in input cause program termination.')
  ENDIF

  RETURN

END SUBROUTINE GetFluidPropertiesData

!**************************************************************************
REAL FUNCTION Pcrit(Refrigerant)
        ! SUBROUTINE INFORMATION:
        !       AUTHOR         Kenneth Tang
        !       DATE WRITTEN   13 December 03
        !       MODIFIED       
        !                      
		!                      
        !       RE-ENGINEERED  na

        ! PURPOSE OF THIS FUNCTION:
        ! This finds the critical pressure (Pa) for a valid refrigerant


        ! METHODOLOGY EMPLOYED:


        ! REFERENCES:
        ! na

        ! USE STATEMENTS:
        ! na

  IMPLICIT NONE           ! Enforce explicit typing of all variables in this routine

        ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN)  :: Refrigerant ! carries in substance name

SELECT CASE (Refrigerant)
CASE ('R22') 
Pcrit=4990E+3
CASE ('R123') 
Pcrit=3661.8E+3
CASE ('R134A') 
Pcrit=4059.28E+3
CASE ('R12') 
Pcrit=4136.1E+3
CASE ('R11') 
Pcrit=4407.638E+3
CASE ('NH3') 
Pcrit=11333.001E+3
CASE ('R410A') 
Pcrit=4769.893E+3
CASE ('R407C') 
Pcrit=4634.046E+3
CASE ('R417A')
Pcrit=4124.385e+3
CASE ('R507A')
Pcrit=3714.739e+3
END SELECT
RETURN
END FUNCTION Pcrit

!**************************************************************************
REAL FUNCTION Tcrit(Refrigerant)
        ! SUBROUTINE INFORMATION:
        !       AUTHOR         Kenneth Tang
        !       DATE WRITTEN   13 December 03
        !       MODIFIED       
        !                      
		!                      
        !       RE-ENGINEERED  na

        ! PURPOSE OF THIS FUNCTION:
        ! This finds the critical temperature (degC) for a valid refrigerant


        ! METHODOLOGY EMPLOYED:


        ! REFERENCES:
        ! na

        ! USE STATEMENTS:
        ! na

  IMPLICIT NONE           ! Enforce explicit typing of all variables in this routine

        ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN)  :: Refrigerant ! carries in substance name

SELECT CASE (Refrigerant)
CASE ('R22') 
Tcrit=96.145
CASE ('R123') 
Tcrit=183.681
CASE ('R134A') 
Tcrit=101.060
CASE ('R12') 
Tcrit=111.970
CASE ('R11') 
Tcrit=197.960
CASE ('NH3') 
Tcrit=132.25
CASE ('R410A') 
Tcrit=70.165
CASE ('R407C') 
Tcrit=86.049
CASE ('R417A')
Tcrit=90.734
CASE ('R507A')
Tcrit=70.745
END SELECT
RETURN
END FUNCTION Tcrit

!**************************************************************************
REAL FUNCTION MW(Refrigerant)
        ! SUBROUTINE INFORMATION:
        !       AUTHOR         Kenneth Tang
        !       DATE WRITTEN   13 December 03
        !       MODIFIED       
        !                      
		!                      
        !       RE-ENGINEERED  na

        ! PURPOSE OF THIS FUNCTION:
        ! This finds the critical molecular weight(kg/mol) for a valid refrigerant


        ! METHODOLOGY EMPLOYED:


        ! REFERENCES:
        ! na

        ! USE STATEMENTS:
        ! na

  IMPLICIT NONE           ! Enforce explicit typing of all variables in this routine

        ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN)  :: Refrigerant ! carries in substance name

SELECT CASE (Refrigerant)
CASE ('R22') 
MW=0.086
CASE ('R123') 
MW=0.153
CASE ('R134A') 
MW=0.102 
CASE ('R12') 
MW=0.121 
CASE ('R11') 
MW=0.137
CASE ('NH3') 
MW=0.017
CASE ('R410A') 
MW=0.073
CASE ('R407C') 
MW=0.086
CASE ('R417A')
MW=0.106
CASE ('R507A')
MW=0.099
END SELECT
RETURN
END FUNCTION MW
!***************************************************************************
REAL FUNCTION PQ(Refrigerant,Pressure,Quality,Property,RefrigIndex,Error)

        ! SUBROUTINE INFORMATION:
        !       AUTHOR         Mike Turner
        !       DATE WRITTEN   10 December 99
        !       MODIFIED       Rick Strand (April 2000, May 2000)
        !                      Simon Rees (May 2002)
		!                      
        !       RE-ENGINEERED  Kenneth Tang (Nov 2003)
		!                      Ipseng Iu (Nov 2006)

        ! PURPOSE OF THIS FUNCTION:
        ! This finds enthalpy for given temperature and a quality under the vapor dome.
        ! This function is only called with a valid refrigerant and quality between 0 and 1.

        ! METHODOLOGY EMPLOYED:
        ! Calls GetInterpolatedSatProp to linearly interpolate between the saturated 
        ! liquid and vapour enthalpies according to the given quality.

        ! REFERENCES:
        ! na

        ! USE STATEMENTS:
        ! na

  IMPLICIT NONE           ! Enforce explicit typing of all variables in this routine

        ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN)  :: Refrigerant ! carries in substance name
  REAL,             INTENT(IN)  :: Pressure    ! actual pressure given as input
  !REAL,             INTENT(IN)  :: Quality     ! actual quality given as input
  REAL                          :: Quality  !RS: Debugging: Allowing quality to be redefined
  CHARACTER(len=*), INTENT(IN)  :: Property
  INTEGER,       INTENT(INOUT)  :: RefrigIndex ! Index to Refrigerant Properties   !RS: Debugging: Removal of plethora of RefrigIndex definitions in the code 
  INTEGER(2),    INTENT(INOUT)  :: Error       ! Error flag: 0-no; otherwise-yes
               
        ! INTERFACE BLOCK SPECIFICATIONS:
        ! na

        ! DERIVED TYPE DEFINITIONS:
        ! na

        ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: RefrigNum    ! index for refrigerant under consideration

  INTEGER    :: LoPressIndex
  INTEGER    :: HiPressIndex
  REAL       :: PressInterpRatio
  REAL       :: SatTemperature             ! calculated saturated temperature
  REAL       :: TempInterpRatio
  REAL       :: LoSatProp
  REAL       :: HiSatProp
  LOGICAL    :: ErrorFlag                  ! error flag for current call

          ! FLOW:

  Error=0 !Initialize - no error
  ErrorFlag = .FALSE.

  IF (GetInput) THEN
    CALL GetFluidPropertiesData
    GetInput = .FALSE.
  END IF

  IF (NumOfRefrigerants == 0) &
    CALL ShowFatalError('GetSatEnthalpyRefrig: No refrigerants found--cannot evaluate fluid enthalpy')

  IF ((Quality < 0.0) .OR. (Quality > 1.0)) THEN
    !CALL ShowFatalError('GetSatEnthalpyRefrig: Saturated refrigerant quality must be between 0 and 1') !RS: Secret Search String
    WRITE(*,*) 'GetSatEnthalpyRefrig: Saturated refrigerant quality must be between 0 and 1'
    IF (Quality <0.0) THEN
        Quality=0.0
    ELSE
        Quality=1.0
    END IF
  END IF

  IF (RefrigIndex > 0) THEN
    RefrigNum=RefrigIndex
  ELSE
    ! Find which refrigerant (index) is being requested
    RefrigNum = FindRefrigerant(Refrigerant)
    RefrigIndex=RefrigNum
  ENDIF

! get the array indices

  !To account for temperature glide - ISI 11/06/06
  IF (Quality < 0.5) THEN
	LoPressIndex = FindArrayIndex(Pressure, RefrigData(RefrigNum)%PsfValues)
	HiPressIndex = LoPressIndex + 1
  ELSE
	LoPressIndex = FindArrayIndex(Pressure, RefrigData(RefrigNum)%PsgValues)
	HiPressIndex = LoPressIndex + 1
  END IF

  ! check for out of data bounds problems

  IF (LoPressIndex == 0) THEN
    LoPressIndex = MAX(1, LoPressIndex)
    WRITE(*,*)'Pressure out of range for saturated region' !- Commented by ISI 01/19/04
  	Error=1
    ErrorFlag = .True.
  END IF
  IF (Quality < 0.5) THEN
	IF(HiPressIndex > Size(RefrigData(RefrigNum)%PsfValues))THEN
		LoPressIndex = MAX(1, LoPressIndex)
		WRITE(*,*)'Pressure out of range for saturated region' !- Commented by ISI 01/19/04
		Error=1
  		ErrorFlag = .True.
    END IF
  ELSE
	IF(HiPressIndex > Size(RefrigData(RefrigNum)%PsgValues))THEN
		LoPressIndex = MAX(1, LoPressIndex)
		WRITE(*,*)'Pressure out of range for saturated region' !- Commented by ISI 01/19/04
		Error=1
  		ErrorFlag = .True.
    END IF
  END IF

  IF(.NOT. ErrorFlag)THEN
    ! find interpolation ratio w.r.t temperature

    IF (Quality < 0.5) THEN
		PressInterpRatio = (Pressure - RefrigData(RefrigNum)%PsfValues(LoPressIndex)) / &
							   (RefrigData(RefrigNum)%PsfValues(HiPressIndex) &
								- RefrigData(RefrigNum)%PsfValues(LoPressIndex))
    ELSE
		PressInterpRatio = (Pressure - RefrigData(RefrigNum)%PsgValues(LoPressIndex)) / &
							   (RefrigData(RefrigNum)%PsgValues(HiPressIndex) &
								- RefrigData(RefrigNum)%PsgValues(LoPressIndex))
	END IF

    ! apply final linear interpolation
    SatTemperature = RefrigData(RefrigNum)%PsTemps(LoPressIndex) + PressInterpRatio * &
                              (RefrigData(RefrigNum)%PsTemps(HiPressIndex) - &
                               RefrigData(RefrigNum)%PsTemps(LoPressIndex))

  ! Apply linear interpolation function
  	SELECT CASE (Property)
	CASE ('surfacetension')

	IF (Quality < 0.5) THEN
		PQ = RefrigData(RefrigNum)%STfValues(LoPressIndex) + PressInterpRatio * &
								  (RefrigData(RefrigNum)%STfValues(HiPressIndex) - &
								   RefrigData(RefrigNum)%STfValues(LoPressIndex))
	ELSE
		PQ = RefrigData(RefrigNum)%STgValues(LoPressIndex) + PressInterpRatio * &
								  (RefrigData(RefrigNum)%STgValues(HiPressIndex) - &
								   RefrigData(RefrigNum)%STgValues(LoPressIndex))
	END IF
	
	CASE ('temperature')
	PQ = SatTemperature
	CASE ('enthalpy')   
	!CALCULATE ENTHAPY OF SATURATED REGION
	PQ = GetInterpolatedSatProp(SatTemperature, RefrigData(RefrigNum)%HTemps, &
                                             RefrigData(RefrigNum)%HfValues, &
                       RefrigData(RefrigNum)%HfgValues, Quality)
	CASE ('entropy')  
	!CALCULATE ENTROPY OF SATURATED REGION	 
	PQ = GetInterpolatedSatProp(SatTemperature, RefrigData(RefrigNum)%STemps, &
                                             RefrigData(RefrigNum)%SfValues, &
                       RefrigData(RefrigNum)%SfgValues, Quality)
	CASE ('density') 	 
	!CALCULATE DENSITY OF SATURATED REGION
  	! find adjacent property values at the given quality
	LoSatProp = 1/RefrigData(RefrigNum)%RhofValues(LoPressIndex) + &
                 Quality*(1/RefrigData(RefrigNum)%RhofgValues(LoPressIndex) - 1/RefrigData(RefrigNum)%RhofValues(LoPressIndex))
	LoSatProp = 1/LoSatProp
	HiSatProp = 1/RefrigData(RefrigNum)%RhofValues(HiPressIndex) + &
                 Quality*(1/RefrigData(RefrigNum)%RhofgValues(HiPressIndex) - 1/RefrigData(RefrigNum)%RhofValues(HiPressIndex))
	HiSatProp= 1/HiSatProp

	! find interpolation ratio in temperature direction
	TempInterpRatio = (SatTemperature - RefrigData(RefrigNum)%RhoTemps(LoPressIndex)) / &
                      (RefrigData(RefrigNum)%RhoTemps(HiPressIndex) - RefrigData(RefrigNum)%RhoTemps(LoPressIndex))
	! apply final linear interpolation
	PQ = 1/LoSatProp + TempInterpRatio*(1/HiSatProp - 1/LoSatProp)
	PQ = 1/PQ   

    CASE ('viscosity') 	 
	!CALCULATE VISCOSITY OF SATURATED REGION
	IF (Quality > 0.5) THEN
    PQ= RefrigData(RefrigNum)%DVfgValues(LoPressIndex) + PressInterpRatio * &
                              (RefrigData(RefrigNum)%DVfgValues(HiPressIndex) - &
                              RefrigData(RefrigNum)%DVfgValues(LoPressIndex))
	ELSEIF (Quality <= 0.5) THEN
	PQ= RefrigData(RefrigNum)%DVfValues(LoPressIndex) + PressInterpRatio * &
                              (RefrigData(RefrigNum)%DVfValues(HiPressIndex) - &
                              RefrigData(RefrigNum)%DVfValues(LoPressIndex))
	ELSE
	PQ=0
	END IF

	CASE ('conductivity') 	 
	!CALCULATE CONDUCTIVITY OF SATURATED REGION
	IF (Quality > 0.5) THEN
	PQ= RefrigData(RefrigNum)%CfgValues(LoPressIndex) + PressInterpRatio * &
                              (RefrigData(RefrigNum)%CfgValues(HiPressIndex) - &
                              RefrigData(RefrigNum)%CfgValues(LoPressIndex))
    
	ELSEIF (Quality <= 0.5) THEN
    PQ= RefrigData(RefrigNum)%CfValues(LoPressIndex) + PressInterpRatio * &
                              (RefrigData(RefrigNum)%CfValues(HiPressIndex) - &
                              RefrigData(RefrigNum)%CfValues(LoPressIndex))
	ELSE
	PQ=0
	END IF

	CASE ('specificheat') 	 
	!CALCULATE SPECIFIC HEAT OF SATURATED REGION
	IF (Quality > 0.5) THEN
	PQ= RefrigData(RefrigNum)%CpfgValues(LoPressIndex) + PressInterpRatio * &
                              (RefrigData(RefrigNum)%CpfgValues(HiPressIndex) - &
                              RefrigData(RefrigNum)%CpfgValues(LoPressIndex))

	ELSEIF (Quality <= 0.5) THEN
    PQ= RefrigData(RefrigNum)%CpfValues(LoPressIndex) + PressInterpRatio * &
                              (RefrigData(RefrigNum)%CpfValues(HiPressIndex) - &
                              RefrigData(RefrigNum)%CpfValues(LoPressIndex))
	ELSE
	PQ=0
	END IF

	END SELECT
  END IF
  RETURN
END FUNCTION PQ

!*******************************************************************************
REAL FUNCTION TQ(Refrigerant,Temperature,Quality,Property,RefrigIndex,Error)

        ! SUBROUTINE INFORMATION:
        !       AUTHOR         Mike Turner
        !       DATE WRITTEN   10 December 99
        !       MODIFIED       Rick Strand (April 2000, May 2000)
        !                      Simon Rees (May 2002)
        !                      
        !       RE-ENGINEERED  Kenneth Tang (Nov 2003)
		!                      Ipseng Iu (Nov 2006)

        ! PURPOSE OF THIS FUNCTION:
        ! This finds enthalpy for given temperature and a quality under the vapor dome.
        ! This function is only called with a valid refrigerant and quality between 0 and 1.

        ! METHODOLOGY EMPLOYED:
        ! Calls GetInterpolatedSatProp to linearly interpolate between the saturated 
        ! liquid  and vapour enthalpies according to the given quality.

        ! REFERENCES:
        ! na

        ! USE STATEMENTS:
        ! na

  IMPLICIT NONE           ! Enforce explicit typing of all variables in this routine

        ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN)  :: Refrigerant ! carries in substance name
  REAL,             INTENT(IN)  :: Temperature ! actual temperature given as input
  REAL,             INTENT(IN)  :: Quality     ! actual quality given as input
  INTEGER,       INTENT(INOUT)  :: RefrigIndex ! Index to Refrigerant Properties   !RS: Debugging: Removal of plethora of RefrigIndex definitions in the code  
  CHARACTER(len=*), INTENT(IN)  :: Property
  INTEGER(2),      INTENT(OUT)  :: Error       !Error flag: 0-no; otherwise-yes

        ! INTERFACE BLOCK SPECIFICATIONS:
        ! na

        ! DERIVED TYPE DEFINITIONS:
        ! na

        ! FUNCTION LOCAL VARIABLE DECLARATIONS:
		
  INTEGER :: RefrigNum    ! index for refrigerant under consideration
  REAL    :: TempInterpRatio
  INTEGER :: LoTempIndex
  INTEGER :: HiTempIndex
  LOGICAL :: ErrorFlag                  ! error flag for current call

          ! FLOW:
  Error=0
  ErrorFlag = .FALSE.

  IF (GetInput) THEN
    CALL GetFluidPropertiesData
    GetInput = .FALSE.
  END IF

  IF (NumOfRefrigerants == 0) &
    CALL ShowFatalError('GetSatEnthalpyRefrig: No refrigerants found--cannot evaluate fluid enthalpy')

  IF ((Quality < 0.0) .OR. (Quality > 1.0)) &
    CALL ShowFatalError('GetSatEnthalpyRefrig: Saturated refrigerant quality must be between 0 and 1')

  IF (RefrigIndex > 0) THEN
    RefrigNum=RefrigIndex
  ELSE
    ! Find which refrigerant (index) is being requested
    RefrigNum = FindRefrigerant(Refrigerant)
    RefrigIndex=RefrigNum
  ENDIF

  LoTempIndex = FindArrayIndex(Temperature, RefrigData(RefrigNum)%PsTemps)
  HiTempIndex = LoTempIndex + 1

  ! check for out of data bounds problems
  IF (LoTempIndex == 0) THEN
    LoTempIndex = MAX(1, LoTempIndex)
	IF (Quality < 0.5) THEN
		TQ = RefrigData(RefrigNum)%PsfValues(LoTempIndex)
	ELSE
		TQ = RefrigData(RefrigNum)%PsgValues(LoTempIndex)
	END IF
    ErrorFlag = .True.
    WRITE(*,*)'Pressure out of bounds' !- Commented by ISI 01/19/04
	Error=1

  ELSE IF(HiTempIndex > Size(RefrigData(RefrigNum)%PsTemps))THEN
    LoTempIndex = MAX(1, LoTempIndex)
	IF (Quality < 0.5) THEN
		TQ = RefrigData(RefrigNum)%PsfValues(LoTempIndex)
	ELSE
		TQ = RefrigData(RefrigNum)%PsgValues(LoTempIndex)
	END IF
    ErrorFlag = .True.
	WRITE(*,*)'Pressure out of bounds' !- Commented by ISI 01/19/04
	Error=1
  END IF
  
  IF (Error .GT. 0) THEN
      RETURN
  END IF

  IF(.NOT. ErrorFlag)THEN
    ! find interpolation ratio w.r.t temperature
    TempInterpRatio = (Temperature - RefrigData(RefrigNum)%PsTemps(LoTempIndex)) / &
                      (RefrigData(RefrigNum)%PsTemps(HiTempIndex) &
                       - RefrigData(RefrigNum)%PsTemps(LoTempIndex))

  END IF

  SELECT CASE (Property)
  CASE ('surfacetension')
  
  IF (Quality < 0.5) THEN
	  TQ = RefrigData(RefrigNum)%STfValues(LoTempIndex) + TempInterpRatio * &
								  (RefrigData(RefrigNum)%STfValues(HiTempIndex) - &
								   RefrigData(RefrigNum)%STfValues(LoTempIndex))
  ELSE
	  TQ = RefrigData(RefrigNum)%STgValues(LoTempIndex) + TempInterpRatio * &
								  (RefrigData(RefrigNum)%STgValues(HiTempIndex) - &
								   RefrigData(RefrigNum)%STgValues(LoTempIndex))
  END IF

  CASE ('enthalpy') 
  !CALCULATE ENTHALPY OF SATURATED REGION GIVEN QUALITY
  TQ = GetInterpolatedSatProp(Temperature, RefrigData(RefrigNum)%HTemps, &
                                             RefrigData(RefrigNum)%HfValues, &
                       RefrigData(RefrigNum)%HfgValues, Quality)
  CASE ('entropy') 
  !CALCULATE ENTROPY OF SATURATED REGION GIVEN QUALITY
  TQ = GetInterpolatedSatProp(Temperature, RefrigData(RefrigNum)%STemps, &
                                             RefrigData(RefrigNum)%SfValues, &
                       RefrigData(RefrigNum)%SfgValues, Quality)
  
  CASE ('pressure') 
  !CALCULATE PRESSURE OF SATURATED REGION 
  ! apply final linear interpolation

  IF (Quality < 0.5) THEN
	TQ = RefrigData(RefrigNum)%PsfValues(LoTempIndex) + TempInterpRatio * &
		                          (RefrigData(RefrigNum)%PsfValues(HiTempIndex) - RefrigData(RefrigNum)%PsfValues(LoTempIndex))
  ELSE
	TQ = RefrigData(RefrigNum)%PsgValues(LoTempIndex) + TempInterpRatio * &
		                          (RefrigData(RefrigNum)%PsgValues(HiTempIndex) - RefrigData(RefrigNum)%PsgValues(LoTempIndex))
  END IF

  END SELECT

  RETURN

END FUNCTION TQ

!*****************************************************************************
REAL FUNCTION TP(Refrigerant,Temperature,Pressure,Property,RefrigIndex,Error)

        ! SUBROUTINE INFORMATION:
        !       AUTHOR         Mike Turner
        !       DATE WRITTEN   10 December 99
        !       MODIFIED       Rick Strand (April 2000, May 2000)
        !                      Simon Rees (May 2002)
        !                      
        !       RE-ENGINEERED  Kenneth Tang (Nov 2003)
		!                      Ipseng Iu (Nov 2006)

        ! PURPOSE OF THIS FUNCTION:
        ! Given temperature and pressure, find enthalpy,density and entropy for 
		! superheated and subcooled region only. Error warning would be shown if 
		! temperature and pressure corresponds to the saturated region.

        ! METHODOLOGY EMPLOYED:
        ! Detemine the state of the properties 

        ! REFERENCES:
        ! na

        ! USE STATEMENTS:
        ! na

 IMPLICIT NONE           ! Enforce explicit typing of all variables in this routine

    ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN)  :: Refrigerant ! carries in substance name
  REAL,				INTENT(IN)  :: Temperature ! actual temperature given as input
  REAL,             INTENT(IN)  :: Pressure    ! pressure(Pa)
  CHARACTER(len=*), INTENT(IN)  :: Property
  INTEGER,       INTENT(INOUT)  :: RefrigIndex ! Index to Refrigerant Properties   !RS: Debugging: Removal of plethora of RefrigIndex definitions in the code
  INTEGER(2),    INTENT(OUT)  :: Error       !Error flag: 0-no; otherwise-yes

        ! INTERFACE BLOCK SPECIFICATIONS:
        ! na

        ! DERIVED TYPE DEFINITIONS:
        ! na

        ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: TempIndex                  ! index value of Temperature from table
  INTEGER :: LoPressIndex               ! index value of next lowest Pressure from table
  INTEGER :: HiPressIndex               ! index value of next highest Pressure from table
  INTEGER :: LoPressIndexLiq            ! index value of next lowest liquid Pressure from table
  INTEGER :: HiPressIndexLiq            ! index value of next highest liquid Pressure from table
  INTEGER :: LoPressIndexVap            ! index value of next lowest vapor Pressure from table
  INTEGER :: HiPressIndexVap            ! index value of next highest vapor Pressure from table
  
  INTEGER :: HiTempIndex                ! index value of next highest Temperature from table
  INTEGER :: LoTempIndex                ! index value of next lowest Temperature from table
  INTEGER :: RefrigNum                  ! index for refrigerant under consideration
  REAL    :: PressInterpRatio           ! ratio to interpolate in Pressure domain
  REAL    :: PressInterpRatioLiq        ! ratio to interpolate in vapor Pressure domain
  REAL    :: PressInterpRatioVap        ! ratio to interpolate in liquid Pressure domain
  REAL    :: TempInterpRatio            ! ratio to interpolate in temperature domain
  REAL    :: LoTempLoProperty           ! low property at low temperature
  REAL    :: LoTempHiProperty           ! high property at low temperature
  REAL    :: HiTempLoProperty           ! low property at high temperature
  REAL    :: HiTempHiProperty           ! high property at high temperature
  REAL    :: TemperatureLow             ! low temperature value
  REAL    :: TemperatureHigh            ! high temperature value
  REAL    :: PropertyLow                ! low property value
  REAL    :: PropertyHigh               ! high property value
  REAL	  :: TSATLoTempHiProperty		! saturated T at high pressure property
  REAL    :: SatTemperatureLiq          ! value for liquid saturated temperature !ISI - 11/16/07
  REAL    :: SatTemperatureVap          ! value for vapor saturated temperature !ISI - 11/16/07 
  
  ! error counters and dummy string
  LOGICAL :: ErrorFlag                  ! error flag for current call

  Error=0
  ErrorFlag = .FALSE.
  
  IF (GetInput) THEN
    CALL GetFluidPropertiesData
    GetInput = .FALSE.
  END IF

  IF (NumOfRefrigerants == 0) &
    CALL ShowFatalError('GetSatPressureRefrig: No refrigerants found--cannot evaluate saturation pressure')

  IF (RefrigIndex > 0) THEN
    RefrigNum=RefrigIndex
  ELSE
    ! Find which refrigerant (index) is being requested
    RefrigNum = FindRefrigerant(Refrigerant)
    RefrigIndex=RefrigNum
  ENDIF

  ! get the array indices
  LoTempIndex = FindArrayIndex(Temperature, RefrigData(RefrigNum)%PsTemps)
  HiTempIndex = LoTempIndex + 1

  LoPressIndexLiq = FindArrayIndex(Pressure, RefrigData(RefrigNum)%PsfValues)
  HiPressIndexLiq = LoPressIndexLiq + 1

  LoPressIndexVap = FindArrayIndex(Pressure, RefrigData(RefrigNum)%PsgValues)
  HiPressIndexVap = LoPressIndexVap + 1

   ! check for out of data bounds problems

  IF (LoPressIndexLiq == 0) THEN
    LoPressIndexLiq = MAX(1, LoPressIndexLiq)
    ErrorFlag = .True. 
	WRITE(*,*) 'Pressure out of range' !- Commented by ISI 01/19/04
	Error=2
  ELSE IF(HiPressIndexLiq > Size(RefrigData(RefrigNum)%PsfValues))THEN
    LoPressIndexLiq = MAX(1, LoPressIndexLiq)
    ErrorFlag = .True.
	WRITE(*,*) 'Pressure out of range' !- Commented by ISI 01/19/04
	Error=2
  END IF

  IF (LoPressIndexVap == 0) THEN
    LoPressIndexVap = MAX(1, LoPressIndexVap)
    ErrorFlag = .True. 
	WRITE(*,*) 'Pressure out of range' !- Commented by ISI 01/19/04
	Error=2
  ELSE IF(HiPressIndexVap > Size(RefrigData(RefrigNum)%PsgValues))THEN
    LoPressIndexVap = MAX(1, LoPressIndexVap)
    ErrorFlag = .True.
	WRITE(*,*) 'Pressure out of range' !- Commented by ISI 01/19/04
	Error=2
  END IF

  IF(.NOT. ErrorFlag)THEN
    ! find interpolation ratio w.r.t temperature

	PressInterpRatioLiq = (Pressure - RefrigData(RefrigNum)%PsfValues(LoPressIndexLiq)) / &
                      (RefrigData(RefrigNum)%PsfValues(HiPressIndexLiq) &
                       - RefrigData(RefrigNum)%PsfValues(LoPressIndexLiq))

	PressInterpRatioVap = (Pressure - RefrigData(RefrigNum)%PsgValues(LoPressIndexVap)) / &
                      (RefrigData(RefrigNum)%PsgValues(HiPressIndexVap) &
                       - RefrigData(RefrigNum)%PsgValues(LoPressIndexVap))

    ! apply final linear interpolation
    
    SatTemperatureLiq = RefrigData(RefrigNum)%PsTemps(LoPressIndexLiq) + PressInterpRatioLiq * &
                              (RefrigData(RefrigNum)%PsTemps(HiPressIndexLiq) - RefrigData(RefrigNum)%PsTemps(LoPressIndexLiq))

    SatTemperatureVap = RefrigData(RefrigNum)%PsTemps(LoPressIndexVap) + PressInterpRatioVap * &
                              (RefrigData(RefrigNum)%PsTemps(HiPressIndexVap) - RefrigData(RefrigNum)%PsTemps(LoPressIndexVap))

   !Check for superheated and subcooled region and calculate states
   IF (Temperature >= SatTemperatureVap) THEN !ISI - 11/16/07
   !CALCULATE PROPERTY FOR SUPERHEATED REGION
      
	  LoPressIndex=LoPressIndexVap
	  HiPressIndex=HiPressIndexVap
	  PressInterpRatio=PressInterpRatioVap

	  TempIndex  = FindArrayIndex(Temperature,RefrigData(RefrigNum)%SHTemps)
	  LoPressIndex = FindArrayIndex(Pressure,RefrigData(RefrigNum)%SHPress)

      ! check temperature data range of superheated region
      IF((TempIndex > 0) .AND. (TempIndex < SIZE(RefrigData(RefrigNum)%SHTemps)) )THEN ! in range
		HiTempIndex   = TempIndex + 1
		TempInterpRatio  = (Temperature - RefrigData(RefrigNum)%SHTemps(TempIndex)) &
                              /(RefrigData(RefrigNum)%SHTemps(HiTempIndex) &
                              - RefrigData(RefrigNum)%SHTemps(TempIndex))
      ELSE
	  ErrorFlag = .True.
	  WRITE(*,*) 'Temperature out of range for superheated region','T=',Temperature !- Commented by ISI 01/19/04
	  Error=3
	  END IF

	  ! check pressure data range and attempt to cap if necessary
	  IF((LoPressIndex > 0) .AND. (LoPressIndex < SIZE(RefrigData(RefrigNum)%SHPress) ))THEN ! in range
	  	HiPressIndex = LoPressIndex + 1
	  	PressInterpRatio = (Pressure - RefrigData(RefrigNum)%SHPress(LoPressIndex)) &
                              /(RefrigData(RefrigNum)%SHPress(HiPressIndex) &
                              - RefrigData(RefrigNum)%SHPress(LoPressIndex))
      ELSE
	  ErrorFlag = .True.
      WRITE(*,*) 'Pressure out of range for superheated region','P=',Pressure !- Commented by ISI 01/19/04
	  Error=3
	  END IF

    	IF(.NOT. ErrorFlag)THEN
		SELECT CASE (Property)
		CASE ('enthalpy')  
		!CALCULATE ENTHAPY OF SUPERHEATED REGION
		LoTempLoProperty = RefrigData(RefrigNum)%HshValues(TempIndex,LoPressIndex)
        LoTempHiProperty = RefrigData(RefrigNum)%HshValues(TempIndex,HiPressIndex)
        HiTempLoProperty = RefrigData(RefrigNum)%HshValues(HiTempIndex,LoPressIndex)
        HiTempHiProperty = RefrigData(RefrigNum)%HshValues(HiTempIndex,HiPressIndex)	
		CASE ('entropy')
	    !CALCULATE ENTROPY OF SUPERHEATED REGION
		LoTempLoProperty = RefrigData(RefrigNum)%SshValues(TempIndex,LoPressIndex)
        LoTempHiProperty = RefrigData(RefrigNum)%SshValues(TempIndex,HiPressIndex)
        HiTempLoProperty = RefrigData(RefrigNum)%SshValues(HiTempIndex,LoPressIndex)
        HiTempHiProperty = RefrigData(RefrigNum)%SshValues(HiTempIndex,HiPressIndex)
		CASE ('density')
        !CALCULATE DENSITY OF SUPERHEATED REGION
        LoTempLoProperty = RefrigData(RefrigNum)%RhoshValues(TempIndex,LoPressIndex)
        LoTempHiProperty = RefrigData(RefrigNum)%RhoshValues(TempIndex,HiPressIndex)
        HiTempLoProperty = RefrigData(RefrigNum)%RhoshValues(HiTempIndex,LoPressIndex)
        HiTempHiProperty = RefrigData(RefrigNum)%RhoshValues(HiTempIndex,HiPressIndex)
		END SELECT				
		IF (LoTempHiProperty.EQ.0) THEN !Estimate the point at saturated region
		LoTempHiProperty= PQ(Refrigerant,RefrigData(RefrigNum)%SHPress(HiPressIndex),1.0,Property,RefrigIndex,Error)
		TSATLoTempHiProperty= PQ(Refrigerant,RefrigData(RefrigNum)%SHPress(HiPressIndex),1.0,'temperature',RefrigIndex,Error)
		!use pressure interpolation to calculate for the low and high temperature
		TemperatureLow=RefrigData(RefrigNum)%SHTemps(TempIndex)+PressInterpRatio*(TSATLoTempHiProperty-RefrigData(RefrigNum)%SHTemps(TempIndex))
		TemperatureHigh=RefrigData(RefrigNum)%SHTemps(HiTempIndex)
		!calculate the new temperature interpolation ratio
		TempInterpRatio=(Temperature-TemperatureLow)/(TemperatureHigh-TemperatureLow)
		END IF	
		PropertyLow = PressInterpRatio*LoTempHiProperty + &
                (1.0-PressInterpRatio)*LoTempLoProperty
        PropertyHigh = PressInterpRatio*HiTempHiProperty + &
                 (1.0-PressInterpRatio)*HiTempLoProperty
        ! interpolate w.r.t. temperature
        TP= TempInterpRatio*PropertyHigh + (1.0-TempInterpRatio)*PropertyLow
        END IF
		                    		
   ELSE IF (Temperature <= SatTemperatureLiq) THEN !ISI - 11/16/07
   !CALCULATE PROPERTY FOR SUBCOOLED REGION

	  LoPressIndex=LoPressIndexLiq
	  HiPressIndex=HiPressIndexLiq
	  PressInterpRatio=PressInterpRatioLiq

       TempIndex  = FindArrayIndexSC(Temperature,RefrigData(RefrigNum)%SCTemps)
	    LoPressIndex = FindArrayIndexSC(Pressure,RefrigData(RefrigNum)%SCPress)

      ! check temperature data range of subcooled region
      IF((TempIndex > 0) .AND. (TempIndex < SIZE(RefrigData(RefrigNum)%SCTemps)) )THEN ! in range
		HiTempIndex   = TempIndex + 1
		TempInterpRatio  = (Temperature - RefrigData(RefrigNum)%SCTemps(TempIndex)) &
                              /(RefrigData(RefrigNum)%SCTemps(HiTempIndex) &
                              - RefrigData(RefrigNum)%SCTemps(TempIndex))
      ELSE
	    ErrorFlag = .True.
	    WRITE(*,*) 'Temperature out of range for subcooled region','T=',Temperature !- Commented by ISI 01/19/04
	    Error=4
	  END IF

	  ! check pressure data range and attempt to cap if necessary
	  IF((LoPressIndex > 0) .AND. (LoPressIndex < SIZE(RefrigData(RefrigNum)%SCPress) ))THEN ! in range
		HiPressIndex = LoPressIndex + 1
		PressInterpRatio = (Pressure - RefrigData(RefrigNum)%SCPress(LoPressIndex)) &
                              /(RefrigData(RefrigNum)%SCPress(HiPressIndex) &
                              - RefrigData(RefrigNum)%SCPress(LoPressIndex))
      ELSE
	    ErrorFlag = .True.
        WRITE(*,*) 'Pressure out of range for subcooled region','P=',Pressure !- Commented by ISI 01/19/04
	    Error=6
	  END IF

		IF(.NOT. ErrorFlag)THEN
		SELECT CASE (Property)
		CASE ('enthalpy')  
    	!CALCULATE ENTHAPY OF SUBCOOLED REGION
		LoTempLoProperty = RefrigData(RefrigNum)%HscValues(TempIndex,LoPressIndex)
        LoTempHiProperty = RefrigData(RefrigNum)%HscValues(TempIndex,HiPressIndex)
        HiTempLoProperty = RefrigData(RefrigNum)%HscValues(HiTempIndex,LoPressIndex)
        HiTempHiProperty = RefrigData(RefrigNum)%HscValues(HiTempIndex,HiPressIndex)
		CASE ('entropy')
	    !CALCULATE ENTROPY OF SUBCOOLED REGION
		LoTempLoProperty = RefrigData(RefrigNum)%SscValues(TempIndex,LoPressIndex)
        LoTempHiProperty = RefrigData(RefrigNum)%SscValues(TempIndex,HiPressIndex)
        HiTempLoProperty = RefrigData(RefrigNum)%SscValues(HiTempIndex,LoPressIndex)
        HiTempHiProperty = RefrigData(RefrigNum)%SscValues(HiTempIndex,HiPressIndex)
		CASE ('density')
        !CALCULATE DENSITY OF SUBCOOLED REGION
        LoTempLoProperty = RefrigData(RefrigNum)%RhoscValues(TempIndex,LoPressIndex)
        LoTempHiProperty = RefrigData(RefrigNum)%RhoscValues(TempIndex,HiPressIndex)
        HiTempLoProperty = RefrigData(RefrigNum)%RhoscValues(HiTempIndex,LoPressIndex)
        HiTempHiProperty = RefrigData(RefrigNum)%RhoscValues(HiTempIndex,HiPressIndex)
		END SELECT
		
		    IF (LoTempHiProperty.EQ.0) THEN !Estimate the point at saturated region
		    LoTempHiProperty= PQ(Refrigerant,RefrigData(RefrigNum)%SCPress(HiPressIndex),0.0,Property,RefrigIndex,Error)
		    TSATLoTempHiProperty= PQ(Refrigerant,RefrigData(RefrigNum)%SCPress(HiPressIndex),0.0,'temperature',RefrigIndex,Error)
		    !use pressure interpolation to calculate for the low and high temperature
		    TemperatureLow=RefrigData(RefrigNum)%SCTemps(TempIndex)+PressInterpRatio*(TSATLoTempHiProperty-RefrigData(RefrigNum)%SCTemps(TempIndex))
		    TemperatureHigh=RefrigData(RefrigNum)%SCTemps(HiTempIndex)

		    !calculate the new temperature interpolation ratio
		    TempInterpRatio=(Temperature-TemperatureLow)/(TemperatureHigh-TemperatureLow)
		    END IF

		PropertyLow = PressInterpRatio*LoTempHiProperty + &
                (1.0-PressInterpRatio)*LoTempLoProperty
        PropertyHigh = PressInterpRatio*HiTempHiProperty + &
                 (1.0-PressInterpRatio)*HiTempLoProperty
        ! interpolate w.r.t. temperature
        TP= TempInterpRatio*PropertyHigh + (1.0-TempInterpRatio)*PropertyLow
		END IF

   ELSE
		SELECT CASE (property)
		CASE ('enthalpy')
		WRITE(*,*) 'Enthalpy in saturated region,Cannot be determined by T & P only' !- Commented by ISI 01/19/04
		CASE ('density')
		WRITE(*,*) 'Density in saturated region,Cannot be determined by T & P only' !- Commented by ISI 01/19/04
		CASE ('entropy')
		WRITE(*,*) 'Entropy in saturated region,Cannot be determined by T & P only' !- Commented by ISI 01/19/04
		END SELECT
   END IF
 END IF
  RETURN
END FUNCTION TP

!*****************************************************************************
REAL FUNCTION PH(Refrigerant,Pressure,Enthalpy,Property,RefrigIndex,Error)

        ! SUBROUTINE INFORMATION:
        !       AUTHOR         Mike Turner
        !       DATE WRITTEN   10 December 99
        !       MODIFIED       Rick Strand (April 2000, May 2000)
        !                      Simon Rees (May 2002)
        !                      
        !       RE-ENGINEERED  Kenneth Tang (Nov 2003)
		!                      Ipseng Iu (Nov 2006)

        ! PURPOSE OF THIS FUNCTION:
        ! Given pressure and enthalpy, find temperature,density,quality,entropy,
		! dynamic viscosity,  thermal conductivity, specific heat.

        ! METHODOLOGY EMPLOYED:
        ! Determine the state of the properties using saturated liquid and vapor 
		! enthalpy and find the property requested. 

        ! REFERENCES:
        ! na

        ! USE STATEMENTS:
        ! na

  IMPLICIT NONE           ! Enforce explicit typing of all variables in this routine

        ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN)  :: Refrigerant ! carries in substance name
  REAL,             INTENT(IN)  :: Pressure ! actual pressure given as input
  REAL,             INTENT(IN)  :: Enthalpy     ! actual enthalpy given as input
  INTEGER,       INTENT(INOUT)  :: RefrigIndex ! Index to Refrigerant Properties    !RS: Debugging: Removal of plethora of RefrigIndex definitions in the code
  CHARACTER(len=*), INTENT(IN)  :: Property
  INTEGER(2),      INTENT(OUT)  :: Error       !Error flag: 0-no; otherwise-yes

        ! INTERFACE BLOCK SPECIFICATIONS:
        ! na

        ! DERIVED TYPE DEFINITIONS:
        ! na

        ! FUNCTION LOCAL VARIABLE DECLARATIONS:
		
  INTEGER :: RefrigNum    ! index for refrigerant under consideration
  REAL    :: PressInterpRatio
  REAL    :: PressInterpRatioLiq
  REAL    :: PressInterpRatioVap
  INTEGER :: LoPressIndex
  INTEGER :: HiPressIndex
  INTEGER :: LoPressIndexLiq
  INTEGER :: HiPressIndexLiq
  INTEGER :: LoPressIndexVap
  INTEGER :: HiPressIndexVap
  INTEGER :: LoPressStart
  INTEGER :: HiPressStart
  INTEGER :: Start
  INTEGER :: Finish
  INTEGER :: PressStart
  INTEGER :: PressFinish
  INTEGER :: LoEnthalpyIndex
  INTEGER :: HiEnthalpyIndex
  REAL    :: EnthInterpRatio
  REAL    :: EnthalpyCheck
  REAL    :: EnthalpyMax
  REAL    :: EnthalpyMin
  REAL    :: SatVapEnthalpy
  REAL    :: SatLiqEnthalpy
  REAL    :: HiEnthLoProperty
  REAL    :: HiEnthHiProperty
  REAL    :: LoEnthLoProperty
  REAL    :: LoEnthHiProperty
  REAL    :: EnthalpyLow
  REAL    :: EnthalpyHigh
  REAL    :: PropertyLow
  REAL    :: PropertyHigh
  REAL	  :: SatTemperature
  REAL    :: Quality
  INTEGER :: Loop
  LOGICAL :: ErrorFlag                  ! error flag for current call
  INTEGER :: PropertyType

  Error=0
  ErrorFlag = .FALSE.
    
  IF (GetInput) THEN
    CALL GetFluidPropertiesData
    GetInput = .FALSE.
  END IF

  IF (NumOfRefrigerants == 0) THEN
    CALL ShowFatalError('GetSatTemperatureRefrig: No refrigerants found--cannot evaluate saturation temperature')
  END IF

  ErrorFlag = .False.

  IF (RefrigIndex > 0) THEN
    RefrigNum=RefrigIndex
  ELSE
    ! Find which refrigerant (index) is being requested
    RefrigNum = FindRefrigerant(Refrigerant)
    RefrigIndex=RefrigNum
  ENDIF

  ! get the array indices

  LoPressIndexLiq = FindArrayIndex(Pressure, RefrigData(RefrigNum)%PsfValues)
  HiPressIndexLiq = LoPressIndexLiq + 1

  LoPressIndexVap = FindArrayIndex(Pressure, RefrigData(RefrigNum)%PsgValues)
  HiPressIndexVap = LoPressIndexVap + 1

  ! check for out of data bounds problems

  IF (LoPressIndexLiq == 0) THEN
    LoPressIndexLiq = MAX(1, LoPressIndexLiq)
    ErrorFlag = .True. 
	WRITE(*,*) 'Pressure out of range' !- Commented by ISI 01/19/04
	Error=1
  ELSE IF(HiPressIndexLiq > Size(RefrigData(RefrigNum)%PsfValues))THEN
    LoPressIndexLiq = MAX(1, LoPressIndexLiq)
    ErrorFlag = .True.
	WRITE(*,*) 'Pressure out of range' !- Commented by ISI 01/19/04
	Error=1
  END IF

  IF (LoPressIndexVap == 0) THEN
    LoPressIndexVap = MAX(1, LoPressIndexVap)
    ErrorFlag = .True. 
	WRITE(*,*) 'Pressure out of range' !- Commented by ISI 01/19/04
	Error=1
  ELSE IF(HiPressIndexVap > Size(RefrigData(RefrigNum)%PsgValues))THEN
    LoPressIndexVap = MAX(1, LoPressIndexVap)
    ErrorFlag = .True.
	WRITE(*,*) 'Pressure out of range' !- Commented by ISI 01/19/04
	Error=1
  END IF

  IF(.NOT. ErrorFlag)THEN

  PressInterpRatioLiq = (Pressure - RefrigData(RefrigNum)%PsfValues(LoPressIndexLiq)) / &
                      (RefrigData(RefrigNum)%PsfValues(HiPressIndexLiq) &
                       - RefrigData(RefrigNum)%PsfValues(LoPressIndexLiq))

  PressInterpRatioVap = (Pressure - RefrigData(RefrigNum)%PsgValues(LoPressIndexVap)) / &
                      (RefrigData(RefrigNum)%PsgValues(HiPressIndexVap) &
                       - RefrigData(RefrigNum)%PsgValues(LoPressIndexVap))

  SatVapEnthalpy = RefrigData(RefrigNum)%HfgValues(LoPressIndexVap) &
                 + PressInterpRatioVap*(RefrigData(RefrigNum)%HfgValues(HiPressIndexVap) &
				    - RefrigData(RefrigNum)%HfgValues(LoPressIndexVap))
 
  SatLiqEnthalpy = RefrigData(RefrigNum)%HfValues(LoPressIndexLiq) &
                 + PressInterpRatioLiq*(RefrigData(RefrigNum)%HfValues(HiPressIndexLiq) &
				    - RefrigData(RefrigNum)%HfValues(LoPressIndexLiq))

  IF ((Enthalpy>=SatLiqEnthalpy).AND. (Enthalpy<=SatVapEnthalpy)) THEN
    Quality=(Enthalpy-SatLiqEnthalpy)/(SatVapEnthalpy-SatLiqEnthalpy)
	
	IF (Quality < 0.5) THEN
		PressInterpRatio=PressInterpRatioLiq
		LoPressIndex=LoPressIndexLiq
		HiPressIndex=HiPressIndexLiq
	ELSE
		PressInterpRatio=PressInterpRatioVap
		LoPressIndex=LoPressIndexVap
		HiPressIndex=HiPressIndexVap
	END IF

	SatTemperature=RefrigData(RefrigNum)%HTemps(LoPressIndex)+PressInterpRatio * &
				(RefrigData(RefrigNum)%HTemps(HiPressIndex)-RefrigData(RefrigNum)%HTemps(LoPressIndex)) 

	SELECT CASE (Property)
	CASE ('quality')
	PH=Quality
	CASE ('temperature')
	PH = SatTemperature
	CASE ('enthalpy')   
	!CALCULATE ENTHAPY OF SATURATED REGION
	PH = GetInterpolatedSatProp(SatTemperature, RefrigData(RefrigNum)%HTemps, &
                                             RefrigData(RefrigNum)%HfValues, &
                       RefrigData(RefrigNum)%HfgValues, Quality)
	CASE ('entropy')  
	!CALCULATE ENTROPY OF SATURATED REGION	 
	PH = GetInterpolatedSatProp(SatTemperature, RefrigData(RefrigNum)%STemps, &
                                             RefrigData(RefrigNum)%SfValues, &
                       RefrigData(RefrigNum)%SfgValues, Quality)
	CASE ('density') 	 
	!CALCULATE DENSITY OF SATURATED REGION
	PH = GetInterpolatedSatProp(SatTemperature, RefrigData(RefrigNum)%RhoTemps, &
                                             1/RefrigData(RefrigNum)%RhofValues, &
                                1/RefrigData(RefrigNum)%RhofgValues, Quality)
	PH=1/PH


    CASE ('viscosity') 	 
	!CALCULATE VISCOSITY OF SATURATED REGION
	IF (Quality > 0.5) THEN
    PH= RefrigData(RefrigNum)%DVfgValues(LoPressIndex) + PressInterpRatio * &
                              (RefrigData(RefrigNum)%DVfgValues(HiPressIndex) - &
                              RefrigData(RefrigNum)%DVfgValues(LoPressIndex))

	ELSEIF (Quality <= 0.5) THEN
	PH= RefrigData(RefrigNum)%DVfValues(LoPressIndex) + PressInterpRatio * &
                              (RefrigData(RefrigNum)%DVfValues(HiPressIndex) - &
                              RefrigData(RefrigNum)%DVfValues(LoPressIndex))
	ELSE
	PH=0
	END IF

	CASE ('conductivity') 	 
	!CALCULATE CONDUCTIVITY OF SATURATED REGION
	IF (Quality > 0.5) THEN
	PH= RefrigData(RefrigNum)%CfgValues(LoPressIndex) + PressInterpRatio * &
                              (RefrigData(RefrigNum)%CfgValues(HiPressIndex) - &
                              RefrigData(RefrigNum)%CfgValues(LoPressIndex))

	ELSEIF (Quality <= 0.5) THEN
    PH= RefrigData(RefrigNum)%CfValues(LoPressIndex) + PressInterpRatio * &
                              (RefrigData(RefrigNum)%CfValues(HiPressIndex) - &
                              RefrigData(RefrigNum)%CfValues(LoPressIndex))
	ELSE
	PH=0
	END IF

	CASE ('specificheat') 	 
	!CALCULATE SPECIFIC HEAT OF SATURATED REGION
   IF (Quality > 0.5) THEN
	PH= RefrigData(RefrigNum)%CpfgValues(LoPressIndex) + PressInterpRatio * &
                              (RefrigData(RefrigNum)%CpfgValues(HiPressIndex) - &
                              RefrigData(RefrigNum)%CpfgValues(LoPressIndex))
    
	ELSEIF (Quality <= 0.5) THEN
    PH= RefrigData(RefrigNum)%CpfValues(LoPressIndex) + PressInterpRatio * &
                              (RefrigData(RefrigNum)%CpfValues(HiPressIndex) - &
                              RefrigData(RefrigNum)%CpfValues(LoPressIndex))
	ELSE
	PH=0
	END IF
	END SELECT
  END IF

  IF (Enthalpy>SatVapEnthalpy) THEN
   !CALCULATE PROPERTY FOR SUPERHEATED REGION
    LoPressIndex = FindArrayIndex(Pressure, RefrigData(RefrigNum)%SHPress)
    HiPressIndex = LoPressIndex + 1

	! check temperature data range and attempt to cap if necessary
	IF((LoPressIndex > 0) .AND. (LoPressIndex < SIZE(RefrigData(RefrigNum)%SHPress)) )THEN ! in range
	HiPressIndex  = LoPressIndex + 1
	ELSE IF (LoPressIndex<1)THEN ! below lower bound
	LoPressIndex = 1
	HiPressIndex = LoPressIndex
	ELSE  ! out of range
    HiPressIndex = LoPressIndex
	END IF
  
  ! check for lowest non-zero value in lower press data
  DO Loop = 1, RefrigData(RefrigNum)%NumSuperTempPts
    IF (RefrigData(RefrigNum)%HshValues(Loop,LoPressIndex) > 0.0) THEN
      LoPressStart = Loop
      EXIT 
    END IF
  END DO

  ! check for lowest non-zero value in high press data
  DO Loop = 1, RefrigData(RefrigNum)%NumSuperTempPts
    IF (RefrigData(RefrigNum)%HshValues(Loop,HiPressIndex) > 0.0) THEN
      HiPressStart = Loop
      EXIT 
    END IF
  END DO

  ! find bounds of both hi and lo pressure data
  PressStart = MAX(LoPressStart, HiPressStart)
  PressFinish = RefrigData(RefrigNum)%NumSuperTempPts
  ! calculate interpolation ratio w.r.t temperature
  ! This ratio is used to find enthalpies at the given temperature
  PressInterpRatio = (Pressure - RefrigData(RefrigNum)%SHPress(LoPressIndex))/ &
                    (RefrigData(RefrigNum)%SHPress(HiPressIndex) - &
                     RefrigData(RefrigNum)%SHPress(LoPressIndex) )

  ! search for array index by bisection
  start = PressStart     ! set the bounds
  finish = PressFinish

  ! find the bounds of the enthalpy data available
  EnthalpyMax = MAX(RefrigData(RefrigNum)%HshValues(PressStart,LoPressIndex), &
                    RefrigData(RefrigNum)%HshValues(PressFinish,HiPressIndex))
  EnthalpyMin = MIN(RefrigData(RefrigNum)%HshValues(PressStart,LoPressIndex), &
                    RefrigData(RefrigNum)%HshValues(PressFinish,HiPressIndex))
   ! go ahead and search

    DO Loop = start, finish

     EnthalpyCheck = RefrigData(RefrigNum)%HshValues(Loop,LoPressIndex) + &
                     PressInterpRatio * (RefrigData(RefrigNum)%HshValues(Loop,HiPressIndex) - &
                     RefrigData(RefrigNum)%HshValues(Loop,LoPressIndex) )  
					  
	IF (Enthalpy < EnthalpyCheck) THEN
    start=loop
	EXIT 
    END IF
    END DO

    LoEnthalpyIndex  = start - 1
    HiEnthalpyIndex = start 

	IF (RefrigData(RefrigNum)%HshValues(LoEnthalpyIndex,HiPressIndex).EQ.0) THEN
	! SscValues(LoEnthalpyIndex,HiPressIndex) is in saturated region
	EnthalpyLow = RefrigData(RefrigNum)%HshValues(LoEnthalpyIndex,LoPressIndex) + &
                  PressInterpRatio * (PQ(Refrigerant,RefrigData(RefrigNum)%SHPress(HiPressIndex),1.0,'enthalpy',RefrigIndex,Error) - &
                  RefrigData(RefrigNum)%HshValues(LoEnthalpyIndex,LoPressIndex) )
	ELSE
    ! calculate enthalpies adjacent specified enthalpy at given temperature
    EnthalpyLow = RefrigData(RefrigNum)%HshValues(LoEnthalpyIndex,LoPressIndex) + &
                  PressInterpRatio * (RefrigData(RefrigNum)%HshValues(LoEnthalpyIndex,HiPressIndex) - &
                  RefrigData(RefrigNum)%HshValues(LoEnthalpyIndex,LoPressIndex) )
	END IF

    EnthalpyHigh =  RefrigData(RefrigNum)%HshValues(HiEnthalpyIndex,LoPressIndex) + &
                    PressInterpRatio * (RefrigData(RefrigNum)%HshValues(HiEnthalpyIndex,HiPressIndex) - &
                    RefrigData(RefrigNum)%HshValues(HiEnthalpyIndex,LoPressIndex) )
    ! calculate an interpolation ratio
    EnthInterpRatio = (Enthalpy - EnthalpyLow) / (EnthalpyHigh - EnthalpyLow)
    ! apply this interpolation ratio to find the final pressure
    IF (property.EQ.'quality') THEN
		!QUALITY FOR SUPERHEATED REGION =100
		PH = 100
	ELSE
		SELECT CASE (property)
		CASE ('temperature')
		LoEnthLoProperty=RefrigData(RefrigNum)%SHTemps(LoEnthalpyIndex)		
		IF (RefrigData(RefrigNum)%HshValues(LoEnthalpyIndex,HiPressIndex).EQ.0) THEN
		LoEnthHiProperty=0
		ELSE
		LoEnthHiProperty=LoEnthLoProperty
		END IF
		HiEnthLoProperty=RefrigData(RefrigNum)%SHTemps(HiEnthalpyIndex)
		HiEnthHiProperty=HiEnthLoProperty
		PropertyType = 1
		CASE ('entropy')
		!CALCULATE ENTROPY OF SUBCOOLED REGION
		LoEnthLoProperty = RefrigData(RefrigNum)%Sshvalues(LoEnthalpyIndex,LoPressIndex)
        LoEnthHiProperty = RefrigData(RefrigNum)%Sshvalues(LoEnthalpyIndex,HiPressIndex)
        HiEnthLoProperty = RefrigData(RefrigNum)%Sshvalues(HiEnthalpyIndex,LoPressIndex)
        HiEnthHiProperty = RefrigData(RefrigNum)%Sshvalues(HiEnthalpyIndex,HiPressIndex)
		PropertyType = 2
		CASE ('density')
		!CALCULATE DENSITY OF SUBCOOLED REGION
		LoEnthLoProperty = RefrigData(RefrigNum)%Rhoshvalues(LoEnthalpyIndex,LoPressIndex)
        LoEnthHiProperty = RefrigData(RefrigNum)%Rhoshvalues(LoEnthalpyIndex,HiPressIndex)
        HiEnthLoProperty = RefrigData(RefrigNum)%Rhoshvalues(HiEnthalpyIndex,LoPressIndex)
        HiEnthHiProperty = RefrigData(RefrigNum)%Rhoshvalues(HiEnthalpyIndex,HiPressIndex)
		PropertyType = 3
		CASE ('viscosity')
        !CALCULATE VISCOSITY OF SUBCOOLED REGION
		LoEnthLoProperty = RefrigData(RefrigNum)%DVshvalues(LoEnthalpyIndex,LoPressIndex)
        LoEnthHiProperty = RefrigData(RefrigNum)%DVshvalues(LoEnthalpyIndex,HiPressIndex)
        HiEnthLoProperty = RefrigData(RefrigNum)%DVshvalues(HiEnthalpyIndex,LoPressIndex)
        HiEnthHiProperty = RefrigData(RefrigNum)%DVshvalues(HiEnthalpyIndex,HiPressIndex)
		PropertyType = 4
		CASE ('specificheat')
        !CALCULATE SPECIFIC HEAT OF SUBCOOLED REGION
		LoEnthLoProperty = RefrigData(RefrigNum)%CPshvalues(LoEnthalpyIndex,LoPressIndex)
        LoEnthHiProperty = RefrigData(RefrigNum)%CPshvalues(LoEnthalpyIndex,HiPressIndex)
        HiEnthLoProperty = RefrigData(RefrigNum)%CPshvalues(HiEnthalpyIndex,LoPressIndex)
        HiEnthHiProperty = RefrigData(RefrigNum)%CPshvalues(HiEnthalpyIndex,HiPressIndex)
		PropertyType = 5
		CASE ('conductivity')
        !CALCULATE CONDUCTIVITY OF SUPERHEATED REGION
		LoEnthLoProperty = RefrigData(RefrigNum)%Cshvalues(LoEnthalpyIndex,LoPressIndex)
        LoEnthHiProperty = RefrigData(RefrigNum)%Cshvalues(LoEnthalpyIndex,HiPressIndex)
        HiEnthLoProperty = RefrigData(RefrigNum)%Cshvalues(HiEnthalpyIndex,LoPressIndex)
        HiEnthHiProperty = RefrigData(RefrigNum)%Cshvalues(HiEnthalpyIndex,HiPressIndex)		
		PropertyType = 6
		END SELECT
!Checks if the property at lower end is zero. Probably this should not be applicable to temperature!! Sankar
		IF (LoEnthHiProperty.EQ.0 .AND. PropertyType .NE. 1) THEN !Estimate the point at saturated region
		LoEnthHiProperty= PQ(Refrigerant,RefrigData(RefrigNum)%SHPress(HiPressIndex),1.0,Property,RefrigIndex,Error)
		END IF
		PropertyLow = LoEnthLoProperty + PressInterpRatio * (LoEnthHiProperty - LoEnthLoProperty )
		PropertyHigh = HiEnthLoProperty + PressInterpRatio * (HiEnthHiProperty - HiEnthLoProperty )
        ! interpolate w.r.t. temperature
        PH = PropertyLow + EnthInterpRatio * (PropertyHigh - PropertyLow)
    END IF
	END IF

	IF (Enthalpy<SatLiqEnthalpy) THEN
	!CALCULATE PROPERTY FOR SUBCOOLED REGION
    LoPressIndex = FindArrayIndexSC(Pressure, RefrigData(RefrigNum)%SCPress)
    HiPressIndex = LoPressIndex + 1

	! check temperature data range and attempt to cap if necessary
	IF((LoPressIndex > 0) .AND. (LoPressIndex < SIZE(RefrigData(RefrigNum)%SCPress)) )THEN ! in range
	HiPressIndex  = LoPressIndex + 1
	ELSE IF (LoPressIndex<1)THEN ! below lower bound
	LoPressIndex = 1
	HiPressIndex = LoPressIndex
	ELSE  ! out of range
    HiPressIndex = LoPressIndex
	END IF
  
  ! check for lowest non-zero value in lower press data
  DO Loop = 1, RefrigData(RefrigNum)%NumSubcoolTempPts
    IF (RefrigData(RefrigNum)%HscValues(Loop,LoPressIndex) > 0.0) THEN
      LoPressStart = Loop
      EXIT 
    END IF
  END DO

  ! check for lowest non-zero value in high press data
  DO Loop = 1, RefrigData(RefrigNum)%NumSubcoolTempPts
    IF (RefrigData(RefrigNum)%HshValues(Loop,HiPressIndex) > 0.0) THEN
      HiPressStart = Loop
      EXIT 
    END IF
  END DO

  ! find bounds of both hi and lo pressure data
  PressStart = MAX(LoPressStart, HiPressStart)
  PressFinish = RefrigData(RefrigNum)%NumSubcoolTempPts
  ! calculate interpolation ratio w.r.t temperature
  ! This ratio is used to find enthalpies at the given temperature
  PressInterpRatio = (Pressure - RefrigData(RefrigNum)%SCPress(LoPressIndex))/ &
                    (RefrigData(RefrigNum)%SCPress(HiPressIndex) - &
                     RefrigData(RefrigNum)%SCPress(LoPressIndex) )

  ! search for array index by bisection
  start = PressStart     ! set the bounds
  finish = PressFinish

  ! find the bounds of the enthalpy data available
  EnthalpyMax = MAX(RefrigData(RefrigNum)%HscValues(PressStart,LoPressIndex), &
                    RefrigData(RefrigNum)%HscValues(PressFinish,HiPressIndex))
  EnthalpyMin = MIN(RefrigData(RefrigNum)%HscValues(PressStart,LoPressIndex), &
                    RefrigData(RefrigNum)%HscValues(PressFinish,HiPressIndex))
   ! go ahead and search

    DO Loop = start, finish
     EnthalpyCheck = RefrigData(RefrigNum)%HscValues(Loop,LoPressIndex) + &
                     PressInterpRatio * (RefrigData(RefrigNum)%HscValues(Loop,HiPressIndex) - &
                     RefrigData(RefrigNum)%HscValues(Loop,LoPressIndex) )   
	IF (Enthalpy > EnthalpyCheck) THEN
    start=loop
	EXIT 
    END IF
    END DO

    LoEnthalpyIndex  = start - 1
    HiEnthalpyIndex = start 

    ! calculate enthalpies adjacent specified enthalpy at given temperature
	IF (RefrigData(RefrigNum)%HscValues(LoEnthalpyIndex,HiPressIndex).EQ.0) THEN
	! SscValues(LoEnthalpyIndex,HiPressIndex) is in saturated region
	EnthalpyLow = RefrigData(RefrigNum)%HscValues(LoEnthalpyIndex,LoPressIndex) + &
                  PressInterpRatio * (PQ(Refrigerant,RefrigData(RefrigNum)%SCPress(HiPressIndex),0.0,'enthalpy',RefrigIndex,Error) - &
                  RefrigData(RefrigNum)%HscValues(LoEnthalpyIndex,LoPressIndex) )
	ELSE
    EnthalpyLow = RefrigData(RefrigNum)%HscValues(LoEnthalpyIndex,LoPressIndex) + &
                  PressInterpRatio * (RefrigData(RefrigNum)%HscValues(LoEnthalpyIndex,HiPressIndex) - &
                  RefrigData(RefrigNum)%HscValues(LoEnthalpyIndex,LoPressIndex) )
	END IF
    EnthalpyHigh =  RefrigData(RefrigNum)%HscValues(HiEnthalpyIndex,LoPressIndex) + &
                    PressInterpRatio * (RefrigData(RefrigNum)%HscValues(HiEnthalpyIndex,HiPressIndex) - &
                    RefrigData(RefrigNum)%HscValues(HiEnthalpyIndex,LoPressIndex) )
    ! calculate an interpolation ratio
    EnthInterpRatio = (Enthalpy - EnthalpyLow) / (EnthalpyHigh - EnthalpyLow)
    ! apply this interpolation ratio to find the final pressure

	IF (property.EQ.'quality') THEN
		!QUALITY FOR SUBCOOLED REGION =-100
		PH = -100
	ELSE
		SELECT CASE (property)
		CASE ('temperature')
		LoEnthLoProperty=RefrigData(RefrigNum)%SCTemps(LoEnthalpyIndex)		
		IF (RefrigData(RefrigNum)%HscValues(LoEnthalpyIndex,HiPressIndex).EQ.0) THEN
		LoEnthHiProperty=0
		ELSE
		LoEnthHiProperty=LoEnthLoProperty
		END IF
		HiEnthLoProperty=RefrigData(RefrigNum)%SCTemps(HiEnthalpyIndex)
		HiEnthHiProperty=HiEnthLoProperty
		CASE ('entropy')
		!CALCULATE ENTROPY OF SUBCOOLED REGION
		LoEnthLoProperty = RefrigData(RefrigNum)%SscValues(LoEnthalpyIndex,LoPressIndex)
        LoEnthHiProperty = RefrigData(RefrigNum)%SscValues(LoEnthalpyIndex,HiPressIndex)
        HiEnthLoProperty = RefrigData(RefrigNum)%SscValues(HiEnthalpyIndex,LoPressIndex)
        HiEnthHiProperty = RefrigData(RefrigNum)%SscValues(HiEnthalpyIndex,HiPressIndex)
		CASE ('density')
		!CALCULATE DENSITY OF SUBCOOLED REGION
		LoEnthLoProperty = RefrigData(RefrigNum)%RhoscValues(LoEnthalpyIndex,LoPressIndex)
        LoEnthHiProperty = RefrigData(RefrigNum)%RhoscValues(LoEnthalpyIndex,HiPressIndex)
        HiEnthLoProperty = RefrigData(RefrigNum)%RhoscValues(HiEnthalpyIndex,LoPressIndex)
        HiEnthHiProperty = RefrigData(RefrigNum)%RhoscValues(HiEnthalpyIndex,HiPressIndex)
		CASE ('viscosity')
        !CALCULATE VISCOSITY OF SUBCOOLED REGION
		LoEnthLoProperty = RefrigData(RefrigNum)%DVscValues(LoEnthalpyIndex,LoPressIndex)
        LoEnthHiProperty = RefrigData(RefrigNum)%DVscValues(LoEnthalpyIndex,HiPressIndex)
        HiEnthLoProperty = RefrigData(RefrigNum)%DVscValues(HiEnthalpyIndex,LoPressIndex)
        HiEnthHiProperty = RefrigData(RefrigNum)%DVscValues(HiEnthalpyIndex,HiPressIndex)
		CASE ('specificheat')
        !CALCULATE SPECIFIC HEAT OF SUBCOOLED REGION
		LoEnthLoProperty = RefrigData(RefrigNum)%CPscValues(LoEnthalpyIndex,LoPressIndex)
        LoEnthHiProperty = RefrigData(RefrigNum)%CPscValues(LoEnthalpyIndex,HiPressIndex)
        HiEnthLoProperty = RefrigData(RefrigNum)%CPscValues(HiEnthalpyIndex,LoPressIndex)
        HiEnthHiProperty = RefrigData(RefrigNum)%CPscValues(HiEnthalpyIndex,HiPressIndex)
		CASE ('conductivity')
        !CALCULATE CONDUCTIVITY OF SUPERHEATED REGION
		LoEnthLoProperty = RefrigData(RefrigNum)%CscValues(LoEnthalpyIndex,LoPressIndex)
        LoEnthHiProperty = RefrigData(RefrigNum)%CscValues(LoEnthalpyIndex,HiPressIndex)
        HiEnthLoProperty = RefrigData(RefrigNum)%CscValues(HiEnthalpyIndex,LoPressIndex)
        HiEnthHiProperty = RefrigData(RefrigNum)%CscValues(HiEnthalpyIndex,HiPressIndex)		
		END SELECT

		IF (LoEnthHiProperty.EQ.0) THEN !Estimate the point at saturated region
		LoEnthHiProperty= PQ(Refrigerant,RefrigData(RefrigNum)%SCPress(HiPressIndex),0.0,Property,RefrigIndex,Error)
		END IF
		PropertyLow = LoEnthLoProperty + PressInterpRatio * (LoEnthHiProperty - LoEnthLoProperty )
		PropertyHigh = HiEnthLoProperty + PressInterpRatio * (HiEnthHiProperty - HiEnthLoProperty )
        ! interpolate w.r.t. temperature
        PH = PropertyLow + EnthInterpRatio * (PropertyHigh - PropertyLow)
	END IF

  END IF
  END IF

END FUNCTION PH

!*****************************************************************************
REAL FUNCTION PS(Refrigerant,Pressure,Entropy,Property,RefrigIndex,Error)

        ! SUBROUTINE INFORMATION:
        !       AUTHOR         Mike Turner
        !       DATE WRITTEN   10 December 99
        !       MODIFIED       Rick Strand (April 2000, May 2000)
        !                      Simon Rees (May 2002)
        !                      
        !       RE-ENGINEERED  Kenneth Tang (Nov 2003)
		!                      Ipseng Iu (Nov 2006)

        ! PURPOSE OF THIS FUNCTION:
        ! This finds enthalpy for given pressure and entropy.

        ! METHODOLOGY EMPLOYED:
        ! Determine the state of the properties using saturated liquid and vapor 
		! entropy and find the property requested.

        ! REFERENCES:
        ! na

        ! USE STATEMENTS:
        ! na

  IMPLICIT NONE           ! Enforce explicit typing of all variables in this routine

        ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN)  :: Refrigerant ! carries in substance name
  REAL,             INTENT(IN)  :: Pressure	   ! actual pressure given as input
  REAL,             INTENT(IN)  :: Entropy     ! actual entropy given as input
  INTEGER,       INTENT(INOUT)  :: RefrigIndex ! Index to Refrigerant Properties   !RS: Debugging: Removal of plethora of RefrigIndex definitions in the code
  CHARACTER(len=*), INTENT(IN)  :: Property
  INTEGER(2),      INTENT(OUT)  :: Error       !Error flag: 0-no; otherwise-yes

        ! INTERFACE BLOCK SPECIFICATIONS:
        ! na

        ! DERIVED TYPE DEFINITIONS:
        ! na

        ! FUNCTION LOCAL VARIABLE DECLARATIONS:
		
  INTEGER :: RefrigNum    ! index for refrigerant under consideration
  REAL    :: PressInterpRatio
  REAL    :: PressInterpRatioLiq
  REAL    :: PressInterpRatioVap
  INTEGER :: LoPressIndex
  INTEGER :: HiPressIndex
  INTEGER :: LoPressIndexLiq
  INTEGER :: HiPressIndexLiq
  INTEGER :: LoPressIndexVap
  INTEGER :: HiPressIndexVap
  INTEGER :: LoPressStart
  INTEGER :: HiPressStart
  INTEGER :: Start
  INTEGER :: Finish
  INTEGER :: PressStart
  INTEGER :: PressFinish
  INTEGER :: LoEntropyIndex
  INTEGER :: HiEntropyIndex
  REAL    :: EntInterpRatio
  REAL    :: EntropyCheck
  REAL    :: EntropyMax
  REAL    :: EntropyMin
  REAL    :: SatVapEntropy
  REAL    :: SatLiqEntropy
  REAL    :: EntropyLow
  REAL    :: EntropyHigh
  REAL    :: PropertyLow
  REAL    :: PropertyHigh
  REAL	  :: SatTemperature
  REAL    :: Quality
  REAL	  :: LoEntLoProperty
  REAL	  :: LoEntHiProperty
  REAL	  :: HiEntLoProperty
  REAL	  :: HiEntHiProperty

  INTEGER :: Loop
  LOGICAL :: ErrorFlag                  ! error flag for current call

  Error=0
  ErrorFlag = .FALSE.
    
  IF (GetInput) THEN
    CALL GetFluidPropertiesData
    GetInput = .FALSE.
  END IF

  IF (NumOfRefrigerants == 0) &
    CALL ShowFatalError('GetSatTemperatureRefrig: No refrigerants found--cannot evaluate saturation temperature')

  ErrorFlag = .False.

  IF (RefrigIndex > 0) THEN
    RefrigNum=RefrigIndex
  ELSE
    ! Find which refrigerant (index) is being requested
    RefrigNum = FindRefrigerant(Refrigerant)
    RefrigIndex=RefrigNum
  ENDIF

  ! get the array indices

  LoPressIndexLiq = FindArrayIndex(Pressure, RefrigData(RefrigNum)%PsfValues)
  HiPressIndexLiq = LoPressIndexLiq + 1

  LoPressIndexVap = FindArrayIndex(Pressure, RefrigData(RefrigNum)%PsgValues)
  HiPressIndexVap = LoPressIndexVap + 1

  ! check for out of data bounds problems

  IF (LoPressIndexLiq == 0) THEN
    LoPressIndexLiq = MAX(1, LoPressIndexLiq)
    ErrorFlag = .True.
	WRITE(*,*) 'Pressure out of range' !- Commented by ISI 01/19/04
	Error=1
  ELSE IF(HiPressIndexLiq > Size(RefrigData(RefrigNum)%PsfValues))THEN
    LoPressIndexLiq = MAX(1, LoPressIndexLiq)
    ErrorFlag = .True.
	WRITE(*,*) 'Pressure out of range' !- Commented by ISI 01/19/04
	Error=1
  END IF

  IF (LoPressIndexVap == 0) THEN
    LoPressIndexVap = MAX(1, LoPressIndexVap)
    ErrorFlag = .True.
	WRITE(*,*) 'Pressure out of range' !- Commented by ISI 01/19/04
	Error=1
  ELSE IF(HiPressIndexVap > Size(RefrigData(RefrigNum)%PsgValues))THEN
    LoPressIndexVap = MAX(1, LoPressIndexVap)
    ErrorFlag = .True.
	WRITE(*,*) 'Pressure out of range' !- Commented by ISI 01/19/04
	Error=1
  END IF

  IF(.NOT. ErrorFlag)THEN

  PressInterpRatioLiq = (Pressure - RefrigData(RefrigNum)%PsfValues(LoPressIndexLiq)) / &
                      (RefrigData(RefrigNum)%PsfValues(HiPressIndexLiq) &
                       - RefrigData(RefrigNum)%PsfValues(LoPressIndexLiq))

  PressInterpRatioVap = (Pressure - RefrigData(RefrigNum)%PsgValues(LoPressIndexVap)) / &
                      (RefrigData(RefrigNum)%PsgValues(HiPressIndexVap) &
                       - RefrigData(RefrigNum)%PsgValues(LoPressIndexVap))

  SatVapEntropy = RefrigData(RefrigNum)%SfgValues(LoPressIndexVap) &
                 + PressInterpRatioVap*(RefrigData(RefrigNum)%SfgValues(HiPressIndexVap) &
				    - RefrigData(RefrigNum)%SfgValues(LoPressIndexVap))

  SatLiqEntropy = RefrigData(RefrigNum)%SfValues(LoPressIndexLiq) &
                 + PressInterpRatioLiq*(RefrigData(RefrigNum)%SfValues(HiPressIndexLiq) &
				    - RefrigData(RefrigNum)%SfValues(LoPressIndexLiq))

!------SATURATED REGION-------- 
  IF ((Entropy>=SatLiqEntropy).AND. (Entropy<=SatVapEntropy)) THEN
    Quality=(Entropy-SatLiqEntropy)/(SatVapEntropy-SatLiqEntropy)

	IF (Quality < 0.5) THEN
		PressInterpRatio=PressInterpRatioLiq
		LoPressIndex=LoPressIndexLiq
		HiPressIndex=HiPressIndexLiq
	ELSE
		PressInterpRatio=PressInterpRatioVap
		LoPressIndex=LoPressIndexVap
		HiPressIndex=HiPressIndexVap
	END IF

	! apply final linear interpolation
	SatTemperature=RefrigData(RefrigNum)%STemps(LoPressIndex)+PressInterpRatio * &
				(RefrigData(RefrigNum)%STemps(HiPressIndex)-RefrigData(RefrigNum)%STemps(LoPressIndex)) 

	SELECT CASE (Property)
	CASE ('enthalpy')   
	!CALCULATE ENTHAPY OF SATURATED REGION
	PS = GetInterpolatedSatProp(SatTemperature, RefrigData(RefrigNum)%HTemps, &
                                             RefrigData(RefrigNum)%HfValues, &
                       RefrigData(RefrigNum)%HfgValues, Quality)

	END SELECT
 END IF
!------SUPERHEATED REGION--------
  IF (Entropy>SatVapEntropy) THEN
    LoPressIndex = FindArrayIndex(Pressure, RefrigData(RefrigNum)%SHPress)
    HiPressIndex = LoPressIndex + 1
	! check temperature data range and attempt to cap if necessary
	IF((LoPressIndex > 0) .AND. (LoPressIndex < SIZE(RefrigData(RefrigNum)%SHPress)) )THEN ! in range
	HiPressIndex  = LoPressIndex + 1
	ELSE IF (LoPressIndex<1)THEN ! below lower bound
	LoPressIndex = 1
	HiPressIndex = LoPressIndex
	ELSE  ! out of range
    HiPressIndex = LoPressIndex
	END IF
  
  ! check for lowest non-zero value in lower press data
  DO Loop = 1, RefrigData(RefrigNum)%NumSuperTempPts
    IF (RefrigData(RefrigNum)%SshValues(Loop,LoPressIndex) > 0.0) THEN
      LoPressStart = Loop
      EXIT 
    END IF
  END DO

  ! check for lowest non-zero value in high press data
  DO Loop = 1, RefrigData(RefrigNum)%NumSuperTempPts
    IF (RefrigData(RefrigNum)%SshValues(Loop,HiPressIndex) > 0.0) THEN
      HiPressStart = Loop
      EXIT 
    END IF
  END DO

  ! find bounds of both hi and lo pressure data
  PressStart = MAX(LoPressStart, HiPressStart)
  PressFinish = RefrigData(RefrigNum)%NumSuperTempPts
  ! calculate interpolation ratio w.r.t temperature
  ! This ratio is used to find enthalpies at the given temperature
  PressInterpRatio = (Pressure - RefrigData(RefrigNum)%SHPress(LoPressIndex))/ &
                    (RefrigData(RefrigNum)%SHPress(HiPressIndex) - &
                     RefrigData(RefrigNum)%SHPress(LoPressIndex) )

  ! search for array index by bisection
  start = PressStart     ! set the bounds
  finish = PressFinish

  ! find the bounds of the enthalpy data available
  EntropyMax = MAX(RefrigData(RefrigNum)%SshValues(PressStart,LoPressIndex), &
                    RefrigData(RefrigNum)%SshValues(PressFinish,HiPressIndex))
  EntropyMin = MIN(RefrigData(RefrigNum)%SshValues(PressStart,LoPressIndex), &
                    RefrigData(RefrigNum)%SshValues(PressFinish,HiPressIndex))
   ! go ahead and search

    DO Loop = start, finish

     EntropyCheck = RefrigData(RefrigNum)%SshValues(Loop,LoPressIndex) + &
                     PressInterpRatio * (RefrigData(RefrigNum)%SshValues(Loop,HiPressIndex) - &
                     RefrigData(RefrigNum)%SshValues(Loop,LoPressIndex) )  
					  
	IF (Entropy < EntropyCheck) THEN
    start=loop
	EXIT 
    END IF
    END DO

    LoEntropyIndex  = start - 1
    HiEntropyIndex = start 

	
    ! calculate enthalpies adjacent specified entropy at given temperature
	IF (RefrigData(RefrigNum)%SshValues(LoEntropyIndex,HiPressIndex).EQ.0) THEN
	!SshValues(LoEntropyIndex,HiPressIndex) is in saturated region
	EntropyLow = RefrigData(RefrigNum)%SshValues(LoEntropyIndex,LoPressIndex) + &
                  PressInterpRatio * (PQ(Refrigerant,RefrigData(RefrigNum)%SHPress(HiPressIndex),1.0,'entropy',RefrigIndex,Error) - &
                  RefrigData(RefrigNum)%SshValues(LoEntropyIndex,LoPressIndex) )
	ELSE
	 EntropyLow = RefrigData(RefrigNum)%SshValues(LoEntropyIndex,LoPressIndex) + &
                  PressInterpRatio * (RefrigData(RefrigNum)%SshValues(LoEntropyIndex,HiPressIndex) - &
                  RefrigData(RefrigNum)%SshValues(LoEntropyIndex,LoPressIndex) )
	END IF

    EntropyHigh =  RefrigData(RefrigNum)%SshValues(HiEntropyIndex,LoPressIndex) + &
                    PressInterpRatio * (RefrigData(RefrigNum)%SshValues(HiEntropyIndex,HiPressIndex) - &
                    RefrigData(RefrigNum)%SshValues(HiEntropyIndex,LoPressIndex) )
    ! calculate an interpolation ratio
    EntInterpRatio = (Entropy - EntropyLow) / (EntropyHigh - EntropyLow)

    ! apply this interpolation ratio to find the final pressure
    SELECT CASE (property)
		CASE ('enthalpy')
        !CALCULATE ENTHALPY OF SUPERHEATED REGION
		LoEntLoProperty = RefrigData(RefrigNum)%HshValues(LoEntropyIndex,LoPressIndex)
        LoEntHiProperty = RefrigData(RefrigNum)%HshValues(LoEntropyIndex,HiPressIndex)
        HiEntLoProperty = RefrigData(RefrigNum)%HshValues(HiEntropyIndex,LoPressIndex)
        HiEntHiProperty = RefrigData(RefrigNum)%HshValues(HiEntropyIndex,HiPressIndex)
		
		IF (LoEntHiProperty.EQ.0) THEN !Estimate the point at saturated region
		LoEntHiProperty = PQ(Refrigerant,RefrigData(RefrigNum)%SHPress(HiPressIndex),1.0,Property,RefrigIndex,Error)
		END IF
		PropertyLow = LoEntLoProperty + PressInterpRatio * (LoEntHiProperty - LoEntLoProperty )
		PropertyHigh = HiEntLoProperty + PressInterpRatio * (HiEntHiProperty - HiEntLoProperty )
        ! interpolate w.r.t. temperature
        PS = PropertyLow + EntInterpRatio * (PropertyHigh - PropertyLow)
	END SELECT
 END IF
!-------SUBCOOLED REGION-------!
 IF (Entropy<SatLiqEntropy) THEN
    LoPressIndex = FindArrayIndexSC(Pressure, RefrigData(RefrigNum)%SCPress)
    HiPressIndex = LoPressIndex + 1

	! check temperature data range and attempt to cap if necessary
	IF((LoPressIndex > 0) .AND. (LoPressIndex < SIZE(RefrigData(RefrigNum)%SCPress)) )THEN ! in range
	HiPressIndex  = LoPressIndex + 1
	ELSE IF (LoPressIndex<1)THEN ! below lower bound
	LoPressIndex = 1
	HiPressIndex = LoPressIndex
	ELSE  ! out of range
    HiPressIndex = LoPressIndex
	END IF
  
  ! check for lowest non-zero value in lower press data
  DO Loop = 1, RefrigData(RefrigNum)%NumSubcoolTempPts
    IF (RefrigData(RefrigNum)%SscValues(Loop,LoPressIndex) > 0.0) THEN
      LoPressStart = Loop
      EXIT 
    END IF
  END DO

  ! check for lowest non-zero value in high press data
  DO Loop = 1, RefrigData(RefrigNum)%NumSubcoolTempPts
    IF (RefrigData(RefrigNum)%SshValues(Loop,HiPressIndex) > 0.0) THEN
      HiPressStart = Loop
      EXIT 
    END IF
  END DO

  ! find bounds of both hi and lo pressure data
  PressStart = MAX(LoPressStart, HiPressStart)
  PressFinish = RefrigData(RefrigNum)%NumSubcoolTempPts
  ! calculate interpolation ratio w.r.t temperature
  ! This ratio is used to find enthalpies at the given temperature
  PressInterpRatio = (Pressure - RefrigData(RefrigNum)%SCPress(LoPressIndex))/ &
                    (RefrigData(RefrigNum)%SCPress(HiPressIndex) - &
                     RefrigData(RefrigNum)%SCPress(LoPressIndex) )

  ! search for array index by bisection
  start = PressStart     ! set the bounds
  finish = PressFinish

  ! find the bounds of the enthalpy data available
  EntropyMax = MAX(RefrigData(RefrigNum)%SscValues(PressStart,LoPressIndex), &
                    RefrigData(RefrigNum)%SscValues(PressFinish,HiPressIndex))
  EntropyMin = MIN(RefrigData(RefrigNum)%SscValues(PressStart,LoPressIndex), &
                    RefrigData(RefrigNum)%SscValues(PressFinish,HiPressIndex))
   ! go ahead and search

    DO Loop = start, finish
     EntropyCheck = RefrigData(RefrigNum)%SscValues(Loop,LoPressIndex) + &
                     PressInterpRatio * (RefrigData(RefrigNum)%SscValues(Loop,HiPressIndex) - &
                     RefrigData(RefrigNum)%SscValues(Loop,LoPressIndex) )   
	IF (Entropy > EntropyCheck) THEN
    start=loop
	EXIT 
    END IF
    END DO

    LoEntropyIndex  = start - 1
    HiEntropyIndex = start 

    ! calculate enthalpies adjacent specified Entropy at given temperature
	IF (RefrigData(RefrigNum)%SscValues(LoEntropyIndex,HiPressIndex).EQ.0) THEN
	!SscValues(LoEntropyIndex,HiPressIndex) is in saturated region
	EntropyLow = RefrigData(RefrigNum)%SscValues(LoEntropyIndex,LoPressIndex) + &
                  PressInterpRatio * (PQ(Refrigerant,RefrigData(RefrigNum)%SCPress(HiPressIndex),0.0,'entropy',RefrigIndex,Error) - &
                  RefrigData(RefrigNum)%SscValues(LoEntropyIndex,LoPressIndex) )
	ELSE
    EntropyLow = RefrigData(RefrigNum)%SscValues(LoEntropyIndex,LoPressIndex) + &
                  PressInterpRatio * (RefrigData(RefrigNum)%SscValues(LoEntropyIndex,HiPressIndex) - &
                  RefrigData(RefrigNum)%SscValues(LoEntropyIndex,LoPressIndex) )
	END IF

    EntropyHigh =  RefrigData(RefrigNum)%SscValues(HiEntropyIndex,LoPressIndex) + &
                    PressInterpRatio * (RefrigData(RefrigNum)%SscValues(HiEntropyIndex,HiPressIndex) - &
                    RefrigData(RefrigNum)%SscValues(HiEntropyIndex,LoPressIndex) )
    ! calculate an interpolation ratio
    EntInterpRatio = (Entropy - EntropyLow) / (EntropyHigh - EntropyLow)
    ! apply this interpolation ratio to find the final pressure
    SELECT CASE (property)
		CASE ('enthalpy')
        !CALCULATE ENTHALPY OF SUPERHEATED REGION
		LoEntLoProperty = RefrigData(RefrigNum)%HscValues(LoEntropyIndex,LoPressIndex)
        LoEntHiProperty = RefrigData(RefrigNum)%HscValues(LoEntropyIndex,HiPressIndex)
        HiEntLoProperty = RefrigData(RefrigNum)%HscValues(HiEntropyIndex,LoPressIndex)
        HiEntHiProperty = RefrigData(RefrigNum)%HscValues(HiEntropyIndex,HiPressIndex)
		
		IF (LoEntHiProperty.EQ.0) THEN !Estimate the point at saturated region
		LoEntHiProperty= PQ(Refrigerant,RefrigData(RefrigNum)%SCPress(HiPressIndex),0.0,Property,RefrigIndex,Error)
		END IF
		PropertyLow = LoEntLoProperty + PressInterpRatio * (LoEntHiProperty - LoEntLoProperty )
		PropertyHigh = HiEntLoProperty + PressInterpRatio * (HiEntHiProperty - HiEntLoProperty )
        ! interpolate w.r.t. temperature
        PS = PropertyLow + EntInterpRatio * (PropertyHigh - PropertyLow)
	END SELECT
	END IF
 END IF
END FUNCTION PS
!*****************************************************************************

INTEGER FUNCTION FindRefrigerant(Refrigerant)

          ! FUNCTION INFORMATION:
          !       AUTHOR         Rick Strand
          !       DATE WRITTEN   May 2000
          !       MODIFIED       Simon Rees (June 2002)
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This function simply determines the index of the refrigerant named
          ! in the input variable to this routine within the derived type.

          ! METHODOLOGY EMPLOYED:
          ! Just checks to see whether or not the refrigerant name coming in can
          ! be found in the refrigerant derived type.  If so, the function is set
          ! to the index within the derived type.  If the input has not been read
          ! yet for some reason, that must be done.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor, ONLY: FindItemInList, MakeUPPERCase

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: Refrigerant ! carries in substance name

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: Found   ! Indicator for found item

          ! FLOW:
          ! Make sure we have already read in the input
  IF (GetInput) THEN
    CALL GetFluidPropertiesData
    GetInput = .FALSE.
  END IF

          ! Check to see if this glycol shows up in the glycol data
  Found=FindItemInList(MakeUPPERCase(Refrigerant),RefrigData%Name,NumOfRefrigerants)

  IF (Found > 0) THEN
    FindRefrigerant = Found
  ELSE ! not found - errors handled in calling proceedure
    FindRefrigerant = 0
  ENDIF

  RETURN

END FUNCTION FindRefrigerant

!*****************************************************************************

INTEGER FUNCTION FindArrayIndex(Value,Array)

          ! FUNCTION INFORMATION:
          !       AUTHOR         Rick Strand
          !       DATE WRITTEN   May 2000
          !       MODIFIED       Simon Rees (May 2002)
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This generic function simply finds the points in an array between
          ! which a single value is found.  The returned value is the index of
          ! the low point.

          ! METHODOLOGY EMPLOYED:
          ! Straight interval halving. It is assumed that the values in the array
          ! appear in ascending order. If the value is below that in the supplied
          ! data array a zero index is returned. If the value is above that in the 
          ! supplied data array, the max index is returned. This allows some error
          ! checking in the calling routine.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
          ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  REAL, INTENT(IN)               :: Value   ! Value to be placed/found within the array of values
  REAL, INTENT(IN), DIMENSION(:) :: Array   ! Array of values in ascending order

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: start                   ! sets low index value
  INTEGER :: finish                  ! sets high index value
  INTEGER :: middle                  ! difference of finish & start

          ! FLOW:
  start  = 1
  finish = SIZE(Array)

  ! check bounds of data and set limiting values of the index
  IF(Value < Array(start)) THEN
    FindArrayIndex = 0
  ELSE IF(Value >Array(finish)) THEN
    FindArrayIndex = finish
  ELSE  ! start searching by bisection method
    DO WHILE ((finish - start) > 1)
      middle = (finish + start) / 2
      IF (Value > Array(middle)) THEN
        start = middle
      ELSE
        finish = middle
      END IF
    END DO
    FindArrayIndex = start
  END IF

  RETURN

END FUNCTION FindArrayIndex

!*****************************************************************************

INTEGER FUNCTION FindArrayIndexSC(Value,Array)

          ! FUNCTION INFORMATION:
          !       AUTHOR         Rick Strand
          !       DATE WRITTEN   May 2000
          !       MODIFIED       Simon Rees (May 2002)
          !       RE-ENGINEERED  Kenneth Tang (Dec 2003)

          ! PURPOSE OF THIS FUNCTION:
          ! This generic function simply finds the points in an array between
          ! which a single value is found.  The returned value is the index of
          ! the low point.

          ! METHODOLOGY EMPLOYED:
          ! Straight interval halving. Assumed array is in (descending) order for 
		  ! (subcooled) region. If the value is below that in the supplied
          ! data array a zero index is returned. If the value is above that in the 
          ! supplied data array, the max index is returned. This allows some error
          ! checking in the calling routine.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
          ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  REAL, INTENT(IN)               :: Value   ! Value to be placed/found within the array of values
  REAL, INTENT(IN), DIMENSION(:) :: Array   ! Array of values in ascending order

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: start                   ! sets high index value
  INTEGER :: finish                  ! sets low index value
  INTEGER :: middle                  ! difference of finish & start

          ! FLOW:
  start  = 1
  finish = SIZE(Array)

  ! check bounds of data and set limiting values of the index
  IF(Value > Array(start)) THEN
    FindArrayIndexSC = 0
  ELSE IF(Value <Array(finish)) THEN
    FindArrayIndexSC = finish
  ELSE  ! start searching by bisection method
    DO WHILE ((finish - start) > 1)
      middle = (finish + start) / 2
      IF (Value < Array(middle)) THEN
        start = middle
      ELSE
        finish = middle
      END IF
    END DO
    FindArrayIndexSC = start
  END IF

  RETURN

END FUNCTION FindArrayIndexSC

!*****************************************************************************

REAL FUNCTION GetInterpolatedSatProp(Temperature, PropTemps, LiqProp, VapProp, Quality)

          ! FUNCTION INFORMATION:
          !       AUTHOR         Simon Rees
          !       DATE WRITTEN   May 2002
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This generic function performs an interpolation on the supplied saturated
          ! liquid and vapor data to find the saturated property value at a given
          ! temperature and quality. This function is used by all the functions that
          ! get saturated property values.

          ! METHODOLOGY EMPLOYED:
          ! Index of arrays either side of given temperature is found using FindArrayIndex.
          ! Double linear interpolation is used to first find property values at the given
          ! quality bounding the required temperature. These values are interpolated in the
          ! temperature domain to find the final value.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
          ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  REAL, INTENT(IN)               :: Temperature   ! Saturation Temp.
  REAL, INTENT(IN), DIMENSION(:) :: PropTemps   ! Array of temperature at which props are available
  REAL, INTENT(IN), DIMENSION(:) :: LiqProp     ! Array of saturated liquid properties
  REAL, INTENT(IN), DIMENSION(:) :: VapProp     ! Array of saturatedv apour properties
  REAL, INTENT(IN)               :: Quality     ! Quality

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:

  INTEGER :: HiTempIndex              ! array index for temp above input temp
  INTEGER :: LoTempIndex              ! array index for temp below input temp
  REAL    :: LoSatProp                  ! Sat. prop. at lower temp & given quality
  REAL    :: HiSatProp                ! Sat. prop. at higher temp & given quality
  REAL    :: TempInterpRatio            ! ratio to interpolate in temperature domain
 ! error counters and dummy string
  LOGICAL :: ErrorFlag                  ! error flag for current call
  INTEGER,SAVE :: TempRangeErrCount=0   ! cumulative error counter
  CHARACTER(len=25) CErrCount

  ErrorFlag = .False.

  LoTempIndex = FindArrayIndex(Temperature, PropTemps)
  HiTempIndex = LoTempIndex + 1

  IF (LoTempIndex == 0) THEN
    LoTempIndex = MAX(1, LoTempIndex)
    GetInterpolatedSatProp = LiqProp(LoTempIndex) + &
               Quality*(VapProp(LoTempIndex) - LiqProp(LoTempIndex))
    ErrorFlag = .True.
    ! Temperature supplied is out of bounds--produce an error message...
  ELSE IF(HiTempIndex > Size(PropTemps))THEN
    LoTempIndex = MAX(1, LoTempIndex)
    GetInterpolatedSatProp = LiqProp(LoTempIndex) + &
               Quality*(VapProp(LoTempIndex) - LiqProp(LoTempIndex))
    ErrorFlag = .True.
  END IF

  IF(.NOT. ErrorFlag)THEN
    ! find adjacent property values at the given quality
    LoSatProp = LiqProp(LoTempIndex) + &
                 Quality*(VapProp(LoTempIndex) - LiqProp(LoTempIndex))

    HiSatProp = LiqProp(HiTempIndex) + &
                 Quality*(VapProp(HiTempIndex) - LiqProp(HiTempIndex))

    ! find interpolation ratio in temperature direction
    TempInterpRatio = (Temperature - PropTemps(LoTempIndex)) / &
                      (PropTemps(HiTempIndex) - PropTemps(LoTempIndex))

    ! apply final linear interpolation
    GetInterpolatedSatProp = LoSatProp + TempInterpRatio*(HiSatProp - LoSatProp)
  ELSE
      TempRangeErrCount = TempRangeErrCount + 1
     ! send warning
      IF (TempRangeErrCount <= 15) THEN
        CALL ShowSevereError('GetInterpolatedSatProp: Saturation temperature for interpolation is out of range '// &
                             'of data supplied: **')
        WRITE(CErrCount,*) Temperature 
        CErrCount=ADJUSTL(CErrCount)
        CALL ShowContinueError('Refrigerant temperature = '//TRIM(CErrCount))
        WRITE(CErrCount,*) GetInterpolatedSatProp 
        CErrCount=ADJUSTL(CErrCount)
        CALL ShowContinueError('Returned saturated property value = '//TRIM(CErrCount))
      ELSEIF (TempRangeErrCount>1 .AND. MOD(TempRangeErrCount,500) == 0) THEN
        WRITE(CErrCount,*) TempRangeErrCount
        CErrCount=ADJUSTL(CErrCount)
        CALL ShowSevereError('GetInterpolatedSatProp: Refrigerant temperature for interpolation out of range error '// &
                              'continues: saturated pressure error count = '//TRIM(CErrCount)//' **')
      ENDIF
  END IF
  RETURN

END FUNCTION GetInterpolatedSatProp

INTEGER FUNCTION CheckFluidPropertyName(NameToCheck)

          ! FUNCTION INFORMATION:
          !       AUTHOR         Linda K. Lawrie
          !       DATE WRITTEN   October 2002
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This function checks on an input fluid property to make sure it is valid.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor, ONLY: FindItemInList

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: NameToCheck  ! Name from input(?) to be checked against valid FluidPropertyNames

          ! FUNCTION PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER Found

  IF (GetInput) THEN
    CALL GetFluidPropertiesData
    GetInput = .FALSE.
  END IF

  ! Item must be either in Refrigerant or Glycol list

  Found=FindItemInList(NameToCheck,RefrigData%Name,NumOfRefrigerants)
  CheckFluidPropertyName=Found

  RETURN

END FUNCTION CheckFluidPropertyName

!     NOTICE
!
!     Copyright © 1996-2003 The Board of Trustees of the University of Illinois
!     and The Regents of the University of California through Ernest Orlando Lawrence
!     Berkeley National Laboratory.  All rights reserved.
!
!     Portions of the EnergyPlus software package have been developed and copyrighted
!     by other individuals, companies and institutions.  These portions have been
!     incorporated into the EnergyPlus software package under license.   For a complete
!     list of contributors, see "Notice" located in EnergyPlus.f90.
!
!     NOTICE: The U.S. Government is granted for itself and others acting on its
!     behalf a paid-up, nonexclusive, irrevocable, worldwide license in this data to
!     reproduce, prepare derivative works, and perform publicly and display publicly.
!     Beginning five (5) years after permission to assert copyright is granted,
!     subject to two possible five year renewals, the U.S. Government is granted for
!     itself and others acting on its behalf a paid-up, non-exclusive, irrevocable
!     worldwide license in this data to reproduce, prepare derivative works,
!     distribute copies to the public, perform publicly and display publicly, and to
!     permit others to do so.
!
!     TRADEMARKS: EnergyPlus is a trademark of the US Department of Energy.
!

END MODULE FluidProperties_HPSim
