MODULE SizingManager

  ! MODULE INFORMATION:
  !       AUTHOR         Fred Buhl
  !       DATE WRITTEN   December 2000
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS MODULE:
  ! This module contains the data and routines relating to managing the sizing
  ! simulations.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES: none

  ! OTHER NOTES: none


  ! USE STATEMENTS:
  USE DataPrecisionGlobals
  USE DataGlobals_HPSimIntegrated, ONLY: DoWeathSim, OutputFileZoneSizing, OutputFileSysSizing, DoZoneSizing, DoPlantSizing, DoDesDaySim, &
  DoOutputReporting, KindOfSim, ZoneSizingCalc, KickOffSizing, KSRunPeriodWeather, BeginEnvrnFlag, &
  EndEnvrnFlag, WarmUpFlag, DayOfSim, DayofSimChr, NumofDayinEnvrn, BeginDay, BeginDayFlag, EndDayFlag, &
  HourofDay, BeginHourFlag, DoSystemSizing, EndHourFlag, Timestep, NumofTimeStepinHour, BeginSimFlag, &
  BeginTimeStepFlag, PreviousHour, Endday, HourofDay, EndZoneSizingCalc, SysSizingCalc, DuringDay, &
  NumofZones, EndSysSizingCalc, OutputFileInits, BeginFullSimFlag
  USE HeatBalanceManager
  USE WeatherManager
  USE DataSizing
  USE DataHVACGlobals, ONLY: NumPrimaryAirSys
  USE DataZoneEquipment, ONLY: ZoneEquipConfig
  USE DataStringGlobals, ONLY: CharTab, CharComma, CharSpace
  USE DataInterfaces

  IMPLICIT NONE    ! Enforce explicit typing of all variables

  PRIVATE

  ! MODULE PARAMETER DEFINITIONS: none

  ! DERIVED TYPE DEFINITIONS: none
  TYPE ZoneListData
    CHARACTER(len=MaxNameLength) :: Name  = ' '
    INTEGER :: NumOfZones=0
    INTEGER,ALLOCATABLE,DIMENSION(:) :: Zones
  END TYPE

  ! INTERFACE BLOCK SPECIFICATIONS: none

  ! MODULE VARIABLE DECLARATIONS:
  INTEGER :: NumAirLoops=0

  ! SUBROUTINE SPECIFICATIONS FOR MODULE SimulationManager

  PUBLIC  ManageSizing
  PRIVATE GetOARequirements
  PRIVATE GetZoneSizingInput
  PRIVATE GetSystemSizingInput
  PRIVATE GetPlantSizingInput
  PRIVATE GetSizingParams
  PRIVATE SetupZoneSizing
  PRIVATE ReportZoneSizing
  PRIVATE ReportSysSizing

CONTAINS

  ! MODULE SUBROUTINES:

  SUBROUTINE  ManageSizing

    ! SUBROUTINE INFORMATION:
    !       AUTHOR         Fred Buhl
    !       DATE WRITTEN   December 2000
    !       MODIFIED       na
    !       RE-ENGINEERED  na

    ! PURPOSE OF THIS SUBROUTINE:
    ! This subroutine manages the sizing simulations (using design day condiions)
    ! for zones, central air systems, and central plants and zone heating and cooling

    ! METHODOLOGY EMPLOYED:
    ! Design day simulations are run with the zones supplied with "Ideal Loads",
    ! yielding zone design supply air flow rates and zone heating and cooling capacities.
    !
    ! Design day simulations are run again with central air systems supplied by
    ! purchased hot and cold water, yielding central heating and cooling capacities.

    ! REFERENCES:
    ! na

    ! USE STATEMENTS:
    USE InputProcessor, ONLY: GetNumRangeCheckErrorsFound
    USE ZoneEquipmentManager, ONLY: UpdateZoneSizing, ManageZoneEquipment
    USE SimAirServingZones, ONLY: ManageAirLoops, UpdateSysSizing
    USE DataEnvironment, ONLY : TotDesDays, OutDryBulbTemp, OutHumRat, OutBaroPress, CurEnvirNum, Month, DayOfMonth, EndMonthFlag
    USE OutputReportPredefined
    USE DataHeatBalance, ONLY: Zone
    USE General, ONLY: TrimSigDigits

    IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

    ! SUBROUTINE ARGUMENT DEFINITIONS:
    ! na

    ! SUBROUTINE PARAMETER DEFINITIONS:
    CHARACTER(len=*), PARAMETER :: RoutineName='ManageSizing: '

    ! INTERFACE BLOCK SPECIFICATIONS: none

    ! DERIVED TYPE DEFINITIONS: none

    ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    LOGICAL, SAVE :: Available=.false. ! an environment is available to process
    LOGICAL, SAVE :: ErrorsFound=.false.
    INTEGER,EXTERNAL :: GetNewUnitNumber
    LOGICAL          :: SimAir=.false.
    LOGICAL          :: SimZoneEquip=.false.
    INTEGER          :: TimeStepInDay=0 ! time step number
    INTEGER          :: LastMonth=0
    INTEGER          :: LastDayOfMonth=0
    INTEGER          :: CtrlZoneNum=0  ! controlled zone index
    INTEGER          :: ZoneNum=0      ! index into the Zone data array for the controlled zone
    REAL(r64)        :: TempAtPeak=0.0       ! Outside temperature at peak cooling/heating for reporting
    REAL(r64)        :: HumRatAtPeak=0.0     ! Outside humidity ratio at peak cooling/heating for reporting
    INTEGER          :: TimeStepAtPeak=0     ! time step number at heat or cool peak
    INTEGER          :: DDNum=0              ! Design Day index
    INTEGER          :: AirLoopNum=0         ! air loop index
    !  EXTERNAL            ReportZoneSizing
    !  EXTERNAL            ReportSysSizing
    CHARACTER(len=MaxNameLength) :: curName
    INTEGER :: NumSizingPeriodsPerformed
    INTEGER :: write_stat

    ! FLOW:

    OutputFileZoneSizing = 0
    OutputFileSysSizing = 0
    TimeStepInDay = 0
    SysSizingRunDone = .FALSE.
    ZoneSizingRunDone = .FALSE.
    curName='Unknown'
    CALL GetOARequirements      ! get the OA requirements object
    CALL GetZoneAirDistribution ! get zone air distribution objects
    CALL GetSizingParams        ! get the building level sizing paramets
    CALL GetZoneSizingInput     ! get the Zone Sizing input
    CALL GetSystemSizingInput   ! get the System Sizing input
    CALL GetPlantSizingInput    ! get the Plant Sizing input

    ! okay, check sizing inputs vs desires vs requirements
    IF (DoZoneSizing .or. DoSystemSizing) THEN
      IF ((NumSysSizInput > 0 .and. NumZoneSizingInput == 0) .or.   &
      (.not. DoZoneSizing .and. DoSystemSizing .and. NumSysSizInput > 0)) THEN
      CALL ShowSevereError(RoutineName//'Requested System Sizing but did not request Zone Sizing.')
      CALL ShowContinueError('System Sizing cannot be done without Zone Sizing')
      CALL ShowFatalError('Program terminates for preceding conditions.')
    ENDIF
  ENDIF

  IF ( (DoZoneSizing) .AND. (NumZoneSizingInput == 0) ) THEN
    CALL ShowWarningError(RoutineName//'For a zone sizing run, there must be at least 1 Sizing:Zone input object.'//  &
    ' SimulationControl Zone Sizing option ignored.')
  END IF

  IF ( (NumZoneSizingInput > 0) .AND. (DoZoneSizing.OR.DoSystemSizing.OR.DoPlantSizing) ) THEN

    IF (DoDesDaySim .OR. DoWeathSim) THEN
      DoOutputReporting = .FALSE.
    END IF
    DoOutputReporting=.false.
    ZoneSizingCalc = .TRUE.
    Available=.true.
    OutputFileZoneSizing = GetNewUnitNumber()
    IF (SizingFileColSep == CharComma) THEN
      OPEN (OutputFileZoneSizing,FILE='epluszsz.csv',Action='write',iostat=write_stat)
      IF (write_stat /= 0) THEN
        CALL ShowFatalError('ManageSizing: Could not open file "epluszsz.csv" for output (write).')
      ENDIF
    ELSEIF (SizingFileColSep == CharTab) THEN
      OPEN (OutputFileZoneSizing,FILE='epluszsz.tab',Action='write',iostat=write_stat)
      IF (write_stat /= 0) THEN
        CALL ShowFatalError('ManageSizing: Could not open file "epluszsz.tab" for output (write).')
      ENDIF
    ELSE
      OPEN (OutputFileZoneSizing,FILE='epluszsz.txt',Action='write',iostat=write_stat)
      IF (write_stat /= 0) THEN
        CALL ShowFatalError('ManageSizing: Could not open file "epluszsz.txt" for output (write).')
      ENDIF
    ENDIF

    !CALL ShowMessage('Beginning Zone Sizing Calculations') !RS: Debugging: Removing the message call

    CALL ResetEnvironmentCounter
    KickOffSizing=.true.
    CALL SetupZoneSizing(ErrorsFound)
    KickOffSizing=.false.

    CALL ResetEnvironmentCounter
    CurOverallSimDay=0
    NumSizingPeriodsPerformed=0
    DO WHILE (Available) ! loop over environments

      CALL GetNextEnvironment(Available,ErrorsFound) ! get an environment

      IF (.not. Available) EXIT
      IF (ErrorsFound) EXIT

      ! check that environment is one of the design days
      IF (KindOfSim == ksRunPeriodWeather) THEN
        CYCLE
      ENDIF

      NumSizingPeriodsPerformed=NumSizingPeriodsPerformed+1

      BeginEnvrnFlag = .TRUE.
      EndEnvrnFlag   = .FALSE.
      EndMonthFlag   = .FALSE.
      WarmupFlag     = .TRUE.
      DayOfSim       =  0
      DayOfSimChr    ='0'
      CurEnvirNumSimDay=1
      CurOverallSimDay=CurOverallSimDay+1
      DO WHILE ((DayOfSim.LT.NumOfDayInEnvrn).OR.(WarmupFlag))  ! Begin day loop ...

        DayOfSim     = DayOfSim + 1
        IF (.not. WarmupFlag .and. DayOfSim > 1) THEN
          CurEnvirNumSimDay=CurEnvirNumSimDay+1
        ENDIF

        WRITE(DayOfSimChr,*) DayOfSim
        DayOfSimChr=ADJUSTL(DayOfSimChr)
        BeginDayFlag = .TRUE.
        EndDayFlag   = .FALSE.

        IF (WarmupFlag) THEN
          CALL DisplayString('Warming up')
        ELSE ! (.NOT.WarmupFlag)
          IF (DayOfSim.EQ.1)   CALL DisplayString('Performing Zone Sizing Simulation')
          CALL UpdateZoneSizing(BeginDay)
        END IF

        DO HourOfDay = 1, 24      ! Begin hour loop ...

          BeginHourFlag = .TRUE.
          EndHourFlag   = .FALSE.

          DO TimeStep = 1, NumOfTimeStepInHour  ! Begin time step (TINC) loop ...

            BeginTimeStepFlag = .TRUE.

            ! Set the End__Flag variables to true if necessary.  Note that
            ! each flag builds on the previous level.  EndDayFlag cannot be
            ! .true. unless EndHourFlag is also .true., etc.  Note that the
            ! EndEnvrnFlag and the EndSimFlag cannot be set during warmup.
            ! Note also that BeginTimeStepFlag, EndTimeStepFlag, and the
            ! SubTimeStepFlags can/will be set/reset in the HVAC Manager.

            IF ((TimeStep.EQ.NumOfTimeStepInHour)) THEN
              EndHourFlag = .TRUE.
              IF (HourOfDay.EQ.24) THEN
                EndDayFlag = .TRUE.
                IF ((.NOT.WarmupFlag).AND.(DayOfSim.EQ.NumOfDayInEnvrn)) THEN
                  EndEnvrnFlag = .TRUE.
                END IF
              END IF
            END IF

            CALL ManageWeather

            IF (.not. WarmupFlag) THEN
              TimeStepInDay = (HourOfDay-1)*NumOfTimeStepInHour + TimeStep
              IF (HourOfDay == 1 .and. TimeStep == 1) THEN
                DesDayWeath(CurOverallSimDay)%DateString=TRIM(TrimSigDigits(Month))//'/'//TRIM(TrimSigDigits(DayOfMonth))
              ENDIF
              DesDayWeath(CurOverallSimDay)%Temp(TimeStepInDay) = OutDryBulbTemp
              DesDayWeath(CurOverallSimDay)%HumRat(TimeStepInDay) = OutHumRat
              DesDayWeath(CurOverallSimDay)%Press(TimeStepInDay) = OutBaroPress
            ENDIF


            CALL ManageHeatBalance

            !  After the first iteration of HeatBalance, all the "input" has been gotten
            IF (BeginSimFlag) THEN
              IF (GetNumRangeCheckErrorsFound() > 0) THEN
                CALL ShowFatalError('Out of "range" values found in input')
              ENDIF
            ENDIF

            BeginHourFlag  = .FALSE.
            BeginDayFlag   = .FALSE.
            BeginEnvrnFlag = .FALSE.
            BeginSimFlag   = .FALSE.

          END DO                              ! ... End time step (TINC) loop.

          PreviousHour=HourOfDay

        END DO                    ! ... End hour loop.

        IF (EndDayFlag) CALL UpdateZoneSizing(EndDay)

        IF (.not. WarmupFlag .and. (DayOfSim > 0) .and. (DayOfSim.LT.NumOfDayInEnvrn)) THEN
          CurOverallSimDay=CurOverallSimDay+1
        ENDIF

      END DO                      ! ... End day loop.

      LastMonth=Month
      LastDayOfMonth=DayOfMonth

    END DO  ! ... End environment loop

    IF (NumSizingPeriodsPerformed > 0) THEN
      CALL UpdateZoneSizing(EndZoneSizingCalc)
      ZoneSizingRunDone = .TRUE.
    ELSE
      CALL ShowSevereError(RoutineName//'No Sizing periods were performed for Zone Sizing.'//  &
      ' No Zone Sizing calculations saved.')
      ErrorsFound=.true.
    ENDIF

  END IF

  ZoneSizingCalc = .FALSE.
  DoOutputReporting = .FALSE.
  Month=LastMonth
  DayOfMonth=LastDayOfMonth

  IF ( (DoSystemSizing)  .AND. (NumSysSizInput == 0) .AND. (NumAirLoops > 0) ) THEN
    CALL ShowWarningError(RoutineName//'For a system sizing run, there must be at least 1 Sizing:System object input.'//  &
    ' SimulationControl System Sizing option ignored.')
  END IF

  IF ( (NumSysSizInput > 0) .AND. (DoSystemSizing.OR.DoPlantSizing) .and. .not. ErrorsFound) THEN

    CALL ShowMessage('Beginning System Sizing Calculations')

    SysSizingCalc = .TRUE.
    Available=.true.
    OutputFileSysSizing = GetNewUnitNumber()
    IF (SizingFileColSep == CharComma) THEN
      OPEN (OutputFileSysSizing,FILE='eplusssz.csv',Action='write',iostat=write_stat)
      IF (write_stat /= 0) THEN
        CALL ShowFatalError('ManageSizing: Could not open file "eplusssz.csv" for output (write).')
      ENDIF
    ELSEIF (SizingFileColSep == CharTab) THEN
      OPEN (OutputFileSysSizing,FILE='eplusssz.tab',Action='write',iostat=write_stat)
      IF (write_stat /= 0) THEN
        CALL ShowFatalError('ManageSizing: Could not open file "eplusssz.tab" for output (write).')
      ENDIF
    ELSE
      OPEN (OutputFileSysSizing,FILE='eplusssz.txt',Action='write',iostat=write_stat)
      IF (write_stat /= 0) THEN
        CALL ShowFatalError('ManageSizing: Could not open file "eplusssz.txt" for output (write).')
      ENDIF
    ENDIF
    SimAir = .TRUE.
    SimZoneEquip = .TRUE.

    CALL ManageZoneEquipment(.TRUE.,SimZoneEquip,SimAir)
    CALL ManageAirLoops(.TRUE.,SimAir,SimZoneEquip)
    IF (GetNumRangeCheckErrorsFound() > 0) THEN
      CALL ShowFatalError(RoutineName//'Out of "range" values found in input')
    ENDIF

    CALL ResetEnvironmentCounter
    CurEnvirNumSimDay=0
    CurOverallSimDay=0
    NumSizingPeriodsPerformed=0
    DO WHILE (Available) ! loop over environments

      CALL GetNextEnvironment(Available,ErrorsFound) ! get an environment

      ! check that environment is one of the design days
      IF (KindOfSim == ksRunPeriodWeather) THEN
        CYCLE
      ENDIF

      IF (.not. Available) EXIT
      IF (ErrorsFound) EXIT

      NumSizingPeriodsPerformed=NumSizingPeriodsPerformed+1

      BeginEnvrnFlag = .TRUE.
      EndEnvrnFlag   = .FALSE.
      WarmupFlag     = .FALSE.
      DayOfSim       =  0
      DayOfSimChr    =  '0'
      CurEnvirNumSimDay=1
      CurOverallSimDay=CurOverallSimDay+1

      DO WHILE ((DayOfSim.LT.NumOfDayInEnvrn).OR.(WarmupFlag))  ! Begin day loop ...

        DayOfSim     = DayOfSim + 1
        IF (.not. WarmupFlag .and. DayOfSim > 1) THEN
          CurEnvirNumSimDay=CurEnvirNumSimDay+1
        ENDIF
        WRITE(DayOfSimChr,*) DayOfSim
        DayOfSimChr=ADJUSTL(DayOfSimChr)
        BeginDayFlag = .TRUE.
        EndDayFlag   = .FALSE.

        IF (WarmupFlag) THEN
          CALL DisplayString('Warming up')
        ELSE ! (.NOT.WarmupFlag)
          IF (DayOfSim.EQ.1)   CALL DisplayString('Calculating System sizing')
          CALL UpdateSysSizing(BeginDay)
        END IF

        DO HourOfDay = 1, 24      ! Begin hour loop ...

          BeginHourFlag = .TRUE.
          EndHourFlag   = .FALSE.

          DO TimeStep = 1, NumOfTimeStepInHour  ! Begin time step (TINC) loop ...

            BeginTimeStepFlag = .TRUE.

            ! Set the End__Flag variables to true if necessary.  Note that
            ! each flag builds on the previous level.  EndDayFlag cannot be
            ! .true. unless EndHourFlag is also .true., etc.  Note that the
            ! EndEnvrnFlag and the EndSimFlag cannot be set during warmup.

            IF ((TimeStep.EQ.NumOfTimeStepInHour)) THEN
              EndHourFlag = .TRUE.
              IF (HourOfDay.EQ.24) THEN
                EndDayFlag = .TRUE.
                IF ((.NOT.WarmupFlag).AND.(DayOfSim.EQ.NumOfDayInEnvrn)) THEN
                  EndEnvrnFlag = .TRUE.
                END IF
              END IF
            END IF

            CALL ManageWeather

            CALL UpdateSysSizing(DuringDay)

            BeginHourFlag  = .FALSE.
            BeginDayFlag   = .FALSE.
            BeginEnvrnFlag = .FALSE.


          END DO                              ! ... End time step (TINC) loop.

          PreviousHour=HourOfDay

        END DO                    ! ... End hour loop.

        IF (EndDayFlag)    CALL UpdateSysSizing(EndDay)

        IF (.not. WarmupFlag .and. (DayOfSim > 0) .and. (DayOfSim.LT.NumOfDayInEnvrn)) THEN
          CurOverallSimDay=CurOverallSimDay+1
        ENDIF

      END DO                      ! ... End day loop.


    END DO  ! ... End environment loop

    IF (NumSizingPeriodsPerformed > 0) THEN
      CALL UpdateSysSizing(EndSysSizingCalc)
      SysSizingRunDone = .TRUE.
    ELSE
      CALL ShowSevereError(RoutineName//'No Sizing periods were performed for System Sizing.'//  &
      ' No System Sizing calculations saved.')
      ErrorsFound=.true.
    ENDIF
  END IF
  SysSizingCalc = .FALSE.

  ! report sizing results to eio file
  IF (ZoneSizingRunDone) THEN
    DO CtrlZoneNum = 1,NumOfZones
      IF (.not. ZoneEquipConfig(CtrlZoneNum)%IsControlled) CYCLE
      ZoneNum = FinalZoneSizing(CtrlZoneNum)%ActualZoneNum
      IF (FinalZoneSizing(CtrlZoneNum)%DesCoolVolFlow > 0.0) THEN
        TimeStepAtPeak = FinalZoneSizing(CtrlZoneNum)%TimeStepNumAtCoolMax
        DDNum = FinalZoneSizing(CtrlZoneNum)%CoolDDNum
        IF (DDNum > 0 .AND. TimeStepAtPeak > 0) THEN
          TempAtPeak = DesDayWeath(DDNum)%Temp(TimeStepAtPeak)
          HumRatAtPeak = DesDayWeath(DDNum)%HumRat(TimeStepAtPeak)
        ELSE
          TempAtPeak = 0.0
          HumRatAtPeak = 0.0
        END IF
        CALL ReportZoneSizing(FinalZoneSizing(CtrlZoneNum)%ZoneName, &
        'Cooling', &
        CalcFinalZoneSizing(CtrlZoneNum)%DesCoolLoad, &
        FinalZoneSizing(CtrlZoneNum)%DesCoolLoad, &
        CalcFinalZoneSizing(CtrlZoneNum)%DesCoolVolFlow, &
        FinalZoneSizing(CtrlZoneNum)%DesCoolVolFlow, &
        FinalZoneSizing(CtrlZoneNum)%CoolDesDay, &
        CoolPeakDateHrMin(CtrlZoneNum), &
        TempAtPeak, &
        HumRatAtPeak, &
        Zone(ZoneNum)%FloorArea, &
        Zone(ZoneNum)%TotOccupants, &
        FinalZoneSizing(CtrlZoneNum)%MinOA)
        curName = FinalZoneSizing(CtrlZoneNum)%ZoneName
        CALL PreDefTableEntry(pdchZnClCalcDesLd,curName,CalcFinalZoneSizing(CtrlZoneNum)%DesCoolLoad)
        CALL PreDefTableEntry(pdchZnClUserDesLd,curName,FinalZoneSizing(CtrlZoneNum)%DesCoolLoad)
        CALL PreDefTableEntry(pdchZnClCalcDesAirFlow,curName,CalcFinalZoneSizing(CtrlZoneNum)%DesCoolVolFlow,3)
        CALL PreDefTableEntry(pdchZnClUserDesAirFlow,curName,FinalZoneSizing(CtrlZoneNum)%DesCoolVolFlow,3)
        CALL PreDefTableEntry(pdchZnClDesDay,curName,FinalZoneSizing(CtrlZoneNum)%CoolDesDay)
        CALL PreDefTableEntry(pdchZnClPkTime,curName,CoolPeakDateHrMin(CtrlZoneNum))
        CALL PreDefTableEntry(pdchZnClPkTemp,curName,TempAtPeak)
        CALL PreDefTableEntry(pdchZnClPkHum,curName,HumRatAtPeak,5)
      END IF
      IF (FinalZoneSizing(CtrlZoneNum)%DesHeatVolFlow > 0.0) THEN
        TimeStepAtPeak = FinalZoneSizing(CtrlZoneNum)%TimeStepNumAtHeatMax
        DDNum = FinalZoneSizing(CtrlZoneNum)%HeatDDNum
        IF (DDNum > 0 .AND. TimeStepAtPeak > 0) THEN
          TempAtPeak = DesDayWeath(DDNum)%Temp(TimeStepAtPeak)
          HumRatAtPeak = DesDayWeath(DDNum)%HumRat(TimeStepAtPeak)
        ELSE
          TempAtPeak = 0.0
          HumRatAtPeak = 0.0
        END IF
        CALL ReportZoneSizing(FinalZoneSizing(CtrlZoneNum)%ZoneName, &
        'Heating', &
        CalcFinalZoneSizing(CtrlZoneNum)%DesHeatLoad, &
        FinalZoneSizing(CtrlZoneNum)%DesHeatLoad, &
        CalcFinalZoneSizing(CtrlZoneNum)%DesHeatVolFlow, &
        FinalZoneSizing(CtrlZoneNum)%DesHeatVolFlow, &
        FinalZoneSizing(CtrlZoneNum)%HeatDesDay, &
        HeatPeakDateHrMin(CtrlZoneNum), &
        TempAtPeak, &
        HumRatAtPeak, &
        Zone(ZoneNum)%FloorArea, &
        Zone(ZoneNum)%TotOccupants, &
        FinalZoneSizing(CtrlZoneNum)%MinOA)
        curName = FinalZoneSizing(CtrlZoneNum)%ZoneName
        CALL PreDefTableEntry(pdchZnHtCalcDesLd,curName,CalcFinalZoneSizing(CtrlZoneNum)%DesHeatLoad)
        CALL PreDefTableEntry(pdchZnHtUserDesLd,curName,FinalZoneSizing(CtrlZoneNum)%DesHeatLoad)
        CALL PreDefTableEntry(pdchZnHtCalcDesAirFlow,curName,CalcFinalZoneSizing(CtrlZoneNum)%DesHeatVolFlow,3)
        CALL PreDefTableEntry(pdchZnHtUserDesAirFlow,curName,FinalZoneSizing(CtrlZoneNum)%DesHeatVolFlow,3)
        CALL PreDefTableEntry(pdchZnHtDesDay,curName,FinalZoneSizing(CtrlZoneNum)%HeatDesDay)
        CALL PreDefTableEntry(pdchZnHtPkTime,curName,HeatPeakDateHrMin(CtrlZoneNum))
        CALL PreDefTableEntry(pdchZnHtPkTemp,curName,TempAtPeak)
        CALL PreDefTableEntry(pdchZnHtPkHum,curName,HumRatAtPeak,5)
      END IF
    END DO
    ! Deallocate arrays no longer needed
    DEALLOCATE(ZoneSizing)
    DEALLOCATE(CalcZoneSizing)
  END IF
  IF (SysSizingRunDone) THEN
    DO AirLoopNum=1,NumPrimaryAirSys
      curName = FinalSysSizing(AirLoopNum)%AirPriLoopName
      CALL ReportSysSizing(curName, &
      'Calculated Cooling Design Air Flow Rate [m3/s]', &
      CalcSysSizing(AirLoopNum)%DesCoolVolFlow)
      CALL PreDefTableEntry(pdchSysSizCalcClAir,curName,CalcSysSizing(AirLoopNum)%DesCoolVolFlow)
      IF (ABS(CalcSysSizing(AirLoopNum)%DesCoolVolFlow) <= 1.E-8) THEN
        CALL ShowWarningError(RoutineName//'Calculated Cooling Design Air Flow Rate for System='//  &
        TRIM(FinalSysSizing(AirLoopNum)%AirPriLoopName)//' is zero.')
        CALL ShowContinueError('Check Sizing:Zone and ZoneControl:Thermostat inputs.')
      ENDIF
      CALL ReportSysSizing(curName, &
      'User Cooling Design Air Flow Rate [m3/s]', &
      FinalSysSizing(AirLoopNum)%DesCoolVolFlow)
      CALL PreDefTableEntry(pdchSysSizUserClAir,curName,FinalSysSizing(AirLoopNum)%DesCoolVolFlow)
      CALL ReportSysSizing(curName, &
      'Calculated Heating Design Air Flow Rate [m3/s]', &
      CalcSysSizing(AirLoopNum)%DesHeatVolFlow)
      CALL PreDefTableEntry(pdchSysSizCalcHtAir,curName,CalcSysSizing(AirLoopNum)%DesHeatVolFlow)
      IF (ABS(CalcSysSizing(AirLoopNum)%DesHeatVolFlow) <= 1.E-8) THEN
        CALL ShowWarningError(RoutineName//'Calculated Heating Design Air Flow Rate for System='//  &
        TRIM(FinalSysSizing(AirLoopNum)%AirPriLoopName)//' is zero.')
        CALL ShowContinueError('Check Sizing:Zone and ZoneControl:Thermostat inputs.')
      ENDIF
      CALL ReportSysSizing(curName, &
      'User Heating Design Air Flow Rate [m3/s]', &
      FinalSysSizing(AirLoopNum)%DesHeatVolFlow)
      CALL PreDefTableEntry(pdchSysSizUserHtAir,curName,FinalSysSizing(AirLoopNum)%DesHeatVolFlow)
    END DO
    ! Deallocate arrays no longer needed
    DEALLOCATE(SysSizing)
    DEALLOCATE(CalcSysSizing)
  END IF

  IF ( (DoPlantSizing)  .AND. (NumPltSizInput == 0)) THEN
    CALL ShowWarningError(RoutineName//'For a plant sizing run, there must be at least 1 Sizing:Plant object input.'//  &
    ' SimulationControl Plant Sizing option ignored.')
  END IF

  IF ( (NumPltSizInput > 0) .AND. (DoPlantSizing) .and. .not. ErrorsFound) THEN

    CALL ShowMessage('Beginning Plant Sizing Calculations')

  END IF

  IF (ErrorsFound) THEN
    CALL ShowFatalError('Program terminates due to preceding conditions.')
  ENDIF

  RETURN
END SUBROUTINE  ManageSizing

SUBROUTINE GetOARequirements

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         R. Raustad - FSEC
  !       DATE WRITTEN   February 2010
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! Obtains input data for the OA Requirements object and stores it in
  ! appropriate data structure.

  ! METHODOLOGY EMPLOYED:
  ! Uses InputProcessor "Get" routines to obtain data.
  ! This object requires only a name where the default values are assumed
  ! if subsequent fields are not entered.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE InputProcessor, ONLY: GetNumObjectsFound, GetObjectDefMaxArgs, GetObjectItem, VerifyName, SameString
  USE ScheduleManager, ONLY: GetScheduleIndex, CheckScheduleValueMinMax, GetScheduleMaxValue
  USE DataIPShortCuts
  USE General, ONLY: RoundSigDigits

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER :: RoutineName='GetOARequirements: ' ! include trailing blank space


  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: NumAlphas        ! Number of Alphas for each GetObjectItem call
  INTEGER :: NumNumbers       ! Number of Numbers for each GetObjectItem call
  INTEGER :: TotalArgs        ! Total number of alpha and numeric arguments (max) for a
  INTEGER :: IOStatus         ! Used in GetObjectItem
  INTEGER :: OAIndex
  LOGICAL :: ErrorsFound = .false.   ! If errors detected in input
  LOGICAL :: IsNotOK                 ! Flag to verify name
  LOGICAL :: IsBlank                 ! Flag for blank name
  !  REAL(r64) :: CalcAmt

  CHARACTER(Len=MaxNameLength)  :: CurrentModuleObject      ! for ease in getting objects
  CHARACTER(Len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: Alphas         ! Alpha input items for object
  CHARACTER(Len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: cAlphaFields   ! Alpha field names
  CHARACTER(Len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: cNumericFields ! Numeric field names
  REAL(r64), ALLOCATABLE, DIMENSION(:) :: Numbers           ! Numeric input items for object
  LOGICAL, ALLOCATABLE, DIMENSION(:)   :: lAlphaBlanks      ! Logical array, alpha field input BLANK = .true.
  LOGICAL, ALLOCATABLE, DIMENSION(:)   :: lNumericBlanks    ! Logical array, numeric field input BLANK = .true.

  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopfully this works.

  OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

  CurrentModuleObject='DesignSpecification:OutdoorAir'
  NumOARequirements = GetNumObjectsFound(TRIM(CurrentModuleObject))
  CALL GetObjectDefMaxArgs(CurrentModuleObject,TotalArgs,NumAlphas,NumNumbers)

  ALLOCATE(Alphas(NumAlphas))
  Alphas=' '
  ALLOCATE(cAlphaFields(NumAlphas))
  cAlphaFields=' '
  ALLOCATE(cNumericFields(NumNumbers))
  cNumericFields=' '
  ALLOCATE(Numbers(NumNumbers))
  Numbers=0.0
  ALLOCATE(lAlphaBlanks(NumAlphas))
  lAlphaBlanks=.true.
  ALLOCATE(lNumericBlanks(NumNumbers))
  lNumericBlanks=.true.

  IF (NumOARequirements .GT. 0) THEN
    ALLOCATE(OARequirements(NumOARequirements))

    !Start Loading the System Input
    DO OAIndex = 1,  NumOARequirements

      CALL GetObjectItem(TRIM(CurrentModuleObject),OAIndex,Alphas,NumAlphas,Numbers,NumNumbers,IOStatus,  &
      AlphaBlank=lAlphaBlanks,NumBlank=lNumericBlanks,  &
      AlphaFieldNames=cAlphaFields,NumericFieldNames=cNumericFields)

      CALL VerifyName(Alphas(1),OARequirements%Name,OAIndex-1,IsNotOK,IsBlank,TRIM(CurrentModuleObject)//' Name')
      IF (IsNotOK) THEN
        ErrorsFound=.true.
        IF (IsBlank) Alphas(1)='xxxxx'
      ENDIF

      OARequirements(OAIndex)%Name     = Alphas(1)

      IF(NumAlphas .GT. 1)THEN
        IF (SameString(Alphas(2) , 'Flow/Person')) THEN
          OARequirements(OAIndex)%OAFlowMethod = OAFlowPPer
        ELSEIF (SameString(Alphas(2), 'Flow/Zone')) THEN
          OARequirements(OAIndex)%OAFlowMethod = OAFlow
        ELSEIF (SameString(Alphas(2), 'Flow/Area')) THEN
          OARequirements(OAIndex)%OAFlowMethod = OAFlowPerArea
        ELSEIF (SameString(Alphas(2), 'AirChanges/Hour')) THEN
          OARequirements(OAIndex)%OAFlowMethod = OAFlowACH
        ELSEIF (SameString(Alphas(2), 'Sum')) THEN
          OARequirements(OAIndex)%OAFlowMethod = OAFlowSum
        ELSEIF (SameString(Alphas(2), 'Maximum')) THEN
          OARequirements(OAIndex)%OAFlowMethod = OAFlowMax
        ELSE
          CAll ShowSevereError(RoutineName//TRIM(CurrentModuleObject)//'="'//TRIM(OARequirements(OAIndex)%Name)//'",')
          CALL ShowContinueError('...Invalid '//TRIM(cAlphaFields(2))//'="'//TRIM(Alphas(2))//'",')
          CALL ShowContinueError('...Valid choices are Flow/Person, Flow/Zone, Flow/Area, AirChanges/Hour, Sum, Maximum.')
          ErrorsFound=.true.
        ENDIF
      ELSE
        ! default value for Outdoor Air Method
        OARequirements(OAIndex)%OAFlowMethod = OAFlowPPer
      END IF
      IF(NumNumbers .GT. 0)THEN
        OARequirements(OAIndex)%OAFlowPerPerson = Numbers(1)
      ELSE
        ! default value for Outdoor Air Flow per Person
        OARequirements(OAIndex)%OAFlowPerPerson = 0.00944D0
      END IF
      ! remaining fields default to 0
      IF(NumNumbers .GT. 1)THEN
        OARequirements(OAIndex)%OAFlowPerArea   = Numbers(2)
      END IF
      IF(NumNumbers .GT. 2)THEN
        OARequirements(OAIndex)%OAFlowPerZone   = Numbers(3)
      END IF
      IF(NumNumbers .GT. 3)THEN
        OARequirements(OAIndex)%OAFlowACH   = Numbers(4)
      END IF
      IF(NumAlphas .GT. 2)THEN
        IF(.NOT. lAlphaBlanks(3))THEN
          OARequirements(OAIndex)%OAFlowFracSchPtr = GetScheduleIndex(Alphas(3))
          IF(OARequirements(OAIndex)%OAFlowFracSchPtr .GT. 0) THEN
            IF (.NOT.CheckScheduleValueMinMax(OARequirements(OAIndex)%OAFlowFracSchPtr,'>=',0.0D0,'<=',1.0D0)) THEN
              CAll ShowSevereError(RoutineName//TRIM(CurrentModuleObject)//'="'//TRIM(OARequirements(OAIndex)%Name)//'",')
              CALL ShowContinueError('Error found in '//TRIM(cAlphaFields(3))//' = '//TRIM(Alphas(3)) )
              CALL ShowContinueError('Schedule values must be (>=0., <=1.)')
              ErrorsFound=.true.
            ELSE
              OARequirements(OAIndex)%MaxOAFractionSchValue = GetScheduleMaxValue(OARequirements(OAIndex)%OAFlowFracSchPtr)
            END IF
          ELSE
            !CAll ShowSevereError(RoutineName//TRIM(CurrentModuleObject)//'="'//TRIM(OARequirements(OAIndex)%Name)//'",')
            !CALL ShowContinueError('...Not Found '//TRIM(cAlphaFields(3))//'="'//TRIM(Alphas(3))//'".') !RS: Secret Search String
            WRITE(DebugFile,*) RoutineName//TRIM(CurrentModuleObject)//'="'//TRIM(OARequirements(OAIndex)%Name)//'",'
            WRITE(DebugFile,*) '...Not Found '//TRIM(cAlphaFields(3))//'="'//TRIM(Alphas(3))//'".'
            ErrorsFound=.TRUE.
          END IF
        END IF
      END IF

    END DO

    DEALLOCATE(Alphas)
    DEALLOCATE(cAlphaFields)
    DEALLOCATE(cNumericFields)
    DEALLOCATE(Numbers)
    DEALLOCATE(lAlphaBlanks)
    DEALLOCATE(lNumericBlanks)

    IF (ErrorsFound) THEN
      !CALL ShowFatalError(RoutineName//'Errors found in input.  Preceding condition(s) cause termination.')    !RS: Secret Search String
      WRITE(DebugFile,*) RoutineName//'Errors found in input. Preceding conditions(s) would like to cause termination.'
    ENDIF

  ENDIF

  RETURN

END SUBROUTINE GetOARequirements

SUBROUTINE GetZoneAirDistribution

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         T. Hong - LBNL
  !       DATE WRITTEN   March 2012
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! Obtains input data for the zone air distribution objects and stores it in
  ! appropriate data structure.

  ! METHODOLOGY EMPLOYED:
  ! Uses InputProcessor "Get" routines to obtain data.
  ! This object requires only a name where the default values are assumed
  ! if subsequent fields are not entered.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE InputProcessor, ONLY: GetNumObjectsFound, GetObjectDefMaxArgs, GetObjectItem, VerifyName, SameString
  USE ScheduleManager, ONLY: GetScheduleIndex, CheckScheduleValueMinMax
  USE DataIPShortCuts

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER :: RoutineName='GetZoneAirDistribution: ' ! include trailing blank space


  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: NumAlphas        ! Number of Alphas for each GetObjectItem call
  INTEGER :: NumNumbers       ! Number of Numbers for each GetObjectItem call
  INTEGER :: TotalArgs        ! Total number of alpha and numeric arguments (max) for a
  INTEGER :: IOStatus         ! Used in GetObjectItem
  INTEGER :: ZADIndex
  LOGICAL :: ErrorsFound = .false.   ! If errors detected in input
  LOGICAL :: IsNotOK                 ! Flag to verify name
  LOGICAL :: IsBlank                 ! Flag for blank name

  CHARACTER(Len=MaxNameLength)  :: CurrentModuleObject      ! for ease in getting objects
  CHARACTER(Len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: Alphas         ! Alpha input items for object
  CHARACTER(Len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: cAlphaFields   ! Alpha field names
  CHARACTER(Len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: cNumericFields ! Numeric field names
  REAL(r64), ALLOCATABLE, DIMENSION(:) :: Numbers           ! Numeric input items for object
  LOGICAL, ALLOCATABLE, DIMENSION(:)   :: lAlphaBlanks      ! Logical array, alpha field input BLANK = .true.
  LOGICAL, ALLOCATABLE, DIMENSION(:)   :: lNumericBlanks    ! Logical array, numeric field input BLANK = .true.

  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopfully this works.

  OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

  CurrentModuleObject='DesignSpecification:ZoneAirDistribution'
  NumZoneAirDistribution = GetNumObjectsFound(TRIM(CurrentModuleObject))
  CALL GetObjectDefMaxArgs(CurrentModuleObject,TotalArgs,NumAlphas,NumNumbers)

  ALLOCATE(Alphas(NumAlphas))
  Alphas=' '
  ALLOCATE(cAlphaFields(NumAlphas))
  cAlphaFields=' '
  ALLOCATE(cNumericFields(NumNumbers))
  cNumericFields=' '
  ALLOCATE(Numbers(NumNumbers))
  Numbers=0.0
  ALLOCATE(lAlphaBlanks(NumAlphas))
  lAlphaBlanks=.true.
  ALLOCATE(lNumericBlanks(NumNumbers))
  lNumericBlanks=.true.

  IF (NumZoneAirDistribution .GT. 0) THEN
    ALLOCATE(ZoneAirDistribution(NumZoneAirDistribution))

    !Start Loading the zone air distribution input
    DO ZADIndex = 1,  NumZoneAirDistribution

      CALL GetObjectItem(TRIM(CurrentModuleObject),ZADIndex,Alphas,NumAlphas,Numbers,NumNumbers,IOStatus,  &
      AlphaBlank=lAlphaBlanks,NumBlank=lNumericBlanks,  &
      AlphaFieldNames=cAlphaFields,NumericFieldNames=cNumericFields)

      CALL VerifyName(Alphas(1),ZoneAirDistribution%Name,ZADIndex-1,IsNotOK,IsBlank,TRIM(CurrentModuleObject)//' Name')
      IF (IsNotOK) THEN
        ErrorsFound=.true.
        IF (IsBlank) Alphas(1)='xxxxx'
      ENDIF

      ZoneAirDistribution(ZADIndex)%Name = Alphas(1)

      ! Zone Air Distribution Effectiveness in Cooling Mode
      IF(NumNumbers .GT. 0)THEN
        ZoneAirDistribution(ZADIndex)%ZoneADEffCooling = Numbers(1)
      ELSE
        ! default value
        ZoneAirDistribution(ZADIndex)%ZoneADEffCooling = 1.0
      END IF

      ! Zone Air Distribution Effectiveness in Heating Mode
      IF(NumNumbers .GT. 1)THEN
        ZoneAirDistribution(ZADIndex)%ZoneADEffHeating = Numbers(2)
      ELSE
        ! default value
        ZoneAirDistribution(ZADIndex)%ZoneADEffHeating = 1.0
      END IF

      ! Zone Secondary Recirculation Fraction
      IF(NumNumbers .GT. 2)THEN
        ZoneAirDistribution(ZADIndex)%ZoneSecondaryRecirculation = Numbers(3)
      ELSE
        ! default value
        ZoneAirDistribution(ZADIndex)%ZoneSecondaryRecirculation = 0.0
      END IF

      IF(NumAlphas .GT. 1)THEN
        IF(.NOT. lAlphaBlanks(2))THEN
          ZoneAirDistribution(ZADIndex)%ZoneADEffSchName = Alphas(2)
          ZoneAirDistribution(ZADIndex)%ZoneADEffSchPtr = GetScheduleIndex(Alphas(2))
          IF(ZoneAirDistribution(ZADIndex)%ZoneADEffSchPtr .GT. 0) THEN
            IF (.NOT.CheckScheduleValueMinMax(ZoneAirDistribution(ZADIndex)%ZoneADEffSchPtr,'>',0.0D0)) THEN
              CAll ShowSevereError(RoutineName//TRIM(CurrentModuleObject)//'="'//TRIM(ZoneAirDistribution(ZADIndex)%Name)//'",')
              CALL ShowContinueError('Error found in '//TRIM(cAlphaFields(2))//' = '//TRIM(Alphas(2)) )
              CALL ShowContinueError('Schedule values must be >0.0)')
              ErrorsFound=.true.
            END IF
          ELSE
            !CAll ShowSevereError(RoutineName//TRIM(CurrentModuleObject)//'="'//TRIM(ZoneAirDistribution(ZADIndex)%Name)//'",')
            !CALL ShowContinueError('...Not Found '//TRIM(cAlphaFields(2))//'="'//TRIM(Alphas(2))//'".') !RS: Secret String Search
            WRITE(DebugFile,*) RoutineName//TRIM(CurrentModuleObject)//'="'//TRIM(ZoneAirDistribution(ZADIndex)%Name)//'",'
            WRITE(DebugFile,*) '...Not Found '//TRIM(cAlphaFields(2))//'="'//TRIM(Alphas(2))//'".'
            ErrorsFound=.TRUE.
          END IF
        END IF
      END IF

    END DO

    DEALLOCATE(Alphas)
    DEALLOCATE(cAlphaFields)
    DEALLOCATE(cNumericFields)
    DEALLOCATE(Numbers)
    DEALLOCATE(lAlphaBlanks)
    DEALLOCATE(lNumericBlanks)

    IF (ErrorsFound) THEN
      !CALL ShowFatalError(RoutineName//'Errors found in input.  Preceding condition(s) cause termination.') !RS: Secret String Search
      WRITE(DebugFile,*) RoutineName//'Errors found in input. Preceding condition(s) would like to cause termination.'
    ENDIF

  ENDIF

  RETURN

END SUBROUTINE GetZoneAirDistribution

SUBROUTINE GetSizingParams

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Fred Buhl
  !       DATE WRITTEN   January 2002
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! Obtains input data for the Sizing Parameters object and stores it in
  ! appropriate data structure.

  ! METHODOLOGY EMPLOYED:
  ! Uses InputProcessor "Get" routines to obtain data.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE InputProcessor, ONLY: GetNumObjectsFound, GetObjectItem, VerifyName
  USE DataIPShortCuts
  USE General, ONLY: RoundSigDigits

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: NumAlphas        ! Number of Alphas for each GetObjectItem call
  INTEGER :: NumNumbers       ! Number of Numbers for each GetObjectItem call
  INTEGER :: IOStatus         ! Used in GetObjectItem
  INTEGER :: NumSizParams
  INTEGER :: Temp

  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopfully this works.

  OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

  cCurrentModuleObject='Sizing:Parameters'
  NumSizParams = GetNumObjectsFound(TRIM(cCurrentModuleObject))

  IF (NumSizParams == 1) THEN
    CALL GetObjectItem(TRIM(cCurrentModuleObject),1,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,IOStatus,  &
    AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
    AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
    IF (lNumericFieldBlanks(1) .or. rNumericArgs(1) < 0.0) THEN
      GlobalHeatSizingFactor = 1.0
    ELSE
      GlobalHeatSizingFactor = rNumericArgs(1)
    ENDIF
    IF (lNumericFieldBlanks(2) .or. rNumericArgs(2) < 0.0) THEN
      GlobalCoolSizingFactor = 1.0
    ELSE
      GlobalCoolSizingFactor = rNumericArgs(2)
    ENDIF
    IF (lNumericFieldBlanks(3) .or. rNumericArgs(3) <= 0.0) THEN
      NumTimeStepsInAvg=NumOfTimeStepInHour
    ELSE
      NumTimeStepsInAvg = INT(rNumericArgs(3))
    ENDIF
  ELSE IF (NumSizParams == 0) THEN
    GlobalHeatSizingFactor = 1.0
    GlobalCoolSizingFactor = 1.0
    NumTimeStepsInAvg = NumOfTimeStepInHour
  ELSE
    CALL ShowFatalError(TRIM(cCurrentModuleObject)//': More than 1 occurence of this object; only 1 allowed')
  END IF

  IF (NumTimeStepsInAvg < NumOfTimeStepInHour) THEN
    !CALL ShowWarningError(TRIM(cCurrentModuleObject)//': note '//TRIM(cNumericFieldNames(3))//' entered value=['//  &
    !   TRIM(RoundSigDigits(NumTimeStepsInAvg))//'] is less than 1 hour (i.e., '//trim(RoundSigDigits(NumOfTimeStepInHour))//  &
    !   ' timesteps).')  !RS: Secret Search String
    WRITE(DebugFile,*) TRIM(cCurrentModuleObject)//': note '//TRIM(cNumericFieldNames(3))//'entered value=['//TRIM(RoundSigDigits(NumTimeStepsInAvg)) &
    //'] is less than 1 hour (i.e., '//TRIM(RoundSigDigits(NumOfTimeStepInHour))//' timesteps).'
  ENDIF

  cCurrentModuleObject='OutputControl:Sizing:Style'
  Temp = GetNumObjectsFound(TRIM(cCurrentModuleObject))

  IF (Temp == 0) THEN
    cAlphaArgs(1)='Comma'
    SizingFileColSep = CharComma !comma
  ELSEIF (Temp == 1) THEN
    CALL GetObjectItem(TRIM(cCurrentModuleObject),1,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,IOStatus,  &
    AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
    AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
    IF (cAlphaArgs(1) == 'COMMA') THEN
      SizingFileColSep = CharComma !comma
      cAlphaArgs(1)='Comma'
    ELSEIF (cAlphaArgs(1) == 'TAB') THEN
      SizingFileColSep = CharTab  !tab
      cAlphaArgs(1)='Tab'
    ELSEIF (cAlphaArgs(1) == 'FIXED' .or. cAlphaArgs(1) == 'SPACE') THEN
      SizingFileColSep = CharSpace ! space
      cAlphaArgs(1)='Space'
    ELSE
      SizingFileColSep = CharComma !comma
      CALL ShowWarningError(TRIM(cCurrentModuleObject)//': invalid '//TRIM(cAlphaFieldNames(1))//' entered value="'//  &
      TRIM(cAlphaArgs(1))//'", Commas will be used to separate fields.')
      cAlphaArgs(1)='Comma'
    ENDIF
    Write(OutputFileInits,'(A)') '! <Sizing Output Files>,Style'
    WRITE(OutputFileInits,"('Sizing Output Files,',A)") TRIM(cAlphaArgs(1))
  ENDIF

  RETURN

END SUBROUTINE GetSizingParams

SUBROUTINE GetZoneSizingInput

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Fred Buhl
  !       DATE WRITTEN   December 2000
  !       MODIFIED       Mangesh Basarkar, 06/2011: Specifying zone outside air based on design specification object
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! Obtains input data for zone sizing objects and stores it in
  ! appropriate data structures.

  ! METHODOLOGY EMPLOYED:
  ! Uses InputProcessor "Get" routines to obtain data.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE InputProcessor, ONLY: GetNumObjectsFound, GetObjectItem, VerifyName, FindItemInList
  USE DataIPShortCuts
  USE General, ONLY: RoundSigDigits

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  TYPE GlobalMiscObject
    CHARACTER(len=MaxNameLength) :: Name  =' '
    INTEGER :: ZoneOrZoneListPtr          =0
    INTEGER :: NumOfZones                 =0
    INTEGER :: StartPtr                   =0
    LOGICAL :: ZoneListActive             =.false.
  END TYPE GlobalMiscObject

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: ZoneSizIndex ! loop index
  INTEGER :: NumAlphas        ! Number of Alphas for each GetObjectItem call
  INTEGER :: NumNumbers       ! Number of Numbers for each GetObjectItem call
  INTEGER :: IOStatus         ! Used in GetObjectItem
  LOGICAL :: ErrorsFound=.false.  ! Set to true if errors in input, fatal at end of routine
  LOGICAL :: IsNotOK                            ! Flag to verify name
  LOGICAL :: IsBlank                            ! Flag for blank name
  INTEGER :: NumDesDays       ! Number of design days in input
  INTEGER :: NumSizingZoneStatements
  INTEGER :: Item
  INTEGER :: Item1
  INTEGER :: ZLItem
  LOGICAL :: ErrFlag
  CHARACTER(len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: ZoneNames
  INTEGER :: NumZones
  TYPE(ZoneListData), ALLOCATABLE, DIMENSION(:) :: ZoneListNames
  INTEGER :: NumZoneLists
  TYPE (GlobalMiscObject), ALLOCATABLE, DIMENSION(:) :: SizingZoneObjects
  INTEGER OAIndex ! Index of design specification object
  INTEGER ObjIndex ! Index of zone air distribution effectiveness object name

  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopfully this works.

  OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

  cCurrentModuleObject='Sizing:Zone'
  NumSizingZoneStatements=GetNumObjectsFound(TRIM(cCurrentModuleObject))
  ALLOCATE(SizingZoneObjects(NumSizingZoneStatements))

  IF (NumSizingZoneStatements > 0) THEN
    Errflag=.false.
    CALL GetZoneAndZoneListNames(ErrFlag,NumZones,ZoneNames,NumZoneLists,ZoneListNames)
  ENDIF

  cCurrentModuleObject='Sizing:Zone'
  NumZoneSizingInput=0
  ErrFlag=.false.
  DO Item=1,NumSizingZoneStatements
    CALL GetObjectItem(TRIM(cCurrentModuleObject),Item,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,IOStatus,  &
    AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
    AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
    IsNotOK = .FALSE.
    IsBlank = .FALSE.
    CALL VerifyName(cAlphaArgs(1),SizingZoneObjects%Name,Item-1,IsNotOK,IsBlank,TRIM(cCurrentModuleObject)//' Name')
    IF (IsNotOK) THEN
      ErrorsFound = .TRUE.
      ErrFlag=.true.
      IF (IsBlank) cAlphaArgs(1) = 'xxxxx'
    END IF
    SizingZoneObjects(Item)%Name = cAlphaArgs(1)

    Item1=FindItemInList(cAlphaArgs(1),ZoneNames,NumZones)
    ZLItem=0
    IF (Item1 == 0 .and. NumZoneLists > 0) &
    ZLItem=FindItemInList(cAlphaArgs(1),ZoneListNames%Name,NumZoneLists)
    IF (Item1 > 0) THEN
      SizingZoneObjects(Item)%StartPtr=NumZoneSizingInput+1
      NumZoneSizingInput=NumZoneSizingInput+1
      SizingZoneObjects(Item)%NumOfZones=1
      SizingZoneObjects(Item)%ZoneListActive=.false.
      SizingZoneObjects(Item)%ZoneOrZoneListPtr=Item1
    ELSEIF (ZLItem > 0) THEN
      SizingZoneObjects(Item)%StartPtr=NumZoneSizingInput+1
      NumZoneSizingInput=NumZoneSizingInput+ZoneListNames(ZLItem)%NumOfZones
      SizingZoneObjects(Item)%NumOfZones=ZoneListNames(ZLItem)%NumOfZones
      SizingZoneObjects(Item)%ZoneListActive=.true.
      SizingZoneObjects(Item)%ZoneOrZoneListPtr=ZLItem
    ELSE
      !CALL ShowSevereError(trim(cCurrentModuleObject)//'="'//trim(cAlphaArgs(1))//'" invalid '//  &
      !     trim(cAlphaFieldNames(1))//' not found.')    !RS: Secret Search String
      WRITE(DebugFile,*) TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//'" invalid '// &
      TRIM(cAlphaFieldNames(1))//' not found.'
      ErrorsFound=.true.
      ErrFlag=.true.
    ENDIF
  ENDDO

  IF (ErrFlag) THEN
    !CALL ShowSevereError('GetZoneSizingInput: Errors with invalid names in '//trim(cCurrentModuleObject)//  &
    !   ' objects.')
    !CALL ShowContinueError('...These will not be read in.  Other errors may occur.')    !RS: Secret Search String
    WRITE(DebugFile,*) 'Get ZoneSizingInput: Errors with invalid names in '//TRIM(cCurrentModuleObject)//' objects.'
    WRITE(DebugFile,*) '...These will not be read in. Other Errors may occur.'
    NumZoneSizingInput=0
  ENDIF

  IF (NumZoneSizingInput > 0) THEN
    NumDesDays = GetNumObjectsFound('SizingPeriod:DesignDay') + GetNumObjectsFound('SizingPeriod:WeatherFileDays') +   &
    GetNumObjectsFound('SizingPeriod:WeatherFileConditionType')
    IF (NumDesDays == 0 .AND. (DoZoneSizing .OR. DoSystemSizing .OR. DoPlantSizing) ) THEN
      CALL ShowSevereError('Zone Sizing calculations need SizingPeriod:* input.')
      ErrorsFound = .TRUE.
    END IF
    ALLOCATE(ZoneSizingInput(NumZoneSizingInput))

    ZoneSizIndex=0
    DO Item = 1, NumSizingZoneStatements

      CALL GetObjectItem(TRIM(cCurrentModuleObject),Item,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,IOStatus,  &
      AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
      AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

      DO Item1=1,SizingZoneObjects(Item)%NumOfZones
        ZoneSizIndex=ZoneSizIndex+1
        IF (.not. SizingZoneObjects(Item)%ZoneListActive) THEN
          IF (SizingZoneObjects(Item)%ZoneOrZoneListPtr > 0) THEN
            ZoneSizingInput(ZoneSizIndex)%ZoneName = ZoneNames(SizingZoneObjects(Item)%ZoneOrZoneListPtr)
          ELSE
            ! Invalid zone, will be caught later
            ZoneSizingInput(ZoneSizIndex)%ZoneName = 'Invalid Zone Name'
          ENDIF
        ELSE  ! Zone list active
          IF (SizingZoneObjects(Item)%ZoneOrZoneListPtr > 0 .and.   &
          ZoneListNames(SizingZoneObjects(Item)%ZoneOrZoneListPtr)%Zones(Item1) > 0) THEN
          ZoneSizingInput(ZoneSizIndex)%ZoneName =   &
          ZoneNames(ZoneListNames(SizingZoneObjects(Item)%ZoneOrZoneListPtr)%Zones(Item1))
        ELSE
          ! Invalid zone, will be caught later
          ZoneSizingInput(ZoneSizIndex)%ZoneName = 'Invalid Zone Name'
        ENDIF
      ENDIF
      IsNotOK=.FALSE.
      IsBlank=.FALSE.
      CALL VerifyName(ZoneSizingInput(ZoneSizIndex)%ZoneName,ZoneSizingInput%ZoneName,  &
      ZoneSizIndex-1,IsNotOK,IsBlank,TRIM(cCurrentModuleObject)//' Name')
      IF (IsNotOK) THEN
        ErrorsFound=.true.
        IF (IsBlank) cAlphaArgs(1)='xxxxx'
      ENDIF
      IF (IsNotOK .and. .not. SizingZoneObjects(Item)%ZoneListActive) THEN
        CALL ShowContinueError('Zone may have been entered in a ZoneList assignment.')
      ENDIF

      ZoneSizingInput(ZoneSizIndex)%DesignSpecOAObjName = cAlphaArgs(2)

      ZoneSizingInput(ZoneSizIndex)%CoolDesTemp = rNumericArgs(1)
      ZoneSizingInput(ZoneSizIndex)%HeatDesTemp = rNumericArgs(2)
      !  N3, \field Zone Cooling Design Supply Air Humidity Ratio
      !      \required-field
      !      \minimum 0.0
      !      \type real
      !      \units kg-H2O/kg-air
      IF (lNumericFieldBlanks(3)) THEN
        ZoneSizingInput(ZoneSizIndex)%CoolDesHumRat = 0.0
      ELSEIF (rNumericArgs(3) < 0.0d0) THEN
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(3))//': '//  &
        TRIM(RoundSigDigits(rNumericArgs(3),2)))
        CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
        ErrorsFound=.true.
      ELSE
        ZoneSizingInput(ZoneSizIndex)%CoolDesHumRat = rNumericArgs(3)
      ENDIF
      !  N4, \field Zone Heating Design Supply Air Humidity Ratio
      !      \required-field
      !      \minimum 0.0
      !      \type real
      !      \units kg-H2O/kg-air
      IF (lNumericFieldBlanks(4)) THEN
        ZoneSizingInput(ZoneSizIndex)%HeatDesHumRat = 0.0
      ELSEIF (rNumericArgs(4) < 0.0d0) THEN
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(4))//': '//  &
        TRIM(RoundSigDigits(rNumericArgs(4),2)))
        CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
        ErrorsFound=.true.
      ELSE
        ZoneSizingInput(ZoneSizIndex)%HeatDesHumRat = rNumericArgs(4)
      ENDIF
      !  A3, \field Design Specification Object Name
      ! Getting zone OA parameters from Design Specification object
      IF (.NOT. lAlphaFieldBlanks(2)) THEN
        OAIndex=FindItemInList(ZoneSizingInput(ZoneSizIndex)%DesignSpecOAObjName,   &
        OARequirements%Name,numOARequirements)
        IF (OAIndex > 0) THEN
          ZoneSizingInput(ZoneSizIndex)%OADesMethod = OARequirements(OAIndex)%OAFlowMethod
          ZoneSizingInput(ZoneSizIndex)%DesOAFlowPPer = OARequirements(OAIndex)%OAFlowPerPerson
          ZoneSizingInput(ZoneSizIndex)%DesOAFlowPerArea = OARequirements(OAIndex)%OAFlowPerArea
          ZoneSizingInput(ZoneSizIndex)%DesOAFlow = OARequirements(OAIndex)%OAFlowPerZone
          ZoneSizingInput(ZoneSizIndex)%ZoneDesignSpecOAIndex = OAIndex
        ELSE
          CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cAlphaFieldNames(2))//'="'//  &
          TRIM(cAlphaArgs(2))//'".')
          ErrorsFound=.true.
        ENDIF
      ELSE ! If no design spec object specified, i.e. no OA, then set OA method to Flow/Person as default but flows to 0
        ZoneSizingInput(ZoneSizIndex)%OADesMethod = 1
        ZoneSizingInput(ZoneSizIndex)%DesOAFlowPPer = 0
        ZoneSizingInput(ZoneSizIndex)%DesOAFlowPerArea = 0
        ZoneSizingInput(ZoneSizIndex)%DesOAFlow = 0
      ENDIF

      !  N5, \field Zone Heating Sizing Factor
      !      \note if blank, global heating sizing factor from Sizing:Parameters is used.
      !      \minimum> 0
      IF (lNumericFieldBlanks(5) .or. rNumericArgs(5) == 0.0d0) THEN
        ZoneSizingInput(ZoneSizIndex)%HeatSizingFactor = GlobalHeatSizingFactor
      ELSEIF (rNumericArgs(5) < 0.0d0) THEN
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(5))//': '//  &
        TRIM(RoundSigDigits(rNumericArgs(5),2)))
        CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
        ErrorsFound=.true.
      ELSE
        ZoneSizingInput(ZoneSizIndex)%HeatSizingFactor = rNumericArgs(5)
      ENDIF
      !  N6, \field Zone Cooling Sizing Factor
      !      \note if blank, global cooling sizing factor from Sizing:Parameters is used.
      !      \minimum> 0
      IF (lNumericFieldBlanks(6) .or. rNumericArgs(6) == 0.0d0) THEN
        ZoneSizingInput(ZoneSizIndex)%CoolSizingFactor = GlobalCoolSizingFactor
      ELSEIF (rNumericArgs(6) < 0.0d0) THEN
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(6))//': '//  &
        TRIM(RoundSigDigits(rNumericArgs(6),2)))
        CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
        ErrorsFound=.true.
      ELSE
        ZoneSizingInput(ZoneSizIndex)%CoolSizingFactor = rNumericArgs(6)
      ENDIF
      !  N7, \field Cooling Design Air Flow Rate
      !      \type real
      !      \units m3/s
      !      \minimum 0
      !      \default 0
      !      \note This input is used if Cooling Design Air Flow Method is Flow/Zone
      !      \note This value will be multiplied by the global or zone sizing factor and
      !      \note by zone multipliers.
      IF (lNumericFieldBlanks(7)) THEN
        ZoneSizingInput(ZoneSizIndex)%DesCoolAirFlow =  0.0
      ELSEIF (rNumericArgs(7) < 0.0d0) THEN
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(7))//': '//  &
        TRIM(RoundSigDigits(rNumericArgs(10),2)))
        CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
        ErrorsFound=.true.
      ELSE
        ZoneSizingInput(ZoneSizIndex)%DesCoolAirFlow = rNumericArgs(7)
      ENDIF
      !  N8,\field Cooling Minimum Air Flow per Zone Floor Area
      !      \type real
      !      \units m3/s-m2
      !      \minimum 0
      !      \default .000762
      !      \note default is .15 cfm/ft2
      !      \note This input is used if Cooling Design Air Flow Method is design day with limit
      IF (lNumericFieldBlanks(8)) THEN
        IF (rNumericArgs(11) <= 0.0) THEN  ! in case someone changes the default in the IDD
          ZoneSizingInput(ZoneSizIndex)%DesCoolMinAirFlowPerArea =  .000762d0
        ELSE
          ZoneSizingInput(ZoneSizIndex)%DesCoolMinAirFlowPerArea = rNumericArgs(8)
        ENDIF
      ELSEIF (rNumericArgs(8) < 0.0d0) THEN
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(8))//': '//  &
        TRIM(RoundSigDigits(rNumericArgs(8),2)))
        CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
        ErrorsFound=.true.
      ELSE
        ZoneSizingInput(ZoneSizIndex)%DesCoolMinAirFlowPerArea = rNumericArgs(8)
      ENDIF
      !  N9,\field Cooling Minimum Air Flow
      !      \type real
      !      \units m3/s
      !      \minimum 0
      !      \default 0
      !      \note This input is used if Cooling Design Air Flow Method is design day with limit
      IF (lNumericFieldBlanks(9)) THEN
        ZoneSizingInput(ZoneSizIndex)%DesCoolMinAirFlow =  0.0
      ELSEIF (rNumericArgs(9) < 0.0d0) THEN
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(9))//': '//  &
        TRIM(RoundSigDigits(rNumericArgs(9),2)))
        CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
        ErrorsFound=.true.
      ELSE
        ZoneSizingInput(ZoneSizIndex)%DesCoolMinAirFlow = rNumericArgs(9)
      ENDIF
      !  N10,\field Cooling Minimum Air Flow Fraction
      !      \note fraction of the Cooling design Air Flow Rate
      !      \type real
      !      \minimum 0
      !      \default 0
      !      \note This input is currently used in sizing the Fan minimum Flow Rate.
      !      \note It does not currently affect other component autosizing.
      IF (lNumericFieldBlanks(10)) THEN
        ZoneSizingInput(ZoneSizIndex)%DesCoolMinAirFlowFrac =   0.0
      ELSEIF (rNumericArgs(10) < 0.0d0) THEN
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(10))//': '//  &
        TRIM(RoundSigDigits(rNumericArgs(10),2)))
        CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
        ErrorsFound=.true.
      ELSE
        ZoneSizingInput(ZoneSizIndex)%DesCoolMinAirFlowFrac = rNumericArgs(10)
      ENDIF
      !  N11,\field Heating Design Air Flow Rate
      !      \type real
      !      \units m3/s
      !      \minimum 0
      !      \default 0
      !      \note This input is used if Heating Design Air Flow Method is Flow/Zone.
      !      \note This value will be multiplied by the global or zone sizing factor and
      !      \note by zone multipliers.
      IF (lNumericFieldBlanks(11)) THEN
        ZoneSizingInput(ZoneSizIndex)%DesHeatAirFlow = 0.0
      ELSEIF (rNumericArgs(11) < 0.0d0) THEN
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(11))//': '//  &
        TRIM(RoundSigDigits(rNumericArgs(11),2)))
        CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
        ErrorsFound=.true.
      ELSE
        ZoneSizingInput(ZoneSizIndex)%DesHeatAirFlow = rNumericArgs(11)
      ENDIF
      !  N12,\field Heating Maximum Air Flow per Zone Floor Area
      !      \type real
      !      \units m3/s-m2
      !      \minimum 0
      !      \default .002032
      !      \note default is .40 cfm/ft2
      !      \note This input is not currently used for autosizing any of the components.
      IF (lNumericFieldBlanks(12)) THEN
        IF (rNumericArgs(12) <= 0.0) THEN  ! in case someone changes the default in the IDD
          ZoneSizingInput(ZoneSizIndex)%DesHeatMaxAirFlowPerArea = 0.002032d0
        ELSE
          ZoneSizingInput(ZoneSizIndex)%DesHeatMaxAirFlowPerArea = rNumericArgs(12)
        ENDIF
      ELSEIF (rNumericArgs(12) < 0.0d0) THEN
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(12))//': '//  &
        TRIM(RoundSigDigits(rNumericArgs(12),2)))
        CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
        ErrorsFound=.true.
      ELSE
        ZoneSizingInput(ZoneSizIndex)%DesHeatMaxAirFlowPerArea = rNumericArgs(12)
      ENDIF
      !  N13,\field Heating Maximum Air Flow
      !      \type real
      !      \units m3/s
      !      \minimum 0
      !      \default .1415762
      !      \note default is 300 cfm
      !      \note This input is not currently used for autosizing any of the components.
      IF (lNumericFieldBlanks(13)) THEN
        IF (rNumericArgs(13) <= 0.0) THEN  ! in case someone changes the default in the IDD
          ZoneSizingInput(ZoneSizIndex)%DesHeatMaxAirFlow = 0.1415762d0
        ELSE
          ZoneSizingInput(ZoneSizIndex)%DesHeatMaxAirFlow = rNumericArgs(13)
        ENDIF
      ELSEIF (rNumericArgs(13) < 0.0d0) THEN
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(13))//': '//  &
        TRIM(RoundSigDigits(rNumericArgs(13),2)))
        CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
        ErrorsFound=.true.
      ELSE
        ZoneSizingInput(ZoneSizIndex)%DesHeatMaxAirFlow = rNumericArgs(13)
      ENDIF
      !  N14;\field Heating Maximum Air Flow Fraction
      !      \note fraction of the Heating Design Air Flow Rate
      !      \note This input is not currently used for autosizing any of the components.
      !      \type real
      !      \minimum 0
      !      \default 0.3
      IF (lNumericFieldBlanks(14)) THEN
        IF (rNumericArgs(14) <= 0.0) THEN  ! in case someone changes the default in the IDD
          ZoneSizingInput(ZoneSizIndex)%DesHeatMaxAirFlowFrac = 0.3d0
        ELSE
          ZoneSizingInput(ZoneSizIndex)%DesHeatMaxAirFlowFrac = rNumericArgs(14)
        ENDIF
      ELSEIF (rNumericArgs(14) < 0.0d0) THEN
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(14))//': '//  &
        TRIM(RoundSigDigits(rNumericArgs(14),2)))
        CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
        ErrorsFound=.true.
      ELSE
        ZoneSizingInput(ZoneSizIndex)%DesHeatMaxAirFlowFrac = rNumericArgs(14)
      ENDIF

      !  A5, \field Zone Air Distribution Object Name
      IF (.NOT. lAlphaFieldBlanks(5)) THEN
        ZoneSizingInput(ZoneSizIndex)%ZoneAirDistEffObjName = cAlphaArgs(5)
        ObjIndex=FindItemInList(ZoneSizingInput(ZoneSizIndex)%ZoneAirDistEffObjName,  &
        ZoneAirDistribution%Name,numZoneAirDistribution)
        IF (ObjIndex > 0) THEN
          ZoneSizingInput(ZoneSizIndex)%ZoneADEffCooling = ZoneAirDistribution(ObjIndex)%ZoneADEffCooling
          ZoneSizingInput(ZoneSizIndex)%ZoneADEffHeating = ZoneAirDistribution(ObjIndex)%ZoneADEffHeating
          ZoneSizingInput(ZoneSizIndex)%ZoneSecondaryRecirculation = ZoneAirDistribution(ObjIndex)%ZoneSecondaryRecirculation
          ZoneSizingInput(ZoneSizIndex)%ZoneAirDistributionIndex = ObjIndex
        ELSE
          ! generate a warning message
          CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cAlphaFieldNames(5))//'="'//  &
          TRIM(cAlphaArgs(5))//'" not found.')
          ErrorsFound=.true.
        ENDIF
      ELSE
        ! assume defaults
        ZoneSizingInput(ZoneSizIndex)%ZoneADEffCooling = 1.0
        ZoneSizingInput(ZoneSizIndex)%ZoneADEffHeating = 1.0
        ZoneSizingInput(ZoneSizIndex)%ZoneSecondaryRecirculation = 0.0
      ENDIF

      SELECT CASE(TRIM(cAlphaArgs(3)))
      CASE('DESIGNDAY')
        ZoneSizingInput(ZoneSizIndex)%CoolAirDesMethod = FromDDCalc
      CASE('FLOW/ZONE')
        ZoneSizingInput(ZoneSizIndex)%CoolAirDesMethod = InpDesAirFlow
      CASE('DESIGNDAYWITHLIMIT')
        ZoneSizingInput(ZoneSizIndex)%CoolAirDesMethod = DesAirFlowWithLim
      CASE DEFAULT
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cAlphaFieldNames(3))//': '//TRIM(cAlphaArgs(3)))
        CALL ShowContinueError('Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
        ErrorsFound=.true.
      END SELECT
      SELECT CASE(TRIM(cAlphaArgs(4)))
      CASE('DESIGNDAY')
        ZoneSizingInput(ZoneSizIndex)%HeatAirDesMethod = FromDDCalc
      CASE('FLOW/ZONE')
        ZoneSizingInput(ZoneSizIndex)%HeatAirDesMethod = InpDesAirFlow
      CASE('DESIGNDAYWITHLIMIT')
        ZoneSizingInput(ZoneSizIndex)%HeatAirDesMethod = DesAirFlowWithLim
      CASE DEFAULT
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cAlphaFieldNames(4))//': '//TRIM(cAlphaArgs(4)))
        CALL ShowContinueError('Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
        ErrorsFound=.true.
      END SELECT
    END DO
  END DO
ENDIF

IF (ErrorsFound) THEN
  !CALL ShowFatalError(TRIM(cCurrentModuleObject)//': Errors found in getting input. Program terminates.') !RS: Secret Search String
  WRITE(DebugFile,*) TRIM(cCurrentModuleObject)//': Errors found in getting input. Program terminates (at least wishes to!).'
END IF


RETURN

END SUBROUTINE GetZoneSizingInput

SUBROUTINE GetZoneAndZoneListNames(ErrorsFound,NumZones,ZoneNames,NumZoneLists,ZoneListNames)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   October 2010
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! Get Zone and ZoneList Names so Sizing:Zone can use global ZoneList.
  ! This is not a full validation of these objects -- only enough to fill
  ! structures for the Sizing:Zone object.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE InputProcessor, ONLY: GetNumObjectsFound, GetObjectItem, VerifyName, FindItemInList
  USE DataIPShortCuts

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT) :: ErrorsFound
  INTEGER, INTENT(INOUT) :: NumZones
  CHARACTER(len=*), ALLOCATABLE, DIMENSION(:) :: ZoneNames
  INTEGER, INTENT(INOUT) :: NumZoneLists
  TYPE(ZoneListData), ALLOCATABLE, DIMENSION(:) :: ZoneListNames

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: Item
  INTEGER :: Found
  INTEGER :: Item1
  INTEGER :: NumAlphas
  INTEGER :: NumNumbers
  INTEGER :: IOStatus
  LOGICAL :: InErrFlag   ! Preserve (no current use) the input status of ErrorsFound

  InErrFlag=ErrorsFound
  cCurrentModuleObject='Zone'
  NumZones=GetNumObjectsFound(cCurrentModuleObject)
  ALLOCATE(ZoneNames(NumZones))

  DO Item=1,NumZones
    CALL GetObjectItem(TRIM(cCurrentModuleObject),Item,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,IOStatus,  &
    AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
    AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
    ! validation, but no error
    Found=FindItemInList(cAlphaArgs(1),ZoneNames,Item-1)
    IF (Found == 0) THEN
      ZoneNames(Item)=cAlphaArgs(1)
    ELSE
      ZoneNames(Item)='xxxxx'
    ENDIF
  ENDDO

  cCurrentModuleObject='ZoneList'
  NumZoneLists=GetNumObjectsFound(cCurrentModuleObject)
  ALLOCATE(ZoneListNames(NumZoneLists))

  DO Item=1,NumZoneLists
    CALL GetObjectItem(TRIM(cCurrentModuleObject),Item,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,IOStatus,  &
    AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
    AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
    ! validation, but no error
    Found=FindItemInList(cAlphaArgs(1),ZoneListNames%Name,Item-1)
    IF (Found == 0) THEN
      ZoneListNames(Item)%Name=cAlphaArgs(1)
    ELSE
      ZoneListNames(Item)%Name='xxxxx'
    ENDIF
    ALLOCATE(ZoneListNames(Item)%Zones(NumAlphas-1))
    ZoneListNames(Item)%NumOfZones=NumAlphas-1
    DO Item1=2,NumAlphas
      Found=FindItemInList(cAlphaArgs(Item1),ZoneNames,NumZones)
      ZoneListNames(Item)%Zones(Item1-1)=Found
    ENDDO
  ENDDO

  RETURN

END SUBROUTINE GetZoneAndZoneListNames

SUBROUTINE GetSystemSizingInput

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Fred Buhl
  !       DATE WRITTEN   January 2001
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! Obtains input data for System Sizing objects and stores it in
  ! appropriate data structures.

  ! METHODOLOGY EMPLOYED:
  ! Uses InputProcessor "Get" routines to obtain data.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataIPShortCuts
  USE InputProcessor, ONLY: GetNumObjectsFound, GetObjectItem, VerifyName
  USE General, ONLY: RoundSigDigits

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: SysSizIndex ! loop index
  INTEGER :: NumAlphas        ! Number of Alphas for each GetObjectItem call
  INTEGER :: NumNumbers       ! Number of Numbers for each GetObjectItem call
  INTEGER :: IOStatus         ! Used in GetObjectItem
  LOGICAL :: ErrorsFound=.false.  ! Set to true if errors in input, fatal at end of routine
  LOGICAL :: IsNotOK                            ! Flag to verify name
  LOGICAL :: IsBlank                            ! Flag for blank name
  INTEGER :: NumDesDays       ! Number of design days in input

  NumAirLoops = GetNumObjectsFound('AirLoopHVAC')
  cCurrentModuleObject='Sizing:System'
  NumSysSizInput = GetNumObjectsFound(TRIM(cCurrentModuleObject))

  IF (NumSysSizInput > 0) THEN
    NumDesDays = GetNumObjectsFound('SizingPeriod:DesignDay') + GetNumObjectsFound('SizingPeriod:WeatherFileDays') +   &
    GetNumObjectsFound('SizingPeriod:WeatherFileConditionType')
    IF (NumDesDays == 0 .AND. (DoSystemSizing .OR. DoPlantSizing) ) THEN
      CALL ShowSevereError('System Sizing calculations need SizingPeriod:* input.')
      ErrorsFound = .TRUE.
    END IF
    ALLOCATE(SysSizInput(NumSysSizInput))
  END IF

  DO SysSizIndex=1,NumSysSizInput
    CALL GetObjectItem(TRIM(cCurrentModuleObject),SysSizIndex,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,IOStatus,  &
    AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
    AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
    IsNotOK=.FALSE.
    IsBlank=.FALSE.
    CALL VerifyName(cAlphaArgs(1),SysSizInput%AirPriLoopName,SysSizIndex-1,IsNotOK,IsBlank,TRIM(cCurrentModuleObject)//' Name')
    IF (IsNotOK) THEN
      ErrorsFound=.true.
      IF (IsBlank) cAlphaArgs(1)='xxxxx'
    ENDIF

    SysSizInput(SysSizIndex)%AirPriLoopName = cAlphaArgs(1)
    SELECT CASE(TRIM(cAlphaArgs(2)))
    CASE('SENSIBLE')
      SysSizInput(SysSizIndex)%LoadSizeType = Sensible
    CASE('LATENT')
      SysSizInput(SysSizIndex)%LoadSizeType = Latent
    CASE('TOTAL')
      SysSizInput(SysSizIndex)%LoadSizeType = Total
    CASE('VENTILATIONREQUIREMENT')
      SysSizInput(SysSizIndex)%LoadSizeType = Ventilation
    CASE DEFAULT
      CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cAlphaFieldNames(2))//': '//TRIM(cAlphaArgs(2)))
      CALL ShowContinueError('Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
      ErrorsFound=.true.
    END SELECT
    SELECT CASE(TRIM(cAlphaArgs(3)))
    CASE('COINCIDENT')
      SysSizInput(SysSizIndex)%SizingOption = 2
    CASE('NONCOINCIDENT')
      SysSizInput(SysSizIndex)%SizingOption = 1
    CASE DEFAULT
      CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cAlphaFieldNames(3))//': '//TRIM(cAlphaArgs(3)))
      CALL ShowContinueError('Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
      ErrorsFound=.true.
    END SELECT
    SELECT CASE(TRIM(cAlphaArgs(4)))
    CASE('YES')
      SysSizInput(SysSizIndex)%CoolOAOption = 1
    CASE('NO')
      SysSizInput(SysSizIndex)%CoolOAOption = 2
    CASE DEFAULT
      CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cAlphaFieldNames(4))//': '//TRIM(cAlphaArgs(4)))
      CALL ShowContinueError('Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
      ErrorsFound=.true.
    END SELECT
    SELECT CASE(TRIM(cAlphaArgs(5)))
    CASE('YES')
      SysSizInput(SysSizIndex)%HeatOAOption = 1
    CASE('NO')
      SysSizInput(SysSizIndex)%HeatOAOption = 2
    CASE DEFAULT
      CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cAlphaFieldNames(5))//': '//TRIM(cAlphaArgs(5)))
      CALL ShowContinueError('Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
      ErrorsFound=.true.
    END SELECT

    !  N1, \field Design Outdoor Air Flow Rate
    !      \type real
    !      \default autosize
    !      \minimum 0.0
    IF (lNumericFieldBlanks(1)) THEN
      SysSizInput(SysSizIndex)%DesOutAirVolFlow = autosize
    ELSEIF (rNumericArgs(1) < 0.0d0 .and. rNumericArgs(1) /= autosize) THEN
      CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(1))//': '//  &
      TRIM(RoundSigDigits(rNumericArgs(1),2)))
      CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
      ErrorsFound=.true.
    ELSE
      SysSizInput(SysSizIndex)%DesOutAirVolFlow = rNumericArgs(1)
    ENDIF
    IF (SysSizInput(SysSizIndex)%DesOutAirVolFlow == autosize) THEN
      SysSizInput(SysSizIndex)%OAAutosized = .TRUE.
    END IF

    !  N2, \field Minimum System Air Flow Ratio
    !      \required-field
    !      \type real
    !      \minimum 0.0
    !      \maximum 1.0
    IF (lNumericFieldBlanks(2)) THEN
      SysSizInput(SysSizIndex)%SysAirMinFlowRat = 0.0
    ELSEIF (rNumericArgs(2) < 0.0) THEN
      CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(2))//': '//  &
      TRIM(RoundSigDigits(rNumericArgs(2),2)))
      CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
      ErrorsFound=.true.
    ELSE
      SysSizInput(SysSizIndex)%SysAirMinFlowRat = rNumericArgs(2)
    ENDIF
    SysSizInput(SysSizIndex)%PreheatTemp = rNumericArgs(3)
    SysSizInput(SysSizIndex)%PreheatHumRat = rNumericArgs(4)
    SysSizInput(SysSizIndex)%PrecoolTemp = rNumericArgs(5)
    SysSizInput(SysSizIndex)%PrecoolHumRat = rNumericArgs(6)
    SysSizInput(SysSizIndex)%CoolSupTemp = rNumericArgs(7)
    SysSizInput(SysSizIndex)%HeatSupTemp = rNumericArgs(8)
    SysSizInput(SysSizIndex)%CoolSupHumRat = rNumericArgs(9)
    SysSizInput(SysSizIndex)%HeatSupHumRat = rNumericArgs(10)
    !  N11, \field Cooling Design Air Flow Rate
    !      \note This input is used if Cooling Design Air Flow Method is Flow/System
    !      \note This value will *not* be multiplied by any sizing factor or by zone multipliers.
    !      \note If using zone multipliers, this value must be large enough to serve the multiplied zones.
    !      \type real
    !      \units m3/s
    !      \minimum 0
    !      \default 0
    IF (lNumericFieldBlanks(11)) THEN
      SysSizInput(SysSizIndex)%DesCoolAirFlow = 0.0
    ELSEIF (rNumericArgs(11) < 0.0) THEN
      CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(11))//': '//  &
      TRIM(RoundSigDigits(rNumericArgs(11),2)))
      CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
      ErrorsFound=.true.
    ELSE
      SysSizInput(SysSizIndex)%DesCoolAirFlow = rNumericArgs(11)
    ENDIF
    !  N12;\field Heating Design Air Flow Rate
    !      \note This input is used if Heating Design Air Flow Method is Flow/System
    !      \note This value will *not* be multiplied by any sizing factor or by zone multipliers.
    !      \note If using zone multipliers, this value must be large enough to serve the multiplied zones.
    !      \type real
    !      \units m3/s
    !      \minimum 0
    !      \default 0
    IF (lNumericFieldBlanks(12)) THEN
      SysSizInput(SysSizIndex)%DesHeatAirFlow = 0.0
    ELSEIF (rNumericArgs(12) < 0.0) THEN
      CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(12))//': '//  &
      TRIM(RoundSigDigits(rNumericArgs(12),2)))
      CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
      ErrorsFound=.true.
    ELSE
      SysSizInput(SysSizIndex)%DesHeatAirFlow = rNumericArgs(12)
    ENDIF
    !  N13;\field Maximum Zone Outdoor Air Fraction
    !      \type real
    !      \default 1.0
    !      \minimum> 0.0
    !      \units dimensionless
    IF (lNumericFieldBlanks(13)) THEN
      SysSizInput(SysSizIndex)%MaxZoneOAFraction = 0.0
    ELSEIF (rNumericArgs(13) < 0.0) THEN
      CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cNumericFieldNames(13))//': '//  &
      TRIM(RoundSigDigits(rNumericArgs(13),2)))
      CALL ShowContinueError('.. value should not be negative. Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
      ErrorsFound=.true.
    ELSE
      SysSizInput(SysSizIndex)%MaxZoneOAFraction = rNumericArgs(13)
    ENDIF
    SELECT CASE(TRIM(cAlphaArgs(6)))
    CASE('DESIGNDAY')
      SysSizInput(SysSizIndex)%CoolAirDesMethod = FromDDCalc
    CASE('FLOW/SYSTEM')
      SysSizInput(SysSizIndex)%CoolAirDesMethod = InpDesAirFlow
    CASE DEFAULT
      CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cAlphaFieldNames(6))//': '//TRIM(cAlphaArgs(6)))
      CALL ShowContinueError('Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
      ErrorsFound=.true.
    END SELECT
    SELECT CASE(TRIM(cAlphaArgs(7)))
    CASE('DESIGNDAY')
      SysSizInput(SysSizIndex)%HeatAirDesMethod = FromDDCalc
    CASE('FLOW/SYSTEM')
      SysSizInput(SysSizIndex)%HeatAirDesMethod = InpDesAirFlow
    CASE DEFAULT
      CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cAlphaFieldNames(7))//': '//TRIM(cAlphaArgs(7)))
      CALL ShowContinueError('Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
      ErrorsFound=.true.
    END SELECT
    SELECT CASE(TRIM(cAlphaArgs(8)))
    CASE('ZONESUM')
      SysSizInput(SysSizIndex)%SystemOAMethod = 1
    CASE('VENTILATIONRATEPROCEDURE')
      SysSizInput(SysSizIndex)%SystemOAMethod = 2
      IF(SysSizInput(SysSizIndex)%DesOutAirVolFlow > 0) THEN
        CALL ShowWarningError(TRIM(cCurrentModuleObject)//  &
        '- Design Outdoor Airflow Rate value will be overridden since SystemOAMethod is set to VRP')
      END IF
    CASE DEFAULT
      CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cAlphaFieldNames(8))//': '//TRIM(cAlphaArgs(8)))
      CALL ShowContinueError('Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
      ErrorsFound=.true.
    END SELECT
  END DO

  IF (ErrorsFound) THEN
    CALL ShowFatalError(TRIM(cCurrentModuleObject)//': Errors found in getting input. Program terminates.')
  END IF

  RETURN

END SUBROUTINE GetSystemSizingInput

SUBROUTINE GetPlantSizingInput

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Fred Buhl
  !       DATE WRITTEN   October 2001
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! Obtains input data for Plant Sizing objects and stores it in
  ! appropriate data structures.

  ! METHODOLOGY EMPLOYED:
  ! Uses InputProcessor "Get" routines to obtain data.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE InputProcessor, ONLY: GetNumObjectsFound, GetObjectItem, VerifyName
  USE DataIPShortCuts

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: PltSizIndex ! loop index
  INTEGER :: NumAlphas        ! Number of Alphas for each GetObjectItem call
  INTEGER :: NumNumbers       ! Number of Numbers for each GetObjectItem call
  INTEGER :: IOStatus         ! Used in GetObjectItem
  LOGICAL :: ErrorsFound=.false.  ! Set to true if errors in input, fatal at end of routine
  LOGICAL :: IsNotOK          ! Flag to verify name
  LOGICAL :: IsBlank          ! Flag for blank name
  INTEGER :: NumDesDays       ! Number of design days in input

  cCurrentModuleObject='Sizing:Plant'
  NumPltSizInput = GetNumObjectsFound(TRIM(cCurrentModuleObject))

  IF (NumPltSizInput > 0) THEN
    NumDesDays = GetNumObjectsFound('SizingPeriod:DesignDay') + GetNumObjectsFound('SizingPeriod:WeatherFileDays') +   &
    GetNumObjectsFound('SizingPeriod:WeatherFileConditionType')
    IF (NumDesDays == 0 .AND. DoPlantSizing ) THEN
      CALL ShowSevereError('Plant Sizing calculations need SizingPeriod:* input')
      ErrorsFound = .TRUE.
    END IF
    ALLOCATE(PlantSizData(NumPltSizInput))
    PlantSizData%PlantLoopName = ' '
    PlantSizData%ExitTemp = 0.0
    PlantSizData%DeltaT = 0.0
    PlantSizData%LoopType = 0
    PlantSizData%DesVolFlowRate = 0.0
  END IF

  DO PltSizIndex=1,NumPltSizInput
    CALL GetObjectItem(TRIM(cCurrentModuleObject),PltSizIndex,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,IOStatus,  &
    AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
    AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
    IsNotOK=.FALSE.
    IsBlank=.FALSE.
    CALL VerifyName(cAlphaArgs(1),PlantSizData%PlantLoopName,PltSizIndex-1,IsNotOK,IsBlank,TRIM(cCurrentModuleObject)//' Name')
    IF (IsNotOK) THEN
      ErrorsFound=.true.
      IF (IsBlank) cAlphaArgs(1)='xxxxx'
    ENDIF
    PlantSizData(PltSizIndex)%PlantLoopName = cAlphaArgs(1)
    PlantSizData(PltSizIndex)%ExitTemp = rNumericArgs(1)
    PlantSizData(PltSizIndex)%DeltaT = rNumericArgs(2)
    SELECT CASE(TRIM(cAlphaArgs(2)))
    CASE('HEATING')
      PlantSizData(PltSizIndex)%LoopType = HeatingLoop
    CASE('COOLING')
      PlantSizData(PltSizIndex)%LoopType = CoolingLoop
    CASE('CONDENSER')
      PlantSizData(PltSizIndex)%LoopType = CondenserLoop
    CASE('STEAM')
      PlantSizData(PltSizIndex)%LoopType = SteamLoop
    CASE DEFAULT
      CALL ShowSevereError(TRIM(cCurrentModuleObject)//': incorrect '//TRIM(cAlphaFieldNames(2))//': '//TRIM(cAlphaArgs(2)))
      CALL ShowContinueError('Occurs in Sizing Object='//TRIM(cAlphaArgs(1)))
      ErrorsFound=.true.
    END SELECT

    CALL SetupEMSInternalVariable('Plant Design Volume Flow Rate',  &
    PlantSizData(PltSizIndex)%PlantLoopName, '[m3/s]', &
    PlantSizData(PltSizIndex)%DesVolFlowRate )
  END DO

  IF (ErrorsFound) THEN
    CALL ShowFatalError(TRIM(cCurrentModuleObject)//': Errors found in getting input. Program terminates.')
  END IF

  RETURN

END SUBROUTINE GetPlantSizingInput

SUBROUTINE SetupZoneSizing(ErrorsFound)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         L. Lawrie/F. Buhl
  !       DATE WRITTEN   March 2010
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  !  execute a few (1) time steps of a simulation to facilitate setting up model for zone sizing
  !  developed to resolve reverse DD problems caused be the differences
  !  that stem from setup and information gathering that occurs during the first pass.

  ! METHODOLOGY EMPLOYED:
  ! Using global flag (kickoff sizing simulation), only a few time steps are executed.
  ! global flag is used in other parts of simulation to terminate quickly.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataEnvironment ,    ONLY: EndMonthFlag
  USE InputProcessor,      ONLY: GetNumRangeCheckErrorsFound
  USE CostEstimateManager, ONLY: SimCostEstimate

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT) :: ErrorsFound

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  LOGICAL :: Available=.false. ! an environment is available to process

  !  return  ! remove comment to do "old way"

  Available = .true.

  CurOverallSimDay=0
  DO WHILE (Available)  ! do for each environment

    CALL GetNextEnvironment(Available,ErrorsFound)

    IF (.not. Available) EXIT
    IF (ErrorsFound) EXIT

    ! check that environment is one of the design days
    IF (KindOfSim == ksRunPeriodWeather) THEN
      CYCLE
    ENDIF

    BeginEnvrnFlag = .TRUE.
    EndEnvrnFlag   = .FALSE.
    EndMonthFlag   = .FALSE.
    WarmupFlag     = .TRUE.
    DayOfSim       =  0

    CurEnvirNumSimDay=1
    CurOverallSimDay=CurOverallSimDay+1

    DayOfSim     = DayOfSim + 1
    BeginDayFlag = .TRUE.
    EndDayFlag   = .FALSE.

    HourOfDay = 1

    BeginHourFlag = .TRUE.
    EndHourFlag   = .FALSE.

    TimeStep = 1

    BeginTimeStepFlag = .TRUE.

    CALL ManageWeather

    CALL ManageHeatBalance

    BeginHourFlag  = .FALSE.
    BeginDayFlag   = .FALSE.
    BeginEnvrnFlag = .FALSE.
    BeginSimFlag   = .FALSE.
    BeginFullSimFlag = .FALSE.

    !          ! do another timestep=1
    CALL ManageWeather

    CALL ManageHeatBalance

    !         do an end of day, end of environment time step

    HourOfDay=24
    TimeStep=NumOfTimeStepInHour
    EndEnvrnFlag   = .True.

    CALL ManageWeather

    CALL ManageHeatBalance

  END DO                        ! ... End environment loop.

  RETURN

END SUBROUTINE SetupZoneSizing

SUBROUTINE ReportZoneSizing(ZoneName,LoadType,CalcDesLoad,UserDesLoad,CalcDesFlow,UserDesFlow,DesDayName,  &
  PeakHrMin,PeakTemp,PeakHumRat,FloorArea, TotOccs, MinOAVolFlow)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Fred Buhl
  !       DATE WRITTEN   Decenber 2001
  !       MODIFIED       August 2008, Greg Stark
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine writes one item of zone sizing data to the "eio" file..

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataPrecisionGlobals
  USE DataGlobals_HPSimIntegrated, ONLY : OutputFileInits
  USE DataStringGlobals, ONLY: VerString
  USE General, ONLY: RoundSigDigits

  ! BSLLC Start
  USE SQLiteProcedures
  ! BSLLC Finish

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: ZoneName     ! the name of the zone
  CHARACTER(len=*), INTENT(IN) :: LoadType     ! the description of the input variable
  REAL(r64), INTENT(IN)        :: CalcDesLoad  ! the value from the sizing calculation [W]
  REAL(r64), INTENT(IN)        :: UserDesLoad  ! the value from the sizing calculation modified by user input [W]
  REAL(r64), INTENT(IN)        :: CalcDesFlow  ! calculated design air flow rate [m3/s]
  REAL(r64), INTENT(IN)        :: UserDesFlow  ! user input or modified design air flow rate [m3/s]
  CHARACTER(len=*), INTENT(IN) :: DesDayName   ! the name of the design day that produced the peak
  CHARACTER(len=*), INTENT(IN) :: PeakHrMin    ! time stamp of the peak
  REAL(r64), INTENT(IN)        :: PeakTemp     ! temperature at peak [C]
  REAL(r64), INTENT(IN)        :: PeakHumRat   ! humidity ratio at peak [kg water/kg dry air]
  REAL(r64), INTENT(IN)        :: FloorArea    ! zone floor area [m2]
  REAL(r64), INTENT(IN)        :: TotOccs      ! design number of occupants for the zone
  REAL(r64), INTENT(IN)        :: MinOAVolFlow ! zone design minimum outside air flow rate [m3/s]

  ! SUBROUTINE PARAMETER DEFINITIONS:

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  LOGICAL, SAVE :: MyOneTimeFlag = .TRUE.

  IF (MyOneTimeFlag) THEN
    WRITE(OutputFileInits, 990)
    MyOneTimeFlag = .FALSE.
  END IF

  WRITE (OutputFileInits, 991) TRIM(ZoneName), TRIM(LoadType), TRIM(RoundSigDigits(CalcDesLoad,5)),  &
  TRIM(RoundSigDigits(UserDesLoad,5)),   &
  TRIM(RoundSigDigits(CalcDesFlow,5)), TRIM(RoundSigDigits(UserDesFlow,5)), TRIM(DesDayName), TRIM(PeakHrMin),  &
  TRIM(RoundSigDigits(PeakTemp,5)), TRIM(RoundSigDigits(PeakHumRat,5)), TRIM(RoundSigDigits(FloorArea,5)),   &
  TRIM(RoundSigDigits(TotOccs,5)), TRIM(RoundSigDigits(MinOAVolFlow,5))

  ! BSLLC Start
  IF (WriteOutputToSQLite) THEN
    CALL AddSQLiteZoneSizingRecord (ZoneName, LoadType, CalcDesLoad, UserDesLoad, CalcDesFlow, UserDesFlow, DesDayName, &
    PeakHrMin, PeakTemp, PeakHumRat, MinOAVolFlow)
  END IF
  ! BSLLC Finish

  990 FORMAT('! <Zone Sizing Information>, Zone Name, Load Type, Calc Des Load {W}, User Des Load {W}, ',  &
  'Calc Des Air Flow Rate {m3/s}, ', &
  'User Des Air Flow Rate {m3/s}, Design Day Name, Date/Time of Peak, Temperature at Peak {C}, ', &
  'Humidity Ratio at Peak {kgWater/kgDryAir}, Floor Area {m2}, # Occupants, Calc Outdoor Air Flow Rate {m3/s}')
  991 FORMAT(' Zone Sizing Information',13(', ',A))

  RETURN

END SUBROUTINE ReportZoneSizing

SUBROUTINE ReportSysSizing(SysName,VarDesc,VarValue)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Fred Buhl
  !       DATE WRITTEN   January 2003
  !       MODIFIED       August 2008, Greg Stark
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine writes one item of system sizing data to the "eio" file..

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataPrecisionGlobals
  USE DataGlobals_HPSimIntegrated, ONLY : OutputFileInits
  USE DataStringGlobals, ONLY: VerString
  USE General, ONLY: RoundSigDigits

  ! BSLLC Start
  USE SQLiteProcedures
  ! BSLLC Finish

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: SysName   ! the name of the zone
  CHARACTER(len=*), INTENT(IN) :: VarDesc   ! the description of the input variable
  REAL(r64), INTENT(IN)        :: VarValue  ! the value from the sizing calculation

  ! SUBROUTINE PARAMETER DEFINITIONS:

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  LOGICAL, SAVE :: MyOneTimeFlag = .TRUE.

  IF (MyOneTimeFlag) THEN
    WRITE(OutputFileInits, 990)
    MyOneTimeFlag = .FALSE.
  END IF

  WRITE (OutputFileInits, 991) TRIM(SysName), TRIM(VarDesc), TRIM(RoundSigDigits(VarValue,5))

  ! BSLLC Start
  IF (WriteOutputToSQLite) CALL AddSQLiteSystemSizingRecord (SysName, VarDesc, VarValue)
  ! BSLLC Finish

  990 FORMAT('! <System Sizing Information>, System Name, ', 'Field Description, Value')
  991 FORMAT(' System Sizing Information',3(', ',A))

  RETURN

END SUBROUTINE ReportSysSizing

!     NOTICE
!
!     Copyright © 1996-2012 The Board of Trustees of the University of Illinois
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

END MODULE SizingManager
