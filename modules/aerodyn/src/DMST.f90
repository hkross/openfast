!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2017  Envision Energy USA, LTD
!
!    This file is part of AeroDyn.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!**********************************************************************************************************************************
module DMST
    
   use NWTC_Library
   use DMST_Types
   use AirfoilInfo

   implicit none
   
   private
   
   type(ProgDesc), parameter  :: DMST_Ver = ProgDesc( 'DMST', '', '' )
   character(*),   parameter  :: DMST_Nickname = 'DMST'
   
   ! ..... Public Subroutines .....................................................................................................

   public :: DMST_Init                           ! Initialization routine
!   public :: DMST_End                            ! Ending routine (includes clean up)
!   public :: DMST_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                  !   continuous states, and updating discrete states
!   public :: DMST_CalcOutput                     ! Routine for computing outputs

   contains

!----------------------------------------------------------------------------------------------------------------------------------   
subroutine DMST_SetParameters( InitInp, p, errStat, errMsg )
! This routine is called from DMST_Init.
! The parameters are set here and not changed during the simulation.
!..................................................................................................................................
   type(DMST_InitInputType),       intent(in   )  :: InitInp     ! Input for initialization routine
   type(DMST_ParameterType),       intent(  out)  :: p           ! Parameters
   integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
   character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables
   integer(IntKi)                                 :: errStat2    ! Temporary error status of the operation
   character(*), parameter                        :: RoutineName = 'DMST_SetParameters'
   integer(IntKi)                                 :: i, j

      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""

   p%numBlades      = InitInp%numBlades   
   p%numBladeNodes  = InitInp%numBladeNodes 
   p%airDens        = InitInp%airDens
   p%kinVisc        = InitInp%kinVisc
   p%Nst            = InitInp%Nst

   allocate ( p%AFindx(p%numBladeNodes, p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%AFindx.', errStat, errMsg, RoutineName )
      return
   end if

   allocate ( p%chord(p%numBladeNodes,1), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%chord.', errStat, errMsg, RoutineName )
      return
   end if 

   allocate ( p%radius(p%numBladeNodes), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%radius.', errStat, errMsg, RoutineName )
      return
   end if 

   allocate ( p%dTheta(p%numBladeNodes), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%dTheta.', errStat, errMsg, RoutineName )
      return
   end if 

   allocate ( p%theta_st(p%numBladeNodes, 2_IntKi*p%Nst), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%theta_st.', errStat, errMsg, RoutineName )
      return
   end if 
   
   p%AFindx = InitInp%AFindx 
   
      ! Compute the radius, total streamtube angle, and azimuthal position of streamtube midpoint
   do i=1,p%numBladeNodes
      p%chord(i,1) = InitInp%chord(i,1)
      p%radius(i) = InitInp%radius(i)
      p%dTheta(i) = pi/p%Nst
      do j=1,p%Nst
         p%theta_st(i,j) = p%dTheta(i)/2.0_ReKi + p%dTheta(i)*(j-1_IntKi)
      end do
      do j=p%Nst+1_IntKi,p%Nst*2_IntKi
         p%theta_st(i,j) = p%theta_st(i,j-p%Nst) + pi
      end do
   end do
   
end subroutine DMST_SetParameters
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DMST_AllocInput( u, p, errStat, errMsg )
! This routine is called from DMST_Init.
!..................................................................................................................................
   type(DMST_InputType),           intent(  out)  :: u           ! Input data
   type(DMST_ParameterType),       intent(in   )  :: p           ! Parameters
   integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
   character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables
   integer(IntKi)                                 :: errStat2    ! Temporary error status of the operation
   character(*), parameter                        :: RoutineName = 'DMST_AllocInput'

      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""

   allocate ( u%Vinf(3_IntKi, p%numBladeNodes), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%Vinf.', errStat, errMsg, RoutineName )
      return
   end if 
   u%Vinf = 0.0_ReKi
   
   allocate ( u%pitch(p%numBladeNodes), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%pitch.', errStat, errMsg, RoutineName )
      return
   end if 
   u%pitch = 0.0_ReKi
   
end subroutine DMST_AllocInput
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DMST_AllocOutput( y, p, errStat, errMsg )
! This routine is called from DMST_Init.
!..................................................................................................................................
   type(DMST_OutputType),          intent(  out)  :: y           ! Output data
   type(DMST_ParameterType),       intent(in   )  :: p           ! Parameters
   integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
   character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables
   character(ErrMsgLen)                           :: errMsg2     ! Temporary error message if ErrStat /= ErrID_None
   integer(IntKi)                                 :: errStat2    ! Temporary error status of the operation
   character(*), parameter                        :: RoutineName = 'DMST_AllocOutput'
   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""

   call allocAry( y%Vind, p%numBladeNodes, p%numBlades, 'y%Vind', errStat2, errMsg2 ); call setErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

   y%Vind = 0.0_ReKi

end subroutine DMST_AllocOutput
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DMST_Init( InitInp, u, p, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
!..................................................................................................................................
   type(DMST_InitInputType),       intent(in   )  :: InitInp     ! Input for initialization routine
   type(DMST_InputType),           intent(  out)  :: u           ! An initial guess for the input; input mesh must be defined
   type(DMST_ParameterType),       intent(  out)  :: p           ! Parameters
   type(DMST_OutputType),          intent(  out)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                 !   only the output mesh is initialized)
   real(DbKi),                     intent(in   )  :: Interval    ! Coupling interval in seconds
   type(DMST_InitOutputType),      intent(  out)  :: InitOut     ! Output for initialization routine
   integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
   character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables
   character(ErrMsgLen)                           :: errMsg2     ! Temporary error message if ErrStat /= ErrID_None
   integer(IntKi)                                 :: errStat2    ! Temporary error status of the operation
   character(*), parameter                        :: RoutineName = 'DMST_Init'

      ! Initialize variables for this routine
   errStat = ErrID_None
   errMsg  = ""

      ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )

      !............................................................................................................................
      ! Define parameters here
      !............................................................................................................................
       
   call DMST_SetParameters( InitInp, p, errStat, errMsg )
      if (errStat >= AbortErrLev) return
   
   p%DT = Interval

      ! allocate all the arrays that store data in the input type:
   call DMST_AllocInput( u, p, errStat2, errMsg2 )      
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat >= AbortErrLev) return

      ! allocate all the arrays that store data in the output type:
   call DMST_AllocOutput(y, p, errStat2, errMsg2)
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat >= AbortErrLev) return
   
   InitOut%Version = DMST_Ver

END SUBROUTINE DMST_Init
! !----------------------------------------------------------------------------------------------------------------------------------
! subroutine BEMT_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! ! This routine is called at the end of the simulation.
! !..................................................................................................................................
!       TYPE(BEMT_InputType),           INTENT(INOUT)  :: u           ! System inputs
!       TYPE(BEMT_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
!       TYPE(BEMT_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
!       TYPE(BEMT_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
!       TYPE(BEMT_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
!       TYPE(BEMT_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other states
!       TYPE(BEMT_OutputType),          INTENT(INOUT)  :: y           ! System outputs
!       INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
!       CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

!          ! Initialize ErrStat

!       ErrStat = ErrID_None
!       ErrMsg  = ""

!          ! Place any last minute operations or calculations here:

!          ! Close files here:
!       if ( p%UA_Flag ) then
!          CALL UA_End(p%UA)
!       end if

!          ! Destroy the input data:

!       CALL BEMT_DestroyInput( u, ErrStat, ErrMsg )
      
!          ! Destroy the parameter data:

!       CALL BEMT_DestroyParam( p, ErrStat, ErrMsg )

!          ! Destroy the state data:

!       CALL BEMT_DestroyContState(   x,           ErrStat, ErrMsg )
!       CALL BEMT_DestroyDiscState(   xd,          ErrStat, ErrMsg )
!       CALL BEMT_DestroyConstrState( z,           ErrStat, ErrMsg )
!       CALL BEMT_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )

!          ! Destroy the output data:

!       CALL BEMT_DestroyOutput( y, ErrStat, ErrMsg )

! END SUBROUTINE BEMT_End
! !----------------------------------------------------------------------------------------------------------------------------------
! subroutine BEMT_UpdateStates( t, n, u1, u2, p, x, xd, z, OtherState, AFInfo, m, errStat, errMsg )
! ! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! ! Continuous, constraint, discrete, and other states are updated for t + Interval
! ! NOTE:  This is a non-standard framework interface!!!!!  GJH
! !..................................................................................................................................
!    real(DbKi),                          intent(in   ) :: t          ! Current simulation time in seconds
!    integer(IntKi),                      intent(in   ) :: n          ! Current simulation time step n = 0,1,...
!    type(BEMT_InputType),                intent(in   ) :: u1,u2      ! Input at t and t+ dt 
!    type(BEMT_ParameterType),            intent(in   ) :: p          ! Parameters   
!    type(BEMT_ContinuousStateType),      intent(inout) :: x          ! Input: Continuous states at t;
!                                                                     !   Output: Continuous states at t + Interval
!    type(BEMT_DiscreteStateType),        intent(inout) :: xd         ! Input: Discrete states at t;
!                                                                     !   Output: Discrete states at t  + Interval
!    type(BEMT_ConstraintStateType),      intent(inout) :: z          ! Input: Constraint states at t;
!                                                                     !   Output: Constraint states at t + Interval
!    type(BEMT_OtherStateType),           intent(inout) :: OtherState ! Input: Other states at t;
!                                                                     !   Output: Other states at t + Interval
!    type(BEMT_MiscVarType),              intent(inout) :: m          ! Misc/optimization variables
!    type(AFI_ParameterType),             intent(in   ) :: AFInfo(:)  ! The airfoil parameter data
!    integer(IntKi),                      intent(  out) :: errStat    ! Error status of the operation
!    character(*),                        intent(  out) :: errMsg     ! Error message if ErrStat /= ErrID_None
      
!    integer(IntKi)                                    :: i,j
   
!    character(ErrMsgLen)                              :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
!    integer(IntKi)                                    :: errStat2    ! temporary Error status of the operation
!    character(*), parameter                           :: RoutineName = 'BEMT_UpdateStates'
!    real(DbKi)                                        :: uTimes(2)
   
!    ErrStat = ErrID_None
!    ErrMsg = ""
   
!    uTimes(1) = t
!    uTimes(2) = t+p%dt
   
!    !...............................................................................................................................
!    ! if we haven't initialized z%phi, we want to get a better guess as to what the actual values of phi at t are:
!    !...............................................................................................................................
!    if (.not. OtherState%nodesInitialized) then
!       call UpdatePhi( u1, p, z%phi, AFInfo, m, OtherState%ValidPhi, errStat2, errMsg2 )
!       OtherState%nodesInitialized = .true.         ! otherState updated to t+dt (i.e., n+1)
!    end if
   
!    !...............................................................................................................................
!    !  compute inputs to DBEMT at step n (also setting inductions--including DBEMT and skewed wake corrections--at time n)
!    !...............................................................................................................................
!    call BEMT_CalcOutput_Inductions( 1, t, .true., .true., z%phi, u1, p, x, xd, z, OtherState, AFInfo, m%axInduction, m%tanInduction, m%chi, m, errStat2, errMsg2 )
!       call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!       if (ErrStat >= AbortErrLev) return

! #ifdef DEBUG_BEMT_RESIDUAL
!       if (p%useInduction) call WriteDEBUGValuesToFile(t, u1, p, x, xd, z, OtherState, m, AFInfo)
! #endif   

!    !...............................................................................................................................
!    !  compute inputs to UA at time n (also setting inductions--including DBEMT and skewed wake corrections--at time n)
!    !...............................................................................................................................
!    if (p%UA_Flag) then
!       m%phi = z%phi
!       call SetInputs_for_UA_AllNodes(u1, p, m%phi, m%axInduction, m%tanInduction, m%u_UA(:,:,1))
!    end if
   
!    !...............................................................................................................................
!    !  update BEMT states to step n+1
!    !...............................................................................................................................
!    call UpdatePhi( u2, p, z%phi, AFInfo, m, OtherState%ValidPhi, errStat2, errMsg2 )
!       call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!       if (errStat >= AbortErrLev) return
  
!    !...............................................................................................................................
!    !  compute inputs to DBEMT at step n+1 (also setting inductions--WITHOUT DBEMT or skewed wake corrections--at step n+1)
!    !...............................................................................................................................
!    call BEMT_CalcOutput_Inductions( 2, t, .true., .false., z%phi, u2, p, x, xd, z, OtherState, AFInfo, m%axInduction, m%tanInduction, m%chi, m, errStat2, errMsg2 )
!       call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!       if (ErrStat >= AbortErrLev) return

!    !...............................................................................................................................
!    !  update DBEMT states to step n+1
!    !...............................................................................................................................
!    if (p%DBEMT_Mod /= DBEMT_none) then

!       !............................................................................................................................
!       ! update DBEMT states to t+dt
!       !............................................................................................................................
!       do j = 1,p%numBlades
!          do i = 1,p%numBladeNodes
!             call DBEMT_UpdateStates(i, j, t, n, m%u_DBEMT, p%DBEMT, x%DBEMT, OtherState%DBEMT, m%DBEMT, errStat2, errMsg2)
!                if (ErrStat2 /= ErrID_None) then
!                   call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//trim(NodeText(i,j)))
!                   if (errStat >= AbortErrLev) return
!                end if
!          end do
!       end do
      
!    end if
   
!    !...............................................................................................................................
!    !  compute inputs to UA at time n+1 (also applying corrections to inductions--including DBEMT and skewed wake corrections)
!    !...............................................................................................................................
!    if (p%UA_Flag) then
!          ! after updating DBEMT states, we can now apply the corrections we omitted on the last call to BEMT_CalcOutput_Inductions()
!       if ( p%useInduction .and. .not. m%UseFrozenWake) then
!             !......................................................................................................................
!             ! apply DBEMT correction to axInduction and tanInduction:
!             !......................................................................................................................
!             if (p%DBEMT_Mod /= DBEMT_none) then
!                call calculate_Inductions_from_DBEMT_AllNodes(2, uTimes(2), u2, p, x, OtherState, m, m%axInduction, m%tanInduction)
!             end if
         
!             call ApplySkewedWakeCorrection_AllNodes(p, u2, m, m%axInduction, m%chi)

!             !......................................................................................................................
!             ! If TSR is too low, (start to) turn off induction
!             !......................................................................................................................
!             call check_turnOffBEMT(p, u2, m%BEM_weight, m%axInduction, m%tanInduction, m%FirstWarn_BEMoff)
            
!       end if
   
!       m%phi = z%phi
!       call SetInputs_for_UA_AllNodes(u2, p, m%phi, m%axInduction, m%tanInduction, m%u_UA(:,:,2))
   
!       !............................................................................................................................
!       !  compute UA states at t+dt
!       !............................................................................................................................
!       do j = 1,p%numBlades
!          do i = 1,p%numBladeNodes

!                ! COMPUTE: x%UA and/or xd%UA, OtherState%UA
!             call UA_UpdateStates( i, j, t, n, m%u_UA(i,j,:), uTimes, p%UA, x%UA, xd%UA, OtherState%UA, AFInfo(p%AFIndx(i,j)), m%UA, errStat2, errMsg2 )
!                if (ErrStat2 /= ErrID_None) then
!                   call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//trim(NodeText(i,j)))
!                   if (errStat >= AbortErrLev) return
!                end if

!          end do
!       end do

!    end if ! is UA used?
    
! end subroutine BEMT_UpdateStates
! !----------------------------------------------------------------------------------------------------------------------------------
! subroutine BEMT_CalcOutput( t, u, p, x, xd, z, OtherState, AFInfo, y, m, errStat, errMsg )
! ! Routine for computing outputs, used in both loose and tight coupling.
! !..................................................................................................................................
!    real(DbKi),                     intent(in   )  :: t           ! Current simulation time in seconds
!    type(BEMT_InputType),           intent(in   )  :: u           ! Inputs at Time t
!    type(BEMT_ParameterType),       intent(in   )  :: p           ! Parameters
!    type(BEMT_ContinuousStateType), intent(in   )  :: x           ! Continuous states at t
!    type(BEMT_DiscreteStateType),   intent(in   )  :: xd          ! Discrete states at t
!    type(BEMT_ConstraintStateType), intent(in   )  :: z           ! Constraint states at t
!    type(BEMT_OtherStateType),      intent(in   )  :: OtherState  ! Other states at t
!    type(BEMT_MiscVarType),         intent(inout)  :: m           ! Misc/optimization variables
!    type(AFI_ParameterType),        intent(in   )  :: AFInfo(:)   ! The airfoil parameter data
!    type(BEMT_OutputType),          intent(inout)  :: y           ! Outputs computed at t (Input only so that mesh con-
!                                                                  !   nectivity information does not have to be recalculated)
!    integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
!    character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None

!       ! Local variables:

!    integer(IntKi)                                 :: i                                               ! Generic index
!    integer(IntKi)                                 :: j                                               ! Loops through nodes / elements
!    integer(IntKi), parameter                      :: InputIndex=1      ! we will always use values at t in this routine
   
!    character(ErrMsgLen)                           :: errMsg2     ! temporary Error message if ErrStat /= ErrID_None
!    integer(IntKi)                                 :: errStat2    ! temporary Error status of the operation
!    character(*), parameter                        :: RoutineName = 'BEMT_CalcOutput'
   
!    type(AFI_OutputType)                           :: AFI_interp

!          ! Initialize some output values
!    errStat = ErrID_None
!    errMsg  = ""

!    y%phi = z%phi ! set this before possibly calling UpdatePhi() because phi is intent(inout) in the solve
   
!    !...............................................................................................................................
!    ! if we haven't initialized z%phi, we want to get a better guess as to what the actual values of phi are:
!    !...............................................................................................................................
!    if (.not. OtherState%nodesInitialized) then
!       call UpdatePhi( u, p, y%phi, AFInfo, m, m%ValidPhi, errStat2, errMsg2 )
!    end if
   
!    !...............................................................................................................................
!    ! calculate inductions using BEMT, applying the DBEMT, and/or skewed wake corrections as applicable:
!    !...............................................................................................................................
!    call BEMT_CalcOutput_Inductions( InputIndex, t, .false., .true., y%phi, u, p, x, xd, z, OtherState, AFInfo, y%axInduction, y%tanInduction, y%chi, m, errStat, errMsg )
   
!    !...............................................................................................................................
!    ! update phi if necessary (consistent with inductions) and calculate inputs to UA (EVEN if UA isn't used, because we use the inputs later):
!    !...............................................................................................................................
!    call SetInputs_for_UA_AllNodes(u, p, y%phi, y%axInduction, y%tanInduction, m%u_UA(:,:,InputIndex))
   
!    !...............................................................................................................................
!    ! unsteady aero and related outputs (cl, cd, cm, AoA, Vrel, Re)
!    !...............................................................................................................................
!    do j = 1,p%numBlades ! Loop through all blades
!       do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements

!             ! set output values:
!          y%AOA( i,j) = m%u_UA(i,j,InputIndex)%alpha
!          y%Vrel(i,j) = m%u_UA(i,j,InputIndex)%U
!          y%Re(  i,j) = m%u_UA(i,j,InputIndex)%Re
         
!       enddo             ! I - Blade nodes / elements
!    enddo          ! J - All blades
   
!       ! Now depending on the option for UA get the airfoil coefs, Cl, Cd, Cm for unsteady or steady implementation
!    if (p%UA_Flag ) then
   
!       do j = 1,p%numBlades ! Loop through all blades
!          do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements

!             call UA_CalcOutput(i, j, m%u_UA(i,j,InputIndex), p%UA, x%UA, xd%UA, OtherState%UA, AFInfo(p%AFindx(i,j)), m%y_UA, m%UA, errStat2, errMsg2 )
!                if (ErrStat2 /= ErrID_None) then
!                   call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//trim(NodeText(i,j)))
!                   if (errStat >= AbortErrLev) return
!                end if
               
!             y%Cl(i,j) = m%y_UA%Cl
!             y%Cd(i,j) = m%y_UA%Cd
!             y%Cm(i,j) = m%y_UA%Cm
!             y%Cpmin(i,j) = 0.0_ReKi !bjj: this isn't set anywhere... ???? 
!          enddo             ! I - Blade nodes / elements
!       enddo          ! J - All blades
   
!       ! if ( mod(REAL(t,ReKi),.1) < p%dt) then
!          call UA_WriteOutputToFile(t, p%UA, m%y_UA)
!       ! end if
      
!    else
!             ! compute steady Airfoil Coefs
!       do j = 1,p%numBlades ! Loop through all blades
!          do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements
         
!             call AFI_ComputeAirfoilCoefs( y%AOA(i,j), y%Re(i,j), u%UserProp(i,j), AFInfo(p%AFindx(i,j)), AFI_interp, errStat2, errMsg2 )
!                if (ErrStat2 /= ErrID_None) then
!                   call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//trim(NodeText(i,j)))
!                   if (errStat >= AbortErrLev) return
!                end if
!             y%Cl(i,j) = AFI_interp%Cl
!             y%Cd(i,j) = AFI_interp%Cd
!             y%Cm(i,j) = AFI_interp%Cm
!             y%Cpmin(i,j) = AFI_interp%Cpmin
         
!          enddo             ! I - Blade nodes / elements
!       enddo          ! J - All blades
      
!    end if

!    !...............................................................................................................................
!    ! Compute Cx, Cy given Cl, Cd and phi
!    !...............................................................................................................................
!    do j = 1,p%numBlades ! Loop through all blades
!       do i = 1,p%numBladeNodes ! Loop through the blade nodes / elements
!             ! NOTE: For these calculations we force the useAIDrag and useTIDrag flags to .TRUE.
!          call Transform_ClCd_to_CxCy( y%phi(i,j), .TRUE., .TRUE., y%Cl(i,j), y%Cd(i,j),y%Cx(i,j), y%Cy(i,j) )
            
!       enddo             ! I - Blade nodes / elements
!    enddo          ! J - All blades

!    return

! end subroutine BEMT_CalcOutput
! !----------------------------------------------------------------------------------------------------------------------------------
end module DMST