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
!   public :: DMST_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
   public :: DMST_CalcOutput                     ! Routine for computing outputs

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
   integer(IntKi)                                 :: lgth

      ! Initialize variables for this routine
   errStat = ErrID_None
   errMsg  = ""

   p%numBlades      = InitInp%numBlades   
   p%numBladeNodes  = InitInp%numBladeNodes 
   p%airDens        = InitInp%airDens
   p%kinVisc        = InitInp%kinVisc
   p%DMSTMod        = InitInp%DMSTMod
   p%Nst            = InitInp%Nst
   p%DMSTRes        = InitInp%DMSTRes

   allocate ( p%AFindx(p%numBladeNodes, p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%AFindx.', errStat, errMsg, RoutineName )
      return
   end if

   allocate ( p%chord(p%numBladeNodes, p%numBlades), STAT = errStat2 )
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

   allocate ( p%theta_st(2_IntKi*p%Nst, p%numBladeNodes), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%theta_st.', errStat, errMsg, RoutineName )
      return
   end if 

   lgth = floor(1.0_ReKi/p%DMSTRes)
   allocate ( p%indf(lgth), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%indf.', errStat, errMsg, RoutineName )
      return
   end if 

   p%AFindx = InitInp%AFindx 
   p%chord = InitInp%chord
   p%radius = InitInp%radius

      ! Compute the radius, total streamtube angle, and azimuthal position of streamtube midpoint
   do i=1,p%numBladeNodes
      p%dTheta(i) = pi/p%Nst
      do j=1,p%Nst
         p%theta_st(j,i) = p%dTheta(i)/2.0_ReKi + p%dTheta(i)*(j-1)
      end do
      do j=p%Nst+1,p%Nst*2
         p%theta_st(j,i) = p%theta_st(j-p%Nst,i) + pi
      end do
      do j = 1,lgth
         p%indf(j) = p%DMSTRes + p%DMSTRes*(j-1)
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

   allocate ( u%PitchAndTwist(p%numBladeNodes,p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%PitchAndTwist.', errStat, errMsg, RoutineName )
      return
   end if 
   u%PitchAndTwist = 0.0_ReKi

   allocate ( u%Vinf(3_IntKi, p%Nst, p%numBladeNodes), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%Vinf.', errStat, errMsg, RoutineName )
      return
   end if 
   u%Vinf = 0.0_ReKi

   allocate ( u%blade_st(p%numBladeNodes,p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%blade_st.', errStat, errMsg, RoutineName )
      return
   end if 
   u%blade_st = 0.0_ReKi

   allocate ( u%UserProp(p%numBladeNodes, p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%UserProp.', errStat, errMsg, RoutineName )
      return
   end if 
   u%UserProp = 0.0_ReKi
   
   u%omega  = 0.0_ReKi
   
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

   call allocAry( y%Vind, 3_IntKi, p%numBladeNodes, p%numBlades, 'y%Vind', errStat2, errMsg2 ); call setErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

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
   type(DMST_OutputType),          intent(  out)  :: y           ! Initial system outputs (outputs are not calculated; only the output mesh is initialized)
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

      ! Define parameters
   call DMST_SetParameters( InitInp, p, errStat, errMsg )
      if (errStat >= AbortErrLev) return
   
   p%DT = Interval

      ! Allocate all the arrays that store data in the input type
   call DMST_AllocInput( u, p, errStat2, errMsg2 )      
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat >= AbortErrLev) return

      ! Allocate all the arrays that store data in the output type
   call DMST_AllocOutput( y, p, errStat2, errMsg2 )
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
!----------------------------------------------------------------------------------------------------------------------------------
subroutine calculate_CTmo( DMSTMod, indf, CTmo )
! This routine is called from DMST_CalcOutput and calculates the thrust coefficient from linear momentum theory.
!..................................................................................................................................
   integer(IntKi),                 intent(in   )  :: DMSTMod     ! Type of momentum theory model
   real(ReKi),                     intent(in   )  :: indf(:)     ! Array of induction factors
   real(ReKi),                     intent(inout)  :: CTmo(:)     ! Thrust coefficient from linear momentum theory
   
      ! Local variables
   integer(IntKi)                                 :: i           ! Loops through induction factors

   if ( DMSTMod == 1 ) then
      do i = 1,size(indf)
         if ( indf(i) < 0.6_ReKi ) then
            CTmo(i) = 0.889_ReKi - ((0.0203_ReKi - (0.857_ReKi - indf(i))**2)/0.6427_ReKi) ! thrust coefficient from linear momentum theory (Glauert's empirical correction)
         else if ( indf(i) >= 0.6_ReKi ) then
            CTmo(i) = 4.0_ReKi*indf(i)*(1.0_ReKi - indf(i)) ! thrust coefficient from linear momentum theory (classic)
         end if
      end do
   else if ( DMSTMod == 2 ) then
      do i = 1,size(indf)
         if ( indf(i) < 0.3_ReKi ) then
            CTmo(i) = 0.889_ReKi - ((0.0203_ReKi - (0.857_ReKi - indf(i))**2)/0.6427_ReKi) ! thrust coefficient from linear momentum theory (Glauert's empirical correction)
         else if ( indf(i) >= 0.3_ReKi ) then
            CTmo(i) = 4.0_ReKi/3.0_ReKi*(1.0_ReKi - indf(i))*(2.0_ReKi + indf(i))/(2.0_ReKi - indf(i)) ! thrust coefficient from linear momentum theory (high load)
         end if
      end do
   end if
   
end subroutine calculate_CTmo
!----------------------------------------------------------------------------------------------------------------------------------
subroutine calculate_CTbe( m, Vinf, indf, p, u, AFinfo, CTbe )
   ! This routine is called from DMST_CalcOutput and calculates the thrust coefficient from blade element theory.
   !..................................................................................................................................
   integer(IntKi),                 intent(in   )                   :: m           ! Loops through sweeps
   real(ReKi),                     intent(in   )                   :: Vinf(:,:,:) ! Free-stream velocity
   real(ReKi),                     intent(in   )                   :: indf(:)     ! Array of induction factors
   type(DMST_ParameterType),       intent(in   )                   :: p           ! Parameters
   type(DMST_InputType),           intent(in   )                   :: u           ! Inputs at time t
   type(AFI_ParameterType),        intent(in   )                   :: AFInfo(:)   ! The airfoil parameter data
   real(ReKi),                     intent(inout)                   :: CTbe(:,:,:) ! Thrust coefficient from blade element theory
      
      ! Local variables
   real(ReKi),     dimension(size(p%indf),p%Nst,p%numBladeNodes)   :: V           ! Free-stream minus induced velocity
   real(ReKi),     dimension(size(p%indf),p%Nst,p%numBladeNodes)   :: lambda      ! Local tip-speed ratio
   real(ReKi),     dimension(size(p%indf),p%Nst,p%numBladeNodes)   :: Vrel        ! Relative velocity
   real(ReKi),     dimension(size(p%indf),p%Nst,p%numBladeNodes)   :: Reb         ! Blade Reynolds number
   real(ReKi),     dimension(size(p%indf),p%Nst,p%numBladeNodes)   :: alpha       ! Angle of attack
   real(ReKi),     dimension(size(p%indf),p%Nst,p%numBladeNodes)   :: phi         ! Inflow angle
   real(ReKi),     dimension(size(p%indf),p%Nst,p%numBladeNodes)   :: Cn          ! Normal force coefficient on the blade
   real(ReKi),     dimension(size(p%indf),p%Nst,p%numBladeNodes)   :: Ct          ! Tangential force coefficient on the blade
   integer(IntKi)                                                  :: i           ! Loops through induction factors
   integer(IntKi)                                                  :: j           ! Loops through streamtubes
   integer(IntKi)                                                  :: k           ! Loops through nodes
   integer(IntKi)                                                  :: n           ! Index for theta_st
   character(ErrMsgLen)                                            :: errMsg2     ! Temporary error message if ErrStat /= ErrID_None
   integer(IntKi)                                                  :: errStat2    ! Temporary error status of the operation
   type(AFI_OutputType)                                            :: AFI_interp  ! Interpolated airfoil coefficients

   do k = 1,p%numBladeNodes
      do j = 1,p%Nst
         n = j + (m-1)*(p%Nst)
         do i = 1,size(p%indf)
            V(i,j,k) = p%indf(i)*Vinf(1,j,k) ! free-stream minus induced velocity
            lambda(i,j,k) = u%omega*p%radius(k)/V(i,j,k) ! local tip-speed ratio
            Vrel(i,j,k) = V(i,j,k)*sqrt(1 + 2.0_ReKi*lambda(i,j,k)*cos(p%theta_st(n,k)) + lambda(i,j,k)**2) ! relative velocity
            Reb(i,j,k) = Vrel(i,j,k)*p%chord(k,1)/p%kinVisc ! blade Reynolds number
            alpha(i,j,k) = -atan2(sin(p%theta_st(n,k)),lambda(i,j,k) + cos(p%theta_st(n,k))) - u%PitchAndTwist(k,1) ! angle of attack
            call AFI_ComputeAirfoilCoefs( alpha(i,j,k), Reb(i,j,k), u%UserProp(k,1), AFInfo(p%AFindx(k,1)), AFI_interp, errStat2, errMsg2 ) ! outputs airfoil coefficients interpolated at given Reb and alpha 
            phi(i,j,k) = alpha(i,j,k) + u%PitchAndTwist(k,1) ! inflow angle
            Cn(i,j,k) = AFI_interp%Cd*sin(phi(i,j,k)) + AFI_interp%Cl*cos(phi(i,j,k)) ! normal force coefficient on the blade
            Ct(i,j,k) = AFI_interp%Cd*cos(phi(i,j,k)) - AFI_interp%Cl*sin(phi(i,j,k)) ! tangential force coefficient on the blade
            CTbe(i,j,k) = Ct(i,j,k)*cos(p%theta_st(n,k)) - Cn(i,j,k)*sin(p%theta_st(n,k)) ! thrust coefficient from blade element theory
         end do
      end do
   end do
      
end subroutine calculate_CTbe
!----------------------------------------------------------------------------------------------------------------------------------
subroutine calculate_Inductions_from_DMST( DMSTMod, indf, tol, crossPts, crossPtsSum, CTmo, CTbe, indf_final, indf_prev )
   ! This routine is called from DMST_CalcOutput and calculates the final induction factor in a streamtube.
   !..................................................................................................................................
   integer(IntKi),                 intent(in   )                   :: DMSTMod         ! Type of momentum theory model
   real(ReKi),                     intent(in   )                   :: indf(:)         ! Array of induction factors
   real(ReKi),                     intent(in   )                   :: tol             ! Tolerance for checking induction factor values
   integer(IntKi),                 intent(in   )                   :: crossPts(:)     ! Crossing points between CTmo and CTbe
   integer(IntKi),                 intent(in   )                   :: crossPtsSum     ! Number of crossing points in a streamtube
   real(ReKi),                     intent(in   )                   :: CTmo(:)         ! Thrust coefficient from linear momentum theory
   real(ReKi),                     intent(in   )                   :: CTbe(:)         ! Thrust coefficient from blade element theory
   real(ReKi),                     intent(inout)                   :: indf_final      ! Final induction factor in a streamtube
   real(ReKi),                     intent(inout)                   :: indf_prev       ! Final induction factor of the current/previous streamtube

      ! Local variables
   real(ReKi),     dimension(crossPtsSum)                          :: CTfinal         ! Array of final CT values
   real(ReKi),     dimension(crossPtsSum)                          :: indf_tmp        ! Temporary array of induction factors
   real(ReKi),     dimension(crossPtsSum)                          :: indf_diff       ! Different between indf_tmp and indf_prev
   integer(IntKi)                                                  :: i               ! Loops through induction factors
   integer(IntKi)                                                  :: m               ! Counter

      ! Initialize some local values
   CTfinal = 0.0
   indf_tmp = 9999.0
   m = 0

   do i = 1,size(indf)-1
      if ( crossPts(i) == 1_IntKi ) then
         if ( CTmo(i) >= CTbe(i) ) then
            m = m + 1
            CTfinal(m) = CTmo(i) + ( ( (CTmo(i) - CTbe(i)) / (CTbe(i+1) - CTbe(i)) ) / ( (1/(CTmo(i+1) - CTmo(i))) - (1/(CTbe(i+1) - CTbe(i))) ) )
            if ( DMSTMod == 1 ) then
               if ( indf(i) < 0.6_ReKi ) then
                  call DMST_QuadSolve_Glauert( tol, CTfinal(m), indf(i), indf(i+1), indf_tmp(m) )
               else if ( indf(i) >= 0.6_ReKi ) then
                  call DMST_QuadSolve_Theoretical( tol, CTfinal(m), indf(i), indf(i+1), indf_tmp(m) )
               end if
            else if ( DMSTMod == 2 ) then
               if ( indf(i) < 0.3_ReKi ) then
                  call DMST_QuadSolve_Glauert( tol, CTfinal(m), indf(i), indf(i+1), indf_tmp(m) )
               else if ( indf(i) >= 0.3_ReKi ) then
                  call DMST_QuadSolve_HighLoad( tol, CTfinal(m), indf(i), indf(i+1), indf_tmp(m) )
               end if
            end if
         else if ( CTmo(i) < CTbe(i) ) then
            m = m + 1
            CTfinal(m) = CTbe(i) + ( ( (CTbe(i) - CTmo(i)) / (CTmo(i+1) - CTmo(i)) ) / ( (1/(CTbe(i+1) - CTbe(i))) - (1/(CTmo(i+1) - CTmo(i))) ) )
            if ( DMSTMod == 1) then
               if ( indf(i) < 0.6_ReKi ) then
                  call DMST_QuadSolve_Glauert( tol, CTfinal(m), indf(i), indf(i+1), indf_tmp(m) )
               else if ( indf(i) >= 0.6_ReKi ) then
                  call DMST_QuadSolve_Theoretical( tol, CTfinal(m), indf(i), indf(i+1), indf_tmp(m) )
               end if
            else if ( DMSTMod == 2 ) then
               if ( indf(i) < 0.3_ReKi ) then
                  call DMST_QuadSolve_Glauert( tol, CTfinal(m), indf(i), indf(i+1), indf_tmp(m) )
               else if ( indf(i) >= 0.3_ReKi ) then
                  call DMST_QuadSolve_HighLoad( tol, CTfinal(m), indf(i), indf(i+1), indf_tmp(m) )
               end if
            end if
         end if
      end if
   end do

   if ( crossPtsSum > 1_IntKi ) then
      do i = 1,size(indf_tmp)
         indf_diff(i) = abs(indf_tmp(i) - indf_prev)
      end do
      indf_final = indf_tmp(minloc(indf_diff,1_IntKi))
   else if ( crossPtsSum == 1_IntKi ) then
      indf_final = indf_tmp(1)
   end if

   indf_prev = indf_final

end subroutine calculate_Inductions_from_DMST
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DMST_QuadSolve_Glauert( tol, CTfinal, indf, indf_plus, indf_tmp )
! This routine is called from calculate_Inductions_from_DMST and calculates induction factors given a thrust coefficient value.
!..................................................................................................................................
   real(ReKi),                     intent(in   )  :: tol          ! Tolerance for checking induction factor values
   real(ReKi),                     intent(in   )  :: CTfinal      ! A possible final CT value
   real(ReKi),                     intent(in   )  :: indf         ! An induction factor
   real(ReKi),                     intent(in   )  :: indf_plus    ! An induction factor
   real(ReKi),                     intent(inout)  :: indf_tmp     ! A possible final induction factor

      ! Local variables
   real(ReKi)                                     :: discriminant ! Discriminant of quadratic solution
   real(ReKi)                                     :: indf_tmp1    ! Temporary induction factor value
   real(ReKi)                                     :: indf_tmp2    ! Temporary induction factor value

   discriminant = 2.6669_ReKi**2 - 4.0_ReKi*1.5559_ReKi*(2.0001_ReKi - CTfinal)
   if ( discriminant >= 0.0_ReKi ) then
      indf_tmp1 = (2.6669_ReKi + sqrt(discriminant))/(2.0_ReKi*1.5559_ReKi)
      indf_tmp2 = (2.6669_ReKi - sqrt(discriminant))/(2.0_ReKi*1.5559_ReKi)
      if ( indf_tmp1 >= indf-tol .and. indf_tmp1 <= indf_plus+tol ) then
         indf_tmp = indf_tmp1
      else if ( indf_tmp2 >= indf-tol .and. indf_tmp2 <= indf_plus+tol ) then
         indf_tmp = indf_tmp2
      end if
   end if
      
end subroutine DMST_QuadSolve_Glauert
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DMST_QuadSolve_Theoretical( tol, CTfinal, indf, indf_plus, indf_tmp )
! This routine is called from calculate_Inductions_from_DMST and calculates induction factors given a thrust coefficient value.
!..................................................................................................................................
   real(ReKi),                     intent(in   )  :: tol          ! Tolerance for checking induction factor values
   real(ReKi),                     intent(in   )  :: CTfinal      ! A possible final CT value
   real(ReKi),                     intent(in   )  :: indf         ! An induction factor
   real(ReKi),                     intent(in   )  :: indf_plus    ! An induction factor
   real(ReKi),                     intent(inout)  :: indf_tmp     ! A possible final induction factor
   
      ! Local variables
   real(ReKi)                                     :: discriminant ! Discriminant of quadratic solution
   real(ReKi)                                     :: indf_tmp1    ! Temporary induction factor value
   real(ReKi)                                     :: indf_tmp2    ! Temporary induction factor value

   discriminant = 4.0_ReKi**2 - 4.0_ReKi*4.0_ReKi*CTfinal
   if ( discriminant >= 0.0_ReKi ) then
      indf_tmp1 = (4.0_ReKi + sqrt(discriminant))/(2.0_ReKi*4.0_ReKi)
      indf_tmp2 = (4.0_ReKi - sqrt(discriminant))/(2.0_ReKi*4.0_ReKi)
      if ( indf_tmp1 >= indf-tol .and. indf_tmp1 <= indf_plus+tol ) then
         indf_tmp = indf_tmp1
      else if ( indf_tmp2 >= indf-tol .and. indf_tmp2 <= indf_plus+tol ) then
         indf_tmp = indf_tmp2
      end if
   end if
         
end subroutine DMST_QuadSolve_Theoretical
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DMST_QuadSolve_HighLoad( tol, CTfinal, indf, indf_plus, indf_tmp )
   ! This routine is called from calculate_Inductions_from_DMST and calculates induction factors given a thrust coefficient value.
   !..................................................................................................................................
      real(ReKi),                     intent(in   )  :: tol          ! Tolerance for checking induction factor values
      real(ReKi),                     intent(in   )  :: CTfinal      ! A possible final CT value
      real(ReKi),                     intent(in   )  :: indf         ! An induction factor
      real(ReKi),                     intent(in   )  :: indf_plus    ! An induction factor
      real(ReKi),                     intent(inout)  :: indf_tmp     ! A possible final induction factor
      
         ! Local variables
      real(ReKi)                                     :: discriminant ! Discriminant of quadratic solution
      real(ReKi)                                     :: indf_tmp1    ! Temporary induction factor value
      real(ReKi)                                     :: indf_tmp2    ! Temporary induction factor value
   
      discriminant = (1.0_ReKi - 3.0_ReKi/4.0_ReKi*CTfinal)**2 - 4.0_ReKi*1.0_ReKi*(3.0_ReKi/2.0_ReKi*CTfinal - 2.0_ReKi)
      if ( discriminant >= 0.0_ReKi ) then
         indf_tmp1 = (3.0_ReKi/4.0_ReKi*CTfinal - 1.0_ReKi + sqrt(discriminant))/(2.0_ReKi)
         indf_tmp2 = (3.0_ReKi/4.0_ReKi*CTfinal - 1.0_ReKi - sqrt(discriminant))/(2.0_ReKi)
         if ( indf_tmp1 >= indf-tol .and. indf_tmp1 <= indf_plus+tol ) then
            indf_tmp = indf_tmp1
         else if ( indf_tmp2 >= indf-tol .and. indf_tmp2 <= indf_plus+tol ) then
            indf_tmp = indf_tmp2
         end if
      end if
            
   end subroutine DMST_QuadSolve_HighLoad
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DMST_CalcOutput( u, p, AFInfo, y, errStat, errMsg )
! Routine for computing outputs.
!..................................................................................................................................
   type(DMST_InputType),           intent(in   )                           :: u              ! Inputs at time t
   type(DMST_ParameterType),       intent(in   )                           :: p              ! Parameters
   type(AFI_ParameterType),        intent(in   )                           :: AFInfo(:)      ! The airfoil parameter data
   type(DMST_OutputType),          intent(inout)                           :: y              ! Outputs computed at t
   integer(IntKi),                 intent(  out)                           :: errStat        ! Error status of the operation
   character(*),                   intent(  out)                           :: errMsg         ! Error message if ErrStat /= ErrID_None

      ! Local variables
   real(ReKi),     dimension(3_IntKi,p%Nst,p%numBladeNodes)                :: Vinf           ! Free-stream velocity
   real(ReKi),     dimension(size(p%indf))                                 :: CTmo           ! Thrust coefficient from linear momentum theory
   real(ReKi),     dimension(size(p%indf),p%Nst,p%numBladeNodes)           :: CTbe           ! Thrust coefficient from blade element theory
   real(ReKi),     dimension(size(p%indf),p%Nst,p%numBladeNodes)           :: CTdiff         ! Difference between CTbe and CTmo
   integer(IntKi), dimension(size(p%indf)-1_IntKi,p%Nst,p%numBladeNodes)   :: crossPts       ! Crossing points between CTmo and CTbe
   integer(IntKi), dimension(p%Nst,p%numBladeNodes)                        :: crossPtsSum    ! Number of crossing points in a streamtube
   integer(IntKi), dimension(p%numBladeNodes)                              :: crossPtsInd    ! Index of the first streamtube with a single crossing point
   real(ReKi),     dimension(p%Nst,p%numBladeNodes)                        :: indf_final     ! Final induction factors
   real(ReKi),     dimension(p%Nst,p%numBladeNodes)                        :: indf_final_u   ! Final induction factors for the upstream sweep
   real(ReKi),     dimension(p%Nst,p%numBladeNodes)                        :: indf_final_d   ! Final induction factors for the downstream sweep
   real(ReKi)                                                              :: indf_prev      ! Final induction factor of the current/previous streamtube
   real(ReKi),     dimension(3_IntKi,p%Nst*2_IntKi,p%numBladeNodes)        :: Vind_st        ! Induced velocity in global coordinates
   integer(IntKi)                                                          :: i              ! Loops through induction factors
   integer(IntKi)                                                          :: j              ! Loops through streamtubes
   integer(IntKi)                                                          :: k              ! Loops through nodes
   integer(IntKi)                                                          :: m              ! Loops through sweeps
   integer(IntKi)                                                          :: n              ! Loops through blades   
   character(*), parameter                                                 :: RoutineName = 'DMST_CalcOutput'

      ! Initialize some output values
   errStat = ErrID_None
   errMsg  = ""

         ! Initialize some local values
   indf_final_u = 0.0
   indf_final_d = 0.0
   Vind_st = 0.0
   CTmo = 0.0

      ! Calculate the thrust coefficient from linear momentum theory
   call calculate_CTmo( p%DMSTMod, p%indf, CTmo )

      ! Loop through upstream and downstream sweeps
   do m = 1,2

         ! Initialize some local values
      CTbe = 0.0
      CTdiff = 0.0
      crossPts = 0
      crossPtsSum = 0
      crossPtsInd = 0
      indf_final = 0.0
      indf_prev = 0.0
      Vinf = 0.0

      if ( m == 1_IntKi ) then
         Vinf = u%Vinf
      else if ( m == 2_IntKi ) then
         do k = 1,p%numBladeNodes
            do j = 1,p%Nst
               if ( p%DMSTMod == 1 ) then
                  Vinf(1,j,k) = (2.0_ReKi*indf_final_u(p%Nst-j+1,k) - 1.0_ReKi)*u%Vinf(1,p%Nst-j+1,k)
               else if ( p%DMSTMod == 2 ) then
                  Vinf(1,j,k) = indf_final_u(p%Nst-j+1,k)/(2.0_ReKi - indf_final_u(p%Nst-j+1,k))*u%Vinf(1,p%Nst-j+1,k)
               end if
            end do 
         end do
      end if 

         ! Calculate the thrust coefficient from blade element theory
      call calculate_CTbe( m, Vinf, p%indf, p, u, AFInfo, CTbe )

         ! Calculate the difference between CTbe and CTmo
      do k = 1,p%numBladeNodes
         do j = 1,p%Nst
            CTdiff(:,j,k) = CTbe(:,j,k) - CTmo
         end do
      end do

         ! Locate crossing points between CTmo and CTbe
      do k = 1,p%numBladeNodes
         do j = 1,p%Nst
            do i = 1,size(p%indf)-1
               if ( CTdiff(i,j,k) < 0.0_ReKi .and. CTdiff(i+1,j,k) >= 0.0_ReKi ) then
                  crossPts(i,j,k) = 1_IntKi
               else if ( CTdiff(i,j,k) >= 0.0_ReKi .and. CTdiff(i+1,j,k) < 0.0_ReKi ) then
                  crossPts(i,j,k) = 1_IntKi
               end if
            end do
            ! Calculate the number of crossing points in a streamtube
            do i = 1,size(p%indf)-1
               crossPtsSum(j,k) = crossPtsSum(j,k) + crossPts(i,j,k)
            end do
         end do
      end do

         ! Identify first streamtube with a single crossing point
      do k = 1,p%numBladeNodes
         do j = 1,p%Nst
            if ( crossPtsSum(j,k) == 1_IntKi ) then
               crossPtsInd(k) = j
               exit
            end if
         end do 
      end do

         ! Calculate final thrust coefficients and induction factors
      do k = 1,p%numBladeNodes
         do j = crossPtsInd(k),p%Nst
            call calculate_Inductions_from_DMST( p%DMSTMod, p%indf, p%DMSTRes, crossPts(:,j,k), crossPtsSum(j,k), CTmo, CTbe(:,j,k), indf_final(j,k), indf_prev )
         end do
         do j = 1,crossPtsInd(k)-1
            call calculate_Inductions_from_DMST( p%DMSTMod, p%indf, p%DMSTRes, crossPts(:,j,k), crossPtsSum(j,k), CTmo, CTbe(:,j,k), indf_final(j,k), indf_prev )
         end do
      end do
      if ( m == 1_IntKi ) then
         indf_final_u = indf_final
      else if ( m == 2_IntKi ) then
         indf_final_d = indf_final
      end if

   end do

   ! Calculate final induced velocities in global coordinates
   do k = 1,p%numBladeNodes
      if ( p%DMSTMod == 1 ) then   
         do j = 1,p%Nst
            Vind_st(1,j,k) = (indf_final_u(j,k) - 1.0_ReKi)*u%Vinf(1,j,k)
         end do
         do j = p%Nst+1,p%Nst*2
            Vind_st(1,j,k) = (2.0_ReKi*indf_final_u(2*p%Nst-j+1,k)*indf_final_d(j-p%Nst,k) - indf_final_d(j-p%Nst,k) - indf_final_u(2*p%Nst-j+1,k))*u%Vinf(1,2*p%Nst-j+1,k) + Vind_st(1,2*p%Nst-j+1,k)
         end do
      else if ( p%DMSTMod == 2 ) then
         do j = 1,p%Nst
            Vind_st(1,j,k) = (indf_final_u(j,k) - 1.0_ReKi)*u%Vinf(1,j,k)
         end do
         do j = p%Nst+1,p%Nst*2
            Vind_st(1,j,k) = (indf_final_d(j-p%Nst,k)/(2.0_ReKi - indf_final_u(2*p%Nst-j+1,k)) - 1.0_ReKi)*indf_final_u(2*p%Nst-j+1,k)*u%Vinf(1,2*p%Nst-j+1,k) + Vind_st(1,2*p%Nst-j+1,k)
         end do
      end if
   end do

   ! Output induced velocity values at blade nodes
   do n = 1,p%numBlades
      do k = 1,p%numBladeNodes
         y%Vind(:,k,n) = Vind_st(:,u%blade_st(k,n),k)
      end do
   end do

end subroutine DMST_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
end module DMST
