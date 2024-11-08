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
   public :: DMST_ReInit                         ! Re-initialization routine
!   public :: DMST_End                            ! Ending routine (includes clean up)
   public :: DMST_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
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

   allocate ( p%radius(p%numBladeNodes,p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%radius.', errStat, errMsg, RoutineName )
      return
   end if 

   allocate ( p%dTheta(p%numBladeNodes), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%dTheta.', errStat, errMsg, RoutineName )
      return
   end if 

   allocate ( p%theta_st(2_IntKi*p%Nst), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%theta_st.', errStat, errMsg, RoutineName )
      return
   end if 

   lgth = floor(2.0_ReKi/p%DMSTRes) - 1_IntKi
   allocate ( p%indf(lgth), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for p%indf.', errStat, errMsg, RoutineName )
      return
   end if 

   p%AFindx = InitInp%AFindx 
   p%chord = InitInp%chord
   p%radius = InitInp%radius

      ! Compute the total streamtube angle and azimuthal position corresponding to each streamtube
   p%dTheta = pi/p%Nst
   do i=1,p%numBladeNodes
      do j=1,p%Nst
         p%theta_st(j) = p%dTheta(i)/2.0_ReKi + p%dTheta(i)*(j-1_IntKi)
      end do
      do j=p%Nst+1,p%Nst*2
         p%theta_st(j) = p%theta_st(j-p%Nst) + pi
      end do
   end do
   do j = 1,lgth
      p%indf(j) = p%DMSTRes + p%DMSTRes*(j-1)
   end do

end subroutine DMST_SetParameters
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DMST_InitOtherStates( OtherState, p, errStat, errMsg )
! This routine is called from DMST_Init.
! The OtherState data is allocated.
!..................................................................................................................................

   type(DMST_OtherStateType),      intent(inout)  :: OtherState  ! OtherState data
   type(DMST_ParameterType),       intent(in   )  :: p           ! Parameters
   integer(IntKi),                 intent(  out)  :: errStat     ! Error status of the operation
   character(*),                   intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables
   integer(IntKi)                                 :: errStat2     ! Temporary error status of the operation
   character(*), parameter                        :: RoutineName = 'DMST_InitOtherStates'
   
      ! Initialize variables for this routine
   errStat = ErrID_None
   errMsg  = ""
     
   allocate ( OtherState%Vstr( 3_IntKi, p%Nst, p%numBladeNodes ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for OtherState%Vstr.', errStat, errMsg, RoutineName )
      return
   end if

   allocate ( OtherState%M_ag( 3_IntKi, 3_IntKi, p%Nst, p%numBladeNodes ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for OtherState%M_ag.', errStat, errMsg, RoutineName )
      return
   end if

   allocate ( OtherState%blade_theta( p%Nst, p%numBladeNodes ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for OtherState%blade_theta.', errStat, errMsg, RoutineName )
      return
   end if

   allocate ( OtherState%indf( 2_IntKi*p%Nst, p%numBladeNodes ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for OtherState%indf.', errStat, errMsg, RoutineName )
      return
   end if
   
   ! Values of the OtherStates are initialized in DMST_ReInit()

end subroutine DMST_InitOtherStates
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

   allocate ( u%Vinf(3_IntKi,p%numBladeNodes,p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%Vinf.', errStat, errMsg, RoutineName )
      return
   end if 
   u%Vinf = 0.0_ReKi

   allocate ( u%blade_theta(p%numBladeNodes,p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%blade_theta.', errStat, errMsg, RoutineName )
      return
   end if 
   u%blade_theta = 0.0_ReKi

   allocate ( u%Vstr(3_IntKi,p%numBladeNodes,p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%Vstr.', errStat, errMsg, RoutineName )
      return
   end if 
   u%Vstr = 0.0_ReKi

   allocate ( u%M_ag(3_IntKi,3_IntKi,p%numBladeNodes,p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%M_ag.', errStat, errMsg, RoutineName )
      return
   end if 
   u%M_ag = 0.0_ReKi

   allocate ( u%blade_st(p%numBladeNodes,p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for u%blade_st.', errStat, errMsg, RoutineName )
      return
   end if 
   u%blade_st = 0_IntKi

   allocate ( u%UserProp(p%numBladeNodes,p%numBlades), STAT = errStat2 )
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
   call allocAry( y%indf, p%numBladeNodes, p%numBlades, 'y%indf', errStat2, errMsg2 ); call setErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

   y%Vind = 0.0_ReKi
   y%indf = 0.0_ReKi

end subroutine DMST_AllocOutput
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DMST_Init( InitInp, u, p, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
!..................................................................................................................................
   type(DMST_InitInputType),       intent(in   )  :: InitInp     ! Input for initialization routine
   type(DMST_InputType),           intent(  out)  :: u           ! An initial guess for the input; input mesh must be defined
   type(DMST_ParameterType),       intent(  out)  :: p           ! Parameters
   type(DMST_OtherStateType),      intent(  out)  :: OtherState  ! Initial other states
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

      ! Initialize other states
   call DMST_InitOtherStates( OtherState, p, errStat, errMsg )
      if (errStat >= AbortErrLev) return

      ! Allocate all the arrays that store data in the input type
   call DMST_AllocInput( u, p, errStat2, errMsg2 )      
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat >= AbortErrLev) return

      ! Allocate all the arrays that store data in the output type
   call DMST_AllocOutput( y, p, errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat >= AbortErrLev) return
   
   InitOut%Version = DMST_Ver

      ! Set initial values for states
   call DMST_ReInit( p, OtherState, errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat >= AbortErrLev) return

end subroutine DMST_Init
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DMST_ReInit(p,OtherState,ErrStat,ErrMsg)

   type(DMST_ParameterType),       intent(in   )  :: p           ! Parameters
   type(DMST_OtherStateType),      intent(inout)  :: OtherState  ! Initial other states
   INTEGER(IntKi),                 intent(  out)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   character(*), parameter                        :: RoutineName = 'DMST_ReInit'

   ErrStat = ErrID_None
   ErrMsg  = ""
   
  
   OtherState%Vstr = 9999.0_ReKi
   OtherState%M_ag = 0.0_ReKi
   OtherState%blade_theta = 0.0_ReKi
   OtherState%indf = 9999.0_ReKi
    
end subroutine DMST_ReInit
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
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DMST_UpdateStates(p, u, y, OtherState, errStat, errMsg)
!
!..................................................................................................................................
   type(DMST_ParameterType),            intent(in   ) :: p          ! Parameters
   type(DMST_InputType),                intent(in   ) :: u(2)       ! Inputs at t and t + dt
   type(DMST_OutputType),               intent(in   ) :: y          ! Outputs at t
   type(DMST_OtherStateType),           intent(inout) :: OtherState ! Input: Other states at t; Output: Other states at t + Interval
   integer(IntKi),                      intent(  out) :: errStat    ! Error status of the operation
   character(*),                        intent(  out) :: errMsg     ! Error message if ErrStat /= ErrID_None

   integer(IntKi)                                     :: j
   integer(IntKi)                                     :: k
      
   ErrStat = ErrID_None
   ErrMsg = ""

   do k = 1,p%numBlades
      do j = 1,p%numBladeNodes
         if (u(2)%blade_st(j,k) <= p%Nst) then
            OtherState%Vstr(:,u(2)%blade_st(j,k),j) = u(2)%Vstr(:,j,k) ! u at t+dt
            OtherState%M_ag(:,:,u(2)%blade_st(j,k),j) = u(2)%M_ag(:,:,j,k) ! u at t+dt
            OtherState%blade_theta(u(2)%blade_st(j,k),j) = u(2)%blade_theta(j,k) ! u at t+dt
         end if
         OtherState%indf(u(1)%blade_st(j,k),j) = y%indf(j,k) ! u,y at t
      end do
   end do

   !do i=1,size(u)
   !   call DMST_DestroyInput( u(i), ErrStat, ErrMsg )
   !enddo

end subroutine DMST_UpdateStates
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
subroutine calculate_CTbe( j, k, Vinf, blade_theta, M_ag, Vstr, p, u, AFinfo, CTbe )
   ! This routine is called from DMST_CalcOutput and calculates the thrust coefficient from blade element theory.
   !..................................................................................................................................
   integer(IntKi),                 intent(in   )                   :: j                ! Blade node number
   integer(IntKi),                 intent(in   )                   :: k                ! Blade number
   real(ReKi),                     intent(in   )                   :: Vinf(:)          ! Free-stream velocity, global coordinates
   real(ReKi),                     intent(in   )                   :: blade_theta      ! Azimuthal angle of each blade node
   real(ReKi),                     intent(in   )                   :: M_ag(:,:)        ! Orientation matrix
   real(ReKi),                     intent(in   )                   :: Vstr(:)          ! Structural velocity, airfoil coordinates
   type(DMST_ParameterType),       intent(in   )                   :: p                ! Parameters
   type(DMST_InputType),           intent(in   )                   :: u                ! Inputs at time t
   type(AFI_ParameterType),        intent(in   )                   :: AFInfo(:)        ! The airfoil parameter data
   real(ReKi),                     intent(inout)                   :: CTbe(:)          ! Thrust coefficient from blade element theory
      
      ! Local variables
   real(ReKi),                     dimension(3)                    :: V           ! Free-stream minus induced velocity, global coordinates
   real(ReKi),                     dimension(3)                    :: Vtot_a      ! Free-stream minus induced velocity minus structural velocity, airfoil coordinates
   real(ReKi)                                                      :: Vrel        ! Relative velocity
   real(ReKi)                                                      :: Reb         ! Blade Reynolds number
   real(ReKi)                                                      :: alpha       ! Angle of attack
   real(ReKi)                                                      :: Cn          ! Normal force coefficient on the blade
   real(ReKi)                                                      :: Ct          ! Tangential force coefficient on the blade
   integer(IntKi)                                                  :: i           ! Loops through induction factors
   character(ErrMsgLen)                                            :: errMsg2     ! Temporary error message if ErrStat /= ErrID_None
   integer(IntKi)                                                  :: errStat2    ! Temporary error status of the operation
   type(AFI_OutputType)                                            :: AFI_interp  ! Interpolated airfoil coefficients

      ! Initialize some local values
   V = 0.0
   Vtot_a = 0.0
   alpha = 0.0
   Vrel = 0.0
   Reb = 0.0
   Cn = 0.0
   Ct = 0.0

   if ( Vstr(1) < 9999.0 ) then
      do i = 1,size(p%indf)
         V(1) = p%indf(i)*Vinf(1) ! global coordinates
         V(2:3) = Vinf(2:3) ! global coordinates
         Vtot_a = matmul(M_ag,V) - Vstr ! airfoil coordinates
         alpha = atan2(Vtot_a(1),Vtot_a(2)) ! angle of attack
         Vrel = sqrt(Vtot_a(1)**2 + Vtot_a(2)**2) ! relative velocity
         Reb = Vrel*p%chord(j,k)/p%kinVisc ! blade Reynolds number
         call AFI_ComputeAirfoilCoefs( alpha, Reb, u%UserProp(j,k), AFInfo(p%AFindx(j,k)), AFI_interp, errStat2, errMsg2 ) ! outputs airfoil coefficients interpolated at given Reb and alpha 
         Cn = AFI_interp%Cd*sin(alpha) + AFI_interp%Cl*cos(alpha) ! normal to chord force coefficient on the blade, positive inwards towards rotation axis
         Ct = AFI_interp%Cd*cos(alpha) - AFI_interp%Cl*sin(alpha) ! tangential to chord force coefficient on the blade, positive towards trailing edge
         CTbe(i) = Ct*cos(blade_theta) - Cn*sin(blade_theta) ! thrust coefficient from blade element theory
      end do
   end if 

end subroutine calculate_CTbe
!----------------------------------------------------------------------------------------------------------------------------------
subroutine calculate_Inductions_from_DMST( DMSTMod, indf, tol, Nst, indf_hist, blade_st, crossPts, crossPtsSum, CTmo, CTbe, indf_final_tmp, indf_final_tmp_store )
   ! This routine is called from DMST_CalcOutput and calculates the final induction factor in a streamtube.
   !..................................................................................................................................
   integer(IntKi),                 intent(in   )                  :: DMSTMod              ! Type of momentum theory model
   real(ReKi),                     intent(in   )                  :: indf(:)              ! Array of induction factors
   real(ReKi),                     intent(in   )                  :: tol                  ! Tolerance for checking induction factor values
   integer(IntKi),                 intent(in   )                  :: Nst                  ! Number of streamtubes   
   real(ReKi),                     intent(in   )                  :: indf_hist(:)         ! Induction factors stored per streamtube
   integer(IntKi),                 intent(in   )                  :: blade_st             ! Streamtube of the current blade node
   integer(IntKi),                 intent(in   )                  :: crossPts(:)          ! Crossing points between CTmo and CTbe
   integer(IntKi),                 intent(in   )                  :: crossPtsSum          ! Number of crossing points for a blade node
   real(ReKi),                     intent(in   )                  :: CTmo(:)              ! Thrust coefficient from linear momentum theory
   real(ReKi),                     intent(in   )                  :: CTbe(:)              ! Thrust coefficient from blade element theory
   real(ReKi),                     intent(inout)                  :: indf_final_tmp       ! Final induction factor for a blade node
   real(ReKi),                     intent(inout)                  :: indf_final_tmp_store ! Stored final induction factor for a blade node

      ! Local variables
   real(ReKi),     dimension(crossPtsSum)                         :: CTfinal          ! Array of final CT values
   real(ReKi),     dimension(crossPtsSum)                         :: indf_tmp         ! Temporary array of induction factors
   real(ReKi),     dimension(crossPtsSum)                         :: indf_diff        ! Difference between indf_tmp and indf_hist
   integer(IntKi)                                                 :: blade_st_ind     ! Index of streamtube just before streamtube of current blade node
   integer(IntKi)                                                 :: i                ! Loops through induction factors
   integer(IntKi)                                                 :: m                ! Counter

      ! Initialize some local values
   CTfinal = 0.0
   indf_tmp = 1.0
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

   if ( blade_st > 1_IntKi ) then
      blade_st_ind = blade_st - 1_IntKi
   elseif ( blade_st == 1_IntKi ) then
      blade_st_ind = 2_IntKi*Nst
   end if

   if ( crossPtsSum > 1_IntKi .and. indf_hist(blade_st_ind) < 9999.0_ReKi ) then
         do i = 1,size(indf_tmp)
            indf_diff(i) = abs(indf_tmp(i) - indf_hist(blade_st_ind))
         end do
         indf_final_tmp = indf_tmp(minloc(indf_diff,1_IntKi))
         indf_final_tmp_store = indf_final_tmp
   else if ( crossPtsSum == 1_IntKi ) then
      indf_final_tmp = indf_tmp(1)
      indf_final_tmp_store = indf_final_tmp
   else
      indf_final_tmp = 1.0_ReKi
      indf_final_tmp_store = 9999.0_ReKi
   end if

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
      else
         indf_tmp = (indf+indf_plus)/2
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
      else
         indf_tmp = (indf+indf_plus)/2
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
      else
         indf_tmp = (indf+indf_plus)/2
      end if
   end if
         
end subroutine DMST_QuadSolve_HighLoad
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DMST_CalcOutput( u, p, OtherState, AFInfo, y, errStat, errMsg )
! Routine for computing outputs.
!..................................................................................................................................
   type(DMST_InputType),           intent(in   )                           :: u              ! Inputs at time t
   type(DMST_ParameterType),       intent(in   )                           :: p              ! Parameters
   type(DMST_OtherStateType),      intent(in   )                           :: OtherState     ! OtherState data
   type(AFI_ParameterType),        intent(in   )                           :: AFInfo(:)      ! The airfoil parameter data
   type(DMST_OutputType),          intent(inout)                           :: y              ! Outputs computed at t
   integer(IntKi),                 intent(  out)                           :: errStat        ! Error status of the operation
   character(*),                   intent(  out)                           :: errMsg         ! Error message if ErrStat /= ErrID_None

      ! Local variables
   real(ReKi)                                                                   :: blade_theta          ! Azimuthal location of each blade node
   real(ReKi),     dimension(3_IntKi,3_IntKi)                                   :: M_ag                 ! Blade orientation matrix
   real(ReKi),     dimension(3_IntKi)                                           :: Vstr                 ! Structural velocity, airfoil coordinates
   real(ReKi),     dimension(3_IntKi)                                           :: Vinf                 ! Inflow velocity, global coordinates
   real(ReKi),     dimension(size(p%indf))                                      :: CTmo                 ! Thrust coefficient from linear momentum theory
   real(ReKi),     dimension(size(p%indf))                                      :: CTbe                 ! Thrust coefficient from blade element theory
   real(ReKi),     dimension(size(p%indf))                                      :: CTdiff               ! Difference between CTbe and CTmo
   integer(IntKi), dimension(size(p%indf)-1_IntKi)                              :: crossPts             ! Crossing points between CTmo and CTbe
   integer(IntKi)                                                               :: crossPtsSum          ! Number of crossing points for a blade node
   real(ReKi)                                                                   :: indf_final_tmp       ! Temporary storage for final induction factors
   real(ReKi)                                                                   :: indf_final_tmp_store ! Temporary storage for stored final induction factors
   real(ReKi),     dimension(p%numBladeNodes,p%numBlades)                       :: indf_final_store     ! Final stored induction factors   
   real(ReKi),     dimension(p%numBladeNodes,p%numBlades)                       :: indf_final_u         ! Final induction factors for the upstream sweep
   real(ReKi),     dimension(p%numBladeNodes,p%numBlades)                       :: indf_final           ! Final induction factors
   real(ReKi),     dimension(3_IntKi,p%numBladeNodes,p%numBlades)               :: Vind_st              ! Induced velocity, global coordinates
   integer(IntKi)                                                               :: i                    ! Loops through induction factors
   integer(IntKi)                                                               :: j                    ! Loops through nodes
   integer(IntKi)                                                               :: k                    ! Loops through blades
   integer(IntKi)                                                               :: m                    ! Loops through sweeps
   character(*), parameter                                                      :: RoutineName = 'DMST_CalcOutput'

      ! Initialize some output values
   errStat = ErrID_None
   errMsg  = ""

      ! Initialize some local values
   CTmo = 0.0
   indf_final_u = 1.0
   indf_final = 0.0
   Vind_st = 0.0

      ! Calculate the thrust coefficient from linear momentum theory
   call calculate_CTmo( p%DMSTMod, p%indf, CTmo )

      ! Loop through upstream and downstream sweeps
   do m = 1,2

      do k = 1,p%numBlades
         do j = 1,p%numBladeNodes

               ! Initialize some local values
            blade_theta = 0.0
            M_ag = 0.0
            Vstr = 0.0
            Vinf = 0.0
            CTbe = 0.0
            CTdiff = 0.0
            crossPts = 0.0
            crossPtsSum = 0.0
            indf_final_tmp = 0.0

            if ( m == 1_IntKi .or. m == 2_IntKi .and. u%blade_st(j,k) > p%Nst ) then

                  ! Set azimuth, orientation matrix, structural velocity
               if ( m == 1_IntKi .and. u%blade_st(j,k) > p%Nst ) then
                  blade_theta = OtherState%blade_theta(p%Nst*2_IntKi+1_IntKi-u%blade_st(j,k),j) ! stored upstream blade azimuth
                  M_ag = OtherState%M_ag(:,:,p%Nst*2_IntKi+1_IntKi-u%blade_st(j,k),j) ! stored upstream blade orientation
                  Vstr = OtherState%Vstr(:,p%Nst*2_IntKi+1_IntKi-u%blade_st(j,k),j) ! stored upstream structural velocity, airfoil coordinates
               else
                  blade_theta = u%blade_theta(j,k)
                  M_ag = u%M_ag(:,:,j,k)
                  Vstr = u%Vstr(:,j,k)
               end if

                  ! Set inflow velocity, accounting for induction from upstream blades
               if ( m == 1_IntKi ) then
                  Vinf = u%Vinf(:,j,k)
               else
                  if ( p%DMSTMod == 1 ) then
                     Vinf(1) = (2.0_ReKi*indf_final_u(j,k) - 1.0_ReKi)*u%Vinf(1,j,k)
                  else if ( p%DMSTMod == 2 ) then
                     Vinf(1) = indf_final_u(j,k)/(2.0_ReKi - indf_final_u(j,k))*u%Vinf(1,j,k)
                  end if
                  Vinf(2:3) = u%Vinf(2:3,j,k)
               end if

                  ! Calculate the thrust coefficient from blade element theory
               call calculate_CTbe( j, k, Vinf, blade_theta, M_ag, Vstr, p, u, AFInfo, CTbe )

               if ( Vstr(1) < 9999.0 ) then
                     ! Calculate difference between CTbe and CTmo
                  CTdiff = CTbe - CTmo

                     ! Locate crossing points between CTmo and CTbe
                  do i = 1,size(p%indf)-1
                     if ( CTdiff(i) < 0.0_ReKi .and. CTdiff(i+1) >= 0.0_ReKi ) then
                        crossPts(i) = 1_IntKi
                     else if ( CTdiff(i) >= 0.0_ReKi .and. CTdiff(i+1) < 0.0_ReKi ) then
                        crossPts(i) = 1_IntKi
                     end if
                  end do

                     ! Calculate the number of crossing points for a blade node
                  do i = 1,size(p%indf)-1
                     crossPtsSum = crossPtsSum + crossPts(i)
                  end do

                     ! Calculate final thrust coefficients and induction factors
                  call calculate_Inductions_from_DMST( p%DMSTMod, p%indf, p%DMSTRes, p%Nst, OtherState%indf(:,j), u%blade_st(j,k), crossPts, crossPtsSum, CTmo, CTbe, indf_final_tmp, indf_final_tmp_store )

               else

                  indf_final_tmp = 1.0_ReKi
                  indf_final_tmp_store = 9999.0_ReKi

               end if

                  ! Store final induction factors
               if ( m == 1_IntKi .and. u%blade_st(j,k) > p%Nst ) then
                  indf_final_u(j,k) = indf_final_tmp
               else
                  indf_final(j,k) = indf_final_tmp
                  indf_final_store(j,k) = indf_final_tmp_store
               end if

            end if
         end do
      end do
   end do

      ! Calculate final induced velocities in global coordinates
   do k = 1,p%numBlades
      do j = 1,p%numBladeNodes
         if ( p%DMSTMod == 1 .and. u%blade_st(j,k) <= p%Nst ) then   
            Vind_st(1,j,k) = (indf_final(j,k) - 1.0_ReKi)*u%Vinf(1,j,k)
         elseif ( p%DMSTMod == 1 .and. u%blade_st(j,k) > p%Nst ) then
            Vind_st(1,j,k) = (2.0_ReKi*indf_final_u(j,k)*indf_final(j,k) - indf_final(j,k) - 1)*u%Vinf(1,j,k)
         elseif ( p%DMSTMod == 2 .and. u%blade_st(j,k) <= p%Nst ) then
            Vind_st(1,j,k) = (indf_final(j,k) - 1.0_ReKi)*u%Vinf(1,j,k)
         elseif ( p%DMSTMod == 2 .and. u%blade_st(j,k) > p%Nst ) then
            Vind_st(1,j,k) = (indf_final_u(j,k)*indf_final(j,k)/(2.0_ReKi - indf_final_u(j,k)) - 1.0_ReKi)*u%Vinf(1,j,k)
         end if
      end do
   end do

      ! Output induced velocity values at blade nodes
   y%Vind = Vind_st
   y%indf = indf_final_store

end subroutine DMST_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
end module DMST
