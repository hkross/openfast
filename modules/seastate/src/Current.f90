!**********************************************************************************************************************************
! The Current and Current_Types modules make up a template for creating user-defined calculations in the FAST Modularization 
! Framework. Currents_Types will be auto-generated based on a description of the variables for the module.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012-2015  National Renewable Energy Laboratory
!
!    This file is part of Current.
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
!    
!**********************************************************************************************************************************
MODULE Current

   USE Current_Types   
   USE NWTC_Library
      
   IMPLICIT NONE
   
   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: Current_ProgDesc = ProgDesc( 'Current', '', '' )

   
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: Current_Init                           ! Initialization routine
      
CONTAINS

!=======================================================================
!JASON: MOVE THIS USER-DEFINED ROUTINE (UserCurrent) TO THE UserSubs.f90 OF HydroDyn WHEN THE PLATFORM LOADING FUNCTIONALITY HAS BEEN DOCUMENTED!!!!!
!> This is a dummy routine for holding the place of a user-specified
!! current profile.  Modify this code to create your own profile.
SUBROUTINE UserCurrent ( zi, EffWtrDpth, DirRoot, CurrVxi, CurrVyi )

      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(SiKi), INTENT(OUT)      :: CurrVxi                                         !< xi-component of the current velocity at elevation zi, m/s.
      REAL(SiKi), INTENT(OUT)      :: CurrVyi                                         !< yi-component of the current velocity at elevation zi, m/s.
      REAL(SiKi), INTENT(IN )      :: EffWtrDpth                                      !< Effective water depth ( EffWtrDpth       >  0 ), meters.
      REAL(SiKi), INTENT(IN )      :: zi                                              !< Elevation   (-EffWtrDpth <= zi <= 0 ), meters.
      
      CHARACTER(*), INTENT(IN )    :: DirRoot                                         !< The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



      CurrVxi = 0.0
      CurrVyi = 0.0



      RETURN
      
END SUBROUTINE UserCurrent
      
      
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine computes the x and y current components for a given water elevation
SUBROUTINE Calc_Current( InitInp, z, h , DirRoot, CurrVxi, CurrVyi )
!----------------------------------------------------------------------------------------------------------------------------------

         ! This routine is used to initialize the variables associated with
         ! current.



      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(SiKi),              INTENT(OUT) :: CurrVxi         !< xi-component of the current velocity at elevation z (m/s)
      REAL(SiKi),              INTENT(OUT) :: CurrVyi         !< yi-component of the current velocity at elevation z (m/s)

      REAL(SiKi),              INTENT(IN ) :: h               !< Effective water depth (meters)  This quantity must be positive-valued
      REAL(SiKi),              INTENT(IN ) :: z               !< Elevation relative to the mean sea level (meters)
      CHARACTER(*),            INTENT(IN ) :: DirRoot         !< The name of the root file including the full path to the current working directory.  
                                                              !! This may be useful if you want this routine to write a permanent record of what it does 
                                                              !! to be stored with the simulation results: the results should be stored in a file whose name 
                                                              !! (including path) is generated by appending any suitable extension to DirRoot.
      TYPE(Current_InitInputType), INTENT(IN ) :: InitInp     !< Initialization data for the current module 


         ! Local Variables:

      REAL(ReKi)                           :: CurrSSV         ! Magnitude of sub -surface current velocity at elevation z (m/s)
      REAL(ReKi)                           :: CurrNSV         ! Magnitude of near-surface current velocity at elevation z (m/s)



         ! If elevation z lies between the seabed and the mean sea level, compute the
         !   xi- and yi-components of the current (which depends on which current
         !   profile model is selected), else set CurrVxi and CurrVyi to zero:

      IF ( ( z < -h ) .OR. ( z > 0.0 ) )  THEN  ! .TRUE. if elevation z lies below the seabed or above mean sea level (exclusive)


            CurrVxi = 0.0  ! Set both the xi- and yi-direction
            CurrVyi = 0.0  ! current velocities to zero


      ELSE                                      ! Elevation z must lie between the seabed and the mean sea level (inclusive)


         SELECT CASE ( InitInp%CurrMod ) ! Which current profile model are we using?

         CASE ( 0 )              ! None!

            CurrVxi = 0.0  ! Set both the xi- and yi-direction
            CurrVyi = 0.0  ! current velocities to zero


         CASE ( 1 )              ! Standard (using inputs from PtfmFile).

            CurrSSV =      InitInp%CurrSSV0*( ( z + h                      )/h                      )**(1.0/7.0)
            CurrNSV = MAX( InitInp%CurrNSV0*( ( z + InitInp%CurrNSRef )/InitInp%CurrNSRef )           , 0.0_SiKi )

            CurrVxi = InitInp%CurrDIV*COS( D2R*InitInp%CurrDIDir ) + CurrSSV*COS( D2R*InitInp%CurrSSDir ) + &
                                   CurrNSV*COS( D2R*InitInp%CurrNSDir )

            CurrVyi = InitInp%CurrDIV*SIN( D2R*InitInp%CurrDIDir ) + CurrSSV*SIN( D2R*InitInp%CurrSSDir ) + &
                                   CurrNSV*SIN( D2R*InitInp%CurrNSDir )


         CASE ( 2 )              ! User-defined current profile model.

            CALL UserCurrent ( z, h, DirRoot, CurrVxi, CurrVyi )

         ENDSELECT

      END IF

   RETURN
   
END SUBROUTINE Calc_Current


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps. 
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE Current_Init( InitInp, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(Current_InitInputType),       INTENT(IN   )  :: InitInp     !< Input data for initialization routine
   TYPE(Current_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


     
      

      ! Local Variables:

   REAL(SiKi)                   :: CurrVxi                          ! xi-component of the current velocity at elevation z (m/s)     
   REAL(SiKi)                   :: CurrVyi                          ! yi-component of the current velocity at elevation z (m/s)
   REAL(SiKi)                   :: CurrVxi0                         ! xi-component of the current velocity at zi =  0.0 meters            (m/s)
   REAL(SiKi)                   :: CurrVyi0                         ! yi-component of the current velocity at zi =  0.0 meters            (m/s)
   REAL(SiKi)                   :: CurrVxiS                         ! xi-component of the current velocity at zi = -SmllNmbr meters       (m/s)
   REAL(SiKi)                   :: CurrVyiS                         ! yi-component of the current velocity at zi = -SmllNmbr meters       (m/s)
   REAL(SiKi), PARAMETER        :: SmllNmbr  = 9.999E-4             ! A small number representing epsilon for taking numerical derivatives.
   
   INTEGER                      :: I                                ! Generic index
      
      
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = ""               
  
   
      ! IF there are Morison elements, then compute the current components at each morison node elevation
      
   IF ( InitInp%NGridPts > 0 ) THEN    
         
      ALLOCATE ( InitOut%CurrVxi( InitInp%NGridPts ) , STAT=ErrStat )
      IF ( ErrStat /= ErrID_None )  THEN
         ErrMsg = ' Error allocating memory for the CurrVxi array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      
      ALLOCATE ( InitOut%CurrVyi( InitInp%NGridPts ) , STAT=ErrStat )
      IF ( ErrStat /= ErrID_None )  THEN
         ErrMsg = ' Error allocating memory for the CurrVyi array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      
      
         ! Loop over all of the points where current information is required
      
      DO I = 1, InitInp%NGridPts
         
         CALL Calc_Current( InitInp, InitInp%WaveKinGridzi(I), InitInp%EffWtrDpth, InitInp%DirRoot, CurrVxi, CurrVyi )        
         
         InitOut%CurrVxi(I) = CurrVxi
         InitOut%CurrVyi(I) = CurrVyi
       
      END DO     
     
   END IF   
      

      ! Compute the partial derivative for wave stretching
   CALL    Calc_Current( InitInp,  0.0_SiKi, InitInp%EffWtrDpth, InitInp%DirRoot, CurrVxi0, CurrVyi0 )
   CALL    Calc_Current( InitInp, -SmllNmbr, InitInp%EffWtrDpth, InitInp%DirRoot, CurrVxiS, CurrVyiS )

   InitOut%PCurrVxiPz0 = ( CurrVxi0 - CurrVxiS )/SmllNmbr                    ! xi-direction
   InitOut%PCurrVyiPz0 = ( CurrVyi0 - CurrVyiS )/SmllNmbr                    ! yi-direction
   
   
END SUBROUTINE Current_Init

   
END MODULE Current
!**********************************************************************************************************************************