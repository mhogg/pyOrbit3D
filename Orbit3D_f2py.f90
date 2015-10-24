
    subroutine ORBIT3D(ALPD,TH0D,R1,R2,IT,JP,H0,H1,PL,M1,M2,XH,YH,ZH, &
                       V1,V2,PHI1,PHI2,Yinit,LOADP,REV,N1,N2,NT,tnd,  &
                       erx,ery,erz,rotx,roty)
    
    ! This program was developed with the aim of determining the
    ! stability of the VentrAssist implantable rotary blood pump.
    ! Orbit3D determines the pressure distribution over the
    ! thrust bearing and conical journal bearing surfaces of the
    ! impeller, calculates the forces and moments on the impeller
    ! resulting from this pressure and calculates the position 
    ! (both translational and rotational) of the impeller within
    ! pump housing. From observation of this 'orbit', the stability
    ! of the impeller can be assessed.
    !
    ! Description of variables:
    !
    ! I          = Index of grid points in the theta direction
    ! J          = Index of grid points in the radial (r and rho) directions
    ! M          = Index of impeller blades
    ! N          = Index of grid points in time
    ! NTH        = Number of grid points in the theta direction
    ! NR         = Number of grid points in the radial (r) direction
    ! NRO        = Number of grid points in the radial (rho) direction
    ! NP         = Number of impeller blades
    ! NT         = Number of temporal grid points
    ! DTH        = Uniform grid spacing in the theta direction (non-dim.)
    ! DR         = Uniform grid spacing in the radial (r) direction (non-dim.)
    ! DRO        = Uniform grid spacing in the radial (rho) direction (non-dim.)
    ! DT         = Uniform grid spacing in time (non-dim.)
    ! TH(,)      = Theta location of grid points on impeller m (non-dim.)
    ! R()        = Radial (r) location of grid points (non-dim.)
    ! RO()       = Radial (rho) location of grid points (non-dim.)
    ! T()        = Time of temporal grid points (non-dim.)
    ! NDT        = Time at the nth time step (non-dim.)
    ! REV        = Number of revolution of the impeller
    ! ALPD       = Half-cone angle (degrees)
    ! ALP        = Half-cone angle (radians)
    ! TH0D       = The angular span of each impeller blade (deg.)
    ! TH0        = The angular span of each impeller blade (rad.)
    ! TH1        = Spacing of impeller blades (rad.) = 2*PI/NP
    ! R1         = The inner radius of the thrust bearing (m)
    ! R2         = The outer radius of the thrust bearing (m)
    ! LR         = R2 - R1 (m)
    ! RO1        = The inner radius of the conical journal bearing in the rho
    !              direction (m)
    ! RO2        = The outer radius of the conical journal bearing in the rho
    !              direction (m)
    ! LRO        = RO2 - RO1 (m)
    ! ROG        = ...
    ! PHI1       = Angle used to specify the gravitational force (rad.)
    ! PHI2       = Angle used to specify the gravitational force (rad.)
    ! LOADP      = Load parameter (non-dim.)
    ! XH,YH,ZH   = Amplitude of the displacement of the pump housing (non-dim.)
    ! V1         = Frequency of "shaking" of the pump housing (non-dim.)
    ! V2         = Frequency of unbalance force (non-dim.)
    ! H0TB       = Minimum (outlet) film thickness of thrust bearing (m) 
    ! H1TB       = Maximum (inlet) film thickness of thrust bearing (m)
    ! KTB        = H1TB/H0TB-1 (non-dim.)
    ! PLTB       = Fraction of thrust bearing that is untapered land
    ! HVC        = Vertical clearance of conical journal bearing = H0TB (m)
    ! H0CJB      = Minimum (outlet) film thickness of conical journal bearing (m) 
    ! H1CJB      = Maximum (inlet) film thickness of conical journal bearing (m)
    ! KCJB       = H1CJB/H0CJB-1 (non-dim.)
    ! PLCJB      = Fraction of conical journal bearing that is untapered land
    ! IT         = Transverse moment of inertia (non-dim.)
    ! JP         = Polar moment of inertia (non-dim.)
    ! M1         = Mass of the impeller (non-dim.)
    ! M2         = Magnitude of the unbalance mass (non-dim.) 
    ! PTB(,,)    = Pressure value at node i,j on the impeller blade m of the 
    !              thrust bearing at time step N (non-dim.) 
    ! PPITB(,,)  = Pressure value at node i,j on the impeller blade m of the 
    !              thrust bearing at time step N-1 (non-dim.)
    ! PCJB(,,)   = Pressure value at node i,j on the impeller blade m of the 
    !              conical journal bearing at time step N (non-dim.)
    ! PPICJB(,,) = Pressure value at node i,j on the impeller blade m of the 
    !              conical journal bearing at time step N-1 (non-dim.)
    ! FXTB       = Force in the X-direction on the impeller from the thrust 
    !              bearing (non-dim.)
    ! MXTB       = Moment about the X-axis on the impeller from the thrust
    !              bearing (non-dim.)
    ! MYTB       = Moment about the Y-axis on the impeller from the thrust 
    !              bearing (non-dim.)
    ! FXCJB      = Force in the X-direction on the impeller from the conical 
    !              journal bearing (non-dim.)
    ! FYCJB      = Force in the Y-direction on the impeller from the conical
    !              journal bearing (non-dim.)
    ! FZCJB      = Force in the Z-direction on the impeller from the conical
    !              journal bearing (non-dim.)
    ! MXCJB      = Moment about the X-axis on the impeller from the conical
    !              journal bearing (non-dim.)
    ! MYCJB      = Moment about the Y-axis on the impeller from the conical
    !              journal bearing (non-dim.)
    ! FX(,)      = Total force in the X-direction on the impeller (non-dim.)
    ! FY(,)      = Total force in the Y-direction on the impeller (non-dim.)
    ! FZ(,)      = Total force in the Z-direction on the impeller (non-dim.)
    ! MX(,)      = Total moment about the X-axis on the impeller (non-dim.)
    ! MY(,)      = Total moment about the Y-axis on the impeller (non-dim.)
    ! K1()       = Variable used in the Runge-Kutta method
    ! K2()       = Variable used in the Runge-Kutta method
    ! K3()       = Variable used in the Runge-Kutta method
    ! K4()       = Variable used in the Runge-Kutta method
    ! Y1()       = Eccentricity ratio in the X-direction
    ! Y2()       = Instantaneous velocity in the X-direction (non-dim.) 
    ! Y3()       = Eccentricity ratio in the Y-direction
    ! Y4()       = Instantaneous velocity in the Y-direction (non-dim.)
    ! Y5()       = Eccentricity ratio in the Z-direction
    ! Y6()       = Instantaneous velocity in the Z-direction (non-dim.)
    ! Y7()       = Rotational position of the impeller about the X-axis (non-dim.)
    ! Y8()       = Rotational velocity of the impeller about the X-axis (non-dim.)
    ! Y9()       = Rotational position of the impeller about the Y-axis (non-dim.)
    ! Y10()      = Rotational velocity of the impeller about the Y-axis (non-dim.)
    ! FLAG1      = Set to true if touchdown of the thrust bearing occurs, otherwise false
    ! FLAG2      = Set to true if touchdown of the conical journal bearing occurs, 
    !              otherwise false

    ! Variable declaration
    ! f2py intent(in)    ALPD,TH0D,R1,R2,IT,JP,M1,M2,Yinit,LOADP,REV
    ! f2py intent(in)    XH,YH,ZH,V1,V2,PHI1,PHI2,H0,H1,PL,N1,N2
    ! f2py intent(inout) tnd,erx,ery,erz,rotx,roty
    ! f2py depends(NT)   tnd,erx,ery,erz,rotx,roty
    
    REAL*4  REV, ALPD, TH0D, R1, R2
    REAL*4  PHI1, PHI2, XH, YH, ZH, V1, V2, LOADP
    REAL*4  H0TB, H1TB, KTB, PLTB, HVC, H0CJB, H1CJB, KCJB, PLCJB
    REAL*4  LR, RO1, RO2, LRO, IT, JP, M1, M2
    REAL*4  FZTB, MXTB, MYTB
    REAL*4  FXCJB, FYCJB, FZCJB, MXCJB, MYCJB
    REAL*8  PI, RTOD, DTOR         
    REAL*8  Yinit(10)
    REAL*8  tnd(NT),erx(NT),ery(NT),erz(NT),rotx(NT),roty(NT)
    REAL*8  K1(10), K2(10), K3(10), K4(10)
    REAL*8  TH0, ALP, TH1, DTH, DR, DRO, DT 
    REAL*8  NDT, SA, CA
    REAL*8, allocatable :: TH(:,:), R(:), RO(:)
    REAL*8, allocatable, dimension(:) :: T,Y1,Y2,Y3,Y4,Y5
    REAL*8, allocatable, dimension(:) :: Y6,Y7,Y8,Y9,Y10
    REAL*8, allocatable, dimension(:,:)   :: FX,FY,FZ,MX,MY
    REAL*8, allocatable, dimension(:,:,:) :: PTB,PPITB,PCJB,PPICJB
    INTEGER N1, N2, NTH, NR, NRO, NP, NT
    INTEGER I, J, M, N, KOUNT
    LOGICAL FLAG1, FLAG2  
    PARAMETER (PI=3.141592654, RTOD=180.0/PI, DTOR=PI/180.0)
    
    ! Initialise variables (do not assign initial values)
    NTH   = N1
    NR    = N1
    NRO   = N1
    NP    = N2
    FLAG1 = .FALSE.
    FLAG2 = .FALSE.
    KOUNT = 1
    
    ! Allocate arrays
    ! NOTE: Y1-Y10 should all be sized NT+1 because of time step loop calculations
    allocate (Y1(NT+1),Y2(NT+1),Y3(NT+1),Y4(NT+1),Y5(NT+1))
    allocate (Y6(NT+1),Y7(NT+1),Y8(NT+1),Y9(NT+1),Y10(NT+1),T(NT))
    allocate (FX(NT,NP),FY(NT,NP),FZ(NT,NP),MX(NT,NP),MY(NT,NP))
    allocate (TH(NTH,NP),R(NR),RO(NRO),PTB(NTH,NR,NP))
    allocate (PPITB(NTH,NR,NP),PCJB(NTH,NRO,NP),PPICJB(NTH,NRO,NP))
    
    ! Open results files
    OPEN (UNIT=1, FILE='plottraj.m', STATUS='UNKNOWN')          
    OPEN (UNIT=2, FILE='plotcjb.m',  STATUS='UNKNOWN')          
    OPEN (UNIT=3, FILE='plottb.m',   STATUS='UNKNOWN')
    
    ! Initialised dependent variables for the
    ! (a) thrust bearing
    H0TB = H0
    H1TB = H1
    PLTB = PL      
    KTB  = H1TB/H0TB-1.0
    
    ! (b) conical journal bearing
    ALP    = ALPD*DTOR
    SA     = DSIN(ALP)
    CA     = DCOS(ALP)
    RO1    = R1/SA
    RO2    = R2/SA
    ROG    = 1.05*RO2           ! Currently undefined. Got to do with RHO value of centre of mass
    ROG    = ROG/(RO2*SA)
    !HVC   = H0TB               ! CHECK THIS - Assigning vertical clearance, not parallel clearance
    !H0CJB = HVC*SA             ! CHECK THIS
    !H1CJB = H0CJB+50E-6        ! CHECK THIS
    !KCJB  = H1CJB/H0CJB-1.0    ! CHECK THIS
    !PLCJB = PLTB               ! CHECK THIS
    H0CJB  = H0
    H1CJB  = H1
    PLCJB  = PL
    KCJB   = H1CJB/H0CJB-1.0
    
    ! Initialise the coordinates of the discrete grid points 
    ! in the circumferential direction for the first impeller blade
    TH0=TH0D*DTOR
    DTH=TH0/(NTH-1)
    TH(1,1)=0.0
    do I=2,NTH
        TH(I,1)=TH(I-1,1)+DTH
    end do
    
    ! Initialise the coordinates of the discrete grid points 
    ! in the circumferential direction for the remaining blades
    TH1=2*PI/NP
    do M=2,NP
        do I=1,NTH
            TH(I,M)=TH(I,1)+(M-1)*TH1
        end do
    end do
    
    ! Initialise the coordinates of the grid points in the radial
    ! direction for the thrust bearing  
    LR=R2-R1
    DR=(LR/R2)/(NR-1)
    R(1)=R1/R2
    do J=2,NR
        R(J)=R(J-1)+DR
    end do
    
    ! Initialise the coordinates of the grid points in the radial
    ! direction for the conical journal bearing
    LRO=RO2-RO1
    DRO=LRO/R2/(NRO-1)
    RO(1)=RO1/R2
    do J=2,NRO
        RO(J)=RO(J-1)+DRO
    end do
    
    ! Initialise the numerical grid in time. The non-dimensional
    ! time step should be in the order of 1 degree, but is typically
    ! 0.5 or even 0.25 degrees depending on the operating conditions. 
    DT=2.0*PI*REV/(NT-1)
    T(1)=0.0D0
    do N=2,NT
        T(N)=T(N-1)+DT
    end do
    
    ! Initialise pressure to ambient (ie. zero) at t=0 for the
    PTB  = 0.0D0; PPITB  = 0.0D0
    PCJB = 0.0D0; PPICJB = 0.0D0
    
    ! Set initial conditions
    Y1(1)  = Yinit(1)
    Y2(1)  = Yinit(2)
    Y3(1)  = Yinit(3)
    Y4(1)  = Yinit(4)
    Y5(1)  = Yinit(5)
    Y6(1)  = Yinit(6)
    Y7(1)  = Yinit(7)
    Y8(1)  = Yinit(8)
    Y9(1)  = Yinit(9)
    Y10(1) = Yinit(10)
    
    ! The following DO WHILE loop uses the classical form of the fourth-order
    ! Runge-Kutta method to solve the five second-order non-linear differential
    ! equations that govern the translational and rotational movement of the
    ! impeller.  This method requires four "estimates" before giving the 
    ! actual position for each time step.  Thus four steps are involved.
    
    ! DO WHILE loop control:
    ! If FLAG1 is becomes true during the course of the run then touchdown 
    ! has occurred in the thrust bearing.  Furthermore, if FLAG2 becomes true 
    ! then touchdown has occurred in the conical journal bearing. If either 
    ! of these occurs then the loop is exited and the program stopped.
    N=1
    
    DO WHILE ((FLAG1 .EQV. (.FALSE.)) .AND. (FLAG2 .EQV. (.FALSE.)) &
                                      .AND. (N .LE. NT))
        
        ! *** FIRST STEP ***
        
        NDT = T(N)       
        
        ! Call Subroutine TB (ie. Thrust Bearing) to calculate the pressure
        ! distribution over the thrust bearing.  This subroutine also integrates
        ! this pressure to determine the forces and moments on the impeller.
        CALL TB(H0TB,H1TB,KTB,PLTB,TH0,TH,DTH,NTH,R1,R2,R,DR,NR,NP,NDT, &   
                N1,N2,Y5(N),Y6(N),Y7(N),Y8(N),Y9(N),Y10(N),PTB,PPITB,   &  
                FZTB,MXTB,MYTB,FLAG1)                                    
        
        ! Call Subroutine CJB (ie. Conical Journal Bearing) to calculate the
        ! pressure distribution over the conical journal bearing.  This subroutine
        ! also integrates this pressure to determine the forces and moments on the
        ! impeller.                                                       
        CALL CJB(H0CJB,H1CJB,KCJB,PLCJB,ALP,TH0,TH,DTH,NTH,RO1,RO2,ROG,  &
                 RO,DRO,NRO,NP,NDT,N1,N2,Y1(N),Y2(N),Y3(N),Y4(N),Y5(N),  &
                 Y6(N),Y7(N),Y8(N),Y9(N),Y10(N),PCJB,PPICJB,FXCJB,FYCJB, &
                 FZCJB,MXCJB,MYCJB,FLAG2)
        
        ! Calculate the net forces and moments on the impeller by summation
        ! of the individual contributions of the thrust bearing and the 
        ! conical journal bearing
        FX(N,1) = FXCJB
        FY(N,1) = FYCJB
        FZ(N,1) = FZTB + FZCJB
        MX(N,1) = MXTB + MXCJB
        MY(N,1) = MYTB + MYCJB
        
        ! Write the pressure field at this time step to file (this is done 
        ! during the first step of the Runge-Kutta method only) at the first
        ! time step and for every twentieth time step thereafter. If you do
        ! not want to print the pressure profile to file, then set 
        ! (N .EQ. 1) to (N .EQ. 0) and (KOUNT .EQ. 20) to (KOUNT .EQ. 0) 
        ! in the if statement below
        IF ( (N .EQ. 0) .OR. (KOUNT .EQ. 0) ) THEN
            !CALL PLTCJB(PCJB,N,NTH,NRO,NP,NT,SA,RO,TH,NDT,N1,N2)        
            !CALL PLTTB(PTB,N,NTH,NR,NP,NT,R,TH,NDT,N1,N2)
            KOUNT = 1
        ELSE
            KOUNT = KOUNT + 1
        END IF
        
        K1(1)  = DT*Y2(N)
        K1(2)  = DT*( (LOADP/M1)*FX(N,1) + (V1**2)*XH*SIN(V1*NDT) &
                 -(M2/M1)*(V2**2)*COS(V2*NDT)+SIN(PHI1)*COS(PHI2)/M1 )
        K1(3)  = DT*Y4(N)
        K1(4)  = DT*( (LOADP/M1)*FY(N,1) + (V1**2)*YH*SIN(V1*NDT) &
                 -(M2/M1)*(V2**2)*SIN(V2*NDT)+SIN(PHI1)*SIN(PHI2)/M1 )
        K1(5)  = DT*Y6(N)
        K1(6)  = DT*( (LOADP/M1)*FZ(N,1) + (V1**2)*ZH*SIN(V1*NDT) + COS(PHI1)/M1 )
        K1(7)  = DT*Y8(N)
        K1(8)  = DT*( JP*Y10(N)+(R2/H0TB)*MX(N,1)/IT )
        K1(9)  = DT*Y10(N)
        K1(10) = DT*( -JP*Y8(N)+(R2/H0TB)*MY(N,1)/IT )
        
        
        ! *** SECOND STEP ***
        
        NDT = T(N)+0.5*DT       
        CALL TB(H0TB,H1TB,KTB,PLTB,TH0,TH,DTH,NTH,R1,R2,R,DR,NR,NP,NDT, &
            N1,N2,Y5(N)+K1(5)/2,Y6(N)+K1(6)/2,Y7(N)+K1(7)/2,Y8(N)+      &
            K1(8)/2,Y9(N)+K1(9)/2,Y10(N)+K1(10)/2,PTB,PPITB,FZTB,MXTB,  &
            MYTB,FLAG1)
        CALL CJB(H0CJB,H1CJB,KCJB,PLCJB,ALP,TH0,TH,DTH,NTH,RO1,RO2,ROG, &
            RO,DRO,NRO,NP,NDT,N1,N2,Y1(N)+K1(1)/2,Y2(N)+K1(2)/2,Y3(N)+  &
            K1(3)/2,Y4(N)+K1(4)/2,Y5(N)+K1(5)/2,Y6(N)+K1(6)/2,Y7(N)+    &
            K1(7)/2,Y8(N)+K1(8)/2,Y9(N)+K1(9)/2,Y10(N)+K1(10)/2,PCJB,   &
            PPICJB,FXCJB,FYCJB,FZCJB,MXCJB,MYCJB,FLAG2)
        FX(N,2) = FXCJB
        FY(N,2) = FYCJB
        FZ(N,2) = FZTB + FZCJB
        MX(N,2) = MXTB + MXCJB
        MY(N,2) = MYTB + MYCJB
        
        K2(1)  = DT*( Y2(N)+0.5*K1(2) )
        K2(2)  = DT*( (LOADP/M1)*FX(N,2) + (V1**2)*XH*SIN(V1*NDT) &
                -(M2/M1)*(V2**2)*COS(V2*NDT)+SIN(PHI1)*COS(PHI2)/M1 )
        K2(3)  = DT*( Y4(N)+0.5*K1(4) )
        K2(4)  = DT*( (LOADP/M1)*FY(N,2) + (V1**2)*YH*SIN(V1*NDT) &
                -(M2/M1)*(V2**2)*SIN(V2*NDT)+SIN(PHI1)*SIN(PHI2)/M1 )
        K2(5)  = DT*( Y6(N)+0.5*K1(6) )
        K2(6)  = DT*( (LOADP/M1)*FZ(N,2) + (V1**2)*ZH*SIN(V1*NDT) + COS(PHI1)/M1 )
        K2(7)  = DT*( Y8(N)+0.5*K1(8) )
        K2(8)  = DT*( JP*(Y10(N)+0.5*K1(10))+(R2/H0TB)*MX(N,2)/IT )
        K2(9)  = DT*( Y10(N)+0.5*K1(10) )
        K2(10) = DT*( -JP*(Y8(N)+0.5*K1(8))+(R2/H0TB)*MY(N,2)/IT )
        
        
        ! *** THIRD STEP ***
        
        NDT = T(N)+0.5*DT       
        CALL TB(H0TB,H1TB,KTB,PLTB,TH0,TH,DTH,NTH,R1,R2,R,DR,NR,NP,NDT, &
             N1,N2,Y5(N)+K2(5)/2,Y6(N)+K2(6)/2,Y7(N)+K2(7)/2,Y8(N)+     &
             K2(8)/2,Y9(N)+K2(9)/2,Y10(N)+K2(10)/2,PTB,PPITB,FZTB,MXTB, &
             MYTB,FLAG1)
        CALL CJB(H0CJB,H1CJB,KCJB,PLCJB,ALP,TH0,TH,DTH,NTH,RO1,RO2,ROG, &
             RO,DRO,NRO,NP,NDT,N1,N2,Y1(N)+K2(1)/2,Y2(N)+K2(2)/2,Y3(N)+ &
             K2(3)/2,Y4(N)+K2(4)/2,Y5(N)+K2(5)/2,Y6(N)+K2(6)/2,Y7(N)+   &
             K2(7)/2,Y8(N)+K2(8)/2,Y9(N)+K2(9)/2,Y10(N)+K2(10)/2,PCJB,  &
             PPICJB,FXCJB,FYCJB,FZCJB,MXCJB,MYCJB,FLAG2)
        FX(N,3) = FXCJB
        FY(N,3) = FYCJB
        FZ(N,3) = FZTB + FZCJB
        MX(N,3) = MXTB + MXCJB
        MY(N,3) = MYTB + MYCJB
        
        K3(1)  = DT*( Y2(N)+0.5*K2(2) )
        K3(2)  = DT*( (LOADP/M1)*FX(N,3) + (V1**2)*XH*SIN(V1*NDT) &
                 -(M2/M1)*(V2**2)*COS(V2*NDT)+SIN(PHI1)*COS(PHI2)/M1 )
        K3(3)  = DT*( Y4(N)+0.5*K2(4) )
        K3(4)  = DT*( (LOADP/M1)*FY(N,3) + (V1**2)*YH*SIN(V1*NDT) &
                 -(M2/M1)*(V2**2)*SIN(V2*NDT)+SIN(PHI1)*SIN(PHI2)/M1 )
        K3(5)  = DT*( Y6(N)+0.5*K2(6) )
        K3(6)  = DT*( (LOADP/M1)*FZ(N,3) + (V1**2)*ZH*SIN(V1*NDT) + COS(PHI1)/M1 )
        K3(7)  = DT*( Y8(N)+0.5*K2(8) )
        K3(8)  = DT*( JP*(Y10(N)+0.5*K2(10))+(R2/H0TB)*MX(N,3)/IT )
        K3(9)  = DT*( Y10(N)+0.5*K2(10) )
        K3(10) = DT*( -JP*(Y8(N)+0.5*K2(8))+(R2/H0TB)*MY(N,3)/IT )
        
        
        ! *** FOURTH STEP ***
        
        NDT = T(N)+DT       
        CALL TB(H0TB,H1TB,KTB,PLTB,TH0,TH,DTH,NTH,R1,R2,R,DR,NR,NP,NDT, &
             N1,N2,Y5(N)+K3(5),Y6(N)+K3(6),Y7(N)+K3(7),Y8(N)+K3(8),     &
             Y9(N)+K3(9),Y10(N)+K3(10),PTB,PPITB,FZTB,MXTB,MYTB,FLAG1)
        CALL CJB(H0CJB,H1CJB,KCJB,PLCJB,ALP,TH0,TH,DTH,NTH,RO1,RO2,ROG, &
             RO,DRO,NRO,NP,NDT,N1,N2,Y1(N)+K3(1),Y2(N)+K3(2),Y3(N)+     &
             K1(3),Y4(N)+K3(4),Y5(N)+K3(5),Y6(N)+K3(6),Y7(N)+K3(7),     &
             Y8(N)+K3(8),Y9(N)+K3(9),Y10(N)+K3(10),PCJB,PPICJB,FXCJB,   &
             FYCJB,FZCJB,MXCJB,MYCJB,FLAG2)
        FX(N,4) = FXCJB
        FY(N,4) = FYCJB
        FZ(N,4) = FZTB + FZCJB
        MX(N,4) = MXTB + MXCJB
        MY(N,4) = MYTB + MYCJB

        K4(1)  = DT*( Y2(N)+K3(2) )
        K4(2)  = DT*( (LOADP/M1)*FX(N,4) + (V1**2)*XH*SIN(V1*NDT) &
                 -(M2/M1)*(V2**2)*COS(V2*NDT)+SIN(PHI1)*COS(PHI2)/M1 )
        K4(3)  = DT*( Y4(N)+K3(4) )
        K4(4)  = DT*( (LOADP/M1)*FY(N,4) + (V1**2)*YH*SIN(V1*NDT) &
                 -(M2/M1)*(V2**2)*SIN(V2*NDT)+SIN(PHI1)*SIN(PHI2)/M1 )
        K4(5)  = DT*( Y6(N)+K3(6) )
        K4(6)  = DT*( (LOADP/M1)*FZ(N,4) + (V1**2)*ZH*SIN(V1*NDT) + COS(PHI1)/M1 )
        K4(7)  = DT*( Y8(N)+K3(8) )
        K4(8)  = DT*( JP*(Y10(N)+K3(10))+(R2/H0TB)*MX(N,4)/IT )
        K4(9)  = DT*( Y10(N)+K3(10) )
        K4(10) = DT*( -JP*(Y8(N)+K3(8))+(R2/H0TB)*MY(N,4)/IT )


        ! Calculate the variables at the next time step based on the
        ! results of these four steps. 
        Y1(N+1)  = Y1(N)+(K1(1)+2.*K2(1)+2.*K3(1)+K4(1))/6.
        Y2(N+1)  = Y2(N)+(K1(2)+2.*K2(2)+2.*K3(2)+K4(2))/6.
        Y3(N+1)  = Y3(N)+(K1(3)+2.*K2(3)+2.*K3(3)+K4(3))/6.
        Y4(N+1)  = Y4(N)+(K1(4)+2.*K2(4)+2.*K3(4)+K4(4))/6.
        Y5(N+1)  = Y5(N)+(K1(5)+2.*K2(5)+2.*K3(5)+K4(5))/6.
        Y6(N+1)  = Y6(N)+(K1(6)+2.*K2(6)+2.*K3(6)+K4(6))/6.
        Y7(N+1)  = Y7(N)+(K1(7)+2.*K2(7)+2.*K3(7)+K4(7))/6.
        Y8(N+1)  = Y8(N)+(K1(8)+2.*K2(8)+2.*K3(8)+K4(8))/6.
        Y9(N+1)  = Y9(N)+(K1(9)+2.*K2(9)+2.*K3(9)+K4(9))/6.
        Y10(N+1) = Y10(N)+(K1(10)+2.*K2(10)+2.*K3(10)+K4(10))/6.

        ! Write results to the screen
        IF (N .EQ. 1) THEN
            WRITE(*,'(A10,5A12)') 't*','erx','ery','erz','xrot*','yrot*'
        END IF
        WRITE(*,'(F10.3,5F12.6)') T(N)*RTOD,Y1(N),Y3(N),Y5(N),Y7(N),Y9(N)

        ! Now redundant, due to f2py i.e. results are passed directly back to Python
        ! Write results of all variables to file 
        ! CALL PLTTRA(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,DT,T,N,NT,FLAG1,FLAG2)

        ! If touchdown occurs, then print message to screen
        IF (FLAG1 .EQV. (.TRUE.)) THEN
            WRITE(*,'(A40)') 'Touchdown reported by SUBROUTINE TB'
        ELSE IF (FLAG2 .EQV. (.TRUE.)) THEN
            WRITE(*,'(A40)') 'Touchdown reported by SUBROUTINE CJB'
        END IF
        
        ! Increment time for next time step
        N = N + 1

    end do
      
    ! Get results for output back to Python
    do i=1,NT
        tnd(i)  = T(i)
        erx(i)  = Y1(i)/(1.0-Y5(i))
        ery(i)  = Y3(i)/(1.0-Y5(i))
        erz(i)  = Y5(i)
        rotx(i) = Y7(i)
        roty(i) = Y9(i)
    end do     

    ! Deallocate arrays. This is required to prevent crash when f2py module
    ! is rerun (otherwise need checks for "if allocated, then...")
    deallocate(T,Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,FX,FY,FZ,MX,MY)  
    deallocate(TH,R,RO,PTB,PPITB,PCJB,PPICJB)

    ! End and close result files
    END FILE (1)
    END FILE (2)
    END FILE (3)
    CLOSE (1)
    CLOSE (2)
    CLOSE (3)
     
    contains
      
! =============================================================================

    SUBROUTINE TB(H0,H1,K,PLAND,TH0,TH,DTH,NTH,R1,R2,R,DR,NR,NP,NDT, &
                  N1,N2,ERZ,ERZDT,XROT,XROTDT,YROT,YROTDT,P,PPI,FZT, &
                  MXT,MYT,FLAG2)

    ! This subroutine first calculates the pressure field over the thrust
    ! bearing surfaces of each impeller blade. This pressure field is then
    ! integrated over the respective bearing surface to calculate the
    ! forces and moments on the impeller from the thrust bearing alone. 

    ! Description of variables (those not already described or different 
    ! to those in the main program):
    ! 
    ! FOR       = Successive over-relaxation factor to accelerate iterative process
    ! DTOR      = Converts degrees to radians 
    ! H0        = Minimum film thickness (m)
    ! H1        = Maximum film thickness (m)
    ! K         = H1/H0-1 (non-dim.)
    ! PLAND     = Fraction of angular span of bearing that is untapered land
    ! X(,,)     = X-coordinate of grid point i,j on impeller blade m (non-dim.)
    ! Y(,,)     = Y-coordinate of grid point i,j on impeller blade m (non-dim.)
    ! FZ()      = Force on the impeller in the Z-direction due to the thrust bearing
    !             surface of impeller blade m (non-dim.)
    ! MX()      = Moment about the X-axis on impeller due to the thrust bearing
    !             surface of impeller blade m (non-dim.)
    ! MY()      = Moment about the Y-axis on impeller due to the thrust bearing
    !             surface of impeller blade m (non-dim.)
    ! FZT       = Total force on impeller in the Z-direction due to the thrust 
    !             bearing (non-dim.)
    ! MXT       = Total moment on the impeller about the X-axis due to the thrust
    !             bearing (non-dim.)
    ! MYT       = Total moment on the impeller about the Y-axis due to the thrust
    !             bearing (non-dim.)
    ! TH0A      = Angular span of the land section (rad.)
    ! TH0B      = Angular span of the tapered section (rad.)
    ! ERZ       = Eccentricity of the impeller in the Z-direction (non-dim.)
    ! ERZDT     = Time derivative of ERZ (ie. Z-velocity) (non-dim.)
    ! XROT      = Rotational position of the impeller about the X-axis (non-dim.)
    ! XROTDT    = Time derivative of XROT (ie. X rotational velocity) (non-dim.)
    ! XROT      = Rotational position of the impeller about the Y-axis (non-dim.)
    ! XROTDT    = Time derivative of YROT (ie. Y rotational velocity) (non-dim.)
    ! A()-D()   = Coefficients for the Thomas algorithm subroutine
    ! H(,,)     = Film thickness at any point (non-dim.)
    ! DHDTH(,,) = Derivative of H with respect to theta (non-dim.)
    ! DHDR(,,)  = Derivative of H with respect to r (non-dim.)
    ! DHDT(,,)  = Derivative of H with respect to time (non-dim.)
    ! P(,,)     = Pressure at node i,j on impeller blade m at time step N (non-dim.)
    ! PPI(,,)   = Pressure at node i,j on impeller blade m at time step N-1 (non-dim.) 
    ! FLAG1     = Set to true if solution converges, otherwise false
    ! FLAG2     = Set to true if touchdown occurs, otherwise false

    ! Declaration of variables
    PARAMETER (PI=3.141592654, DTOR=PI/180.0, FOR=1.5)
    REAL H0, H1, K, R1, R2, LR, NDR, PLAND, TH0D 
    REAL X(N1,N1,N2), Y(N1,N1,N2)
    REAL PR(N1,N1,N2), PRX(N1,N1,N2), PRY(N1,N1,N2) 
    REAL FZ(N2), MX(N2), MY(N2), FZT, MXT, MYT
    REAL*8 TH(N1,N2), TH0, TH0A, TH0B, TH1, R(N1)
    REAL*8 DTH, DTH2, DR, DR2, NDT
    REAL*8 ERZ, ERZDT, XROT, XROTDT, YROT, YROTDT 
    REAL*8 A(N1), B(N1), C(N1), D(N1)
    REAL*8 H(N1,N1,N2), DHDTH(N1,N1,N2), DHDR(N1,N1,N2)
    REAL*8 DHDT(N1,N1,N2), P(N1,N1,N2), PPI(N1,N1,N2)
    INTEGER NTH, NTHM1, NTHM2, NR, NRM1, NRM2
    INTEGER I, J, M, N, NP, N1, N2
    LOGICAL FLAG1, FLAG2

    ! Initialise dependent variables
    NTHM1 = NTH-1
    NTHM2 = NTH-2
    NRM1  = NR-1
    NRM2  = NR-2
    TH0A  = TH0*PLAND
    TH0B  = TH0-TH0A
    TH1   = 2*PI/NP
    LR    = R2-R1
    NDR   = (R2/LR)**2
    DTH2  = DTH**2
    DR2   = DR**2

    ! Initialise variables
    FZT = 0.0
    MXT = 0.0
    MYT = 0.0

    ! Calculate the pressure distribution and subsequently the forces and moments 
    ! on the impeller as a result of this pressure over each pad m=1,2,...,Np                
    do M=1,NP
               
        ! Calculate the local film thickness and spatial and temporal derivatives
        do J=1,NR
        
            do I=1,NTH
                
                ! Film thickness for untapered land
                H(I,J,M)= 1.0 + ERZ + XROT*R(J)*SIN(TH(I,M)+NDT) &
                                    - YROT*R(J)*COS(TH(I,M)+NDT)  
                
                ! Derivative with respect to (w.r.t.) theta for untapered land
                DHDTH(I,J,M) = XROT*R(J)*COS(TH(I,M)+NDT) &
                              +YROT*R(J)*SIN(TH(I,M)+NDT)
                
                ! Film thickess and derivative w.r.t. theta for taper
                IF ((TH(I,M)-(M-1)*TH1) .GE. TH0A) THEN
                    H(I,J,M) = H(I,J,M)+K*(TH(I,M)-(TH0A+(M-1)*TH1))/TH0B
                    DHDTH(I,J,M) = DHDTH(I,J,M) + K/TH0B
                END IF           
                
                ! If the film thickness is less than or equal to zero, then
                ! touchdown has occurred. Set FLAG2 to true.
                IF (H(I,J,M) .LE. 0.0D0) THEN
                    FLAG2 = .TRUE.
                    RETURN
                END IF
                
                ! Radial and time derivatives
                DHDR(I,J,M) = XROT*SIN(TH(I,M)+NDT)-YROT*COS(TH(I,M)+NDT)
                DHDT(I,J,M) = ERZDT + (XROT-YROTDT)*R(J)*COS(TH(I,M)+NDT) &
                                    + (XROTDT+YROT)*R(J)*SIN(TH(I,M)+NDT)
                
            end do
        end do        

        ! Initialise iterative process loop control variable FLAG1
        FLAG1 = .FALSE.

        ! Iterative process to solve for the pressure distribution

        !  --------------- BEGIN ITERATIVE PROCESS ---------------
        
        DO WHILE(FLAG1 .EQV. (.FALSE.))
            
            ! Sweep all rows j=2,3,...,J-1
            do J=2,NRM1
                
                do I=2,NTHM1
                    A(I-1) = 1./(R(J)*DTH2)-1.5*DHDTH(I,J,M)/(R(J)*H(I,J,M)*DTH)
                    B(I-1) = -2.*(R(J)/DR2+1./(R(J)*DTH2))
                    C(I-1) = 1./(R(J)*DTH2)+1.5*DHDTH(I,J,M)/(R(J)*H(I,J,M)*DTH)
                    D(I-1) = -(NDR*R(J))/(H(I,J,M)**3)*DHDTH(I,J,M)   &
                             +(2.*NDR*R(J))/(H(I,J,M)**3)*DHDT(I,J,M) &
                             -(R(J)/DR2-(1.5*R(J))/(H(I,J,M)*DR)      &
                             *DHDR(I,J,M)-0.5/DR)*P(I,J-1,M)          &
                             -(R(J)/DR2+(1.5*R(J))/(H(I,J,M)*DR)      &
                             *DHDR(I,J,M)+0.5/DR)*P(I,J+1,M)       
                end do
            
                ! Use Thomas algorithm to solve tridiagonal system of algebraic equations
                CALL THOMAS(NTHM2,A,B,C,D,N1)
                
                ! Capture solution from Thomas and apply SOR factor, For
                do I=2,NTHM1
                    P(I,J,M)   = PPI(I,J,M)+FOR*(D(I-1)-PPI(I,J,M))             
                    PPI(I,J,M) = P(I,J,M)
                end do
         
            end do
         
            ! Sweep all columns i=2,3,...,I-1
            do I=2,NTHM1
                
                do J=2,NRM1
                    A(J-1) = (R(J)/DR-(1.5*R(J)*DHDR(I,J,M))/H(I,J,M)-0.5)/DR
                    B(J-1) = -2.*(R(J)/DR2+1./(R(J)*DTH2))
                    C(J-1) = (R(J)/DR+(1.5*R(J)*DHDR(I,J,M))/H(I,J,M)+0.5)/DR
                    D(J-1) = -(NDR*R(J))/(H(I,J,M)**3)*DHDTH(I,J,M)   &
                             +(2.*NDR*R(J))/(H(I,J,M)**3)*DHDT(I,J,M) &
                             -(1./(R(J)*DTH2)-1.5/(R(J)*H(I,J,M)*DTH) &
                             *DHDTH(I,J,M))*P(I-1,J,M)                &
                             -(1./(R(J)*DTH2)+1.5/(R(J)*H(I,J,M)*DTH) &
                             *DHDTH(I,J,M))*P(I+1,J,M)
                end do
                
                ! Use Thomas algorithm to solve tridiagonal system of algebraic equations
                CALL THOMAS(NRM2,A,B,C,D,N1)
                
                ! Capture solution from Thomas and apply SOR factor, For
                do J=2,NRM1
                    P(I,J,M) = PPI(I,J,M)+FOR*(D(J-1)-PPI(I,J,M))
                end do
            
            end do
         
            ! Convergence test
            CALL CVERGE(NTH,NR,P,PPI,FLAG1,M,N1,N2)
          
        end do
        
        !  --------------- END ITERATIVE PROCESS ---------------
        
        ! Transform the position of each grid point on the lubricating surfaces of 
        ! the rotating impeller in Cartesian x-y coordinates relative to the pump 
        ! housing (which is fixed)
        do J=1,NR
            do I=1,NTH
                X(I,J,M) = R(J)*COS(TH(I,M)+NDT)
                Y(I,J,M) = R(J)*SIN(TH(I,M)+NDT)
            end do
        end do
        
        ! Calculate the force (Fz) and moments (Mx and My) on the impeller 
        ! due to each separate blade
        do J=1,NR
            do I=1,NTH
                PR(I,J,M)  =  P(I,J,M)*R(J)
                PRX(I,J,M) = -P(I,J,M)*R(J)*X(I,J,M)
                PRY(I,J,M) =  P(I,J,M)*R(J)*Y(I,J,M)
            end do
        end do
        FZ(M) = DBLINT(PR, DTH,DR,NTHM1,NRM1,N1,N2,M)
        MX(M) = DBLINT(PRY,DTH,DR,NTHM1,NRM1,N1,N2,M) 
        MY(M) = DBLINT(PRX,DTH,DR,NTHM1,NRM1,N1,N2,M)
        
        ! Calculate the total force (Fz) and moment (Mx and My) on the centre
        ! of mass of the impeller
        FZT = FZT + FZ(M)  
        MXT = MXT + MX(M)
        MYT = MYT + MY(M)
        
    end do
    
    return
    end
    
    
    ! =========================================================================

    SUBROUTINE CJB(H0,H1,K,PLAND,ALP,TH0,TH,DTH,NTH,RO1,RO2,ROG,RO,    &
                   DRO,NRO,NP,NDT,N1,N2,ERX,ERXDT,ERY,ERYDT,ERZ,ERZDT, &
                   XROT,XROTDT,YROT,YROTDT,P,PPI,FXT,FYT,FZT,MXT,MYT,  &
                   FLAG2)

    ! This subroutine first calculates the pressure field over the conical
    ! journal bearing surfaces of each impeller blade. This pressure field 
    ! is then integrated over the respective bearing surface to calculate
    ! the forces and moments on the impeller from the conical journal bearing 
    ! alone. 

    ! Description of variables (those not already described or different 
    ! to those in the main program):
    ! 
    ! FOR       = Successive over-relaxation factor to accelerate iterative process
    ! DTOR      = Converts degrees to radians 
    ! H0        = Minimum film thickness (m)
    ! H1        = Maximum film thickness (m)
    ! K         = H1/H0-1 (non-dim.)
    ! PLAND     = Fraction of angular span of bearing that is untapered land
    ! FX()      = Force on the impeller in the X-direction due to the conical 
    !             journal bearing surface of impeller blade m (non-dim.)
    ! FY()      = Force on the impeller in the Y-direction due to the conical
    !             journal bearing surface of impeller blade m (non-dim.)
    ! FZ()      = Force on the impeller in the Z-direction due to the conical
    !             journal bearing surface of impeller blade m (non-dim.)
    ! MX()      = Moment about the X-axis on impeller due to the conical journal
    !             bearing surface of impeller blade m (non-dim.)
    ! MY()      = Moment about the Y-axis on impeller due to the conical journal
    !             bearing surface of impeller blade m (non-dim.)
    ! FXT       = Total force on impeller in the X-direction due to the conical 
    !             journal bearing (non-dim.)
    ! FYT       = Total force on impeller in the Y-direction due to the conical
    !             journal bearing (non-dim.)
    ! FZT       = Total force on impeller in the Z-direction due to the conical
    !             journal bearing (non-dim.)
    ! MXT       = Total moment on the impeller about the X-axis due to the conical
    !             journal bearing (non-dim.)
    ! MYT       = Total moment on the impeller about the Y-axis due to the conical
    !             journal bearing (non-dim.)
    ! TH0A      = Angular span of the land section (rad.)
    ! TH0B      = Angular span of the tapered section (rad.)
    ! ERX       = Eccentricity of the impeller in the X-direction (non-dim.)
    ! ERXDT     = Time derivative of ERX (ie. X-velocity) (non-dim.)
    ! ERY       = Eccentricity of the impeller in the Y-direction (non-dim.)
    ! ERYDT     = Time derivative of ERY (ie. Y-velocity) (non-dim.)
    ! ERZ       = Eccentricity of the impeller in the Z-direction (non-dim.)
    ! ERZDT     = Time derivative of ERZ (ie. Z-velocity) (non-dim.)
    ! XROT      = Rotational position of the impeller about the X-axis (non-dim.)
    ! XROTDT    = Time derivative of XROT (ie. X rotational velocity) (non-dim.)
    ! XROT      = Rotational position of the impeller about the Y-axis (non-dim.)
    ! XROTDT    = Time derivative of YROT (ie. Y rotational velocity) (non-dim.)
    ! A()-D()   = Coefficients for the Thomas algorithm subroutine
    ! H(,,)     = Film thickness at any point (non-dim.)
    ! DHDTH(,,) = Derivative of H with respect to theta (non-dim.)
    ! DHDRO(,,) = Derivative of H with respect to r (non-dim.)
    ! DHDT(,,)  = Derivative of H with respect to time (non-dim.)
    ! P(,,)     = Pressure at node i,j on impeller blade m at time step N (non-dim.)
    ! PPI(,,)   = Pressure at node i,j on impeller blade m at time step N-1 (non-dim.) 
    ! FLAG1     = Set to true if solution converges, otherwise false
    ! FLAG2     = Set to true if touchdown occurs, otherwise false

    ! Declaration of variables
    PARAMETER (PI=3.141592654, DTOR=PI/180.0, FOR = 1.5)
    REAL H0, H1, K, RO1, RO2, ROG, LRO, NDR, PLAND
    REAL TH0D, ALPD
    REAL FX(N2), FY(N2), FZ(N2), MX(N2), MY(N2)
    REAL FXT, FYT, FZT, MXT, MYT
    REAL DFX(N1,N1,N2), DFY(N1,N1,N2), DFZ(N1,N1,N2)
    REAL DMX(N1,N1,N2), DMY(N1,N1,N2)
    REAL*8 TH(N1,N2), TH0, TH0A, TH0B, TH1, RO(N1)
    REAL*8 DTH, DTH2, DRO, DRO2
    REAL*8 ERX, ERXDT, ERY, ERYDT, ERZ, ERZDT
    REAL*8 XROT, XROTDT, YROT, YROTDT
    REAL*8 A(N1), B(N1), C(N1), D(N1), NDT, ALP, SA, CA
    REAL*8 H(N1,N1,N2), DHDRO(N1,N1,N2), DHDTH(N1,N1,N2)
    REAL*8 DHDT(N1,N1,N2), P(N1,N1,N2), PPI(N1,N1,N2)
    INTEGER NTH, NRO, NTHM1, NTHM2, NROM1, NROM2
    INTEGER I, J, M, NP, N1, N2
    LOGICAL FLAG1, FLAG2

    ! Independent variables
    SA   = DSIN(ALP)
    CA   = DCOS(ALP)

    ! Dependent variables
    NTHM1 = NTH-1
    NTHM2 = NTH-2
    NROM1 = NRO-1
    NROM2 = NRO-2
    TH0A  = TH0*PLAND
    TH0B  = TH0-TH0A
    TH1   = 2*PI/NP
    LRO   = RO2-RO1
    NDR   = (RO2/LRO)**2
    DTH2  = DTH**2
    DRO2  = DRO**2

    ! Initialise variables
    FXT = 0.0
    FYT = 0.0
    FZT = 0.0
    MXT = 0.0
    MYT = 0.0

    ! Calculate the pressure distribution and subsequently the forces and moments 
    ! on the impeller as a result of this pressure over each pad m=1,2,...,Np       
    DO M=1,NP
      
        ! Calculate local film thickness and the spatial and temporal derivatives 
        DO J=1,NRO
            DO I=1,NTH
                H(I,J,M) = (1.0-ERZ                                &
                            -ERX*COS(TH(I,M)+NDT)                  &
                            -ERY*SIN(TH(I,M)+NDT) )*SA             &            
                            -XROT*(RO(J)-ROG)*SIN(TH(I,M)+NDT)     &
                            +YROT*(RO(J)-ROG)*COS(TH(I,M)+NDT) 
                DHDTH(I,J,M) = ( ERX*SIN(TH(I,M)+NDT)              &               
                                -ERY*COS(TH(I,M)+NDT) )*SA         &                
                                -XROT*(RO(J)-ROG)*COS(TH(I,M)+NDT) &
                                -YROT*(RO(J)-ROG)*SIN(TH(I,M)+NDT) 
                IF ((TH(I,M)-(M-1)*TH1) .GE. TH0A) THEN           
                    H(I,J,M) = H(I,J,M) + (K*(TH(I,M)-(TH0A+(M-1)*TH1))/TH0B)*SA                  
                    DHDTH(I,J,M) = DHDTH(I,J,M) + ( K/TH0B )*SA                                    
                END IF                                                       
                IF (H(I,J,M) .LE. 0.0D0) THEN
                    FLAG2 = .TRUE.
                    RETURN
                END IF
                DHDRO(I,J,M) = -XROT*SIN(TH(I,M)+NDT)+YROT*COS(TH(I,M)+NDT)            
                DHDT(I,J,M)  = (-ERZDT                                      &
                                +(ERX-ERYDT)*SIN(TH(I,M)+NDT)               &
                                -(ERXDT+ERY)*COS(TH(I,M)+NDT) )*SA          &               
                                -(XROT-YROTDT)*(RO(J)-ROG)*COS(TH(I,M)+NDT) &
                                -(XROTDT+YROT)*(RO(J)-ROG)*SIN(TH(I,M)+NDT)
            end do
        end do
        
        ! Initialise iterative process loop control variable FLAG1
        FLAG1 = .FALSE.


        ! Iterative process to solve for the pressure distribution

        ! --------------- BEGIN ITERATIVE PROCESS ---------------
        
        DO WHILE(FLAG1 .EQV. (.FALSE.))
        
            ! Sweep all rows j=2,3,...,J-1
            do J=2,NROM1
            
                do I=2,NTHM1
                    A(I-1) = (1./DTH-1.5/H(I,J,M)*DHDTH(I,J,M))/(DTH*(RO(J)*SA)**2)
                    B(I-1) = -2*(1./DRO2+1./(DTH*RO(J)*SA)**2)
                    C(I-1) = (1/DTH+1.5/H(I,J,M)*DHDTH(I,J,M))/(DTH*(RO(J)*SA)**2)
                    D(I-1) = (2.0*NDR*DHDT(I,J,M))/H(I,J,M)**3           &
                             -NDR*DHDTH(I,J,M)/H(I,J,M)**3               &
                             -(1./DRO)*(1./DRO-1.5/H(I,J,M)*DHDRO(I,J,M) &
                             -0.5/RO(J))*P(I,J-1,M)                      &
                             -(1./DRO)*(1./DRO+1.5/H(I,J,M)*DHDRO(I,J,M) &
                             +0.5/RO(J))*P(I,J+1,M)
                end do

                ! Use Thomas algorithm to solve tridiagonal system of algebraic equations  
                CALL THOMAS(NTHM2,A,B,C,D,N1)
            
                ! Capture solution from Thomas and apply SOR factor, For
                do I=2,NTHM1
                    P(I,J,M)   = PPI(I,J,M)+FOR*(D(I-1)-PPI(I,J,M))
                    PPI(I,J,M) = P(I,J,M)
                end do
          
            end do
          
            ! Sweep all columns i=2,3,...,I-1
            DO I=2,NTHM1
                
                DO J=2,NROM1
                    A(J-1) = (1./DRO)*(1./DRO-1.5/H(I,J,M)*DHDRO(I,J,M)-0.5/RO(J))  
                    B(J-1) = -2*(1./DRO2+1./(DTH*RO(J)*SA)**2)
                    C(J-1) = (1./DRO)*(1./DRO+1.5/H(I,J,M)*DHDRO(I,J,M)+0.5/RO(J))
                    D(J-1) = (2.0*NDR*DHDT(I,J,M))/H(I,J,M)**3   &
                             -NDR*DHDTH(I,J,M)/H(I,J,M)**3       &
                             -(1./DTH-1.5/H(I,J,M)*DHDTH(I,J,M)) &
                             /(DTH*(RO(J)*SA)**2)*P(I-1,J,M)     &
                             -(1./DTH+1.5/H(I,J,M)*DHDTH(I,J,M)) &
                             /(DTH*(RO(J)*SA)**2)*P(I+1,J,M)
                end do
                
                ! Use Thomas algorithm to solve tridiagonal system of algebraic equations
                CALL THOMAS(NROM2,A,B,C,D,N1)
                
                ! Capture solution from Thomas and apply SOR factor, For
                do J=2,NROM1
                    P(I,J,M) = PPI(I,J,M)+FOR*(D(J-1)-PPI(I,J,M))
                end do
                
            end do
            
            ! Convergence test
            CALL CVERGE(NTH,NRO,P,PPI,FLAG1,M,N1,N2)
            
        end do       
        
        
        !  --------------- END ITERATIVE PROCESS ---------------
        
        
        ! Calculate forces (Fx, Fy and Fz) and moments (Mx and My) on 
        ! the impeller due to each separate blade
        do J=1,NRO
            do I=1,NTH
                DFX(I,J,M)=P(I,J,M)*RO(J)*COS(TH(I,M)+NDT)
                DFY(I,J,M)=P(I,J,M)*RO(J)*SIN(TH(I,M)+NDT) 
                DFZ(I,J,M)=P(I,J,M)*RO(J)
                DMX(I,J,M)=P(I,J,M)*RO(J)*(RO(J)-ROG)*SIN(TH(I,M)+NDT)
                DMY(I,J,M)=P(I,J,M)*RO(J)*(RO(J)-ROG)*COS(TH(I,M)+NDT)
            end do
        end do
        FX(M) = -SA*CA*DBLINT(DFX,DRO,DTH,NROM1,NTHM1,N1,N2,M)
        FY(M) = -SA*CA*DBLINT(DFY,DRO,DTH,NROM1,NTHM1,N1,N2,M)
        FZ(M) = -SA*SA*DBLINT(DFZ,DRO,DTH,NROM1,NTHM1,N1,N2,M)
        MX(M) = -SA*DBLINT(DMX,DRO,DTH,NROM1,NTHM1,N1,N2,M)
        MY(M) =  SA*DBLINT(DMY,DRO,DTH,NROM1,NTHM1,N1,N2,M)
        
        ! Calculate the total force (Fx, Fy and Fz) and moment (Mx and My) 
        ! on the centre of mass of the impeller
        FXT = FXT + FX(M)  
        FYT = FYT + FY(M)
        FZT = FZT + FZ(M)
        MXT = MXT + MX(M)
        MYT = MYT + MY(M)
        
    end do
    
    return
    end


    ! =========================================================================

      SUBROUTINE THOMAS(NXM2,A,B,C,D,N1)

!      This subroutine uses the Thomas algorithm for solving tri-diagonal
!      systems of algebraic equations with Dirichlet boundary conditions.  
!      It is used in iterative line-by-line sweeps to determine the pressure
!      distribution over the bearing surfaces. 

      REAL*8 A(N1),B(N1),C(N1),D(N1)
      INTEGER I, N1, NXM2

      B(1)=1.0/B(1)

      DO 5,I=2,NXM2
        B(I)=1.0/(B(I)-A(I)*B(I-1)*C(I-1))
        D(I)=D(I)-A(I)*B(I-1)*D(I-1)
    5 CONTINUE

      D(NXM2)=D(NXM2)*B(NXM2)

      DO 10,I=NXM2-1,1,-1
        D(I)=(D(I)-D(I+1)*C(I))*B(I)
   10 CONTINUE

      RETURN
      END


    ! =========================================================================
    
    SUBROUTINE CVERGE(NX,NY,P,PPI,FLAG,M,N1,N2)
    
    ! This subroutine tests for convergence of the solution. If the 
    ! residual is less than or equal to RESLIM then the solution is 
    ! converged.  If the residual is greater than RESMAX, the solution 
    ! is diverging and the program is stopped.  Otherwise, an over -
    ! relaxation factor is applied to speed up the iterative process 
    ! (This is performed in subroutines TB and CJB)
    
    REAL*8 P(N1,N1,N2), PPI(N1,N1,N2), RESLIM, RESID
    PARAMETER (RESLIM=1.0D-4, RESMAX=1.0D4)
    LOGICAL FLAG
    
    RESID=0.0D0
    DO J=1,NY
        DO I=1,NX
            RESID = MAX(RESID,ABS(P(I,J,M)-PPI(I,J,M)))
        end do
    end do
    
    IF (RESID .LE. RESLIM) THEN
        FLAG = .TRUE.
    ELSE IF (RESID .GE. RESMAX) THEN
        STOP 'Solution not converging'
    ELSE
        DO J=2,NY-1
            DO I=2,NX-1
                PPI(I,J,M) = P(I,J,M)          
            end do
        end do
    END IF
    
    RETURN
    END


    ! =========================================================================

    REAL*8 FUNCTION DBLINT(FXY,DX,DY,NXM1,NYM1,N1,N2,M)

    ! This subroutine performs a double integration using the
    ! trapezoidal rule of integration.

    REAL SUMM, FXY(N1,N1,N2)
    REAL*8 DX, DY
    INTEGER I, J, M, NXM1, NYM1

    SUMM = 0.0
    do J=1,NYM1
        do I=1,NXM1
            SUMM = SUMM+FXY(I,J,M)+FXY(I,J+1,M)+FXY(I+1,J,M)+FXY(I+1,J+1,M)
        end do
    end do

    DBLINT = 0.25*DX*DY*SUMM

    END


!~    ! =========================================================================
!~
!~      SUBROUTINE PLTTRA(ERX,ERXDT,ERY,ERYDT,ERZ,ERZDT,XROT,XROTDT,YROT,
!~     +                  YROTDT,DT,T,N,NT,FLAG1,FLAG2)
!~
!~!      This subroutine creates a Matlab script file that plots the
!~!      trajectory of the impeller.  It doubles as a results file.
!~
!~      PARAMETER(PI=3.141592654, RTOD=180.0/PI)
!~      REAL*8 ERX(NT), ERY(NT), ERZ(NT), XROT(NT), YROT(NT)
!~      REAL*8 ERXDT(NT), ERYDT(NT), ERZDT(NT)
!~      REAL*8 XROTDT(NT), YROTDT(NT), DT, T(NT)
!~      INTEGER N, NT
!~      LOGICAL FLAG1, FLAG2
!~
!~      IF (N .EQ. 1) THEN 
!~        WRITE(1,*)'% Where C refers to the column no. of matrix data:'
!~        WRITE(1,*)'% C1:  Non-dim. time'
!~        WRITE(1,*)'% C2:  Eccentricity ratio in the X-direction'
!~        WRITE(1,*)'% C3:  Translational velocity in the X-direction'
!~        WRITE(1,*)'% C4:  Eccentricity ratio in the Y-direction'
!~        WRITE(1,*)'% C5:  Translational velocity in the Y-direction'
!~        WRITE(1,*)'% C6:  Eccentricity ratio in the Z-direction'
!~        WRITE(1,*)'% C7:  Translational velocity in the X-direction'
!~        WRITE(1,*)'% C8:  Non-dim. rotation about the X-axis'
!~        WRITE(1,*)'% C9:  Non-dim. rotational velocity about the X-axis'
!~        WRITE(1,*)'% C10: Non-dim. rotation about the X-axis'
!~        WRITE(1,*)'% C11: Non-dim. rotational velocity about the Y-axis'
!~        WRITE(1,*)
!~        WRITE(1,*) 'data=['
!~      END IF
!~
!~      WRITE(1,5) T(N)*RTOD,ERX(N),ERXDT(N),ERY(N),ERYDT(N),ERZ(N), 
!~     +           ERZDT(N),XROT(N),XROTDT(N),YROT(N),YROTDT(N)
!~
!~      IF ((N .EQ. NT) .OR. (FLAG1 .EQV. (.TRUE.))
!~     +                  .OR. (FLAG2 .EQV. (.TRUE.))) THEN
!~        WRITE(1,*) '];'
!~        WRITE(1,*)
!~        WRITE(1,*) 'sf=', 2*PI/DT, ';'
!~        WRITE(1,*) 'K=', K, ';'
!~        WRITE(1,*) 'erx=data(:,2);' 
!~        WRITE(1,*) 'ery=data(:,4);'
!~        WRITE(1,*) 'erz=data(:,6);'
!~        WRITE(1,*) 'xrot=data(:,8);'
!~        WRITE(1,*) 'yrot=data(:,10);'
!~!        A correction is made to the eccentricities in the X and Y
!~!        directions where the erz term below accounts for increased / 
!~!        decreased clearance as the impeller moves down/up inside the
!~!        cone shaped upper pump housing.
!~        WRITE(1,*) 'erx=erx./(1.0-erz);'
!~        WRITE(1,*) 'ery=ery./(1.0-erz);'
!~      END IF
!~    5 FORMAT(F10.3,10F12.6)
!~
!~      RETURN
!~      END
!~
!~      
!~!  ======================================================================
!~
!~      SUBROUTINE PLTCJB(P,N,NTH,NRO,NP,NT,SA,RO,TH,NDT,N1,N2)
!~
!~!      This subroutine plots the pressure distribution over the conical
!~!      journal bearing of the impeller at each time step.  A Matlab movie
!~!      is created by playing these plots in the correct sequence.  To run
!~!      this movie load and run the Matlab script file "Orbit3D.m" and
!~!      click the appropriate selection on the pop-up menu.
!~
!~      REAL*8 P(N1,N1,N2), TH(N1,N2), RO(N1), NDT, SA
!~      INTEGER I, J, M, N, NTH, NRO, NP, NT 
!~
!~      DO 65,M=1,NP
!~
!~        IF ((N .EQ. 1) .AND. (M .EQ. 1)) THEN
!~          WRITE(2,*)  'dth=2*pi/100;'
!~          WRITE(2,*)  'theta=0:dth:2*pi;'
!~          WRITE(2,40) 'xro1=', RO(1),   '*cos(theta);'
!~          WRITE(2,40) 'yro1=', RO(1),   '*sin(theta);'
!~          WRITE(2,40) 'xro2=', RO(NRO), '*cos(theta);'
!~          WRITE(2,40) 'yro2=', RO(NRO), '*sin(theta);'
!~          WRITE(2,*)  'MOV1=moviein(', NT, ');'
!~          WRITE(2,*)
!~        END IF
!~
!~        IF (N .LT. 10) THEN
!~          WRITE(2,10) 'xcjb', N, M, '=['
!~        ELSE IF (N .LT. 100) THEN
!~          WRITE(2,15) 'xcjb', N, M, '=['
!~        ELSE
!~          WRITE(2,20) 'xcjb', N, M, '=['
!~        END IF
!~        DO 5,J=1,NRO
!~          WRITE(2,25) ((RO(J)*COS(TH(I,M)+NDT)),I=1,NTH)
!~    5   CONTINUE
!~   10   FORMAT(A4,I1,I1,A2)
!~   15   FORMAT(A4,I2,I1,A2)
!~   20   FORMAT(A4,I3,I1,A2)
!~   25   FORMAT(101F10.5)
!~        WRITE(2,*) '];'
!~        WRITE(2,*)
!~
!~        IF (N .LT. 10) THEN
!~          WRITE(2,10) 'ycjb', N, M, '=['
!~        ELSE IF (N .LT. 100) THEN
!~          WRITE(2,15) 'ycjb', N, M, '=['
!~        ELSE
!~          WRITE(2,20) 'ycjb', N, M, '=['
!~        END IF
!~        DO 30,J=1,NRO
!~          WRITE(2,25) ((RO(J)*SIN(TH(I,M)+NDT)),I=1,NTH)
!~   30   CONTINUE
!~        WRITE(2,*) '];'
!~        WRITE(2,*)
!~
!~        IF (N .LT. 10) THEN
!~          WRITE(2,10) 'pcjb', N, M, '=['
!~        ELSE IF (N .LT. 100) THEN
!~          WRITE(2,15) 'pcjb', N, M, '=['
!~        ELSE
!~          WRITE(2,20) 'pcjb', N, M, '=['
!~        END IF
!~        DO 35,J=1,NRO
!~          WRITE(2,25) (P(I,J,M),I=1,NTH)
!~   35   CONTINUE
!~        WRITE(2,*) '];'
!~        WRITE(2,*)     
!~
!~        IF (N .LT. 10) THEN
!~          WRITE(2,45) 'surf(xcjb',N,M,',ycjb',N,M,',pcjb',N,M,');'
!~        ELSE IF (N .LT. 100) THEN
!~          WRITE(2,50) 'surf(xcjb',N,M,',ycjb',N,M,',pcjb',N,M,');' 
!~        ELSE
!~          WRITE(2,55) 'surf(xcjb',N,M,',ycjb',N,M,',pcjb',N,M,');' 
!~        END IF
!~        IF (M .EQ. 1) THEN           
!~          WRITE(2,*) 'hold on;'
!~!          WRITE(2,60) 'axis([',-1./SA,1./SA,-1./SA,1./SA,0,0.5,']);'
!~          WRITE(2,*) 'axis([',-1./SA,1./SA,-1./SA,1./SA,0,0.5,']);'
!~          WRITE(2,*) 'set(gca,''PlotBoxAspectRatio'',[2 2 1]);'
!~          WRITE(2,*) 'view(25,30);'
!~          WRITE(2,*) 'plot(xro1,yro1,''--'');'
!~          WRITE(2,*) 'plot(xro2,yro2,''--'');' 
!~          WRITE(2,*) 'xlabel(''X / sin\ alpha'');'   
!~          WRITE(2,*) 'ylabel(''Y / sin\ alpha'');'
!~          WRITE(2,*) 'zlabel(''Pressure, p^*'');'
!~          WRITE(2,*)
!~        ELSE IF (M .EQ. NP) THEN
!~          WRITE(2,*) 'shading interp;'
!~          WRITE(2,*) 'hold off;'
!~          WRITE(2,*) 'MOV1(', N ,') = getframe;'
!~          WRITE(2,*)
!~        END IF
!~   40   FORMAT(A5,F6.4,A12)
!~   45   FORMAT(A9,2I1,A5,2I1,A5,2I1,A2)
!~   50   FORMAT(A9,I2,I1,A5,I2,I1,A5,I2,I1,A2)
!~   55   FORMAT(A9,I3,I1,A5,I3,I1,A5,I3,I1,A2)
!~!   60   FORMAT(A6,6F10.6,A3)
!~
!~   65 CONTINUE
!~
!~      RETURN    
!~      END 
!~
!~
!~!  ======================================================================
!~
!~      SUBROUTINE PLTTB(P,N,NTH,NR,NP,NT,R,TH,NDT,N1,N2)
!~
!~!      This subroutine plots the pressure distribution over the thrust
!~!      bearing of the VentrAssist impeller at each time step.  A Matlab 
!~!      movie is created by playing these plots in the correct sequence.
!~!      To run this movie load and run the Matlab script file "Orbit3D.m" 
!~!      and click the appropriate selection on the pop-up menu.
!~
!~      REAL*8 P(N1,N1,N2), TH(N1,N2), R(N1), NDT
!~      INTEGER I, J, M, N, NTH, NR, NP, NT 
!~
!~      DO 60,M=1,NP
!~
!~        IF ((N .EQ. 1) .AND. (M .EQ. 1)) THEN
!~          WRITE(3,*) 'dth=2*pi/100;'
!~          WRITE(3,*) 'theta=0:dth:2*pi;'
!~          WRITE(3,40) 'xr1=', R(1),  '*cos(theta);'
!~          WRITE(3,40) 'yr1=', R(1),  '*sin(theta);'
!~          WRITE(3,40) 'xr2=', R(NR), '*cos(theta);'
!~          WRITE(3,40) 'yr2=', R(NR), '*sin(theta);'
!~          WRITE(3,*)  'MOV2=moviein(', NT, ');'
!~          WRITE(3,*)
!~        END IF
!~
!~        IF (N .LT. 10) THEN
!~          WRITE(3,10) 'xtb', N, M, '=['
!~        ELSE IF (N .LT. 100) THEN
!~          WRITE(3,15) 'xtb', N, M, '=['
!~        ELSE
!~          WRITE(3,20) 'xtb', N, M, '=['
!~        END IF
!~        DO 5,J=1,NR
!~          WRITE(3,25) ((R(J)*COS(TH(I,M)+NDT)),I=1,NTH)
!~    5   CONTINUE
!~   10   FORMAT(A4,I1,I1,A2)
!~   15   FORMAT(A4,I2,I1,A2)
!~   20   FORMAT(A4,I3,I1,A2)
!~   25   FORMAT(101F10.5)
!~        WRITE(3,*) '];'
!~        WRITE(3,*)
!~
!~        IF (N .LT. 10) THEN
!~          WRITE(3,10) 'ytb', N, M, '=['
!~        ELSE IF (N .LT. 100) THEN
!~          WRITE(3,15) 'ytb', N, M, '=['
!~        ELSE
!~          WRITE(3,20) 'ytb', N, M, '=['
!~        END IF
!~        DO 30,J=1,NR
!~          WRITE(3,25) ((R(J)*SIN(TH(I,M)+NDT)),I=1,NTH)
!~   30   CONTINUE
!~        WRITE(3,*) '];'
!~        WRITE(3,*)
!~
!~        IF (N .LT. 10) THEN
!~          WRITE(3,10) 'ptb', N, M, '=['
!~        ELSE IF (N .LT. 100) THEN
!~          WRITE(3,15) 'ptb', N, M, '=['
!~        ELSE
!~          WRITE(3,20) 'ptb', N, M, '=['
!~        END IF
!~        DO 35,J=1,NR
!~          WRITE(3,25) (P(I,J,M),I=1,NTH)
!~   35   CONTINUE
!~        WRITE(3,*) '];'
!~        WRITE(3,*)     
!~     
!~        IF (N .LT. 10) THEN
!~          WRITE(3,45) 'surf(xtb',N,M,',ytb',N,M,',ptb',N,M,');'
!~        ELSE IF (N .LT. 100) THEN
!~          WRITE(3,50) 'surf(xtb',N,M,',ytb',N,M,',ptb',N,M,');' 
!~        ELSE
!~          WRITE(3,55) 'surf(xtb',N,M,',ytb',N,M,',ptb',N,M,');' 
!~        END IF
!~        IF (M .EQ. 1) THEN
!~          WRITE(3,*) 'hold on;'
!~          WRITE(3,*) 'axis([-1 1 -1 1 0 0.5]);'
!~          WRITE(3,*) 'set(gca,''PlotBoxAspectRatio'',[2 2 1]);'
!~          WRITE(3,*) 'view(25,30);'
!~          WRITE(3,*) 'plot(xr1,yr1,''--'');'
!~          WRITE(3,*) 'plot(xr2,yr2,''--'');'
!~          WRITE(3,*) 'xlabel(''X'');'
!~          WRITE(3,*) 'ylabel(''Y'');'
!~          WRITE(3,*) 'ylabel(''Y'');'
!~          WRITE(3,*) 'zlabel(''Pressure, p^*'');'
!~          WRITE(3,*)
!~        ELSE IF (M .EQ. NP) THEN
!~          WRITE(3,*) 'shading interp;'
!~          WRITE(3,*) 'hold off;'
!~          WRITE(3,*) 'MOV2(', N ,') = getframe;'
!~          WRITE(3,*)
!~        END IF
!~   40   FORMAT(A4,F6.4,A12)
!~   45   FORMAT(A9,2I1,A5,2I1,A5,2I1,A2)
!~   50   FORMAT(A9,I2,I1,A5,I2,I1,A5,I2,I1,A2)
!~   55   FORMAT(A9,I3,I1,A5,I3,I1,A5,I3,I1,A2)
!~
!~   60 CONTINUE
!~
!~      RETURN    
!~      END
      
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      
      END
      