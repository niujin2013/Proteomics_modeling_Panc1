C*********************************************************************
C                           ADAPT                                     *
C                         Version 5                                   *
C**********************************************************************
C                                                                     *
C                      MODEL - 1COMPCL                                *
C                                                                     *
C    This file contains Fortran subroutines into which the user       *
C    must enter the relevant model equations and constants.           *
C    Consult the User's Guide for details concerning the format for   *
C    entered equations and definition of symbols.                     *
C                                                                     *
C       1. Symbol-  Parameter symbols and model constants             *
C       2. DiffEq-  System differential equations                     *
C       3. Output-  System output equations                           *
C       4. Varmod-  Error variance model equations                    *
C       5. Covmod-  Covariate model equations (ITS,MLEM)              *
C       6. Popinit- Population parameter initial values (ITS,MLEM)    *
C       7. Prior -  Parameter mean and covariance values (ID,NPD,STS) *
C       8. Sparam-  Secondary parameters                              *
C       9. Amat  -  System state matrix                               *
C                                                                     *
C**********************************************************************

C######################################################################C
        
        Subroutine SYMBOL
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

CC
C----------------------------------------------------------------------C
C       Enter as Indicated                                             C
C----c-----------------------------------------------------------------C

      NDEqs   =  6      ! Enter # of Diff. Eqs.
      NSParam =  15     ! Enter # of System Parameters.
      NVparam =  3      ! Enter # of Variance Model Parameters.
      NSecPar =  0      ! Enter # of Secondary Parameters.
      NSecOut =  0      ! Enter # of Secondary Outputs (not used).
      Ieqsol  =  1 ! Model type: 1 - DIFFEQ, 2 - AMAT, 3 - OUTPUT only.
      Descr = 'cell cycle,total apoptosis% and total cell#'
 
c fit  G0G1%, S%, G2M%, PL%, apop% and normalized cell# 
c k10 * (I0) 
c G2 Ma togather 
c polyploid: kpl only 
CC
C----------------------------------------------------------------------C
C     Enter Symbol for Each System Parameter (eg. Psym(1)='Kel')       C
C----c-----------------------------------------------------------------C
        PSym(1)  = 'k12'  
        PSym(2)  = 'k23'                  
        PSym(3)  = 'k31'
        PSym(4)  = 'LNmax'     ! 0). Gpz growth 
        PSym(5)  = 'kap'      ! natural apop rate 
        PSym(6)  = 'DT'       ! transition time before prone to death 
        PSym(7)  = 'kpl'  	   !   Polyploidy           
c BRP parameters  
        PSym(8)  = 'Smax_b_ap' ! 1). BRP induced apoptosis 
        PSym(9)  = 'SC50_b'
        PSym(10) = 'psi' ! no 
c PTX parameters        
        PSym(11) = 'Smax_p_ap' ! 3). PTX to induce apoptosis 
        PSym(12) = 'SC50_p'
        PSym(13) = 'gamma_p' ! no 
        PSym(14) = 'Kmax_p'    ! 5). PTX to induce mitotic arrest
        PSym(15) = 'KC50_p'  
        PSym(16) = ''


CC
C----------------------------------------------------------------------C
C    Enter Symbol for Each Variance Parameter {eg: PVsym(1)='Sigma'}   C
C----c-----------------------------------------------------------------C

       PVsym(1) = 'int_%'
       PVsym(2) = 'slp_#'
       PVsym(3) = 'int_#'
       
       
C----------------------------------------------------------------------C
C    Enter Symbol for Each Secondary Parameter {eg: PSsym(1)='CLt'}    C
C----c-----------------------------------------------------------------C
        PSsym(1) = 'SC50_b*psi'
        PSsym(2) = 'SC50_p*psi'

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine DIFFEQ(T,X,XP)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 T,X(MaxNDE),XP(MaxNDE)
        Real*8 k12,k23,k31,LNmax,kap,knc,kat,Smax_b_ap,SC50_b,Smax_b_nc
        Real*8 Smax_p_ap,SC50_p,Smax_p_nc,Kmax_p,KC50_p,gamma_p,kdec_ma
        Real*8 kpl,psi,live,live0,S_b_ap,S_b_nc,S_p_ap,S_p_nc
        Real*8 kar1,kma,I0
		Real*8 kap2,PTX,BRP,DT,kapm


CC
C----------------------------------------------------------------------C
C     Enter Differential Equations Below  {e.g.  XP(1) = -P(1)*X(1) }  C
C----c-----------------------------------------------------------------C
		k12	=	P(1)
		k23	=	P(2)
		k31	=	P(3)
		LNmax	=	P(4)
		kap		=	P(5)
		DT   	= 	P(6) ! transition time before prone to death   
		kpl		=	P(7) ! transition  to PL 		
		Smax_b_ap	=	P(8)
		SC50_b	=	P(9)
		psi 	=	P(10)			
		Smax_p_ap	=	P(11)
		SC50_p	=	P(12)
		gamma_p	=	P(13)
		Kmax_p	=	P(14)
		KC50_p	=	P(15)
		      
        live   = x(1) +x(2) +x(3)+ x(4)+ x(5)  !total live cells 
        live0  = IC(1)+IC(2)+IC(3)             !initial live cells 
		
		BRP = R(1) 
		PTX = R(2) 
      IF (PTX.EQ.0.AND.BRP.EQ.0) THEN
          S_b_ap = 0.0         
          S_p_ap = 0.0		 		  
		  
      ELSE IF (PTX.NE.0.AND.BRP.EQ.0) THEN
          S_b_ap = 0.0	
          S_p_ap = Smax_p_ap *PTX /( SC50_p+ PTX )         
	  
      ELSE IF (PTX.EQ.0.AND.BRP.NE.0) THEN
        S_b_ap = Smax_b_ap* BRP /( SC50_b+ BRP )      
        S_p_ap = 0.0  
        
      ELSE IF (PTX.NE.0.AND.BRP.NE.0) THEN
        S_b_ap = Smax_b_ap *BRP /( (SC50_b*psi)+ BRP ) 
        S_p_ap = Smax_p_ap *PTX /( (SC50_p)+ PTX )        
      ENDIF	  
  
      kar1   = 0		  
      IF (R(2).GT.(5) ) THEN           
          kar1     = Kmax_p *PTX** gamma_p /
     $          ( PTX** gamma_p +KC50_p** gamma_p)       
      ENDIF	  

      kma= kar1 ! * exp( -kdec_ma * T )
	  kapm = kap * ( 1+S_b_ap+(T>DT)*S_p_ap) 

c Growth inhibition      
      I0 = dlog( LNmax * live0 ) -dlog( live )        
c G1
      xp(1) = 2* k31 * x(3)
     $       - k12 *I0 * x(1)    
     $       -kap * (1 + S_b_ap + S_p_ap) * x(1) 
c S      
      xp(2) = k12 *x(1) *I0 - k23 *x(2) 
     $       -kap * (1 + S_b_ap + S_p_ap)   * x(2)
cG2M
      xp(3) = k23* x(2)
     $       - k31 *x(3)  
     $       - kma * x(3)
     $       -kap * (1 + S_b_ap + S_p_ap)   * x(3)   
c arrested M
      xp(4) = kma * x(3) - kpl* x(4) 
     $       -kapm * x(4)    
	 
c polyploidy      
      xp(5) = kpl * x(4)  -kpl * x(5)
      
c Apoptosis       
      xp(6) = kap * (1 + S_b_ap + S_p_ap) *( x(1)+x(2)+x(3))
     $       +kapm *x(4) -kap * x(6) 
      



C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine OUTPUT(Y,T,X)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 Y(MaxNOE),T,X(MaxNDE)
        Real*8 live0, kdec, I0, live,Nmax, Lmax,total 

CC
C----------------------------------------------------------------------C
C     Enter Output Equations Below   {e.g.  Y(1) = X(1)/P(2) }         C
C----c-----------------------------------------------------------------C
      total = x(1)+x(2)+x(3)+x(4)+x(5)+x(6)     !total  
      live = x(1)+x(2)+x(3) +x(4) +x(5)
      
      y(1) = (x(1) ) /live*100          !G0G1%
      y(2) = (x(2) ) /live*100          !S%
      y(3) = (x(3) +x(4)) /live*100     !G2M%
      y(4) = x(5)    /live*100          !poly %
	  y(5) = x(6) /total *100           !apop%
	  y(6) = total                      !total cell#
 
        !live0= IC(1)+IC(2)+IC(3)+IC(4)+IC(7)+IC(8)
        !Nmax = P(6)
        !Lmax = Nmax * live0
      
        !  I0 =log(Lmax)-log(live) ! Growth inhibition


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine VARMOD(V,T,X,Y)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 V(MaxNOE),T,X(MaxNDE),Y(MaxNOE)

CC
C----------------------------------------------------------------------C
C       Enter Variance Model Equations Below                           C
C        {e.g. V(1) = (PV(1) + PV(2)*Y(1))**2 }                        C
C----c-----------------------------------------------------------------C

        V(1)  = PV(1)
        V(2)  = PV(1)
        V(3)  = PV(1)    
        V(4)  = PV(1)
        V(5)  = PV(1)
        V(6)  = (PV(3) + PV(2)*Y(6))**2

       

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End
C######################################################################C

        Subroutine COVMOD(Pmean, ICmean, PC)
C  Defines any covariate model equations (MLEM, ITS)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 PC(MaxNCP)
        Real*8 Pmean(MaxNSP+MaxNDE), ICmean(MaxNDE)

CC
C----------------------------------------------------------------------C
C     Enter # of Covariate Parameters                                  C
C----c-----------------------------------------------------------------C

        NCparam = 0    ! Enter # of Covariate Parameters.

CC
C----------------------------------------------------------------------C
C   Enter Symbol for Covariate Params {eg: PCsym(1)='CLRenal'}         C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   For the Model Params. that Depend on Covariates Enter the Equation C
C         {e.g. Pmean(1) =  PC(1)*R(2) }                               C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine POPINIT(PmeanI,ICmeanI,PcovI,ICcovI, PCI)
C  Initial parameter values for population program parameters (ITS, MLEM)

        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Integer I,J
        Real*8 PmeanI(MaxNSP+MaxNDE), ICmeanI(MaxNDE)
        Real*8 PcovI(MaxNSP+MaxNDE,MaxNSP+MaxNDE), ICcovI(MaxNDE,MaxNDE)
        Real*8 PCI(MaxNCP)

CC
C----------------------------------------------------------------------C
C  Enter Initial Values for Population Means                           C
C          {  e.g. PmeanI(1) = 10.0    }                               C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   Enter Initial Values for Pop. Covariance Matrix (Lower Triang.)    C
C         {  e.g. PcovI(2,1) = 0.25    }                               C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   Enter Values for Covariate Model Parameters                        C
C         {  e.g. PCI(1) = 2.0    }                                    C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine PRIOR(Pmean,Pcov,ICmean,ICcov)
C  Parameter mean and covariance values for MAP estimation (ID,NPD,STS)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Integer I,J
        Real*8 Pmean(MaxNSP+MaxNDE), ICmean(MaxNDE)
        Real*8 Pcov(MaxNSP+MaxNDE,MaxNSP+MaxNDE), ICcov(MaxNDE,MaxNDE)

CC
C----------------------------------------------------------------------C
C  Enter Nonzero Elements of Prior Mean Vector                         C
C          {  e.g. Pmean(1) = 10.0    }                                C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   Enter Nonzero Elements of Covariance Matrix (Lower Triang.)       C
C         {  e.g. Pcov(2,1) = 0.25    }                                C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine SPARAM(PS,P,IC)
        Implicit None

        Include 'globals.inc'

        Real*8 PS(MaxNSECP), P(MaxNSP+MaxNDE), IC(MaxNDE) 

CC
C----------------------------------------------------------------------C
C       Enter Equations Defining Secondary Paramters                   C
C           {  e.g.  PS(1) = P(1)*P(2)   }                             C
C----c-----------------------------------------------------------------C
        PS(1) = P(19) * p(9)
        PS(2) = P(19) * P(12)  
              
        
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End
        
C######################################################################C
 
        Subroutine AMAT(A)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Integer I,J
        Real*8 A(MaxNDE,MaxNDE)

        DO I=1,Ndeqs
           Do J=1,Ndeqs
              A(I,J)=0.0D0
           End Do
        End Do

CC
C----------------------------------------------------------------------C
C    Enter non zero elements of state matrix  {e.g.  A(1,1) = -P(1) }  C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C
