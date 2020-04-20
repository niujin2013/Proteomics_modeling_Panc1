**********************************************************************
C                           ADAPT                                     *
C                         Version 5                                   *
C**********************************************************************
C                                                                     *
C                           MODEL                                     *
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
C   Enter as Indicated                                                 C
C----c-----------------------------------------------------------------C

      NDEqs   =  9   ! Enter # of Diff. Eqs.
      NSParam =  12   ! Enter # of System Parameters.
      NVparam =  2  ! Enter # of Variance Parameters.
      NSecPar =  2   ! Enter # of Secondary Parameters.
      NSecOut =  0  ! Enter # of Secondary Outputs (not used).
      Ieqsol  =  1  ! Model type: 1 - DIFFEQ, 2 - AMAT, 3 - OUTPUT only.
      Descr   = ' both PTX and BRP on killing cells
     $            psi_b on p'

CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each System Parameter (eg. Psym(1)='Kel')         C
C----c-----------------------------------------------------------------C
      Psym(1) = 'R0'
      Psym(2) = 'Kg'
      psym(3) = 'Rss'
      Psym(4) = 'Emax_p'
      Psym(5) = 'EC50_p'
      Psym(6) = 'M_p'   
      Psym(7) = 'Psi_bp'      
      Psym(8) = 'Emax_b'
      Psym(9) = 'EC50_b'
      Psym(10) = 'M_b'
      psym(11)= 'ktau_p'
      psym(12)= 'ktau_b'

      
CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each Variance Parameter {eg: PVsym(1)='Sigma'}    C
C----c-----------------------------------------------------------------C
      PVsym(1)='Int'
      PVsym(2)='Slo'


CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each Secondary Parameter {eg: PSsym(1)='CLt'}     C
C----c-----------------------------------------------------------------C

      pssym(1) = 'tau_p'
      pssym(2) = 'tau_b'

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
        Real*8 R0,Kg,rss,Emax_P,EC50_P,M_P,ktau_P,Emax_b
        real*8 EC50_b,m_b, psi_bp,psi_pb,Ek_P,Ek_B,ktau_b
 

CC
C----------------------------------------------------------------------C
C   Enter Differential Equations Below  {e.g.  XP(1) = -P(1)*X(1) }    C
C----c-----------------------------------------------------------------C
      R0= P(1)
      Kg = P(2)
      Rss = P(3)
          
      Emax_p = P(4)
      EC50_p = P(5)
      M_p = P(6)
      
      psi_bp = P(7)
           
      Emax_b = P(8)
      EC50_b =P(9)
      M_b = P(10)
      
      ktau_p= p(11)
      ktau_b= p(12)
   
      
      IF (R(1).EQ.0.AND.R(2).EQ.0) THEN
      Ek_p = 0
      Ek_b = 0
      ELSE IF (R(1).NE.0.AND.R(2).EQ.0) THEN
      Ek_p = Emax_p * R(1)**M_p/ (EC50_p**M_p + R(1)**M_p)
      Ek_b = 0
      ELSE IF (R(1).EQ.0.AND.R(2).NE.0) THEN
      Ek_p = 0
      Ek_b = Emax_b * R(2)**M_b/ (EC50_b**M_b + R(2)**M_b)
      ELSE IF (R(1).NE.0.AND.R(2).NE.0) THEN
      Ek_p  =  Emax_p * R(1)**M_p/ ((psi_bp*EC50_p)**M_p + R(1)**M_p)
      Ek_b  =  Emax_b * R(2)**M_b/ ((psi_bp*EC50_b)**M_b + R(2)**M_b)
      ENDIF
      
C comb
      XP(1) = Kg * (X(1)+R0) !*(1-(R0+x(1))/Rss) 
     $  - (X(1)+R0) * (x(4) +x(8))

      XP(2) = Ktau_p * (Ek_p -X(2) )    !kp1
      XP(3) = Ktau_p * (X(2) - X(3) )    !k2
      XP(4) = Ktau_p * (X(3) - X(4) )    !k3
      !XP(5) = Ktau_p * (X(4) - X(5) )    !k4
      
      
      xp(6)=ktau_b*(Ek_b-x(6))   !k1_b
      xp(7)=ktau_b*(x(6)-x(7))     
      xp(8)=ktau_b*(x(7)-x(8))     
      !xp(9)=ktau_b*(x(8)-x(9))
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
        Real*8 R0

CC
C----------------------------------------------------------------------C
C   Enter Output Equations Below   {e.g.  Y(1) = X(1)/P(2) }           C
C----c-----------------------------------------------------------------C

       R0= P(1)
             
       Y(1) = X(1) + R0
    
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
C   Enter Variance Model Equations Below                               C
C         {e.g. V(1) = (PV(1) + PV(2)*Y(1))**2 }                       C
C----c-----------------------------------------------------------------C

      V(1) = (PV(1) + PV(2)*Y(1))**2
      

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
C   Enter Equations Defining Secondary Paramters                       C
C           {  e.g.  PS(1) = P(1)*P(2)   }                             C
C----c-----------------------------------------------------------------C
      
      PS(1) = 1/p(11)
      ps(2) = 1/p(12)
        
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
C   Enter non zero elements of state matrix  {e.g.  A(1,1) = -P(1) }   C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C