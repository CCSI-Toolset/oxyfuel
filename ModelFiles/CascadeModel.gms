****************************
***** CascadeModel.gms *****
****************************

* This file is home for the Edmister-based group method distillation cascade
* model. Flowsheet topology is defined in a separate file.

*--------------------------------**
* Variables for Cascade
*--------------------------------**
Positive Variables
   PhiAEd1(Ed1,Comp)          Phi for abosrption
   PhiSEd1(Ed1,Comp)          Phi for Stripping
   VNEd1(Ed1)                 Vn stream flowrate
   AeEd1(Ed1,Comp)            Average absorption factor
   SeEd1(Ed1,Comp)            Average stripping factore
   L1Ed1(Ed1)                 L1 stream flowrate
   ANEd1(Ed1,Comp)
   A1Ed1(Ed1,Comp)
   SNEd1(Ed1,Comp)
   S1Ed1(Ed1,Comp)
   StgEd1(Ed1)                Number of trays in each cascade
   DummyAe(Ed1,Comp)
   DummySe(Ed1,Comp)
   CascInter(Ed1,Comp)        Intermediate (prevents log of zero or negative)
   CascInterA(Ed1,Comp)
   CascInterS(Ed1,Comp)
   PresDrop(Ed1)
;

Positive Variables
   RefRatio(Clmn)
   RebRatio(Clmn);

Variables
   CscSlack(Ed1,Comp)        Slack variable for performance equations;

*----------------------------------------------------*
*        Equations for Cascade
*----------------------------------------------------*

Equations
EqTotCompMolBalEd1(Ed1,Comp) Overall component mass balance
EqTotEgyBalEd1(Ed1) Overall energy balance
EqCascEd1(Ed1,Comp) Cascade (group) equation
EqCascEd1WithSlack(Ed1,Comp) Cascade (group) equation

EqApprx1Ed1(Ed1,Str,Str2) Approximation for Vn vapor flow
EqApprx2Ed1(Ed1,Comp,Str,Str2,Str3) Approximation for Vn vapor flow

* Obsolete for now
EqCascEdRevised(Ed1,Comp) Revised cascade equation
EqCascInter(Ed1,Comp,Str,Str2,Str3) Intermediate for revised cascade equation
EqPhiAApprox(Ed1,Comp)
EqPhiSApprox(Ed1,Comp)

* Original Model
EqPhiAEd1Orig(Ed1,Comp)  Phi for absorption
EqPhiSEd1Orig(Ed1,Comp)  Phi for stripping

* Alternate 1
EqPhiAEd1Alt1(Ed1,Comp)  Phi for absorption
EqPhiSEd1Alt1(Ed1,Comp)  Phi for stripping
EqCascInterA(Ed1,Comp) Cascade intermediate
EqCascInterS(Ed1,Comp) Cascade intermediate

* Alternate 2
EqPhiAEd1Alt2(Ed1,Comp)  Phi for absorption
EqPhiSEd1Alt2(Ed1,Comp)  Phi for stripping
EqAeEd1Extra(Ed1,Comp)
EqSeEd1Extra(Ed1,Comp)

* Alternate 3
EqCascEd1Alt3(Ed1,Comp)  Alternate performance eqn
EqPhiAEd1Alt3(Ed1,Comp)  Phi* for absorption
EqPhiSEd1Alt3(Ed1,Comp)  Phi* for stripping

EqAbsorbingSafe(Ed1,Comp)
EqStrippingSafe(Ed1,Comp)
EqAeEd1(Ed1,Comp)    Effective Abs
EqSeEd1(Ed1,Comp)    Effective Stp
EqA1Ed1Smp(Ed1,Comp,Str) Absorption factor at top
EqANEd1Smp(Ed1,Comp,Str) Absorption factor at bottom


EqA1Ed1(Ed1,Comp,Str,Str2,Str3) Absorption factor at top
EqANEd1(Ed1,Comp,Str,Str2,Str3) Absorption factor at bottom
EqS1Ed1(Ed1,Comp)
EqSNEd1(Ed1,Comp)


EqCscNoPresDropV(Ed1,Str,Str2)
EqCscNoPresDropL(Ed1,Str,Str2)
EqCscConstPres(Ed1,Str,Str2)

EqCalcPresDrop(Ed1)
EqCscPresDropV(Ed1,Str,Str2)
EqCscPresDropL(Ed1,Str,Str2)
EqCscPresDropOverall(Ed1,Str,Str2)

EqRefRatio(Clmn,Str,Str2)
EqRebRatio(Clmn,Str,Str2)
;


EqTotCompMolBalEd1(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  SUM(Str $InVEd1(Ed1,Str), Fc(Str,Comp)) + SUM(Str $InLEd1(Ed1,Str), Fc(Str,Comp)) =E=
  SUM(Str $OutVEd1(Ed1,Str), Fc(Str,Comp)) + SUM(Str $OutLEd1(Ed1,Str), Fc(Str,Comp));

EqTotEgyBalEd1(Ed1)$(NOT InactiveCsc(Ed1))..
    SUM(Str $InVEd1(Ed1,Str), F(Str)*H(Str))
  + SUM(Str $InLEd1(Ed1,Str), F(Str)*H(Str))  =E=
    SUM(Str $OutVEd1(Ed1,Str), F(Str)*H(Str))
  + SUM(Str $OutLEd1(Ed1,Str), F(Str)*H(Str)) ;

EqCascEd1(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  SUM(Str $OutVEd1(Ed1,Str), Fc(Str,Comp)) =E= SUM(Str $InVEd1(Ed1,Str), Fc(Str,Comp))*PhiAEd1(Ed1,Comp)
  + SUM(Str $InLEd1(Ed1,Str), Fc(Str,Comp))*(1- PhiSEd1(Ed1,Comp)) ;

EqCascEd1WithSlack(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  SUM(Str $OutVEd1(Ed1,Str), Fc(Str,Comp)) =E= SUM(Str $InVEd1(Ed1,Str), Fc(Str,Comp))*PhiAEd1(Ed1,Comp)
  + SUM(Str $InLEd1(Ed1,Str), Fc(Str,Comp))*(1- PhiSEd1(Ed1,Comp)) + CscSlack(Ed1,Comp) ;

EqApprx1Ed1(Ed1,Str,Str2)$(OutLEd1(Ed1,Str) AND OutVEd1(Ed1,Str2) AND NOT InactiveCsc(Ed1))..
   L1Ed1(Ed1) - F(Str) =E=  F(Str2) - VNEd1(Ed1) ;

Sets
  CfA(Comp)      Component for inner cascade approximation             /O2/;

EqApprx2Ed1(Ed1,CfA,Str,Str2,Str3)$(OutLEd1(Ed1,Str) AND OutVEd1(Ed1,Str2) AND InLEd1(Ed1,Str3))..
  VNEd1(Ed1)*Xc(Str,CfA)*(PhiAEd1(Ed1,CfA)/(1 - PhiAEd1(Ed1,CfA)*(DummyAe(Ed1,CfA)*(1/AeEd1(Ed1,CfA) + 1/Power(AeEd1(Ed1,CfA),2)))))
  + L1Ed1(Ed1)*Xc(Str2,CfA)*(1 - PhiSEd1(Ed1,CfA)/(1 - PhiSEd1(Ed1,CfA)*(DummySe(Ed1,CfA)*(1/SeEd1(Ed1,CfA) + 1/Power(SeEd1(Ed1,CfA),2)))))
  =e= L1Ed1(Ed1)*Xc(Str2,CfA) + Fc(Str2,CfA) - Fc(Str2,CfA);

*EqApprx1Ed1(Ed1,Str,Str2)$(InVEd1(Ed1,Str) AND OutVEd1(Ed1,Str2))..
*   VNEd1(Ed1) =E= F(Str)*(F(Str2)/F(Str))**(1/StgEd1(Ed1)) ;


***** Obsolete Equations *****

EqCascEdRevised(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  StgEd1(Ed1)*log(AeEd1(Ed1,Comp)) =e= log(CascInter(Ed1,Comp));

EqCascInter(Ed1,Comp,Str,Str2,Str3)$(InVEd1(Ed1,Str) AND OutVEd1(Ed1,Str2) AND InLEd1(Ed1,Str3) AND NOT InactiveCsc(Ed1))..
  CascInter(Ed1,Comp)*(Fc(Str3,Comp) - Fc(Str2,Comp)* AeEd1(Ed1,Comp)) =e=
         Fc(Str,Comp)*(1 - AeEd1(Ed1,Comp)) - Fc(Str2,Comp)*Fc(Str3,Comp) ;

*EqCascEdRevised(Ed1,Comp,Str,Str2,Str3)$(InVEd1(Ed1,Str) AND OutVEd1(Ed1,Str2) AND InLEd1(Ed1,Str3))..
*  exp(CascInter(Ed1,Comp))*(Fc(Str3,Comp) - Fc(Str2,Comp)* AeEd1(Ed1,Comp))
*         =e= Fc(Str,Comp)*(1 - AeEd1(Ed1,Comp)) - Fc(Str2,Comp)*Fc(Str3,Comp);

*EqCascInter(Ed1,Comp,Str,Str2,Str3)$(InVEd1(Ed1,Str) AND OutVEd1(Ed1,Str2) AND InLEd1(Ed1,Str3))..
*  CascInter(Ed1,Comp) =e= StgEd1(Ed1)*log(AeEd1(Ed1,Comp)) ;

CascInter.lo(Ed1,Comp) = 1E-8;
CascInter.up(Ed1,Comp) = 1E6;

EqPhiAApprox(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  2*log(PhiAEd1(Ed1,Comp)) + 2*(StgEd1(Ed1)+1)*log(AeEd1(Ed1,Comp)) =E= log(Power(AeEd1(Ed1,Comp) - 1,2)) ;

EqPhiSApprox(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  2*log(PhiSEd1(Ed1,Comp)) + 2*(StgEd1(Ed1)+1)*log(SeEd1(Ed1,Comp)) =E= log(Power(SeEd1(Ed1,Comp) - 1,2)) ;

***** Original Model *****

EqPhiAEd1Orig(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  PhiAEd1(Ed1,Comp)*(AeEd1(Ed1,Comp)**(StgEd1(Ed1)+1) - 1) =E= AeEd1(Ed1,Comp) - 1 ;
*  PhiAEd1(Ed1,Comp) =E= (AeEd1(Ed1,Comp) - 1)/(AeEd1(Ed1,Comp)**(StgEd1(Ed1)+1) - 1) ;

EqPhiSEd1Orig(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  PhiSEd1(Ed1,Comp)*(SeEd1(Ed1,Comp)**(StgEd1(Ed1)+1) - 1) =E= SeEd1(Ed1,Comp) - 1 ;
*  PhiSEd1(Ed1,Comp) =E= (SeEd1(Ed1,Comp) - 1)/(SeEd1(Ed1,Comp)**(StgEd1(Ed1)+1) - 1) ;

Model CascOrigEq /EqCascEd1, EqPhiAEd1Orig, EqPhiSEd1Orig/;

***** Alternate 1 *****

EqPhiAEd1Alt1(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  (StgEd1(Ed1)+1)*log(AeEd1(Ed1,Comp)) =e= log(CascInterA(Ed1,Comp)) - log(PhiAEd1(Ed1,Comp));

EqPhiSEd1Alt1(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  (StgEd1(Ed1)+1)*log(SeEd1(Ed1,Comp)) =e= log(CascInterS(Ed1,Comp)) - log(PhiSEd1(Ed1,Comp));

EqCascInterA(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  CascInterA(Ed1,Comp) =e= AeEd1(Ed1,Comp) - 1 + PhiAEd1(Ed1,Comp);

CascInterA.lo(Ed1,Comp) = 1E-6;
CascInterA.up(Ed1,Comp) = 1E8;

EqCascInterS(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  CascInterS(Ed1,Comp) =e= SeEd1(Ed1,Comp) - 1 + PhiSEd1(Ed1,Comp);

CascInterS.lo(Ed1,Comp) = 1E-6;
CascInterS.up(Ed1,Comp) = 1E8;

Model CascAlt1Eq /EqCascEd1, EqPhiAEd1Alt1, EqPhiSEd1Alt1, EqCascInterA, EqCascInterS/;

***** Alternate 2 *****

EqPhiAEd1Alt2(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  PhiAEd1(Ed1,Comp)*(DummyAe(Ed1,Comp) - 1) =E= AeEd1(Ed1,Comp) - 1 ;

EqPhiSEd1Alt2(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  PhiSEd1(Ed1,Comp)*(DummySe(Ed1,Comp) - 1) =E= SeEd1(Ed1,Comp) - 1 ;

EqAeEd1Extra(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  log(DummyAe(Ed1,Comp)) =e= (StgEd1(Ed1)+1)*log(AeEd1(Ed1,Comp));

DummyAe.lo(Ed1,Comp) = 1E-6;

EqSeEd1Extra(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  log(DummySe(Ed1,Comp)) =e= (StgEd1(Ed1)+1)*log(SeEd1(Ed1,Comp));

DummySe.lo(Ed1,Comp) = 1E-6;

Model CascAlt2Eq /EqCascEd1, EqPhiAEd1Alt2, EqPhiSEd1Alt2, EqAeEd1Extra, EqSeEd1Extra/;

***** Alternate 3 *****
EqPhiAEd1Alt3(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  (StgEd1(Ed1)+1)*log(AeEd1(Ed1,Comp)) =e= log(PhiAEd1(Ed1,Comp)*(AeEd1(Ed1,Comp) - 1) + 1);

EqPhiSEd1Alt3(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  (StgEd1(Ed1)+1)*log(SeEd1(Ed1,Comp)) =e= log(PhiSEd1(Ed1,Comp)*(SeEd1(Ed1,Comp) - 1) + 1);

EqCascEd1Alt3(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  SUM(Str $OutVEd1(Ed1,Str), Fc(Str,Comp)) =E= SUM(Str $InVEd1(Ed1,Str), Fc(Str,Comp))/PhiAEd1(Ed1,Comp)
  + SUM(Str $InLEd1(Ed1,Str), Fc(Str,Comp))/(1 - PhiSEd1(Ed1,Comp)) ;

Model CascAlt3Eq /EqCascEd1Alt3, EqPhiAEd1Alt3, EqPhiSEd1Alt3/;

EqAbsorbingSafe(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  (StgEd1(Ed1)+1)*log(AeEd1(Ed1,Comp)) =l= 25;

EqStrippingSafe(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  (StgEd1(Ed1)+1)*log(SeEd1(Ed1,Comp)) =l= 25;

EqS1Ed1(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  S1Ed1(Ed1,Comp)*A1Ed1(Ed1,Comp) =e= 1;

EqSNEd1(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  SNEd1(Ed1,Comp)*ANEd1(Ed1,Comp) =e= 1;

EqAeEd1(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  AeEd1(Ed1,Comp) =E=  (ANEd1(Ed1,Comp)*(A1Ed1(Ed1,Comp)+1) + 0.25)**0.5 - 0.5 ;

EqSeEd1(Ed1,Comp)$(NOT InactiveCsc(Ed1))..
  SeEd1(Ed1,Comp) =E=  (S1Ed1(Ed1,Comp)*(SNEd1(Ed1,Comp)+1) + 0.25)**0.5 - 0.5 ;

***** Simple Thermo Model Specific Equations *****

EqA1Ed1Smp(Ed1,Comp,Str)$(OutVEd1(Ed1,Str) AND NOT InactiveCsc(Ed1))..
  A1Ed1(Ed1,Comp)* F(Str) * Pvap(Str,Comp) =E=  L1Ed1(Ed1) * P(Str) ;

EqANEd1Smp(Ed1,Comp,Str)$(OutLEd1(Ed1,Str) AND NOT InactiveCsc(Ed1))..
  ANEd1(Ed1,Comp)* VNEd1(Ed1) * Pvap(Str,Comp) =E=  F(Str) * P(Str) ;

***** CEOS Thermo Model Specific Equations *****

EqA1Ed1(Ed1,j,Str,Str2,Str3)$(OutVEd1(Ed1,Str) and ComMap(Str,Str3,Str2))..
  A1Ed1(Ed1,j)* F(Str)*PhiEOS(Str2,j) =E=  L1Ed1(Ed1) * PhiEOS(Str,j) ;

EqANEd1(Ed1,j,Str,Str2,Str3)$(OutLEd1(Ed1,Str) and ComMap(Str,Str2,Str3))..
  ANEd1(Ed1,j)* VNEd1(Ed1)*PhiEOS(Str,j) =E=  F(Str)*PhiEOS(Str2,j)    ;


***** No Pressure Drop *****

EqCscNoPresDropV(Ed1,Str,Str2)$(OutVEd1(Ed1,Str) AND InVEd1(Ed1,Str2))..
  P(Str) =e= P(Str2);

EqCscNoPresDropL(Ed1,Str,Str2)$(OutLEd1(Ed1,Str) AND InLEd1(Ed1,Str2))..
  P(Str) =e= P(Str2);

EqCscConstPres(Ed1,Str,Str2)$(InVEd1(Ed1,Str) AND InLEd1(Ed1,Str2))..
  P(Str) =e= P(Str2);

***** Constant Pressure Drop per Tray *****
EqCalcPresDrop(Ed1)..
  PresDrop(Ed1) =g= StgEd1(Ed1)*PresDropCnst(Ed1);

EqCscPresDropV(Ed1,Str,Str2)$(OutVEd1(Ed1,Str) AND InVEd1(Ed1,Str2))..
  P(Str) + PresDrop(Ed1) =e= P(Str2);

EqCscPresDropL(Ed1,Str,Str2)$(OutLEd1(Ed1,Str) AND InLEd1(Ed1,Str2))..
  P(Str) =e= P(Str2) + PresDrop(Ed1);

EqCscPresDropOverall(Ed1,Str,Str2)$(InVEd1(Ed1,Str) AND InLEd1(Ed1,Str2))..
  P(Str) =e= P(Str2) + PresDrop(Ed1);

***** Reboil and Reflux Ration Calculations *****

EqRefRatio(Clmn,Str,Str2)$(ClmnRflx(Clmn,Str) AND ClmnDst(Clmn,Str2))..
  RefRatio(Clmn) * F(Str2) =e= F(Str);

EqRebRatio(Clmn,Str,Str2)$(ClmnRb(Clmn,Str) AND ClmnBtm(Clmn,Str2))..
  RebRatio(Clmn) * F(Str2) =e= F(Str);

Model CscPresDropEqns /EqCscNoPresDropV, EqCscNoPresDropL, EqCscConstPres/;
*Model CscPresDropEqns/ EqCalcPresDrop, EqCscPresDropV, EqCscPresDropL, EqCscPresDropOverall/;

Model EqCascade /EqTotCompMolBalEd1, EqTotEgyBalEd1, EqApprx1Ed1, CascAlt1Eq,
         EqS1Ed1, EqSNEd1, EqAeEd1, EqSeEd1, CscPresDropEqns  /;

* EqCscNoPresDropV, EqCscNoPresDropL, EqCscConstPres

* EqApprx1Ed1

* CascAlt1Eq

*Model EqCascade /EqTotCompMolBalEd1, EqTotEgyBalEd1, EqCascEdRevised, EqCascInter,
*      EqApprx1Ed1, EqAeEd1, EqCscNoPresDropV, EqCscNoPresDropL, EqCscConstPres /;

* Model EqCascade /EqTotCompMolBalEd1, EqTotEgyBalEd1, EqCascEd1, EqApprx1Ed1, EqCascInterA, EqCascInterS,
*         EqPhiAEd1, EqPhiSEd1, EqS1Ed1, EqSNEd1, EqAeEd1, EqSeEd1, EqCscNoPresDropV, EqCscNoPresDropL, EqCscConstPres /;
*         EqAeEd1Extra, EqSeEd1Extra/;

*Model EqCascade /EqTotCompMolBalEd1, EqTotEgyBalEd1, EqCascEd1, EqApprx1Ed1,
*         EqPhiAApprox, EqPhiSApprox, EqS1Ed1, EqSNEd1, EqAeEd1, EqSeEd1,
*         EqCscNoPresDropV, EqCscNoPresDropL, EqCscConstPres /;

Model EqCsdSimple /EqA1Ed1Smp, EqANEd1Smp/;

Model EqCsdCEOS /EqA1Ed1,EqANEd1/;

Model EqReduceCsdSimpe /EqTotCompMolBalEd1, EqTotEgyBalEd1, CscPresDropEqns/;

Model EqCscRegressionSimple /EqCascade - EqTotCompMolBalEd1 - EqTotEgyBalEd1 -  EqCascEd1, EqCascEd1WithSlack, EqA1Ed1Smp, EqANEd1Smp/;

Model EqCscRegressionCEOS /EqCascade - EqTotCompMolBalEd1 - EqTotEgyBalEd1 -  EqCascEd1, EqCascEd1WithSlack, EqA1Ed1,EqANEd1/;
