*******************************
***** CEOSDewBubModel.gms *****
*******************************

* TO DO: Rewrite this introduction

* This file contains the equations required for dew and bubble point
* calculations using the cubic EOS framework. The Edmister cascade models
* require the exiting streams are at dew/bubble points.

* Shadow streams are required for the calculations. This is because the real
* stream is in thermodynamic equilibrium with the shadow stream of the opposite
* phase. This is practically implemented by copying the temperature, pressure
* and composition to the shadow stream of the same phase as the real stream.
* The two shadow streams are then constrained to be in thermodynamic equilibrium.

* Format reminder
* ComMap(Str,Str2,Str3)
* Str1 - actual process stream
* Str2 - vapor shadow
* Str3 - liquid shadow

*----------------------------------------------------*
*       Equations for Dew / bubble point calculations
*----------------------------------------------------*

***** Type 1: Actual stream is at the dew or bubble point *****

* Str: vapor, Str2: liquid, Str3: unused
$macro type1cond(junk) ((BubPoint(Str2) AND ComMap(Str2,Str,Str3) AND NOT InactiveStr(Str2)) OR (DewPoint(Str) AND ComMap(Str,Str3,Str2) AND NOT InactiveStr(Str)))

Equations
  EqT1SumMF                      Type 1 - sum of mass fractions
  EqT1cpyT                       Type 1 - copy temperature
  EqT1cpyP                       Type 1 - copy pressure
  EqT1VLE                        Type 1 - equilibrium
;

* Rachford-Rice formulation. Recall F and Fc are not used for shadow streams.
EqT1SumMF(Str,Str2,Str3)$type1cond(0)..
  SUM(Comp, Xc(Str, Comp) - Xc(Str2,Comp)) =e= 0;

EqT1cpyP(Str,Str2,Str3)$type1cond(0)..
  P(Str) =E= P(Str2);

EqT1cpyT(Str,Str2,Str3)$type1cond(0)..
  Tscaled(Str) =E= Tscaled(Str2) ;

EqT1VLE(Str,Str2,Str3,Comp)$type1cond(0)..
  Xc(Str,Comp)*phiEOS(Str,Comp) =E= Xc(Str2,Comp)*phiEOS(Str2,Comp) ;

***** Type 2: Calculate the dew or bubble point for the actual stream *****

* Str: Actual process stream, Str2: vapor shadow, Str3: liquid shadow
$macro type2cond(junk) (PhaseStability(Str) AND ComMap(Str,Str2,Str3) AND NOT InactiveStr(Str))

Equations
  EqT2CopyX                      Type 2 - copy composition
  EqT2SumMF                      Type 2 - sum mass fraction
  EqT2EqlT                       Type 2 - temperature equilibrium
  EqT2EqlP                       Type 2 - pressure equilibrium
  EqT2CopyP                      Type 2 - copy pressure
  EqT2VLE                        Type 2 - vapor-liquid equilibrium
;

* Str: Actual process stream, Str2: shadow stream with same phase as Str, Str3: unused
EqT2CopyX(Str,Str2,Str3,Comp)$(PhaseStability(Str) AND NOT InactiveStr(Str) AND
         ( (ComMap(Str,Str2,Str3) AND (FlashVap(Str) OR VapStr(Str))) OR
           (ComMap(Str,Str3,Str2) AND (FlashLiq(Str) OR LiqStr(Str))) ) )..
  Xc(Str,Comp) =e= Xc(Str2,Comp);

EqT2SumMF(Str,Str2,Str3)$type2cond(0)..
  SUM(Comp, Xc(Str3, Comp) - Xc(Str2,Comp)) =e= 0;

EqT2EqlT(Str,Str2,Str3)$type2cond(0)..
  Tscaled(Str2) =e= Tscaled(Str3);

EqT2EqlP(Str,Str2,Str3)$type2cond(0)..
  P(Str2) =e= P(Str3);

EqT2CopyP(Str,Str2,Str3)$type2cond(0)..
  P(Str) =e= P(Str2);

EqT2VLE(Str,Str2,Str3,Comp)$type2cond(0)..
  Xc(Str3,Comp)*phiEOS(Str3,Comp) =E= Xc(Str2,Comp)*phiEOS(Str2,Comp) ;


Model BubbleDewEqns / EqT1SumMF, EqT1cpyT, EqT1cpyP, EqT1VLE /;

Model PhaseStabilityEqns /EqT2CopyX, EqT2SumMF, EqT2EqlT, EqT2EqlP, EqT2CopyP,
         EqT2VLE, EqPhaseStabLiq, EqPhaseStabFlashLiq, EqPhaseStabVap, EqPhaseStabFlashVap/;
