*******************************
***** CEOSDewBubModel.gms *****
*******************************

* This file contains the equations required for dew and bubble point
* calculations using the cubic EOS framework. The Edmister cascade models
* require the exiting streams are at dew/bubble points.

* Streams involved in BP/DP calculations are specified in the following sets:
*   BPStr(Str,Str2,Str3)
*   DPStr(Str,Str2,Str3)
* "Str" is the real process stream at either its bubble or dew point.
* "Str2" is the vapor shadow stream
* "Str3" is the liquid shadow stream

* Shadow streams are required for the calculations. This is because the real
* stream is in thermodynamic equilibrium with the shadow stream of the opposite
* phase. This is practically implemented by copying the temperature, pressure
* and composition to the shadow stream of the same phase as the real stream.
* The two shadow streams are then constrained to be in thermodynamic equilibrium.

*----------------------------------------------------*
*       Equations for Dew / bubble point calculations
*----------------------------------------------------*

****** Common Equipments for both Dew & Bubble Point Calcs *****
Equations
EqSumMFBDP(Str,Str2)  VLE equation for DEW BUB point calculation
EqtrP(Str,Str2)   "Transfer P of liq = P of Vap"
EqtrT(Str,Str2)   "Transfer T of liq = T of Vap"
;

EqSumMFBDP(Str,Str2)$((BPStr(Str,Str2)+ DPStr(Str,Str2))*RealStr(Str))..
  SUM(j,Xc(Str,j)) - SUM(j,Xc(Str2,j)) =E= 0  ;

EqtrP(Str,Str2)$((BPStr(Str,Str2)+ DPStr(Str,Str2))*RealStr(Str))..
  P(Str) =E= P(Str2) ;

EqtrT(Str,Str2)$((BPStr(Str,Str2)+ DPStr(Str,Str2))*RealStr(Str))..
  Tscaled(Str) =E= Tscaled(Str2) ;


***** Equations only for bubble point calc *****
Equations
EqVLEBP(Str,Str2,j)
;

EqVLEBP(Str,Str2,j)$(BPStr(Str,Str2)*RealStr(Str))..
  Xc(Str2,j)*phiEOS(Str2,j) =E= Xc(Str,j)*phiEOS(Str,j)  ;


***** Equations only for dew point calcs *****
Equations
EqVLEDP(Str,Str2,j)
;

EqVLEDP(Str,Str2,j)$(DPStr(Str,Str2)*RealStr(Str))..
  Xc(Str,j)*phiEOS(Str,j) =E= Xc(Str2,j)*phiEOS(Str2,j)  ;


Model BubbleDewEqns /EqSumMFBDP, EqtrP, EqtrT, EqVLEBP,
                         EqVLEDP/;
