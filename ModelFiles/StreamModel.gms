****************************
***** StreamModels.gms *****
****************************

*----------------------------**
** Variables for STREAMS
*----------------------------**

Positive Variables
   F(Str)        Flowrate of Stream (mol per time)
   Fc(Str,Comp)  Component flowrate of stream (mol per time)
   Xc(Str,Comp)  Component split fraction
   T(Str)        Temperature of stream (K)
   P(Str)        Pressure of stream (bar)
*   SptrFrac(Str)  Split fraction for stream (streams that are output of splitter)
;

Variables
H(Str)           Enthalpy (kJ per mole)
S(Str)           Entropy (J per mole-K)
V(Str)           Specific Volume (m^3 per kmol)
HIG(Str)         IG Enthalpy (J per mole)
SIG(Str)         IG Entropy (kJ per mole-K)
MW(Str)          Molecular weight
phiEOS(Str,Comp) Fugacity
Tscaled(Str)     Scaled temperature
Pd(Str)          Stream dew point (bar)
Pb(Str)          Stream bubble point (bar)
Pvap(Str,Comp)   Vapor pressure (bar);

Positive Variables
  sL(Str)          Slack variable for disappearing liquid streams
  sV(Str)          Slack variable for disappearing vapor streams;

* Note: The thermodynamic properties - H, S, V, Pvap, etc - are
* calculated in the thermodynamic model files.

*----------------------------------------------------*
*        Equations for streams
*----------------------------------------------------*
Equations
EqTotMolStr(Str)         Total mass balance for stream (Total Mole flow = sum of component Mole flow)
EqSplitFrac(Str,Comp)    Component split fraction
EqTscaled(Str)           Dimensionless temperature
EqSumFrac(Str) ;

EqTotMolStr(Str)$(RealStr(Str))..
    F(Str) =E=  Sum(Comp, Fc(Str,Comp)) ;

EqSplitFrac(Str,Comp)$(RealStr(Str) AND ord(Comp) < card(Comp))..
*EqSplitFrac(Str,Comp)$(RealStr(Str))..
    Xc(Str,Comp)*F(Str) =e= Fc(Str,Comp);

EqTscaled(Str)$(HCalc(Str) OR SCalc(Str))..
  Tscaled(Str)*Tref =e= T(Str);

* EqSumFrac is used during CEOS initialization. It is important this equation
* is written for all streams, not just RealStr.


* EqSumFrac(Str)$(RealStr(Str) AND NOT OutThrmEStr(Str))..
*EqSumFrac(Str)$( RealStr(Str) AND StrModCalc(Str) )..
EqSumFrac(Str)$RealStr(Str)..
  SUM(Comp, Xc(Str,Comp)) =e= 1;

Model EqStreams /EqTotMolStr, EqSplitFrac, EqTscaled, EqSumFrac/;
* EqSumFrac
* EqTotMolStr
