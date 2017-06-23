**************************
***** ValveModel.gms *****
**************************

*---------------------------------*
* Variables for Valve
*---------------------------------*
*Declared with the equations

Positive Variables
TinV(Valve)              Temperature in for each valve
ToutV(Valve)             Temperature out for each valve
PinV(Valve)              Pressure in for each valve
PoutV(Valve)             Pressure out for each valve;

*----------------------------------------*
*   Equations for valve                  *
*----------------------------------------*
Equations
EqPRel1Valve(Valve,Str,Str2)     Feed pressure greater than or equal to exit pressure
EqPinValve(Valve,Str)            Calculate Pin
EqPoutValve(Valve,Str)           Calculate Pout
EqTinValve(Valve,Str)            Calculate Tin
EqToutValve(Valve,Str)           Calculate Tout;

* Fix pressure
Qout.fx(Valve) = 0;
Qin.fx(Valve) = 0;

EqPRel1Valve(Valve,Str,Str2)$(InThrmEPres(Valve,Str) and OutLValve(Valve,Str2) AND NOT InactiveGnrlE(Valve))..
  P(Str) =g= P(Str2);

* Require for now multiple stream pairs entering a valve must be at the same
* pressure and temperature.

EqPinValve(Valve,Str)$(InThrmEPres(Valve,Str) AND NOT InactiveGnrlE(Valve))..
  PinV(Valve) =e= P(Str);

EqPoutValve(Valve,Str)$(OutLValve(Valve,Str) AND NOT InactiveGnrlE(Valve))..
  PoutV(Valve) =e= P(Str);

EqTinValve(Valve,Str)$(InThrmEPres(Valve,Str) AND NOT InactiveGnrlE(Valve))..
  TinV(Valve) =e= Tscaled(Str)*Tref;

EqToutValve(Valve,Str)$(OutLValve(Valve,Str) AND NOT InactiveGnrlE(Valve))..
  ToutV(Valve) =e= Tscaled(Str)*Tref;

Model EqValve /EqPRel1Valve/;

* These extra equations are used for the Simple thermo model only
Model EqValveExtra /EqPinValve,EqPoutValve,EqTinValve,EqToutValve/;
