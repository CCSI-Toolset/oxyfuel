**************************
***** FlashModel.gms *****
**************************

* This file contains supplemental equations for the flash unit models. Recall
* flash vessels are considered ThrmE units. They inherit equations for
* ThrmEModel.gms.

*--------------------------------**
* Variables for Flash
*--------------------------------**

Positive Variables
TFlsh(Flsh)
PFlsh(Flsh)
;

*----------------------------------------------------*
*        Equations for Flash
*----------------------------------------------------*
Equations

EqPRel1Flsh(Flsh)        Calculate PFlsh based on outlet pressure
EqTRel1Flsh(Flsh)        Calculate TFlsh based on outlet temperature
EqPFlashLgc(Flsh,Str)    Inlet pressure is greater than or equal to outlet pres
EqPFlashLgc2(Flsh,Str) Flash operates at constant pressure - for simple thermo model
;

EqPRel1Flsh(Flsh)$(NOT InactiveGnrlE(Flsh))..
        Pflsh(Flsh) =E= SUM(Str$OutVFlsh(Flsh,Str), P(Str)) ;

EqTRel1Flsh(Flsh)$(NOT InactiveGnrlE(Flsh))..
        Tflsh(Flsh) =E= SUM(Str$OutVFlsh(Flsh,Str), Tscaled(Str)) ;

EqPFlashLgc(Flsh,Str)$(InThrmEPres(Flsh,Str) AND NOT InactiveGnrlE(Flsh))..
  P(Str) =g= PFlsh(Flsh);
*  P(Str) =e= PFlsh(Flsh);

EqPFlashLgc2(Flsh,Str)$(InThrmEPres(Flsh,Str) AND NOT InactiveGnrlE(Flsh))..
  P(Str) =l= PFlsh(Flsh) + 5;

***** Flashes are adiabatic *****
Qin.fx(Flsh) = 0;
Qout.fx(Flsh) = 0;

Model EqFlsh /EqPRel1Flsh, EqPFlashLgc/
* EqPFlashLgc2
* EqTRel1Flsh
