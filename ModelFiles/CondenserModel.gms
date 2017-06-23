******************************
***** CondenserModel.gms *****
******************************

* This file contains equations for the total condenser model. Recall the
* condensers are not considered ThrmE units; they only have one exit stream (liquid).

*---------------------------------*
* Variables for total Condenser
*---------------------------------*
positive variables
  PTCond(Tcond)      Pressure of total Condensor
;

*----------------------------------------------------*
* Equations for Total Condensor
*----------------------------------------------------*

Equations
EqCompMolBalTCond(TCond,Comp)            Component mole balance
*EqBubTTCond(TCond,Str)                   Bubble point eq.
EqQTCond(Tcond)                          Energy balance
EqPRel1TCond(Str,TCond)                  Calculate PTCond & ensure const pres.
EqBubTCond(TCond,Str)                    Ensure exiting stream is a liquid - obsolete
EqTRelTCond(TCond,Str,Str2)              Temperature can only decrease in a condenser
;

EqCompMolBalTcond(TCond,Comp)$(NOT InactiveGnrlE(TCond))..
   sum(Str $InTCond(TCond,Str), Fc(Str,Comp)) =e=
   sum(Str $OutTCond(TCond,Str), Fc(Str,Comp));

EqQTCond(TCond)$(NOT InactiveGnrlE(TCond))..
    Qout(TCond) =E=  SUM(Str $InTCond(TCond,Str), F(Str)*H(Str))
                           - SUM(Str $OutTCond(TCond,Str), F(Str)*H(Str))  ;

EqPRel1TCond(Str,TCond)$(OutTCond(TCond,Str) OR InTCond(TCond,Str) AND NOT InactiveGnrlE(TCond))..
   P(Str) =e= PTCond(TCond);

EqBubTCond(TCond,Str)$OutTCond(TCond,Str)..
  Pb(Str) =l= P(Str);

EqTRelTCond(TCond,Str,Str2)$(InTCond(TCond,Str) and OutTCond(TCond,Str2) AND NOT InactiveGnrlE(TCond))..
  Tscaled(Str) =g= Tscaled(Str2);

Model EqTCond /EqPRel1TCond, EqTRelTCond/;

Qin.fx(TCond) = 0;

* Removing temporarily to avoid pressure over/under specification in cascade
* OR InTCond(TCond,Str) from EqPRel1TCond
