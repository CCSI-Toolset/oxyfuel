*****************************
***** ReboilerModel.gms *****
*****************************

* This file contains the partial reboiler model. Reboilers are ThrmE unit and
* thus inherit equations from ThrmEModel.gms.

*---------------------------------*
* Variables for Partial Reboiler
*---------------------------------*
Positive Variables
   PPReb(PReb)        operating pressure of partial reboiler
;

*----------------------------------------------------*
*        Equations for Partial Reboiler
*----------------------------------------------------*
Equations
EqPRel1PReb(PReb)  Pv = PPReb
EqPRel3PReb(PReb,Str)  Pin = PPReb
EqTRelPReb(PReb,Str,Str2) ;

EqPRel1PReb(PReb)$(NOT InactiveGnrlE(PReb))..
    SUM(Str $OutVPReb(PReb,Str), P(Str)) =E= PPReb(PReb) ;

EqPRel3PReb(PReb,Str)$(NOT InactiveGnrlE(PReb) AND InThrmEPres(PReb,Str) )..
  P(Str) =e= PPReb(PReb);

EqTRelPReb(PReb,Str,Str2)$(InPReb(PReb,Str) and OutVPReb(PReb,Str2) AND NOT InactiveGnrlE(PReb))..
  Tscaled(Str) =l= Tscaled(Str2);


Model EqPReb /EqPRel1PReb, EqTRelPReb, EqPRel3PReb/;

* Removing temporarily to avoid pressure over/under specification in cascade
* EqPRel3PReb
