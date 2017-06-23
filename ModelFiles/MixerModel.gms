*****************************
***** MixerModel.gms *****
*****************************

* This file includes the mixer model. This model structure assumes there
* is only one outlet stream for each mixer

*----------------------------------------*
*  Equation for mixer
*----------------------------------------*
Equations
EqCompMolBalMxr(Mxr,Comp)                Component balance for mixer
EqEnrgBalMxr(Mxr)                        Energy balance for mixer
EqPRel1Mxr(Mxr,Str)                      Pressure relation
;

EqCompMolBalMxr(Mxr,Comp)..
  SUM(Str$InMxr(Mxr,Str), Fc(Str,Comp)) =e= SUM(Str2$OutMxr(Mxr,Str2),Fc(Str2,Comp));

EqEnrgBalMxr(Mxr)..
  SUM(Str$InMxr(Mxr,Str), F(Str)*H(Str)) =e= SUM(Str2$OutMxr(Mxr,Str2),F(Str2)*H(Str2));

EqPRel1Mxr(Mxr,Str)$(InMxr(Mxr,Str) AND InMxrPres(Str))..
  SUM(Str2$OutMxr(Mxr,Str2), P(Str2)) =E= P(Str) ;


Model EqMxr /EqCompMolBalMxr, EqEnrgBalMxr, EqPRel1Mxr/;
