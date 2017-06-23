***********************
***** ObjFncs.gms *****
***********************

* This file contains equations required for the objective function. There are
* several objective functions available in this file.

***** Complementarity Penalties *****

Positive Variables
  liqPen
  vapPen;

Equations
  EqLiqPen               Calculate liquid penalty - used for complementarities
  EqVapPen               Calculate vapor penalty - used for complementarities;

EqVapPen..
  vapPen =e= sum(Str$(FlashVap(Str) AND ( (NOT InactiveStr(Str) AND RealStr(Str)) OR ActShdStr(Str)) ), sV(Str)*F(Str));


EqLiqPen..
  liqPen =e= sum(Str$(FlashLiq(Str) AND ( (NOT InactiveStr(Str) AND RealStr(Str)) OR ActShdStr(Str)) ), sL(Str)*F(Str));


Model EqGenericObj /EqLiqPen, EqVapPen/;
