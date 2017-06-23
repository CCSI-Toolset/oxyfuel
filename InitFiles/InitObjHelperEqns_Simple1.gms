*********************************
***** InitObjHelperEqns.gms *****
*********************************

* This file initializes variables used in the object function "helper" equations

***** Complementarity Slacks *****
  liqPen.l = sum(Str$FlashLiq(Str), sL.l(Str)*F.l(Str));
  vapPen.l = sum(Str$FlashVap(Str), sV.l(Str)*F.l(Str));

***** Oxygen Recovery and Purity Calcs *****
  o2pure.l = Xc.l('S10','O2');
  o2rec.l = Fc.l('S10','O2')/max(SUM(FeedStr, Fc.l(FeedStr,'O2') ),epsi);
