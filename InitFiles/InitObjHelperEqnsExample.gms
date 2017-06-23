*********************************
***** InitObjHelperEqns.gms *****
*********************************

* This file initializes variables used in the object function "helper" equations

***** Complementarity Slacks *****
  liqPen.l = sum(Str$FlashLiq(Str), sL.l(Str)*F.l(Str));
  vapPen.l = sum(Str$FlashVap(Str), sV.l(Str)*F.l(Str));

***** Compressor Work *****
  CmprPwr.l(FeedStr) = Ncmpr*gamma/(gamma - 1)*Rsi*300*( P.l(FeedStr)**((gamma - 1)/(gamma*Ncmpr)) - 1);
