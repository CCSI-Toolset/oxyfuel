*********************************
***** InitObjHelperEqns.gms *****
*********************************

* This file initializes variables used in the object function "helper" equations

***** Complementarity Slacks *****
  liqPen.l = sum(Str$FlashLiq(Str), sL.l(Str)*F.l(Str));
  vapPen.l = sum(Str$FlashVap(Str), sV.l(Str)*F.l(Str));

***** Oxygen Recovery and Purity Calcs *****
  o2pure.l = (Fc.l('S15','O2')+Fc.l('S17','O2'))/max(sum(Comp, Fc.l('S15',Comp) + Fc.l('S17',Comp)),epsi);
  o2rec.l = (Fc.l('S15','O2')+Fc.l('S17','O2'))/max(SUM(FeedStr, Fc.l(FeedStr,'O2') ),epsi);

***** Compressor Work *****
  CmprPwr.l(FeedStr) = F.l(FeedStr)*Ncmpr*gamma/(gamma - 1)*Rsi*300*( P.l(FeedStr)**((gamma - 1)/(gamma*Ncmpr)) - 1)
                         /max(Fc.l('S15','O2') + Fc.l('S17','O2'),epsi)/(1000*115.2);
