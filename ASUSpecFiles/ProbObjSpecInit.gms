*******************************
***** ProbObjSpecInit.gms *****
*******************************
* Problem specific objective function initialization

***** Complementarity Slacks *****
  liqPen.l = sum(Str$FlashLiq(Str), sL.l(Str)*F.l(Str));
  vapPen.l = sum(Str$FlashVap(Str), sV.l(Str)*F.l(Str));

***** Oxygen Recovery and Purity Calcs *****
  o2pure.l = (Fc.l('S15','O2')+Fc.l('S17','O2'))/max(sum(Comp, Fc.l('S15',Comp) + Fc.l('S17',Comp)),epsi);
  o2rec.l = (Fc.l('S15','O2')+Fc.l('S17','O2'))/max(SUM(FeedStr, Fc.l(FeedStr,'O2') ),epsi);

***** Compressor Work *****
  CmprPwr.l(FeedStr) = F.l(FeedStr)*Ncmpr*gamma/(gamma - 1)*Rsi*300*( P.l(FeedStr)**((gamma - 1)/(gamma*Ncmpr)) - 1)
                         /max(Fc.l('S15','O2') + Fc.l('S17','O2'),epsi)/(1000*115.2);

$offlisting

***** Section 0 *****

if(SectionSwitch eq 0,
$onlisting

$offlisting

***** Section 1 *****
elseif SectionSwitch eq 1,
$onlisting

$offlisting

***** Section 2 *****
elseif SectionSwitch eq 2,
$onlisting

$offlisting

***** Section 3 *****
elseif SectionSwitch eq 3,
$onlisting

Z9.l = ComplPen*vapPen.l + ComplPen*liqPen.l + SUM(FeedStr, CmprPwr.l(FeedStr))
                 + QwWeight*Qw.l + Qs.l;

$offlisting

***** Section 4 *****
elseif SectionSwitch eq 4,
$onlisting

Z9.l = ComplPen*vapPen.l + ComplPen*liqPen.l + SUM(FeedStr, CmprPwr.l(FeedStr))
                 + QwWeight*Qw.l + Qs.l;

$offlisting

***** Section 5 *****
elseif SectionSwitch eq 5,
$onlisting

$offlisting

***** Section 6 *****
elseif SectionSwitch eq 6,
$onlisting

Z9.l = ComplPen*vapPen.l + ComplPen*liqPen.l + SUM(FeedStr, CmprPwr.l(FeedStr))
                 + QwWeight*Qw.l + Qs.l;

$offlisting

);
