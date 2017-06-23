*******************************
***** ProbObjSpecInit.gms *****
*******************************
* Problem specific objective function initialization

***** Complementarity Slacks *****
  liqPen.l = sum(Str$FlashLiq(Str), sL.l(Str)*F.l(Str));
  vapPen.l = sum(Str$FlashVap(Str), sV.l(Str)*F.l(Str));

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

$offlisting

***** Section 4 *****
elseif SectionSwitch eq 4,
$onlisting


$offlisting

***** Section 5 *****
elseif SectionSwitch eq 5,
$onlisting

$offlisting

***** Section 6 *****
elseif SectionSwitch eq 6,
$onlisting

$offlisting

);
