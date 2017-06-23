*******************************
***** ProbSecSpecInit.gms *****
*******************************

HRAT('Z2') = 0.4;

$offlisting

***** Section 0 *****

if(SectionSwitch eq 0,
$onlisting

HRAT(HIZone) = HRATSimple;

* Turn off GnrlE model for valves - only for ASU with some simple thermo
* enthalpy models. Otherwise ensure GnrlE model (energy balance) is on
CalcGnrlE(Valve) = yes;

o2pure.l = 0.5;
o2rec.l = 0.5;

$offlisting

***** Section 1 *****
elseif SectionSwitch eq 1,
$onlisting

o2rec.lo = O2RecLowCEOS + O2RecLowDelta;
o2pure.lo = o2pureSpec;
o2rec.up = 0.97;

P.fx(FeedStr) = 3;

T.up(ColumnFeedStr) = 200;

Qin.up(PReb) = 30;

$offlisting

***** Section 2 *****
elseif SectionSwitch eq 2,
$onlisting

$offlisting

***** Section 3 *****
elseif SectionSwitch eq 3,
$onlisting

epsi = 1e-6;

ComplPen = 100;

QwWeight = 5;

o2pure.lo = o2pureSpec;
o2rec.lo = O2RecLowCEOS;


Xc.lo(Str,Comp) = 0;

T.lo(CscStr) = Tmin;

T.up(LiqStr) = smax(Comp, Tc(Comp)) + 5;

$offlisting

***** Section 4 *****
elseif SectionSwitch eq 4,
$onlisting

EpsilonZ = 1E-6;
Inter1.lo(Str) = 1E-6;

QwWeight = 5;

$offlisting

***** Section 5 *****
elseif SectionSwitch eq 5,
$onlisting

$offlisting

***** Section 6 *****
elseif SectionSwitch eq 6,
$onlisting

o2pure.lo = o2pureSpec;

if(DualReboilerASU eq 1,
  HRAT(HIZone) = HRATCEOS;
else
  HRAT('Z1') = HRATCEOS;
);

$offlisting

);
