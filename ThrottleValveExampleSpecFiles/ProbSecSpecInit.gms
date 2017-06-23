*******************************
***** ProbSecSpecInit.gms *****
*******************************


$offlisting

***** Section 0 *****

if(SectionSwitch eq 0,
$onlisting

HRAT(HIZone) = HRATSimple;

* Turn off GnrlE model for valves - only for ASU with some simple thermo
* enthalpy models. Otherwise ensure GnrlE model (energy balance) is on
CalcGnrlE(Valve) = yes;

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
