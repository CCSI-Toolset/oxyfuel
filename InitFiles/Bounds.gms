*------------------------**
*      Bounds
*------------------------**

* Flows
F.lo(Str) = 0;
F.up(Str) = 100 ;
Fc.lo(Str,Comp) = 0;
Fc.up(Str,Comp) = 100 ;
F.lo(FeedStr) = 0.1;
*F.lo(FeedStr) = 0;
F.up(FeedStr) = 2;
F.lo(CscStr) = 0.001;

* Mole fraction
*Xc.lo(Str,Comp) = 1E-8;
Xc.lo(Str,Comp) = 0;
Xc.up(Str,Comp) = 1;

* Temperature
T.lo(Str) = Tmin ;
T.up(Str) = Tmax ;
*T.up(LiqStr) = 126.2;
*T.up(LiqStr) = smin(Comp, Tc(Comp)) - 0.1;
*T.lo('S26') = 60;
*T.lo('S27') = 60;

*T.lo(FeedStr) = 300 ;
*T.up(FeedStr) = 350;

Tscaled.lo(Str) = T.lo(Str)/Tref;
Tscaled.up(Str) = T.up(Str)/Tref;


P.lo(Str) = Pmin ;
P.up(Str) = Pmax ;

* variables of cascade
StgEd1.lo(Ed1) = LBStgEd1 ;
StgEd1.up(Ed1) = UBStgEd1 ;

* upper bound on Pressure (PEd1) ;
*PEd1.up(Ed1) = Pmax ;

* bounds on absorption and stripping factors
ANEd1.up(Ed1,Comp) = UBASFact ;
A1Ed1.up(Ed1,Comp) = UBASFact ;

ANEd1.lo(Ed1,Comp) = LBASFact ;
A1Ed1.lo(Ed1,Comp) = LBASFact ;

AeEd1.lo(Ed1,Comp) = LBASFact ;
SeEd1.lo(Ed1,Comp) = LBASFact ;

AeEd1.up(Ed1,Comp) = UBASFact ;
SeEd1.up(Ed1,Comp) = UBASFact ;

* bounds on recovery factor phi
PhiAEd1.up(Ed1,Comp) = 1 ;
PhiSEd1.up(Ed1,Comp) = 1 ;

PhiAEd1.lo(Ed1,Comp) = 1E-6;
PhiSEd1.lo(Ed1,Comp) = 1E-6;

DummyAe.lo(Ed1,Comp) = 1E-10;
DummySe.lo(Ed1,Comp) = 1E-10;


* bounds on vapor and liquid flowrate
L1Ed1.up(Ed1) = 50 ;
VNEd1.up(Ed1) = 50 ;

loop(Clmn,
*  P.up(Str)$ClmnFeed(Clmn,Str) = min(PmaxC(Clmn), Pmax);
*  P.lo(Str)$ClmnFeed(Clmn,Str) = max(PminC(Clmn), Pmin);
  P.up(Str)$ClmnRflx(Clmn,Str) = min(PmaxC(Clmn), Pmax);
  P.lo(Str)$ClmnRflx(Clmn,Str) = max(PminC(Clmn), Pmin);
);

* bounds on reboil ratio of reboiler
RebRatio.lo(Clmn) = 0.05 ;
RebRatio.up(Clmn) = 30 ;

* bounds on reflux ratio
RefRatio.lo(Clmn) = 0.1;
RefRatio.up(Clmn) = 30;

* Valves
Qin.fx(Valve) = 0;
Qout.fx(Valve) = 0;

* upper bound on heat duty of reboiler
Qin.up(PReb) = 100 ;
Qout.up(TCond) = 150;

* upper bound on PPReb ;
PPReb.up(Preb) = Pmax ;

* Temperature of flash
Tflsh.lo(Flsh) = Tmin ;
Tflsh.up(Flsh) = Tmax ;

* P of flash
Pflsh.lo(Flsh) = Pmin ;
Pflsh.up(Flsh) = Pmax ;

beta.lo(ThrmE) = -1E3;
beta.up(ThrmE) = 1E3;

***** Thermo - CEOS *****

phiEOS.lo(fEOSStr,j) = epsi;
phiEOS.up(fEOSStr,j) = 25;

* Set bounds... potential for divide by zero
bmEOS.lo(fEOSStr) = epsi;
bbEOS.lo(fEOSStr) = epsi;
amEOS.lo(fEOSStr) = epsi;

* Potential sqrt of a negative number
aEOS.lo(fEOSStr,Comp) = 0;

* Avoid sqrt of a negative number
aEOS.lo(fEOSStr,Comp) = epsi;

ZEOS.lo(Str) = epsi;
ZEOS.up(Str) = 2;
*Tscaled.lo(Str) = epsi;

sL.fx(LiqStr) = 0;
sV.fx(VapStr) = 0;

sL.fx(CscStr) = 0;
sV.fx(CscStr) = 0;

sL.up(Str)$(NOT sL.up(Str) eq 0) = 1000;
sV.up(Str)$(NOT sV.up(Str) eq 0) = 1000;

F.fx(Str)$ManualRemove(Str) = 0;
Fc.fx(Str,Comp)$ManualRemove(Str) = 0;
