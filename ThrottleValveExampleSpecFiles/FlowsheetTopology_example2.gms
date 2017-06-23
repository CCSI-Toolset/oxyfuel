*********************************
***** FlowsheetTopology.gms *****
*********************************

* This file contains all of the flowsheet topology information. In principle
* only this file needs to be changed to consider a new flowsheet. (One may
* also find tweaking the initialization procedures improves performance of
* the NLP solvers.)

* Notes: This version (FlowsheetTopology_23.gms) includes extra heat exchangers.

**********************************
***** Streams and Components *****
**********************************

Sets
* Note AllComp must include all of the components listed in ThermoData
AllComp          All available components        / N2 Nitrogen
                                                   O2 Oxygen
                                                   Ar Argon/

Comp(AllComp)    Components in flowsheet         / N2, O2, Ar /

* Streams
Str              All available streams           /S1*S7
                                                  V1*V500,L1*L500,
                                                  Vb1*Vb350, Lb1*Lb350,
                                                  Ve1*Ve350, Le1*Le350,
                                                  Vin1*Vin350, Lin1*Lin350,
                                                  SdV1*SdV200, SdL1*SdL200/
* Available Shadow Streams
VapShd(Str)      Vapor shadow streams            /SdV1*SdV200/
LiqShd(Str)      Liquid shadow streams           /SdL1*SdL200/ ;

Alias         (Comp , Comp2), (Comp, Comp3), (Str , Str2), (Str, Str3), (Str, Str4)  ;

*j and j2 are used an indicies for components in the thermo model files
Alias (Comp, j);
Alias (Comp2, j2);


******************************
***** Flowsheet Topology *****
******************************

Sets

* Actual process streams. Shadow streams are not included.
RealStr(Str)     actual process streams          / /,

* Streams pruned from the flowsheet are added to this set. The actually pruning
* is done by removing all members of this set from RealStr.
InactiveStr(Str) inactive streams... from flowsheet pruning / /,

* Active Shadow Streams.
ActShdStr(Str)  Active shadow streams           / /,

***** Equipment *****
* Thermodynamic Equipment. This is a generalization of all equipment with at
* least one feed stream but only two possible exit streams: one liquid and one
* vapor. The exit streams are in thermodynamic equilibrium, hence the name.
GnrlE         General equipment              /R1, C1, F1, Vlv1, HX1*HX2, SHX1*SHX250/
ThrmE(GnrlE)  Equip with thermo calcs        /R1, F1, Vlv1, HX1*HX2, SHX1*SHX250/

Ed1           Cascade                        /1/
PReb(ThrmE)   Partial Reboiler               /R1/
Tcond(GnrlE)  Total Condenser                /C1/
Flsh(ThrmE)   Flash                          /F1/
Sptr          Splitter                       /1/
Mxr           Mixer                          /M1/
Valve(ThrmE)  isenthalpic valve              /Vlv1/

HtEx(ThrmE)    Heat exchangers                /HX1*HX2/
CoolHtEx(HtEx) HX providing cooling          /HX1/
HeatHtEx(HtEx) HX providing heating          /HX2/
HtExGrps       HX groupings                  /1/

SHtEx(ThrmE)   Sub heat exchangers           /SHX1*SHX250/

*-------------------------------------
* Flowsheet topology
*-------------------------------------

***** Cascades *****
* Each cascade phase inlet vapor, outlet vapor, inlet liquid and outlet streams.
* Connectivity is managed by the following sets.
InVEd1(Ed1,Str)       Inlet vapor stream to cascade          / /
InLEd1(Ed1,Str)       Inlet liquid stream to cascade         / /
OutVEd1(Ed1,Str)       Outlet vapor stream to cascade        / /
OutLEd1(Ed1,Str)       Outlet liquid stream to cascade       / /

* This set contains all of the cascade inlet or outlet streams. In principle
* population of the set should be automated.
CscStr(Str)            Inlet and outlet streams for cascades     / /

***** Bubble and Dew Point Calculations *****
* For every bubble/dew point calculation with the CEOS thermo there are
* two extra shadow streams. These extra streams are required to ensure the
* real process stream is in the two phase region. The fake (shadow) stream
* is used for equilibrium calculations.
BPStr(Str,Str2)  "streams for which bubble point calculation are done"      / /
DPStr(Str,Str2)  "streams for which dew point calculation are done"         / /

***** Partial Reboilers *****
* Partial reboilers are thermodynamic equipment, thus they have two
* outlet streams: one liquid and one vapor.
InPReb(PReb,Str)       Inlet liquid stream to partial reboiler   / /
OutVPReb(PReb,Str)     Outlet Vapor stream to partial reboiler   / /
OutLPReb(PReb,Str)     Outlet Liquid stream to partial reboiler  / /

***** Total Condenser *****
InTCond(TCond,Str)       Inlet vapor stream to total condensor     / /
OutTCond(TCond,Str)      Outlet liquid stream at bubble point      / /

***** Flash Vessel *****
InFlsh(Flsh,Str)    Inlet stream to flash           / /
OutVFlsh(Flsh,Str)  Outlet vapor stream of flash    / /
OutLFlsh(Flsh,Str)  Outlet liquid stream of flash   / /

ExtraClmnFeeds(Str) Streams that feed into first or last stage of a column / /

***** Heat Exchangers *****
InHtEx(HtEx,Str)       Inlet Stream to valve           / HX1.S1, HX2.(S4,S5)/
OutVHtEx(HtEx,Str)     Vapor Outlet Stream of valve    / HX1.S2, HX2.S6 /
OutLHtEx(HtEx,Str)     Liquid Outlet Stream of valve   / HX1.S3, HX2.S7 /

* To better approximate nonlinearities in heat capacity, each large heat
* exchanger is broken into multiple units.
* OBSOLETE - needs to be removed
* HtExBig(HtExGrps,HtEx)                                   /1.HX1/


* These streams maybe constrained to be in the two phase region.
* OBSOLETE - needs to be removed
* TwoPhaseLiq(Str)                                        / /
* TwoPhaseVap(Str)                                        / /

***** Splitter *****
InSptr(Sptr,Str)       Inlet stream to splitter / /
OutSptr(Sptr,Str)  Outlet streams of splitter / /

OutSptrAll(Str)    All outlet streams for splitter / /

***** InMxr *****
InMxr(Mxr, Str)          Inlet streams for mixers                / /
OutMxr(Mxr,Str)          Outlet streams for mixers               / /
InMxrPres(Str)           Inlets to consider for pressure "balance" / /

***** Joule-Thomson Valves *****
InValve(Valve,Str)       Inlet Stream to valve           /Vlv1.(S2,S3) /
OutVValve(Valve,Str)     Vapor Outlet Stream of valve    /Vlv1.S4 /
OutLValve(Valve,Str)     Liquid Outlet Stream of valve   /Vlv1.S5 /

***** Thermo Equipment *****
ThrmEOutDewPoint(ThrmE)  ThrmE equipment whose outlet should be at its dew point /  /

ThrmEOutBubPoint(ThrmE)  ThrmE equipment whose outlet should be at its bubble point /  /;

loop(Sptr,
  OutSptrAll(Str)$OutSptr(Sptr,Str) = yes;
);

loop(Ed1,
  CscStr(Str)$(InVEd1(Ed1,Str) OR InLEd1(Ed1,Str) OR OutVEd1(Ed1,Str) OR OutLEd1(Ed1,Str)) = yes;
);

Variable
  PresDropCnst(Ed1);

****************************
***** Heat Integration *****
****************************

Sets
HIZone                          /Z1 /
HIMap(GnrlE,HIZone)      / /
AlreadyMapped(GnrlE)             / /;

* By default map all other GnrlE into Zone 1

loop(HIZone,
  AlreadyMapped(GnrlE)$HIMap(GnrlE, HIZone) = yes;
);

HIMap(PReb,'Z1')$(NOT AlreadyMapped(PReb)) = yes;
HIMap(TCond,'Z1')$(NOT AlreadyMapped(TCond)) = yes;
HIMap(HtEx,'Z1')$(NOT AlreadyMapped(HtEx)) = yes;

display HIMap;

Parameter
  HRAT(HIZone)   Mininimium dTemp for heat integration (K)

*********************************************************
***** Prevent pressure constraint overspecification *****
*********************************************************

* This set is used to keep track of the stream pairs for every ThrmE unit.
StrPairs(Str,Str2)       Outlet stream pairs     / / ;

****************************************
***** Distillation Column Groupings*****
****************************************

* Two or more cascades can be linked together into a distillation column. The
* following sets allow for this, which enables the calculation of recoveries,
* reflux ratios, reboil ratios, bottom-to-feed ratios and distillate-to-feed
* ratios. Not only can these quantities be added to the optimization problem,
* they also enable easier validation of designs in Aspen.

Set
  Clmn                   Specify and name columns                / DummyClmn /
  ClmnFeed(Clmn,Str)     Feed streams for each column            / /
  ClmnBtm(Clmn,Str)      Bottom streams for each column          / /
  ClmnDst(Clmn,Str)      Distillate streams for each column      / /
  ClmnRflx(Clmn,Str)     Reflux streams for each column          / /
  ClmnRb(Clmn,Str)       Reboil streams for each column          / /
  ClmnLghtK(Clmn,Comp)   Light key for each column               / /
  ClmnHvyK(Clmn,Comp)    Heavy key for each column               / /;

* Each column may also have pressure bounds specified for it. These parameters
* are used in Bounds.gms.
Parameters
  PmaxC(Clmn)    Max pressure for column    / /
  PminC(Clmn)    Min pressure for column    / /;

*********************************************
***** Distillation Column Pressure Drop *****
*********************************************

Variable
  PresDropCnst(Ed1);

PresDropCnst.fx('1') = 0;

***********************************************
***** Thermodynamic Sets and Calculations *****
***********************************************

* Primarily used for the CEOS, these sets manage which thermodynamic
* equations should be evaulated for each stream.

Sets

fEOSStr(Str)          "Stream for which EOS should be evaluated"                 / /
CPStr(Str)            Streams to calculate heat capacity of gas                  / /
LiqStr(Str)           Streams that are liquid and phase won't disappear          / /
VapStr(Str)           Streams that are vapor and phase won't disappear           /S1/

FlashLiq(Str)         Liquid streams that may disappear                          /  /
FlashVap(Str)         Vapor streams that may disappear                           /  /;

*********************************
***** Manual Stream Removal *****
*********************************

Set
  ManualRemove(Str) / /;

***************************************************************
***** Maintain Compatibility with Post Processing Scripts *****
***************************************************************

Alias (ActShdStr, ShdwStr), (RealStr, NotShdwStr);


