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
Str              All available streams           /S1*S38,S6V,S21V,
                                                  SV2,SL2,SV7,SL7,SV14,SL14,SV22,SL22,
                                                  SV5,SL5,SV9,SL9,SV20,SL20,SV25,SL25,
                                                  SF0,SF10,SFV1*SFV3,SFL1*SFL3,
                                                  SFV11*SFV13,SFL11*SFL13,
                                                  SNV0*SNV2,SNL0*SNL2,
                                                  SNV11*SNV12,SNL11*SNL12,
                                                  SOV0*SOV2,SOL0*SOL2,
                                                  SOV10*SOV12,SOL10*SOL12,
                                                  V1*V500,L1*L500,
                                                  Vb1*Vb350, Lb1*Lb350,
                                                  Ve1*Ve350, Le1*Le350,
                                                  Vin1*Vin350, Lin1*Lin350,
                                                  SSVt1*SSVt3, SSLt1*SSLt3, SSVm1*SSVm3,
                                                  SSLm1*SSLm3, SSVb1*SSVb3, SSLb1*SSLb3,
                                                  SFlL3*SFlL6, SFlV3*SFlV6,
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
GnrlE         General equipment              /R1*R2, C1*C2, F1*F6, Vlv1*Vlv8, HX1*HX15, SHX1*SHX500/
ThrmE(GnrlE)  Equip with thermo calcs        /R1*R2, F1*F6, Vlv1*Vlv8, HX1*HX15, SHX1*SHX500/

Ed1           Cascade                        /1*4/
PReb(ThrmE)   Partial Reboiler               /R1*R2/
Tcond(GnrlE)  Total Condenser                /C1*C2/
Flsh(ThrmE)   Flash                          /F1*F6/
Sptr          Splitter                       /1*10/
Mxr           Mixer                          /M1/
Valve(ThrmE)  isenthalpic valve              /Vlv1*Vlv8/

HtEx(ThrmE)    Heat exchangers                /HX1*HX15/
CoolHtEx(HtEx) HX providing cooling          /HX1*HX4,HX13*HX14/
HeatHtEx(HtEx) HX providing heating          /HX5*HX12, HX15/
HtExGrps       HX groupings                  /1*6/

SHtEx(ThrmE)   Sub heat exchangers           /SHX1*SHX500/

*-------------------------------------
* Flowsheet topology
*-------------------------------------

***** Cascades *****
* Each cascade phase inlet vapor, outlet vapor, inlet liquid and outlet streams.
* Connectivity is managed by the following sets.
InVEd1(Ed1,Str)       Inlet vapor stream to cascade          /1.S16, 2.S8,  3.SFlV5, 4.S23/
InLEd1(Ed1,Str)       Inlet liquid stream to cascade         /1.S1, 2.SFlL4, 3.S19, 4.SFlL6/
OutVEd1(Ed1,Str)       Outlet vapor stream to cascade        /1.S5, 2.S9,  3.S20, 4.S25/
OutLEd1(Ed1,Str)       Outlet liquid stream to cascade       /1.S2, 2.S7,  3.S14, 4.S22/

* This set contains all of the cascade inlet or outlet streams. In principle
* population of the set should be automated.
CscStr(Str)            Inlet and outlet streams for cascades     / /


***** Partial Reboilers *****
* Partial reboilers are thermodynamic equipment, thus they have two
* outlet streams: one liquid and one vapor.
InPReb(PReb,Str)       Inlet liquid stream to partial reboiler   /R1.SFlL3, R2.SFlL5/
OutVPReb(PReb,Str)     Outlet Vapor stream to partial reboiler   /R1.S4, R2.S18/
OutLPReb(PReb,Str)     Outlet Liquid stream to partial reboiler  /R1.S3, R2.S15/

***** Total Condenser *****
InTCond(TCond,Str)       Inlet vapor stream to total condensor     /C1.S35, C2.S31/
OutTCond(TCond,Str)      Outlet liquid stream at bubble point      /C1.S11, C2.S24/

***** Flash Vessel *****
InFlsh(Flsh,Str)    Inlet stream to flash           /F1.(S5,S7,SSVm1, SSLm1), F2.(S20,S22,SSVm2, SSVm3, SSLm2, SSLm3),
                                                     F3.(S2, S4, SSVb1, SSLb1), F4.(S9, S10, SSVt1, SSLt1),
                                                     F5.(S14, S18, SSVb2, SSLb2, SSVb3, SSLb3), F6.(S24, S25, SSVt2, SSLt2, SSVt3, SSLt3, S26, S27)/
OutVFlsh(Flsh,Str)  Outlet vapor stream of flash    /F1.S8, F2.S23, F3.SFlV3, F4.SFlV4, F5.SFlV5, F6.SFlV6/
OutLFlsh(Flsh,Str)  Outlet liquid stream of flash   /F1.S1, F2.S19, F3.SFlL3, F4.SFlL4, F5.SFlL5, F6.SFlL6/

ExtraClmnFeeds(Str) Streams that feed into first or last stage of a column /SSVb1*SSVb3, SSLb1*SSLb3, SSVt1*SSVt3, SSLt1*SSLt3/

***** Heat Exchangers *****
InHtEx(HtEx,Str)       Inlet Stream to valve           /HX1.SF0, HX2.(SFV1, SFL1), HX3.SF10, HX4.(SFV11, SFL11),
                                                         HX5.(SNV0, SNL0), HX6.(SNV1, SNL1), HX7.(S30, S34), HX8.(SNV11, SNL11),
                                                         HX9.(SOV0, SOL0), HX10.(SOV1, SOL1),
                                                         HX11.(SOV10, SOL10), HX12.(SOV11, SOL11),
                                                         HX13.(SFV2, SFL2), HX14.(SFV12, SFL12), HX15.(S29, S33) /
OutVHtEx(HtEx,Str)     Vapor Outlet Stream of valve    /HX1.SFV1, HX2.SFV2, HX3.SFV11, HX4.SFV12,
                                                         HX5.SNV1, HX6.SNV2, HX7.SNV11, HX8.SNV12,
                                                         HX9.SOV1, HX10.SOV2, HX11.SOV11, HX12.SOV12,
                                                         HX13.SFV3, HX14.SFV13, HX15.S37 /
OutLHtEx(HtEx,Str)     Liquid Outlet Stream of valve   /HX1.SFL1, HX2.SFL2, HX3.SFL11, HX4.SFL12,
                                                         HX5.SNL1, HX6.SNL2, HX7.SNL11, HX8.SNL12,
                                                         HX9.SOL1, HX10.SOL2, HX11.SOL11, HX12.SOL12,
                                                         HX13.SFL3, HX14.SFL13, HX15.S38 /

* To better approximate nonlinearities in heat capacity, each large heat
* exchanger is broken into multiple units.
HtExBig(HtExGrps,HtEx)                                   /1.(HX1,HX2,HX13), 2.(HX3,HX4,HX14), 3.(HX5,HX6), 4.(HX7,HX8), 5.(HX9,HX10), 6.(HX11,HX12)/


* These streams maybe constrained to be in the two phase region.
TwoPhaseLiq(Str)                                        /SFL1,SFL11/
TwoPhaseVap(Str)                                        /SFV1,SFV11/

***** Splitter *****
InSptr(Sptr,Str)       Inlet stream to splitter /1.S11, 2.SFlV4, 3.SFlV3, 4.SFlV6, 5.S6, 6.S6V, 7.S21, 8.S21V, 9.S37, 10.S38/
OutSptr(Sptr,Str)  Outlet streams of splitter /1.(S10,S12,S13), 2.(S35, S36), 3.(S16,S32), 4.(S28,S31),
                                               5.(SSLt1, SSLm1, SSLb1), 6.(SSVt1, SSVm1, SSVb1),
                                               7.(SSLt2, SSLm2, SSLb2), 8.(SSVt2, SSVm2, SSVb2),
                                               10.(SSLt3, SSLm3, SSLb3, S17), 9.(SSVt3, SSVm3, SSVb3)/

OutSptrAll(Str)    All outlet streams for splitter / /

***** InMxr *****
InMxr(Mxr, Str)          Inlet streams for mixers                / /
OutMxr(Mxr,Str)          Outlet streams for mixers               / /
InMxrPres(Str)           Inlets to consider for pressure "balance" / /

***** Joule-Thomson Valves *****
InValve(Valve,Str)       Inlet Stream to valve           /Vlv1.(S3,S32), Vlv2.S12, Vlv3.(SFV3, SFL3), Vlv4.(SFV13, SFL13), Vlv5.S28, Vlv6.S15, Vlv7.S17, Vlv8.(S13, S36)  /
OutVValve(Valve,Str)     Vapor Outlet Stream of valve    /Vlv1.S33, Vlv2.S34, Vlv3.S6V, Vlv4.S21V, Vlv5.SNV0, Vlv6.SOV0, Vlv7.SOV10, Vlv8.S26 /
OutLValve(Valve,Str)     Liquid Outlet Stream of valve   /Vlv1.S29, Vlv2.S30, Vlv3.S6, Vlv4.S21, Vlv5.SNL0, Vlv6.SOL0, Vlv7.SOL10, Vlv8.S27 /

***** Thermo Equipment *****
ThrmEOutDewPoint(ThrmE)  ThrmE equipment whose outlet should be at its dew point /HX1, HX3 /
* HX1

ThrmEOutBubPoint(ThrmE)  ThrmE equipment whose outlet should be at its bubble point /HX2, HX4 /;
* HX2

loop(Sptr,
  OutSptrAll(Str)$OutSptr(Sptr,Str) = yes;
);

loop(Ed1,
  CscStr(Str)$(InVEd1(Ed1,Str) OR InLEd1(Ed1,Str) OR OutVEd1(Ed1,Str) OR OutLEd1(Ed1,Str)) = yes;
);

Variable
  PresDropCnst(Ed1);

PresDropCnst.fx('1') = 0.0082;
PresDropCnst.fx('2') = 0.0082;
PresDropCnst.fx('3') = 0.0024;
PresDropCnst.fx('4') = 0.0024;

****************************
***** Heat Integration *****
****************************

Sets
HIZone                          /Z1, Z2/
HIMap(GnrlE,HIZone)      / /
AlreadyMapped(GnrlE)             / /;


* Specify C1 and R2 into Zones 2

HIMap('C1','Z2') = yes;
HIMap('R2','Z2') = yes;

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

StrPairs('S13','S36') = yes;
StrPairs('S3','S32') = yes;

* Added
StrPairs('SSLb1','SSVb1') = yes;
StrPairs('SSLm1','SSVm1') = yes;
StrPairs('SSLt1','SSVt1') = yes;
StrPairs('SSLb2','SSVb2') = yes;
StrPairs('SSLm2','SSVm2') = yes;
StrPairs('SSLt2','SSVt2') = yes;
StrPairs('SSLb3','SSVb3') = yes;
StrPairs('SSLm3','SSVm3') = yes;
StrPairs('SSLt3','SSVt3') = yes;

* Usage: Outlet, Inlet

* Prevents recycle loop issues?
StrPairs('SFlV6','S24') = yes;
StrPairs('SFlV4','S10') = yes;

* Added on May 15th, 2014
StrPairs('SFlL5','S18') = yes;
StrPairs('SFlL3','S4') = yes;

$ontext
* Added on May 16th, 2014
StrPairs('SFlL6','SSLt2') = yes;
StrPairs('SFlL6','SSLt3') = yes;

StrPairs('SFlV5','SSLb2') = yes;
StrPairs('SFlV5','SSLb3') = yes;

StrPairs('SFlL4','SSLt1') = yes;
StrPairs('SFlV3','SSLb1') = yes;
$offtext

****************************************
***** Distillation Column Groupings*****
****************************************

* Two or more cascades can be linked together into a distillation column. The
* following sets allow for this, which enables the calculation of recoveries,
* reflux ratios, reboil ratios, bottom-to-feed ratios and distillate-to-feed
* ratios. Not only can these quantities be added to the optimization problem,
* they also enable easier validation of designs in Aspen.

Set
  Clmn                   Specify and name columns                /HPC, LPC/
  ClmnFeed(Clmn,Str)     Feed streams for each column            /HPC.(S6V,S6), LPC.(S21V, S21, S13, S36, S37, S38)/
  ClmnBtm(Clmn,Str)      Bottom streams for each column          /HPC.(S3), LPC.(S15)/
  ClmnDst(Clmn,Str)      Distillate streams for each column      /HPC.(S12, S36, S13), LPC.S28/
  ClmnRflx(Clmn,Str)     Reflux streams for each column          /HPC.SFlL4, LPC.SFlL6/
  ClmnRb(Clmn,Str)       Reboil streams for each column          /HPC.S16, LPC.SFlV5/
  ClmnLghtK(Clmn,Comp)   Light key for each column               /HPC.N2, LPC.N2/
  ClmnHvyK(Clmn,Comp)    Heavy key for each column               /HPC.O2, LPC.O2/;

* Each column may also have pressure bounds specified for it. These parameters
* are used in Bounds.gms.
Parameters
  PmaxC(Clmn)    Max pressure for column    /HPC = 12, LPC = 8/
  PminC(Clmn)    Min pressure for column    /HPC = 1, LPC = 1/;

***********************************************
***** Thermodynamic Sets and Calculations *****
***********************************************

* Primarily used for the CEOS, these sets manage which thermodynamic
* equations should be evaulated for each stream.

Sets

fEOSStr(Str)          "Stream for which EOS should be evaluated"                 / /
CPStr(Str)            Streams to calculate heat capacity of gas                  / /
LiqStr(Str)           Streams that are liquid and phase won't disappear          / /
VapStr(Str)           Streams that are vapor and phase won't disappear           /SF0, SF10/

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


