*******************************************
***** FlowsheetTopologyProcessing.gms *****
*******************************************

* This file contains all of the flowsheet topology information. In principle
* only this file needs to be changed to consider a new flowsheet. (One may
* also find tweaking the initialization procedures improves performance of
* the NLP solvers.)

* Notes: This version (FlowsheetTopology_23.gms) includes extra heat exchangers.

Sets

*-------------------------------------
* Flowsheet topology
*-------------------------------------

***** Cascades *****

***** Bubble and Dew Point Calculations *****
***** Bubble and Dew Point Calculations *****
* For every bubble/dew point calculation with the CEOS thermo there are
* one or two extra shadow streams. These extra streams are required to ensure the
* real process stream is in the two phase region. The fake (shadow) stream(s)
* is used for equilibrium calculations.

*** Common Mapping for DP and BP calculations
* This is a common mapping between actual, shadow vapor and shadow liquid streams.
* ComMap is used to populate BPStr/DPStr and BPMap/DPMap
ComMap(Str,Str2,Str3)    Common mapping for shadow streams                       / /

*** Type 1: Process streams at either at their bubble or dew point
BubPoint(Str)            Streams at their bubble point                           / /
DewPoint(Str)            Stream at their dew point                               / /

*** Type 2: Constraint process stream using bubble or dew point
* Phase stability analysis is only considered for certain streams to maintain a small
* problem size
PhaseStability(Str)      Consider phase stability analysis for these streams             / /

***** Thermo and General Equipment *****
* All of these streams are automatically populated below.
InGnrlE(GnrlE, Str)                                  / /
OutVGnrlE(GnrlE, Str)                                / /
OutLGnrlE(GnrlE, Str)                                / /
InOneGnrlE(GnrlE, Str)                               / /

* These two set are used to prevent redundant equations. If the outlet streams from
* one piece of thermo equipment are both inlets for a different piece of
* equipment, pressure (and some temperature) constraints should only be
* written for one of the streams. Both sets manage this.
InThrmEPres(ThrmE,Str)   Inlet streams to write pressure constraints for / /

* This set is used to keep track of the stream pairs for every ThrmE unit.
* StrPairs(Str,Str2)       Outlet stream pairs     / /

* Energy balances are only considered for ThrmE in this set. In the simple thermo
* model (used for initialization) a correlation is used to calculate Joule-Thomson
* cooling in valves. As a result an energy balance should not be considered
* around the valves initially.
CalcGnrlE(GnrlE)                             / /

* This set is used to prune inactive equipment.
InactiveGnrlE(GnrlE) Used for pruning        / /
InactiveCsc(Ed1)     Used for pruning        / /;

************************************
***** Automated Set Population *****
************************************

CalcGnrlE(PReb) = yes;
CalcGnrlE(Flsh) = yes;
CalcGnrlE(HtEx) = yes;
CalcGnrlE(TCond) = yes;

* Thermo Equip Connectivity
loop(ThrmE,
* Reboilers
  InGnrlE(PReb,Str)$(InPReb(PReb,Str)) = yes;
  OutVGnrlE(PReb,Str)$(OutVPReb(PReb,Str)) = yes;
  OutLGnrlE(PReb,Str)$(OutLPReb(PReb,Str)) = yes;

* Flashes
  InGnrlE(Flsh,Str)$(InFlsh(Flsh,Str)) = yes;
  OutVGnrlE(Flsh,Str)$(OutVFlsh(Flsh,Str)) = yes;
  OutLGnrlE(Flsh,Str)$(OutLFlsh(Flsh,Str)) = yes;

* Valves
  InGnrlE(Valve,Str)$(InValve(Valve,Str)) = yes;
  OutVGnrlE(Valve,Str)$(OutVValve(Valve,Str)) = yes;
  OutLGnrlE(Valve,Str)$(OutLValve(Valve,Str)) = yes;

* HtEx
  InGnrlE(HtEx,Str)$(InHtEx(HtEx,Str)) = yes;
  OutVGnrlE(HtEx,Str)$(OutVHtEx(HtEx,Str)) = yes;
  OutLGnrlE(HtEx,Str)$(OutLHtEx(HtEx,Str)) = yes;

);

loop(GnrlE,
* Condensers
  InGnrlE(TCond,Str)$InTCond(TCond,Str) = yes;
  OutLGnrlE(TCond,Str)$OutTCond(TCond,Str) = yes;

* Set stream pairs
  StrPairs(Str,Str2)$(OutVGnrlE(GnrlE,Str) AND OutLGnrlE(GnrlE,Str2)) = yes;
);

* Initially populate InThrmEPres with all of the inlet streams for all thermo
* equipment.
InThrmEPres(ThrmE,Str)$(InGnrlE(ThrmE,Str)) = yes;

* Remove the liquid stream in each pair (leaving ThrmE A) if the
* vapor and liquid stream are both feeds for a different unit (feeds for ThrmE B)
loop(ThrmE,
  loop(Str$InGnrlE(ThrmE,Str),
    loop(Str2$InGnrlE(ThrmE,Str2),
      if(StrPairs(Str,Str2),
        InThrmEPres(ThrmE,Str2) = no;
      );
    );
  );
);

*loop(ThrmE,
*  loop(Str$(OutVGnrlE(ThrmE,Str) OR OutLGnrlE(ThrmE,Str)),
*   Collect all of the outlet thrme streams
*    OutThrmEStr(Str) = yes;
*  );
*);

* Usage to remove recycles. Consider StrPairs(1,2). If 1 is an outlet for unit A,
* don't consider 2 as an inlet.

loop(ThrmE,
  loop(Str$(OutLGnrlE(ThrmE,Str) OR OutVGnrlE(ThrmE,Str)),
    loop(Str2$StrPairs(Str,Str2),
      InThrmEPres(ThrmE,Str2) = no;
    );
  );
);

* Display the automatically populated sets for debugging purposes.
display InGnrlE, OutVGnrlE, OutLGnrlE, StrPairs, InThrmEPres;

***** General Equipment *****

loop(TCond,
  InGnrlE(TCond,Str)$InTCond(TCond,Str) = yes;
  OutLGnrlE(TCond,Str)$OutTCond(TCond,Str) = yes;
);

loop(GnrlE,
  loop(Str$InGnrlE(GnrlE, Str),
    if(SUM(Str2$InOneGnrlE(GnrlE, Str2), 1) < 1,
      InOneGnrlE(GnrlE, Str) = yes;
    );
  );
);

*display InGnrlE, OutVGnrlE, OutLGnrlE, InOneGnrlE;

***** Populate Real Stream *****

* General equipment
loop(GnrlE,
  RealStr(Str)$(InGnrlE(GnrlE,Str) OR OutLGnrlE(GnrlE, Str) OR OutVGnrlE(GnrlE, Str)) = yes;
);

* Cascade
RealStr(CscStr) = yes;

* Splitters
loop(Sptr,
  RealStr(Str)$(InSptr(Sptr,Str) OR OutSptr(Sptr,Str)) = yes;
);

***** Shadow Streams *****

Scalar StrCount          /1/;

* Assign every RealStr corresponding shadow streams
loop(Str$RealStr(Str),
  loop(VapShd$(ord(VapShd) eq StrCount),
    loop(LiqShd$(ord(LiqShd) eq StrCount),
      ComMap(Str,VapShd,LiqShd) = yes;
    );
  );
  StrCount = StrCount + 1;
);

* Initially constrain all cascade outlets to be at their bubble or dew point
loop(Ed1,
  DewPoint(Str)$OutVEd1(Ed1,Str) = yes;
  BubPoint(Str)$OutLEd1(Ed1,Str) = yes;
);

* Propagate DewPoint and BubPoint through splitters
loop(Str$DewPoint(Str),
  loop(Sptr$InSptr(Sptr,Str),
    DewPoint(Str2)$OutSptr(Sptr,Str2) = yes;
  );
);

loop(Str$BubPoint(Str),
  loop(Sptr$InSptr(Sptr,Str),
    BubPoint(Str2)$OutSptr(Sptr,Str2) = yes;
  );
);

* Initially consider phase stability for all single phase streams (LiqStr, VapStr)
* that are NOT already constrained to be at their bubble or dew point

PhaseStability(Str)$(LiqStr(Str) AND NOT BubPoint(Str)) = yes;
PhaseStability(Str)$(VapStr(Str) AND NOT DewPoint(Str)) = yes;

* Also consider all total condenser outputs for phase stability
loop(TCond,
  PhaseStability(Str)$(OutTCond(TCond,Str) AND NOT (BubPoint(Str) OR DewPoint(Str))) = yes;
);



* Last Step - Run ProcessShadowStreams.gms
* This is done at the end. PhiCalc must be declared first

************************************************
***** Pinch Location Heat Integration Sets *****
************************************************

* In this section the sets required for pinch location heat integration
* are defined and managed.

Sets
PStr(Str,HIZone)                   pinch candidate streams        / /
InPStr(Str,HIZone)                 inlet streams in heat integration / /
OutPStr(Str,HIZone)                outlet streams in heat integration / /;

Sets
  CEqp(GnrlE)            Equipment that NEEDS cooling (hot streams)         / /
  HEqp(GnrlE)             Equipment that NEEDS heating (cold streams)        / /;

$include ../ProcsFiles/HeatIntegrationProcessing.gms

***********************************************
***** Thermodynamic Sets and Calculations *****
***********************************************

* Automatically assign phases to streams based on flowsheet connectivity. Here
* are the rules:
*
* 1. Except for the feed streams, all other streams are
* considered "vanishable". This allows for units to disappear in the superstructure.
*
* 2. Outlet streams for splitters "inhert" their stream classification for the
* splitters single inlet stream (not yet)

loop(GnrlE,
  FlashVap(Str)$OutVGnrlE(GnrlE, Str) = yes;
  FlashLiq(Str)$OutLGnrlE(GnrlE, Str) = yes;
);

loop(Ed1,
  FlashLiq(Str)$(InLEd1(Ed1,Str) OR OutLEd1(Ed1,Str)) = yes;
  FlashVap(Str)$(InVEd1(Ed1,Str) OR OutVEd1(Ed1,Str)) = yes;
);

* Primarily used for the CEOS, these sets manage which thermodynamic
* equations should be evaulated for each stream.

Sets
HCalc(Str)            Streams to calculate enthalpy                              / /
SCalc(Str)            Streams to calculate entropy                               / /
VCalc(Str)            Streams to calculate specific volume                       / /
PhiCalc(Str)          Streams to calculate fugacity                              / /
StrModCalc(Str)       Streams to evaluate stream model                           / /;

* Heat capacity is only required for valves (and only for the simple thermo model)
loop(Valve,
  CpStr(Str)$OutLValve(Valve,Str) = yes;
  CpStr(Str)$InValve(Valve,Str) = yes;
);

* Populate FlashVap and FlashLiq with ThrmE outlets
* Calculate enthalpy for GnrlE inlets and outlets (energy balance)
loop(GnrlE,
  HCalc(Str)$(InGnrlE(GnrlE,Str) OR OutVGnrlE(GnrlE,Str) OR OutLGnrlE(GnrlE,Str)) = yes;
);

* Calculate fugacity for ThrmE outlets (equilibrium)
loop(ThrmE,
  PhiCalc(Str)$(OutVGnrlE(ThrmE,Str) OR OutLGnrlE(ThrmE,Str)) = yes;
);

* Calculate fugacity for cascade outlets (equilibrium - BP/DP calcs)
loop(Ed1,
  PhiCalc(Str)$(OutVEd1(Ed1,Str) OR OutLEd1(Ed1,Str)) = yes;
);

* Calculate enthalpy for cascade streams (energy balance)
HCalc(Str)$CscStr(Str) = yes;

* Calculate fugacity for shadow streams (equilibrium - BP/DP calcs)
PhiCalc(Str)$ActShdStr(Str) = yes;

* Copy all RealStr into StrModCalc
StrModCalc(RealStr) = yes;

Set CurrentSptr(Sptr)    / /;

* Propogate set membership from outlet splitter streams to inlet streams
loop(Sptr,
  loop(Str$OutSptr(Sptr,Str),
    if(HCalc(Str),
      HCalc(Str2)$InSptr(Sptr,Str2) = yes;
    );
    if(SCalc(Str),
      SCalc(Str2)$InSptr(Sptr,Str2) = yes;
    );
    if(PhiCalc(Str),
      PhiCalc(Str2)$InSptr(Sptr,Str2) = yes;
    );
  );
  HCalc(Str)$OutSptr(Sptr,Str) = no;
  PhiCalc(Str)$OutSptr(Sptr,Str) = no;
  SCalc(Str)$OutSptr(Sptr,Str) = no;
  StrModCalc(Str)$OutSptr(Sptr,Str) = no;
);

$Offwarning

* Propogate phase through splitter
loop(Sptr,

* Step 1: Propogate phases from outlet to inlet. Be mindful of the following rules:
* A. Prexisting phase designation for inlet trumps any inherited phase designation from outlet
* B. When inheriting, VapStr or LiqStr trumps FlashVap/FlashLiq (in the event outlet streams have multiple phase designations)
  loop(Str$InSptr(Sptr,Str),
    if(LiqStr(Str) OR FlashLiq(Str) OR VapStr(Str) OR FlashVap(Str),
* Do nothing. Inlet has a pre-existing phase specification
    else
      loop(Str2$OutSptr(Sptr,Str2),
* An outlet stream is liquid-only
        if(LiqStr(Str2),
          FlashLiq(Str) = no;
          LiqStr(Str) = yes;
* An outlet stream is a liquid that may disappear, and the inlet stream is not liquid-only
        elseif FlashLiq(Str2) AND NOT LiqStr(Str),
          FlashLiq(Str) = yes;
        );
* An outlet stream is vapor-only
       if(VapStr(Str2),
          FlashVap(Str) = no;
          VapStr(Str) = yes;
* An outlet stream is a vapor that may disappear, and the inlet stream is not vapor only
       elseif FlashVap(Str2) AND NOT VapStr(Str),
          FlashVap(Str) = yes;
        );
      );
    );
  );

* Step 2: Check for pre-existing inconsistencies
  CurrentSptr(Sptr) = no;
  loop(Str$OutSptr(Sptr,Str),
    loop(Str2$InSptr(Sptr,Str),
      if( (LiqStr(Str) OR FlashLiq(Str)) AND (FlashVap(Str2) OR VapStr(Str2)),
        CurrentSptr(Sptr) = yes;
        display CurrentSptr;
        display "Check splitters. Case of LiqStr/FlashLiq inlet and VapStr/FlashVap outlet detected";
      elseif (LiqStr(Str2) OR FlashLiq(Str2)) AND (FlashVap(Str) OR VapStr(Str)),
        CurrentSptr(Sptr) = yes;
        display CurrentSptr;
        display "Check splitters. Case of LiqStr/FlashLiq outlet and VapStr/FlashVap inlet detected";
      );
    );
  );

* Step 3: Propogate inlet phase specification to outlets. Important rules:
* A: Non-disappearing phases trump phases that may disappear
* B: There is no priority for pre-existing outlet stream phase specifications

  loop(Str$InSptr(Sptr,Str),
    if(LiqStr(Str),
      LiqStr(Str2)$OutSptr(Sptr,Str2) = yes;
      FlashLiq(Str2)$OutSptr(Sptr,Str2) = no;
    elseif FlashLiq(Str),
      FlashLiq(Str2)$OutSptr(Sptr,Str2) = yes;
    );
    if(VapStr(Str),
      VapStr(Str2)$OutSptr(Sptr,Str2) = yes;
      FlashVap(Str2)$OutSptr(Sptr,Str2) = no;
    elseif FlashVap(Str),
      FlashVap(Str2)$OutSptr(Sptr,Str2) = yes;
    );
  );

);

display LiqStr, VapStr, FlashLiq, FlashVap, HCalc, SCalc, VCalc, PhiCalc;

* Populate fEOSStr
fEOSStr(FlashVap) = yes;
fEOSStr(FlashLiq) = yes;
fEOSStr(LiqStr) = yes;
fEOSStr(VapStr) = yes;

display fEOSStr;

* For some reason these commented out lines result in ASU_CEOS9b2 solving to
* one of the higher objective function value local points with the
* ASU_Simple6_good1.gdx starting point.
*fEOSStr(HCalc) = yes;
*fEOSStr(PhiCalc) = yes;

* Remove exit casacde stream from InThrmEPres(Flsh,Str)
* This prevents a redundent pressure equation
loop(Ed1,
  loop(Flsh,
    InThrmEPres(Flsh,Str)$(OutLEd1(Ed1,Str) OR OutVEd1(Ed1,Str)) = no;
  );
);

* Remove splitter outlet streams from fEOSStr
loop(Sptr,
  loop(Str$OutSptr(Sptr,Str),
    fEOSStr(Str) = no;
  );
);

* Finish handling shadow streams
$include ../ProcsFiles/ProcessShadowStreams.gms

********************************************
***** Determine if a cascade is unused *****
********************************************

* If a cascade doesn't have any inlet or outlet streams, remove it
InactiveCsc(Ed1)$(SUM(Str$(InVEd1(Ed1,Str) OR InLEd1(Ed1,Str) OR OutVEd1(Ed1,Str) OR OutLEd1(Ed1,Str)), 1) + 0 < 1) = yes;

