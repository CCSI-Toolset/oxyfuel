************************************
***** ProcessShadowStreams.gms *****
************************************
* This script activates shadow streams based on the real streams in
* BubPoint, DewPoint and PhaseStability

***** Step 1 *****
* Remove unnecessary streams from PhaseStability

* "Anti-propagate" PhaseStability through splitters
* Transfer membership in PhaseStability from outlet splitter streams to inlet
* splitter streams
loop(Str$PhaseStability(Str),
  loop(Sptr$OutSptr(Sptr,Str),
    PhaseStability(Str) = no;
    PhaseStability(Str2)$InSptr(Sptr,Str2) = yes;
  );
);

* Ensure PhaseStability is not considered for any outlet splitter streams
loop(Str$PhaseStability(Str),
  loop(Sptr$InSptr(Sptr,Str),
    PhaseStability(Str2)$OutSptr(Sptr,Str2) = no;
  );
);

* Remove any streams in BubPoint or DewPoint from PhaseStability. The streams
* in PhaseStability may have changed due to splitter processing
PhaseStability(Str)$(BubPoint(Str) OR DewPoint(Str)) = no;

**** Step 2 *****
* Populate ActShdStr

* Start by deactivating all shadow streams
ActShdStr(Str) = no;
PhiCalc(LiqShd) = no;
PhiCalc(VapShd) = no;

* Loop over DewPoint and BubPoint sets
loop(Str$(DewPoint(Str) OR BubPoint(Str)),
  loop(Str2,
    loop(Str3$((ComMap(Str,Str2,Str3) AND BubPoint(Str)) OR (ComMap(Str,Str3,Str2) AND DewPoint(Str))),
      ActShdStr(Str2) = yes;
    );
  );
);

* Loop over PhaseStability
loop(Str$PhaseStability(Str),
  loop(LiqShd,
    loop(VapShd$ComMap(Str,VapShd,LiqShd),
      ActShdStr(LiqShd) = yes;
      ActShdStr(VapShd) = yes;
    );
  );
);

***** Step 3 *****
* Turn on shadow streams for CEOS model

* Activate members of ActShdStr
fEOSStr(ActShdStr) = yes;
LiqStr(Str)$(ActShdStr(Str) AND LiqShd(Str)) = yes;
VapStr(Str)$(ActShdStr(Str) AND VapShd(Str)) = yes;
PhiCalc(ActShdStr) = yes;

display ActShdStr, DewPoint, BubPoint, PhaseStability;
