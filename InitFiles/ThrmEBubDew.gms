***************************
***** ThrmEBubDew.gms *****
***************************

* This file fixes slack and flow variables to restrict specified ThrmE outlet
* streams to be at their bubble or dew point. This is useful for locating
* phaes changes for heat integration.

***** Restrict dew point outlets *****
loop(ThrmE$ThrmEOutDewPoint(ThrmE),
  loop(Str$OutLGnrlE(ThrmE,Str),
    sL.fx(Str) = 0;
    F.fx(Str) = 0;
    Fc.l(Str,Comp) = 0;
  );
);
