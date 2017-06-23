* Empty Sets
PStr(Str,HIZone) = no;
InPStr(Str,HIZone) = no;
OutPStr(Str,HIZone) = no;

* These loops automatically populate the pinch candidate streams by iterating
* through the connectivity sets for each type of equipment considered for
* heat integration.
loop(HIZone,
  loop(GnrlE$HIMap(GnrlE, HIZone),
    loop(Str$(InOneGnrlE(GnrlE,Str)),
      InPStr(Str,HIZone) = yes;
    );
  );
);

loop(HIZone,
  loop(GnrlE$HIMap(GnrlE, HIZone),
    loop(Str$(OutLGnrlE(GnrlE,Str) AND NOT InPStr(Str, HIZone)),
      OutPStr(Str,HIZone) = yes;
    );
  );
);

PStr(Str,HIZone)$InPStr(Str,HIZone) = yes;
PStr(Str,HIZone)$OutPStr(Str,HIZone) = yes;

display PStr, InPStr, OutPStr;

CEqp(CoolHtEx) = yes;
CEqp(TCond) = yes;
HEqp(HeatHtEx) = yes;
HEqp(PReb) = yes;