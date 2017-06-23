* Bounds for global optimization

*** Bounds for ASU Simple ***

H.lo(Str) = -20;
H.up(Str) = 150;

Pd.lo(Str) = 0;
Pd.up(Str) = 15000;

Pb.lo(Str) = 0;
Pd.up(Str) = 15000;

Pvap.up(Str,Comp) = 15000;

K.lo(ThrmE,Comp) = 0;
K.up(ThrmE,Comp) = 15000;

SNEd1.up(Ed1,Comp) = 10;
S1Ed1.up(Ed1,Comp) = 10;

DummyAe.up(Ed1,Comp) = 10E8;
DummySe.up(Ed1,Comp) = 10E8;

Tin.lo(GnrlE) = Tmin - 2;
Tout.lo(GnrlE) = Tmin - 2;
TPcnd.lo(Str, HIZone) = Tmin - 2;
Tout.up(GnrlE) = 350 + 2;
Tin.up(GnrlE) = 350 + 2;
TPcnd.up(Str, HIZone) = 350 + 2;

QsZ.up(HIZone) = 100;
QwZ.up(HIZone) = 100;
Qs.up = 100;
Qw.up = 100;

QAhU.lo(Str, HIZone) = 0;
QAhU.up(Str, HIZone) = 100;

QAcU.lo(Str, HIZone) = 0;
QAcU.up(Str, HIZone) = 100;

FCp.up(GnrlE) = 10;

liqPen.up = 1000;
vapPen.up = 1000;

CmprPwr.up(FeedStr) = 10;

*** Bounds for ASU CEOS ***

*HIG.lo(Str) = -10;
*HIG.up(Str) = 10;

*bmEOS.up(Str) = 1;
*bbEOS.up(Str) = 1;
*aEOS.up(Str,Comp) = 5;
*amEOS.up(Str) = 5;
*aaEOS.up(Str) = 1;

*dadT.lo(Str) = -0.1;
*dadT.up(Str) = 0.1;

*delta.lo(Str,Comp) = 0;
*delta.up(Str,Comp) = 3;

*bRatio.lo(Str,Comp) = 0;
*bRatio.up(Str,Comp) = 2;

*Inter1.up(Str) = 1.5;
*Inter2.up(Str) = 1.5;
*Inter3.up(Str) = 10;
