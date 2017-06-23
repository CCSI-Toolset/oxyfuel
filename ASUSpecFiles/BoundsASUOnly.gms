*o2rec.lo = 0.05;
o2rec.up = 1;

o2pure.lo = o2pureSpec;

o2pure.up = 1;

*if(ThermoSwitch eq 1,
*  F.lo('S12') = 0.1;
  F.lo('S33') = 0.1;
*else
*  F.lo('S12') = 0;
*  F.lo('S3') = 0;
*  F.lo('SF0') = 0;
*  F.lo('SF10') = 0;
*);
  F.lo('S13') = 0.1;

Tscaled.lo(FeedStr) = 300/Tref;


