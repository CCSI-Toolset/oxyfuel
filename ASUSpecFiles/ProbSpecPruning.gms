**********************************************
***** Problem Specific Flowsheet Pruning *****
**********************************************

Sets
  Config0(Str)   /S17, SOV10*SOV12, SOL10*SOL12/
  Config1(Str)   /SSVt1*SSVt3, SSLt1*SSLt3, SSVb1*SSVb3, SSLb1*SSLb3/
  Config2(Str)   /SSVt1*SSVt2, SSLt1*SSLt2, SSVb2*SSVb3, SSLb2*SSLb3, SSVm1, SSVm3, SSLm1, SSLm3, S4, S13, S36, S26, S27/
  Config3(Str)   /SSVt1*SSVt2, SSLt1*SSLt2, SSVb2*SSVb3, SSLb2*SSLb3, SSVm1, SSVm3, SSLm1, SSLm3, S4, S12, S30, S34, SNV11, SNL11, SNV12, SNL12/
  Config4(Str)   Remove C2       /S24, S31/
  Config5(Str)   Possible ideal configuration            /SSVt1*SSVt3, SSLt1*SSLt3, SSLb2, SSVb2, SSLm3, SSVm3, S24, S31/;


*Set
*  ManualRemove(Str)      This set contains streams that should be manually removed. This forces certain flowsheet configurations        /  /;

if(PruneConfig eq 0,
  ManualRemove(Str)$Config0(Str) = yes;
elseif PruneConfig eq 1,
  ManualRemove(Str)$Config0(Str) = yes;
  ManualRemove(Str)$Config1(Str) = yes;
elseif PruneConfig eq 2,
  ManualRemove(Str)$Config0(Str) = yes;
  ManualRemove(Str)$Config2(Str) = yes;
elseif PruneConfig eq 3,
  ManualRemove(Str)$Config0(Str) = yes;
  ManualRemove(Str)$Config3(Str) = yes;
elseif PruneConfig eq 4,
  ManualRemove(Str)$Config0(Str) = yes;
  ManualRemove(Str)$Config4(Str) = yes;
elseif PruneConfig eq 5,
  ManualRemove(Str)$Config0(Str) = yes;
  ManualRemove(Str)$Config5(Str) = yes;
);

***** End Problem Specific Flowsheet Pruning *****