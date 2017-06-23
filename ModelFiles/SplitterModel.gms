*****************************
***** SplitterModel.gms *****
*****************************

* This file includes the splitter model. This model structure assumes there
* is only one inlet stream for the splitter.

*----------------------------------------*
*  Equation for splitter
*----------------------------------------*
Equations
*EqCompMolBalSptr(Sptr,Str,Comp)           Component balance in splitter (not used)
EqMolBalSptr(Sptr)                        Splitter mole balance
EqMolFracSptr(Sptr,Str,Str2,Comp)         Splitter split propogation
EqPRel1Sptr(Sptr,Str)                     Pressure relation
EqTRel1Sptr(Sptr,Str)                     Temperature relation
EqCopyHSptr(Sptr,Str,Str2)                Copy enthalpy from inlet to outlet
EqCopyPhiSptr(Sptr,Str,Str2,Comp)         Copy fugacity from inlet to outlet
EqCopyZSptr(Sptr,Str,Str2)                Copy ZEOS from inlet to outlet
EqCopyPvapSptr(Sptr,Str,Str2,Comp)        Copy vapor pressure from inlet to outlet
EqCopySSptr(Sptr,Str,Str2)                Copy entrop from inlet to outlet
;

*EqCompMolBalSptr(Sptr,Str,Comp)$OutSptr(Sptr,Str)..
*     Fc(Str,Comp)  =E= Sum(Str2 $InSptr(Sptr,Str2),Fc(Str2,Comp)*SptrFrac(Str));

EqMolBalSptr(Sptr)..
     SUM(Str$InSptr(Sptr,Str), F(Str)) =e= SUM(Str$OutSptr(Sptr,Str), F(Str));

EqMolFracSptr(Sptr,Str,Str2,Comp)$(InSptr(Sptr,Str) AND OutSptr(Sptr,Str2) AND ord(Comp) < card(Comp))..
     Xc(Str,Comp) =e= Xc(Str2,Comp);

EqPRel1Sptr(Sptr,Str)$OutSptr(Sptr,Str)..
        SUM(Str2 $InSptr(Sptr,Str2), P(Str2)) =E= P(Str) ;

EqTRel1Sptr(Sptr,Str)$OutSptr(Sptr,Str)..
        SUM(Str2 $InSptr(Sptr,Str2), Tscaled(Str2)) =E= Tscaled(Str) ;

EqCopyHSptr(Sptr,Str,Str2)$(InSptr(Sptr,Str) AND OutSptr(Sptr,Str2))..
  H(Str) =e= H(Str2);

EqCopyPhiSptr(Sptr,Str,Str2,Comp)$(InSptr(Sptr,Str) AND OutSptr(Sptr,Str2))..
  phiEOS(Str,Comp) =e= phiEOS(Str2,Comp);

EqCopyZSptr(Sptr,Str,Str2)$(InSptr(Sptr,Str) AND OutSptr(Sptr,Str2))..
  ZEOS(Str) =e= ZEOS(Str2);

EqCopyPvapSptr(Sptr,Str,Str2,Comp)$(InSptr(Sptr,Str) AND OutSptr(Sptr,Str2))..
  Pvap(Str,Comp) =e= Pvap(Str2,Comp);

EqCopySSptr(Sptr,Str,Str2)$(InSptr(Sptr,Str) AND OutSptr(Sptr,Str2) AND SCalc(Str))..
  S(Str) =e= S(Str2);

Model EqSpltr /EqMolBalSptr, EqMolFracSptr, EqPRel1Sptr, EqTRel1Sptr/;

* Should this include EqCopyPvapSptr?
Model EqSpltrSimpleThermo /EqCopyHSptr/;

Model EqSpltrCEOSThermo /EqCopyHSptr, EqCopyPhiSptr, EqCopyZSptr, EqCopySSptr/;
