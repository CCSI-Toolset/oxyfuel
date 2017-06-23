*****************************
***** Init_Cascade2.gms *****
*****************************

* This script attempts to initialize
* the approximate cascade in an iterative fashion

*StgEd1.l(Ed1) = 10;

* Assume guesses for PhiAEd1 and PhiSEd1 are already in memory

count = 0;

while(count < 3,

***** Calculate vapor outlet flows from performance equation
if(count > 0,
loop(Ed1,
  Fc.l(Str,Comp)$OutVEd1(Ed1,Str) =
         SUM(Str2$InVEd1(Ed1,Str2), Fc.l(Str2,Comp)*PhiAEd1.l(Ed1,Comp))
         + SUM(Str2$InLEd1(Ed1,Str2), Fc.l(Str2,Comp)*(1 - PhiSEd1.l(Ed1,Comp)));

***** Calculate liquid outlet flows from mass balance
  Fc.l(Str,Comp)$OutLEd1(Ed1,Str) =
         SUM(Str2$(InVEd1(Ed1,Str2) OR InLEd1(Ed1,Str2)), Fc.l(Str2,Comp))
         - SUM(Str2$OutVEd1(Ed1,Str2), Fc.l(Str2,Comp));
);
);

***** Guess/Calculate L1 and VN *****

  VNEd1.l(Ed1) = SUM(Str $InVEd1(Ed1,Str), F.l(Str));
  L1Ed1.l(Ed1) =  SUM(Str $InLEd1(Ed1,Str), F.l(Str));

***** Calculate A1Ed1, ANEd1 *****
if( ThermoSwitch eq 1,
  A1Ed1.l(Ed1,Comp) =  Max(L1Ed1.l(Ed1)* SUM(Str $OutVEd1(Ed1,Str), P.l(Str))
         /( SUM(Str $OutVEd1(Ed1,Str), Pvap.l(Str,Comp)*(Max(F.l(Str),epsi)))), epsi) ;

  ANEd1.l(Ed1,Comp) =  Max(SUM(Str $OutLEd1(Ed1,Str), F.l(Str) * P.l(Str))
         /( SUM(Str $OutLEd1(Ed1,Str),  Pvap.l(Str,Comp))*(Max(VNEd1.l(Ed1),epsi))),epsi) ;

elseif ThermoSwitch eq 2,
  A1Ed1.l(Ed1,Comp) =  SUM(Str$OutVEd1(Ed1,Str), SUM(Str2$DPStr(Str,Str2),
         Max(L1Ed1.l(Ed1)*PhiEOS.l(Str,Comp)/( Max(PhiEOS.l(Str2,Comp)*F.l(Str),epsi)), epsi ))) ;

  ANEd1.l(Ed1,Comp) =  SUM(Str$OutLEd1(Ed1,Str), SUM(Str2$BPStr(Str,Str2),
         Max(F.l(Str)*PhiEOS.l(Str2,Comp)/( Max(VNEd1.l(Ed1)*PhiEOS.l(Str,Comp),epsi)), epsi ))) ;
else
  A1Ed1.l(Ed1,Comp) = 1;
  ANEd1.l(Ed1,Comp) = 1;
);

***** Calculate S1Ed1, SNEd1 *****
S1Ed1.l(Ed1,Comp) = 1/A1Ed1.l(Ed1,Comp);
SNEd1.l(Ed1,Comp) = 1/ANEd1.l(Ed1,Comp);

***** Calculate AeEd1, SeEd1 *****
AeEd1.l(Ed1,Comp) =  sqrt(ANEd1.l(Ed1,Comp)*(A1Ed1.l(Ed1,Comp)+1) + 0.25) - 0.5 ;
SeEd1.l(Ed1,Comp) = sqrt(S1Ed1.l(Ed1,Comp)*(SNEd1.l(Ed1,Comp)+1) + 0.25) - 0.5 ;

***** Update guess for PhiAEd1, PhiSEd1 *****
PhiAEd1.l(Ed1,Comp) = (AeEd1.l(Ed1,Comp) - 1)/(AeEd1.l(Ed1,Comp)**(StgEd1.l(Ed1)+1) - 1) ;
PhiSEd1.l(Ed1,Comp) = (SeEd1.l(Ed1,Comp) - 1)/(SeEd1.l(Ed1,Comp)**(StgEd1.l(Ed1)+1) - 1) ;

count = count + 1;

);

***** Initialize intermediates used in alternate models
DummyAe.l(Ed1,Comp) = min(AeEd1.l(Ed1,Comp)**(StgEd1.l(Ed1)+1), 1E6);
DummySe.l(Ed1,Comp) = min(SeEd1.l(Ed1,Comp)**(StgEd1.l(Ed1)+1), 1E6);

loop(Ed1,
  loop(Str$InVEd1(Ed1,Str),
    loop(Str2$OutVEd1(Ed1,Str2),
      loop(Str3$InLEd1(Ed1,Str3),
        CascInter.l(Ed1,Comp) = (Fc.l(Str,Comp)*(1 - AeEd1.l(Ed1,Comp)) - Fc.l(Str2,Comp)*Fc.l(Str3,Comp))/
          max(Fc.l(Str3,Comp) - Fc.l(Str2,Comp)* AeEd1.l(Ed1,Comp),1E-4);
*        CascInter.l(Ed1,Comp) = Power(Fc.l(Str3,Comp) - Fc.l(Str2,Comp)* AeEd1.l(Ed1,Comp),2);
*        CascInter.l(Ed1,Comp) = StgEd1.l(Ed1)*log(AeEd1.l(Ed1,Comp))
      );
    );
  );
);

CascInterA.l(Ed1,Comp) = AeEd1.l(Ed1,Comp) - 1 + PhiAEd1.l(Ed1,Comp);
CascInterS.l(Ed1,Comp) = SeEd1.l(Ed1,Comp) - 1 + PhiSEd1.l(Ed1,Comp);

***** Initialize pressure drop *****
PresDrop.l(Ed1) = StgEd1.l(Ed1)*PresDropCnst.l(Ed1);

***** Ensure cascade streams are not relaxed *****
sL.fx(Str)$CscStr(Str) = 0;
sV.fx(Str)$CscStr(Str) = 0;
