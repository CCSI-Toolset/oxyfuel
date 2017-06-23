*StgEd1.l(Ed1) = 10;

if( ThermoSwitch eq 1,
  A1Ed1.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) =  Max(L1Ed1.l(Ed1)* SUM(Str $OutVEd1(Ed1,Str), P.l(Str))
         /( SUM(Str $OutVEd1(Ed1,Str), Pvap.l(Str,Comp)*(Max(F.l(Str),epsi)))), epsi) ;

  ANEd1.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) =  Max(SUM(Str $OutLEd1(Ed1,Str), F.l(Str) * P.l(Str))
         /( SUM(Str $OutLEd1(Ed1,Str),  Pvap.l(Str,Comp))*(Max(VNEd1.l(Ed1),epsi))),epsi) ;

*** Adding this for ASU23 ***


elseif ThermoSwitch eq 2,
  A1Ed1.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) =  SUM(Str$OutVEd1(Ed1,Str), SUM(Str2, SUM(Str3$ComMap(Str,Str3,Str2),
         Max(L1Ed1.l(Ed1)*PhiEOS.l(Str,Comp)/( Max(PhiEOS.l(Str2,Comp)*F.l(Str),epsi)), epsi )))) ;

  ANEd1.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) =  SUM(Str$OutLEd1(Ed1,Str), SUM(Str2, SUM(Str3$ComMap(Str,Str2,Str3),
         Max(F.l(Str)*PhiEOS.l(Str2,Comp)/( Max(VNEd1.l(Ed1)*PhiEOS.l(Str,Comp),epsi)), epsi ))) ) ;
else
  A1Ed1.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) = 1;
  ANEd1.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) = 1;
);


S1Ed1.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) = 1/A1Ed1.l(Ed1,Comp);
SNEd1.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) = 1/ANEd1.l(Ed1,Comp);

*display A1Ed1.l, ANEd1.l, S1Ed1.l, SNEd1.l;

AeEd1.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) =  sqrt(ANEd1.l(Ed1,Comp)*(A1Ed1.l(Ed1,Comp)+1) + 0.25) - 0.5 ;
SeEd1.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) = sqrt(S1Ed1.l(Ed1,Comp)*(SNEd1.l(Ed1,Comp)+1) + 0.25) - 0.5 ;

*display AeEd1.l, SeEd1.l;

PhiAEd1.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) = (AeEd1.l(Ed1,Comp) - 1)/(AeEd1.l(Ed1,Comp)**(StgEd1.l(Ed1)+1) - 1) ;
PhiSEd1.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) = (SeEd1.l(Ed1,Comp) - 1)/(SeEd1.l(Ed1,Comp)**(StgEd1.l(Ed1)+1) - 1) ;

*display PhiAEd1.l, PhiSEd1.l;

* Removed in an attempt to get InitCsc working
*loop(Ed1,
*  loop(Comp,
*    Fc.l(Str,Comp)$OutVEd1(Ed1,Str) = SUM(Str2$InVEd1(Ed1,Str2), Fc.l(Str2,Comp)*PhiAEd1.l(Ed1,Comp))
*                               + SUM(Str3$InLEd1(Ed1,Str3), Fc.l(Str3,Comp)*PhiSEd1.l(Ed1,Comp));
*    Fc.l(Str,Comp)$OutLEd1(Ed1,Str) = SUM(Str2$(InVEd1(Ed1,Str2) OR InLEd1(Ed1,Str2)), Fc.l(Str,Comp)) - SUM(Str3$OutVEd1(Ed1,Str3), Fc.l(Str3,Comp));
*  );
*);

DummyAe.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) = min(AeEd1.l(Ed1,Comp)**(StgEd1.l(Ed1)+1), 1E6);
DummySe.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) = min(SeEd1.l(Ed1,Comp)**(StgEd1.l(Ed1)+1), 1E6);



loop(Ed1,
  loop(Str$InVEd1(Ed1,Str),
    loop(Str2$OutVEd1(Ed1,Str2),
      loop(Str3$InLEd1(Ed1,Str3),
        CascInter.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) = (Fc.l(Str,Comp)*(1 - AeEd1.l(Ed1,Comp)) - Fc.l(Str2,Comp)*Fc.l(Str3,Comp))/
          max(Fc.l(Str3,Comp) - Fc.l(Str2,Comp)* AeEd1.l(Ed1,Comp),1E-4);
*        CascInter.l(Ed1,Comp) = Power(Fc.l(Str3,Comp) - Fc.l(Str2,Comp)* AeEd1.l(Ed1,Comp),2);
*        CascInter.l(Ed1,Comp) = StgEd1.l(Ed1)*log(AeEd1.l(Ed1,Comp))
      );
    );
  );
);

CascInterA.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) = AeEd1.l(Ed1,Comp) - 1 + PhiAEd1.l(Ed1,Comp);
CascInterS.l(Ed1,Comp)$(NOT InactiveCsc(Ed1)) = SeEd1.l(Ed1,Comp) - 1 + PhiSEd1.l(Ed1,Comp);

* Initialize pressure drop
PresDrop.l(Ed1)$(NOT InactiveCsc(Ed1)) = StgEd1.l(Ed1)*PresDropCnst.l(Ed1);

* For the Edmister model one must assume the inlet and outlet streams exist
sL.fx(Str)$CscStr(Str) = 0;
sV.fx(Str)$CscStr(Str) = 0;
