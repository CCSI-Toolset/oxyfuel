***************************
***** Thermo_Data.gms *****
***************************

* This file contains data for both the simple and CEOS thermodynamic models.

Parameters
* Critical temperature
  Tc(AllComp)                   Critical temperature  /O2 = 154.58, N2 = 126.2, Ar = 150.86/ ;
*  SelectCEOSModel            Toggle between CEOS models          /1/;

****** Toggle Between CEOS Models *****
* SelectCEOSModel = 1    ----->  Soave (SRK)
* SelectCEOSModel = 2    ----->  Soave modified by Graboski and Daubert (SRK)
* SelectCEOSModel = 3    ----->  Peng-Robinson (PR)

Scalars
Pref      Reference pressure (in bar) for thermo calculations    /1.0/
Tref      Reference temperature (K) for Enthalpy calculations    /298.15/
hscale    Enthalpy scaling                                       /1.0/
R         Ideal gas constant (m^3 bar per K-kmol)                /8.314E-2/
Rsi       Ideal gas constant (J per K-mol)                       /8.314/
bigM      Big M for 2nd derivative of CEOS                       /10/;

Alias (AllComp, AllComp2);

Table kEOS(AllComp,AllComp2)  Binary interaction parameters for PR-EOS
                 N2              O2              Ar
N2               0               -0.0119         -0.0026
O2               -0.0119         0               0.0104
Ar               -0.0026         0.0104          0;

*Table kEOS(Comp,Comp2)  Binary interaction parameters for SRK-EOS
*                 N2              O2              Ar
*N2               0               -0.0078         0
*O2               -0.0078         0               0.0178
*Ar               0               0.0178          0;


* CpIG units: J/(gmol-K)
Table CPIG (*,AllComp)  Constants for Specific heat capacity for ideal gas (from Prop. Gases & Liquids)
                 N2              O2              Ar
1                31.12896        28.087192       20.790296
2                -1.356E-02      -3.678E-06      -3.209E-05
3                2.678E-05       1.745E-05       5.163E-08
4                -1.167E-08      -1.064E-08      0 ;

* From Yannic's file
*Table CPIG(*,Comp)
*                 O2                      N2                      Ar
*5                1.00321110903601e-011   6.4863934487993e-013    0.0
*4                -3.6329587575635e-008   -6.96486792964864e-009  0.0
*3                4.20502644993578e-005   1.50552864518556e-005   0.0
*2                -0.0105472703429208     -0.00566985412879804    0.0
*1                29.7351909214832        29.5810529469761        20.786;

* Antoine constants for vapor pressure
* Properties of Gases and Liquids, Reid, Prausnitz and Sherwood, McGraw-Hill, New York, 3rd edition (1977)
Table AntConst(*,AllComp) (T in K and P in bar (after conversion))
                 N2              O2              Ar
A                14.9342         15.4075         15.2330
B                588.72          734.55          700.51
C                -6.60           -6.45           -5.84;

Table HVapSurf(*,AllComp)
                 N2                      O2                      Ar
1                -9.419638082            -9.281345395            -6.764363316
2                -0.134383294            -0.119208846            -0.117390304
3                -6.65834E-05            -6.81689E-05            -6.87623E-05
4                11.87900068             10.8997866              8.531344943
5                -6.179092207            -4.420145046            -4.307560989
6                2.969998132             2.050844214             2.002033931
7                0.266669566             0.217390761             0.216158842
8                -0.138705834            -0.106580082            -0.106878162
;

Table HLiqSurf(*,AllComp)
                 N2                      O2                      Ar
1                -23.46909382            -19.29553654            -16.47459134
2                0                       0                       0
3                0                       0                       0
4                85.6988809              32.13683952             31.42021147
5                -234.7677232            -60.50977368            -67.12027812
6                261.776955              67.02456925             74.8956398
7                0.03868993              0.017744641             0.018784856
8                -0.113036951            -0.040616087            -0.044294344
;

* liquid molar volume (m^3/kmol) (assume incompressible)
Parameter
  Vm(AllComp)       Molar Volume    /O2 = 2.8e-2, N2 = 3.7e-2, Ar = 3e-2 /;

Parameters
* Critical pressure [bar] (from Prop. Gases & Liquids)
  Pc(AllComp)      Critical pressure  /N2 = 33.943875, O2 = 50.45985, Ar = 48.737325/

* Pitzer's   accentric factor omega (from Prop. Gases & Liquids)
  omega(AllComp)  accentric factor    /N2 = 0.040, O2 = 0.021, Ar = -0.004/

* fw factor calculated from omega for CEOS
  fw(AllComp)       s factor (unique value for each CEOS)
  uEOS           constant that specifies unique CEOS
  wEOS           constant that specifies unique CEOS
  Inter0         intermediate for CEOS calculations
  omegaA         capital omega (unique value for each CEOS)

* b for individual components for PR-EOS
  bEOS(AllComp) b for component

* Molecular weight
  MWcomp(AllComp)      Molecular weight /N2 = 28.013, O2 = 31.999, Ar = 39.948/;

* Calculation of other parameters from existing parameters

* R = 8.314E-2 bar m3 / (kmol K)

* Soave Equation (SRK)
if(SelectCEOSModel eq 1,
  bEOS(AllComp) = 0.08664*8.314E-2*Tc(AllComp)/Pc(AllComp);
  uEOS = 1;
  wEOS = 0;
  omegaA = 0.42748;
  fw(AllComp) = 0.480 + 1.574*omega(AllComp) - 0.176*Power(omega(AllComp),2);
);

* Soave Equation Modified by Graboski and Daubert (SRK)
* Note: This is NOT complete
if(SelectCEOSModel eq 2,
  bEOS(AllComp) = 0.08664*8.314E-2*Tc(AllComp)/Pc(AllComp);
  uEOS = 1;
  wEOS = 0;
  omegaA = 0.42747;
  fw(AllComp) = 0.48508 + 1.55171*omega(AllComp) - 0.15613*Power(omega(AllComp),2);
);

* Peng-Robinson Equation (PR)
if(SelectCEOSModel eq 3,
  bEOS(AllComp) = 0.07780*8.314E-2*Tc(AllComp)/Pc(AllComp);
  uEOS = 2;
  wEOS = -1;
  omegaA = 0.45724;
  fw(AllComp) = 0.37464 + 1.54226*omega(AllComp) - 0.26992*Power(omega(AllComp),2);
);

* Calculate remaining parameters
Inter0 = sqrt(uEOS*uEOS - 4*wEOS);
