*********************************
***** GenericParameters.gms *****
*********************************

* This file contains generic flowsheet parameters. Do not edit the names
* of scalars or parameters, as they are used in the generic initialization
* procedure.

Scalars
Tmax             Max temperature (K) for the flowsheet                   /320/
* Note: Tmin is part of the initialization procedure

Pmax             Max pressure (bar) for flowsheet                        /40/
Pmin             Min pressure (bar) for flowsheet                        /1.01/

* for zero flow in cascade
epsi             small number for zero situation                         /1E-4/


LBStgEd1         Lower bound on stages in cascade                       /2.0/
UBStgEd1         Upper bound on stages in cascade                       /50.0/

UBFlowEd1        Upper bound on flow in cascade                         /50/
LBFlowEd1        Lower bound on flow in cascade                         /0.0001/

* flash
UBFlowFlsh       Upper bound on outlet flow to flsh                      /50/
UBInFlsh         Upper bound on inlet feed to flsh                       /50/
LBInFlsh         Lower bound on inlet feed to flsh                       /0.0001/

* flow through side draw
LBFlowVSD        Lower bound on flow of side-draw                        /0.0001/
UBFlowVSD        Upper bound on flow of side-draw                        /50/

* for single input mizers
BMTStr           Big M for T in single input mxr                         /130/
BMPStr           Big M for P in single input mxr                         /20/
UBFStr           Upper bound on flow through single input mxr or sptr    /50/

MaxDT            Maximium temperature difference for heat integration    /30/

ComplPen         Complementarity penalty (initial value)                 /100/

epsiSmax         Small epsilon used in smoothed max (initial value)      /1E-4/

StartCPPower     ComplPen = 10^(StartCPPower + 0  1  or  2) in loops     /1/

* Debug mode settings
DebugMode        Amount of debug info to display... 0 = none ... 3 = all /3/
* Mode 0: No debug info is displayed
* Mode 1: Only sets are displayed in .lst file (recommended)
* Mode 2: Add some intermediate variable values
* Mode 3: Add some initialization parameters (such as iterations, etc)
;
