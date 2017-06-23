$Title Optimal Synthesis of Sorbent Carbon Capture Processes

$ONTEXT
================================================================================
Base version: Murthy Konda, NETL (May 19, 2011)
Modifications:
 - YoungJung Chang, CMU (August, 2011)
 - Nick Sahinidis, CMU (September, 2011; December, 2011)
 - Murthy Konda, NETL (December, 2011)
 - Zhihong Yuan, CMU (March, 2012 - August, 2014)

This model identifies an optimal process configuration and corresponding optimal technology and optimal
design/operation levels for bubbling fluidized bed based CO2 capture process

ADDITIONAL NOTES -
* flue gas from FGD is saturated at 1 atm
* 650 MWe (Gross) Pulverized Coal-fired power plant
    P:    14.689 psia
    T:    149.82 F
    H*:   18.47 Btu/lb
    Mgas: 1673 lb/s ~ 26.14336 kgmol/s
    Mash: 0.006 lb/s
    V:    25654.3 ft^3/s
    Mole Composition %
      N2   70.123
      O2   4.547
      CO2  12.272
      H2O  12.21
      Ar   0.843
      SO2  0.005
* Minimize Total Annualized Cost or Levelized Cost of Electricity

MAJOR ASSUMPTIONS
*   Each stage is a single stage operation (stacked BFB will be considered later)
*   The steam mixed with the CO2 feed is always at the higher pressure
*   No pressure change for liquid and solid flow
*   Very simple linear heat transfer for sorbent HX
*   Utility cost for sorbent HX is negligible
*   LMTD based heat transfer for flueHX

STREAMS NOMENCLATURE
*   flueIn, flueOut: flue gas around the flueHX
*   utilIn, utilOut: coolant water for flueHX
*   gasIn, gasOut: process gas around Separator
*   solidOut: sorbent output from Separator
*   steam: steam input to Desorption Mixer
*   pureCO2: compressed CO2 input to Desorption Mixer
*   feedCO2: CO2 input to Desorption Compressor
*   coldIn, coldOut: coolant water for Adsorber
*   hotIn, hotOut: heating steam for Desorber
*   solidLean: sorbent output from leanHX
*   solidRich: sorbent output from richHX
*   coolIn, coolOut: coolant water for leanHX
*   warmIn, warmOut: heating steam for richHX

BALANCE EQUATIONS SUFFIXING SCHEME
*   flueHX:      0g,     0u   (only in the adsorption train)
*   separator:   1g, 1s, 1u
*   blower:      2g
*   sorbentHX:   3s, 3u
*   mixer:       4g           (only in the regeneration train)

DESIGN VARIABLES
*   y(s)
*   sorbentF
*   steamF, feedCO2F
*   coldInF(a)
*   hotInF(d)
*   unitD(s), unitL(s)
*   unitHXArea(s) or {unitDx(s) and unitLx(s)}
*   deltaPa
*   deltaPd
*   richHXArea, leanHXArea
*   utilInF

================================================================================
$OFFTEXT


* Maximum number of serial separators in a train
$SET  noADS  4
$SET  noDES  4


$EOLCOM  //
$OFFSYMLIST
$OFFSYMXREF
*$OFFLISTING

SETS
  s      Separator stages                / a1*a%noADS%, d1*d%noDES% /
  a(s)   Stages in adsorption section    / a1*a%noADS% /
  d(s)   Stages in regeneration section  / d1*d%noDES% /
  t      Process technologies            / BOF, BUF /
  c      Components superset             / CO2, N2, H2O, HCO3, NH2COO /
  fc(c)  Components in fluid             / CO2, N2, H2O /
  sc(c)  Components in solid             / H2O, HCO3, NH2COO /
  v      Thermodynamic properties        / P, T /;
ALIAS (fc,fcp);
ALIAS (a, aa);
ALIAS (d, dd);

SCALAR  modelInput  To be substituted by surrogate models  / 0 /;

SCALARS
  pi             Ratio of circle's circumference to diameter    / 3.14159265 /
  R              Universal ideal gas constant (J per mol*K)     / 8.314472 /
* Individual gas constants: H2O(8.3162), CO2(8.3135), N2(8.3163), Air(8.3115)
  Rp             Conversion of R to (psi*ft^3 per kgmol*K)
  fgF            Total flow rate of flue gas (kgmol per sec)    / 27.882089/



  CRF            Capital Recovery Factor                        / 0.124 /
$ONTEXT
CRF determines the annual expense due to purchase of equipment including debt
and interest payments.  It may also be affected by the tax consequences of
property taxes and depreciation allowances.  According to NETL's Cost Estimation
Quality Guidelines DOE/NETL-2010/1455, the capital charge factor for a high
risk, 5-year project is 0.124.
$OFFTEXT
  OCF            Overnight Cost Factor                          / 5.82 /
$ONTEXT
OCF converts from FOB cost to total capital outlay.  It is the product of the
Lang Factor (5.04 for fluids handling), the CE Index (550/500), and the Delivery
Factor (1.05).
$OFFTEXT
  maintOverhead  TDCF*(% Capital on Maintenacne)*O&M Overhead   / 0.0828 /
$ONTEXT
TDCF(Total Direct Cost Factor) is applied to the Overnight Cost to "back-out"
the total direct cost of equipment for the purpose of evaluating O&M costs.
From SSL Table 22.17: 0.6
$OFFTEXT
  opWage         Wages & benefits to operators ($ per hour)     / 40 /
  opOverhead     Overhead on Operator Labor*O&M Overhead        / 1.5 /
  rho            Carbon Steel Density (lbs per ft^3)            / 490 /
  fm             Carbon Steel Fabrication Factor                / 1 /
  CF             Capacity Factor                                / 0.85 /
  k              Ratio of specific heats of diatomic gas        / 1.4 /
* Heat capacity ratio (Cp/Cv): H2O(1.33), CO2(1.28), N2&Air(1.40)
*                              monatomic(5/3), linear(7/5), nonlinear(8/6)
  nameplate      Nameplate Capacity before capture (MW)         / 650 /
  plateSpace     Space between MB reactor plates (ft)           / 0.656/
  elevWidth      Elevator Width (inches)                        / 12 /
  unitT          Wall Thickness (ft) [1 in - 1.25 in]           / 0.104 /
  unitEta        Blower Efficiency                              / 0.85 /
  vtr           Bound on maximum superficial gas velocity       / 4 /
  solidCp        Heat capacity of sorbent (kJ per K*kg)         / 1.13 /
  solidU         Heat transfer coeff. for solid (kW per K*m^2)  / 0.3 /
  cndU           Condenser heat transfer coeff. (kW per K*m^2)  / 0.1 /;
* cndU may go up to 1 kW/m^2/K
Rp = 14.696/101.325*100/2.8317*R;

PARAMETERS
  cP(fc)        Specific heat capacity (kJ per kgmol.K)    / CO2 = 41.3, N2 = 29.2, H2O = 34.2 /
* Isobar specific heat capacity data (kJ per kg.K)
*  T(K)   H2O    CO2     N2    Air
*  375  1.890  0.918  1.042  1.009 (100 C)
*  400  1.901  0.939  1.044  1.013 (120-140 C)
*  450  1.926  0.978  1.049  1.017 (160 C) ~ 1.022 (180 C)
  fgC(fc)       Mole fractions of flue gas                 / CO2 = 0.12, N2 = 0.74, H2O = 0.14 /
  fgV(v)        P and T of flue gas from power plant       / P = 1.01, T = 54.18 /
  flueInV(v)    P and T of flue gas treated by one train
  steamC(fc)    Mole fractions of steam injection          / CO2 = 0, H2O = 1 /
  steamV(v)     P and T of the steam to regenerators       / P = 6.8, T = 180/
  feedCO2C(fc)  Mole fractions of feedCO2                  / CO2 = 0.42, H2O = 0.58 /
  feedCO2V(v)   P and T of the CO2 stream to regenerators  / P = 1, T = 138.07 /
  utilInC(fc)   Mole fractions of coolant water to flueHX  / CO2 = 0, N2 = 0, H2O = 1 /
  utilInT       T of the coolant water to flueHX           / 32.22 /
  coldInC(fc)   Mole fractions of coolant water
  coldInT       T of the coolant water to adsorber
  hotInC(fc)    Mole fractions of heating steam
  hotInT        T of the heating steam to desorber
  coolInT       T of the solid cooling water
  warmInT       T of the solid heating steam
  unitNp(s)                                                /a1=1,a2=1,a3=1,a4=1,d1=1,d2=1,d3=1,d4=1/



;

flueInV(v) = fgV(v);
coldInC(fc) = utilInC(fc);       coldInT = utilInT;
hotInC(fc) = steamC(fc);         hotInT = steamV('T');
coolInT = utilInT;
warmInT = steamV('T');


positive variable  COE  Adder to the cost of electricity due to capture ($ per MWh);
positive variable  soi  sum of constraint infeasibilities;
variable fplus objective;


BINARY VARIABLES
  x(s,t)                 1 iff technology t is used for stage s
  y(s)                   1 iff stage s is ever utilized;


INTEGER VARIABLES
  Nu                     Number of parallel trains

POSITIVE VARIABLES
  capEX                  Capital overnight cost
  unitCpa(a)             FOB cost of Adsorber stage s ($)
  unitCpd(d)             FOB cost of Regenerate stage s ($)
  richHXArea             Area of the rich sorbent heater (ft^2)
  leanHXArea             Area of the lean sorbent cooler (ft^2)
  flueHXArea             Area of the flue gas condenser (ft^2)
  derate                 Derating of the plant due to steam takeoff (MW)
  CaptureTarget          CO2 fraction to be captured from flue gas
 
* FLOW RATES (default unit: kgmol/sec)
  sorbentF               Sorbent flow rate (kg per hr)
  flueInC(fc), flueOutC(fc), flueOut
  utilInF
  gasInC(s,fc), gasOutC(s,fc), gasOut(s)
  solidOutC(s,sc)
  steamF
  feedCO2F
  coldInF(a)
  hotInF(d)
* PRESSURE/TEMPERATURE Variables (Pressure in atm, Temperature in C)
  flueOutV(v)
  utilOutT
  gasInV(s,v), gasOutV(s,v)
  solidOutT(s)
  pureCO2V(v)
  coldOutT(a)
  hotOutT(d)
  deltaPa              The pressure increase through the compressor in adsorber
  deltaPd              The pressure increase through the compressor in regenerator
  solidLeanT
  solidRichT
* PROCESS VARIABLES
  unitD(s)               Diameter (ft)
  unitL(s)               Length (ft)
  unitPi(s)              Inlet Absolute Pressure (psi)
  unitPo(s)              Outlet Absolute Pressure (psi)
  unitF(s)               Input Gas Flow Rate to Unit (ft^3 per min)
  unitFb(s)              Input Gas Flow Rate to Blower (ft^3 per min)
  unitW(s)               Vessel weight (lbs)
  unitPb(s)              Blower Power (hp)
  unitPe(s)              Elevator Power (hp)
  unitHXArea(s)          Area of HXs for Separators (ft^2)
  unitDx(s)              HX tube diameter (m)
  unitNx(s)              Number of HX tubes
  unitDLX(s)             HX tube spacing (m)
  steamFlow              Steam Take-off amount from power plant (kg per sec)
  gasIn(s)               Molar flow rate of the gas inlet to stage s
  gasInX(s,fc)           Fractional molar composition of fc in the inlet gas stream to stage s
  gasOutX(s,fc)          Fractional molar composition of fc in the outlet gas stream to stage s
  solidOutC(s,sc)        Solid component Loading of the solid stream out of stage s

* SLACK VARIABLES
  P(*,s), N(*,s), P1(*), N1(*), P2(a,fc), N2(a,fc), P3(s,fc), N3(s,fc), P4(s,fc), N4(s,fc), P5(*), N5(*), P6(*),N6(*)
;

* Bounds



set i /1*19/;  //flags to turn on/off model equations
parameters ads(i), des(i),econ(i); ads(i)=0; des(i)=0; econ(i)=0;
set kinf /BFBads1*BFBads14, BFBDes1*BFBDes14, extra1/;
parameter scale(kinf)   Scale for equation ;

EQUATIONS
  objective, SumOfInfeasibilities, combo
  eqDerate, eqSteam, eqNu(fc)
  eqCapEX
  eqUnitCpa(s),eqUnitCpd(s), eqUnitW(s), eqUnitPb(s), eqUnitPe(s),
  eqUnitFbAds(a), eqUnitFbDes(d), eqUnitFAds(a), eqUnitFDes(d)


  eqUnitPiAds(a), eqUnitPoAds(a), eqUnitPiDes(d), eqUnitPoDes(d)
  eqUnitHXArea(s), eqFlueHXArea
  eqEnrgAds0g(v)
  eqStageAds(a), eqStageDes(d)
  eqstage(s)
  eqNumAds(a), eqNumDes(d)

  eqMassAds0g(fc), eqEnrgAds0u, eqEnrgAds2g(a,v), eqEnrgAds3s
  eqEnrgDes2g(v), eqEnrgDes3s, eqEnrgDes4g(d,v)

  eqCaptureTarget
  eqRegenerationbalance
  eqAdsorbebalance

  cont1
  cont2


  extra1  // Depth of the bed should be no larger than the diameter
  BFBads1
  BFBads2
  BFBads3
  BFBads4
  BFBads6
  BFBads7
  BFBads8
  BFBads9
  BFBads10
  BFBads11
  BFBads12
  BFBDes1
  BFBDes2
  BFBDes3
  BFBDes4
  BFBDes6
  BFBDes7
  BFBDes8
  BFBDes9
  BFBDes10
  BFBDes11
  gasOutConvads(a,fc)
  gasOutConvdes(d,fc)
  gasInadConv(a,fc)
  gasIndeConv(d,fc)
  gasInCS(s,fc)
  gasOutSum(a)
  gasInSum(a)
  gasOutSumd(d)
  gasInSumd(d)
  flueOutSum
;

*ECONOMIC MODULE

objective..

  CF*365*24*(nameplate/7446-derate/7446)*(COE-7.53) =G= (0.124*(1316.19E6+capEX)+(1.25*(6.80E6+0.0053343*(1316.19E6+capEX))+0.0163*(1316.19E6+capEX))+0.008*(1316.19E6+capEX)+14.41E6+70.08E6)/7446;

SumOfInfeasibilities ..
  soi =G= sum((kinf, a), (P(kinf,a)+N(kinf,a))/scale(kinf))+sum((kinf, d), (P(kinf,d)+N(kinf,d))/scale(kinf))
        + P1('flueOut') + N1('flueOut')
        + sum((a,fc), P2(a,fc)+N2(a,fc))
        + sum((s,fc), P3(s,fc)+P4(s,fc)+N3(s,fc)+N4(s,fc));

combo .. fplus =e= 0.2*COE+ 1*soi;

eqDerate $econ('1') ..
  100*derate =G= 100*0.924*steamFlow //420.42*steamFlow;
  + 100*7.46E-4*SUM(s,Nu*y(s)*unitPb(s)+Nu*y(s)*unitPe(s));

eqSteam $econ('2') ..
  steamFlow =G= 18*Nu*(steamF+SUM(d,y(d)*hotInF(d)));

eqNu(fc) $econ('3') ..
  flueInC(fc) =E= fgF*fgC(fc)/Nu;

eqCapEX $econ('4') ..
   (capEX/(OCF))/1e5 =G= SUM(a,Nu*y(a)*unitCpa(a)/1e5) + SUM(d,Nu*y(d)*unitCpd(d)/1e5)
  + 1.10 *Nu* (EXP(11.9052-0.8709*LOG(richHXArea)+0.09005*LOG(richHXArea)**2-LOG(1e5)))
  + 1.10 *Nu* (EXP(11.9052-0.8709*LOG(leanHXArea)+0.09005*LOG(leanHXArea)**2-LOG(1e5)))
  //+ 1.03*1.12*4.25*EXP(11.9052-0.8709*LOG(flueHXArea)+0.09005*LOG(flueHXArea)*LOG(flueHXArea));
  + 1.12*0.98*Nu*((7*EXP(11.9052-0.8709*LOG(12000)+0.09005*LOG(12000)*LOG(12000)-LOG(1e5))+EXP(11.9052-0.8709*LOG(flueHXArea-7*12000)+0.09005*LOG(flueHXArea-7*12000)*LOG(flueHXArea-7*12000)-LOG(1e5))));

eqUnitCpa(a) $econ('5') ..
  unitCpa(a)/1e5 =G= fm*EXP(7.0132+0.18255*LOG(unitW(a))+0.02297*LOG(unitW(a))**2-LOG(1e5))   // Cost of vessels
    + EXP(7.59176+0.7932*LOG(unitPb(a))-0.0129*LOG(unitPb(a))**2-LOG(1e5))$(ORD(a) eq 1)+0* EXP(7.59176+0.7932*LOG(unitPb(a))-0.0129*LOG(unitPb(a))**2-LOG(1e5))$(ORD(a) ne 1)                // Blower (compressor) cost
    + 1.10 * (EXP(11.9052-0.8709*LOG(unitHXArea(a))+0.09005*LOG(unitHXArea(a))**2-LOG(1e5)))      // Cost of in-reactor heat exchanger
    + unitNp(a)*734.96*1.87*(0.3048*unitD(a))**2.0049/1e5                                // Cost of plates =?= unitNp(s,t)*468.0*EXP(0.1739*unitD(s))
    + 361.8*(unitD(a)**0.7396)*(unitL(a)**0.70684)/1e5 //300.9*(unitD(s)**0.63316)*(unitL(s)**0.80161)                                                 // Cost of platforms and ladders
    + 610*SQRT(elevWidth)*unitL(a)**0.57/1e5                                                           // Cost of elevator
    + 1.3*EXP(5.8259+0.13141*LOG(unitPe(a))+0.053255*LOG(unitPe(a))**2+0.028628*LOG(unitPe(a))**3-0.0035549*LOG(unitPe(a))**4-LOG(1e5));  // Cost of elevator motor

eqUnitCpd(d) $econ('6') ..
  unitCpd(d)/1e5 =G= fm*EXP(7.0132+0.18255*LOG(unitW(d))+0.02297*LOG(unitW(d))**2-LOG(1e5))   // Cost of vessels
    + EXP(7.5800+0.8*LOG(unitPb(d))-LOG(1e5))$(ORD(d) eq 1)+ 0*EXP(6.8929+0.7900*LOG(unitPb(d))-LOG(1e5))$(ORD(d) ne 1)                                       // Blower (compressor) cost
    + 1.10 * (EXP(11.9052-0.8709*LOG(unitHXArea(d))+0.09005*LOG(unitHXArea(d))**2-LOG(1e5)))      // Cost of in-reactor heat exchanger
    + unitNp(d)*734.96*1.87*(0.3048*unitD(d))**2.0049/1e5                             // Cost of plates =?= unitNp(s,t)*468.0*EXP(0.1739*unitD(s))
    + 361.8*(unitD(d)**0.7396)*(unitL(d)**0.70684)/1e5 //300.9*(unitD(d)**0.63316)*(unitL(d)**0.80161)                                                 // Cost of platforms and ladders
    + 610*SQRT(elevWidth)*unitL(d)**0.57/1e5                                                            // Cost of elevator
    + 1.3*EXP(5.8259+0.13141*LOG(unitPe(d))+0.053255*LOG(unitPe(d))**2+0.028628*LOG(unitPe(d))**3-0.0035549*LOG(unitPe(d))**4-LOG(1e5));  // Cost of elevator motor

eqUnitW(s) $econ('7') ..
 //  unitW(s) =G= pi*(unitD(s)+unitT)*(unitL(s)+5+0.8*unitD(s))*rho*unitT;
     unitW(s)/500 =G= pi*(unitD(s)+unitT)*(unitL(s)+0.8*unitD(s))*rho*unitT/500;

eqUnitPb(s) $econ('8') ..
     unitPb(s) =G= 0.00436*(k/(k-1))*(unitFb(s)*unitPi(s)/unitEta)*((unitPo(s)/unitPi(s))**(1-1/k)-1);

eqUnitPe(s) $econ('9') ..  // basic power to drive buckets + power required to lift sorbent
     3600*unitPe(s) =G= 0.020*2.2046*sorbentF*unitL(s)**0.63 + 0.00182*2.2046*sorbentF*unitL(s);

* LINK between ECONOMICS MODULE and PROCESS MODULE

eqUnitFbAds(a) $econ('10') ..
     unitPo(a)*unitFb(a) =E= Rp*60*SUM(fc,gasOutC(a-1,fc)+flueOutC(fc)$(ORD(a)=1))*(gasInV(a,'T')+273.15);

eqUnitFbDes(d) $econ('11') ..
     unitPo(d)*unitFb(d) =E= Rp*60*SUM(fc,feedCO2F*feedCO2C(fc)$(ORD(d)=1)+gasOutC(d-1,fc))*((pureCO2V('T')$(ORD(d) eq 1)+gasInV(d,'T')$(ORD(d) ne 1))+273.15);

eqUnitFAds(a) $econ('12') ..
     unitPo(a)*unitF(a) =E= Rp*60*SUM(fc,gasOutC(a-1,fc)+flueOutC(fc)$(ORD(a)=1))*(gasInV(a,'T')+273.15);

eqUnitFDes(d) $econ('13') ..
     unitPo(d)*unitF(d) =E= Rp*60*SUM(fc,gasInC(d,fc))*(gasInV(d,'T')+273.15);

eqUnitPiAds(a) $econ('14') ..
  unitPi(a) =E= 14.696*(gasOutV(a-1,'P')+flueOutV('P')$(ORD(a)=1));

eqUnitPiDes(d) $econ('15') ..
  unitPi(d) =E= 14.696*feedCO2V('P')$(ORD(d) eq 1)+unitPo(d)$(ORD(d) ne 1);

eqUnitPoAds(a) $econ('16') ..
  unitPo(a) =E= 14.696*gasInV(a,'P');

eqUnitPoDes(d) $econ('17') ..
  unitPo(d) =E= 14.696*(pureCO2V('P')$(ORD(d) eq 1)+gasInV(d,'P')$(ORD(d) ne 1));

eqUnitHXArea(s) $econ('18') ..
  unitHXArea(s)/9000 =G= unitNx(s)*3.28084*unitDx(s)*pi*unitL(s)/9000;

eqFlueHXArea $econ('19') ..
  cndU*0.0929*flueHXarea*((flueOutV('T')-utilInT)-(flueInV('T')-utilOutT))/LOG((flueOutV('T')-utilInT)/(flueInV('T')-utilOutT)) =E= cP('H2O')*utilInF*(utilOutT-utilInT);



* LOGICAL ASSIGNMENT CONSTRAINTS
eqStageAds(a)$(ORD(a)<CARD(a))..  y(a) =G= y(a+1);
eqStageDes(d)$(ORD(d)<CARD(d))..  y(d) =G= y(d+1);
eqNumAds(a).. SUM(aa,y(aa)) =G= 1;
eqNumDes(d).. SUM(dd,y(dd)) =G= 1;
eqstage(s).. SUM(t,x(s,t)) =E= y(s);


* MATERIAL BALANCES
eqMassAds0g(fc)..  // gasC change around flueHX for H2O vs noncondensible

  (1$SAMEAS(fc,'H2O')+1$(NOT SAMEAS(fc,'H2O')))*flueOutC(fc) =E= flueInC(fc)$(NOT SAMEAS(fc,'H2O'))+0.58*flueInC(fc)$SAMEAS(fc,'H2O');


* HYDRODYNAMIC/ENERGY BALANCES

eqEnrgAds0g(v)..  // flash calculation
   flueOutV(v) =E= flueInV(v)$SAMEAS(v,'P')+(5.9405*rpower(flueInV(v),1/2))$SAMEAS(v,'T');

eqEnrgAds0u..  // HX model assuming counter-current operation
   utilOutT*utilInF =E=
        utilInT*utilInF
        + (40683*((374-flueOutV('T'))/274)**0.38*(flueInC('H2O')-flueOutC('H2O'))
           -SUM(fc,cP(fc)*flueInC(fc)*(flueOutV('T')-flueInV('T')))
          )/cP('H2O');

eqEnrgAds2g(a,v)..  // gasV change around adsorber

   gasInV(a,v) =E= ((((flueOutV(v)+273.15)*(gasInV(a,'P')/(gasInV(a,'P')-deltaPa))**(1-1/k)-273.15)$(ORD(a)=1)+gasOutV(a-1,v)$(ORD(a) ne 1))$SAMEAS(v,'T')+((flueOutV(v)+deltaPa)$(ORD(a)=1)+gasOutV(a-1,v)$(ORD(a) ne 1))$SAMEAS(v,'P'));

eqEnrgAds3s..  // solidT change around leanHX
  solidCp*sorbentF*(solidLeanT-solidOutT('d1'))/3600
        =E= solidU*0.0929*leanHXArea*(coolInT-0.5*solidLeanT-0.5*solidOutT('d1'));

eqEnrgDes2g(v)..  // gasV change around Blower(d) assuming adiabatic compression
  pureCO2V(v) =E= ((feedCO2V(v)+273.15)*(pureCO2V('P')/feedCO2V('P'))**(1-1/k)-273.15)$SAMEAS(v,'T')+(feedCO2V(v)+deltaPd)$SAMEAS(v,'P');

eqEnrgDes3s..  // solidT change around richHX (to use LMTD, include equation for warmOutT)
  solidCp*sorbentF*(solidRichT-solidOutT('a1'))/3600 =E= solidU*0.0929*richHXArea*(warmInT-0.5*solidRichT-0.5*solidOutT('a1'));

eqEnrgDes4g(d,v)..  // gasV change around regenerator
  gasInV(d,v) =E= (pureCO2V(v)$SAMEAS(v,'P') + (feedCO2F*pureCO2V('T')+steamF*steamV('T'))$SAMEAS(v,'T')/(1$SAMEAS(v,'P')+(feedCO2F+steamF)$SAMEAS(v,'T')))$(ORD(d) eq 1) + (gasOutV(d-1,v)$SAMEAS(v,'P')+gasOutV(d-1,v)$SAMEAS(v,'T'))$(ORD(d) ne 1);

* CAPTURE TARGET must be met

eqCaptureTarget..
  flueInC('CO2') - gasOutC('a%noADS%','CO2') =G= CaptureTarget*flueInC('CO2');

eqRegenerationbalance..
  1000*flueInC('CO2')+solidOutC('d1','NH2COO')*sorbentF+solidOutC('d1','HCO3')*sorbentF =E= 1000*gasOutC('a%noADS%','CO2')+solidOutC('a1','NH2COO')*sorbentF+solidOutC('a1','HCO3')*sorbentF;

eqAdsorbebalance..
  1000*feedCO2F*feedCO2C('CO2')+solidOutC('a1','NH2COO')*sorbentF+solidOutC('a1','HCO3')*sorbentF =E= 1000*gasOutC('d%noADS%','CO2')+solidOutC('d1','NH2COO')*sorbentF+solidOutC('d1','HCO3')*sorbentF ;

* Additional process related conditions must be also added
* For example, the input variable bounds for each technology and its model

gasOutConvads(a,fc)..   gasOutC(a,fc) =E= gasOutX(a,fc)*gasOut(a);
gasOutConvdes(d,fc)..   gasOutC(d,fc) =E= gasOutX(d,fc)*gasOut(d);

gasInadConv(a,fc)..    gasInC(a,fc) =E= gasOutC(a-1,fc)$(ORD(a) ne 1)+flueOutC(fc)$(ORD(a) eq 1);
gasIndeConv(d,fc)..    gasInC(d,fc) =E= (feedCO2F*feedCO2C(fc) + steamF*steamC(fc))$(ORD(d) eq 1)+gasOutX(d-1,fc)*gasIn(d)$(ORD(d) ne 1);
                 
gasOutSum(a)..         SUM(fc,gasOutX(a,fc)) =E= 1;
gasInSum(a)..          SUM(fc,gasInX(a,fc)) =E= 1;
gasInSumd(d)..         gasInX(d,'CO2')+gasInX(d,'H2O') =E= 1;
gasOutSumd(d)..        gasOutX(d,'CO2')+gasOutX(d,'H2O') =E= 1;


flueOutSum..           flueOut =E= Sum(fc,flueOutC(fc))+ P1('flueOut') - N1('flueOut');

cont1(a) ..            gasIn(a) =E= gasOut(a-1) $(ORD(a) ne 1) + flueOut $(ORD(a) eq 1);
cont2(d) ..            gasIn(d) =E= (feedCO2F+steamF)$(ORD(d) eq 1) + gasOut(d-1)$(ORD(d) ne 1);

gasInCS(s,fc)..        100*gasInC(s,fc) =E= 100*gasInX(s,fc)*gasIn(s) + P4(s,fc)-N4(s,fc);

extra1(a) ..           unitL(a) =L= unitD(a);

* BFB adsorber surrogate model

BFBads1(a) $ads('1') .. unitNx(a) =E=

      (7074.4424 * EXP(unitDx(a)) + 0.4017 * unitDLX(a)**(-3.0) - 0.04366 * unitDLX(a)**(-4.0) - 2.9778908 * (unitDx(a)*unitDLX(a))**(-1.0) + 0.004788605 * (unitDx(a)*unitDLX(a))**(-2.0) - 7493.5386083 * (unitDx(a)*unitDLX(a))**(1/2) + 373215.0257197 * (unitDx(a)*unitDLX(a))**2 - 1439.07841395 * ((unitD(a)/3.28084)/unitDLX(a)/18)**(-1.0) - 7646.46344045 * ((unitD(a)/3.28084)/unitDLX(a)/18)**(1/2) + 3022.081222 * ((unitD(a)/3.28084)/unitDLX(a)/18) + 52.712785 * ((unitD(a)/3.28084)/unitDLX(a)/18.)**2.0 - 25704.67341 * (unitDx(a)/unitDLX(a))**2.0)*x(a,'BOF')+(7036.4564618 * EXP(unitDx(a))+ 0.41037566762 * unitDLX(a)**(-3.0) - 0.44133E-001 * unitDLX(a)**(-4.0) -2.970380967 * (unitDx(a)*unitDLX(a))**(-1.0) + 0.4790502 * (unitDx(a)*unitDLX(a))**(-2.0) - 7437.8045041 * (unitDx(a)*unitDLX(a))**(1/2) +372823.88653688 * (unitDx(a)*unitDLX(a))**2.0 - 1429.896542801 * ((unitD(a)/3.28084)/unitDLX(a)/18.)**(-1.0) - 7612.463174627 * ((unitD(a)/3.28084)/unitDLX(a)/18)**(1/2) + 3011.680210 * ((unitD(a)/3.28084)/unitDLX(a)/18) + 52.92446264 * ((unitD(a)/3.28084)/unitDLX(a)/18)**2.0 - 25692.63828924 * (unitDx(a)/unitDLX(a))**2.0)*x(a,'BUF')
        +  P('BFBads1',a)  -  N('BFBads1',a)    ;

BFBads2(a) $ads('2') .. gasInV(a,'P') =E=
        (0.5217168512 * EXP(gasOutV(a,'P')/1.4) + 0.1766880 * ((unitL(a)/3.28084)*gasOutV(a,'P')/11)+ P('BFBads2',a)  -  N('BFBads2',a))*x(a,'BOF')+(0.521638756 * EXP(gasOutV(a,'P')/1.4) + 0.177059343 * ((unitL(a)/3.28084)*gasOutV(a,'P')/11) +  P('BFBads2',a)  -  N('BFBads2',a) )*x(a,'BUF') + gasOutV(a,'P')*(1-y(a))    ;

BFBads3(a) $ads('3') .. 3600*gasOut(a) =E=
       (8491.7202136 * gasIn(a)*3600/0.90E+04 + 0.462851E-001 * unitDLX(a)**(-3.0) + 35.1941308 * ((unitL(a)/3.28084)*unitDx(a)/8.0)**(1/2) - 0.20365412E-002 * ((unitL(a)/3.28084)*unitDx(a)/8.0)**(-2.0) + 0.9440E-001 * (unitDLX(a)*sorbentF/0.90E+06)**(-2.0) + 446.036301 * sorbentF*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a))) + 1040.71614 * (unitDLX(a)*gasIn(a)*3600/0.90E+04)**2.0 + 33572.3618004 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a))))**2.0 + 7584.4753563 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**2.0 + 21.97001280 * ((unitD(a)/3.28084)/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/15)**(-1.0) - 207.739495 * (unitDLX(a)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(-1.0) - 650.793883 * (unitDLX(a)/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(-1.0) - 25.9726522 * ((unitD(a)/3.28084)/unitDLX(a)/18.)**(-2.0) + 84.0033494 * ((unitD(a)/3.28084)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.16)**(-2.0) + 30.749462 * ((unitD(a)/3.28084)/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/15)**(-2.0) + 130.20905 * (unitDLX(a)/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(-2.0) - 19.6528809 * (gasInV(a,'T')/sorbentF/0.11E-03)**(-2.0) - 165.35111 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/0.83)**0.50 - 1355.55513 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/gasOutV(a,'P')/0.71) - 512.98702 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02) - 141.691318 * (unitDLX(a)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**2.0 +  P('BFBads3',a)  -  N('BFBads3',a))*x(a,'BOF')+(9004.315343 * gasIn(a)*3600/0.90E+04 + 20.8513750 * ((unitL(a)/3.28084)*unitDx(a)/8.0)**0.50 - 5.02669976310 * ((unitL(a)/3.28084)*unitDx(a)/8.0)**(-1.0) - 10.3808750 * (unitDLX(a)*sorbentF/0.90E+06)**(-1.0) + 0.6538200 * (unitDLX(a)*gasIn(a)*3600/0.90E+04)**(-2.0) - 0.48166E-003 * (unitDLX(a)*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(-2.0) + 0.309281771 * (unitDLX(a)*sorbentF/0.90E+06)**(-2.0) + 0.3320919 * (unitDLX(a)*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03)**(-2.0) - 0.11599926E-001 * (gasIn(a)*3600*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/0.90E+04)**(-2.0) + 6.5313460 * (sorbentF*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.99E+08)**(-2.0) - 296.32391 * ((unitD(a)/3.28084)*gasIn(a)*3600/0.16E+06) + 326.6966211 * (unitDLX(a)*sorbentF/0.90E+06) + 2547.41914 * (gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))  - 5752.706080 * ((unitD(a)/3.28084)*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/18.)**2.0 - 93774.50375 * ((unitL(a)/3.28084)*unitDx(a)/8.0)**2.0 + 1091.29369079 * (unitDLX(a)*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03)**2.0 + 366.9711630 * ((solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))/0.11E+03)**2.0 - 460.27407660 * (gasInV(a,'T')/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/0.10E+03)**(1/2) - 610.25247 * (unitDLX(a)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(-1.0) - 377.9922391 * (unitDLX(a)/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(-1.0) + 9.9873412 * (gasInV(a,'T')/sorbentF/0.11E-03)**(-1.0)  - 240.45051 * ((unitD(a)/3.28084)/unitDLX(a)/18)**(-2.0) + 34.84502 * ((unitD(a)/3.28084)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.16)**(-2.0) +  22.049716*((unitD(a)/3.28084)/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/15)**(-2.0)+ 53.28543* ((unitD(a)/3.28084)/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/15)**(-2.0) - 52.418372* ((unitL(a)/3.28084)/unitDLX(a)/8.0)**(-2.0) +  436.562863* ((unitL(a)/3.28084)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/8.0)**(-2.0) + 16.691110* ((unitL(a)/3.28084)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.73E-01)**(-2.0) + 155.96485*(unitDLX(a)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(-2.0) - 12.358223 * (unitDLX(a)/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(-2.0) - 0.47748682 * (unitDLX(a)/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/0.83)**(-2.0) + 0.5923E-001 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/sorbentF/0.11E-05)**(-2.0) - 200.996877 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/0.83)**(1/2) - 498.115 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02) - 245.815954 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.83)**2.0  + 296.3867*((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/sorbentF/0.11E-05)**2.0 - 1109.695550 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**2.0 - 57.3707835*(sorbentF/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.82E+04)**2.0 +  P('BFBads3',a)  -  N('BFBads3',a))*x(a,'BUF')+3600*(flueOut$(ORD(a) eq 1)+gasOut(a-1)$(ORD(a) ne 1))*(1-y(a));

BFBads4(a) $ads('4') .. gasOutV(a, 'T') =E=

       (56.899125 * exp(unitDLX(a))- 76.53910 * unitDLX(a)**4.0 - 0.051 * (gasInV(a,'T')/0.10E+03)**(-4.0) - 0.4545029* (unitDLX(a)*gasIn(a)*3600/0.90E+04)**(-1.0) - 0.383082 * ((unitD(a)/3.28084)*gasInV(a,'T')/0.18E+04)**(-2.0) + 0.7301286430E-002 * ((unitD(a)/3.28084)*gasInV(a,'T')/0.10E+03)**(-2.0) - 115.68605 * ((unitL(a)/3.28084)*unitDx(a)/8.0)**(1/2) + 13.1366895 * ((unitL(a)/3.28084)*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.88E+03)**(1/2) + 37.48765819 * (unitDLX(a)*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(1/2) + 9.2633212 * (gasInV(a,'T')*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/0.10E+03)**(1/2)- 18.98608 * unitDLX(a)*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))+40.71591 *(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*gasOutV(a,'P')/1.4- 11.9825 * ((unitD(a)/3.28084)*(unitL(a)/3.28084)/0.14E+03)**2.0 + 24.5376 * ((unitL(a)/3.28084)*unitDLX(a)/8.0)**2.0 + 19.419642 * (unitDLX(a)*gasInV(a,'T')/0.10E+03)**2.0 - 9.47471350 * (unitDLX(a)*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**2.0 - 20.737519 * (unitDLX(a)*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.2)**2.0 - 230.2949157 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**2.0 + 56.9774721 * ((unitD(a)/3.28084)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/18)**0.50 + 2.11012 * ((unitD(a)/3.28084)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.16)**(1/2) - 17.9382567 * (unitDLX(a)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(1/2) - 3.79470 * (sorbentF/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.75E+06)**(1/2) - 29.41709 * ((unitD(a)/3.28084)/unitDLX(a)/18)**(-1.0) + 3.85982 * ((unitL(a)/3.28084)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.73E-01)**(-1.0) + 1099.23568* ((unitD(a)/3.28084)/unitDx(a)/18)**(-2.0) - 4.1319 * ((unitD(a)/3.28084)/unitDLX(a)/18)**(-2.0) + 4.9333* ((unitD(a)/3.28084)/gasInV(a,'T')/0.18)**(-2.0) - 2.168173 * ((unitL(a)/3.28084)/unitDLX(a)/8.0)**(-2.0) + 0.1330E-001 * (unitDLX(a)/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.83)**(-2.0) - 1.6867476 * (gasIn(a)*3600/gasInV(a,'T')/90)**(-2.0)+ 1.63342 * (gasIn(a)*3600/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/82)**(-2.0) + 9.9955 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/0.83)**(1/2) - 0.241057E-001 * ((unitL(a)/3.28084)/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/8.0) - 81.04505940 * unitDx(a)/unitDLX(a) + 79.0135 * (unitDx(a)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02) - 1.716881 * (gasIn(a)*3600/sorbentF/0.10E-01) + 0.7256483 * ((unitL(a)/3.28084)/gasInV(a,'T')/0.80E-01)**2.0 + 91.879495 * (unitDx(a)/unitDLX(a))**2.0 - 9.66646644 * (unitDLX(a)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**2.0 - 15.662110 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/sorbentF/0.11E-05)**2.0 +  P('BFBads4',a)  -  N('BFBads4',a))*x(a,'BOF')+(- 8.317256 * exp(unitDLX(a)) - 2.7186437 * ((solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**3.0 + 5.7411315 * (gasInV(a,'T')/0.10E+03)**4.0  + 1.193360 * (sorbentF/0.90E+06)**(-2.0) - 0.5927090 * (gasIn(a)*3600/0.90E+04)**(-3.0)  - 0.236503 * ((unitL(a)/3.28084)*gasOutV(a,'P')/11)**(-2.0) + 0.338248E-004 * (unitDLX(a)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(-2.0) - 36.9403061 * ((unitL(a)/3.28084)*unitDx(a)/8.0)**(1/2)  + 58.04638 * (unitDLX(a)*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03)**(1/2) + 40.945864 * (gasIn(a)*3600*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/0.90E+04)**(1/2) + 15.346660 * (gasInV(a,'T')*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+05)**(1/2) + 6.0739635 * ((unitL(a)/3.28084)*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/9.6) + 18.57295 * (unitDLX(a)*gasInV(a,'T')/0.10E+03) - 16.176403 * (unitDLX(a)*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.2) + 27.625987 * (gasIn(a)*3600*gasInV(a,'T')/0.90E+06) - 16.670232 * (gasInV(a,'T')*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+05)- 6.12261620 * ((solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/0.13E+03) - 4488.81204 * ((unitD(a)/3.28084)*unitDx(a)/18.)**2.0  + 36.780581 * ((unitD(a)/3.28084)*unitDLX(a)/18.)**2.0 - 19.1880702 * ((unitD(a)/3.28084)*gasInV(a,'T')/0.18E+04)**2.0 - 3.30340282 *((unitL(a)/3.28084)*gasInV(a,'T')/0.80E+03)**2.0 + 1378.45764 *(unitDx(a)*gasIn(a)*3600/0.90E+04)**2.0 - 33.933113 * (unitDLX(a)*gasIn(a)*3600/0.90E+04)**2.0 - 71.422885 * (unitDLX(a)*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a))))**2.0 - 1.921273 * (gasInV(a,'T')*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))/0.10E+03)**2.0 + 3.291738 * (gasInV(a,'T')*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.12E+03)**2.0 + 10.37838 * (unitDx(a)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(1/2) + 30.929636 * (gasIn(a)*3600/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/82)**(1/2)  + 2.78630* ((unitL(a)/3.28084)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.73E-01)**(-1.0) - 2.0276802 * (sorbentF/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.75E+06)**(-1.0) - 4.377932 * ((unitL(a)/3.28084)/unitDLX(a)/8.0)**(-2.0) - 0.1707645 * ((unitL(a)/3.28084)/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/6.7)**(-2.0) - 0.3550E-002 * (unitDx(a)/unitDLX(a))**(-2.0)  - 0.488416E-001 * (unitDx(a)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(-2.0) - 0.2055E-001 * (unitDLX(a)/gasInV(a,'T')/0.10E-01)**(-2.0) + 0.4019E-001 * (unitDLX(a)/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.83)**(-2.0) - 1.7032874 * (sorbentF/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.82E+04)**(-2.0) + 4.2564 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/0.83)**(1/2) + 6.91731 * ((unitD(a)/3.28084)/(unitL(a)/3.28084)/2.3) - 4.6106136 * ((unitD(a)/3.28084)/unitDLX(a)/18) - 3.1824 * (gasIn(a)*3600/sorbentF/0.10E-01) + 0.14662 * ((unitD(a)/3.28084)/unitDLX(a)/18.)**2.0 + 1178.2448 * (unitDx(a)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**2.0 + 1.916236 * (unitDLX(a)/sorbentF/0.11E-05)**2.0 - 2.556486 * (unitDLX(a)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**2.0 +  P('BFBads4',a)  -  N('BFBads4',a))*x(a,'BUF') + gasInV(a,'T')*(1-y(a));



BFBads6(a) $ads('6') .. coldOutT(a) =E=
        (0.3024012E-002 * (unitL(a)/3.28084)/8.0 + 37.205355)*x(a,'BOF')+(0.302653E-002 * (unitL(a)/3.28084)/8.0 + 36.6582)*x(a,'BUF')
        +  P('BFBads6',a)  -  N('BFBads6',a);

BFBads7(a) $ads('7') .. solidOutT(a) =E=
        (- 1.15518 * EXP(unitDLX(a)) + 31.518829 * EXP(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))+ 9.898139 *((solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03)**3.0 - 36.22825 * unitDLX(a)**4.0 + 2.153780 * ((unitD(a)/3.28084)*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/22)**(1/2)  + 0.58171813 * ((unitL(a)/3.28084)*unitDx(a)/8.0)**(1/2) - 0.8213056E-001 * (unitDLX(a)*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**(-1.0) - 0.11436E-003 * ((unitL(a)/3.28084)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/8.0)**(-2.0) - 0.1485948E-006 * (unitDx(a)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(-2.0)  + 0.14179E-004 * (unitDLX(a)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(-2.0)  - 0.40372216 * (gasIn(a)*3600*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.99E+06)**(-2.0) + 71.35173 * (unitDLX(a)*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03)**(1/2) - 20.08101 * (gasIn(a)*3600*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.99E+06)**(1/2) - 4.3491254 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(1/2) - 519.8079 * ((unitL(a)/3.28084)*unitDx(a)/8.0) + 15.29240 * (unitDLX(a)*gasIn(a)*3600/0.90E+04)  + 110.25239 * unitDLX(a)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)) - 21.403203 * (unitDLX(a)*SorbentF/0.90E+06) + 20.3821293 * (sorbentF*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.99E+08)  - 1.374 * ((solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)  + 12.9239 * ((unitD(a)/3.28084)*unitDLX(a)/18)**2.0 - 5.62224 * ((unitD(a)/3.28084)*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.20E+04)**2.0  + 3874.533356 * ((unitL(a)/3.28084)*unitDx(a)/8.0)**2.0 + 18.462012 * ((unitL(a)/3.28084)*unitDLX(a)/8.0)**2.0 - 6.61325 * ((unitL(a)/3.28084)*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.88E+03)**2.0  + 2.287 * ((unitL(a)/3.28084)*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/9.6)**2.0  + 1022.67995 * (unitDx(a)*sorbentF/0.90E+06)**2.0 - 55.84266* (unitDLX(a)*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a))))**2.0 - 13.701694 * (unitDLX(a)*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**2.0 - 21.62892 * (unitDLX(a)*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.2)**2.0 + 13.28621 * ((unitD(a)/3.28084)/gasInV(a,'T')/0.18)**(1/2) - 2.7560 * (unitDLX(a)/sorbentF/0.11E-05)**(1/2) + 0.11906 * (unitDx(a)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**(-1.0) - 0.36450 * (unitDLX(a)/gasInV(a,'T')/0.10E-01)**(-1.0)  + 2.20296 * (gasInV(a,'T')/sorbentF/0.11E-03)**(-1.0)- 2.26983 * ((unitL(a)/3.28084)/unitDLX(a)/8.0)**(-2.0)  - 0.4208E-002 * (unitDx(a)/unitDLX(a))**(-2.0)  + 0.1402259E-002 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**(-2.0)  - 1.5913 * (gasOutV(a,'P')/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**(-2.0) - 0.927368 * (sorbentF/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.82E+04)**(-2.0) + 0.4185304 * ((solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/92.)**(-2.0)  - 1.99786 * (unitDLX(a)/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(1/2) + 6.907554 * (unitDLX(a)/sorbentF/0.11E-05)**(1/2) + 14.79475 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.83)**(1/2)  - 1.911873 * (sorbentF/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.75E+06)**(1/2)  - 1.47872 * ((unitD(a)/3.28084)/unitDLX(a)/18)  + 14.64508 * unitDx(a)/unitDLX(a) - 1.845 * (unitDLX(a)/gasInV(a,'T')/0.10E-01) + 0.472703 * ((unitD(a)/3.28084)/(unitL(a)/3.28084)/2.3)**2.0 + 0.60163E-001 * ((unitD(a)/3.28084)/unitDLX(a)/18)**2.0 +  P('BFBads7',a)  -  N('BFBads7',a))*x(a,'BOF')+(25.959 * exp(unitDLX(a)) - 7.385462 * unitDLX(a)**(1/2) + 18.3336 * ((unitL(a)/3.28084)*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/9.6)**(1/2) - 3.1209 * ((unitL(a)/3.28084)*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/9.6)**(-1.0) - 0.303548 * (unitDLX(a)*gasOutV(a,'P')/1.4)**(-1.0) - 0.65271 * (unitDLX(a)*sorbentF/0.90E+06)**(-1.0) - 0.14598 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**(-1.0) - 0.34918 * (sorbentF*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/0.11E+07)**(-1.0) - 0.7650E-003 * ((unitD(a)/3.28084)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/18)**(-2.0) - 0.6412989E-001 * ((unitL(a)/3.28084)*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.88E+03)**(-2.0) + 0.520695E-001 * ((unitL(a)/3.28084)*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/9.6)**(-2.0) + 0.6166247E-002 * (unitDLX(a)*sorbentF/0.90E+06)**(-2.0) + 0.4030313E-003 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*sorbentF/0.90E+06)**(-2.0) + 0.321337E-003 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**(-2.0) - 220.609 * (unitDx(a)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(1/2)  + 142.207 * (unitDLX(a)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(1/2) + 191.6102 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(1/2) - 633.476 * (gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)) + 23.32949 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*sorbentF/0.90E+06) - 11.74624 * ((solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))/0.11E+03) - 14.279892 * (unitDLX(a)*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**2.0 - 20.55699 * (unitDLX(a)*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.2)**2.0 - 286.5945418 * (gasIn(a)*3600*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/0.90E+04)**2.0 + 9.5612 * (gasIn(a)*3600/gasInV(a,'T')/90.)**(1/2) + 1.735135 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(-1.0) - 1.283934 * ((unitD(a)/3.28084)/sorbentF/0.20E-04)**(-2.0) + 4.43777 * ((unitD(a)/3.28084)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.16)**(-2.0) + 133.9396 * ((unitL(a)/3.28084)/unitDx(a)/8.0)**(-2.0) - 1.2228475 * ((unitL(a)/3.28084)/unitDLX(a)/8.0)**(-2.0) - 0.2325795E-002 * (unitDx(a)/unitDLX(a))**(-2.0) + 1.1802 * (unitDLX(a)/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(-2.0) - 0.278887E-001 * (unitDLX(a)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**(-2.0) - 0.60720 * (sorbentF/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.82E+04)**(-2.0) - 39.89957 * (unitDx(a)/gasOutV(a,'P')/0.71)**(1/2) - 10.44641 * (unitDx(a)/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/0.83)**(1/2) + 10.231 * ((solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/92.)**(1/2) + 0.640852 * ((unitD(a)/3.28084)/unitDLX(a)/18) + 0.2748* unitDLX(a)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)) + 0.1486785E-001 * ((unitL(a)/3.28084)/unitDLX(a)/8.0)**2.0 + 0.1761E-002 * ((unitL(a)/3.28084)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/8.0)**2.0 - 1.583679 * (unitDx(a)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**2.0 + 696.950* (unitDx(a)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**2.0 - 5.881363 * (unitDLX(a)/(gasIn(a)*3600)/0.11E-03)**2.0 + 0.59877E-002 * (unitDLX(a)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**2.0 - 10.34750 * (unitDLX(a)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**2.0 - 24.45096 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/sorbentF/0.11E-05)**2.0 + 16.28740 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**2.0 - 0.9060990E-001 * ((solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a))))**2.0 +  P('BFBads7',a)  -  N('BFBads7',a))*x(a,'BUF')+(solidOutT(a+1)$(ORD(a) ne CARD(a)) + solidLeanT$(ORD(a) eq CARD(a)))*(1-y(a))
        ;


BFBads8(a) $ads('8') .. solidOutC(a,'H2O') =E=
        (0.6686E-003 * ((solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03)**(-5.0) - 0.10838E-002 * ((unitD(a)/3.28084)*unitDx(a)/18)**(-1.0) - 0.1413486E-002 * (unitDLX(a)*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(-1.0)  + 0.104891 * (gasIn(a)*3600*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.11E+05)  + 1.27717 * ((unitD(a)/3.28084)*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))/18)**2.0 - 0.15432E-004 * (unitDLX(a)*gasIn(a)*3600/0.90E+04)**(-3.0) + 0.126924E-004 * (unitDLX(a)*sorbentF/0.90E+06)**(-3.0) + 0.55399 * (unitDLX(a)/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/0.83)**(1/2) - 0.11600E-001 * (unitDx(a)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(-1.0)  - 0.15206E-001 * ((unitL(a)/3.28084)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.73E-01)**(-2.0)  - 0.138491E-002 * (unitDLX(a)/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/0.83)**(-2.0) + 0.290645E-001 * (gasInV(a,'T')/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/83)**(-2.0)  + 0.73234E-002 * ((unitL(a)/3.28084)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/8.0)  + 2.860378 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/sorbentF/0.11E-05)**2.0  + 3.48062 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**2.0 - 0.1986E-001 * ((unitD(a)/3.28084)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.16)**(-4.0)+  P('BFBads8',a)  -  N('BFBads8',a))*x(a,'BOF')+(0.17630 * log((unitD(a)/3.28084)/18) + 0.180 * log((unitL(a)/3.28084)/8.0) + 0.95201E-002 * (gasInV(a,'T')*sorbentF/0.90E+08)**(-2.0) + 0.298E-004 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*sorbentF/0.90E+06)**(-2.0) - 0.49365 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*gasOutV(a,'p')/1.4) + 0.4783518 * ((solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2) + 0.1268 * ((solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.4)**2.0 + 1.2382 * (unitDLX(a)/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(1/2) - 0.3843E-001 * (unitDx(a)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(-1.0) + 0.3718 * ((solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/92)**(-1.0) - 0.9343E-001 * (unitDLX(a)/sorbentF/0.11E-05)**2.0 + 0.3687E-001 * (gasIn(a)*3600/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/82)**2.0 +  P('BFBads8',a)  -  N('BFBads8',a))*x(a,'BUF')+(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a)) + solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))*(1-y(a))
        ;

BFBads9(a) $ads('9') .. solidOutC(a,'HCO3') =E=
        (0.71235E-001 * log((unitD(a)/3.28084)/18) + 0.21968E-001 * log((unitL(a)/3.28084)/8.0) - 0.1821E-003 * (unitDLX(a)*3600*gasIn(a)/0.90E+04)**(-2.0) + 0.38000 * ((solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.2)**0.50 + 0.82535 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2) - 0.81131E-001 * (sorbentF*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.99E+08)+ 22.6356 * (unitDx(a)*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a))))**2.0  + 0.46177E-001 * ((solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.4)**2.0 + 0.2558 * (unitDLX(a)/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**0.50 + 0.25580 * (unitDLX(a)/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**0.50 - 0.92364E-002 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(-1.0)  + 0.18241E-002 * (gasIn(a)*3600/sorbentF/0.10E-01)**4.0 +  P('BFBads9',a)  -  N('BFBads9',a))*x(a,'BOF')+(- 0.340231E-001 * log((unitD(a)/3.28084)/18) + 0.90930E-002 * log((unitL(a)/3.28084)/8.0) + 0.2947 * (solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))+ 0.34733E-002 * (unitDLX(a)*sorbentF/0.90E+06)**(-1.0) - 0.259236E-002 * (unitDLX(a)*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03)**(-1.0) + 0.9216E-006 * ((unitL(a)/3.28084)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/8.0)**(-2.0) - 0.12323936E-003 * (unitDLX(a)*gasIn(a)*3600/0.90E+04)**(-2.0) - 0.84792E-005 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*sorbentF/0.90E+06)**(-2.0)  + 0.1040369E-005 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.2)**(-2.0) - 0.525494 * (unitDLX(a)*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a))))**0.50 + 4.86616 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**0.50 + 0.32215 * ((solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.2)**0.50 + 0.553477533 * ((unitL(a)/3.28084)*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/8.0) + 4.170886 * unitDx(a)*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a))) - 0.180667 * (unitDLX(a)*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03) + 0.156414 * (unitDLX(a)*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2) - 0.816652 * (gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a))) - 1.94865099 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))*sorbentF/0.90E+06) + 0.81967E-001 * ((solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2) + 0.5458E-001 * ((solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.4) + 3.3669 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**2.0 - 1.04687 * ((unitD(a)/3.28084)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/18)**0.50 + 0.354E-001 * ((unitD(a)/3.28084)/sorbentF/0.20E-04)**(-2.0) - 0.3335E-002 * (sorbentF/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.75E+06)**(-2.0) + 0.8465E-001 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**0.50 + 0.44723E-002 *unitDLX(a)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)) + 0.13490E-001 * (gasIn(a)*3600/sorbentF/0.10E-01)**2.0 + 1.030842 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**2.0 +  P('BFBads9',a)  -  N('BFBads9',a))*x(a,'BUF')+(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a)) + solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))*(1-y(a))
        ;

BFBads10(a) $ads('10') .. solidOutC(a,'NH2COO') =E=
        (1.0659 * ((solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**0.50 - 21.1084 * (gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))**3.0  - 0.213E-001 * ((solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03)**4.0 - 0.1922E-002 * (gasIn(a)*3600/0.90E+04)**(-4.0) - 0.90852E-003 * ((unitD(a)/3.28084)*unitDx(a)/18.)**(-1.0) - 0.4524E-002 * (unitDLX(a)*gasIn(a)*3600/0.90E+04)**(-1.0) - 0.398E-003 * (unitDLX(a)*gasOutV(a,'P')/1.4)**(-1.0) - 0.727E-002 * (unitDLX(a)*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03)**(-1.0) - 0.162E-002 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03)**(-1.0) - 0.2690E-003 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**(-1.0) + 0.2936E-005 * ((unitL(a)/3.28084)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/8.0)**(-2.0) - 0.271878E-005 * ((unitL(a)/3.28084)*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/8.0)**(-2.0)  + 0.526911E-008 * (unitDx(a)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(-2.0) + 0.240787E-007 * (unitDLX(a)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**(-2.0) + 0.19237E-003 * (unitDLX(a)*gasOutV(a,'P')/1.4)**(-2.0)  + 0.1004E-004 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))*sorbentF/0.90E+06)**(-2.0) + 0.33673 * ((unitL(a)/3.28084)*unitDx(a)/8.0)**0.50 + 0.63478E-001 * ((unitL(a)/3.28084)*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))/8.0) + 0.646825 * unitDLX(a)*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a))) - 0.413834 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.2)- 0.19234 * (gasOutV(a,'P')*sorbentF/0.13E+07) + 0.16743 * (sorbentF*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))/0.90E+06) + 0.38727 * (sorbentF*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.11E+07) - 0.174251 * ((solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.13E+03)- 8.9684 * (unitDLX(a)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**2.0 - 0.38183 * (gasIn(a)*3600*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))/0.90E+04)**2.0 - 0.12677 * (gasIn(a)*3600*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.11E+05)**2.0  + 0.404908E-001 * (gasInV(a,'T')*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))/0.10E+03)**2.0  - 4.056 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**2.0 + 0.11447 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.2)**2.0  + 0.3141 * (unitDLX(a)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**0.50 + 0.677E-001 * (unitDLX(a)/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.83)**0.50 - 1.2337 * (gasIn(a)*3600/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/0.90E+04)**0.50 + 0.9635E-003 * (unitDx(a)/sorbentF/0.11E-05)**(-1.0)  - 2.812 * (gasIn(a)*3600/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/0.90E+04)**(-1.0) + 1.2824 * (gasIn(a)*3600/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/0.90E+04)**(-1.0) + 0.11977E-001 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(-1.0)- 0.341E-002 * ((solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a))))**(-1.0)  - 0.1877E-001 * ((unitD(a)/3.28084)/gasInV(a,'T')/0.18)**(-2.0) - 0.52227E-002 * ((unitL(a)/3.28084)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.73E-01)**(-2.0)  + 3.747558 * (gasIn(a)*3600/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/0.90E+04)**(-2.0) + 0.328E-001 * (gasIn(a)*3600/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.75E+04)**(-2.0) + 0.25988E-004 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.83)**(-2.0)  - 0.51926E-001 * (gasOutV(a,'P')/sorbentF/0.16E-05)**(-2.0)  - 0.14971E-001 * (sorbentF/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.75E+06)**(-2.0) - 0.99875E-001 * (unitDLX(a)/(gasIn(a)*3600)/0.11E-03)**0.50  - 0.206E-001 * (unitDLX(a)/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**0.50  + 2.392* ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/sorbentF/0.11E-05)**0.50 + 0.57171E-001 * (gasIn(a)*3600/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/82.) + 0.762 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.83) - 0.4476E-003 * ((unitD(a)/3.28084)/unitDLX(a)/18.)**2.0 + 0.831E-002 * ((unitL(a)/3.28084)/sorbentF/0.89E-05)**2.0 + 0.833E-001 * (unitDLX(a)/(gasIn(a)*3600)/0.11E-03)**2.0 - 0.6510E-001 * (unitDLX(a)/sorbentF/0.11E-05)**2.0 + 0.278E-002 * (gasIn(a)*3600/sorbentF/0.10E-01)**2.0 - 1.4128 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/sorbentF/0.11E-05)**2.0 - 0.3527 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.83)**2.0 +  P('BFBads10',a)  -  N('BFBads10',a))*x(a,'BOF')+(0.3998 * (gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))**0.50 + 1.15707 * ((solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**2.0 - 40.9811 * (gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))**3.0 - 0.5635 * ((solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**3.0 - 0.3497E-001 * ((unitD(a)/3.28084)*(unitL(a)/3.28084)/0.14E+03)**(-1.0) - 0.1665E-001 * ((solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.13E+03)**(-1.0) -0.204E-002 * ((unitD(a)/3.28084)*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/22.)**(-2.0) + 0.1143E-005 * ((unitL(a)/3.28084)*unitDx(a)/8.0)**(-2.0) - 0.195E-005 * (unitDx(a)*sorbentF/0.90E+06)**(-2.0) + 0.8539 * (unitDx(a)*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a))))**0.50 - 0.3860 * (gasInV(a,'T')*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/0.10E+03)**0.50 + 0.431072 * unitDLX(a)*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a))) - 2.1696 *(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))- 0.7325 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2) + 0.3197 * (sorbentF*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))/0.90E+06) + 0.3928 * (sorbentF*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.11E+07)- 0.256 * ((solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2) + 0.6625E-001 * ((unitL(a)/3.28084)*gasOutV(a,'P')/11.)**2.0 - 0.5601 * (unitDLX(a)*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03)**2.0- 0.291490 * (gasIn(a)*3600*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.11E+05)**2.0 - 5.42990 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**2.0 + 0.8086 * (unitDLX(a)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**0.50 - 0.7915 * (gasIn(a)*3600/sorbentF/0.10E-01)**0.50 + 1.96842 * (sorbentF/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.75E+06)**0.50 - 0.1223E-002 * (unitDx(a)/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.83)**(-1.0) - 0.286E-001 * (unitDLX(a)/sorbentF/0.11E-05)**(-1.0) - 0.9397 * (sorbentF/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.75E+06)**(-1.0)- 0.06783 * ((unitD(a)/3.28084)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.16)**(-2.0) + 0.373 * (gasIn(a)*3600/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/0.90E+04)**(-2.0) + 0.288E-001 * (gasIn(a)*3600/sorbentF/0.10E-01)**(-2.0) - 0.644E-002 * (gasIn(a)*3600/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/0.75E+04)**(-2.0) + 0.126E-003 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/sorbentF/0.11E-05)**(-2.0) + 0.737E-001 * (sorbentF/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.75E+06)**(-2.0) + 1.4450* ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/0.83)**0.50 - 0.62E-004 * ((unitD(a)/3.28084)/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/18)**2.0 + 0.59E-001 * (unitDLX(a)/(gasIn(a)*3600)/0.11E-03)**2.0 - 0.61003 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/sorbentF/0.11E-05)**2.0 +  P('BFBads10',a)  -  N('BFBads10',a))*x(a,'BUF')+(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a)) + solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))*(1-y(a))
        ;


BFBads11(a) $ads('11') .. gasOutX(a,'CO2') =E=
        (0.1402E-002 * log((unitD(a)/3.28084)/18.)-0.360E-002 * log((unitL(a)/3.28084)/8.0)+0.497739 * unitDLX(a)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))+ 0.4721411 * (gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))+0.3819 * (gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2+3.602 * (gasIn(a)*3600*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/0.90E+04)**3.0 -0.404E-002 * (gasIn(a)*3600/sorbentF/0.10E-01)**(-1.0)+0.289E-004 * (unitDLX(a)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**(-2.0)+0.342E-002 * ((solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a)))/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/0.83)- 0.9894 * (unitDx(a)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**2.0 + 0.6206 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/sorbentF/0.11E-05)**3.0+ 0.132E-002 * ((unitD(a)/3.28084)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.16)**(-4.0) +  P('BFBads11',a)  -  N('BFBads11',a))*x(a,'BOF')+(- 0.99392E-002 * log((unitD(a)/3.28084)/18.) - 0.1340E-002 * log((unitL(a)/3.28084)/8.0) -0.1375E-003 * (unitDLX(a)*gasIn(a)*3600/0.90E+04)**(-1.0)  + 0.3520 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2) + 2.5926 * (unitDLX(a)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)))**2.0  + 1.0774 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03)**2.0 - 0.38920 * (gasIn(a)*3600/(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/0.90E+04)**(-2.0) + 0.3608 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))/sorbentF/0.11E-05)**2.0 +  P('BFBads11',a)  -  N('BFBads11',a))*x(a,'BUF')+(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(1-y(a))
        ;

BFBads12(a) $ads('12') .. gasOutX(a,'H2O') =E=
        (0.9974E-002 * ((solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03)**4.0 + 0.5113E-004 * ((solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**(-3.0) - 0.4655E-003 * (gasIn(a)*3600/0.90E+04)**(-4.0) - 0.180822E-002 * (gasInV(a,'T')*sorbentF/0.90E+08)**(-1.0) - 0.155368E-003 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.11E+03)**(-1.0) + 0.1877E-005 * (gasIn(a)*3600*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/0.90E+04)**(-2.0) - 0.77853E-001 * ((unitD(a)/3.28084)*(unitL(a)/3.28084)/0.14E+03)**0.50 + 0.5526 * (unitDLX(a)*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**0.50 + 0.305560E-001 * ((solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/0.13E+03)**0.50 + 0.192519* unitDLX(a)*(gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)) + 0.2691E-001 * (unitDLX(a)*gasInV(a,'T')/0.10E+03)**2.0 + 0.2034 * (unitDLX(a)*(solidOutC(a+1,'HCO3')$(ORD(a) ne CARD(a))+solidOutC('d1','HCO3')$(ORD(a) eq CARD(a))))**2.0 - 0.1260* (unitDLX(a)*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**3.0  + 0.2608E-009 * (unitDLX(a)*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.2)**(-4.0) + 0.2715E-002 * (unitDx(a)/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1)))**(-1.0) + 0.30626E-002 * ((unitD(a)/3.28084)/(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/15.)**(-2.0)  + 0.37285E-003 * ((unitL(a)/3.28084)/(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/8.0) + 0.866E-002 * ((unitL(a)/3.28084)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.73E-01) - 0.1914E-001 * (unitDLX(a)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**2.0 + 0.104 * ((gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/sorbentF/0.11E-05)**2.0 - 0.2254E-002 * ((unitL(a)/3.28084)/unitDLX(a)/8.0)**(-3.0) - 0.268E-001 * ((unitD(a)/3.28084)/unitDLX(a)/18.)**(-4.0) + 0.53449991E-003 * (gasIn(a)*3600/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/82.)**(-4.0) - 0.1942E-003 * (sorbentF/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.82E+04)**(-4.0)+  P('BFBads12',a) - N('BFBads12',a))*x(a,'BOF')+(- 0.186E-001 * log((unitD(a)/3.28084)/18.) - 0.88691E-002 * log((unitL(a)/3.28084)/8.0) - 0.6612E-003 * (unitDLX(a)*(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/1.2)**(-1.0)  - 0.1096E-003 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.2)**(-1.0) + 1.366 * unitDLX(a)*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))  + 0.1533 * (gasInV(a,'T')*(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))/0.10E+03) + 1.3377 * ((gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1))*(solidOutC(a+1,'H2O')$(ORD(a) ne CARD(a))+solidOutC('d1','H2O')$(ORD(a) eq CARD(a)))/1.2)**2.0 + 0.4972E-005 * (unitDx(a)/(solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/0.91E-02)**(-2.0) + 0.956436E-002 * ((solidLeanT$(ORD(a) eq CARD(a))+solidOutT(a+1)$(ORD(a) ne CARD(a)))/(solidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a))+solidOutC('d1','NH2COO')$(ORD(a) eq CARD(a)))/92.)+  P('BFBads12',a) - N('BFBads12',a))*x(a,'BUF')+(gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))*(1-y(a))
        ;


* BFB regenerator surrogate model

BFBDes1(d) $des('1') .. unitNx(d) =E=
        (- 0.2852 * unitDx(d)**(-2.0) - 0.140689E-001 * unitDLX(d)**(-4.0) - 363.718 * ((unitD(d)/3.28084)/16.)**4.0 + 0.2090E-002 * (unitDx(d)*unitDLX(d))**(-2.0) + 21.872 * (unitDx(d)/unitDLX(d))**(-1.0) + 6779.0566 * ((unitD(d)/3.28084)/unitDLX(d)/16.)**(-2.0) + 0.11833 * (unitDx(d)/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**(-2.0) - 1034.1328 * ((unitD(d)/3.28084)/unitDx(d)/16.)**0.50 + 103.037 * ((unitD(d)/3.28084)/unitDx(d)/16.) + 1111.12949 * ((unitD(d)/3.28084)/unitDLX(d)/16.) + 55.24623 * ((unitD(d)/3.28084)/unitDLX(d)/16.)**2.0 - 20302.8303 * (unitDx(d)/unitDLX(d))**2.0)*x(d,'BOF')+(0.11222E-002 * (unitDx(d)*unitDLX(d))**(-2.0) - 0.1042 * (unitDx(d)*(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**(-2.0) + 27.2640 * (unitDx(d)/unitDLX(d))**(-1.0) + 7748.8387 * ((unitD(d)/3.28084)/unitDLX(d)/16.)**(-2.0) - 1161.5500 * ((unitD(d)/3.28084)/unitDx(d)/16.)**0.50 + 95.13755 * ((unitD(d)/3.28084)/unitDx(d)/16.) + 1284.499 * ((unitD(d)/3.28084)/unitDLX(d)/16.) + 44.835090 * ((unitD(d)/3.28084)/unitDLX(d)/16.)**2.0 - 23227.2978 * (unitDx(d)/unitDLX(d))**2.0)*x(d,'BUF')
        +  P('BFBDes1',d)  -  N('BFBDes1',d);

BFBDes2(d) $des('2') .. gasInV(d,'P') =E=
        (0.485193 * exp(gasOutV(d,'P')/1.3) + 0.18412 * ((unitL(d)/3.28084)*gasOutV(d,'P')/10) +  P('BFBDes2',d)  -  N('BFBDes2',d))*x(d,'BOF')+(0.48509 * exp(gasOutV(d,'P')/1.3) + 0.18346 * ((unitL(d)/3.28084)*gasOutV(d,'P')/10.) +  P('BFBDes2',d)  -  N('BFBDes2',d))*x(d,'BUF') + gasOutV(d,'P')*(1-y(d));

BFBDes3(d) $des('3') .. 3600*gasOut(d) =E=
        (701.7595 * log((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03) + 140.0799 * log((solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3) - 8.774868 * (sorbentF/0.90E+06)**(-3.0) + 3158.840 * (gasIn(d)*3600)/0.30E+04 + 153.18180 * ((unitD(d)/3.28084)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/21.)**0.50 + 100.991 * ((unitL(d)/3.28084)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/10.)**0.50 - 81.1949 * (unitDLX(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**0.50 + 95.0851 * ((unitL(d)/3.28084)*sorbentF/0.72E+07)**(-1.0) + 0.815219E-001 * (unitDx(d)*unitDLX(d))**(-1.0) - 0.149963 * (unitDx(d)*(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**(-1.0) + 0.28081 * (unitDx(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**(-1.0) + 2.8331 * (gasIn(d)*3600*(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/0.30E+04)**(-1.0) - 3.7096 * ((unitL(d)/3.28084)*sorbentF/0.72E+07)**(-2.0) - 0.84035E-003 * (unitDx(d)*(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**(-2.0) + 0.3529E-002 * (unitDLX(d)*(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1)))**(-2.0) - 0.229294 * (unitDLX(d)*sorbentF/0.90E+06)**(-2.0) + 29.2433 * (sorbentF*(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.14E+09)**(-2.0) + 0.3273 * (sorbentF*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.12E+07)**(-2.0)- 0.4462E-001 * ((solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**(-2.0) - 1.3932 * ((solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/2.3)**(-2.0) + 148.046 * (gasIn(d)*3600*sorbentF/0.27E+10) + 573.1971 * (gasOutV(d,'P')*sorbentF/0.12E+07) + 1042.032 * (sorbentF*(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/0.90E+06) + 218.6349 * (sorbentF*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.16E+07) + 311.415 * (sorbentF*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.12E+07) + 82.12791 * ((unitD(d)/3.28084)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/29.)**2.0 + 97.4603 *((unitL(d)/3.28084)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/10.)**2.0 - 103581.5106 * (unitDx(d)*(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**2.0 + 103581.510 * (unitDx(d)*(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**2.0 + 90198.1519 * (unitDx(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**2.0 + 54544.2724 * (unitDx(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**2.0 - 1359.4401 * (unitDLX(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**2.0 + 540.61946 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**2.0 - 80.12754 * ((unitL(d)/3.28084)/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/8.0)**0.50 - 796.3093 * ((unitL(d)/3.28084)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/4.4)**0.50 - 52.502 * (unitDx(d)/sorbentF/0.11E-05)**0.50 - 140.5114 * (unitDLX(d)/gasOutV(d,'P')/0.77)**0.50 + 110.7108 * (unitDLX(d)/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**0.50 + 329.7823 * (unitDLX(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**0.50 - 288.78823 * ((unitD(d)/3.28084)/unitDLX(d)/16.)**(-1.0) + 8.56917 * (unitDLX(d)/(gasIn(d)*3600)/0.33E-03)**(-1.0) - 28.74580 * ((unitD(d)/3.28084)/(gasIn(d)*3600)/0.53E-02)**(-2.0) - 16.644 * ((unitD(d)/3.28084)/sorbentF/0.18E-04)**(-2.0) + 68.054 * ((unitL(d)/3.28084)/unitDLX(d)/8.0)**(-2.0) + 10.0093 * ((unitL(d)/3.28084)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/4.4)**(-2.0) - 11.32355 * (gasIn(d)*3600/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.23E+04)**(-2.0) - 0.1365 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.77)**(-2.0) + 207.4231 * (gasOutV(d,'P')/sorbentF/0.14E-05)**(-2.0) + 20.921* ((solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.4)**(-2.0) - 90.8309 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/sorbentF/0.11E-05) + 0.54 * ((unitD(d)/3.28084)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/12.)**2.0 + 228.30422 * (unitDx(d)/(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1)))**2.0 + 1441.510 * (unitDx(d)/sorbentF/0.11E-05)**2.0 - 26.04843 * (gasIn(d)*3600/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.17E+04)**2.0 - 539.3837 * (sorbentF/(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.60E+04)**2.0 +  P('BFBDes3',d)  -  N('BFBDes3',d))*x(d,'BOF')+(- 0.1880E-003 * unitDx(d)**(-3.0) + 4465.72943 * gasIn(d)*3600/0.30E+04+ 163.671 * ((solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**4.0 - 65.57033 * (unitDLX(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**0.50 - 47.795058 * ((unitD(d)/3.28084)*sorbentF/0.14E+08)**(-1.0) - 14.1722 * ((solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**(-1.0) + 0.80188E-002 * (unitDx(d)*sorbentF/0.90E+06)**(-2.0) - 0.36060 * (unitDLX(d)*sorbentF/0.90E+06)**(-2.0) + 3.0436 * (sorbentF*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.12E+07)**(-2.0) +258.949 * (unitL(d)/3.28084)*(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/8.0 + 169.616 * ((unitD(d)/3.28084)*gasIn(d)*3600/0.48E+05)**2.0 + 339.143 * ((unitL(d)/3.28084)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/14.)**2.0 + 108583.52805 * (unitDx(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**2.0 - 3777.5195 * (unitDLX(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**2.0 + 729.571 * (gasInV(d,'T')*(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/0.17E+03)**2.0 + 1237.839 * ((solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**2.0 + 65.9550 * (unitDLX(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**(-1.0) - 61.8721 * ((unitD(d)/3.28084)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/8.9)**(-2.0) - 12.145 * ((unitL(d)/3.28084)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/6.2)**(-2.0) - 0.2990E-001 * (unitDx(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**(-2.0) - 5.047 * (unitDLX(d)/(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1)))**(-2.0) + 0.5947 * (unitDLX(d)/sorbentF/0.11E-05)**(-2.0) + 18.584 * ((solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.4)**(-2.0) + 308.8718 * ((unitL(d)/3.28084)/unitDLX(d)/8.0)**0.50 - 253.705 * ((unitL(d)/3.28084)/sorbentF/0.89E-05) - 1399.3043 * (gasIn(d)*3600/(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/20.) - 1.8656 * ((unitL(d)/3.28084)/unitDLX(d)/8.0)**2.0 + 1539.409 * (unitDx(d)/(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1)))**2.0 + 21.14924 * (unitDLX(d)/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**2.0 + 185.387 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/sorbentF/0.11E-05)**2.0 +  P('BFBDes3',d)  -  N('BFBDes3',d))*x(d,'BUF')+3600*((feedCO2F+steamF)$(ORD(d) eq 1)+gasOut(d-1)$(ORD(d) ne 1))*(1-y(d))
       ;

BFBDes4(d) $des('4') .. gasOutV(d,'T') =E=
        (202.864 * exp(unitDx(d)) - 21.383 * ((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**(-1.0) + 0.312E-002 * (unitDx(d)*unitDLX(d))**(-1.0) - 0.1059 * (sorbentF*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.16E+07)**(-2.0) - 78.869 * (unitDLX(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**0.50 + 34.111 * ((unitD(d)/3.28084)*unitDLX(d)/16.) - 8.944 * (sorbentF*(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/0.90E+06) + 4.8397 * ((unitL(d)/3.28084)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/14.)**2.0 - 94.3711 * (unitDLX(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**2.0 + 8.607 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**2.0 - 10.848 * ((unitL(d)/3.28084)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/6.2)**0.50 - 6.4174 * (gasInV(d,'T')/sorbentF/0.19E-03)**0.50 - 0.8035 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**0.50 - 1.1125 * ((unitD(d)/3.28084)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/8.9)**(-2.0) + 2.46734 * ((unitL(d)/3.28084)/unitDLX(d)/8.0)**(-2.0) - 0.5383E-003 * ((unitD(d)/3.28084)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**(-2.0) - 0.4016E-003 * (unitDx(d)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.77)**(-2.0) - 0.354E-001 * (unitDLX(d)/(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1)))**(-2.0) + 0.482 * (sorbentF/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.69E+06)**(-2.0) + 3.101 * (gasIn(d)*3600/sorbentF/0.33E-02)**0.50 + 5.949* ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56) - 0.6109 * ((unitD(d)/3.28084)/(gasIn(d)*3600)/0.53E-02)**2.0 - 419.1247 * (unitDx(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**2.0 + 5.634 * (unitDLX(d)/(gasIn(d)*3600)/0.33E-03)**2.0+  P('BFBDes4',d)  -  N('BFBDes4',d))*x(d,'BOF')+(126.99* exp(unitDx(d)) + 0.14E-001 * ((unitL(d)/3.28084)/8.0)**(-4.0) + 55.147 * ((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**0.50 + 1.47 * ((unitL(d)/3.28084)/8.0)**4.0 + 0.26E-001 * ((unitL(d)/3.28084)*unitDx(d)/8.0)**(-1.0) - 0.360E-001 * (unitDx(d)*sorbentF/0.90E+06)**(-1.0) - 0.6627 * ((unitD(d)/3.28084)*gasOutV(d,'P')/21.)**(-2.0) - 103.306 * (unitDLX(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**0.50 + 7.4981 * ((unitL(d)/3.28084)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/14.) + 35.229 * (unitDLX(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8) - 7.7828 * (gasOutV(d,'P')*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/2.3) - 7.72593 * (sorbentF*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.12E+07) - 139.0107 * ((unitL(d)/3.28084)*unitDLX(d)/8.0)**2.0 + 138.36018 * (unitDLX(d)*(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**2.0 + 5.26508 * (unitDLX(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**2.0 - 11.7906777 * ((unitL(d)/3.28084)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/6.2)**0.50 - 7.224004 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**0.50 - 4.09384 * ((unitD(d)/3.28084)/(gasIn(d)*3600)/0.53E-02)**(-1.0) + 1.848 * (gasOutV(d,'P')/(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.87E-02)**(-1.0) + 0.303 * ((unitL(d)/3.28084)/sorbentF/0.89E-05)**(-2.0) + 0.54254E-002 * (unitDx(d)/unitDLX(d))**(-2.0) - 0.7695E-003 * (unitDx(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**(-2.0) + 0.1352E-001 * (unitDLX(d)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.77)**(-2.0) + 0.257886 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**(-2.0) - 0.36333E-001 * ((solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.77)**(-2.0) + 18.58 * (gasIn(d)*3600/sorbentF/0.33E-02)**0.50 - 3.222 * ((unitD(d)/3.28084)/(unitL(d)/3.28084)/2.0) + 0.1646E-001 * ((unitD(d)/3.28084)/(gasIn(d)*3600)/0.53E-02)**2.0 + 0.1583E-001 * ((unitD(d)/3.28084)/(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/16.)**2.0 + 0.45710 * ((unitL(d)/3.28084)/sorbentF/0.89E-05)**2.0 - 11.9169 * (unitDLX(d)/sorbentF/0.11E-05)**2.0 + 3.39178 * (unitDLX(d)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.77)**2.0 - 0.2567 * ((solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.4)**2.0+  P('BFBDes4',d)  -  N('BFBDes4',d))*x(d,'BUF')+ gasInV(d,'T')*(1-y(d))
        ;



BFBDes6(d) $des('6') .. hotOutT(d) =E=
        (- 0.5771 * (unitL(d)/3.28084)/8.0 + 164.7497)*x(d,'BOF')+(- 0.61194 * (unitL(d)/3.28084)/8.0 - 0.1387E-002 * ((unitL(d)/3.28084)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/14.)**(-2.0) - 13.258 * (unitDx(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**2.0 - 0.2250E-003 * ((unitL(d)/3.28084)/unitDLX(d)/8.0)**2.0 + 164.7303)*x(d,'BUF')
        +  P('BFBDes6',d)  -  N('BFBDes6',d);

BFBDes7(d) $des('7') .. solidOutT(d) =E=
  (161.094* exp(unitDx(d))- 11.2092 * ((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**(-2.0)+ 3.1605 * ((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**4.0  + 0.1356E-001 * (unitDx(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**(-1.0)+ 0.111E-005 * (unitDx(d)*(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**(-2.0)  - 0.1523E-002 * (unitDLx(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**(-2.0) - 25.8999 * ((solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/2.3)**0.50- 78.468 * (unitDLX(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8) - 7.27 * ((solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8) + 3.47399 * ((unitL(d)/3.28084)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/10.)**2.0 + 111.8711 * (unitDLX(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**2.0 - 3.73830 * ((unitD(d)/3.28084)/(gasIn(d)*3600)/0.53E-02)**0.50 - 2.42566 * ((unitD(d)/3.28084)/(unitL(d)/3.28084)/2.0)**(-1.0) - 1.0058 * ((unitL(d)/3.28084)/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/8.0)**(-1.0) + 1.7600 * ((unitL(d)/3.28084)/unitDLX(d)/8.0)**(-2.0) - 0.64602E-003 * (unitDx(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**(-2.0) + 0.24275E-001 * (unitDLX(d)/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**(-2.0) - 0.26853 * (gasIn(d)*3600/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.23E+04)**(-2.0) + 0.4851 * (gasOutV(d,'P')/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d))))**(-2.0) + 0.51140 * (sorbentF/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.69E+06)**(-2.0) - 1.7026 * ((solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.4)**(-2.0) + 2.678 * ((unitL(d)/3.28084)/unitDLX(d)/8.0)**0.50 + 9.5538 * (gasIn(d)*3600/sorbentF/0.33E-02)**0.50 + 5.8721 * ((unitL(d)/3.28084)/gasOutV(d,'P')/6.2)+ 1.701305 * ((unitL(d)/3.28084)/sorbentF/0.89E-05) - 40.32617 * (unitDx(d)/(gasIn(d)*3600)/0.33E-03) - 0.1848 * ((unitD(d)/3.28084)/(gasIn(d)*3600)/0.53E-02)**2.0 - 0.2092E-003 * ((unitL(d)/3.28084)/unitDx(d)/8.0)**2.0 - 0.7331 * ((unitL(d)/3.28084)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/4.4)**2.0 + 108.984 * (unitDx(d)/sorbentF/0.11E-05)**2.0 + 8.488277 * (unitDLX(d)/(gasIn(d)*3600)/0.33E-03)**2.0 - 9.6179* (unitDLX(d)/sorbentF/0.11E-05)**2.0 + 0.499E-001 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/gasOutV(d,'P')/0.77)**2.0 + 3.55250 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**2.0 +  P('BFBDes7',d)  -  N('BFBDes7',d))*x(d,'BOF')+(264.054 * exp(unitDx(d)) - 45.5614 * ((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**0.50 - 0.709E-003 * (unitDLX(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**(-2.0) - 0.502E-002 * ((gasIn(d)*3600)*(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/0.30E+04)**(-2.0) - 117.865 * (unitDLX(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**0.50 + 10.774 * ((unitL(d)/3.28084)*(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/8.0)**2.0 + 5.2069 * ((unitL(d)/3.28084)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/14.)**2.0 + 2122.0283 * (unitDx(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**2.0 + 68.734 * (unitDLX(d)*(gasIn(d)*3600)/0.30E+04)**2.0 - 22.733 * ((unitD(d)/3.28084)/sorbentF/0.18E-04)**0.50 - 10.089 * ((unitL(d)/3.28084)/unitDx(d)/8.0)**0.50 - 8.771 * ((unitL(d)/3.28084)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/6.2)**0.50 - 3.67025 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**0.50 - 0.8641E-001 * (unitDx(d)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.77)**(-1.0) + 0.78716E-002 * (unitDx(d)/unitDLX(d))**(-2.0) - 0.598399E-003 * (unitDx(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**(-2.0) + 0.2084E-001 * (unitDLX(d)/sorbentF/0.11E-05)**(-2.0) - 9.8210 * (unitDx(d)/unitDLX(d))**2.0 - 634.0885 * (unitDx(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**2.0 + 7.4380 * (unitDLX(d)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.77)**2.0 - 0.2346414 * (gasOutV(d,'P')/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d))))**2.0 +  P('BFBDes7',d)  -  N('BFBDes7',d))*x(d,'BUF')+(solidOutT(d+1)$(ORD(d) ne CARD(d)) + solidRichT$(ORD(d) eq CARD(d)))*(1-y(d))
       ;

BFBDes8(d) $des('8') .. solidOutC(d,'H2O') =E=
       (- 0.104 * ((unitL(d)/3.28084)/8.0)**(-1.0) - 0.351E-005 * (gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))**(-4.0) - 0.990E-002 * (unitDx(d)*gasIn(d)*3600/0.30E+04)**0.50 + 0.1010 * (gasIn(d)*3600*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.54E+04)**0.50 + 0.1062 * ((unitL(d)/3.28084)*(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.12E+04)**(-1.0) - 0.4210E-006 * (unitDx(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**(-2.0) - 0.1106E-004 * (unitDLX(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**(-2.0) + 0.4600E-005 * (unitDLX(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**(-2.0) - 0.5529 * ((unitL(d)/3.28084)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/10.)**0.50 - 7.50893 * unitDx(d)*unitDLX(d) - 0.20708 * (unitDLX(d)*gasInV(d,'T')/0.17E+03) + 0.99534 * unitDLX(d)*(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))) + 0.352 * ((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.20E+03) - 0.514E-001 * ((unitD(d)/3.28084)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/29.)**2.0 + 29.18627 * (unitDx(d)*(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**2.0 - 48.3422 * (unitDx(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**2.0 + 0.5230 * (unitDLX(d)*gasIn(d)*3600/0.30E+04)**2.0 - 3.4373 * (unitDLX(d)*(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**2.0 - 2.8390 * (unitDLX(d)*(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**2.0 + 1.6314 * (unitDLX(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**2.0 + 1.673314 * (unitDLX(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**2.0 - 0.25318 * (gasIn(d)*3600*(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/0.30E+04)**2.0 - 0.9143E-001 * (gasIn(d)*3600*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.39E+04)**2.0 - 0.1197 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**2.0 - 0.21312 * ((unitL(d)/3.28084)/unitDLX(d)/8.0)**0.50 + 0.311E-001 * ((unitD(d)/3.28084)/(gasIn(d)*3600)/0.53E-02)**(-1.0) + 0.7385E-001 * ((unitD(d)/3.28084)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/12.)**(-1.0) + 0.1702E-002 * (unitDx(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**(-1.0) - 0.307320E-002 * (unitDLX(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**(-1.0) - 0.8362E-002 * (unitDLX(d)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.77)**(-1.0) - 0.110* (gasOutV(d,'P')/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d))))**(-1.0) - 0.120 * (sorbentF/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.69E+06)**(-1.0)+ 0.6300 * ((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.12E+03)**(-1.0) - 0.26521E-001 * ((unitL(d)/3.28084)/(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/8.0)**(-2.0) + 0.5056E-002 * ((unitL(d)/3.28084)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/4.4)**(-2.0) + 0.8488E-002 * (gasIn(d)*3600/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.23E+04)**(-2.0) + 0.819E-002 * (sorbentF/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.50E+06)**(-2.0) - 0.12015E-001 * ((solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.4)**(-2.0) + 0.1531E-001 * ((unitD(d)/3.28084)/unitDx(d)/16.)**0.50 + 0.145 * (unitDLX(d)/(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1)))**0.50 + 1.165 * (unitDLX(d)/(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.67E-02) - 0.228E-001 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/sorbentF/0.11E-05) + 0.25573E-003 * ((unitD(d)/3.28084)/unitDLX(d)/16.)**2.0 + 0.2184E-003 * ((unitL(d)/3.28084)/unitDLX(d)/8.0)**2.0 + 0.6542E-002 * ((unitL(d)/3.28084)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/4.4)**2.0 + 0.18670 * (unitDx(d)/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**2.0 - 0.4373E-001 * ((solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**2.0+  P('BFBDes8',d)  -  N('BFBDes8',d))*x(d,'BOF')+(- 0.2893 * log((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03) - 0.1019E-005 * (unitDx(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**(-2.0) - 0.297E-003 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))*sorbentF/0.90E+06)**(-2.0) + 0.457 * (unitDLX(d)*(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**0.50 - 1.6813 * ((unitL(d)/3.28084)*unitDx(d)/8.0) - 0.586E-001 * ((unitD(d)/3.28084)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/29.)**2.0 - 1.3087 * (unitDLX(d)*(gasIn(d)*3600)/0.30E+04)**2.0 - 4.823 * (unitDLX(d)*(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**2.0 + 5.4446 * (unitDLX(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**2.0 + 3.105297 * (unitDLX(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**2.0 - 0.51150 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**2.0 + 0.5773E-001 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.77)**0.50 + 0.217661 * ((unitD(d)/3.28084)/unitDLX(d)/16.)**(-1.0) + 0.867255E-001 * ((unitD(d)/3.28084)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/12.)**(-1.0) - 0.84131E-001 * ((solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.4)**(-1.0) - 0.49392E-001 * ((unitL(d)/3.28084)/(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/8.0)**(-2.0) + 0.93E-002 * ((unitL(d)/3.28084)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/4.4)**(-2.0) + 0.741E-002 * ((unitL(d)/3.28084)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/6.2)**(-2.0) + 0.283E-004 * (unitDx(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**(-2.0) - 0.1638E-001 * (sorbentF/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.69E+06)**(-2.0) + 0.4275 * (unitDLX(d)/(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1)))**0.50 - 0.9146E-001 * ((unitL(d)/3.28084)/sorbentF/0.89E-05) + 0.137E-001 * unitDLX(d)/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))) - 0.731E-001 * (unitDLX(d)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.77) + 0.2170E-002 * ((unitL(d)/3.28084)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/6.2)**2.0 - 1.255* (unitDx(d)/(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1)))**2.0 + 0.3278E-001 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.77)**2.0+  P('BFBDes8',d)  -  N('BFBDes8',d))*x(d,'BUF')+(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d)) + solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))*(1-y(d))
        ;

BFBDes9(d) $des('9') .. solidOutC(d,'HCO3') =E=
        (- 0.3352 * log((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03) + 0.346 * exp(unitDLX(d)) - 0.474E-001 * ((solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/2.3)**0.50 - 26.67 * (unitDx(d)*unitDLX(d))**2.0 + 0.108 * ((solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/2.3)**2.0 + 0.2215 * ((unitD(d)/3.28084)/unitDLX(d)/16.)**0.50 + 0.9713E-002 * (gasIn(d)*3600/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.17E+04)**(-1.0) + 0.145E-002 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.67E-02)**(-1.0) - 0.1160 * ((unitL(d)/3.28084)/sorbentF/0.89E-05)**0.50 - 0.4740 * (unitDx(d)/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**0.50 + 0.360 * unitDx(d)/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))) + 0.693E-002 * ((unitD(d)/3.28084)/(gasIn(d)*3600)/0.53E-02)**2.0 - 0.327 * (unitDLX(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**2.0 +  P('BFBDes9',d)  -  N('BFBDes9',d))*x(d,'BOF')+(0.6977E-001 * log(sorbentF/0.90E+06) - 0.9808 * (gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))**4.0 - 0.7735E-001 * ((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**4.0 + 0.54127E-002 * ((unitD(d)/3.28084)*(unitL(d)/3.28084)/0.13E+03)**(-1.0) + 0.844867 * (unitDLX(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**0.50 - 0.3273E-001 * ((unitD(d)/3.28084)*(gasIn(d)*3600)/0.48E+05) - 2.1116 * (unitL(d)/3.28084*unitDx(d)/8.0) - 22.78 * (unitDx(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**2.0 + 0.976E-001 * ((unitD(d)/3.28084)/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/16.)**0.50 + 0.165E-001 * ((unitD(d)/3.28084)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/12.)**(-2.0) + 0.370E-002 * ((unitL(d)/3.28084)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/4.4)**(-2.0) + 0.110E-004 * (unitDx(d)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.77)**(-2.0) - 0.1720E-002 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**(-2.0) + 0.1618E-002 * (sorbentF/(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.60E+04)**(-2.0) + 0.3186E-001 * (unitDx(d)/unitDLX(d))**2.0 - 0.6466E-001 * (unitDLX(d)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.77)**2.0 +  P('BFBDes9',d)  -  N('BFBDes9',d))*x(d,'BUF')+(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d)) + solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))*(1-y(d))
       ;

BFBDes10(d) $des('10') .. solidOutC(d,'NH2COO') =E=
      (- 0.2211 * log(gasOutV(d,'P')/1.3) + 0.3733E-002 * ((solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**(-4.0) + 0.790E-002 * ((unitD(d)/3.28084)*unitDx(d)/16.)**0.50  - 0.189E-001 * (unitDLX(d)*sorbentF/0.90E+06)**0.50 + 0.571E-002 * (gasIn(d)*3600*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.39E+04)**(-1.0)- 0.262E-001 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))*sorbentF/0.90E+06)**(-1.0) + 0.5799E-003 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))*sorbentF/0.90E+06)**(-2.0) + 0.869 * (gasOutV(d,'P')*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/2.3)**0.50 - 0.6588E-001 * unitDLX(d)*(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1)) + 1.0949 * (unitDLX(d)*gasIn(d)*3600/0.30E+04)**2.0 + 0.5516 * (gasIn(d)*3600*(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/0.30E+04)**2.0 - 0.200 * (gasIn(d)*3600*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.54E+04)**2.0 + 0.709 * ((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/83.)**0.50 - 0.1713E-003 * (unitDLX(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**(-2.0) + 0.278E-002 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**(-2.0) + 0.905E-002 * (sorbentF/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.69E+06)**(-2.0) - 0.343E-001 * ((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.12E+03)**(-2.0) - 0.6649E-001 * ((unitL(d)/3.28084)/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/8.0)**0.50 - 0.2672E-001 * (sorbentF/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.50E+06)**2.0+  P('BFBDes10',d)  -  N('BFBDes10',d))*x(d,'BOF')+(- 0.2573473 * log((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03) + 0.6300 * exp((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))) + 0.1261E-002 * (unitDLX(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**0.50 - 0.133E-002 * ((solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**(-1.0) - 0.7563E-003 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))*gasOutV(d,'P')/1.3)**(-2.0) + 0.5010 * (unitDLX(d)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**0.50 + 0.2487* ((unitD(d)/3.28084)*(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/16.) - 38.317 * (unitDx(d)*(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1)))**2.0 - 38.317 * (unitDx(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**2.0 + 0.765 * (unitDLX(d)*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**2.0 + 0.11274 * (gasOutV(d,'P')*(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/1.3)**2.0 + 0.2767 * ((unitD(d)/3.28084)/sorbentF/0.18E-04)**(-1.0) - 0.682E-001 * ((unitD(d)/3.28084)/sorbentF/0.18E-04)**(-2.0) + 0.2318E-001 * ((unitD(d)/3.28084)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/8.9)**(-2.0) + 0.539E-005 * (unitDx(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**(-2.0) - 0.354E-003 * (unitDLX(d)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**(-2.0) - 0.2945 * (unitDx(d)/unitDLX(d))**0.50 - 0.2104 * ((unitL(d)/3.28084)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/4.4) + 0.6166E-001 * ((unitL(d)/3.28084)/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/4.4)**2.0 + 0.2543 * (unitDx(d)/unitDLX(d))**2.0 + 1.0333 * (unitDx(d)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.77)**2.0 - 0.730E-003 * ((gasIn(d)*3600)/(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/0.30E+04)**2.0 - 0.1365E-001 * (sorbentF/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.50E+06)**2.0+  P('BFBDes10',d)  -  N('BFBDes10',d))*x(d,'BUF')+(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d)) + solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))*(1-y(d))
        ;

BFBDes11(d) $des('11') .. gasOutX(d,'CO2') =E=

          (- 0.12665E-002 * (unitDx(d)*(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1)))**0.50 + 0.4524E-002 * (gasIn(d)*3600*(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/0.30E+04)**(-1.0) + 0.53E-006 * (unitDx(d)*gasOutV(d,'P')/1.3)**(-2.0) + 0.946E-005 * (unitDLX(d)*gasIn(d)*3600/0.30E+04)**(-2.0) - 0.129E-002 * (gasIn(d)*3600*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.54E+04)**(-2.0) - 0.55E-003 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))*(solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/0.15E+03)**(-2.0) + 0.233E-001 * ((unitL(d)/3.28084)*(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/8.0)**0.50- 0.544E-002 * (unitDx(d)*gasIn(d)*3600/0.30E+04)**0.50 - 0.5471E-001 * (unitDLX(d)*gasIn(d)*3600/0.30E+04)**0.50 + 0.871 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**0.50 - 0.2418 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.3)**0.50 + 0.460E-001 * ((unitD(d)/3.28084)*(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/16.) - 0.1289302E-001 * ((unitD(d)/3.28084)*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/21.) + 0.632E-001 * (sorbentF*(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/0.90E+06) - 0.47908E-001 * (gasIn(d)*3600*(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)))/0.30E+04)**2.0 - 0.21249E-001 * ((solidRichT$(ORD(d) eq CARD(d))+solidOutT(d+1)$(ORD(d) ne CARD(d)))/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/83.)**0.50 + 0.3141E-001 * ((unitD(d)/3.28084)/(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/16.)**(-2.0) - 0.393E-002 * ((unitD(d)/3.28084)/gasOutV(d,'P')/12.)**(-2.0) - 0.379E-003 * (unitDLX(d)/(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1)))**(-2.0) - 0.365E-002 * (gasIn(d)*3600/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.23E+04)**(-2.0) + 0.896E-003 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**(-2.0) + 0.620E-003 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**(-2.0) - 0.922E-002 * (sorbentF/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.50E+06)**(-2.0) + 0.5686E-002 * (sorbentF/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.69E+06)**(-2.0) + 0.7525E-002 * ((unitL(d)/3.28084)/(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/8.0)**0.50 - 0.982E-002 * ((unitL(d)/3.28084)/sorbentF/0.89E-05)**0.50  + 0.6044E-003 * ((unitL(d)/3.28084)/unitDLX(d)/8.0) - 0.85E-002 * (sorbentF/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.50E+06) + 0.950E-002 * (gasIn(d)*3600/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.17E+04)**2.0 - 0.1450E-002 * (gasIn(d)*3600/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.23E+04)**2.0 + 0.399E-001 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/sorbentF/0.11E-05)**2.0 + 0.113 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**2.0+ P('BFBDes11',d)  -  N('BFBDes11',d))*x(d,'BOF')+(0.18832E-001 * (gasIn(d)*3600*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.39E+04)**(-1.0) + 0.22E-003 * (gasIn(d)*3600*(gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/0.30E+04)**(-2.0) - 0.49E-002 * (gasIn(d)*3600*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.54E+04)**(-2.0) - 0.18E-003 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))*sorbentF/0.90E+06)**(-2.0) + 0.40E-002 * (sorbentF*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.16E+07)**(-2.0) - 0.20E-002 * (sorbentF*(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/0.12E+07)**(-2.0) + 0.73 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))*(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/1.8)**0.50 - 0.108E-001 * ((unitL(d)/3.28084)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/6.2)**(-1.0) - 0.109E-001 * ((unitD(d)/3.28084)/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/12.)**(-2.0) + 0.98E-003 * ((gasOutX(d-1,'CO2')$(ORD(d) ne 1)+((feedCO2F*feedCO2C('CO2')+steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1))/(solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/0.56)**(-2.0) - 0.588E-001 * (unitDLX(d)/(solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d))+solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))))**0.50 + 0.11E-001 * ((solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d))+solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)))/(solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d))+solidOutC('a1','H2O')$(ORD(d) eq CARD(d)))/1.4)+ P('BFBDes11',d)  -  N('BFBDes11',d))*x(d,'BUF')+(((feedCO2F*feedCO2C('CO2') + steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1)+gasOutX(d-1,'CO2')$(ORD(d) ne 1))*(1-y(d))
        ;


MODEL  ProSyn  / ALL /;
MaxExecError = 10000000;
option sys12 = 1;
option domlim=100;

OPTIONS
  LIMROW = 1000, LIMCOL = 0
* PROFILE = 0
  SOLPRINT = ON;
*ProSyn.scaleopt=1;

scale(kinf) = 1;
scale('BFBads1') = 6000;
scale('BFBads2') = 0.29;
scale('BFBads3') = 1.33;
scale('BFBads4') = 1.705;

scale('BFBads6') = 35.75;
scale('BFBads7') = 67.89;
scale('BFBads8') = 0.83;
scale('BFBads9') = 0.29;
scale('BFBads10') = 1.09;
scale('BFBads11') = 0.04;
scale('BFBads12') = 0.07;

scale('BFBDes1') = 6000;
scale('BFBDes2') = 0.16;
scale('BFBDes3') = 1.28;
scale('BFBDes4') = 0.7;

scale('BFBDes6') = 164.46;
scale('BFBDes7') = 147.56;
scale('BFBDes8') = 0.42;
scale('BFBDes9') = 0.22;
scale('BFBDes10') = 0.95;
scale('BFBDes11') = 0.3;

gasIn.lo(a) = 1.11 ;
gasIn.up(a) = 2.5 ;
gasIn.scale(a) = 1/360;
gasIn.lo(d) = 0.26 ;
gasIn.up(d) = 0.85 ;
gasIn.scale(d) = 1/3600;
gasIn.l(d) = gasIn.lo(d);
gasOut.lo(a) = 1;
gasOut.up(a) = 2.5;

gasOut.lo(d) = 0.2 ;
gasOut.up(d) = 1;

gasInV.lo(a,'T') = 40;
gasInV.up(a,'T') = 100;
gasInV.lo(d,'T') = 140;
gasInV.up(d,'T') = 170;
gasInV.l(d,'T')= gasInV.lo(d,'T');
gasInV.lo(s,'P') = 1 ;
gasInV.up(s,'P') = 2.1 ;
gasOutV.lo(a,'P') = 1 ;
gasOutV.up(a,'P') = 1.4;
gasOutV.lo(a,'T') = 40;
gasOutV.up(a,'T') = 100;
gasOutV.lo(d,'P') = 1 ;
gasOutV.up(d,'P') = 1.3;
gasOutV.lo(d,'T') = 140;
gasOutV.up(d,'T') = 170;

feedCO2F.lo = 0.2 ;
feedCO2F.up = 0.42;
feedCO2F.l = feedCO2F.lo;
steamF.lo = 0.2;
steamF.up = 0.42;
steamF.l = steamF.lo;
hotInF.lo(d) = 0.3;
hotInF.up(d) = 0.4;

deltaPa.lo=0.01;
deltaPa.up=1;
deltaPd.lo=0.01;
deltaPd.up=1;

flueOut.lo = 1;
flueOut.up = 2.5;
unitD.lo(a) = 32;
unitD.up(a) = 60;
unitD.lo(d) = 26;
unitD.up(d) = 53;
//sorbentF.lo = 300000;
sorbentF.lo = 400000;
sorbentF.up = 900000;
//sorbentF.scale = 5e5 ;
//unitDx.lo(s) = 0.01;
//unitDx.up(s) = 0.04;
unitDx.lo(a) = 0.0175;
unitDx.up(a) = 0.03;
unitDx.lo(d) = 0.014;
unitDx.up(d) = 0.026;
unitL.lo(s) = 6.56168;
unitL.up(s) = 26.24672;
unitDLX.lo(a) = 0.05;
unitDLX.up(a) = 0.55;
unitDLX.lo(d) = 0.05;
unitDLX.up(d) = 0.25;
unitDLX.l(d) = unitDLX.lo(d);
solidOutT.lo(a)= 50;
solidOutT.up(a)= 110;
solidOutT.lo(d)= 120;
solidOutT.up(d)= 150;
solidOutT.l(d) = solidOutT.lo(d);
solidLeanT.lo = 110;
solidLeanT.up = 130;
solidRichT.lo = 90;
solidRichT.up = 120;
solidRichT.l = solidRichT.lo;
utilOutT.lo= 40;
utilOutT.up= 60;
utilInF.lo = 5;
utilInF.up = 12;
flueOutV.lo('T') = utilInT+5;
flueOutV.up('T') = fgV('T');
flueHXArea.lo = 85000 ;
flueHXArea.up = 100000 ;
// flueHXArea.scale = 1e4 ;
flueHXArea.l = flueHXArea.lo;
richHXArea.lo = 1000 ;
richHXArea.up = 100000 ;
// richHXArea.scale = 1e3;
leanHXArea.lo = 1000 ;
leanHXArea.up = 100000 ;
// leanHXArea.scale = 1e3;

coldOutT.lo(a) = 30;
coldOutT.up(a) = 50;
hotOutT.lo(d) = 150;
hotOutT.up(d) = 180;
coldInF.lo(a) = 1.5;
coldInF.up(a) = 100;


gasOutX.lo(a,'CO2') = 0.001;
gasOutX.up(a,'CO2') = 0.14;
gasOutX.lo(a,'H2O')= 0.02;
gasOutX.up(a,'H2O')= 0.14;
gasOutX.lo(d,'CO2') = 0.1;
gasOutX.up(d,'CO2') = 0.4;
gasOutX.l(d, 'CO2') = gasOutX.lo(d,'CO2');
gasInX.lo(d,'CO2') = 0.1;
gasInX.up(d,'CO2') = 0.4;

flueOutC.lo('CO2') = 0.01;
flueOutC.up('CO2') = 0.5;
flueOutC.lo('H2O') = 0.01;
flueOutC.up('H2O') = 0.5;

solidOutC.lo(a,'H2O') =0.2;
solidOutC.up(a,'H2O') =1.2;
solidOutC.lo(a,'HCO3') =0.001;
solidOutC.up(a,'HCO3') =0.6;
solidOutC.lo(a,'NH2COO') =0.2;
solidOutC.up(a,'NH2COO') =1.4;
solidOutC.lo(d,'H2O') =0.3;
solidOutC.up(d,'H2O') =1.3;
solidOutC.lo(d,'HCO3') =0.1;
solidOutC.up(d,'HCO3') =0.5;
solidOutC.lo(d,'NH2COO') =0.8;
solidOutC.up(d,'NH2COO') =1.4;

CaptureTarget.lo = 0.90;
CaptureTarget.up = 0.95;


capEX.lo = 1e8;
capEX.up = 1e10;
capEX.scale = 1e8;
COE.up = 1000;
derate.up = nameplate;
Nu.lo = 12;
Nu.up = 16;
Nu.l = Nu.lo;

unitPi.lo(s) = 14.696;
unitPi.up(s) = 45 ;
unitPo.lo(s) = 14.696;
unitPo.up(s) = 45 ;

unitNx.lo(s) = 600;
unitNx.up(s) = 12000;
//unitNx.scale(s) = 3000 ;
unitHXArea.lo(s) = 6000 ;
unitHXArea.up(s) = 170000 ;
unitHXArea.scale(s) = 6000 ;
unitW.lo(s) = 50000 ;
unitW.up(s) = 500000 ;
unitW.scale(s) = 2e5;
unitCpa.lo(a) = 100000 ;
unitCpa.up(a) = 5000000;
unitCpa.scale(a) = 1e5 ;
unitCpd.lo(d) = 100000 ;
unitCpd.up(d) = 5000000;
unitCpd.scale(d) = 1e5 ;
unitPb.lo(s) = 50 ;
//unitPb.up(s) = 1600;
unitPe.lo(s) = 10 ;
unitPe.up(s) = 500 ;

econ(i)=1; ads(i)=1; des(i)=1;

option reslim=720000;
option minlp=dicopt;
SOLVE ProSyn USING MINLP MINIMIZING fplus;
display y.l,Nu.l, CaptureTarget.l;
