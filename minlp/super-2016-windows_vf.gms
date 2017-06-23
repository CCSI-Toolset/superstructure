$Title Optimal Synthesis of Sorbent Carbon Capture Processes

$ONTEXT
================================================================================
Base version: Murthy Konda, NETL (May 19, 2011)
Modifications:
 - Miguel Zamarripa, NETL (August, 2015 - date)
 - YoungJung Chang, CMU (August, 2011)
 - Nick Sahinidis, CMU (September, 2011; December, 2011)
 - Murthy Konda, NETL (December, 2011)
 - Zhihong Yuan, CMU (March, 2012 - September, 2014)

This model identifies an optimal process configuration and corresponding optimal technology and optimal
design/operation levels for bubbling fluidized bed based CO2 capture process
*** New changes for June 2016 release
- Update of the algebraic models with the last version of the BFB model
- Using 2000 samples for each technology the models are fitted with ALAMO.
- Once the structure of the surrogate models is obtained
- The coefficients are optimized to improve the R2 metric (not used anymore)
- Changes in the superstructure code
         * Nu - Number of parallel corresponds for the adsorbers
         * the regenerators are intalled for half of Nu (so the cost are divided by 2 - only regenerators -)
         * Then the Sorbent Flow in the regenerator units is 2 times the sorbent flow in the adsorbers
- COE calculations remains unchanged
- new capital cost for the units is listed at the end of the document and calculated to compare with the previous costs
         * The next release must optimize based on the new costs calculations
- upper and lower bounds for the slack variables have been updated

CASE STUDY NOTES -
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
$OFFLISTING

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
  steamV(v)     P and T of the steam to regenerators       / P = 6.8, T = 170/   //changed T = 180
  feedCO2C(fc)  Mole fractions of feedCO2                  / CO2 = 0.42, H2O = 0.58 /
  feedCO2V(v)   P and T of the CO2 stream to regenerators  / P = 1, T = 138.07 /
  utilInC(fc)   Mole fractions of coolant water to flueHX  / CO2 = 0, N2 = 0, H2O = 1 /
  utilInT       T of the coolant water to flueHX           / 32.22 /  //32.22
  coldInC(fc)   Mole fractions of coolant water  ( =E= TO UTILINC(FC) )
  coldInT       T of the coolant water to adsorber   (=E= TO COLDINT)
  hotInC(fc)    Mole fractions of heating steam      (=E= STEAMC(FC))
  hotInT        T of the heating steam to desorber
  coolInT       T of the solid cooling water
  warmInT       T of the solid heating steam
  unitNp(s)                                                /a1=1,a2=1,a3=1,a4=1,d1=1,d2=1,d3=1,d4=1/

*new parameters Miguel (12-08-2015)
 cbt Bubble Region Gas Total Concentration kmol per m3 Value: 0.04 /0.04/


;

flueInV(v) = fgV(v);
coldInC(fc) = utilInC(fc);       coldInT = utilInT;
hotInC(fc) = steamC(fc);         hotInT = steamV('T');
coolInT = utilInT;
warmInT = steamV('T');


positive variable  COE  Adder to the cost of electricity due to capture ($ per MWh);
variable  soi  sum of constraint infeasibilities;
variable fplus objective;


binary VARIABLES
  x(s,t)                 1 iff technology t is used for stage s
  y(s)                   1 iff stage s is ever utilized;


integer VARIABLES
  Nu                     Number of parallel trains (even number)
  Nup                    Number of parallel trains (used to make Nu an even)

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
  sorbentFDes            Sorbent flow rate (kg per hr input for regenerators)  __Miguel 02_05_2016
  flueInC(fc), flueOutC(fc), flueOut
  utilInF                Coolant water for flue HX
  gasInC(s,fc)           Flue gas molar composition at the inlet of stage s
  gasOutC(s,fc)          Flue gas molar composition at the outlet of stage s
  gasOut(s)              flue gas molar flow at the outlet of stage s
  solidOutC(s,sc)        Solid component Loading of the solid stream outlet of stage s
  steamF                 Steam flow to the overall sytem (desorbers)
  feedCO2F
  coldInF(a)            design variable (not used in the model)must be cold inlet flow
  hotInF(d)             Heating steam to desorber d(s)
* PRESSURE/TEMPERATURE Variables (Pressure in atm, Temperature in C)
  flueOutV(v)
  utilOutT               Temperature of colant water for HX
  gasInV(s,v), gasOutV(s,v)
  solidOutT(s)           Temperature of Sorbent (solids) coming out from s
  pureCO2V(v)
  coldOutT(a)            Temperature for colant water adsorbers
  hotOutT(d)             Temperature steam out of the desorber d(s)
  deltaPa                The pressure increase through the compressor in adsorber
  deltaPd                The pressure increase through the compressor in regenerator
  solidLeanT             Temperature after the heat exchanger (solidLeanHX)
  solidRichT             Temperature after the heat exchanger (solidRichHX)
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
  unitLb(s) depth of the solid bed in regenerator d
  unitAx(d) area of the vessels
  unitvg(d) fluid velocity
* SLACK VARIABLES
  P(*,s), N(*,s), P1(*), N1(*), P2(a,fc), N2(a,fc), P3(s,fc), N3(s,fc), P4(s,fc), N4(s,fc), P5(*), N5(*), P6(*),N6(*)
;


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

  cont1          Material balance around the first heat exchanger (flue gas HX)
  cont2


  extra1         Depth of the bed should be no larger than the diameter
  extra2         CO2 outlet of stage a must be lower than xfactor *gasInX(CO2) xfactor obtained from the data sets
  extra2b        CO2 outlet of stage a must be greater than xfactor *gasInX(CO2)
  extra3         unitL is 2 times the unitLb (from BFB model -ACM model developed by Andrew Lee)
  extra4         gasinX eq gasOutX
  extra5a        car regeneration rules similar to extra 2 and 2b (not used in this release)
  extra5b        car regeneration rules similar to extra 2 and 2b (not used in this release)
  extra6         Nu = number of adsorbers _ Nu over 2 = number of regenerators then SorbentFDes = 2*SorbentF
  extra7         used to make Nu an even number
  extra8
  BFBads1        Surrogate model adsorbers Nx
  BFBads2        Surrogate model adsorbers flue Gas outlet Pressure
  BFBads3        Surrogate model adsorbers flue gas outlet flow
  BFBads4        Surrogate model adsorbers flue gas outlet temperature
  BFBads6        Surrogate model adsorbers Cold out temperature (coolant HX adsorbers)
  BFBads7        Surrogate model adsorbers solid outlet temperature
  BFBads8        Surrogate model adsorbers solid outlet load H2O
  BFBads9        Surrogate model adsorbers solid outlet load HCO3
  BFBads10       Surrogate model adsorbers solid outlet load NH2COO
  BFBads11       Surrogate model adsorbers gas outlet CO2 (fraction)
  BFBads12       Surrogate model adsorbers gas outlet H2O (fraction)
  BFBDes1        Surrogate model regenerator Nx
  BFBDes2        Surrogate model regenerator gas outlet pressure
  BFBDes3        Surrogate model regenerator gas outlet flow
  BFBDes4        Surrogate model regenerator gas outlet temperature
  BFBDes6        Surrogate model regenerator Hot outlet temperature (steam flow)
  BFBDes7        Surrogate model regenerator solid outlet temperature
  BFBDes8        Surrogate model regenerator solid outlet load H2O
  BFBDes9        Surrogate model regenerator solid outlet load HCO3
  BFBDes10       Surrogate model regenerator solid Outlet load NH2COO
  BFBDes11       Surrogate model regenerator gas outlet CO2 (fraction)
  gasOutConvads(a,fc)    Material balances gas inlet and outlet (adsorbers)
  gasOutConvdes(d,fc)    Material balances gas inlet and outlet (regenerators)
  gasInadConv(a,fc)      Material balances around single units (adsorbers)
  gasIndeConv(d,fc)      Material balances around single units (regenerators)
  gasInCS(s,fc)          Material balances per components (adsorbers and regenerators)
  gasOutSum(a)           Material balances around the system by component(adsorbers)
  gasInSum(a)            Material balances around the system by component (adsorbers)
  gasOutSumd(d)          Material balances around the system by component (regenerators)
  gasInSumd(d)           Material balances around the system by component (regenerators)
  flueOutSum             Material balances flue gas (system)
;

*ECONOMIC MODULE

objective..

  CF*365*24*(nameplate/7446-derate/7446)*(COE-7.53) =G= (0.124*(1316.19E6+capEX)+(1.25*(6.80E6+0.0053343*(1316.19E6+capEX))+0.0163*(1316.19E6+capEX))+0.008*(1316.19E6+capEX)+14.41E6+70.08E6)/7446;

SumOfInfeasibilities ..
  soi =G= sum((kinf, a), (P(kinf,a)+N(kinf,a))/scale(kinf))+sum((kinf, d), (P(kinf,d)+N(kinf,d))/scale(kinf))
        + P1('flueOut') + N1('flueOut')
        + sum((a,fc), P2(a,fc)+N2(a,fc))
        + sum((s,fc), P3(s,fc)+P4(s,fc)+N3(s,fc)+N4(s,fc))
         ;

** ORIGINAL combo .. fplus =e= 0.2*COE+ 1*soi;
combo .. fplus =e= 0.2*COE+ 100*soi;

eqDerate $econ('1') ..
  100*derate =G= 100*0.924*steamFlow //420.42*steamFlow;
  + 100*7.46E-4*SUM(s,Nu*y(s)*unitPb(s)+Nu*y(s)*unitPe(s));

eqSteam $econ('2') ..
  steamFlow =G= 18*Nu/2*(steamF+SUM(d,y(d)*hotInF(d)));
                   // regenerator trains are half of the adsorbers trains
eqNu(fc) $econ('3') ..
  flueInC(fc) =E= fgF*fgC(fc)/Nu;

eqCapEX $econ('4') ..                                           // the unit cost is half of the trains for adsorbers
   (capEX/(OCF))/1e5 =G= SUM(a,Nu*y(a)*unitCpa(a)/1e5) + SUM(d,Nu/2*y(d)*unitCpd(d)/1e5)
  + 1.10 *Nu* (EXP(11.9052-0.8709*LOG(richHXArea)+0.09005*LOG(richHXArea)**2-LOG(1e5)))
  + 1.10 *Nu* (EXP(11.9052-0.8709*LOG(leanHXArea)+0.09005*LOG(leanHXArea)**2-LOG(1e5)))
  //+ 1.03*1.12*4.25*EXP(11.9052-0.8709*LOG(flueHXArea)+0.09005*LOG(flueHXArea)*LOG(flueHXArea));
  + 1.12*0.98*Nu*((7*EXP(11.9052-0.8709*LOG(12000)+0.09005*LOG(12000)*LOG(12000)-LOG(1e5))+EXP(11.9052-0.8709*LOG(flueHXArea-7*12000)+0.09005*LOG(flueHXArea-7*12000)*LOG(flueHXArea-7*12000)-LOG(1e5))));

eqUnitCpa(a) $econ('5') ..
  unitCpa(a)/1e5 =G= fm*EXP(7.0132+0.18255*LOG(unitW(a))+0.02297*LOG(unitW(a))**2-LOG(1e5))   // Cost of vessels
    + EXP(7.59176+0.7932*LOG(unitPb(a))-0.0129*LOG(unitPb(a))**2-LOG(1e5))$(ORD(a) eq 1)+0* EXP(7.59176+0.7932*LOG(unitPb(a))-0.0129*LOG(unitPb(a))**2-LOG(1e5))$(ORD(a) ne 1)                // Blower (compressor) cost
    + 1.10 * (EXP(11.9052-0.8709*LOG(unitHXArea(a))+0.09005*LOG(unitHXArea(a))**2-LOG(1e5)))      // Cost of in-reactor heat exchanger
    + unitNp(a)*734.96*1.87*(0.3048*unitD(a))**2.0049/1e5                                // Cost of plates =?= unitNp(s,t)*468.0*EXP(0.1739*unitD(s))
    + 361.8*(unitD(a)**0.7396)*(unitL(a)**0.70684)/1e5 //300.9*(unitD(s)**0.63316)*(unitL(s)**0.80161)  // Cost of platforms and ladders
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
     3600*unitPe(s) =G= 0.020*2.2046*(sorbentF$(ord(s) lt 5) + sorbentFDes$(ord(s) gt 4))*unitL(s)**0.63 + 0.00182*2.2046*(sorbentF$(ord(s) lt 5) + sorbentFDes$(ord(s) gt 4))*unitL(s);

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

  (1$SAMEAS(fc,'H2O')+1$(NOT SAMEAS(fc,'H2O')))*flueOutC(fc)*100 =E= 100*(flueInC(fc)$(NOT SAMEAS(fc,'H2O'))+0.58*flueInC(fc)$SAMEAS(fc,'H2O'));


* HYDRODYNAMIC/ENERGY BALANCES

eqEnrgAds0g(v)..  // flash calculation
   flueOutV(v) =E= flueInV(v)$SAMEAS(v,'P')+(5.9405*rpower(flueInV(v),1/2))$SAMEAS(v,'T');

***FLUE GAS HEAT EXCHANGER EQUATION
eqEnrgAds0u..  // HX model assuming counter-current operation
   utilOutT*utilInF  =E=
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

*INFEASIBLE BY 0.00000004
variable slackz2; slackz2.lo = -0.001; slackz2.up=0.001;
eqRegenerationbalance..
  1000*flueInC('CO2')+solidOutC('d1','NH2COO')*sorbentFDes+solidOutC('d1','HCO3')*sorbentFDes
*          =E= 1000*gasOutC('a%noADS%','CO2') + solidOutC('a1','NH2COO')*sorbentFDes + solidOutC('a1','HCO3')*sorbentFDes;
        =E= 1000*gasOutC('a%noADS%','CO2') + solidOutC('a1','NH2COO')*sorbentF + solidOutC('a1','HCO3')*sorbentF
                                           + solidOutC('d1','NH2COO')*(sorbentFDes-sorbentF) + solidOutC('d1','HCO3')*(sorbentFDes-sorbentF)  ;
eqAdsorbebalance..
  1000*feedCO2F*feedCO2C('CO2')+solidOutC('a1','NH2COO')*sorbentF+solidOutC('a1','HCO3')*sorbentF
         =E= 1000*gasOutC('d%noADS%','CO2')+solidOutC('d1','NH2COO')*sorbentF+solidOutC('d1','HCO3')*sorbentF ;


* Additional process related conditions must be also added
* For example, the input variable bounds for each technology and its model

gasOutConvads(a,fc)..   gasOutC(a,fc) =E= gasOutX(a,fc)*gasOut(a);
gasOutConvdes(d,fc)..   gasOutC(d,fc) =E= gasOutX(d,fc)*gasOut(d);

gasInadConv(a,fc)..    gasInC(a,fc) =E= gasOutC(a-1,fc)$(ORD(a) ne 1)+flueOutC(fc)$(ORD(a) eq 1);
gasIndeConv(d,fc)..    gasInC(d,fc) =E= (feedCO2F*feedCO2C(fc) + steamF*steamC(fc))$(ORD(d) eq 1)+gasOutX(d-1,fc)*gasIn(d)$(ORD(d) ne 1);

gasOutSum(a)..         SUM(fc,gasOutX(a,fc)) =E= 1;
gasInSum(a)..          SUM(fc,gasInX(a,fc)) =E= 1;
gasInSumd(d)..         gasInX(d,'CO2')+gasInX(d,'H2O')  =E= 1;
gasOutSumd(d)..        gasOutX(d,'CO2')+gasOutX(d,'H2O')  =E= 1;

variable slackz slack variable similar than P and N;
slackz.lo = -1; slackz.up = 1;

flueOutSum..           flueOut =E= Sum(fc,flueOutC(fc)) + P1('flueOut')- N1('flueOut');

cont1(a) ..            gasIn(a) =e= gasOut(a-1) $(ORD(a) ne 1) + flueOut $(ORD(a) eq 1);
cont2(d) ..            gasIn(d) =E= (feedCO2F+steamF)$(ORD(d) eq 1) + gasOut(d-1)$(ORD(d) ne 1);


gasInCS(s,fc)..        100*gasInC(s,fc) =E= 100*gasInX(s,fc)*gasIn(s)+ P4(s,fc)-N4(s,fc);
P4.lo(a,fc) = 0; P4.up(a,fc) = 0.1; N4.lo(a,fc) = 0; N4.up(a,fc)=0.1;  P4.l(a,fc) = 0; N4.l(a,fc)=0;

extra1(a) ..           unitL(a) =L= unitD(a);

extra2(a)..            gasOutX(a,'CO2') =l= gasInX(a,'CO2')*1*(1-y(a)) + gasInX(a,'CO2')*0.85*(y(a));
extra2b(a)..           gasOutX(a,'CO2') =g= gasInX(a,'CO2')*1*(1-y(a)) + gasInX(a,'CO2')*0.70*(y(a));

extra3(a) ..           unitL(a) =e= 2*unitLb(a)*3.24084 + slackz ;

extra4(a,fc)$(ord(a) gt 1)..     gasInX(a,fc) =e= gasOutX(a-1,fc);

** from data samples solidOut Car is 2 times the solidsIn for CO2 in the adsorbers (not used in this release version)
extra5a(a)..     SolidOutC(a,'NH2COO') =e=
                 (2*SolidOutC('d1','NH2COO')$(ord(a) eq card(a)) + 2*SolidOutC(a+1,'NH2COO')$(ord(a) ne card(a)))*(y(a))
                         + (1*SolidOutC('d1','NH2COO')$(ord(a) eq card(a)) + 1*SolidOutC(a+1,'NH2COO')$(ord(a) ne card(a)))*(1-y(a));
extra5b(a)..     SolidOutC(a,'NH2COO') =l= (2.1*SolidOutC('d1','NH2COO')$(ord(a) eq card(a)) + 2.1*SolidOutC(a,'NH2COO')$(ord(a) ne card(a)))$(y(a));

extra6..         sorbentFDes =g= 2*sorbentF ;
**** New variables to estimate the number of parallel trains using 2 integer variables
extra7..         Nu =e= 2 * Nup;
Nup.lo = 6;
Nup.up = 8;
extra8(a)..      ColdOutT(a) =g= UtilInT;

*** $include function, calls a gams document to be included in this document. (Use in the same directory or call it using the path C:users/....)
$include ADSORBER_vf.gms
$include REGENERATOR_vf.gms

MODEL  ProSyn  / ALL /;

MaxExecError = 10000000;
option sys12 = 1;
option domlim=100;

OPTIONS
  LIMROW = 1000, LIMCOL = 1000
* PROFILE = 0
  SOLPRINT = ON;
  ProSyn.scaleopt=1;

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

gasIn.lo(a) = 0.833;
gasIn.up(a) = 3.333;
gasIn.scale(a) = 1/3600;
unitvg.lo(d) =0.1;
unitvg.up(d) =1;
unitLb.lo(d) = 1;
unitLb.up(d) = 10;

gasIn.lo(d) = 0.05 ;
gasIn.up(d) = 8 ;
gasIn.scale(d) = 1/3600;
gasIn.l(d) = 3;
*gasOut.lo(a) = 1;
*gasOut.up(a) = 2.5;

gasOut.lo(d) = 0.2 ;
gasOut.up(d) = 1;

gasInV.lo(a,'T') = 35;
gasInV.up(a,'T') = 80;
gasInV.lo(d,'T') = 120;
gasInV.up(d,'T') = 170;
gasInV.l(d,'T')= gasInV.lo(d,'T');
gasInV.lo(s,'P') = 1 ;
gasInV.up(s,'P') = 2.4 ;
gasInV.l('a1','P')=1.178;
gasOutV.lo(a,'P') = 1.2 ;
gasOutV.up(a,'P') = 2.4;
gasOutV.lo(a,'T') = 40;
gasOutV.up(a,'T') = 100;
gasOutV.lo(d,'P') = 1 ;
gasOutV.up(d,'P') = 1.94;
gasOutV.lo(d,'T') = 93;
gasOutV.up(d,'T') = 168;

feedCO2F.lo = 0.2 ;
feedCO2F.up = 0.42;
feedCO2F.l = feedCO2F.lo;
steamF.lo = 0.2;
steamF.up = 0.42;
steamF.l = steamF.lo;
***kmol/hr 6.98 - 13000  -->/3600
hotInF.lo(d) = 0.3;
hotInF.up(d) = 0.4;

deltaPa.lo=0.01;
deltaPa.up=1;
deltaPa.fx=0.178;
deltaPd.lo=0.01;
deltaPd.up=1;

flueOut.lo = 1;
flueOut.up = 2.5;
unitD.lo(a) =13;
*32;
unitD.up(a) =49;
*60;
unitD.lo(d) = 13;
unitD.up(d) = 49;

//sorbentF.lo = 300000;
sorbentF.lo = 600000;
*400000;
sorbentF.up = 800000;
*900000;
sorbentFDes.lo = 1000000;
sorbentFDes.up = 2000000;
sorbentFDes.l = 1850000.0;
sorbentFDes.scale = 1E5;
sorbentF.scale = 5e5 ;
//unitDx.lo(s) = 0.01;
//unitDx.up(s) = 0.04;
unitDx.lo(a)= 0.01;
*unitvg.fx(d) = 0.22;
unitLb.l(s)=3.6;
unitLb.lo(s)=1;
unitLb.up(s)=10;
*unitD.fx(d) = 29; unitDLX.fx(d) = 0.24; unitDx.fx(a)=0.02;
*unitdx.fx(d) = 0.075;
unitDx.up(a) =0.1;

unitDx.lo(d) = 0.01;
unitDx.up(d) = 0.15;

unitL.lo(s) = 13;
unitL.up(s) = 49;
unitDLX.lo(a) =0.15;
*0.05;
unitDLX.up(a) =0.5;
*0.55;
unitDLX.lo(d) =0.015;
unitDLX.up(d) =0.5;

solidOutT.lo(a)=34;
solidOutT.up(a)= 102;
SolidOutT.l(a)=73;
*solidOutT.fx('a2')=73;
solidOutT.lo(d)= 89;
solidOutT.up(d)= 173;
solidOutT.l(d) = solidOutT.lo(d);
solidLeanT.lo = 70;
solidLeanT.up = 102;
***solidinTemperature to regenerators
solidRichT.lo =120;
*90;
solidRichT.up = 150;
solidRichT.l = solidRichT.lo;
utilOutT.lo= 30;
*40;
utilOutT.up= 60;
utilInF.lo = 5;      //kg/s?
utilInF.up = 12;
flueOutV.lo('T') = utilInT+5;
flueOutV.up('T') = fgV('T')+30;     // it was only = fgV('T') changed miguel 2/03/2016
flueHXArea.lo = 85000 ;
flueHXArea.up = 100000 ;
flueHXArea.scale = 1e4 ;
flueHXArea.l = flueHXArea.lo;
richHXArea.lo = 1000 ;
richHXArea.up = 100000 ;
richHXArea.scale = 1e3;
leanHXArea.lo = 1000 ;
leanHXArea.up = 100000 ;
leanHXArea.scale = 1e3;

coldOutT.lo(a) = 30;
coldOutT.up(a) = 65;
hotOutT.lo(d) = 89;
hotOutT.up(d) = 173;
coldInF.lo(a) = 1.5;
coldInF.up(a) = 100;


gasOutX.lo(a,'CO2') = 0.0001;
gasOutX.up(a,'CO2') = 0.16;
*gasOutX.l('a3','CO2')= 0.093;
*gasOutX.l('a4','CO2')= 0.093;
*gasOutX.fx('a3','H2O')= 0.180;
*gasOutX.fx('a4','H2O')= 0.180;
gasOutX.lo(a,'H2O')= 0.0092;
gasOutX.up(a,'H2O')= 0.18;
gasOutX.lo(d,'CO2') = 0.22;
gasOutX.up(d,'CO2') = 0.88;
gasOutX.l(d, 'CO2') = gasOutX.up(d,'CO2');
gasInX.lo(d,'CO2') = 0.5;
gasInX.up(d,'CO2') = 0.899;

flueOutC.lo('CO2') = 0.01;
flueOutC.up('CO2') = 0.5;
flueOutC.lo('H2O') = 0.01;
flueOutC.up('H2O') = 0.5;
flueOutC.lo('N2') = 0.01;
flueOutC.up('N2') = 10;

solidOutC.lo(a,'H2O') =0.155;
solidOutC.up(a,'H2O') =2.7286;
solidOutC.lo(a,'HCO3') =0.01089;
solidOutC.up(a,'HCO3') =0.58;
solidOutC.lo(a,'NH2COO') =0.434;
solidOutC.up(a,'NH2COO') =2.2293;
solidOutC.lo(d,'H2O') =0.020;
solidOutC.up(d,'H2O') =1.02;
solidOutC.lo(d,'HCO3') =0.00051;
solidOutC.up(d,'HCO3') =0.02152;
solidOutC.lo(d,'NH2COO') =0.4306;
solidOutC.up(d,'NH2COO') =1.8827;

CaptureTarget.lo = 0.90;
CaptureTarget.up = 0.98;


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

unitNx.lo(a) = 40;
*600;
unitNx.up(a) = 6240;
unitNx.up(d) = 37;
unitNx.up(d) = 43500;
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
unitPe.scale(s) = 100;
econ(i)=1; ads(i)=1; des(i)=1;

*unitFb.lo(a) =( Rp*60*0.0001*(gasInV.lo(a,'T')+273.15)  )/unitPo.lo(a);
*unitFb.up(a) =( Rp*60*0.0001*(gasInV.up(a,'T')+273.15)  )/unitPo.up(a);
unitFb.scale(s) = 1E5;
*unitF.lo(d) = (Rp*60*1*(gasInV.lo(d,'T')+273.15)  )/unitPo.lo(d);
*unitF.up(d) = (Rp*60*1*(gasInV.up(d,'T')+273.15)  )/unitPo.up(d);
unitF.scale(s) = 1E5;
UnitPb.scale(s)= 100;

$ontext
y.fx(s)=0;
x.fx(s,t)=0;
y.fx('a1')= 1.000;
y.fx('a2')= 0.000;
x.fx('a1','bof')=1;
x.fx('a2','bof')=0;
x.fx('d1','bof')=1;
y.fx('d1')= 1.000;
y.l('d2')= 0;
$offtext

y.l(s)=0;
x.l(s,t)=0;
y.l('a1')= 1.000;
y.l('a2')= 1.000;
y.l('d1')= 1.000;
y.l('d2')= 0;
x.l('a1','BOF')=1;
x.l('a2','BOF')=1;
x.l('d1','BOF')=1;


unitAx.lo(d) = unitNx.lo(d)*( (unitdlx.lo(d) + unitdx.lo(d))**2 + ((3.14159265/4)*(unitdx.lo(d))**2));
unitAx.up(d) = unitNx.up(d)*( (unitdlx.up(d) + unitdx.up(d))**2 + ((3.14159265/4)*(unitdx.up(d))**2));


$ontext
**Initial Conditions test case
GasIn.fx('a1')        =        1.219230083;
GasOut.fx('a1')        =       3631/3600;
*FlueOut.fx=1.219;
*gasInX.fx('a1','CO2')        =        0.115602897;
*gasInX.fx('a1','H2O')        =        0.131244412;
unitD.fx('a1')        =        29.04844257     ;
*unitL.fx('a1')        =        29.04844257     ;
unitLb.fx('a1')        =        6.655483315   ;
gasInV.l('a1','P')        =        1.680997632 ;
*gasInV.fx('a1','T')        =        74.1480437  ;
unitDx.fx('a1')=        0.0819   ;
unitNx.fx('a1')=        659.30  ;
unitDlx.fx('a1')=        0.270;
Sorbentf.fx=461142.3165;
SolidOutT.fx('a1')= 60;
gasOutX.lo(a,'CO2')=0.004;
gasOutX.up(a,'CO2')=0.15;
*gasOutX.lo(a,'CO2')=GasInC;
Nu.fx=12;

*gasOut.lo(a) =     1073.138/3600;
*gasOut.lo(a)=    11000.059/3600;

GasIn.l('a1')        =        1.219230083;
GasOut.l('a1')        =       3631/3600;
$offtext
ColdOutT.l('a1') = 34.15;
*FlueOut.l              = 2.187;
gasInX.l('a1','CO2')   = 0.1275;
gasInX.l('a1','H2O')   = 0.08642;
*gasOutX.fx('a1','CO2')  =0.85*0.1275;
unitD.l('a1')           = 29.04844257     ;
*unitL.fx('a1')         = 29.04844257     ;
unitLb.l('a1')          = 3.655483315   ;
*gasInV.fx('a1','P')      = 1.680997632 ;
*gasInV.fx('a1','T')    = 74.1480437  ;
unitDx.l('a1')          = 0.0819   ;
unitNx.l('a1')          = 659.30  ;
unitDlx.l('a1')         = 0.270;
Sorbentf.l              = 461142.3165;
SolidOutT.l('a1')       = 60;
solidOutC.l('a1','NH2COO') =1.1382;
solidOutC.l('a1','HCO3') =0.1333;
solidOutC.fx('a1','H2O') =0.63227;
*gasOutX.lo(a,'CO2')    = 0.004;
*gasOutX.up(a,'CO2')    = 0.15;
Nu.l                   =16;
N.up('BFBads8',a)=0.1;
*N.up('BFBads9',a)=0.1;
N.up('BFBads10',a)=0.1;
N.up('BFBdes8',a)=0.1;
*N.up('BFBdes9',a)=0.1;
N.up('BFBdes10',a)=0.1;

model prosyn2 /

SumOfInfeasibilities,
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
  BFBDes12
*  UAx
  objective
  cont1
  cont2

  extra2
* extra2b
  extra3
  extra1
  extra4
*  extra5a
*  extra5b
  extra6
  extra7
  extra8
  combo,
  eqDerate, eqSteam, eqNu
  eqCapEX
  eqUnitCpa,eqUnitCpd, eqUnitW, eqUnitPb, eqUnitPe,
  eqUnitFbAds, eqUnitFbDes, eqUnitFAds, eqUnitFDes


  eqUnitPiAds, eqUnitPoAds, eqUnitPiDes, eqUnitPoDes
  eqUnitHXArea, eqFlueHXArea
  eqEnrgAds0g
  eqStageAds, eqStageDes
  eqstage
  eqNumAds, eqNumDes


  eqMassAds0g, eqEnrgAds0u, eqEnrgAds2g, eqEnrgAds3s
  eqEnrgDes2g, eqEnrgDes3s, eqEnrgDes4g

  eqCaptureTarget
  eqRegenerationbalance
  eqAdsorbebalance


  gasOutConvads
  gasOutConvdes
  gasInadConv
  gasIndeConv
 gasInCS
  gasOutSum
  gasInSum
  gasOutSumd
  gasInSumd
  flueOutSum
$ontext
$offtext
/

option nlp=conopt;
option mip=cplex;
option reslim=720000;
option MINLP=dicopt;
option RMINLP=conopt;
option sysout=on;
*option MINLP=baron;

SOLVE ProSyn2 USING MINLP MINIMIZING fplus;
display COE.l, soi.l,y.l,Nu.l,x.l, CaptureTarget.l, derate.l, steamFlow.l,     capEX.l,      unitCpa.l,  unitCpd.l ;

* ** Display results
parameters
         par1COE, par1soi, par1y, par1CaptureTarget, par1derate,
         par1steamFlow, par1CapEx, par1unitCpa(a), par1unitCpd(d),
         par1X(s,t)            Solution if installed or not each stage s using technology t in each iteration lo
         par1Nu               Nu fixed values to evaluate in the loop /1 12, 2 14, 3 16/
         cputime1             cpu time in seconds
         solution1(*)          mod 1=optimal 2=local optimal 4=infeasible stat eq 1=normal termination 2=out iterations 10 to 13 failure
;
par1COE             = COE.l;
par1soi             = soi.l;
par1y(s)            = y.l(s);
par1CaptureTarget   = CaptureTarget.l;
par1derate          = derate.l;
par1steamFlow       = steamFlow.l;
par1CapEx       = CapEx.l;
par1unitCpa(a)      = unitCpa.l(a);
par1unitCpd(d)      = unitCpd.l(d);
par1X(s,t)          = X.l(s,t);
cputime1            = ProSyn2.resusd;
solution1('mod')    =ProSyn2.modelstat;
solution1('stat')   =ProSyn2.solvestat;

** Print  values
Execute_unload "Superstructure_results.gdx" FgF, FgC,fgV,COE.l, soi.l, y.l, CaptureTarget.l, derate.l, steamFlow.l
                                                  CapEx.l, unitCpa.l, unitCpd.l, X.l, Nu.l, cputime1, solution1
*        Gas PRoperties
         FlueInC.l, gasInC.l,FlueOutC.l, gasInV.l, FlueOutV.l,gasOutV.l,gasOutC.l,
*        Coolant and steam
         utilInF.l,ColdOutT.l, HotOutT.l,UtilInT, steamV
*        STEAMF FEEDCO2F FEEDCO2C
         steamF.l, FeedCO2F.l, FeedCO2C, FeedCO2V
*        Solids
         SolidOutC.l, SolidOutT.l,SolidRichT.l, SolidLeanT.l, sorbentF.l, sorbentFDes.l
*        HX
         FlueHXArea.l, LeanHXArea.l, RichHXArea.l, unitHXArea.l
* Unit specs
         unitD.l, unitL.l, unitPi.l, unitPo.l, unitF.l, unitFb.l, unitW.l, unitPb.l, unitPe.l, unitDx.l, unitNx.l, unitDLX.l, unitLb.l
;

*rng = R2 creates a new sheet called "R2" in excel (but, you provide the cells and results can be displayed in a specific location)
**economics
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=COE.l rng=Output!B3'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=soi.l rng=Output!B4'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=derate.l rng=Output!B5'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=CapEx.l rng=Output!B6'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=steamFlow.l rng=Output!B8'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=CaptureTarget.l rng=Output!B7'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=Nu.l rng=Output!B9'

Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitCpa.l rng=Calc!B1:E2'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitCpd.l rng=Calc!B4:E5'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=y.l rng=Calc!B7:I8'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=X.l rng=Calc!B10'

*Gas input and temperatures
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=gasInC.l rng=Calc!B20'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=FlueInC.l rng=Calc!G20'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=FlueOutC.l rng=Calc!G24'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=gasInV.l rng=Calc!K20'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=gasOutV.l rng=Calc!O20'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=gasOutC.l rng=Calc!U20'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=FlueOutV.l rng=Calc!G27'
*coolant
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=UtilinF.l rng=Calc!B31'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=ColdOutT.l rng=Calc!D31'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=HotOutT.l rng=Calc!J31'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=steamV.l rng=Calc!B34'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=UtilinT.l rng=Calc!B32'

*fluegas in
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=FgF rng=Calc!B38'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=FgC rng=Calc!D38'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=FgV rng=Calc!H38'
*steamF and FeedCO2F
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=SteamF.l rng=Calc!B42'

Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=FeedCO2F rng=Calc!B43'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=FeedCO2V rng=Calc!D43'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=FeedCO2C rng=Calc!H43'

*Solids
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=SolidOutC rng=Calc!B47'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=SolidOutT rng=Calc!G47'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=SolidRichT rng=Calc!G50'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=SolidLeanT rng=Calc!G52'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=SorbentF rng=Calc!G54'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=SorbentFDes rng=Calc!G55'

*Heat exchangers area
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=FlueHXArea rng=Calc!B59'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=LeanHXArea rng=Calc!B60'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=RichHXArea rng=Calc!B61'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitHXArea rng=Calc!B63'

*unit specs
*unitD.l, unitL.l, unitPi.l, unitPo.l, unitF.l, unitFb.l, unitW.l, unitPb.l, unitPe.l, unitDx.l, unitNx.l, unitDLX.l
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitD rng=Calc!B65'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitL rng=Calc!B67'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitPi rng=Calc!B69'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitPo rng=Calc!B71'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitF rng=Calc!B73'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitFb rng=Calc!B75'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitW rng=Calc!B77'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitPb rng=Calc!B79'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitPe rng=Calc!B81'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitDx rng=Calc!B83'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitNx rng=Calc!B85'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitDLX rng=Calc!B87'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=unitLb rng=Calc!B89'

*model stats
*Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N var=NewCostUnits.l rng=NewCostUnits!'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=cputime1.l rng=Output!B10'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=solution1.l rng=solution!'

$ontext
$offtext

********************************************************************************
************** Print results to FOQUS GamsOutput.txt file***********************
parameter
*Inlet parameters
parGasInF,parGasInP,parGasInT,parGasInCO2,parGasInH2O,parSolidInFm,parSolidInT,parSolidInP,parSolidInBic,parSolidInCar,parSolidInH2O
*Outlet parameters
parGasOutF,parGasOutP,parGasOutT,parGasOutCO2,parGasOutH2O,parSolidOutFm,parSolidOutT,parSolidOutP,parSolidOutBic,parSolidOutCar,parSolidOutH2O
parDt;

parGasInF=sum(fc,gasInC.l('a1',fc))*3600;
parGasInP=GasInV.l('a1','P');
parGasInT=GasInV.l('a1','T');
parGasInCO2=gasInC.l('a1','CO2')/(parGasInF/3600);
parGasInH2O=gasInC.l('a1','H2O')/(parGasInF/3600);
parSolidInFm=SorbentFDes.l/2;
parSolidInT=solidOutT.l('a1');
parSolidInP=1.4;
parSolidInBic=solidOutC.l('a1','HCO3');
parSolidInCar=solidOutC.l('a1','NH2COO');
parSolidInH2O=solidOutC.l('a1','H2O');

parGasOutF=sum(fc,gasInC.l('a2',fc))*3600;
parGasOutP=GasInV.l('a2','P');
parGasOutT=GasInV.l('a2','T');
parGasOutCO2=gasInC.l('a2','CO2')/(parGasOutF/3600);
parGasOutH2O=gasInC.l('a2','H2O')/(parGasOutF/3600);
parSolidOutFm=SorbentFDes.l/2;
parSolidOutT=solidOutT.l('a2');
parSolidOutP=1.4;
parSolidOutBic=solidOutC.l('a2','HCO3');
parSolidOutCar=solidOutC.l('a2','NH2COO');
parSolidOutH2O=solidOutC.l('a2','H2O');
parDt=unitD.l('a1')/3.28084;

display gasInV.l,
parGasInF,parGasInP,parGasInT,parGasInCO2,parGasInH2O,parSolidInT,parSolidInP,parSolidInBic,parSolidInCar,parSolidInH2O
parGasOutF,parGasOutP,parGasOutT,parGasOutCO2,parGasOutH2O,parSolidOutT,parSolidOutP,parSolidOutBic,parSolidOutCar,parSolidOutH2O;

*Output variables
display gasoutV.l;

file GamsOutput /GamsOutput.txt/;
put GamsOutput;
*< justified to the left :6 number of decimals to be displayed
put
*Inlet to the adsorber
parGasInF:<:6    /
parGasInP:<:6    /
parGasInT:<:6    /
parGasInCO2:<:6  /
parGasInH2O:<:6  /
parSolidInFm:<:6 /
parSolidInT:<:6  /
parSolidInBic:<:6  /
parSolidInCar:<:6  /
parSolidInH2O:<:6  /

*Outlet of the adsorber
parGasOutF:<:6    /
parGasOutP:<:6    /
parGasOutT:<:6    /
parGasOutCO2:<:6  /
parGasOutH2O:<:6  /
parSolidOutFm:<:6 /
parSolidOutT:<:6  /
parSolidOutBic:<:6  /
parSolidOutCar:<:6  /
parSolidOutH2O:<:6  /
parDt:<:6/
unitDx.l('a1'):<:6/
unitLb.l('a1'):<:6/
unitDLX.l('a1'):<:6/
UtilInT:<:6/
ColdOutT.l('a1'):<:6
;
putclose GamsOutput;


** This section evaluates the optimization problem using the new cost calculations (only for the adsorbers and regenerators)
**Next release should replace the previous cost with this calculations.


variable cost_height(s), cost_Pdesing(s), Cost_tp_calc(s), Cost_Pdesign(s), Cost_tw(s),
         Cost_ts_est(s), Cost_Shell_w(s),  Cost_Shell(s),
         Cost_pl(s), Cost_plates(s), Cost_hx_area(s), Cost_unit(s), Cost_hx(s),
Total_cost;
Parameter
Cost_MaxStress        /948/
Cost_weldEff          /0.85/
Cost_tc               /0.003175/
Cost_ts               /0.01905/
Cost_Steel_dens       /7850/
Cost_shell_mFact      /1/
Cost_hx_fl               /1.12/
Cost_hx_oversize      /1.15/
Cost_tp_adj           /0.0127/ ;
Cost_shell_w.lo(s)=10;

Equations Cost1, Cost2, Cost3, Cost4, Cost5, Cost6, Cost7, Cost8, Cost9, Cost10, Cost11, Cost12, eqCapEX_new ,Tot_cost ;
        //Shell height, add 2 m to bed depth for internal stuff
Cost1(s)..      Cost_height(s) =e= unitLb(s) + 2;

        //Design pressure (SSL 2009 eq 22.61), for low pressures use minimum 10 psig
cost2(s)..      Cost_Pdesign(s) =e= 1.69;

        //Thickness based on design pressure (SSL 2009 eq 22.60) (Probably much lower than min thickness, calcuated as a check)
Cost3(s)..      Cost_tp_calc(s)*(2*Cost_MaxStress*Cost_WeldEff - 1.2*Cost_Pdesign(s)) =e= Cost_Pdesign(s)*unitD(s);

        //Additional thickness needed for wind load (SSL 2009 eq 22.62) Diameter should be outside, but not much different than inside diameter
Cost4(s)..      Cost_tw(s)*39.3701*Cost_MaxStress*14.504*(unitD(s)*39.3701)**2 * 1e-8 =e= 0.22*(unitD(s)*39.3701 + 18)*(Cost_height(s)*39.3701)**2*1e-8;

        //Thickness is windload + min thickness + corrosion allowance + rounding adjustment
Cost5(s)..      Cost_ts_est(s) =e= Cost_tw(s) + Cost_tp_adj + Cost_tc;

        //Weight of shell and ends
Cost6(s)..      Cost_shell_w(s) =e= 3.14159*(unitD(s) + Cost_ts)*(Cost_height(s) + 0.8*unitD(s))*Cost_ts*Cost_steel_dens;

        // Cost of the shell for tower (SSL 2009 eq 22.57)
Cost7(s)..      Cost_shell(s) =e= Cost_shell_mFact*exp(7.0132 + 0.18255*log(Cost_shell_w(s)*2.204623)
                                                                         + 0.02297*log10(Cost_shell_w(s)*2.204623)**2);

        // Cost of the platforms and ladders
Cost8(s)..      Cost_pl(s) =e= 361.8*(unitD(s)*3.28084)**0.73960 *(Cost_height(s)*3.28084)**0.70684;

        // Cost of plates (Ulrich 2004 Fig 5.48), diameter in m, (SSL 2009) 1.87 factor to adjust to bubble cap, carbon steel
Cost9(s)..      Cost_plates(s) =e= 1.87*587.97*(unitD(s)**2.0049)*500/400;

        // Cost of internal heat exchanger, no pressure correction needed, low pressure, carbon steel
Cost10(s)..      Cost_hx_area(s) =e= 3.14159*unitLb(s)*unitNx(s)*unitdx(s);
Cost11(s)..      Cost_hx(s) =e= Cost_hx_fl*exp(11.9052 - 0.8709*log(Cost_hx_area(s)*Cost_hx_oversize*10.7639)
                                                 + 0.09005*log10(Cost_hx_area(s)*Cost_hx_oversize*10.7639)**2);

        // Total cost
Cost12(s)..      Cost_unit(s) =e= Cost_shell(s) + Cost_plates(s) + Cost_hx(s) + Cost_pl(s);

Tot_cost.. Total_cost =e= sum(s,Cost_unit(s));
eqCapEX_new $econ('4') ..
   (capEX/(OCF))/1e5 =G= SUM(a,Nu*y(a)*Cost_unit(a)/1e5) + SUM(d,Nu*y(d)*Cost_unit(d)/1e5)
  + 1.10 *Nu* (EXP(11.9052-0.8709*LOG(richHXArea)+0.09005*LOG(richHXArea)**2-LOG(1e5)))
  + 1.10 *Nu* (EXP(11.9052-0.8709*LOG(leanHXArea)+0.09005*LOG(leanHXArea)**2-LOG(1e5)))
  //+ 1.03*1.12*4.25*EXP(11.9052-0.8709*LOG(flueHXArea)+0.09005*LOG(flueHXArea)*LOG(flueHXArea));
  + 1.12*0.98*Nu*((7*EXP(11.9052-0.8709*LOG(12000)+0.09005*LOG(12000)*LOG(12000)-LOG(1e5))+EXP(11.9052-0.8709*LOG(flueHXArea-7*12000)+0.09005*LOG(flueHXArea-7*12000)*LOG(flueHXArea-7*12000)-LOG(1e5))));

* Declare new model considering the new cost calculation

Cost_hx_area.l(s)= 3.14159*unitLb.l(s)*unitNx.l(s)*unitdx.l(s);
model prosyn3 /

** New cost*************
  eqCapEX_new
  Cost1, Cost2, Cost3, Cost4, Cost5, Cost6, Cost7, Cost8, Cost9, Cost10, Cost11, Cost12, tot_cost         /
************************
$ontext

SumOfInfeasibilities, BFBads1,  BFBads2,  BFBads3,  BFBads4,  BFBads6,  BFBads7
  BFBads8,  BFBads9,  BFBads10,  BFBads11,  BFBads12,
  BFBDes1,  BFBDes2,  BFBDes3,  BFBDes4,  BFBDes6,  BFBDes7,  BFBDes8,  BFBDes9
  BFBDes10, BFBDes11,  BFBDes12,  UAx
  objective
  cont1
  cont2
  extra2, extra3,  extra1,  extra4
*  extra5a
*  extra5b

  combo,
  eqDerate, eqSteam, eqNu

  eqUnitCpa,eqUnitCpd, eqUnitW, eqUnitPb, eqUnitPe,
  eqUnitFbAds, eqUnitFbDes, eqUnitFAds, eqUnitFDes


  eqUnitPiAds, eqUnitPoAds, eqUnitPiDes, eqUnitPoDes
  eqUnitHXArea, eqFlueHXArea
  eqEnrgAds0g
  eqStageAds, eqStageDes
  eqstage
  eqNumAds, eqNumDes

  eqMassAds0g, eqEnrgAds0u, eqEnrgAds2g, eqEnrgAds3s
  eqEnrgDes2g, eqEnrgDes3s, eqEnrgDes4g

  eqCaptureTarget
  eqRegenerationbalance
  eqAdsorbebalance

  gasOutConvads
  gasOutConvdes
  gasInadConv
  gasIndeConv
  gasInCS
  gasOutSum
  gasInSum
  gasOutSumd
  gasInSumd
  flueOutSum
/
$offtext
;
Nu.fx= Nu.l; unitLb.fx(s)=unitLb.l(s);   unitD.fx(s)=unitD.l(s);  unitNx.fx(s)=unitNx.l(s); unitdx.fx(s)= unitDx.l(s);
richHXArea.fx= richHXArea.l; leanHXArea.fx=leanHXArea.l;
* The new cost calculations can be used to optimize the problem by including all the other equations in the prosyn3 model
*SOLVE ProSyn3 USING MINLP MINIMIZING Total_cost;
*display COE.l, soi.l,y.l,Nu.l,x.l, CaptureTarget.l, derate.l, steamFlow.l,     capEX.l,      unitCpa.l,  unitCpd.l ;

$ontext
******** loop of initial points evaluations
* Dicopt (cplex and conopt) solver depends on the initial solution and the local optima obtained could be very near from this point
* So, an option to explore more solution could be to fix decisions in diferent inital points and evaluate the solutions.
* The optimization code can be used then to explore more than one solution (since it runs very fast).
* We decided to create a loop to evaluate more initial solutions and obtain the COE for different Nu points

*** Loop to explore more solutions Nu = 12, 14, 16
set lo /1,2,3/;
* these parameters are used to store the solutions from each iteration.
parameters
         parCOE(lo), parsoi(lo), pary(s,lo), parCaptureTarget(lo), parderate(lo),
         parsteamFlow(lo), parCapEx(lo), parunitCpa(a,lo), parunitCpd(d,lo),
         parX(s,t,lo)            Solution if installed or not each stage s using technology t in each iteration lo
         parNu(lo)               Nu fixed values to evaluate in the loop /1 12, 2 14, 3 16/
         NewCostUnits(lo,s)      new costs evaluations
         cputime(lo)             cpu time in seconds
         solution(*,lo)          mod 1=optimal 2=local optimal 4=infeasible stat eq 1=normal termination 2=out iterations 10 to 13 failure
;

loop(lo$(ord(lo) = 10),
Nu.fx=parNu(lo);
SOLVE ProSyn2 USING RMINLP MINIMIZING fplus or soi
******SOLVE ProSyn2 USING MINLP MINIMIZING fplus;
display COE.l, soi.l,y.l,Nu.l,x.l, CaptureTarget.l, derate.l, steamFlow.l,     capEX.l,      unitCpa.l,  unitCpd.l ;

parCOE(lo)$(ord(lo))             = COE.l;
parsoi(lo)$(ord(lo))             = soi.l;
pary(s,lo)$(ord(lo))             = y.l(s);
parCaptureTarget(lo)$(ord(lo))   = CaptureTarget.l;
parderate(lo)$(ord(lo))          = derate.l;
parsteamFlow(lo)$(ord(lo))       = steamFlow.l;
parCapEx(lo)$(ord(lo))         = CapEx.l;
parunitCpa(a,lo)$(ord(lo))       = unitCpa.l(a);
parunitCpd(d,lo)$(ord(lo))       = unitCpd.l(d);
parX(s,t,lo)$(ord(lo))           = X.l(s,t);
cputime(lo)$(ord(lo))                = ProSyn2.resusd;
solution('mod',lo)$(ord(lo))           =ProSyn2.modelstat;
solution('stat',lo)$(ord(lo))           =ProSyn2.solvestat;

Nu.fx= Nu.l; unitLb.fx(s)=unitLb.l(s);   unitD.fx(s)=unitD.l(s);  unitNx.fx(s)=unitNx.l(s); unitdx.fx(s)= unitDx.l(s);
richHXArea.fx= richHXArea.l; leanHXArea.fx=leanHXArea.l;

SOLVE ProSyn3 USING MINLP MINIMIZING Total_cost;
NewCostunits(lo,s)$(ord(lo))     = Cost_unit.l(s);

*);
*close loop
display    parCOE, parsoi, pary, parCaptureTarget, parderate, parsteamFlow
                                                  parCapEx, parunitCpa, parunitCpd, parX, parNu,   NewCostUnits, cputime,solution;


** Print  values
Execute_unload "Superstructure_results.gdx" parCOE, parsoi, pary, parCaptureTarget, parderate, parsteamFlow
                                                  parCapEx, parunitCpa, parunitCpd, parX, parNu,   NewCostUnits, cputime, solution;

*rng = R2 creates a new sheet called "R2" in excel (but, you provide the cells and results can be displayed in a specific location)
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=parCOE.l rng=parCOE!'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=parsoi.l rng=parsoi!'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=pary.l rng=pary!'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=parCaptureTarget.l rng=parCaptureTarget!'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=parderate.l rng=parderate!'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=parsteamFlow.l rng=parsteamFlow!'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=parCapEx.l rng=parCapEx!'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=parunitCpa.l rng=parunitCpa!'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=parunitCpd.l rng=parunitCpd!'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=parX.l rng=parX!'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=parNu.l rng=parNu!'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=NewCostUnits.l rng=NewCostUnits!'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=cputime.l rng=cputime!'
Execute 'gdxxrw.exe Superstructure_results.gdx Squeeze=N par=solution.l rng=solution!'


$offtext
