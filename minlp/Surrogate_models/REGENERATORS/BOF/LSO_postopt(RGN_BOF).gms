*least square optimization problem (Post optimization), solved to find the best coefficients
*Jan 11, 2016
*By Miguel Zamarripa
*Project: Superstructure optimization

set i(*);  alias (i,ip); set j /1*12/; alias(Z,j);
parameter
*** GasInput variables
vgValue0, GasInP, GasInT, GasInCO2,
*GasInH2O,
*Solids input variables
SolidInFm, SolidInT, SolidInBic, SolidInCar, SolidInH2O,
** Design input variables
Lb, Dx, Dt, Lhx,
** cooling input variables
HXInT,
***Gas output vars (possible surrogate models)
GasInF, GasOutF, GasOutP, GasOutT, GasOutCO2, GasOutH2O, GasOutN2,
*** Solids output vars (possible surrogate models
SolidInP, SolidOutFm, SolidOutP, SolidOutT, SolidOutBic, SolidOutCar, SolidOutH2O,
*** Design output variables
*Ahx, Dte, not used in Regenerator models
vhx, Nx, HXOutT, HXInP, HXInF
** min, max values for the surrogate models
zminval(j), zmaxval(j)
;
;

$CALL  gdxxrw.exe C:\Users\mzamarripa\Documents\Superstructure_project\01_Project_files\SurrogateModels\Regenerator\BOF_surrogates\RGN_BOF_Variables.xlsx @C:\Users\mzamarripa\Documents\Superstructure_project\01_Project_files\SurrogateModels\Regenerator\BOF_surrogates\surrogate_models\parameters_RGN_BOF.txt
$GDXIN RGN_BOF_Variables.gdx
$LOAD i
$LOAD vgValue0 GasInF GasInP GasInT GasInCO2
$LOAD SolidInFm SolidInT SolidInBic SolidInCar SolidInH2O
$LOAD Lb Dx Dt Lhx
$LOAD HXInT HxInF HxInP

*Outputs
$LOAD GasOutF GasOutP GasOutT GasOutCO2 GasOutN2
$LOAD SolidOutFm SolidOutP SolidOutT SolidOutBic SolidOutCar SolidOutH2O
$LOAD Nx vhx
$LOAD HXOutT zminval zmaxval

*C:\Users\mzamarripa\Documents\Superstructure_project\01_Project_files\SurrogateModels\BOF_surrogates\parameters_ADS_BOF.txt

Display i, GasInF, GasInP, GasInT, GasInCO2,
SolidInFm, SolidInT, SolidInBic, SolidInCar, SolidInH2O,
Lb, Dx, Dt, Lhx,
HXInT,
GasOutF, GasOutP, GasOutT, GasOutCO2,
SolidOutFm, SolidOutP, SolidOutT, SolidOutBic, SolidOutCar, SolidOutH2O,
Nx, vhx, HXOutT
zminval, zmaxval
;
*$ontext
parameter valmed aritmetic mean value from sample points ,dim scale of the lst square equation;
parameter minv, maxv,xf(j),xx,estim(i);
variables fx(i),val(i),error(i),lst2,
a(j), b(j), c(j), d(j), a1(j),a2(j),a3(j),a4(j),slack(i)
lst, err(i), SE(i), SEM(i), R2;

equation obj least square error objective function
         c1  Squared error meassured with the predicted parameters
         c2  Squared error the meassured parameters against the mean of the sample
         c3  R2 calculation
         c4  error                      ;

slack.lo(i)=-0.001;Slack.up(i)=.001;
xx=1; xf(j)=0;
dim=1;

equation zmin, zmax, Value,c5,obj2,func;

c1(i)..  SE(i)   =e= (GasOutP(i) - (1.2*GasInP(i)-0.11*(GasInP(i)**2) -0.0012*Dt(i)*Lb(i)));
c2(i)..  SEM(i)  =e= (GasOutP(i) - (sum(ip,GasOutP(ip))/card(i) )) ;
c3..     R2      =e= 1- ( sum(i,power(SE(i),2))/0.00001+sum(i,power(SEM(i),2)) );
c4(i)..  err(i)  =e= GasOutP(i) - (1.2*GasInP(i)-0.11*power(GasInP(i),2) - 0.0012*Dt(i)*Lb(i))   ;
obj..    lst     =e= sum(i, power(err(i),2));

model Least /c4,obj/;
*a.lo=0.0001; b.lo=0.0001;c.lo=0.0001;d.lo=0.0001;
*d.fx=1;
*a.fx=1;
*b.fx=1;
*c.fx=1;

d.l(j)=1;
a.l(j)=2;
b.l(j)=1;
c.l(j)=1;
GasOutH2O(i)=(1 - (GasOutCO2(i) + GasOutN2(i)));
display GasOutH2O;
;
*Nx
func(i).. fx(i) =e=
*- 0.92E+03 * log(Dt(i)) - 0.49E+04 * log(lhx(i)) - 0.11E+05 * exp(dx(i)) + 0.52E+03 * Dt(i) + 0.20E+05 * lhx(i) + 0.26E+05 * dx(i)**2 - 0.10E+05 * lhx(i)**3 - 0.76E+03 * Dt(i)*dx(i) - 0.73E+03 * Dt(i)*lhx(i) + 0.37E+05 * dx(i)*lhx(i)
*xf('1')* ( 458.090398 * (Dt(i)/3.28084) - 1114.34410773 ) Old model (bad)
*+xf('1')* ( 302.75675 * (Lhx(i)**(-1)))   new model (very bad)
+xf('1')*(  (- a1('1') * log(Dt(i)) - 0.49E+04 * log(lhx(i)) - 0.11E+05 * exp(dx(i)) + 0.52E+03 * (Dt(i)**a('1')) + 0.20E+05 * (lhx(i)**b('1')) + 0.26E+05 * dx(i)**2.0 - 0.10E+05 * lhx(i)**3.0 - 0.76E+03 * Dt(i)*dx(i) - 0.73E+03 * Dt(i)*lhx(i) + a2('1') * dx(i)*lhx(i))   )

*GasOutPressure
*-0.8385165779E-1 * Log( unitL(a)/3.28084 ) + 0.9990489 * gasInV(a,'P') )+ P('BFBads2',a) -  N('BFBads2',a) )* x(a,'BOF')
*+ xf('2') * ( -0.8385165779E-1 * Log( dt(i) ) + 0.9990489 * gasInP(i) )
+ xf('2') * ( a1('2') * Log( dt(i)**a('2') ) + a2('2') * gasInP(i)**b('2') )

*GasOutF
*+ xf('3') * ( 2293.84724167 * (Dt(i)) * vgValue0(i) - 58.3781187 * ( GasInT(i)) * vgValue0(i) )
*+ xf('3') * ( (a1('3')*gasInF(i)**a('3')) - (a2('3')*(gasInCO2(i)**b('3'))*gasInCO2(i)**c('3')) )
+ xf('3') * ( a1('3')* Dt(i)**a('3') * vgValue0(i)**b('3') - a2('3') * ( GasInT(i) * vgValue0(i) )**c('3'))

*GasOutT(i) =
*+xf('4')* ( 28.61461625 * log( GasInT(i)) - 34.0536064 * lhx(i) * SolidInCar(i) )
*+ xf('4') * (a1('4')+(a2('4')*log(GasInT(i))) + a3('4')*(lhx(i)**a('4'))*(GasInCO2(i)**b('4')) + a4('4')*(SolidInCar(i)**c('4'))*(lhx(i)**d('4'))   )
+ xf('4') * (a1('4')+(a2('4')*log(GasInT(i))) - a3('4')*(lhx(i)**a('4'))*(SolidInCar(i)**c('4'))   )

*SolidOutT
*+ (26.53430 * log( solidOutT(d+1)$(ORD(d) ne CARD(d)) + solidRichT$(ORD(d) eq CARD(d)) ) - 30.03665 * unitDLX(d) * (SolidoutC(d+1,'NH2COO')$(ORD(d) ne Card(d))+ SolidoutC('a1','NH2COO')$(ORD(d) eq Card(d))) )*x(d,'BOF')
*+ xf('5')*(26.53430 * log( SolidInT(i) ) - 30.03665 * LhX(i) * SolidInCar(i) )
+xf('5')* ((a1('5')*log(SolidInT(i)**a('5'))) - (a2('5') * (SolidInCar(i)**b('5')) * (Lhx(i)**c('5'))) )

*SolidOutH2O
*+( - 0.39933E-001 * unitD(d)/3.28084 * (solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d)) + solidOutC('a1','H2O')$(ORD(d) eq CARD(d)) ) + 1019.8970 * ((solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d)) + solidOutC('a1','HCO3')$(ORD(d) eq CARD(d)) ) * (solidOutC(d+1,'H2O')$(ORD(d) ne CARD(d)) + solidOutC('a1','H2O')$(ORD(d) eq CARD(d)) )) )*x(d,'BOF')
*+xf('6')*( - 0.39933E-001 * Dt(i) * ( SolidInH2O(i) ) + 1019.8970 * (SolidInBic(i)) * (SolidInH2O(i)) )
+xf('6')* (  - a1('6') * Dt(i)**c('6') * solidInH2O(i)  + a2('6') * (SolidInBic(i)**a('6')) * solidInH2O(i)**b('6')    )

*SolidOutBic
*+xf('7')* ( 0.450207E-003 * log(Dt(i)) )
*+ ( 0.82889168 * unitD(a)/3.28084 * (gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)) + 0.780967811E-4 * GasIn(a) * GasInX(a,'H2O')  )
*+xf('7')* (  +(0.82889168 * Dt(i) * gasOutCO2(i)) + (0.780967811E-4 * GasInF(i) * GasInH2O(i))  )
+xf('7') * ( a1('7') * SolidInBic(i) + a3('7') * log( Dt(i)**a('7')) + a2('7')*log(SolidInFm(i)/1000) + a4('7') * GasInCO2(i)**c('7') + 0.00012*Lb(i)**0.100)

*SolidOutCar
*+(0.6002879 * exp(unitDlx(d)) + 0.177458 * ( solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d)) + solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)) )**2.0 )*x(d,'BOF')
*+xf('8') * (0.6002879 * exp(Lhx(i)) + 0.177458 * ( SolidInCar(i) )**2.0 )
+xf('8')* ( a1('8') * exp(Lhx(i)**a('8')) + a2('8') * Dt(i)**c('8') * SolidInCar(i)**b('8') )

*GasOutCO2
*(0.65720137E-002 * unitD(d)/3.28084* ( solidOutC(d+1,'NH2COO')$(ORD(d) ne CARD(d)) + solidOutC('a1','NH2COO')$(ORD(d) eq CARD(d)) ) + 1003.431984 * ((((feedCO2F*feedCO2C('CO2') + steamF*steamC('CO2'))/(feedCO2F+steamF))$(ORD(d) eq 1)+gasOutX(d-1,'CO2')$(ORD(d) ne 1))) * (solidOutC(d+1,'HCO3')$(ORD(d) ne CARD(d)) + solidOutC('a1','HCO3')$(ORD(d) eq CARD(d))) )*x(d,'BOF')
*+xf('9') * (0.65720137E-002 * Dt(i)* ( SolidInCar(i) ) + 1003.431984 * (GasInCO2(i)) * (SolidInBic(i)))
+xf('9')*( (a1('9') * (Dt(i)**a('9'))* SolidInCar(i)**b('9') ) + a2('9')*(GasInCO2(i)**c('9'))*SolidInBic(i)**d('9') )

*GasOutH2O
*( 0.74292231 * (gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))  + P('BFBads12',a) - N('BFBads12',a) )*x(a,'BOF')
*+ xf('10')*( (0.74292231 * GaSInH2O(i)) );
+xf('10')*( 1 - (GasOutCO2(i) +GasOutN2(i)) )

*ColdOut (HXOutT)
+xf('11')*(    HXInT(i) )

;

*Initialization of coefficients
a1.l('9')=0.6572E-002;  a2.l('9')=1003.43; a.l('9')=1; b.l('9')=1; c.l('9')=1; d.l('9')=1;
a1.l('1')=1241.6;  a2.l('1')=35637; a.l('1')=1.014; b.l('1')=0.988; c.l('1')=3;
a1.l('8')=0.6002879 ;  a2.l('8')=0.177458; a.l('8')=1; b.l('8')=2; c.l('8')=1;
a1.l('7')=0.450207E-003;  a2.l('7')=0.780967811E-4; a.l('7')=1; b.l('7')=1;
a.l('5')=1;b.l('5')=1; a2.l('5')=369.25; a1.l('5')=52.9266;
a.l('6')=0.235; b.l('6')=1.002; c.l('6')= -0.095; a1.l('6')=-4.471; a2.l('6')=-17.651;
a1.l('4')=-883.918;  a2.l('4')=14.713; a3.l('4')=22.899; a4.l('4')=902.821; a.l('4')=0.372; b.l('4')=0.548; c.l('4')=0.005; d.l('4')=0.01;


Value(i)..                       val(i) =e= fx(i);
*+ slack(i);
zmin(i,j)$(ord(j) eq xx)..       val(i) =g= zminval(j);
zmax(i,j)$(ord(j) eq xx)..       val(i) =l= zmaxval(j);
**Change GASOUTT or GASOUTP

c5(i)..                          error(i)=e= power((estim(i)- val(i)),2);
obj2..                           lst2=e= sum(i,error(i))/(1+dim);
*a.fx=1.014;b.fx=0.988;a1.fx=1241.6;a2.fx=35637;
option NLP=CONOPT;
*model leastconst /func,Value,zmin,zmax,c5,obj2/
Parameters ParR2(j), ParMinv(j), ParMaxv(j), ParA(j), ParB(j), ParC(j), ParD(j), ParA1(j), ParA2(j), ParA3(j), ParA4(j), Time(j);
model leastconst /func,Value,c5,obj2,zmin,zmax/
*Solve leastconst using nlp minimizing lst2;
*Display val.l, lst2.l, a.l,b.l,c.l,d.l;
*xf selects the j surrogate model to be solved
*xx selects the surrogate model to be used (1 to 10)
loop(Z$(ord(Z)
eq 3),
*gt 0 and ord(z) lt 12),
*initialize parameters to select the variable to be optimized (1 Nx, 2=GasOutP, 3=GasOutF, 4=GasOutT, 5=SolidOutT, 6=SolidOutH2O,
*7=SolidOutBic, 8=SolidOutCar, 9=GasOutCO2, 10= GasOutH2O, 11=HXOutT)
xf(j)=0;xx=0;
xx=(ord(Z));
xf(j)$(ord(j) eq ord(Z))=1;
dim=1; dim=1000$(ord(Z) lt 5);
         loop(i,
         estim(i)=
         +xf('1')*Nx(i) +xf('2')*GasOutP(i) +xf('3')*GasOutF(i) +xf('4')*GasOutT(i) +xf('5')*SolidOutT(i) +xf('6')*SolidOutH2O(i)
         +xf('7')*SolidOutBic(i) +xf('8')*SolidOutCar(i) +xf('9')*GasOutCO2(i) +xf('10')*GasOutH2O(i) +xf('11')*HxOutT(i) );
         display xf, estim,xx;
*to calculate SEM squared error mean value
valmed=sum(ip,estim(ip))/card(i);
*Initialization of the variables
*fx.l(i)=

Solve leastconst using nlp minimizing lst2;
Display val.l, lst2.l, a.l,b.l,c.l,d.l,zminval,zmaxval,valmed;

**R2 calculations
SE.l(i) = power(abs( estim(i) - fx.l(i) ),2);
SEM.l(i)= power((estim(i) - valmed),2) ;
R2.l=  1 - ( sum(i,error.l(i))/sum(i,SEM.l(i)) );
minv=smin(i,val.l(i));
maxv=smax(i,val.l(i));
display SE.l,SEM.l,R2.l,minv,maxv,xf,estim;

**
ParR2(j)$(ord(j) eq ord(Z)) =1 - ( sum(i,error.l(i))/sum(i,SEM.l(i)) );
Parminv(j)$(ord(j) eq ord(Z)) = smin(i,val.l(i));
Parmaxv(j)$(ord(j) eq ord(Z)) = smax(i,val.l(i));


ParR2(j)$(ord(j) eq ord(Z)) =1 - ( sum(i,error.l(i))/sum(i,SEM.l(i)) );
Parminv(j)$(ord(j) eq ord(Z)) = smin(i,val.l(i));
Parmaxv(j)$(ord(j) eq ord(Z)) = smax(i,val.l(i));
Para(j)$(ord(j) eq ord(Z)) = a.l(j);
ParB(j)$(ord(j) eq ord(Z)) = b.l(j);
ParC(j)$(ord(j) eq ord(Z)) = c.l(j);
ParD(j)$(ord(j) eq ord(Z)) = d.l(j);
ParA1(j)$(ord(j) eq ord(Z)) = a1.l(j);
ParA2(j)$(ord(j) eq ord(Z)) = a2.l(j);
ParA3(j)$(ord(j) eq ord(Z)) = a3.l(j);
ParA4(j)$(ord(j) eq ord(Z)) = a4.l(j);
time(j)$(ord(j) eq ord(Z)) = leastconst.resusd;

);
*close loop
display ParR2, Parminv, Parmaxv, ParA, ParB, ParC, ParD, ParA1, ParA2, ParA3, ParA4, time;


parameter sempar;
sempar=sum(i,SEM.l(i));
display sempar;
$ontext
correcto
SE.l(i) = abs( estim(i)*100 - fx.l(i)*100 );
SEM.l(i)= abs( estim(i)*100 - (sum(ip,estim(ip)*100)/card(i)) ) ;
R2.l=  1- ( sum(i,power(SE.l(i),2))/sum(i,power(SEM.l(i),2)) );
minv=smin(i,val.l(i));
maxv=smax(i,val.l(i));

*A simplier option is
Execute_unload "RGN_BOF_Results.gdx" ParR2, Parminv, Parmaxv, ParA, ParB, ParC, ParD, ParA1, ParA2, ParA3, ParA4, time;
*rng = R2 creates a new sheet called "R2" in excel (but, you provide the cells and results can be displayed in a specific location)
Execute 'gdxxrw.exe RGN_BOF_Results.gdx Squeeze=N par=ParR2.l rng=R2!'
Execute 'gdxxrw.exe RGN_BOF_Results.gdx Squeeze=N par=Parminv.l rng=Parminv!'
Execute 'gdxxrw.exe RGN_BOF_Results.gdx Squeeze=N par=Parmaxv.l rng=Parmaxv!'
Execute 'gdxxrw.exe RGN_BOF_Results.gdx Squeeze=N par=ParA.l rng=ParA!'
Execute 'gdxxrw.exe RGN_BOF_Results.gdx Squeeze=N par=ParB.l rng=ParB!'
Execute 'gdxxrw.exe RGN_BOF_Results.gdx Squeeze=N par=ParC.l rng=ParC!'
Execute 'gdxxrw.exe RGN_BOF_Results.gdx Squeeze=N par=ParD.l rng=ParD!'
Execute 'gdxxrw.exe RGN_BOF_Results.gdx Squeeze=N par=ParA1.l rng=ParA1!'
Execute 'gdxxrw.exe RGN_BOF_Results.gdx Squeeze=N par=ParA2.l rng=ParA2!'
Execute 'gdxxrw.exe RGN_BOF_Results.gdx Squeeze=N par=ParA3.l rng=ParA3!'
Execute 'gdxxrw.exe RGN_BOF_Results.gdx Squeeze=N par=ParA4.l rng=ParA4!'

$offtext



