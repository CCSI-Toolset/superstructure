*least square optimization problem (Post optimization), solved to find the best coefficients
*Jan 11, 2016
*By Miguel Zamarripa
*Project: Superstructure optimization

set i(*);  alias (i,ip); set j /1*11/; alias(Z,j);
parameter
GasInF, GasInP, GasInT, GasInCO2, GasInH2O,
SolidInFm, SolidInT, SolidInBic, SolidInCar, SolidInH2O,
Lb, Dx, Dt, Lhx,
HXInT,
GasOutF, GasOutP, GasOutT, GasOutCO2, GasOutH2O,
SolidOutFm, SolidOutP, SolidOutT, SolidOutBic, SolidOutCar, SolidOutH2O,
Ahx, Dte, Nx, HXOutT
zminval(j), zmaxval(j)
;

*$CALL  gdxxrw.exe C:\Users\mzamarripa\Documents\Superstructure_project\01_Project_files\SurrogateModels\BOF_surrogates\ADS_BOF_Variables.xlsx @C:\Users\mzamarripa\Documents\Superstructure_project\01_Project_files\SurrogateModels\BOF_surrogates\parameters_ADS_BOF.txt
$CALL  gdxxrw.exe ADS_BOF_Variables.xlsx @parameters_ADS_BOF.txt
$GDXIN ADS_BOF_Variables.gdx
$LOAD i
$LOAD GasInF GasInP GasInT GasInCO2 GasInH2O
$LOAD SolidInFm SolidInT SolidInBic SolidInCar SolidInH2O
$LOAD Lb Dx Dt Lhx
$LOAD HXInT

*Outputs
$LOAD GasOutF GasOutP GasOutT GasOutCO2 GasOutH2O
$LOAD SolidOutFm SolidOutP SolidOutT SolidOutBic SolidOutCar SolidOutH2O
$LOAD Ahx Dte Nx
$LOAD HXOutT zminval zmaxval

*C:\Users\mzamarripa\Documents\Superstructure_project\01_Project_files\SurrogateModels\BOF_surrogates\parameters_ADS_BOF.txt

Display i, GasInF, GasInP, GasInT, GasInCO2, GasInH2O,
SolidInFm, SolidInT, SolidInBic, SolidInCar, SolidInH2O,
Lb, Dx, Dt, Lhx,
HXInT,
GasOutF, GasOutP, GasOutT, GasOutCO2, GasOutH2O,
SolidOutFm, SolidOutP, SolidOutT, SolidOutBic, SolidOutCar, SolidOutH2O,
Ahx, Dte, Nx, HXOutT
zminval, zmaxval
;

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


;
*Nx
func(i).. fx(i) =e=
*- 0.92E+03 * log(Dt(i)) - 0.49E+04 * log(lhx(i)) - 0.11E+05 * exp(dx(i)) + 0.52E+03 * Dt(i) + 0.20E+05 * lhx(i) + 0.26E+05 * dx(i)**2 - 0.10E+05 * lhx(i)**3 - 0.76E+03 * Dt(i)*dx(i) - 0.73E+03 * Dt(i)*lhx(i) + 0.37E+05 * dx(i)*lhx(i)
+xf('1')*(  (- a1('1') * log(Dt(i)) - 0.49E+04 * log(lhx(i)) - 0.11E+05 * exp(dx(i)) + 0.52E+03 * (Dt(i)**a('1')) + 0.20E+05 * (lhx(i)**b('1')) + 0.26E+05 * dx(i)**2.0 - 0.10E+05 * lhx(i)**3.0 - 0.76E+03 * Dt(i)*dx(i) - 0.73E+03 * Dt(i)*lhx(i) + a2('1') * dx(i)*lhx(i))   )

*GasOutPressure
*-0.8385165779E-1 * Log( unitL(a)/3.28084 ) + 0.9990489 * gasInV(a,'P') )+ P('BFBads2',a) -  N('BFBads2',a) )* x(a,'BOF')
*+ xf('2') * ( -0.8385165779E-1 * Log( dt(i) ) + 0.9990489 * gasInP(i) )
+ xf('2') * ( a1('2') * Log( dt(i)**a('2') ) + a2('2') * gasInP(i)**b('2') )

*GasOutF
*+ xf('3') * ( (gasInF(i)*0.9469509058) - (23302.38325895*gasInH2O(i)*gasInCO2(i)) )
+ xf('3') * ( (a1('3')*gasInF(i)**a('3')) - (a2('3')*(gasInH2O(i)**b('3'))*gasInCO2(i)**c('3')) )

*GasOutT(i) =
*+((  (15.304891940*Log(gasInV(a,'T')) ) + (262.55799341*gasInX(a,'CO2')*unitDLX(a))  ) +  P('BFBads4',a)  -  N('BFBads4',a) )*x(a,'BOF')
****15*log(GasInT(i)) + 5.5*Dt(i)*GasInCO2(i) + 0.13E+03 * GasInH2O(i)*lhx(i) ;
*+ xf('4') * (15.304891940*(log(GasInT(i))) + (262.55799*GasInCO2(i) *lhx(i))   )
*+ xf('4') * (15.304891940*(log(GasInT(i))) + 5.5*(Dt(i))*(GasInCO2(i)) + 0.13E3*(GasInH2O(i))*(lhx(i))   )
+ xf('4') * ( a1('4')*log(GasInT(i)) + a2('4')*SolidInFm(i)/1E6 + a3('4')*(Dt(i)**a('4'))*(GasInCO2(i)**b('4')) + a4('4')*(GasInH2O(i)**c('4'))*(lhx(i)**d('4')) )

*SolidOutT
*+xf('5')* (52.926674*EXP(gasInH2O(i)) + 369.2543375 * gasInCO2(i) * Lhx(i) )
* 52.926674*EXP(gasInX(a,'H2O')) + 369.2543375 * gasInX(a,'CO2') * unitDLX(a)
*+xf('5') *( ( 16.718826*log(SolidInT(i)) + 463.4708590 * gasInCO2(i) * gasInH2O(i)   ) )
+xf('5')* ((a1('5')*EXP(gasInH2O(i)**a('5'))) + (a2('5') * (gasInCO2(i)**b('5')) * (Lhx(i)**c('5'))) )

*SolidOutH2O
*+xf('6')* (  -0.182760 * log(gasInCO2(i)) + (6.6747347 * gasInH2O(i)) + solidInH2O(i) )
+xf('6')* (  -a1('6') * log(gasInCO2(i)) + a2('6') * gasInH2O(i)**a('6') + a3('6') * solidInH2O(i)**b('6')  + a4('6')*SolidInFm(i)**c('6')  )
*+xf('6')* (  -0.1827 * log(gasInCO2(i)) + 6.6747 * gasInH2O(i) + solidInH2O(i)  + a1('6')*SolidInFm(i)**a('6') +a2('6')*Lb(i)*a3('6')*Dx(i)+ a4('6')*Dt(i)  )
*SolidOutBic
*+ ( 0.82889168 * unitD(a)/3.28084 * (gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)) + 0.780967811E-4 * GasIn(a) * GasInX(a,'H2O')  )
*+xf('7')* (  +(0.82889168 * Dt(i) * gasOutCO2(i)) + (0.780967811E-4 * GasInF(i) * GasInH2O(i))  )
+xf('7')* (  +(a1('7') * Dt(i)**a('7') * gasOutCO2(i)) + ( a2('7')* GasInF(i) * GasInH2O(i)**b('7') ) )

*SolidOutCar
*+( 0.389431104 * exp(SolidOutC(a+1,'NH2COO')$(ORD(a) ne CARD(a)) *exp(SolidOutC('d1','NH2COO'))$(ORD(a) eq CARD(a)) + 0.4755337466 * unitD(a)/3.28084 * GasInX(a,'CO2')  )
*+xf('8')* ( 0.389431104 * exp(SolidInCar(i)) + 0.4755337466*Dt(i) * GasInCO2(i) )
+xf('8')* ( a1('8') * exp(SolidInCar(i)**a('8')) + a2('8') * Dt(i)**c('8') * GasInCO2(i)**b('8') )

*GasOutCO2
*+( 0.48976547 * (gasOutX(a-1,'CO2')$(ORD(a) ne 1)+(flueOutC('CO2')/flueOut)$(ORD(a) eq 1)) + P('BFBads11',a)  -  N('BFBads11',a) )*x(a,'BOF')
*+xf('9')*( 0.48976547 * GaSInCO2(i) )
+xf('9')*( (a1('9') * GaSInCO2(i)**a('9')) + a2('9')*Dt(i)**b('9') +a3('9') * Lb(i)**c('9') * SolidInCar(i)+ a4('9') * Dx(i) * LhX(i)**d('9') )

*GasOutH2O
*( 0.74292231 * (gasOutX(a-1,'H2O')$(ORD(a) ne 1)+(flueOutC('H2O')/flueOut)$(ORD(a) eq 1))  + P('BFBads12',a) - N('BFBads12',a) )*x(a,'BOF')
*+ xf('10')*( (0.74292231 * GaSInH2O(i)) );
+xf('10')*( (a1('10') * GaSInH2O(i)**a('10')) + (a2('10')*Dt(i)**b('10'))+ (a3('10')*GasInCO2(i)*c('10')) )

*ColdOut (HXOutT)
*(-2.4 * Log(unitDx(a)) + 84 * ColdInT *SolidOutC(a+1,'HCO3')$(ord(a) ne card(a)) *SolidOutC('d1','HCO3')$(ord(a) eq card(a)) )*x(a,'BOF')
*+xf('11')*(    +(-2.4*Log(Dx(i)))+ ( 84 * HXInT(i) * SolidInBic(i) )  )
*+xf('11')*( +( a1('11') * Log(Dx(i)**a('11')) )+ ( a2('11') * HXInT(i)**b('11') * SolidInBic(i)**c('11') )  )
+xf('11')*(a1('11') + a2('11')*HXInT(i) + a3('11') * Lb(i)**b('11') )
* + a3('11') * Dx(i)**c('11')  )
;

*Initialization of coefficients
a1.l('1')=1241.6;  a2.l('1')=35637; a.l('1')=1.014; b.l('1')=0.988; c.l('1')=3;
a1.l('8')=0.38943;  a2.l('8')=0.47553; a.l('8')=1; b.l('8')=1;
a1.l('7')=0.82889168;  a2.l('7')=0.780967811E-4; a.l('7')=1; b.l('7')=1;
a.l('5')=1;b.l('5')=1; a2.l('5')=369.25; a1.l('5')=52.9266;
a.l('6')=0.196;b.l('6')=-0.137; a1.l('6')=0.919; a2.l('6')=-1.236;
a1.l('4')=15.304891940;  a2.l('4')=1E-5; a3.l('4')=5.5; a4.l('4')=0.13E3; a.l('4')=0.372; b.l('4')=0.548; c.l('4')=0.005; d.l('4')=0.01;
a1.l('11')=1;  a2.l('11')=1; a3.l('11')=1; a.l('11')=1.014; b.l('11')=0.988; c.l('11')=1;

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
PArameters ParR2(j), ParMinv(j), ParMaxv(j), ParA(j), ParB(j), ParC(j), ParD(j), ParA1(j), ParA2(j), ParA3(j), ParA4(j), Time(j);
model leastconst /func,Value,c5,obj2,zmin,zmax/
*Solve leastconst using nlp minimizing lst2;
*Display val.l, lst2.l, a.l,b.l,c.l,d.l;
*xf selects the j surrogate model to be solved
*xx selects the surrogate model to be used (1 to 10)
loop(Z$(ord(Z) eq 9),
*loop(Z$(ord(Z) gt 0),
*initialize parameters to select the variable to be optimized (1 Nx, 2=GasOutP, 3=GasOutF, 4=GasOutT, 5=SolidOutT, 6=SolidOutH2O,
*7=SolidOutBic, 8=SolidOutCar, 9=GasOutCO2, 10= GasOutH2O, 11=HXOutT)
xf(j)=0;xx=0;
xx=(ord(Z));
xf(j)$(ord(j) eq ord(Z))=1;
dim=1; dim=1000$(ord(Z) lt 5 and Ord(Z) eq 11);
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
*Display val.l, lst2.l, a.l,b.l,c.l,d.l,zminval,zmaxval,valmed;



**R2 calculations
SE.l(i) = power(abs( estim(i) - fx.l(i) ),2);
SEM.l(i)= power((estim(i) - valmed),2) ;
R2.l=  1 - ( sum(i,error.l(i))/sum(i,SEM.l(i)) );
minv=smin(i,val.l(i));
maxv=smax(i,val.l(i));

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
display ParR2, Parminv, Parmaxv,ParA,ParB,ParC,ParA1,ParA2,ParA3, zminval,zmaxval;



*display SE.l,SEM.l,R2.l,minv,maxv,xf,estim;
parameter sempar;
sempar=sum(i,SEM.l(i));
*display sempar;
$ontext
correcto
SE.l(i) = abs( estim(i)*100 - fx.l(i)*100 );
SEM.l(i)= abs( estim(i)*100 - (sum(ip,estim(ip)*100)/card(i)) ) ;
R2.l=  1- ( sum(i,power(SE.l(i),2))/sum(i,power(SEM.l(i),2)) );
minv=smin(i,val.l(i));
maxv=smax(i,val.l(i));


** Export results from GAMS to excel (using GDX files)
*this option neds the correct path and parameters file
*execute_unload 'C:\Users\mzamarripa\Documents\Superstructure_project\01_Project_files\SurrogateModels\BOF_surrogates\ADS_BOF_Results.gdx' ParR2, Parminv, Parmaxv, ParA, ParB, ParC, ParA1, ParA2, ParA3, time;
*Execute 'GDXXRW.EXE "C:\Users\mzamarripa\Documents\Superstructure_project\01_Project_files\SurrogateModels\BOF_surrogates\ADS_BOF_Results.gdx" O="PAth\ADS_BOF_Results.xlsx" @"Output_parameters"';


*A simplier option is
Execute_unload "ADS_BOF_Results.gdx" ParR2, Parminv, Parmaxv, ParA, ParB, ParC, ParD, ParA1, ParA2, ParA3, ParA4, time;
*rng = R2 creates a new sheet called "R2" in excel (but, you provide the cells and results can be displayed in a specific location)
Execute 'gdxxrw.exe ADS_BOF_Results.gdx Squeeze=N par=ParR2.l rng=R2!'
Execute 'gdxxrw.exe ADS_BOF_Results.gdx Squeeze=N par=Parminv.l rng=Parminv!'
Execute 'gdxxrw.exe ADS_BOF_Results.gdx Squeeze=N par=Parmaxv.l rng=Parmaxv!'
Execute 'gdxxrw.exe ADS_BOF_Results.gdx Squeeze=N par=ParA.l rng=ParA!'
Execute 'gdxxrw.exe ADS_BOF_Results.gdx Squeeze=N par=ParB.l rng=ParB!'
Execute 'gdxxrw.exe ADS_BOF_Results.gdx Squeeze=N par=ParC.l rng=ParC!'
Execute 'gdxxrw.exe ADS_BOF_Results.gdx Squeeze=N par=ParD.l rng=ParD!'
Execute 'gdxxrw.exe ADS_BOF_Results.gdx Squeeze=N par=ParA1.l rng=ParA1!'
Execute 'gdxxrw.exe ADS_BOF_Results.gdx Squeeze=N par=ParA2.l rng=ParA2!'
Execute 'gdxxrw.exe ADS_BOF_Results.gdx Squeeze=N par=ParA3.l rng=ParA3!'
Execute 'gdxxrw.exe ADS_BOF_Results.gdx Squeeze=N par=ParA4.l rng=ParA4!'

$offtext
