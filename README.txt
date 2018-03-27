    = Superstructure Formulation Project =

The superstructure optimization code is based on GAMS modeling system (.gms), 
and is an update from the previous release version (2015 and 2014). 

The project files are distributed as follows: 
• DOCs folder (located in Superstructure Bundle) 
	o Superstructure Formulation User Manual.pdf

• Minlp folder
	o 1_Supr_final_proj.gpr
	o 2_Super-2016-windows.gms
	o 3_Super-2016-windows.lst
	o 4_ADSORBER.gms
	o 5_REGENERATOR.gms
	o 6_Superstructure_results.xlsx
	o Surrogate Models folder
		 ADSORBERS folder
			•BOF folder	
				o 7_LSO_postopt(ADS_BOF).gms
				o 8_ADS_BOF_Variables.xlsx
				o 9_ADS_BOF_Results.xlsx
			•BUF folder
				o LSO_postopt(ADS_BUF).gms
				o ADS_BOF_Variables.xlsx
				o ADS_BOF_Results.xlsx
		 REGENERATORS folder
			• BOF folder
				o LSO_postopt(RGN_BOF).gms
				o ADS_BOF_Variables.xlsx
				o ADS_BOF_Results.xlsx
			• BUF folder
				o LSO_postopt(RGN_BUF).gms
				o ADS_BOF_Variables.xlsx
				o ADS_BOF_Results.xlsx

First initialize GAMS IDLE and open the project file (file # 1, “Supr_final_proj.gpr”), 
by doing this the user will be running all the files in the same directory (minlp). 

The file “super-2016-windows.gms” is the main optimization code, 
which includes/calls the ADSORBER.gms and REGENERATOR.gms surrogate models using “$include” function of GAMS. 

Files 4 and 5, correspond to the surrogate models created with the data obtained with FOQUS session. 

The folders called “surrogate models” include all the data sets for each technology. 
The ADS_BOF_Variables.xlsx includes: 
	i) 	the input and output variables required for the surrogate models; 
	ii) 	the maximum and minimum values used to simulate and run all the samples (the sampling scheme used was Latin Hypercube); 
	iii) 	the data samplings input and output variables are included in this file. The data sets are then used by ALAMO to obtain the surrogate models. 

The file LSO_postopt(Technology).gms runs a least square optimization problem to improve the fitting of surrogate models obtained with ALAMO and the results are included in the Technology_Results.xlsx file.

