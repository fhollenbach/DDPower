clear all

timer clear
timer on 1

set more off
set graphics off

cd "//Users//fho.egb//Documents//GitHub//SimData_firms//"
matrix results_imp = (1, 1, 1, 1, 1)

	
set seed 12345
set matsize 11000

local i = 4001

while `i' < 8001 {
	clear
	import delimited "Data_firms_`i'.csv"
	sort gvkey time
	by gvkey: egen group = min(time) if post_period == 1
	capture  did_imputation y gvkey time group
	matrix results_imp = ( results_imp \ `i' , 0, -100,  e(b) , e(V) )
	capture  did_imputation y gvkey time group,  leaveout
	matrix results_imp = ( results_imp \ `i', 1, -100,  e(b) , e(V) )
	display `i'
	local i = `i'+1
}

matsave results_imp, dropall replace
cd "//Users//fho.egb//Documents//GitHub//StataSim_firms//"
export delimited using "sim_imputation_firm_stata2.csv", replace 
timer off 1
timer list 1
clear
