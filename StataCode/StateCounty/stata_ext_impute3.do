clear all

timer clear
timer on 1

set more off
set graphics off

cd "//Users//fho.egb//Documents//GitHub//SimData_State_opioids//"

matrix results_imp = (1, 1, 1, 1, 1)

	
set seed 12345
set matsize 11000

local i = 3001

while `i' < 4501 {
	clear
	import delimited "Data_state_`i'.csv"
	sort statecode year
	by statecode: egen group = min(year) if post_period == 1
	capture  did_imputation y statecode year group, cluster(statecode)
	matrix results_imp = ( results_imp \ `i' , 0, -100,  e(b) , e(V) )
	capture   did_imputation y statecode year group,  leaveout cluster(statecode)
	matrix results_imp = ( results_imp \ `i', 1, -100,  e(b) , e(V) )
	display `i'
	local i = `i'+1
}

matsave results_imp, dropall replace
cd "//Users//fho.egb//Documents//GitHub//StataSim_state//"
export delimited using "sim_imputation_stata3.csv", replace 
timer off 1
timer list 1
clear
