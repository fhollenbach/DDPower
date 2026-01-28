clear all

timer clear
timer on 1

set more off
set graphics off

cd "//Users//fho.egb//Documents//GitHub//SimData_State_modelbased//"
matrix results_CdH = (1, 1, 1)

	
set seed 12345
set matsize 11000

local i = 1

while `i' < 1501 {
	clear
	import delimited "Data_state_`i'.csv"
	capture did_multiplegt_dyn y statecode year post_period, cluster(statecode) effects(15)
	matrix results_CdH = ( results_CdH \ `i',  e(Av_tot_effect) , e(se_avg_total_effect))
	display `i'
	local i = `i'+1
}

matsave results_CdH, dropall replace
cd "//Users//fho.egb//Documents//GitHub//StataSim_state_modelbased//"
export delimited using "sim_CdH_stata1.csv", replace 
timer off 1
timer list 1
clear
