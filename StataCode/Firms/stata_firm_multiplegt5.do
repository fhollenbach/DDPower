clear all

timer clear
timer on 1

set more off
set graphics off

cd "//Users//fho.egb//Documents//GitHub//SimData_firms//"
matrix results_CdH = (1, 1, 1, 1)

	
set seed 12345
set matsize 11000

local i = 32001

while `i' < 36001 {
	clear
	import delimited "Data_firms_`i'.csv"
	capture did_multiplegt_dyn y_log firm_id time post_period, cluster(firm_id) effects(11) graph_off
    matrix results_CdH = ( results_CdH \ `i', -100, e(Av_tot_effect) , e(se_avg_total_effect) )
	display `i'
	local i = `i'+1
}

matsave results_CdH, dropall replace
cd "//Users//fho.egb//Documents//GitHub//StataSim_firms//"
export delimited using "sim_CdH_firm_stata5.csv", replace 
timer off 1
timer list 1
clear
