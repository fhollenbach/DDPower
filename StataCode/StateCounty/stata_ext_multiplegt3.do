clear all

timer clear
timer on 1

set more off
set graphics off

cd "//Users//fho.egb//Documents//GitHub//SimData_State_opioids//"
matrix results_CdH = (1, 1, 1, 1)

	
set seed 12345
set matsize 11000

local i = 3001

while `i' < 4501 {
	clear
	import delimited "Data_state_`i'.csv"
	capture did_multiplegt_dyn y statecode year post_period, cluster(statecode) effects(15)
	matrix results_CdH = ( results_CdH \ `i', -100,  e(Av_tot_effect) , e(se_avg_total_effect) \  `i', 0, e(Effect_1) , e(se_effect_1) \  `i', 1, e(Effect_2) , e(se_effect_2)  \  `i', 2, e(Effect_3) , e(se_effect_3)  \  `i', 3, e(Effect_4) , e(se_effect_4)  \  `i', 4, e(Effect_5) , e(se_effect_5) )
	display `i'
	local i = `i'+1
}

matsave results_CdH, dropall replace
cd "//Users//fho.egb//Documents//GitHub//StataSim_state//"
export delimited using "sim_CdH_stata3.csv", replace 
timer off 1
timer list 1
clear
