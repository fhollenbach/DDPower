clear all

timer clear
timer on 1

set more off
set graphics off

cd "//Users//fho.egb//Documents//GitHub//SimData_State_opioids//"
matrix results_step = (1, 1, 1, 1)

	
set seed 12345
set matsize 11000

local i = 3001

while `i' < 4501 {
	clear
	import delimited "Data_state_`i'.csv"
	sort statecode year
	by statecode: egen group = min(year) if post_period == 1
	capture  did_stepwise y statecode year group, agg cluster(statecode) minn(0)
	matrix results_step = ( results_step \ `i' , -100,  e(b)[1,"tau_average"] , e(V)["tau_average", "tau_average"] )
	display `i'
	local i = `i'+1
}

matsave results_step, dropall replace
cd "//Users//fho.egb//Documents//GitHub//StataSim_state//"
export delimited using "sim_stepwise_stata3.csv", replace 
timer off 1
timer list 1
clear
