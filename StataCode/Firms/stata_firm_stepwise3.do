clear all

timer clear
timer on 1

set more off
set graphics off

cd "//Users//fho.egb//Documents//GitHub//SimData_firms//"
matrix results_step = (1, 1, 1)

	
set seed 12345
set matsize 11000

local i = 10001

while `i' < 15001 {
	clear
	import delimited "Data_firms_`i'.csv"
	sort firm_id time
	by firm_id: egen group = min(time) if post_period == 1
	capture did_stepwise y_log firm_id time group, agg cluster(firm_id) minn(0)
	matrix results_step = ( results_step \ `i' ,  e(b)[1,"tau_average"] , e(V)["tau_average", "tau_average"] )
	display `i'
	local i = `i'+1
}

matsave results_step, dropall replace
cd "//Users//fho.egb//Documents//GitHub//StataSim_firms//"
export delimited using "sim_stepwise_firm_stata3.csv", replace 
timer off 1
timer list 1
clear
