clear all

timer clear
timer on 1

set more off
set graphics off

cd "//Users//fho.egb//Documents//GitHub//SimData_firms//"
matrix results_step = (1, 1, 1, 1)

	
set seed 12345
set matsize 11000

local i = 12001

while `i' < 16001 {
	clear
	import delimited "Data_firms_`i'.csv"
	sort gvkey time
	by gvkey: egen group = min(time) if post_period == 1
	capture did_stepwise y gvkey time group, agg cluster(gvkey) minn(0)
	matrix results_step = ( results_step \ `i' , -100,  e(b)[1,"tau_average"] , e(V)["tau_average", "tau_average"] )
	display `i'
	local i = `i'+1
}

matsave results_step, dropall replace
cd "//Users//fho.egb//Documents//GitHub//StataSim_firms//"
export delimited using "sim_stepwise_firm_stata4.csv", replace 
timer off 1
timer list 1
clear
