clear all

timer clear
timer on 1

set more off
set graphics off

cd "//Users//fho.egb//Documents//GitHub//SimData_County_opioids//"
matrix results_step = (1, 1, 1, 1)

	
set seed 12345
set matsize 11000

local i = 1501

while `i' < 3001 {
	clear
	import delimited "Data_county_`i'.csv"
	sort countycode year
	by countycode: egen group = min(year) if post_period == 1
	capture  did_stepwise y countycode year group, agg cluster(statecode) minn(0)
	matrix results_step = ( results_step \ `i' , -100,  e(b)[1,"tau_average"] , e(V)["tau_average", "tau_average"])
	display `i'
	local i = `i'+1
}

matsave results_step, dropall replace
cd "//Users//fho.egb//Documents//GitHub//StataSim_county//"
export delimited using "sim_stepwise_county_stata2.csv", replace 
timer off 1
timer list 1
clear

