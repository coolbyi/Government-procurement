clear all
use tdata.dta, clear
//// Table 1 
rdrobust cash dis, c(0) p(1) kernel(triangular) bwselect(mserd)  all
estimates store m1
rdrobust cash dis, c(0) p(1) kernel(uniform) bwselect(mserd)  all
estimates store m2
rdrobust cash dis, c(0) p(1) kernel(epanechnikov) bwselect(mserd)  all
estimates store m3
rdrobust pm dis, c(0) p(1) kernel(triangular) bwselect(mserd)  all
estimates store m4
rdrobust pm dis, c(0) p(1) kernel(uniform) bwselect(mserd)  all
estimates store m5
rdrobust pm dis, c(0) p(1) kernel(epanechnikov) bwselect(mserd)  all
estimates store m6
esttab m1 m2 m3 m4 m5 m6 using Table.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f)

//// Table C.2
clear all
use tdata.dta, clear
eststo m1: reg cash same##c.x , vce(robust)
eststo m2: reg cash same##c.(x x2) , vce(robust)
eststo m3: reg cash same##c.(x x2 x3) , vce(robust)
eststo m4: reg cash same##c.(x x2 x3 x4) , vce(robust)
esttab m1 m2 m3 m4  using Table1_rdrobust.rtf, star(* 0.10 ** 0.05 *** 0.01 ) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))

use tdata.dta, clear
eststo m1: reg pm same##c.x , vce(robust)
eststo m2: reg pm same##c.(x x2) , vce(robust)
eststo m3: reg pm same##c.(x x2 x3) , vce(robust)
eststo m4: reg pm same##c.(x x2 x3 x4) , vce(robust)
esttab m1 m2 m3 m4  using Table.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))



//// Table 2 
clear all
use tdata.dta, clear
global controls "pop gdp gov light  scale"
rdrobust cash dis, c(0) p(1) kernel(triangular) bwselect(mserd)  covs($controls)
estimates store m1
rdrobust pm dis, c(0) p(1) kernel(triangular) bwselect(mserd)  covs($controls) 
estimates store m2
eststo m3: reg cash same##c.(x x2), vce(robust)
eststo m4: reg pm same##c.(x x2 x3 x4), vce(robust)
eststo m5: reg cash same##c.(x x2) $controls, vce(robust)
eststo m6: reg pm same##c.(x x2 x3 x4) $controls, vce(robust)
esttab m1 m2 m3 m4 m5 m6  using Table.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))

//// Figure 1 
rdplot cash dis, c(0) p(2) graph_options(legend(off) title(""))
graph save a,replace
rdplot pm dis, c(0) p(4) graph_options(legend(off) title(""))
graph save b,replace
graph combine a.gph b.gph  

//// Figure B.1
clear all
use tdata.dta, clear
quietly sum dis
local N  = r(N)
local lo = r(min)
local hi = r(max)
local J = 100
local w = (`hi' - `lo') / `J'
preserve
gen double mid = `lo' + (`w'/2) + `w'*floor((dis - `lo')/`w')
collapse (count) n = dis, by(mid)
gen double lndens = ln((n/`N')/`w')
local bw = 0.15*(`hi' - `lo')
gen double kw = max(0, 1 - abs(mid)/`bw')
reg lndens mid [aw=kw] if mid < 0
local bl = _b[_cons]
local sl = _se[_cons]
reg lndens mid [aw=kw] if mid > 0 & mid < .
local br = _b[_cons]
local sr = _se[_cons]
local diff = `br' - `bl'
local sed  = sqrt(`sl'^2 + `sr'^2)
di "McCrary test: log-density diff at cutoff = " %6.4f `diff' "  se = " %6.4f `sed'
di "z = " %6.3f `diff'/`sed' "   p-value = " %6.4f 2*(1-normal(abs(`diff'/`sed')))
cap rdplot lndens mid, c(0) p(0) graph_options(title("McCrary-type density plot"))
if _rc rdplot lndens mid, c(0) p(1) graph_options(title("McCrary-type density plot"))
restore


//// Figure B.2
clear all
use tdata.dta, clear
rdplot pop dis, c(0) p(4) graph_options(xtitle("(a)distance") ytitle("pop") legend(off))
graph save a.gph, replace
rdplot gdp dis, c(0) p(4) graph_options(xtitle("(b)distance") ytitle("gdp") legend(off))
graph save b.gph, replace
rdplot gov dis, c(0) p(4)  graph_options(xtitle("(c)distance") ytitle("gov") legend(off))
graph save c.gph, replace
rdplot light dis, c(0) p(4)  graph_options(xtitle("(d)distance") ytitle("light") legend(off))
graph save d.gph, replace
rdplot scale dis, c(0) p(4)  graph_options(xtitle("(e)distance") ytitle("scale") legend(off))
graph save e.gph, replace
graph combine a.gph b.gph c.gph d.gph e.gph, cols(3) graphregion(color(white))


clear all
use tdata.dta, clear
gen absdis_km = abs(dis)
eststo m1: reg cash same##c.(x x2) absdis_km , vce(robust)
eststo m2: reg pm same##c.(x x2 x3 x4) absdis_km , vce(robust)
esttab m1 m2  using Table.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))

//// Table C.1

clear all
use tdata.dta, clear
quietly rdrobust cash dis, kernel(triangular) bwselect(mserd) 
local h_opt_cash = e(h_l)
quietly rdrobust pm dis, kernel(triangular) bwselect(mserd)
local h_opt_pm = e(h_l)

eststo clear
local mults "1.5 2 2.5 3"
local names "15 20 25 30"
local i = 0
foreach m of local mults {
    local i = `i' + 1
    local nm : word `i' of `names'
    local hc = `h_opt_cash' * `m'
    local hp = `h_opt_pm' * `m'
    quietly rdrobust cash dis, kernel(triangular) h(`hc')  
    estimates store bw_cash_`nm'
    quietly rdrobust pm dis, kernel(triangular) h(`hp') 
    estimates store bw_pm_`nm'
}
esttab bw_cash_15 bw_cash_20 bw_cash_25 bw_cash_30 bw_pm_15 bw_pm_20 bw_pm_25 bw_pm_30   using "Table_C1_Bandwidth.rtf", replace star(* 0.10 ** 0.05 *** 0.01) nogap b(%5.3f) se(%5.3f)  stats(N, fmt(0) labels("N")) mtitle("BW = 1.5xMSE" "BW = 2.0xMSE" "BW = 2.5xMSE" "BW = 3.0xMSE") title("Table C.1: Bandwidth robustness")




//// Table C.3
clear all
use tdata.dta, clear
reg cash same##c.(x x2 ) trend, vce(robust)   
estimates store m1
reg pm same##c.(x x2 x3 x4) trend, vce(robust)
estimates store m2
reghdfe cash same##c.(x x2), absorb(year)  vce(robust)
estimates store m3
reghdfe pm same##c.(x x2 x3 x4), absorb(year)  vce(robust)
estimates store m4
esttab m1 m2 m3 m4 using Table.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))


//// Table C.4

clear all
use tdata.dta, clear
gen s  = lon
gen s2 = s^2
gen s3 = s^3
gen s4 = s^4
eststo clear
eststo m1: reg cash same##c.(s s2), vce(robust)
eststo m2: reg pm same##c.(s s2 s3 s4), vce(robust)
gen t  = lat
gen t2 = t^2
gen t3 = t^3
gen t4 = t^4
eststo m3: reg cash same##c.(t t2), vce(robust)
eststo m4: reg pm same##c.(t t2 t3 t4), vce(robust)
esttab m1 m2 m3 m4 using Table.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))


//// Figure C.5

clear all
use tdata.dta, clear
cap confirm variable cash
cap confirm variable dit
if !_rc {
    cap drop dis
    gen dis = dit/1000
}
cap which reghdfe

else {
 
    quietly reg cash same , vce(robust)
    scalar beta_cash = _b[same]
    scalar pval_cash = 2*ttail(e(df_r), abs(_b[same]/_se[same]))

    set seed 1
    permute same beta=_b[same] se=_se[same] df=e(df_r), reps(500) rseed(123) saving("simpm.dta", replace): reghdfe cash same , vce(r)

    preserve
    use "simpm.dta", clear
    gen t_value = beta/se
    gen p_value = 2*ttail(df, abs(beta/se))
    twoway (scatter p_value beta, mcolor(navy%60) msize(small)) (kdensity beta, yaxis(2) lcolor(cranberry)), xline(`=beta_cash', lc(black*0.5) lp(dash)) xline(0, lc(black*0.5) lp(solid)) yline(`=pval_cash', lc(black*0.5) lp(dash)) xlabel(-0.4(0.1)0.4) xtitle("Estimator", size(*0.8)) xlabel(, format(%4.1f) labsize(small)) ytitle("Density", size(*0.8)) ylabel(, nogrid format(%4.1f) labsize(small)) ytitle("P Value", size(*0.8) axis(2)) ylabel(, nogrid format(%4.1f) labsize(small) axis(2)) legend(r(1) order(1 "P Value" 2 "Estimator")) graphregion(color(white)) name(c5a, replace)
    restore

    quietly reg pm same , vce(robust)
    scalar beta_pm = _b[same]
    scalar pval_pm = 2*ttail(e(df_r), abs(_b[same]/_se[same]))

    set seed 1
    permute same beta=_b[same] se=_se[same] df=e(df_r), reps(500) rseed(123) saving("simpm_pm.dta", replace): reghdfe pm same , vce(r)

    preserve
    use "simpm_pm.dta", clear
    gen t_value = beta/se
    gen p_value = 2*ttail(df, abs(beta/se))
    twoway (scatter p_value beta, mcolor(navy%60) msize(small)) (kdensity beta, yaxis(2) lcolor(cranberry)), xline(`=beta_pm', lc(black*0.5) lp(dash)) xline(0, lc(black*0.5) lp(solid)) yline(`=pval_pm', lc(black*0.5) lp(dash)) xlabel(-0.4(0.1)0.4) xtitle("Estimator", size(*0.8)) xlabel(, format(%4.1f) labsize(small)) ytitle("Density", size(*0.8)) ylabel(, nogrid format(%4.1f) labsize(small)) ytitle("P Value", size(*0.8) axis(2)) ylabel(, nogrid format(%4.1f) labsize(small) axis(2)) legend(r(1) order(1 "P Value" 2 "Estimator")) graphregion(color(white)) name(c5b, replace)
    restore

    graph combine c5a c5b , rows(1) graphregion(color(white))
}



//// Table C.6-1
clear all
use tdata.dta, clear
eststo clear
eststo m1: reg cash same##c.(x x2)        if local == 0, vce(robust)
eststo m2: reg pm  same##c.(x x2 x3 x4)  if local == 0, vce(robust)
esttab m1 m2 using Table.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))



//// Table C.6-2
eststo clear
eststo m1: reg cash same##c.(x x2)        if time > 0, vce(robust)
eststo m2: reg cash same##c.(x x2)        if time > 1, vce(robust)
eststo m3: reg cash same##c.(x x2)        if time > 2, vce(robust)
eststo m4: reg pm same##c.(x x2 x3 x4)  if time > 0, vce(robust)
eststo m5: reg pm same##c.(x x2 x3 x4)  if time > 1, vce(robust)
eststo m6: reg pm same##c.(x x2 x3 x4)  if time > 2, vce(robust)
esttab m1 m2 m3 m4 m5 m6 using Table.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))


//// Table C.6-3
eststo clear
eststo m1: reg cash same##c.(x x2)        if cancel == 0, vce(robust)
eststo m2: reg pm  same##c.(x x2 x3 x4)  if cancel == 0, vce(robust)
esttab m1 m2  using Tablec4.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))


//// Table C.6-4
clear all
use tdata.dta, clear
drop if owned > 0 & owned < 3
eststo m1: reg cash same##c.(x x2) , vce(robust)
eststo m2: reg pm   same##c.(x x2 x3 x4) , vce(robust)
esttab m1 m2  using Table.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))


//// Table C.7
clear all
use tdata.dta, clear
eststo m1: reg pm10 same##c.(x x2 x3 x4), vce(robust)
eststo m2: reg o3 same##c.(x x2 x3 x4), vce(robust)
esttab m1 m2  using Table.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))


//// Table 3 Panel A 
clear all
use tdata.dta, clear
eststo m1: reg trade_entity same##c.(x x2), vce(robust)
eststo m2: reg cash  trade_entity same##c.(x x2 ), vce(robust)
esttab m1 m2  using Table.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))

//// Table 3 Panel B 
clear all
use tdata.dta, clear
eststo clear
eststo m1: reg trade_service 1.same##c.(x x2 x3 x4), vce(robust)
eststo m2: reg pm trade_service 1.same##c.(x x2 x3 x4), vce(robust)
esttab m1 m2  using Table.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))



//// Heterogeneity analysis 
clear all
use tdata.dta, clear
//// Heterogeneity: performance assessment ////
rdrobust cash dis if assess == 1, c(0) p(1) kernel(triangular) bwselect(mserd)
estimates store m1
rdrobust cash dis if assess == 0, c(0) p(1) kernel(triangular) bwselect(mserd) 
estimates store m2
rdrobust pm  dis if assess == 1, c(0) p(1) kernel(triangular) bwselect(mserd)
estimates store m3
rdrobust pm  dis if assess == 0, c(0) p(1) kernel(triangular) bwselect(mserd) 
estimates store m4
esttab m1 m2 m3 m4 using Table.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))


//// Heterogeneity: pilot policy ////
clear all
use tdata.dta, clear
rdrobust cash dis if pilot == 1, c(0) p(1) kernel(triangular) bwselect(mserd)
estimates store m1
rdrobust cash dis if pilot == 0, c(0) p(1) kernel(triangular) bwselect(mserd)
estimates store m2
rdrobust pm  dis if pilot == 1, c(0) p(1) kernel(triangular) bwselect(mserd) 
estimates store m3
rdrobust pm  dis if pilot == 0, c(0) p(1) kernel(triangular) bwselect(mserd) 
estimates store m4
esttab m1 m2 m3 m4 using Table.rtf, star(* 0.10 ** 0.05 *** 0.01) nogap replace  se(%5.3f) stats(N r2_a, fmt(0 3))



