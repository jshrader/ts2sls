clear all 
sysuse auto.dta
expand 2, gen(group)
replace group = 2 if group == 0
replace mpg = . if group == 1
// Make sure you are using the newest version of ts2sls
capture program drop ts2sls
copy ~/Dropbox/bin/stata/ts2sls/ts2sls.ado ~/Dropbox/bin/stata/ado/t/ts2sls.ado, replace
// Test baseline: one instrument, one exogenous regressor
ts2sls price (mpg = headroom) weight, group(group)
ivregress 2sls price (mpg = headroom) weight if group == 2
// Test multiple instruments
ts2sls price (mpg = c.weight##c.weight) headroom, group(group)
ivregress 2sls price (mpg = c.weight##c.weight) headroom if group == 2
// Test no constant
ts2sls price (mpg = weight) headroom, group(group) noc
ivregress 2sls price (mpg = weight) headroom if group == 2, noc
// Test subsetting (if)
ts2sls price (mpg = weight) headroom if foreign == 0, group(group)
ivregress 2sls price (mpg = weight) headroom if group == 2 & foreign == 0
// Test large matrices
clear
set obs 200000
set seed 1924
gen rand = runiform()
sort rand
gen group = (_n<100001)
replace group = group + 1
gen eps = 50*rnormal()
gen z = rnormal()
gen x1 = z + eps
gen x2 = rnormal()
gen y = 1.5*x1 + .6*x2 + eps
ts2sls y (x1 = z) x2, group(group)
ivregress 2sls y (x1 = z) x2 if group == 2
