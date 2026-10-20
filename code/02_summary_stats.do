*===============================================================================
* 02_summary_stats.do - summary statistics and sample comparisons of
*   "Natural capital restoration and household wealth in sub-Saharan Africa"
*   (Rizzi, Xu, Heal). Replication package, github.com/claudiorizzi/natural_capital_restoration.
*
* Defines four programs; run them from code/run_demo.do (which sets the path
* globals data_dir, out_tab, out_fig) or interactively after `do` of this file.
*   sum_stats_farms          summary_stats.tex             (Supp. Table 3, household rows)
*   sum_stats_evi            summary_stats_EVI.tex         (Supp. Table 3, EVI rows)
*   wb_comparison            comparison_table.tex          (Supp. Table 11)
*   single_multi_comparison  comparison_single_multi.tex   (Supp. Table 12)
*===============================================================================

*** Preamble ***
clear all
set more off
set matsize 11000
set maxvar 32700

//data and output folders are set by code/run_demo.do (globals data_dir, out_tab, out_fig)
cd "$data_dir"

global controls i.wb_adult_mostedu i.wb_gender wb_w_npeople wb_w_size wb_w_participantage
global gcpc w1ch_gcpc_2017usd

* =========================================
* SUMMARY STATS - household sample
*   -> $out_tab/summary_stats.tex  (Supp. Table 3)
* =========================================
capture program drop sum_stats_farms
program define sum_stats_farms

eststo clear
use "final/TREES_data_pseudo", clear

//farm size and household members of 0 (or below) are non-response coded as
//zero, not landless farms or empty households: the per-ha variables are already
//missing for zero-size farms (division by 0), so treat the size itself as
//missing too. Summary stats for these two variables are
//therefore computed only over observations with values strictly greater than 0.
replace size=. if size<=0
replace npeople=. if npeople<=0
 
label variable havelivestock "1(Own livestock)"
label variable tot_animal "N. livestock"
label variable gender "1(Male)"
label variable npeople "Tot. household members"
label variable females_tot "\% of female members"
label variable hfiascategory "Household ood Insecurity Access Survey (HFIAS)"
label variable worriedfood "1(Worried about food)"
label variable reliantothers "1(Reliant on others)"
label variable gcp_2017usd "Gross Cell Product (USD mil.)"
//label variable pop_cell "Population per Cell"
label variable w_foodcrops "N. food crops"
label variable w_marketableproducts "N. marketable products"
label variable nofoodid "1(No food)"
label variable hfiascatid4 "1(Food insecure)"
label variable hdds "Food Diversity"
label variable hdds_id3 "Low food diversity"
label variable w1ch_gcpc_2017usd "$\Delta$ GCP per capita"
label variable w1_gcpc_2017usd "GCP per capita (USD)"
label var size "Farm size (ha)"
label var w_trees "N. trees"


//simple summary stats
local summary_vars havelivestock tot_animal w_trees w_marketableproducts w_foodcrops nofoodid hfiascatid4 hdds_id3 size npeople w1_gcpc_2017usd 


eststo clear
estpost sum `summary_vars' if gap_ok==1 & gap_ok_fl==1 , detail
est store a
esttab a using "$out_tab/summary_stats.tex", nonumbers booktabs replace compress nogaps b(2) /// ///
collabels({Mean} {Median} {St.Dev.} {Obs.} {Min.} {Max.}) ///
cells("mean(fmt(2)) p50(fmt(2)) sd(fmt(2)) count(fmt(%9.0fc)) min(fmt(2)) max(fmt(2))") label noobs

end

* =========================================
* SUMMARY STATS - EVI
*   -> $out_tab/summary_stats_EVI.tex  (Supp. Table 3)
* =========================================
capture program drop sum_stats_evi
program define sum_stats_evi

//EVI summary stats
eststo clear
use "final/EVI_farms_pseudo", clear
label var w1_evi_mean "EVI mean ($\times$ 100)"
label var w1_evi_median "EVI median ($\times$ 100)"
gen event_window=year-calyear_treat if calyear_treat!=0
keep if (event_window>=-3 & event_window<=4) | calyear_treat==0
keep if !missing(w1_evi_median) & !missing(gcpc_2017usd_w1)

local summary_vars w1_evi_median w1_evi_mean

eststo clear
estpost sum `summary_vars', detail
est store a
esttab a using "$out_tab/summary_stats_EVI.tex", nonumbers booktabs replace compress nogaps b(2) /// ///
collabels({Mean} {Median} {St. Dev.} {Obs.} {Min.} {Max.}) ///
cells("mean(fmt(2)) p50(fmt(2)) sd(fmt(2)) count(fmt(%9.0fc)) min(fmt(2)) max(fmt(2))") label noobs

end

* =========================================
* WORLD BANK COMPARISON - TREES farms vs LSMS-ISA Tanzania households
*   -> $out_tab/comparison_table.tex  (Supp. Table 11)
* =========================================
capture program drop wb_comparison
program define wb_comparison

eststo clear
use "final/TREES_data_pseudo", clear
	//(the package ships one household file; the paper used the earlier
	// household extract farm_income.dta for this table, with year_hhs = year_survey)

//farm size of exactly 0 is non-response coded as zero: treat as missing,
//consistent with the other summary-stats programs
replace size=. if size==0

keep if inlist(year_survey,2019, 2020, 2021, 2022) & iso_trees=="TZA" & treat_year==0

gen farm_main_income=1 if farmingmainsource=="Yes"
replace farm_main_income=0 if farmingmainsource=="No"
gen reduce_n_meals=1 if fewermeals=="Yes"
replace reduce_n_meals=0 if fewermeals=="No"
gen child_in_school=1 if childreninschool=="Yes"
replace child_in_school=0 if childreninschool=="No"
gen children=1 if householdchildren=="Yes"
replace children=0 if householdchildren=="No"

gen hfias_category=0 if strpos(hfiascategory, "Secure")>0
replace hfias_category=1 if strpos(hfiascategory, "Mild")>0
replace hfias_category=2 if strpos(hfiascategory, "Moderate")>0 
replace hfias_category=3 if strpos(hfiascategory, "Severe")>0 
drop hfiascategory

//adjust householdhighestschoollevel to an indicator
gen highest_sch_level=0 if householdhighestschoollevel=="never attended"
replace highest_sch_level=1 if householdhighestschoollevel=="some primary" | householdhighestschoollevel=="primary"
replace highest_sch_level=2 if householdhighestschoollevel=="some secondary/tech" | householdhighestschoollevel=="secondary/tech"
replace highest_sch_level=3 if householdhighestschoollevel=="higher"

//number of crop types
gen n_crops=marketableproducts+foodcrops
rename size farm_area 

//females %
gen p_females=householdfemales/npeople

//children %
gen p_child5yrs=n_childrenunder5/npeople

keep forestgardenid year_survey project proj_area iso_trees npeople householdfemales p_females n_childrennumber n_childrenunder5 p_child5yrs highest_sch_level farm_main_income farm_area trees trees_per_ha havelivestock tot_animal tot_animal_ha n_crops hfias hfias_category hfiascatid4 worriedfood nofoodid  reduce_n_meals nofoodid 
 
gen dataset=1 //1 = TREES
append using "final/tanzania_merged"
replace dataset=0 if missing(dataset) //World Bank Survey
//children number from TREES and world bank might differ due to not specified age limit for children in TREES

*-------
//only same regions as trees. proj_area= OSU is in Marwa village (Manyara region)
keep if inlist(region, "manyara", "singida", "tabora") | missing(region)

label variable npeople              "Household size"
label variable householdfemales     "Number of females in household"
label variable p_females            "Share of females in household"

label variable n_childrennumber     "Number of children"
label variable n_childrenunder5     "Number of children under 5"
label variable p_child5yrs          "Share of children under 5"

label variable highest_sch_level    "Highest education level"

label variable farm_main_income     "Main income from farming"

label variable farm_area            "Farm size (hectares)"

label variable trees                "Number of trees"
label variable trees_per_ha         "Trees per hectare"

label variable tot_animal           "Number of livestock"
label variable tot_animal_ha        "Livestock per hectare"

label variable n_crops              "Number of crops"

label variable hfias                "HFIAS score"
label variable hfias_category       "HFIAS category"
label variable hfiascatid4 "1(Food Insecure)"

label variable worriedfood          "Worried about food"
label variable nofoodid             "1(No Food)"
label variable reduce_n_meals       "Reduced number of meals"

tempfile temp
save `temp', replace

* ------------------------------------------
* Variables to compare
* ------------------------------------------
use `temp', clear


local vars npeople n_childrennumber n_childrenunder5 p_child5yrs highest_sch_level farm_main_income farm_area trees trees_per_ha tot_animal tot_animal_ha n_crops  nofoodid

  
* Create empty results dataset
clear
set obs 0
gen varname = ""
gen mean_trees = .
gen mean_wb = .
gen diff = .
gen N_trees = .
gen N_wb = .

tempfile results
save `results', replace

* ------------------------------------------
* Loop over variables
* ------------------------------------------
foreach v of local vars {

    use `temp', clear
	
	* Grab the variable label BEFORE doing anything else
    local vlabel : variable label `v'
    if "`vlabel'" == "" local vlabel "`v'"   // fallback to name if no label
	
    * Means by group. For livestock calculation, make it conditional on owning livestock (WB survey is conditional on livestock own.)
	if inlist("`v'", "tot_animal" , "tot_animal_ha") {
		quietly summarize `v' if dataset == 1 & havelivestock==1 & `v' != .
		local mean1 = r(mean)
		local N1 = r(N)
	}
	else {
		quietly summarize `v' if dataset == 1 & `v' != .
		local mean1 = r(mean)
		local N1 = r(N)
	}
	
    quietly summarize `v' if dataset == 0 & `v' != .
    local mean0 = r(mean)
	local N0 = r(N)

    * Regression
	if inlist("`v'", "tot_animal" , "tot_animal_ha") {
		quietly reg `v' dataset if `v' != . & havelivestock==1, robust
	} 
	else {
		quietly reg `v' dataset if `v' != ., robust
		}	
	
	
    local b = _b[dataset]
    local t = _b[dataset] / _se[dataset]

    * Append to results
    use `results', clear
    set obs `=_N+1'

	replace varname = "`v'" in L
    replace mean_trees = `mean1' in L
    replace mean_wb = `mean0' in L
    replace diff = `b' in L
	replace N_trees    = `N1' in L
    replace N_wb       = `N0' in L

    save `results', replace
}

use `results', clear
gen diff_star = string(diff, "%9.3f") 

* Keep only the columns you want, in order
keep varname mean_trees mean_wb diff_star N_trees N_wb

* Format numbers
format mean_trees mean_wb %9.3f

* Export to LaTeX
//close a handle left open by an interrupted earlier run before reopening
capture file close tex
file open tex using "$out_tab/comparison_table.tex", write replace
file write tex "\begin{tabular}{lccccc}" _n
file write tex "\hline" _n
file write tex "Variable & Mean  & Mean & Difference & Obs. & Obs. \\" _n
file write tex "& (TREES) & (WB) & & (TREES) & (WB) \\" _n
file write tex "\hline" _n

* Replace variable names with descriptive labels
replace varname = "Household size"                if varname == "npeople"
replace varname = "Number of females in household" if varname == "householdfemales"
replace varname = "Share of females in household"  if varname == "p_females"
replace varname = "Number of children"             if varname == "n_childrennumber"
replace varname = "Number of children under 5"     if varname == "n_childrenunder5"
replace varname = "Share of children under 5"      if varname == "p_child5yrs"
replace varname = "Highest education level"        if varname == "highest_sch_level"
replace varname = "Main income from farming"       if varname == "farm_main_income"
replace varname = "Farm size (ha)"           if varname == "farm_area"
replace varname = "Number of trees"                if varname == "trees"
replace varname = "Trees per hectare"              if varname == "trees_per_ha"
replace varname = "Number of livestock"            if varname == "tot_animal"
replace varname = "Livestock per hectare"          if varname == "tot_animal_ha"
replace varname = "Number of crops"                if varname == "n_crops"
replace varname = "HFIAS score"                    if varname == "hfias"
replace varname = "HFIAS category"                 if varname == "hfias_category"
replace varname = "Food insecure"               if varname == "hfiascatid4"
replace varname = "Worried about food"             if varname == "worriedfood"
replace varname = "No food in household"           if varname == "nofoodid"
replace varname = "Reduced number of meals"        if varname == "reduce_n_meals"


quietly {
    forvalues i = 1/`=_N' {
        file write tex ///
        "`=varname[`i']' & " ///
        "`=string(mean_trees[`i'],"%9.3f")' & " ///
        "`=string(mean_wb[`i'],"%9.3f")' & " ///
        "`=diff_star[`i']' & " ///
        "`=N_trees[`i']' & " ///
        "`=N_wb[`i']' \\" _n
    }
}

file write tex "\hline" _n
file write tex "\end{tabular}" _n
file close tex

end

* =========================================
* SINGLE- VS MULTI-SURVEY FARMS - baseline characteristics comparison
*   Compares the first-survey characteristics of farms surveyed once with
*   farms surveyed more than once (the latter restricted to the estimation
*   sample). Difference = coefficient on the multi-survey dummy, SE clustered
*   by project area; stars: * p<0.1, ** p<0.05, *** p<0.01.
*   -> $out_tab/comparison_single_multi.tex  (Supp. Table 12)
*   Note: in the file-write loop the row terminator `\\" _n` is written on the
*   difference column (the two observation-count columns are not written).
* =========================================
capture program drop single_multi_comparison
program define single_multi_comparison
eststo clear
	use "final/TREES_data_pseudo", clear

	//farm size of exactly 0 is non-response coded as zero, not a landless
	//farm: the per-ha variables are already missing for these farms, so
	//treat the size itself as missing too
	replace size=. if size==0
	replace trees=0 if missing(trees)
 

	//farm-level flag: a farm enters the baseline REGRESSION sample only if
	//BOTH its first and last survey pass the row filter. Farms with a single
	//valid row (174 farms whose last survey fails gap_ok) are dropped by
	//reghdfe as fe_farm singletons, so checking the filter on the first-survey
	//row alone yields 674 farms, not the 500 in the regressions.
	gen byte rowpass = insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1
	bysort fe_farm: egen n_valid = total(rowpass)

	//baseline (first survey) observation of every farm
	keep if first_obs==1

	//drop multi-survey farms that are NOT in the regression sample, so they
	//do not fall into the single-survey column
	keep if n_obs==1 | n_valid>=2

	//multi-survey indicator: 1 = farm in the baseline estimation sample
	gen byte multi=(n_obs>1)

	//keep only single-survey farms surveyed in the same years and the same
	//countries as the multi-survey farms' baseline surveys (all rows here are
	//first surveys, so year_survey is the baseline year)
	bysort year_survey: egen byte has_multi_yr = max(multi)
	bysort iso_trees:   egen byte has_multi_ct = max(multi)
	keep if multi==1 | (has_multi_yr==1 & has_multi_ct==1)
	drop has_multi_yr has_multi_ct


	//number of crop types (as in wb_comparison)
	gen n_crops=marketableproducts+foodcrops
	replace n_crops=0 if missing(n_crops)

	tempfile base
	save `base', replace

	local vars npeople size trees w1_trees_per_ha havelivestock tot_animal wtot_animal_ha n_crops hdds_id3 nofoodid hfiascatid4

	* Create empty results dataset
	clear
	set obs 0
	gen varname = ""
	gen mean_single = .
	gen mean_multi = .
	gen diff_star = ""
	gen N_single = .
	gen N_multi = .
	tempfile results
	save `results', replace

	* ------------------------------------------
	* Loop over variables
	* ------------------------------------------
	foreach v of local vars {

		use `base', clear

		quietly summarize `v' if multi==0 & `v' != .
		local mean0 = r(mean)
		local N0 = r(N)

		quietly summarize `v' if multi==1 & `v' != .
		local mean1 = r(mean)
		local N1 = r(N)

		//difference: coefficient on the multi dummy, SE clustered by project area
		quietly reg `v' multi if `v' != ., vce(cluster fe_proj_area)
		local b = _b[multi]
		local p = 2*ttail(e(df_r), abs(_b[multi]/_se[multi]))
		local stars ""
		if `p'<0.10 local stars "*"
		if `p'<0.05 local stars "**"
		if `p'<0.01 local stars "***"

		* Append to results
		use `results', clear
		set obs `=_N+1'
		replace varname     = "`v'" in L
		replace mean_single = `mean0' in L
		replace mean_multi  = `mean1' in L
		replace diff_star   = string(`b', "%9.3f") + "`stars'" in L
		replace N_single    = `N0' in L
		replace N_multi     = `N1' in L
		save `results', replace
	}

	use `results', clear

	* Replace variable names with descriptive labels
	replace varname = "Household Size"        if varname=="npeople"
	replace varname = "Farm Size (ha)"        if varname=="size"
	replace varname = "Number of Trees"       if varname=="trees"
	replace varname = "Trees per Hectare"     if varname=="w1_trees_per_ha"
	replace varname = "1(Own Livestock)"      if varname=="havelivestock"
	replace varname = "Number of Livestock"   if varname=="tot_animal"
	replace varname = "Livestock per Hectare" if varname=="wtot_animal_ha"
	replace varname = "N. crops"       if varname=="n_crops"	
	replace varname = "Low food diversity"    if varname=="hdds_id3"
	replace varname = "No Food"            if varname=="nofoodid"
	replace varname = "Food Insecure"      if varname=="hfiascatid4"

	* Export to LaTeX
	//close a handle left open by an interrupted earlier run before reopening
	capture file close tex
	file open tex using "$out_tab/comparison_single_multi.tex", write replace
	file write tex "\begin{tabular}{lccc}" _n
	file write tex "\hline" _n
	file write tex "Variable & Mean & Mean & Difference  \\" _n
	file write tex "& (1 survey) & (2+ surveys)  \\" _n
	file write tex "\hline" _n

	quietly {
		forvalues i = 1/`=_N' {
			file write tex ///
			"`=varname[`i']' & " ///
			"`=string(mean_single[`i'],"%9.3f")' & " ///
			"`=string(mean_multi[`i'],"%9.3f")' & " ///
			"`=diff_star[`i']' \\" _n
		}
	}

	file write tex "\hline" _n
	file write tex "\end{tabular}" _n
	file close tex

end
