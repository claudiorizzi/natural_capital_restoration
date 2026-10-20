*===============================================================================
* 03_main_estimates.do - estimation programs of
*   "Natural capital restoration and household wealth in sub-Saharan Africa"
*   (Rizzi, Xu, Heal). Replication package, github.com/claudiorizzi/natural_capital_restoration.
*
* This file only DEFINES the programs listed in `main` at the bottom; run them
* from code/run_demo.do (which sets the three path globals) or interactively:
*     global data_dir "<folder with final/TREES_data_pseudo.dta and final/EVI_farms_pseudo.dta>"
*     global out_tab  "<folder for LaTeX tables>"
*     global out_fig  "<folder for PNG figures>"
*     do code/03_main_estimates.do
*     baseline
*
* Programs and the outputs they produce (manuscript / Supplementary Information)
*   baseline        own_n_livestock.tex (Supp. Table 4), rob_alternative_fes.tex (Supp. Table 8)
*   leave_one_out   leave_one_out_animal.png (Supp. Fig. 5), leave_one_out.png (Supp. Fig. 4)
*   median_reg      rob_median_reg.tex (Supp. Table 10)
*   plots_wealth    wealth_post_coefficients.png (Fig. 2)
*   evi_analysis    evi_median_mean.tex (Supp. Table 7), evi_median_combined.png and
*                   evi_mean_combined.png (Fig. 1, panels a and b)
*   evi_robustness  evi_median_matchedcontrols.png, evi_mean_matchedcontrols.png
*                   (Supp. Fig. 2), evi_median_rand_inference.png (Supp. Fig. 3)
*   trees_crops     trees_crops.tex (Supp. Table 5), trees_crops_post_coefficients.png (Fig. 3)
*   channel         nlivestock_trees_crops.tex (Supp. Table 13)
*   health          health.tex (Supp. Table 6), health_post_coefficients.png (Fig. 4)
*   add_rob         rob_size.tex (Supp. Table 9), cross_farm_char_livestock.tex (Supp. Table 14)
*
* Household estimates: within-farm OLS / PPML with farm and survey-month fixed
* effects on the estimation sample (insample & n_obs>1 & gap_ok & gap_ok_fl:
* farms with a first survey within 12 months of project start and a later
* survey at least 12 months after it), standard errors clustered by project
* area (fe_proj_area).
* Satellite (EVI) estimates: stacked DiD (Wooldridge 2021 extended TWFE, jwdid)
* and Sun & Abraham (2021), never-treated control locations, event window
* [-3,+4] with t-1 as the reference period, standard errors clustered by
* project area (fe_proj_area) throughout. esttab reports t statistics in
* parentheses.
*
* Required packages (SSC): ftools reghdfe hdfe ppmlhdfe jwdid eventstudyinteract
*   avar estout erepost distinct winsor2 qreg2 psmatch2
*===============================================================================

*** Preamble ***
clear all
set more off
set matsize 11000
set maxvar 32700

//data and output folders are set by code/run_demo.do (globals data_dir, out_tab, out_fig)
cd "$data_dir"

global controls i.wb_adult_mostedu i.wb_gender wb_w_npeople wb_w_size wb_w_participantage
global gcpc w5_gcpc_2017usd
//EVI event-time indicators: treat_yearidN = event year N-8 (id5 = -3 ... id12 = +4);
//id7 (t = -1) is the omitted reference period
global treat_years treat_yearid5 treat_yearid6 treat_yearid8 treat_yearid9 treat_yearid10 treat_yearid11 treat_yearid12

*-------------------------------------------------------------------------------
* Helpers to harvest jwdid event-study output (used by the EVI programs below)
*-------------------------------------------------------------------------------
* jwdid_evmat: after `jwdid ... ; estat event`, capture `matrix __ETAB = r(table)`
* and call jwdid_evmat. It builds matrix EVMAT with columns (event_year b1 V1),
* one row per relative period returned by jwdid. Leaves the data in memory
* untouched; the caller then does `clear' / `svmat EVMAT, names(col)'.
* Column names of r(table) look like "r2vs1._at@2.__event__"; stripping the
* prefixes leaves the __event__ level, whose value label is the relative period
* (text, parsed with real(); columns that do not parse are skipped and reported).
capture program drop jwdid_evmat
program define jwdid_evmat
    capture matrix drop EVMAT
    local coln : colname __ETAB
    local k = 0
    foreach name of local coln {
        local ++k
        //pull the __event__ level integer, ignoring base/omitted (#b/#o/#bn) factor markers
        if !regexm("`name'", "([0-9]+)[bno]*\.__event__") {
            display as text "jwdid_evmat: skipped column `k' (|`name'|) - not an __event__ level"
            continue
        }
        local i = regexs(1)
        local et : label (__event__) `i'
        //the label is TEXT: parse it with real() so a non-event label is skipped
        //instead of raising a type mismatch inside the matrix expression below
        local ev = real("`et'")
        if missing(`ev') {
            display as text "jwdid_evmat: skipped column `k' (|`i'|, label |`et'|) - not an event period"
            continue
        }
        matrix EVMAT = nullmat(EVMAT) \ (`ev', __ETAB[1,`k'], (__ETAB[2,`k'])^2)
    }
    matrix colnames EVMAT = event_year b1 V1
end

* jwdid_ev_relabel: run right after `estat event, post`. Renames the posted
* event-study ATTs to treat_yearid# so they line up, by event time, with the
* Sun & Abraham columns (which use $treat_years) in the combined esttab table.
* POSITION-BASED: the estimation sample is the contiguous [-3,+4] window, so
* e(b) column k is event time k-4 (col 3 = t=-1 reference, coef 0). treat_yearid
* = event+8 (id5=-3 ... id8=0 ... id12=+4). The t=-1 column (treat_yearid7) is
* dropped so t=-1 is the reference, matching SA ($treat_years skips id7).
capture program drop jwdid_ev_relabel
program define jwdid_ev_relabel
    tempname b V b2 V2
    matrix `b' = e(b)
    matrix `V' = e(V)
    //hard guard: the position->event mapping below is only valid for the full
    //8-column [-3,+4] window. Fail loudly rather than mislabel.
    local K = colsof(`b')
    if `K' != 8 {
        display as error "jwdid_ev_relabel: expected 8 event columns ([-3,+4] window), got `K'"
        exit 498
    }
    local newnames ""
    local keepidx ""
    forvalues k = 1/`K' {
        local ev = `k' - 4
        if `ev' != -1 {
            local newnames `newnames' treat_yearid`=`ev'+8'
            local keepidx `keepidx' `k'
        }
    }
    * Keep the selected columns/rows. Stata's own matrix subscripts do NOT accept a
    * comma-separated index list, so do the non-contiguous selection in Mata.
    mata: idx = strtoreal(tokens(st_local("keepidx"))); st_matrix(st_local("b2"), st_matrix(st_local("b"))[1, idx]); st_matrix(st_local("V2"), st_matrix(st_local("V"))[idx, idx])
    matrix colnames `b2' = `newnames'
    matrix colnames `V2' = `newnames'
    matrix rownames `V2' = `newnames'
    * Post via -erepost- (Ben Jann): a fresh -ereturn post- that keeps the other
    * e() results, so the b/V dimension may shrink from 8 to 7 (official
    * -ereturn repost- forbids dimension changes).
    erepost b = `b2' V = `V2', rename
end

* =========================================
* Household outcomes of the intervention: livestock ownership (OLS) and number
* of livestock (PPML), within-farm design
*   -> $out_tab/own_n_livestock.tex      (Supp. Table 4)
*   -> $out_tab/rob_alternative_fes.tex  (Supp. Table 8: alternative fixed effects,
*      single-survey farms, livestock per hectare, herd value)
* =========================================
capture program drop baseline
program define baseline

	eststo clear
	use "final/TREES_data_pseudo", clear

	//GCP in thousands
	replace w5_gcpc_2017usd=w5_gcpc_2017usd/1000

	//farm size of 0 is non-response coded as zero: treat as missing in the
	//wb_w_size control of the alternative-FEs table
	replace size=. if size==0
	replace wb_w_size=. if wb_w_size==0
		
	eststo m1: quietly reghdfe havelivestock post if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local estimator "OLS"
			 estadd local contr "Y"
			 estadd local farmfe "Y" 	
			 estadd local monthfe "Y"
			 estadd local post2019 "N"
			 estadd scalar myr2 e(r2)
			 
	eststo m2: quietly reghdfe havelivestock post $gcpc if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1 , absorb( month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local estimator "OLS"
			 estadd local contr "Y"
			 estadd local farmfe "Y" 	
			 estadd local monthfe "Y"
			 estadd local post2019 "Y"
			 estadd scalar myr2 e(r2)
			 			 
	*-------------
	//Number of animals PPML
	eststo m3: quietly ppmlhdfe tot_animal post if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1 , absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local estimator "PPML"
			 estadd local contr "Y"
			 estadd local farmfe "Y" 	
			 estadd local monthfe "Y"
			 estadd local myr2 "-"
			 
	eststo m4: quietly ppmlhdfe tot_animal post $gcpc  if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey  fe_farm) vce(cluster fe_proj_area)
			 estadd local estimator "PPML"
			 estadd local contr "Y"
			 estadd local farmfe "Y" 	
			 estadd local monthfe "Y"
			 estadd local myr2 "-"
	
	esttab m* using  "$out_tab/own_n_livestock.tex" , nonotes nobaselevels noconstant noomitted width(\hsize) compress booktabs label interaction(" $\times$ ") replace  nogaps b(3) r2(2) nomtitles ///
		hlinechar(`=char(151)')   ///
		prehead("\begin{tabular}{l*{4}{c}}" ///
            "\toprule" ///
            "Dep. var.: & \multicolumn{2}{c}{Own livestock} & \multicolumn{2}{c}{N. livestock}  \\" ///
            "\cmidrule(lr){2-3}\cmidrule(lr){4-5}")	///
			 drop(_cons) /// 
		 coeflabels(w5_gcpc_2017usd "GCP per capita (\textdollar 1,000)" post "1(Post 1st survey)") ///
		star(* 0.1 ** 0.05 *** .01)    ///
		stats(estimator farmfe monthfe N myr2, fmt(0 0 0 %12.0gc 3) labels("Estimator" "Farm F.E." "Month F.E." "Observations" "R\$^2\$") ) 


	*------------
	//Robustness: alternative F.E.s specs
	eststo clear
	global controls wb_adult_mostedu wb_gender wb_w_npeople wb_w_size wb_w_participantage
	
	//define treatment from non-multi-observation farms has a late survey, i.e., after 1 year.
	gen late_365=(time_from_start>365)
	label var late_365 "1(Later Survey)"

	//own livestock	 	
	eststo m1: quietly reghdfe havelivestock late_365 $gcpc $controls if gap_ok==1 & gap_ok_fl==1, absorb(fe_year fe_country#month_survey fe_proj_area) vce(cluster fe_proj_area)
			 estadd local sample "All"
			 estadd local estimator "OLS"
			 estadd local contr "Y"
			 estadd local farmfe "N"
			 estadd local monthfe "N"
			 estadd local cityfe "Y" 
			 estadd local ctrymofe "Y"
			 estadd local yearfe "Y" 
			 estadd scalar myr2 e(r2)

	//own livestock, no double survey 	
	eststo m2: quietly reghdfe havelivestock late_365 $gcpc $controls if gap_ok==1 & gap_ok_fl==1 & n_obs==1, absorb(fe_year fe_country#month_survey fe_proj_area) vce(cluster fe_proj_area)
			 estadd local sample "1 survey"
			 estadd local estimator "OLS"
			 estadd local contr "Y"
			 estadd local farmfe "N"
			 estadd local monthfe "N"
			 estadd local cityfe "Y" 
			 estadd local ctrymofe "Y"
			 estadd local yearfe "Y" 
			 estadd scalar myr2 e(r2)
			 
	//N. livestock (PPML)	 
	eststo m3: quietly 	ppmlhdfe tot_animal late_365 $gcpc $controls if gap_ok==1 & gap_ok_fl==1, absorb( fe_year fe_country#month_survey fe_proj_area) vce(cluster fe_proj_area)
			 estadd local sample "All"
			 estadd local estimator "PPML"
			 estadd local contr "Y"
			 estadd local farmfe "N"
			 estadd local monthfe "N"
			 estadd local cityfe "Y" 
			 estadd local ctrymofe "Y"
			 estadd local yearfe "Y" 
			 estadd local myr2 "-"

	//N. livestock (PPML), no double survey	 
	eststo m4: quietly ppmlhdfe tot_animal late_365 $gcpc $controls if gap_ok==1 & gap_ok_fl==1 & n_obs==1, absorb( fe_year fe_country#month_survey fe_proj_area) vce(cluster fe_proj_area)
			 estadd local sample "1 survey"
			 estadd local estimator "PPML"
			 estadd local contr "Y"
			 estadd local farmfe "N"
			 estadd local monthfe "N"
			 estadd local cityfe "Y" 
			 estadd local ctrymofe "Y"
			 estadd local yearfe "Y" 
			 estadd local myr2 "-"
				 
	//total animal per hectare
	eststo m5: quietly ppmlhdfe wtot_animal_ha post $gcpc  if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb( month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local sample "$>$ 2 surveys"
			 estadd local estimator "PPML"
			 estadd local contr "Y"
			 estadd local farmfe "Y"
			 estadd local monthfe "Y"
			 estadd local cityfe "N" 
			 estadd local ctrymofe "N"
			 estadd local yearfe "N" 
			 estadd local myr2 "-"

	// Animal value	 
	eststo m6: quietly ppmlhdfe w_animals_value post $gcpc if iso_trees!="SEN" & insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb( month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local sample "$>$ 2 surveys"
			 estadd local estimator "PPML"
			 estadd local contr "Y"
			 estadd local farmfe "Y"
			 estadd local monthfe "Y"
			 estadd local cityfe "N" 
			 estadd local ctrymofe "N"
			 estadd local yearfe "N" 
			 estadd local myr2 "-"

			 
	esttab m* using  "$out_tab/rob_alternative_fes.tex" , nonotes nobaselevels noconstant noomitted width(\hsize) compress booktabs label interaction(" $\times$ ") replace nogaps b(3) r2(2) nomtitles ///
		hlinechar(`=char(151)')   ///
	   prehead("\begin{tabular}{l*{6}{c}}" ///
            "\toprule" ///
            "Dep. var.: & \multicolumn{2}{c}{1(Own Livestock)} & \multicolumn{2}{c}{N. livestock} & N. Livestock/ha & Herd Value \\" "\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-6}\cmidrule(lr){7-7}")	///
		 keep(post late_365 $gcpc $controls) /// 
		 coeflabels(w5_gcpc_2017usd "GCP per capita (\textdollar 1,000)" post "1(Post 1st survey)") ///
		star(* 0.1 ** 0.05 *** .01)    ///
		stats(estimator sample contr cityfe farmfe monthfe ctrymofe yearfe N myr2, fmt(0 0 0 0 0 0 0 0 %12.0gc 3) labels("Estimator" "Sample" "Controls" "City/Province F.E." "Farm F.E." "Month F.E." "Country-month F.E." "Year F.E."  "Observations" "R\$^2\$") ) 
end

* =========================================
* Robustness: leave one animal type out, leave one country out
*   -> $out_fig/leave_one_out_animal.png  (Supp. Fig. 5)
*   -> $out_fig/leave_one_out.png         (Supp. Fig. 4)
* Step 1 builds country-specific unit values for each animal type (herd value).
* =========================================
capture program drop leave_one_out
program define leave_one_out
	eststo clear
	use "final/TREES_data_pseudo", clear
	//GCP in thousands (same rescaling as in baseline)
	replace w5_gcpc_2017usd=w5_gcpc_2017usd/1000

	foreach animal in chicken cattle mule goat sheep {
		gen price_`animal' = .
		replace price_`animal' = 5     if iso_trees == "KEN" & "`animal'" == "chicken"
		replace price_`animal' = 619   if iso_trees == "KEN" & "`animal'" == "cattle"
		replace price_`animal' = 116   if iso_trees == "KEN" & "`animal'" == "mule"
		replace price_`animal' = 154   if iso_trees == "KEN" & "`animal'" == "goat"
		replace price_`animal' = 46    if iso_trees == "KEN" & "`animal'" == "sheep"

		replace price_`animal' = 10    if iso_trees == "TZA" & "`animal'" == "chicken"
		replace price_`animal' = 440   if iso_trees == "TZA" & "`animal'" == "cattle"
		replace price_`animal' = 220   if iso_trees == "TZA" & "`animal'" == "mule"
		replace price_`animal' = 52    if iso_trees == "TZA" & "`animal'" == "goat"
		replace price_`animal' = 48    if iso_trees == "TZA" & "`animal'" == "sheep"

		replace price_`animal' = 10    if iso_trees == "UGA" & "`animal'" == "chicken"
		replace price_`animal' = 279   if iso_trees == "UGA" & "`animal'" == "cattle"
		replace price_`animal' = 935   if iso_trees == "UGA" & "`animal'" == "mule"
		replace price_`animal' = 58    if iso_trees == "UGA" & "`animal'" == "goat"
		replace price_`animal' = 52    if iso_trees == "UGA" & "`animal'" == "sheep"
		
	}

	* --- Step 2: Full (non-excluded) totals for reference ---

	gen tot_animals = n_chickencount + n_cattlecount + n_mulecount + n_goatcount + n_sheepcount

	gen animals_value_full = price_chicken * n_chickencount + price_cattle * n_cattlecount ///
		+ price_mule * n_mulecount + price_goat * n_goatcount + price_sheep * n_sheepcount

	* --- Step 3: Leave-one-out variables and regressions ---

	local animals "chicken cattle mule goat sheep"
	local counts  "n_chickencount n_cattlecount n_mulecount n_goatcount n_sheepcount"

	tempname results
	tempfile results_file
	postfile `results' str20 outcome str10 excluded double(b se lo hi) using `results_file'

	local n_animals : word count `animals'

	forvalues i = 1/`n_animals' {
		local animal : word `i' of `animals'
		local count  : word `i' of `counts'

		* Leave-one-out total count
		cap drop loo_tot_animals
		gen loo_tot_animals = tot_animals - `count'

		* Leave-one-out total count per hectare)
		cap drop loo_tot_animals_ha
		gen temp = loo_tot_animals / size
		winsor2 temp, cuts(1 99)
		rename temp_w loo_tot_animals_ha
		drop temp 

		* Leave-one-out value (country-specific)
		cap drop loo_animals_value
		gen loo_animals_value = animals_value_full - price_`animal' * `count'

		* --- Regressions ---
		foreach y in loo_tot_animals loo_tot_animals_ha loo_animals_value {
			qui ppmlhdfe `y' post $gcpc if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, ///
				absorb(month_survey fe_farm) vce(cluster fe_proj_area)

			local b  = _b[post]
			local se = _se[post]
			local lo = `b' - 1.645 * `se'
			local hi = `b' + 1.645 * `se'

			post `results' ("`y'") ("`animal'") (`b') (`se') (`lo') (`hi')
		}
	}

	postclose `results'

	* --- Step 4: Plot ---

	preserve
	use `results_file', clear

	* Numeric outcome axis
	gen outcome_n = 1 if outcome == "loo_tot_animals"
	replace outcome_n = 2 if outcome == "loo_tot_animals_ha"
	replace outcome_n = 3 if outcome == "loo_animals_value"

	* Jitter by excluded animal for visibility
	gen     jitter = -0.20 if excluded == "chicken"
	replace jitter = -0.10 if excluded == "cattle"
	replace jitter =  0.00 if excluded == "mule"
	replace jitter =  0.10 if excluded == "goat"
	replace jitter =  0.20 if excluded == "sheep"

	gen x = outcome_n + jitter

	* Marker symbols by excluded animal
	twoway (rcap hi lo x if excluded == "chicken", lcolor(navy))       ///
		   (scatter b x     if excluded == "chicken", mc(navy)    ms(O))  ///
		   (rcap hi lo x if excluded == "cattle",  lcolor(cranberry))  ///
		   (scatter b x     if excluded == "cattle",  mc(cranberry) ms(D))  ///
		   (rcap hi lo x if excluded == "mule",    lcolor(forest_green)) ///
		   (scatter b x     if excluded == "mule",    mc(forest_green) ms(T))  ///
		   (rcap hi lo x if excluded == "goat",    lcolor(dkorange))   ///
		   (scatter b x     if excluded == "goat",    mc(dkorange) ms(S))  ///
		   (rcap hi lo x if excluded == "sheep",   lcolor(purple))     ///
		   (scatter b x     if excluded == "sheep",   mc(purple)  ms(+))  ///
		   , xscale(range(0.5 3.5))                                    ///
		   xlabel(1 "N. livestock" 2 "N. livestock/ha" 3 "Herd value", labsize(small) nogrid) ///
		   ytitle("Change relative to first survey")           ///
		   xtitle("")                                                   ///
		   yline(0, lpattern(dash) lcolor(gs8))                        ///
		   legend(order(2 "Excl. chicken" 4 "Excl. cattle" 6 "Excl. mule" 8 "Excl. goat" 10 "Excl. sheep") ///
				  rows(1) size(vsmall) position(6))                     ///
		   graphregion(color(white))
	graph export "$out_fig/leave_one_out_animal.png", replace

	restore

	* --- Cleanup ---
	cap drop price_* tot_animals animals_value_full
	cap drop loo_tot_animals loo_tot_animals_ha loo_animals_value
 	
	
	*-------------
	//robustness: leave one country out
	* Set up a postfile to store results
	tempname results
	tempfile results_file
	postfile `results' str20 outcome str10 excl_country double(b se lo hi) using `results_file'

	* Get list of countries to exclude one at a time
	levelsof iso_trees, local(countries)

	* Loop over outcomes and excluded countries
	foreach y in havelivestock tot_animal wtot_animal_ha { 
		foreach c of local countries {
			if "`y'"=="havelivestock" {
				qui reghdfe `y' post $gcpc if iso_trees != "`c'" & insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, ///
				absorb(month_survey fe_farm) vce(cluster fe_proj_area)
				}
				else {
					if "`c'"=="TZA" {
					ppmlhdfe wtot_animal post $gcpc if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1 & iso_trees!="TZA", absorb(month_survey  fe_farm) vce(cluster fe_proj_area)
					}
					
					else {
						qui ppmlhdfe `y' post $gcpc if iso_trees != "`c'" & insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)

					}
				}
			
			local b  = _b[post]
			local se = _se[post]
			local lo = `b' - 1.645*`se'   // 90% CI
			local hi = `b' + 1.645*`se'
			
			post `results' ("`y'") ("`c'") (`b') (`se') (`lo') (`hi')
		}
	}
	postclose `results'

	* Load results and plot
	preserve
	use `results_file', clear

	* Create numeric x-axis for outcomes
	encode outcome, gen(outcome_n)

	* one point per (outcome, excluded country)

	twoway (rcap hi lo outcome_n, lcolor(navy)) ///
		   (scatter b outcome_n, mcolor(navy) msymbol(O)), ///
		   xscale(range(0.5 3.5)) ///
	xlabel(1 "1(Own Livestock)" 2 "N. livestock" 3 "N. Livestock per ha", nogrid)    ytitle("Change relative to first survey") ///
		xtitle("Outcome variable") ///
		yline(0, lpattern(dash) lcolor(gs8)) ///
		legend(off) 
	graph export "$out_fig/leave_one_out.png", replace
	restore
end

* =========================================
* Robustness: median regressions of first-to-last-survey changes
*   -> $out_tab/rob_median_reg.tex  (Supp. Table 10)
* =========================================
capture program drop median_reg
program define median_reg
	eststo clear
	use "final/TREES_data_pseudo", clear
	replace ch_w5_gcpc_2017usd=ch_w5_gcpc_2017usd*100

	//first-to-last changes (last survey - first survey, same convention as the
	//ch_* variables built in 3_put_files_together); missing for single-survey farms
	foreach var of varlist w_trees w_marketableproducts w_foodcrops {
		bysort fe_farm (time_from_start): gen ch_`var'=`var'[_N]-`var'[1] if _N>1
	}

	eststo m1: quietly qreg2 ch_tot_animal ch_w5_gcpc_2017usd if last_obs==1 & insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, quantile(0.5) cluster(fe_proj_area)
	
	eststo m2: quietly qreg2 wch_tot_animal_ha ch_w5_gcpc_2017usd if last_obs==1 & insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, quantile(0.5) cluster(fe_proj_area)
	
	eststo m3: quietly qreg2 ch_w_trees ch_w5_gcpc_2017usd if last_obs==1 & insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, quantile(0.5) cluster(fe_proj_area)

	eststo m4: quietly qreg2 wch_trees_per_ha ch_w5_gcpc_2017usd if last_obs==1 & insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, quantile(0.5) cluster(fe_proj_area)
	


	
	esttab m* using  "$out_tab/rob_median_reg.tex" , nonotes  noomitted width(\hsize) compress booktabs label interaction(" $\times$ ") replace nogaps b(3) r2(2) nomtitles ///
		hlinechar(`=char(151)')   ///
		keep( ch_w5_gcpc_2017usd _cons) ///
		prehead("\begin{tabular}{l*{4}{c}}" ///
            "\toprule" ///
            "Dep. var.: $\Delta$ &  N. livestock & N. livestock/ha & N. trees & N. trees/ha  \\" ///
			"\cmidrule(lr){2-2}\cmidrule(lr){3-3}\cmidrule(lr){4-4}\cmidrule(lr){5-5}\cmidrule(lr){6-6}")	///
		 coeflabels( _cons "Median $\Delta$, last vs. 1st survey" ch_w5_gcpc_2017usd "$\Delta$ GCP per capita (\%)") ///
		star(* 0.1 ** 0.05 *** .01)    ///
		stats( N, fmt(0 %12.0gc) labels("Observations" ) ) 
	

end

* =========================================
* Post coefficients across the wealth outcomes (own livestock, N. livestock,
* N. livestock/ha, herd value)
*   -> $out_fig/wealth_post_coefficients.png  (Fig. 2)
* =========================================
capture program drop plots_wealth
program define plots_wealth
	*-----------------------
	// Post-indicator wealth effects: single coefficient plot across outcomes.
	// Reports the coefficient on the treatment "post" indicator from:
	//   havelivestock    reghdfe
	//   tot_animal       ppmlhdfe
	//   wtot_animal_ha   ppmlhdfe
	//   w_animals_value  ppmlhdfe  (drops SEN, no price data)
	// NOTE: havelivestock (LPM) is a probability-point effect; the ppmlhdfe
	// outcomes are proportional (log-point) effects, so the y-axis mixes scales.
	// To make that explicit, each coefficient (x100) is written above its CI
	// with its unit: "pp" (percentage points) for the LPM, "log p." for PPML.
	eststo clear
	use "final/TREES_data_pseudo", clear
	replace w5_gcpc_2017usd=w5_gcpc_2017usd/1000
	tempname results
	tempfile results_file
	postfile `results' str20 outcome double(b se lo hi) using `results_file'

	local outcomes "havelivestock tot_animal wtot_animal_ha w_animals_value"

	foreach y of local outcomes {

		if "`y'"=="havelivestock" {
			reghdfe `y' post $gcpc if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
		}
		else if "`y'"=="wtot_animal_ha" {
			ppmlhdfe `y' post $gcpc if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
		}
		else if "`y'"=="w_animals_value" {
			ppmlhdfe `y' post $gcpc if iso_trees!="SEN" & insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
		}
		else {
			ppmlhdfe `y' post $gcpc  if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
		}

		local b  = _b[post]
		local se = _se[post]
		local lo = `b' - 1.645*`se'
		local hi = `b' + 1.645*`se'
		post `results' ("`y'") (`b') (`se') (`lo') (`hi')
	}
	postclose `results'

	* ---- Plot: coefficient on post, one point per outcome ----
	preserve
	use `results_file', clear

	gen outcome_n = .
	replace outcome_n = 1 if outcome=="havelivestock"
	replace outcome_n = 2 if outcome=="tot_animal"
	replace outcome_n = 3 if outcome=="wtot_animal_ha"
	replace outcome_n = 4 if outcome=="w_animals_value"
	sort outcome_n

	//coefficient label shown above each CI: signed coefficient x100 with its unit
	//(pp for the LPM column, log points for the PPML columns). The explicit "+"
	//prevents a quick reader from mistaking sub-1 PPML coefficients for declines.
	gen blabel = cond(b>=0, "+", "") + string(100*b, "%9.1f") + cond(outcome=="havelivestock", " pp", " log p.")

	twoway ///
		(rcap hi lo outcome_n, lcolor(navy)) ///
		(scatter b outcome_n, mcolor(navy) msymbol(O)) ///
		(scatter hi outcome_n, msymbol(none) mlabel(blabel) mlabposition(12) mlabcolor(navy) mlabsize(small)), ///
		xscale(range(0.5 4.5)) ///
		xlabel(1 "Own livestock" 2 "N. livestock" 3 "N. livestock/ha" 4 "Herd value", labsize(small) nogrid) ///
		xtitle("") ///
		ytitle("Change relative to first survey") ///
		yline(0, lpattern(dash) lcolor(gs8)) ///
		legend(off) ///
		graphregion(color(white)) plotregion(color(white))
	graph export "$out_fig/wealth_post_coefficients.png", replace
	restore
end

* =========================================
* Satellite analysis: EVI event studies, treated farms vs never-treated
* control locations, window [-3,+4], reference t-1
*   stacked DiD (jwdid, Wooldridge 2021 ETWFE) and Sun & Abraham (2021);
*   farm and country-year (median) / cell-year (mean) fixed effects, GCP control,
*   standard errors clustered by project area
*   -> $out_tab/evi_median_mean.tex       (Supp. Table 7)
*   -> $out_fig/evi_median_combined.png   (Fig. 1a)
*   -> $out_fig/evi_mean_combined.png     (Fig. 1b)
* =========================================
capture program drop evi_analysis
program define evi_analysis

	eststo clear
	use "final/EVI_farms_pseudo", clear
	gen event_window=year-calyear_treat if calyear_treat!=0
	keep if (event_window>=-3 & event_window<=4) | calyear_treat==0
	egen fe_cellyr=group(fe_cell year)


	distinct forestgardenid if !missing(w1_evi_median) & !missing(gcpc_2017usd_w1) 
	
	distinct forestgardenid if !missing(w1_evi_median) & !missing(gcpc_2017usd_w1) & !missing(treat_year)

	distinct forestgardenid if !missing(w1_evi_median) & !missing(gcpc_2017usd_w1) & missing(treat_year)
		
	//Median stacked DID (Wooldridge 2021 ETWFE via jwdid)
	egen fe_ctry_yrcoh=group(projunit_calyear_treat iso_trees)
	//reghdfe helper: defines the estimation sample used by jwdid and Sun & Abraham
	quietly reghdfe w1_evi_median $treat_years gcpc_2017usd_w1, absorb(fe_farm fe_ctryyr) vce(cluster fe_proj_area)
			 cap drop insample_evi_med
			 gen byte insample_evi_med = e(sample)
			 local model_r2=e(r2)
			 local model_n=e(N)

	//jwdid stacked DiD; estat event, post -> event-study ATTs; relabel to treat_yearid# for the table
	
	quietly jwdid w1_evi_median gcpc_2017usd_w5 if insample_evi_med==1, ivar(fe_farm) tvar(year) gvar(calyear_treat) fevar(fe_ctryyr) method(reghdfe) cluster(fe_proj_area) never	
	quietly estat event, post
	jwdid_ev_relabel
	eststo m1
			 estadd local stacked_did "Y"
			 estadd local Sun_Abraham "N"
			 estadd local contr "Y"
			 estadd local farmfe "Y"
			 estadd local cyfe "Y"
			 estadd scalar myr2 `model_r2'
			 estadd scalar myn `model_n'
			 
	//Median Sun and Abrahams
	xtset fe_farm year
	quietly	eventstudyinteract w1_evi_median $treat_years if insample_evi_med==1, cohort(calyear_treat) control_cohort(never_treat) absorb(fe_farm fe_cellyr) vce(cluster fe_proj_area)
	matrix b_iw = e(b_iw)
	matrix V_iw = e(V_iw)
	erepost b = b_iw V = V_iw 
	estimates store m2
	estadd local stacked_did "N"
	estadd local Sun_Abraham "Y"
	estadd local contr "Y"
	estadd local farmfe "Y" 	
	estadd local cyfe "Y" 
			 estadd scalar myr2 `model_r2'
			 estadd scalar myn `model_n'

	///EVI MEAN
	//reghdfe helper: defines the estimation sample used by jwdid and Sun & Abraham
	quietly reghdfe w1_evi_mean $treat_years gcpc_2017usd_w5, absorb(fe_farm fe_ctryyr ) vce(cluster fe_proj_area)
			 cap drop insample_evi_mean
			 gen byte insample_evi_mean = e(sample)
			local model_r2=e(r2)
			local model_n=e(N)

	//jwdid stacked DiD (Wooldridge 2021 ETWFE)
	quietly jwdid w1_evi_mean gcpc_2017usd_w5 if insample_evi_mean==1 , ivar(fe_farm) tvar(year) gvar(calyear_treat)  method(reghdfe) fevar(fe_cellyr) cluster(fe_proj_area) never
	quietly 	estat event, post
	jwdid_ev_relabel
	eststo m3
			 estadd local stacked_did "Y"
			 estadd local Sun_Abraham "N"
			 estadd local contr "Y"
			 estadd local farmfe "Y"
			 estadd local cyfe "Y"
			 estadd scalar myr2 `model_r2'
			 estadd scalar myn `model_n'			 
	//Sun and Abraham DID
	replace w5_evi_median=w5_evi_median*100
	replace w1_evi_mean=w5_evi_median
	quietly	eventstudyinteract w1_evi_mean $treat_years if insample_evi_mean==1  , cohort(calyear_treat) control_cohort(never_treat) absorb(fe_farm fe_cellyr) vce(cluster fe_proj_area)
	matrix b_iw = e(b_iw)
	matrix V_iw = e(V_iw)
	erepost b = b_iw V = V_iw 
	estimates store m4
			 estadd local stacked_did "N"
			 estadd local Sun_Abraham "Y"
			 estadd local contr "Y"
			 estadd local farmfe "Y" 	
			 estadd local cyfe "Y" 		
			 estadd scalar myr2 `model_r2'
			 estadd scalar myn `model_n'			
	//coef names: treat_yearidN = event N-8 (t=-1/id7 omitted as reference in all 4 models)
	esttab m* using  "$out_tab/evi_median_mean.tex" , nonotes nobaselevels noconstant noomitted width(\hsize) compress booktabs label interaction(" $\times$ ") replace  nogaps b(3) r2(2) nomtitles ///
		hlinechar(`=char(151)')   ///
		    prehead("\begin{tabular}{l*{4}{c}}" ///
            "\toprule" ///
            "Dep. var.: & \multicolumn{2}{c}{EVI Median} & \multicolumn{2}{c}{EVI Mean}  \\" ///
            "\cmidrule(lr){2-3}\cmidrule(lr){4-5}")	///
			coeflabels(treat_yearid5 "Treat. $\times$ 1(Year -3)" treat_yearid6 "Treat. $\times$ 1(Year -2)" treat_yearid8 "Treat. $\times$ 1(Year 0)" treat_yearid9 "Treat. $\times$ 1(Year 1)" treat_yearid10 "Treat. $\times$ 1(Year 2)" treat_yearid11 "Treat. $\times$ 1(Year 3)" treat_yearid12 "Treat. $\times$ 1(Year 4)"		) ///
		star(* 0.1 ** 0.05 *** .01)    ///
		stats(stacked_did Sun_Abraham contr farmfe cyfe myn myr2, fmt(0 0 0 0 0 %12.0gc 3) labels("Stacked DiD" "Sun and Abraham (2021)" "Controls" "Farm F.E." "Country-Year F.E." "Observations" "R\$^2\$") ) 
			 

	*===============================================================
	* Fig. 1: EVI event study - two panels
	*   Panel a: Median EVI - Stacked DiD (jwdid) + Sun & Abraham overlaid
	*   Panel b: Mean   EVI - Stacked DiD (jwdid) + Sun & Abraham overlaid
	*   90% CIs; reference period t = -1; treatment onset marked at 0
	*===============================================================
	local z = 1.645

	foreach m in median mean {

		* ---------- (1) Stacked DiD (jwdid, Wooldridge 2021 ETWFE) ----------
		use "final/EVI_farms_pseudo", clear
		gen event_window=year-calyear_treat if calyear_treat!=0
		egen fe_cellyr=group(fe_cell year)
		keep if (event_window>=-3 & event_window<=4) | calyear_treat==0   // same window as the table
		if "`m'"=="median" {
			quietly jwdid w1_evi_`m' gcpc_2017usd_w5, ivar(fe_farm) tvar(year) gvar(calyear_treat) fevar(fe_ctryyr) method(reghdfe) cluster(fe_proj_area) never
			} 
			else {
				quietly jwdid w1_evi_`m' gcpc_2017usd_w5, ivar(fe_farm) tvar(year) gvar(calyear_treat)  method(reghdfe) fevar(fe_cellyr) cluster(fe_proj_area) never
			}
		quietly estat event, post
		* Event-study ATTs straight from e(b)/e(V). The kept window is the contiguous
		* [-3,+4], so column/row j -> event time j-4 (row 3 = t=-1, the normalized
		* reference: coef 0, variance 0 -> plots as a dot with no CI).
		matrix b = e(b)'
		matrix V = e(V)
		local K = rowsof(b)
		//hard guard: the position->event mapping (_n - 4) is only valid for the
		//full 8-column [-3,+4] window. Fail loudly rather than mislabel the plot.
		if `K' != 8 {
			display as error "expected 8 event columns ([-3,+4] window), got `K'"
			exit 498
		}
		clear
		set obs `K'
		gen event_year = _n - 4
		gen b1 = .
		gen V1 = .
		forvalues i = 1/`K' {
			replace b1 = b[`i',1]   in `i'
			replace V1 = V[`i',`i'] in `i'
		}
		replace b1 = 0 if event_year==-1
		replace V1 = 0 if event_year==-1
		sort event_year
		gen estid = 1
		tempfile sd_`m'
		save `sd_`m''

		* ---------- (2) Sun & Abraham (2021) ----------
		eststo clear
		use "final/EVI_farms_pseudo", clear
		gen event_window=year-calyear_treat if calyear_treat!=0
		replace w5_evi_median=w5_evi_median*100
		replace w1_evi_mean=w5_evi_median

		keep if (event_window>=-3 & event_window<=4) | calyear_treat==0   // same window as the table
		cap drop fe_cellyr
		egen fe_cellyr=group(fe_cell year)
		//recover the stacked-DiD estimation sample so SA uses the exact same observations
		qui reghdfe w1_evi_`m' ib0.($treat_years), absorb(fe_farm fe_cellyr) vce(cluster fe_proj_area)
		cap drop insample_evi
		gen byte insample_evi = e(sample)
		xtset fe_farm year
		if "`m'"=="median" {
		eventstudyinteract w1_evi_`m' $treat_years if insample_evi==1, cohort(calyear_treat) control_cohort(never_treat) absorb(fe_farm fe_cellyr) vce(cluster fe_proj_area)
		} 
		else {
			quietly	eventstudyinteract w1_evi_`m' $treat_years if insample_evi==1, cohort(calyear_treat) control_cohort(never_treat) absorb(fe_farm fe_cellyr) vce(cluster fe_proj_area)
		}
		* Sun & Abraham IW estimates straight from e(b_iw)/e(V_iw). Each column is a
		* treat_yearidN indicator with event time N-8 (id7 = -1 is the omitted
		* reference, added back below as a 0 dot). Mapping is by column NAME, so it
		* stays aligned with the jwdid panel's true event times over [-3,+4].
		matrix b = e(b_iw)
		matrix V = e(V_iw)
		local names : colnames b
		local K = colsof(b)
		clear
		set obs `K'
		gen event_year = .
		gen b1 = .
		gen V1 = .
		local j = 0
		foreach nm of local names {
			local ++j
			if regexm("`nm'", "treat_yearid([0-9]+)") {
				replace event_year = real(regexs(1)) - 8 in `j'
				replace b1 = b[1,`j']   in `j'
				replace V1 = V[`j',`j'] in `j'
			}
		}
		drop if missing(event_year)
		* add the omitted reference period t=-1 as an exact dot (no CI)
		set obs `=_N + 1'
		replace event_year = -1 in L
		replace b1 = 0 in L
		replace V1 = 0 in L
		sort event_year
		gen estid = 2

		* ---------- combine the two estimators ----------
		append using `sd_`m''
		sort estid event_year
		gen ci_upper = b1 + `z'*sqrt(V1)
		gen ci_lower = b1 - `z'*sqrt(V1)

		* horizontal dodge so the CIs do not overlap
		gen event_plot = event_year
		replace event_plot = event_year - 0.10 if estid==1
		replace event_plot = event_year + 0.10 if estid==2

		tempfile evi_`m'
		save `evi_`m''
	}

	* ---------------- Panel a: Median EVI ----------------
	use `evi_median', clear
	twoway ///
		(rcap ci_upper ci_lower event_plot if estid==1, lcolor(navy)) ///
		(scatter b1 event_plot if estid==1, mcolor(navy) msymbol(O)) ///
		(rcap ci_upper ci_lower event_plot if estid==2, lcolor(cranberry)) ///
		(scatter b1 event_plot if estid==2, mcolor(cranberry) msymbol(S) mfcolor(white)) ///
		, xlabel(-3(1)4, nogrid labsize(medsmall)) ylabel(, grid angle(0) labsize(medsmall)) ///
		xtitle("Years relative to treatment") ytitle("Effect on median EVI (x 100)") ///
		yline(0, lcolor(gs10) lpattern(dash)) xline(-0, lcolor(gs12)) ///
		legend(order(2 "Stacked DiD" 4 "Sun & Abraham (2021)") rows(1) position(6) size(small) region(lwidth(none))) ///
		graphregion(color(white)) plotregion(color(white)) 
	graph export "$out_fig/evi_median_combined.png", replace

	* ---------------- Panel b: Mean EVI ----------------
	use `evi_mean', clear
	twoway ///
		(rcap ci_upper ci_lower event_plot if estid==1, lcolor(navy)) ///
		(scatter b1 event_plot if estid==1, mcolor(navy) msymbol(O)) ///
		(rcap ci_upper ci_lower event_plot if estid==2, lcolor(cranberry)) ///
		(scatter b1 event_plot if estid==2, mcolor(cranberry) msymbol(S) mfcolor(white)) ///
		, xlabel(-3(1)4, nogrid labsize(medsmall)) ylabel(, grid angle(0) labsize(medsmall)) ///
		xtitle("Years relative to treatment") ytitle("Effect on mean EVI (x 100)") ///
		yline(0, lcolor(gs10) lpattern(dash)) xline(-0, lcolor(gs12)) ///
		legend(order(2 "Stacked DiD" 4 "Sun & Abraham (2021)") rows(1) position(6) size(small) region(lwidth(none))) ///
		graphregion(color(white)) plotregion(color(white)) 
	graph export "$out_fig/evi_mean_combined.png", replace
		

end

* =========================================
* EVI robustness (1): matched control locations
*   random control locations matched to treated farms on EVI two years before
*   treatment, within country x cohort (psmatch2, logit, 1 nearest neighbour);
*   same estimators and clustering as evi_analysis
*   -> $out_fig/evi_median_matchedcontrols.png, evi_mean_matchedcontrols.png (Supp. Fig. 2)
* =========================================
capture program drop evi_rob_matched
program define evi_rob_matched

	* ===================================================================
	* MATCHED CONTROLS - controls matched to treated farms on baseline EVI
	*     within COUNTRY x COHORT (psmatch2, logit, 1 nearest neighbor).
	*
	*     Baseline = each cohort's OWN t-2 (calendar year g-2): the panel keeps
	*     event times [-3,+5] only, so a common calendar baseline year would
	*     leave the later cohorts without a baseline observation.
	*
	*     For each country c and treated cohort g: treated farms of cohort g
	*     are matched to never-treated farms of the same country on EVI in
	*     calendar year g-2 (the same calendar year for both groups, so the
	*     comparison is not contaminated by year effects). Treated farms with
	*     NO EVI reading at t-2 cannot be matched and are dropped from the
	*     sample (the final merge keeps matched farms only). Cells whose
	*     baseline year falls outside the panel are skipped automatically:
	*     with a t-2 baseline the panel starts at 2017, so the 2018 cohort
	*     (baseline 2016) drops out and cohorts 2019+ are retained. Matched
	*     controls enter unweighted (a control used by several treated farms
	*     still appears once).
	* ===================================================================
	local z = 1.645
	use "final/EVI_farms_pseudo", clear
	gen event_window=year-calyear_treat if calyear_treat!=0
	keep if (event_window>=-3 & event_window<=4) | calyear_treat==0
	drop if calyear_treat>0 & calyear_treat<2018

	//farm-year extract used to build the match (kept small on purpose)
	tempfile fy
	preserve
		keep forestgardenid iso_trees calyear_treat year w1_evi_median
		drop if missing(w1_evi_median)
		save `fy'
	restore

	set seed 20260811
	tempfile matchedfarms
	local haveany = 0
	quietly levelsof iso_trees, local(ctries)
	quietly levelsof calyear_treat if calyear_treat>0, local(cohorts)

	foreach c of local ctries {
		foreach g of local cohorts {
			local base = `g' - 2
			preserve
				quietly use `fy', clear
				quietly keep if iso_trees=="`c'" & year==`base'
				//treated farms of THIS cohort vs never-treated farms of the same country
				quietly keep if calyear_treat==`g' | calyear_treat==0
				quietly bysort forestgardenid: keep if _n==1
				gen byte treatgrp=(calyear_treat==`g')
				quietly count if treatgrp==1
				local ntr=r(N)
				quietly count if treatgrp==0
				local nct=r(N)
				if `ntr'>0 & `nct'>0 {
					gen double u=runiform()
					sort u                       // random order for tie-breaking
					//the cell is already restricted to country `c' by the keep
					//above; the if-condition is kept explicit so the
					//within-country restriction is visible at the call site
					capture psmatch2 treatgrp w1_evi_median if iso_trees=="`c'", neighbor(1) logit
					//store the return code immediately: _rc is reset by later
					//commands, and `_rc' as a macro is empty (it is a scalar)
					local prc = _rc
					if `prc'==0 {
						//_weight is non-missing for treated units AND for the
						//control units actually used as matches
						quietly keep if _weight!=.
						keep forestgardenid
						if `haveany'==1 {
							quietly append using `matchedfarms'
						}
						quietly duplicates drop forestgardenid, force
						quietly save `matchedfarms', replace
						local haveany = 1
					}
					else {
						display as text "evi_rob_matched: psmatch2 skipped `c' cohort `g' (rc=`prc')"
					}
				}
			restore
		}
	}

	if `haveany'==0 {
		display as error "evi_rob_matched: no country x cohort cell could be matched"
		exit 498
	}

	merge m:1 forestgardenid using `matchedfarms', keep(match) nogen

	//sample-size report for the paper (log only, no table written)
	display as text _n "{hline 60}"
	display as text "evi_rob_matched: MATCHED SAMPLE"
	display as text "{hline 60}"
	quietly count
	display as text "  farm-year observations (matched panel) = " r(N)
	quietly distinct forestgardenid
	display as text "  farms (total)                          = " r(ndistinct)
	quietly distinct forestgardenid if calyear_treat>0
	display as text "  farms (treated)                        = " r(ndistinct)
	quietly distinct forestgardenid if calyear_treat==0
	display as text "  farms (matched controls)               = " r(ndistinct)
	quietly count if calyear_treat>0
	display as text "  observations (treated farms)           = " r(N)
	quietly count if calyear_treat==0
	display as text "  observations (control farms)           = " r(N)
	quietly levelsof calyear_treat if calyear_treat>0, local(kept)
	display as text "  cohorts retained                       = `kept'"

	* -------------------------------------------------------------------
	* Estimators on the SAME matched sample, one regression each, run for
	* BOTH outcomes (w1_evi_median and w1_evi_mean) -> one figure per outcome.
	* Specifications follow evi_analysis:
	*   (a) reghdfe (TWFE):  $treat_years gcpc_2017usd_w1,
	*                        absorb(fe_farm) vce(cluster fe_proj_area)
	*                        - HELPER ONLY, not plotted: its e(sample)
	*                          defines the sample jwdid and SA run on
	*   (b) jwdid (ETWFE):   gcpc_2017usd_w5, ivar(fe_farm) tvar(year)
	*                        gvar(calyear_treat) fevar(fe_ctryyr)
	*                        method(reghdfe) cluster(fe_proj_area) never
	*   (c) Sun & Abraham:   $treat_years ch_gcpc_2017usd_w1,
	*                        cohort(calyear_treat) control_cohort(never_treat)
	*                        absorb(fe_farm fe_cellyr) vce(cluster fe_proj_area)
	* Coefficients are mapped by column NAME (treat_yearidN -> event N-8),
	* with the omitted t=-1 (id7) added back as a 0 dot.
	* NOTE: the matched sample itself is built ONCE, on baseline MEDIAN EVI at
	* t-2, so both figures use the same farms and are directly comparable.
	* -------------------------------------------------------------------
	cap drop fe_cellyr
	egen fe_cellyr=group(fe_cell year)
	tempfile matchedsample
	save `matchedsample'

	foreach m in median mean {

		display as text _n "  --- outcome: w1_evi_`m' ---"

		* ---------- (a) reghdfe / TWFE: sample-definition helper only ----------
		* NOT plotted. Kept (as in evi_analysis) so its e(sample) defines the
		* estimation sample that jwdid and Sun & Abraham both run on.
		use `matchedsample', clear
		quietly reghdfe w1_evi_`m' $treat_years gcpc_2017usd_w1, absorb(fe_farm ) vce(cluster fe_proj_area)
		display as text "  N (reghdfe helper, not plotted)        = " e(N)
		cap drop insample_matched
		gen byte insample_matched = e(sample)
		tempfile matchedsample2
		save `matchedsample2'

		* ---------- (b) jwdid / ETWFE ----------
		use `matchedsample2', clear
		quietly jwdid w1_evi_`m' gcpc_2017usd_w5 if insample_matched==1, ivar(fe_farm) tvar(year) gvar(calyear_treat) fevar(fe_ctryyr) method(reghdfe) cluster(fe_proj_area) never
		display as text "  N (jwdid / ETWFE)                      = " e(N)
		quietly estat event
		matrix __ETAB = r(table)
		jwdid_evmat
		clear
		svmat EVMAT, names(col)
		sort event_year
		gen estid = 1
		tempfile est_jwdid
		save `est_jwdid'

		* ---------- (c) Sun & Abraham (2021) ----------
		use `matchedsample2', clear
		xtset fe_farm year
			quietly eventstudyinteract w1_evi_`m' $treat_years ch_gcpc_2017usd_w1 if insample_matched==1 , cohort(calyear_treat) control_cohort(never_treat) absorb(fe_farm fe_cellyr) vce(cluster fe_proj_area)
		
		display as text "  N (Sun & Abraham)                      = " e(N)
		matrix b = e(b_iw)
		matrix V = e(V_iw)
		local names : colnames b
		local K = colsof(b)
		clear
		set obs `K'
		gen event_year = .
		gen b1 = .
		gen V1 = .
		local j = 0
		foreach nm of local names {
			local ++j
			if regexm("`nm'", "treat_yearid([0-9]+)") {
				replace event_year = real(regexs(1)) - 8 in `j'
				replace b1 = b[1,`j']   in `j'
				replace V1 = V[`j',`j'] in `j'
			}
		}
		drop if missing(event_year)
		set obs `=_N + 1'
		replace event_year = -1 in L
		replace b1 = 0 in L
		replace V1 = 0 in L
		sort event_year
		gen estid = 3

		* ---------- combine the two plotted estimators and plot ----------
		append using `est_jwdid'
		sort estid event_year
		gen ci_upper = b1 + `z'*sqrt(V1)
		gen ci_lower = b1 - `z'*sqrt(V1)

		* horizontal dodge so the CIs do not overlap
		gen event_plot = event_year
		replace event_plot = event_year - 0.10 if estid==1
		replace event_plot = event_year + 0.10 if estid==3

		local mlab = proper("`m'")
		twoway ///
			(rcap ci_upper ci_lower event_plot if estid==1, lcolor(navy)) ///
			(scatter b1 event_plot if estid==1, mcolor(navy) msymbol(O)) ///
			(rcap ci_upper ci_lower event_plot if estid==3, lcolor(cranberry)) ///
			(scatter b1 event_plot if estid==3, mcolor(cranberry) msymbol(S) mfcolor(white)) ///
			, xlabel(-3(1)4, nogrid labsize(medsmall)) ylabel(, grid angle(0) labsize(medsmall)) ///
			xtitle("Years relative to treatment") ytitle("Effect on `mlab' EVI (x 100)") ///
			yline(0, lcolor(gs10) lpattern(dash)) xline(-0, lcolor(gs12)) ///
			legend(order(2 "Stacked DiD" 4 "Sun & Abraham (2021)") rows(1) position(6) size(small) region(lwidth(none))) ///
			graphregion(color(white)) plotregion(color(white))
		graph export "$out_fig/evi_`m'_matchedcontrols.png", replace
	}
	display as text "{hline 60}"

end

* =========================================
* EVI robustness (2): placebo-timing (randomization inference) test
*   pseudo treatment cohorts assigned to never-treated locations; distribution
*   of the placebo average of the t+1..t+4 ATTs vs the actual estimate
*   -> $out_fig/evi_median_rand_inference.png (Supp. Fig. 3)
* =========================================
capture program drop evi_rob_randinf
program define evi_rob_randinf

	* ===================================================================
	* RANDOMIZATION INFERENCE - placebo cohorts among never-treated farms
	*     Actual estimate: the jwdid (Wooldridge ETWFE) spec of evi_analysis;
	*     the test statistic is the AVERAGE of the event-study ATTs from t+1
	*     to t+4 (harvested via estat event / jwdid_evmat; averaged over the
	*     available horizons, so draws stay comparable if one is missing).
	*     Placebos: never-treated farms only; each rep assigns pseudo cohorts
	*     (drawn from the empirical treated cohort distribution, same treated
	*     share), trims pseudo-treated farms to the same -3..+4 event window,
	*     re-runs the SAME jwdid with gvar(pseudo cohort), and stores the
	*     average of the t+1..t+4 ATTs. The figure is the placebo distribution
	*     with the actual marked; the empirical two-sided p-value is in the note.
	*     R = 5 draws in this package so that the demo runs quickly; the paper
	*     uses R = 1,000 (set local R = 1000 below to reproduce Supp. Fig. 3).
	* ===================================================================
	local R = 5
	set seed 20260717

	use "final/EVI_farms_pseudo", clear
	gen event_window=year-calyear_treat if calyear_treat!=0
	keep if (event_window>=-3 & event_window<=4) | calyear_treat==0

	//estimation-sample tag, exactly as in evi_analysis
	quietly reghdfe w1_evi_median $treat_years gcpc_2017usd_w1, absorb(fe_farm fe_ctryyr) vce(cluster fe_proj_area)
	cap drop insample_evi_med
	gen byte insample_evi_med = e(sample)

	//actual: jwdid as in evi_analysis; statistic = average of the ATTs t+1..t+4
	quietly jwdid w1_evi_median gcpc_2017usd_w5 if insample_evi_med==1, ivar(fe_farm) tvar(year) gvar(calyear_treat) fevar(fe_ctryyr) method(reghdfe) cluster(fe_proj_area) never
	quietly estat event
	matrix __ETAB = r(table)
	jwdid_evmat
	//average of the ATTs over the AVAILABLE horizons in t+1..t+4 (an average,
	//unlike a sum, stays comparable if a horizon is missing)
	local b_actual = 0
	local n_ev = 0
	forvalues i=1/`=rowsof(EVMAT)' {
		if inrange(EVMAT[`i',1],1,4) {
			local b_actual = `b_actual' + EVMAT[`i',2]
			local ++n_ev
		}
	}
	local b_actual = `b_actual'/`n_ev'
	display as text "actual average of event-study ATTs (t+1..t+4, `n_ev' horizons) = " as result %8.4f `b_actual'

	//empirical cohort distribution (farm level) and overall treated share,
	//within the actual estimation sample
	quietly {
		preserve
			keep if insample_evi_med==1
			bysort forestgardenid: keep if _n==1
			count
			local n_farms=r(N)
			count if calyear_treat>0
			local n_treat=r(N)
			local p_treat=`n_treat'/`n_farms'
			keep if calyear_treat>0
			contract calyear_treat, freq(nf)
			gen double share=nf/`n_treat'
			gen double cumshare=sum(share)
			local ncoh=_N
			forvalues i=1/`ncoh' {
				local coh`i'=calyear_treat[`i']
				local cum`i'=cumshare[`i']
			}
		restore
	}

	//placebo loop on never-treated farms only (same sample conditions)
	keep if calyear_treat==0 & insample_evi_med==1
	tempname ri
	tempfile rifile
	postfile `ri' double(b_placebo) using `rifile', replace
	forvalues r=1/`R' {
		quietly {
			preserve
				//one uniform draw per farm, spread to all its rows
				bysort forestgardenid (year): gen double u1=runiform() if _n==1
				bysort forestgardenid (year): replace u1=u1[1]
				bysort forestgardenid (year): gen double u2=runiform() if _n==1
				bysort forestgardenid (year): replace u2=u2[1]
				//pseudo cohort: treated with prob p_treat, year from empirical
				//shares; 0 = pseudo never-treated (matches gvar==0 + `never')
				gen pseudo_year=0
				forvalues i=`ncoh'(-1)1 {
					replace pseudo_year=`coh`i'' if u1<`p_treat' & u2<=`cum`i''
				}
				//mirror the actual sample: pseudo-treated farms trimmed to the
				//-3..+4 pseudo event window
				gen pseudo_event=year-pseudo_year if pseudo_year>0
				drop if pseudo_year>0 & !inrange(pseudo_event,-3,4)
				//same jwdid, with the pseudo cohort as gvar; a failed draw is
				//skipped (counted in the note via the posted draw count)
				capture jwdid w1_evi_median gcpc_2017usd_w5, ivar(fe_farm) tvar(year) gvar(pseudo_year) fevar(fe_ctryyr) method(reghdfe) cluster(fe_proj_area) never
				if _rc==0 {
					capture estat event
					if _rc==0 {
						matrix __ETAB = r(table)
						jwdid_evmat
						//a draw whose event table holds no event-time column
						//leaves no EVMAT: skip it like a failed draw
						capture confirm matrix EVMAT
						if _rc==0 {
							local bsum = 0
							local n_ev = 0
							forvalues i=1/`=rowsof(EVMAT)' {
								if inrange(EVMAT[`i',1],1,4) {
									local bsum = `bsum' + EVMAT[`i',2]
									local ++n_ev
								}
							}
							//average over available horizons (comparable across
							//draws even when a horizon is missing)
							if `n_ev'>0 post `ri' (`=`bsum'/`n_ev'')
						}
						else noisily display as text "  randomization inference: draw `r' skipped (no event-time columns)"
					}
				}
			restore
		}
		if mod(`r',25)==0 display as text "  randomization inference: `r'/`R' done"
	}
	postclose `ri'

	use `rifile', clear
	quietly count
	local Rdone=r(N)
	quietly count if abs(b_placebo)>=abs(`b_actual')
	local pval=r(N)/`Rdone'
	quietly summarize b_placebo
	//widen the x-axis so it always spans BOTH the placebo distribution and the
	//actual estimate; otherwise histogram auto-scales to the placebos only and
	//the actual (which lies far in the tail) falls off-plot and the line vanishes
	local lo = r(min)
	local hi = r(max)
	if `b_actual' < `lo' local lo = `b_actual'
	if `b_actual' > `hi' local hi = `b_actual'
	local pad = 0.10*(`hi'-`lo')
	if `pad'==0 local pad = 0.1
	local lo = `lo' - `pad'
	local hi = `hi' + `pad'
	//empirical two-sided p-value for the figure note (0 of R draws -> p < 1/R,
	//i.e. p < 0.001 with the paper's 1,000 draws)
	if `pval' == 0 local ptxt "p < `: display %5.3f 1/`Rdone''"
	else local ptxt "p = `: display %5.3f `pval''"
	display as text "randomization inference: `Rdone' of `R' draws completed; `ptxt'"
	histogram b_placebo, bin(30) fcolor(navy%40) lcolor(navy) ///
		xline(`b_actual', lcolor(cranberry) lwidth(medthick)) ///
		xscale(range(`lo' `hi')) xlabel(#8, nogrid) ///
		xtitle("Placebo average of event-study ATTs, t+1 to t+4 (median EVI x 100)") ///
		ytitle("Density") ///
		text(0 `b_actual' "Actual", place(e) color(cranberry) size(small)) ///
		note("Actual average of ATTs t+1..t+4 (red line) = `: display %6.3f `b_actual''; randomization-inference `ptxt'") ///
		graphregion(color(white)) plotregion(color(white))
	graph export "$out_fig/evi_median_rand_inference.png", replace

end

* wrapper: both EVI robustness tests in sequence
capture program drop evi_robustness
program define evi_robustness
	evi_rob_matched
	evi_rob_randinf
end

* =========================================
* Validation of the intervention in the survey data: numbers of trees,
* marketable crops and food crops (levels and per hectare), PPML
*   -> $out_tab/trees_crops.tex                  (Supp. Table 5)
*   -> $out_fig/trees_crops_post_coefficients.png (Fig. 3)
* =========================================
capture program drop trees_crops
program define trees_crops
	eststo clear
	use "final/TREES_data_pseudo", clear
	replace w5_gcpc_2017usd=w5_gcpc_2017usd/1000
	//N trees
	eststo m1: quietly ppmlhdfe w_trees post $gcpc if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y" 	
			 estadd local month "Y"
			 estadd local year_cohort "Y"
	*----------------
	//N marketable 
	eststo m2: quietly ppmlhdfe w_marketableproducts post $gcpc if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y" 
			 estadd local month "Y"
			 estadd local year_cohort "Y" 
			 
	*-----------------
	//N Food crops
	eststo m3: quietly ppmlhdfe w_foodcrops post $gcpc if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y" 	
			 estadd local month "Y"
			 estadd local year_cohort "Y"

	*-----------------
	//N trees per hectare 
	eststo m4: quietly ppmlhdfe w1_trees_per_ha post $gcpc if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y" 	
			 estadd local month "Y"
			 estadd local year_cohort "Y"
			 
	*-----------------
	//N mkt crops per hectare 
	eststo m5: quietly ppmlhdfe w1_mark_prod_ha post $gcpc if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y" 	
			 estadd local month "Y"
			 estadd local year_cohort "Y"
			 
	*-----------------
	//N food crops per hectare 
	eststo m6: quietly ppmlhdfe w1_foodcrop_ha post $gcpc if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y" 	
			 estadd local month "Y"
			 estadd local year_cohort "Y"

	esttab m* using  "$out_tab/trees_crops.tex" , nonotes nobaselevels noconstant noomitted width(\hsize) compress booktabs label interaction(" $\times$ ") replace  nogaps b(3) r2(2) nomtitles ///
		hlinechar(`=char(151)')   ///
		    prehead("\begin{tabular}{l*{6}{c}}" ///
            "\toprule" ///
            "Dep. var.: Number of & Trees & Mkt. crops & Food crops & Trees/ha & Mkt. crops/ha & Food crops/ha \\" "\cmidrule(lr){2-2}\cmidrule(lr){3-3}\cmidrule(lr){4-4}\cmidrule(lr){5-5}\cmidrule(lr){6-6}\cmidrule(lr){7-7}")	///
		 keep(post $gcpc) ///
		 coeflabels(w5_gcpc_2017usd "GPC per capita (\textdollar 1,000)" post "1(Post 1st survey)") ///
		star(* 0.1 ** 0.05 *** .01)    ///
		stats(contr farmfe month N , fmt(0 0 0 %12.0gc ) labels("Controls" "Farm F.E." "Month F.E." "Observations" ) ) 
	
	*---------------
	//plots
	//Post-indicator effects across tree/crop outcomes (single coefficient plot)
	eststo clear
	use "final/TREES_data_pseudo", clear
	tempname _res
	tempfile _res_file
	postfile `_res' str24 outcome double(b se lo hi) using `_res_file'
	foreach y in w_trees w_marketableproducts w_foodcrops {
		ppmlhdfe `y' post $gcpc if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
		local b  = _b[post]
		local se = _se[post]
		local lo = `b' - 1.645*`se'
		local hi = `b' + 1.645*`se'
		post `_res' ("`y'") (`b') (`se') (`lo') (`hi')
	}
	postclose `_res'
	preserve
	use `_res_file', clear
	gen outcome_n = .
	replace outcome_n = 1 if outcome=="w_trees"
	replace outcome_n = 2 if outcome=="w_marketableproducts"
	replace outcome_n = 3 if outcome=="w_foodcrops"
	sort outcome_n

	//coefficient label shown above each CI: signed coefficient x100 in log points
	//(all three outcomes are PPML). The explicit "+" prevents a quick reader
	//from mistaking the coefficients for declines.
	gen blabel = cond(b>=0, "+", "") + string(100*b, "%9.1f") + " log p."

	twoway ///
		(rcap hi lo outcome_n, lcolor(navy)) ///
		(scatter b outcome_n, mcolor(navy) msymbol(O)) ///
		(scatter hi outcome_n, msymbol(none) mlabel(blabel) mlabposition(12) mlabcolor(navy) mlabsize(small)), ///
		xscale(range(0.5 3.5)) ///
		xlabel(1 "N. trees" 2 "N. mkt. crops" 3 "N. food crops", labsize(small) nogrid) ///
		xtitle("") ytitle("Change relative to first survey") ///
		yline(0, lpattern(dash) lcolor(gs8)) legend(off) ///
		graphregion(color(white)) plotregion(color(white))
	graph export "$out_fig/trees_crops_post_coefficients.png", replace
	restore
end

* =========================================
* Natural capital intensity and livestock: first-difference regressions of the
* change in livestock per hectare on the change in trees / crops per hectare
*   -> $out_tab/nlivestock_trees_crops.tex (Supp. Table 13)
* =========================================
capture program drop channel
program define channel
	eststo clear
	use "final/TREES_data_pseudo", clear

	//survey months of the first and last survey (constant within farm; computed
	//before duplicates drop, which keeps an arbitrary row's month_survey)
	bysort fe_farm (time_from_start): gen beg_month  = month_survey[1]
	bysort fe_farm (time_from_start): gen last_month = month_survey[_N]
	label var beg_month "Month of first survey"
	label var last_month "Month of last survey"

	//one observation per farm
	keep if last_obs==1 & insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1

	eststo m1: quietly 	reghdfe wch_tot_animal_ha wch_trees1000_per_ha w5ch_gcpc_2017usd , absorb(beg_month last_month) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local month "Y" 
			 estadd local yearfe "Y" 

	eststo m2: quietly reghdfe wch_tot_animal_ha wch_mark_prod_ha w5ch_gcpc_2017usd , absorb(beg_month last_month year_cohort) vce(cluster fe_proj_area)	
			 estadd local contr "Y"
			 estadd local month "Y" 
			 estadd local yearfe "Y"  	

	eststo m3: quietly reghdfe wch_tot_animal_ha wch_foodcrop_ha  w5ch_gcpc_2017usd , absorb(beg_month last_month year_cohort) vce(cluster fe_proj_area)	
			 estadd local contr "Y"
			 estadd local month "Y" 
			 estadd local yearfe "Y" 

	eststo m4: quietly reghdfe wch_tot_animal_ha wch_trees1000_per_ha wch_foodcrop_ha  w5ch_gcpc_2017usd , absorb(beg_month last_month year_cohort) vce(cluster fe_proj_area)	
			 estadd local contr "Y"
			 estadd local month "Y" 
			 estadd local yearfe "Y" 

	eststo m5: quietly reghdfe wch_tot_animal_ha wch_trees1000_per_ha wch_mark_prod_ha w5ch_gcpc_2017usd , absorb(beg_month last_month year_cohort) vce(cluster fe_proj_area)	
			 estadd local contr "Y"
			 estadd local month "Y" 
			 estadd local yearfe "Y" 

	esttab m* using  "$out_tab/nlivestock_trees_crops.tex" , nonotes nobaselevels noomitted noconstant width(\hsize) compress booktabs label interaction(" $\times$ ") replace  nogaps b(3) r2(2) nomtitles nonumbers ///
		hlinechar(`=char(151)')   ///
		mlabels("(1)" "(2)" "(3)" "(4)" "(5)" , lhs("Dep. Var.: $\Delta$ N. Livestock/ha")) ///
		drop(_cons w5ch_gcpc_2017usd) ///
		coeflabels(wch_trees1000_per_ha "$\Delta$ Trees/ha (1,000s/ha)" wch_mark_prod_ha "$\Delta$ Market Crops/ha" wch_foodcrop_ha "$\Delta$ Food Crops/ha" w5ch_gcpc_2017usd "$\Delta$ GCP per capita (\textdollar 1,000)") ///
		star(* 0.1 ** 0.05 *** .01)    ///
		stats(month N r2, fmt(0 %12.0gc 3) labels( "First/Last Month F.E."  "Observations" "R\$^2\$") )
end

* =========================================
* Food security and dietary diversity outcomes (OLS)
*   -> $out_tab/health.tex                  (Supp. Table 6)
*   -> $out_fig/health_post_coefficients.png (Fig. 4)
* =========================================
capture program drop health
program define health
	eststo clear
	use "final/TREES_data_pseudo", clear
	replace w5_gcpc_2017usd=w5_gcpc_2017usd/1000


	foreach var of varlist slepthungry noeating24hrs guestsstaylastnight childreninschool {
		gen `var'id=1 if `var'=="Yes"
		replace `var'id=0 if `var'=="No"
	}

	*-------------
	//No food
	//Food-security tables: effect of the post indicator
	gen post2019=(year_survey>2019)
	eststo m1: quietly reghdfe nofoodid post $gcpc post2019 if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y"
			 estadd local monthfe "Y"

	eststo m2: quietly		reghdfe hfiascatid4 post $gcpc post2019 if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y"
			 estadd local monthfe "Y"

	
	eststo m3: quietly reghdfe hdds_id3 post $gcpc if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y"
			 estadd local monthfe "Y"

	esttab m* using  "$out_tab/health.tex" , nonotes nobaselevels noconstant noomitted width(\hsize) compress booktabs label interaction(" $\times$ ") replace  nogaps nomtitles b(3) r2(2) ///
hlinechar(`=char(151)')   ///
		    prehead("\begin{tabular}{l*{3}{c}}" ///
            "\toprule" ///
            "Dep. var.: & No food & Food insecurity & Low food diversity \\" "\cmidrule(lr){2-2}\cmidrule(lr){3-3}\cmidrule(lr){4-4}")	///
		 keep(post $gcpc) ///
		 coeflabels(post "1(Post 1st survey)" $gcpc "GCP per capita (\textdollar 1,000)") ///
		star(* 0.1 ** 0.05 *** .01)    ///
		stats(contr farmfe monthfe N r2, fmt(0 0 0 %12.0gc 3) labels("Controls" "Farm F.E." "Month F.E." "Observations" "R\$^2\$") )
	
	//Post-indicator effects across food-security outcomes (single coefficient plot)
	tempname _res
	tempfile _res_file
	postfile `_res' str24 outcome double(b se lo hi) using `_res_file'
	foreach y in nofoodid hfiascatid4 hdds_id3 {
		if "`y'"=="hdds_id3" {
			reghdfe `y' post $gcpc if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)

		} 
		else {
			reghdfe `y' post $gcpc post2019 if insample & n_obs>1 & gap_ok==1 & gap_ok_fl==1, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
		}
			
		local b  = _b[post]
		local se = _se[post]
		local lo = `b' - 1.645*`se'
		local hi = `b' + 1.645*`se'
		post `_res' ("`y'") (`b') (`se') (`lo') (`hi')
	}
	postclose `_res'
	preserve
	use `_res_file', clear
	gen outcome_n = .
	replace outcome_n = 1 if outcome=="nofoodid"
	replace outcome_n = 2 if outcome=="hfiascatid4"
	replace outcome_n = 3 if outcome=="hdds_id3"
	sort outcome_n

	//coefficient label shown above each CI: coefficient x100 in percentage
	//points (all three outcomes are LPMs)
	gen blabel = string(100*b, "%9.1f") + " pp"

	twoway ///
		(rcap hi lo outcome_n, lcolor(navy)) ///
		(scatter b outcome_n, mcolor(navy) msymbol(O)) ///
		(scatter hi outcome_n, msymbol(none) mlabel(blabel) mlabposition(12) mlabcolor(navy) mlabsize(small)), ///
		xscale(range(0.5 3.5)) ///
		xlabel(1 "No food" 2 "Food insecurity" 3 "Low food diversity", labsize(small) nogrid) ///
		xtitle("") ytitle("Change relative to first survey") ///
		yscale(range(0)) yline(0, lpattern(dash) lcolor(gs8)) legend(off) ///
		graphregion(color(white)) plotregion(color(white))
	graph export "$out_fig/health_post_coefficients.png", replace
	restore
end

* =========================================
* Additional robustness of the livestock results
*   farm-size trimming  -> $out_tab/rob_size.tex                  (Supp. Table 9)
*   heterogeneity by baseline household characteristics
*                       -> $out_tab/cross_farm_char_livestock.tex (Supp. Table 14)
* =========================================
capture program drop add_rob
program define add_rob

	eststo clear
	use "final/TREES_data_pseudo", clear
	replace w5_gcpc_2017usd=w5_gcpc_2017usd/1000
	replace l1_gcpc_2017usd=l1_gcpc_2017usd/1000000

	//farm size of 0 is non-response coded as zero: treat as missing in the
	//wb_w_size control of the characteristics table. Raw `size` is not recoded
	//here because the percentile cuts below handle zeros explicitly.
	replace wb_w_size=. if wb_w_size==0

	//estimation sample (as in baseline); the median splits below are computed on it
	keep if insample==1 & n_obs>1 & gap_ok==1 & gap_ok_fl==1
	egen p95_size = pctile(size), p(95)
	//p5 computed excluding farms with zero size (zeros are non-response and
	//would disturb the percentile calculation)
	summarize size if size!=0, detail
	gen p5_size = r(p5)
	
	
	//removing very large and very small farms farms
	eststo m1: quietly 	reghdfe havelivestock post $gcpc if baselinesize<p95_size & baselinesize!=0, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y" 
			 estadd local monthfe "Y"
			 estadd local yearfe "Y"
			 estadd local rob_sample "Size $< 95^{th}$" 
			 estadd scalar myr2 e(r2)
			 
	eststo m2: quietly 	reghdfe havelivestock post $gcpc if baselinesize>p5_size & !missing(baselinesize), absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y" 
			 estadd local monthfe "Y"
			 estadd local yearfe "Y"
			 estadd local rob_sample "Size $> 5^{th}$" 
			 estadd scalar myr2 e(r2)
	
	eststo m3: quietly reghdfe havelivestock post $gcpc if baselinesize>p5_size & baselinesize<p95_size, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y" 
			 estadd local monthfe "Y"
			 estadd local yearfe "Y"
			 estadd local rob_sample "5$^{th} <$ Size $< 95^{th}$" 
			 estadd scalar myr2 e(r2)
	
	eststo m4: quietly ppmlhdfe tot_animal post $gcpc if baselinesize<p95_size & baselinesize!=0, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y" 
			 estadd local monthfe "Y"
			 estadd local yearfe "Y"
			 estadd local rob_sample "Size $< 95^{th}$" 
			 estadd local myr2 "-"
			 
	eststo m5: quietly ppmlhdfe tot_animal post $gcpc if baselinesize>p5_size & !missing(baselinesize), absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y" 
			 estadd local monthfe "Y"
			 estadd local yearfe "Y"
			 estadd local rob_sample "Size $> 5^{th}$" 
			 estadd local myr2 "-"

	eststo m6: quietly ppmlhdfe tot_animal post $gcpc if baselinesize>p5_size & baselinesize<p95_size, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			 estadd local contr "Y"
			 estadd local farmfe "Y" 
			 estadd local monthfe "Y"
			 estadd local yearfe "Y"
			 estadd local rob_sample "5$^{th} <$ Size $< 95^{th}$" 
			 estadd local myr2 "-"
		esttab m* using  "$out_tab/rob_size.tex" , nonotes nobaselevels noconstant noomitted width(\hsize) compress booktabs label interaction(" $\times$ ") replace  nogaps b(3) r2(2) nomtitles ///
		hlinechar(`=char(151)')   ///
		prehead("\begin{tabular}{l*{6}{c}}" ///
            "\toprule" ///
            "Dep. var.: &  \multicolumn{3}{c}{Own livestock} & \multicolumn{3}{c}{N. livestock}  \\" ///
"\cmidrule(lr){2-4}\cmidrule(lr){5-7}") ///
		 keep(post $gcpc) ///
		 coeflabels($gcpc "GCP per capita (\textdollar 1,000)" post "1(Post 1st survey)") ///
		star(* 0.1 ** 0.05 *** .01)    ///
		stats(rob_sample contr farmfe monthfe N myr2, fmt(0 0 0 0 %12.0gc 3) labels("Sample" "Controls" "Farm F.E." "Month F.E." "Observations" "R\$^2\$") ) 
		
	*-------
	//Heterogeneity: does the treatment effect on livestock vary with baseline
	//household characteristics? Panel PPML with the post indicator interacted
	//with an indicator for the characteristic (above-median split for continuous
	//characteristics, the raw dummy for binary ones), farm + month FE, standard
	//errors clustered by project area. The main effect of the
	//characteristic is absorbed by the farm FE; the reported coefficients are
	//post and post x characteristic.
	eststo clear

	//baseline household/cell characteristics (measured at the baseline year)
	label variable wb_gender           "1(House. Head Male)"
	label variable wb_adult_mostedu    "1(Adult Most Edu.)"
	label variable wb_worriedfood      "1(Worried Food)"
	label variable wb_lackfood         "1(Lacked Food)"
	label variable post                "1(Post 1st survey)"

	//children relative to household size, at baseline
	gen b_n_child_tot=b_n_childrennumber/wb_w_npeople if !missing(b_n_childrennumber) & !missing(wb_w_npeople)

	//above-median indicators for the continuous characteristics (xtile nq(2)
	//over the loaded sample, as specified; the if-guard keeps true missings
	//missing rather than coding them as below-median — wb_w_size has missings
	//from the ==0 non-response recode at the top of this program)
	xtile q_age     = wb_w_participantage, nq(2)
	xtile q_npeople = wb_w_npeople,        nq(2)
	xtile q_child   = b_n_child_tot,       nq(2)
	xtile q_femtot  = wb_females_tot,      nq(2)
	xtile q_size    = wb_w_size,           nq(2)

	gen above_age     = (q_age==2)     if !missing(q_age)
	gen above_npeople = (q_npeople==2) if !missing(q_npeople)
	gen above_child   = (q_child==2)   if !missing(q_child)
	gen above_femtot  = (q_femtot==2)  if !missing(q_femtot)
	gen above_size    = (q_size==2)    if !missing(q_size)

	label var above_age     "1(Age $>$ Med.)"
	label var above_npeople "1(N. House. Members $>$ Med.)"
	label var above_child   "1(Child. to Tot. $>$ Med.)"
	label var above_femtot  "1(Females to Tot. $>$ Med.)"
	label var above_size    "1(Farm Size $>$ Med.)"

	//one characteristic at a time
	local hetvars above_age wb_gender wb_adult_mostedu above_child wb_worriedfood wb_lackfood above_size
	foreach var of local hetvars {

			eststo m_`var': quietly ppmlhdfe tot_animal c.post##c.`var' $gcpc, absorb(month_survey fe_farm) vce(cluster fe_proj_area)
			estadd local farmfe "Y"
			 estadd local monthfe "Y"

	}

	esttab m_* using  "$out_tab/cross_farm_char_livestock.tex" , nonotes nobaselevels noconstant noomitted width(\hsize) compress booktabs label interaction(" $\times$ ") replace nogaps nomtitles nonumbers  nogaps b(3) ///
		hlinechar(`=char(151)')   ///
		mlabels("(1)" "(2)" "(3)" "(4)" "(5)" "(6)" "(7)" , lhs("Dep. Var.: N. Livestock")) ///
		keep(post *#*) ///
		star(* 0.1 ** 0.05 *** .01)    ///
		stats(farmfe monthfe N, fmt(0 0 %12.0gc) labels("Farm F.E." "Month F.E." "Observations") )
end

* =========================================
* MAIN: run every program in order (called by code/run_demo.do)
* =========================================
capture program drop main
program define main
	baseline
	leave_one_out
	median_reg
	plots_wealth
	evi_analysis
	evi_robustness
	trees_crops
	channel
	health
	add_rob
end
