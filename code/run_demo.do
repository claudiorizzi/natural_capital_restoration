*===============================================================================
* run_demo.do - run the paper's Stata programs on the pseudo (synthetic) data
*
* Usage (from the repository root, the folder holding README.md, so that c(pwd)
* is the package root):
*   Windows:   "C:\Program Files\Stata19\StataMP-64.exe" /e do code\run_demo.do
*   Mac/Linux: stata-mp -b do code/run_demo.do
*   or open Stata, `cd` to the repository root and type:  do code/run_demo.do
*
* What it does
*   1. sets the three path globals the scripts rely on
*        $data_dir  folder holding the final .dta files   (data)
*        $out_tab   LaTeX tables written by esttab         (output/tables)
*        $out_fig   PNG figures written by graph export    (output/figures)
*   2. runs code/02_summary_stats.do, which only DEFINES its programs, then
*      calls them one by one
*   3. runs code/03_main_estimates.do and calls its programs in the order of
*      that script's `main` program
*   Every program is wrapped in `capture noisily`, so one failure does not stop
*   the demo; the exit code and run time of each program are listed at the end
*   (and in output/run_demo.log).
*
* Required user-written packages (install once with ssc install <pkg>):
*   reghdfe ftools hdfe ppmlhdfe jwdid eventstudyinteract avar erepost estout
*   distinct winsor2 qreg2 psmatch2
*===============================================================================
clear all
set more off
set rmsg off

* ---- locate the package root ------------------------------------------------
global root "`c(pwd)'"
capture confirm file "$root/data/final/TREES_data_pseudo.dta"
if _rc {
    * allow launching from inside code/
    global root "`c(pwd)'/.."
    capture confirm file "$root/data/final/TREES_data_pseudo.dta"
    if _rc {
        display as error "run_demo.do: cannot find data/final/TREES_data_pseudo.dta"
        display as error "Run this file from the repository root (the folder holding README.md), or create the pseudo data first:"
        display as error "    uv run python code/01_make_pseudo_data.py"
        exit 601
    }
}
global data_dir "$root/data"
global out_tab  "$root/output/tables"
global out_fig  "$root/output/figures"
capture mkdir "$root/output"
capture mkdir "$out_tab"
capture mkdir "$out_fig"

capture log close _all
log using "$root/output/run_demo.log", replace text

display as text _n "{hline 78}"
display as text "TREES replication demo - pseudo data"
display as text "  package root : $root"
display as text "  Stata        : `c(stata_version)' (`c(edition_real)'), `c(os)' `c(machine_type)'"
display as text "  started      : `c(current_date)' `c(current_time)'"
display as text "{hline 78}"

* ---- programs and their sources -----------------------------------------------
* 02_summary_stats.do
local progs_sum   sum_stats_farms sum_stats_evi wb_comparison single_multi_comparison
* 03_main_estimates.do (order of its `main` program)
local progs_main  baseline leave_one_out median_reg plots_wealth evi_analysis ///
                  evi_robustness trees_crops channel health add_rob

* results are accumulated in local macros (the analysis scripts start with
* `clear all`, which would close a postfile opened here)
local log_prog ""
local log_script ""
local log_rc ""
local log_secs ""

* ---- 1. summary statistics --------------------------------------------------------
display as text _n ">>> do code/02_summary_stats.do (defines the programs)"
do "$root/code/02_summary_stats.do"
foreach p of local progs_sum {
    display as text _n "{hline 78}" _n ">>> `p'" _n "{hline 78}"
    timer clear 1
    timer on 1
    capture noisily `p'
    local rc = _rc
    timer off 1
    quietly timer list 1
    local secs = r(t1)
    if `rc' display as error ">>> `p' FAILED with rc=`rc' after `=round(`secs',0.1)' s"
    else    display as result ">>> `p' done in `=round(`secs',0.1)' s"
    local log_prog   `log_prog' `p'
    local log_script `log_script' 02
    local log_rc     `log_rc' `rc'
    local log_secs   `log_secs' `=round(`secs',0.1)'
}

* ---- 2. main estimates ----------------------------------------------------------
display as text _n ">>> do code/03_main_estimates.do (defines the programs)"
do "$root/code/03_main_estimates.do"
foreach p of local progs_main {
    display as text _n "{hline 78}" _n ">>> `p'" _n "{hline 78}"
    timer clear 1
    timer on 1
    capture noisily `p'
    local rc = _rc
    timer off 1
    quietly timer list 1
    local secs = r(t1)
    if `rc' display as error ">>> `p' FAILED with rc=`rc' after `=round(`secs',0.1)' s"
    else    display as result ">>> `p' done in `=round(`secs',0.1)' s"
    local log_prog   `log_prog' `p'
    local log_script `log_script' 03
    local log_rc     `log_rc' `rc'
    local log_secs   `log_secs' `=round(`secs',0.1)'
}

* ---- summary -------------------------------------------------------------------------
clear
local n : word count `log_prog'
set obs `n'
gen str40 program = ""
gen str12 script  = ""
gen double rc = .
gen double seconds = .
forvalues i = 1/`n' {
    replace program = "`: word `i' of `log_prog''"   in `i'
    replace script  = "`: word `i' of `log_script''" in `i'
    replace rc      = `: word `i' of `log_rc''        in `i'
    replace seconds = `: word `i' of `log_secs''      in `i'
}
format seconds %9.1f
display as text _n "{hline 78}" _n "RUN SUMMARY (rc = Stata return code, 0 = ok)" _n "{hline 78}"
list program script rc seconds, noobs clean
quietly summarize seconds
display as text "total run time: " %6.1f r(sum) " seconds"
export delimited using "$root/output/run_demo_summary.csv", replace
quietly count if rc != 0
if r(N) display as error "`r(N)' program(s) failed - see output/run_demo.log"
else display as result "all programs ran without error"
display as text "finished: `c(current_date)' `c(current_time)'"
log close
