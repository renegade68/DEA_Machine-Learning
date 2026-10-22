*! _dea_ml_frontier_graph.ado v0.1.0
*! Internal postestimation graph used by dea_ml, frontiergraph
*! Author: Choonjoo Lee, Korea National Defense University
*!         Nonsan-si, Republic of Korea; bloom.rampike@gmail.com
*! Method references: This helper visualizes results produced by the estimation
*!   commands; see the corresponding standalone command ado/help files.

program define _dea_ml_frontier_graph
    version 16.1
    syntax [, options(string)]

    local xvars `e(x)'
    local yvars `e(y)'
    local nx : word count `xvars'
    local ny : word count `yvars'

    if `ny' != 1 | !inlist(`nx',1,2) {
        di as err "frontiergraph requires one output and one or two inputs"
        di as err "the graph is available only when the displayed frontier has two or three dimensions"
        exit 198
    }
    capture confirm matrix e(frontier)
    if _rc {
        di as err "frontiergraph is not available for method `e(cmd)'"
        exit 498
    }

    tempvar touse
    quietly generate byte `touse' = 1
    foreach v of varlist `xvars' `yvars' {
        quietly replace `touse' = 0 if missing(`v')
    }
    quietly count if `touse'
    if r(N) != rowsof(e(frontier)) {
        di as err "frontiergraph could not align e(frontier) with the estimation observations"
        exit 498
    }
    tempname F
    matrix `F' = e(frontier)
    preserve
    quietly keep if `touse'
    capture drop zzdeamlfit
    quietly generate double zzdeamlfit = .
    mata: st_store(., "zzdeamlfit", st_matrix("`F'")[.,1])

    local yvar : word 1 of `yvars'
    local ylab : variable label `yvar'
    if `"`ylab'"' == "" local ylab "`yvar'"

    if `nx' == 1 {
        local x1 : word 1 of `xvars'
        local xlab : variable label `x1'
        if `"`xlab'"' == "" local xlab "`x1'"
        local connect ""
        if inlist("`e(cmd)'","eatree","rfeat","beateat","svfeat") ///
            local connect "connect(stairstep)"
        twoway ///
            (scatter `yvar' `x1', msymbol(Oh) mcolor(gs10) msize(vsmall)) ///
            (line zzdeamlfit `x1', sort `connect' lcolor(black) lpattern(solid) lwidth(medthick)), ///
            xtitle(`"`xlab'"') ytitle(`"`ylab'"') ///
            legend(order(1 "Observed output" 2 "Fitted frontier") cols(1)) ///
            graphregion(color(white)) plotregion(color(white)) `options'
    }
    else {
        local x1 : word 1 of `xvars'
        local x2 : word 2 of `xvars'
        local xlab1 : variable label `x1'
        local xlab2 : variable label `x2'
        if `"`xlab1'"' == "" local xlab1 "`x1'"
        if `"`xlab2'"' == "" local xlab2 "`x2'"
        twoway contour zzdeamlfit `x2' `x1', ///
            interp(thinplatespline) levels(8) ///
            crule(linear) scolor(gs15) ecolor(gs4) ///
            xtitle(`"`xlab1'"') ytitle(`"`xlab2'"') ///
            ztitle("Fitted output") ///
            graphregion(color(white)) plotregion(color(white)) `options'
    }
    restore
end
