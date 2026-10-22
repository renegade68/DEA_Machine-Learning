*! csvf.ado v0.1.2
*! Convexified Support Vector Frontier, Valero-Carreras et al. (2021), model (10)
*! Author: Choonjoo Lee, Korea National Defense University
*!         Nonsan-si, Republic of Korea; bloom.rampike@gmail.com
*! Method reference: Valero-Carreras, Aparicio, and Guerrero (2021),
*!   "Support vector frontiers: A new approach for estimating production
*!   functions through support vector machines", Omega 104:102490.
*!   DOI: 10.1016/j.omega.2021.102490

program define csvf, eclass
    version 16.1
    tempname previous
    _estimates hold `previous', copy restore nullok
    capture noisily _csvf_fit `0'
    local rc = _rc
    if `rc' {
        _estimates unhold `previous'
        exit `rc'
    }
    _estimates drop `previous'
end

program define _csvf_fit, eclass
    version 16.1
    syntax [if] [in] , X(varlist numeric) Y(varname numeric) [GEN(name) D(integer 3) ///
        C(real 1) EPS(real 0) SEED(integer 12345) TUNE ///
        CGRID(numlist >0) EPSGRID(numlist >=0) DGRID(numlist integer >0) ///
        FOLDS(integer 5)]

    if "`gen'"!="" {
        confirm new variable `gen'
    }
    if `d'<1 | `c'<=0 | `eps'<0 | `folds'<2 {
        display as error "require d()>=1, c()>0, eps()>=0, and folds()>=2"
        exit 198
    }

    tempvar touse
    marksample touse, novarlist
    markout `touse' `x' `y'
    quietly count if `touse'
    if r(N)<2 {
        display as error "Too few complete observations."
        exit 2001
    }
    foreach v of varlist `x' {
        quietly count if `touse' & `v'<0
        if r(N) {
            display as error "input variables must be nonnegative"
            exit 459
        }
    }
    quietly count if `touse' & `y'<=0
    if r(N) {
        display as error "the output variable must be strictly positive for radial efficiency"
        exit 459
    }

    local opts d(`d') c(`c') eps(`eps') seed(`seed') folds(`folds')
    if "`tune'"!="" local opts `opts' tune
    if "`cgrid'"!="" local opts `opts' cgrid(`cgrid')
    if "`epsgrid'"!="" local opts `opts' epsgrid(`epsgrid')
    if "`dgrid'"!="" local opts `opts' dgrid(`dgrid')

    quietly svfeat if `touse', x(`x') y(`y') `opts'
    local N=e(N)
    local selected_C=e(c_pen)
    local selected_eps=e(eps)
    local selected_d=e(d)
    local tuned=e(tuned)
    local folds_used=e(folds)
    tempname GX GF Xm Ym PHI LAM FT
    matrix `GX'=e(grid_inputs)
    matrix `GF'=e(grid_frontier)
    local gridN=rowsof(`GX')
    mkmat `x' if `touse', matrix(`Xm')
    mkmat `y' if `touse', matrix(`Ym')

    mata: _csvf_all=csvf_compute(st_matrix("`GX'"),st_matrix("`GF'"),st_matrix("`Xm'"),st_matrix("`Ym'"))
    mata: _csvf_phi=_csvf_all[.,1]
    mata: _csvf_lambda=_csvf_all[.,2..cols(_csvf_all)]
    mata: _csvf_failed=sum(_csvf_phi:>=.)
    mata: st_numscalar("r(_csvf_failed_all)", _csvf_failed)
    if r(_csvf_failed_all) == `N' {
        display as error "CSVF convexification failed for every observation; no results were stored"
        exit 498
    }
    mata: st_matrix("`PHI'",_csvf_phi)
    mata: st_matrix("`LAM'",_csvf_lambda)
    mata: st_matrix("`FT'",_csvf_phi:*st_matrix("`Ym'"))

    if "`gen'"!="" {
        generate double `gen'=.
        mata: st_store(.,"`gen'","`touse'",_csvf_phi)
        label variable `gen' "CSVF output-oriented efficiency"
    }

    mata: st_numscalar("r(csvf_mean)",mean(_csvf_phi))
    mata: st_numscalar("r(csvf_failed)",_csvf_failed)
    ereturn clear
    ereturn post, esample(`touse')
    ereturn scalar N=`N'
    ereturn scalar mean_phi=r(csvf_mean)
    ereturn scalar n_failed=r(csvf_failed)
    ereturn scalar c_pen=`selected_C'
    ereturn scalar eps=`selected_eps'
    ereturn scalar d=`selected_d'
    ereturn scalar tuned=`tuned'
    ereturn scalar folds=`folds_used'
    ereturn scalar n_grid=`gridN'
    ereturn matrix frontier=`FT'
    ereturn matrix lambda=`LAM'
    ereturn matrix grid_inputs=`GX'
    ereturn matrix grid_frontier=`GF'
    ereturn local x "`x'"
    ereturn local y "`y'"
    ereturn local orientation "output"
    ereturn local returns_to_scale "vrs"
    ereturn local solver "LinearProgram"
    ereturn local cmd "csvf"
    ereturn local cmdline `"csvf `0'"'

    display as text _newline "Convexified Support Vector Frontier (CSVF)"
    display as text "  Observations: " as result `N'
    display as text "  Grid points : " as result `gridN'
    display as text "  Mean phi    : " as result %9.4f e(mean_phi)
    display as text "  LP failures : " as result e(n_failed)
end

mata:
mata set matastrict off
real rowvector csvf_one(real matrix lv_a, real matrix lv_yhat,
    real colvector x0, real colvector y0)
{
    real scalar rc,nX,nY,nJ
    real rowvector c,lb,ub
    real rowvector sol
    real matrix A,Aeq
    real colvector b,beq
    class LinearProgram scalar L
    nX=rows(x0); nY=rows(y0); nJ=rows(lv_a)
    c=(J(1,nJ,0),-1); lb=(J(1,nJ,0),1); ub=J(1,nJ+1,.)
    A=(lv_a',J(nX,1,0))\(-lv_yhat',y0)
    b=x0\J(nY,1,0)
    Aeq=(J(1,nJ,1),0); beq=1
    L=LinearProgram(); L.setMaxOrMin("min"); L.setCoefficients(c)
    L.setInequality(A,b); L.setEquality(Aeq,beq); L.setBounds(lb,ub)
    L.setTrace("off"); rc=L.optimize()
    if (L.converged() & L.errorcode()==0) {
        sol=L.parameters()
        return((-L.value(),sol[1..nJ]))
    }
    return(J(1,nJ+1,.))
}

real matrix csvf_compute(real matrix aa, real matrix ff, real matrix xx, real matrix yy)
{
    real scalar n,q,i
    real matrix ans
    real colvector x0,y0
    n=rows(yy); q=rows(aa); ans=J(n,q+1,.)
    for (i=1;i<=n;i++) {
        x0=xx[i,.]'; y0=yy[i,.]'
        ans[i,.]=csvf_one(aa,ff,x0,y0)
    }
    return(ans)
}
end
