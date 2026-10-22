*! lsbmafs.ado  v0.1.1
*! Least Squares Boosting with MARS (LSB-MAFS) Production Frontier
*! Author: Choonjoo Lee, Korea National Defense University
*!         Nonsan-si, Republic of Korea; bloom.rampike@gmail.com
*! Method reference: Guillen, Aparicio, Zofio, and Espana (2024), "Improving
*!   the predictive accuracy of production frontier models for efficiency
*!   measurement using machine learning: The LSB-MAFS method",
*!   Computers & Operations Research 171:106793. DOI: 10.1016/j.cor.2024.106793
*!
*! Estimates a production frontier using gradient boosting (L2 loss)
*! with Multivariate Adaptive Regression Splines (MARS) as weak learners.
*!
*! Implementation of published Algorithm 2:
*!   F_0 = max(y); r = y-F; fit a constrained reflected-pair MAFS
*!   learner directly to r; F = F + lr*h.
*!
*!   phi_i = F(x_i) / y_i  (output-oriented efficiency)
*!
*! Each reflected pair uses (x_j-t)+ and (t-x_j)+ with LP constraints
*! alpha>=0, beta<=0, and alpha+beta<=0, plus envelopment slacks.
*!
*! Syntax:
*!   lsbmafs, x(varlist) y(varlist) [gen(name) maxterms(#) iter(#)
*!             lr(real) numstop(#) seed(#)]
*!
*!   x(varlist)    variable inputs
*!   y(varname)    one good output (as in the published estimator)
*!   gen(name)     store efficiency phi in dataset (>= 1)
*!   maxterms(#)   max MARS basis terms per weak learner (default 5)
*!   iter(#)       boosting iterations (default 100)
*!   lr(real)      shrinkage / learning rate (default 0.1)
*!   numstop(#)    optional MAFS minimum-span override (default automatic)
*!   seed(#)       random seed (default 12345)

program define lsbmafs, eclass
    version 16.1
    syntax [if] [in] , ///
        X(varlist numeric)      ///
        Y(varlist numeric)      ///
        [ GEN(name)             ///
          MAXTerms(integer 5)   ///
          ITER(integer 100)     ///
          LR(real 0.1)          ///
          NUMStop(integer 0)    ///
          SEED(integer 12345)   ]

    if "`gen'" != "" confirm new variable `gen'

    local xvars `x'
    local yvars `y'
    local nX : word count `xvars'
    local nY : word count `yvars'

    if `nY' != 1 {
        di as err "lsbmafs requires exactly one output variable"
        di as err "The published LSB-MAFS estimator and its reference implementation are single-output."
        exit 198
    }

    if `maxterms' < 1 {
        di as err "maxterms() must be >= 1"
        exit 198
    }
    if `iter' < 1 {
        di as err "iter() must be >= 1"
        exit 198
    }
    if `lr' <= 0 | `lr' > 1 {
        di as err "lr() must be in (0,1]"
        exit 198
    }
    if `numstop' < 0 {
        di as err "numstop() must be >= 0"
        exit 198
    }

    tempvar touse
    marksample touse, novarlist
    markout `touse' `xvars' `yvars'
    qui count if `touse'
    local N = r(N)
    if `N' < 2 {
        di as err "Too few observations."
        exit 2001
    }
    foreach v of varlist `xvars' {
        quietly count if `touse' & `v' < 0
        if r(N) {
            di as err "input variables must be nonnegative"
            exit 459
        }
    }
    quietly count if `touse' & `yvars' <= 0
    if r(N) {
        di as err "the output variable must be strictly positive for radial efficiency"
        exit 459
    }

    tempname Xm Ym
    mkmat `xvars' if `touse', matrix(`Xm')
    mkmat `yvars' if `touse', matrix(`Ym')

    mata: _lsbmafs_res = lsbmafs_main( ///
        st_matrix("`Xm'"), st_matrix("`Ym'"), ///
        `maxterms', `iter', `lr', `numstop')
    mata: st_numscalar("r(_lsb_final_failed)", _lsbmafs_res.n_final_failed)
    if r(_lsb_final_failed) > 0 {
        di as err "LSB-MAFS selected-basis optimization failed; no results were stored"
        exit 498
    }

    if "`gen'" != "" {
        quietly gen double `gen' = . if `touse'
        mata: st_store(., "`gen'", "`touse'", _lsbmafs_res.phi)
        label var `gen' "LSB-MAFS efficiency (phi >= 1)"
    }

    mata: lsbmafs_display(_lsbmafs_res, "`xvars'", "`yvars'", ///
        `maxterms', `iter', `lr')

    ereturn clear
    ereturn post, esample(`touse')
    ereturn scalar N          = `N'
    mata: st_numscalar("r(mean_phi)", lsbmafs_valid_mean(_lsbmafs_res.phi))
    mata: st_numscalar("r(n_lp_failed)", _lsbmafs_res.n_lp_failed)
    mata: st_numscalar("r(n_final_failed)", _lsbmafs_res.n_final_failed)
    ereturn scalar mean_phi  = r(mean_phi)
    ereturn scalar n_lp_failed = r(n_lp_failed)
    ereturn scalar n_final_failed = r(n_final_failed)
    ereturn scalar n_iter    = `iter'
    ereturn scalar maxterms  = `maxterms'
    ereturn scalar lr        = `lr'
    tempname Fmat
    mata: st_matrix("`Fmat'", _lsbmafs_res.F)
    ereturn matrix frontier = `Fmat'
    ereturn local  x         "`xvars'"
    ereturn local  y         "`yvars'"
    ereturn local  solver    "LinearProgram"
    ereturn local  cmd       "lsbmafs"
    ereturn local  cmdline   `"lsbmafs `0'"'
end


mata:
mata set matastrict off

struct LSBMAFSResult {
    real colvector phi       // n x 1 efficiency (>= 1)
    real matrix    F         // n x nY frontier predictions
    real scalar    n_iter    // boosting iterations run
    real scalar    maxterms  // max MARS terms per weak learner
    real scalar    n_lp_failed     // all failed candidate/final LP solves
    real scalar    n_final_failed  // failed selected-basis LP solves
}

void lsbmafs_lp_fit(real matrix B, real colvector y,
                    real colvector coef, real colvector pred,
                    real scalar ok)
{
    real scalar N,p,q,g,rc,nv
    real rowvector c,lb,ub
    real matrix Aeq,Aineq
    real colvector beq,bineq,sol
    class LinearProgram scalar S
    N=rows(B); p=cols(B); q=(p-1)/2; nv=p+N
    c=(J(1,p,0),J(1,N,1))
    Aeq=(B,-I(N)); beq=y
    Aineq=J(3*q,nv,0); bineq=J(3*q,1,0)
    for (g=1;g<=q;g++) {
        // columns 2g and 2g+1 are right (alpha) and left (beta).
        Aineq[3*g-2,2*g]=1; Aineq[3*g-2,2*g+1]=1 // alpha+beta<=0
        Aineq[3*g-1,2*g]=-1                       // alpha>=0
        Aineq[3*g,2*g+1]=1                        // beta<=0
    }
    lb=(J(1,p,.),J(1,N,0)); ub=J(1,nv,.)
    S=LinearProgram(); S.setMaxOrMin("min"); S.setCoefficients(c)
    S.setEquality(Aeq,beq)
    if (q>0) S.setInequality(Aineq,bineq)
    S.setBounds(lb,ub); S.setTrace("off"); rc=S.optimize()
    ok=(S.converged() & S.errorcode()==0)
    if (!ok) {
        coef=J(p,1,.)
        pred=J(N,1,.)
        return
    }
    sol=S.parameters()'; coef=sol[1..p]; pred=B*coef
}

real colvector lsbmafs_mars_fit(
    real matrix X,
    real colvector r,
    real scalar maxterms,
    real scalar numstop,
    real scalar n_lp_failed,
    real scalar final_failed
    )
{
    real scalar N,nX,j,k,best_j,best_k,err,err_try,best_err,ok,Le,L,idx0
    real matrix B, Btry
    real colvector coef,pred,pred_try,right,left,index,eligible,knots,used

    N = rows(X); nX = cols(X)
    B = J(N, 1, 1)
    coef=max(r); pred=J(N,1,coef[1]); err=round(mean((r-pred):^2),.0001)
    used=J(N,nX,0)
    n_lp_failed=0; final_failed=0
    Le=ceil(3-ln(.05/nX)/ln(2)); if (Le<1) Le=1

    // Match boostingDEA's forward stopping rule and total-basis convention.
    while (cols(B)+2 < maxterms) {
        best_err=err; best_j=0; best_k=0
        index=J(N,1,1)
        if (2*Le>=N) index=J(N,1,0)
        else {
            index[1..Le]=J(Le,1,0)
            index[(N-Le+1)..N]=J(Le,1,0)
        }
        for (j=1;j<=nX;j++) {
            eligible=index
            if (numstop>0) L=numstop
            else L=ceil(-(ln(-(1/(nX*sum(index)))*ln(1-.05))/ln(2))/2.5)
            for (idx0=1;idx0<=N;idx0++) if (used[idx0,j]) {
                eligible[max((1,idx0-L))..min((N,idx0+L))]=J(min((N,idx0+L))-max((1,idx0-L))+1,1,0)
            }
            knots=uniqrows(sort(select(X[.,j],eligible),1))
            if (rows(knots)<=1) continue
            for (k=1;k<=rows(knots);k++) {
                right=(X[.,j]:-knots[k]):*(X[.,j]:>knots[k])
                left =(knots[k]:-X[.,j]):*(X[.,j]:<knots[k])
                Btry=(B,right,left)
                lsbmafs_lp_fit(Btry,r,coef,pred_try,ok)
                if (!ok) {
                    n_lp_failed++
                    continue
                }
                err_try=round(mean((r-pred_try):^2),.0001)
                if (err_try<best_err) {
                    best_err=err_try; best_j=j
                    best_k=selectindex(X[.,j]:==knots[k])[1]
                }
            }
        }
        if (best_j==0 | best_err>=err*(1-.01)) break
        right=(X[.,best_j]:-X[best_k,best_j]):*(X[.,best_j]:>X[best_k,best_j])
        left =(X[best_k,best_j]:-X[.,best_j]):*(X[.,best_j]:<X[best_k,best_j])
        B=(B,right,left); used[best_k,best_j]=1; err=best_err
    }
    lsbmafs_lp_fit(B,r,coef,pred,ok)
    if (!ok) {
        n_lp_failed++
        final_failed=1
        pred=J(N,1,.)
    }
    return(pred)
}

// -----------------------------------------------------------------------
// Main LSB-MAFS function
// -----------------------------------------------------------------------
struct LSBMAFSResult scalar lsbmafs_main(
    real matrix X,
    real matrix Y,
    real scalar maxterms,
    real scalar iter,
    real scalar lr,
    real scalar numstop
    )
{
    struct LSBMAFSResult scalar res
    real scalar N, nX, nY, m, r, k
    real scalar phi_k, phi_kr
    real scalar lp_failed_m, final_failed_m, n_lp_failed, n_final_failed
    real matrix F, r_mat
    real colvector h, y_r

    N = rows(X); nX = cols(X); nY = cols(Y)

    // F_0 = max(y) for each output (tight constant upper bound)
    F = J(N, nY, 0)
    for (r = 1; r <= nY; r++) {
        F[.,r] = J(N, 1, max(Y[.,r]))
    }
    n_lp_failed = 0
    n_final_failed = 0

    // Boosting loop
    for (m = 1; m <= iter; m++) {
        for (r = 1; r <= nY; r++) {
            if (any(missing(F[.,r]))) continue
            y_r   = Y[.,r]
            r_mat = y_r - F[.,r]
            h = lsbmafs_mars_fit(X, r_mat, maxterms, numstop,
                                 lp_failed_m, final_failed_m)
            n_lp_failed = n_lp_failed + lp_failed_m
            n_final_failed = n_final_failed + final_failed_m
            if (final_failed_m) F[.,r] = J(N,1,.)
            else F[.,r] = F[.,r] :+ lr :* h
        }
    }

    // Efficiency: phi_i = max(1, min_r{ F[i,r] / Y[i,r] })
    res.phi = J(N, 1, .)
    for (k = 1; k <= N; k++) {
        if (any(missing(F[k,.]))) continue
        phi_k = 1e15
        for (r = 1; r <= nY; r++) {
            if (Y[k,r] > 1e-12) {
                phi_kr = F[k,r] / Y[k,r]
                if (phi_kr < phi_k) phi_k = phi_kr
            }
        }
        if (phi_k == 1e15) phi_k = 1
        if (phi_k < 1)     phi_k = 1
        res.phi[k] = phi_k
    }

    res.F        = F
    res.n_iter   = iter
    res.maxterms = maxterms
    res.n_lp_failed = n_lp_failed
    res.n_final_failed = n_final_failed
    return(res)
}

real scalar lsbmafs_valid_mean(real colvector v)
{
    real colvector w
    w = select(v, v :!= .)
    return(length(w) > 0 ? mean(w) : .)
}

void lsbmafs_display(struct LSBMAFSResult scalar res,
                     string scalar xnames, string scalar ynames,
                     real scalar maxterms, real scalar iter, real scalar lr)
{
    real scalar N
    real colvector v, eff1
    N = rows(res.phi)
    v = select(res.phi, res.phi :!= .)

    printf("{txt}\n")
    printf("{txt}LSB-MAFS (Guillen, Aparicio, Zofio, España 2024, C&OR)\n")
    printf("{hline 60}\n")
    printf("{txt}  Variable inputs (x)     : {res}%s{txt}\n", xnames)
    printf("{txt}  Good outputs (y)        : {res}%s{txt}\n", ynames)
    printf("{txt}  Observations            : {res}%g{txt}\n", N)
    printf("{txt}  MARS max terms          : {res}%g{txt}\n", maxterms)
    printf("{txt}  Boosting iterations     : {res}%g{txt}\n", iter)
    printf("{txt}  Learning rate           : {res}%g{txt}\n", lr)
    if (res.n_lp_failed > 0) {
        printf("{err}  LP solve failures        : %g (selected-basis: %g){txt}\n",
               res.n_lp_failed, res.n_final_failed)
    }
    printf("{hline 60}\n")
    if (length(v) > 0) {
        eff1 = select(res.phi, abs(res.phi :- 1) :<= 1e-6)
        printf("{txt}  Output-oriented efficiency (phi, >= 1):\n")
        printf("{txt}    Mean= %7.4f  Min= %7.4f  Max= %7.4f\n",
               mean(v), min(v), max(v))
        printf("{txt}    Efficient DMUs (phi~=1): %g\n", length(eff1))
    }
    printf("{hline 60}\n")
}

end
