*! svfspline.ado  v0.3.0
*! Support Vector Frontier with Splines (SVF-Spline)
*! Author: Choonjoo Lee, Korea National Defense University
*!         Nonsan-si, Republic of Korea; bloom.rampike@gmail.com
*! Method reference: Guerrero, Moragues, Aparicio, and Valero-Carreras (2024),
*!   "Support Vector Frontiers with kernel splines", Omega 128:103130.
*!   DOI: 10.1016/j.omega.2024.103130
*!
*! Estimates a concave, non-decreasing production frontier using the additive
*! linear-spline transformation and squared-L2 regularization.  For one output
*! this is model (14); for multiple outputs it implements the separable model
*! (19) with a common input grid.
*!
*! Frontier model (additive, per output r):
*!   f_r(x) = sum_j [ beta0_jr + beta_jr*x_j + sum_k beta_jkr*(x_j-t_jk)+ ]
*!
*!   Monotonicity: cumulative segment slopes are nonnegative
*!   Concavity:     hinge coefficients are nonpositive
*!
*! Convex formulation (per output r):
*!   min  beta'beta + C*sum(xi_i)
*!   s.t. y_i <= Phi[i,.]*beta <= y_i + eps + xi_i
*!        cumulative slopes >= 0; hinges <= 0; xi_i >= 0
*! Solved by Mata optimize() with a decreasing logarithmic barrier.
*!
*! Efficiency (output-oriented, >= 1):
*!   phi_i = max(1, min_r{ f_r(x_i) / y_ir })
*!
*! Syntax:
*!   svfspline, x(varlist) y(varlist) [gen(name) nknots(#) c(real) eps(real)
*!              lambda(real) seed(#) tune cgrid(numlist) epsgrid(numlist)
*!              dgrid(numlist) folds(#) solver(native|highs)
*!              pythonexec(string)]
*!
*!   x(varlist)   variable inputs
*!   y(varlist)   good outputs
*!   gen(name)    store efficiency phi in dataset (>= 1)
*!   nknots(#)    equal-width intervals d per input (default 3; d+1 nodes)
*!   c(real)      slack penalty (default 10.0)
*!   eps(real)    epsilon insensitivity (default 0.0)
*!   lambda(real) retained only as a deprecated compatibility option; ignored
*!   seed(#)      random seed (default 12345)

program define svfspline, eclass
    version 16.1
    syntax [if] [in], ///
        X(varlist numeric)     ///
        Y(varlist numeric)     ///
        [ GEN(name)            ///
          NKNots(integer 3)    ///
          C(real 10.0)         ///
          EPS(real 0.0)        ///
          LAMbda(real 0.01)    ///
          SEED(integer 12345)  ///
          TUNE                 ///
          CGRID(numlist >0)    ///
          EPSGRID(numlist >=0) ///
          DGRID(numlist integer >0) ///
          FOLDS(integer 5)     ///
          SOLVER(string)       ///
          PYTHONEXEC(string) ]

    local xvars `x'
    local yvars `y'
    local nX : word count `xvars'
    local nY : word count `yvars'
    if "`solver'" == "" local solver native
    if !inlist("`solver'", "native", "highs") {
        di as err "solver() must be native or highs"
        exit 198
    }
    if "`solver'" == "highs" & "`tune'" != "" {
        di as err "solver(highs) currently supports fixed c(), eps(), and nknots() only; do not specify tune"
        exit 198
    }
    if "`solver'" == "highs" & `"`pythonexec'"' == "" {
        di as err "solver(highs) requires pythonexec() pointing to Python with cvxpy and highspy installed"
        exit 198
    }
    if "`solver'" == "highs" {
        capture confirm file `"`pythonexec'"'
        if _rc {
            di as err "pythonexec() does not identify a readable Python executable: `pythonexec'"
            exit 601
        }
    }

    if `nknots' < 1 {
        di as err "nknots() must be >= 1"
        exit 198
    }
    if `c' <= 0 {
        di as err "c() must be > 0"
        exit 198
    }
    if `eps' < 0 {
        di as err "eps() must be >= 0"
        exit 198
    }
    if `lambda' != 0.01 {
        di as txt "note: lambda() is deprecated and ignored by the current SVF-Spline estimator"
    }
    if `folds' < 2 {
        di as err "folds() must be >= 2"
        exit 198
    }
    if "`cgrid'" == "" local cgrid ".1 .5 1 2 5"
    if "`epsgrid'" == "" local epsgrid "0 .001 .01 .1 .2"
    if "`dgrid'" == "" local dgrid "1 2 3"

    tempvar touse
    quietly mark `touse' `if' `in'
    markout `touse' `xvars' `yvars'

    if "`gen'" != "" {
        confirm new variable `gen'
    }
    foreach v of varlist `xvars' {
        quietly count if `touse' & `v' < 0
        if r(N) {
            di as err "input variable `v' contains negative values in the estimation sample"
            exit 459
        }
    }
    foreach v of varlist `yvars' {
        quietly count if `touse' & `v' <= 0
        if r(N) {
            di as err "output variable `v' must be strictly positive in the estimation sample"
            exit 459
        }
    }
    qui count if `touse'
    local N = r(N)
    if `N' < 2 {
        di as err "Too few observations."
        exit 2001
    }
    if "`tune'" != "" & `folds' > `N' {
        di as err "folds() may not exceed the number of complete observations (`N')"
        exit 198
    }

    tempname Xm Ym Bmat Fmat Rmat
    mkmat `xvars' if `touse', matrix(`Xm')
    mkmat `yvars' if `touse', matrix(`Ym')

    local solver_used optimize_log_barrier
    if "`solver'" == "highs" {
        tempfile hx hy hbeta hfront hranges hdiag
        quietly findfile svfspline_highs_bridge.py
        if _rc {
            di as err "svfspline_highs_bridge.py was not found on the ado path"
            exit 499
        }
        local hbridge `"`r(fn)'"'
        preserve
            keep if `touse'
            keep `xvars'
            export delimited using `"`hx'"', replace
        restore
        preserve
            keep if `touse'
            keep `yvars'
            export delimited using `"`hy'"', replace
        restore
        shell `pythonexec' `hbridge' --x `hx' --y `hy' --beta `hbeta' ///
            --frontier `hfront' --ranges `hranges' --diagnostics `hdiag' ///
            --d `nknots' --c `c' --eps `eps'
        foreach f in `hbeta' `hfront' `hranges' `hdiag' {
            capture confirm file `"`f'"'
            if _rc {
                di as err "HiGHS did not return a certified SVF-Spline solution; no efficiency scores were stored"
                exit 498
            }
        }
        preserve
            import delimited using `"`hbeta'"', clear varnames(1)
            mkmat _all, matrix(`Bmat')
        restore
        preserve
            import delimited using `"`hfront'"', clear varnames(1)
            mkmat _all, matrix(`Fmat')
        restore
        preserve
            import delimited using `"`hranges'"', clear varnames(1)
            mkmat _all, matrix(`Rmat')
        restore
        preserve
            import delimited using `"`hdiag'"', clear varnames(1)
            local highs_objective = objective[1]
            local highs_violation = max_violation[1]
            local highs_stationarity = max_stationarity[1]
            local highs_complementarity = max_complementarity[1]
            local highs_iterations = iterations[1]
        restore
        mata: _svfspline_res = svfspline_highs_result( ///
            st_matrix("`Bmat'"), st_matrix("`Fmat'"), st_matrix("`Rmat'"), ///
            st_matrix("`Ym'"), `c', `eps', `nknots', ///
            `highs_objective', `highs_violation', `highs_stationarity', ///
            `highs_complementarity', `highs_iterations')
        local solver_used highs_qp
    }
    else if "`tune'" != "" {
        set seed `seed'
        mata: _svfspline_res = svfspline_cv_main( ///
            st_matrix("`Xm'"), st_matrix("`Ym'"), ///
            strtoreal(tokens("`cgrid'")), strtoreal(tokens("`epsgrid'")), ///
            strtoreal(tokens("`dgrid'")), `folds', `seed')
    }
    else {
        mata: _svfspline_res = svfspline_main( ///
            st_matrix("`Xm'"), st_matrix("`Ym'"), ///
            `nknots', `c', `eps', `lambda')
    }

    mata: st_numscalar("r(_svfspline_fit_failed)", _svfspline_res.n_failed)
    if r(_svfspline_fit_failed) > 0 {
        di as err "SVF-Spline optimization failed; no efficiency scores were stored"
        exit 498
    }

    if "`gen'" != "" {
        gen double `gen' = .
        mata: st_store(., "`gen'", "`touse'", _svfspline_res.phi)
        label var `gen' "SVF-Spline efficiency (phi >= 1)"
    }

    mata: svfspline_display(_svfspline_res, "`xvars'", "`yvars'", `nknots', `c', `eps')

    ereturn clear
    ereturn post, esample(`touse')
    ereturn scalar N          = `N'
    mata: st_numscalar("r(mean_phi)",  svfspline_valid_mean(_svfspline_res.phi))
    mata: st_numscalar("r(n_failed)",  _svfspline_res.n_failed)
    mata: st_numscalar("r(n_basis)",   _svfspline_res.n_basis)
    mata: st_numscalar("r(max_violation)", _svfspline_res.max_violation)
    mata: st_numscalar("r(objective)", _svfspline_res.objective)
    mata: st_numscalar("r(max_stationarity)", _svfspline_res.max_stationarity)
    mata: st_numscalar("r(max_barrier_gap)", _svfspline_res.max_barrier_gap)
    mata: st_numscalar("r(all_stages_converged)", _svfspline_res.all_stages_converged)
    mata: st_numscalar("r(barrier_stage_warning)", _svfspline_res.barrier_stage_warning)
    mata: st_numscalar("r(numerical_warning)", _svfspline_res.numerical_warning)
    mata: st_numscalar("r(max_complementarity)", _svfspline_res.max_complementarity)
    mata: st_numscalar("r(solver_iterations)", _svfspline_res.solver_iterations)
    mata: st_numscalar("r(selected_C)", _svfspline_res.selected_C)
    mata: st_numscalar("r(selected_eps)", _svfspline_res.selected_eps)
    mata: st_numscalar("r(selected_d)", _svfspline_res.selected_d)
    mata: st_numscalar("r(n_folds)", _svfspline_res.n_folds)
    mata: st_numscalar("r(tuned)", _svfspline_res.tuned)
    mata: st_numscalar("r(cv_failures)", _svfspline_res.cv_failures)
    mata: st_matrix("`Bmat'", _svfspline_res.B)
    mata: st_matrix("`Fmat'", _svfspline_res.F)
    mata: st_matrix("`Rmat'", _svfspline_res.ranges)
    if "`tune'" != "" {
        tempname CVmat CVFoldmat CVPredmat Foldmat
        mata: st_matrix("`CVmat'", _svfspline_res.cv_results)
        mata: st_matrix("`CVFoldmat'", _svfspline_res.cv_fold_results)
        mata: st_matrix("`CVPredmat'", _svfspline_res.cv_predictions)
        mata: st_matrix("`Foldmat'", _svfspline_res.fold_id)
    }
    ereturn scalar mean_phi  = r(mean_phi)
    ereturn scalar n_failed  = r(n_failed)
    ereturn scalar n_basis   = r(n_basis)
    ereturn scalar max_violation = r(max_violation)
    ereturn scalar objective = r(objective)
    ereturn scalar max_stationarity = r(max_stationarity)
    ereturn scalar max_barrier_gap = r(max_barrier_gap)
    ereturn scalar all_stages_converged = r(all_stages_converged)
    ereturn scalar barrier_stage_warning = r(barrier_stage_warning)
    ereturn scalar numerical_warning = r(numerical_warning)
    ereturn scalar max_complementarity = r(max_complementarity)
    ereturn scalar solver_iterations = r(solver_iterations)
    ereturn matrix coefficients = `Bmat'
    ereturn matrix frontier  = `Fmat'
    matrix colnames `Rmat' = min max
    matrix rownames `Rmat' = `xvars'
    ereturn matrix xranges = `Rmat'
    if "`tune'" != "" {
        matrix colnames `CVmat' = C eps d mse
        matrix colnames `CVFoldmat' = C eps d fold mse
        matrix colnames `CVPredmat' = C eps d fold obs output prediction
        matrix colnames `Foldmat' = fold
        ereturn matrix cv_results = `CVmat'
        ereturn matrix cv_fold_results = `CVFoldmat'
        ereturn matrix cv_predictions = `CVPredmat'
        ereturn matrix fold_id    = `Foldmat'
    }
    ereturn scalar nknots    = r(selected_d)
    ereturn scalar c_pen     = r(selected_C)
    ereturn scalar eps       = r(selected_eps)
    ereturn scalar folds     = r(n_folds)
    ereturn scalar tuned     = r(tuned)
    ereturn scalar cv_failures = r(cv_failures)
    ereturn scalar n_inputs   = `nX'
    ereturn scalar n_outputs  = `nY'
    ereturn local  solver    "`solver_used'"
    ereturn local  predict   "svfspline_p"
    ereturn local  x         "`xvars'"
    ereturn local  y         "`yvars'"
    ereturn local  cmd       "svfspline"
    ereturn local  cmdline   `"svfspline `0'"'
end

mata:
mata set matastrict off

struct SVFSplineResult {
    real colvector phi       // n x 1 efficiency (>= 1)
    real matrix    F         // n x nY frontier predictions
    real matrix    B         // spline coefficients, basis x outputs
    real scalar    n_failed  // output-specific QP failures
    real scalar    n_basis   // number of basis functions used
    real scalar    max_violation
    real scalar    objective
    real matrix    ranges       // training min/max for prediction
    real scalar    max_stationarity
    real scalar    max_barrier_gap
    real scalar    all_stages_converged
    real scalar    barrier_stage_warning
    real scalar    numerical_warning
    real scalar    max_complementarity
    real scalar    solver_iterations
    real matrix    cv_results // C, eps, d, mean validation MSE
    real matrix    cv_fold_results // C, eps, d, fold, validation MSE
    real matrix    cv_predictions // C, eps, d, fold, obs, output, prediction
    real colvector fold_id
    real scalar    selected_C
    real scalar    selected_eps
    real scalar    selected_d
    real scalar    n_folds
    real scalar    tuned
    real scalar    cv_failures
}

// -----------------------------------------------------------------------
// Build spline basis matrix (N x p) for ONE input dimension
// Paper block [1,x,(x-t_1)+,...,(x-t_(d+1))+].  The d+1 nodes are
// equally spaced endpoints of the training range, as in the Omega article.
// -----------------------------------------------------------------------
real matrix svfspline_basis1_at(real colvector x, real scalar d,
                               real scalar xmin, real scalar xmax)
{
    real scalar N, k
    real matrix B
    real colvector knots, col_k

    N = rows(x)
    knots = J(d+1,1,xmin)
    for (k=1; k<=d+1; k++) {
        knots[k] = xmin + (k-1)*(xmax-xmin)/d
    }

    // Paper block: [1, x, (x-t1)+, ..., (x-t_(d+1))+].
    B = J(N,d+3,0)
    B[.,1] = J(N,1,1)
    B[.,2] = x
    for (k=1; k<=d+1; k++) {
        col_k = x :- knots[k]
        B[.,k+2] = col_k :* (col_k :> 0)
    }
    return(B)
}

real matrix svfspline_basis1(real colvector x, real scalar d)
{
    return(svfspline_basis1_at(x,d,min(x),max(x)))
}

real matrix svfspline_design_at(real matrix X, real scalar d,
                                real matrix ranges)
{
    real scalar j
    real matrix Phi
    Phi=J(rows(X),0,0)
    for (j=1;j<=cols(X);j++) {
        Phi=(Phi,svfspline_basis1_at(X[.,j],d,ranges[j,1],ranges[j,2]))
    }
    return(Phi)
}
// -----------------------------------------------------------------------
// Published squared-L2 specification solved with optimize() and a log barrier.
// Inequalities are A*p <= b and remain strictly feasible at every stage.
// -----------------------------------------------------------------------
void svfspline_barrier_eval(
    real scalar evaltype, real rowvector par,
    real matrix A, real colvector b, real scalar C,
    real scalar q, real scalar mu,
    real scalar value, real rowvector gradient, real matrix Hessian)
{
    real scalar nv
    real colvector gap, invgap
    real matrix H0
    nv = cols(par)
    gap = b-A*par'
    if (min(gap)<=0 | any(gap:>=.)) {
        value = .
        return
    }
    value = sum(par[1..q]:^2) + C*sum(par[(q+1)..nv]) - mu*sum(ln(gap))
    if (evaltype>=1) {
        invgap = 1:/gap
        gradient = (2*par[1..q],J(1,nv-q,C)) + mu*invgap'*A
    }
    if (evaltype==2) {
        H0 = J(nv,nv,0)
        H0[|1,1\q,q|] = 2*I(q)
        Hessian = H0 + mu*A'*diag((1:/gap):^2)*A
    }
}

real rowvector svfspline_barrier_optimize(
    real matrix A, real colvector b, real rowvector p0,
    real scalar C, real scalar q,
    real scalar stages_converged_out)
{
    transmorphic S
    real rowvector p, mus
    real scalar k
    p = p0
    stages_converged_out = 1
    // Stop at 1e-7: below the command's 1e-6 original-constraint tolerance.
    // Smaller barriers materially worsen conditioning on empirical data.
    mus = (1e-2,1e-3,1e-4,1e-5,1e-6,1e-7)
    for (k=1;k<=cols(mus);k++) {
        S = optimize_init()
        optimize_init_which(S,"min")
        optimize_init_evaluator(S,&svfspline_barrier_eval())
        optimize_init_evaluatortype(S,"d2")
        optimize_init_params(S,p)
        optimize_init_narguments(S,5)
        optimize_init_argument(S,1,A)
        optimize_init_argument(S,2,b)
        optimize_init_argument(S,3,C)
        optimize_init_argument(S,4,q)
        optimize_init_argument(S,5,mus[k])
        optimize_init_technique(S,"nr")
        optimize_init_conv_ptol(S,1e-10)
        optimize_init_conv_vtol(S,1e-12)
        optimize_init_conv_nrtol(S,1e-10)
        // Warm starts across six barrier levels make a large per-level cap
        // unnecessary; 100 prevents numerically flat candidates from making
        // a CV grid effectively nonterminating.
        optimize_init_conv_maxiter(S,100)
        optimize_init_conv_warning(S,"off")
        optimize_init_tracelevel(S,"none")
        p = optimize(S)
        // Most stages converge well within 100 iterations.  Retry only a
        // difficult stage from its latest strict point, avoiding a global
        // runtime penalty while preventing a merely feasible result from
        // masquerading as an optimum.
        if (!optimize_result_converged(S) & !any(p:>=.) & min(b-A*p')>0) {
            optimize_init_params(S,p)
            optimize_init_conv_maxiter(S,1000)
            p=optimize(S)
        }
        if (!optimize_result_converged(S)) stages_converged_out = 0
        // Near-singular barrier Hessians can leave optimize()'s convergence
        // flag unset even when the warm-start iterate is finite and strictly
        // feasible.  Continue the barrier path in that case; the public
        // command accepts a solution only after checking the original
        // (non-barrier) constraint residuals below.
        if (any(p:>=.) | min(b-A*p')<=0) return(J(1,cols(p0),.))
    }
    return(p)
}

void svfspline_qp_solve(
    real matrix Phi, real colvector y, real scalar C, real scalar eps,
    real scalar nX, real scalar d,
    real colvector beta_out, real colvector F_out,
    real colvector xi_out, real scalar objective_out,
    real scalar violation_out, real scalar failed_out,
    real scalar stationarity_out, real scalar barrier_gap_out,
    real scalar stages_converged_out)
{
    real scalar N,q,nv,j,k,off,shift,mu
    real matrix A,row
    real colvector b,xi,gap,gradient0,dual,dual_gradient,stationarity
    real rowvector w,p0,p

    N=rows(Phi); q=cols(Phi); nv=q+N
    A=J(0,nv,0); b=J(0,1,0)
    A=A\(-Phi,J(N,N,0)); b=b\(-y)
    A=A\( Phi,-I(N)); b=b\(y:+eps)
    A=A\(J(N,q,0),-I(N)); b=b\J(N,1,0)

    for (j=1;j<=nX;j++) {
        off=(j-1)*(d+3)
        row=J(1,nv,0); row[off+1]=-1; A=A\row; b=b\0
        row=J(1,nv,0); row[off+2]=-1; A=A\row; b=b\0
        for (k=3;k<=d+3;k++) {
            row=J(1,nv,0); row[off+k]=1; A=A\row; b=b\0
        }
        for (k=2;k<=d+3;k++) {
            row=J(1,nv,0)
            row[(off+2)..(off+k)]=J(1,k-1,-1)
            A=A\row; b=b\0
        }
    }

    // Construct a strict interior point for arbitrary input scaling.
    w=J(1,q,0)
    for (j=1;j<=nX;j++) {
        off=(j-1)*(d+3)
        w[off+1]=1
        w[off+2]=1
        w[(off+3)..(off+d+3)]=J(1,d+1,-.1/(d+1))
    }
    F_out=Phi*w'
    shift=max((0,max(y :+ 1 :- F_out)))
    if (shift>0) {
        for (j=1;j<=nX;j++) w[(j-1)*(d+3)+1]=w[(j-1)*(d+3)+1]+shift/nX
    }
    F_out=Phi*w'
    xi=rowmax((F_out-y:-eps:+1,J(N,1,1)))
    p0=(w,xi')
    if (min(b-A*p0')<=0) {
        failed_out=1; beta_out=J(q,1,.); xi_out=J(N,1,.)
        objective_out=.; violation_out=.; stationarity_out=.
        barrier_gap_out=.; stages_converged_out=0; return
    }

    p=svfspline_barrier_optimize(A,b,p0,C,q,stages_converged_out)
    if (any(p:>=.)) {
        failed_out=1; beta_out=J(q,1,.); xi_out=J(N,1,.)
        objective_out=.; violation_out=.; stationarity_out=.
        barrier_gap_out=.; stages_converged_out=0; return
    }
    w=p[1..q]; xi=p[(q+1)..nv]'; F_out=Phi*w'
    gap=b-A*p'
    violation_out=max((0,-min(gap)))
    objective_out=sum(w:^2)+C*sum(xi)
    // Approximate KKT diagnostics at the final logarithmic-barrier iterate.
    // For A*p<=b, lambda=mu/(b-A*p) and grad(f)+A'*lambda should vanish.
    mu=1e-7
    dual=mu:/gap
    gradient0=(2*w,J(1,N,C))'
    dual_gradient=A'*dual
    stationarity=gradient0+dual_gradient
    stationarity_out=max(abs(stationarity))/max((1,max(abs(gradient0)),max(abs(dual_gradient))))
    barrier_gap_out=(rows(A)*mu)/max((1,abs(objective_out)))
    beta_out=w'; xi_out=xi
    // Primal infeasibility is fatal.  The approximate barrier KKT quantities
    // are reported separately because an intermediate iteration cap can make
    // them conservative even when an independent QP oracle agrees.
    failed_out=(violation_out>1e-6)
}

// -----------------------------------------------------------------------
// Main SVF-Spline function
// -----------------------------------------------------------------------
struct SVFSplineResult scalar svfspline_main(
    real matrix X,
    real matrix Y,
    real scalar nknots,
    real scalar C,
    real scalar eps,
    real scalar lambda
    )
{
    struct SVFSplineResult scalar res
    real scalar N, nX, nY, r, k, j, n_failed_r, p, n_fail_r_temp
    real scalar phi_k, phi_kr, obj_r, viol_r, maxviol, objsum
    real scalar stat_r, gap_r, conv_r, maxstat, maxgap, allconv
    real matrix Phi, F, Bcoef, ranges
    real colvector beta_r, F_r, xi_r

    N = rows(X); nX = cols(X); nY = cols(Y)

    // Store full-sample ranges so predict uses the fitted nodes unchanged.
    ranges=J(nX,2,.)
    for (j=1; j<=nX; j++) {
        ranges[j,1]=min(X[.,j]); ranges[j,2]=max(X[.,j])
    }
    Phi=svfspline_design_at(X,nknots,ranges)
    p = cols(Phi)

    F          = J(N, nY, .)
    Bcoef      = J(p,nY,.)
    n_failed_r = 0
    maxviol    = 0
    maxstat    = 0
    maxgap     = 0
    allconv    = 1
    objsum     = 0

    for (r=1; r<=nY; r++) {
        svfspline_qp_solve(Phi,Y[.,r],C,eps,nX,nknots,
                           beta_r,F_r,xi_r,obj_r,viol_r,n_fail_r_temp,
                           stat_r,gap_r,conv_r)
        n_failed_r = n_failed_r + n_fail_r_temp
        if (n_fail_r_temp) F[.,r]=J(N,1,.)
        else {
            F[.,r]=F_r
            Bcoef[.,r]=beta_r
            maxviol=max((maxviol,viol_r))
            maxstat=max((maxstat,stat_r))
            maxgap=max((maxgap,gap_r))
            allconv=allconv & conv_r
            objsum=objsum+obj_r
        }
    }

    // Efficiency: phi_i = max(1, min_r{ F[i,r] / Y[i,r] })
    res.phi = J(N, 1, 1)
    for (k=1; k<=N; k++) {
        phi_k = 1e15
        for (r=1; r<=nY; r++) {
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
    res.B        = Bcoef
    res.n_failed = n_failed_r
    res.n_basis  = p
    res.max_violation = (n_failed_r ? . : maxviol)
    res.objective = (n_failed_r ? . : objsum)
    res.ranges = ranges
    res.max_stationarity = (n_failed_r ? . : maxstat)
    res.max_barrier_gap = (n_failed_r ? . : maxgap)
    res.all_stages_converged = (n_failed_r ? 0 : allconv)
    res.barrier_stage_warning = (n_failed_r ? 1 : (allconv==0))
    res.numerical_warning = (n_failed_r ? 1 : ///
        (maxviol>1e-6 | maxstat>1e-2 | maxgap>1e-3))
    res.max_complementarity = .
    res.solver_iterations = .
    res.cv_results = J(0,4,.)
    res.cv_fold_results = J(0,5,.)
    res.cv_predictions = J(0,7,.)
    res.fold_id = J(N,1,.)
    res.selected_C = C
    res.selected_eps = eps
    res.selected_d = nknots
    res.n_folds = 0
    res.tuned = 0
    res.cv_failures = 0
    return(res)
}

// Certified result returned by the optional external HiGHS QP bridge.
struct SVFSplineResult scalar svfspline_highs_result(
    real matrix B, real matrix F, real matrix ranges, real matrix Y,
    real scalar C, real scalar eps, real scalar d,
    real scalar objective, real scalar violation, real scalar stationarity,
    real scalar complementarity, real scalar iterations)
{
    struct SVFSplineResult scalar res
    real scalar N,nY,k,r,phi_k,phi_kr

    N=rows(Y); nY=cols(Y)
    res.phi=J(N,1,1)
    for (k=1;k<=N;k++) {
        phi_k=1e15
        for (r=1;r<=nY;r++) {
            phi_kr=F[k,r]/Y[k,r]
            if (phi_kr<phi_k) phi_k=phi_kr
        }
        res.phi[k]=max((1,phi_k))
    }
    res.F=F; res.B=B; res.ranges=ranges
    res.n_failed=0; res.n_basis=rows(B)
    res.objective=objective; res.max_violation=violation
    res.max_stationarity=stationarity
    res.max_barrier_gap=.
    res.all_stages_converged=1
    res.barrier_stage_warning=0
    res.max_complementarity=complementarity
    res.solver_iterations=iterations
    res.numerical_warning=(violation>1e-7 | stationarity>1e-7 | complementarity>1e-7)
    res.cv_results=J(0,4,.); res.cv_fold_results=J(0,5,.)
    res.cv_predictions=J(0,7,.); res.fold_id=J(N,1,.)
    res.selected_C=C; res.selected_eps=eps; res.selected_d=d
    res.n_folds=0; res.tuned=0; res.cv_failures=0
    return(res)
}

struct SVFSplineResult scalar svfspline_cv_main(
    real matrix X, real matrix Y,
    real rowvector Cgrid, real rowvector Egrid, real rowvector Dgrid,
    real scalar Kfold, real scalar seed)
{
    struct SVFSplineResult scalar res
    real scalar N,nX,nY,ic,ie,id,f,r,j,rownum,mse,err,nval,failed
    real scalar obj,viol,fail,stat,gap,stageconv
    real scalar bestrow,bestmse,Cval,Eval,dval
    real colvector foldid,ord,tr,va,beta,Ftr,xi,pred
    real matrix ranges,Ptr,Pva,cv,cvfold,cvpred

    N=rows(X); nX=cols(X); nY=cols(Y)
    if (Kfold>N) Kfold=N
    foldid=J(N,1,.)
    ord=order(runiform(N,1),1)
    for (j=1;j<=N;j++) foldid[ord[j]]=mod(j-1,Kfold)+1

    cv=J(cols(Cgrid)*cols(Egrid)*cols(Dgrid),4,.)
    cvfold=J(0,5,.)
    cvpred=J(0,7,.)
    rownum=0
    for (id=1;id<=cols(Dgrid);id++) {
        dval=Dgrid[id]
        for (ic=1;ic<=cols(Cgrid);ic++) {
            Cval=Cgrid[ic]
            for (ie=1;ie<=cols(Egrid);ie++) {
                Eval=Egrid[ie]
                rownum++
                mse=0; failed=0
                for (f=1;f<=Kfold;f++) {
                    tr=selectindex(foldid:!=f); va=selectindex(foldid:==f)
                    ranges=J(nX,2,.)
                    for (j=1;j<=nX;j++) {
                        ranges[j,1]=min(X[tr,j]); ranges[j,2]=max(X[tr,j])
                    }
                    Ptr=svfspline_design_at(X[tr,.],dval,ranges)
                    Pva=svfspline_design_at(X[va,.],dval,ranges)
                    err=0; nval=rows(va)*nY
                    for (r=1;r<=nY;r++) {
                        svfspline_qp_solve(Ptr,Y[tr,r],Cval,Eval,nX,dval,
                            beta,Ftr,xi,obj,viol,fail,stat,gap,stageconv)
                        if (fail) {
                            failed=1
                            break
                        }
                        pred=Pva*beta
                        for (j=1;j<=rows(va);j++) {
                            cvpred=cvpred\(Cval,Eval,dval,f,va[j],r,pred[j])
                        }
                        err=err+sum((Y[va,r]-pred):^2)
                    }
                    if (failed) {
                        cvfold=cvfold\(Cval,Eval,dval,f,.)
                    }
                    else {
                        cvfold=cvfold\(Cval,Eval,dval,f,err/nval)
                        mse=mse+err/nval
                    }
                }
                cv[rownum,.]=(Cval,Eval,dval,(failed ? . : mse/Kfold))
            }
        }
    }

    bestrow=.; bestmse=.
    for (j=1;j<=rows(cv);j++) {
        if (cv[j,4]<. & (bestmse>=. | cv[j,4]<bestmse-1e-14)) {
            bestrow=j; bestmse=cv[j,4]
        }
    }
    if (bestrow>=.) {
        res=svfspline_main(X,Y,Dgrid[1],Cgrid[1],Egrid[1],.01)
        res.n_failed=res.n_failed+1
        res.cv_results=cv; res.cv_fold_results=cvfold
        res.cv_predictions=cvpred; res.fold_id=foldid
        res.n_folds=Kfold; res.tuned=1
        res.cv_failures=sum(cv[.,4]:>=.)
        return(res)
    }

    res=svfspline_main(X,Y,cv[bestrow,3],cv[bestrow,1],cv[bestrow,2],.01)
    res.cv_results=cv
    res.cv_fold_results=cvfold
    res.cv_predictions=cvpred
    res.fold_id=foldid
    res.selected_C=cv[bestrow,1]
    res.selected_eps=cv[bestrow,2]
    res.selected_d=cv[bestrow,3]
    res.n_folds=Kfold
    res.tuned=1
    res.cv_failures=sum(cv[.,4]:>=.)
    return(res)
}

real scalar svfspline_valid_mean(real colvector v)
{
    real colvector w
    w = select(v, v :!= .)
    return(length(w) > 0 ? mean(w) : .)
}

void svfspline_display(struct SVFSplineResult scalar res,
                       string scalar xnames, string scalar ynames,
                       real scalar nknots, real scalar C, real scalar eps)
{
    real scalar N
    real colvector v, eff1
    N = rows(res.phi)
    v = select(res.phi, res.phi :!= .)

    printf("{txt}\n")
    printf("{txt}SVF-Spline (Guerrero et al. 2024, Omega)\n")
    printf("{hline 60}\n")
    printf("{txt}  Variable inputs (x)     : {res}%s{txt}\n", xnames)
    printf("{txt}  Good outputs (y)        : {res}%s{txt}\n", ynames)
    printf("{txt}  Observations            : {res}%g{txt}\n", N)
    printf("{txt}  Equal-width intervals d : {res}%g{txt}\n", res.selected_d)
    printf("{txt}  Basis functions         : {res}%g{txt}\n", res.n_basis)
    printf("{txt}  Parameters: C=%g  eps=%g\n", res.selected_C, res.selected_eps)
    if (res.tuned) printf("{txt}  Selection: %g-fold cross-validation\n",res.n_folds)
    if (res.n_failed==0) {
        printf("{txt}  Max primal violation    : {res}%9.2e{txt}\n",res.max_violation)
        printf("{txt}  Max KKT stationarity    : {res}%9.2e{txt}\n",res.max_stationarity)
        if (res.max_barrier_gap<.) {
            printf("{txt}  Max barrier duality gap : {res}%9.2e{txt}\n",res.max_barrier_gap)
        }
        else {
            printf("{txt}  Max KKT complementarity: {res}%9.2e{txt}\n",res.max_complementarity)
            printf("{txt}  HiGHS QP iterations     : {res}%g{txt}\n",res.solver_iterations)
        }
        if (res.barrier_stage_warning) {
            printf("{txt}  Note: at least one intermediate barrier stage reached its iteration cap;\n")
            printf("{txt}        see e(barrier_stage_warning) and inspect final diagnostics.\n")
        }
        if (res.numerical_warning) {
            printf("{err}  Warning: numerical diagnostics require review; see e(numerical_warning).{txt}\n")
        }
    }
    if (res.n_failed > 0) {
        printf("{err}  Optimization failures   : %g{txt}\n", res.n_failed)
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
