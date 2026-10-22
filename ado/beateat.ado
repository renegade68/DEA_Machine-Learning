*! beateat.ado  v0.1.1
*! Boosted EAT: Gradient Boosting with Efficiency Analysis Trees as base learners
*! Author: Choonjoo Lee, Korea National Defense University
*!         Nonsan-si, Republic of Korea; bloom.rampike@gmail.com
*! Method reference: Guillen, Aparicio, and Esteve (2023), "Gradient tree
*!   boosting and the estimation of production frontiers", Expert Systems with
*!   Applications 214:119134. DOI: 10.1016/j.eswa.2022.119134
*!
*! Gradient boosting algorithm with EAT as the weak learner.
*! Each round fits an EAT (with fixed numleaves) to the pseudo-residuals Y - F_prev,
*! then updates the frontier estimate F:
*!
*!   F_0 = colmax(Y) for all obs   (conservative upper bound)
*!   For m=1..iter:
*!     r     = Y - F_{m-1}            (residuals, <= 0)
*!     r_min = colmin(r)              (minimum per output, <= 0)
*!     r_pos = r - r_min             (shift to >= 0 for EAT)
*!     T_m   = EAT(X, r_pos, maxleaves=numleaves)
*!     h     = T_m.predict(X)        (predictions >= 0)
*!     F_m   = F_{m-1} + lr*(h + r_min)   [F decreases toward true frontier]
*!
*!   phi_i = max(1, min_r{F_iter(x_i) / y_ir})   [output-oriented, >= 1]
*!
*! Syntax:
*!   beateat, x(varlist) y(varlist) [gen(name) numleaves(#) iter(#) lr(real)
*!            numstop(#) seed(#)]
*!
*!   x(varlist)     variable inputs
*!   y(varlist)     good outputs
*!   gen(name)      store efficiency phi in dataset (>= 1)
*!   numleaves(#)   leaves per weak learner (default 5)
*!   iter(#)        boosting iterations (default 100)
*!   lr(real)       learning rate / shrinkage (default 0.3)
*!   numstop(#)     min obs per leaf in each weak learner (default 5)
*!   seed(#)        random seed (default 12345)

program define beateat, eclass
    version 16.1
    syntax [if] [in] , ///
        X(varlist numeric)    ///
        Y(varlist numeric)    ///
        [ GEN(name)           ///
          NUMLeaves(integer 5) ///
          ITER(integer 100)   ///
          LR(real 0.3)        ///
          NUMStop(integer 5)  ///
          SEED(integer 12345) ]

    if "`gen'" != "" confirm new variable `gen'

    local xvars `x'
    local yvars `y'
    local nX : word count `xvars'
    local nY : word count `yvars'

    if `numleaves' < 1 {
        di as err "numleaves() must be >= 1"
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
    if `N' < 4 {
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
    tempvar ypos
    generate byte `ypos' = 0
    foreach v of varlist `yvars' {
        quietly count if `touse' & `v' < 0
        if r(N) {
            di as err "output variables must be nonnegative"
            exit 459
        }
        quietly replace `ypos' = 1 if `touse' & `v' > 0
    }
    quietly count if `touse' & !`ypos'
    if r(N) {
        di as err "each observation must have at least one strictly positive output"
        exit 459
    }

    tempname Xm Ym
    mkmat `xvars' if `touse', matrix(`Xm')
    mkmat `yvars' if `touse', matrix(`Ym')

    mata: _beateat_res = beateat_main( ///
        st_matrix("`Xm'"), st_matrix("`Ym'"), ///
        `numstop', `numleaves', `iter', `lr')

    if "`gen'" != "" {
        quietly gen double `gen' = . if `touse'
        mata: st_store(., "`gen'", "`touse'", _beateat_res.phi)
        label var `gen' "Boosted EAT efficiency (phi >= 1)"
    }

    mata: beateat_display(_beateat_res, "`xvars'", "`yvars'", `numleaves', `iter', `lr')

    ereturn clear
    ereturn post, esample(`touse')
    ereturn scalar N          = `N'
    mata: st_numscalar("r(mean_phi)", beateat_valid_mean(_beateat_res.phi))
    mata: st_numscalar("r(n_iter)",   _beateat_res.n_iter)
    mata: st_numscalar("r(n_leaves)", _beateat_res.n_leaves)
    ereturn scalar mean_phi  = r(mean_phi)
    ereturn scalar n_iter    = r(n_iter)
    ereturn scalar n_leaves  = r(n_leaves)
    ereturn scalar lr        = `lr'
    tempname Fmat
    mata: st_matrix("`Fmat'", _beateat_res.F_final)
    ereturn matrix frontier = `Fmat'
    ereturn local  x         "`xvars'"
    ereturn local  y         "`yvars'"
    ereturn local  cmd       "beateat"
    ereturn local  cmdline   `"beateat `0'"'
end


mata:
mata set matastrict off

// ======================================================================
// EAT functions (prefixed beateat_eat_)
// ======================================================================

real scalar beateat_eat_posNode(real matrix nm, real scalar id)
{
    real colvector v
    v = selectindex(nm[.,1] :== id)
    return(length(v) ? v[1] : 0)
}

real scalar beateat_eat_isFinal(real matrix Xa, real colvector obs, real scalar ns)
{
    real matrix Xo
    real scalar j
    if (length(obs) <= ns) return(1)
    Xo = Xa[obs,.]
    for (j=1; j<=cols(Xo); j++) {
        if (length(uniqrows(Xo[.,j])) > 1) return(0)
    }
    return(1)
}

real scalar beateat_eat_mse(real matrix Y, real colvector obs,
                              real rowvector yhat, real scalar N, real scalar nY)
{
    real scalar err, r
    if (length(obs)==0) return(.)
    err = 0
    for (r=1; r<=nY; r++) err = err + sum((yhat[r] :- Y[obs,r]):^2)
    return(err/(N*nY))
}

real scalar beateat_eat_comparePareto(real rowvector a1, real rowvector b1,
                                       real rowvector a2, real rowvector b2)
{
    real rowvector b1f, b2f
    b1f = editvalue(b1, ., 1e15)
    b2f = editvalue(b2, ., 1e15)
    if (mreldif(a1,a2)<1e-12 & mreldif(b1f,b2f)<1e-12) return(0)
    if (all(a1 :< b2f)) return(-1)
    if (all(a2 :< b1f)) return(1)
    return(0)
}

void beateat_eat_estimNode(
    real matrix  X,
    real matrix  Y,
    real matrix  oleaves,
    real matrix  prow,
    real scalar  xi_sp, real scalar s_sp,
    real scalar  N, real scalar nX, real scalar nY,
    real colvector pobs,
    real matrix    tRes,
    real colvector oL, real colvector oR
    )
{
    real scalar  ca, cb, cy, j, r, maxL_r, n_ol
    real colvector obsL, obsR
    real rowvector a, b, yh, aL, bL, aR, bR, yhL, yhR, yInf, la, lb

    ca = 8; cb = 8+nX; cy = 8+2*nX

    obsL = select(pobs, X[pobs,xi_sp] :<  s_sp)
    obsR = select(pobs, X[pobs,xi_sp] :>= s_sp)

    tRes = prow \ prow
    oL   = J(0,1,.); oR = J(0,1,.)

    if (length(obsL)==0 | length(obsR)==0) {
        tRes[1,7] = .; tRes[2,7] = .
        return
    }

    a  = prow[1, ca..ca+nX-1]
    b  = prow[1, cb..cb+nX-1]
    yh = prow[1, cy..cy+nY-1]

    aL = a; bL = b; bL[xi_sp] = s_sp
    aR = a; bR = b; aR[xi_sp] = s_sp

    yInf = J(1,nY,0)
    n_ol = rows(oleaves)
    for (j=1; j<=n_ol; j++) {
        la = oleaves[j, ca..ca+nX-1]
        lb = oleaves[j, cb..cb+nX-1]
        if (beateat_eat_comparePareto(aL,bL,la,lb)==1) {
            for (r=1; r<=nY; r++) {
                if (yInf[r] < oleaves[j,cy+r-1]) yInf[r] = oleaves[j,cy+r-1]
            }
        }
    }

    yhL = J(1,nY,0)
    for (r=1; r<=nY; r++) {
        maxL_r = max(Y[obsL,r])
        yhL[r] = (maxL_r >= yInf[r]) ? maxL_r : yInf[r]
    }
    yhR = yh

    tRes[1,1]=.; tRes[1,2]=prow[1,1]; tRes[1,3..6]=J(1,4,-1)
    tRes[1, ca..ca+nX-1]=aL; tRes[1, cb..cb+nX-1]=bL
    tRes[1, cy..cy+nY-1]=yhL; tRes[1,7]=beateat_eat_mse(Y,obsL,yhL,N,nY)

    tRes[2,1]=.; tRes[2,2]=prow[1,1]; tRes[2,3..6]=J(1,4,-1)
    tRes[2, ca..ca+nX-1]=aR; tRes[2, cb..cb+nX-1]=bR
    tRes[2, cy..cy+nY-1]=yhR; tRes[2,7]=beateat_eat_mse(Y,obsR,yhR,N,nY)

    oL = obsL; oR = obsR
}

void beateat_eat_grow(
    real matrix X, real matrix Y,
    real scalar numStop, real scalar maxLeaves,
    real scalar N, real scalar nX, real scalar nY,
    real matrix final_nm
    )
{
    real scalar  K, ca, cb, cy
    real scalar  next_id, n_nodes, leaf_count, qlen
    real scalar  xi, i, err, err_min, best_xi, best_s
    real scalar  id_L, id_R, pidx
    real scalar  iL, iR
    real matrix  nm, prow, tRes, oleaves, btRes
    real colvector pobs, obsL, obsR, boL, boR, queue, aobs, Sx
    transmorphic obs_db

    K = 7+2*nX+nY; ca = 8; cb = 8+nX; cy = 8+2*nX

    nm = J(1,K,.)
    nm[1,1]=1; nm[1,2]=-1; nm[1,3]=-1; nm[1,4]=-1
    nm[1,5]=-1; nm[1,6]=-1
    nm[1, ca..ca+nX-1] = colmin(X)
    nm[1, cb..cb+nX-1] = J(1,nX,.)
    nm[1, cy..cy+nY-1] = colmax(Y)

    obs_db = asarray_create("real")
    aobs   = range(1,N,1)
    asarray(obs_db, 1, aobs)
    nm[1,7] = beateat_eat_mse(Y, aobs, nm[1,cy..cy+nY-1], N, nY)

    n_nodes=1; next_id=2; leaf_count=1
    queue = J(1,1,1)

    while (1) {
        qlen = length(queue)
        if (qlen==0) break
        if (maxLeaves>0 & leaf_count>=maxLeaves) break

        pidx  = queue[qlen]
        if (qlen>1) {
            queue = queue[1..qlen-1]
        } else {
            queue = J(0,1,.)
        }
        prow  = nm[pidx,.]
        pobs  = asarray(obs_db, pidx)

        if (beateat_eat_isFinal(X, pobs, numStop)) continue

        oleaves = select(nm, (nm[.,3]:==-1) :& (nm[.,1]:!=prow[1,1]))

        err_min=1e300; best_xi=-1; best_s=.
        btRes=J(2,K,.); boL=J(0,1,.); boR=J(0,1,.)

        for (xi=1; xi<=nX; xi++) {
            Sx = sort(uniqrows(X[pobs,xi]), 1)
            if (length(Sx)<=1) continue
            for (i=2; i<=length(Sx); i++) {
                beateat_eat_estimNode(X,Y,oleaves,prow,xi,Sx[i],N,nX,nY,pobs,
                                      tRes,obsL,obsR)
                if (tRes[1,7]==. | tRes[2,7]==.) continue
                err = tRes[1,7]+tRes[2,7]
                if (err < err_min) {
                    err_min=err; best_xi=xi; best_s=Sx[i]
                    btRes=tRes; boL=obsL; boR=obsR
                }
            }
        }

        if (best_xi==-1) continue

        id_L=next_id; id_R=next_id+1; next_id=next_id+2
        nm[pidx,3]=id_L; nm[pidx,4]=id_R; nm[pidx,5]=best_xi; nm[pidx,6]=best_s

        btRes[1,1]=id_L; btRes[1,2]=prow[1,1]
        btRes[2,1]=id_R; btRes[2,2]=prow[1,1]
        nm = nm \ btRes
        n_nodes = n_nodes+2

        asarray(obs_db, n_nodes-1, boL)
        asarray(obs_db, n_nodes,   boR)

        iL=n_nodes-1; iR=n_nodes

        if (beateat_eat_isFinal(X,boR,numStop)) {
            queue=(iR\queue)
        } else {
            queue=(queue\iR)
        }
        if (beateat_eat_isFinal(X,boL,numStop)) {
            queue=(iL\queue)
        } else {
            queue=(queue\iL)
        }
        leaf_count++
    }

    final_nm = nm
}

real rowvector beateat_eat_predictor(real matrix nm, real rowvector xobs,
                                      real scalar nX, real scalar nY)
{
    real scalar ti, next_ti
    ti = 1
    while (nm[ti,3] != -1) {
        if (xobs[nm[ti,5]] < nm[ti,6]) {
            next_ti = beateat_eat_posNode(nm, nm[ti,3])
        } else {
            next_ti = beateat_eat_posNode(nm, nm[ti,4])
        }
        ti = next_ti
    }
    return(nm[ti, (8+2*nX)..(7+2*nX+nY)])
}

real matrix beateat_predict_all(real matrix nm, real matrix X,
                                 real scalar nX, real scalar nY)
{
    real scalar  k, N
    real matrix  pred
    N = rows(X)
    pred = J(N, nY, 0)
    for (k=1; k<=N; k++) {
        pred[k,.] = beateat_eat_predictor(nm, X[k,.], nX, nY)
    }
    return(pred)
}

// ======================================================================
// Boosted EAT specific Mata code
// ======================================================================

struct BEATResult {
    real matrix    F_final   // n x nY final frontier predictions
    real colvector phi       // n x 1 output-oriented efficiency (>= 1)
    real scalar    n_iter    // boosting rounds
    real scalar    n_leaves  // leaves per weak learner
}

struct BEATResult scalar beateat_main(
    real matrix X,        // n x nX
    real matrix Y,        // n x nY
    real scalar numStop,
    real scalar numLeaves,
    real scalar iter,
    real scalar lr
    )
{
    struct BEATResult scalar res
    real scalar N, nX, nY, m, k, r
    real matrix F, r_mat, r_pos, h, nm
    real rowvector r_min
    real colvector phi
    real scalar phi_k, ratio

    N = rows(X); nX = cols(X); nY = cols(Y)

    // Step 1: Initialize F = colmax(Y) for all observations
    F = J(N, nY, 0)
    for (r=1; r<=nY; r++) F[.,r] = J(N,1,max(Y[.,r]))

    // Step 2: Boosting rounds
    for (m=1; m<=iter; m++) {
        r_mat = Y - F                      // N x nY residuals (<= 0)
        r_min = colmin(r_mat)              // 1 x nY (<= 0)
        r_pos = r_mat :- r_min            // N x nY (>= 0)

        // Fit EAT with maxLeaves = numLeaves to (X, r_pos)
        beateat_eat_grow(X, r_pos, numStop, numLeaves, N, nX, nY, nm)

        // Predict for all obs
        h = beateat_predict_all(nm, X, nX, nY)  // N x nY

        // Update frontier: F += lr * (h + r_min)
        F = F + lr :* (h :+ r_min)
    }

    // Step 3: Compute output-oriented efficiency
    // phi_i = max(1, min_r{ F_r(x_i) / Y_ir })
    phi = J(N, 1, 1)
    for (k=1; k<=N; k++) {
        phi_k = 1e15
        for (r=1; r<=nY; r++) {
            if (Y[k,r] > 1e-12) {
                ratio = F[k,r] / Y[k,r]
                if (ratio < phi_k) phi_k = ratio
            }
        }
        if (phi_k == 1e15 | phi_k < 1) phi_k = 1
        phi[k] = phi_k
    }

    res.F_final  = F
    res.phi      = phi
    res.n_iter   = iter
    res.n_leaves = numLeaves
    return(res)
}

real scalar beateat_valid_mean(real colvector v)
{
    real colvector w
    w = select(v, v :!= .)
    return(length(w) > 0 ? mean(w) : .)
}

void beateat_display(struct BEATResult scalar res,
                     string scalar xnames, string scalar ynames,
                     real scalar numLeaves, real scalar iter, real scalar lr)
{
    real scalar N
    real colvector v, eff1
    N = rows(res.phi)
    v = select(res.phi, res.phi :!= .)

    printf("{txt}\n")
    printf("{txt}Boosted EAT (Guillen, Aparicio, Esteve 2023, ESWA)\n")
    printf("{hline 60}\n")
    printf("{txt}  Variable inputs (x)     : {res}%s{txt}\n", xnames)
    printf("{txt}  Good outputs (y)        : {res}%s{txt}\n", ynames)
    printf("{txt}  Observations            : {res}%g{txt}\n", N)
    printf("{txt}  Boosting iterations     : {res}%g{txt}\n", iter)
    printf("{txt}  Leaves per weak learner : {res}%g{txt}\n", numLeaves)
    printf("{txt}  Learning rate           : {res}%g{txt}\n", lr)
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
