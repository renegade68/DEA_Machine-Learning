*! rfeat.ado  v0.1.3
*! Random Forest Efficiency Analysis Trees (RFEAT) for Stata
*! Author: Choonjoo Lee, Korea National Defense University
*!         Nonsan-si, Republic of Korea; bloom.rampike@gmail.com
*! Method reference: Esteve M, Aparicio J, Rodriguez-Sala JJ, Zhu J (2023),
*!   "Random Forests and the measurement of super-efficiency in the
*!   context of Free Disposal Hull", EJOR 304(2):729-744.
*!   DOI: 10.1016/j.ejor.2022.04.024
*! Replicates R package 'eat' rfeat() v0.1.4
*!
*! Syntax:
*!   rfeat , x(varlist) y(varlist) [ntrees(#) numstop(#) s(#)
*!           maxdepth(#) maxleaves(#) gen(newvar) gensd(newvar) seed(#)]

program define rfeat, eclass
    version 16.1
    syntax [if] [in] , ///
        X(varlist numeric) ///
        Y(varlist numeric) ///
        [ NTrees(integer 200)   ///  number of bootstrap trees
          NUMStop(integer 5)    ///  min obs to attempt a split
          S(integer 0)          ///  features per split (0 = all = pure bagging)
          MAXDepth(integer 0)   ///  max tree depth (0 = unlimited)
          MAXLeaves(integer 0)  ///  max leaves (0 = unlimited)
          GENerate(name)        ///  variable name for mean efficiency scores
          GENsd(name)           ///  variable name for SD of efficiency scores
          SEED(integer 12345)   ]

    if "`generate'" != "" confirm new variable `generate'
    if "`gensd'" != "" confirm new variable `gensd'
    if "`generate'" != "" & "`generate'" == "`gensd'" {
        di as err "generate() and gensd() must name distinct new variables"
        exit 198
    }

    local xvars `x'
    local yvars `y'
    local nX : word count `xvars'
    local nY : word count `yvars'

    // -- validation --------------------------------------------------
    if `ntrees' < 1 {
        di as err "ntrees() must be >= 1"
        exit 198
    }
    if `numstop' < 0 {
        di as err "numstop() must be >= 0"
        exit 198
    }
    if `maxdepth' < 0 | `maxleaves' < 0 {
        di as err "maxdepth() and maxleaves() must be >= 0"
        exit 198
    }
    if `s' < 0 | `s' > `nX' {
        di as err "s() must be between 0 (all features) and `nX'"
        exit 198
    }
    if `s' == 0 local s `nX'   // 0 means all features (pure bagging)

    // -- sample ------------------------------------------------------
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

    // -- matrices ----------------------------------------------------
    tempname Xm Ym
    mkmat `xvars' if `touse', matrix(`Xm')
    mkmat `yvars' if `touse', matrix(`Ym')

    // -- run RFEAT in Mata --------------------------------------------
    set seed `seed'
    mata: _rfeat_res = rfeat_main( ///
        st_matrix("`Xm'"), st_matrix("`Ym'"), ///
        `ntrees', `numstop', `s', `maxdepth', `maxleaves')

    // -- store results in dataset -------------------------------------
    if "`generate'" != "" {
        qui gen double `generate' = . if `touse'
        mata: st_store(., "`generate'", "`touse'", _rfeat_res.eff_mean)
        label var `generate' "RFEAT mean efficiency"
    }
    if "`gensd'" != "" {
        qui gen double `gensd' = . if `touse'
        mata: st_store(., "`gensd'", "`touse'", _rfeat_res.eff_sd)
        label var `gensd' "RFEAT SD of efficiency"
    }

    // -- display ------------------------------------------------------
    mata: rfeat_display(_rfeat_res, "`xvars'", "`yvars'", `ntrees', `s')

    // -- ereturn ------------------------------------------------------
    ereturn clear
    ereturn post, esample(`touse')
    ereturn scalar N       = `N'
    ereturn scalar ntrees  = `ntrees'
    ereturn scalar s       = `s'
    mata: st_numscalar("r(mean_phi)", mean(_rfeat_res.eff_mean))
    ereturn scalar mean_phi = r(mean_phi)
    ereturn local  x   "`xvars'"
    ereturn local  y   "`yvars'"
    ereturn local  cmd "rfeat"
    ereturn local  cmdline `"rfeat `0'"'
    tempname Fmat
    mata: st_matrix("`Fmat'", _rfeat_res.frontier)
    ereturn matrix frontier = `Fmat'
end


// ======================================================================
mata:
mata set matastrict off

// -- result struct -----------------------------------------------------
struct RFEATResult {
    real colvector eff_mean   // mean efficiency across B trees
    real colvector eff_sd     // SD of efficiency across B trees
    real matrix    frontier   // componentwise mean forest prediction
    real scalar    ntrees     // number of trees actually fitted
    real scalar    s_feat     // features per split used
}

// -- shared helpers (identical to eatree.ado; safe to redefine) ---------

real scalar rfeat_posNode(real matrix nm, real scalar id)
{
    real colvector v
    v = selectindex(nm[.,1] :== id)
    return(length(v) ? v[1] : 0)
}

real scalar rfeat_isFinal(real matrix Xa, real colvector obs, real scalar ns)
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

real scalar rfeat_mse(real matrix Y, real colvector obs,
                      real rowvector yhat, real scalar N, real scalar nY)
{
    real scalar err, r
    if (length(obs)==0) return(.)
    err = 0
    for (r=1; r<=nY; r++) err = err + sum((yhat[r] :- Y[obs,r]):^2)
    return(err/(N*nY))
}

real scalar rfeat_comparePareto(real rowvector a1, real rowvector b1,
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

void rfeat_estimNode(
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
        if (rfeat_comparePareto(aL,bL,la,lb)==1) {
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
    tRes[1, cy..cy+nY-1]=yhL; tRes[1,7]=rfeat_mse(Y,obsL,yhL,N,nY)

    tRes[2,1]=.; tRes[2,2]=prow[1,1]; tRes[2,3..6]=J(1,4,-1)
    tRes[2, ca..ca+nX-1]=aR; tRes[2, cb..cb+nX-1]=bR
    tRes[2, cy..cy+nY-1]=yhR; tRes[2,7]=rfeat_mse(Y,obsR,yhR,N,nY)

    oL = obsL; oR = obsR
}

// -- RFEAT predictor: traverse nm to leaf, return yhat -----------------
real rowvector rfeat_predictor(real matrix nm, real rowvector xobs,
                               real scalar nX, real scalar nY)
{
    real scalar ti
    ti = 1
    while (nm[ti,3] != -1) {
        ti = (xobs[nm[ti,5]] < nm[ti,6]) ?
             rfeat_posNode(nm, nm[ti,3]) : rfeat_posNode(nm, nm[ti,4])
    }
    return(nm[ti, (8+2*nX)..(7+2*nX+nY)])
}

// -- RFEAT efficiency: same logic as eat_efficiency --------------------
real colvector rfeat_efficiency(real matrix X, real matrix Y,
                                real matrix nm,
                                real scalar nX, real scalar nY)
{
    real scalar  ca, cy, N, j, k, r, n_leaves
    real scalar  theta_j, theta_best, ratio
    real matrix  leaves
    real colvector leaf_rows, eff
    real rowvector xk, yk, a_j, yhat_j

    ca = 8; cy = 8+2*nX
    N = rows(X)
    leaf_rows = selectindex(nm[.,3]:==-1)
    leaves    = nm[leaf_rows,.]
    n_leaves  = rows(leaves)
    eff       = J(N,1,.)

    for (k=1; k<=N; k++) {
        xk = X[k,.]; yk = Y[k,.]
        theta_best = .
        for (j=1; j<=n_leaves; j++) {
            a_j    = leaves[j, ca..ca+nX-1]
            yhat_j = leaves[j, cy..cy+nY-1]
            if (!all(a_j :<= xk)) continue
            theta_j = 1e15
            for (r=1; r<=nY; r++) {
                if (yk[r]<=0) continue
                ratio = yhat_j[r]/yk[r]
                if (ratio < theta_j) theta_j = ratio
            }
            if (theta_j >= 1e14) theta_j = 1
            if (theta_best==. | theta_j>theta_best) theta_best = theta_j
        }
        if (theta_best == .) theta_best = 1   // no dominating leaf -> on frontier
        eff[k] = theta_best   // may be < 1 (super-efficient w.r.t. bootstrap tree)
    }
    return(eff)
}

// -- RFEAT tree grower (no pruning states; feature subsampling) ---------
// Returns the fully grown node matrix.
// s_feat: number of features to consider at each split (0 = all)
real matrix rfeat_grow_tree(
    real matrix X, real matrix Y,
    real scalar numStop, real scalar maxDepth, real scalar maxLeaves,
    real scalar N, real scalar nX, real scalar nY,
    real scalar s_feat
    )
{
    real scalar  K, ca, cb, cy
    real scalar  next_id, n_nodes, leaf_count, qlen
    real scalar  xi, xi_pos, i, err, err_min, best_xi, best_s
    real scalar  id_L, id_R, pidx, anc_id, d_cnt
    real scalar  iL, iR
    real matrix  nm, prow, tRes, oleaves, btRes
    real colvector pobs, obsL, obsR, boL, boR, queue, aobs, Sx, feat_idx
    transmorphic obs_db

    K = 7+2*nX+nY; ca = 8; cb = 8+nX; cy = 8+2*nX

    // initialise root node
    nm = J(1,K,.)
    nm[1,1]=1; nm[1,2]=-1; nm[1,3]=-1; nm[1,4]=-1
    nm[1,5]=-1; nm[1,6]=-1
    nm[1, ca..ca+nX-1] = colmin(X)
    nm[1, cb..cb+nX-1] = J(1,nX,.)
    nm[1, cy..cy+nY-1] = colmax(Y)

    obs_db = asarray_create("real")
    aobs   = range(1,N,1)
    asarray(obs_db, 1, aobs)
    nm[1,7] = rfeat_mse(Y, aobs, nm[1,cy..cy+nY-1], N, nY)

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

        if (rfeat_isFinal(X, pobs, numStop)) continue

        if (maxDepth > 0) {
            anc_id = prow[1,2]; d_cnt = 0
            while (anc_id != -1 & d_cnt < maxDepth+3) {
                anc_id = nm[rfeat_posNode(nm,anc_id), 2]; d_cnt++
            }
            if (d_cnt >= maxDepth) continue
        }

        oleaves = select(nm, (nm[.,3]:==-1) :& (nm[.,1]:!=prow[1,1]))

        err_min=1e300; best_xi=-1; best_s=.
        btRes=J(2,K,.); boL=J(0,1,.); boR=J(0,1,.)

        // -- feature subsampling -------------------------------------
        // range(1,nX,1) is nXx1 column; runiform(nX,1) is nXx1 column
        if (s_feat < nX) {
            // shuffle feature indices, take first s_feat
            feat_idx = sort((runiform(nX,1), range(1,nX,1)), 1)[1..s_feat, 2]
        } else {
            feat_idx = range(1,nX,1)
        }

        for (xi_pos=1; xi_pos<=rows(feat_idx); xi_pos++) {
            xi = feat_idx[xi_pos,1]
            Sx = sort(uniqrows(X[pobs,xi]), 1)
            if (length(Sx)<=1) continue
            for (i=2; i<=length(Sx); i++) {
                rfeat_estimNode(X,Y,oleaves,prow,xi,Sx[i],N,nX,nY,pobs,
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

        if (rfeat_isFinal(X,boR,numStop)) {
            queue=(iR\queue)
        } else {
            queue=(queue\iR)
        }
        if (rfeat_isFinal(X,boL,numStop)) {
            queue=(iL\queue)
        } else {
            queue=(queue\iL)
        }
        leaf_count++
    }
    return(nm)
}

// -- main RFEAT function ------------------------------------------------
struct RFEATResult scalar rfeat_main(
    real matrix X, real matrix Y,
    real scalar ntrees, real scalar numStop, real scalar s_feat,
    real scalar maxDepth, real scalar maxLeaves
    )
{
    struct RFEATResult scalar res
    real scalar  N, nX, nY, b, Nb, k, r, phi_k, ratio
    real colvector boot_idx, eff_b
    real matrix  Xb, Yb, nm_b, eff_mat, pred_b, pred_sum, pred_mean

    N=rows(X); nX=cols(X); nY=cols(Y)

    // accumulate efficiency scores: N x ntrees matrix
    eff_mat = J(N, ntrees, .)
    pred_sum = J(N, nY, 0)

    for (b=1; b<=ntrees; b++) {
        // -- bootstrap sample with replacement ---------------------
        boot_idx = ceil(runiform(N,1) :* N)  // N indices in 1..N
        Nb = N
        Xb = X[boot_idx,.]
        Yb = Y[boot_idx,.]

        // -- grow fully (no CV pruning) -----------------------------
        nm_b = rfeat_grow_tree(Xb, Yb, numStop, maxDepth, maxLeaves,
                               Nb, nX, nY, s_feat)

        // Predict every output component for the ORIGINAL observations.
        // Paper eq. (9) averages these vectors before radial scoring.
        pred_b = J(N, nY, .)
        for (k=1; k<=N; k++) {
            pred_b[k,.] = rfeat_predictor(nm_b, X[k,.], nX, nY)
        }
        pred_sum = pred_sum + pred_b

        // Retain per-tree scores only for the reported ensemble SD.
        eff_b = J(N,1,.)
        for (k=1; k<=N; k++) {
            phi_k = 1e15
            for (r=1; r<=nY; r++) if (Y[k,r]>0) {
                ratio = pred_b[k,r]/Y[k,r]
                if (ratio<phi_k) phi_k=ratio
            }
            if (phi_k==1e15) phi_k=.
            eff_b[k]=phi_k
        }
        eff_mat[., b] = eff_b
    }

    // Aggregate output predictions first, then apply Proposition 2.
    pred_mean = pred_sum :/ ntrees
    res.frontier = pred_mean
    res.eff_mean = J(N,1,.)
    for (k=1; k<=N; k++) {
        phi_k = 1e15
        for (r=1; r<=nY; r++) if (Y[k,r]>0) {
            ratio = pred_mean[k,r]/Y[k,r]
            if (ratio<phi_k) phi_k=ratio
        }
        if (phi_k==1e15) phi_k=.
        res.eff_mean[k]=phi_k
    }
    res.eff_sd   = J(N, 1, .)
    if (ntrees > 1) {
        real rowvector row_k
        for (k=1; k<=N; k++) {
            row_k = eff_mat[k,.]
            res.eff_sd[k] = rfeat_sample_sd(row_k)
        }
    }
    res.ntrees = ntrees
    res.s_feat = s_feat
    return(res)
}

real scalar rfeat_sample_sd(real rowvector v)
{
    if (cols(v) < 2) return(.)
    // mean() computes column means; for a row vector it returns v itself.
    return(sqrt(sum((v :- sum(v)/cols(v)):^2) / (cols(v)-1)))
}

// -- display -----------------------------------------------------------
void rfeat_display(struct RFEATResult scalar res,
                   string scalar xnames, string scalar ynames,
                   real scalar ntrees, real scalar s_feat)
{
    real scalar N, eff_min, eff_max, eff_mean, n_eff1

    N        = rows(res.eff_mean)
    eff_mean = mean(res.eff_mean)
    eff_min  = min(res.eff_mean)
    eff_max  = max(res.eff_mean)
    n_eff1   = sum(res.eff_mean :>= 1 - 1e-8)

    printf("{txt}\n")
    printf("{txt}Random Forest Efficiency Analysis Trees (RFEAT)\n")
    printf("{hline 55}\n")
    printf("{txt}  Inputs      : {res}%s{txt}\n", xnames)
    printf("{txt}  Outputs     : {res}%s{txt}\n", ynames)
    printf("{txt}  Observations: {res}%g{txt}\n", N)
    printf("{txt}  Trees (B)   : {res}%g{txt}\n", ntrees)
    printf("{txt}  Features/split (s): {res}%g{txt}\n", s_feat)
    printf("{hline 55}\n")
    printf("{txt}  Mean efficiency     : {res}%9.4f{txt}\n", eff_mean)
    printf("{txt}  Min  efficiency     : {res}%9.4f{txt}\n", eff_min)
    printf("{txt}  Max  efficiency     : {res}%9.4f{txt}\n", eff_max)
    printf("{txt}  Observations with phi>=1: {res}%g{txt}\n", n_eff1)
    printf("{hline 55}\n")
}

end
// ======================================================================
