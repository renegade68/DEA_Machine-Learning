*! dyneat.ado  v0.2.1
*! Dynamic Efficiency Analysis Trees (dynamic EAT)
*! Author: Choonjoo Lee, Korea National Defense University
*!         Nonsan-si, Republic of Korea; bloom.rampike@gmail.com
*! Method reference: Aparicio, Esteve, and Kapelko (2023), "Measuring dynamic
*!   inefficiency through machine learning techniques", Expert Systems with
*!   Applications 228:120417. DOI: 10.1016/j.eswa.2023.120417
*!
*! EAT step-function version of dynamic CEAT (no LP required; leaf enumeration).
*! Measures dynamic inefficiency via DDF:
*!   beta = max{b: x - b*gx achievable, y + b*gy feasible, z + b*gI feasible}
*!
*! One EAT tree T* is grown with:
*!   inputs  = variable inputs x
*!   outputs = [y (good outputs), z_h = I_h - delta_h*k_h (net investment)]
*! Then for each DMU, iterate over all leaves and compute the maximum feasible beta.
*!
*! For leaf t and DMU k:
*!   beta_x = min_j{ (x0_j - a_j_t) / gx0_j }  [gx0_j > 0]
*!   beta_y = min_r{ (yhat_t_r - y0_r) / gy0_r } [gy0_r > 0]
*!   beta_z = min_h{ (zhat_t_h - z0_h) / gI_h }  [gI_h > 0]
*!   beta_t = min(beta_x, beta_y, beta_z)
*! DMU beta = max over leaves with beta_t >= 0
*!
*! Default DDF directions (translation-invariant):
*!   gx_j0  = x_j0 - min_all(x_j)
*!   gy_r0  = max_all(y_r) - y_r0
*!   gI_h0  = 0.2 * (max_all(k_h) - k_h0)
*!
*! Syntax:
*!   dyneat, x(varlist) y(varlist) inv(varlist) kfix(varlist) dep(numlist)
*!           [gen(name) numstop(#) fold(#) seed(#)]

program define dyneat, eclass
    version 16.1
    syntax [if] [in] , ///
        X(varlist numeric)    ///
        Y(varlist numeric)    ///
        INV(varlist numeric)  ///
        KFix(varlist numeric) ///
        DEP(string)           ///
        [ GEN(name)           ///
          NUMStop(integer 5)  ///
          FOLD(integer 5)     ///
          SEED(integer 12345) ]

    if "`gen'" != "" confirm new variable `gen'

    local xvars  `x'
    local yvars  `y'
    local invars `inv'
    local kvars  `kfix'

    local nX  : word count `xvars'
    local nY  : word count `yvars'
    local nI  : word count `invars'
    local nK  : word count `kvars'

    if `numstop' < 0 {
        di as err "numstop() must be >= 0"
        exit 198
    }
    if `fold' == 1 | `fold' < 0 {
        di as err "fold() must be 0 or >= 2"
        exit 198
    }

    if `nI' != `nK' {
        di as err "inv() and kfix() must have the same number of variables"
        exit 198
    }

    // Parse depreciation rates
    local nDep : word count `dep'
    if `nDep' != `nK' {
        di as err "dep() must have exactly `nK' depreciation rate(s)"
        exit 198
    }
    tempname depmat
    matrix `depmat' = J(1, `nK', 0)
    forvalues h = 1/`nK' {
        local dh : word `h' of `dep'
        confirm number `dh'
        if `dh' < 0 | `dh' > 1 {
            di as err "Depreciation rate `dh' out of [0,1] range"
            exit 198
        }
        matrix `depmat'[1, `h'] = `dh'
    }

    tempvar touse
    marksample touse, novarlist
    markout `touse' `xvars' `yvars' `invars' `kvars'
    qui count if `touse'
    local N = r(N)
    if `N' < 4 {
        di as err "Too few observations."
        exit 2001
    }
    if `fold' > `N' {
        di as err "fold() may not exceed the estimation-sample size"
        exit 198
    }
    foreach v of varlist `xvars' `kvars' {
        quietly count if `touse' & `v' < 0
        if r(N) {
            di as err "variable and quasi-fixed inputs must be nonnegative"
            exit 459
        }
    }
    tempvar ypos
    generate byte `ypos' = 0
    foreach v of varlist `yvars' {
        quietly count if `touse' & `v' < 0
        if r(N) {
            di as err "good outputs must be nonnegative"
            exit 459
        }
        quietly replace `ypos' = 1 if `touse' & `v' > 0
    }
    quietly count if `touse' & !`ypos'
    if r(N) {
        di as err "each observation must have at least one strictly positive good output"
        exit 459
    }

    tempname Xm Ym Im Km
    mkmat `xvars'  if `touse', matrix(`Xm')
    mkmat `yvars'  if `touse', matrix(`Ym')
    mkmat `invars' if `touse', matrix(`Im')
    mkmat `kvars'  if `touse', matrix(`Km')

    if `fold' > 0 set seed `seed'

    mata: _dyneat_res = dyneat_main( ///
        st_matrix("`Xm'"), st_matrix("`Ym'"), ///
        st_matrix("`Im'"), st_matrix("`Km'"), ///
        st_matrix("`depmat'"),                ///
        `numstop', `fold')

    if "`gen'" != "" {
        quietly gen double `gen' = . if `touse'
        mata: st_store(., "`gen'", "`touse'", _dyneat_res.beta)
        label var `gen' "Dynamic EAT DDF (beta)"
    }

    mata: dyneat_display(_dyneat_res, "`xvars'", "`yvars'", "`invars'", "`kvars'", "`dep'")

    ereturn clear
    ereturn post, esample(`touse')
    ereturn scalar N          = `N'
    mata: st_numscalar("r(n_leaves)",  _dyneat_res.n_leaves)
    mata: st_numscalar("r(n_nodes)",   _dyneat_res.n_nodes)
    mata: st_numscalar("r(n_zero_direction)", _dyneat_res.n_zero_direction)
    mata: st_numscalar("r(mean_beta)", dyneat_valid_mean(_dyneat_res.beta))
    ereturn scalar n_leaves  = r(n_leaves)
    ereturn scalar n_nodes   = r(n_nodes)
    ereturn scalar n_zero_direction = r(n_zero_direction)
    ereturn scalar mean_beta = r(mean_beta)
    ereturn matrix dep_rates = `depmat'
    if `nK' == 1 ereturn scalar dep = `dep'
    ereturn local  x       "`xvars'"
    ereturn local  y       "`yvars'"
    ereturn local  inv     "`invars'"
    ereturn local  kfix    "`kvars'"
    ereturn local  cmd     "dyneat"
    ereturn local  cmdline `"dyneat `0'"'
end


mata:
mata set matastrict off

// ======================================================================
// EAT functions (prefixed dyneat_eat_)
// Identical algorithm to dynceat.ado, different prefix
// ======================================================================

real scalar dyneat_eat_posNode(real matrix nm, real scalar id)
{
    real colvector v
    v = selectindex(nm[.,1] :== id)
    return(length(v) ? v[1] : 0)
}

real scalar dyneat_eat_isFinal(real matrix Xa, real colvector obs, real scalar ns)
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

real scalar dyneat_eat_mse(real matrix Y, real colvector obs,
                             real rowvector yhat, real scalar N, real scalar nY)
{
    real scalar err, r
    if (length(obs)==0) return(.)
    err = 0
    for (r=1; r<=nY; r++) err = err + sum((yhat[r] :- Y[obs,r]):^2)
    return(round(err/(N*nY), .0001))
}

real scalar dyneat_eat_comparePareto(real rowvector a1, real rowvector b1,
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

void dyneat_eat_estimNode(
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

    // Net investment may be negative. Zero is not a technology bound.
    // The sample minimum is a harmless finite lower bound and translates
    // with the responses, preserving the dynamic model's origin invariance.
    yInf = colmin(Y)
    n_ol = rows(oleaves)
    for (j=1; j<=n_ol; j++) {
        la = oleaves[j, ca..ca+nX-1]
        lb = oleaves[j, cb..cb+nX-1]
        if (dyneat_eat_comparePareto(aL,bL,la,lb)==1) {
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
    tRes[1, cy..cy+nY-1]=yhL; tRes[1,7]=dyneat_eat_mse(Y,obsL,yhL,N,nY)

    tRes[2,1]=.; tRes[2,2]=prow[1,1]; tRes[2,3..6]=J(1,4,-1)
    tRes[2, ca..ca+nX-1]=aR; tRes[2, cb..cb+nX-1]=bR
    tRes[2, cy..cy+nY-1]=yhR; tRes[2,7]=dyneat_eat_mse(Y,obsR,yhR,N,nY)

    oL = obsL; oR = obsR
}

real rowvector dyneat_eat_rbranch(real matrix nm, real scalar id)
{
    real scalar pos
    real rowvector rL, rR
    pos = dyneat_eat_posNode(nm, id)
    if (pos==0 | nm[pos,3]==-1) return((nm[pos,7], 1))
    rL = dyneat_eat_rbranch(nm, nm[pos,3])
    rR = dyneat_eat_rbranch(nm, nm[pos,4])
    return((rL[1]+rR[1], rL[2]+rR[2]))
}

real scalar dyneat_eat_alpha(real matrix nm)
{
    real scalar amin, alpha_i, i
    real rowvector rb
    amin = .
    for (i=1; i<=rows(nm); i++) {
        if (nm[i,3]==-1) continue
        rb = dyneat_eat_rbranch(nm, nm[i,1])
        if (rb[2]<=1) continue
        alpha_i = (nm[i,7]-rb[1])/(rb[2]-1)
        if (alpha_i<0) alpha_i = 0
        if (amin==. | alpha_i<amin) amin = alpha_i
        if (amin==0) break
    }
    return(amin==. ? 1e300 : amin)
}

void dyneat_eat_grow(
    real matrix X, real matrix Y,
    real scalar numStop, real scalar maxDepth, real scalar maxLeaves,
    real scalar N, real scalar nX, real scalar nY,
    transmorphic tree_db
    )
{
    real scalar  K, ca, cb, cy
    real scalar  next_id, n_nodes, leaf_count, qlen
    real scalar  xi, i, err, err_min, best_xi, best_s
    real scalar  id_L, id_R, pidx, anc_id, d_cnt
    real scalar  iL, iR, ns
    real matrix  nm, prow, tRes, oleaves, btRes
    real colvector pobs, obsL, obsR, boL, boR, queue, aobs, Sx
    real vector  av
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
    nm[1,7] = dyneat_eat_mse(Y, aobs, nm[1,cy..cy+nY-1], N, nY)

    n_nodes=1; next_id=2; leaf_count=1
    ns = 1
    av = J(1, 2*N+4, .)
    asarray(tree_db, 1, nm)
    av[1] = dyneat_eat_alpha(nm)
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

        if (dyneat_eat_isFinal(X, pobs, numStop)) continue

        if (maxDepth > 0) {
            anc_id = prow[1,2]; d_cnt = 0
            while (anc_id != -1 & d_cnt < maxDepth+3) {
                anc_id = nm[dyneat_eat_posNode(nm,anc_id), 2]; d_cnt++
            }
            if (d_cnt >= maxDepth) continue
        }

        oleaves = select(nm, (nm[.,3]:==-1) :& (nm[.,1]:!=prow[1,1]))

        err_min=1e300; best_xi=-1; best_s=.
        btRes=J(2,K,.); boL=J(0,1,.); boR=J(0,1,.)

        for (xi=1; xi<=nX; xi++) {
            Sx = sort(uniqrows(X[pobs,xi]), 1)
            if (length(Sx)<=1) continue
            for (i=2; i<=length(Sx); i++) {
                dyneat_eat_estimNode(X,Y,oleaves,prow,xi,Sx[i],N,nX,nY,pobs,
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

        if (dyneat_eat_isFinal(X,boR,numStop)) {
            queue=(iR\queue)
        } else {
            queue=(queue\iR)
        }
        if (dyneat_eat_isFinal(X,boL,numStop)) {
            queue=(iL\queue)
        } else {
            queue=(queue\iL)
        }
        leaf_count++

        ns++
        asarray(tree_db, ns, nm)
        av[ns] = dyneat_eat_alpha(nm)
    }

    asarray(tree_db, 0, J(1,1,ns))
    asarray(tree_db, -1, av[1..ns])
}

real matrix dyneat_eat_pickTree(transmorphic db, real vector av,
                                 real scalar ns, real scalar amid)
{
    real scalar best, j
    best = ns
    for (j=ns; j>=2; j--) {
        if (av[j]!=. & av[j] <= amid) best = j
    }
    return(asarray(db, best))
}

real rowvector dyneat_eat_predictor(real matrix nm, real rowvector xobs,
                                     real scalar nX, real scalar nY)
{
    real scalar ti, next_ti
    ti = 1
    while (nm[ti,3] != -1) {
        if (xobs[nm[ti,5]] < nm[ti,6]) {
            next_ti = dyneat_eat_posNode(nm, nm[ti,3])
        } else {
            next_ti = dyneat_eat_posNode(nm, nm[ti,4])
        }
        ti = next_ti
    }
    return(nm[ti, (8+2*nX)..(7+2*nX+nY)])
}

real matrix dyneat_eat_cv_select(
    real matrix X, real matrix Y,
    real scalar numStop, real scalar fold,
    real scalar N, real scalar nX, real scalar nY,
    transmorphic tree_db_full,
    real scalar  ns_full
    )
{
    real scalar  v, i, r, reg, rr
    real scalar  a1, a2, amid, Rcv
    real scalar  a1b, a2b, amid2, s2, SE, Tk_score
    real scalar  min_cv, best_state, final_state, margin
    real scalar  ns_v, dif1, nfold
    real vector  fold_id, cv_scores, av_v, perm, ixbest
    real colvector tr, val
    real matrix  nm_v, Xt, Yt, Xv, Yv, nm_b
    real rowvector pred
    transmorphic fold_dbs, db_v
    real vector av_full

    if (ns_full<=1) return(asarray(tree_db_full,1))

    av_full = asarray(tree_db_full, -1)

    fold_id = J(N,1,0)
    perm = order(runiform(N,1),1)
    nfold = floor(N/fold)
    for (v=1; v<=fold; v++) {
        fold_id[perm[((v-1)*nfold+1)..(v*nfold)]] = J(nfold,1,v)
    }

    fold_dbs = asarray_create("real")

    for (v=1; v<=fold; v++) {
        tr = selectindex(fold_id :!= v)
        Xt = X[tr,.]; Yt = Y[tr,.]
        db_v = asarray_create("real")
        dyneat_eat_grow(Xt,Yt,numStop,0,0,rows(Xt),nX,nY, db_v)
        ns_v = asarray(db_v, 0)[1,1]
        av_v = asarray(db_v, -1)
        asarray(fold_dbs, v*3-2, db_v)
        asarray(fold_dbs, v*3-1, av_v)
        asarray(fold_dbs, v*3,   J(1,1,ns_v))
    }

    cv_scores = J(1, ns_full, 1e15)
    for (i=2; i<=ns_full; i++) {
        a1 = av_full[i]
        a2 = av_full[i-1]
        if (a1==.) a1=0
        if (a2==.) a2=0
        amid = sqrt(a1*a2); Rcv = 0
        for (v=1; v<=fold; v++) {
            val  = selectindex(fold_id :== v)
            db_v = asarray(fold_dbs, v*3-2)
            av_v = asarray(fold_dbs, v*3-1)
            ns_v = asarray(fold_dbs, v*3)[1,1]
            nm_v = dyneat_eat_pickTree(db_v, av_v, ns_v, amid)
            Xv   = X[val,.]; Yv = Y[val,.]
            for (reg=1; reg<=rows(Xv); reg++) {
                pred = dyneat_eat_predictor(nm_v, Xv[reg,.], nX, nY)
                for (r=1; r<=nY; r++) Rcv = Rcv+(Yv[reg,r]-pred[r])^2
            }
        }
        cv_scores[i] = Rcv/(N*nY)
    }

    min_cv     = min(cv_scores)
    ixbest      = selectindex(cv_scores :== min_cv)
    best_state  = max(ixbest)
    Tk_score   = cv_scores[best_state]

    a1b  = av_full[best_state]
    a2b  = av_full[best_state-1]
    if (a1b==.) a1b=0
    if (a2b==.) a2b=0
    amid2 = sqrt(a1b*a2b)

    s2 = 0
    for (v=1; v<=fold; v++) {
        val  = selectindex(fold_id :== v)
        db_v = asarray(fold_dbs, v*3-2)
        av_v = asarray(fold_dbs, v*3-1)
        ns_v = asarray(fold_dbs, v*3)[1,1]
        nm_b = dyneat_eat_pickTree(db_v, av_v, ns_v, amid2)
        Xv   = X[val,.]; Yv = Y[val,.]
        for (reg=1; reg<=rows(Xv); reg++) {
            pred = dyneat_eat_predictor(nm_b, Xv[reg,.], nX, nY)
            for (rr=1; rr<=nY; rr++) {
                dif1 = (Yv[reg,rr]-pred[rr])^2
                s2   = s2+(dif1-Tk_score)^2
            }
        }
    }
    SE = sqrt(s2/(N*nY)/(N*nY))

    final_state = best_state
    margin = Tk_score + SE
    for (i=2; i<best_state; i++) {
        if (cv_scores[i]!=1e15 & cv_scores[i]<=margin) {
            final_state = i
            break
        }
    }

    return(asarray(tree_db_full, final_state))
}

// ======================================================================
// Dynamic EAT specific Mata code
// ======================================================================

struct DynEATResult {
    real matrix    node_mat
    real colvector beta
    real scalar    n_leaves
    real scalar    n_nodes
    real scalar    n_zero_direction
}

// -----------------------------------------------------------------------
// Compute directional vectors (same as dynceat)
// -----------------------------------------------------------------------
void dyneat_directions(
    real matrix X,
    real matrix Y,
    real matrix K,
    real matrix GX,
    real matrix GY,
    real matrix GI
    )
{
    real scalar nX, nY, nZ
    real rowvector minX, maxY, maxK

    nX = cols(X); nY = cols(Y); nZ = cols(K)
    minX = colmin(X)
    maxY = colmax(Y)
    maxK = colmax(K)

    GX = X :- minX
    GY = maxY :- Y
    GI = 0.2 :* (maxK :- K)
}

// -----------------------------------------------------------------------
// Leaf enumeration for one DMU
// For each leaf t:
//   beta_x = min_j{ (x0_j - a_j_t) / gx0_j }  (gx0_j > 0)
//   beta_y = min_r{ (yhat_t_r - y0_r) / gy0_r } (gy0_r > 0)
//   beta_z = min_h{ (zhat_t_h - z0_h) / gI_h }  (gI_h > 0)
//   beta_t = min(beta_x, beta_y, beta_z)
// DMU beta = max{beta_t : beta_t >= 0}
// -----------------------------------------------------------------------
real scalar dyneat_beta_leaf(
    real matrix    lv_a,    // nJ x nX
    real matrix    lv_yhat, // nJ x nY
    real matrix    lv_zhat, // nJ x nZ
    real colvector x0,      // nX x 1
    real colvector y0,      // nY x 1
    real colvector z0,      // nZ x 1
    real colvector gx0,     // nX x 1
    real colvector gy0,     // nY x 1
    real colvector gI,      // nZ x 1
    real scalar nX, real scalar nY, real scalar nZ, real scalar nJ
    )
{
    real scalar t, j, r, h, beta_t, beta_best, bval
    real scalar infeasible

    beta_best = .

    for (t=1; t<=nJ; t++) {
        beta_t  = 1e15
        infeasible = 0

        // Input constraints: a_t + beta*gx0 <= x0
        for (j=1; j<=nX; j++) {
            if (gx0[j] > 1e-12) {
                bval = (x0[j] - lv_a[t,j]) / gx0[j]
                if (bval < beta_t) beta_t = bval
            } else {
                // gx0_j = 0: constraint is a_t_j <= x0_j (no beta limit)
                if (lv_a[t,j] > x0[j] + 1e-10) {
                    infeasible = 1
                    break
                }
            }
        }
        if (infeasible) continue

        // Output constraints: yhat_t - beta*gy0 >= y0  ->  beta <= (yhat - y0)/gy0
        for (r=1; r<=nY; r++) {
            if (gy0[r] > 1e-12) {
                bval = (lv_yhat[t,r] - y0[r]) / gy0[r]
                if (bval < beta_t) beta_t = bval
            } else {
                // gy0_r = 0: constraint is yhat_t_r >= y0_r
                if (lv_yhat[t,r] < y0[r] - 1e-10) {
                    infeasible = 1
                    break
                }
            }
        }
        if (infeasible) continue

        // Investment constraints: zhat_t - beta*gI >= z0  ->  beta <= (zhat - z0)/gI
        for (h=1; h<=nZ; h++) {
            if (gI[h] > 1e-12) {
                bval = (lv_zhat[t,h] - z0[h]) / gI[h]
                if (bval < beta_t) beta_t = bval
            } else {
                // gI_h = 0: constraint is zhat_t_h >= z0_h
                if (lv_zhat[t,h] < z0[h] - 1e-10) {
                    infeasible = 1
                    break
                }
            }
        }
        if (infeasible) continue

        // Leaf t is feasible if beta_t >= 0
        if (beta_t >= -1e-10) {
            if (beta_t < 0) beta_t = 0
            if (beta_best==. | beta_t > beta_best) beta_best = beta_t
        }
    }

    return(beta_best==. ? 0 : beta_best)
}

// -----------------------------------------------------------------------
// Main dynamic EAT function
// -----------------------------------------------------------------------
struct DynEATResult scalar dyneat_main(
    real matrix X,
    real matrix Y,
    real matrix Inv,
    real matrix K,
    real rowvector delta,
    real scalar numStop,
    real scalar fold
    )
{
    struct DynEATResult scalar res
    real scalar N, nX, nY, nZ, nYZ, nJ, k, n_zero_direction
    real scalar ca, cy, ns
    real matrix Z_comp, Y_combined, final_tree, leaves
    real matrix lv_a, lv_yhat, lv_zhat
    real matrix GX, GY, GI
    real colvector x0, y0, z0, gx0, gy0, gI0, beta
    transmorphic tree_db

    N  = rows(X); nX = cols(X); nY = cols(Y); nZ = cols(K)
    nYZ = nY + nZ

    // Aparicio, Esteve, and Kapelko (2023), model (11): I_h - delta_h*k_h.
    Z_comp     = Inv - (K :* delta)
    Y_combined = (Y, Z_comp)

    ca = 8; cy = 8+2*nX

    tree_db = asarray_create("real")
    dyneat_eat_grow(X, Y_combined, numStop, 0, 0, N, nX, nYZ, tree_db)
    ns = asarray(tree_db, 0)[1,1]

    if (ns<=1 | fold<=0) {
        final_tree = asarray(tree_db, ns)
    } else {
        final_tree = dyneat_eat_cv_select(X, Y_combined, numStop, fold, N, nX, nYZ,
                                          tree_db, ns)
    }

    leaves   = select(final_tree, final_tree[.,3]:==-1)
    nJ       = rows(leaves)
    lv_a     = leaves[., ca..ca+nX-1]
    lv_yhat  = leaves[., cy..cy+nY-1]
    lv_zhat  = leaves[., cy+nY..cy+nY+nZ-1]

    dyneat_directions(X, Y, K, GX, GY, GI)

    beta = J(N, 1, .)
    n_zero_direction = 0

    for (k=1; k<=N; k++) {
        x0  = X[k,.]'
        y0  = Y[k,.]'
        z0  = Z_comp[k,.]'
        gx0 = GX[k,.]'
        gy0 = GY[k,.]'
        gI0 = GI[k,.]'

        // The DDF is undefined when every component of the direction is zero.
        if (max(abs(gx0 \ gy0 \ gI0)) <= 1e-12) {
            n_zero_direction++
            continue
        }

        beta[k] = dyneat_beta_leaf(lv_a, lv_yhat, lv_zhat,
                                    x0, y0, z0, gx0, gy0, gI0,
                                    nX, nY, nZ, nJ)
    }

    res.node_mat = final_tree
    res.beta     = beta
    res.n_leaves = nJ
    res.n_nodes  = rows(final_tree)
    res.n_zero_direction = n_zero_direction
    return(res)
}

real scalar dyneat_valid_mean(real colvector v)
{
    real colvector w
    w = select(v, v :!= .)
    return(length(w) > 0 ? mean(w) : .)
}

void dyneat_display(struct DynEATResult scalar res,
                    string scalar xnames, string scalar ynames,
                    string scalar invnames, string scalar knames,
                    string scalar dep)
{
    real scalar N
    real colvector v, eff0
    N = rows(res.beta)
    v = select(res.beta, res.beta :!= .)

    printf("{txt}\n")
    printf("{txt}Dynamic EAT (Aparicio, Esteve, Kapelko 2023, ESWA)\n")
    printf("{hline 60}\n")
    printf("{txt}  Variable inputs (x)     : {res}%s{txt}\n", xnames)
    printf("{txt}  Good outputs (y)        : {res}%s{txt}\n", ynames)
    printf("{txt}  Gross investments (I)   : {res}%s{txt}\n", invnames)
    printf("{txt}  Quasi-fixed inputs (k)  : {res}%s{txt}\n", knames)
    printf("{txt}  Depreciation rates (d)  : {res}%s{txt}\n", dep)
    printf("{txt}  Observations            : {res}%g{txt}\n", N)
    printf("{txt}  EAT leaves              : {res}%g{txt}  (nodes: %g)\n",
           res.n_leaves, res.n_nodes)
    if (res.n_zero_direction > 0) {
        printf("{err}  Zero direction vectors : %g (scores set to missing){txt}\n",
               res.n_zero_direction)
    }
    printf("{hline 60}\n")
    if (length(v) > 0) {
        eff0 = select(res.beta, abs(res.beta) :<= 1e-6)
        printf("{txt}  Dynamic DDF (beta, TI direction):\n")
        printf("{txt}    Mean= %7.4f  Min= %7.4f  Max= %7.4f\n",
               mean(v), min(v), max(v))
        printf("{txt}    Efficient DMUs (beta~=0): %g\n", length(eff0))
    }
    printf("{hline 60}\n")
}

end
