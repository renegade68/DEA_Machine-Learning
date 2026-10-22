*! bpeat.ado  v0.2.1
*! By-production Efficiency Analysis Trees (by-production EAT)
*! Author: Choonjoo Lee, Korea National Defense University
*!         Nonsan-si, Republic of Korea; bloom.rampike@gmail.com
*! Method reference: Guillen, Aparicio, Kapelko, and Esteve (2025),
*!   "Measuring environmental inefficiency through machine learning: An
*!   approach based on efficiency analysis trees and by-production technology",
*!   EJOR 321(2):529-542. DOI: 10.1016/j.ejor.2024.10.003
*!
*! Measures environmental inefficiency by leaf enumeration (no LP required):
*!   - T1*: good output technology  (inputs = x1+x2, outputs = y)
*!   - T2*: bad output technology   (inputs = z, outputs = x2)
*!   - beta_Psi1: max DDF over leaves dominating DMU inputs
*!   - beta_Psi2: max DDF over leaves with predicted p >= actual p
*!   - Combined: B = delta1*beta_Psi1/(1+beta_Psi1) + delta2*beta_Psi2
*!
*! Syntax:
*!   bpeat, x1(varlist) x2(varlist) y(varlist) z(varlist)
*!          [gen(name) genpsi1(name) genpsi2(name)
*!           delta1(#) delta2(#) numstop(#) fold(#) seed(#)]

program define bpeat, eclass
    version 16.1
    syntax [if] [in] , ///
        X1(varlist numeric) ///
        X2(varlist numeric) ///
        Y(varlist numeric) ///
        Z(varlist numeric) ///
        [ GEN(name)           ///
          GENpsi1(name)        ///
          GENpsi2(name)        ///
          DELta1(real 0.5)     ///
          DELta2(real 0.5)     ///
          NUMStop(integer 5)   ///
          FOLD(integer 5)      ///
          SEED(integer 12345)  ]

    foreach out in `gen' `genpsi1' `genpsi2' {
        confirm new variable `out'
    }
    if ("`gen'" != "" & "`gen'" == "`genpsi1'") | ///
       ("`gen'" != "" & "`gen'" == "`genpsi2'") | ///
       ("`genpsi1'" != "" & "`genpsi1'" == "`genpsi2'") {
        di as err "gen(), genpsi1(), and genpsi2() must name distinct new variables"
        exit 198
    }

    if `delta1' < 0 | `delta2' < 0 | abs(`delta1' + `delta2' - 1) > 1e-10 {
        di as err "delta1() and delta2() must be nonnegative and sum to 1"
        exit 198
    }
    if `numstop' < 0 {
        di as err "numstop() must be >= 0"
        exit 198
    }
    if `fold' == 1 | `fold' < 0 {
        di as err "fold() must be 0 or >= 2"
        exit 198
    }

    tempvar touse
    marksample touse, novarlist
    markout `touse' `x1' `x2' `y' `z'
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
    foreach v of varlist `x1' `x2' `y' `z' {
        quietly count if `touse' & `v' < 0
        if r(N) {
            di as err "inputs and good/bad outputs must be nonnegative"
            exit 459
        }
    }
    tempvar ypos zpos
    generate byte `ypos' = 0
    generate byte `zpos' = 0
    foreach v of varlist `y' {
        quietly replace `ypos' = 1 if `touse' & `v' > 0
    }
    foreach v of varlist `z' {
        quietly replace `zpos' = 1 if `touse' & `v' > 0
    }
    quietly count if `touse' & !`ypos'
    if r(N) {
        di as err "each observation must have at least one strictly positive good output"
        exit 459
    }
    quietly count if `touse' & !`zpos'
    if r(N) {
        di as err "each observation must have at least one strictly positive bad output"
        exit 459
    }

    tempname X1m X2m Ym Zm
    mkmat `x1' if `touse', matrix(`X1m')
    mkmat `x2' if `touse', matrix(`X2m')
    mkmat `y'  if `touse', matrix(`Ym')
    mkmat `z'  if `touse', matrix(`Zm')

    if `fold' > 0 set seed `seed'

    mata: _bpeat_res = bpeat_main( ///
        st_matrix("`X1m'"), st_matrix("`X2m'"), ///
        st_matrix("`Ym'"),  st_matrix("`Zm'"),  ///
        `numstop', `fold', `delta1', `delta2')

    if "`gen'" != "" {
        quietly gen double `gen' = . if `touse'
        mata: st_store(., "`gen'", "`touse'", _bpeat_res.beta_comb)
        label var `gen' "BP-EAT combined environmental inefficiency"
    }
    if "`genpsi1'" != "" {
        quietly gen double `genpsi1' = . if `touse'
        mata: st_store(., "`genpsi1'", "`touse'", _bpeat_res.beta_psi1)
        label var `genpsi1' "BP-EAT Psi1 DDF (good output)"
    }
    if "`genpsi2'" != "" {
        quietly gen double `genpsi2' = . if `touse'
        mata: st_store(., "`genpsi2'", "`touse'", _bpeat_res.beta_psi2)
        label var `genpsi2' "BP-EAT Psi2 DDF (bad output)"
    }

    mata: bpeat_display(_bpeat_res, "`x1'", "`x2'", "`y'", "`z'", `delta1', `delta2')

    ereturn clear
    ereturn post, esample(`touse')
    ereturn scalar N         = `N'
    mata: st_numscalar("r(nJ1)", _bpeat_res.nJ1)
    mata: st_numscalar("r(nJ2)", _bpeat_res.nJ2)
    mata: st_numscalar("r(mean_beta1)", bpeat_valid_mean(_bpeat_res.beta_psi1))
    mata: st_numscalar("r(mean_beta2)", bpeat_valid_mean(_bpeat_res.beta_psi2))
    mata: st_numscalar("r(mean_b)",  bpeat_valid_mean(_bpeat_res.beta_comb))
    ereturn scalar n_leaves1  = r(nJ1)
    ereturn scalar n_leaves2  = r(nJ2)
    ereturn scalar mean_beta1 = r(mean_beta1)
    ereturn scalar mean_beta2 = r(mean_beta2)
    ereturn scalar mean_b  = r(mean_b)
    ereturn scalar delta1     = `delta1'
    ereturn scalar delta2     = `delta2'
    ereturn local  x1      "`x1'"
    ereturn local  x2      "`x2'"
    ereturn local  y       "`y'"
    ereturn local  z       "`z'"
    ereturn local  cmd     "bpeat"
    ereturn local  cmdline `"bpeat `0'"'
end


mata:
mata set matastrict off

// ======================================================================
// EAT Mata functions (embedded; prefixed bpeat_eat_)
// Same algorithm as eatree.ado / ceat.ado
// ======================================================================

real scalar bpeat_eat_posNode(real matrix nm, real scalar id)
{
    real colvector v
    v = selectindex(nm[.,1] :== id)
    return(length(v) ? v[1] : 0)
}

real scalar bpeat_eat_isFinal(real matrix Xa, real colvector obs, real scalar ns)
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

real scalar bpeat_eat_mse(real matrix Y, real colvector obs,
                           real rowvector yhat, real scalar N, real scalar nY)
{
    real scalar err, r
    if (length(obs)==0) return(.)
    err = 0
    for (r=1; r<=nY; r++) err = err + sum((yhat[r] :- Y[obs,r]):^2)
    return(round(err/(N*nY), .0001))
}

real scalar bpeat_eat_comparePareto(real rowvector a1, real rowvector b1,
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

void bpeat_eat_estimNode(
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
        if (bpeat_eat_comparePareto(aL,bL,la,lb)==1) {
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
    tRes[1, cy..cy+nY-1]=yhL; tRes[1,7]=bpeat_eat_mse(Y,obsL,yhL,N,nY)

    tRes[2,1]=.; tRes[2,2]=prow[1,1]; tRes[2,3..6]=J(1,4,-1)
    tRes[2, ca..ca+nX-1]=aR; tRes[2, cb..cb+nX-1]=bR
    tRes[2, cy..cy+nY-1]=yhR; tRes[2,7]=bpeat_eat_mse(Y,obsR,yhR,N,nY)

    oL = obsL; oR = obsR
}

real rowvector bpeat_eat_rbranch(real matrix nm, real scalar id)
{
    real scalar pos
    real rowvector rL, rR
    pos = bpeat_eat_posNode(nm, id)
    if (pos==0 | nm[pos,3]==-1) return((nm[pos,7], 1))
    rL = bpeat_eat_rbranch(nm, nm[pos,3])
    rR = bpeat_eat_rbranch(nm, nm[pos,4])
    return((rL[1]+rR[1], rL[2]+rR[2]))
}

real scalar bpeat_eat_alpha(real matrix nm)
{
    real scalar amin, alpha_i, i
    real rowvector rb
    amin = .
    for (i=1; i<=rows(nm); i++) {
        if (nm[i,3]==-1) continue
        rb = bpeat_eat_rbranch(nm, nm[i,1])
        if (rb[2]<=1) continue
        alpha_i = (nm[i,7]-rb[1])/(rb[2]-1)
        if (alpha_i<0) alpha_i = 0
        if (amin==. | alpha_i<amin) amin = alpha_i
        if (amin==0) break
    }
    return(amin==. ? 1e300 : amin)
}

void bpeat_eat_grow(
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
    nm[1,7] = bpeat_eat_mse(Y, aobs, nm[1,cy..cy+nY-1], N, nY)

    n_nodes=1; next_id=2; leaf_count=1
    ns = 1
    av = J(1, 2*N+4, .)
    asarray(tree_db, 1, nm)
    av[1] = bpeat_eat_alpha(nm)
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

        if (bpeat_eat_isFinal(X, pobs, numStop)) continue

        if (maxDepth > 0) {
            anc_id = prow[1,2]; d_cnt = 0
            while (anc_id != -1 & d_cnt < maxDepth+3) {
                anc_id = nm[bpeat_eat_posNode(nm,anc_id), 2]; d_cnt++
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
                bpeat_eat_estimNode(X,Y,oleaves,prow,xi,Sx[i],N,nX,nY,pobs,
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

        if (bpeat_eat_isFinal(X,boR,numStop)) {
            queue=(iR\queue)
        } else {
            queue=(queue\iR)
        }
        if (bpeat_eat_isFinal(X,boL,numStop)) {
            queue=(iL\queue)
        } else {
            queue=(queue\iL)
        }
        leaf_count++

        ns++
        asarray(tree_db, ns, nm)
        av[ns] = bpeat_eat_alpha(nm)
    }

    asarray(tree_db, 0, J(1,1,ns))
    asarray(tree_db, -1, av[1..ns])
}

real matrix bpeat_eat_pickTree(transmorphic db, real vector av,
                                real scalar ns, real scalar amid)
{
    real scalar best, j
    best = ns
    for (j=ns; j>=2; j--) {
        if (av[j]!=. & av[j] <= amid) best = j
    }
    return(asarray(db, best))
}

real rowvector bpeat_eat_predictor(real matrix nm, real rowvector xobs,
                                    real scalar nX, real scalar nY)
{
    real scalar ti, next_ti
    ti = 1
    while (nm[ti,3] != -1) {
        if (xobs[nm[ti,5]] < nm[ti,6]) {
            next_ti = bpeat_eat_posNode(nm, nm[ti,3])
        } else {
            next_ti = bpeat_eat_posNode(nm, nm[ti,4])
        }
        ti = next_ti
    }
    return(nm[ti, (8+2*nX)..(7+2*nX+nY)])
}

real matrix bpeat_eat_cv_select(
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
        bpeat_eat_grow(Xt,Yt,numStop,0,0,rows(Xt),nX,nY, db_v)
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
            nm_v = bpeat_eat_pickTree(db_v, av_v, ns_v, amid)
            Xv   = X[val,.]; Yv = Y[val,.]
            for (reg=1; reg<=rows(Xv); reg++) {
                pred = bpeat_eat_predictor(nm_v, Xv[reg,.], nX, nY)
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
        nm_b = bpeat_eat_pickTree(db_v, av_v, ns_v, amid2)
        Xv   = X[val,.]; Yv = Y[val,.]
        for (reg=1; reg<=rows(Xv); reg++) {
            pred = bpeat_eat_predictor(nm_b, Xv[reg,.], nX, nY)
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
// By-production EAT specific Mata code (leaf enumeration, no LP)
// ======================================================================

struct BPEATResult {
    real matrix    nm1
    real matrix    nm2
    real colvector beta_psi1
    real colvector beta_psi2
    real colvector beta_comb
    real scalar    nJ1
    real scalar    nJ2
    real scalar    n_nodes1
    real scalar    n_nodes2
}

// -----------------------------------------------------------------------
// Compute beta_Psi1 for one DMU by leaf enumeration (MILP equivalent)
//
// Select ONE leaf t from T1* that:
//  - dominates x120 in all input dimensions: lv_a[t,j] <= x120[j] for all j
//  - maximizes min_r( lv_yhat[t,r] / y0[r] ) - 1
//
// Returns the max DDF value (>= 0), or 0 if no dominating leaf found.
// -----------------------------------------------------------------------
real scalar bpeat_beta_psi1(
    real matrix lv_a,     // nJ1 x nX12
    real matrix lv_yhat,  // nJ1 x nY
    real rowvector x120,  // 1 x nX12
    real rowvector y0,    // 1 x nY
    real scalar nX12, real scalar nY, real scalar nJ1
    )
{
    real scalar t, j, r, best_beta, beta_t, ratio_r
    real scalar feasible

    best_beta = 0
    for (t=1; t<=nJ1; t++) {
        // Check input domination
        feasible = 1
        for (j=1; j<=nX12; j++) {
            if (lv_a[t,j] > x120[j] + 1e-10) {
                feasible = 0
                break
            }
        }
        if (!feasible) continue

        // Compute DDF = min_r(yhat_t,r / y0_r) - 1
        beta_t = .
        for (r=1; r<=nY; r++) {
            if (y0[r] <= 0) continue
            ratio_r = lv_yhat[t,r] / y0[r] - 1
            if (beta_t==. | ratio_r < beta_t) beta_t = ratio_r
        }
        if (beta_t == .) beta_t = 0
        if (beta_t < 0) beta_t = 0

        if (beta_t > best_beta) best_beta = beta_t
    }
    return(best_beta)
}

// -----------------------------------------------------------------------
// Compute beta_Psi2 for one DMU by leaf enumeration (MILP equivalent)
//
// Select ONE leaf t from T2* that:
//  - has sufficient predicted polluting input: lv_phat[t,j] >= x20[j] all j
//  - maximizes min_j( 1 - lv_az[t,j] / z0[j] ) over dims with z0[j]>0
//
// Returns the max beta (>= 0), or 0 if no feasible leaf found.
// -----------------------------------------------------------------------
real scalar bpeat_beta_psi2(
    real matrix lv_az,    // nJ2 x nZ
    real matrix lv_phat,  // nJ2 x nX2
    real rowvector x20,   // 1 x nX2
    real rowvector z0,    // 1 x nZ
    real scalar nX2, real scalar nZ, real scalar nJ2
    )
{
    real scalar t, j, best_beta, beta_t, feasible

    best_beta = 0
    for (t=1; t<=nJ2; t++) {
        // Check: phat[t,j] >= x20[j] for all j (polluting input constraint)
        feasible = 1
        for (j=1; j<=nX2; j++) {
            if (lv_phat[t,j] < x20[j] - 1e-10) {
                feasible = 0
                break
            }
        }
        if (!feasible) continue

        // Compute beta_t = min_j (1 - lv_az[t,j]/z0[j]) for z0[j]>0
        // = min reduction in z such that leaf lower bound still satisfied
        beta_t = .
        for (j=1; j<=nZ; j++) {
            if (z0[j] <= 1e-12) {
                // Zero direction still imposes az_j <= z0_j. A positive
                // support cannot produce a zero observed bad output.
                if (lv_az[t,j] > z0[j] + 1e-10) feasible = 0
                continue
            }
            real scalar bj
            bj = 1 - lv_az[t,j] / z0[j]
            if (beta_t==. | bj < beta_t) beta_t = bj
        }
        if (!feasible) continue
        if (beta_t == .) beta_t = 1     // no z0>0: unconstrained -> cap at 1
        if (beta_t < 0) beta_t = 0
        if (beta_t > 1) beta_t = 1

        if (beta_t > best_beta) best_beta = beta_t
    }
    return(best_beta)
}

// -----------------------------------------------------------------------
// Main by-production EAT function
// -----------------------------------------------------------------------
struct BPEATResult scalar bpeat_main(
    real matrix X1,
    real matrix X2,
    real matrix Y,
    real matrix Z,
    real scalar numStop,
    real scalar fold,
    real scalar delta1,
    real scalar delta2
    )
{
    struct BPEATResult scalar res
    real scalar N, nX1, nX2, nY, nZ, nX12, nJ1, nJ2, k
    real scalar ca1, cy1, ca2, cy2, ns1, ns2
    real matrix X12, leaves1, leaves2, lv_a, lv_yhat, lv_az, lv_phat
    real colvector b1, b2, bc
    transmorphic tree_db1, tree_db2

    N    = rows(X1)
    nX1  = cols(X1)
    nX2  = cols(X2)
    nY   = cols(Y)
    nZ   = cols(Z)
    nX12 = nX1 + nX2

    X12 = (X1, X2)
    ca1 = 8; cy1 = 8+2*nX12
    ca2 = 8; cy2 = 8+2*nZ

    // Grow T1*: inputs=(X1,X2), outputs=Y
    tree_db1 = asarray_create("real")
    bpeat_eat_grow(X12, Y, numStop, 0, 0, N, nX12, nY, tree_db1)
    ns1 = asarray(tree_db1, 0)[1,1]

    if (ns1<=1 | fold<=0) {
        res.nm1 = asarray(tree_db1, ns1)
    } else {
        res.nm1 = bpeat_eat_cv_select(X12, Y, numStop, fold, N, nX12, nY,
                                       tree_db1, ns1)
    }

    leaves1  = select(res.nm1, res.nm1[.,3]:==-1)
    nJ1      = rows(leaves1)
    lv_a     = leaves1[., ca1..ca1+nX12-1]
    lv_yhat  = leaves1[., cy1..cy1+nY-1]

    // Grow T2*: inputs=Z, outputs=X2
    tree_db2 = asarray_create("real")
    bpeat_eat_grow(Z, X2, numStop, 0, 0, N, nZ, nX2, tree_db2)
    ns2 = asarray(tree_db2, 0)[1,1]

    if (ns2<=1 | fold<=0) {
        res.nm2 = asarray(tree_db2, ns2)
    } else {
        res.nm2 = bpeat_eat_cv_select(Z, X2, numStop, fold, N, nZ, nX2,
                                       tree_db2, ns2)
    }

    leaves2  = select(res.nm2, res.nm2[.,3]:==-1)
    nJ2      = rows(leaves2)
    lv_az    = leaves2[., ca2..ca2+nZ-1]
    lv_phat  = leaves2[., cy2..cy2+nX2-1]

    // Compute efficiency for each DMU
    b1 = J(N, 1, .)
    b2 = J(N, 1, .)
    bc = J(N, 1, .)

    for (k=1; k<=N; k++) {
        b1[k] = bpeat_beta_psi1(lv_a, lv_yhat, X12[k,.], Y[k,.], nX12, nY, nJ1)
        b2[k] = bpeat_beta_psi2(lv_az, lv_phat, X2[k,.], Z[k,.], nX2, nZ, nJ2)
        // Guillen et al. (2025): convert the good-output expansion DDF
        // to 1 - 1/(1+beta_Psi1) before combining the two technologies.
        bc[k] = delta1 * (b1[k] / (1 + b1[k])) + delta2 * b2[k]
    }

    res.beta_psi1 = b1
    res.beta_psi2 = b2
    res.beta_comb = bc
    res.nJ1       = nJ1
    res.nJ2       = nJ2
    res.n_nodes1  = rows(res.nm1)
    res.n_nodes2  = rows(res.nm2)
    return(res)
}

real scalar bpeat_valid_mean(real colvector v)
{
    real colvector w
    w = select(v, v :!= .)
    return(length(w) > 0 ? mean(w) : .)
}

void bpeat_display(struct BPEATResult scalar res,
                   string scalar x1names, string scalar x2names,
                   string scalar ynames,  string scalar znames,
                   real scalar delta1, real scalar delta2)
{
    real scalar N
    real colvector v1, v2, vc
    N  = rows(res.beta_psi1)
    v1 = select(res.beta_psi1, res.beta_psi1 :!= .)
    v2 = select(res.beta_psi2, res.beta_psi2 :!= .)
    vc = select(res.beta_comb,  res.beta_comb  :!= .)

    printf("{txt}\n")
    printf("{txt}By-production EAT (Guillen et al. 2025, EJOR)\n")
    printf("{hline 60}\n")
    printf("{txt}  Non-polluting inputs (x1) : {res}%s{txt}\n", x1names)
    printf("{txt}  Polluting inputs (x2)     : {res}%s{txt}\n", x2names)
    printf("{txt}  Good outputs (y)          : {res}%s{txt}\n", ynames)
    printf("{txt}  Bad outputs (z)           : {res}%s{txt}\n", znames)
    printf("{txt}  Observations              : {res}%g{txt}\n", N)
    printf("{txt}  delta1 / delta2           : {res}%g / %g{txt}\n", delta1, delta2)
    printf("{txt}  T1* leaves (Psi1)         : {res}%g{txt}  (nodes: %g)\n",
           res.nJ1, res.n_nodes1)
    printf("{txt}  T2* leaves (Psi2)         : {res}%g{txt}  (nodes: %g)\n",
           res.nJ2, res.n_nodes2)
    printf("{hline 60}\n")
    if (length(v1) > 0) {
        printf("{txt}  beta_Psi1 (good output DDF):\n")
        printf("{txt}    Mean= %7.4f  Min= %7.4f  Max= %7.4f\n",
               mean(v1), min(v1), max(v1))
    }
    if (length(v2) > 0) {
        printf("{txt}  beta_Psi2 (bad output DDF):\n")
        printf("{txt}    Mean= %7.4f  Min= %7.4f  Max= %7.4f\n",
               mean(v2), min(v2), max(v2))
    }
    if (length(vc) > 0) {
        printf("{txt}  Combined B = d1*b1/(1+b1) + d2*b2:\n")
        printf("{txt}    Mean= %7.4f  Min= %7.4f  Max= %7.4f\n",
               mean(vc), min(vc), max(vc))
    }
    printf("{hline 60}\n")
}

end
