*! ceat.ado  v0.2.1
*! Convexified Efficiency Analysis Trees (CEAT) for Stata
*! Author: Choonjoo Lee, Korea National Defense University
*!         Nonsan-si, Republic of Korea; bloom.rampike@gmail.com
*! Self-contained: EAT Mata functions embedded directly.
*! Method reference: Aparicio, Esteve, Rodriguez-Sala, and Zofio (2021),
*!   "The Estimation of Productive Efficiency Through Machine Learning
*!   Techniques: Efficiency Analysis Trees", Data-Enabled Analytics: DEA for
*!   Big Data, pp. 51-92. DOI: 10.1007/978-3-030-75162-3_3
*!
*! CEAT takes the EAT step-function frontier and convexifies it via LP:
*!   For each DMU k, solve output-oriented VRS DEA using EAT leaf
*!   reference points (lower-bound inputs, frontier outputs) as the
*!   reference set.  Uses Stata 16+ LinearProgram() class.
*!
*! Syntax:
*!   ceat , x(varlist) y(varlist) [numstop(#) fold(#) maxdepth(#)
*!          maxleaves(#) rts(vrs|crs|drs|irs) gen(newvar)
*!          genlambda(stub) foldid(varname) seed(#)]

program define ceat, eclass
    version 16.1
    syntax [if] [in] , ///
        X(varlist numeric) ///
        Y(varlist numeric) ///
        [ NUMStop(integer 5)    ///
          FOLD(integer 5)       ///
          MAXDepth(integer 0)   ///
          MAXLeaves(integer 0)  ///
          RTS(string)           ///
          GENerate(name)        ///
          GENlambda(string)     ///
          FOLDID(varname numeric) ///
          SEED(integer 12345)   ]

    if "`generate'" != "" confirm new variable `generate'

    local xvars `x'
    local yvars `y'
    local nX : word count `xvars'
    local nY : word count `yvars'

    if `numstop'<0 {
        di as err "numstop() must be >= 0"
        exit 198
    }
    if `fold'==1 | `fold'<0 {
        di as err "fold() must be 0 or >= 2"
        exit 198
    }
    if `maxdepth'<0 | `maxleaves'<0 {
        di as err "maxdepth() and maxleaves() must be >= 0"
        exit 198
    }

    if "`rts'" == "" local rts "vrs"
    local rts = lower("`rts'")
    if !inlist("`rts'", "vrs","crs","drs","irs") {
        di as err "rts() must be vrs, crs, drs, or irs"
        exit 198
    }

    tempvar touse
    marksample touse, novarlist
    markout `touse' `xvars' `yvars' `foldid'
    qui count if `touse'
    local N = r(N)
    if `N' < 4 {
        di as err "Too few observations."
        exit 2001
    }
    if `fold'>`N' {
        di as err "fold() may not exceed the estimation-sample size"
        exit 198
    }
    foreach v of varlist `xvars' {
        quietly count if `touse' & `v'<0
        if r(N) {
            di as err "input variables must be nonnegative"
            exit 459
        }
    }
    foreach v of varlist `yvars' {
        quietly count if `touse' & `v'<=0
        if r(N) {
            di as err "output variables must be strictly positive for radial efficiency scores"
            exit 459
        }
    }

    tempname Xm Ym
    mkmat `xvars' if `touse', matrix(`Xm')
    mkmat `yvars' if `touse', matrix(`Ym')
    tempname Fm
    if "`foldid'" != "" {
        mkmat `foldid' if `touse', matrix(`Fm')
        quietly assert `foldid'==floor(`foldid') & inrange(`foldid',1,`fold') if `touse'
    }
    else matrix `Fm' = J(1,1,.)

    if `fold' > 0 & "`foldid'" == "" set seed `seed'
    mata: _ceat_res = ceat_main( ///
        st_matrix("`Xm'"), st_matrix("`Ym'"), ///
        `numstop', `fold', `maxdepth', `maxleaves', "`rts'", st_matrix("`Fm'"))
    mata: st_numscalar("r(_ceat_failed)", _ceat_res.n_failed)
    if r(_ceat_failed) == `N' {
        di as err "CEAT optimization failed for every observation; no results were stored"
        exit 498
    }

    // Preflight the entire lambda output set before writing any result.
    if "`genlambda'" != "" {
        mata: st_numscalar("r(_ceat_nJ)", cols(_ceat_res.lambda_mat))
        local ceat_nJ = r(_ceat_nJ)
        forvalues j=1/`ceat_nJ' {
            confirm new variable `genlambda'`j'
            if "`generate'" != "" & "`generate'" == "`genlambda'`j'" {
                di as err "generate() and genlambda() output names must be distinct"
                exit 198
            }
        }
    }

    if "`generate'" != "" {
        qui gen double `generate' = . if `touse'
        mata: st_store(., "`generate'", "`touse'", _ceat_res.eff)
        label var `generate' "CEAT efficiency (output-oriented)"
    }
    if "`genlambda'" != "" {
        mata: ceat_store_lambda(_ceat_res, "`genlambda'", "`touse'")
    }

    mata: ceat_display(_ceat_res, "`xvars'", "`yvars'", "`rts'")

    ereturn clear
    ereturn post, esample(`touse')
    ereturn scalar N          = `N'
    mata: st_numscalar("r(n_leaves)", _ceat_res.n_leaves)
    mata: st_numscalar("r(n_nodes)",  _ceat_res.n_nodes)
    mata: st_numscalar("r(n_failed)", _ceat_res.n_failed)
    ereturn scalar n_leaves   = r(n_leaves)
    ereturn scalar n_nodes    = r(n_nodes)
    ereturn scalar n_failed   = r(n_failed)
    mata: st_numscalar("r(mean_phi)", mean(select(_ceat_res.eff, _ceat_res.eff :!= .)))
    ereturn scalar mean_phi   = r(mean_phi)
    ereturn local  rts     "`rts'"
    ereturn local  solver  "LinearProgram"
    ereturn local  x       "`xvars'"
    ereturn local  y       "`yvars'"
    ereturn local  cmd     "ceat"
    ereturn local  cmdline `"ceat `0'"'
    tempname nodes_tmp
    mata: st_matrix("`nodes_tmp'", _ceat_res.node_mat)
    ereturn matrix nodes = `nodes_tmp'
    tempname Fmat
    mata: st_matrix("`Fmat'", _ceat_res.eff :* st_matrix("`Ym'"))
    ereturn matrix frontier = `Fmat'
end


mata:
mata set matastrict off

// ======================================================================
// EAT Mata functions (embedded from eatree.ado v0.1.4)
// node_mat columns:
//  1=id  2=F(parent)  3=SL  4=SR  5=xi  6=s  7=R(MSE)
//  8..7+nX    = a[1..nX]  (lower bounds)
//  8+nX..7+2nX= b[1..nX]  (upper bounds, .=+Inf)
//  8+2nX..end = yhat[1..nY] (frontier outputs)
// ======================================================================

real scalar ceat_eat_posNode(real matrix nm, real scalar id)
{
    real colvector v
    v = selectindex(nm[.,1] :== id)
    return(length(v) ? v[1] : 0)
}

real scalar ceat_eat_isFinal(real matrix Xa, real colvector obs, real scalar ns)
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

real scalar ceat_eat_mse(real matrix Y, real colvector obs,
                         real rowvector yhat, real scalar N, real scalar nY)
{
    real scalar err, r
    if (length(obs)==0) return(.)
    err = 0
    for (r=1; r<=nY; r++) err = err + sum((yhat[r] :- Y[obs,r]):^2)
    return(round(err/(N*nY), .0001))
}

real scalar ceat_eat_comparePareto(real rowvector a1, real rowvector b1,
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

void ceat_eat_estimNode(
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
        if (ceat_eat_comparePareto(aL,bL,la,lb)==1) {
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
    tRes[1, cy..cy+nY-1]=yhL; tRes[1,7]=ceat_eat_mse(Y,obsL,yhL,N,nY)

    tRes[2,1]=.; tRes[2,2]=prow[1,1]; tRes[2,3..6]=J(1,4,-1)
    tRes[2, ca..ca+nX-1]=aR; tRes[2, cb..cb+nX-1]=bR
    tRes[2, cy..cy+nY-1]=yhR; tRes[2,7]=ceat_eat_mse(Y,obsR,yhR,N,nY)

    oL = obsL; oR = obsR
}

real rowvector ceat_eat_rbranch(real matrix nm, real scalar id)
{
    real scalar pos
    real rowvector rL, rR
    pos = ceat_eat_posNode(nm, id)
    if (pos==0 | nm[pos,3]==-1) return((nm[pos,7], 1))
    rL = ceat_eat_rbranch(nm, nm[pos,3])
    rR = ceat_eat_rbranch(nm, nm[pos,4])
    return((rL[1]+rR[1], rL[2]+rR[2]))
}

real scalar ceat_eat_alpha(real matrix nm)
{
    real scalar amin, alpha_i, i
    real rowvector rb
    amin = .
    for (i=1; i<=rows(nm); i++) {
        if (nm[i,3]==-1) continue
        rb = ceat_eat_rbranch(nm, nm[i,1])
        if (rb[2]<=1) continue
        alpha_i = (nm[i,7]-rb[1])/(rb[2]-1)
        if (alpha_i<0) alpha_i = 0
        if (amin==. | alpha_i<amin) amin = alpha_i
        if (amin==0) break
    }
    return(amin==. ? 1e300 : amin)
}

void ceat_eat_grow(
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
    nm[1,7] = ceat_eat_mse(Y, aobs, nm[1,cy..cy+nY-1], N, nY)

    n_nodes=1; next_id=2; leaf_count=1
    ns = 1
    av = J(1, 2*N+4, .)
    asarray(tree_db, 1, nm)
    av[1] = ceat_eat_alpha(nm)
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

        if (ceat_eat_isFinal(X, pobs, numStop)) continue

        if (maxDepth > 0) {
            anc_id = prow[1,2]; d_cnt = 0
            while (anc_id != -1 & d_cnt < maxDepth+3) {
                anc_id = nm[ceat_eat_posNode(nm,anc_id), 2]; d_cnt++
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
                ceat_eat_estimNode(X,Y,oleaves,prow,xi,Sx[i],N,nX,nY,pobs,
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

        if (ceat_eat_isFinal(X,boR,numStop)) {
            queue=(iR\queue)
        } else {
            queue=(queue\iR)
        }
        if (ceat_eat_isFinal(X,boL,numStop)) {
            queue=(iL\queue)
        } else {
            queue=(queue\iL)
        }
        leaf_count++

        ns++
        asarray(tree_db, ns, nm)
        av[ns] = ceat_eat_alpha(nm)
    }

    asarray(tree_db, 0, J(1,1,ns))
    asarray(tree_db, -1, av[1..ns])
}

real matrix ceat_eat_pickTree(transmorphic db, real vector av,
                               real scalar ns, real scalar amid)
{
    real scalar best, j
    best = ns
    for (j=ns; j>=2; j--) {
        if (av[j]!=. & av[j] <= amid) {
            best = j
        }
    }
    return(asarray(db, best))
}

real rowvector ceat_eat_predictor(real matrix nm, real rowvector xobs,
                                   real scalar nX, real scalar nY)
{
    real scalar ti, next_ti
    ti = 1
    while (nm[ti,3] != -1) {
        if (xobs[nm[ti,5]] < nm[ti,6]) {
            next_ti = ceat_eat_posNode(nm, nm[ti,3])
        } else {
            next_ti = ceat_eat_posNode(nm, nm[ti,4])
        }
        ti = next_ti
    }
    return(nm[ti, (8+2*nX)..(7+2*nX+nY)])
}

real matrix ceat_eat_cv_select(
    real matrix X, real matrix Y,
    real scalar numStop, real scalar fold,
    real scalar N, real scalar nX, real scalar nY,
    transmorphic tree_db_full,
    real scalar  ns_full,
    real colvector fold_id_user
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

    if (rows(fold_id_user)==N) {
        fold_id=fold_id_user
    }
    else {
        fold_id = J(N,1,0)
        perm = order(runiform(N,1),1)
        nfold = floor(N/fold)
        for (v=1; v<=fold; v++) {
            fold_id[perm[((v-1)*nfold+1)..(v*nfold)]] = J(nfold,1,v)
        }
    }

    fold_dbs = asarray_create("real")

    for (v=1; v<=fold; v++) {
        tr = selectindex(fold_id :!= v)
        Xt = X[tr,.]; Yt = Y[tr,.]
        db_v = asarray_create("real")
        ceat_eat_grow(Xt,Yt,numStop,0,0,rows(Xt),nX,nY, db_v)
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
            nm_v = ceat_eat_pickTree(db_v, av_v, ns_v, amid)
            Xv   = X[val,.]; Yv = Y[val,.]
            for (reg=1; reg<=rows(Xv); reg++) {
                pred = ceat_eat_predictor(nm_v, Xv[reg,.], nX, nY)
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
        nm_b = ceat_eat_pickTree(db_v, av_v, ns_v, amid2)
        Xv   = X[val,.]; Yv = Y[val,.]
        for (reg=1; reg<=rows(Xv); reg++) {
            pred = ceat_eat_predictor(nm_b, Xv[reg,.], nX, nY)
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
// CEAT-specific Mata code
// ======================================================================

struct CEATResult {
    real matrix    node_mat
    real colvector eff
    real matrix    lambda_mat
    real scalar    n_leaves
    real scalar    n_nodes
    real scalar    n_failed
}

// LP solver for one DMU
// Variables: [lambda_1..lambda_nJ, phi]  (nJ+1 total)
// Objective: min -phi  (equivalent to max phi)
// Constraints:
//   lv_a' * lambda <= x_k          [nX input constraints]
//   -lv_yhat' * lambda + phi*y_k <= 0  [nY output constraints]
//   sum(lambda) = 1 (VRS)
//   lambda >= 0, phi >= 1
void ceat_lp_solve(
    real matrix    lv_a,      // nJ x nX leaf lower bounds
    real matrix    lv_yhat,   // nJ x nY leaf frontier outputs
    real colvector x0,        // nX x 1 DMU inputs
    real colvector y0,        // nY x 1 DMU outputs
    real scalar    nX,
    real scalar    nY,
    real scalar    nJ,
    string scalar  rts,
    real scalar    phi_out,
    real rowvector lam_out
    )
{
    real rowvector c, lb, ub, sol
    real matrix    A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar    rc
    class LinearProgram scalar S

    // Objective: min -phi (last variable)
    c  = (J(1, nJ, 0), -1)
    lb = (J(1, nJ, 0),  1)    // lambda>=0, phi>=1
    ub =  J(1, nJ+1, .)       // unbounded above

    // Input constraints: lv_a' * lambda <= x0
    A_ineq = (lv_a', J(nX, 1, 0))        // nX x (nJ+1)
    b_ineq = x0                            // nX x 1

    // Output constraints: -lv_yhat' * lambda + phi*y0 <= 0
    A_ineq = A_ineq \ (-lv_yhat', y0)    // (nX+nY) x (nJ+1)
    b_ineq = b_ineq \ J(nY, 1, 0)        // (nX+nY) x 1

    // Convexity constraint (RTS)
    A_eq = J(0, nJ+1, .)
    b_eq = J(0, 1, .)

    if (rts == "vrs") {
        A_eq = (J(1, nJ, 1), 0)
        b_eq = 1
    } else if (rts == "drs") {
        A_ineq = A_ineq \ (J(1, nJ, 1), 0)
        b_ineq = b_ineq \ 1
    } else if (rts == "irs") {
        A_ineq = A_ineq \ (J(1, nJ, -1), 0)
        b_ineq = b_ineq \ (-1)
    }
    // crs: no convexity constraint

    S = LinearProgram()
    S.setMaxOrMin("min")
    S.setCoefficients(c)
    S.setInequality(A_ineq, b_ineq)
    if (rows(A_eq) > 0) {
        S.setEquality(A_eq, b_eq)
    }
    S.setBounds(lb, ub)
    S.setTrace("off")
    rc = S.optimize()

    if (S.converged() & S.errorcode()==0) {
        sol     = S.parameters()
        phi_out = -S.value()
        lam_out = sol[1..nJ]
    } else {
        phi_out = .
        lam_out = J(1, nJ, .)
    }
}

// Main CEAT function
struct CEATResult scalar ceat_main(
    real matrix   X,
    real matrix   Y,
    real scalar   numStop,
    real scalar   fold,
    real scalar   maxDepth,
    real scalar   maxLeaves,
    string scalar rts,
    real colvector fold_id_user)
{
    struct CEATResult scalar res
    real scalar  N, nX, nY, nJ, k, n_failed, ca, cy, ns
    real matrix  leaves, final_tree, lv_a, lv_yhat
    real colvector eff
    real matrix  lam_mat
    real scalar  phi_k
    real rowvector lam_k
    real colvector x0, y0
    transmorphic tree_db

    N=rows(X); nX=cols(X); nY=cols(Y)
    ca=8; cy=8+2*nX

    // Step 1: grow EAT tree
    tree_db = asarray_create("real")
    ceat_eat_grow(X, Y, numStop, maxDepth, maxLeaves, N, nX, nY, tree_db)
    ns = asarray(tree_db, 0)[1,1]

    if (ns <= 1 | fold <= 0) {
        final_tree = asarray(tree_db, ns)
    } else {
        final_tree = ceat_eat_cv_select(X, Y, numStop, fold, N, nX, nY, tree_db, ns,
                                        fold_id_user)
    }

    // Step 2: extract leaf reference points
    leaves  = select(final_tree, final_tree[.,3]:==-1)
    nJ      = rows(leaves)
    lv_a    = leaves[., ca..ca+nX-1]      // nJ x nX lower bounds
    lv_yhat = leaves[., cy..cy+nY-1]      // nJ x nY frontier outputs

    // Step 3: solve CEAT LP for each DMU
    eff      = J(N, 1, .)
    lam_mat  = J(N, nJ, .)
    n_failed = 0

    for (k=1; k<=N; k++) {
        x0 = X[k,.]'
        y0 = Y[k,.]'
        ceat_lp_solve(lv_a, lv_yhat, x0, y0, nX, nY, nJ, rts,
                      phi_k, lam_k)
        if (phi_k == .) {
            n_failed++
        }
        eff[k]       = phi_k
        lam_mat[k,.] = lam_k
    }

    res.node_mat   = final_tree
    res.eff        = eff
    res.lambda_mat = lam_mat
    res.n_leaves   = nJ
    res.n_nodes    = rows(final_tree)
    res.n_failed   = n_failed
    return(res)
}

void ceat_store_lambda(struct CEATResult scalar res,
                       string scalar stub, string scalar touse)
{
    real scalar j, nJ
    string scalar vname
    nJ = cols(res.lambda_mat)
    for (j=1; j<=nJ; j++) {
        vname = stub + strtrim(strofreal(j))
        stata("quietly generate double " + vname + " = .")
        st_store(., vname, touse, res.lambda_mat[., j])
    }
}

void ceat_display(struct CEATResult scalar res,
                  string scalar xnames, string scalar ynames,
                  string scalar rts)
{
    real scalar N, eff_min, eff_max, eff_mean, n_eff1
    real colvector valid_eff
    N = rows(res.eff)
    valid_eff = select(res.eff, res.eff :!= .)
    if (length(valid_eff) == 0) {
        printf("{err}All LP solvers failed.\n")
        return
    }
    eff_mean = mean(valid_eff)
    eff_min  = min(valid_eff)
    eff_max  = max(valid_eff)
    n_eff1   = sum(abs(valid_eff :- 1) :<= 1e-6)

    printf("{txt}\n")
    printf("{txt}Convexified Efficiency Analysis Trees (CEAT)\n")
    printf("{hline 55}\n")
    printf("{txt}  Inputs      : {res}%s{txt}\n", xnames)
    printf("{txt}  Outputs     : {res}%s{txt}\n", ynames)
    printf("{txt}  Observations: {res}%g{txt}\n", N)
    printf("{txt}  RTS         : {res}%s{txt}\n", strupper(rts))
    printf("{txt}  EAT leaves  : {res}%g{txt}\n", res.n_leaves)
    printf("{txt}  EAT nodes   : {res}%g{txt}\n", res.n_nodes)
    if (res.n_failed > 0) {
        printf("{txt}  LP failures : {err}%g{txt}\n", res.n_failed)
    }
    printf("{hline 55}\n")
    printf("{txt}  Efficiency (output-oriented %s):\n", strupper(rts))
    printf("{txt}    Mean= %7.4f  Min= %7.4f  Max= %7.4f\n",
           eff_mean, eff_min, eff_max)
    printf("{txt}    Efficient DMUs (score ~= 1): %g\n", n_eff1)
    printf("{hline 55}\n")
}

end
