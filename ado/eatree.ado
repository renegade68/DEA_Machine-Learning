*! eatree.ado  v0.1.5
*! Efficiency Analysis Trees (EAT) for Stata
*! Author: Choonjoo Lee, Korea National Defense University
*!         Nonsan-si, Republic of Korea; bloom.rampike@gmail.com
*! Method reference: Esteve, Aparicio, Rabasa, and Rodriguez-Sala (2020),
*!   "Efficiency analysis trees: A new methodology for estimating production
*!   frontiers through decision trees", Expert Systems with Applications 162:113783.
*!   DOI: 10.1016/j.eswa.2020.113783
*! Replicates R package 'eat' v0.1.4
*!
*! Syntax:
*!   eatree , x(varlist) y(varlist) [numstop(#) fold(#) maxdepth(#)
*!         maxleaves(#) gen(newvar) seed(#)]

program define eatree, eclass
    version 16.1
    syntax [if] [in] , ///
        X(varlist numeric) ///
        Y(varlist numeric) ///
        [ NUMStop(integer 5)   ///
          FOLD(integer 5)      ///
          MAXDepth(integer 0)  ///
          MAXLeaves(integer 0) ///
          FOLDID(varname numeric) ///
          GENerate(name)       ///
          SEED(integer 12345)  ]

    if "`generate'" != "" confirm new variable `generate'

    local xvars `x'
    local yvars `y'

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

    tempvar touse
    marksample touse, novarlist
    markout `touse' `xvars' `yvars' `foldid'
    qui count if `touse'
    local N = r(N)
    if `N' < 4 {
        di as err "Too few obs."
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

    tempname Xm Ym Fm
    mkmat `xvars' if `touse', matrix(`Xm')
    mkmat `yvars' if `touse', matrix(`Ym')
    if "`foldid'" != "" {
        mkmat `foldid' if `touse', matrix(`Fm')
        quietly assert `foldid'==floor(`foldid') & inrange(`foldid',1,`fold') if `touse'
    }
    else matrix `Fm'=J(1,1,.)

    if `fold' > 0 & "`foldid'" == "" set seed `seed'
    mata: _eatree_res = eatree_main(st_matrix("`Xm'"), st_matrix("`Ym'"), ///
        `numstop', `fold', `maxdepth', `maxleaves', st_matrix("`Fm'"))

    if "`generate'" != "" {
        qui gen double `generate' = . if `touse'
        mata: st_store(., "`generate'", "`touse'", _eatree_res.eff_scores)
    }

    mata: eatree_display(_eatree_res, "`xvars'", "`yvars'")

    ereturn clear
    ereturn post, esample(`touse')
    ereturn scalar N = `N'
    mata: st_numscalar("r(n_leaves)", _eatree_res.n_leaves)
    mata: st_numscalar("r(n_nodes)",  _eatree_res.n_nodes)
    mata: st_numscalar("r(mean_phi)", mean(_eatree_res.eff_scores))
    ereturn scalar n_leaves = r(n_leaves)
    ereturn scalar n_nodes  = r(n_nodes)
    ereturn scalar mean_phi = r(mean_phi)
    ereturn local  x "`xvars'"
    ereturn local  y "`yvars'"
    ereturn local  cmd "eatree"
    ereturn local  cmdline `"eatree `0'"'
    tempname nodes_tmp
    mata: st_matrix("`nodes_tmp'", _eatree_res.node_mat)
    ereturn matrix nodes = `nodes_tmp'
    tempname Fmat
    mata: st_matrix("`Fmat'", _eatree_res.eff_scores :* st_matrix("`Ym'"))
    ereturn matrix frontier = `Fmat'
end


// ==================================================================
mata:
mata set matastrict off

struct EATResult {
    real matrix    node_mat
    real colvector eff_scores
    real scalar    n_leaves
    real scalar    n_nodes
}

// node_mat columns:
//  1=id  2=F  3=SL  4=SR  5=xi  6=s  7=R
//  8..7+nX    = a[1..nX]
//  8+nX..7+2nX= b[1..nX]  (.=+Inf)
//  8+2nX..end = yhat[1..nY]
//
// tree_db structure (asarray with real keys):
//   key 1..ns  : tree state (nm matrix) at each pruning step
//   key 0      : 1x1 matrix holding ns
//   key -1     : row vector holding av[1..ns]

real scalar eatree_posNode(real matrix nm, real scalar id)
{
    real colvector v
    v = selectindex(nm[.,1] :== id)
    return(length(v) ? v[1] : 0)
}

real scalar eatree_isFinal(real matrix Xa, real colvector obs, real scalar ns)
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

real scalar eatree_mse(real matrix Y, real colvector obs,
                    real rowvector yhat, real scalar N, real scalar nY)
{
    real scalar err, r
    if (length(obs)==0) return(.)
    err = 0
    for (r=1; r<=nY; r++) err = err + sum((yhat[r] :- Y[obs,r]):^2)
    // The reference implementation rounds every node error to four decimals
    // before comparing candidate splits and computing pruning alphas.
    return(round(err/(N*nY), .0001))
}

real scalar eatree_comparePareto(real rowvector a1, real rowvector b1,
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

void eatree_estimNode(
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
        if (eatree_comparePareto(aL,bL,la,lb)==1) {
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
    tRes[1, cy..cy+nY-1]=yhL; tRes[1,7]=eatree_mse(Y,obsL,yhL,N,nY)

    tRes[2,1]=.; tRes[2,2]=prow[1,1]; tRes[2,3..6]=J(1,4,-1)
    tRes[2, ca..ca+nX-1]=aR; tRes[2, cb..cb+nX-1]=bR
    tRes[2, cy..cy+nY-1]=yhR; tRes[2,7]=eatree_mse(Y,obsR,yhR,N,nY)

    oL = obsL; oR = obsR
}

real rowvector eatree_rbranch(real matrix nm, real scalar id)
{
    real scalar pos
    real rowvector rL, rR
    pos = eatree_posNode(nm, id)
    if (pos==0 | nm[pos,3]==-1) return((nm[pos,7], 1))
    rL = eatree_rbranch(nm, nm[pos,3])
    rR = eatree_rbranch(nm, nm[pos,4])
    return((rL[1]+rR[1], rL[2]+rR[2]))
}

real scalar eatree_alpha(real matrix nm)
{
    real scalar amin, alpha_i, i
    real rowvector rb
    amin = .
    for (i=1; i<=rows(nm); i++) {
        if (nm[i,3]==-1) continue
        rb = eatree_rbranch(nm, nm[i,1])
        if (rb[2]<=1) continue
        alpha_i = (nm[i,7]-rb[1])/(rb[2]-1)
        if (alpha_i<0) alpha_i = 0
        if (amin==. | alpha_i<amin) amin = alpha_i
        if (amin==0) break
    }
    // With no internal node, alpha() in the R reference returns Inf.
    return(amin==. ? 1e300 : amin)
}

real rowvector eatree_predictor(real matrix nm, real rowvector xobs,
                              real scalar nX, real scalar nY)
{
    real scalar ti
    ti = 1
    while (nm[ti,3] != -1) {
        ti = (xobs[nm[ti,5]] < nm[ti,6]) ?
             eatree_posNode(nm, nm[ti,3]) : eatree_posNode(nm, nm[ti,4])
    }
    return(nm[ti, (8+2*nX)..(7+2*nX+nY)])
}

// -- grow full tree; store states in tree_db ----------------------
// tree_db is pre-allocated transmorphic (asarray) by caller
// On return: tree_db[1..ns] = states, tree_db[0] = ns, tree_db[-1] = av
void eatree_grow(
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
    nm[1,7] = eatree_mse(Y, aobs, nm[1,cy..cy+nY-1], N, nY)

    n_nodes=1; next_id=2; leaf_count=1
    ns = 1
    av = J(1, 2*N+4, .)
    asarray(tree_db, 1, nm)
    av[1] = eatree_alpha(nm)
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

        if (eatree_isFinal(X, pobs, numStop)) continue

        if (maxDepth > 0) {
            anc_id = prow[1,2]; d_cnt = 0
            while (anc_id != -1 & d_cnt < maxDepth+3) {
                anc_id = nm[eatree_posNode(nm,anc_id), 2]; d_cnt++
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
                eatree_estimNode(X,Y,oleaves,prow,xi,Sx[i],N,nX,nY,pobs,
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

        if (eatree_isFinal(X,boR,numStop)) {
            queue=(iR\queue)
        } else {
            queue=(queue\iR)
        }
        if (eatree_isFinal(X,boL,numStop)) {
            queue=(iL\queue)
        } else {
            queue=(queue\iL)
        }
        leaf_count++

        ns++
        asarray(tree_db, ns, nm)
        av[ns] = eatree_alpha(nm)
    }

    // store ns and av inside tree_db for caller to retrieve
    asarray(tree_db, 0, J(1,1,ns))
    asarray(tree_db, -1, av[1..ns])
}

real matrix eatree_pickTree(transmorphic db, real vector av,
                          real scalar ns, real scalar amid)
{
    // R reference order is full tree -> root.  Starting from the full tree,
    // retain each encountered state whose alpha is no larger than the
    // threshold; the last retained state is the simplest eligible subtree.
    real scalar best, j
    best = ns
    for (j=ns; j>=2; j--) {
        if (av[j]!=. & av[j] <= amid) {
            best = j
        }
    }
    return(asarray(db, best))
}

real matrix eatree_cv_select(
    real matrix X, real matrix Y,
    real scalar numStop, real scalar fold,
    real scalar N, real scalar nX, real scalar nY,
    transmorphic tree_db_full,
    real scalar  ns_full,
    real colvector fold_id_user
    )
{
    // Pruning-order CV with alpha midpoints:
    //   Growth order: state 1=root (simplest), state ns=full tree (most complex)
    //   Pruning order: state ns, ns-1, ..., 2, 1 (decreasing complexity)
    //   av[j] = transition alpha to prune growth state j to growth state j-1
    //   In pruning order: av[ns] < av[ns-1] < ... < av[2]  (INCREASING)
    //
    //   For cv_scores[i] (growth state i, i=2..ns):
    //     a1 = av_full[i]     (alpha to prune state i -> state i-1)
    //     a2 = av_full[i-1]   (alpha to prune state i-1 -> state i-2)
    //     amid = sqrt(a1 * a2)
    //     Use eatree_pickTree(fold_db, av_v, ns_v, amid) -> fold state
    //
    //   1-SE rule: find SIMPLEST state (smallest i) within min_cv + SE

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

    // Match generateLv(): randomly permute, then take equally sized
    // contiguous validation blocks.  Any remainder stays in every training
    // sample, as in eat 0.1.4.
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
        eatree_grow(Xt,Yt,numStop,0,0,rows(Xt),nX,nY, db_v)
        ns_v = asarray(db_v, 0)[1,1]
        av_v = asarray(db_v, -1)
        asarray(fold_dbs, v*3-2, db_v)
        asarray(fold_dbs, v*3-1, av_v)
        asarray(fold_dbs, v*3,   J(1,1,ns_v))
    }

    // The R list is the reverse of growth order.  State i receives the score
    // at the geometric midpoint of alpha(state i) and alpha(state i-1).
    cv_scores = J(1, ns_full, 1e15)
    for (i=2; i<=ns_full; i++) {
        a1 = av_full[i]
        a2 = av_full[i-1]
        if (a1==.) a1=0
        if (a2==.) a2=0
        amid = sqrt(a1*a2)
        Rcv = 0
        for (v=1; v<=fold; v++) {
            val  = selectindex(fold_id :== v)
            db_v = asarray(fold_dbs, v*3-2)
            av_v = asarray(fold_dbs, v*3-1)
            ns_v = asarray(fold_dbs, v*3)[1,1]
            nm_v = eatree_pickTree(db_v, av_v, ns_v, amid)
            Xv   = X[val,.]; Yv = Y[val,.]
            for (reg=1; reg<=rows(Xv); reg++) {
                pred = eatree_predictor(nm_v, Xv[reg,.], nX, nY)
                for (r=1; r<=nY; r++) Rcv = Rcv+(Yv[reg,r]-pred[r])^2
            }
        }
        cv_scores[i] = Rcv/(N*nY)
    }

    min_cv     = min(cv_scores)
    ixbest     = selectindex(cv_scores :== min_cv)
    // scores() traverses full -> simple and retains the first strict minimum.
    best_state = max(ixbest)
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
        nm_b = eatree_pickTree(db_v, av_v, ns_v, amid2)
        Xv   = X[val,.]; Yv = Y[val,.]
        for (reg=1; reg<=rows(Xv); reg++) {
            pred = eatree_predictor(nm_b, Xv[reg,.], nX, nY)
            for (rr=1; rr<=nY; rr++) {
                dif1 = (Yv[reg,rr]-pred[rr])^2
                s2   = s2+(dif1-Tk_score)^2
            }
        }
    }
    SE = sqrt(s2/(N*nY)/(N*nY))

    // 1-SE rule: find SIMPLEST state (smallest i = fewest leaves) within margin
    // Growth order: smaller i = simpler. Search from i=2 upward.
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

real colvector eatree_efficiency(real matrix X, real matrix Y,
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
        eff[k] = (theta_best==.) ? 1 : theta_best
    }
    return(eff)
}

struct EATResult scalar eatree_main(
    real matrix X, real matrix Y,
    real scalar numStop, real scalar fold,
    real scalar maxDepth, real scalar maxLeaves,
    real colvector fold_id_user)
{
    struct EATResult scalar res
    real scalar N, nX, nY, ns
    real vector  av
    transmorphic tree_db
    real matrix  final_tree

    N=rows(X); nX=cols(X); nY=cols(Y)
    tree_db = asarray_create("real")
    eatree_grow(X,Y,numStop,maxDepth,maxLeaves,N,nX,nY, tree_db)
    ns = asarray(tree_db, 0)[1,1]
    av = asarray(tree_db, -1)

    if (ns<=1 | fold<=0) {
        final_tree = asarray(tree_db, ns)
    } else {
        final_tree = eatree_cv_select(X,Y,numStop,fold,N,nX,nY,tree_db,ns,fold_id_user)
    }

    res.node_mat   = final_tree
    res.eff_scores = eatree_efficiency(X,Y,final_tree,nX,nY)
    res.n_leaves   = sum(final_tree[.,3]:==-1)
    res.n_nodes    = rows(final_tree)
    return(res)
}

void eatree_display(struct EATResult scalar res,
                 string scalar xnames, string scalar ynames)
{
    printf("{txt}\n")
    printf("{txt}Efficiency Analysis Trees (EAT)\n")
    printf("{hline 50}\n")
    printf("{txt}  Inputs  : {res}%s{txt}\n", xnames)
    printf("{txt}  Outputs : {res}%s{txt}\n", ynames)
    printf("{hline 50}\n")
    printf("{txt}  Total nodes    : {res}%g{txt}\n", res.n_nodes)
    printf("{txt}  Internal nodes : {res}%g{txt}\n", res.n_nodes-res.n_leaves)
    printf("{txt}  Leaf nodes     : {res}%g{txt}\n", res.n_leaves)
    printf("{hline 50}\n")
    printf("{txt}  Efficiency (BCC-OUT):\n")
    printf("{txt}    Min=%7.4f  Mean=%7.4f  Max=%7.4f\n",
           min(res.eff_scores), mean(res.eff_scores), max(res.eff_scores))
    printf("{txt}    (>=1.000 ; efficient DMU = 1.000)\n")
    printf("{hline 50}\n")
}

end
