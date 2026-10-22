*! svfeat.ado  v0.2.1
*! Support Vector Frontier (SVF), published model (9)
*! Author: Choonjoo Lee, Korea National Defense University
*!         Nonsan-si, Republic of Korea; bloom.rampike@gmail.com
*! Method reference: Valero-Carreras, Aparicio, and Guerrero (2021),
*!   "Support vector frontiers: A new approach for estimating production
*!   functions through support vector machines", Omega 104:102490.
*!   DOI: 10.1016/j.omega.2021.102490
*!
*! Estimates a non-decreasing step-function production frontier using the
*! equal-width dominated-cell grid transformation in published model (9).
*! 
*! Let phi_g(x) indicate whether grid cell g is weakly dominated by x.  With
*! free grid weights w_g, the LP formulation (per output) is:
*!   min  sum_g |w_g| + C*sum_i xi_i
*!   s.t. y_i <= sum_g w_g*phi_g(x_i) <= y_i + eps + xi_i
*!        adjacent grid values are nondecreasing; xi_i >= 0.
*! Free weights are split into positive and negative parts to represent the
*! exact L1 objective in Mata's LinearProgram class.
*!
*! Efficiency (output-oriented, >= 1):
*!   phi_i = max(1, min_r{ f_r(x_i) / y_ir })
*!   epsilon-insensitive: DMU efficient if f(x_i) <= y_i + eps
*!
*! Syntax:
*!   svfeat, x(varlist) y(varlist) [gen(name) d(#) c(real) eps(real) seed(#)]
*!
*!   x(varlist)   variable inputs
*!   y(varlist)   good outputs
*!   gen(name)    store efficiency phi in dataset (>= 1)
*!   d(#)         equal-width intervals per input (default 3)
*!   nu(real)     deprecated compatibility option; ignored
*!   c(real)      slack penalty (default 1.0, larger = tighter frontier)
*!   eps(real)    epsilon insensitivity band (default 0.0)
*!   seed(#)      random seed (default 12345)

program define svfeat, eclass
    version 16.1
    syntax [if] [in] , ///
        X(varlist numeric)     ///
        Y(varlist numeric)     ///
        [ GEN(name)            ///
          NU(real 0.5)         ///
          D(integer 3)         ///
          C(real 1.0)          ///
          EPS(real 0.0)        ///
          SEED(integer 12345)  ///
          TUNE                 ///
          CGRID(numlist >0)    ///
          EPSGRID(numlist >=0) ///
          DGRID(numlist integer >0) ///
          FOLDS(integer 5)     ]

    if "`gen'" != "" confirm new variable `gen'

    local xvars `x'
    local yvars `y'
    local nX : word count `xvars'
    local nY : word count `yvars'

    if `nu' < 0 {
        di as err "nu() must be >= 0"
        exit 198
    }
    if `d' < 1 {
        di as err "d() must be >= 1"
        exit 198
    }
    if `nu' != 0.5 {
        di as txt "note: nu() is deprecated and ignored by the current SVF estimator"
    }
    if `c' <= 0 {
        di as err "c() must be > 0"
        exit 198
    }
    if `eps' < 0 {
        di as err "eps() must be >= 0"
        exit 198
    }
    if `folds' < 2 {
        di as err "folds() must be >= 2"
        exit 198
    }
    if "`cgrid'" == "" local cgrid ".1 .5 1 2 5"
    if "`epsgrid'" == "" local epsgrid "0 .001 .01 .1 .2"
    if "`dgrid'" == "" local dgrid "1 2 3"
    local maxd `d'
    if "`tune'" != "" {
        foreach dd of numlist `dgrid' {
            if `dd' > `maxd' local maxd `dd'
        }
    }
    local maxcells = (`maxd'+1)^`nX'
    if `maxcells' > 2500 {
        di as err "SVF grid would contain `maxcells' cells; maximum is 2500"
        di as err "reduce d()/dgrid() or the number of input variables"
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
    if "`tune'" != "" & `folds' > `N' {
        di as err "folds() may not exceed the number of complete observations (`N')"
        exit 198
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

    if "`tune'" != "" set seed `seed'

    if "`tune'" != "" {
        mata: _svfeat_res = svfeat_cv_main( ///
            st_matrix("`Xm'"), st_matrix("`Ym'"), ///
            strtoreal(tokens("`cgrid'")), strtoreal(tokens("`epsgrid'")), ///
            strtoreal(tokens("`dgrid'")), `folds', `seed')
    }
    else {
        mata: _svfeat_res = svfeat_main( ///
            st_matrix("`Xm'"), st_matrix("`Ym'"), ///
            `d', `c', `eps')
    }

    mata: st_numscalar("r(_svfeat_fit_failed)", _svfeat_res.n_failed)
    if r(_svfeat_fit_failed) > 0 {
        di as err "SVF optimization failed; no efficiency scores were stored"
        exit 498
    }

    if "`gen'" != "" {
        quietly gen double `gen' = . if `touse'
        mata: st_store(., "`gen'", "`touse'", _svfeat_res.phi)
        label var `gen' "SVF efficiency (phi >= 1)"
    }

    mata: svfeat_display(_svfeat_res, "`xvars'", "`yvars'", `d', `c', `eps')

    ereturn clear
    ereturn post, esample(`touse')
    ereturn scalar N          = `N'
    mata: st_numscalar("r(mean_phi)",  svfeat_valid_mean(_svfeat_res.phi))
    mata: st_numscalar("r(n_failed)",  _svfeat_res.n_failed)
    mata: st_numscalar("r(n_svec)",    _svfeat_res.n_svec)
    mata: st_numscalar("r(max_violation)", _svfeat_res.max_violation)
    mata: st_numscalar("r(objective)", _svfeat_res.objective)
    mata: st_numscalar("r(selected_C)", _svfeat_res.selected_C)
    mata: st_numscalar("r(selected_eps)", _svfeat_res.selected_eps)
    mata: st_numscalar("r(selected_d)", _svfeat_res.selected_d)
    mata: st_numscalar("r(n_folds)", _svfeat_res.n_folds)
    mata: st_numscalar("r(tuned)", _svfeat_res.tuned)
    mata: st_numscalar("r(cv_failures)", _svfeat_res.cv_failures)
    tempname Wmat Fmat GXmat GFmat
    mata: st_matrix("`Wmat'", _svfeat_res.W)
    mata: st_matrix("`Fmat'", _svfeat_res.F)
    mata: st_matrix("`GXmat'", _svfeat_res.grid_X)
    mata: st_matrix("`GFmat'", _svfeat_res.grid_F)
    if "`tune'" != "" {
        tempname CVmat CVFoldmat CVPredmat Foldmat
        mata: st_matrix("`CVmat'", _svfeat_res.cv_results)
        mata: st_matrix("`CVFoldmat'", _svfeat_res.cv_fold_results)
        mata: st_matrix("`CVPredmat'", _svfeat_res.cv_predictions)
        mata: st_matrix("`Foldmat'", _svfeat_res.fold_id)
    }
    ereturn scalar mean_phi  = r(mean_phi)
    ereturn scalar n_failed  = r(n_failed)
    ereturn scalar n_svec    = r(n_svec)
    ereturn scalar max_violation = r(max_violation)
    ereturn scalar objective = r(objective)
    ereturn matrix coefficients = `Wmat'
    ereturn matrix frontier = `Fmat'
    ereturn matrix grid_inputs = `GXmat'
    ereturn matrix grid_frontier = `GFmat'
    if "`tune'" != "" {
        matrix colnames `CVmat' = C eps d mse
        matrix colnames `CVFoldmat' = C eps d fold mse
        matrix colnames `CVPredmat' = C eps d fold obs output prediction
        matrix colnames `Foldmat' = fold
        ereturn matrix cv_results = `CVmat'
        ereturn matrix cv_fold_results = `CVFoldmat'
        ereturn matrix cv_predictions = `CVPredmat'
        ereturn matrix fold_id = `Foldmat'
    }
    ereturn scalar d         = r(selected_d)
    ereturn scalar c_pen     = r(selected_C)
    ereturn scalar eps       = r(selected_eps)
    ereturn scalar folds     = r(n_folds)
    ereturn scalar tuned     = r(tuned)
    ereturn scalar cv_failures = r(cv_failures)
    ereturn local  x         "`xvars'"
    ereturn local  y         "`yvars'"
    ereturn local  solver    "LinearProgram"
    ereturn local  cmd       "svfeat"
    ereturn local  cmdline   `"svfeat `0'"'
end


mata:
mata set matastrict off

struct SVFResult {
    real colvector phi       // n x 1  output-oriented efficiency (>= 1)
    real matrix    F         // n x nY frontier predictions
    real scalar    n_failed  // LP failures
    real scalar    n_svec    // support vectors (alpha > threshold)
    real matrix    W
    real matrix    grid_X
    real matrix    grid_F
    real scalar    max_violation
    real scalar    objective
    real matrix    cv_results
    real matrix    cv_fold_results
    real matrix    cv_predictions
    real colvector fold_id
    real scalar    selected_C
    real scalar    selected_eps
    real scalar    selected_d
    real scalar    n_folds
    real scalar    tuned
    real scalar    cv_failures
}

real matrix svfeat_grid_index(real scalar m, real scalar d)
{
    real scalar q,g,j,base,den
    real matrix idx
    base=d+1; q=base^m; idx=J(q,m,0)
    for (g=1;g<=q;g++) {
        for (j=1;j<=m;j++) {
            den=base^(m-j)
            idx[g,j]=mod(floor((g-1)/den),base)
        }
    }
    return(idx)
}

real matrix svfeat_grid_phi(real matrix X, real scalar d, real matrix ranges,
                            real matrix idx)
{
    real scalar N,m,q,i,j,g,width,ok
    real rowvector cell
    real matrix Phi
    N=rows(X); m=cols(X); q=rows(idx); Phi=J(N,q,0)
    for (i=1;i<=N;i++) {
        cell=J(1,m,0)
        for (j=1;j<=m;j++) {
            width=(ranges[j,2]-ranges[j,1])/d
            if (X[i,j]<ranges[j,1]) cell[j]=-1
            else if (width<=1e-14) cell[j]=d
            else cell[j]=min((d,floor((X[i,j]-ranges[j,1])/width+1e-12)))
        }
        for (g=1;g<=q;g++) {
            ok=1
            for (j=1;j<=m;j++) {
                if (idx[g,j]>cell[j]) {
                    ok=0
                    break
                }
            }
            Phi[i,g]=ok
        }
    }
    return(Phi)
}

void svfeat_model9_solve(real matrix Phi, real colvector y,
                         real matrix idx, real scalar C, real scalar eps,
                         real colvector w_out, real colvector F_out,
                         real scalar obj_out, real scalar viol_out,
                         real scalar conv_out)
{
    real scalar N,q,m,nv,i,g,j,prev,base,rc
    real rowvector c,lb,ub,sol,w
    real matrix A,row,Gphi
    real colvector b,gap
    class LinearProgram scalar S
    N=rows(Phi); q=cols(Phi); m=cols(idx); nv=2*q+N; base=round(q^(1/m))
    c=(J(1,q,1),J(1,q,1),J(1,N,C)); lb=J(1,nv,0); ub=J(1,nv,.)
    A=(-Phi,Phi,J(N,N,0)); b=-y
    A=A\(Phi,-Phi,-I(N)); b=b\(y:+eps)

    // Grid-cell cumulative values and monotonic adjacent-cell restrictions.
    Gphi=J(q,q,0)
    for (g=1;g<=q;g++) {
        for (i=1;i<=q;i++) Gphi[g,i]=all(idx[i,.]:<=idx[g,.])
    }
    row=J(1,nv,0); row[1..q]=-Gphi[1,.]; row[(q+1)..(2*q)]=Gphi[1,.]
    A=A\row; b=b\0
    for (g=1;g<=q;g++) for (j=1;j<=m;j++) if (idx[g,j]>0) {
        prev=g-base^(m-j)
        row=J(1,nv,0)
        row[1..q]=Gphi[prev,.]-Gphi[g,.]
        row[(q+1)..(2*q)]=-row[1..q]
        A=A\row; b=b\0
    }
    S=LinearProgram(); S.setMaxOrMin("min"); S.setCoefficients(c)
    S.setInequality(A,b); S.setBounds(lb,ub); S.setTrace("off"); rc=S.optimize()
    if (!S.converged() | S.errorcode()!=0) {
        conv_out=0; w_out=J(q,1,.); F_out=J(N,1,.); obj_out=.; viol_out=.; return
    }
    sol=S.parameters(); w=sol[1..q]-sol[(q+1)..(2*q)]
    w_out=w'; F_out=Phi*w'; gap=b-A*sol'
    obj_out=S.value(); viol_out=max((0,-min(gap))); conv_out=(viol_out<=1e-7)
}

// -----------------------------------------------------------------------
// Main SVF function
// -----------------------------------------------------------------------
struct SVFResult scalar svfeat_main(
    real matrix X,
    real matrix Y,
    real scalar d,
    real scalar C,
    real scalar eps
    )
{
    struct SVFResult scalar res
    real scalar N, nX, nY, r, k, n_failed, n_svec_r, conv, obj,viol
    real scalar phi_k, phi_kr, maxviol,objsum,j
    real matrix Phi,F,ranges,idx,W,Gphi,gridX
    real colvector w_r,y_r,F_r

    N = rows(X); nX = cols(X); nY = cols(Y)

    ranges=J(nX,2,.)
    for (j=1;j<=nX;j++) {
        ranges[j,1]=min(X[.,j])
        ranges[j,2]=max(X[.,j])
    }
    idx=svfeat_grid_index(nX,d)
    Phi=svfeat_grid_phi(X,d,ranges,idx)
    gridX=J(rows(idx),nX,.)
    for (j=1;j<=nX;j++) gridX[.,j]=ranges[j,1]:+idx[.,j]:*((ranges[j,2]-ranges[j,1])/d)

    F         = J(N, nY, .)
    W         = J(cols(Phi),nY,.)
    n_failed  = 0
    n_svec_r  = 0
    maxviol=0; objsum=0

    for (r=1; r<=nY; r++) {
        y_r = Y[.,r]
        svfeat_model9_solve(Phi,y_r,idx,C,eps,w_r,F_r,obj,viol,conv)
        if (!conv) {
            n_failed++
            F[.,r] = J(N,1,.)
        } else {
            W[.,r]=w_r; F[.,r]=F_r
            n_svec_r=n_svec_r+sum(abs(w_r):>1e-8)
            maxviol=max((maxviol,viol)); objsum=objsum+obj
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
    res.n_failed = n_failed
    res.n_svec   = n_svec_r
    res.W        = W
    Gphi=J(rows(idx),rows(idx),0)
    for (k=1;k<=rows(idx);k++) for (j=1;j<=rows(idx);j++) Gphi[k,j]=all(idx[j,.]:<=idx[k,.])
    res.grid_X=gridX
    res.grid_F=Gphi*W
    res.max_violation=(n_failed ? . : maxviol)
    res.objective=(n_failed ? . : objsum)
    res.cv_results=J(0,4,.)
    res.cv_fold_results=J(0,5,.)
    res.cv_predictions=J(0,7,.)
    res.fold_id=J(N,1,.)
    res.selected_C=C
    res.selected_eps=eps
    res.selected_d=d
    res.n_folds=0
    res.tuned=0
    res.cv_failures=0
    return(res)
}

struct SVFResult scalar svfeat_cv_main(
    real matrix X, real matrix Y,
    real rowvector Cgrid, real rowvector Egrid, real rowvector Dgrid,
    real scalar Kfold, real scalar seed)
{
    struct SVFResult scalar res
    real scalar N,nX,nY,ic,ie,id,f,r,j,rownum,mse,err,nval,failed
    real scalar obj,viol,conv,bestrow,bestmse,Cval,Eval,dval
    real colvector foldid,ord,tr,va,w,Ftr,pred
    real matrix ranges,idx,Ptr,Pva,cv,cvfold,cvpred

    N=rows(X); nX=cols(X); nY=cols(Y)
    if (Kfold>N) Kfold=N
    foldid=J(N,1,.)
    ord=order(runiform(N,1),1)
    for (j=1;j<=N;j++) foldid[ord[j]]=mod(j-1,Kfold)+1

    cv=J(cols(Cgrid)*cols(Egrid)*cols(Dgrid),4,.)
    cvfold=J(0,5,.); cvpred=J(0,7,.); rownum=0
    for (id=1;id<=cols(Dgrid);id++) {
        dval=Dgrid[id]
        idx=svfeat_grid_index(nX,dval)
        for (ic=1;ic<=cols(Cgrid);ic++) {
            Cval=Cgrid[ic]
            for (ie=1;ie<=cols(Egrid);ie++) {
                Eval=Egrid[ie]; rownum++; mse=0; failed=0
                for (f=1;f<=Kfold;f++) {
                    tr=selectindex(foldid:!=f); va=selectindex(foldid:==f)
                    ranges=J(nX,2,.)
                    for (j=1;j<=nX;j++) {
                        ranges[j,1]=min(X[tr,j]); ranges[j,2]=max(X[tr,j])
                    }
                    Ptr=svfeat_grid_phi(X[tr,.],dval,ranges,idx)
                    Pva=svfeat_grid_phi(X[va,.],dval,ranges,idx)
                    err=0; nval=rows(va)*nY
                    for (r=1;r<=nY;r++) {
                        svfeat_model9_solve(Ptr,Y[tr,r],idx,Cval,Eval,
                                            w,Ftr,obj,viol,conv)
                        if (!conv) {
                            failed=1
                            break
                        }
                        pred=Pva*w
                        for (j=1;j<=rows(va);j++) {
                            cvpred=cvpred\(Cval,Eval,dval,f,va[j],r,pred[j])
                        }
                        err=err+sum((Y[va,r]-pred):^2)
                    }
                    if (failed) cvfold=cvfold\(Cval,Eval,dval,f,.)
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
        res=svfeat_main(X,Y,Dgrid[1],Cgrid[1],Egrid[1])
        res.n_failed=res.n_failed+1
        res.cv_results=cv; res.cv_fold_results=cvfold
        res.cv_predictions=cvpred; res.fold_id=foldid
        res.n_folds=Kfold; res.tuned=1
        res.cv_failures=sum(cv[.,4]:>=.)
        return(res)
    }
    res=svfeat_main(X,Y,cv[bestrow,3],cv[bestrow,1],cv[bestrow,2])
    res.cv_results=cv; res.cv_fold_results=cvfold
    res.cv_predictions=cvpred; res.fold_id=foldid
    res.selected_C=cv[bestrow,1]
    res.selected_eps=cv[bestrow,2]
    res.selected_d=cv[bestrow,3]
    res.n_folds=Kfold; res.tuned=1
    res.cv_failures=sum(cv[.,4]:>=.)
    return(res)
}

real scalar svfeat_valid_mean(real colvector v)
{
    real colvector w
    w = select(v, v :!= .)
    return(length(w) > 0 ? mean(w) : .)
}

void svfeat_display(struct SVFResult scalar res,
                    string scalar xnames, string scalar ynames,
                    real scalar d, real scalar C, real scalar eps)
{
    real scalar N
    real colvector v, eff1
    N = rows(res.phi)
    v = select(res.phi, res.phi :!= .)

    printf("{txt}\n")
    printf("{txt}Support Vector Frontier (Valero-Carreras, Aparicio, Guerrero 2021, Omega)\n")
    printf("{hline 60}\n")
    printf("{txt}  Variable inputs (x)     : {res}%s{txt}\n", xnames)
    printf("{txt}  Good outputs (y)        : {res}%s{txt}\n", ynames)
    printf("{txt}  Observations            : {res}%g{txt}\n", N)
    printf("{txt}  Parameters: d=%g  C=%g  eps=%g\n",
           res.selected_d,res.selected_C,res.selected_eps)
    if (res.tuned) printf("{txt}  Selection: %g-fold cross-validation\n",res.n_folds)
    printf("{txt}  Nonzero grid weights     : {res}%g{txt}\n", res.n_svec)
    if (res.n_failed > 0) {
        printf("{err}  LP failures             : %g{txt}\n", res.n_failed)
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
