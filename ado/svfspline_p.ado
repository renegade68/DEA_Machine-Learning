*! svfspline_p.ado  v0.2.1
program define svfspline_p
    version 16.1
    if "`e(cmd)'" != "svfspline" {
        error 301
    }
    syntax newvarname [if] [in] [, OUTput(integer 1)]
    if `output' < 1 | `output' > e(n_outputs) {
        di as err "output() must be between 1 and " e(n_outputs)
        exit 198
    }

    tempvar touse
    quietly mark `touse' `if' `in'
    local xvars `e(x)'
    foreach v of local xvars {
        confirm numeric variable `v'
        markout `touse' `v'
        quietly count if `touse' & `v' < 0
        if r(N) {
            di as err "input variable `v' contains negative values in the prediction sample"
            exit 459
        }
    }

    quietly count if `touse'
    if r(N) == 0 {
        gen `typlist' `varlist' = .
        exit
    }

    tempname Xm Bmat Rmat Pmat
    mkmat `xvars' if `touse', matrix(`Xm')
    matrix `Bmat' = e(coefficients)
    matrix `Rmat' = e(xranges)
    mata: st_matrix("`Pmat'", svfspline_p_design_at( ///
        st_matrix("`Xm'"), st_numscalar("e(nknots)"), ///
        st_matrix("`Rmat'")) * st_matrix("`Bmat'"))
    gen `typlist' `varlist' = .
    mata: st_store(., "`varlist'", "`touse'", st_matrix("`Pmat'")[.,`output'])
    label variable `varlist' "SVF-Spline frontier prediction, output `output'"
end

mata:
real matrix svfspline_p_basis1_at(real colvector x, real scalar d,
                                  real scalar xmin, real scalar xmax)
{
    real scalar N,k
    real matrix B
    real colvector knots,col_k
    N=rows(x)
    knots=J(d+1,1,xmin)
    for (k=1;k<=d+1;k++) knots[k]=xmin+(k-1)*(xmax-xmin)/d
    B=J(N,d+3,0)
    B[.,1]=J(N,1,1)
    B[.,2]=x
    for (k=1;k<=d+1;k++) {
        col_k=x:-knots[k]
        B[.,k+2]=col_k:*(col_k:>0)
    }
    return(B)
}

real matrix svfspline_p_design_at(real matrix X, real scalar d,
                                   real matrix ranges)
{
    real scalar j
    real matrix Phi
    Phi=J(rows(X),0,0)
    for (j=1;j<=cols(X);j++) {
        Phi=(Phi,svfspline_p_basis1_at(X[.,j],d,ranges[j,1],ranges[j,2]))
    }
    return(Phi)
}
end
