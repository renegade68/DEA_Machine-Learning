*! dea_ml.ado  v0.2.2
*! Unified machine-learning frontier estimation command
*! Author: Choonjoo Lee, Korea National Defense University
*!         Nonsan-si, Republic of Korea; bloom.rampike@gmail.com
*! Method references: See the corresponding standalone command ado/help files.
*!
*! Syntax:
*!   dea_ml, method(name) [method-specific options]
*!
*! method() dispatches to the corresponding sub-command, passing all
*! remaining options through unchanged. Sub-commands are also available
*! as standalone Stata commands (e.g., ceat, beateat).
*!
*! Methods:
*!   eatree    EAT (Efficiency Analysis Trees)
*!   ceat      Convex EAT
*!   rfeat     Random Forest EAT
*!   csvf      Convexified Support Vector Frontier
*!   bpeat     By-production EAT
*!   bpceat    By-production Convex EAT
*!   dyneat    Dynamic EAT
*!   dynceat   Dynamic Convex EAT
*!   beateat   Boosted EAT
*!   svfeat    Support Vector Frontier (equal-width dominated-cell grid)
*!   svfspline SVF with Splines
*!   lsbmafs   Least Squares Boosting with MARS
*!
*! Examples:
*!   dea_ml, method(ceat) x(x1 x2) y(y) gen(phi) rts(vrs)
*!   dea_ml, method(beateat) x(x) y(y) gen(phi) iter(100) lr(0.3)
*!   dea_ml, method(lsbmafs) x(x1 x2) y(y) gen(phi) maxterms(5) iter(100)
*!   dea_ml, method(svfeat) x(x) y(y) gen(phi) d(3) c(1.0)

program define dea_ml, eclass
    version 16.1

    // [*] captures all remaining options (x, y, gen, etc.) into `options'
    syntax [if] [in] , method(string) [FRONTIERGraph FGROptions(string) *]

    // Valid method list -- add new methods here
    local valid "eatree ceat rfeat csvf bpeat bpceat dyneat dynceat"
    local valid "`valid' beateat svfeat svfspline lsbmafs"

    local found 0
    foreach m of local valid {
        if `"`method'"' == `"`m'"' local found 1
    }

    if !`found' {
        di as err `"Unknown method: "`method'""'
        di as err "Valid methods: `valid'"
        di as err "Type {cmd:help dea_ml} for details."
        exit 198
    }

    // Dispatch: sub-command receives all options except method()
    `method' `if' `in', `options'

    if "`frontiergraph'" != "" {
        _dea_ml_frontier_graph, options(`"`fgroptions'"')
    }
end
