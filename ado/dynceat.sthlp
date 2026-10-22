{smcl}
{* Help revision: 2026-10-22.1; command version unchanged.}{...}
{* *! dynceat  v0.2.1  08sep2026}{...}
{viewerjumpto "Syntax" "dynceat##syntax"}{...}
{viewerjumpto "Description" "dynceat##description"}{...}
{viewerjumpto "Options" "dynceat##options"}{...}
{viewerjumpto "Stored results" "dynceat##results"}{...}
{viewerjumpto "Examples" "dynceat##examples"}{...}
{viewerjumpto "References" "dynceat##references"}{...}

{title:Title}
{p2colset 5 18 20 2}{...}
{p2col:{cmd:dynceat} {hline 2}}Dynamic Convex EAT (Dynamic CEAT){p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}
{p 8 16 2}
{cmd:dynceat} {ifin}{cmd:,} {opt x(varlist)} {opt y(varlist)} {opt inv(varlist)} {opt kfix(varlist)} {opt dep(numlist)}
[{opt gen(name)} {opt numstop(#)} {opt fold(#)} {opt seed(#)}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt x(varlist)}}nonnegative variable inputs{p_end}
{synopt:{opt y(varlist)}}nonnegative good outputs; at least one positive per observation{p_end}
{synopt:{opt inv(varlist)}}investment variables{p_end}
{synopt:{opt kfix(varlist)}}nonnegative fixed capital stock variables{p_end}
{synopt:{opt dep(numlist)}}one depreciation rate in [0,1] for each variable in {cmd:kfix()}, e.g. {cmd:dep(0.15)}{p_end}
{syntab:Optional}
{synopt:{opt gen(name)}}new variable to store the DDF efficiency score beta (>= 0; 0 = efficient){p_end}
{synopt:{opt numstop(#)}}minimum observations in a node before a split is attempted; default {cmd:numstop(5)}{p_end}
{synopt:{opt fold(#)}}number of CV folds for pruning; default {cmd:fold(5)}{p_end}
{synopt:{opt seed(#)}}random seed for CV fold assignment; unused with {cmd:fold(0)}{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:dynceat} estimates dynamic productive efficiency using the Convex Efficiency
Analysis Tree (CEAT) under the directional distance function (DDF) framework. It
follows the same two-stage structure as {helpb dyneat}: a single EAT tree is grown
and pruned with cross-validation. Following Aparicio, Esteve, and Kapelko
(2023, model 11), net-investment responses are constructed as
Z_comp = inv - delta * kfix using the rate vector in {opt dep()}; one rate is
required for each paired {opt inv()} and {opt kfix()} variable.{p_end}

{pstd}
The key distinction from {helpb dyneat} is that {cmd:dynceat} solves a linear program
(LP) over convex combinations of leaf nodes (requiring Stata 16+ LinearProgram class)
rather than taking the best single leaf. Because CEAT permits convex
combinations of the EAT leaf references, the resulting efficiency score is
weakly greater than or equal to the DynEAT score for observations with a
defined direction and a successfully solved LP, up to numerical tolerance.{p_end}

{pstd}
Direction vectors are identical to {helpb dyneat}:{p_end}
{phang2}gx = x0 - min(x){p_end}
{phang2}gy = max(y) - y0{p_end}
{phang2}gI_i = 0.2 * [max(k) - k_i]{p_end}

{pstd}
A value of beta = 0 indicates full efficiency under the convex frontier; larger
values indicate greater inefficiency.{p_end}

{marker options}{...}
{title:Options}

{phang}{opt x(varlist)} specifies the variable input variables. Required.{p_end}

{phang}{opt y(varlist)} specifies the good (desirable) output variables. Required.{p_end}

{phang}{opt inv(varlist)} specifies the investment variables used to form compound
outputs. Required.{p_end}

{phang}{opt kfix(varlist)} specifies the fixed capital stock variables used to form
compound outputs. Required.{p_end}

{phang}{opt dep(numlist)} specifies one depreciation rate in [0,1] for each
paired investment and fixed-capital variable. For example, {cmd:dep(0.15)}
sets the rate for one pair to 0.15. Required.{p_end}

{phang}{opt gen(name)} stores the DDF efficiency score (beta >= 0) in the new
variable {it:name}. If omitted, results are stored in {cmd:e()} only.{p_end}

{phang}{opt numstop(#)} sets the minimum number of observations in a node before
a split is attempted. Default is 5.{p_end}

{phang}{opt fold(#)} sets the number of folds for cross-validated cost-complexity
pruning. Default is 5.{p_end}

{phang}{opt seed(#)} sets the random seed for CV fold assignment. It is unused
with {cmd:fold(0)}. Default is 12345.{p_end}

{pstd}If all components of an observation's translation-invariant direction are
zero, its DDF is undefined. The score is missing and the observation is counted
in {cmd:e(n_zero_direction)}, separately from LP failures.{p_end}

{pstd}Net investment may be negative. The embedded tree uses response-specific
sample minima as initial lower bounds, so negative responses are not raised
to zero and translations of net investment are preserved.{p_end}

{marker results}{...}
{title:Stored results}

{pstd}{cmd:dynceat} stores the following in {cmd:e()}:{p_end}
{synoptset 18 tabbed}{...}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sample)}}estimation-sample indicator{p_end}
{synopt:{cmd:e(mean_beta)}}mean DDF efficiency score across all DMUs{p_end}
{synopt:{cmd:e(n_leaves)}}number of leaves in the pruned tree{p_end}
{synopt:{cmd:e(dep)}}depreciation rate when one quasi-fixed input is used{p_end}
{synopt:{cmd:e(dep_rates)}}row vector of all depreciation rates{p_end}
{synopt:{cmd:e(n_nodes)}}total number of nodes in the selected tree{p_end}
{synopt:{cmd:e(n_failed)}}number of convexification LP failures{p_end}
{synopt:{cmd:e(n_zero_direction)}}observations with an undefined all-zero DDF direction{p_end}
{synopt:{cmd:e(solver)}}LP solver identifier{p_end}
{synopt:{cmd:e(x)}}names of input variables{p_end}
{synopt:{cmd:e(y)}}names of output variables{p_end}
{synopt:{cmd:e(inv)}}names of investment variables{p_end}
{synopt:{cmd:e(kfix)}}names of fixed capital stock variables{p_end}
{synopt:{cmd:e(cmd)}}"dynceat"{p_end}
{synopt:{cmd:e(cmdline)}}executed sub-command invocation{p_end}

{marker examples}{...}
{title:Examples}

{pstd}Run the complete example block in the Do-file Editor, omitting leading dots.
The examples create their own data; {cmd:preserve} and {cmd:restore} retain your
original dataset on successful completion. Settings and estimation results may change.
Small samples and tuning grids demonstrate syntax rather than a research design.{p_end}
{* DEA_ML_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}

{pstd}Create reproducible dynamic-production data:{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 100}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. generate double x = 1 + 4*runiform()}{p_end}
{phang2}{cmd:. generate double y = x^0.6*exp(-0.2*runiform())}{p_end}
{phang2}{cmd:. generate double inv = 0.2 + 0.8*runiform()}{p_end}
{phang2}{cmd:. generate double k = 5 + 5*runiform()}{p_end}

{pstd}Basic dynamic CEAT with one input, one output, one investment, one capital
stock variable, and a 15% depreciation rate:{p_end}
{phang2}{cmd:. dynceat, x(x) y(y) inv(inv) kfix(k) dep(0.15) gen(beta_basic)}{p_end}

{pstd}Compare DynEAT and DynCEAT scores (weak dominance, up to solver tolerance):{p_end}
{phang2}{cmd:. dyneat,  x(x) y(y) inv(inv) kfix(k) dep(0.15) gen(beta_eat)}{p_end}
{phang2}{cmd:. dynceat, x(x) y(y) inv(inv) kfix(k) dep(0.15) gen(beta_ceat)}{p_end}
{phang2}{cmd:. summarize beta_eat beta_ceat}{p_end}

{pstd}Display mean efficiency after estimation:{p_end}
{phang2}{cmd:. dynceat, x(x) y(y) inv(inv) kfix(k) dep(0.15) gen(beta_display)}{p_end}
{phang2}{cmd:. display e(mean_beta)}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEA_ML_HELP_EXAMPLES_END}{...}


{marker references}{...}
{title:References}

{phang}Aparicio, J., M. Esteve, and M. Kapelko. 2023. Measuring dynamic
inefficiency through machine learning techniques.
{it:Expert Systems with Applications} 228: 120417.
{browse "https://doi.org/10.1016/j.eswa.2023.120417"}
{p_end}

{marker author}{...}
{title:Author}

{pstd}Choonjoo Lee{break}
Korea National Defense University{break}
Nonsan-si, Republic of Korea{break}
bloom.rampike@gmail.com{p_end}

{title:Publication and license}

{pstd}Lee, C. Forthcoming. {it:dea_ml: Machine-learning frontier estimators for efficiency analysis in Stata}.
Accepted for publication in the {browse "https://www.stata-journal.com/":Stata Journal}.{p_end}
{pstd}Author-maintained software under the Academic Free License 3.0 (AFL-3.0).
See LICENSE and LICENSE_CODE.txt in the repository. The acceptance notice refers
to the article; later software revisions are maintained separately.{p_end}
