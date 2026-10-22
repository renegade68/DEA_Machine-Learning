{smcl}
{* Help revision: 2026-10-22.1; command version unchanged.}{...}
{* *! csvf  v0.1.2  08sep2026}{...}
{viewerjumpto "Syntax" "csvf##syntax"}{...}
{viewerjumpto "Description" "csvf##description"}{...}
{viewerjumpto "Options" "csvf##options"}{...}
{viewerjumpto "Stored results" "csvf##results"}{...}
{viewerjumpto "Examples" "csvf##examples"}{...}
{viewerjumpto "References" "csvf##references"}{...}

{title:Title}
{p2colset 5 18 20 2}{...}
{p2col:{cmd:csvf} {hline 2}}Convexified Support Vector Frontier (CSVF){p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}
{p 8 16 2}{cmd:csvf} {ifin}{cmd:,} {opt x(varlist)} {opt y(varname)}
[{opt gen(name)} {opt d(#)} {opt c(#)} {opt eps(#)} {opt tune}
 {opt cgrid(numlist)} {opt epsgrid(numlist)} {opt dgrid(numlist)}
 {opt folds(#)} {opt seed(#)}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt x(varlist)}}nonnegative input variables{p_end}
{synopt:{opt y(varname)}}one strictly positive output variable{p_end}
{syntab:Optional}
{synopt:{opt gen(name)}}new variable for output-oriented efficiency phi (>=1){p_end}
{synopt:{opt d(#)}}equal-width intervals per input; default 3{p_end}
{synopt:{opt c(#)}}SVF slack penalty; default 1{p_end}
{synopt:{opt eps(#)}}SVF epsilon-insensitivity band; default 0{p_end}
{synopt:{opt tune}}select C, epsilon, and d by cross-validation{p_end}
{synopt:{opt cgrid(numlist)}}candidate C values{p_end}
{synopt:{opt epsgrid(numlist)}}candidate epsilon values{p_end}
{synopt:{opt dgrid(numlist)}}candidate interval counts{p_end}
{synopt:{opt folds(#)}}number of cross-validation folds; default 5{p_end}
{synopt:{opt seed(#)}}random seed for fold assignment; default 12345{p_end}
{synoptline}

{marker description}{...}
{title:Description}
{pstd}{cmd:csvf} implements the convexification in model (10) of
Valero-Carreras, Aparicio, and Guerrero (2021). It first fits their model (9)
Support Vector Frontier (SVF). Each equal-width grid cell supplies a virtual
point consisting of its lower knot and fitted SVF output. A variable-returns-to-
scale DEA-type envelopment of these virtual points yields a concave production
frontier and a convex production possibility set.{p_end}

{pstd}For observation i, the command maximizes phi subject to the virtual-point
convex combination using no more than x_i, producing at least phi*y_i, and
having weights that are nonnegative and sum to one. Thus phi>=1 and larger
values indicate greater output-oriented inefficiency. The 2021 model is a
production-function estimator and therefore {opt y()} accepts one output.{p_end}

{pstd}Requires Stata 16 or later because it uses Mata's LinearProgram class.{p_end}

{marker options}{...}
{title:Options}
{phang}{opt x(varlist)} specifies nonnegative inputs. Required.{p_end}
{phang}{opt y(varname)} specifies one strictly positive output. Required.{p_end}
{phang}{opt gen(name)} stores phi in a new variable.{p_end}
{phang}{opt d(#)}, {opt c(#)}, and {opt eps(#)} specify the underlying SVF
grid, slack penalty, and epsilon band. Their defaults are 3, 1, and 0.{p_end}
{phang}{opt tune}, {opt cgrid()}, {opt epsgrid()}, {opt dgrid()},
{opt folds()}, and {opt seed()} have the same leakage-free tuning definitions
as in {help svfeat}.{p_end}

{pstd}Failed second-stage LP scores are missing and counted in {cmd:e(n_failed)}.
If every score fails, the command exits with error 498 without creating the
requested output variable. Any prior estimation results, including {cmd:e(sample)},
are restored on error; first-stage SVF results are not left as a successful CSVF
estimate. Summary scores exclude missing observations.{p_end}

{marker results}{...}
{title:Stored results}
{pstd}{cmd:csvf} stores the following in {cmd:e()}:{p_end}
{synoptset 20 tabbed}{...}
{synopt:{cmd:e(N)}}number of complete observations{p_end}
{synopt:{cmd:e(sample)}}estimation-sample indicator{p_end}
{synopt:{cmd:e(mean_phi)}}mean output-oriented efficiency score{p_end}
{synopt:{cmd:e(n_failed)}}number of failed envelopment LPs{p_end}
{synopt:{cmd:e(solver)}}LP solver identifier{p_end}
{synopt:{cmd:e(n_grid)}}number of virtual grid points{p_end}
{synopt:{cmd:e(d)}}selected interval count{p_end}
{synopt:{cmd:e(c_pen)}}selected SVF slack penalty{p_end}
{synopt:{cmd:e(eps)}}selected epsilon band{p_end}
{synopt:{cmd:e(tuned)}}1 if cross-validation was requested{p_end}
{synopt:{cmd:e(folds)}}number of folds, or 0 for a fixed fit{p_end}
{synopt:{cmd:e(frontier)}}radial fitted-output targets{p_end}
{synopt:{cmd:e(lambda)}}virtual-point intensity weights by observation{p_end}
{synopt:{cmd:e(grid_inputs)}}virtual input points{p_end}
{synopt:{cmd:e(grid_frontier)}}SVF outputs at virtual points{p_end}
{synopt:{cmd:e(x)}}input variable names{p_end}
{synopt:{cmd:e(y)}}output variable name{p_end}
{synopt:{cmd:e(orientation)}}{cmd:output}{p_end}
{synopt:{cmd:e(returns_to_scale)}}{cmd:vrs}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:csvf}{p_end}
{synopt:{cmd:e(cmdline)}}executed sub-command invocation{p_end}

{marker examples}{...}
{title:Examples}
{pstd}Run the complete example block in the Do-file Editor, omitting leading dots.
The examples create their own data; {cmd:preserve} and {cmd:restore} retain your
original dataset on successful completion. Settings and estimation results may change.
Small samples and tuning grids demonstrate syntax rather than a research design.{p_end}
{* DEA_ML_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}
{pstd}Create reproducible example data:{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 80}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. generate double x1 = 1 + 4*runiform()}{p_end}
{phang2}{cmd:. generate double x2 = 2 + 3*runiform()}{p_end}
{phang2}{cmd:. generate double x = x1}{p_end}
{phang2}{cmd:. generate double y = x1^0.4*x2^0.3*exp(-0.2*runiform())}{p_end}
{phang2}{cmd:. csvf, x(x) y(y) gen(phi_basic)}{p_end}
{phang2}{cmd:. csvf, x(x1 x2) y(y) gen(phi_grid) d(2) c(5) eps(.05)}{p_end}
{phang2}{cmd:. dea_ml, method(csvf) x(x) y(y) gen(phi_dispatch) d(2)}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEA_ML_HELP_EXAMPLES_END}{...}


{marker references}{...}
{title:References}
{phang}Valero-Carreras, D., J. Aparicio, and N. M. Guerrero. 2021.
Support vector frontiers: A new approach for estimating production functions
through support vector machines. {it:Omega} 104: 102490.
{browse "https://doi.org/10.1016/j.omega.2021.102490"}{p_end}

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
