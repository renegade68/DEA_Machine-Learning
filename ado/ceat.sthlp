{smcl}
{* Help revision: 2026-10-22.1; command version unchanged.}{...}
{* *! ceat  v0.2.1  08sep2026}{...}
{viewerjumpto "Syntax" "ceat##syntax"}{...}
{viewerjumpto "Description" "ceat##description"}{...}
{viewerjumpto "Options" "ceat##options"}{...}
{viewerjumpto "Stored results" "ceat##results"}{...}
{viewerjumpto "Examples" "ceat##examples"}{...}
{viewerjumpto "References" "ceat##references"}{...}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col:{cmd:ceat} {hline 2}}Convex Efficiency Analysis Trees (CEAT){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:ceat} {ifin}{cmd:,}
{cmdab:x(}{it:varlist}{cmd:)}
{cmdab:y(}{it:varlist}{cmd:)}
[{it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt x(varlist)}}nonnegative variable inputs{p_end}
{synopt:{opt y(varlist)}}strictly positive good outputs{p_end}
{syntab:Optional}
{synopt:{opt gen:erate(name)}}store output-oriented efficiency phi (>= 1) in dataset{p_end}
{synopt:{opt genlambda(string)}}prefix for lambda weight variables{p_end}
{synopt:{opt rts(string)}}returns to scale: {cmd:vrs} (default), {cmd:crs}, {cmd:drs}, or {cmd:irs}{p_end}
{synopt:{opt numstop(#)}}minimum observations in a node before a split is attempted (default 5){p_end}
{synopt:{opt fold(#)}}cross-validation folds for EAT pruning (default 5){p_end}
{synopt:{opt foldid(varname)}}use a supplied fold assignment (mainly for exact replication){p_end}
{synopt:{opt maxdepth(#)}}maximum EAT tree depth; 0 = unlimited (default 0){p_end}
{synopt:{opt maxleaves(#)}}maximum EAT leaves; 0 = unlimited (default 0){p_end}
{synopt:{opt seed(#)}}random seed for internally generated CV folds; unused with {cmd:fold(0)} or {cmd:foldid()}{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:ceat} estimates a production frontier by solving a linear program (LP)
over the convex combination of leaves of an Efficiency Analysis Tree (EAT).
The CEAT frontier is output-oriented, monotone non-decreasing in inputs, and
satisfies variable returns to scale (VRS) by default.

{pstd}
EAT partitions the input space into rectangular regions using recursive
binary splitting and pruned by cost-complexity cross-validation.
CEAT then applies a linear program to convex combinations of the EAT leaf
reference points, yielding a piecewise-linear convexification of the selected
EAT frontier.

{pstd}
Output-oriented efficiency: phi_i = max{phi : phi*y_i <= F(x_i)}, phi >= 1.
Efficient DMUs have phi = 1.


{marker options}{...}
{title:Options}

{phang}
{opt x(varlist)} specifies the input variables. All must be non-negative.

{phang}
{opt y(varlist)} specifies the output variables. All must be strictly positive.

{phang}
{opt generate(name)} saves the output-oriented efficiency score phi (>= 1) as
a new variable named {it:name}.

{phang}
{opt genlambda(string)} saves the DEA lambda weights (convex combination
coefficients over EAT leaves) as variables {it:string}1, {it:string}2, ...

{phang}
{opt rts(string)} specifies returns to scale. {cmd:vrs} (default) imposes
sum(lambda) = 1; {cmd:crs} drops this constraint; {cmd:drs} imposes
sum(lambda) <= 1; and {cmd:irs} imposes sum(lambda) >= 1.

{phang}
{opt numstop(#)} sets the minimum number of observations required in a node
before a split is attempted. Default is 5.

{phang}
{opt fold(#)} sets the number of cross-validation folds used for EAT
cost-complexity pruning. Default is 5.

{phang}
{opt foldid(varname)} supplies integer fold identifiers from 1 through
{cmd:fold()}. This option is intended primarily for exact replication and
external implementation checks.

{phang}
{opt maxdepth(#)} sets the maximum depth of the EAT tree. 0 (default) imposes
no limit.

{phang}
{opt maxleaves(#)} sets the maximum number of leaves. 0 (default) imposes no
limit.

{phang}
{opt seed(#)} sets the random seed only when CV folds are generated internally.
It is unused with {cmd:fold(0)} or {cmd:foldid()}. Default is 12345.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:ceat} stores the following in {cmd:e()}:

{synoptset 18 tabbed}{...}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sample)}}estimation-sample indicator{p_end}
{synopt:{cmd:e(mean_phi)}}mean output-oriented efficiency{p_end}
{synopt:{cmd:e(n_failed)}}number of LP failures (should be 0){p_end}
{synopt:{cmd:e(solver)}}LP solver identifier{p_end}
{synopt:{cmd:e(n_leaves)}}number of EAT leaves used{p_end}
{synopt:{cmd:e(n_nodes)}}number of EAT nodes used{p_end}
{synopt:{cmd:e(rts)}}returns to scale ({cmd:vrs}, {cmd:crs}, {cmd:drs}, or {cmd:irs}){p_end}
{synopt:{cmd:e(x)}}input variable names{p_end}
{synopt:{cmd:e(y)}}output variable names{p_end}
{synopt:{cmd:e(cmd)}}{cmd:ceat}{p_end}
{synopt:{cmd:e(cmdline)}}executed sub-command invocation{p_end}
{synopt:{cmd:e(nodes)}}selected EAT node database{p_end}
{synopt:{cmd:e(frontier)}}radial fitted-output targets at the estimation observations{p_end}


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

{pstd}Basic 1×1 frontier with VRS:{p_end}
{phang2}{cmd:. ceat, x(x) y(y) gen(phi_vrs)}{p_end}

{pstd}CRS frontier with efficiency stored:{p_end}
{phang2}{cmd:. ceat, x(x1 x2) y(y) gen(phi_crs) rts(crs)}{p_end}

{pstd}Via unified dispatcher:{p_end}
{phang2}{cmd:. dea_ml, method(ceat) x(x) y(y) gen(phi_dispatch) rts(vrs) seed(99)}{p_end}


{phang2}{cmd:. restore}{p_end}
{* DEA_ML_HELP_EXAMPLES_END}{...}


{marker references}{...}
{title:References}

{phang}
Aparicio, J., M. Esteve, W. Rodriguez-Sala, and J. L. Zofio. 2021.
The estimation of productive efficiency through machine learning techniques:
Efficiency Analysis Trees.
In {it:Data-Enabled Analytics: DEA for Big Data}, 51{c -}92.
{browse "https://doi.org/10.1007/978-3-030-75162-3_3"}

{phang}
Banker, R. D., A. Charnes, and W. W. Cooper. 1984.
Some models for estimating technical and scale inefficiencies in
data envelopment analysis.
{it:Management Science} 30: 1078–1092.
{browse "https://doi.org/10.1287/mnsc.30.9.1078"}

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
