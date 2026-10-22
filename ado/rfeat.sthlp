{smcl}
{* Help revision: 2026-10-22.1; command version unchanged.}{...}
{* *! rfeat  v0.1.3  08sep2026}{...}
{viewerjumpto "Syntax" "rfeat##syntax"}{...}
{viewerjumpto "Description" "rfeat##description"}{...}
{viewerjumpto "Options" "rfeat##options"}{...}
{viewerjumpto "Stored results" "rfeat##results"}{...}
{viewerjumpto "Examples" "rfeat##examples"}{...}
{viewerjumpto "References" "rfeat##references"}{...}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col:{cmd:rfeat} {hline 2}}Random Forest Efficiency Analysis Trees (RFEAT){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:rfeat} {ifin}{cmd:,}
{cmdab:x(}{it:varlist}{cmd:)}
{cmdab:y(}{it:varlist}{cmd:)}
[{it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt x(varlist)}}nonnegative variable inputs{p_end}
{synopt:{opt y(varlist)}}nonnegative good outputs; each observation must have at least one positive output{p_end}
{syntab:Optional}
{synopt:{opt gen:erate(name)}}store mean efficiency phi in dataset (typically >= 1; values below 1 indicate super-efficiency w.r.t. bootstrap trees){p_end}
{synopt:{opt gensd(name)}}store standard deviation of phi across trees{p_end}
{synopt:{opt ntrees(#)}}number of bootstrap EAT trees (default 200){p_end}
{synopt:{opt numstop(#)}}minimum observations in a node before a split is attempted (default 5){p_end}
{synopt:{opt s(#)}}features considered at each split; 0 = all (default 0){p_end}
{synopt:{opt maxdepth(#)}}maximum EAT tree depth; 0 = unlimited (default 0){p_end}
{synopt:{opt maxleaves(#)}}maximum EAT leaves; 0 = unlimited (default 0){p_end}
{synopt:{opt seed(#)}}random seed (default 12345){p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:rfeat} estimates a production frontier from an ensemble of B bootstrap EAT
trees (Random Forest EAT). For each tree, a bootstrap sample is drawn with
replacement and an EAT is fitted to that sample. The fitted output vector is
predicted for every estimation observation and averaged component by component
over all trees. The radial score is then computed from this mean output vector.

{pstd}
Efficiency is output-oriented: phi_i = 1 indicates an efficient DMU;
phi_i > 1 indicates the proportional output expansion needed to reach the frontier.
Because the trees are grown on bootstrap samples, highly productive DMUs may lie
above the averaged bootstrap frontier, yielding phi_i < 1 (super-efficiency).
This is expected under the published construction.

{pstd}
When {opt s(#)} is set to a value greater than 0, only a random subset of s
input features is considered at each candidate split, giving a proper random
forest. Setting s(0) (the default) uses all features at every split, which
corresponds to bagging.

{pstd}
RFEAT aggregates EAT frontiers while preserving monotonicity
and the output-orientation property of the underlying EAT trees.


{marker options}{...}
{title:Options}

{phang}
{opt x(varlist)} specifies the input variables. All must be non-negative.

{phang}
{opt y(varlist)} specifies the output variables. All must be non-negative.

{phang}
{opt generate(name)} saves the mean output-oriented efficiency score phi
(typically >= 1; values below 1 are possible for super-efficient DMUs),
computed from the componentwise average of all fitted-tree predictions, as a
new variable named {it:name}.

{phang}
{opt gensd(name)} saves the standard deviation of phi across trees as a new
variable named {it:name}. It uses the sample denominator B-1 and centers on
the mean of the tree-specific radial scores, not on the forest's final radial
score. It is missing when {cmd:ntrees(1)} is requested. Larger values indicate greater dispersion among the
tree-specific radial scores for a given observation; this is not a sampling
standard error.

{phang}
{opt ntrees(#)} sets the number of bootstrap EAT trees to grow. Default is 200.
Larger values reduce variance at the cost of computation time.

{phang}
{opt numstop(#)} sets the minimum number of observations required in a node
before a split is attempted. Default is 5.

{phang}
{opt s(#)} sets the number of input features randomly selected as candidates at
each split node. Default is 0, which means all features are considered (bagging).
A common choice for a proper random forest is floor(sqrt(p)) where p is the total
number of inputs.

{pstd}
Candidate features are drawn from all input columns, in random draw order;
constant candidates offer no split. Split losses retain full precision and
exact ties keep the first encountered candidate. The R {cmd:eat} package
uses different candidate ordering/filtering and rounds node losses, so even
matched samples do not guarantee identical native-R trees.

{phang}
{opt maxdepth(#)} sets the maximum depth of each EAT tree. 0 (default) imposes
no limit.

{phang}
{opt maxleaves(#)} sets the maximum number of leaves per EAT tree. 0 (default)
imposes no limit.

{phang}
{opt seed(#)} sets the random seed for bootstrap sampling. Default is 12345.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:rfeat} stores the following in {cmd:e()}:

{synoptset 18 tabbed}{...}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sample)}}estimation-sample indicator{p_end}
{synopt:{cmd:e(mean_phi)}}mean output-oriented efficiency across all observations{p_end}
{synopt:{cmd:e(ntrees)}}number of bootstrap trees requested{p_end}
{synopt:{cmd:e(s)}}number of candidate input features per split{p_end}
{synopt:{cmd:e(x)}}input variable names{p_end}
{synopt:{cmd:e(y)}}output variable names{p_end}
{synopt:{cmd:e(cmd)}}{cmd:rfeat}{p_end}
{synopt:{cmd:e(cmdline)}}executed sub-command invocation{p_end}
{synopt:{cmd:e(frontier)}}componentwise mean forest predictions at the estimation observations{p_end}


{marker examples}{...}
{title:Examples}

{pstd}Run the complete example block in the Do-file Editor, omitting leading dots.
The examples create their own data; {cmd:preserve} and {cmd:restore} retain your
original dataset on successful completion. Settings and estimation results may change.
Small samples and tuning grids demonstrate syntax rather than a research design.{p_end}
{* DEA_ML_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}

{pstd}Create reproducible data supporting all examples below:{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 80}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. forvalues j = 1/4 {c -(}}{p_end}
{phang2}{cmd:.     generate double x`j' = 1 + 2*runiform()}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. generate double y = (x1*x2)^0.3*exp(-0.2*runiform())}{p_end}
{phang2}{cmd:. generate double y1 = y}{p_end}
{phang2}{cmd:. generate double y2 = (x3*x4)^0.3*exp(-0.2*runiform())}{p_end}

{pstd}Basic RFEAT with default settings:{p_end}
{phang2}{cmd:. rfeat, x(x1 x2) y(y) gen(phi_basic)}{p_end}

{pstd}RFEAT with 200 trees, fixed seed, and standard deviation stored:{p_end}
{phang2}{cmd:. rfeat, x(x1 x2) y(y) gen(phi_200) gensd(phi_sd) ntrees(200) seed(42)}{p_end}

{pstd}Random forest (not bagging) using sqrt(p) features per split:{p_end}
{phang2}{cmd:. rfeat, x(x1 x2 x3 x4) y(y1 y2) gen(phi_rf) s(2) ntrees(500) seed(1)}{p_end}

{pstd}Via unified dispatcher:{p_end}
{phang2}{cmd:. dea_ml, method(rfeat) x(x1 x2) y(y) gen(phi_dispatch) ntrees(200) seed(99)}{p_end}


{phang2}{cmd:. restore}{p_end}
{* DEA_ML_HELP_EXAMPLES_END}{...}


{marker references}{...}
{title:References}

{phang}
Esteve, M., J. Aparicio, J. J. Rodriguez-Sala, and J. Zhu. 2023.
Random Forests and the measurement of super-efficiency in the context of Free
Disposal Hull.
{it:European Journal of Operational Research} 304(2): 729{c -}744.
{browse "https://doi.org/10.1016/j.ejor.2022.04.024"}

{phang}
Breiman, L. 2001.
Random forests.
{it:Machine Learning} 45(1): 5{c -}32.
{browse "https://doi.org/10.1023/A:1010933404324"}

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
