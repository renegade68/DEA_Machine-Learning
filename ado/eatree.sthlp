{smcl}
{* Help revision: 2026-10-22.1; command version unchanged.}{...}
{* *! eatree  v0.1.5  08sep2026}{...}
{viewerjumpto "Syntax"         "eatree##syntax"}{...}
{viewerjumpto "Description"    "eatree##description"}{...}
{viewerjumpto "Options"        "eatree##options"}{...}
{viewerjumpto "Examples"       "eatree##examples"}{...}
{viewerjumpto "Stored results" "eatree##results"}{...}
{viewerjumpto "References"     "eatree##references"}{...}
{viewerjumpto "Author"         "eatree##author"}{...}

{title:Title}

{phang}
{bf:eatree} {hline 2} Efficiency analysis trees


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:eatree} {ifin}{cmd:,}
{cmd:x(}{varlist}{cmd:)}
{cmd:y(}{varlist}{cmd:)}
[{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt x(varlist)}}input (predictor) variables; must be numeric and non-negative{p_end}
{synopt:{opt y(varlist)}}output variables; must be numeric and strictly positive{p_end}
{syntab:Tree complexity}
{synopt:{opt numstop(#)}}minimum number of observations in a node for a split to occur;
  default {cmd:numstop(5)}{p_end}
{synopt:{opt fold(#)}}number of cross-validation folds for tree selection;
  default {cmd:fold(5)}; use {cmd:fold(0)} to skip CV and return the fully grown tree{p_end}
{synopt:{opt maxdepth(#)}}maximum depth of the tree (0 = unlimited);
  default {cmd:maxdepth(0)}{p_end}
{synopt:{opt maxleaves(#)}}maximum number of leaf nodes (0 = unlimited);
  default {cmd:maxleaves(0)}{p_end}
{syntab:Output}
{synopt:{opt gen:erate(newvar)}}name of new variable to store BCC-output efficiency scores{p_end}
{synopt:{opt seed(#)}}random seed for fold assignment in cross-validation;
  default {cmd:seed(12345)}{p_end}
{synopt:{opt foldid(varname)}}use a supplied integer fold assignment; primarily for exact replication{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
All variables in {it:varlist} must be numeric.
Observations with any missing value in x or y are automatically excluded.


{marker description}{...}
{title:Description}

{pstd}
{cmd:eatree} estimates an Efficiency Analysis Tree (EAT) frontier and computes
BCC-output efficiency scores for each observation.
EAT combines the CART decision-tree framework (Breiman et al., 1984) with
the isotonicity constraint of variable-returns-to-scale DEA (Banker et al., 1984).

{pstd}
The algorithm (i) grows a full binary tree by recursively partitioning the
input space to minimise the mean squared prediction error, enforcing at each
step the monotonicity constraint that no leaf may predict less output than
any input-dominating leaf; (ii) generates a nested sequence of pruned subtrees
via cost-complexity pruning; and (iii) selects the optimal subtree by
{it:k}-fold cross-validation using the 1-SE rule.

{pstd}
The resulting piecewise-constant frontier satisfies free disposability by
construction and avoids the overfitting that arises when every observation
anchors the frontier (as in FDH) or when the frontier is the convex hull of
all observations (as in DEA).  Efficiency scores are scaled so that a value
of 1 indicates operation on the frontier and values above 1 indicate the
factor by which outputs could be expanded.

{pstd}
The Stata implementation follows the pruning and 1-SE selection sequence in
the R package {bf:eat} v0.1.4. With identical fold assignments, single- and
multiple-output reference checks reproduce the selected tree and agree with
the R package's efficiency scores within its three-decimal reporting precision.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt x(varlist)} specifies the input variables. At least one variable is
required. Inputs must be nonnegative.

{phang}
{opt y(varlist)} specifies the output (response) variables.  At least one
variable is required.  All variables must be strictly positive numeric.

{dlgtab:Tree complexity}

{phang}
{opt numstop(#)} sets the minimum number of observations in a node
required for a split to take place.  Smaller values allow more splits and
produce larger trees; larger values restrict growth.  Default is 5.

{phang}
{opt fold(#)} sets the number of folds for cross-validation.
Setting {cmd:fold(0)} skips cross-validation entirely and returns
the fully grown tree, which may overfit.  Default is 5.

{phang}
{opt foldid(varname)} supplies an integer fold identifier from 1 through
{opt fold()}. This option makes an externally defined cross-validation
partition exactly reproducible. If omitted, the command randomly permutes the
observations and forms approximately equal-sized folds using {opt seed()}.

{phang}
{opt maxdepth(#)} limits the maximum number of edges on any path from
the root to a leaf.  Setting {cmd:maxdepth(0)} imposes no limit.
Default is 0.

{phang}
{opt maxleaves(#)} limits the maximum number of leaf nodes in the selected
tree.  Setting {cmd:maxleaves(0)} imposes no limit.  Default is 0.

{dlgtab:Output}

{phang}
{opt generate(newvar)} stores the BCC-output efficiency score for each
observation in {it:newvar}.  Excluded observations (missing values in x or y)
receive a missing value.

{phang}
{opt seed(#)} sets the random seed before assigning observations to
cross-validation folds.  Setting the same seed across runs guarantees
reproducible fold assignments.  Default is 12345.


{marker examples}{...}
{title:Examples}

{pstd}Run the complete example block in the Do-file Editor, omitting leading dots.
The examples create their own data; {cmd:preserve} and {cmd:restore} retain your
original dataset on successful completion. Settings and estimation results may change.
Small samples and tuning grids demonstrate syntax rather than a research design.{p_end}
{* DEA_ML_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}

{pstd}Create reproducible one-input, one-output data:{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 80}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. generate double x = 1 + 9*runiform()}{p_end}
{phang2}{cmd:. generate double y = x^0.6*exp(-0.3*runiform())}{p_end}

{pstd}Basic usage with default options{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. eatree, x(x) y(y) numstop(5) fold(5) gen(eff)}{p_end}

{pstd}Summarize efficiency scores{p_end}
{phang2}{cmd:. sum eff, detail}{p_end}

{pstd}Multiple inputs and outputs{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 80}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. forvalues j = 1/3 {c -(}}{p_end}
{phang2}{cmd:.     generate double x`j' = 1 + 2*runiform()}{p_end}
{phang2}{cmd:. {c )-}}{p_end}
{phang2}{cmd:. generate double y1 = (x1*x2)^0.3*exp(-0.2*runiform())}{p_end}
{phang2}{cmd:. generate double y2 = (x2*x3)^0.3*exp(-0.2*runiform())}{p_end}
{phang2}{cmd:. generate double y3 = (x1*x3)^0.3*exp(-0.2*runiform())}{p_end}
{phang2}{cmd:. eatree, x(x1 x2 x3) y(y1 y2 y3) numstop(5) fold(5) gen(eff3)}{p_end}

{pstd}Full tree without cross-validation{p_end}
{phang2}{cmd:. eatree, x(x1) y(y1) numstop(5) fold(0) gen(eff_full)}{p_end}

{pstd}Restrict tree depth to 2 levels{p_end}
{phang2}{cmd:. eatree, x(x1) y(y1) numstop(5) fold(5) maxdepth(2) gen(eff_d2)}{p_end}

{pstd}Longitudinal use: create panel-shaped data and run EAT separately by period{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 200}{p_end}
{phang2}{cmd:. generate int year = 2001 + mod(_n-1,10)}{p_end}
{phang2}{cmd:. generate double capital = 5 + 5*runiform()}{p_end}
{phang2}{cmd:. generate double labor = 2 + 3*runiform()}{p_end}
{phang2}{cmd:. generate double output = capital^0.4*labor^0.3*exp(-0.2*runiform())}{p_end}
{phang2}{cmd:. forvalues yr = 2001/2010 {c -(}}{p_end}
{phang2}{cmd:.     quietly eatree if year == `yr', x(capital labor) y(output) numstop(2) fold(0) gen(eff_`yr')}{p_end}
{phang2}{cmd:. {c )-}}{p_end}


{phang2}{cmd:. restore}{p_end}
{* DEA_ML_HELP_EXAMPLES_END}{...}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:eatree} stores the following in {cmd:e()}:

{synoptset 18 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:e(N)}}number of observations used{p_end}
{synopt:{cmd:e(sample)}}estimation-sample indicator{p_end}
{synopt:{cmd:e(mean_phi)}}mean output-oriented efficiency score{p_end}
{synopt:{cmd:e(n_leaves)}}number of leaf nodes in the selected tree{p_end}
{synopt:{cmd:e(n_nodes)}}total number of nodes (internal + leaf) in the selected tree{p_end}

{syntab:Macros}
{synopt:{cmd:e(cmd)}}{cmd:eatree}{p_end}
{synopt:{cmd:e(cmdline)}}executed sub-command invocation{p_end}
{synopt:{cmd:e(x)}}names of input variables{p_end}
{synopt:{cmd:e(y)}}names of output variables{p_end}
{syntab:Matrices}
{synopt:{cmd:e(nodes)}}selected tree node database{p_end}
{synopt:{cmd:e(frontier)}}radial fitted-output targets at the estimation observations{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{phang}
Esteve, M., J. Aparicio, A. Rabasa, and J. J. Rodr{c i'}guez-Sala. 2020.
Efficiency analysis trees: A new methodology for estimating production
frontiers through decision trees.
{it:Expert Systems with Applications} 162: 113783.
{browse "https://doi.org/10.1016/j.eswa.2020.113783"}

{phang}
Aparicio, J., M. Esteve, J. J. Rodr{c i'}guez-Sala, and J. L. Zofio. 2021.
The estimation of productive efficiency through machine learning techniques:
Efficiency Analysis Trees.
In {it:Data-Enabled Analytics: DEA for Big Data}, 51{c -}92.
{browse "https://doi.org/10.1007/978-3-030-75162-3_3"}

{phang}
Breiman, L., J. H. Friedman, R. A. Olshen, and C. J. Stone. 1984.
{it:Classification and Regression Trees}. Belmont, CA: Wadsworth.

{phang}
Charnes, A., W. W. Cooper, and E. Rhodes. 1978.
Measuring the efficiency of decision making units.
{it:European Journal of Operational Research} 2(6): 429{c -}444.
{browse "https://doi.org/10.1016/0377-2217(78)90138-8"}

{phang}
Banker, R. D., A. Charnes, and W. W. Cooper. 1984.
Some models for estimating technical and scale inefficiencies in data
envelopment analysis.
{it:Management Science} 30(9): 1078{c -}1092.
{browse "https://doi.org/10.1287/mnsc.30.9.1078"}


{marker author}{...}
{title:Author}

{pstd}
Choonjoo Lee{break}
Korea National Defense University{break}
Nonsan-si, Republic of Korea{break}
bloom.rampike@gmail.com{p_end}

{title:Publication and license}

{pstd}Lee, C. Forthcoming. {it:dea_ml: Machine-learning frontier estimators for efficiency analysis in Stata}.
Accepted for publication in the {browse "https://www.stata-journal.com/":Stata Journal}.{p_end}
{pstd}Author-maintained software under the Academic Free License 3.0 (AFL-3.0).
See LICENSE and LICENSE_CODE.txt in the repository. The acceptance notice refers
to the article; later software revisions are maintained separately.{p_end}
