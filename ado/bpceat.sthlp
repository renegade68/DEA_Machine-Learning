{smcl}
{* Help revision: 2026-10-22.1; command version unchanged.}{...}
{* *! bpceat  v0.2.2  08sep2026}{...}
{viewerjumpto "Syntax" "bpceat##syntax"}{...}
{viewerjumpto "Description" "bpceat##description"}{...}
{viewerjumpto "Options" "bpceat##options"}{...}
{viewerjumpto "Stored results" "bpceat##results"}{...}
{viewerjumpto "Examples" "bpceat##examples"}{...}
{viewerjumpto "References" "bpceat##references"}{...}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col:{cmd:bpceat} {hline 2}}By-production Convex Efficiency Analysis Trees (BP-CEAT){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:bpceat} {ifin}{cmd:,}
{cmdab:x1(}{it:varlist}{cmd:)}
{cmdab:x2(}{it:varlist}{cmd:)}
{cmdab:y(}{it:varlist}{cmd:)}
{cmdab:z(}{it:varlist}{cmd:)}
[{it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt x1(varlist)}}nonnegative non-polluting inputs{p_end}
{synopt:{opt x2(varlist)}}nonnegative pollution-generating inputs{p_end}
{synopt:{opt y(varlist)}}nonnegative good outputs; at least one positive per observation{p_end}
{synopt:{opt z(varlist)}}nonnegative bad outputs; at least one positive per observation{p_end}
{syntab:Optional}
{synopt:{opt gen(name)}}store B = delta1*beta_Psi1/(1+beta_Psi1) + delta2*beta_Psi2{p_end}
{synopt:{opt genpsi1(name)}}store phase-1 component beta_Psi1{p_end}
{synopt:{opt genpsi2(name)}}store phase-2 component beta_Psi2{p_end}
{synopt:{opt delta1(#)}}nonnegative weight on beta_Psi1 (default 0.5){p_end}
{synopt:{opt delta2(#)}}nonnegative weight on beta_Psi2; delta1 + delta2 must equal 1 (default 0.5){p_end}
{synopt:{opt numstop(#)}}minimum observations in a node before a split is attempted (default 5){p_end}
{synopt:{opt fold(#)}}cross-validation folds for pruning (default 5){p_end}
{synopt:{opt seed(#)}}random seed for CV fold assignment; unused with {cmd:fold(0)}{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:bpceat} estimates a by-production model using convex frontiers and is the
convex analogue of {help bpeat}. The model separates a good-output technology
from a bad-output technology but replaces the EAT step-function frontiers with
CEAT (Convex EAT) linear-program frontiers that operate over the convex hull of
EAT leaves.

{pstd}
{bf:Psi1 (good-output technology):} A CEAT LP is solved with inputs (X1, X2)
and good outputs Y.
The convex hull of EAT leaves provides a convex, piecewise-linear frontier.
The component score beta_Psi1 measures good-output
inefficiency relative to the CEAT frontier.

{pstd}
{bf:Psi2 (bad-output technology):} A CEAT LP is solved for the bad-output
component using bad outputs Z and pollution-generating inputs X2. The component
score beta_Psi2 measures bad-output inefficiency relative to this convex
frontier.

{pstd}
The combined efficiency score is:

{pstd}
{it:B = delta1 * [beta_Psi1/(1 + beta_Psi1)] + delta2 * beta_Psi2}

{pstd}
The transformation beta_Psi1/(1 + beta_Psi1), equivalently
1 - 1/(1 + beta_Psi1), places the good-output expansion component on
the bounded scale used in the published by-production measure before it is
combined with beta_Psi2.

{pstd}
where delta1 and delta2 are user-specified weights (default 0.5 each). All
scores are bounded below by 0.

{pstd}
Because CEAT permits convex combinations of the EAT leaf references, its
feasible reference set contains the single-leaf EAT references. For matched
fitted leaf reference sets and successfully solved component programs,
bpceat scores are therefore weakly greater than the corresponding bpeat
scores, up to LP numerical tolerance: beta_Psi1(CEAT) >= beta_Psi1(EAT)
and beta_Psi2(CEAT) >= beta_Psi2(EAT).

{pstd}
The CEAT LP is solved internally; there is no need to specify a convexity
parameter explicitly.
The internal LP convergence tolerance is 1e-10. This is a solver stopping
criterion, not a universal bound on score error for ill-conditioned data.


{marker options}{...}
{title:Options}

{phang}
{opt x1(varlist)} specifies non-polluting inputs. These enter the good-output
technology together with X2.

{phang}
{opt x2(varlist)} specifies pollution-generating inputs. These enter both
components of the by-production model.

{phang}
{opt y(varlist)} specifies the final output variables. All must be non-negative.

{phang}
{opt z(varlist)} specifies bad outputs, or undesirable outputs. All must be
non-negative.

{phang}
{opt gen(name)} saves the combined efficiency score B as a new variable named
{it:name}.

{phang}
{opt genpsi1(name)} saves the phase-1 component score beta_Psi1 as a new
variable named {it:name}.

{phang}
{opt genpsi2(name)} saves the phase-2 component score beta_Psi2 as a new
variable named {it:name}.

{phang}
{opt delta1(#)} sets the weight applied to the transformed good-output
component beta_Psi1/(1 + beta_Psi1) in the combined score.
Default is 0.5.

{phang}
{opt delta2(#)} sets the weight applied to beta_Psi2 in the combined score.
Default is 0.5. The published measure requires delta1 + delta2 = 1; the
command rejects weights that do not satisfy this condition.

{phang}
{opt numstop(#)} sets the minimum number of observations required in a node
before a split is attempted in the underlying EAT trees. Default is 5.

{phang}
{opt fold(#)} sets the number of cross-validation folds used for EAT
cost-complexity pruning. Default is 5.

{phang}
{opt seed(#)} sets the random seed for cross-validation. It is unused with
{cmd:fold(0)}. Default is 12345.

{pstd}A failed component LP is stored as missing, never as zero.  The associated
failure counter identifies the affected component; a genuine boundary optimum
remains zero. If no observation has both component LPs solved, the command
exits with error 498 before creating output variables or replacing estimation
results. This includes failures that affect different observations in the two
technologies. Combined-score summaries exclude missing scores.{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:bpceat} stores the following in {cmd:e()}:

{synoptset 18 tabbed}{...}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sample)}}estimation-sample indicator{p_end}
{synopt:{cmd:e(mean_b)}}mean combined efficiency score B{p_end}
{synopt:{cmd:e(mean_beta1)}}mean first-component directional score{p_end}
{synopt:{cmd:e(mean_beta2)}}mean second-component directional score{p_end}
{synopt:{cmd:e(n_leaves1)}}number of leaves in the first component tree{p_end}
{synopt:{cmd:e(n_leaves2)}}number of leaves in the second component tree{p_end}
{synopt:{cmd:e(n_failed1)}}first-component convexification LP failures{p_end}
{synopt:{cmd:e(n_failed2)}}second-component convexification LP failures{p_end}
{synopt:{cmd:e(solver)}}LP solver identifier{p_end}
{synopt:{cmd:e(delta1)}}weight on beta_Psi1{p_end}
{synopt:{cmd:e(delta2)}}weight on beta_Psi2{p_end}
{synopt:{cmd:e(x1)}}fixed input variable names{p_end}
{synopt:{cmd:e(x2)}}discretionary input variable names{p_end}
{synopt:{cmd:e(y)}}output variable names{p_end}
{synopt:{cmd:e(z)}}intermediate output variable names{p_end}
{synopt:{cmd:e(cmd)}}{cmd:bpceat}{p_end}
{synopt:{cmd:e(cmdline)}}executed sub-command invocation{p_end}


{marker examples}{...}
{title:Examples}

{pstd}Run the complete example block in the Do-file Editor, omitting leading dots.
The examples create their own data; {cmd:preserve} and {cmd:restore} retain your
original dataset on successful completion. Settings and estimation results may change.
Small samples and tuning grids demonstrate syntax rather than a research design.{p_end}
{* DEA_ML_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}

{pstd}Create reproducible by-production data:{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 100}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. generate double x1 = 1 + 4*runiform()}{p_end}
{phang2}{cmd:. generate double x2 = 1 + 3*runiform()}{p_end}
{phang2}{cmd:. generate double y = x1^0.3*x2^0.4*exp(-0.2*runiform())}{p_end}
{phang2}{cmd:. generate double z = x2^0.6*exp(-0.2*runiform())}{p_end}

{pstd}Basic BP-CEAT with equal weights:{p_end}
{phang2}{cmd:. bpceat, x1(x1) x2(x2) y(y) z(z) gen(b_equal)}{p_end}

{pstd}BP-CEAT saving all component scores:{p_end}
{phang2}{cmd:. bpceat, x1(x1) x2(x2) y(y) z(z) gen(b_parts) genpsi1(b1) genpsi2(b2)}{p_end}

{pstd}Comparing BP-EAT and BP-CEAT scores (bpceat >= bpeat):{p_end}
{phang2}{cmd:. bpeat,  x1(x1) x2(x2) y(y) z(z) gen(b_eat)}{p_end}
{phang2}{cmd:. bpceat, x1(x1) x2(x2) y(y) z(z) gen(b_ceat)}{p_end}
{phang2}{cmd:. assert b_ceat >= b_eat - 1e-8}{p_end}

{pstd}BP-CEAT emphasising the good-output component with custom weights:{p_end}
{phang2}{cmd:. bpceat, x1(x1) x2(x2) y(y) z(z) gen(b_weighted) delta1(0.7) delta2(0.3)}{p_end}

{pstd}Via unified dispatcher:{p_end}
{phang2}{cmd:. dea_ml, method(bpceat) x1(x1) x2(x2) y(y) z(z) gen(b_dispatch) seed(99)}{p_end}


{phang2}{cmd:. restore}{p_end}
{* DEA_ML_HELP_EXAMPLES_END}{...}


{marker references}{...}
{title:References}

{phang}
Guillen, M., J. Aparicio, M. Kapelko, and M. Esteve. 2025.
Measuring environmental inefficiency through machine learning: An approach
based on efficiency analysis trees and by-production technology.
{it:European Journal of Operational Research} 321(2): 529{c -}542.
{browse "https://doi.org/10.1016/j.ejor.2024.10.003"}

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
