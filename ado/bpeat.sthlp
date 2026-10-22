{smcl}
{* Help revision: 2026-10-22.1; command version unchanged.}{...}
{* *! bpeat  v0.2.1  08sep2026}{...}
{viewerjumpto "Syntax" "bpeat##syntax"}{...}
{viewerjumpto "Description" "bpeat##description"}{...}
{viewerjumpto "Options" "bpeat##options"}{...}
{viewerjumpto "Stored results" "bpeat##results"}{...}
{viewerjumpto "Examples" "bpeat##examples"}{...}
{viewerjumpto "References" "bpeat##references"}{...}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col:{cmd:bpeat} {hline 2}}By-production Efficiency Analysis Trees (BP-EAT){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:bpeat} {ifin}{cmd:,}
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
{synopt:{opt fold(#)}}cross-validation folds for EAT pruning (default 5){p_end}
{synopt:{opt seed(#)}}random seed for CV fold assignment; unused with {cmd:fold(0)}{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:bpeat} estimates a by-production model for good and bad outputs. The
model separates a good-output technology from a bad-output technology and
represents each component with an Efficiency Analysis Tree (EAT):

{pstd}
{bf:Psi1 (good-output technology):} An EAT is fitted with inputs (X1, X2) and
good outputs Y. This tree characterises the production frontier for desirable
outputs. The component score beta_Psi1 measures good-output inefficiency
relative to this frontier.

{pstd}
{bf:Psi2 (bad-output technology):} An EAT is fitted with bad outputs Z and
pollution-generating inputs X2. The component score beta_Psi2 measures
bad-output inefficiency relative to this technology.

{pstd}
A zero component of Z has a zero contraction direction, but its technology
constraint remains active: the selected leaf must also have a zero lower
bound in that component. Such a component is not simply discarded.

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
where delta1 and delta2 are user-specified weights (default 0.5 each, summing
to 1). All scores are bounded below by 0. A higher B indicates greater
environmental inefficiency in the by-production system.

{pstd}
EAT trees are pruned by cost-complexity cross-validation, yielding monotone
non-decreasing step-function frontiers.


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
before a split is attempted in either EAT. Default is 5.

{phang}
{opt fold(#)} sets the number of cross-validation folds used for EAT
cost-complexity pruning in both phases. Default is 5.

{phang}
{opt seed(#)} sets the random seed for cross-validation. It is unused with
{cmd:fold(0)}. Default is 12345.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:bpeat} stores the following in {cmd:e()}:

{synoptset 18 tabbed}{...}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sample)}}estimation-sample indicator{p_end}
{synopt:{cmd:e(mean_b)}}mean combined efficiency score B{p_end}
{synopt:{cmd:e(mean_beta1)}}mean first-component directional score{p_end}
{synopt:{cmd:e(mean_beta2)}}mean second-component directional score{p_end}
{synopt:{cmd:e(n_leaves1)}}number of leaves in the first component tree{p_end}
{synopt:{cmd:e(n_leaves2)}}number of leaves in the second component tree{p_end}
{synopt:{cmd:e(delta1)}}weight on beta_Psi1{p_end}
{synopt:{cmd:e(delta2)}}weight on beta_Psi2{p_end}
{synopt:{cmd:e(x1)}}fixed input variable names{p_end}
{synopt:{cmd:e(x2)}}discretionary input variable names{p_end}
{synopt:{cmd:e(y)}}output variable names{p_end}
{synopt:{cmd:e(z)}}intermediate output variable names{p_end}
{synopt:{cmd:e(cmd)}}{cmd:bpeat}{p_end}
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

{pstd}Basic BP-EAT with equal weights:{p_end}
{phang2}{cmd:. bpeat, x1(x1) x2(x2) y(y) z(z) gen(b_equal) delta1(0.5) delta2(0.5)}{p_end}

{pstd}BP-EAT saving all component scores:{p_end}
{phang2}{cmd:. bpeat, x1(x1) x2(x2) y(y) z(z) gen(b_parts) genpsi1(b1) genpsi2(b2)}{p_end}

{pstd}BP-EAT emphasising the good-output component with custom weights:{p_end}
{phang2}{cmd:. bpeat, x1(x1) x2(x2) y(y) z(z) gen(b_weighted) delta1(0.7) delta2(0.3)}{p_end}

{pstd}Via unified dispatcher:{p_end}
{phang2}{cmd:. dea_ml, method(bpeat) x1(x1) x2(x2) y(y) z(z) gen(b_dispatch) seed(99)}{p_end}


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
