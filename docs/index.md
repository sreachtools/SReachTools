---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
title: "Stochastic Reachability Toolbox"
---

SReachTools is an open-source MATLAB Toolbox for performing stochastic
verification and reachability analysis. We had our first stable release of
this toolbox on [October, 2018](./jekyll/update/2018/10/22/release-of-v1.html).

- Can you show me some examples of SReachTools working? 
    - We have cataloged a number of
      [examples](https://unm-hscl.github.io/SReachTools/examples/) implemented
      using SReachTools. These examples are also available as part of the source
      code of SReachTools, see `examples/*.m`. 
- How do I install this? What are the dependencies?
    - Our [quick start guide](#quick-start-guide), described further
      down this page, walks through the installation process.
- Where do I get the source code from?
    - See our [Github repository](https://github.com/unm-hscl/SReachTools), or
      our [release page](https://github.com/unm-hscl/SReachTools/releases) for
      zip files. 
- How can I use this toolbox? What are the terms and conditions to follow to use
  SReachTools?
    - SReachTools is licensed under [GNU General Public License v3]
      (https://www.gnu.org/licenses/), or (at your option) any later version.
      See our [License](license/).
- Can I contribute to this toolbox?
    - Of course, we welcome pull requests. See [Contributing guidelines](contributing/). 
- Where do I ask questions or give feedback? 
    - Use our [Google groups](https://groups.google.com/d/forum/sreachtools) or
      the [Github issues](https://github.com/unm-hscl/SReachTools/issues) page.

The authors of this toolbox are [Abraham P.
Vinod](http://www.unm.edu/~abyvinod/) and [Joseph D.
Gleason](http://www.unm.edu/~gleasonj/). Please cite their [relevant
papers](https://scholar.google.com/citations?user=yb5Z7AwAAAAJ&hl=en) when using
the toolbox. The authors are PhD advisees of [Prof. Meeko
Oishi](http://www.unm.edu/~oishi/).

We have submitted a tool paper describing the features of SReachTools to the
*22nd ACM International Conference on Hybrid Systems: Computation and Control
summarizing the features of SReachTools*. A copy of this submission is
[available in the
repository](https://github.com/unm-hscl/SReachTools/raw/master/SReachTools.pdf).

{% include news.html %}

## What does SReachTools do?

SReachTools focuses on the problem of stochastic reachability of a target
tube[^1] --- Construct **controllers** and characterize the **set of initial
states** such that 
1. the controller satisfies the specified control bounds,
1. the stochastic system stays within a time-varying target tube with a
   probability above a given threshold. 

For example, a typical **reach-avoid** constraint is to stay within a *safe set*
to stay within the time horizon and reach a *target set* at the time horizon
when starting from an initial state \\(\overline{x}\_0\\), as shown in the
figure below.
<div class="desc-figure">
    <img src="{{ "/assets/StochReachAvoidCartoon.jpeg" | absolute_url }}" alt="A
    cartoon depicting the stochastic reach-avoid problem"/>
</div>
Here, we would like to pick the *green* controller over the *red* controller and
compute the collection, the *orange set*, of all initial states such that the
probability of success (reach-avoid) \\(\mathbb{P}\\) is above a given threshold
\\(\theta\\).

This problem appears in a wide range of applications --- space applications
(spacecraft rendezvous and docking problem), transport
(semi-autonomous/fully-autonomous cars and airplanes), biomedical applications
(automated anesthesia delivery system), to name a few.  Some of these examples
have been analyzed using SReachTools. See the `examples` folder.

Our approaches rely on **convex optimization**, **stochastic programming**,
**Fourier transforms**, and **computational geometry** to provide *scalable*,
*grid-free*, and anytime algorithms for stochastic reachability analysis of
linear systems.  Specifically, SReachTools tackles the **stochastic reachability
of a target tube problem** [^TAC2018_verification] [^HSCC2018_cvxcmpt] for
Gaussian-perturbed linear (time-varying or time-invariant) systems.  SReachTools
can construct polytopic (over- and under-) approximations and (open-loop and
affine) controller synthesis for this problem.  Our solution techniques include:
1. chance-constrained approaches[^CDC2013_Lesser] [^HSCC2019_chance],
1. Fourier transforms[^CSSL2017_genzps], 
1. particle filter control [^CDC2013_Lesser] [^ACC2019_Voronoi], 
1. Lagrangian (set-operations)[^CDC2017_Lagrangian], and
1. dynamic programming[^Automatica_Abate] [^Automatica_Summers]. 
<div class="desc-figure">
    <img src="{{ "/assets/scalability.png" | absolute_url }}" alt="Scalability
    of the various set computation techniques"/>
</div>
The above figure shows how `SReachSet` scales in the stochastic reach set
computation for a chain of integrator dynamics, in comparison with the dynamic
programming approach.

SReachTools also provides APIs to analyze the forward stochastic reachability
problem[^HSCC2017_Fwd] using Genz's algorithm [^GenzAlgorithm]. 

The following
table[^table_ack] summarizes the features of SReachTools.

|    Function   |   method-str  |                                                       Utility                                                       | Notes                                      |
|:-------------:|:---------------:|:-------------------------------------------------------------------------------------------------------------------:|--------------------------------------------|
| `SReachPoint` |                 |          **Approximation of the maximal reach  probability for a target tube from  a given initial state** [^TAC2018_verification]         | **Synthesize open-loop or affine disturbance feedback controllers** |
|               |  `chance-open`  |                                            Guaranteed underapproximation [^CDC2013_Lesser] [^HSCC2019_chance]                                             | Open-loop                                  |
|               |  `genzps-open`  |                                  Approximate up to \\( \\epsilon\_\\mathrm{genz}\\), a user-specified quadrature error tolerance [^CSSL2017_genzps]                                 | Open-loop                                  |
|               | `particle-open` |                        Approximate with quality proportional  to the number of particles used [^CDC2013_Lesser]                     | Open-loop                                  |
|               |  `voronoi-open` |                          Probabilistically enforced upper  bound on overapproximation error  [^ACC2019_Voronoi]                          | Open-loop                                  |
|               | `chance-affine` |                                            Guaranteed underapproximation  [^HSCC2019_chance]                                            | Affine   disturbance-feedback              |
|  `SReachSet`  |                 |  **Polytopic approximation of the stochastic  reach sets for the stochastic reachabilty  of a target tube problem**[^TAC2018_verification] [^HSCC2018_cvxcmpt] | **Synthesize open-loop controllers in some cases** |
|               |  `chance-open`  |                                            Guaranteed underapproximation  [^TAC2018_verification]                                           | Optimal  open-loop controllers at vertices |
|               |  `genzps-open`  |                                 Approximation up to \\( \\epsilon\_\\mathrm{genz}\\), a user-specified quadrature error tolerance  [^TAC2018_verification] [^HSCC2018_cvxcmpt]                                | Optimal  open-loop controllers at vertices |
|               |   `lag-under`   |                                            Guaranteed underapproximation [^CDC2017_Lagrangian]                                            |                                            |
|               |    `lag-over`   |                                             Guaranteed overapproximation [^CDC2017_Lagrangian]                                            |                                            |
|  `SReachFwd`  |                 |      **Forward stochastic reachability analysis of an uncontrolled LTI/LTV system from a given initial state** [^HSCC2017_Fwd] [^GenzAlgorithm]    |                                            |
|               |  `state-stoch`  |                                     Stochasticity of the state at a future time                                     |                                            |
|               |  `concat-stoch` |                     Stochasticity of the concatenated state vector up to a specified future time                    |                                            |
|               | `state-prob`    | Probability that the concatenated state vector (trajectory) up to a future time will lie in a given target tube set [^GenzAlgorithm] |                                            |
|               | `concat-prob`   | Probability that the state at a future time will lie in a given target set [^GenzAlgorithm]                                          |                                            |
|  `SReachDyn`  |                 |              **Dynamic programming approximation of the maximal reach probability and the reach  set**              |  **Analyze 2D and 3D LTI/LTV systems**     |


## Quick start guide

### Dependencies

You can skip installing the dependencies marked **optional**.
This will disable some of the features of SReachTools or hamper performance.

1. MATLAB (>2017a)
    1. Toolboxes
        1. MATLAB's Statistics and Machine Learning Toolbox
        1. MATLAB's Global Optimization Toolbox (**Optional**)
1. MPT3 ([https://www.mpt3.org/](https://www.mpt3.org/))
    1. Copy the MATLAB script [install_mpt3.m](https://www.mpt3.org/Main/Installation?action=download&upname=install_mpt3.m) provided by MPT3 from the browser, and run it in MATLAB to automatically download MPT3 and its dependencies.
1. CVX ([http://cvxr.com/cvx/](http://cvxr.com/cvx/))
    1. Install the CVX (Standard bundle, including Gurobi and/or MOSEK)
    1. Installation instructions are given in [http://cvxr.com/cvx/download/](http://cvxr.com/cvx/download/).
1. (**Optional**) We recommend using Gurobi as the backend solver for the convex programs
   formulated by SReachTools. In practice, we find both CVX and MPT3 perform
   much better with Gurobi. See
   [http://www.gurobi.com/registration/download-reg](http://www.gurobi.com/registration/download-reg)
   for more details. Note that Gurobi offers free academic license.

### Installation

1. Install the necessary dependencies listed above
1. Clone the SReachTools repository (or download the latest zip file from
   [Releases](https://github.com/unm-hscl/SReachTools/releases))
1. Change the MATLAB current working directory to where SReachTools was
   downloaded
1. Run `srtinit` in MATLAB to add the toolbox to the paths and ensure all
   must-have dependencies are properly installed.
   - You can add `cd <path_to_sreachtools_repo>;srtinit` to your MATLAB's
     `startup.m` to automatically have this done in future.
   - Additional steps (optional):
       - Run `srtinit -t` to run all the unit tests.
       - Run `srtinit -v` to visualize the steps the changes to the path and
         check for recommended dependencies.  
       - Run `srtinit -x` to remove functions of SReachTools from MATLAB's path
         after use.  

------
[^TAC2018_verification]: A. Vinod and M. Oishi, "[Stochastic reachability of a target tube:  Theory and computation](https://arxiv.org/pdf/1810.05217.pdf)", submitted to IEEE Transactions of Automatic Control, 2018 (submitted).
[^HSCC2018_cvxcmpt]: A. Vinod and M. Oishi, "[Scalable Underapproximative Verification of Stochastic LTI Systems using Convexity and Compactness](https://doi.org/10.1145/3178126.3178148)", in Proceedings of Hybrid Systems: Computation and Control, pp. 1--10, 2018.
[^HSCC2019_chance]: A. Vinod and M. Oishi, "[Affine controller synthesis for stochastic reachability via difference of convex programming](https://hscl.unm.edu/affinecontrollersynthesis/)", in Proceedings of Hybrid Systems: Computation and Control, 2019 (submitted).
[^CSSL2017_genzps]: A. Vinod and M. Oishi, "[Scalable Underapproximation for Stochastic Reach-Avoid Problem for High-Dimensional LTI Systems using Fourier Transforms](https://ieeexplore.ieee.org/document/7950904/)", in IEEE Control Systems Letters (CSS-L), pp. 316--321, 2017. 
[^CDC2013_Lesser]: K. Lesser, M. Oishi, and R. S. Erwin, "[Stochastic reachability for control of spacecraft relative motion ](https://doi.org/10.1109/CDC.2013.6760626)," in Proceedings of the IEEE Conference on Decision and Control, pp. 4705-4712, 2013.
[^ACC2019_Voronoi]: H. Sartipizadeh, A. Vinod,  B. Acikmese, and M. Oishi, "Voronoi Partition-based Scenario Reduction for Fast Sampling-based Stochastic Reachability Computation of LTI Systems", In Proceedings of American Control Conference, 2019 (submitted).
[^CDC2017_Lagrangian]: J. Gleason, A. Vinod, and M. Oishi, "[Underapproximation of Reach-Avoid Sets for Discrete-Time Stochastic Systems via Lagrangian Methods](https://doi.org/10.1109/CDC.2017.8264291)," in Proceedings of the IEEE Conference on Decision and Control, pp. 4283-4290, 2017.
[^Automatica_Summers]: S. Summers and J. Lygeros, "[Verification of discrete time stochastic hybrid systems: A stochastic reach-avoid decision problem](https://doi.org/10.1016/j.automatica.2010.08.006)," Automatica, 2010.  
[^Automatica_Abate]: A. Abate, M. Prandini, J. Lygeros, and S. Sastry, "[Probabilistic reachability and safety for controlled discrete time stochastic hybrid systems](https://doi.org/10.1016/j.automatica.2008.03.027)," Automatica, 2008.
[^HSCC2017_Fwd]:  A. Vinod, B. HomChaudhuri, and M. Oishi, "[Forward Stochastic Reachability Analysis for Uncontrolled Linear Systems using Fourier Transforms](https://dl.acm.org/citation.cfm?doid=3049797.3049818)", in Proceedings of the 20th International Conference on Hybrid Systems: Computation and Control (HSCC), pp. 35-44, 2017. 
[^GenzAlgorithm]: A. Genz, "[Quadrature of a multivariate normal distribution over a region specified by linear inequalities: QSCMVNV](http://www.math.wsu.edu/faculty/genz/software/matlab/qscmvnv.m)", 2014. 
[^table_ack]: This table was generated using [https://www.tablesgenerator.com/markdown_tables#](https://www.tablesgenerator.com/markdown_tables#).
