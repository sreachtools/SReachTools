---
layout: home
title: "Stochastic Reachability Toolbox"
---

SReachTools focuses on the problem of stochastic reachability of a target
tube[^TAC2018_verification] --- Construct **controllers** and characterize the
**set of initial states** such that 
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
of a target tube problem**
[^TAC2018_verification]<sup>,</sup>[^HSCC2018_cvxcmpt] for stochastic linear
(time-varying or time-invariant) systems.  SReachTools can construct polytopic
(over- and under-) approximations and (open-loop and affine) controller
synthesis for this problem.  Our solution techniques include:
1. chance-constrained approaches[^CDC2013_Lesser]<sup>,</sup>[^CDC2019_chance],
1. Fourier transforms[^CSSL2017_genzps], 
1. particle control (and Voronoi partition-based undersampling)
   [^CDC2013_Lesser]<sup>,</sup>[^ACC2019_Voronoi], 
1. Lagrangian (set-operations)[^CDC2017_Lagrangian], and
1. dynamic programming[^Automatica_Abate]<sup>,</sup>[^Automatica_Summers]. 
<div class="desc-figure">
    <img src="{{ "/assets/scalability.png" | absolute_url }}" alt="Scalability
    of the various set computation techniques"/>
</div>
The above figure shows how `SReachSet` scales in the stochastic reach set
computation for a chain of integrator dynamics, in comparison with the dynamic
programming approach. Among these techniques, Lagrangian and particle control
(along with Voronoi-based extensions) can handle arbitrary disturbances.

SReachTools also provides APIs to analyze the forward stochastic reachability
problem[^HSCC2017_Fwd] using Genz's algorithm [^GenzAlgorithm].

------
[^TAC2018_verification]: A. Vinod and M. Oishi, "[Stochastic reachability of a target tube:  Theory and computation](https://arxiv.org/pdf/1810.05217.pdf)", submitted to IEEE Transactions of Automatic Control, 2018 (submitted).
[^HSCC2018_cvxcmpt]: A. Vinod and M. Oishi, "[Scalable Underapproximative Verification of Stochastic LTI Systems using Convexity and Compactness](https://doi.org/10.1145/3178126.3178148)", in Proceedings of Hybrid Systems: Computation and Control, pp. 1--10, 2018.
[^CDC2019_chance]: A. Vinod and M. Oishi, "[Affine controller synthesis for stochastic reachability via difference of convex programming](https://hscl.unm.edu/affinecontrollersynthesis/)", in Proceedings of Conference on Decision and Control, 2019 (submitted).
[^CSSL2017_genzps]: A. Vinod and M. Oishi, "[Scalable Underapproximation for Stochastic Reach-Avoid Problem for High-Dimensional LTI Systems using Fourier Transforms](https://ieeexplore.ieee.org/document/7950904/)", in IEEE Control Systems Letters (CSS-L), pp. 316--321, 2017. 
[^CDC2013_Lesser]: K. Lesser, M. Oishi, and R. S. Erwin, "[Stochastic reachability for control of spacecraft relative motion](https://doi.org/10.1109/CDC.2013.6760626)," in Proceedings of the IEEE Conference on Decision and Control, pp. 4705-4712, 2013.
[^ACC2019_Voronoi]: H. Sartipizadeh, A. Vinod,  B. Acikmese, and M. Oishi, "[Voronoi Partition-based Scenario Reduction for Fast Sampling-based Stochastic Reachability Computation of LTI Systems](https://arxiv.org/abs/1811.03643)", In Proceedings of American Control Conference, 2019 (accepted).
[^CDC2017_Lagrangian]: J. Gleason, A. Vinod, and M. Oishi, "[Underapproximation of Reach-Avoid Sets for Discrete-Time Stochastic Systems via Lagrangian Methods](https://doi.org/10.1109/CDC.2017.8264291)," in Proceedings of the IEEE Conference on Decision and Control, pp. 4283-4290, 2017.
[^Automatica_Summers]: S. Summers and J. Lygeros, "[Verification of discrete time stochastic hybrid systems: A stochastic reach-avoid decision problem](https://doi.org/10.1016/j.automatica.2010.08.006)," Automatica, 2010.  
[^Automatica_Abate]: A. Abate, M. Prandini, J. Lygeros, and S. Sastry, "[Probabilistic reachability and safety for controlled discrete time stochastic hybrid systems](https://doi.org/10.1016/j.automatica.2008.03.027)," Automatica, 2008.
[^HSCC2017_Fwd]:  A. Vinod, B. HomChaudhuri, and M. Oishi, "[Forward Stochastic Reachability Analysis for Uncontrolled Linear Systems using Fourier Transforms](https://dl.acm.org/citation.cfm?doid=3049797.3049818)", in Proceedings of the 20th International Conference on Hybrid Systems: Computation and Control (HSCC), pp. 35-44, 2017. 
[^GenzAlgorithm]: A. Genz, "[Quadrature of a multivariate normal distribution over a region specified by linear inequalities: QSCMVNV](http://www.math.wsu.edu/faculty/genz/software/matlab/qscmvnv.m)", 2014. 

