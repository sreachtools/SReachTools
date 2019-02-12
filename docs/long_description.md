---
layout: home
title: "Stochastic Reachability Toolbox: An overview"
---

SReachTools is an open-source MATLAB Toolbox for performing stochastic
verification and controller synthesis.  The authors of this toolbox are [Abraham
P.  Vinod](http://www.unm.edu/~abyvinod/) and [Joseph D.
Gleason](http://www.unm.edu/~gleasonj/). The authors are PhD advisees of [Prof.
Meeko Oishi](http://www.unm.edu/~oishi/). We had our first stable release of
this toolbox on [October, 2018](./jekyll/update/2018/10/22/release-of-v1.html).

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


