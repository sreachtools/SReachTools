---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
title: "Stochastic Reachability Toolbox"
---

SReachTools is an open-source MATLAB Toolbox for performing stochastic verification and reachability analysis.  

- Where do I get the code from? See our [Github repository](https://github.com/unm-hscl/SReachTools).
- How do I install this? What do I need to get it working? See our [quick start guide](#quick-start-guide).
- Can you show me some examples? See our [examples page](./examples.md).
- Where do I ask/give questions/feedback/comments? Use our [Google groups](https://groups.google.com/d/forum/sreachtools) or the [Github issues](https://github.com/unm-hscl/SReachTools/issues) page.
- I would like to contribute to this toolbox. See [Contributing guidelines](contributing/).
- How can I use this toolbox? See our [License](license/).

The authors of this toolbox are [Abraham P.  Vinod](http://www.unm.edu/~abyvinod/) and [Joseph D.  Gleason](http://www.unm.edu/~gleasonj/). Please cite their [relevant papers](https://scholar.google.com/citations?user=yb5Z7AwAAAAJ&hl=en) when using the toolbox. The authors are PhD advisees of [Prof. Meeko Oishi](http://www.unm.edu/~oishi/).

{% include news.html %}

## What does SReachTools do?

SReachTools focuses on the following problem --- Construct **controllers** and characterize the **set of initial states** such that 
1. the controller satisfies the specified control bounds,
1. the stochastic system stays within a time-varying **target tube** with a probability above a given threshold? \\
For example, a typical **reach-avoid** constraint is to stay within a *safe set* to stay within the time horizon and reach a *target set* at the time horizon when starting from an initial state \\(\overline{x}\_0\\), as shown in the figure below.
<div class="desc-figure">
    <img src="{{ "/assets/StochReachAvoidCartoon.jpeg" | absolute_url }}" alt="A cartoon depicting the stochastic reach-avoid problem"/>
    <!-- [A cartoon depicting the stochastic reach-avoid problem]({{ "/assets/StochReachAvoidCartoon.jpeg" | absolute_url }})\\ -->
</div>
Here, we would like to pick the *green* controller over the *red* controller and compute the collection, the *orange set*, of all initial states such that the probability of success (reach-avoid) \\(\mathbb{P}\\) is above a given threshold \\(\theta\\).

This problem appears in a wide range of applications --- space applications ([spacecraft rendezvous and docking problem](./examples/FtCVXUnderapproxVerifyCWH.html)), transport (semi-autonomous/fully-autonomous cars and airplanes), biomedical applications (automated anesthesia delivery system), to name a few.

This toolbox provides MATLAB APIs to tackle this problem for Gaussian-perturbed linear time-invariant systems using [Fourier transforms](./examples/FtCVXUnderapproxVerifyCWH.html) [^1], [Lagrangian (set-operations)](./examples/lagrangianApproximations.html) [^2], and [dynamic programming](./examples/doubleIntegratorDynamicProgramming.html) [^3] [^4] methods.
We currently provide polytopic underapproximation and open-loop controller synthesis for this problem.
In future, we will provide extensions to linear time-varying systems, closed-loop controller synthesis, and non-Gaussian disturbances.

[^1]: A. P. Vinod and M. M. K. Oishi, "[Scalable Underapproximative Verification of Stochastic LTI Systems using Convexity and Compactness](https://doi.org/10.1145/3178126.3178148)", in Proceedings of Hybrid Systems: Computation and Control, 2018
[^2]: J. D. Gleason, A. P. Vinod, M. M. K. Oishi, "[Underapproximation of Reach-Avoid Sets for Discrete-Time Stochastic Systems via Lagrangian Methods](https://doi-org/10.1109/CDC.2017.8264291)," in Proceedings of the IEEE Conference on Decision and Control, 2017
[^3]: S. Summers and J. Lygeros, "[Verification of discrete time stochastic hybrid systems: A stochastic reach-avoid decision problem](https://doi.org/10.1016/j.automatica.2010.08.006)," Automatica, 2010.
[^4]: A. Abate, M. Prandini, J. Lygeros, S. Sastry, "[Probabilistic reachability and safety for controlled discrete time stochastic hybrid systems](https://doi.org/10.1016/j.automatica.2008.03.027)," Automatica, 2008.

## Quick start guide

### Dependencies

You can skip installing the dependencies marked **optional**.
This will disable some of the features of SReachTools.

1. MATLAB (2017a or newer)
    * Toolboxes
        * MATLAB's Global Optimization Toolbox (**Optional**)
        * MATLAB's Statistics and Machine Learning Toolbox (**Optional**)
        * MATLAB's Control System Toolbox (**Optional**)
1. MPT3 ([http://people.ee.ethz.ch/~mpt/3/](http://people.ee.ethz.ch/~mpt/3/))
    * Do an automatic install using a MATLAB script [install_mpt3.m](http://control.ee.ethz.ch/~mpt/3/Main/Installation?action=download&upname=install_mpt3.m) provided by MPT3.
1. CVX ([http://cvxr.com/cvx/](http://cvxr.com/cvx/)) (**Optional**)

### Installation

1. Install the necessary dependencies (MATLAB and MPT3 are a must)
1. Clone the *SReachTools* repository (or download the zip file)
1. Run `srtinit -v -t` in MATLAB to add the toolbox to the paths, visualize the steps, and test the installation.  
   - You can add `p=pwd();cd /path/to/SReachTools/folder;srtinit;cd(p);` to your MATLAB's `startup.m` to automatically have this done in future.

------
