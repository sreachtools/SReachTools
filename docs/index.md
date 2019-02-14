---
layout: home
title: "Stochastic Reachability Toolbox for MATLAB"
---

{% include news.html %}


## What does SReachTools do? 

Given a discrete-time linear system with an additive stochastic disturbance,
SReachTools constructs a set of safe initial states that satisfy some
reachability/safety specification with at least a desired likelihood
(*verification*). It can also *synthesize* controllers for these specifications
under soft/hard control bounds. 
<div class="desc-figure">
    <img src="{{ "/assets/stochTubeCartoon.png" | absolute_url }}"
    alt="Illustration of stochastic reachability of a target tube problem"/>
    <b>SReachTools can design controllers that maximize the probability of
    staying within the target tube, and characterize the set of initial states
    \(\mathcal{L}_\mathrm{SR}(\cdot)\) that satisfy a minimum reach
    probability. Image licensed under a <a rel="license"
    href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons
    Attribution-ShareAlike 4.0 International License</a>.</b>
</div>
<br>
SReachTools exploits convex optimization, Fourier transforms, and computational
geometry to obtain *scalable* results. 
<div class="desc-figure">
    <img src="{{ "/assets/scalability.png" | absolute_url }}" alt="Illustration
    of the scalability of SReachTools"/>
    <b>Scalability of SReachTools for verification of a chain of
    integrators.</b>
</div>
<br>
See [this page](long_description) for a detailed description of SReachTools.

## Can you show me some examples of SReachTools working? 

See the [examples page](https://unm-hscl.github.io/SReachTools/examples/). They
are also available in SReachTools code base, see the `examples/` folder. 

## How do I install this toolbox? What are the dependencies?

See the [installation page](installation). 
<!--Where do I get the source code from?-->
<!--See our [Github repository](https://github.com/unm-hscl/SReachTools), or our-->
<!--[release page](https://github.com/unm-hscl/SReachTools/releases) for zip files. -->

## How can I use this toolbox?

SReachTools is licensed under [GNU General Public License
v3](https://www.gnu.org/licenses/), or (at your option) any later version.  See
our [License](license/).  Please cite our
[toolpaper](https://github.com/unm-hscl/SReachTools/raw/master/SReachTools.pdf),
if it helps you in your research. 
- IEEE citation style

  A. P. Vinod, J. D. Gleason, and M. M. K. Oishi. "SReachTools: A MATLAB
  Stochastic Reachability Toolbox", In _Proceedings of the  International
  Conference on Hybrid Systems: Computation and Control_, Montreal, Canada,
  April 16--18, 2019. https://unm-hscl.github.io/SReachTools/ (accepted).

- BibTeX entry for use in LaTeX with `\usepackage{url}`: 
```
@misc{SReachTools,
    author    = {Vinod, Abraham P. and Gleason, Joseph D. and Oishi, Meeko M. K.},
    title     = {{ '{' }}{S}{R}each{T}ools: A {MATLAB} {S}tochastic {R}eachability {T}oolbox},
    booktitle = {Proceedings of the International Conference on Hybrid Systems: Computation and Control},
    year      = {2019},
    address   = {Montreal, Canada},
    month     = {April 16--18},
    note      = {\url{https://unm-hscl.github.io/SReachTools/} (accepted)}
}
```

## Where do I ask questions or give feedback? 

For better documentation, use our [Github issues
page](https://github.com/unm-hscl/SReachTools/issues).  Alternatively, see our
[Google groups page](https://groups.google.com/d/forum/sreachtools).


## Can I contribute to this toolbox?

Of course, we welcome contributions. See [Contributing guidelines](contributing/). 

## Credits

The authors of this toolbox are [Abraham P. Vinod](https://abyvinod.github.io/)
and [Joseph D.  Gleason](http://www.unm.edu/~gleasonj/). The authors are PhD
advisees of [Prof. Meeko Oishi](http://www.unm.edu/~oishi/).  SReachTools
leverages several existing toolboxes and third-party codes:
1. [MPT3](https://www.mpt3.org/) developed by M. Herceg, M. Kvasnica, C.N.
   Jones, and M. Morari, along with their dependencies.
2. [CVX](http://cvxr.com/cvx/) developed by Michael Grant and Stephen Boyd.
3. [GeoCalcLib](http://worc4021.github.io/GeoCalcLib/) developed by Rainer
   Schaich.
4. [An algorithm for numerical computation of multivariate normal distribution values](http://www.math.wsu.edu/faculty/genz/software/matlab/qscmvnv.m) developed by Alan Genz (Distributed with SReachTools under the [license](https://unm-hscl.github.io/SReachTools/docs/src/helperFunctions/qscmvnv/)).
5. [allcomb.m](https://www.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin)
developed by Jos van der Geest (Distributed with SReachTools under the [license](https://unm-hscl.github.io/SReachTools/docs/src/helperFunctions/allcomb/#license)).

When available, we also use commercial toolboxes from MATLAB and Gurobi.
