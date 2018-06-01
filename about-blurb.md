---
layout: default
title: "About Blurb"
math: true
---

SReachTools is an open-source MATLAB Toolbox for performing stochastic verification and reachability analysis.  

- [Github repository](https://github.com/abyvinod/SReachTools)
- [Quick start guide](#quick-start)
- [Google groups forum](https://groups.google.com/d/forum/sreachtools)
- [License](license/)
- [Contributing guidelines](contributing/)

This toolbox currently focuses on the following problem --- Construct **controllers** and characterize the **set of initial states** such that 
1. the controller satisfies the specified control bounds,
1. the stochastic system stays within a **safe set** to stay within the time horizon and/or reaches the **target set** at the time horizon, or stays within a time-varying **target tube** with a probability above a given threshold?\\
![A cartoon depicting the stochastic reach-avoid problem]({{ "/assets/StochReachAvoidCartoon.jpeg" | absolute_url }})\\
Here, we would like to pick the *green* controller over the *red* controller and compute the collection, the *orange set*, of all initial states \\(\overline{x}\_0\\) such that the probability of success (reachability) \\(\mathbb{P}\\) is above a given threshold \\(\theta\\).
<!--1. **Prediction problem**: -->


The authors of this toolbox are [Abraham P.  Vinod](http://www.unm.edu/~abyvinod/) and [Joseph D.  Gleason](http://www.unm.edu/~gleasonj/). Please cite their [relevant papers](https://scholar.google.com/citations?user=yb5Z7AwAAAAJ&hl=en) when using the toolbox. The authors are PhD advisees of [Prof. Meeko Oishi](http://www.unm.edu/~oishi/).

