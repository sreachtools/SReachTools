---
layout: post
title:  "Release of v1 (Stable)"
date:   2018-10-22 19:27:00 -0600
categories: jekyll update
author: "Abraham P. Vinod  & Joseph D. Gleason"
---

SReachTools Toolbox [v1](https://github.com/unm-hscl/SReachTools/tree/v1.0.0) is
now out! 

{% include important-note.html content="This is the first stable release of the
toolbox." %}

The SReachTools toolbox is a set of MATLAB codes to facilitate stochastic
reachability computations. The toolbox currently supports reachability
computations with several methods:
* Dynamic programming computation (implemented in `SReachDynProg`)
* Stochastic reach set computation (implemented in `SReachSet`) and controller
  synthesis (implemented in `SReachPoint`) using:
    * Lagrangian (set-based) over- and under-approximations
    * Fourier transformation approximations, with support for open-loop
      controller synthesis
    * Chance-constrained under-approximations, with support for open-loop and
      affine controller synthesis
    * Particle filter-based approximation for open-loop controller synthesis
* Forward stochastic reachability analysis (implemented in `SReachFwd`)

We can analyze *stochastic linear time-invariant and time-varying systems*.

We have submitted a tool paper describing the features of SReachTools to the
*22nd ACM International Conference on Hybrid Systems: Computation and Control
summarizing the features of SReachTools*. A copy of this submission is
[available in the
repository](https://github.com/unm-hscl/SReachTools/raw/master/SReachTools.pdf).

See the quick-start guide below for installation instructions.

## Important notes on this release

- Updated License to GNU GPLv3, or (at your option) a later version.
  See [LICENSE](https://unm-hscl.github.io/SReachTools/license/).
- A complete overhaul of the APIs. For more details, see our [submitted tool
  paper](https://github.com/unm-hscl/SReachTools/raw/master/SReachTools.pdf).
- All promises made in our previous [blog post](./2018-05-31-release-of-v0x2.md)
  have been met:
    - Extension to LTV systems
    - Support for chance-constrained and particle filter-based verification
    - Added online API documentation. See [online documentation](https://unm-hscl.github.io/SReachTools/docs/index.html).
- Affine control synthesis technique based on chance constraints and difference
  of convex programs have been added. See our [submitted
  paper (arXiv)](https://hscl.unm.edu/affinecontrollersynthesis/) for more
  details.
- Lagrangian over-approximation has been added into the toolbox. See our
  [submitted paper (arXiv)](https://arxiv.org/pdf/1810.07118) for more details.

## Quick start guide: installation and examples

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
