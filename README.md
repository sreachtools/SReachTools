# Stochastic reachability toolbox (SReachTools)

SReachTools is a MATLAB toolbox to tackle various problems in stochastic
reachability.

This is an area of **active research**, and this toolbox will attempt to cater
certain classes of problems.  

We aim to support the following problems:
 - **Stochastic reachability of a target tube** (guaranteeing safety for stochastic
   systems to lying in a collection of time-varying safe sets while satisfying
   input bounds)
    - This problem subsumes existing work on terminal hitting stochastic reach-avoid
      problems as well as stochastic viability problems.
    - **Open-loop controller synthesis** (admissible controller satisfying hard
      control bounds with maximum safety probability):
        - `chance-open`: Chance constraint formulation solved via linear
          programming
        - `genzps-open`: Fourier transforms-based compuation (Genz's algorithm +
          patternsearch)
        - `particle-open`: Particle filter approach (mixed-integer linear
          program approach)
    - **Affine controller synthesis** (admissible controller with chance constrained
      input bounds with maximum safety probability)
        - `chance-affine`: Chance constraint formulation solved via
          difference-of-convex programming
    - **Stochastic reach set computation** (set of initial states from which an 
      admissible controller exists such that the probability of safety is above a 
      given threshold)
        - `chance-open`: Chance constraint-based under-approximation
        - `genzps-open`: Fourier transforms-based under-approximation
        - `lag-over/lag-under`: Lagrangian methods-based over- and
          under-approximation
 - **Forward stochastic reachability** (characterizing the stochasticity of the
      state at a future time of interest)
      - `state-stoch/concat-stoch`: Stochasticity of the state or the
          concatenated state vector
      - `state-prob/concat-prob`: Probability of the state or the concatenated
          state vector lying in a target set or a tube respectively

Do check our [project blog](https://unm-hscl.github.io/SReachTools/blog/) for
updates!

## Installation, documentation, and examples

### Dependencies

You can skip installing the dependencies marked **optional**.
This will disable some of the features of SReachTools.

1. MATLAB (>2017a)
    1. Toolboxes
        1. MATLAB's Statistics and Machine Learning Toolbox
        1. MATLAB's Global Optimization Toolbox (**Optional**)
1. MPT3 ([https://www.mpt3.org/](https://www.mpt3.org/))
    1. Do an automatic install using a MATLAB script
       [install_mpt3.m](https://www.mpt3.org/Main/Installation?action=download&upname=install_mpt3.m)
       provided by MPT3.
1. CVX ([http://cvxr.com/cvx/](http://cvxr.com/cvx/))
    1. Install the CVX (Standard bundle, including Gurobi and/or MOSEK)
    1. Installation instructions are given in [http://cvxr.com/cvx/download/](http://cvxr.com/cvx/download/).
1. We recommend using Gurobi as the backend solver for the convex programs
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

### Examples

For easy start, we have cataloged in our [project
webpage](https://unm-hscl.github.io/SReachTools/examples/) a number of relevant,
easy-to-follow examples. These are also part of the repository (see
`examples/*.m`). 

## Contributions

SReachTools is an open-source MATLAB toolboxes. We welcome feedback, issues,
bug-fixes, and other enhancements.  Please see [our contribution
guidelines](./CONTRIBUTING.md) of this project to know more on how you can help
us.

## License

See [LICENSE](./LICENSE).

## Credits

The authors of this toolbox are [Abraham P.
Vinod](https://abyvinod.github.io/) (primary maintainer) and [Joseph D.
Gleason](http://www.unm.edu/~gleasonj/).  Please cite their [relevant
papers](https://scholar.google.com/citations?user=yb5Z7AwAAAAJ&hl=en) when using
the toolbox.  The authors are PhD advisees of [Prof. Meeko
Oishi](http://www.unm.edu/~oishi/).
