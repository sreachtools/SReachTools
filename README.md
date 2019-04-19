# Stochastic reachability toolbox (SReachTools)

SReachTools is a MATLAB toolbox to tackle various problems in stochastic
reachability. Currently, the toolbox can perform verification and (open-loop or
affine disturbance-feedback) controller synthesis for linear
(time-varying/time-invariant) systems with additive (Gaussian/non-Gaussian)
disturbance. By verification, we are referring to the problem of [*stochastic
reachability of a target tube*](https://arxiv.org/abs/1810.05217). Our project
website is at [https://sreachtools.github.io](https://sreachtools.github.io).

This is an area of **active research**, and this toolbox will attempt to cater
certain classes of problems.  

We aim to support the following problems:
 - **Stochastic reachability of a target tube** (guaranteeing safety for
   stochastic systems to lying in a collection of time-varying safe sets while
   satisfying input bounds):
    - This problem subsumes existing work on terminal hitting stochastic
      reach-avoid problems as well as stochastic viability problems. We
      implement a dynamic programming solution, limited to 3-dimensional LTI
      systems using `SReachDynProg`.
    - **Open-loop controller synthesis** using `SReachPoint` (admissible
      controller satisfying hard control bounds with maximum safety
      probability):
        - `chance-open`: Chance constraint formulation solved via linear
          programming
        - `genzps-open`: Fourier transforms-based compuation (Genz's algorithm +
          patternsearch)
        - `particle-open`: Particle control approach (mixed-integer linear
          program approach)
        - `voronoi-open`: Voronoi partition-based undersampled particle control
          approach (mixed-integer linear program approach)
    - **Affine controller synthesis** using `SReachPoint` (admissible controller
      with chance constrained input bounds with maximum safety probability):
        - `chance-affine`: Chance constraint formulation solved via
          difference-of-convex programming (risk allocation and controller
          synthesis performed simultaneously)
        - `chance-affine-uni`: Chance constraint formulation solved via
          bisection for uniform risk allocation and second order cone programs
          for controller synthesis (risk allocation and controller synthesis
          performed separately)
    - **Stochastic reach set computation** using `SReachSet` (set of initial
      states from which an admissible controller exists such that the
      probability of safety is above a given threshold):
        - `chance-open`: Chance constraint-based under-approximation
        - `genzps-open`: Fourier transforms-based under-approximation
        - `lag-over/lag-under`: Lagrangian methods-based over- and
          under-approximation
 - **Forward stochastic reachability** using `SReachFwd` (characterizing the
   stochasticity of the state at a future time of interest):
      - `state-stoch/concat-stoch`: Stochasticity of the state or the
          concatenated state vector
      - `state-prob/concat-prob`: Probability of the state or the concatenated
          state vector lying in a target set or a tube respectively

Do check our [project blog](https://sreachtools.github.io/blog/) for
updates!

## Examples

For easy start, we have cataloged in our [project
webpage](https://sreachtools.github.io/examples/) a number of relevant,
easy-to-follow examples. These are also part of the repository (see
`examples/*.m`). 

Further, you can see SReachTools in action at Code Ocean. Check out
[https://codeocean.com/explore/capsules/?query=SReachTools](https://codeocean.com/explore/capsules/?query=SReachTools).


## Installation

### Dependencies

You can skip installing the dependencies marked **optional**.
This will disable some of the features of SReachTools or hamper performance.
You can skip installing the dependencies marked **optional**.
This will disable some of the features of SReachTools or hamper performance.
1. MATLAB (>2017a)
    1. Toolboxes
        1. MATLAB's Statistics and Machine Learning Toolbox
        1. (**Optional**) MATLAB's Global Optimization Toolbox --- required for
           `genzps-open` options in `SReachPoint` and `SReachSet`
        1. (**Optional**) MATLAB's Optimization Toolbox --- recommended
           installation for MATLAB's Global Optimization Toolbox
1. MPT3 ([https://www.mpt3.org/](https://www.mpt3.org/)) --- for polytopic
   computational geometry
    1. Copy the MATLAB script [install_mpt3.m](https://www.mpt3.org/Main/Installation?action=download&upname=install_mpt3.m)
       provided by MPT3 from the browser, and run it in MATLAB to automatically
       download MPT3 and its dependencies.
1. CVX v2.1 ([http://cvxr.com/cvx/](http://cvxr.com/cvx/)) --- for
       parsing convex and mixed-integer programs
    1. Install the CVX (Standard bundle, including Gurobi and/or MOSEK)
    1. Installation instructions are given in
       [http://cvxr.com/cvx/download/](http://cvxr.com/cvx/download/).
    1. No license is required, if you do not plan on using Gurobi (see next step). See [http://web.cvxr.com/cvx/doc/intro.html#licensing](http://web.cvxr.com/cvx/doc/intro.html#licensing) for more details.
1. (**Optional**) Gurobi --- recommended backend solver for the convex programs
   formulated by SReachTools and required for all particle-based approaches
   in `SReachPoint`. We also find both CVX and MPT3 perform much better with
   Gurobi.
    1. CVX works best with Gurobi v7.5.2. See installation details:
       http://web.cvxr.com/cvx/doc/gurobi.html#gurobi        
    1. To use Gurobi with CVX, we requires two licenses (one for CVX and one for
       Gurobi). Both of these license are free for non-commercial academic
       research.
        1. Gurobi offers free academic license, which can be requested at
           [http://www.gurobi.com/registration/download-reg](http://www.gurobi.com/registration/download-reg).
        1. CVX provides free academic license, which can be requested at
           [http://cvxr.com/cvx/academic/](http://cvxr.com/cvx/academic/).
    1. MPT3 will automatically update its backend solver to Gurobi, when Gurobi
       is installed as a standalone and the license is found.
1. (**Optional**) [GeoCalcLib](https://github.com/worc4021/GeoCalcLib) --- a
   MATLAB interface to Avis's [LRS vertex-facet enumeration
   library](http://cgm.cs.mcgill.ca/~avis/C/lrs.html), an alternative to MPT's
   preferred approach for vertex-facet enumeration,
   [CDD](https://www.inf.ethz.ch/personal/fukudak/cdd_home/index.html).

    1. **NOTE**: GeoCalcLib currently works only Unix
    and MAC OS.  SReachTools will gracefully switch back to CDD, if installation
    of GeoCalcLib is not correct.
    1. Install [GMP](https://gmplib.org/)
        1. Get the tar ball from [https://gmplib.org/#DOWNLOAD](https://gmplib.org/#DOWNLOAD)
        1. Follow the installation instructions
           [https://gmplib.org/manual/Installing-GMP.html#Installing-GMP](https://gmplib.org/manual/Installing-GMP.html#Installing-GMP)
    1. Download the zip file from
       [https://github.com/worc4021/GeoCalcLib/archive/master.zip](https://github.com/worc4021/GeoCalcLib/archive/master.zip).
    1. Extract the contents of this zip file to a desired location, whose full path is referred to here as `/path/to/GeoCalcLib`
    1. Open a terminal and change directory to GeoCalcLib by `$cd /path/to/GeoCalcLib`. We will refer to this location as the GeoCalcLib root folder.
    1. Create a folder `mexfiles` in GeoCalcLib root folder. 
    1. Create a file named `User.make` in GeoCalcLib root folder using your
       favorite editor with the following contents, and save it. See
       [https://www.mathworks.com/matlabcentral/answers/66570-what-is-the-default-installation-path-for-matlab-on-architecture-x#answer_78163](https://www.mathworks.com/matlabcentral/answers/66570-what-is-the-default-installation-path-for-matlab-on-architecture-x#answer_78163)
       for hints on how to identify your matlab root folder for your OS.
        ```
        # Specify the absolute path to the root folder of your Matlab
        # installation where <FULL-PATH-TO-YOUR-MATLAB-INSTALLATION>/bin/mex
        # exists
        MATLABROOT = <FULL-PATH-TO-YOUR-MATLAB-INSTALLATION>
        
        # Path to which everything should be installed
        INSTALLDIR = ../mexfiles/
        ```
    1. In the command prompt in GeoCalcLib root folder, execute `$ make`.
    1. Add `/path/to/GeoCalcLib/mexfiles` to MATLAB path. If you want to use
       this across sessions, we recommend adding the following command to
       your MATLAB startup.
        ```
        addpath('/path/to/GeoCalcLib/mexfiles');
        ```

### Installation

1. Install the necessary dependencies listed above
1. Clone the SReachTools repository (or download the latest zip file from
   [Releases](https://github.com/sreachtools/SReachTools/releases))
1. Change the MATLAB current working directory to where SReachTools was
   downloaded. **WARNING**: Please do not add the folder to the MATLAB path
   manually.
1. Run `srtinit` in MATLAB to add the toolbox to the paths and ensure all
   must-have dependencies are properly installed.
   - You can add `cd <path_to_sreachtools_repo>;srtinit` to your MATLAB's
     `startup.m` to automatically have this done in future.
   - (**Optional**) Additional steps:
       - Run `srtinit -t` to run all the unit tests.
       - Run `srtinit -v` to visualize the steps the changes to the path and
         check for recommended dependencies.  
       - Run `srtinit -x` to remove functions of SReachTools from MATLAB's path
         after use.  

## Contributions

SReachTools is an open-source MATLAB toolboxes. We welcome feedback, issues,
bug-fixes, and other enhancements.  Please see [our contribution
guidelines](./CONTRIBUTING.md) of this project to know more on how you can help
us.

## License

The Stochastic Reachability Toolbox (SReachTools) is free software: you can
redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.  You
should have received a copy of the GNU General Public License along with this
toolbox (see [LICENSE](./LICENSE)).  If not, see
<https://www.gnu.org/licenses/>.

It is the user's responsibility in assessing the correctness of the theory and
software implementation before putting it to use in their own research or
exploiting the results commercially. We are, however, very happy to answer any
questions and investigate any bug reports.

## Credits

This toolbox was developed by [Abraham P. Vinod](https://abyvinod.github.io/)
and [Joseph D. Gleason](http://www.unm.edu/~gleasonj/), under the supervision
of [Prof. Meeko Oishi](http://www.unm.edu/~oishi/).

If this toolbox comes handy in your research, please consider citing our
work. A copy of this paper is [available in the
repository](https://github.com/sreachtools/SReachTools/raw/master/SReachTools.pdf).

IEEE citation style:

A. P. Vinod, J. D. Gleason, and M. M. K. Oishi. "SReachTools: A MATLAB
Stochastic Reachability Toolbox," In Proceedings of the  International
Conference on Hybrid Systems: Computation and Control, Montreal, Canada,
April 16--18, 2019.  https://sreachtools.github.io.
    
BibTeX entry for use in LaTeX with `\usepackage{url}`: 
```
@misc{SReachTools,
    author    = {Vinod, Abraham P. and Gleason, Joseph D. and Oishi, Meeko M. K.},
    title     = {{ '{' }}{S}{R}each{T}ools: A {MATLAB} {S}tochastic {R}eachability {T}oolbox},
    booktitle = {Proceedings of the International Conference on Hybrid Systems: Computation and Control},
    year      = {2019},
    address   = {Montreal, Canada},
    month     = {April 16--18},
    note      = {\url{https://sreachtools.github.io}}
}
```
As seen from the installation instructions, this toolbox leverages several
existing toolboxes and third-party codes:
1. [MPT3](https://www.mpt3.org/) developed by M. Herceg, M. Kvasnica, C.N. Jones, and M. Morari, along with their dependencies.
2. [CVX](http://cvxr.com/cvx/) developed by Michael Grant and Stephen Boyd
3. [GeoCalcLib](http://worc4021.github.io/GeoCalcLib/) developed by Rainer Schaich.
4. [An algorithm for numerical computation of multivariate normal distribution values](http://www.math.wsu.edu/faculty/genz/software/matlab/qscmvnv.m) developed by Alan Genz (Distributed with SReachTools under its appropriate [license](https://sreachtools.github.io/docs/src/helperFunctions/qscmvnv/)).
5. [allcomb.m](https://www.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin)
developed by Jos van der Geest (Distributed with SReachTools under its appropriate [license](https://sreachtools.github.io/docs/src/helperFunctions/allcomb/#license)).

When available, we also use commercial toolboxes from MATLAB and Gurobi.
