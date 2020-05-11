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

We will denote MATLAB's command prompt by `>>`, while the system command prompt
by `$ `.

### Dependencies

External dependencies of SReachTools are:
1. MATLAB's Statistics and Machine Learning Toolbox
1. MPT3 ([https://www.mpt3.org/](https://www.mpt3.org/))
1. CVX ([http://cvxr.com/cvx/](http://cvxr.com/cvx/))
1. MATLAB's Global Optimization Toolbox (**optional**)
1. MATLAB's Optimization Toolbox (**optional**)
1. GeoCalcLib
   ([https://github.com/sreachtools/GeoCalcLib](https://github.com/sreachtools/GeoCalcLib))
   (**optional**)
1. GUROBI (**optional**)
1. MOSEK (**optional**)

#### Essential dependencies

These dependencies are platform-independent, and have been tested in Windows,
MacOS, and Linux.

1. **MATLAB (>2017a)** with MATLAB's Statistics and Machine Learning Toolbox
1. **MPT3** ([https://www.mpt3.org/](https://www.mpt3.org/)) --- for polytopic
   computational geometry
    1. Copy the MATLAB script [install_mpt3.m](https://www.mpt3.org/Main/Installation?action=download&upname=install_mpt3.m)
       provided by MPT3 from the browser to your local computer.
    1. Download and install MPT3 by running the following command in 
       MATLAB's command prompt (after changing
       directory to the folder containing `install_mpt3.m`)
       ```
       >> install_mpt3
       ```
    1. Add `cd <PATH-TO-TBXMANAGER>;tbxmanager restorepath`
       to your MATLAB `startup` script for the MPT3
       installation to persist across MATLAB runs.
1. **CVX** ([http://cvxr.com/cvx/](http://cvxr.com/cvx/)) --- for parsing convex
   and mixed-integer programs. Use the *Standard bundles, including Gurobi
   and/or MOSEK*,  even if you do not plan to use Gurobi or MOSEK. CVX does not
   require any additional licenses to work with GUROBI or MOSEK in an academic
   setting.
   1. Download the zip file from
      [http://cvxr.com/cvx/download/](http://cvxr.com/cvx/download/), and
      extract it to the `cvx` folder.
   1. Install CVX by running the following command in MATLAB's command prompt 
      (after changing the current working directory to the `cvx` folder)
      ```
      >> cvx_setup
      ```
   1. Add `cd <PATH-TO-CVX>;cvx_setup` to your MATLAB `startup` script for the 
      CVX installation to persist across MATLAB runs.
   1. Other notes:
      - Detailed installation instructions are given in
        [http://cvxr.com/cvx/download/](http://cvxr.com/cvx/download/).
      - SDPT3 (the default backend solver of CVX) performs reasonably well with
        CVX, when compared to MOSEK, and significantly poorly when compared to
        GUROBI in the tested examples and CVX v2.1 version. See Step 5 for
        instructions in installing external solvers for SReachTools.


#### Optional dependencies

These dependencies will allow non-essential features of SReachTools or enable
superior performance.

1. MATLAB's Global Optimization Toolbox and Optimization Toolbox
   - Affects `genzps-open` options in `SReachPoint` and `SReachSet` computation.
1. [GeoCalcLib](https://github.com/worc4021/GeoCalcLib) --- a MATLAB interface
   to Avis's [LRS vertex-facet enumeration
   library](http://cgm.cs.mcgill.ca/~avis/C/lrs.html). In empirical tests, we
   found LRS to be a superior alternative to MPT's preferred approach for
   vertex-facet enumeration,
   [CDD](https://www.inf.ethz.ch/personal/fukudak/cdd_home/index.html).  
   - Affects `lag-over` and `lag-under` options in `SReachSet` computation.
   - See
   [https://github.com/sreachtools/GeoCalcLib](https://github.com/sreachtools/GeoCalcLib)
   for detailed installation instructions.
   - :warning: GeoCalcLib currently works only in Unix and MAC OS. 
1. External solvers --- **GUROBI** and/or **MOSEK**.
    1. Why do we need to install external solvers?
        - Mixed-integer programming enabled by GUROBI or MOSEK is required for
          particle-based approaches in `SReachPoint`.
            - Affects `voronoi-open` and `particle-open` options in
              `SReachPoint` computation.
        - External solvers are typically more numerically robust and
          computationally faster than free solvers like SDPT3 that come with
          CVX.
            - Affects overall computation speed and solution quality.
    1. **GUROBI** ([http://www.gurobi.com/](http://www.gurobi.com/)):  A backend
       solver for CVX that also enables mixed-integer programming associated
       with particle-based approaches, apart from convex programming. MPT3 +
       GUROBI also provides robust polyhedral computation. (:warning: Currently
       facing issues).  
       1. GUROBI offers free academic license, which can be requested at
          [http://www.gurobi.com/registration/download-reg](http://www.gurobi.com/registration/download-reg).  
       1. **MPT3 + GUROBI**: Requires the external installation.
           1. Obtain a copy of GUROBI Optimizer from
              [http://www.gurobi.com/](http://www.gurobi.com/) (Requires signing
              up)
           1. Unzip the installation to desired folder.
           1. Generate the license file by running the following command in the 
              Unix command prompt (after changing the current working directory 
              to the `<PATH-TO-GUROBI-HOME>/gurobi902/<OS>/bin/`)
              ```
              $ grbgetkey OUTPUT_OF_GUROBI_LICENSE_REQUEST
              ```               
           1. Setup GUROBI by running the following command in MATLAB's command 
              prompt (after changing the current working directory to 
              `<PATH-TO-GUROBI-HOME>/gurobi902/<OS>/matlab/`)
              ```
              >> gurobi_setup
              ```
           1. Add `GRB_LICENSE_FILE` environment variable that has the location
              of the `gurobi.lic` file for MPT3 to detect GUROBI. Alternatively,
              add the following command in `startup.m`
              ```
              >> setenv('GRB_LICENSE_FILE','/home/ubuntu/gurobi.lic')
              ```
           1. Update MPT3 with GUROBI by running the following command in MATLAB's 
              command prompt
              ```
              >> mpt_init
              ```
       1. **CVX + GUROBI**: The current build of CVX v2.2 does not play well
            with GUROBI v9.0.2. (:warning: Currently facing issues)
          - Previously, CVX v2.1 worked with GUROBI v7.5.2 successfully. We will
            update these instructions once CVX + GUROBI becomes a viable option
            again. 
          - See the following CVX forum posts for more details:
            1. [http://ask.cvxr.com/t/cvx-with-gurobi-error-warning/7072](http://ask.cvxr.com/t/cvx-with-gurobi-error-warning/7072)
            1. [http://ask.cvxr.com/t/mismangement-of-matlab-path-in-cvxv2-2-and-external-gurobi-v9-0-2/7346/3](http://ask.cvxr.com/t/mismangement-of-matlab-path-in-cvxv2-2-and-external-gurobi-v9-0-2/7346/3)
          - We found CVX + GUROBI to be faster than SDPT3 and MOSEK 
    1. **MOSEK**  ([https://www.mosek.com/](https://www.mosek.com/)): A backend
       solver for CVX that also enables mixed-integer programming associated
       with particle-based approaches, apart from convex programming.
       1. MOSEK offers free academic license, which can be requested at
          [https://www.mosek.com/license/request/](https://www.mosek.com/license/request/).
       1. **CVX + MOSEK**: 
          1. Save the license file obtained via email in your home folder in a
             folder named `mosek`. See
             [https://docs.mosek.com/9.2/licensing/quickstart.html](https://docs.mosek.com/9.2/licensing/quickstart.html)
             for more details.
          1. To enable MOSEK bundled with CVX, run the following command in
             MATLAB command prompt
             ```
             >> cvx_setup
             ```
          1. See
             [http://web.cvxr.com/cvx/doc/mosek.html](http://web.cvxr.com/cvx/doc/mosek.html)
             for more details.
       1. **MPT3 + MOSEK**: Unfortunately, MPT3 does not support MOSEK. See
          [https://www.mpt3.org/Main/FAQ](https://www.mpt3.org/Main/FAQ).
       1. In empirical tests, SDPT3 and MOSEK have demonstrated similar
          performance in computation time and solution quality. We found MOSEK
          to be slower than GUROBI.

### Installation

1. Install the necessary dependencies listed above
1. Clone the SReachTools repository (or download the latest zip file from
   [Releases](https://github.com/sreachtools/SReachTools/releases))
   ```
   $ git clone https://github.com/sreachtools/SReachTools
   ```
   > :warning: Please do not add the folder to the MATLAB path manually.
1. Install SReachTools by running the following command in MATLAB command prompt 
   (after changing the current working directory in MATLAB to `SReachTools` 
   folder)
   ```
   >> srtinit
   ```
   - (**Optional**) Additional steps:
       - Run `srtinit -t` to run all the unit tests.
       - Run `srtinit -v` to visualize the steps the changes to the path and
         check for recommended dependencies.  
       - Run `srtinit -x` to remove functions of SReachTools from MATLAB's path
         after use.  
1. You can add `cd <PATH-TO-SREACHTOOLS>;srtinit` to your MATLAB's `startup.m`
   to automatically have this done in future.
   
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

> A. P. Vinod, J. D. Gleason, and M. M. K. Oishi.  "SReachTools: A MATLAB Stochastic Reachability Toolbox," In Proceedings of the International Conference on Hybrid Systems: Computation and Control, Montreal, Canada, April 16--18, 2019.  https://sreachtools.github.io.
    
BibTeX entry for use in LaTeX with `\usepackage{url}`: 
```
@misc{SReachTools,
    author    = {Vinod, Abraham P. and Gleason, Joseph D. and Oishi, Meeko M. K.},
    title     = {{S}{R}each{T}ools: A {MATLAB} {S}tochastic {R}eachability {T}oolbox},
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
