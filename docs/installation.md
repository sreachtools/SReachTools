---
layout: page
title: "Dependencies & Installation"
---

<a name="posttop"></a>

- [Dependencies](#dependencies)
- [Installation](#installation)

### Dependencies

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

    {% include important-note.html content="GeoCalcLib currently works only Unix
    and MAC OS.  SReachTools will gracefully switch back to CDD, if installation
    of GeoCalcLib is not correct." %}

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

[Go to top](#posttop)

### Installation

1. Install the necessary dependencies listed above
1. Clone the SReachTools repository (or download the latest zip file from
   [Releases](https://github.com/unm-hscl/SReachTools/releases))
1. Change the MATLAB current working directory to where SReachTools was
   downloaded. 
   {% include important-note.html content="Do not add the SReachTools folder to the path manually." %}
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

[Go to top](#posttop)
