# Stochastic reachability toolbox (SReachTools)

SReachTools is a MATLAB toolbox to tackle various problems in stochastic reachability.

This is an area of **active research**, and this toolbox will attempt to cater certain classes of problems. 
We aim to support the following problems:
 - **Stochastic reach-avoid problem** (guaranteeing safety for stochastic
   systems)
    - Lagrangian methods-based underapproximation
    - Fourier transforms-based underapproximation
 - **Forward stochastic reachability** (characterizing the stochasticity of the
   state at a future time of interest)

We will add more features like stochastic obstacle avoidance using reachability and closed-loop controller synthesis.
Do check our [project blog](https://abyvinod.github.io/SReachTools/blog/) for updates!

## Installation, documentation, and examples

### Dependencies

You can skip installing the dependencies marked **optional**.
This will disable some of the features of SReachTools.

1. MATLAB (>2017a)
    1. Toolboxes
        1. MATLAB's Global Optimization Toolbox (**Optional**)
        1. MATLAB's Statistics and Machine Learning Toolbox (**Optional**)
        1. MATLAB's Control System Toolbox (**Optional**)
1. MPT3 ([http://people.ee.ethz.ch/~mpt/3/](http://people.ee.ethz.ch/~mpt/3/))
    1. Do an automatic install using a MATLAB script [install_mpt3.m](http://control.ee.ethz.ch/~mpt/3/Main/Installation?action=download&upname=install_mpt3.m) provided by MPT3.
1. CVX ([http://cvxr.com/cvx/](http://cvxr.com/cvx/)) (**Optional**)

### Installation

1. Install the necessary dependencies (MATLAB and MPT3 are a must)
1. Clone the *SReachTools* repository (or download the zip file)
1. Run `srtinit -v -t` in MATLAB to add the toolbox to the paths, visualize the steps, and test the installation.  
   - You can add `cd <path_to_sreachtools_repo>;srtinit` to your MATLAB's `startup.m` to automatically have this done in future.

### Examples

See `examples/*.pdf` for the PDF version of various examples run using SReachTools.
These are also catalouged in our [project webpage](https://abyvinod.github.io/SReachTools/examples/). 

## Contributions

SReachTools is an open-source MATLAB toolboxes. We welcome feedback, issues, bug-fixes, and other enhancements. 
Please see [our contribution guidelines](./CONTRIBUTING.md) of this project to know more on how you can help us.

## License

See [LICENSE](./LICENSE).

## Credits

The authors of this toolbox are [Abraham P. Vinod](http://www.unm.edu/~abyvinod/) and [Joseph D.  Gleason](http://www.unm.edu/~gleasonj/). 
Please cite their [relevant papers](https://scholar.google.com/citations?user=yb5Z7AwAAAAJ&hl=en) when using the toolbox. 
The authors are PhD advisees of [Prof. Meeko Oishi](http://www.unm.edu/~oishi/).
