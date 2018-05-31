---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
title: "Stochastic Reachability Toolbox"
---

## Quick start guide: installation, documentation, and examples

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
These are also catalouged in our [project page](https://abyvinod.github.io/SReachTools/examples/).
