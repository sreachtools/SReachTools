---
layout: post
title:  "Release of v0.2"
date:   2018-05-31 13:40:00 -0600
categories: jekyll update
author: "Joseph D. Gleason & Abraham P. Vinod"
---

SReachTools Toolbox [v0.2](https://github.com/unm-hscl/SReachTools/tree/v0.2) is now out! The SReachTools toolbox is a set of MATLAB codes to facilitate stochastic reachability computations. The toolbox currently supports reachability computations with several methods:
* Dynamic programming computation
* Fourier transformation approximations
* Lagrangian (set-based) approximations
We currently support forward and backward stochastic reachability of *stochastic LTI systems*.

See the quick-start guide below for installation instructions.

## Quick start guide: installation and examples

### Dependencies

You can skip installing the dependencies marked **optional**.
This will disable some of the features of SReachTools.

1. MATLAB (>2017a)
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

### Examples

See our [examples page](https://unm-hscl.github.io/SReachTools/examples/).

## What are we working on next?

- Extension to LTV systems
- Support for chance-constrained and particle filter-based verification
- Adding online API documentation
