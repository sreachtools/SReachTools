---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
title:
---

### Dependencies

* MATLAB (>2017a)
* MPT3 (http://people.ee.ethz.ch/~mpt/3/)
* MATLAB's Global Optimization Toolbox (Optional)
* MATLAB's Statistics and Machine Learning Toolbox (Optional)
* MATLAB's Control System Toolbox (Optional)
* CVX (http://cvxr.com/cvx/) (Optional)

### Installation

1. Install the necessary dependencies
1. Clone the SReachTools repository (or download the [zip file](https://github.com/abyvinod/SReachTools/archive/master.zip))
        
        git clone https://github.com/abyvinod/SReachTools.git

1. Run `sreachinit -v -t` in MATLAB to add the toolbox to the paths, visualize the steps, and test the installation.
* You can add `p = pwd;cd '/path/to/sreachtools';srtinit;cd(p);` to your MATLAB's `startup.m` to automatically have this done in future.
