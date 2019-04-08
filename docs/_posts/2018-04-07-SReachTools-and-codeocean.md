---
layout: post
title:  "Using SReachTools on CodeOcean"
date:   2019-04-07 08:00:00 -0600
categories: repeatability
author: "Abraham P. Vinod"
---

CodeOcean is a  cloud-based computational reproducibility platform. This post
shows how you can utilize SReachTools in your CodeOcean capsule for
repeatability. We list few examples that show SReachTools setup for CodeOcean.
Then, we detail the procedure to set up a new capsule.

### Examples

1. [Dubins' car: `SReachSet` demonstration](https://codeocean.com/capsule/9849812/tree) --- <2 min runtime.
1. [ARCH 2019: Automated anesthesia delivery](https://codeocean.com/capsule/81dbdd67-83a4-4b7a-87eb-fbfc1fe72ef2/tree?ID=1763713383dc4aadaf77c4a8c9085b7f) --- <1 min runtime.
1. [ARCH 2019: Building automation system](https://codeocean.com/capsule/ddfa8988-7061-4b9c-8cc8-0d3393a6bf02/tree?ID=64794d9f443f4fbc9ba1710cbd0e72c7) --- <1 min runtime.
<!--1. Spacecraft rendezvous problem: `SReachPoint` demonstration-->

### Detailed procedure

We choose MATLAB 2017b and the follow the installation instructions given in
[the installation page](/installation/) to setup `GMP`, `GeoCalcLib`, `YALMIP`,
`CVX`, `MPT3`, and finally `SReachTools`.

1. Create a new capsule
1. Choose MATLAB 2017b environment
    1. At the time of writing, this was the latest MATLAB environment available.
1. In additional packages, add the following:
    1. `bzip2` (tested with 1.0.6-8) --- for unzipping MATLAB toolboxes
    1. `build-essential` (tested with 12.1ubuntu2) --- for
       [GeoCalcLib](https://github.com/worc4021/GeoCalcLib) installation (`gcc`
       and `make`)
    1. `m4` (tested with 1.4.17-5) --- for
       [GeoCalcLib](https://github.com/worc4021/GeoCalcLib) installation
1. Next, copy the following text into `post-install script` section:
    ``` bash
    #!/usr/bin/env bash
    
    # tells script to exit, if any line in the script fails (CodeOcean default)
    set -e
    
    ### Install GMP
    # Requires packages build-essentials for gcc and make, and package m4 (don't
    # really know why)
    curl -sL https://gmplib.org/download/gmp/gmp-6.1.2.tar.bz2 | tar xj
    echo 'Downloaded GMP'
    cd gmp-6.1.2
    # Setup gmp
    ./configure 
    make 
    make install
    cd ..
    echo 'Installed GMP'
    
    ### GeoCalcLib (Download and setup)
    # Requires installation of gmp and packages build-essentials for gcc and make 
    curl -sL https://github.com/worc4021/GeoCalcLib/archive/master.zip --output GeoCalcLib.zip
    unzip -qq GeoCalcLib.zip
    mv GeoCalcLib-master GeoCalcLib
    rm GeoCalcLib.zip
    echo 'Downloaded GeoCalcLib'
    # Setup GeoCalcLib
    cd GeoCalcLib
    mkdir mexfiles
    echo 'MATLABROOT = /MATLAB' > User.make
    echo 'INSTALLDIR = /GeoCalcLib/mexfiles' >> User.make
    make
    cd ..
    echo 'Installed GeoCalcLib'
    
    ### YALMIP (Download only)
    # MATLAB simply adds the files to the path to complete the setup
    YALMIP_RELEASE=R20181012
    curl -sL https://github.com/yalmip/YALMIP/archive/$YALMIP_RELEASE.tar.gz | tar xz
    mv YALMIP-$YALMIP_RELEASE YALMIP
    echo 'Downloaded YALMIP'
    
    ### CVX (Download only)
    # We call cvx_setup in MATLAB to complete the setup
    curl -sL http://web.cvxr.com/cvx/cvx-a64.tar.gz | tar zx
    echo 'Downloaded CVX'
    
    ### SReachTools (Download only)
    # Call srtinit in MATLAB to complete the setup
    ## OPTION 1: Fetch a tagged release
    #
    # SREACHTOOLS_RELEASE=1.1
    # curl -sL https://abyvinod.github.io/code/SReachTools-v$SREACHTOOLS_RELEASE.tar.gz | tar xz
    # echo 'Downloaded SReachTools (stable)'
    # mv SReachTools-$SReachTools SReachTools
    #
    ## OPTION 2: Fetch nightly
    #
    curl -sL https://github.com/unm-hscl/SReachTools/archive/master.zip --output SReachTools.zip
    echo 'Downloaded SReachTools (nightly)'
    unzip -qq SReachTools.zip
    mv SReachTools-master SReachTools
    rm SReachTools.zip
    
    ### Setup MATLAB env for GeoCalcLib, YALMIP, CVX, MPT3, and SReachTools
    # Trailing \ implies newline in bash. MATLAB executes the commands in the quotes
    # Make sure each line ends with ;\ to avoid MATLAB throwing errors
    # Use GPLK instead of LCP to avoid wierd shifting | Don't do mpt_init again
    # 1. Add GeoCalcLib to path
    # 2. Add YALMIP to the path
    # 3. Setup CVX
    # 4. Fetch and install MPT3 using tbxmanager; Initialize MPT3
    # 5. Setup SReachTools
    # Even though, we do not recommend using savepath, we have to do it here due to
    # CodeOcean's setup
    matlab -nodisplay -r "\
    addpath('/GeoCalcLib/mexfiles');\
    disp('Installed GeoCalcLib');\
    addpath(genpath('/YALMIP'));\
    disp('Installed YALMIP');\
    cd('/cvx');\
    evalc('cvx_setup()');\
    disp('Installed CVX (Standard bundle)');\
    mkdir('/tbxmanager');\
    cd('/tbxmanager');\
    urlwrite('http://www.tbxmanager.com/tbxmanager.m', 'tbxmanager.m');\
    a=evalc('tbxmanager');\
    disp('Installed tbxmanager');\
    evalc('tbxmanager install mpt mptdoc cddmex fourier glpkmex hysdel lcp sedumi espresso');\
    evalc('mpt_init');\
    a=mptopt('lpsolver','glpk','qpsolver','quadprog');\
    disp('Installed MPT3');\
    cd('/SReachTools');\
    srtinit;\
    disp('Installed SReachTools');\
    savepath;"
    ```
1. Create a file `run.sh` in the `code` folder. 
    ```
    #!/bin/bash
    matlab -nodisplay -nosoftwareopengl -r "dubinsSReachSetGauss;"
    ```
    This runs a default example in `SReachTools`. If you wish to run a separate
    MATLAB script `main.m`, you can upload it the `code` folder, and overwrite
    `run.sh` with
    ```
    #!/bin/bash
    matlab -nodisplay -nosoftwareopengl -r "main;"
    ```

#### Notes

- The first run will take about 3 minutes for setting up of the environment.
  However, if the `post-install script` is left untouched, then subsequent runs
  will be significantly faster.
- CodeOcean, by default, does not save figures. Use the following lines to save
  a figure as a `png` or `fig` file.
  ```
  saveas(gcf, '../results/FILENAME.png');
  saveas(gcf, '../results/FILENAME.fig');
  ```

