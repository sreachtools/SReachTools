---
layout: docs
title: CwhSystemParameters.m
---

<ul class="doc-list">
    <li class="doc-list"><a href="#CwhSystemParameters">CwhSystemParameters</a></li>
    <ul class="doc-list">
        <li><a href="#CwhSystemParameters-CwhSystemParameters">Constructor</a></li>
        <li>Properties</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#CwhSystemParameters-prop-sampling_period">sampling_period</a></li>
            <li class="doc-list"><a href="#CwhSystemParameters-prop-orbital_radius">orbital_radius</a></li>
            <li class="doc-list"><a href="#CwhSystemParameters-prop-grav_constant">grav_constant</a></li>
            <li class="doc-list"><a href="#CwhSystemParameters-prop-celes_mass">celes_mass</a></li>
            <li class="doc-list"><a href="#CwhSystemParameters-prop-chief_mass">chief_mass</a></li>
            <li class="doc-list"><a href="#CwhSystemParameters-prop-grav_body">grav_body</a></li>
            <li class="doc-list"><a href="#CwhSystemParameters-prop-orbit_ang_vel">orbit_ang_vel</a></li>
            <li class="doc-list"><a href="#CwhSystemParameters-prop-disc_orbit_dist">disc_orbit_dist</a></li>
        </ul>
        <li>Methods</li>
        <ul class="doc-list">
        </ul>
    </ul>
</ul>

{:#CwhSystemParameters}
### CwhSystemParameters
```
  A MATLAB class to store/retrive the default parameters used in CWH dynamics
  with modifications, if any
  =============================================================================
 
  Notes:
  ------
  * This code and the parameters were obtained from Lesser's repeatability code
    for the 2013 CDC paper.
  * The default parameters (with their Names to specify changes) are:
        sampling_period : sampling period        = 20 s
        orbital_radius  : orbital radius         = 850 + 6378.1 m
        grav_constant   : gravitational constant = 6.673e-11
        celes_mass      : celestial body mass    = 5.9472e24 kg
        chief_mass      : chief mass             = 300 kg
  * Along with these parameters, the class provides these parameters that are
    computed using the above parameters
        grav_body       : gravitational body           = grav_constant *
                                                            celes_mass / 1e6
        orbit_ang_vel   : orbital angular velocity     = sqrt(grav_body /
                                                               orbital_radius^3)
        disc_orbit_dist : discretized orbital distance = orbit_ang_vel *
                                                            sampling_period rad
 
  =============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  

    Reference page in Doc Center
       doc CwhSystemParameters

```

{:#CwhSystemParameters-CwhSystemParameters}
### Constructor
```
  A MATLAB class to store/retrive the default parameters used in CWH dynamics
  with modifications, if any
  =============================================================================
 
  Notes:
  ------
  * This code and the parameters were obtained from Lesser's repeatability code
    for the 2013 CDC paper.
  * The default parameters (with their Names to specify changes) are:
        sampling_period : sampling period        = 20 s
        orbital_radius  : orbital radius         = 850 + 6378.1 m
        grav_constant   : gravitational constant = 6.673e-11
        celes_mass      : celestial body mass    = 5.9472e24 kg
        chief_mass      : chief mass             = 300 kg
  * Along with these parameters, the class provides these parameters that are
    computed using the above parameters
        grav_body       : gravitational body           = grav_constant *
                                                            celes_mass / 1e6
        orbit_ang_vel   : orbital angular velocity     = sqrt(grav_body /
                                                               orbital_radius^3)
        disc_orbit_dist : discretized orbital distance = orbit_ang_vel *
                                                            sampling_period rad
 
  =============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  

    Reference page in Doc Center
       doc CwhSystemParameters

```

### Property: sampling_period
{:#CwhSystemParameters-prop-sampling_period}
```
  sampling period in sec
  Default value: 20
```

### Property: orbital_radius
{:#CwhSystemParameters-prop-orbital_radius}
```
  Orbital radius in m
  Default value: 7228.1
```

### Property: grav_constant
{:#CwhSystemParameters-prop-grav_constant}
```
  Universal gravitational constant  
  Default value: 6.673e-11
```

### Property: celes_mass
{:#CwhSystemParameters-prop-celes_mass}
```
  Mass of the celestial body (default is Earth)
  Default value: 5.9472e24
```

### Property: chief_mass
{:#CwhSystemParameters-prop-chief_mass}
```
  Mass of the chief kg
  Default value: 300
```

### Property: grav_body
{:#CwhSystemParameters-prop-grav_body}
```
  Gravitation constant for the pull of the celestial body 
  (default Earth)
  Set via the equation
    grav_body = grav_constant * celes_mass / 1000^3;
```

### Property: orbit_ang_vel
{:#CwhSystemParameters-prop-orbit_ang_vel}
```
  Angular velocity in the orbit
  Set via the equatoin
    orbit_ang_vel = sqrt(grav_body / orbital_radius^3);
```

### Property: disc_orbit_dist
{:#CwhSystemParameters-prop-disc_orbit_dist}
```
  Discretized orbital distance
  Set via the equation
     disc_orbit_dist = orbit_ang_vel * sampling_period;   
```

