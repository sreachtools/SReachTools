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
  input handling

    Reference page in Doc Center
       doc CwhSystemParameters

```

{:#CwhSystemParameters-CwhSystemParameters}
### Constructor
```
  input handling

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
  Gravitation constant for the pull of the celestial body (default Earth)
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

