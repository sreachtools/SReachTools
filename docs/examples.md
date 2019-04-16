---
layout: page
permalink: /examples/
title: Examples
---

We apply SReachTools to perform the following analysis on various stochastic LTI
and LTV systems:
1. [backward stochastic reachability](#backward-stochastic-reachability)
   (verification and controller synthesis)
1. [forward stochastic reachability](#forward-stochastic-reachability)
   (prediction of the stochasticity of the state in future) 

Please feel free to add requests for more examples at
[https://github.com/unm-hscl/SReachTools/issues](https://github.com/unm-hscl/SReachTools/issues).
The code for `XYZ.html` is available in `examples/XYZ.m`. These html pages were
created using MATLAB's publish command (see `examples/publish_examples.m`).

{% include important-note.html content="HTML examples are not mobile friendly. 
We're working on it." %}

## Backward stochastic reachability

### (Under)approximative verification and controller synthesis
Using convex optimization, Fourier transforms, mixed-integer linear programming,
and computation geometry, we can perform backward stochastic reachability. Here,
we define safety by requiring the state to stay within a pre-specified target
tube (a collection of time-varying safe sets).

- *Verification problem*: We can use `SReachSet` to compute the set of
  initial states and their associated admissible controllers such that the
  probability of safety of the system is above a user-specified threshold. 
- *Controller synthesis problem*: We can use `SReachPoint` to compute the
  optimal (open-loop or affine) controller that maximizes the safety
  probability, the likelihood of staying with the target tube.
- See [the table in the documentation page]({{ '/docs/#features-of-sreachtools' | relative_url }}) 
  for how SReachTools can solve these problems and the available approximation
      guarantees.
- In all of these examples, we will validate the optimal controllers using
  Monte-Carlo simulations.

We provide the following demonstration examples:
1. **Spacecraft rendezvous problem**: In a satellite rendezvous problem, we wish
   to steer a spacecraft (deputy) towards another spacecraft on the same
   elliptical orbit (chief), while staying within the line-of-sight cone and
   respecting actuation limits. The relative dynamics of the deputy with respect
   to the chief spacecraft is a LTI system ([Clohessy-Wiltshire-Hill
   dynamics](https://en.wikipedia.org/wiki/Clohessy-Wiltshire_equations)) with
   additive Gaussian disturbance. 
    - Application of `SReachPoint` from a given state
      [[cwhSReachPoint.html](./publish/cwhSReachPoint.html)]
    - Application of `SReachSet` for a probability threshold
      [[cwhSReachSet.html](./publish/cwhSReachSet.html)]
1. **Dubin's vehicle with a known turn rate**: Given a [Dubin's car](https://en.wikipedia.org/wiki/Dubins_path) with a known turning rate and additive disturbance in position, we wish to ascertain the set of initial states and its associated controllers from which the car can be steered to lie within a target tube. 
    - Application of `SReachPoint` from a given state
      [[dubinsSReachPointGauss.html](./publish/dubinsSReachPointGauss.html)]
    - Application of `SReachSet` for a probability threshold [![Open in Code
      Ocean](https://codeocean.com/codeocean-assets/badge/open-in-code-ocean.svg)](https://doi.org/10.24433/CO.9849812.v1)
      [[dubinsSReachSetGauss.html](./publish/dubinsSReachSetGauss.html)] 
1. **Automated anesthesia delivery system**: Given the dynamics for the sedation
   levels of a patient, we wish to design an automated anesthesia delivery
   system that maintains the probability of the depth of hypnosis to stay within
   pre-specified bounds, above 0.99. 
   [![Open in Code
      Ocean](https://codeocean.com/codeocean-assets/badge/open-in-code-ocean.svg)](https://doi.org/10.24433/CO.3325937.v1)
    - Application of `SReachSet` for a probability threshold. [[automatedAnesthesiaDeliveryARCH2018.html](./publish/automatedAnesthesiaDeliveryARCH2018.html)]
   <!--This script serves as a repeatability-->
   <!--evaluation of SReachTools code in [Abate et. al, ARCH-->
   <!--2018](https://doi.org/10.29007/7ks7).-->
1. **Lagrangian approximations for the stochastic double integrator**: For a
   stochastic double integrator, we compute the over- and under-approximations
   for the stochastic viability problem, and compare it with the dynamic
       programming. 
    - Application of `SReachSet` for a probability threshold via Lagrangian options. [[dIntSReachSetLag.html](./publish/dIntSReachSetLag.html)]
1. **Building Automation System**: For assuring probabilistic safety in a
   multi-room building automation system. [![Open in Code
   Ocean](https://codeocean.com/codeocean-assets/badge/open-in-code-ocean.svg)](
https://doi.org/10.24433/CO.8093142.v1)

{% include important-note.html content="Examples 1,2, and 3 can not be
implemented using dynamic programming! <\br> SReachTools provides verification of
these systems for the first time." %}
### Dynamic programming

**Stochastic double integrator**: For a stochastic double integrator, we compute
a grid-based dynamic programming solution to the stochastic reachability of
different target tubes listed below.  

- Stochastic viability problem
- Stochastic reach-avoid problem (terminal hitting time)
- Stochastic reachability of a target tube (a collection of time-varying safe sets)

The application of `SReachDyn` is given in
[[dIntSReachDyn.html](./publish/dIntSReachDyn.html)].

## Forward stochastic reachability

1. **Spacecraft rendezvous problem**: Computes the probability of the state of
   the spacecraft lying in a set at future time of interest with/without safety
   requirements. We also demonstrate how SReachTools can handle a deterministic
   and a stochastic initial state. The dynamics of the spacecraft is assumed to
   be  a LTI system ([Clohessy-Wiltshire-Hill
   dynamics](https://en.wikipedia.org/wiki/Clohessy-Wiltshire_equations)) with
   additive Gaussian disturbance, and the controller is given by a LQR
   controller.
    - Application of `SReachFwd`
      [[cwhSReachFwd.html](./publish/cwhSReachFwd.html)].

<!-- Add {:target="_blank"} if it is desired that the page opens in a new window.-->
