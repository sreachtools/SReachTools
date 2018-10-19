---
layout: page
permalink: /examples/
title: Examples
---

In these set of examples, we apply SReachTools to perform forward (prediction of the stochasticity of the state in future) and backward (verification and controller synthesis) stochastic reachability in various stochastic LTI and LTV systems.

Please feel free to add requests for more examples at [https://github.com/unm-hscl/SReachTools/issues](https://github.com/unm-hscl/SReachTools/issues).
They are available in the `examples` folder of the project. 

{% include important-note.html content="HTML examples are not mobile friendly. We're working on it." %}

1. **Stochastic reachability of a target tube using dynamic programming**:  Demonstrates the use of `SReachDynProg`.
    1. **Stochastic double integrator**: For a stochastic double integrator, we compute a grid-based dynamic programming solution to the stochastic reachability of different target tubes listed below. [[PDF](https://github.com/unm-hscl/SReachTools/raw/master/examples/doubleIntegratorDynamicProgramming.pdf)] [[HTML](doubleIntegratorDynamicProgramming.html)]
        - Stochastic viability problem
        - Stochastic reach-avoid problem (terminal hitting time)
        - Stochastic reachability of a target tube (a collection of time-varying safe sets)
1. **Underapproximative verification**: Demonstrates the use of `SReachSet` and `SReachPoint` on a host of examples. 
    - *Verification problem*: We can use `SReachSet` to compute the set of initial states and their associated admissible controllers such that the probability of the system evolving in a pre-specified target tube (a collection of time-varying sets) is above a user-specified threshold. 
    - *Controller synthesis problem*: We can use `SReachPoint` to compute the optimal (open-loop or affine) controller that maximizes the safety probability, the likelihood of staying with the target tube.
    - For both of these problems, we will validate the optimal controllers using \\(10^5\\) Monte-Carlo simulations.
    - Examples: 
        1. **Spacecraft rendezvous and docking problem**: In a satellite rendezvous problem, we wish to steer a spacecraft (deputy) towards another spacecraft on the same elliptical orbit (chief), while staying within the line-of-sight cone and respecting actuation limits. The relative dynamics of the deputy with respect to the chief spacecraft is a LTI system ([Clohessy-Wiltshire-Hill dynamics](https://en.wikipedia.org/wiki/Clohessy-Wiltshire_equations)) with additive stochastic disturbance. 
            - Application of `SReachPoint` from a given state. [[PDF](https://github.com/unm-hscl/SReachTools/raw/master/examples/publish/cwhSReachPointDemo.pdf)] [[HTML](cwhSReachPointDemo.html)]
            - Application of `SReachSet` for a probability threshold. [[PDF](https://github.com/unm-hscl/SReachTools/raw/master/examples/publish/cwhSReachPointDemo.pdf)] [[HTML](cwhSReachPointDemo.html)]
        1. **Dubin's vehicle with a known turn rate**: Given a [Dubin's car](https://en.wikipedia.org/wiki/Dubins_path) with a known turning rate and additive disturbance in position, we wish to ascertain the set of initial states and its associated controllers from which the car can be steered to lie within a target tube. 
            - Application of `SReachPoint` from a given state. [[PDF](https://github.com/unm-hscl/SReachTools/raw/master/examples/publish/cwhSReachPointDemo.pdf)] [[HTML](cwhSReachPointDemo.html)]
            - Application of `SReachSet` for a probability threshold. [[PDF](https://github.com/unm-hscl/SReachTools/raw/master/examples/publish/cwhSReachPointDemo.pdf)] [[HTML](cwhSReachPointDemo.html)]
        1. **Automated anesthesia delivery system**: Given the dynamics for the sedation levels of a patient, we wish to design an automated anesthesia delivery system that maintains the probability of the depth of hypnosis to stay within pre-specified bounds, above 0.99. 
            - Application of `SReachSet` for a probability threshold. [[PDF](https://github.com/unm-hscl/SReachTools/raw/master/examples/publish/cwhSReachPointDemo.pdf)] [[HTML](cwhSReachPointDemo.html)]
1. **Forward stochastic reachability**: Computes the probability of the state of the spacecraft lying in a set at future time of interest with/without safety requirements. We also demonstrate how SReachTools can handle a deterministic and a stochastic initial state. [[PDF](https://github.com/unm-hscl/SReachTools/raw/master/examples/publish/forwardStochasticReachCWH.pdf)] [[HTML](forwardStochasticReachCWH.html)] 
<!-- Add {:target="_blank"} if it is desired that the page opens in a new window.-->


{% comment %}

like:

* the relative dynamics of a spacecraft in space to a docking station,
* the stochastic double integrator, and
* the automated anesthesia delivery system.

We used MATLAB's LiveEditor to write the following examples. 

{% endcomment %}
