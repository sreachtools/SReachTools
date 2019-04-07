---
layout: docs
title: Documentation
permalink: /docs/
---

<a name="posttop"></a>

- [List of functions in SReachTools](#list-of-functions-in-sreachtools)
- [Features of SReachTools](#features-of-sreachtools)

Do check out the [examples](../examples) for more detailed explanations of
various functionalities of the toolbox.

## List of functions in SReachTools

In MATLAB, you can use:

- `help FUNCTION_NAME` to understand the details of a function
- `methods(CLASS_NAME)` to understand the details of a class
- `<Ctrl+F1>` to get function hints for a given function

Click on any of the SReachTools functions listed below to learn more about it.
This information is obtained from their docstrings. 

<!-- DO NOT REMOVE: docs2md START FILE DUMP HERE -->

<ul class="doc-list">
    <li>src/</li>
    <ul class="doc-list">
        <li>classes/</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="src/classes/@LtiSystem/LtiSystem">@LtiSystem/</a></li>
            <li class="doc-list"><a href="src/classes/@LtvSystem/LtvSystem">@LtvSystem/</a></li>
            <li class="doc-list"><a href="src/classes/RandomVector">RandomVector.m</a></li>
            <li class="doc-list"><a href="src/classes/SReachEllipsoid">SReachEllipsoid.m</a></li>
            <li class="doc-list"><a href="src/classes/SReachLagController">SReachLagController.m</a></li>
            <li class="doc-list"><a href="src/classes/Tube">Tube.m</a></li>
        </ul>
        <li>exceptions/</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="src/exceptions/SrtBaseException">SrtBaseException.m</a></li>
            <li class="doc-list"><a href="src/exceptions/SrtDevError">SrtDevError.m</a></li>
            <li class="doc-list"><a href="src/exceptions/SrtInternalError">SrtInternalError.m</a></li>
            <li class="doc-list"><a href="src/exceptions/SrtInvalidArgsError">SrtInvalidArgsError.m</a></li>
            <li class="doc-list"><a href="src/exceptions/SrtRuntimeError">SrtRuntimeError.m</a></li>
            <li class="doc-list"><a href="src/exceptions/SrtSetupError">SrtSetupError.m</a></li>
            <li class="doc-list"><a href="src/exceptions/SrtTestError">SrtTestError.m</a></li>
        </ul>
        <li>helperFunctions/</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="src/helperFunctions/allcomb">allcomb.m </a> (Third-party code distributed with license information)</li>
            <li class="doc-list"><a href="src/helperFunctions/ellipsoidsFromMonteCarloSims">ellipsoidsFromMonteCarloSims.m</a></li>
            <li class="doc-list"><a href="src/helperFunctions/generateMonteCarloSims">generateMonteCarloSims.m</a></li>
            <li class="doc-list"><a href="src/helperFunctions/getBsetWithProb">getBsetWithProb.m</a></li>
            <li class="doc-list"><a href="src/helperFunctions/getSrtWarning">getSrtWarning.m</a></li>
            <li class="doc-list"><a href="src/helperFunctions/normalizeForParticleControl">normalizeForParticleControl.m</a></li>
            <li class="doc-list"><a href="src/helperFunctions/polytopesFromMonteCarloSims">polytopesFromMonteCarloSims.m</a></li>
            <li class="doc-list"><a href="src/helperFunctions/qscmvnv">qscmvnv.m </a> (Third-party code distributed with license information)</li>
            <li class="doc-list"><a href="src/helperFunctions/setSrtWarning">setSrtWarning.m</a></li>
            <li class="doc-list"><a href="src/helperFunctions/spreadPointsOnUnitSphere">spreadPointsOnUnitSphere.m</a></li>
        </ul>
        <li>modules/</li>
        <ul class="doc-list">
            <li>fwdStochReach/</li>
            <ul class="doc-list">
                <li class="doc-list"><a href="src/modules/fwdStochReach/SReachFwd">SReachFwd.m</a></li>
            </ul>
            <li>nonStochReach/</li>
            <ul class="doc-list">
                <li class="doc-list"><a href="src/modules/nonStochReach/getSReachLagOverapprox">getSReachLagOverapprox.m</a></li>
                <li class="doc-list"><a href="src/modules/nonStochReach/getSReachLagUnderapprox">getSReachLagUnderapprox.m</a></li>
            </ul>
            <li>stochReach/</li>
            <ul class="doc-list">
                <li>dynProg/</li>
                <ul class="doc-list">
                    <li class="doc-list"><a href="src/modules/stochReach/dynProg/SReachDynProg">SReachDynProg.m</a></li>
                    <li class="doc-list"><a href="src/modules/stochReach/dynProg/getDynProgLevelSets2D">getDynProgLevelSets2D.m</a></li>
                </ul>
                <li>pointBased/</li>
                <ul class="doc-list">
                    <li class="doc-list"><a href="src/modules/stochReach/pointBased/SReachPoint">SReachPoint.m</a></li>
                    <li class="doc-list"><a href="src/modules/stochReach/pointBased/SReachPointCcA">SReachPointCcA.m</a></li>
                    <li class="doc-list"><a href="src/modules/stochReach/pointBased/SReachPointCcAu">SReachPointCcAu.m</a></li>
                    <li class="doc-list"><a href="src/modules/stochReach/pointBased/SReachPointCcO">SReachPointCcO.m</a></li>
                    <li class="doc-list"><a href="src/modules/stochReach/pointBased/SReachPointGpO">SReachPointGpO.m</a></li>
                    <li class="doc-list"><a href="src/modules/stochReach/pointBased/SReachPointOptions">SReachPointOptions.m</a></li>
                    <li class="doc-list"><a href="src/modules/stochReach/pointBased/SReachPointPaO">SReachPointPaO.m</a></li>
                    <li class="doc-list"><a href="src/modules/stochReach/pointBased/SReachPointVoO">SReachPointVoO.m</a></li>
                    <li class="doc-list"><a href="src/modules/stochReach/pointBased/computeNormCdfInvOverApprox">computeNormCdfInvOverApprox.m</a></li>
                </ul>
                <li>setBased/</li>
                <ul class="doc-list">
                    <li class="doc-list"><a href="src/modules/stochReach/setBased/SReachSet">SReachSet.m</a></li>
                    <li class="doc-list"><a href="src/modules/stochReach/setBased/SReachSetCcO">SReachSetCcO.m</a></li>
                    <li class="doc-list"><a href="src/modules/stochReach/setBased/SReachSetGpO">SReachSetGpO.m</a></li>
                    <li class="doc-list"><a href="src/modules/stochReach/setBased/SReachSetLag">SReachSetLag.m</a></li>
                    <li class="doc-list"><a href="src/modules/stochReach/setBased/SReachSetLagBset">SReachSetLagBset.m</a></li>
                    <li class="doc-list"><a href="src/modules/stochReach/setBased/SReachSetOptions">SReachSetOptions.m</a></li>
                </ul>
            </ul>
        </ul>
        <li>systems/</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="src/systems/CwhSystemParameters">CwhSystemParameters.m</a></li>
            <li class="doc-list"><a href="src/systems/getChainOfIntegLtiSystem">getChainOfIntegLtiSystem.m</a></li>
            <li class="doc-list"><a href="src/systems/getCwhLtiSystem">getCwhLtiSystem.m</a></li>
            <li class="doc-list"><a href="src/systems/getDubinsCarLtv">getDubinsCarLtv.m</a></li>
        </ul>
    </ul>
</ul>

<!-- DO NOT REMOVE: docs2md END FILE DUMP HERE -->

[Go to top](#posttop)

## Features of SReachTools

The following table[^table_ack] summarizes the features in SReachTools.


|    Function   |   method-str  |                                                       Utility                                                       | Notes                                      |
|:-------------:|:---------------:|:-------------------------------------------------------------------------------------------------------------------:|--------------------------------------------|
| `SReachPoint` |                 |          **Approximation of the maximal reach  probability for a target tube from  a given initial state** [^TAC2018_verification]         | **Synthesize open-loop or affine disturbance feedback controllers** |
|               |  `chance-open`  |                                            Guaranteed underapproximation [^CDC2013_Lesser]<sup>,</sup>  [^CDC2019_chance]                                             | Open-loop                                  |
|               |  `genzps-open`  |                                  Approximate up to \\( \\epsilon\_\\mathrm{genz}\\), a user-specified quadrature error tolerance [^CSSL2017_genzps]                                 | Open-loop                                  |
|               | `particle-open` |                        Approximate with quality proportional  to the number of particles used [^CDC2013_Lesser]                     | Open-loop                                  |
|               |  `voronoi-open` |                          Probabilistically enforced upper  bound on overapproximation error  [^ACC2019_Voronoi]                          | Open-loop                                  |
|               | `chance-affine` |                                            Guaranteed underapproximation  [^CDC2019_chance]                                            | Affine   disturbance-feedback              |
|  `SReachSet`  |                 |  **Polytopic approximation of the stochastic  reach sets for the stochastic reachabilty  of a target tube problem**[^TAC2018_verification]<sup>,</sup>[^HSCC2018_cvxcmpt] | **Synthesize open-loop controllers in some cases** |
|               |  `chance-open`  |                                            Guaranteed underapproximation  [^TAC2018_verification]                                           | Optimal  open-loop controllers at vertices |
|               |  `genzps-open`  |                                 Approximation up to \\( \\epsilon\_\\mathrm{genz}\\), a user-specified quadrature error tolerance  [^TAC2018_verification]<sup>,</sup>[^HSCC2018_cvxcmpt]                                | Optimal  open-loop controllers at vertices |
|               |   `lag-under`   |                                            Guaranteed underapproximation [^CDC2017_Lagrangian]                                            |   Set-based feedback controller for all points within the set |
|               |    `lag-over`   |                                             Guaranteed overapproximation [^CDC2017_Lagrangian]                                            |                                            |
|  `SReachFwd`  |                 |      **Forward stochastic reachability analysis of an uncontrolled LTI/LTV system from a given initial state** [^HSCC2017_Fwd]<sup>,</sup>[^GenzAlgorithm]    |                                            |
|               |  `state-stoch`  |                                     Stochasticity of the state at a future time                                     |                                            |
|               |  `concat-stoch` |                     Stochasticity of the concatenated state vector up to a specified future time                    |                                            |
|               | `state-prob`    | Probability that the concatenated state vector (trajectory) up to a future time will lie in a given target tube set [^GenzAlgorithm] |                                            |
|               | `concat-prob`   | Probability that the state at a future time will lie in a given target set [^GenzAlgorithm]                                          |                                            |
|  `SReachDyn`  |                 |              **Dynamic programming approximation of the maximal reach probability and the reach  set**              |  **Analyze 2D and 3D LTI/LTV systems**     |

[Go to top](#posttop)

------
[^TAC2018_verification]: A. Vinod and M. Oishi, "[Stochastic reachability of a target tube:  Theory and computation](https://arxiv.org/pdf/1810.05217.pdf)", submitted to IEEE Transactions of Automatic Control, 2018 (submitted).
[^HSCC2018_cvxcmpt]: A. Vinod and M. Oishi, "[Scalable Underapproximative Verification of Stochastic LTI Systems using Convexity and Compactness](https://doi.org/10.1145/3178126.3178148)", in Proceedings of Hybrid Systems: Computation and Control, pp. 1--10, 2018.
[^CDC2019_chance]: A. Vinod and M. Oishi, "[Affine controller synthesis for stochastic reachability via difference of convex programming](https://hscl.unm.edu/affinecontrollersynthesis/)", in Proceedings of Conference on Decision and Control, 2019 (submitted).
[^CSSL2017_genzps]: A. Vinod and M. Oishi, "[Scalable Underapproximation for Stochastic Reach-Avoid Problem for High-Dimensional LTI Systems using Fourier Transforms](https://ieeexplore.ieee.org/document/7950904/)", in IEEE Control Systems Letters (CSS-L), pp. 316--321, 2017. 
[^CDC2013_Lesser]: K. Lesser, M. Oishi, and R. S. Erwin, "[Stochastic reachability for control of spacecraft relative motion](https://doi.org/10.1109/CDC.2013.6760626)," in Proceedings of the IEEE Conference on Decision and Control, pp. 4705-4712, 2013.
[^ACC2019_Voronoi]: H. Sartipizadeh, A. Vinod,  B. Acikmese, and M. Oishi, "[Voronoi Partition-based Scenario Reduction for Fast Sampling-based Stochastic Reachability Computation of LTI Systems](https://arxiv.org/abs/1811.03643)", In Proceedings of American Control Conference, 2019 (accepted).
[^CDC2017_Lagrangian]: J. Gleason, A. Vinod, and M. Oishi, "[Underapproximation of Reach-Avoid Sets for Discrete-Time Stochastic Systems via Lagrangian Methods](https://doi.org/10.1109/CDC.2017.8264291)," in Proceedings of the IEEE Conference on Decision and Control, pp. 4283-4290, 2017.
[^Automatica_Summers]: S. Summers and J. Lygeros, "[Verification of discrete time stochastic hybrid systems: A stochastic reach-avoid decision problem](https://doi.org/10.1016/j.automatica.2010.08.006)," Automatica, 2010.  
[^Automatica_Abate]: A. Abate, M. Prandini, J. Lygeros, and S. Sastry, "[Probabilistic reachability and safety for controlled discrete time stochastic hybrid systems](https://doi.org/10.1016/j.automatica.2008.03.027)," Automatica, 2008.
[^HSCC2017_Fwd]:  A. Vinod, B. HomChaudhuri, and M. Oishi, "[Forward Stochastic Reachability Analysis for Uncontrolled Linear Systems using Fourier Transforms](https://dl.acm.org/citation.cfm?doid=3049797.3049818)", in Proceedings of the 20th International Conference on Hybrid Systems: Computation and Control (HSCC), pp. 35-44, 2017. 
[^GenzAlgorithm]: A. Genz, "[Quadrature of a multivariate normal distribution over a region specified by linear inequalities: QSCMVNV](http://www.math.wsu.edu/faculty/genz/software/matlab/qscmvnv.m)", 2014. 
[^table_ack]: This table was generated using [https://www.tablesgenerator.com/markdown_tables#](https://www.tablesgenerator.com/markdown_tables#)
