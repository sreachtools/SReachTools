---
layout: docs
title: Documentation
permalink: /docs/
---

This page lists all the functions included as part of SReachTools. This information is also available in the docstrings. Please use:

- `help FUNCTION_NAME` to understand the details of a function
- `methods(CLASS_NAME)` to understand the details of a class
- `<Ctrl+F1>` to get function hints for a given function



Do check out the [examples](../examples) for more detailed explanations of various functionalities of the toolbox.

## Function list

<ul class="doc-list">
    <li>src/</li>
    <ul class="doc-list">
        <li>classes/</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="src/classes/@LtiSystem/LtiSystem">@LtiSystem/</a></li>
            <li class="doc-list"><a href="src/classes/@LtvSystem/LtvSystem">@LtvSystem/</a></li>
            <li class="doc-list"><a href="src/classes/RandomVector">RandomVector.m</a></li>
            <li class="doc-list"><a href="src/classes/SReachEllipsoid">SReachEllipsoid.m</a></li>
            <li class="doc-list"><a href="src/classes/Tube">Tube.m</a></li>
        </ul>
        <li>exceptions/</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="src/exceptions/SrtBaseException">SrtBaseException.m</a></li>
            <li class="doc-list"><a href="src/exceptions/SrtDevError">SrtDevError.m</a></li>
            <li class="doc-list"><a href="src/exceptions/SrtInternalError">SrtInternalError.m</a></li>
            <li class="doc-list"><a href="src/exceptions/SrtInvalidArgsError">SrtInvalidArgsError.m</a></li>
            <li class="doc-list"><a href="src/exceptions/SrtSetupError">SrtSetupError.m</a></li>
            <li class="doc-list"><a href="src/exceptions/SrtTestError">SrtTestError.m</a></li>
        </ul>
        <li>helperFunctions/</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="src/helperFunctions/allcomb">allcomb.m</a></li>
            <li class="doc-list"><a href="src/helperFunctions/ellipsoidsFromMonteCarloSims">ellipsoidsFromMonteCarloSims.m</a></li>
            <li class="doc-list"><a href="src/helperFunctions/generateMonteCarloSims">generateMonteCarloSims.m</a></li>
            <li class="doc-list"><a href="src/helperFunctions/getSrtWarning">getSrtWarning.m</a></li>
            <li class="doc-list"><a href="src/helperFunctions/iteratedQscmvnv">iteratedQscmvnv.m</a></li>
            <li class="doc-list"><a href="src/helperFunctions/minkSumInner">minkSumInner.m</a></li>
            <li class="doc-list"><a href="src/helperFunctions/qscmvnv">qscmvnv.m</a></li>
            <li class="doc-list"><a href="src/helperFunctions/setSrtWarning">setSrtWarning.m</a></li>
        </ul>
        <li>modules/</li>
        <ul class="doc-list">
            <li>fwdStochReach/</li>
            <ul class="doc-list">
                <li class="doc-list"><a href="src/modules/fwdStochReach/SReachFwd">SReachFwd.m</a></li>
                <li class="doc-list"><a href="src/modules/fwdStochReach/computeReachProb">computeReachProb.m</a></li>
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
