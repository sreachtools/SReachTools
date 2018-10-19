clear
clc
close all

srtinitHome = which('srtinit');
SReachToolsHome = strsplit(srtinitHome,'srtinit.m');
cd(SReachToolsHome{1});

set(0,'DefaultFigureWindowStyle','normal')
options = struct('format','pdf','outputDir','./examples/publish');
publish('AutomatedAnesthesiaDelivery.m', options);
options = struct('format','pdf','outputDir','./examples/publish');
publish('cwhSReachPointDemo.m', options);
options = struct('format','pdf','outputDir','./examples/publish');
publish('doubleIntegratorDynamicProgramming.m', options);
options = struct('format','pdf','outputDir','./examples/publish');
publish('forwardStochasticReachCWH.m', options);

% options = struct('format','html','outputDir','./examples/publish');
% publish('AutomatedAnesthesiaDelivery.m', options);
% options = struct('format','html','outputDir','./examples/publish');
% publish('cwhSReachPointDemo.m', options);
% options = struct('format','html','outputDir','./examples/publish');
% publish('doubleIntegratorDynamicProgramming.m', options);
% options = struct('format','html','outputDir','./examples/publish');
% publish('forwardStochasticReachCWH.m', options);

close all;