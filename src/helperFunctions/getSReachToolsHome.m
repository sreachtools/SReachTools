function SReachToolsHome = getSReachToolsHome()
    % Provides a string containing the home folder location of SReachTools
    % ==========================================================================
    srtinitHome = which('srtinit');
    SReachToolsHomeSplits = strsplit(srtinitHome,'srtinit.m');
    SReachToolsHome = SReachToolsHomeSplits{1};    
end
