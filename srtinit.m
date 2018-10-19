function varargout = srtinit(varargin)
%  Initialization function
% ======================================================================
%
%   Function to initialize and add source functions of the SReachTools toolbox
%   to the path.
%
%   Usage:
%   ------
%   srtinit
%   srtinit('--options');
%
% ======================================================================
% 
% srtinit options
% srtinit('options');
%
% Inputs:
% -------
%   Available options:
%       -v, --verbose    Have initalization function explicitly print to
%                        console which folders are being added to the path
%       -x, --deinit     Remove SReachTools toolbox folders from the path
%       -t, --test       Perform unit testing after initialization or deinit
%       -T               Perform unit testing without initialization or deinit,
%                        will cancel out any other parameters, e.g. '-x', '-v'
%
% Outputs:
% --------
%   None
%
% Notes:
% ------
% * Performing a deinit and testing '-x -t' will deinit the SReachTools toolbox 
%   and then perform unit testing, causing all unit tests to fail.
% 
% =========================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
% 
    
    SRTINIT_PATH = fileparts(which('srtinit.m'));

    verbose   = false;
    deinit    = false;
    run_tests = false;
    
    for lv = 1:length(varargin)
        if strcmp(varargin{lv}, '-v') || strcmp(varargin{lv}, '--verbose')
            verbose = true;
        elseif strcmp(varargin{lv}, '-x') || strcmp(varargin{lv}, '--deinit')
            deinit = true;
        elseif strcmp(varargin{lv}, '-t') || strcmp(varargin{lv}, '--test')
            run_tests = true;
        elseif strcmp(varargin{lv}, '-T')
            test_results = srttest();
            if nargout == 1
                varargout{1} = test_results;
            end
            return;
        elseif strcmp(varargin{lv}, '--update')
            update_sreachtools();
            return;
        elseif strcmp(varargin{lv}, '--version')
            dets = ver(SRTINIT_PATH);
            fprintf('SReachTools version %s\n', dets.Version);
            return;
        else
            throwAsCaller(SrtInvalidArgsError(['Invalid input option, ', ...
                'see help srtinit']));
        end
    end
    
    % get the parent dir of this function
    script_path = fileparts(mfilename('fullpath'));
    
    % absolute path to the SReachTools src directory
    src_path = fullfile(script_path, 'src');
    ex_path = fullfile(script_path, 'examples');
    
    if ismac                            % MAC (TODO: Untested)
        % new paths to add
        new_paths = strsplit([script_path ':' genpath(src_path), ...
            genpath(ex_path)], ':');
   
        % current MATLAB folders on path
        p = path;
        path_cell = split(p, ':');
    elseif ispc                         % WINDOWS
        % new paths to add
        new_paths = strsplit([script_path ';' genpath(src_path), ...
            genpath(ex_path)], ';');
   
        % current MATLAB folders on path
        p = path;
        path_cell = split(p, ';');
    elseif isunix                       % UNIX
        % new paths to add
        new_paths = strsplit([script_path ':' genpath(src_path), ...
            genpath(ex_path)], ':');
   
        % current MATLAB folders on path
        p = path;
        path_cell = split(p, ':');
    else
        disp('Platform not supported')
    end 
    if ~deinit
        % add paths
        for i = 1:length(new_paths) - 1
            if ~ismember(new_paths{i}, path_cell)
                if verbose
                    fprintf('Adding to path: %s\n', new_paths{i});
                end
                addpath(new_paths{i});
            end
        end
        
        %% Dependency checks                
        check_must_have_dependencies();
        prev_warn_state = getSrtWarning('SReachTools:setup');
        if verbose
            setSrtWarning('SReachTools:setup','on');
        else
            setSrtWarning('SReachTools:setup','off');
        end
        check_recommended_dependencies();
        setSrtWarning('SReachTools:setup', prev_warn_state);        
    else
        % remove paths
        for i = 1:length(new_paths) - 1
            if ismember(new_paths{i}, path_cell)
                if verbose
                    fprintf('Removing from path: %s\n', new_paths{i});
                end
                rmpath(new_paths{i});
            end
        end
    end
    
    if run_tests
        srttest();
    end
end

function test_results = srttest()
    % get the parent dir of this function
    script_path = fileparts(mfilename('fullpath'));
    
    current_dir = pwd;
    
    cd(fullfile(script_path, 'tests'));
    test_results = runtests();
    disp(test_results)
    cd(current_dir);
end

function update_sreachtools()

    SRTINIT_PATH = fileparts(which('srtinit.m'));

    % First check current version and newest version from repository
    % get repository tags
    fprintf('Fetching repository data...\n')
    tag_data = webread(['https://api.github.com/repos/unm-hscl/', ...
        'SReachTools/tags']);

    local_version = ver(SRTINIT_PATH);
    fprintf('Current SReachTools version  :  %s\n', local_version.Version);

    fprintf('Newest SReachTools version   :  %s\n', tag_data(1).name(2:end));

    if strcmp(tag_data(1).name(2:end), local_version.Version)
        fprintf('No newer version available.\n')
    else
        % download and insall newest version
        fprintf('Updating...\n')
    end

end

function check_must_have_dependencies()
    v = ver;    
    %% Check for Gaussian pdf/cdf etc in Statistics and Machine Learning Toolbox
    has_normcdf = any(strcmp(cellstr(char(v.Name)),...
        'Statistics and Machine Learning Toolbox'));
    if ~has_normcdf
        exc = SrtSetupError(['SReachTools needs MATLAB''s ', ...
            'Statistics and Machine Learning Toolbox.']);
        throw(exc);
    end
    %% Check for MPT3
    has_mpt3 = any(strcmp(cellstr(char(v.Name)),'Multi-Parametric Toolbox'));
    if ~has_mpt3
        exc = SrtSetupError(['SReachTools needs Multi-Parameteric Toolbox ',...
            '(MPT3). See https://www.mpt3.org/Main/Installation for ',...
            'installation instructions.']);
        throw(exc);
    end
    %% Check for CVX
    has_cvx = (exist('cvx_begin', 'file') == 2) &&...
        (exist('cvx_end', 'file') == 2) && (exist('cvx_setup', 'file') == 2); 
    if ~has_cvx
        exc = SrtSetupError(['SReachTools needs CVX. See ',...
            'http://cvxr.com/cvx/download/ for installation instructions.']);
        throw(exc);
    end
end

function check_recommended_dependencies()
    v = ver;    
    %% Check for Symbolic Math Toolbox for PWA
    has_syms = any(strcmp(cellstr(char(v.Name)), 'Symbolic Math Toolbox'));
    if ~has_syms
        warning('SReachTools:setup',['Piecewise approximation ',...
            'construction getPWAOverAndUnderApprox() in SReachTools ',...
            'requires MATLAB''s Symbolic Math Toolbox.']);
    end    
    %% Check for Global Optimization Toolbox (patternsearch)
    has_patternsearch = any(strcmp(cellstr(char(v.Name)),... 
        'Global Optimization Toolbox'));
    if ~has_patternsearch
        warning('SReachTools:setup',['''genzps-open'' option of ',...
            'SReachPoint() function in SReachTools requires MATLAB''s ',...
            'Global Optimization Toolbox.']);
    end    
    %% Check for Optimization Toolbox (fmincon)
    has_fmincon = any(strcmp(cellstr(char(v.Name)), 'Optimization Toolbox'));
    if ~has_fmincon
        warning('SReachTools:setup',['''lag-over/under'' option of ',...
            'SReachSet() function in SReachTools requires MATLAB''s ',...
            'Optimization Toolbox in one of its techniques: ''optim-box''.']);
    end
    %% Check for Gurobi
    [default_solver, solvers_cvx] = cvx_solver;
    options_mpt = mptopt;
    if ~(contains(default_solver,'Gurobi') &&...
            any(contains(options_mpt.solvers_list.MIQP,'GUROBI')))    
        warning('SReachTools:setup',['Gurobi is the recommended backend ',...
            'solver for MPT3 and CVX when using SReachTools.']);
    end
end