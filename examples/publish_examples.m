function publish_examples(varargin)
    %  Publish all examples in the folder to update their HTML and/or PDF
    %  versions
    % ======================================================================
    %
    %   publish_examples works like a unix command. It accepts options that
    %   modify its behaviour. By default, the command publishes all the
    %   examples in their HTML and PDF format.
    %
    %   Usage:
    %   ------
    %   publish_examples
    %   publish_examples --no-pdf
    %   publish_examples --no-html
    %   publish_examples --no-html
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


    publish_pdf = 1;
    publish_html = 1;

    exclusions = {'publish_examples.m'};

    if ~isempty(varargin)
        switch lower(varargin{1})
            case '--exclude'
                exclusions = [exclusions, varargin(2:end)];
            case '--no-pdf'
                publish_pdf = 0;
                if length(varargin) > 1
                    switch lower(varargin{2})
                        case '--exclude'
                            exclusions = [exclusions, varargin(3:end)];
                        case '--no-html'
                            disp('Nothing left to do');
                            return;
                        otherwise
                            error('Unhandled option');
                    end
                end
            case '--no-html'
                publish_html = 0;
                if length(varargin) > 1
                    switch lower(varargin{2})
                        case '--exclude'
                            exclusions = [exclusions, varargin(3:end)];
                        case '--no-pdf'
                            disp('Nothing left to do');
                            return;
                        otherwise
                            error('Unhandled option');
                    end
                end
            otherwise
                error('Unhandled option');
        end
    end

    % clc
    close all;

    set(0,'DefaultFigureWindowStyle','normal')

    % Get SReachTools path
    srt_rootpath = srtinit('--rootpath');
    examples_path = fullfile(srt_rootpath, 'examples');
    publish_path = fullfile(srt_rootpath, 'examples', 'publish');

    % Get the examples
    dl = dir(examples_path);

    for lv = 1:length(dl)
        if strcmp(dl(lv).name, '.') || strcmp(dl(lv).name, '..')
            continue;
        end

        [fpath, fname, ext] = fileparts(fullfile(examples_path, dl(lv).name));
        if ~strcmp(ext, '.m')
            fprintf('Skipping %s\n', dl(lv).name);
        else
            isExcepted = false;
            for lexc = 1:length(exclusions)
                if regexpi(dl(lv).name, ['^' exclusions{lexc}])
                    isExcepted = true;
                    break;
                end
            end

            if isExcepted
                fprintf('File ''%s'' matches exclusion filter. Skipping...\n', ...
                    dl(lv).name);
                continue;
            end

            save('PUBLISH_WORKSPACE.mat');

            % pdf
            fprintf('Publishing %s\n', dl(lv).name);
            if publish_pdf
                fprintf('    pdf... ');
                publish(dl(lv).name, 'format', 'pdf', 'outputDir',publish_path);
                fprintf('Done!\n');
            end

            % html
            if publish_html
                fprintf('    html... ');
                publish(dl(lv).name, 'format', 'html','outputDir',publish_path);
                fprintf('Done!\n')
            end

            load('PUBLISH_WORKSPACE.mat');
        end

        close all;
    end

    delete('PUBLISH_WORKSPACE.mat');
    % Copy the publish contents to doc for website
    current_wd = pwd;
    cd(srtinit('--rootpath'));
    copyfile ./examples/publish/ ./docs/examples/publish;
    cd(current_wd);
end


