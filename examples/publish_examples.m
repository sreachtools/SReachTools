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
%   All types of function calls
%
%   1. Publish all examples
%       publish_examples
%   2a. Publish all examples with no pdf only html
%       publish_examples --no-pdf
%   2b. Publish certain examples with no pdf only html
%       publish_examples --no-pdf --publish-only SPACE_SEPARATED_FILENAMES
%   2b. Publish after skipping certain examples with no pdf only html
%       publish_examples --no-pdf --exclude SPACE_SEPARATED_FILENAMES
%   3a. Publish all examples with no html only pdf
%       publish_examples --no-html
%   3b. Publish certain examples with no html only pdf
%       publish_examples --no-html --publish-only SPACE_SEPARATED_FILENAMES
%   3b. Publish after skipping certain examples with no html only pdf
%       publish_examples --no-html --exclude SPACE_SEPARATED_FILENAMES
%   4. Publish only certain examples
%       publish_examples --publish-only SPACE_SEPARATED_FILENAMES
%   5. Publish after skipping certain examples
%       publish_examples --exclude SPACE_SEPARATED_FILENAMES
%
% =========================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://sreachtools.github.io/license/
% 
% 


    publish_pdf = 1;
    publish_html = 1;

    % Get SReachTools path
    srt_rootpath = srtinit('--rootpath');
    examples_path = fullfile(srt_rootpath, 'examples');
    publish_path = fullfile(srt_rootpath, 'examples', 'publish');

    % Get all the MATLAB script files
    dl = dir(strcat(examples_path,'/*.m'));
    
    % Collect the file names into
    cell_dl_filenames = {dl(:).name};
    
    % Makes sure to skip this script
    index_of_this_script = find(contains(cell_dl_filenames, ...
        'publish_examples.m')==1);
    
    if ~isempty(varargin)
        switch lower(varargin{1})
            case '--publish-only'
                % Remove all the scripts not provided in varargin
                indices_to_remove = find(contains(cell_dl_filenames, ...
                    {varargin{2:end}})==0);
                dl([index_of_this_script,indices_to_remove])=[];
            case '--exclude'
                % Remove all the scripts provided in varargin
                indices_to_remove = find(contains(cell_dl_filenames, ...
                    {varargin{2:end}})==1);
                dl([index_of_this_script,indices_to_remove])=[];
            case '--no-pdf'
                publish_pdf = 0;
                if length(varargin) == 1
                    dl(index_of_this_script)=[];                            
                else
                    switch lower(varargin{2})
                        case '--publish-only'
                            % Remove all the scripts not provided
                            indices_to_remove = find( ...
                                contains(cell_dl_filenames, ...
                                    {varargin{3:end}})==0);
                            dl([index_of_this_script, ...
                                indices_to_remove])=[];
                            if length(indices_to_remove) == ...
                                    length(cell_dl_filenames)
                                error('Nothing to publish? Check filenames.');
                            end
                        case '--exclude'
                            % Remove all the scripts provided in varargin
                            indices_to_remove = find( ...
                                contains(cell_dl_filenames, ...
                                    {varargin{3:end}})==1);
                            dl([index_of_this_script, ...
                                indices_to_remove])=[];
                            if isempty(indices_to_remove)
                                error('Nothing to exclude? Check filenames.');
                            end
                        case '--no-pdf'
                            disp('Nothing left to do');
                            return;
                        otherwise
                            error('Unhandled option');
                    end
                end
            case '--no-html'
                publish_html = 0;
                if length(varargin) == 1
                    dl(index_of_this_script)=[];                            
                else
                    switch lower(varargin{2})
                        case '--publish-only'
                            % Remove all the scripts not provided
                            indices_to_remove = find( ...
                                contains(cell_dl_filenames, ...
                                    {varargin{3:end}})==0);
                            dl([index_of_this_script, ...
                                indices_to_remove])=[];
                            if length(indices_to_remove) == ...
                                    length(cell_dl_filenames)
                                error('Nothing to publish? Check filenames.');
                            end
                        case '--exclude'
                            % Remove all the scripts provided in varargin
                            indices_to_remove = find( ...
                                contains(cell_dl_filenames, ...
                                    {varargin{3:end}})==1);
                            dl([index_of_this_script, ...
                                indices_to_remove])=[];
                            if isempty(indices_to_remove)
                                error('Nothing to exclude? Check filenames.');
                            end
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

    for lv = 1:length(dl)
        save('PUBLISH_WORKSPACE.mat');

        % pdf
        fprintf('Publishing %s\n', dl(lv).name);
        if publish_pdf
            fprintf('    pdf... ');
            publish(dl(lv).name, 'format', 'pdf','outputDir',publish_path);
            fprintf('Done!\n');
        end

        % html
        if publish_html
            fprintf('    html... ');
            publish(dl(lv).name, 'format','html','outputDir',publish_path);
            fprintf('Done!\n')
        end

        load('PUBLISH_WORKSPACE.mat');
        
        close all;
    end

    delete('PUBLISH_WORKSPACE.mat');
    fprintf(['\nRemember to copy over the ./examples/publish/ over to ', ...
        'the webpage folder, for publishing.\n']);
end


