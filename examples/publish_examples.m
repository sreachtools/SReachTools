function publish_examples(varargin)

    if strcmp(varargin{1}, '--exclude')
        exclusions = varargin(2:end);
    else
        error('Unhandled option');
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
            fprintf('    pdf... ');
            publish(dl(lv).name, 'format', 'pdf', 'outputDir', publish_path);
            fprintf('Done!\n');

            % html
            fprintf('    html... ');
            publish(dl(lv).name, 'format', 'html', 'outputDir', publish_path);
            fprintf('Done!\n')

            load('PUBLISH_WORKSPACE.mat');
        end

        close all;
    end

    delete('PUBLISH_WORKSPACE.mat');
end


