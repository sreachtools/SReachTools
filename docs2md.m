%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of the programmer: Joseph Gleason %
% Date: 2018-10-19                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Purpose
% Generate Markdown files from the docstrings

function docs2md()

    global tabstr FUNCTION_BASE_PATH DOCS_ROOT_PATH
    tabstr = '';
    FUNCTION_BASE_PATH = fileparts(mfilename('fullpath'));
    DOCS_ROOT_PATH = fullfile(FUNCTION_BASE_PATH, 'docs/_docs/');

    rmdir(DOCS_ROOT_PATH, 's');
    checkAndMakeFolder(DOCS_ROOT_PATH);

    % write the base docs index file
    fid = fopen(fullfile(DOCS_ROOT_PATH, 'index.md'), 'w+');
    fprintf(fid, '---\n');
    fprintf(fid, 'layout: docs\n');
    fprintf(fid, 'title: Function List\n');
    fprintf(fid, 'permalink: /docs/\n');
    fprintf(fid, '---\n\n');
    fclose(fid);

    checkAndMakeFolder(fullfile(DOCS_ROOT_PATH, 'src'));

    fid = fopen(fullfile(DOCS_ROOT_PATH, 'index.md'), 'a');
    fprintf(fid, sprintf('%s<ul class="doc-list">\n', tabstr));
    fclose(fid);

    tabstr = [tabstr, '    '];
    listAllFilesRecursive('src', FUNCTION_BASE_PATH);
    tabstr = '';

    fid = fopen(fullfile(DOCS_ROOT_PATH, 'index.md'), 'a');
    fprintf(fid, sprintf('%s</ul>\n', tabstr));
    fclose(fid);
end

function checkAndMakeFolder(fpath)
    if exist(fpath, 'dir') ~= 7
        % fprintf('Creating %s\n', fpath)
        try
            mkdir(fpath)
        catch err
            throw(SrtInternalError('Could not create directory %s', ...
                fpath));
        end
    end
end

function listAllFilesRecursive(folder_name, base_dir)

    global tabstr DOCS_ROOT_PATH

    len_folder_name = length(folder_name);
    old_tabstr = tabstr;

    IGNORE_DIRS = {'.', '..', '.git'};

    % check if the folder is a class folder; currently not going to print
    % the contents of a class folder because this will get its own more
    % dedicated page

    if regexp(folder_name, '^@')
        const_name = fullfile(base_dir, folder_name, [folder_name(2:end) '.m']);
        if exist(const_name, 'file') == 2
            % found constructor, write class
            [cpath, cname, cext] = fileparts(const_name);
            
            fid = fopen(fullfile(DOCS_ROOT_PATH, 'index.md'), 'a');
            md_file_path = functionPath2DocsPath(const_name);
            fprintf(fid, sprintf(['%s<li class="doc-list"><a href="%s">', ...
                '%s/</a></li>\n'], ...
                tabstr, getLocalPath(const_name, 'toolbox'), folder_name));
            fclose(fid);

            wrtieClassPageDoc(const_name);
        else
            warning('Class file %s not found. Skipping...', const_name);
        end

        % fid = fopen(fullfile(DOCS_ROOT_PATH, 'index.md'), 'a');
        % fprintf(fid, sprintf(['%s<li class="doc-list"><a href="#%s">%s/</a></li>\n'], tabstr, ...
        %     folder_name(2:end), folder_name));

        % fclose(fid);

        % checkAndMakeFolder(functionPath2DocsPath(fullfile(base_dir, ...
        %     folder_name(2:end))));

        % if exist(fullfile(base_dir, folder_name, ...
        %     [folder_name(2:end) '.m']), 'file') == 2

        %     wrtieClassPageDoc(fullfile(base_dir, folder_name, ...
        %         [folder_name(2:end) '.m']));
        % else
        %     warning('Class file %s not found. Skipping...', ...
        %         fullfile(base_dir, folder_name, [folder_name(2:end) '.m']));
        % end

        return;
    end

    dl = dir(fullfile(base_dir, folder_name));

    fid = fopen(fullfile(DOCS_ROOT_PATH, 'index.md'), 'a');
    fprintf(fid, sprintf('%s<li>%s/</li>\n', tabstr, folder_name));
    fprintf(fid, sprintf('%s<ul class="doc-list">\n', tabstr));
    fclose(fid);
    tabstr = [tabstr, '    '];

    for lv = 1:length(dl)
        if any(strcmp(dl(lv).name, IGNORE_DIRS))
            continue;
        end

        if dl(lv).isdir()
            listAllFilesRecursive(dl(lv).name, dl(lv).folder);
        else
            fname = fullfile(dl(lv).folder, dl(lv).name);
            [~, name, ext] = fileparts(fname);
            if any(strcmp(ext, {'.m', '.matlab'}))
                fid = fopen(fullfile(DOCS_ROOT_PATH, 'index.md'), 'a');
                fprintf(fid, sprintf(['%s<li class="doc-list"><a href="%s">%s</a></li>\n'], ...
                    tabstr, getLocalPath(fname, 'toolbox'), [name, ext]));
                fclose(fid);

                writeDocStringToFile(fname);
            end
        end
    end

    tabstr = old_tabstr;
    fid = fopen(fullfile(DOCS_ROOT_PATH, 'index.md'), 'a');
    fprintf(fid, sprintf('%s</ul>\n', tabstr));
    fclose(fid);
%     movefile(strcat(DOCS_ROOT_PATH','index.md'),strcat(DOCS_ROOT_PATH','../docs.md'));
end

function bool = isclassfile(f)

    bool = false;
    return;

    fid = fopen(f);

    while ~eof(fid)
        ln = fgetl(fid);
        if regexp(ln, '^\s*%')
            continue;
        else
            if regexp(ln, '^\s*classdef')
                bool = true;
            else
                bool = false;
            end
            break;
        end
    end
end

function localpath = getLocalPath(file_name, opt)

    global FUNCTION_BASE_PATH DOCS_ROOT_PATH

    switch(opt)
        case 'toolbox'
            file_name = functionPath2DocsPath(file_name);
            localpath = getLocalPath(file_name, 'doc');
        case 'doc'
            ln = length(DOCS_ROOT_PATH);

            % eliminate the preamble
            localpath = file_name(ln+1:end);

            % get rid of the .md
            localpath = localpath(1:end-3);

            % add index.md
            % localpath = fullfile(localpath, 'index.md');
        otherwise
            error('Unhandled option')
    end
end

function writeDocStringToFile(file_name)

    global FUNCTION_BASE_PATH DOCS_ROOT_PATH

    [fpath, fname, fext] = fileparts(file_name);

    metaobj = meta.class.fromName(fname);

    if ~isempty(metaobj);
        wrtieClassPageDoc(file_name);
    else
        writeFunctionPageDoc(file_name);
        % hstr = help(file_name);

        % [~, fun_name, ext] = fileparts(file_name);

        % file_name = functionPath2DocsPath(file_name);

        % % remove the file extension '.md'
        % file_name = file_name(1:end-3);

        % % create the base folder so that can write to /some/path/index.md
        % if exist(file_name, 'dir') ~= 7
        %     fprintf('Writing %s ...\n', fullfile(file_name, 'index.md'));
        %     mkdir(file_name);
        % end

        % % write index.md
        % fid = fopen(fullfile(file_name, 'index.md'), 'w+');
        % fprintf(fid, '---\n');
        % fprintf(fid, 'layout: docs\n');

        % % function name with extension for title
        % fprintf(fid, 'title: %s\n', [fun_name, ext]);
        % fprintf(fid, '---\n\n');

        % fprintf(fid, '```\n%s```\n', hstr);
        % fclose(fid);
    end
end

function docsPath = functionPath2DocsPath(file_name)

    global FUNCTION_BASE_PATH DOCS_ROOT_PATH

    % need to remove the function path and replace it with the docs path
    si = length(FUNCTION_BASE_PATH) + 2;

    [p, ~, file_ext] = fileparts(file_name);

    switch(file_ext)
        case '.m'
            % remove the .m and replace with .md
            file_name = fullfile(DOCS_ROOT_PATH, file_name(si:end-2));
            file_name = [file_name, '.md'];
        case '.matlab'
            % remove the .matlab and replace with .md
            file_name = fullfile(DOCS_ROOT_PATH, file_name(si:end-7));
            file_name = [file_name, '.md'];
        case ''
            file_name = fullfile(DOCS_ROOT_PATH, file_name(si:end));
        otherwise
            error('Unhandled file option')
    end

    docsPath = file_name;
end

function writeFunctionPageDoc(file_name)
    hstr = help(file_name);

    [fun_path, fun_name, ext] = fileparts(file_name);

    [md_path, md_name, ~] = fileparts(functionPath2DocsPath(file_name));

    % create the base folder so that can write to /some/path/index.md
    checkAndMakeFolder(fullfile(md_path, md_name))

    fprintf('Writing %s ...\n', fullfile(md_path, md_name, 'index.md'));
    % write index.md
    fid = fopen(fullfile(md_path, md_name, 'index.md'), 'w+');
    fprintf(fid, '---\n');
    fprintf(fid, 'layout: docs\n');

    % function name with extension for title
    fprintf(fid, 'title: %s\n', [fun_name, ext]);
    fprintf(fid, '---\n\n');

    fprintf(fid, '```\n%s```\n', hstr);
    fclose(fid);
end

function wrtieClassPageDoc(ffpath)

    [class_path, class_name, ext] = fileparts(ffpath);

    atinds = regexp(class_path, '@');
    offset = 0;
    for lv = 1:length(atinds)
        class_path = class_path(1:atinds(lv)-offset-1);
        offset = offset + 1;
    end

    tabstr = '';

    md_name = functionPath2DocsPath(ffpath);
    md_path = md_name(1:end-3);

    fprintf('Writing %s ...\n', fullfile(md_path, 'index.md'));
    checkAndMakeFolder(md_path);

    fid = fopen(fullfile(md_path, 'index.md'), 'w+');
    fprintf(fid, '---\n');
    fprintf(fid, 'layout: docs\n');

    metaobj = meta.class.fromName(class_name);
    disp_property_count = 0;
    disp_method_count = 0;

    % for lv = 1:length(metaobj.PropertyList)
    %     if metaobj.PropertyList(lv).Hidden || ...
    %         strcmp(metaobj.PropertyList(lv).GetAccess, 'private')

    %     else
    %         disp_property_count

    % function name with extension for title
    fprintf(fid, 'title: %s\n', [class_name, ext]);
    fprintf(fid, '---\n\n');

    fprintf(fid, '%s<ul %s>\n', tabstr, getUlClassStr());
    tabstr = incrementTab(tabstr);
    % +1 tab inc
    fprintf(fid, '%s<li %s><a href="%s">%s</a></li>\n', tabstr, ...
        getLiClassStr(), ['#', metaobj.Name], metaobj.Name);
    fprintf(fid, '%s<ul %s>\n', tabstr, getUlClassStr());
    tabstr = incrementTab(tabstr);
    % +2 tab inc

    fprintf(fid, '%s<li><a href="%s">Constructor</a></li>\n', tabstr, ...
        ['#' metaobj.Name '-' metaobj.Name]);

    fprintf(fid, '%s<li>Properties</li>\n', tabstr);
    fprintf(fid, '%s<ul %s>\n', tabstr, getUlClassStr());
    tabstr = incrementTab(tabstr);
    % +3 tab inc
    for lv = 1:length(metaobj.PropertyList)
        if isDispProperty(metaobj, lv)

            fprintf(fid, '%s<li %s><a href="%s">%s</a></li>\n', ...
                tabstr, ...
                getLiClassStr(), ...
                ['#' metaobj.Name '-prop-' metaobj.PropertyList(lv).Name], ...
                metaobj.PropertyList(lv).Name);
        end
    end
    tabstr = decrementTab(tabstr);
    % +2 tab inc
    fprintf(fid, '%s</ul>\n', tabstr);

    fprintf(fid, '%s<li>Methods</li>\n', tabstr);
    fprintf(fid, '%s<ul %s>\n', tabstr, getUlClassStr());
    tabstr = incrementTab(tabstr);
    % +3 tab inc
    for lv = 1:length(metaobj.MethodList)
        if isDispMethod(metaobj, lv)

            fprintf(fid, '%s<li %s><a href="%s">%s</a></li>\n', ...
                tabstr, ...
                getLiClassStr(), ...
                ['#' metaobj.Name '-method-' metaobj.MethodList(lv).Name], ...
                metaobj.MethodList(lv).Name);
        end
    end
    tabstr = decrementTab(tabstr);
    % +2 tab inc
    fprintf(fid, '%s</ul>\n', tabstr);


    tabstr = decrementTab(tabstr);
    % +1 tab inc
    fprintf(fid, '%s</ul>\n', tabstr);

    tabstr = decrementTab(tabstr);
    % +0 tab inc
    fprintf(fid, '%s</ul>\n', tabstr);


    fprintf(fid, '\n');

    % print the main help
    fprintf(fid, '{:#%s}\n', class_name);
    fprintf(fid, '### %s\n', class_name);
    fprintf(fid, '```\n%s```\n\n', help(class_name));

    % print constructor help
    fprintf(fid, '{:#%1$s-%1$s}\n', class_name);
    fprintf(fid, '### Constructor\n');
    fprintf(fid, '```\n%s```\n\n', help(class_name));

    % print the help for properties
    for lv = 1:length(metaobj.PropertyList)
        if isDispProperty(metaobj, lv)

            fprintf(fid, '### %s\n', ...
                ['Property: ' metaobj.PropertyList(lv).Name]);
            fprintf(fid, '{:#%s-prop-%s}\n', ...
                class_name, metaobj.PropertyList(lv).Name);
            hstr = help([class_name, '.', metaobj.PropertyList(lv).Name]);
            if isempty(hstr)
                fprintf(fid, ['{%% include important-note.html ', ...
                    'content="Property currently has no help ', ...
                    'documentation." %%}\n']);
            else
                fprintf(fid, '```\n%s```\n', help([class_name, '.', ...
                    metaobj.PropertyList(lv).Name]));
            end
            fprintf(fid, '\n');
        end
    end

    % print the help for methods
    for lv = 1:length(metaobj.MethodList)
        if isDispMethod(metaobj, lv)

            fprintf(fid, '### %s\n', ...
                ['Method: ' metaobj.MethodList(lv).Name]);
            fprintf(fid, '{:#%s-method-%s}\n', ...
                class_name, metaobj.MethodList(lv).Name);
            hstr = help([class_name, '.', metaobj.MethodList(lv).Name]);
            if isempty(hstr)
                fprintf(fid, ['{%% include important-note.html ', ...
                    'content="Method currently has no help ', ...
                    'documentation." %%}\n']);
            else
                fprintf(fid, '```\n%s```\n', help([class_name, '.', ...
                    metaobj.MethodList(lv).Name]));
            end
            fprintf(fid, '\n');
        end
    end

    fclose(fid);
end

function tabstr = incrementTab(tabstr)
    tabstr = [tabstr, '    '];
end

function tabstr = decrementTab(tabstr)
    tabstr = tabstr(1:end-4);
end

function str = getLiClassStr()
    str = 'class="doc-list"';
end

function str = getUlClassStr()
    str = 'class="doc-list"';
end

function bool = isDispMethod(mo, ind)
    ol_methods = {'subsref', 'subsassn', 'length', 'size', 'disp', ...
        'end'};

    bool = ~mo.MethodList(ind).Hidden && ...
            strcmp(mo.MethodList(ind).Access, 'public') && ...
            ~any(strcmp(mo.MethodList(ind).Name, ol_methods)) && ...
            ~strcmp(mo.MethodList(ind).Name, mo.Name) && ...
            mo.MethodList(ind).DefiningClass == mo;
end

function bool = isDispProperty(mo, ind)
    bool = ~mo.PropertyList(ind).Hidden && ...
            strcmp(mo.PropertyList(ind).GetAccess, 'public') && ...
            mo.PropertyList(ind).DefiningClass == mo;
end
