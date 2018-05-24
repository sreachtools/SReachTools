function new_nvcell = removeNameValuePairs(orig_nvcell, names)

    inpar = inputParser();
    inpar.addRequired('orig_nvcell', @(x) validateattributes(x, ...
        {'cell'}, {'row', 'nonempty'}));
    inpar.addRequired('names', @(x) validateattributes(x, ...
        {'cell', 'char'}, {'row', 'nonempty'}));

    inpar.parse(orig_nvcell, names)

    for lv = 1:2:length(orig_nvcell)-1
        if ~ischar(orig_nvcell{lv})
            error('Original name-value cell array must be name-value pairs')
        end
    end

    if ischar(names)
        names = {names};
    end

    new_nvcell = cell(1, length(orig_nvcell) - 2*length(names));
    nvk = 1;
    for lv = 1:2:length(orig_nvcell)-1
        if ~any(strcmp(orig_nvcell{lv}, names))
            new_nvcell{nvk}   = orig_nvcell{lv};
            new_nvcell{nvk+1} = orig_nvcell{lv+1};

            nvk = nvk + 2;
        end
    end

end