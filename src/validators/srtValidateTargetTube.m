function srtValidateTargetTube(tt)
    validateattributes(tt, {'cell'}, {'nonempty'});

    for lv = 1:length(target_tube)
        validateattributes(target_tube{lv}, {'Polyhedron'}, {'nonempty'});
    end
end