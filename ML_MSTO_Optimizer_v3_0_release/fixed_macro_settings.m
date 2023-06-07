% LIMITS CERTAIN MACRO OPTIMIZATION SETTINGS SO ERRORS ARE NOT THROWN
function macro = fixed_macro_settings(macro)
    
    if macro.cont ~= 0
        macro.cont = 0;
    end
    
    if macro.flt_den_min ~= 1.0
        macro.flt_den_min = 1.0;
    end

    if macro.x_lb(2) ~= 0
        macro.x_lb(2) = 0;
    end
    if macro.x_ub(2) ~= 1
        macro.x_ub(2) = 1;
    end
    
    if macro.h_fd ~= 0
        macro.h_fd = 0;
    end

    if macro.vf_cutoff ~= 0.1
        macro.vf_cutoff = 0.1;
    end

    if macro.alg ~= 3
        macro.alg = 3;
    end

    if macro.dehom ~= 1
        macro.dehom = 1;
    end

    if macro.scaler ~= 1
        macro.scaler = 1;
    end

    if macro.epi ~= 1
        macro.epi = 1;
    end

    if macro.den_threshold ~= 0.5
        macro.den_threshold = 0.5;
    end

    if macro.tri_infill ~= 0
        macro.tri_infill = 0;
    end

    if macro.rm_skel ~= 0
        macro.rm_skel = 0;
    end

    if macro.rm_low_SE ~= 0
        macro.rm_low_SE = 0;
    end

    if macro.final_comp ~= 0
        macro.final_comp = 0;
    end
    
    if macro.phys ~= 1
        macro.phys = 1;
    end

    if macro.obj ~= 1
        macro.obj = 1;
    end

    if macro.bc_type ~= 1
        macro.bc_type = 1;
    end

    if macro.design_ini ~= 1
        macro.design_ini = 1;
    end
end