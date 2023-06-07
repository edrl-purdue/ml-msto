% DISPLAY MACROSCALE OPTIMIZATION SETTINGS
function display_macro(macro)
    fprintf(' Macroscale Optimization Settings:\n');
    fprintf('==========================================================\n');
    fprintf(['     Problem name: ' macro.macro_bc_name '\n']);
    if macro.dim == 2
        fprintf('     Mesh Resolution: nelx = %3d, nely = %3d\n',macro.nelx,macro.nely);
    elseif macro.dim == 3
        fprintf('     Mesh Resolution: nelx = %3d, nely = %3d, nelz = %3d\n',macro.nelx,macro.nely,macro.nelz);
    end
    if macro.phys == 1
        fprintf('     Physics mode: Mechanical\n');
    elseif macro.phys == 2
        fprintf('     Physics mode: Thermal\n');
    elseif all(macro.phys == [1 2])
        fprintf('     Physics mode: Thermomechanical\n');
    end
    fprintf('     Volume Fraction: vf = %1.2f\n',macro.volfrac);
    fprintf('     Density Filter Radius: rmin = %1.2f\n',macro.flt_den_min);
    fprintf('     Sensitivity Filter Radius: rmin = %1.2f\n',macro.flt_sen_min);
    fprintf('     Density Grading limit: delta = %1.2f\n',macro.delta);
    fprintf('     Theta bounds: [%1.2f %1.2f]\n',macro.x_lb(1),macro.x_ub(1));
    fprintf('     Weight bounds: [%1.2f %1.2f]\n',macro.x_lb(2),macro.x_ub(2));
    fprintf('     Finite Difference Pertubation size: h = %1.4f\n',macro.h_fd(1));
    fprintf('     Volume Fraction Cutoff for Interpolation: vf = %1.2f\n',macro.vf_cutoff);
    if macro.alg == 1
        fprintf('     Algorithm: MATLAB''s fmincon\n');
    elseif macro.alg == 2
        fprintf('     Algorithm: Optimality Criterion\n');
    elseif macro.alg == 3
        fprintf('     Algorithm: Generalized Optimality Criterion\n');
    elseif macro.alg == 4
        fprintf('     Algorithm: Method of Moving Asymptotes\n');
    end
    fprintf('     Loops: %3d, Tolerance: %1.1e\n\n',macro.maxloop,macro.tolx);
    
    fprintf(' Macroscale De-homogenization Settings:\n');
    fprintf('==========================================================\n');
    if macro.dehom == 0
        fprintf('     No de-homogenization. \n');
    elseif macro.dehom == 1
        fprintf('     TOM de-homogenization. \n');
        fprintf('     TOM mesh resolution scale factor: %1.2f\n',macro.scaler);
        fprintf('     Unit cell size: epi = %1.2f\n',macro.epi);
        fprintf('     Density threshold for post-processing: %1.2f\n',macro.den_threshold);
        if macro.tri_infill
            fprintf('     Apply triangular infill connectivity post process: True\n');
        else
            fprintf('     Apply triangular infill connectivity post process: False\n');
        end
        if macro.rm_skel
            fprintf('     Apply skeletonization post process: True\n');
        else
            fprintf('     Apply skeletonization post process: False\n');
        end
    elseif macro.dehom == 2
        fprintf('     Projection de-homogenization with principal stress direction. \n');
    elseif macro.dehom == 3
        fprintf('     Projection de-homogenization via mesh. \n');
    end

    if macro.dehom ~= 0
        fprintf('     Number of remove low strain energy solid iterations: %1d\n',macro.rm_low_SE);
        if macro.rm_skel
            fprintf('     Calculate final compliance and compatibility: True\n');
        else
            fprintf('     Calculate final compliance and compatibility: False\n');
        end
    end
    fprintf('\n');
end