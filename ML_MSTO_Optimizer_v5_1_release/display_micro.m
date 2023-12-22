% DISPLAY MICROSCALE OPTIMIZATION SETTINGS
function display_micro(micro)
    fprintf(' Microscale Optimization Settings:\n');
    fprintf('==========================================================\n');
    if micro.bc_type == 1
        fprintf('     Boundary Condition: Force-driven\n');
    elseif micro.bc_type == 2
        fprintf('     Boundary Condition: Displacement-driven\n');
    end
    if micro.dim == 2
        fprintf('     Mesh Resolution: nelx = %3d, nely = %3d\n',micro.nelx,micro.nely);
    elseif micro.dim == 3
        fprintf('     Mesh Resolution: nelx = %3d, nely = %3d, nelz = %3d\n',micro.nelx,micro.nely,micro.nelz);
    end
    if micro.phys == 1
        fprintf('     Physics mode: Mechanical\n');
    elseif micro.phys == 2
        fprintf('     Physics mode: Thermal\n');
    elseif all(micro.phys == [1 2])
        fprintf('     Physics mode: Thermomechanical\n');
    end
    if micro.cont
        fprintf('     Penalization: p = %1d, Continuation: On\n',micro.penal);
    else
        fprintf('     Penalization: p = %1d, Continuation: Off\n',micro.penal);
    end
    fprintf('     Density Filter Radius: rmin = %1.2f\n',micro.flt_den_min);
    fprintf('     Sensitivity Filter Radius: rmin = %1.2f\n',micro.flt_sen_min);
    fprintf('     Theta bounds: [%1.2f %1.2f]\n',micro.x_lb(1),micro.x_ub(1));
    fprintf('     Weight bounds: [%1.2f %1.2f]\n',micro.x_lb(2),micro.x_ub(2));
    fprintf('     Algorithm: Optimality Criterion\n');
    fprintf('     Loops: %3d, Tolerance: %1.1e\n\n',micro.maxloop,micro.tolx);
end