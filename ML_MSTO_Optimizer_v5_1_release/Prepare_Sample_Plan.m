% GENERATE FEATURES AND SAMPLE PLAN
function [X,scale] = Prepare_Sample_Plan(scale)
    %% GENERATE FEATURES AND SAMPLE PLAN FOR HOMOGENIZATION
    % Dimension Counting
    scale.spdim = 0; %[Dimension]
    scale.spdim = scale.spdim + 1; % volume fraction input
    scale.spdim = scale.spdim + 1; % objective function weight input
    % scale.spdim = scale.spdim + 1; % penalization parameter input
    scale.spdim_nums = zeros(2,3); % rows = start dim num, end dim num for Fx/Ux_e. cols = number of phyics models
    if any(scale.phys == 1)
        scale.spdim_nums(1,1) = scale.spdim;
        scale.spdim = scale.spdim + scale.nen*scale.dim;
        scale.spdim_nums(2,1) = scale.spdim;
    end
    if any(scale.phys == 2)
        scale.spdim_nums(1,2) = scale.spdim;
        scale.spdim = scale.spdim + scale.nen;
        scale.spdim_nums(2,2) = scale.spdim;
    end

    % Feature Bounds
    lb = ones(1,scale.spdim)*-1; ub = ones(1,scale.spdim); % boundary condition input range
    lb(1) = scale.x_lb(1); ub(1) = scale.x_ub(1); % volume fraction input range
    if length(scale.phys) == 1; lb(2) = 1; ub(2) = 1; else; lb(2) = scale.x_lb(2); ub(2) = scale.x_ub(2); end % objective function weight input range
    % if scale.cont == 1; lb(3) = 1; ub(3) = scale.penal; else; lb(3) = scale.penal; ub(3) = scale.penal; end % penalization parameter input range

    % Latin Hypercube Sampling
    X_unnorm = lhsdesign(ceil(scale.N*scale.rr),scale.spdim-2,'iterations',10); %Volume fraction and Nodal force/disp Input
    X_unnorm = [X_unnorm(:,1) lhsdesign(ceil(scale.N*scale.rr),2,'iterations',10) X_unnorm(:,2:end)]; %Objective function Weight and penalization parameter Input
    for i = 1:ceil(scale.N*scale.rr)
        X_unnorm(i,:) = X_unnorm(i,:).*(ub-lb) + lb;
    end

    % Normalize
    X = X_unnorm;
    for i = 1:ceil(scale.N*scale.rr)
        if any(scale.phys == 1)
            X(i,(scale.spdim_nums(1,1)+1):scale.spdim_nums(2,1)) = Normalize_TOM_BC(X_unnorm(i,(scale.spdim_nums(1,1)+1):scale.spdim_nums(2,1)),scale); % normalize subproblem boundary condition values
        end
        if any(scale.phys == 2)
            X(i,(scale.spdim_nums(1,2)+1):scale.spdim_nums(2,2)) = Normalize_TOM_BC(X_unnorm(i,(scale.spdim_nums(1,2)+1):scale.spdim_nums(2,2)),scale); % normalize subproblem boundary condition values
        end
    end
end