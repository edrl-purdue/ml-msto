% SENSITIVITY COEFFICIENT APPROXIMATION - CENTRAL FINITE DIFFERENCES
function [dc, satime] = Approx_SA_CFD(Xin,macro,micro,mat,Ue_int_Macro,Fe_int_Macro,ML_model,U_Macro)
    %% INITIALIZE TOM DATA
    if all(macro.phys == [1 2])
        xMacro = Xin(1:length(Xin)/2);
        xWeight = Xin(length(Xin)/2+1:end);
        num_dv = 2;
        macro.h_fd = [(macro.x_ub(1) - macro.x_lb(1))*macro.h_fd(1), (macro.x_ub(2) - macro.x_lb(2))*macro.h_fd(2)];
    else
        xMacro = Xin;
        xWeight = ones(size(xMacro));
        num_dv = 1;
        macro.h_fd(1) = (macro.x_ub(1) - macro.x_lb(1))*macro.h_fd(1);
    end
    x_lb = macro.x_lb; x_ub = macro.x_ub;
    [~, ch_index, phys_index] = Indexing('full',micro,1);
    [ml_index, ~, ~] = Indexing('part',micro,1);
    TOM_e = [xMacro, xWeight];
    for p = 1:length(micro.phys)
        if micro.bc_type == 1
            TOM_e = [TOM_e, Normalize_TOM_BC(Fe_int_Macro{micro.phys(p)},micro)]; % TOM features
        elseif micro.bc_type == 2
            TOM_e = [TOM_e, Normalize_TOM_BC(Ue_int_Macro{micro.phys(p)},micro)]; % TOM features
        end
    end
    
    %% PERTURB VOLUME FRACTION
    TOM_e_ph = []; TOM_e_mh = [];
    for v = 1:num_dv
        TOM_temp = TOM_e; TOM_temp(:,v) = TOM_temp(:,v) + macro.h_fd(v); % plus H pertubation
        TOM_temp(:,v) = min(TOM_temp(:,v),macro.x_ub(v)); % if perturbed volume fraction is outside of macroscale design variable bounds then set volume fraction to bounds
        TOM_e_ph = [TOM_e_ph; TOM_temp];
        TOM_temp = TOM_e; TOM_temp(:,v) = TOM_temp(:,v) - macro.h_fd(v); % minus H pertubation
        TOM_temp(:,v) = max(TOM_temp(:,v),macro.x_lb(v)); % if perturbed volume fraction is outside of macroscale design variable bounds then set volume fraction to bounds
        TOM_e_mh = [TOM_e_mh; TOM_temp];
    end
    if micro.ml
        [C_lower_HS, C_upper_HS] = Prepare_HS_bounds(micro.mat,micro.dim,0);
        C_Homvec_ph = Cap_ANN(Interp_CH(ML_model, TOM_e_ph, micro),C_lower_HS(TOM_e_ph(:,1)'),C_upper_HS(TOM_e_ph(:,1)'))';
        C_Homvec_mh = Cap_ANN(Interp_CH(ML_model, TOM_e_mh, micro),C_lower_HS(TOM_e_mh(:,1)'),C_upper_HS(TOM_e_mh(:,1)'))';
    end
    
    %% CONTINUATION INITIALIZATION FOR SUBPROBLEMS
    if micro.ml == 0 && micro.cont == 1 && macro.cont == 0 % if full TO is used and TOM continutation is on and Macro continuation is off
        penal_k = 1; % p0
        del_penal = 0.3; % del_p
        penal_k = penal_k:del_penal:micro.penal; % penalization values
        if penal_k(end) ~= micro.penal; penal_k(end+1) = micro.penal; end % add last penalization value is max
        tolx_k = 1e-2; % w0
        tolx_mult = (micro.tolx/tolx_k)^(1/(length(penal_k)-1)); % w multiplier
        for i = 1:(length(penal_k)-1)
            tolx_k = [tolx_k min(tolx_k)*tolx_mult]; %#ok<AGROW>
        end
        loop_k = ones(length(tolx_k))*micro.maxloop/10; loop_k(end) = micro.maxloop;
    else
        penal_k = micro.penal;
        tolx_k = micro.tolx;
        loop_k = micro.maxloop;
    end

    %% APPROXIMATE SENSITIVITY COEFFICIENTS
    dKE_Hom = cell(numel(xMacro),1);
    WaitMessage = parfor_waitbar(numel(xMacro)*num_dv); %Add ('Waitbar', true) as arguments to turn on a waitbar popup
    fprintf('     |        Approximate Sensitivity Coefficients       |\n'); tic
%     fprintf('     |0---10---20---30---40---50---60---70---80---90-100%%|\n'); tic
    parfor e = 1:numel(xMacro)
        if mod(e,numel(xMacro)) == 0
            eind = numel(xMacro); % element index for multiple design variables num_dv>=2
        else
            eind = mod(e,numel(xMacro)); % element index for multiple design variables num_dv>=2
        end
        
        %% DERIVATIVE OF HOMOGENIZED STIFFNESS TENSORS WITH RESPECT TO VF
        dC_Hom = cell(3,3);
        if xMacro(eind) < micro.vf_cutoff
            if any(micro.phys == 1)
                dC_Hom{1,1} = (mat.Emin + 3*(mat.Emax - mat.Emin)*xMacro(e)^2) * Prepare_C(mat.nu,1,micro.dim);
            end
            if any(micro.phys == 2)
                dC_Hom{2,2} = (mat.kmin + 3*(mat.kmax - mat.kmin)*xMacro(e)^2) * Prepare_C(mat.nu,2,micro.dim);
            end
        else
            dC_Homvec = (C_Homvec_ph(:,e) - C_Homvec_mh(:,e))/(2*macro.h_fd(ceil(e/numel(xMacro))));
            for p = 1:size(phys_index,1)
                dC_Hom{phys_index(p,1),phys_index(p,2)} = zeros(max(ch_index{phys_index(p,1),phys_index(p,2)}(:,1)),max(ch_index{phys_index(p,1),phys_index(p,2)}(:,2)));
                for ind = 1:size(ch_index{phys_index(p,1),phys_index(p,2)},1)
                    % Derivative of homogenized constitutive tensor with respect to VF
                    dC_Hom{phys_index(p,1),phys_index(p,2)}(ch_index{phys_index(p,1),phys_index(p,2)}(ind,1),ch_index{phys_index(p,1),phys_index(p,2)}(ind,2)) = dC_Homvec(all(sort(ch_index{1,1}(ind,:),'descend') == ml_index{1,1},2));
                end
            end
        end
        
        %% DERIVATIVE OF HOMOGENIZED ELEMENTAL STIFFNESS MATRIX WITH RESPECT TO VF
        dKE_Hom{e} = cell(3,3); % # of rows/cols is the number of unique physics modes || KE_ij: coupling stiffness matrices of physics i to physics j
        if any(micro.phys == 1)
            [dKE_Hom{e}{1,1}, ~, ~] = Prepare_KE(dC_Hom{1,1},1,micro); % Mechanical
            dKE_Hom{e}{1,1} = reshape(dKE_Hom{e}{1,1},[micro.nen*micro.dim,micro.nen*micro.dim]);
        end
        if any(micro.phys == 2)
            [dKE_Hom{e}{2,2}, ~, ~] = Prepare_KE(dC_Hom{2,2},2,micro); % Thermal
            dKE_Hom{e}{2,2} = reshape(dKE_Hom{e}{2,2},[micro.nen,micro.nen]);
        end
        if all(micro.phys == [1 2])
            [dKE_Hom{e}{2,1}, ~, ~] = Prepare_KE(dC_Hom{1,1},[1,2],micro); % Thermomechanical
            dKE_Hom{e}{2,1} = reshape(dKE_Hom{e}{2,1},[micro.nen*micro.dim,micro.nen]);
        end

        WaitMessage.Send;
    end
    WaitMessage.Destroy;
    fprintf('     |         Done. Elapsed Time: %4.2f seconds.         |\n\n',toc); satime = toc;

    % Calculate analytical derivative of homogenized macroscale using dKE_Hom
    [~, dc, ~] = Obj_Fns([],macro,U_Macro,[],[],[],[],dKE_Hom,0.5); % Case 3: Calling Obj_Fns([],macro,U_Macro,[],[],[],[],dKE_Hom,weight_e) is dc for Homogenized stiffness matrix
    dc(dc>0) = -dc(dc>0);
end