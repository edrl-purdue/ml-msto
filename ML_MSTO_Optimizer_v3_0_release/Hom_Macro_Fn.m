% HOMOGENIZED MACROSCALE OBJECTIVE FUNCTION
function [c, dc] = Hom_Macro_Fn(Xin,macro,micro,mat,KE,FE,F0_Macro,U0_Macro,N_Macro,ML_model)
    %% LOAD TEMP FILES AND ASSEMBLE OPT_HISTORY
    temp_files = dir(['temp_files/temp_' num2str(macro.nelx) 'x' num2str(macro.nely) 'x' num2str(macro.nelz) '_' num2str(macro.rand_id) '_*']);
    load(['temp_files/' temp_files(end).name],'opt_history')
    
    %% INITIALIZE DESIGN VECTORS
    if all(macro.phys == [1 2])
        xMacro = Xin(1:length(Xin)/2);
        xWeight = Xin(length(Xin)/2+1:end);
        xMacro(:) = (macro.den_H*xMacro(:))./macro.den_Hs; % apply density filter
        Xin = [xMacro(:); xWeight(:)];
    else
        xMacro = Xin;
        xWeight = ones(size(xMacro));
        xMacro(:) = (macro.den_H*xMacro(:))./macro.den_Hs; % apply density filter
        Xin = xMacro(:);
    end
    xMacro(:) = (macro.den_H*xMacro(:))./macro.den_Hs; % apply density filter
    
    %% SAVE ITERATION INFORMATION
    opt_history{end+1,2} = reshape(xMacro,[macro.nelx, macro.nely, max(macro.nelz,1)]);
    opt_history{end,3} = reshape(xWeight,[macro.nelx, macro.nely, max(macro.nelz,1)]);
    
    %% INITIALIZE TOM DATA
    [~, ch_index, phys_index] = Indexing('full',micro,micro.ortho);
    [ml_index, ~, ~] = Indexing('part',micro,micro.ortho);
    Ue_int_Macro = opt_history{end-1,9};
    Fe_int_Macro = opt_history{end-1,10};
    U_Macro = opt_history{end-1,11};
    TOM_e = [xMacro, xWeight];
    for p = 1:length(micro.phys)
        if micro.bc_type == 1
            TOM_e = [TOM_e, Normalize_TOM_BC(Fe_int_Macro{micro.phys(p)},micro)]; % TOM features
        elseif micro.bc_type == 2
            TOM_e = [TOM_e, Normalize_TOM_BC(Ue_int_Macro{micro.phys(p)},micro)]; % TOM features
        end
    end
    if micro.ml
        [C_lower_HS, C_upper_HS] = Prepare_HS_bounds(micro.mat,micro.dim,micro.ortho);
        C_Homvec = Cap_ANN(Interp_CH(ML_model, TOM_e, micro),C_lower_HS(TOM_e(:,1)'),C_upper_HS(TOM_e(:,1)'))'; % interpolate and cap ANN predictions
    else
        C_Homvec = [];
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
        loop_k = ones(size(tolx_k))*micro.maxloop/10; loop_k(end) = micro.maxloop;
    else
        penal_k = micro.penal;
        tolx_k = micro.tolx;
        loop_k = micro.maxloop;
    end
    
    %% GENERATE SUBPROBLEMS
    xMicro = cell(numel(xMacro),1);
    C_Hom = cell(numel(xMacro),1);
    KE_Hom = cell(numel(xMacro),1);
    WaitMessage = parfor_waitbar(numel(xMacro)); %Add ('Waitbar', true) as arguments to turn on a waitbar popup
    fprintf('     |It.:%5i                                          |\n',opt_history{end-1,1}+1);
    fprintf('     |   Approximate Inner Optimization Problems (TOMs)  |\n'); tic
    parfor e = 1:numel(xMacro)
        %% HOMOGENIZED CONSTITUTIVE TENSORS
        C_Hom{e} = cell(3,3);
        if micro.ml % ANN Prediction
            for p = 1:size(phys_index,1)
                C_Hom{e}{phys_index(p,1),phys_index(p,2)} = zeros(max(ch_index{phys_index(p,1),phys_index(p,2)}(:,1)),max(ch_index{phys_index(p,1),phys_index(p,2)}(:,2)));
                for ind = 1:size(ch_index{phys_index(p,1),phys_index(p,2)},1)
                    C_Hom{e}{phys_index(p,1),phys_index(p,2)}(ch_index{phys_index(p,1),phys_index(p,2)}(ind,1),ch_index{phys_index(p,1),phys_index(p,2)}(ind,2)) = C_Homvec(all(sort(ch_index{1,1}(ind,:),'descend') == ml_index{1,1},2),e);
                end
            end
        else % No ML
            if xMacro(e) < micro.nelx/micro.nele
                if micro.ml == 0; xMicro{e} = zeros(micro.nelx,micro.nely,max(1,micro.nelz)); end
                if any(micro.phys == 1)
                    C_Hom{e}{1,1} = mat.Emin * Prepare_C(mat.nu,1,micro.dim);
                end
                if any(micro.phys == 2)
                    C_Hom{e}{2,2} = mat.kmin * Prepare_C(mat.nu,2,micro.dim);
                end
            else
                xMicro{e} = ini_design(micro,xMacro(e),macro.x_lb(1),macro.x_ub(1));
                for pk = 1:numel(penal_k)
                    xMicro{e} = Micro_TO(micro,TOM_e(e,:),xMicro{e},penal_k(pk),loop_k(pk),tolx_k(pk),mat,KE);
                end
                C_Hom{e} = Num_Hom(xMicro{e},micro,mat,KE,FE); % returns C_Hom{2,1} but this is not needed to build any stiffness matrix for Hom_FEA
            end
        end
        [~, C_Hom{e}] = Check_CH(C_Hom{e},micro,mat,xMacro(e)); % returns 1 then CHi is valid, returns 0 then CHi is not valid
        
        %% HOMOGENIZED ELEMENTAL STIFFNESS MATRIX
        KE_Hom{e} = cell(3,3); % # of rows/cols is the number of unique physics modes || KE_ij: coupling stiffness matrices of physics i to physics j
        if any(micro.phys == 1)
            [KE_Hom{e}{1,1}, ~, ~] = Prepare_KE(C_Hom{e}{1,1},1,micro); % Mechanical
            KE_Hom{e}{1,1} = reshape(KE_Hom{e}{1,1},[micro.nen*micro.dim,micro.nen*micro.dim]);
        end
        if any(micro.phys == 2)
            [KE_Hom{e}{2,2}, ~, ~] = Prepare_KE(C_Hom{e}{2,2},2,micro); % Thermal
            KE_Hom{e}{2,2} = reshape(KE_Hom{e}{2,2},[micro.nen,micro.nen]);
        end
        if all(micro.phys == [1 2])
            [~ , ~, KE_Hom{e}{2,1}] = Prepare_KE(C_Hom{e}{1,1},[1,2],micro); % Thermomechanical
            KE_Hom{e}{2,1} = reshape(KE_Hom{e}{2,1},[micro.nen*micro.dim,micro.nen]);
        end
        WaitMessage.Send;
    end
    WaitMessage.Destroy;
    fprintf('     |         Done. Elapsed Time: %4.2f seconds.         |\n',toc); ptomtime = toc;
    opt_history{end,7} = C_Hom; opt_history{end,8} = KE_Hom;
    
    %% SENSITIVITY ANALYSIS
    if macro.h_fd == 0
        [dc, satime] = Approx_SA_ANN(Xin,macro,micro,Ue_int_Macro,Fe_int_Macro,ML_model,U_Macro);
    else
        [dc, satime] = Approx_SA_CFD(Xin,macro,micro,mat,Ue_int_Macro,Fe_int_Macro,ML_model,U_Macro);
    end
    if all(macro.phys == [1 2])
        dc(1:length(Xin)/2) = (macro.den_H*dc(1:length(Xin)/2))./macro.den_Hs; % apply density filter
        dc(1:length(Xin)/2) = macro.sen_H*(xMacro(:).*dc(1:length(Xin)/2))./macro.sen_Hs./max(1e-3,xMacro(:)); % filter Macroscale sensitivities
    else
        dc(:) = (macro.den_H*dc(:))./macro.den_Hs; % apply density filter
        dc(:) = macro.sen_H*(xMacro(:).*dc(:))./macro.sen_Hs./max(1e-3,xMacro(:)); % filter Macroscale sensitivities
    end
    opt_history{end,5} = dc; opt_history{end,6} = [ptomtime, satime];

    %% HOMOGENIZED FINITE ELEMENT ANALYSIS
    [U_Macro,Ue_int_Macro,Fe_int_Macro] = Hom_FEA(KE_Hom,macro,F0_Macro,U0_Macro,N_Macro);
    opt_history{end,9} = Ue_int_Macro; opt_history{end,10} = Fe_int_Macro; opt_history{end,11} = U_Macro;
    
    %% OBJECTIVE FUNCTION EVALUATION
    [c, ~, ~] = Obj_Fns([],macro,U_Macro,[],[],[],KE_Hom,[],0.5); % Case 2: Calling Obj_Fns([],macro,U_Macro,[],[],[],KE_Hom,[],weight_e) is c for Homogenized stiffness matrix
    opt_history{end,4} = c;

    %% PLOT MACROSCALE ITERATION
    opt_history{end,1} = opt_history{end-1,1} + 1; % Macroscale iteration number
    output_plot(xMacro,xWeight,macro,xMicro,c,opt_history{end,1},max(abs(opt_history{end,2}(:) - opt_history{end-1,2}(:))));

    %% SAVE LAST TEMP FILE AND OPT_HISTORY PART
    save_size = 10; % number of iterations to save per temp file (keeps file size low)
    if size(opt_history,1) > save_size
        fnum = str2double(temp_files(end).name(end-6:end-4))+1;
        if length(int2str(fnum)) == 1
            zeros_num = '00';
        elseif length(int2str(fnum)) == 2
            zeros_num = '0';
        elseif length(int2str(fnum)) == 3
            zeros_num = [];
        end
        temp_files(end + 1).name = [temp_files(end).name(1:end-7) zeros_num num2str(fnum) '.mat'];
        opt_history(1:save_size,:) = [];
    end
    save(['temp_files/' temp_files(end).name],'opt_history')
end