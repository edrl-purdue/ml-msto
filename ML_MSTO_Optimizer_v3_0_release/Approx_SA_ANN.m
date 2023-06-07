% SENSITIVITY COEFFICIENT APPROXIMATION - ANALYTICAL DERIVATIVE THROUGH ANN
function [dc, satime] = Approx_SA_ANN(Xin,macro,micro,Ue_int_Macro,Fe_int_Macro,ML_model,U_Macro)
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
    [~, ch_index, phys_index] = Indexing('full',micro,micro.ortho);
    [ml_index, ~, ~] = Indexing('part',micro,micro.ortho);
    TOM_e = [xMacro, xWeight];
    for p = 1:length(micro.phys)
        if micro.bc_type == 1
            TOM_e = [TOM_e, Normalize_TOM_BC(Fe_int_Macro{micro.phys(p)},micro)]; % TOM features
        elseif micro.bc_type == 2
            TOM_e = [TOM_e, Normalize_TOM_BC(Ue_int_Macro{micro.phys(p)},micro)]; % TOM features
        end
    end
    Cminvec = micro.Cminvec; Cmaxvec = micro.Cmaxvec;

    %% ANN TRANSFER FUNCTIONS
    tansig_n = @(yn_1,wn,bn) 2./(1+exp(-2*(yn_1*wn+bn))) - 1; % inputs are: y_(n-1), w_n, b_n; output: F_n for Fn = tansig
    purelin_n = @(yn_1,wn,bn) yn_1*wn+bn; % inputs are: y_(n-1), w_n, b_n; output: F_n for Fn = purelin
    
    %% DERIVATIVES OF ANN TRANSFER FUNCTIONS
    dtansig_n = @(yn_1,wn,bn) (4*wn.*exp(-2*(yn_1*wn+bn)))./(1+exp(-2*(yn_1*wn+bn))).^2; % inputs are: y_(n-1), w_n, b_n; output: dF_n/dy_(n-1) for Fn = tansig
    dpurelin_n = @(yn_1,wn,bn) wn; % inputs are: y_(n-1), w_n, b_n; output: dF_n/dy_(n-1) for Fn = purelin

    %% INITIALIZE ANN FOR DERIVATIVE CALCULATION
    % Initialize ANN parameters (ANN derivative code adapted from mathworks.com/matlabcentral/answers/398531-how-to-manually-calculate-a-neural-network-output)
    w = cellfun(@transpose,[ML_model.IW{1},ML_model.LW(2:size(ML_model.LW,1)+1:end)],'UniformOutput',false); % weight matrices
    b = cellfun(@transpose,ML_model.b','UniformOutput',false); % bias vectors
    tf = cellfun(@(x)x.transferFcn,ML_model.layers','UniformOutput',false); % transfer functions
    dtf = cell(size(tf)); % derivative of transfer functions
    for f = 1:numel(dtf)
        if strcmp(tf{f},'tansig')
            tf{f} = 'tansig_n'; dtf{f} = 'dtansig_n';
        elseif strcmp(tf{f},'purelin')
            tf{f} = 'purelin_n'; dtf{f} = 'dpurelin_n';
        end
    end
    % Prepare mapminmax on TOM_in (pre-process)
    if any(cellfun(@(c) strcmp(c,'mapminmax'),ML_model.Inputs{1}.processFcns))
        gain_in = ML_model.Inputs{1}.processSettings{cellfun(@(c) strcmp(c,'mapminmax'),ML_model.Inputs{1}.processFcns)}.gain; % gain of input mapminmax processSetting
        xmin_in = ML_model.Inputs{1}.processSettings{cellfun(@(c) strcmp(c,'mapminmax'),ML_model.Inputs{1}.processFcns)}.xmin; % xmin of input mapminmax processSetting
        ymin_in = ML_model.Inputs{1}.processSettings{cellfun(@(c) strcmp(c,'mapminmax'),ML_model.Inputs{1}.processFcns)}.ymin; % ymin of input mapminmax processSetting
    else
        gain_in = []; xmin_in = []; ymin_in = [];
    end
    % Prepare inverse mapminmax on CH_out (post-process)
    if any(cellfun(@(c) strcmp(c,'mapminmax'),ML_model.Outputs{end}.processFcns))
        gain_in = ML_model.Inputs{1}.processSettings{cellfun(@(c) strcmp(c,'mapminmax'),ML_model.Inputs{1}.processFcns)}.gain; % gain of input mapminmax processSetting
        gain_out = ML_model.outputs{end}.processSettings{cellfun(@(c) strcmp(c,'mapminmax'),ML_model.Outputs{end}.processFcns)}.gain; % gain of output mapminmax processSetting
        ymin_out = ML_model.outputs{end}.processSettings{cellfun(@(c) strcmp(c,'mapminmax'),ML_model.Outputs{end}.processFcns)}.ymin; % ymin of input mapminmax processSetting
        xmin_out = ML_model.outputs{end}.processSettings{cellfun(@(c) strcmp(c,'mapminmax'),ML_model.Outputs{end}.processFcns)}.xmin; % xmin of input mapminmax processSetting
    else
        gain_out = []; xmin_out = []; ymin_out = [];
    end
    
    %% APPROXIMATE SENSITIVITY COEFFICIENTS
    dKE_Hom = cell(numel(xMacro),1);
    WaitMessage = parfor_waitbar(numel(xMacro)); %Add ('Waitbar', true) as arguments to turn on a waitbar popup
    fprintf('     |        Approximate Sensitivity Coefficients       |\n'); tic
%     fprintf('     |0---10---20---30---40---50---60---70---80---90-100%%|\n'); tic
    parfor e = 1:numel(xMacro)
        %% DEFINE TOM INPUT PARAMETERS
        if all(micro.phys == [1 2])
            TOM_in = TOM_e(e,:)'; % includes VF, obj weight, and TOM BCs
        else
            TOM_in = [TOM_e(e,1)'; TOM_e(e,3:end)']; % includes only VF and TOM BCs
        end

        TOM_in0 = TOM_in;
        if micro.vf_cutoff > 0
            TOM_in(TOM_in0(1)<=micro.vf_cutoff,1) = micro.vf_cutoff;
        end

        %% DERIVATIVE OF HOMOGENIZED CONSTITUTIVE TENSORS WITH RESPECT TO VF
        % Perform mapminmax on TOM_in (pre-process)
        if isempty(gain_in)
            y0 = TOM_in; % input to 1st hidden layer, no pre-processing
        else
            y0 = bsxfun(@plus,bsxfun(@times,bsxfun(@minus,TOM_in,xmin_in),gain_in),ymin_in); % input to 1st hidden layer
        end

        % Loop through hidden layers + output layer
        y = cell(1,length(w)); % output of layer n
        dy = y; % derivative of y_n with respect to y_(n-1)
        y{1} = tansig_n(transpose(y0),w{1},b{1}); % 1st hidden layer output
        dy{1} = dtansig_n(transpose(y0),w{1},b{1}); dy_product = dy{1}; % 1st hidden layer output's derivative
        for n=2:length(w)
            if strcmp(tf{f},'tansig_n')
                y{n} = tansig_n(y{n-1},w{n},b{n}); % hidden layer n output
                dy{n} = dtansig_n(y{n-1},w{n},b{n}); % derivative of hidden layer n output
            elseif strcmp(tf{f},'purelin_n')
                y{n} = purelin_n(y{n-1},w{n},b{n}); % hidden layer n output
                dy{n} = dpurelin_n(y{n-1},w{n},b{n}); % derivative of hidden layer n output
            end
            dy_product = dy_product*dy{n}; % product ANN derivative via chain rule
        end
        
        % Perform inverse mapminmax on CH_out (post-process)
        if isempty(gain_out)
            CH_out = y{end}; % output of output layer, no post-processing
            dCH_out = dy_product; % derivative of output layer with respect to y0
        else
            CH_out = transpose((y{end} - ymin_out)./transpose(gain_out) + transpose(xmin_out)); % output of output layer
            dCH_out = transpose(bsxfun(@rdivide,bsxfun(@times,dy_product,gain_in),transpose(gain_out))); % derivative of output layer with respect to y0
        end

        dC_Hom = cell(3,3);
        dC_Homvec = dCH_out(:,1)'; % no capping of dc
        if TOM_in0(1) <= micro.vf_cutoff
            dC_Homvec = (CH_out' - Cminvec).*(1/micro.vf_cutoff).^1; % vf_cutoff derivative
        end
        for p = 1:size(phys_index,1)
            dC_Hom{phys_index(p,1),phys_index(p,2)} = zeros(max(ch_index{phys_index(p,1),phys_index(p,2)}(:,1)),max(ch_index{phys_index(p,1),phys_index(p,2)}(:,2)));
            for ind = 1:size(ch_index{phys_index(p,1),phys_index(p,2)},1)
                % Derivative of homogenized constitutive tensor with respect to VF
                dC_Hom{phys_index(p,1),phys_index(p,2)}(ch_index{phys_index(p,1),phys_index(p,2)}(ind,1),ch_index{phys_index(p,1),phys_index(p,2)}(ind,2)) = dC_Homvec(all(sort(ch_index{1,1}(ind,:),'descend') == ml_index{1,1},2));
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
    fprintf('     |         Done. Elapsed Time: %4.2f seconds.         |\n',toc); satime = toc;
    
    % Calculate analytical derivative of homogenized macroscale using dKE_Hom
    [~, dc, ~] = Obj_Fns([],macro,U_Macro,[],[],[],[],dKE_Hom,0.5); % Case 3: Calling Obj_Fns([],macro,U_Macro,[],[],[],[],dKE_Hom,weight_e) is dc for Homogenized stiffness matrix
    dc(dc>0) = -dc(dc>0);
end