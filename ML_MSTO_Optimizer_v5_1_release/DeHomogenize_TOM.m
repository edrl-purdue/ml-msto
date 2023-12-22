% DE-HOMOGENIZE MACROSCALE WITH FINAL SET OF TOMs
function xMicro = DeHomogenize_TOM(xMacro,xWeight,macro,micro,mat,KE,opt_history)
    %% INITIALIZE TOM DATA
    Ue_int_Macro = opt_history{end,9};
    Fe_int_Macro = opt_history{end,10};
    TOM_e = [xMacro, xWeight];
    for p = 1:length(micro.phys)
        if micro.bc_type == 1
            TOM_e = [TOM_e, Normalize_TOM_BC(Fe_int_Macro{micro.phys(p)},micro)]; % TOM features
        elseif micro.bc_type == 2
            TOM_e = [TOM_e, Normalize_TOM_BC(Ue_int_Macro{micro.phys(p)},micro)]; % TOM features
        end
    end

    %% CONTINUATION INITIALIZATION FOR SUBPROBLEMS
    if micro.cont == 1 && macro.cont == 0 % if TOM continutation is on and Macro continuation is off
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

    %% MICROSCALE DE-HOMOGENIZATION
    xMicro = cell(numel(xMacro),1); cMicro = cell(numel(xMacro),1);
    WaitMessage = parfor_waitbar(numel(xMacro)); %Add ('Waitbar', true) as arguments to turn on a waitbar popup
    fprintf('     |   Solve Final Inner Optimization Problems (TOMs)  |\n'); tic
    parfor e = 1:numel(xMacro)
        if xMacro(e) < 0.02
            xMicro{e} = zeros(micro.nelx,micro.nely,max(1,micro.nelz));
        else
            xMicro{e} = ini_design(micro,xMacro(e),macro.x_lb(1),macro.x_ub(1));
            for pk = 1:numel(penal_k)
                [xMicro{e}, cMicro{e}] = Micro_TO(micro,TOM_e(e,:),xMicro{e},penal_k(pk),loop_k(pk),tolx_k(pk),mat,KE);
            end
            if ~isreal(sum(sum(sum(xMicro{e}))))
                xMicro{e} = zeros(micro.nelx,micro.nely,max(1,micro.nelz));
            end
        end
        WaitMessage.Send;
    end
    WaitMessage.Destroy;
    fprintf('     |         Done. Elapsed Time: %4.2f seconds.         |\n',toc);
    xMicro = reshape(xMicro,[macro.nelx,macro.nely,max(macro.nelz,1)]);
end