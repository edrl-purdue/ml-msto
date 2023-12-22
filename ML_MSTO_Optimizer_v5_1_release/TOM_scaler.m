% SCALE MICROSCALE MESH RESOLUTION (NOT COMPLETED YET, PLACE-HOLDER CODE)
function xMicro_scaled = TOM_scaler(macro,micro,xMicro)
    %% SCALING TO SPECIFIED EPI SIZE
    [Ye_c, Xe_c, Ze_c] = meshgrid((-macro.nely/2 + 1/2):1:(macro.nely/2 - 1/2),(-macro.nelx/2 + 1/2):1:(macro.nelx/2 - 1/2),min(-macro.nelz/2 + 1/2,0):1:max(macro.nelz/2 - 1/2,0)); % middle of coarse, unscaled element coordinates on point-wise coordinates
    [ye, xe, ze] = meshgrid((-0.5 + 0.5/micro.nely):(1/micro.nely):(0.5 - 0.5/micro.nely),(-0.5 + 0.5/micro.nelx):(1/micro.nelx):(0.5 - 0.5/micro.nelx),min(-0.5 + 0.5/micro.nelz,0):(1/micro.nelz):max(0.5 - 0.5/micro.nelz,0)); % middle of resolution-scaled microscale element coordinates on global reference

    xMicro_data = cell(size(xMicro));
    for e = 1:macro.nele
        xMicro_data{e} = zeros(micro.nele,4);
        xMicro_data{e}(:,1) = Xe_c(e) + xe(:);
        xMicro_data{e}(:,2) = Ye_c(e) + ye(:);
        xMicro_data{e}(:,3) = Ze_c(e) + ze(:);
        xMicro_data{e}(:,4) = xMicro{e}(:);
    end

    macro_scaled = macro;
    [Ye_f, Xe_f, Ze_f] = meshgrid((-macro.nely/2 + macro.epi/2):macro.epi:(macro.nely/2 - macro.epi/2),(-macro.nelx/2 + macro.epi/2):macro.epi:(macro.nelx/2 - macro.epi/2),min(-macro.nelz/2 + macro.epi/2,0):macro.epi:max(macro.nelz/2 - macro.epi/2,0)); % middle of coarse, unscaled element coordinates on point-wise coordinates
    macro_scaled.nelx = macro.nelx*(1/macro.epi); macro_scaled.nely = macro.nely*(1/macro.epi); macro_scaled.nelz = macro.nelz*(1/macro.epi);
    [macro_scaled.epi_necon, ~] = Prepare_Mesh(macro_scaled); % Macro mesh connectivity
    macro_scaled = mesh_info(macro_scaled);

    [ye, xe, ze] = meshgrid((-0.5 + 0.5/micro.nely)*macro.epi:(1/micro.nely)*macro.epi:(0.5 - 0.5/micro.nely)*macro.epi,(-0.5 + 0.5/micro.nelx)*macro.epi:(1/micro.nelx)*macro.epi:(0.5 - 0.5/micro.nelx)*macro.epi,min(-0.5 + 0.5/micro.nelz,0)*macro.epi:(1/micro.nelz)*macro.epi:max(0.5 - 0.5/micro.nelz,0)*macro.epi); % middle of resolution-scaled microscale element coordinates on global reference

    xMicro_scaled = cell(size(xMicro)*(1/macro.epi));
    xMicro_scaled_data = cell(size(xMicro)*(1/macro.epi));
    for e = 1:macro_scaled.nele
        xMicro_scaled_data{e} = zeros(micro.nele,3);
        xMicro_scaled_data{e}(:,1) = Xe_f(e) + xe(:);
        xMicro_scaled_data{e}(:,2) = Ye_f(e) + ye(:);
        xMicro_scaled_data{e}(:,3) = Ze_f(e) + ze(:);
    end
    
    xMicro_data_x = cell(micro.nele,1); xMicro_data_y = cell(micro.nele,1); xMicro_data_z = cell(micro.nele,1); xMicro_data_s = cell(micro.nele,1);
    xMicro_scaled_data_x = cell(micro.nele,1); xMicro_scaled_data_y = cell(micro.nele,1); xMicro_scaled_data_z = cell(micro.nele,1); xMicro_scaled_data_s = cell(micro.nele,1);
    parfor s = 1:micro.nele
        xMicro_data_x{s} = zeros(size(xMicro_data)); xMicro_data_y{s} = zeros(size(xMicro_data)); xMicro_data_z{s} = zeros(size(xMicro_data)); xMicro_data_s{s} = zeros(size(xMicro_data));
        for e = 1:macro.nele
            xMicro_data_x{s}(e) = xMicro_data{e}(s,1);
            xMicro_data_y{s}(e) = xMicro_data{e}(s,2);
            xMicro_data_z{s}(e) = xMicro_data{e}(s,3);
            xMicro_data_s{s}(e) = xMicro_data{e}(s,4);
        end
        xMicro_scaled_data_x{s} = zeros(size(xMicro_scaled_data)); xMicro_scaled_data_y{s} = zeros(size(xMicro_scaled_data)); xMicro_scaled_data_z{s} = zeros(size(xMicro_scaled_data)); xMicro_scaled_data_s{s} = zeros(size(xMicro_scaled_data));
        for e = 1:macro_scaled.nele
            xMicro_scaled_data_x{s}(e) = xMicro_scaled_data{e}(s,1);
            xMicro_scaled_data_y{s}(e) = xMicro_scaled_data{e}(s,2);
            xMicro_scaled_data_z{s}(e) = xMicro_scaled_data{e}(s,3);
        end
        if macro.dim == 2
            xMicro_scaled_data_s{s} = interp2(xMicro_data_y{s},xMicro_data_x{s},xMicro_data_s{s},max(min(xMicro_scaled_data_y{s},max(xMicro_data_y{s}(:))),min(xMicro_data_y{s}(:))),max(min(xMicro_scaled_data_x{s},max(xMicro_data_x{s}(:))),min(xMicro_data_x{s}(:))),'linear'); % X_scaled and Y_scaled coordinates are rounded to be on boundary inside this interp function call so that a nearest neighbor boundary condition is effectively applied.
        elseif macro.dim == 3
            xMicro_scaled_data_s{s} = interp3(xMicro_data_y{s},xMicro_data_x{s},xMicro_data_z{s},xMicro_data_s{s},max(min(xMicro_scaled_data_y{s},max(xMicro_data_y{s}(:))),min(xMicro_data_y{s}(:))),max(min(xMicro_scaled_data_x{s},max(xMicro_data_x{s}(:))),min(xMicro_data_x{s}(:))),max(min(xMicro_scaled_data_z{s},max(xMicro_data_z{s}(:))),min(xMicro_data_z{s}(:))),'linear'); % X_scaled, Y_scaled, and Z_scaled coordinates are rounded to be on boundary inside this interp function call so that a nearest neighbor boundary condition is effectively applied.
        end
    end

    
    parfor e = 1:macro_scaled.nele
        xMicro_scaled{e} = ini_design(micro,0,micro.x_lb(1),micro.x_ub(1));
        for s = 1:micro.nele
            xMicro_scaled{e}(s) = xMicro_scaled_data_s{s}(e);
        end
    end
    
end