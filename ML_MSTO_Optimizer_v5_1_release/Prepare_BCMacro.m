% PREPARE MACRO BOUNDARY CONDITIONS
function [scale,F,U,N] = Prepare_BCMacro(scale,phys)
    %% DEFINE MACROSCALE BOUNDARY CONDITIONS
    % scale.bc{phys = (1,2,3),1} = {'name', ##, [magx magy magz], optional} % define loads
    % scale.bc{phys = (1,2,3),2} = {'name', ##, dofs, optional} % define supports
    if scale.dim == 2
        if phys == 1 % 2D Mechanical
            scale.bc{1,1} = {'corn', 22, [0 -1 0], ''}; scale.bc{1,2} = {'edge', 13, 12, ''}; scale.macro_bc_name = 'MBB Beam';
            % scale.bc{1,1} = {'corn', 33, [0 -1 0], ''}; scale.bc{1,2} = {'edge', 13, 1, ''; 'corn', 22, 2, ''}; scale.macro_bc_name = 'Half MBB Beam';
            % scale.bc{1,1} = {'mide', 24, [0 -1 0], ''}; scale.bc{1,2} = {'edge', 13, 12, ''}; scale.macro_bc_name = 'Cantilever Beam';
            % scale.bc{1,1} = {'edge', 24, [0 -1 0], ''}; scale.bc{1,2} = {'edge', 13, 12, ''}; scale.macro_bc_name = 'Sheared Beam';
            % scale.bc{1,1} = {'midf', 11, [0 -1 0], ''}; scale.bc{1,2} = {'edge', 13, 12, ''; 'edge', 24, 12, ''}; scale.macro_bc_name = 'Clamped Beam with center load';
            % scale.bc{1,1} = {'edge', 34, [0 -1 0], ''; 'corn', 33, [0 -0.5 0], ''; 'corn', 44, [0 -0.5 0], ''}; scale.bc{1,2} = {'edge', 13, 12, ''; 'edge', 24, 12, ''}; scale.macro_bc_name = 'Clamped Beam with distributed top edge load';
            % scale.bc{1,1} = {'mide', 34, [0 -1 0], ''}; scale.bc{1,2} = {'corn', 11, 2, ''; 'corn', 22, 2, ''; 'midy', 11, 1, ''}; scale.macro_bc_name = 'Three-point Bending Beam';
            % scale.bc{1,1} = {'edge', 34, [0 -1 0], ''; 'corn', 33, [0 -0.5 0], ''; 'corn', 44, [0 -0.5 0], ''}; scale.bc{1,2} = {'edge', 12, 12, ''}; scale.macro_bc_name = 'Compressed Block (top edge load, bottom edge support)';
            % scale.bc{1,1} = {'mide', 12, [0 -1 0], ''}; scale.bc{1,2} = {'corn', 11, 2, ''; 'corn', 22, 2, ''; 'midy', 11, 1, ''}; scale.macro_bc_name = 'Bridge';
            % scale.bc{1,1} = {'corn', 11, [0 -1 0], ''}; scale.bc{1,2} = {'edge', 13, 1, ''; 'corn', 22, 2, ''}; scale.macro_bc_name = 'Wheel';
            % scale.bc{1,1} = {'mide', 24, [0 -1 0], ''}; scale.bc{1,2} = {'edge', 34, 12, '33_40'}; scale.macro_bc_name = 'Lbracket';
        end
    elseif scale.dim == 3
        if phys == 1 % 3D Mechanical
            % scale.bc{1,1} = {'edge', 24, [0 0 -1], ''}; scale.bc{1,2} = {'face', 55, 123, ''}; scale.macro_bc_name = 'MBB Beam';
            % scale.bc{1,1} = {'midf', 66, [0 0 -1], ''}; scale.bc{1,2} = {'face', 55, 123, ''}; scale.macro_bc_name = 'Cantilever Beam';
            scale.bc{1,1} = {'midf', 22, [0 0 -1], ''}; scale.bc{1,2} = {'corn', 11, 3, ''; 'corn', 22, 3, ''; 'corn', 33, 3, ''; 'corn', 44, 3, ''; 'mids', 34, 2, ''; 'mids', 56, 1, ''}; scale.macro_bc_name = 'Five-point Bending Plate';
        end
    end

    %% MESH INFORMATION
    [nnodes, dx, dy, dz, xyzc, lcoorp] = mesh_info_detailed(scale);
    
    %% INITIALIZE PHYSICS MODELS
    % Mechanical
    if phys == 1
        F = zeros(nnodes*3,1);
        U = zeros(nnodes*3,1);
        N = ones(nnodes*3,1);
    elseif phys == 2
        F = zeros(nnodes,1);
        U = zeros(nnodes,1);
        N = ones(nnodes,1);
    elseif phys == 3
        
    end
    
    %% BOUNDARY CONDITIONS
    %                                                                             NAME    FEATURE                                  ##
    %                 DOMAIN NODE AND FACE NUMBERINGS                            'midp' = Midpoint of domain (point)/(area/volume) 00
    %      2D Nodes               3D Nodes                  3D Faces             'corn' = Corner (point)/(area/volume)             2 x Corner #
    %    3 ________ 4         5 ________ 7             ________                  'edge' = Edge (line)/(partial line)               Corner # to Corner #
    %     |        |  ^+Y      |\       |\    ^+Z     |\   2   |\    1 = bottom  'face' = Face (surface)/(volume)                  2 x Face #
    %     |  1 =   |  |        |6\ _____|_\ 8 |       | \ _____|_\   2 = top     'mide' = Midpoint of edge (point)/(line)          Corner # to Corner #
    %     |  face  |  |        |  |     |  |  |       |3 | 5   |  |  3 = front   'midx' = x-Midline (line)/(area)                  2 x Face #
    %    1|________|2 +---->  1|__|_____|3 |  +---->  |__|_____| 4|  4 = back    'midy' = y-Midline (line)/(area)                  2 x Face #
    %                    +X     \ |      \ |   \  +Y   \ |   6  \ |  5 = left    'midz' = z-Midline (line)/(area)                  2 x Face #
    %                           2\|_______\|4   \       \|_______\|  6 = right   'midf' = Midpoint of face (point)/(area)          2 x Face #
    %                                            v +X        1                   'mids' = Surface btwn faces (surface)/(volume)    Face # to Face #
    %                                                                                     (feature loaded/constrained)/(optional feature loaded/constrained)
    %                                                                            OPTIONAL
    %                                                                            '##_@@' where # is the feature number and @ is the percent of domain to cover. Usually ## is the same feature numbers specified in the second column. The exception to this is 'edge'
    %                                                                            For example {'edge', 13, 12, '13_20'} as a mechanical support means that the middle 20% of the 13 edge should be fixed in the x and y directions
    %                                                                                        {'edge', 13, 12, '11_20'} as a mechanical support means that the left (node 1 side) 20% of the 13 edge should be fixed in the x and y directions
    % Mechanical Loads and Supports
    % scale.bc{phys = (1,2,3),1} = {'name', ##, [magx magy magz], optional} % define loads
    % scale.bc{phys = (1,2,3),2} = {'name', ##, dofs, optional} % define supports
    for i = 1:3:length(lcoorp)
        for k = 1:size(scale.bc{phys,1},1)
            if i == 1
                verify_BC(scale.bc{phys,1}{k,1},scale.bc{phys,1}{k,2},scale.bc{phys,1}{k,3},scale.bc{phys,1}{k,4},1,phys,scale.dim); % verify load BC are leagal
            end
            if isPoint(scale.bc{phys,1}{k,1},scale.bc{phys,1}{k,2},scale.bc{phys,1}{k,4},i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,scale.bc{phys,1}{k,3},F,scale.dim,phys);
            end
        end
        for k = 1:size(scale.bc{phys,2},1)
            if i == 1
                verify_BC(scale.bc{phys,2}{k,1},scale.bc{phys,2}{k,2},scale.bc{phys,2}{k,3},scale.bc{phys,2}{k,4},2,phys,scale.dim); % verify support BC are leagal
            end
            if isPoint(scale.bc{phys,2}{k,1},scale.bc{phys,2}{k,2},scale.bc{phys,2}{k,4},i,lcoorp,xyzc,dx,dy,dz)
                N = assign_supp(i,N,scale.bc{phys,2}{k,3},phys);
            end
        end
    end
    
    if scale.dim == 2 && phys == 1
        F(3:3:end) = [];
        U(3:3:end) = [];
        N(3:3:end) = [];
    end
end

% MESH INFORMATION
function [nnodes, dx, dy, dz, xyzc, lcoorp] = mesh_info_detailed(scale)
    nnodes = (scale.nelx+1)*(scale.nely+1)*(scale.nelz+1);

    dx = 1/scale.nelx;
    dy = 1/scale.nely;
    if scale.dim == 2
        dz = 1;
    elseif scale.dim == 3
        dz = 1/scale.nelz;
    end

    xyzc(1) = 0;
    xyzc(2) = scale.nelx*dx;
    xyzc(3) = 0;
    xyzc(4) = scale.nely*dy;
    xyzc(5) = 0;
    xyzc(6) = scale.nelz*dz;
    
    x = 0:dx:1;
    y = 0:dy:1;
    z = 0:dz:1;

    lcoorp = zeros(nnodes,3);
    r = 1;
    for k = 1:(scale.nelz+1)
        for j = 1:(scale.nely+1)
            for i = 1:(scale.nelx+1)
                lcoorp(r,1) = x(i);
                lcoorp(r,2) = y(j);
                lcoorp(r,3) = z(k);
                r = r + 1;
            end
        end
    end
    lcoorp = lcoorp'; lcoorp = lcoorp(:);
end

% ASSIGN NODAL LOADS
function LOADout = assign_load(i,LOADin,LOADout,dim,phys)
    if phys == 1
        for d = 1:dim
            LOADout(i+d-1) = LOADin(d);
        end
    elseif phys == 2
        LOADout((i-1)/3 + 1) = LOADin;
    end
end

% ASSIGN NODAL SUPPORTS
function DOFout = assign_supp(i,DOFout,loadDOFs,phys)
    if phys == 1
        if loadDOFs == 1
            DOFout(i) = 0;
        elseif loadDOFs == 2
            DOFout(i+1) = 0;
        elseif loadDOFs == 3
            DOFout(i+2) = 0;
        elseif loadDOFs == 12 || loadDOFs == 21
            DOFout(i) = 0; DOFout(i+1) = 0;
        elseif loadDOFs == 23 || loadDOFs == 32
            DOFout(i+1) = 0; DOFout(i+2) = 0;
        elseif loadDOFs == 13 || loadDOFs == 31
            DOFout(i) = 0; DOFout(i+2) = 0;
        elseif loadDOFs == 123 || loadDOFs == 132 || loadDOFs == 213 || loadDOFs == 231 || loadDOFs == 312 || loadDOFs == 321
            DOFout(i) = 0; DOFout(i+1) = 0; DOFout(i+2) = 0;
        end
    elseif phys == 2
        DOFout((i-1)/3 + 1) = 0;
    end
end

% VERIFY INPUT BOUNDARY CONDITIONS
function verify_BC(name,idnums,mag_dofs,optional_ids,bctype,phys,dim)
    if length(name) ~= 4 % Name
        error(['Error: ' name ' is not a valid boundary condition name'])
    end
    
    if length(num2str(idnums)) ~= 2 % ID numbers
        error(['Error: ' num2str(idnums) ' is not a valid boundary condition ##'])
    end
    
    if phys == 1 && bctype == 1 && length(mag_dofs) ~= 3 % Mechanical, Load magnitude
        error(['Error: [' num2str(mag_dofs) '] is not a valid mechanical load magnitude'])
    end
    
    if phys == 2 && bctype == 1 && length(mag_dofs) ~= 1 % Thermal, Load magnitude
        error(['Error: [' num2str(mag_dofs) '] is not a valid thermal load magnitude'])
    end
    
    if phys == 1 && dim == 2 && bctype == 2 &&... % Mechanical, 2D Support Dofs
       (mag_dofs ~= 1 && mag_dofs ~= 2 && mag_dofs ~= 12 && mag_dofs ~= 21)
        error(['Error: [' num2str(mag_dofs) '] is not a valid mechanical DOF constraint'])
    end
    
    if phys == 1 && dim == 3 && bctype == 2 && (mag_dofs ~= 1 && mag_dofs ~= 2 && mag_dofs ~= 3 && mag_dofs ~= 12 && mag_dofs ~= 21 && mag_dofs ~= 23 && mag_dofs ~= 32 && mag_dofs ~= 13 && mag_dofs ~= 31 &&...
       mag_dofs ~= 123 && mag_dofs ~= 132 && mag_dofs ~= 213 && mag_dofs ~= 231 && mag_dofs ~= 321 && mag_dofs ~= 312) % Mechanical, 3D Support Dofs
        error(['Error: [' num2str(mag_dofs) '] is not a valid mechanical DOF constraint'])
    end
    
    if phys == 2 && bctype == 2 &&  mag_dofs ~= 1 % Thermal, Support Dofs
        error(['Error: [' num2str(mag_dofs) '] is not a valid thermal DOF constraint'])
    end
    
    if strcmp(name,'midp') % Midpoint of domain
        if idnums ~= 00
            error(['Error: Incorrect ID specfified for ' name])
        end
    end
    
    if strcmp(name,'corn') % Corner
        if dim == 2
            if idnums ~= 11 && idnums ~= 22 && idnums ~= 33 && idnums ~= 44
                error(['Error: Incorrect ID specfified for ' name])
            end
        elseif dim == 3
            if idnums ~= 11 && idnums ~= 22 && idnums ~= 33 && idnums ~= 44 && idnums ~= 55 && idnums ~= 66 && idnums ~= 77 && idnums ~= 88
                error(['Error: Incorrect ID specfified for ' name])
            end
        end
    end
    
    if strcmp(name,'edge') || strcmp(name,'mide') % Edge, Midpoint of edge
        if dim == 2
            if idnums ~= 12 && idnums ~= 21 && idnums ~= 24 && idnums ~= 42 && idnums ~= 34 && idnums ~= 43 && idnums ~= 13 && idnums ~= 31
                error(['Error: Incorrect ID specfified for ' name])
            end
        elseif dim == 3
            if idnums ~= 12 && idnums ~= 21 && idnums ~= 24 && idnums ~= 42 && idnums ~= 34 && idnums ~= 43 && idnums ~= 13 && idnums ~= 31 &&...
               idnums ~= 56 && idnums ~= 65 && idnums ~= 68 && idnums ~= 86 && idnums ~= 78 && idnums ~= 87 && idnums ~= 57 && idnums ~= 75 &&...
               idnums ~= 15 && idnums ~= 51 && idnums ~= 26 && idnums ~= 62 && idnums ~= 48 && idnums ~= 84 && idnums ~= 37 && idnums ~= 73
                error(['Error: Incorrect ID specfified for ' name])
            end
        end
    end
    
    if strcmp(name,'face') || strcmp(name,'midx') || strcmp(name,'midy') || strcmp(name,'midz') || strcmp(name,'midf') % Face, x-Midline, y-Midline, z-Midline, Midpoint of face
        if dim == 2
            if idnums ~= 11
                error(['Error: Incorrect ID specfified for ' name])
            end
        elseif dim == 3
            if idnums ~= 11 && idnums ~= 22 && idnums ~= 33 && idnums ~= 44 && idnums ~= 55 && idnums ~= 66
                error(['Error: Incorrect ID specfified for ' name])
            end
        end
    end
    
    if strcmp(name,'mids') % Middle surface between two faces
        if dim == 2
            error(['Error: Incorrect ID specfified for ' name])
        elseif dim == 3
            if idnums ~= 12 && idnums ~= 21 && idnums ~= 34 && idnums ~= 43 && idnums ~= 56 && idnums ~= 65
                error(['Error: Incorrect ID specfified for ' name])
            end
        end
    end

    % If optional ids is not empty then verify input
    if ~isempty(optional_ids)
        % Midpoint of domain (add 0.5*percentage in the x, y, and z-directions)
        % Corner (adds percentage in the x, y, and z-directions)
        % Edge (adds percentage in the x, y, or z-direction)
        % Face (adds percentage in the x, y, or z-direction for 3D problems)
        % Midpoint of edge (adds 0.5*percentage in the x, y, or z-direction)
        % x-Midline (adds 0.5*percentage in the y and|or z direction)
        % y-Midline (adds 0.5*percentage in the x and|or z direction)
        % z-Midline (adds 0.5*percentage in the x and|or y direction)
        % Midpoint of face (add 0.5*percentage in the x and y, y and z, or x and z-directions)
        % Middle surface between two faces (adds 0.5*percentage|percentage in the x, y, or z-direction)
        
        if ~(isfloat(str2double(optional_ids(1:2))) && ~isempty(str2double(optional_ids(1:2))) && length(optional_ids(1:2)) == 2) % Check for correct format in optional area ID
            error(['Error: ' optional_ids(1:2) ' is not a valid optional area ID'])
        end

        area_ids = mod(floor(str2double(optional_ids(1:2)) ./ 10 .^ (floor(log10(str2double(optional_ids(1:2)))):-1:0)), 10); % Check for a valid optional area ID for the selected boundary condition
        idnums_vec = num2str(idnums); idnums_vec = [str2num(idnums_vec(1)), str2num(idnums_vec(2))];
        if any(area_ids ~= [idnums_vec(1), idnums_vec(1)]) && any(area_ids ~= [idnums_vec(2), idnums_vec(2)]) && any(area_ids ~= [idnums_vec(1), idnums_vec(2)]) && any(area_ids ~= [idnums_vec(2), idnums_vec(1)])
        % if idnums ~= str2double([num2str(area_ids(1)) num2str(area_ids(1))]) && idnums ~= str2double([num2str(area_ids(2)) num2str(area_ids(2))]) && idnums ~= str2double([num2str(area_ids(1)) num2str(area_ids(2))]) && idnums ~= str2double([num2str(area_ids(2)) num2str(area_ids(1))])
            error(['Error: The optional area ID ' optional_ids(1:2) ' is invalid for boundary condition ' num2str(idnums)])
        end
        
        if isnan(str2double(optional_ids(4:end))) || str2double(optional_ids(4:end)) <= 0 || str2double(optional_ids(4:end)) >= 100 % Check for a valid percentage of area to be added
            error(['Error: ' optional_ids(4:end) ' is not a valid optional area size'])
        end
        
        % Usually ## is the same feature numbers specified in the second column. The exception to this is 'edge'
        if strcmp(name,'midp') || strcmp(name,'corn') || strcmp(name,'face') || strcmp(name,'mide') || strcmp(name,'midx') || strcmp(name,'midy') || strcmp(name,'midz') || strcmp(name,'midf') || strcmp(name,'mids')
            if str2double(optional_ids(1:2)) ~= idnums
                error(['Error: Optional area ID should match feature ## for ' name])
            end
        end
        
        if strcmp(name,'edge')
            if dim == 2
                if str2double(optional_ids(1:2)) ~= 11 && str2double(optional_ids(1:2)) ~= 22 && str2double(optional_ids(1:2)) ~= 33 && str2double(optional_ids(1:2)) ~= 44
                    error(['Error: Incorrect optional area ID specfified for ' name])
                end
            elseif dim == 3
                if str2double(optional_ids(1:2)) ~= 11 && str2double(optional_ids(1:2)) ~= 22 && str2double(optional_ids(1:2)) ~= 33 && str2double(optional_ids(1:2)) ~= 44 && str2double(optional_ids(1:2)) ~= 55 && str2double(optional_ids(1:2)) ~= 66 && str2double(optional_ids(1:2)) ~= 77 && str2double(optional_ids(1:2)) ~= 88
                    error(['Error: Incorrect optional area ID specfified for ' name])
                end
            end
        end
    end
end