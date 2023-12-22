% PREPARE MICRO BOUNDARY CONDITIONS
function [F,U,N] = Prepare_BCMicro(scale,PTOM_e,phys)
    %% MESH INFORMATION
    [nnodes, dx, dy, dz, xyzc, lcoorp] = mesh_info(scale);
    
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
    if scale.bc_type == 1 % Internal Force BC
        Fe_int = PTOM_e(scale.spdim_nums(1,phys)+1:scale.spdim_nums(2,phys));
        Ue_int = zeros(size(Fe_int));
    elseif scale.bc_type == 2 % Internal Displacement BC
        Ue_int = PTOM_e(scale.spdim_nums(1,phys)+1:scale.spdim_nums(2,phys));
        Fe_int = zeros(size(Ue_int));
    end
    
    %% BOUNDARY CONDITIONS
%                                                                            NAME    FEATURE                      ##
%                 DOMAIN NODE AND FACE NUMBERINGS                            'midp' = Midpoint of domain(point)    00
%      2D Nodes               3D Nodes                  3D Faces             'corn' = Corner (point)               2 x Corner #
%    3 ________ 4         5 ________ 7             ________                  'edge' = Edge (line)                  Corner # to Corner #
%     |        |  ^+Y      |\       |\    ^+Z     |\   2   |\    1 = bottom  'face' = Face (surface)               2 x Face #
%     |  1 =   |  |        |6\ _____|_\ 8 |       | \ _____|_\   2 = top     'mide' = Midpoint of edge (point)     Corner # to Corner #
%     |  face  |  |        |  |     |  |  |       |3 | 5   |  |  3 = front   'midx' = x-Midline (line)             2 x Face #
%    1|________|2 +---->  1|__|_____|3 |  +---->  |__|_____| 4|  4 = back    'midy' = y-Midline (line)             2 x Face #
%                    +X     \ |      \ |   \  +Y   \ |   6  \ |  5 = left    'midz' = z-Midline (line)             2 x Face #
%                           2\|_______\|4   \       \|_______\|  6 = right   'midf' = Midpoint of face (point)     2 x Face #
%                                            v +X        1                   'mids' = Surface btwn faces (surface) Face # to Face #
%   Mechanical Loads and Supports
%   macro.bc{phys = (1,2,3),1} = {'name', ##, [magx magy magz]} % define loads
%   macro.bc{phys = (1,2,3),2} = {'name', ##, dofs} % define supports
    for i = 1:3:length(lcoorp)
        % Calculate Natural Coordinates
        if scale.dim == 2
            coord = [lcoorp(i)*2 - 1, lcoorp(i+1)*2 - 1];
        elseif scale.dim == 3
            coord = [lcoorp(i)*2 - 1, lcoorp(i+1)*2 - 1, lcoorp(i+2)*2 - 1];
        end
        
        %% NODAL LOADS
        if any(scale.load_bc == 1) % Corner Loads
            if isPoint('corn',11,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('corn',22,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('corn',33,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('corn',44,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('corn',55,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('corn',66,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('corn',77,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('corn',88,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
        end
        
        if any(scale.load_bc == 2) % Midpoint of edge Loads
            if isPoint('mide',12,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('mide',24,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('mide',13,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('mide',34,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('mide',56,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('mide',68,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('mide',57,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('mide',78,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('mide',15,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('mide',26,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('mide',37,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('mide',48,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
        end
        
        if any(scale.load_bc == 3) % Edge Loads
            if isPoint('edge',12,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('edge',24,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('edge',13,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('edge',34,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('edge',56,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('edge',68,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('edge',57,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('edge',78,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('edge',15,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('edge',26,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('edge',37,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('edge',48,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
        end
        
        if any(scale.load_bc == 4) % Face Loads
            if isPoint('face',11,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('face',22,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('face',33,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('face',44,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('face',55,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('face',66,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
        end
        
        if any(scale.load_bc == 5) % Midpoint of face Loads
            if isPoint('midf',11,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('midf',22,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('midf',33,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('midf',44,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('midf',55,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
            if isPoint('midf',66,'',i,lcoorp,xyzc,dx,dy,dz)
                F = assign_load(i,coord,F,Fe_int,scale.dim,phys);
            end
        end
        
        %% NODAL SUPPORTS
        if any(scale.supp_bc == 1) % Corner Supports
            if isPoint('corn',11,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('corn',22,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('corn',33,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('corn',44,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('corn',55,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('corn',66,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('corn',77,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('corn',88,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
        end

        if any(scale.supp_bc == 2) % Midpoint of edge Supports
            if isPoint('mide',12,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('mide',24,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('mide',13,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('mide',34,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('mide',56,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('mide',68,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('mide',57,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('mide',78,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('mide',15,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('mide',26,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('mide',37,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('mide',48,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
        end
        
        if any(scale.supp_bc == 3) % Edge Supports
            if isPoint('edge',12,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('edge',24,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('edge',13,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('edge',34,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('edge',56,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('edge',68,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('edge',57,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('edge',78,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('edge',15,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('edge',26,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('edge',37,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('edge',48,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
        end
        
        if any(scale.supp_bc == 4) % Face Supports
            if isPoint('face',11,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('face',22,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('face',33,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('face',44,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('face',55,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('face',66,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
        end
        
        if any(scale.supp_bc == 5) % Midpoint of face Supports
            if isPoint('midf',11,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('midf',22,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('midf',33,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('midf',44,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('midf',55,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
            if isPoint('midf',66,'',i,lcoorp,xyzc,dx,dy,dz)
                [U, N] = assign_supp(i,coord,U,Ue_int,N,scale.dim,phys);
            end
        end

        % Observations:
        % Using only load_bc with scale.supp_bc = [], should not work because the stiffness matrix should be singular det(K_m) = 0, however it does work 
        % All of the scale.supp_bc options limit the deformation of the subproblem for force driven fea
    end
    
    if scale.dim == 2 && phys == 1
        F(3:3:end) = [];
        U(3:3:end) = [];
        N(3:3:end) = [];
    end
end

% MESH INFORMATION
function [nnodes, dx, dy, dz, xyzc, lcoorp] = mesh_info(scale)
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
function F = assign_load(i,coord,F,Fe_int,dim,phys)
    if phys == 1
        for d = 1:dim
            F(i+d-1) = N_interp(coord,Fe_int(d:dim:end));
        end
    elseif phys == 2
        F((i-1)/3 + 1) = N_interp(coord,Fe_int);
    end
end

% ASSIGN NODAL SUPPORTS
function [U, N] = assign_supp(i,coord,U,Ue_int,N,dim,phys)
    if phys == 1
        for d = 1:dim
            U(i+d-1) = N_interp(coord,Ue_int(d:dim:end));
            N(i+d-1) = 0;
        end
    elseif phys == 2
        U((i-1)/3 + 1) = N_interp(coord,Ue_int);
        N((i-1)/3 + 1) = 0;
    end
end

% SHAPE FUNCTION INTERPOLATION
function Xi = N_interp(coord,efield)
    Xi = dot(N_dN_calc(coord),efield);
end