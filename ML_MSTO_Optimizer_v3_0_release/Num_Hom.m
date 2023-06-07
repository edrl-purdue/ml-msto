% NUMERICAL HOMOGENIZATION
function CH = Num_Hom(xPhys,scale,mat,KE,FE)
    % This code was adapted from numerical homogenization code from "How to
    % determine composite material properties using numerical
    % homogenization." by Erik Andreassen and Casper Schousboe Andreasen,
    % published in Computational Materials Science, 2014.
    % https://doi.org/10.1016/j.commatsci.2013.09.006

    %% ELEMENT CONNECTIVITY MATRICES
    edof_per = cell(1,3);
    if any(scale.phys == 1) % 1 = Mechanical
        edof_per{1} = zeros(scale.nele,scale.nen*scale.dim);
    end
    if any(scale.phys == 2) % 2 = Thermal
        edof_per{2} = zeros(scale.nele,scale.nen);
    end
    
    for i = 1:scale.nele
        for j = 1:scale.nen
            for k = 1:scale.dim
                if any(scale.phys == 1) % 1 = Mechanical
                    edof_per{1}(i,(j-1)*scale.dim + k) = scale.dim * scale.necon_per((i-1)*scale.nen + j) + k;
                end
            end
            if any(scale.phys == 2) % 2 = Thermal
                edof_per{2}(i,j) = scale.necon_per((i-1)*scale.nen + j) + 1;
            end
        end
    end
    num_LC(:,:,1) = [3, 1; 1, 2]; num_LC(:,:,2) = [6, 1; 1, 3]; % number of load cases for FE. Rows/Cols are physics modes. Pages are dimension, 2D/3D
     
    %% ASSEMBLE STIFFNESS MATRICES
    if scale.dim == 3
        KE{1,1} = KE{1,1}./scale.nelx; % make KE mesh dependent to be compatible with FE
        KE{2,2} = KE{2,2}./scale.nelx; % make KE mesh dependent to be compatible with FE
    end
    [K, ~, F] = Prepare_K(xPhys,scale,scale.penal,mat,KE,[],FE,edof_per,num_LC); % Case 2: Calling Prepare_K(xPhys,scale,scale.penal,mat,KE,[],FE,edof_per,num_LC) is SIMP stiffness and loading matrices
    
    x_E = mat.Emin + (xPhys.^scale.penal)*(mat.Emax-mat.Emin); % assumed the same penalization parameter is used for all 3 material properties
    x_k = mat.kmin + (xPhys.^scale.penal)*(mat.kmax-mat.kmin);
    x_A = mat.Amin + (xPhys.^scale.penal)*(mat.Amax-mat.Amin);

    %% COMPUTE FIELD VARIABLES RESULTING FROM PRESCRIBED UNIT STRAINS
    chi = cell(3,3);
    dof = cell(3,3);
    if any(scale.phys == 1) % 1 = Mechanical
        chi{1,1} = zeros(scale.dim*scale.nele,num_LC(1,1,scale.dim-1)); %ndof number of rows
        dof{1,1} = (scale.dim + 1):scale.dim*scale.nele;
        try
            R = chol(K{1,1}(dof{1,1},dof{1,1}));
            chi{1,1}(dof{1,1},:) = R\(R'\F{1,1}(dof{1,1},:)); %full calculations % Solve Cholesky (remember to constrain one node)
        catch
            chi{1,1}(dof{1,1},:) = K{1,1}(dof{1,1},dof{1,1})\F{1,1}(dof{1,1},:); %full calculations % Solve standard (remember to constrain one node)
        end
    end
    if any(scale.phys == 2) % 2 = Thermal
        chi{2,2} = zeros(scale.nele,num_LC(2,2,scale.dim-1));
        dof{2,2} = 2:scale.nele;
        try
            R = chol(K{2,2}(dof{2,2},dof{2,2}));
            chi{2,2}(dof{2,2},:) = R\(R'\F{2,2}(dof{2,2},:)); %full calculations % Solve Cholesky (remember to constrain one node)
        catch
            chi{2,2}(dof{2,2},:) = K{2,2}(dof{2,2},dof{2,2})\F{2,2}(dof{2,2},:); %full calculations % Solve standard (remember to constrain one node)
        end
    end
    if all(scale.phys == [1 2]) % 1,2 = Thermomechanical
        chi{2,1} = zeros(scale.dim*scale.nele,num_LC(2,1,scale.dim-1));
        dof{2,1} = (scale.dim + 1):scale.dim*scale.nele;
        try
            R = chol(K{1,1}(dof{2,1},dof{2,1}));
            chi{2,1}(dof{2,1},1) = R\(R'\F{2,1}(dof{2,1},1)); %full calculations % Solve Cholesky (remember to constrain one node)
        catch
            chi{2,1}(dof{2,1},1) = K{1,1}(dof{2,1},dof{2,1})\F{2,1}(dof{2,1},1); %full calculations % Solve standard (remember to constrain one node)
        end
    end
    
    %% COMPUTE PRESCRIBED FIELD VARIABLES
    chi0 = cell(3,3); chi0e = cell(3,3);
    dofe = cell(3,3);
    if any(scale.phys == 1) % 1 = Mechanical
        chi0{1,1} = zeros(scale.nele, scale.dim*scale.nen, num_LC(1,1,scale.dim-1)); % The displacement vectors corresponding to the unit strain cases
        chi0e{1,1} = zeros(scale.dim*scale.nen, num_LC(1,1,scale.dim-1)); % The element displacements for the three unit strains
        if scale.dim == 2; dofe{1,1} = [3, 5:8]; elseif scale.dim == 3; dofe{1,1} = [4, 6:8, 10:12, 14:24]; end
        KE{1,1} = reshape(KE{1,1}, [scale.dim*scale.nen,scale.dim*scale.nen]); % Here the exact ratio does not matter, because
        FE{1,1} = reshape(FE{1,1}, [scale.dim*scale.nen,num_LC(1,1,scale.dim-1)]); % it is reflected in the load vector
        try
            R = chol(KE{1,1}(dofe{1,1},dofe{1,1}));
            chi0e{1,1}(dofe{1,1},:) = R\(R'\FE{1,1}(dofe{1,1},:)); % solve Cholesky
        catch
            chi0e{1,1}(dofe{1,1},:) = KE{1,1}(dofe{1,1},dofe{1,1})\FE{1,1}(dofe{1,1},:); % solve standard
        end
        for lc = 1:num_LC(1,1,scale.dim-1)
            chi0{1,1}(:,:,lc) = kron(chi0e{1,1}(:,lc)', ones(scale.nele,1)); % epsilon0_lclc = (0, ..., 1, ..., 0)
        end
    end
    if any(scale.phys == 2) % 2 = Thermal
        chi0{2,2} = zeros(scale.nele, scale.nen, num_LC(2,2,scale.dim-1)); % The displacement vectors corresponding to the unit strain cases
        chi0e{2,2} = zeros(scale.nen, num_LC(2,2,scale.dim-1)); % The element displacements for the three unit strains
        dofe{2,2} = 2:scale.nen;
        KE{2,2} = reshape(KE{2,2}, [scale.nen,scale.nen]); % Here the exact ratio does not matter, because
        FE{2,2} = reshape(FE{2,2}, [scale.nen,num_LC(2,2,scale.dim-1)]); % it is reflected in the load vector
        try
            R = chol(KE{2,2}(dofe{2,2},dofe{2,2}));
            chi0e{2,2}(dofe{2,2},:) = R\(R'\FE{2,2}(dofe{2,2},:)); % solve Cholesky
        catch
            chi0e{2,2}(dofe{2,2},:) = KE{2,2}(dofe{2,2},dofe{2,2})\FE{2,2}(dofe{2,2},:); % solve standard
        end
        for lc = 1:num_LC(2,2,scale.dim-1)
            chi0{2,2}(:,:,lc) = kron(chi0e{2,2}(:,lc)', ones(scale.nele,1)); % epsilon0_lclc = (0, ..., 1, ..., 0)
        end
    end
    if all(scale.phys == [1 2]) % 1,2 = Thermomechanical
        chi0{2,1} = zeros(scale.nele, scale.dim*scale.nen, num_LC(2,1,scale.dim-1)); % The displacement vectors corresponding to the unit strain cases
        chi0e{2,1} = zeros(scale.dim*scale.nen, num_LC(2,1,scale.dim-1)); % The element displacements for the three unit strains
        if scale.dim == 2; dofe{2,1} = [3, 5:8]; elseif scale.dim == 3; dofe{2,1} = [4, 6:8, 10:12, 14:24]; end
        FE{2,1} = reshape(FE{2,1}, [scale.dim*scale.nen,num_LC(2,1,scale.dim-1)]); % it is reflected in the load vector
        try
            R = chol(KE{1,1}(dofe{2,1},dofe{2,1}));
            chi0e{2,1}(dofe{2,1},1) = R\(R'\FE{2,1}(dofe{2,1},1)); % solve Cholesky
        catch
            chi0e{2,1}(dofe{2,1},1) = KE{1,1}(dofe{2,1},dofe{2,1})\FE{2,1}(dofe{2,1},1); % solve standard
        end
        for lc = 1:num_LC(2,1,scale.dim-1)
            chi0{2,1}(:,:,lc) = x_A(:)*chi0e{2,1}(:,lc)';
        end
    end
    
    %% HOMOGENIZATION
    CH = cell(3,3); % # of rows/cols is the number of unique physics modes
    cellVolume = 1*1*1; % lx, ly, lz = 1: Unit cell length in x-direction, y-direction, and z-direction.
    if any(scale.phys == 1) % 1 = Mechanical
        for i = 1:num_LC(1,1,scale.dim-1)
            for j = 1:num_LC(1,1,scale.dim-1)
                sum_KE = ((chi0{1,1}(:,:,i) - chi{1,1}(edof_per{1}+(i-1)*max(dof{1,1})))*KE{1,1}).*(chi0{1,1}(:,:,j) - chi{1,1}(edof_per{1}+(j-1)*max(dof{1,1})));
                sum_KE = reshape(sum(sum_KE,2), size(xPhys));
                CH{1,1}(i,j) = 1/cellVolume*sum(sum(sum(x_E.*sum_KE))); % Homogenized elasticity tensor
            end
        end
    end
    if any(scale.phys == 2) % 2 = Thermal
        for i = 1:num_LC(2,2,scale.dim-1)
            for j = 1:num_LC(2,2,scale.dim-1)
                sum_KE = ((chi0{2,2}(:,:,i) - chi{2,2}(edof_per{2}+(i-1)*max(dof{2,2})))*KE{2,2}).*(chi0{2,2}(:,:,j) - chi{2,2}(edof_per{2}+(j-1)*max(dof{2,2})));
                sum_KE = reshape(sum(sum_KE,2), size(xPhys));
                CH{2,2}(i,j) = 1/cellVolume*sum(sum(sum(x_k.*sum_KE))); % Homogenized conductivity tensor
            end
        end
    end
    if all(scale.phys == [1 2]) % 1,2 = Thermomechanical
        for i = 1:num_LC(1,1,scale.dim-1)
            sum_KE = ((chi0{2,1}(:,:) - chi{2,1}(edof_per{1}))*KE{1,1}).*(chi0{1,1}(:,:,i) - chi{1,1}(edof_per{1}+(i-1)*max(dof{1,1})));
            sum_KE = reshape(sum(sum_KE,2), size(xPhys));
            CH{2,1}(i,1) = 1/cellVolume*sum(sum(sum(x_E.*sum_KE))); % Homogenized thermal stress vector
        end
%         CH{2,1} = CH{1,1}/CH{2,1}'; % Homogenized thermal strain vector
    end
    CH{1,1} = real(0.5*(CH{1,1}+CH{1,1}'));
    CH{2,2} = real(0.5*(CH{2,2}+CH{2,2}'));
    CH{2,1} = real(CH{2,1});
end