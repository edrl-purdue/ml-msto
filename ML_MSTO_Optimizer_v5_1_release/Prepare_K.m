% ASSEMBLE STIFFNESS MATRIX - MULTIPHYSICS
function [K, K_nele, F] = Prepare_K(xPhys,scale,penal,mat,KE,KE_Hom,FE,edof,num_LC)
    % Case 1: Calling Prepare_K(xPhys,scale,penal,mat,KE,[],[],edof,[]) is SIMP stiffness matrix only
    % Case 2: Calling Prepare_K(xPhys,scale,scale.penal,mat,KE,[],FE,edof_per,num_LC) is SIMP stiffness and loading matrices
    % Case 3: Calling Prepare_K([],scale,[],[],[],KE_Hom,[],edof,[]) is Homogenized stiffness matrix only
    
    K = cell(3,3);
    K_nele = cell(3,3);
    K_i = cell(3,3);
    edof_k_i = cell(3,3); edof_k_j = cell(3,3);
    if ~isempty(FE) % if FE is not empty
        F = cell(3,3);
        F_nele = cell(3,3);
        F_i = cell(3,3);
        edof_f_i = cell(3,3); edof_f_j = cell(3,3); % input edof is edof_per for FE case
    else
        F = [];
        F_nele = [];
    end
    
    %% INITIALIZE STIFFNESS MATRICES
    if any(scale.phys == 1) % 1 = Mechanical
        K_nele{1,1} = zeros(scale.nele, (scale.nen*scale.dim)*(scale.nen*scale.dim));
        edof_k_i{1,1} = zeros(scale.nele*(scale.nen*scale.dim)*(scale.nen*scale.dim),1); edof_k_j{1,1} = zeros(scale.nele*(scale.nen*scale.dim)*(scale.nen*scale.dim),1);
        K_i{1,1} = zeros(scale.nele*(scale.nen*scale.dim)*(scale.nen*scale.dim), 1);
        if ~isempty(FE) % if FE is not empty
            F_nele{1,1} = zeros(scale.nele, (scale.nen*scale.dim)*num_LC(1,1,scale.dim-1));
            edof_f_i{1,1} = zeros(scale.nele*(scale.nen*scale.dim)*num_LC(1,1,scale.dim-1),1); edof_f_j{1,1} = zeros(scale.nele*(scale.nen*scale.dim)*num_LC(1,1,scale.dim-1),1);
            F_i{1,1} = zeros(scale.nele*(scale.nen*scale.dim)*num_LC(1,1,scale.dim-1), 1);
        end
    end
    if any(scale.phys == 2) % 2 = Thermal
        K_nele{2,2} = zeros(scale.nele, scale.nen*scale.nen);
        edof_k_i{2,2} = zeros(scale.nele*scale.nen*scale.nen, 1);   edof_k_j{2,2} = zeros(scale.nele*scale.nen*scale.nen, 1);
        K_i{2,2} = zeros(scale.nele*scale.nen*scale.nen, 1);
        if ~isempty(FE) % if FE is not empty
            F_nele{2,2} = zeros(scale.nele, scale.nen*num_LC(2,2,scale.dim-1));
            edof_f_i{2,2} = zeros(scale.nele*scale.nen*num_LC(2,2,scale.dim-1),1); edof_f_j{2,2} = zeros(scale.nele*scale.nen*num_LC(2,2,scale.dim-1),1);
            F_i{2,2} = zeros(scale.nele*scale.nen*num_LC(2,2,scale.dim-1), 1);
        end
    end
    if all(scale.phys == [1 2]) % 1,2 = Thermomechanical
        K_nele{2,1} = zeros(scale.nele, (scale.nen*scale.dim)*scale.nen);
        edof_k_i{2,1} = zeros(scale.nele*(scale.nen*scale.dim)*scale.nen, 1); edof_k_j{2,1} = zeros(scale.nele*(scale.nen*scale.dim)*scale.nen, 1);
        K_i{2,1} = zeros(scale.nele*(scale.nen*scale.dim)*scale.nen, 1);
        if ~isempty(FE) % if FE is not empty
            F_nele{2,1} = zeros(scale.nele, (scale.nen*scale.dim)*num_LC(2,1,scale.dim-1));
            edof_f_i{2,1} = zeros(scale.nele*(scale.nen*scale.dim)*num_LC(2,1,scale.dim-1),1); edof_f_j{2,1} = zeros(scale.nele*(scale.nen*scale.dim)*num_LC(2,1,scale.dim-1),1);
            F_i{2,1} = zeros(scale.nele*(scale.nen*scale.dim)*num_LC(2,1,scale.dim-1), 1);
        end
    end
    
    %% INTERPOLATE MATERIAL PROPERTIES
    for i = 1:scale.nele
        if any(scale.phys == 1) % 1 = Mechanical
            if ~isempty(KE)
                dens = mat.Emin + (xPhys(i)^penal)*(mat.Emax-mat.Emin);
                K_nele{1,1}(i,:) = KE{1,1}*dens;
                if ~isempty(FE) % if FE is not empty
                    F_nele{1,1}(i,:) = FE{1,1}*dens;
                end
            elseif ~isempty(KE_Hom)
                K_nele{1,1}(i,:) = KE_Hom{i}{1,1}(:);
            end
        end
        if any(scale.phys == 2) % 2 = Thermal
            if ~isempty(KE)
                dens = mat.kmin + (xPhys(i)^penal)*(mat.kmax-mat.kmin);
                K_nele{2,2}(i,:) = KE{2,2}*dens;
                if ~isempty(FE) % if FE is not empty
                    F_nele{2,2}(i,:) = FE{2,2}*dens;
                end
            elseif ~isempty(KE_Hom)
                K_nele{2,2}(i,:) = KE_Hom{i}{2,2}(:);
            end
        end
        if all(scale.phys == [1 2]) % 1,2 = Thermomechanical
            if ~isempty(KE)
                dens = mat.Amin + (xPhys(i)^penal)*(mat.Amax-mat.Amin);
                K_nele{2,1}(i,:) = KE{2,1}*dens;
                if ~isempty(FE) % if FE is not empty
                    F_nele{2,1}(i,:) = FE{2,1}*dens;
                end
            elseif ~isempty(KE_Hom)
                K_nele{2,1}(i,:) = KE_Hom{i}{2,1}(:);
            end
        end
    end
    
    %% ASSEMBLE STIFFNESS MATRICES
    for i = 1:scale.nele
        if any(scale.phys == 1) % 1 = Mechanical
            L1 = 1; L2 = 1;
            for j = 1:(scale.nen*scale.dim)
                for k = 1:(scale.nen*scale.dim)
                    edof_k_i{1,1}((i-1)*(scale.nen*scale.dim)*(scale.nen*scale.dim) + L1) = edof{1}(i,k);
                    edof_k_j{1,1}((i-1)*(scale.nen*scale.dim)*(scale.nen*scale.dim) + L1) = edof{1}(i,j);
                    L1 = L1 + 1;
                    if ~isempty(FE) % if FE is not empty
                        if j <= num_LC(1,1,scale.dim-1)
                            edof_f_i{1,1}((i-1)*(scale.nen*scale.dim)*num_LC(1,1,scale.dim-1) + L2) = edof{1}(i,k);
                            edof_f_j{1,1}((i-1)*(scale.nen*scale.dim)*num_LC(1,1,scale.dim-1) + L2) = j;
                            L2 = L2 + 1;
                        end
                    end
                end
            end
            K_i{1,1}(((i-1)*(scale.nen*scale.dim)*(scale.nen*scale.dim):(i*(scale.nen*scale.dim)*(scale.nen*scale.dim) - 1)) + 1) = K_nele{1,1}(i,:);
            if ~isempty(FE) % if FE is not empty
                F_i{1,1}(((i-1)*(scale.nen*scale.dim)*num_LC(1,1,scale.dim-1):(i*(scale.nen*scale.dim)*num_LC(1,1,scale.dim-1) - 1)) + 1) = F_nele{1,1}(i,:);
            end
        end
        if any(scale.phys == 2) % 2 = Thermal
            L1 = 1; L2 = 1;
            for j = 1:scale.nen
                for k = 1:scale.nen
                    edof_k_i{2,2}((i-1)*scale.nen*scale.nen + L1) = edof{2}(i,k);
                    edof_k_j{2,2}((i-1)*scale.nen*scale.nen + L1) = edof{2}(i,j);
                    L1 = L1 + 1;
                    if ~isempty(FE) % if FE is not empty
                        if j <= num_LC(2,2,scale.dim-1)
                            edof_f_i{2,2}((i-1)*scale.nen*num_LC(2,2,scale.dim-1) + L2) = edof{2}(i,k);
                            edof_f_j{2,2}((i-1)*scale.nen*num_LC(2,2,scale.dim-1) + L2) = j;
                            L2 = L2 + 1;
                        end
                    end
                end
            end
            K_i{2,2}(((i-1)*scale.nen*scale.nen:(i*scale.nen*scale.nen - 1)) + 1) = K_nele{2,2}(i,:);
            if ~isempty(FE) % if FE is not empty
                F_i{2,2}(((i-1)*scale.nen*num_LC(2,2,scale.dim-1):(i*scale.nen*num_LC(2,2,scale.dim-1) - 1)) + 1) = F_nele{2,2}(i,:);
            end
        end
        if all(scale.phys == [1 2]) % 1,2 = Thermomechanical
            L1 = 1;
            for j = 1:scale.nen
                for k = 1:(scale.nen*scale.dim)
                    edof_k_i{2,1}((i-1)*(scale.nen*scale.dim)*scale.nen + L1) = edof{1}(i,k);
                    edof_k_j{2,1}((i-1)*(scale.nen*scale.dim)*scale.nen + L1) = edof{2}(i,j);
                    L1 = L1 + 1;
                end
            end
            K_i{2,1}(((i-1)*(scale.nen*scale.dim)*scale.nen:(i*(scale.nen*scale.dim)*scale.nen - 1)) + 1) = K_nele{2,1}(i,:);
            if ~isempty(FE) % if FE is not empty
                L2 = 1;
                for j = 1:num_LC(2,1,scale.dim-1)
                    for k = 1:(scale.nen*scale.dim)
                        edof_f_i{2,1}((i-1)*(scale.nen*scale.dim)*num_LC(2,1,scale.dim-1) + L2) = edof{1}(i,k);
                        edof_f_j{2,1}((i-1)*(scale.nen*scale.dim)*num_LC(2,1,scale.dim-1) + L2) = j;
                        L2 = L2 + 1;
                    end
                end
                F_i{2,1}(((i-1)*(scale.nen*scale.dim)*num_LC(2,1,scale.dim-1):(i*(scale.nen*scale.dim)*num_LC(2,1,scale.dim-1) - 1)) + 1) = F_nele{2,1}(i,:);
            end
        end
    end
    
    if any(scale.phys == 1) % 1 = Mechanical
        if ~isempty(FE) % if FE is not empty
            K{1,1} = sparse(edof_k_i{1,1},edof_k_j{1,1},K_i{1,1},scale.nele*scale.dim,scale.nele*scale.dim); K{1,1} = 0.5*(K{1,1}+K{1,1}'); % perodic sizing
            F{1,1} = sparse(edof_f_i{1,1},edof_f_j{1,1},F_i{1,1},scale.nele*scale.dim,num_LC(1,1,scale.dim-1)); %three load cases for the three strain cases
        else
            K{1,1} = sparse(edof_k_i{1,1},edof_k_j{1,1},K_i{1,1},(scale.nelx+1)*(scale.nely+1)*(scale.nelz+1)*scale.dim,(scale.nelx+1)*(scale.nely+1)*(scale.nelz+1)*scale.dim); K{1,1} = 0.5*(K{1,1}+K{1,1}'); % nonperiodic sizing
        end
    end
    if any(scale.phys == 2) % 2 = Thermal
        if ~isempty(FE) % if FE is not empty
            K{2,2} = sparse(edof_k_i{2,2},edof_k_j{2,2},K_i{2,2},scale.nele,scale.nele); K{2,2} = 0.5*(K{2,2}+K{2,2}'); % periodic sizing
            F{2,2} = sparse(edof_f_i{2,2},edof_f_j{2,2},F_i{2,2},scale.nele,num_LC(2,2,scale.dim-1)); %two load cases
        else
            K{2,2} = sparse(edof_k_i{2,2},edof_k_j{2,2},K_i{2,2},(scale.nelx+1)*(scale.nely+1)*(scale.nelz+1),(scale.nelx+1)*(scale.nely+1)*(scale.nelz+1)); K{2,2} = 0.5*(K{2,2}+K{2,2}'); % nonperiodic sizing
        end
    end
    if all(scale.phys == [1 2]) % 1,2 = Thermomechanical
        if ~isempty(FE) % if FE is not empty
            K{2,1} = sparse(edof_k_i{2,1},edof_k_j{2,1},K_i{2,1},scale.nele*scale.dim,scale.nele); % periodic sizing
            F{2,1} = sparse(edof_f_i{2,1},edof_f_j{2,1},F_i{2,1},scale.nele*scale.dim,num_LC(2,1,scale.dim-1)); %one load cases
        else
            K{2,1} = sparse(edof_k_i{2,1},edof_k_j{2,1},K_i{2,1},(scale.nelx+1)*(scale.nely+1)*(scale.nelz+1)*scale.dim,(scale.nelx+1)*(scale.nely+1)*(scale.nelz+1)); % nonperiodic sizing
        end
    end
end