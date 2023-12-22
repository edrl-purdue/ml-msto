% FINITE ELEMENT ANALYSIS FOR ANY SCALE
function [U,Ue_int,Fe_int] = FEA(xPhys,scale,penal,mat,KE,F,U,N)
    %% ELEMENT CONNECTIVITY MATRICES
    edof = cell(1,3);
    if any(scale.phys == 1) % 1 = Mechanical
        edof{1} = zeros(scale.nele,scale.nen*scale.dim);
    end
    if any(scale.phys == 2) % 2 = Thermal
        edof{2} = zeros(scale.nele,scale.nen);
    end
    
    for i = 1:scale.nele
        for j = 1:scale.nen
            for k = 1:scale.dim
                if any(scale.phys == 1) % 1 = Mechanical
                    edof{1}(i,(j-1)*scale.dim + k) = scale.dim * scale.necon((i-1)*scale.nen + j) + k;
                end
            end
            if any(scale.phys == 2) % 2 = Thermal
                edof{2}(i,j) = scale.necon((i-1)*scale.nen + j) + 1;
            end
        end
    end
    
    %% ASSEMBLE STIFFNESS MATRICES
    [K, K_nele, ~] = Prepare_K(xPhys,scale,penal,mat,KE,[],[],edof,[]); % Case 1: Calling Prepare_K(xPhys,scale,penal,mat,KE,[],[],edof,[]) is SIMP stiffness matrix only
    
    %% SOLVE FEA
    if any(scale.phys == 2) % 2 = Thermal
        F{2}(N{2}==1) = F{2}(N{2}==1) - K{2,2}(N{2}==1,N{2}==0)*U{2}(N{2}==0);
        for LC = 1:length(U{2}(1,:))
            try
                R = chol(K{2,2}(N{2}(:,LC)==1,N{2}(:,LC)==1));
                U{2}(N{2}(:,LC)==1,LC) = R\(R'\F{2}(N{2}(:,LC)==1,LC)); % solve Cholesky thermal problem
            catch
                U{2}(N{2}(:,LC)==1,LC) = K{2,2}(N{2}(:,LC)==1,N{2}(:,LC)==1)\F{2}(N{2}(:,LC)==1,LC); % solve standard thermal problem
            end
        end
    end
    
    if all(scale.phys == [1 2]) % 1,2 = Thermomechanical
        Ft = zeros(size(F{1}));
        Ft(N{1}==1) = K{2,1}(N{1}==1,N{2}==1)*U{2}(N{2}==1); % solve thermomechanical problem
        F{1} = F{1} + Ft; % add thermally-induced forces to mechanical load
    end
    if any(scale.phys == 1) % 1 = Mechanical
        F{1}(N{1}==1) = F{1}(N{1}==1) - K{1,1}(N{1}==1,N{1}==0)*U{1}(N{1}==0);
        for LC = 1:length(U{1}(1,:))
            try
                R = chol(K{1,1}(N{1}(:,LC)==1,N{1}(:,LC)==1));
                U{1}(N{1}(:,LC)==1,LC) = R\(R'\F{1}(N{1}(:,LC)==1,LC)); % solve Cholesky mechanical problem
            catch
                U{1}(N{1}(:,LC)==1,LC) = K{1,1}(N{1}(:,LC)==1,N{1}(:,LC)==1)\F{1}(N{1}(:,LC)==1,LC); % solve standard mechanical problem
            end
        end
    end
    
    Ue_int = cell(1,3);
    Fe_int = cell(1,3);
    if any(scale.phys == 1) % 1 = Mechanical
        Ue_int{1} = zeros(scale.nele,scale.nen*scale.dim,length(U{1}(1,:)));
        Fe_int{1} = zeros(scale.nele,scale.nen*scale.dim,length(U{1}(1,:)));
        if length(F{1}(1,:)) == 1
            for i = 1:scale.nele
                Ue_int{1}(i,:) = U{1}(edof{1}(i,:)); % local elemental displacements
                Fe_int{1}(i,:) = reshape(K_nele{1,1}(i,:),[scale.nen*scale.dim,scale.nen*scale.dim])*Ue_int{1}(i,:)'; % solve for internal forces
            end
        end
    end
    if any(scale.phys == 2) % 2 = Thermal
        Ue_int{2} = zeros(scale.nele,scale.nen,length(U{2}(1,:)));
        Fe_int{2} = zeros(scale.nele,scale.nen,length(U{2}(1,:)));
        if length(F{2}(1,:)) == 1
            for i = 1:scale.nele
                Ue_int{2}(i,:) = U{2}(edof{2}(i,:)); % local elemental temperatures
                Fe_int{2}(i,:) = reshape(K_nele{2,2}(i,:),[scale.nen,scale.nen])*Ue_int{2}(i,:)'; % solve for internal fluxes
            end
        end
    end
    if any(scale.phys == 3) % 3 = Fluid
        Ue_int{3} = zeros(scale.nele,scale.nen*scale.dim,length(U{3}(1,:)));
        Fe_int{3} = zeros(scale.nele,scale.nen*scale.dim,length(U{3}(1,:)));
    end
end