% OBJECTIVE FUNCTION AND ANALYTICAL DERIVATIVE
function [c, dc, c_vec] = Obj_Fns(xTheta,scale,U,penal,mat,KE,KE_Hom,dKE_Hom,weight)
    % Case 1: Calling [c, dc, c_vec] = Obj_Fns(xMicro,scale,U_Micro,penal,mat,KE,[],[],weight_e) is c and dc for SIMP stiffness matrix
    % Case 2: Calling [c, ~, c_vec] = Obj_Fns([],macro,U_Macro,[],[],[],KE_Hom,[],weight_e) is c for Homogenized stiffness matrix
    % Case 3: Calling [~, dc, ~] = Obj_Fns([],macro,U_Macro,[],[],[],[],dKE_Hom,weight_e) is dc for Homogenized stiffness matrix
    
    %% OBJECTIVE FUNCTION WEIGHTS
    if scale.phys == 1
        w1 = 1.0; %mechanical
        w2 = 0.0; %thermal
    elseif scale.phys == 2
        w1 = 0.0; %mechanical
        w2 = 1.0; %thermal
    elseif all(scale.phys == [1 2])
        w1 = weight; %mechanical
        w2 = (1-weight); %thermal
    end
    
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
    
    %% COMPUTE OBJECTIVE FUNCTION AND ANALYTICAL SENSITIVITIES
    if ~isempty(KE) % Case 1
        c = zeros(scale.nele,2);
        dc = zeros(scale.nele,3);
%         xTheta(:) = (scale.H*xTheta(:))./scale.Hs;
        for i = 1:scale.nele
            uKu_m = 0;
            uKu_t = 0;
            uKu_mt = 0;
            if any(scale.phys == 1)
                for k = 1:scale.nen*scale.dim
                    for h = 1:scale.nen*scale.dim
                        uKu_m = uKu_m + U{1}(edof{1}(i,k))*KE{1,1}((k-1)*scale.nen*scale.dim+h)*U{1}(edof{1}(i,h)); % multiplied by [simp E] or [derivative of simp E]
                    end
                end
            end
            if any(scale.phys == 2)
                for k = 1:scale.nen
                    for h = 1:scale.nen
                        uKu_t = uKu_t + U{2}(edof{2}(i,k))*KE{2,2}((k-1)*scale.nen+h)*U{2}(edof{2}(i,h)); % multiplied by [simp k] or [derivative of simp k]
                    end
                end
            end
            if all(scale.phys == [1 2])
                for k = 1:scale.nen*scale.dim
                    for h = 1:scale.nen
                        uKu_mt = uKu_mt + U{1}(edof{1}(i,k))*KE{2,1}((k-1)*scale.nen+h)*U{2}(edof{2}(i,h)); % multiplied by [derivative of simp alpha - (simp alpha/simp k)*derivative of simp k]
                    end
                end
            end
            
            % objective function
            c(i,1) = w1*(mat.Emin + (xTheta(i)^penal)*(mat.Emax-mat.Emin))*uKu_m;
            c(i,2) = w2*(mat.kmin + (xTheta(i)^penal)*(mat.kmax-mat.kmin))*uKu_t;
            if all(scale.phys == [1 2])
                dc(i,1) = -w1*penal*(xTheta(i)^(penal-1))*(mat.Emax-mat.Emin)*uKu_m;
                dc(i,2) = -w2*penal*(xTheta(i)^(penal-1))*(mat.kmax-mat.kmin)*uKu_t;
                dc(i,3) = w1*(penal*(xTheta(i)^(penal-1))*(mat.Amax-mat.Amin) - ((mat.Amin + (xTheta(i)^penal)*(mat.Amax-mat.Amin))/(mat.kmin + (xTheta(i)^penal)*(mat.kmax-mat.kmin)))*(penal*(xTheta(i)^(penal-1))*(mat.kmax-mat.kmin)))*uKu_mt;
            else
                dc(i,1) = -w1*penal*(xTheta(i)^(penal-1))*(mat.Emax-mat.Emin)*uKu_m; % no dKmt/dTheta derivative term
                dc(i,2) = -w2*penal*(xTheta(i)^(penal-1))*(mat.kmax-mat.kmin)*uKu_t;
            end
        end
        c_vec = sum(c,2);
        c = sum(c(:));
        dc = sum(dc,2);
    elseif ~isempty(KE_Hom) % Case 2
        c = zeros(numel(KE_Hom),2);
        for i = 1:numel(KE_Hom)
            uKu_m = 0;
            uKu_t = 0;
            uKu_mt = 0;
            if any(scale.phys == 1)
                for k = 1:scale.nen*scale.dim
                    for h = 1:scale.nen*scale.dim
                        uKu_m = uKu_m + U{1}(edof{1}(i,k))*KE_Hom{i}{1,1}((k-1)*scale.nen*scale.dim+h)*U{1}(edof{1}(i,h)); % multiplied by [simp E] or [derivative of simp E]
                    end
                end
            end
            if any(scale.phys == 2)
                for k = 1:scale.nen
                    for h = 1:scale.nen
                        uKu_t = uKu_t + U{2}(edof{2}(i,k))*KE_Hom{i}{2,2}((k-1)*scale.nen+h)*U{2}(edof{2}(i,h)); % multiplied by [simp k] or [derivative of simp k]
                    end
                end
            end
            if all(scale.phys == [1 2])
                for k = 1:scale.nen*scale.dim
                    for h = 1:scale.nen
                        uKu_mt = uKu_mt + U{1}(edof{1}(i,k))*KE_Hom{i}{2,1}((k-1)*scale.nen+h)*U{2}(edof{2}(i,h)); % multiplied by [derivative of simp alpha - (simp alpha/simp k)*derivative of simp k]
                    end
                end
            end
            
            % objective function
            c(i,1) = w1*uKu_m;
            c(i,2) = w2*uKu_t;
        end
        c_vec = sum(c,2);
        c = sum(c(:));
        dc = [];
    elseif ~isempty(dKE_Hom) % Case 3 NOT EXTENDED TO THERMOMECHANICAL YET
        dc = zeros(numel(dKE_Hom),3);
        for i = 1:numel(dKE_Hom)
            udKu_m = 0;
            udKu_t = 0;
            if any(scale.phys == 1)
                for k = 1:scale.nen*scale.dim
                    for h = 1:scale.nen*scale.dim
                        udKu_m = udKu_m + U{1}(edof{1}(i,k))*dKE_Hom{i}{1,1}((k-1)*scale.nen*scale.dim+h)*U{1}(edof{1}(i,h)); % multiplied by [simp E] or [derivative of simp E]
                    end
                end
            end
            if any(scale.phys == 2)
                for k = 1:scale.nen
                    for h = 1:scale.nen
                        udKu_t = udKu_t + U{2}(edof{2}(i,k))*dKE_Hom{i}{2,2}((k-1)*scale.nen+h)*U{2}(edof{2}(i,h)); % multiplied by [simp E] or [derivative of simp E]
                    end
                end
            end
            if all(scale.phys == [1 2])
                for k = 1:scale.nen*scale.dim
                    for h = 1:scale.nen
                         % NOT EXTENDED TO THERMOMECHANICAL YET
%                         udKu_mt = udKu_mt + U{1}(edof{1}(i,k))*KE_Hom{i}{2,1}((k-1)*scale.nen+h)*U{2}(edof{2}(i,h)); % multiplied by [derivative of simp alpha - (simp alpha/simp k)*derivative of simp k]
                    end
                end
            end
            
            if all(scale.phys == [1 2]) % NOT EXTENDED TO THERMOMECHANICAL YET
                dc(i,1) = w1*udKu_m;
                dc(i,2) = w2*udKu_t;
%                 dc(i,3) = ???
%                 dc(i,3) = w1*(penal*(xTheta(i)^(penal-1))*(mat.Amax-mat.Amin) - ((mat.Amin + (xTheta(i)^penal)*(mat.Amax-mat.Amin))/(mat.kmin + (xTheta(i)^penal)*(mat.kmax-mat.kmin)))*(penal*(xTheta(i)^(penal-1))*(mat.kmax-mat.kmin)))*uKu_mt;
            else
                dc(i,1) = w1*udKu_m; % no udKu_mt
                dc(i,2) = w2*udKu_t;
            end
        end
        c_vec = [];
        c = [];
        dc = sum(dc,2);
    end
end