% GENERATE ELEMENT STIFFNESS MATRIX AND HOMOGENIZATION LOAD VECTOR - MULTIPHYSICS
function [KE, FE, GE] = Prepare_KE(C,phys,scale)
    % input C is stiffness/conductivity tensor for KE and FE generation (thermal, mechanical)
    % input C is the thermal stress tensor for FE generation (thermomechanical)
    % no input required for coupling tensors GE
    
    %nels is the number of side elements for a microscale that is homogenized
    nels = scale.nelx; % assumes that nelx = nely = nelz
    
    if scale.dim == 2
%         a = (1/nels)/2; b = (1/nels)/2;
        a = (scale.DDx/scale.nelx)/2; b = (scale.DDy/scale.nely)/2;
        X = [-a a a -a]; % a and b matter for FE. They do not affect 2D KE, but affect 3D KE. This is adjusted for below (line 137).
        Y = [-b -b b b];
        if phys == 1 % Mechanical
            ndofn = 2;
            epsilon = diag([1 1 1]);
            alpha = zeros(size(epsilon,1),ndofn,scale.dim);
            alpha(1,1,1) = 1; alpha(3,2,1) = 1;
            alpha(2,2,2) = 1; alpha(3,1,2) = 1;
        elseif phys == 2 % Thermal
            ndofn = 1;
            epsilon = diag([1 1]);
            alpha = zeros(size(epsilon,1),ndofn,scale.dim);
            alpha(1,1,1) = 1;
            alpha(2,1,2) = 1;
        elseif all(phys == [1 2]) % Thermomechanical
            ndofn = 2;
            epsilon = [1; 1; 0];
            alpha = zeros(size(epsilon,1),ndofn,scale.dim);
            alpha(1,1,1) = 1; alpha(3,2,1) = 1;
            alpha(2,2,2) = 1; alpha(3,1,2) = 1;
        elseif phys == 3 % Fluid
            ndofn = 3;
        end
    elseif scale.dim == 3
%         a = (1/nels)/2; b = (1/nels)/2; c = (1/nels)/2;
        a = (scale.DDx/scale.nelx)/2; b = (scale.DDy/scale.nely)/2; c = (scale.DDz/scale.nelz)/2;
        X = [-a a a -a -a a a -a];
        Y = [-b -b b b -b -b b b];
        Z = [-c -c -c -c c c c c];
        if phys == 1 % Mechanical
            ndofn = 3;
            epsilon = diag([1 1 1 1 1 1]);
            alpha = zeros(size(epsilon,1),ndofn,scale.dim);
            alpha(1,1,1) = 1; alpha(4,2,1) = 1; alpha(6,3,1) = 1;
            alpha(2,2,2) = 1; alpha(4,1,2) = 1; alpha(5,3,2) = 1;
            alpha(3,3,3) = 1; alpha(5,2,3) = 1; alpha(6,1,3) = 1;
        elseif phys == 2 % Thermal
            ndofn = 1;
            epsilon = diag([1 1 1]);
            alpha = zeros(size(epsilon,1),ndofn,scale.dim);
            alpha(1,1,1) = 1;
            alpha(2,1,2) = 1;
            alpha(3,1,3) = 1;
        elseif all(phys == [1 2]) % Thermomechanical
            ndofn = 3;
            epsilon = [1; 1; 1; 0; 0; 0];
            alpha = zeros(size(epsilon,1),ndofn,scale.dim);
            alpha(1,1,1) = 1; alpha(4,2,1) = 1; alpha(6,3,1) = 1;
            alpha(2,2,2) = 1; alpha(4,1,2) = 1; alpha(5,3,2) = 1;
            alpha(3,3,3) = 1; alpha(5,2,3) = 1; alpha(6,1,3) = 1;
        elseif phys == 3 % Fluid
            ndofn = 4;
        end
    end

    GP = [-1/sqrt(3), 1/sqrt(3)];
    W = [1 1];

    KE = zeros(scale.nen*ndofn,scale.nen*ndofn);
    FE = zeros(scale.nen*ndofn,size(epsilon,2));
    GE = zeros(scale.nen*ndofn,scale.nen);
    beta = zeros(size(epsilon,1),ndofn);
    for ii = 1:2
        for jj = 1:2
            if scale.dim == 2
                [N, dN] = N_dN_calc([GP(ii), GP(jj)]);
                J(1,1) = dot(dN(:,1),X);
                J(1,2) = dot(dN(:,1),Y);
                J(2,1) = dot(dN(:,2),X);
                J(2,2) = dot(dN(:,2),Y);
                invJ = inv(J);
                weight = W(ii)*W(jj)*det(J)*scale.DDz;
                B = zeros(size(epsilon,1),scale.nen*ndofn);

                for ll = 1:scale.dim
                    for i = 1:size(epsilon,1)
                        for j = 1:ndofn
                            beta(i,j) = invJ(ll,ll)*alpha(i,j,ll);
                        end
                        for j = 1:(scale.nen*ndofn)
                            B(i,j) = B(i,j) + beta(i,mod(j-1,ndofn)+1)*dN(floor((j-1)/ndofn)+1,ll);
                        end
                    end
                end

                if all(phys == 1) || all(phys == 2)
                    KE = KE + weight*B'*C*B; % this 2D KE calculation is already mesh-independent
                    FE = FE + weight*B'*C*epsilon; % FE is inherently mesh-dependent
                    GE = [];
                elseif all(phys == [1 2])
                    KE = [];
                    FE = FE + weight*B'*C*epsilon;
                    GE = GE + weight*B'*epsilon*N;
                end
            elseif scale.dim == 3
                for kk = 1:2
                    [N, dN] = N_dN_calc([GP(ii), GP(jj), GP(kk)]);
                    J(1,1) = dot(dN(:,1),X);
                    J(1,2) = dot(dN(:,1),Y);
                    J(1,3) = dot(dN(:,1),Z);
                    J(2,1) = dot(dN(:,2),X);
                    J(2,2) = dot(dN(:,2),Y);
                    J(2,3) = dot(dN(:,2),Z);
                    J(3,1) = dot(dN(:,3),X);
                    J(3,2) = dot(dN(:,3),Y);
                    J(3,3) = dot(dN(:,3),Z);
                    invJ = inv(J);
                    weight = W(ii)*W(jj)*W(kk)*det(J);
                    B = zeros(size(epsilon,1),scale.nen*ndofn);

                    for ll = 1:scale.dim
                        for i = 1:size(epsilon,1)
                            for j = 1:ndofn
                                beta(i,j) = invJ(ll,ll)*alpha(i,j,ll);
                            end
                        end
                        for i = 1:size(epsilon,1)
                            for j = 1:(scale.nen*ndofn)
                                B(i,j) = B(i,j) + beta(i,mod(j-1,ndofn)+1)*dN(floor((j-1)/ndofn)+1,ll);
                            end
                        end
                    end
                    
                    if all(phys == 1) || all(phys == 2)
                        KE = KE + weight*B'*C*B*nels; % nels factor added to make 3D KE mesh-independent (for fea, but not for homogenization)
                        FE = FE + weight*B'*C*epsilon; % FE is inherently mesh-dependent
                        GE = [];
                    elseif all(phys == [1 2])
                        KE = [];
                        FE = FE + weight*B'*C*epsilon;
                        GE = GE + weight*B'*epsilon*N;
                    end
                end
            end
        end
    end
    KE = 0.5*(KE+KE'); KE = KE(:);
    FE = FE(:);
    GE = GE(:);
end