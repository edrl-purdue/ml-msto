% GENERATE STIFFNESS TENSOR WITHOUT E, k, etc
function [C] = Prepare_C(nu,phys,dim) %without material properties (E, k, etc)
    if phys == 1 % Mechanical
        L1 = nu*1/((1+nu)*(1 - 2*nu));
        L2 = 1/(2*(1+nu));
        if dim == 2 % 2D
            C_L1 = zeros(3,3); C_L1(1:2,1:2) = L1;
            C_L2 = diag([2 2 1])*L2;
            C = C_L1 + C_L2; % plane strain assumption (usually go with this, better for 2d extruded designs)
%             C = (1/(1-nu*nu))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]; % plane stress assumption (better for 2d thin plate designs)
        elseif dim == 3 % 3D
            C_L1 = zeros(6,6); C_L1(1:3,1:3) = L1;
            C_L2 = diag([2 2 2 1 1 1])*L2;
            C = C_L1 + C_L2;
        end
    elseif phys == 2 % Thermal
        if dim == 2 % 2D
            C = eye(2);
        elseif dim == 3 % 3D
            C = eye(3);
        end
    elseif phys == 3 % Fluid
        if dim == 2 % 2D
            C = 0;
        elseif dim == 3 % 3D
            C = 0;
        end
    end
end