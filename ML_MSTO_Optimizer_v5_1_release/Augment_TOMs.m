% AUGMENT TOPOLOGY-OPTIMIZED MICROSTRUCTURES
function [X,Y_C] = Augment_TOMs(scale,X0,Y_C0)
    %% CREATE COORDINATE DATA FOR OPERATIONS: X0
    X_data = zeros(scale.nen,7); %[Fx/Ux_e, Fy/Uy_e, Fz/Uz_e, Q/T_e, x_e, y_e, z_e]
    if scale.dim == 2
        max_aug = 8;
        if any(scale.phys == 1)
            X_data(:,1) = X0((scale.spdim_nums(1,1)+1):scale.dim:scale.spdim_nums(2,1)); % Fx/Ux of X
            X_data(:,2) = X0((scale.spdim_nums(1,1)+2):scale.dim:scale.spdim_nums(2,1)); % Fy/Uy of X
            X_data(:,3) = 0;
        end
        if any(scale.phys == 2)
            X_data(:,4) = X0((scale.spdim_nums(1,2)+1):scale.spdim_nums(2,2)); % Q/T of X
        end
        X_data(:,5:7) = [-1 -1 0; % node 1
                          1 -1 0; % node 2
                          1  1 0; % node 4
                         -1  1 0]*scale.nelx/2; %node 3
    elseif scale.dim == 3
        max_aug = 48;
        if any(scale.phys == 1)
            X_data(:,1) = X0((scale.spdim_nums(1,1)+1):scale.dim:scale.spdim_nums(2,1)); % Fx/Ux of X
            X_data(:,2) = X0((scale.spdim_nums(1,1)+2):scale.dim:scale.spdim_nums(2,1)); % Fy/Uy of X
            X_data(:,3) = X0((scale.spdim_nums(1,1)+3):scale.dim:scale.spdim_nums(2,1)); % Fz/Uz of Z
        end
        if any(scale.phys == 2)
            X_data(:,4) = X0((scale.spdim_nums(1,2)+1):scale.spdim_nums(2,2)); % Q/T of X
        end
        X_data(:,5:7) = [-1 -1 -1; % node 1
                          1 -1 -1; % node 2
                          1  1 -1; % node 4
                         -1  1 -1; % node 3
                         -1 -1  1; % node 5
                          1 -1  1; % node 6
                          1  1  1; % node 8
                         -1  1  1]*scale.nelx/2; % node 7
    end
    
    %% INITIALIZE OUTPUTS
    X = zeros(max_aug,scale.spdim);
    X(:,1:min(scale.spdim_nums(1,:))) = X0(1:min(scale.spdim_nums(1,:)));
    Y_C = cell(max_aug,1); % Homogenized Tensor Output
    
    %% DEFINE 2D/3D OPERATIONS
    ops = cell(max_aug,1);
    if scale.dim == 2
        ops{2} =  {'Rz', 90};
        ops{3} =  {'Rz', 180};
        ops{4} =  {'Rz', 270};
        ops{5} =  {'Fx', []};
        ops{6} =  {'Rz', 90; 'Fx', []};
        ops{7} =  {'Fx', []; 'Rz', 90};
        ops{8} =  {'Rz', 180; 'Fx', []};
    elseif scale.dim == 3
        ops{2} =  {'Rx', 90};
        ops{3} =  {'Rx', 180};
        ops{4} =  {'Rx', 270};
        ops{5} =  {'Ry', 90};
        ops{6} =  {'Rx', 90; 'Ry', 90};
        ops{7} =  {'Ry', 90; 'Rx', 90};
        ops{8} =  {'Rx', 180; 'Ry', 90};
        ops{9} =  {'Ry', 90; 'Rx', 180};
        ops{10} = {'Rx', 270; 'Ry', 90};
        ops{11} = {'Ry', 90; 'Rx', 270};
        ops{12} = {'Ry', 180};
        ops{13} = {'Rx', 90; 'Ry', 180};
        ops{14} = {'Ry', 180; 'Rx', 90};
        ops{15} = {'Rx', 180; 'Ry', 180};
        ops{16} = {'Ry', 270};
        ops{17} = {'Rx', 90; 'Ry', 270};
        ops{18} = {'Ry', 270; 'Rx', 90};
        ops{19} = {'Rx', 270; 'Ry', 270};
        ops{20} = {'Ry', 270; 'Rx', 270};
        ops{21} = {'Rz', 90};
        ops{22} = {'Rx', 180; 'Rz', 90};
        ops{23} = {'Rz', 90; 'Rx', 180};
        ops{24} = {'Rx', 180; 'Ry', 180; 'Rz', 90};
        ops{25} = {'Fx', []};
        ops{26} = {'Rx', 90; 'Fx', []};
        ops{27} = {'Rx', 180; 'Fx', []};
        ops{28} = {'Rx', 270; 'Fx', []};
        ops{29} = {'Ry', 90; 'Fx', []};
        ops{30} = {'Fx', []; 'Ry', 90};
        ops{31} = {'Rx', 90; 'Ry', 90; 'Fx', []};
        ops{32} = {'Rx', 90; 'Fx', []; 'Ry', 90};
        ops{33} = {'Ry', 90; 'Rx', 90; 'Fx', []};
        ops{34} = {'Fx', []; 'Ry', 90; 'Rx', 90};
        ops{35} = {'Rx', 180; 'Ry', 90; 'Fx', []};
        ops{36} = {'Rx', 180; 'Fx', []; 'Ry', 90};
        ops{37} = {'Rx', 270; 'Ry', 90; 'Fx', []};
        ops{38} = {'Rx', 270; 'Fx', []; 'Ry', 90};
        ops{39} = {'Ry', 90; 'Rx', 270; 'Fx', []};
        ops{40} = {'Fx', []; 'Ry', 90; 'Rx', 270};
        ops{41} = {'Ry', 180; 'Fx', []};
        ops{42} = {'Rx', 90; 'Ry', 180; 'Fx', []};
        ops{43} = {'Ry', 180; 'Rx', 90; 'Fx', []};
        ops{44} = {'Rx', 180; 'Ry', 180; 'Fx', []};
        ops{45} = {'Rz', 90; 'Fx', []};
        ops{46} = {'Fx', []; 'Rz', 90};
        ops{47} = {'Rx', 180; 'Rz', 90; 'Fx', []};
        ops{48} = {'Rx', 180; 'Fx', []; 'Rz', 90};
    end
    
    %% PERFORM AUGMENTATIONS
    X(1,:) = X0(1:scale.spdim);
    Y_C{1} = Y_C0;
    for c = 2:max_aug
        X_data(:,:,c) = X_data(:,:,1);
        [X_data(:,:,c), A] = compute_A(X_data(:,:,c),ops{c});
        X = transform_X(X,X_data,c,scale.spdim_nums,scale);
        Y_C = transform_C(Y_C,A,c,scale);
    end
end

%% COMPUTE THE TRANSFORMATION MATRIX A FOR THE GIVEN OPERATIONS
function [X_data_c, A] = compute_A(X_data_c,ops)
    A = eye(3);
    for op = 1:size(ops,1)
        if strcmp(ops{op,1},'Rx')
            [X_data_c, A] = rot_x(X_data_c,ops{op,2},A);
        elseif strcmp(ops{op,1},'Ry')
            [X_data_c, A] = rot_y(X_data_c,ops{op,2},A);
        elseif strcmp(ops{op,1},'Rz')
            [X_data_c, A] = rot_z(X_data_c,ops{op,2},A);
        elseif strcmp(ops{op,1},'Fx')
            [X_data_c, A] = flip_x(X_data_c,A);
        elseif strcmp(ops{op,1},'Fy')
            [X_data_c, A] = flip_y(X_data_c,A);
        elseif strcmp(ops{op,1},'Fz')
            [X_data_c, A] = flip_z(X_data_c,A);
        end
    end
end

%% TRANSFORM COORDINATE AND BOUNDARY CONDITION DATA
function X = transform_X(X,X_data,c,spdim_nums,scale)
    for e1 = 1:scale.nen %original coord
        for e2 = 1:scale.nen %transformed coord
            if any(X_data(e1,5,1) == X_data(e2,5,c)) && any(X_data(e1,6,1) == X_data(e2,6,c)) && any(X_data(e1,7,1) == X_data(e2,7,c))
                if any(scale.phys == 1)
                    X(c,spdim_nums(1,1)+e1*scale.dim-(scale.dim-1)) = X_data(e2,1,c)'; % Fx/Ux of X
                    X(c,spdim_nums(1,1)+e1*scale.dim-(scale.dim-2)) = X_data(e2,2,c)'; % Fy/Uy of X
                    if scale.dim == 3; X(c,spdim_nums(1,1)+e1*scale.dim) = X_data(e2,3,c)'; end % Fz/Uz of X
                end
                if any(scale.phys == 2)
                    X(c,spdim_nums(1,2)+e1) = X_data(e2,4,c)'; % Q/T of X
                end
            end
        end
    end
    X(c,1:min(spdim_nums(spdim_nums>0))) = X(1,1:min(spdim_nums(spdim_nums>0)));
end

%% TRANSFORM HOMOGENIZED CONSTITUTIVE TENSORS
function Y_C = transform_C(Y_C,A,c,scale)
    if scale.dim == 2
        Abar = [A(1,1)*A(1,1), A(1,2)*A(1,2),             2*A(1,1)*A(1,2);
                A(2,1)*A(2,1), A(2,2)*A(2,2),             2*A(2,1)*A(2,2);
                A(1,1)*A(2,1), A(1,2)*A(2,2), A(1,1)*A(2,2)+A(1,2)*A(2,1)];
        F = eye(3); F(3,3) = 2;
    elseif scale.dim == 3
        Abar = [A(1,1)*A(1,1), A(1,2)*A(1,2), A(1,3)*A(1,3),               2*A(1,2)*A(1,3),               2*A(1,1)*A(1,3),               2*A(1,1)*A(1,2);
                A(2,1)*A(2,1), A(2,2)*A(2,2), A(2,3)*A(2,3),               2*A(2,2)*A(2,3),               2*A(2,1)*A(2,3),               2*A(2,1)*A(2,2);
                A(3,1)*A(3,1), A(3,2)*A(3,2), A(3,3)*A(3,3),               2*A(3,2)*A(3,3),               2*A(3,1)*A(3,3),               2*A(3,1)*A(3,2);
                A(2,1)*A(3,1), A(2,2)*A(3,2), A(2,3)*A(3,3), A(2,2)*A(3,3) + A(2,3)*A(3,2), A(2,1)*A(3,3) + A(2,3)*A(3,1), A(2,1)*A(3,2) + A(2,2)*A(3,1);
                A(1,1)*A(3,1), A(1,2)*A(3,2), A(1,3)*A(3,3), A(1,2)*A(3,3) + A(1,3)*A(3,2), A(1,1)*A(3,3) + A(1,3)*A(3,1), A(1,1)*A(3,2) + A(1,2)*A(3,1);
                A(1,1)*A(2,1), A(1,2)*A(2,2), A(1,3)*A(2,3), A(1,2)*A(2,3) + A(1,3)*A(2,2), A(1,1)*A(2,3) + A(1,3)*A(2,1), A(1,1)*A(2,2) + A(1,2)*A(2,1)];
        F = eye(6); F(4,4) = 2; F(5,5) = 2; F(6,6) = 2;
    end
    Ma = F*(Abar/F);
    Y_C{c} = cell(3,3); % # of rows/cols is the number of unique physics modes
    if any(scale.phys == 1)
        Y_C{c}{1,1} = Ma'*Y_C{1}{1,1}*Ma; % Mechanical
    end
    if any(scale.phys == 2)
        Y_C{c}{2,2} = A(1:scale.dim,1:scale.dim)'*Y_C{1}{2,2}*A(1:scale.dim,1:scale.dim); % Thermal
    end
    if all(scale.phys == [1 2])
        Y_C{c}{2,1} = Ma*Y_C{1}{2,1}; % Thermomechanical
    end
end

%% ROTATE ABOUT THE X, Y, OR Z AXIS
function [cd, R] = rot_x(cd0,a,A)
    R = [1 0 0; 0 cosd(a) -sind(a); 0 sind(a) cosd(a)]; %Rx
    RR = [R zeros(3,4); 0 0 0 1 0 0 0; zeros(3,4) R]; % expanded rotation matrix for operation on the 1,2,3,5,6,7 columns of X_data (the columns with data that is orientation-relevant)
    
    if size(cd0,2) == 3 % if coordinates only
        cd = cd0*R';
    else % operate on the expanded cd
        cd = cd0*RR';
    end
    R = R*A;
end
function [cd, R] = rot_y(cd0,a,A)
    R = [cosd(a) 0 sind(a); 0 1 0; -sind(a) 0 cosd(a)]; %Ry
    RR = [R zeros(3,4); 0 0 0 1 0 0 0; zeros(3,4) R]; % expanded rotation matrix for operation on the 1,2,3,5,6,7 columns of X_data (the columns with data that is orientation-relevant)
    
    if size(cd0,2) == 3 % if coordinates only
        cd = cd0*R';
    else % operate on the expanded cd
        cd = cd0*RR';
    end
    R = R*A;
end
function [cd, R] = rot_z(cd0,a,A)
    R = [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1]; %Rz
    RR = [R zeros(3,4); 0 0 0 1 0 0 0; zeros(3,4) R]; % expanded rotation matrix for operation on the 1,2,3,5,6,7 columns of X_data (the columns with data that is orientation-relevant)
    
    if size(cd0,2) == 3 % if coordinates only
        cd = cd0*R';
    else % operate on the expanded cd
        cd = cd0*RR';
    end
    R = R*A;
end

%% FLIP COORDINATES ALONG THE X, Y, OR Z DIRECTION
function [cd, R] = flip_x(cd0,A)
    R = eye(3);
    R(1,1) = -R(1,1);
    RR = [R zeros(3,4); 0 0 0 1 0 0 0; zeros(3,4) R]; % expanded rotation matrix for operation on the 1,2,3,5,6,7 columns of X_data (the columns with data that is orientation-relevant)
    
    if size(cd0,2) == 3 % if coordinates only
        cd = cd0*R';
    else
        cd = cd0*RR';
    end
    R = R*A;
end
function [cd, R] = flip_y(cd0,A)
    R = eye(3);
    R(2,2) = -R(2,2);
    RR = [R zeros(3,4); 0 0 0 1 0 0 0; zeros(3,4) R]; % expanded rotation matrix for operation on the 1,2,3,5,6,7 columns of X_data (the columns with data that is orientation-relevant)
    
    if size(cd0,2) == 3 % if coordinates only
        cd = cd0*R';
    else
        cd = cd0*RR';
    end
    R = R*A;
end
function [cd, R] = flip_z(cd0,A)
    R = eye(3);
    R(3,3) = -R(3,3);
    RR = [R zeros(3,4); 0 0 0 1 0 0 0; zeros(3,4) R]; % expanded rotation matrix for operation on the 1,2,3,5,6,7 columns of X_data (the columns with data that is orientation-relevant)
    
    if size(cd0,2) == 3 % if coordinates only
        cd = cd0*R';
    else
        cd = cd0*RR';
    end
    R = R*A;
end