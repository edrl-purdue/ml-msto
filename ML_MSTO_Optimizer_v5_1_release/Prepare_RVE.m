% NUMERICAL HOMOGENIZATION - MULTIPHYSICS
function [line_data, point_data, RVE, node_conn, elem_conn, norm_scale, renums, CH_rot, ops] = Prepare_RVE(input_line_data, CH_unrot, scale)
    
    coord_tol = 1e-3;

    %% Create Point Data
    input_point_data = zeros(4, 1); ii = 1;
    for L = 1:size(input_line_data,1)
        for r = 1:2
            if all(input_line_data(L,r) ~= input_point_data(1,:))
                input_point_data(1,ii) = input_line_data(L,r);
                input_point_data(2:4,ii) = input_line_data(L, 3*r:3*r+2)';
                ii = ii + 1;
            end
        end
    end

    %% Shift Input Coordinates
    centerx = mean(input_point_data(2,any(unique(input_line_data(:,1:2)) == input_point_data(1,:),1)));
    centery = mean(input_point_data(3,any(unique(input_line_data(:,1:2)) == input_point_data(1,:),1)));
    centerz = mean(input_point_data(4,any(unique(input_line_data(:,1:2)) == input_point_data(1,:),1)));
    input_line_data(:,3) = input_line_data(:,3) - centerx; input_line_data(:,6) = input_line_data(:,6) - centerx;
    input_line_data(:,4) = input_line_data(:,4) - centery; input_line_data(:,7) = input_line_data(:,7) - centery;
    input_line_data(:,5) = input_line_data(:,5) - centerz; input_line_data(:,8) = input_line_data(:,8) - centerz;
    input_point_data(2,:) = input_point_data(2,:) - centerx;
    input_point_data(3,:) = input_point_data(3,:) - centery;
    input_point_data(4,:) = input_point_data(4,:) - centerz;

    %% Analyze Beam Elements
    for L = 1:size(input_line_data,1)
        input_line_data(L,9) = sqrt(sum((input_line_data(L,3:5) - input_line_data(L,6:8)).^2)); % length
        input_line_data(L,10) = 0; % initialize thickness
        input_line_data(L,11) = atan2d(input_line_data(L,7)-input_line_data(L,4),input_line_data(L,6)-input_line_data(L,3)); % element angle
        % input_line_data(L,12) = 0; % initialize alpha
    end

    %% Rotate RVE and CH_target
    if size(input_line_data,1) == 3 % Triangle
        [~, longest_ind] = max(input_line_data(:,9));
    elseif size(input_line_data,1) == 4 % Quadrilateral
        vec1 = input_line_data(1,6:8) - input_line_data(1,3:5);
        vec2 = input_line_data(2,6:8) - input_line_data(2,3:5);
        vec3 = input_line_data(3,6:8) - input_line_data(3,3:5);
        vec4 = input_line_data(4,6:8) - input_line_data(4,3:5);
        LR_angle_diff = 180 - acosd(dot(vec3,vec1)/(norm(vec3)*norm(vec1)));
        TB_angle_diff = 180 - acosd(dot(vec4,vec2)/(norm(vec4)*norm(vec2)));
        if LR_angle_diff == TB_angle_diff
            [~, longest_ind] = max(input_line_data(:,9));
        elseif LR_angle_diff < TB_angle_diff
            [~, longest_ind] = max(input_line_data(:,9).*[1 0 1 0]');
        elseif LR_angle_diff > TB_angle_diff
            [~, longest_ind] = max(input_line_data(:,9).*[0 1 0 1]');
        end
        
    end
    
    forward_vec = input_line_data(longest_ind, 6:8) - input_line_data(longest_ind, 3:5);
    forward_angle = atan2d(forward_vec(2),forward_vec(1));
    Rz = [cosd(-forward_angle), -sind(-forward_angle), 0; sind(-forward_angle), cosd(-forward_angle), 0; 0, 0, 1];
    
    rot_line_data = input_line_data;
    rot_line_data(:,3:5) = transpose(Rz*rot_line_data(:,3:5)');
    rot_line_data(:,6:8) = transpose(Rz*rot_line_data(:,6:8)');
    rot_point_data = input_point_data;
    rot_point_data(2:4,:) = Rz*input_point_data(2:4,:);

    % check to see if the forward rotation makes the 3rd point above (y+) the longest line
    % [~, off_ind] = max(all(input_line_data(longest_ind, 1:2)' ~= input_point_data(1,:),1));
    if rot_point_data(3,all(input_line_data(longest_ind, 1:2)' ~= input_point_data(1,:),1)) > rot_point_data(3,input_line_data(longest_ind, 1)' == input_point_data(1,:))
        ops = cell(1,2); ops{1} = 'Rz'; ops{2} = -forward_angle;
        RVE.origin = rot_line_data(longest_ind,1); % record origin (bottom left node) for rigid body constraint and RVE equilibrium
    else % use backward rotation
        backward_vec = forward_vec*-1;
        backward_angle = atan2d(backward_vec(2),backward_vec(1));
        Rz = [cosd(-backward_angle), -sind(-backward_angle), 0; sind(-backward_angle), cosd(-backward_angle), 0; 0, 0, 1];

        rot_line_data = input_line_data;
        rot_line_data(:,3:5) = transpose(Rz*rot_line_data(:,3:5)');
        rot_line_data(:,6:8) = transpose(Rz*rot_line_data(:,6:8)');
        rot_point_data = input_point_data;
        rot_point_data(2:4,:) = Rz*input_point_data(2:4,:);
        ops = cell(1,2); ops{1} = 'Rz'; ops{2} = -backward_angle;
        RVE.origin = rot_line_data(longest_ind,2); % record origin (bottom left node) for rigid body constraint and RVE equilibrium
    end

    CH_rot = Transform_CH({CH_unrot},scale,ops); CH_rot = CH_rot{1,1};

    %% Adjust Quadrilateral Geometry
    left_lean = zeros(size(rot_line_data,1),1);
    if size(input_line_data,1) == 4 && LR_angle_diff ~= TB_angle_diff
        % Make Quadrilateral have Parallel Sides via Rotation
        opp_inds = [3, 4, 1, 2]; % opposite side indexes
        centerx = mean(rot_line_data(opp_inds(longest_ind),[3,6]));
        centery = mean(rot_line_data(opp_inds(longest_ind),[4,7]));
        centerz = mean(rot_line_data(opp_inds(longest_ind),[5,8]));
        rot_angle = -min(LR_angle_diff,TB_angle_diff);
        Rz = [cosd(rot_angle), -sind(rot_angle), 0; sind(rot_angle), cosd(rot_angle), 0; 0, 0, 1];
        rot_line_data0 = rot_line_data;
        rot_line_data(rot_line_data(:,1) == rot_line_data(opp_inds(longest_ind),1),3:5) = transpose(Rz*(rot_line_data(rot_line_data(:,1) == rot_line_data(opp_inds(longest_ind),1),3:5)-[centerx,centery,centerz])') + [centerx,centery,centerz];
        rot_line_data(rot_line_data(:,1) == rot_line_data(opp_inds(longest_ind),2),3:5) = transpose(Rz*(rot_line_data(rot_line_data(:,1) == rot_line_data(opp_inds(longest_ind),2),3:5)-[centerx,centery,centerz])') + [centerx,centery,centerz];
        rot_line_data(rot_line_data(:,2) == rot_line_data(opp_inds(longest_ind),1),6:8) = transpose(Rz*(rot_line_data(rot_line_data(:,2) == rot_line_data(opp_inds(longest_ind),1),6:8)-[centerx,centery,centerz])') + [centerx,centery,centerz];
        rot_line_data(rot_line_data(:,2) == rot_line_data(opp_inds(longest_ind),2),6:8) = transpose(Rz*(rot_line_data(rot_line_data(:,2) == rot_line_data(opp_inds(longest_ind),2),6:8)-[centerx,centery,centerz])') + [centerx,centery,centerz];
        if ~(abs(rot_line_data(opp_inds(longest_ind),4) - rot_line_data(opp_inds(longest_ind),7)) < coord_tol)
            rot_line_data = rot_line_data0;
            Rz = [cosd(-rot_angle), -sind(-rot_angle), 0; sind(-rot_angle), cosd(-rot_angle), 0; 0, 0, 1];
            rot_line_data(rot_line_data(:,1) == rot_line_data(opp_inds(longest_ind),1),3:5) = transpose(Rz*(rot_line_data(rot_line_data(:,1) == rot_line_data(opp_inds(longest_ind),1),3:5)-[centerx,centery,centerz])') + [centerx,centery,centerz];
            rot_line_data(rot_line_data(:,1) == rot_line_data(opp_inds(longest_ind),2),3:5) = transpose(Rz*(rot_line_data(rot_line_data(:,1) == rot_line_data(opp_inds(longest_ind),2),3:5)-[centerx,centery,centerz])') + [centerx,centery,centerz];
            rot_line_data(rot_line_data(:,2) == rot_line_data(opp_inds(longest_ind),1),6:8) = transpose(Rz*(rot_line_data(rot_line_data(:,2) == rot_line_data(opp_inds(longest_ind),1),6:8)-[centerx,centery,centerz])') + [centerx,centery,centerz];
            rot_line_data(rot_line_data(:,2) == rot_line_data(opp_inds(longest_ind),2),6:8) = transpose(Rz*(rot_line_data(rot_line_data(:,2) == rot_line_data(opp_inds(longest_ind),2),6:8)-[centerx,centery,centerz])') + [centerx,centery,centerz];
        end

        % figure(23); clf; display_lines(rot_line_data, '-r'); display_lines(rot_line_data(longest_ind,:), '-g'); plot(0,0,'*b'); plot(centerx,centery,'or')
        
        % Rebuild Point Data
        rot_point_data = zeros(4, 1); ii = 1;
        for L = 1:size(rot_line_data,1)
            for r = 1:2
                if all(rot_line_data(L,r) ~= rot_point_data(1,:))
                    rot_point_data(1,ii) = rot_line_data(L,r);
                    rot_point_data(2:4,ii) = rot_line_data(L, 3*r:3*r+2)';
                    ii = ii + 1;
                end
            end
        end
        
        % Identify leaning vertical walls
        for L = 1:size(rot_line_data,1)
            if L ~= longest_ind && L ~= opp_inds(longest_ind)
                vec_data = sortrows([[rot_line_data(L,3:4), 1]; [rot_line_data(L,6:7), 2]],2);
                if vec_data(2,1) <= vec_data(1,1) % is left leaning
                    left_lean(L) = 1;
                end % otherwise is right leaning
            end
        end

        % Flip double right leaning quads
        if all(left_lean == 0)
            rot_line_data(:,3) = rot_line_data(:,3)*-1;
            rot_line_data(:,6) = rot_line_data(:,6)*-1;
            rot_point_data(2,:) = rot_point_data(2,:)*-1;
        end
        % figure(24); clf; display_lines(rot_line_data, '-r'); display_lines(rot_line_data(longest_ind,:), '-g'); plot(0,0,'*b');

        % align nodes if slightly off
        rot_point_data(3,abs(min(rot_point_data(3,:)) - rot_point_data(3,:)) < coord_tol) = mean(rot_point_data(3,abs(min(rot_point_data(3,:)) - rot_point_data(3,:)) < coord_tol)); % bottom edge nodes are consistent
        rot_point_data(3,abs(max(rot_point_data(3,:)) - rot_point_data(3,:)) < coord_tol) = mean(rot_point_data(3,abs(max(rot_point_data(3,:)) - rot_point_data(3,:)) < coord_tol)); % top edge nodes are consistent
        if sum(abs(min(rot_point_data(2,:)) - rot_point_data(2,:)) < coord_tol) == 2
            rot_point_data(2,abs(min(rot_point_data(2,:)) - rot_point_data(2,:)) < coord_tol) = mean(rot_point_data(2,abs(min(rot_point_data(2,:)) - rot_point_data(2,:)) < coord_tol)); % left edge nodes are consistent
        end
        if sum(abs(max(rot_point_data(2,:)) - rot_point_data(2,:)) < coord_tol) == 2
            rot_point_data(2,abs(max(rot_point_data(2,:)) - rot_point_data(2,:)) < coord_tol) = mean(rot_point_data(2,abs(max(rot_point_data(2,:)) - rot_point_data(2,:)) < coord_tol)); % right edge nodes are consistent
        end
        
        for L = 1:size(rot_line_data,1)
            rot_line_data(L, 3:5) = rot_point_data(2:4,rot_point_data(1,:)==rot_line_data(L, 1))';
            rot_line_data(L, 6:8) = rot_point_data(2:4,rot_point_data(1,:)==rot_line_data(L, 2))';
        end
        
        % Recheck left leaning walls
        left_lean = zeros(size(rot_line_data,1),1);
        for L = 1:size(rot_line_data,1)
            if L ~= longest_ind && L ~= opp_inds(longest_ind)
                vec_data = sortrows([[rot_line_data(L,3:4), 1]; [rot_line_data(L,6:7), 2]],2);
                if vec_data(2,1) <= vec_data(1,1) % is left leaning
                    left_lean(L) = 1;
                end % otherwise is right leaning
            end
        end
    end

    %% Reanalyze Beam Elements
    for L = 1:size(rot_line_data,1)
        rot_line_data(L,9) = sqrt(sum((rot_line_data(L,3:5) - rot_line_data(L,6:8)).^2)); % length
        rot_line_data(L,10) = 0; % initialize thickness
        rot_line_data(L,11) = atan2d(rot_line_data(L,7)-rot_line_data(L,4),rot_line_data(L,6)-rot_line_data(L,3)); % element angle
        % rot_line_data(L,12) = 0; % initialize alpha
    end

    %% Add Extra Elements and Nodes
    full_point_data = rot_point_data; % initial point data
    full_line_data = [];
    if size(input_line_data,1) == 3 % Triangle
        % fully shifted triangle periodicity
        % initial element data
        full_line_data(:,[1,2,12]) = rot_line_data(:,[1,2,12]); % add elements 1, 2, and 3 and their indexes

        % point data
        full_point_data(1,end+1) = max(full_point_data(1,:))+1; % add 4th point
        [~, ii] = max(full_point_data(3,:));
        full_point_data(2,end) = full_point_data(2,ii) - rot_line_data(longest_ind,9); % left aligned RVE
        full_point_data(3:4,end) = full_point_data(3:4,ii);
    elseif size(input_line_data,1) == 4 % Quadrilateral
        opp_inds = [3, 4, 1, 2]; % opposite side indexes
        BT_ind = [longest_ind, opp_inds(longest_ind)]; % bottom to top element inds
        LR_ind = setdiff(opp_inds,BT_ind)';
        LR_ind = sortrows([LR_ind, mean(rot_line_data(setdiff(opp_inds,BT_ind),[3,6]),2)],2); % sort element inds from left to right
        LR_ind = LR_ind(:,1)';
        
        % initial element data
        full_line_data(end+1,[1,2,12]) = rot_line_data(LR_ind(1),[1,2,12]); % add element 1
        full_line_data(end+1,[1,2,12]) = rot_line_data(LR_ind(2),[1,2,12]); % add element 2

        % point data
        full_point_data(1,end+1) = max(full_point_data(1,:))+1; % add 7th point
        full_point_data(2:4,end) = full_point_data(2:4,intersect(rot_line_data(BT_ind(1),1:2),rot_line_data(LR_ind(2),1:2))==full_point_data(1,:));
        full_point_data(2,end) = full_point_data(2,end) + rot_line_data(BT_ind(2),9);
        % element data
        full_line_data(end+1,1:2) = [full_point_data(1,intersect(rot_line_data(BT_ind(1),1:2),rot_line_data(LR_ind(2),1:2))==full_point_data(1,:)), full_point_data(1,end)]; % add element 3
        full_line_data(end,12) = rot_line_data(opp_inds(longest_ind),12); % add index for element 3

        if sum(left_lean) == 1 % No shift trapezoid (smaller top than bottom base)
            % point data
            full_point_data(1,end+1) = max(full_point_data(1,:))+1; % add 8th point
            full_point_data(2:4,end) = full_point_data(2:4,intersect(rot_line_data(BT_ind(1),1:2),rot_line_data(LR_ind(1),1:2))==full_point_data(1,:));
            full_point_data(3,end) = full_point_data(3,intersect(rot_line_data(BT_ind(2),1:2),rot_line_data(LR_ind(1),1:2))==full_point_data(1,:));
            
            % point data
            full_point_data(1,end+1) = max(full_point_data(1,:))+1; % add 10th point
            full_point_data(2:4,end) = full_point_data(2:4,intersect(rot_line_data(BT_ind(1),1:2),rot_line_data(LR_ind(2),1:2))==full_point_data(1,:));
            full_point_data(2,end) = full_point_data(2,end) + rot_line_data(BT_ind(2),9);
            full_point_data(3,end) = full_point_data(3,intersect(rot_line_data(BT_ind(2),1:2),rot_line_data(LR_ind(2),1:2))==full_point_data(1,:));
            
            % point data
            full_point_data(1,end+1) = max(full_point_data(1,:))+1; % add 5th point
            full_point_data(2:4,end) = full_point_data(2:4,intersect(rot_line_data(BT_ind(2),1:2),rot_line_data(LR_ind(1),1:2))==full_point_data(1,:));
            full_point_data(3,end) = full_point_data(3,intersect(rot_line_data(BT_ind(1),1:2),rot_line_data(LR_ind(1),1:2))==full_point_data(1,:));
            % element data
            full_line_data(end+1,1:2) = [full_point_data(1,intersect(rot_line_data(BT_ind(1),1:2),rot_line_data(LR_ind(1),1:2))==full_point_data(1,:)), full_point_data(1,end)]; % add element 4
            full_line_data(end,12) = rot_line_data(longest_ind,12); % add index for element 4

            % point data
            full_point_data(1,end+1) = max(full_point_data(1,:))+1; % add 6th point
            full_point_data(2:4,end) = full_point_data(2:4,intersect(rot_line_data(BT_ind(2),1:2),rot_line_data(LR_ind(2),1:2))==full_point_data(1,:));
            full_point_data(3,end) = full_point_data(3,intersect(rot_line_data(BT_ind(1),1:2),rot_line_data(LR_ind(2),1:2))==full_point_data(1,:));
            % element data
            full_line_data(end+1,1:2) = [full_point_data(1,end), full_point_data(1,intersect(rot_line_data(BT_ind(1),1:2),rot_line_data(LR_ind(2),1:2))==full_point_data(1,:))]; % add element 5
            full_line_data(end,12) = rot_line_data(longest_ind,12); % add index for element 5
            full_line_data(end+1,1:2) = full_point_data(1,end-1:end); % add element 6
            full_line_data(end,12) = rot_line_data(longest_ind,12); % add index for element 6
            
            % point data
            full_point_data(1,end+1) = max(full_point_data(1,:))+1; % add 9th point
            full_point_data(2:4,end) = full_point_data(2:4,intersect(rot_line_data(BT_ind(1),1:2),rot_line_data(LR_ind(2),1:2))==full_point_data(1,:));
            full_point_data(3,end) = full_point_data(3,intersect(rot_line_data(BT_ind(2),1:2),rot_line_data(LR_ind(2),1:2))==full_point_data(1,:));
        elseif sum(left_lean) == 0 || sum(left_lean) == 2 % left leaning quadrilateral or square (shifted by the degree of lean)
            % point data
            full_point_data(1,end+1) = max(full_point_data(1,:))+1; % add 10th point
            full_point_data(2:4,end) = full_point_data(2:4,intersect(rot_line_data(BT_ind(1),1:2),rot_line_data(LR_ind(2),1:2))==full_point_data(1,:));
            full_point_data(2,end) = full_point_data(2,end) + rot_line_data(BT_ind(2),9) - abs(rot_line_data(LR_ind(1),6) - rot_line_data(LR_ind(1),3));
            full_point_data(3,end) = full_point_data(3,intersect(rot_line_data(BT_ind(2),1:2),rot_line_data(LR_ind(2),1:2))==full_point_data(1,:));
            
            if rot_line_data(BT_ind(1),9) ~= rot_line_data(BT_ind(2),9) % left leaning quadrilateral
                % point data
                full_point_data(1,end+1) = max(full_point_data(1,:))+1; % add 6th point
                full_point_data(2:4,end) = full_point_data(2:4,intersect(rot_line_data(BT_ind(1),1:2),rot_line_data(LR_ind(1),1:2))==full_point_data(1,:));
                full_point_data(2,end) = full_point_data(2,end) + rot_line_data(BT_ind(2),9);
                % element data
                full_line_data(end+1,1:2) = [full_point_data(1,intersect(rot_line_data(BT_ind(1),1:2),rot_line_data(LR_ind(1),1:2))==full_point_data(1,:)), full_point_data(1,end)]; % add element 4
                full_line_data(end,12) = rot_line_data(longest_ind,12); % add index for element 4
                full_line_data(end+1,1:2) = [full_point_data(1,end), full_point_data(1,intersect(rot_line_data(BT_ind(1),1:2),rot_line_data(LR_ind(2),1:2))==full_point_data(1,:))]; % add element 5
                full_line_data(end,12) = rot_line_data(longest_ind,12); % add index for element 5
                
                % point data
                full_point_data(1,end+1) = max(full_point_data(1,:))+1; % add 9th point
                full_point_data(2:4,end) = full_point_data(2:4,intersect(rot_line_data(BT_ind(2),1:2),rot_line_data(LR_ind(2),1:2))==full_point_data(1,:));
                full_point_data(2,end) = full_point_data(2,end) + (rot_line_data(BT_ind(1),9) - rot_line_data(BT_ind(2),9));
            else % square
                % element data
                full_line_data(end+1,1:2) = rot_line_data(BT_ind(1),1:2); % add element 4
                full_line_data(end,12) = rot_line_data(longest_ind,12); % add index for element 4
            end
        end
    end
    full_point_data(3,abs(min(full_point_data(3,:)) - full_point_data(3,:)) < coord_tol) = mean(full_point_data(3,abs(min(full_point_data(3,:)) - full_point_data(3,:)) < coord_tol)); % bottom edge nodes are consistent
    full_point_data(3,abs(max(full_point_data(3,:)) - full_point_data(3,:)) < coord_tol) = mean(full_point_data(3,abs(max(full_point_data(3,:)) - full_point_data(3,:)) < coord_tol)); % top edge nodes are consistent

    %% Assemble Full Line Data
    for L = 1:size(full_line_data,1)
        full_line_data(L, 3:5) = full_point_data(2:4,full_point_data(1,:)==full_line_data(L, 1))';
        full_line_data(L, 6:8) = full_point_data(2:4,full_point_data(1,:)==full_line_data(L, 2))';
    end

    %% Normalize RVE Length Scale
    norm_scale = range(full_point_data(2,abs(min(full_point_data(3,:)) - full_point_data(3,:)) < coord_tol));
    full_point_data(2:4,:) = full_point_data(2:4,:)*(1/norm_scale);
    full_line_data(:,3:8) = full_line_data(:,3:8)*(1/norm_scale);
    
    %% Renumber Points
    line_data = full_line_data;
    renums = [transpose(1:numel(unique(full_point_data(1,:)))), transpose(unique(full_point_data(1,:)))];
    for L = 1:size(line_data,1)
        line_data(L,1) = renums(full_line_data(L,1) == renums(:,2),1);
        line_data(L,2) = renums(full_line_data(L,2) == renums(:,2),1);
    end
    line_data = sortrows(line_data,[1,2]);

    point_data = full_point_data;
    for p = 1:size(full_point_data,2)
        point_data(1,p) = renums(full_point_data(1,p)==renums(:,2),1);
    end
    renums0 = renums; renums = cell(1,2);
    renums{1,1} = renums0; renums{1,2} = line_data(:,[1, 2, 12]);
    
    %% RVE Data Points
    % RVE_points code works for RVEs where the top and bottom boundaries are parallel with the x-axis. The left and right boundaries can be slanted.
    sort_point_data = sortrows(point_data',[3, 2],'ascend')';
    RVE.bot_points = sort_point_data(1,min(sort_point_data(3,:)) == sort_point_data(3,:)); % bottom edge (order sensitive; +x direction left to right)
    RVE.bot_mid_points = RVE.bot_points(2:end-1); % middle points of bottom edge (order sensitive; +x direction left to right)
    RVE.top_points = sort_point_data(1,max(sort_point_data(3,:)) == sort_point_data(3,:)); % top edge (order sensitive; +x direction left to right)
    RVE.top_mid_points = RVE.top_points(2:end-1); % middle points of top edge (order sensitive; +x direction left to right)
    RVE.left_points = [RVE.bot_points(1), RVE.top_points(1)]; % left corner points from the top and bottom edge (order sensitive; +y direction left to right)
    RVE.left_mid_points = []; % middle points of left edge (order sensitive; +y direction left to right)
    for n = 1:size(sort_point_data,2)
        if ~any(sort_point_data(1,n) == RVE.left_points) && all(abs(cross(sort_point_data(2:4,sort_point_data(1,:) == RVE.left_points(2)) - sort_point_data(2:4,sort_point_data(1,:) == RVE.left_points(1)), sort_point_data(2:4,n)-sort_point_data(2:4,sort_point_data(1,:) == RVE.left_points(1)))) < coord_tol)
            RVE.left_mid_points = [RVE.left_mid_points, sort_point_data(1,n)];
        end
    end
    RVE.left_points = [RVE.left_points(1), RVE.left_mid_points, RVE.left_points(end)]; % left edge (order sensitive; +y direction left to right)
    RVE.right_points = [RVE.bot_points(end), RVE.top_points(end)]; % right corner points from the top and bottom edge (order sensitive; +y direction left to right)
    RVE.right_mid_points = []; % middle points of right edge (order sensitive; +y direction left to right)
    for n = 1:size(sort_point_data,2)
        if ~any(sort_point_data(1,n) == RVE.right_points) && all(abs(cross(sort_point_data(2:4,sort_point_data(1,:) == RVE.right_points(2)) - sort_point_data(2:4,sort_point_data(1,:) == RVE.right_points(1)), sort_point_data(2:4,n)-sort_point_data(2:4,sort_point_data(1,:) == RVE.right_points(1)))) < coord_tol)
            RVE.right_mid_points = [RVE.right_mid_points, sort_point_data(1,n)];
        end
    end
    RVE.right_points = [RVE.right_points(1), RVE.right_mid_points, RVE.right_points(end)]; % right edge (order sensitive; +y direction left to right)
    RVE.corner_points = [RVE.bot_points([1,end]), RVE.top_points([1,end])]; % corner points from the top and bottom edge (not order sensitive)
    RVE.int_points = setdiff(point_data(1,:), unique([RVE.corner_points, RVE.bot_points, RVE.top_points, RVE.left_points, RVE.right_points])); % internal points (not order sensitive)
    RVE.eq_point = RVE.bot_points(1); % bottom left corner is point for RVE equilibrium equation (MATTERS FOR EQ 4, NO RIGID BODY MOTION) 
    if numel(RVE.top_mid_points) ~= numel(RVE.bot_mid_points)
        if numel(RVE.top_mid_points) == 0 || numel(RVE.bot_mid_points) == 0
            % warning('RVE boundary points in the middle of the top and bottom edges do not match! Correcting...')
            RVE.top_mid_points = []; RVE.bot_mid_points = [];
        else
            error('RVE boundary points in the middle of the top and bottom edges do not match!')
        end
    end
    if numel(RVE.left_mid_points) ~= numel(RVE.right_mid_points)
        if numel(RVE.left_mid_points) == 0 || numel(RVE.right_mid_points) == 0
            % warning('RVE boundary points in the middle of the left and right edges do not match! Correcting...')
            RVE.left_mid_points = []; RVE.right_mid_points = [];
        else
            error('RVE boundary points in the middle of the left and right edges do not match!')
        end
    end
    
    %% FEA Connectivity Data
    nnodes = numel(unique(full_point_data(1,:)));
    elem_conn = line_data(:,1:2);
    node_conn = cell(nnodes,1);
    for n = 1:nnodes
        % node_conn{n} = sort([elem_conn(elem_conn(:,2)==n,1)', elem_conn(elem_conn(:,1)==n,2)']); % global node numbers that are connected to node n
        node_conn{point_data(1,n)} = sort([elem_conn(elem_conn(:,2)==point_data(1,n),1)', elem_conn(elem_conn(:,1)==point_data(1,n),2)']); % global node numbers that are connected to node n
    end

    %% Reanalyze Beam Elements
    for L = 1:size(line_data,1)
        line_data(L,9) = sqrt(sum((line_data(L,3:5) - line_data(L,6:8)).^2)); % length
        line_data(L,10) = 0; % initialize thickness
        line_data(L,11) = atan2d(line_data(L,7)-line_data(L,4),line_data(L,6)-line_data(L,3)); % element angle
        % line_data(L,12) = 0; % initialize alpha
    end
end






