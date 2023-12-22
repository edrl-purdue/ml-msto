% TRIANGULAR INFILL CONNECTIVITY POST-PROCESSING
function xMicro = Tri_Infill_Method(macro,micro,xMicro)
    %% LINE CONNECTIVITY
    if macro.dim == 2
        % indexing for the microscale corners (counting +x then +y then +z order)
        mask_indexing = [micro.nelx+0.5, micro.nely+0.5, 0.5, 0, 0, 0, -1; % corner 4 of elem 1 [x, y, z, x_size, y_size, z_size, rad_density]
                         0.5, micro.nely+0.5, 0.5, 0, 0, 0, -1; % corner 3 of elem 2 [x, y, z, x_size, y_size, z_size, rad_density]
                         micro.nelx+0.5, 0.5, 0.5, 0, 0, 0, -1; % corner 2 of elem 3 [x, y, z, x_size, y_size, z_size, rad_density]
                         0.5, 0.5, 0.5, 0, 0, 0, -1]; % corner 1 of elem 4 [x, y, z, x_size, y_size, z_size, rad_density]
    elseif macro.dim == 3
        % indexing for the microscale corners (counting +x then +y then +z order)
        mask_indexing = [micro.nelx+0.5, micro.nely+0.5, micro.nelz+0.5, 0, 0, 0, -1; % corner 8 of elem 1 [x, y, z, x_size, y_size, z_size, rad_density]
                         0.5, micro.nely+0.5, micro.nelz+0.5, 0, 0, 0, -1; % corner 7 of elem 2 [x, y, z, x_size, y_size, z_size, rad_density]
                         micro.nelx+0.5, 0.5, micro.nelz+0.5, 0, 0, 0, -1; % corner 6 of elem 3 [x, y, z, x_size, y_size, z_size, rad_density]
                         0.5, 0.5, micro.nelz+0.5, 0, 0, 0, -1; % corner 5 of elem 4 [x, y, z, x_size, y_size, z_size, rad_density]
                         micro.nelx+0.5, micro.nely+0.5, 0.5, 0, 0, 0, -1; % corner 4 of elem 5 [x, y, z, x_size, y_size, z_size, rad_density]
                         0.5, micro.nely+0.5, 0.5, 0, 0, 0, -1; % corner 3 of elem 6 [x, y, z, x_size, y_size, z_size, rad_density]
                         micro.nelx+0.5, 0.5, 0.5, 0, 0, 0, -1; % corner 2 of elem 7 [x, y, z, x_size, y_size, z_size, rad_density]
                         0.5, 0.5, 0.5, 0, 0, 0, -1]; % corner 1 of elem 8 [x, y, z, x_size, y_size, z_size, rad_density]
    end
    
    %% PERFORM LINE CONNECTIVITY METHOD AT ALL OF THE INNER NODES n
    elem_nums = reshape(1:macro.nele,size(xMicro)); elem_nums = elem_nums(:); % list of element numbers
    node_nums = reshape(1:(macro.nelx+1)*(macro.nely+1)*(macro.nelz+1),size(xMicro)+1); % array of nodal numbers
    if macro.dim == 2
        inner_nums = node_nums(2:end-1,2:end-1); % array of inner node numbers
    elseif macro.dim == 3
        inner_nums = node_nums(2:end-1,2:end-1,2:end-1); % array of inner node numbers
    end
    macro_necon = reshape(macro.necon+1,[macro.nen,macro.nele])'; % macroscale element connectivity matrix
    step = [-1, 1]; % used for vector indexing for counting up or down when finding edge
    eind = [2, 1]; % used for element indexing for selecting adjacent elements
    for n = transpose(inner_nums(:)) % loop through inner nodes
        elem_n = [elem_nums(any(macro_necon == n,2)), zeros(macro.nen,7)]; % Locate the elements connected at inner node n; preallocate extra space for relevant nodal data

        %% CALCULATE OPEN END SIZE AT INNER NODE n FOR EACH CONNECTED ELEMENT e
        conn_data = repmat(mask_indexing,[1,1,size(elem_n,1)]); % duplicate the mask index array for calculating open ends
        for e = 1:size(elem_n,1) % loop through the elements touching
            corn_ind = zeros(size(xMicro{1})); % corner index for the corner of e that is connected to node n
            corn_ind(interp1([1, micro.nelx], [1, micro.nelx], mask_indexing(e,1), 'nearest', 'extrap'),...
                     interp1([1, micro.nely], [1, micro.nely], mask_indexing(e,2), 'nearest', 'extrap'),...
                     interp1([1, micro.nelz], [1, micro.nelz], mask_indexing(e,3), 'nearest', 'extrap')) = 1;
            mask = bwdist(corn_ind) <= micro.rmin; % mask for microscale elements within rmin radius of corner that is connected to node n
            conn_data(e,7,e) = mean(xMicro{elem_n(e,1)}(mask)); % only the corner under consideration will have a non -1 value for its page e
            
            %% CALCULATE OPEN END SIZE AT ALL NODES OF ELEMENT e
            for en = 1:macro.nen
                % Find edges/faces that are connected to node en in the x-direction
                xen = interp1([1, micro.nelx], [1, micro.nelx], conn_data(all(conn_data(en,[2,3],e) == conn_data(:,[2,3],e),2),1:3,e), 'nearest', 'extrap'); % coordinates of x line to node en
                if ~all(xen(1,:) == interp1([1, micro.nelx], [1, micro.nelx], conn_data(en,1:3,e), 'nearest', 'extrap'))
                    xen = flip(xen,1); % flip order of the two points if the en node is not listed first
                end
                x_edge = reshape(xMicro{elem_n(e,1)}(xen(1,1):step([xen(1,1)>xen(2,1),xen(1,1)<=xen(2,1)]):xen(2,1),...
                                                     xen(1,2):step([xen(1,2)>xen(2,2),xen(1,2)<=xen(2,2)]):xen(2,2),...
                                                     xen(1,3):step([xen(1,3)>xen(2,3),xen(1,3)<=xen(2,3)]):xen(2,3)),[micro.nelx,1]) > macro.den_threshold;
                if all(x_edge)
                    conn_data(en,4,e) = micro.nelx-1; % edge size in x-direction
                else
                    conn_data(en,4,e) = (find(x_edge==0,1,'first') - 1); % edge size in x-direction
                end

                % Find edges that are connected to node en in the y-direction
                yen = interp1([1, micro.nely], [1, micro.nely], conn_data(all(conn_data(en,[1,3],e) == conn_data(:,[1,3],e),2),1:3,e), 'nearest', 'extrap'); % coordinates of y line to node en
                if ~all(yen(1,:) == interp1([1, micro.nely], [1, micro.nely], conn_data(en,1:3,e), 'nearest', 'extrap'))
                    yen = flip(yen,1); % flip order of the two points if the en node is not listed first
                end
                y_edge = reshape(xMicro{elem_n(e,1)}(yen(1,1):step([yen(1,1)>yen(2,1),yen(1,1)<=yen(2,1)]):yen(2,1),...
                                                     yen(1,2):step([yen(1,2)>yen(2,2),yen(1,2)<=yen(2,2)]):yen(2,2),...
                                                     yen(1,3):step([yen(1,3)>yen(2,3),yen(1,3)<=yen(2,3)]):yen(2,3)),[micro.nely,1]) > macro.den_threshold;
                if all(y_edge)
                    conn_data(en,5,e) = micro.nely-1; % edge size in y-direction
                else
                    conn_data(en,5,e) = (find(y_edge==0,1,'first') - 1); % edge size in y-direction
                end

                % Find edges that are connected to node en in the z-direction
                if macro.dim == 3
                    zen = interp1([1, micro.nelz], [1, micro.nelz], conn_data(all(conn_data(en,[1,2],e) == conn_data(:,[1,2],e),2),1:3,e), 'nearest', 'extrap'); % coordinates of z line to node en
                    if ~all(zen(1,:) == interp1([1, micro.nelz], [1, micro.nelz], conn_data(en,1:3,e), 'nearest', 'extrap'))
                        zen = flip(zen,1); % flip order of the two points if the en node is not listed first
                    end
                    z_edge = reshape(xMicro{elem_n(e,1)}(zen(1,1):step([zen(1,1)>zen(2,1),zen(1,1)<=zen(2,1)]):zen(2,1),...
                                                         zen(1,2):step([zen(1,2)>zen(2,2),zen(1,2)<=zen(2,2)]):zen(2,2),...
                                                         zen(1,3):step([zen(1,3)>zen(2,3),zen(1,3)<=zen(2,3)]):zen(2,3)),[micro.nelz,1]) > macro.den_threshold;
                    if all(z_edge)
                        conn_data(en,6,e) = micro.nelz-1; % edge size in z-direction
                    else
                        conn_data(en,6,e) = (find(z_edge==0,1,'first') - 1); % edge size in z-direction
                    end
                end
            end
            
            % Extract relevant nodal data en from each element e that is connected to inner node n
            elem_n(e,2:end) = conn_data(conn_data(:,7,e) >= 0,:,e); % [elem #, n_corner_x, n_corner_y, n_corner_z, size_x, size_y, size_z, rad_density]
        end
        
        %% FILL VOID AREAS USING TRIANGULAR POINTS
        elem_grid = reshape(elem_n(:,1),2*ones(1,macro.dim)); % element numbers of node in 2x2 or 2x2x2 grid
        void_e = zeros(size(elem_grid)); void_e_nums = elem_n(~(elem_n(:,8) > macro.den_threshold),1);
        for i = 1:numel(void_e_nums)
            void_e = or(void_e,elem_grid == void_e_nums(i));
        end
        [row, col, page] = ind2sub(size(void_e),find(void_e)); % find elem_grid indices for element with void corner
        for i = 1:numel(row) % loop through the number of void corners at inner node n
            if macro.dim == 2 &&...
                elem_n(elem_n(:,1)==elem_grid(eind(row(i)),col(i)),8) > macro.den_threshold &&... % check if neighboring element in x-direction has a full corner
                elem_n(elem_n(:,1)==elem_grid(row(i),eind(col(i))),8) > macro.den_threshold       % check if neighboring element in y-direction has a full corner
                % x points for triangular mask
                x_points = [elem_n(elem_n(:,1)==elem_grid(row(i),col(i)),2),... % corner point x
                            elem_n(elem_n(:,1)==elem_grid(row(i),col(i)),2),... % y edge point x
                            plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i)),2),elem_n(elem_n(:,1)==elem_grid(row(i),eind(col(i))),5),[0.5, micro.nelx+0.5])]; % x edge point x
                % y points for triangular mask
                y_points = [elem_n(elem_n(:,1)==elem_grid(row(i),col(i)),3),... % corner point y
                            plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i)),3),elem_n(elem_n(:,1)==elem_grid(eind(row(i)),col(i)),6),[0.5, micro.nely+0.5]),... % y edge point y
                            elem_n(elem_n(:,1)==elem_grid(row(i),col(i)),3)]; % x edge point y
                % Modify triangular points if the triangle spans over half of a edge's length
                if range(x_points) > micro.nely/2 && range(y_points) > micro.nely/2
                    x_points(3) = plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i)),2),micro.nelx/2,[0.5, micro.nelx+0.5]); % use half of the edge as the x edge
                    y_points(2) = plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i)),3),micro.nely/2,[0.5, micro.nely+0.5]); % use half of the edge as the y edge
                elseif range(x_points) > micro.nely/2
                    x_points(3) = plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i)),2),elem_n(elem_n(:,1)==elem_grid(eind(row(i)),col(i)),6),[0.5, micro.nely+0.5]); % use the y size for the x edge
                elseif range(y_points) > micro.nely/2
                    y_points(2) = plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i)),3),elem_n(elem_n(:,1)==elem_grid(row(i),eind(col(i))),5),[0.5, micro.nelx+0.5]); % use the x size for the y edge
                end
                tri_points{1} = [y_points', x_points', ones(size(x_points'))]; % xy plane mask
                tri_points{2} = [y_points', x_points', zeros(size(x_points'))]; % xy plane mask

                % Fill in triangular mask
                xMicro{elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),1)}(blendedPolymask(tri_points,1:micro.nely,1:micro.nelx,1:max(1,micro.nelz))) = 1;
            elseif macro.dim == 3 &&...
                elem_n(elem_n(:,1)==elem_grid(eind(row(i)),col(i),page(i)),8) > macro.den_threshold &&... % check if neighboring element in x-direction has a full corner
                elem_n(elem_n(:,1)==elem_grid(row(i),eind(col(i)),page(i)),8) > macro.den_threshold &&... % check if neighboring element in y-direction has a full corner
                elem_n(elem_n(:,1)==elem_grid(row(i),col(i),eind(page(i))),8) > macro.den_threshold       % check if neighboring element in z-direction has a full corner
                % x points for triangular mask
                x_points = [elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),2),... % corner point, x
                            elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),2),... % y edge point, x
                            plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),2),min(elem_n(or(elem_n(:,1)==elem_grid(row(i),eind(col(i)),page(i)),elem_n(:,1)==elem_grid(row(i),col(i),eind(page(i)))),5)),[0.5, micro.nelx+0.5])]; % x edge point, x
                % y points for triangular mask
                y_points = [elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),3),... % corner point, y
                            plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),3),min(elem_n(or(elem_n(:,1)==elem_grid(eind(row(i)),col(i),page(i)),elem_n(:,1)==elem_grid(row(i),col(i),eind(page(i)))),6)),[0.5, micro.nely+0.5]),... % y edge point, y
                            elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),3)]; % x edge point, y
                % z points for triangular mask
                z_points = ones(size(x_points))*elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),4); % [corner point, y edge point, x edge point], z
                % Modify triangular points if the triangle spans over half of a edge's length
                if range(x_points) > micro.nely/2 && range(y_points) > micro.nely/2
                    x_points(3) = plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),2),micro.nelx/2,[0.5, micro.nelx+0.5]); % use half of the edge as the x edge
                    y_points(2) = plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),3),micro.nely/2,[0.5, micro.nely+0.5]); % use half of the edge as the y edge
                elseif range(x_points) > micro.nely/2
                    x_points(3) = plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),2),elem_n(elem_n(:,1)==elem_grid(eind(row(i)),col(i),page(i)),6),[0.5, micro.nely+0.5]); % use the y size for the x edge
                elseif range(y_points) > micro.nely/2
                    y_points(2) = plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),3),elem_n(elem_n(:,1)==elem_grid(row(i),eind(col(i)),page(i)),5),[0.5, micro.nelx+0.5]); % use the x size for the y edge
                end
                tri_points{1} = [y_points', x_points', z_points']; % xy plane mask (large)

                % x points for triangular mask
                x_points = [elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),2),... % corner point, x
                            elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),2),... % y edge point, x
                            plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),2),2,[0.5, micro.nelx+0.5])]; % x edge point, x
                % y points for triangular mask
                y_points = [elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),3),... % corner point, y
                            plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),3),2,[0.5, micro.nely+0.5]),... % y edge point, y
                            elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),3)]; % x edge point, y
                % z points for triangular mask
                z_points = ones(size(x_points))*plusminus(elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),3),min(elem_n(or(elem_n(:,1)==elem_grid(eind(row(i)),col(i),page(i)),elem_n(:,1)==elem_grid(row(i),eind(col(i)),page(i))),7)),[0.5, micro.nelz+0.5]); % [corner point, y edge point, x edge point], z
                tri_points{2} = [y_points', x_points', z_points']; % xy plane mask (point)
                
                % Fill in triangular mask
                xMicro{elem_n(elem_n(:,1)==elem_grid(row(i),col(i),page(i)),1)}(blendedPolymask(tri_points,1:micro.nely,1:micro.nelx,1:max(1,micro.nelz))) = 1;
            end
        end
    end

    %% PLUS/MINUS FUNCTION WITH LIMIT CHECKING
    function results = plusminus(A,B,limits) % results A + B or A - B depending on which is inside the limits of [lower, upper]; A ~= B
        if A + B > limits(2)
            results = A - B;
        elseif A - B < limits(1)
            results = A + B;
        end
    end
end