% SKELETONIZATION POST PROCESSING
function xMulti = PostProcess_Skel(xMulti)
    xMulti0 = zeros(size(xMulti));

    %% ITERATE THROUGH THE SKELETONIZATION POST PROCESS METHOD
    fprintf('     |       Apply Skeletonization Post Processing       |\n'); tic
    k = 1;
    while any(xMulti(:) ~= xMulti0(:)) && k <= 10
        % Save past iteration
        xMulti0 = xMulti;

        % Remove Islands
        xMulti_rm = remove_islands(xMulti);
    
        % Calculate Distance Transform and Stencils
        [D, Dlist, stencils, stencil_size] = calc_D_stencils(xMulti_rm);
    
        % Generate Skeletons
        [~, xSkel_good, xSkel_bad] = gen_skel(xMulti_rm);
        
        % Remove Skeletons
        [xMulti_er, xSkel_rebuild] = rm_skel(xMulti_rm,xSkel_good,xSkel_bad,D,Dlist,stencils,stencil_size);
        
        % Rebuild Skeletons
        xMulti = rb_skel(xMulti_er,xSkel_rebuild,D,Dlist,stencils,stencil_size);

        k = k + 1;
        clear xMulti_er xMulti_rm xSkel_bad xSkel_good xSkel_rebuild stencils stencil_size D Dlist
    end
    fprintf('     |         Done. Elapsed Time: %4.2f seconds.         |\n',toc); ptomtime = toc;

    %% REMOVE ISLANDS
    function xMulti_rm = remove_islands(xMulti)
        if size(xMulti,3) == 1 % 2D
            filter_neighborhood = 4; % 2d conn
        else % 3D
            filter_neighborhood = 6; % 3d conn
        end
        stats = regionprops(bwconncomp(xMulti,filter_neighborhood),'Area','Image','Centroid','BoundingBox','PixelIdxList');
        area = zeros(numel(stats),1);
        for s = 1:numel(stats)
            area(s) = stats(s).Area;
        end
        [~, ind] = max(area);
        xMulti_rm = zeros(size(xMulti));
        xMulti_rm(stats(ind).PixelIdxList) = 1;
        xMulti_rm = logical(xMulti_rm);
        
        %%%%%%% POSSIBLE IMPROVEMENT
        % could make this part more robust by checking loads/supports of
        % each mass. Keep masses (even if distinct) that are statically
        % determinate with loads and supports. By removing islands based
        % off this criteria and not just connectivity this subroutine is
        % more robust to a variety of topologies and problems (e.g. two
        % distinct columns in a 2D top-distributed load and bottom-distributed support problem)
    end

    %% CALCULATE DISTANCE TRANSFORM AND STENCILS ON XMULTI_RM
    function [D, Dlist, stencils, stencil_size] = calc_D_stencils(xMulti)
        D = bwdist(~xMulti,'euclidean'); % calculate distance transform on input top with islands removed
        Dlist = unique(D(:)); % sort list of all the unique element distances
        template = zeros(2*ceil(max(Dlist)) + 1); % create a square template to calculate the discrete set of circle stencils on. Use the max D value to set the size of the template
        template(ceil(max(Dlist)) + 1,ceil(max(Dlist)) + 1) = 1; % Set center element to 1
        template = bwdist(template); % calculate distance transform to every point on the template
        
        if size(xMulti,3) == 1 % 2D
            filter_neighborhood = 4; % 2d conn
        else % 3D
            filter_neighborhood = 6; % 3d conn
        end

        stencil_size = cell(numel(Dlist),1);
        parfor s = 1:numel(Dlist)
            stats = regionprops(bwconncomp(template <= Dlist(s),filter_neighborhood),'BoundingBox');
            stencil_size{s} = stats.BoundingBox(3:4);
        end

        uni_sizes = (unique(cell2mat(stencil_size))-1)/2;
        stencils = cell(numel(uni_sizes),1);
        parfor s = 1:numel(uni_sizes)
            stencils{s} = double(template((ceil(max(Dlist)) + 1 - uni_sizes(s)):(ceil(max(Dlist)) + 1 + uni_sizes(s)),(ceil(max(Dlist)) + 1 - uni_sizes(s)):(ceil(max(Dlist)) + 1 + uni_sizes(s))));
        end
    end

    %% GENERATE SKELETONS
    function [xSkel_full, xSkel_good, xSkel_bad] = gen_skel(xin)
        xSkel_full = bwskel(xin); % Create full skeleton
        xSkel_good = bwskel(bwskel(bwskel(bwskel(bwskel(bwskel(bwskel(bwskel(bwskel(xin,'MinBranchLength',1/eps),'MinBranchLength',1/eps),'MinBranchLength',1/eps),'MinBranchLength',1/eps),'MinBranchLength',1/eps),'MinBranchLength',1/eps),'MinBranchLength',1/eps),'MinBranchLength',1/eps),'MinBranchLength',1/eps); % Prune full skeleton to arrive at the good skeleton to keep
        xSkel_bad = xSkel_full - xSkel_good; % Subtract skeletons to fine the bad skeleton to remove
    end
    
    %% REMOVE SKELETONS
    function [xMulti_er, xSkel_rebuild] = rm_skel(xMulti_er,xSkel_good,xSkel_bad,D,Dlist,stencils,stencil_size)
        % Eroded topology (xMulti_er) starts as the full topology with removed islands
        xSkel_good_er = xSkel_good;
        xnums = reshape(1:numel(xMulti_er),size(xMulti_er));
        
        data = [xnums(:), D(:), xSkel_bad(:), xMulti_er(:)];
        % data = sortrows(data,[3,2,1],{'descend','descend','ascend'});
        data = sortrows(data,[3,1,2],{'descend','ascend','descend'});
        comp = cell(sum(xSkel_bad(:)),1); % distance comparison cell array
        uni_sizes = unique(cell2mat(stencil_size));
        for L = 1:sum(xSkel_bad(:))
            % Form Indexing Arrays for translate D stencils for distance comparison
            [row,col,page] = ind2sub(size(xMulti_er),data(L,1)); % array index for element L of xSkel_bad
            r_size = (stencil_size{Dlist==data(L,2)} - 1)/2; % get the radius size for element L according to the D stencil
            rows = round((row - r_size(1)):(row + r_size(1)),0); % stencil rows
            row_ind = and(rows >= 1, rows <= size(xMulti_er,1)); % stencil rows that are within bounds
            cols = round((col - r_size(2)):(col + r_size(2)),0); % stencil cols
            col_ind = and(cols >= 1, cols <= size(xMulti_er,2)); % stencil cols that are within bounds
            
            % Compare element distance array to D
            try
                ind = combvec(rows(row_ind),cols(col_ind))';
                comp{L} = sparse(ind(:,1), ind(:,2), stencils{mean(stencil_size{Dlist==data(L,2)})==uni_sizes}(row_ind,col_ind), size(xMulti_er,1), size(xMulti_er,2)); % located elements that are within stencil size of D_L and less than the distance function D at element L.
            catch
                C = cols(col_ind);
                comp{L} = sparse(size(xMulti_er,1), size(xMulti_er,2));
                sten = stencils{mean(stencil_size{Dlist==data(L,2)})==uni_sizes}(row_ind,col_ind);
                for idx=1:sum(col_ind) % loop to avoid large memory usage from combvec
                   comp{L} = comp{L} + sparse(rows(row_ind)', repmat(C(:,idx),1,length(rows(row_ind)))', sten(:,idx), size(xMulti_er,1), size(xMulti_er,2));
                end
            end
        end
        for L = 1:sum(xSkel_bad(:))
            xMulti_er(logical(comp{L})) = 0;
            xSkel_good_er(logical(comp{L})) = 0;
        end
        xSkel_rebuild = xSkel_good - xSkel_good_er;

        %%%%%%% POSSIBLE IMPROVEMENT
        % could also check removals in location to external loads and
        % supports. Structure material should not be eroded if it is
        % connected to a boundary load or support.
    end

    %% REBUILD SKELETONS
    function xMulti_rb = rb_skel(xMulti_rb,xSkel_rebuild,D,Dlist,stencils,stencil_size)
        % Projected topology starts as the full topology
        xnums = reshape(1:numel(xMulti_rb),size(xMulti_rb));
        
        data = [xnums(:), D(:), xSkel_rebuild(:), xMulti_rb(:)];
        % data = sortrows(data,[3,2,1],{'descend','descend','ascend'});
        data = sortrows(data,[3,1,2],{'descend','ascend','descend'});
        comp = cell(sum(xSkel_rebuild(:)),1); % distance comparison cell array
        uni_sizes = unique(cell2mat(stencil_size));
        for L = 1:sum(xSkel_rebuild(:))
            % Form Indexing Arrays for translate D stencils for distance comparison
            [row,col,page] = ind2sub(size(xMulti_rb),data(L,1)); % array index for element L of xSkel_bad
            r_size = (stencil_size{Dlist==data(L,2)} - 1)/2; % get the radius size for element L according to the D stencil
            rows = round((row - r_size(1)):(row + r_size(1)),0); % stencil rows
            row_ind = and(rows >= 1, rows <= size(xMulti_rb,1)); % stencil rows that are within bounds
            cols = round((col - r_size(2)):(col + r_size(2)),0); % stencil cols
            col_ind = and(cols >= 1, cols <= size(xMulti_rb,2)); % stencil cols that are within bounds
        
            % Compare element distance array to D
            ind = combvec(rows(row_ind),cols(col_ind))';
            comp{L} = sparse(ind(:,1), ind(:,2), stencils{mean(stencil_size{Dlist==data(L,2)})==uni_sizes}(row_ind,col_ind), size(xMulti_rb,1), size(xMulti_rb,2)); % located elements that are within stencil size of D_L and less than the distance function D at element L.

            try
                ind = combvec(rows(row_ind),cols(col_ind))';
                comp{L} = sparse(ind(:,1), ind(:,2), stencils{mean(stencil_size{Dlist==data(L,2)})==uni_sizes}(row_ind,col_ind), size(xMulti_rb,1), size(xMulti_rb,2)); % located elements that are within stencil size of D_L and less than the distance function D at element L.
            catch
                C = cols(col_ind);
                comp{L} = sparse(size(xMulti_rb,1), size(xMulti_rb,2));
                sten = stencils{mean(stencil_size{Dlist==data(L,2)})==uni_sizes}(row_ind,col_ind);
                for idx=1:sum(col_ind) % loop to avoid large memory usage from combvec
                   comp{L} = comp{L} + sparse(rows(row_ind)', repmat(C(:,idx),1,length(rows(row_ind)))', sten(:,idx), size(xMulti_rb,1), size(xMulti_rb,2));
                end
            end
        end
        for L = 1:sum(xSkel_rebuild(:))
            xMulti_rb(logical(comp{L})) = 1;
        end
    end
end