% PREPARE MESH CONNECTIVITY MATRICES
function [necon, necon_per] = Prepare_Mesh(scale)
    % Node numbers for non-periodic mesh
    [Y, X, Z] = meshgrid(0:scale.nely,0:scale.nelx,0:scale.nelz); % nonperiodic mesh
    nnums = X + Y*(scale.nelx+1) + Z*(scale.nely+1)*(scale.nelx+1);

    % Node numbers for periodic mesh
    nelx_1 = scale.nelx-1; nely_1 = scale.nely-1; nelz_1 = max(scale.nelz-1,0);
    [Y_per, X_per, Z_per] = meshgrid(0:nely_1,0:nelx_1,0:nelz_1); % sub periodic mesh
    nnums_per = X_per + Y_per*(nelx_1+1) + Z_per*(nely_1+1)*(nelx_1+1); % nnums created with nel-1 (sub periodic mesh)
    nnums_per(end+1,:,:) = nnums_per(1,:,:); nnums_per(:,end+1,:) = nnums_per(:,1,:); nnums_per(:,:,end+1) = nnums_per(:,:,1); % periodic mesh

    % Create connectivity matrices for non-periodic and periodic meshes
    if scale.dim == 2 % 2D
        nnums_per = nnums_per(:,:,1);
        nnums1 = nnums(1:end-1,1:end-1); nnums1 = nnums1(:);
        nnums2 = nnums(2:end,1:end-1);   nnums2 = nnums2(:);
        nnums3 = nnums(2:end,2:end);     nnums3 = nnums3(:);
        nnums4 = nnums(1:end-1,2:end);   nnums4 = nnums4(:);
        necon = [nnums1, nnums2, nnums3, nnums4]'; necon = necon(:);
        nnums_per1 = nnums_per(1:end-1,1:end-1); nnums_per1 = nnums_per1(:);
        nnums_per2 = nnums_per(2:end,1:end-1);   nnums_per2 = nnums_per2(:);
        nnums_per3 = nnums_per(2:end,2:end);     nnums_per3 = nnums_per3(:);
        nnums_per4 = nnums_per(1:end-1,2:end);   nnums_per4 = nnums_per4(:); 
        necon_per = [nnums_per1, nnums_per2, nnums_per3, nnums_per4]'; necon_per = necon_per(:);
    elseif scale.dim == 3 % 3D
        nnums1 = nnums(1:end-1,1:end-1,1:end-1); nnums1 = nnums1(:);
        nnums2 = nnums(2:end,1:end-1,1:end-1);   nnums2 = nnums2(:);
        nnums3 = nnums(2:end,2:end,1:end-1);     nnums3 = nnums3(:);
        nnums4 = nnums(1:end-1,2:end,1:end-1);   nnums4 = nnums4(:);
        nnums5 = nnums(1:end-1,1:end-1,2:end);   nnums5 = nnums5(:);
        nnums6 = nnums(2:end,1:end-1,2:end);     nnums6 = nnums6(:);
        nnums7 = nnums(2:end,2:end,2:end);       nnums7 = nnums7(:);
        nnums8 = nnums(1:end-1,2:end,2:end);     nnums8 = nnums8(:);
        necon = [nnums1, nnums2, nnums3, nnums4, nnums5, nnums6, nnums7, nnums8]'; necon = necon(:);
        nnums_per1 = nnums_per(1:end-1,1:end-1,1:end-1); nnums_per1 = nnums_per1(:);
        nnums_per2 = nnums_per(2:end,1:end-1,1:end-1);   nnums_per2 = nnums_per2(:);
        nnums_per3 = nnums_per(2:end,2:end,1:end-1);     nnums_per3 = nnums_per3(:);
        nnums_per4 = nnums_per(1:end-1,2:end,1:end-1);   nnums_per4 = nnums_per4(:);
        nnums_per5 = nnums_per(1:end-1,1:end-1,2:end);   nnums_per5 = nnums_per5(:);
        nnums_per6 = nnums_per(2:end,1:end-1,2:end);     nnums_per6 = nnums_per6(:);
        nnums_per7 = nnums_per(2:end,2:end,2:end);       nnums_per7 = nnums_per7(:);
        nnums_per8 = nnums_per(1:end-1,2:end,2:end);     nnums_per8 = nnums_per8(:);
        necon_per = [nnums_per1, nnums_per2, nnums_per3, nnums_per4, nnums_per5, nnums_per6, nnums_per7, nnums_per8]'; necon_per = necon_per(:);
    end
end