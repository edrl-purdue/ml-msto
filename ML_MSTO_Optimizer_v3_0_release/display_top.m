% DISPLAY 2D/3D TOPOLOGIES
function display_top(rho,threshold,view_name)
    % rho: input 2D/3D array for topologies [0, 1]
    % threshold: plotting threshold for rho elements (always used for 3D plots, has to be user-defined for 2D plots)
    % view_name: camera view for 3D only (iso, top, bottom, front, back, left, right)

    %% MANAGE INPUT ARGUMENTS
    if  nargin == 1
        threshold = 0.5; % default 50% plotting threshold
        view_name = 'iso';
    elseif nargin == 2
        if isstring(threshold) || ischar(threshold) % 2nd input is view_name, threshold not defined
            view_name = threshold;
            if size(rho,3) == 1
                threshold = 0.5; % default 50% plotting threshold for 2D
            else
                threshold = mean(rho(:)); % default 50% plotting threshold for 3D
            end
        else % 2nd input is not view_name and assumed to be the threshold
            view_name = 'iso';
        end
    end
    
    %% PLOT 2D AND 3D TOPOLOGIES
    if size(rho,3) == 1 % 2D plot
        if nargin > 1 % threshold is only used for 2D plots that specify it. Otherwise no threshold is used for 2D plots
            rho = rho > threshold;
        end
        colormap(gray); imagesc(1-rho,[0, 1]); caxis([0 1]);
        axis equal; axis off; axis tight; box on; view([-90,90]); drawnow;
    elseif sum(sum(sum(rho > threshold))) ~= 0 % 3D plot (rho > threshold cannot be all zeros)
        patch = PATCH_3Darray(rho > threshold); % from https://www.mathworks.com/matlabcentral/fileexchange/28497-plot-a-3d-array-using-patch?s_tid=ta_fx_results\
        patch.FaceColor = ones(1,3)*0.85; % gray faces
        patch.EdgeColor = ones(1,3)*0.1; % black edges
        camlight('headlight')
        axis equal; axis off; axis tight; box on;
        if strcmp(view_name ,'top')
            view([0,90]);
        elseif strcmp(view_name ,'bottom')
            view([0,-90]);
        elseif strcmp(view_name ,'front')
            view([0,0]);
        elseif strcmp(view_name ,'back')
            view([180,0]);
        elseif strcmp(view_name ,'right')
            view([90,0]);
        elseif strcmp(view_name ,'left')
            view([-90,0]);
        else % default to iso view
            view([20,25]);
        end
        drawnow;
    end
end