% LOAD MICROSCALE DATA STRUCTURE AND ANN MODEL
function [micro, ML_model] = Prepare_Micro_Model(useANN,dim,phys)
    if useANN % LOAD MICROSCALE SURROGATE MODEL
        if all(phys == [1])
            if dim == 2
                ann_name = 'Homog2D_001_Mecha_N1000000_nels50';
            elseif dim == 3
                ann_name = 'Homog3D_001_Mecha_N1000000_nels10';
            end
            load(['ML_Models/Phys001_Mecha/' ann_name '.mat'],'ML_model','micro'); % Load Machine Learning Model
        elseif all(phys == [2])
            if dim == 2
                ann_name = 'Homog2D_011_Therm_N1000000_nels50';
            elseif dim == 3
                ann_name = 'Homog3D_011_Therm_N1000000_nels50';
            end
            load(['ML_Models/Phys002_Therm/' ann_name '.mat'],'ML_model','micro'); % Load Machine Learning Model
        elseif all(phys == [1 2])
            if dim == 2
                ann_name = 'Homog2D_011_ThMec_N1000000_nels50';
            elseif dim == 3
                ann_name = 'Homog3D_011_ThMec_N1000000_nels50';
            end
            load(['ML_Models/Phys012_ThMec/' ann_name '.mat'],'ML_model','micro'); % Load Machine Learning Model
        end
        micro.ann_name = ann_name;
    else % DEFINE MICROSCALE OPTIMIZATION PARAMETERS
        ML_model = []; micro.ann_name = []; % Machine Learning Model
        if dim == 2
            micro.nelx = 50;
        elseif dim == 3
            micro.nelx = 10;
        end
        if rem(micro.nelx,2) ~= 0; micro.nelx = micro.nelx + 1; end % nels should be an even number
        micro.dim = dim; micro.phys = phys; micro.nely = micro.nelx; if micro.dim == 2; micro.nelz = 0; elseif micro.dim == 3; micro.nelz = micro.nelx; end
        micro.penal = 3; micro.cont = 1; % Continuation: 1 = Yes, = 0 No % CURRENTLY IT IS ASSUMED THAT ALL PHYS MODES USE THE SAME PENALIZATION SCHEME
        micro.flt_den_min = micro.nelx/20; %     density filter minimum radius (flt_den_min <= 1 is off)
        micro.flt_sen_min = 1.0; % sensitivity filter minimum radius (flt_sen_min <= 1 is off)
        micro.xMicro_lb = 0; micro.xMicro_ub = 1; % bounds for the initial design (final designs are always binary bounds)
        
        micro.maxloop = 50;    % Maximum number of iterations
        micro.tolx = 1e-4;    % Termination criterion
        micro.displayflag = 0; % Display structure flag
    
        micro.bc_type = 2; % TOM Type of Boundary Conditions: 1 = Internal Force BC, 2 = Internal Displacement BC
        micro.load_bc = []; % TOM Load Boundary Conditions: 1 = Corners, 2 = Midpoints of Edges, 3 = Edges, 4 = Face, 5 = Midpoints of Faces
        micro.supp_bc  = [1]; % TOM Support Boundary Conditions: 1 = Corners, 2 = Midpoints of Edges, 3 = Edges, 4 = Face, 5 = Midpoints of Faces
        micro.design_ini = 1; % TOM Initial Design: 1 = Uniform, 2 = Normal, 3 = Random Distributions
        micro.ML_model_type = 1;
        if all(phys == [1])
            if dim == 2
                micro.spdim_nums = [2, 0, 0; 10, 0, 0]; % TOM Dimension Indexing
                micro.spdim = 10; % Number of TOM Dimensions
            elseif dim == 3
                micro.spdim_nums = [2, 0, 0; 26, 0, 0]; % TOM Dimension Indexing
                micro.spdim = 26; % Number of TOM Dimensions
            end
        elseif all(phys == [2])
            if dim == 2
                micro.spdim_nums = [0, 2, 0; 0, 6, 0]; % TOM Dimension Indexing
                micro.spdim = 6; % Number of TOM Dimensions
            elseif dim == 3
                
            end
        elseif all(phys == [1 2])
            if dim == 2
                micro.spdim_nums = [2, 10, 0; 10, 14, 0]; % TOM Dimension Indexing
                micro.spdim = 14; % Number of TOM Dimensions
            elseif dim == 3
                
            end
        end
        if micro.bc_type == 2; micro.load_bc = []; end
    
        [micro.necon, micro.necon_per] = Prepare_Mesh(micro); % Micro mesh connectivity
        micro = mesh_info(micro);
        micro.DDx = 1; micro.DDy = 1; micro.DDz = 1; % Design Domain size. DDz is thickness in 2D
        [micro.den_H, micro.den_Hs] = Prepare_Filter(micro,micro.flt_den_min); % Micro density filter
        [micro.sen_H, micro.sen_Hs] = Prepare_Filter(micro,micro.flt_sen_min); % Micro sensitivity filter
    end
    micro.ml = useANN; % record if an ANN was loaded
end