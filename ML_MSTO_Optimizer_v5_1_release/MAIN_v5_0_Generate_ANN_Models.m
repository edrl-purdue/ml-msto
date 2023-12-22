%% GENERATE SAMPLE PLAN AND TRAIN ANN MODELS
% By Joel Najmon Ph.D.
% Last date updated: December, 2023
clear; clc; close all; 

%% SELECT SAMPLE PLAN SIZE
micro.N = 1000; %[Number of Sample Points]
micro.rr = 1.05; %[N Multiplier for Redundant Samples]

%% SELECT PHYSICAL DIMENSION
micro.dim = 2; % 2D or 3D

%% SURROGATE MODEL SETTINGS
n_uni_train = micro.N; % number of unique samples for training, N = max
if micro.dim == 2
    micro.n_aug = 8;
elseif micro.dim == 3
    micro.n_aug = 48; % Number of augmentations to use in ANN training. Integer between 1-8 (2D) or 1-48 (3D).
end
NL = 1; % number of hidden layers
NN = 2; % number of neurons per hidden layer
reps = 8; % number of times to generate an ANN to get the best one

%% SELECT PHYSICS MODELS
% 1 = Mechanical
% 2 = Thermal
% 3 = Fluid
micro.phys = [1];

%% MICROSCALE OPTIMIZATION PARAMETERS
micro.nelx = 100; if rem(micro.nelx,2) ~= 0; micro.nelx = micro.nelx + 1; end % nels should be an even number
micro.nely = micro.nelx; if micro.dim == 2; micro.nelz = 0; elseif micro.dim == 3; micro.nelz = micro.nelx; end
micro.penal = 3; micro.cont = 1; % Continuation: 1 = Yes, = 0 No % CURRENTLY IT IS ASSUMED THAT ALL PHYS MODES USE THE SAME PENALIZATION VALUE
micro.flt_den_min = micro.nelx/20; %     density filter minimum radius (flt_den_min <= 1 is off)
micro.flt_sen_min = 1.0; % sensitivity filter minimum radius (flt_sen_min <= 1 is off)
micro.x_lb(1) = 0; micro.x_ub(1) = 1; % bounds for the initial design (final designs are always binary bounds)
micro.x_lb(2) = 0; micro.x_ub(2) = 1; % bounds for the weight design
micro.ortho = 0; % ortho = 1: orthotropic, ortho = 0: anisotropic

micro.maxloop = 50;    % Maximum number of iterations
micro.tolx = 1e-4;    % Termination criterion
micro.displayflag = 0; % Display structure flag

% SELECT MICROSCALE BOUNDARY CONDITIONS
micro.bc_type = 2; % TOM Type of Boundary Conditions: 1 = Internal Force BC, 2 = Internal Displacement BC
micro.load_bc = []; % TOM Load Boundary Conditions: 1 = Corners, 2 = Midpoints of Edges, 3 = Edges, 4 = Face, 5 = Midpoints of Faces
micro.supp_bc = [1]; % TOM Support Boundary Conditions: 1 = Corners, 2 = Midpoints of Edges, 3 = Edges, 4 = Face, 5 = Midpoints of Faces
micro.design_ini = 1; % TOM Initial Design: 1 = Uniform, 2 = Normal, 3 = Random Distributions
if micro.bc_type == 2; micro.load_bc = []; end

%% TYPE OF PRE-DEFINED PARAMETERIZED MICROSTRUCTURE
% 1 = Unit cell with rectangular void (2 - parameters + rotation)
% 2 = Unit cell with cross (2 - parameters + rotation)
% 3 = Unit cell with anisotruss box/cross (6 - parameters + rotation)
micro.unit_cell_type = 1;
if micro.unit_cell_type == 1 || micro.unit_cell_type == 2; micro.para_num = 3; elseif micro.unit_cell_type == 3; micro.para_num = 7; end

%% MATERIAL PROPERTIES
mat.Emax = 1.0;
mat.Emin = 1e-9;
mat.nu = 0.3;
mat.kmax = 0.5;  % Good thermal conductivity
mat.kmin = 1e-9; % Poor thermal conductivity
mat.Amax = 1e-3;
mat.Amin = 1e-9;
micro.mat = mat; % save mat into micro

%% PREPARE MESH AND FILTER
[micro.necon, micro.necon_per] = Prepare_Mesh(micro); % Micro mesh connectivity
micro = mesh_info(micro);
micro.DDx = 1; micro.DDy = 1; micro.DDz = 1; % Design Domain size. DDz is thickness in 2D
[micro.den_H, micro.den_Hs] = Prepare_Filter(micro,micro.flt_den_min); % Micro density filter
[micro.sen_H, micro.sen_Hs] = Prepare_Filter(micro,micro.flt_sen_min); % Micro sensitivity filter

%% PREPARE FEA
KE = cell(3,3); FE = cell(3,3); Cmin = cell(3,3); Cmax = cell(3,3); % # of rows/cols is the number of unique physics modes || KE_ij: coupling stiffness matrices of physics i to physics j
if any(micro.phys == 1)
    [KE{1,1}, FE{1,1},       ~] = Prepare_KE(Prepare_C(mat.nu,1,micro.dim),1,micro); % Mechanical
    Cmin{1,1} = micro.mat.Emin * Prepare_C(micro.mat.nu,1,micro.dim);
    Cmax{1,1} = micro.mat.Emax * Prepare_C(micro.mat.nu,1,micro.dim);
end
if any(micro.phys == 2)
    [KE{2,2}, FE{2,2},       ~] = Prepare_KE(Prepare_C(mat.nu,2,micro.dim),2,micro); % Thermal
    Cmin{2,2} = micro.mat.kmin * Prepare_C(micro.mat.nu,2,micro.dim);
    Cmax{2,2} = micro.mat.kmax * Prepare_C(micro.mat.nu,2,micro.dim);
end
if all(micro.phys == [1 2])
    [      ~, FE{2,1}, KE{2,1}] = Prepare_KE(Prepare_C(mat.nu,1,micro.dim),[1,2],micro); % Thermomechanical
    Cmin{2,1} = micro.mat.Amin * Prepare_C(micro.mat.nu,1,micro.dim); % might not be finished yet
    Cmax{2,1} = micro.mat.Amax * Prepare_C(micro.mat.nu,1,micro.dim); % might not be finished yet
end

%% RECORD C RANGE
[ml_index, ~, phys_index] = Indexing('part',micro,micro.ortho);
outnum = 0; for p = 1:size(phys_index,1); for ind = 1:size(ml_index{phys_index(p,1),phys_index(p,2)},1); outnum = outnum + 1; end; end
an = 1; Cminvec = zeros(1,outnum); Cmaxvec = zeros(1,outnum);
for p = 1:size(phys_index,1)
    for ind = 1:size(ml_index{phys_index(p,1),phys_index(p,2)},1)
        Cminvec(1,an) = Cmin{phys_index(p,1),phys_index(p,2)}(ml_index{phys_index(p,1),phys_index(p,2)}(ind,1),ml_index{phys_index(p,1),phys_index(p,2)}(ind,2));
        Cmaxvec(1,an) = Cmax{phys_index(p,1),phys_index(p,2)}(ml_index{phys_index(p,1),phys_index(p,2)}(ind,1),ml_index{phys_index(p,1),phys_index(p,2)}(ind,2));
        an = an + 1;
    end
end
micro.Cminvec = Cminvec; micro.Cmaxvec = Cmaxvec;

%% GENERATE SAMPLE PLAN INPUTS
[SP_para, micro] = Prepare_Sample_Plan(micro);

%% CONTINUATION INITIALIZATION
if micro.cont
    penal_k = 1; % p0
    del_penal = 0.3; % del_p
    penal_k = penal_k:del_penal:micro.penal; % penalization values
    if penal_k(end) ~= micro.penal; penal_k(end+1) = micro.penal; end % add last penalization value is max
    tolx_k = 1e-2; % w0
    tolx_mult = (micro.tolx/tolx_k)^(1/(length(penal_k)-1)); % w multiplier
    for i = 1:(length(penal_k)-1)
        tolx_k = [tolx_k min(tolx_k)*tolx_mult]; %#ok<AGROW>
    end
    loop_k = ones(size(tolx_k))*micro.maxloop/10; loop_k(end) = micro.maxloop;
else
    penal_k = micro.penal;
    tolx_k = micro.tolx;
    loop_k = micro.maxloop;
end

%% PRINT GENERATE SAMPLE PLAN INFORMATION
fprintf('==========================================================\n');
fprintf(' Generate Sample Plan and Train ANNs\n');
fprintf('==========================================================\n\n');
fprintf(' ANN Settings:\n');
fprintf('==========================================================\n');
fprintf('     Number of Sample Points: N = %1.2e\n',micro.N);
fprintf('     Number of Augmentations per TOM: %2d\n',micro.n_aug);
fprintf('     Number of Hidden Layers: %3d\n',NL);
fprintf('     Number of Nodes per Hidden Layer: %3d\n',NN);
fprintf('     Number of Training Repetitions: %3d\n\n',reps);
display_micro(micro)
Initialize_Program();

%% GENERATE SAMPLE PLAN OUTPUTS
fprintf(' Generate Sample Plan, Augment TOMs, and Train ANNs:\n');
fprintf('==========================================================\n');
SP_theta = cell(ceil(micro.N*micro.rr),1); % Microstructure Output
SP_C = cell(ceil(micro.N*micro.rr),1); % Homogenized Tensor Output
SP_solve_time = zeros(ceil(micro.N*micro.rr),1); % Average TOM Solution Time Output
SP_C_check = cell(ceil(micro.N*micro.rr),1); % Check if sample is valid
WaitMessage = parfor_waitbar(ceil(micro.N*micro.rr)); % Add ('Waitbar', true) as arguments to turn on a waitbar popup
fprintf('     |      Solve Inner Optimization Problems (TOMs)     |\n'); tic
for i = 1:ceil(micro.N*micro.rr)
    TOM_time = tic;
    if SP_para(i,1) < micro.nelx/micro.nele
        % Design
        SP_theta{i} = zeros(micro.nelx,micro.nely,max(1,micro.nelz));
        % Homogenized Stiffness Tensors
        SP_C{i} = cell(3,3);
        if any(micro.phys == 1)
            SP_C{i}{1,1} = mat.Emin * Prepare_C(mat.nu,1,micro.dim);
        end
        if any(micro.phys == 2)
            SP_C{i}{2,2} = mat.kmin * Prepare_C(mat.nu,2,micro.dim);
        end
        if 2 == sum(micro.phys == [1 2])
            SP_C{i}{2,1} = mat.Amin * Prepare_C(mat.nu,1,micro.dim); % might not be finished yet
        end
        SP_C_check{i} = 1;
    else
        % Design
        SP_para(i,1) = 0.5;
        micro.design_ini = 3;
        SP_theta{i} = ini_design(micro,SP_para(i,1),micro.x_lb(1),micro.x_ub(1));
        for pk = 1:numel(penal_k)
            SP_theta{i} = Micro_TO(micro,SP_para(i,:),SP_theta{i},penal_k(pk),loop_k(pk),tolx_k(pk),mat,KE);
        end
        % Homogenized Stiffness Tensors
        SP_C{i} = Num_Hom(SP_theta{i}>0.48,micro,mat,KE,FE);
        SP_C_check{i} = Check_CH(SP_C{i},micro,mat,SP_para(i,1)); % returns 1 then CHi is valid, returns 0 then CHi is not valid
    end
    SP_solve_time(i) = toc(TOM_time);
    WaitMessage.Send;
end
WaitMessage.Destroy;
fprintf('     |         Done. Elapsed Time: %4.2f seconds.         |\n\n',toc/micro.rr);
sptime = toc/micro.rr; % Total time to solve all TOMs
avgTOMtime = mean(SP_solve_time); % Average time to solve one TOM

%% REMOVE INVALID SAMPLES
SP_C_check = cell2mat(SP_C_check);
% fprintf('\n NUMBER OF FAILED SAMPLES: %d\n',ceil(micro.N*micro.rr)-sum(SP_C_check));
SP_para(SP_C_check==0,:) = [];
SP_theta(SP_C_check==0) = [];
SP_C(SP_C_check==0) = [];

if sum(SP_C_check) < micro.N
    error('\n TOO MANY INVALID SAMPLES. INCREASE rr VALUE.')
else
    del_nums = sum(SP_C_check) - micro.N;
    del_nums = randperm(sum(SP_C_check),del_nums);
end
SP_para(del_nums,:) = [];
SP_theta(del_nums) = [];
SP_C(del_nums) = [];

%% SAVE SAMPLE PLAN
curr_path = pwd;
model_name = 'Homog';
if all(micro.phys == 1)
    physname = 'Phys001_Mecha';
    physshort = 'Mecha';
elseif all(micro.phys == 2)
    physname = 'Phys002_Therm';
    physshort = 'Therm';
elseif all(micro.phys == 3)
    physname = 'Phys003_Fluid';
    physshort = 'Fluid';
elseif all(micro.phys == [1 2])
    physname = 'Phys012_ThMec';
    physshort = 'ThMec';
elseif all(micro.phys == [2 3])
    physname = 'Phys023_ThFlu';
    physshort = 'ThFlu';
end
% Output File Name
existing_file = dir([curr_path '/ML_Models/' physname '/' model_name num2str(micro.dim) 'D_*']);
if numel(existing_file) == 0
    fnum = 1; zeros_num = '00';
else
    file_num = cell(numel(existing_file),1);
    for i = 1:numel(existing_file)
        file_num{i} = existing_file(i).name;
        file_num{i} = str2double(file_num{i}(9:11));
    end
    fnum = max(cell2mat(file_num))+1;
    if length(int2str(fnum)) == 1
        zeros_num = '00';
    elseif length(int2str(fnum)) == 2
        zeros_num = '0';
    elseif length(int2str(fnum)) == 3
        zeros_num = [];
    end
end
fname = [model_name num2str(micro.dim) 'D_' zeros_num num2str(fnum) '_' physshort '_N' num2str(micro.N) '_nels' num2str(micro.nelx) '.mat'];
% Save Sample Plan
save([curr_path '/ML_Models/' physname '/' fname],'SP_para','SP_C','sptime','avgTOMtime','micro','-v7.3');

%% GENERATE SURROGATE MODEL
Train_ANNs([curr_path '/ML_Models/' physname '/' fname],SP_para,SP_C,micro,n_uni_train,NN,NL,reps);

fprintf('==========================================================\n');
fprintf(' End of Generate Sample Plans and Train ANNs Program\n');
fprintf('==========================================================\n');