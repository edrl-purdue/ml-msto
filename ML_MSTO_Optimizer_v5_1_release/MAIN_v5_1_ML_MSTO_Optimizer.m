%% MACHINE LEARNING-BASED MULTISCALE TOPOLOGY OPTIMIZATION v5_0
% By Joel Najmon Ph.D.
% Last date updated: December, 2023
clear; clc; close all; 

%% MACROSCALE OPTIMIZATION SETTINGS
macro.nelx = 6;
macro.nely = 2;
macro.nelz = 0; if macro.nelz == 0; macro.dim = 2; else; macro.dim = 3; end % set to 0 if 2D
macro.volfrac = 0.25;
macro.cont = 0; % Continuation: 1 = Yes, 0 = No (Untested)
macro.flt_den_min = 1.0; %     density filter minimum radius (flt_den_min <= 1 is off)
macro.flt_sen_min = 1.5; % sensitivity filter minimum radius (flt_sen_min <= 1 is off)
macro.delta = 0.5; % Maximum allowable limit in density between adjacent macroscale elements. [0 1] (1 = no constraint), (0 = fixed design)
macro.x_lb(1) = 0; macro.x_ub(1) = 1; % bounds for the density design
macro.x_lb(2) = 0; macro.x_ub(2) = 1; % bounds for the weight design
macro.h_fd = 0; % SA via ANN: h_fd = 0; SA via CFD: h_fd = [0.001, 0.001] (pertubation size for density and weight design variables)
macro.vf_cutoff = 0.1; % cutoff volume fraction for CH interpolation

macro.alg = 3; % alg=1: fmincon, alg=2: OC, alg=3: GOCM
macro.maxloop = 10;   % Maximum number of iterations
macro.tolx = 1e-3;    % Termination criterion
macro.displayflag = 1; % Display structure actively flag

%% SELECT DE-HOMOGENIZATION APPROACH AND SETTINGS
% 0 = No De-homogenization
% 1 = Local topology optimization (TOMs)
macro.dehom = 1;
if macro.dehom == 1
    macro.Micro_scaler = 1; % TOM mesh resolution scaler (scales mesh to the nearest integer)
    macro.epi = 1; % unit cell size (should be 1/integer for dehom==1)
    macro.PP_tri_infill = 0; % Triangular Infilling post processing method. 0 = off, 1 = on
    macro.den_threshold = 0.5; % this density threshold parameter controls the limit used for labelling a corner radii or edge element as solid or void (used for logically identifying solid or void elements)
    macro.PP_rm_skel = 0; % Skeletization post processing method. 0 = off, 1 = on
end
if macro.dehom ~= 0
    macro.PP_rm_low_SE = 0; % Number of low strain energy solid elements removal iterations (requires high-fidelity fea). 0 = off, Non-zero = on
    macro.final_comp = 0; % Calculate final de-homogenization structure's compliance and compatibilty (requires one high-fidelity fea).
end

%% MACROSCALE PHYSICS MODEL
macro.phys = 1; % 1 = Mechanical

%% MACROSCALE OBJECTIVE FUNCTION
macro.obj = 1; % 1 = Compliance

%% DEFINE MACROSCALE BOUNDARY CONDITIONS
macro.bc_type = 1; % TOM Type of Boundary Conditions: 1 = Homogeneous Dirichlet BC, 2 = Non-homogeneous Dirichlet BC
macro.design_ini = 1; % Initial Design: 1 = Uniform, 2 = Normal, 3 = Random Distributions
macro.bc = cell(3,2); % rows = physics modes, columns = [loads, supports]
F0_Macro = cell(1,3); % KU = F for potentially 3 physic models
U0_Macro = cell(1,3);
N_Macro = cell(1,3);
for p = 1:length(macro.phys)
    [macro,F0_Macro{macro.phys(p)},U0_Macro{macro.phys(p)},N_Macro{macro.phys(p)}] = Prepare_BCMacro(macro,macro.phys(p));
end

%% DEFINE MICROSCALE MODEL
[micro, ML_model] = Prepare_Micro_Model(true,macro.dim,macro.phys); % true = use ANNs, false = do not use ANNs (i.e., solve inner optimization problems explicitly).
% Specify saved ann mat file name or micro properties in Prepare_Micro_Model
micro.vf_cutoff = macro.vf_cutoff;

%% MATERIAL PROPERTIES
mat.Emax = 1.0;
mat.Emin = 1e-9;
mat.kmax = 0.5;
mat.kmin = 1e-9;
mat.Amax = 1e-3;
mat.Amin = 1e-9;
mat.nu = 0.3;
macro.mat = mat;

%% PREPARE MESH AND FILTER
[macro.necon, ~] = Prepare_Mesh(macro); % Macro mesh connectivity
macro = mesh_info(macro);
macro.DDx = macro.nelx; macro.DDy = macro.nely; if macro.nelz == 0; macro.DDz = 1; else; macro.DDz = macro.nelz; end % Design Domain size. DDz is thickness in 2D
[macro.den_H, macro.den_Hs] = Prepare_Filter(macro,macro.flt_den_min); % Macro density filter
[macro.sen_H, macro.sen_Hs] = Prepare_Filter(macro,macro.flt_sen_min); % Macro sensitivity filter

%% PREPARE FEA
KE = cell(3,3); % # of rows/cols is the number of unique physics modes || KE_ij: coupling stiffness matrices of physics i to physics j
FE = cell(3,3);
if any(macro.phys == 1)
    [KE{1,1}, FE{1,1},       ~] = Prepare_KE(Prepare_C(mat.nu,1,macro.dim),1,micro); % Mechanical
end
if any(macro.phys == 2)
    [KE{2,2}, FE{2,2},       ~] = Prepare_KE(Prepare_C(mat.nu,2,macro.dim),2,micro); % Thermal
end
if all(macro.phys == [1 2])
    [      ~, FE{2,1}, KE{2,1}] = Prepare_KE(Prepare_C(mat.nu,1,macro.dim),[1,2],micro); % Thermomechanical
end

%% PRINT GENERATE SAMPLE PLAN INFORMATION
fprintf('==========================================================\n');
fprintf(' ML-MSTO Optimizer Program\n');
fprintf('==========================================================\n\n');
display_macro(macro)
t_start = Initialize_Program();

%% PREPARE ALGORITHM
xMacro = ini_design(macro,macro.volfrac,macro.x_lb(1),macro.x_ub(1));
if all(macro.phys == [1 2])
    xWeight = ini_design(macro,mean([macro.x_lb(2) macro.x_ub(2)]),macro.x_lb(2),macro.x_ub(2));
    Xin = [xMacro(:); xWeight(:)];
    [A, B] = density_grad(macro); % density grading inequality constraint
    A = [A, zeros(size(A))];
else
    xWeight = [];
    Xin = xMacro(:);
    [A, B] = density_grad(macro); % density grading inequality constraint
end
Aeq = [ones(1,numel(xMacro)), zeros(1,numel(xWeight))]; % volume fraction equality constraint
Beq = macro.nele*macro.volfrac; % volume fraction equality constraint
LB = [ones(numel(xMacro),1)*macro.x_lb(1); ones(numel(xWeight),1)*macro.x_lb(2)]; % lower design variable bounds
UB = [ones(numel(xMacro),1)*macro.x_ub(1); ones(numel(xWeight),1)*macro.x_ub(2)]; % upper design variable bounds
nonlcon = [];
if macro.delta == 0 % no macroscale optimization, fixed macroscale design
    macro.maxloop = 1;
end

%% INITIAL DESIGN FIELD RESPONSE
c0_time = tic;
[U_Macro,Ue_int_Macro,Fe_int_Macro] = FEA(xMacro,macro,3,mat,KE,F0_Macro,U0_Macro,N_Macro); % SIMP-based initial design field response of macroscale
[c0, ~, ~] = Obj_Fns(xMacro,macro,U_Macro,3,mat,KE,[],[],0.5); % Case 1: Calling Obj_Fns(xMicro,scale,U_Micro,penal,mat,KE,[],[],weight_e) is c and dc for SIMP stiffness matrix
c0_time = toc(c0_time);

%% INITIALIZE OPTIMIZATION HISTORY
opt_history = cell(1,11);
opt_history{end,1} = 0; opt_history{end,2} = xMacro; opt_history{end,3} = xWeight; opt_history{end,4} = c0; opt_history{end,9} = Ue_int_Macro; opt_history{end,10} = Fe_int_Macro; opt_history{end,11} = U_Macro; % save initial macroscale data
macro.rand_id = round(rand*10000);
save(['temp_files/temp_' num2str(macro.nelx) 'x' num2str(macro.nely) 'x' num2str(macro.nelz) '_' num2str(macro.rand_id) '_001' '.mat'],'opt_history'); % created temp files are delete at program's end

%% MULITSCALE TOPOLOGY OPTIMIZATION LOOP
% fprintf('     Solve Outer Optimization Problem:\n');
% fprintf('     -----------------------------------------------------\n');
fprintf(' Solve Outer Optimization Problem:\n');
fprintf('==========================================================\n');
% PREPARE CONTINUATION SCHEME (NOT TESTED YET FOR OUTER OPT PROBLEM)
if macro.cont
    penal_k = 1; % p0
    del_penal = 0.3; % del_p
    penal_k = penal_k:del_penal:micro.penal; % penalization values
    if penal_k(end) ~= micro.penal; penal_k(end+1) = micro.penal; end % add last penalization value is max
    tolx_k = 1e-2; % w0
    tolx_mult = (macro.tolx/tolx_k)^(1/(length(penal_k)-1)); % w multiplier
    for i = 1:(length(penal_k)-1)
        tolx_k = [tolx_k min(tolx_k)*tolx_mult]; %#ok<AGROW>
    end
    loop_k = ones(1,size(tolx_k))*macro.maxloop/10; loop_k(end) = macro.maxloop;
else
    penal_k = micro.penal; % uses microscale penal
    tolx_k = macro.tolx;
    loop_k = macro.maxloop;
end

% SOLVE MACROSCALE OPTIMIZATION PROBLEM
opt_t = tic;
for pk = 1:numel(penal_k)
    macro.penal = penal_k(pk);
    macro.tolx = tolx_k(pk);
    macro.maxloop = loop_k(pk);

    % Define Homogenized Macroscale Function
    fun = @(Xin) Hom_Macro_Fn(Xin,macro,micro,mat,KE,FE,F0_Macro,U0_Macro,N_Macro,ML_model);

    % FMINCON
    if macro.alg == 1
        % Define options
        options = optimoptions('fmincon',...
                               'Algorithm','interior-point',...
                               'MaxIterations',macro.maxloop, ...
                               'SpecifyObjectiveGradient',true,...
                               'SpecifyConstraintGradient',false,...
                               'MaxFunctionEvaluations',(macro.nele+1)*macro.maxloop,...
                               'Display','final-detailed',...
                               'StepTolerance',macro.tolx,... % add termination criteria
                               'PlotFcns',@optimplotfval);
        
        % Solve macroscale optimization problem with fmincon
        Xin = fmincon(fun,Xin,A,B,Aeq,Beq,LB,UB,nonlcon,options);

    % OPTIMALITY CRITERION
    elseif macro.alg ==2
        % Define options
        options.lmin = 0; % lagrange multiplier lower bound
        options.lmax = 1e9; % lagrange multiplier upper bound
        options.move = 0.2; % positive move limit per iteration for design variables
        options.eta = 0.5; % numerical damping coefficient
        options.maxloop = macro.maxloop; % maximum number of iterations
        options.tolx = macro.tolx; % converage tolerance for design variables
        
        % Solve macroscale optimization problem with Optimality Criterion
        Xin = alg_OCM(fun,Xin,A,B,Aeq,Beq,LB,UB,nonlcon,options);

    % GENERALIZED OPTIMALITY CRITERION
    elseif macro.alg == 3
        % Define options
        options.lmin = 1e-3; % lagrange multiplier lower bound
        options.lmax = 1e9; % lagrange multiplier upper bound
        options.move = 0.2; % positive move limit per iteration for design variables
        options.eta = 0.5; % numerical damping coefficient
        options.tolp = 0.05; % tolerance for checking sign of g_change
        options.maxloop = macro.maxloop; % maximum number of iterations
        options.tolx = macro.tolx; % converage tolerance for design variables

        % Solve macroscale optimization problem with Generalized Optimality Criterion
        Xin = alg_GOCM(fun,Xin,A,B,Aeq,Beq,LB,UB,nonlcon,options,macro);
    end
end
fprintf('     Outer Optimization Problem Finished:\n');

%% LOAD OPT_HISTORY
temp_files = dir(['temp_files/temp_' num2str(macro.nelx) 'x' num2str(macro.nely) 'x' num2str(macro.nelz) '_' num2str(macro.rand_id) '_*']);
for f = 1:numel(temp_files)
    load(['temp_files/' temp_files(f).name],'opt_history')
    temp_files(f).opt_history = opt_history;
end
opt_history = cell(1,11);
opt_history{1,1} = 'loop'; opt_history{1,2} = 'xMacro'; opt_history{1,3} = 'xWeight'; opt_history{1,4} = 'c_Macro'; opt_history{1,5} = 'dc_Macro'; opt_history{1,6} = 'ptomtime, satime'; opt_history{1,7} = 'C_Hom'; opt_history{1,8} = 'KE_Hom'; opt_history{1,9} = 'Ue_int_Macro'; opt_history{1,10} = 'Fe_int_Macro'; opt_history{1,11} = 'U_Macro';
for f = 1:numel(temp_files)
    opt_history = [opt_history; temp_files(f).opt_history];
end
output_plot(opt_history{end,2},opt_history{end,3},macro,cell(numel(opt_history{end,2}),1),opt_history{end,4},opt_history{end,1},0);
Xin = opt_history{end,2}(:);

%% NORMALIZE DESIGN VARIABLES IF THEY ARE NEAR BOUNDS
foutput = save_results([],[],[],macro,micro); % Initialize output folder and save macro and micro variables
save_results(foutput,Xin,'Macro_obj_plot',[],[]); % Saves unnormalized Xin and saves fmincon obj fn plot

%% SAVE FINAL MACROSCALE DESIGN DATA
final_history = cell(1,7); final_history{1,1} = 'name'; final_history{1,2} = 'topology'; final_history{1,3} = 'weights'; final_history{1,4} = 'objective'; final_history{1,5} = 'time'; final_history{1,6} = 'compliance ratio'; final_history{1,7} = 'compatibility';
final_history{end+1,1} = 'macroscale_ini'; final_history{end,2} = opt_history{2,2}; final_history{end,3} = opt_history{2,3}; final_history{end,4} = c0; final_history{end,5} = c0_time; final_history{end,6} = c0/c0; % initial macroscale design
final_history{end+1,1} = 'macroscale_final'; final_history{end,2} = reshape(Xin,size(xMacro)); final_history{end,3} = xWeight; final_history{end,4} = opt_history{end,4}; final_history{end,5} = toc(opt_t); final_history{end,6} = final_history{end,4}/c0; % final macroscale design
save_results(foutput,opt_history,[],macro,[]); % Saves opt_history and vtk files

%% DE-HOMOGENIZE FINAL MACROSCALE DESIGN
[dehom_history, final_history] = DeHomogenize_Macro(Xin,macro,micro,mat,KE,opt_history,final_history,foutput);
save_results(foutput,dehom_history,[],[],[]);

%% SAVE FINAL MULTISCALE DESIGN DATA
final_history{end+1,1} = 'total_time'; final_history{end,5} = toc(t_start); % total time it took the ML-MSTO code to run
save_results(foutput,final_history,[],[],[]); % Saves final_history
if macro.displayflag
    figure(1); clf;
    display_top(final_history{3,2}); title('Final Macroscale Design'); drawnow;
    if macro.dehom ~= 0
       figure(3); clf;
       display_top(final_history{5,2}); title('Final Multiscale Design'); drawnow;
    end
end
fprintf('\n');
fprintf('==========================================================\n');
fprintf(' End of ML-MSTO Optimizer Program\n');
fprintf('==========================================================\n');

%% DELETE TEMP FILES
for f = 1:numel(temp_files)
    if exist(['temp_files/' temp_files(f).name])==2 %#ok<*EXIST> 
        delete(['temp_files/' temp_files(f).name]);
    end
end
