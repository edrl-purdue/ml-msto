# ml-msto
Machine learning-based multiscale topology optimization

By: Joel Najmon and Andres Tovar

This is a MATLAB code repository for the machine learning-based multiscale topology optimization (ML-MSTO) method proposed in the 2023 journal paper: "Multiscale Topology Optimization via Artificial Neural Netowrks and Displacement-driven Toplogy-optimized Microstructures" submitted to Structural and Multidisciplinary Optimization (SMO) in June 2023. The paper is currently under review.

Only lines 7-23 of the MAIN_v3_0_ML_MSTO_Optimizer.m file need to be modified to run the examples. If you wish to skip the computationally expensive de-homogenization step then change line 30 from 'macro.dehom = 1;' to 'macro.dehom = 0;'. Results from the program are found in the 'Results' folder.

Note that there is a typo in the manuscript regarding the delta value for the 3D density-graded and infilled design cases. The delta value for these 3D design cases is delta = 0.9 (i.e., it is not the same delta value of delta = 0.5 used in the 2D density-graded and infilled design cases). This correction is reflected in the 3D example instructions below.

Below are instructions on how to run reduced fidelity versions of the 2D and 3D examples found in the submitted manuscript (simply copy and paste the following lines over lines 7-23).
\
\
2D Cantilever Beam (fixed design case): %\
macro.nelx = 60; %\
macro.nely = 30; %\
macro.nelz = 0; if macro.nelz == 0; macro.dim = 2; else; macro.dim = 3; end % set to 0 if 2D %\
macro.volfrac = 0.50; %\
macro.cont = 0; % Continuation: 1 = Yes, 0 = No (Untested) %\
macro.flt_den_min = 1.0; %     density filter minimum radius (flt_den_min <= 1 is off) %\
macro.flt_sen_min = 3.0; % sensitivity filter minimum radius (flt_sen_min <= 1 is off) %\
macro.delta = 0; % Maximum allowable limit in density between adjacent macroscale elements. [0 1] (1 = no constraint), (0 = fixed design) %\
macro.x_lb(1) = 0; macro.x_ub(1) = 1; % bounds for the density design %\
macro.x_lb(2) = 0; macro.x_ub(2) = 1; % bounds for the weight design %\
macro.h_fd = 0; % SA via ANN: h_fd = 0; SA via CFD: h_fd = [0.001, 0.001] (pertubation size for density and weight design variables) %\
macro.vf_cutoff = 0.1; % cutoff volume fraction for CH interpolation %\

macro.alg = 3; % alg=1: fmincon, alg=2: OC, alg=3: GOC, alg=4: MMA %\
macro.maxloop = 100;   % Maximum number of iterations %\
macro.tolx = 1e-3;    % Termination criterion %\
macro.displayflag = 1; % Display structure actively flag %\
\
\
2D Cantilever Beam (density-graded design case): %\
macro.nelx = 60; %\
macro.nely = 30; %\
macro.nelz = 0; if macro.nelz == 0; macro.dim = 2; else; macro.dim = 3; end % set to 0 if 2D %\
macro.volfrac = 0.50; %\
macro.cont = 0; % Continuation: 1 = Yes, 0 = No (Untested) %\
macro.flt_den_min = 1.0; %     density filter minimum radius (flt_den_min <= 1 is off) %\
macro.flt_sen_min = 3.0; % sensitivity filter minimum radius (flt_sen_min <= 1 is off) %\
macro.delta = 0.5; % Maximum allowable limit in density between adjacent macroscale elements. [0 1] (1 = no constraint), (0 = fixed design) %\
macro.x_lb(1) = 0; macro.x_ub(1) = 1; % bounds for the density design %\
macro.x_lb(2) = 0; macro.x_ub(2) = 1; % bounds for the weight design %\
macro.h_fd = 0; % SA via ANN: h_fd = 0; SA via CFD: h_fd = [0.001, 0.001] (pertubation size for density and weight design variables) %\
macro.vf_cutoff = 0.1; % cutoff volume fraction for CH interpolation %\

macro.alg = 3; % alg=1: fmincon, alg=2: OC, alg=3: GOC, alg=4: MMA %\
macro.maxloop = 100;   % Maximum number of iterations %\
macro.tolx = 1e-3;    % Termination criterion %\
macro.displayflag = 1; % Display structure actively flag %\
\
\
2D Cantilever Beam (infilled design case): %\
macro.nelx = 60; %\
macro.nely = 30; %\
macro.nelz = 0; if macro.nelz == 0; macro.dim = 2; else; macro.dim = 3; end % set to 0 if 2D %\
macro.volfrac = 0.50; %\
macro.cont = 0; % Continuation: 1 = Yes, 0 = No (Untested) %\
macro.flt_den_min = 1.0; %     density filter minimum radius (flt_den_min <= 1 is off) %\
macro.flt_sen_min = 3.0; % sensitivity filter minimum radius (flt_sen_min <= 1 is off) %\
macro.delta = 0.5; % Maximum allowable limit in density between adjacent macroscale elements. [0 1] (1 = no constraint), (0 = fixed design) %\
macro.x_lb(1) = 0; macro.x_ub(1) = 0.8; % bounds for the density design %\
macro.x_lb(2) = 0; macro.x_ub(2) = 1; % bounds for the weight design %\
macro.h_fd = 0; % SA via ANN: h_fd = 0; SA via CFD: h_fd = [0.001, 0.001] (pertubation size for density and weight design variables) %\
macro.vf_cutoff = 0.1; % cutoff volume fraction for CH interpolation %\
 
macro.alg = 3; % alg=1: fmincon, alg=2: OC, alg=3: GOC, alg=4: MMA %\
macro.maxloop = 100;   % Maximum number of iterations %\
macro.tolx = 1e-3;    % Termination criterion %\
macro.displayflag = 1; % Display structure actively flag %\
\
\
3D Five-point Bending Plate (infilled design case): %\
macro.nelx = 24; %\
macro.nely = 24; %\
macro.nelz = 4; if macro.nelz == 0; macro.dim = 2; else; macro.dim = 3; end % set to 0 if 2D %\
macro.volfrac = 0.25; %\
macro.cont = 0; % Continuation: 1 = Yes, 0 = No (Untested) %\
macro.flt_den_min = 1.0; %     density filter minimum radius (flt_den_min <= 1 is off) %\
macro.flt_sen_min = 1.2; % sensitivity filter minimum radius (flt_sen_min <= 1 is off) %\
macro.delta = 0; % Maximum allowable limit in density between adjacent macroscale elements. [0 1] (1 = no constraint), (0 = fixed design) %\
macro.x_lb(1) = 0; macro.x_ub(1) = 1; % bounds for the density design %\
macro.x_lb(2) = 0; macro.x_ub(2) = 1; % bounds for the weight design %\
macro.h_fd = 0; % SA via ANN: h_fd = 0; SA via CFD: h_fd = [0.001, 0.001] (pertubation size for density and weight design variables) %\
macro.vf_cutoff = 0.1; % cutoff volume fraction for CH interpolation %\

macro.alg = 3; % alg=1: fmincon, alg=2: OC, alg=3: GOC, alg=4: MMA %\
macro.maxloop = 50;   % Maximum number of iterations %\
macro.tolx = 1e-3;    % Termination criterion %\
macro.displayflag = 1; % Display structure actively flag %\
\
\
3D Five-point Bending Plate (density-graded design case): %\
macro.nelx = 24; %\
macro.nely = 24; %\
macro.nelz = 4; if macro.nelz == 0; macro.dim = 2; else; macro.dim = 3; end % set to 0 if 2D %\
macro.volfrac = 0.25; %\
macro.cont = 0; % Continuation: 1 = Yes, 0 = No (Untested) %\
macro.flt_den_min = 1.0; %     density filter minimum radius (flt_den_min <= 1 is off) %\
macro.flt_sen_min = 1.2; % sensitivity filter minimum radius (flt_sen_min <= 1 is off) %\
macro.delta = 0.9; % Maximum allowable limit in density between adjacent macroscale elements. [0 1] (1 = no constraint), (0 = fixed design) %\
macro.x_lb(1) = 0; macro.x_ub(1) = 1; % bounds for the density design %\
macro.x_lb(2) = 0; macro.x_ub(2) = 1; % bounds for the weight design %\
macro.h_fd = 0; % SA via ANN: h_fd = 0; SA via CFD: h_fd = [0.001, 0.001] (pertubation size for density and weight design variables) %\
macro.vf_cutoff = 0.1; % cutoff volume fraction for CH interpolation %\

macro.alg = 3; % alg=1: fmincon, alg=2: OC, alg=3: GOC, alg=4: MMA %\
macro.maxloop = 50;   % Maximum number of iterations %\
macro.tolx = 1e-3;    % Termination criterion %\
macro.displayflag = 1; % Display structure actively flag %\
\
\
3D Five-point Bending Plate (infilled design case): %\
macro.nelx = 24; %\
macro.nely = 24; %\
macro.nelz = 4; if macro.nelz == 0; macro.dim = 2; else; macro.dim = 3; end % set to 0 if 2D %\
macro.volfrac = 0.25; %\
macro.cont = 0; % Continuation: 1 = Yes, 0 = No (Untested) %\
macro.flt_den_min = 1.0; %     density filter minimum radius (flt_den_min <= 1 is off) %\
macro.flt_sen_min = 1.2; % sensitivity filter minimum radius (flt_sen_min <= 1 is off) %\
macro.delta = 0.9; % Maximum allowable limit in density between adjacent macroscale elements. [0 1] (1 = no constraint), (0 = fixed design) %\
macro.x_lb(1) = 0; macro.x_ub(1) = 0.6; % bounds for the density design %\
macro.x_lb(2) = 0; macro.x_ub(2) = 1; % bounds for the weight design %\
macro.h_fd = 0; % SA via ANN: h_fd = 0; SA via CFD: h_fd = [0.001, 0.001] (pertubation size for density and weight design variables) %\
macro.vf_cutoff = 0.1; % cutoff volume fraction for CH interpolation %\

macro.alg = 3; % alg=1: fmincon, alg=2: OC, alg=3: GOC, alg=4: MMA %\
macro.maxloop = 50;   % Maximum number of iterations %\
macro.tolx = 1e-3;    % Termination criterion %\
macro.displayflag = 1; % Display structure actively flag %\
