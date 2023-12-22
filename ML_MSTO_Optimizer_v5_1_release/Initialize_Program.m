% INITIALIZE ML-MSTO PROGRAM
function t_start = Initialize_Program()
    t_start = tic; rng shuffle

    %% ADD PROGRAM PATHS
    curr_path = pwd;
    addpath(genpath(curr_path));
    
    %% SETUP PARALLEL COMPUTING POOL (SUPERCOMPUTERS ONLY)
    patchJobStorageLocation();
    setup_parpool();
end