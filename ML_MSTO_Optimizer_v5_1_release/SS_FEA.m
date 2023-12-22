% SINGLE-SCALE FINITE ELEMENT ANALYSIS
function [c, c_vec] = SS_FEA(proj_rho,scale,mat)
    try
        % PREPARE MESH AND FILTER
        scale.DDx = scale.nelx; scale.DDy = scale.nely; if scale.nelz == 0; scale.DDz = 1; else; scale.DDz = scale.nelz; end % Design Domain size. DDz is thickness in 2D
        [scale.H,scale.Hs] = Prepare_Filter(scale,0.1); % No filter needed. placer holder H, Hs
        [scale.necon, scale.necon_per] = Prepare_Mesh(scale); % Micro mesh connectivity
        
        % PREPARE FEA
        KE = cell(3,3); % # of rows/cols is the number of unique physics modes || KE_ij: coupling stiffness matrices of physics i to physics j
        if any(scale.phys == 1)
            [KE{1,1}, ~,       ~] = Prepare_KE(Prepare_C(mat.nu,1,scale.dim),1,scale); % Mechanical
        end
        if any(scale.phys == 2)
            [KE{2,2}, ~,       ~] = Prepare_KE(Prepare_C(mat.nu,2,scale.dim),2,scale); % Thermal
        end
        if all(scale.phys == [1 2])
            [      ~, ~, KE{2,1}] = Prepare_KE(Prepare_C(mat.nu,1,scale.dim),[1,2],scale); % Thermomechanical
        end
        
        % PREPARE MACROSCALE BOUNDARY CONDITIONS
        scale.bc = cell(3,2); % rows = physics modes, columns = [loads, supports]
        F0_Macro = cell(1,3); % KU = F for 3 physics
        U0_Macro = cell(1,3);
        N_Macro = cell(1,3);
        for p = 1:length(scale.phys)
            [scale,F0_Macro{scale.phys(p)},U0_Macro{scale.phys(p)},N_Macro{scale.phys(p)}] = Prepare_BCMacro(scale,scale.phys(p));
        end

        % FINITE ELEMENT ANALYSIS
        disp('Solving single-scale FEA '); tic
        [U_Macro,~,~] = FEA(proj_rho,scale,3,mat,KE,F0_Macro,U0_Macro,N_Macro);
        [c, ~, c_vec] = Obj_Fns(proj_rho,scale,U_Macro,3,mat,KE,[],[],0.5); % 0.5 argument is weight between multiphysics if on
        disp(['Finished Solving single-scale FEA. Time: ' num2str(toc) ' seconds.']);
    catch ME
        c_vec = NaN; c = [];
        disp(['Skipping single-scale FEA due to error:' ME.message]); % try-catch is used here in case the system runs-out-of-memory
    end
end