% DE-HOMOGENIZE FINAL MACROSCALE DESIGN
function [dehom_history, final_history] = DeHomogenize_Macro(Xin,macro,micro,mat,KE,opt_history,final_history,foutput)
    if all(macro.phys == [1 2])
        xMacro = Xin(1:length(Xin)/2);
        xWeight = Xin(length(Xin)/2+1:end);
    else
        xMacro = Xin;
        xWeight = ones(size(xMacro));
    end

    %% INITIALIZE DEHOMOGENIZATION HISTORY
    dehom_history = cell(1,7); dehom_history{1,1} = 'Description'; dehom_history{1,2} = 'xMacro'; dehom_history{1,3} = 'xMicro'; dehom_history{1,4} = 'xMulti'; dehom_history{1,5} = 'xWeight'; dehom_history{1,6} = 'Time';  dehom_history{1,7} = 'c_vec';
    dehom_history{end+1,1} = ['Initial, Perform DH_' num2str(macro.dehom)]; dehom_history{end,2} = xMacro; dehom_history{end,5} = xWeight; dhtime = tic;

    %% DE-HOMOGENIZE FINAL MACROSCALE DESIGN
    fprintf(' De-homogenize Macroscale Design:\n');
    fprintf('==========================================================\n');
    if macro.dehom == 0
        fprintf('     |               Skip De-homogenization              |\n');
        xMicro = xMacro;
        xMulti = xMacro;
    else
        if macro.dehom == 1 % Local topology optimization (TOM)

            % 1.a De-homogenize TOMs
            subtime = tic;
            xMicro = DeHomogenize_TOM(xMacro,xWeight,macro,micro,mat,KE,opt_history); % Solve local TOMs
            dehom_history{end+1,1} = 'DH_1a_TOM_Multi'; dehom_history{end,2} = xMacro; dehom_history{end,3} = xMicro; dehom_history{end,4} = cell2mat(xMicro); dehom_history{end,5} = xWeight; dehom_history{end,6} = toc(subtime);
    
            % 1.b Connect TOMs
            subtime = tic;
            if macro.PP_tri_infill
                xMicro = Tri_Infill_Method(macro,micro,xMicro);
            end
            dehom_history{end+1,1} = 'DH_1b_Connected_Multi'; dehom_history{end,2} = xMacro; dehom_history{end,3} = xMicro; dehom_history{end,4} = cell2mat(xMicro); dehom_history{end,5} = xWeight; dehom_history{end,6} = toc(subtime);
            
            % 1.c Scale TOMs
            subtime = tic;
            xMicro = TOM_scaler(macro,micro,xMicro);
            dehom_history{end+1,1} = 'DH_1c_Scaled_Multi'; dehom_history{end,2} = xMacro; dehom_history{end,3} = xMicro; dehom_history{end,4} = cell2mat(xMicro); dehom_history{end,5} = xWeight; dehom_history{end,6} = toc(subtime);

            % Define Multi
            xMulti = cell2mat(xMicro);
            multi = macro;
            multi.nelx = macro.nelx*micro.nelx; multi.nely = macro.nely*micro.nely; multi.nelz = macro.nelz*micro.nelz;
            multi.nele = multi.nelx*multi.nely*max(multi.nelz,1); % number of elements
            multi.DDx = multi.nelx; multi.DDy = multi.nely; if multi.nelz == 0; multi.DDz = 1; else; macro.DDz = macro.nelz; end % Design Domain size. DDz is thickness in 2D
            
        end
        save_results(foutput,dehom_history,[],[],[]);
        final_history{end+1,1} = 'multiscale_ini'; final_history{end,2} = ones(size(xMulti))*macro.volfrac;
        final_history{end+1,1} = 'multiscale_final'; final_history{end,2} = xMulti;
        
        % Post-process final de-homogenized topology
        if (macro.final_comp || macro.PP_rm_low_SE > 0 || macro.PP_rm_skel)
            if macro.final_comp
                subtime = tic;
                [c0, ~] = SS_FEA(final_history{end-1,2},multi,mat); % add command line output eventually for final_comp operations
                final_history{end-1,4} = c0; final_history{end-1,5} = toc(subtime); final_history{end-1,6} = c0/c0;
            end
            
            % Skeletonization post processing
            if macro.PP_rm_skel
                subtime = tic;
                xMulti = PostProcess_Skel(xMulti>macro.den_threshold);
                dehom_history{end+1,1} = 'DH_PP0_Skel_Multi'; dehom_history{end,2} = xMacro; dehom_history{end,3} = xMicro; dehom_history{end,4} = xMulti; dehom_history{end,5} = xWeight; dehom_history{end,6} = toc(subtime);
            end
            
            % xMulti FEA analysis after PPskel, needed for PP_rm_low_SE or final_comp
            if (macro.final_comp || macro.PP_rm_low_SE > 0)
                subtime = tic;
                [c, c_vec] = SS_FEA(xMulti,multi,mat); % initial FEA analysis needed for PP_rm_low_SE or final_comp
            else
                c = [];
            end

            % Remove low strain energy solids
            if macro.PP_rm_low_SE > 0
                filter_neighborhood = [0 1 0; 1 1 1; 0 1 0];
                xMulti(reshape(c_vec(:,1),size(xMulti)) < mean(c_vec(:,1))*0.001) = 0;
                xMulti = imopen(xMulti,filter_neighborhood); % open
                xMulti = imclose(xMulti,filter_neighborhood); % close
                dehom_history{end+1,1} = 'DH_PP1_Post_Multi'; dehom_history{end,2} = xMacro; dehom_history{end,3} = xMicro; dehom_history{end,4} = xMulti; dehom_history{end,5} = xWeight; dehom_history{end,6} = toc(subtime); dehom_history{end,7} = c_vec;

                for f = 2:macro.PP_rm_low_SE % repeat PP_rm_low_SE iterations
                    subtime = tic;
                    [c, c_vec] = SS_FEA(xMulti,multi,mat);
                    xMulti(reshape(c_vec(:,1),size(xMulti)) < mean(c_vec(:,1))*0.001) = 0;
                    xMulti = imopen(xMulti,filter_neighborhood); % open
                    xMulti = imclose(xMulti,filter_neighborhood); % close
                    dehom_history{end+1,1} = ['DH_PP' num2str(f) '_Post_Multi']; dehom_history{end,2} = xMacro; dehom_history{end,3} = xMicro; dehom_history{end,4} = xMulti; dehom_history{end,5} = xWeight; dehom_history{end,6} = toc(subtime); dehom_history{end,7} = c_vec;
                end
            end
            final_history{end,4} = c; % save final multiscale topology's compliance
        end
        final_history{end,2} = xMulti; % save final multiscale topology

        % Calculate Compatibility
        final_history{end,5} = toc(dhtime);
        if macro.final_comp
            final_history{end,6} = c/c0; % Multiscale Compliance Ratio
            final_history{end,7} = c/opt_history{end,4}; % DH Compatibility
        end
    end
    dehom_history{end+1,1} = 'DH_Final_Multi'; dehom_history{end,2} = xMacro; dehom_history{end,3} = xMicro; dehom_history{end,4} = xMulti; dehom_history{end,5} = xWeight; dehom_history{end,6} = toc(dhtime);
end