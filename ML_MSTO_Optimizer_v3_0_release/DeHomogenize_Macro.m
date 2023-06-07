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
%             if macro.displayflag
%                 figure(10); clf;
%                 display_top(cell2mat(xMicro)); title('DH 1a TOM Multi'); drawnow;
%                 save_results(foutput,[],'DH_1a_TOM_Multi',[],[]);
%             end
    
            % 1.b Connect TOMs
            subtime = tic;
            if macro.tri_infill
                xMicro = Tri_Infill_Method(macro,micro,xMicro);
            end
            dehom_history{end+1,1} = 'DH_1b_Connected_Multi'; dehom_history{end,2} = xMacro; dehom_history{end,3} = xMicro; dehom_history{end,4} = cell2mat(xMicro); dehom_history{end,5} = xWeight; dehom_history{end,6} = toc(subtime);
%             if macro.displayflag
%                 figure(10); clf;
%                 display_top(cell2mat(xMicro)); title('DH 1b Connected Multi'); drawnow;
%                 save_results(foutput,[],'DH_1b_Connected_Multi',[],[]);
%             end
            
            % 1.c Scale TOMs
            subtime = tic;
            xMicro = TOM_scaler(macro,micro,xMicro);
            dehom_history{end+1,1} = 'DH_1c_Scaled_Multi'; dehom_history{end,2} = xMacro; dehom_history{end,3} = xMicro; dehom_history{end,4} = cell2mat(xMicro); dehom_history{end,5} = xWeight; dehom_history{end,6} = toc(subtime);
%             if macro.displayflag
%                 figure(10); clf;
%                 display_top(cell2mat(xMicro)); title('DH 1c Scaled Multi'); drawnow;
%                 save_results(foutput,[],'DH_1c_Scaled_Multi',[],[]);
%             end

            % Define Multi
            xMulti = cell2mat(xMicro);
            multi = macro;
            multi.nelx = macro.nelx*micro.nelx; multi.nely = macro.nely*micro.nely; multi.nelz = macro.nelz*micro.nelz;
            multi.nele = multi.nelx*multi.nely*max(multi.nelz,1); % number of elements
            multi.DDx = multi.nelx; multi.DDy = multi.nely; if multi.nelz == 0; multi.DDz = 1; else; macro.DDz = macro.nelz; end % Design Domain size. DDz is thickness in 2D
        elseif macro.dehom == 2 % Projection Method with Principal Stress Direction (Mechanical only)
            [xMulti, phi_f, eta_f, multi, phase_shift_f] = DeHomogenize_Proj(xMacro,macro,micro,cell(numel(xMacro),1),opt_history);
        elseif macro.dehom == 3 % Projection Method with Topology Approximation (anisotruss)
    %         xMicro = DeHomogenize_TOM(xMacro,xWeight,macro,micro,mat,KE,opt_history);
    %         xMulti = DeHomogenize_SDF(Xin,macro,micro,xMicro);
    
    %         [xMulti, phi_f, eta_f, multi, phase_shift_f] = DeHomogenize_Proj(Xin,macro,micro,xMicro,opt_history);
    
            [xMulti, phi_f, eta_f, multi, phase_shift_f] = DeHomogenize_Proj_SDF(Xin,macro,micro,opt_history);


%             if exist('phi_f','var')
%                 final_history{end+1,1} = 'phi_f, eta_f, phase_shift_f'; final_history{end,2} = phi_f; final_history{end,3} = eta_f; final_history{end,4} = phase_shift_f;
%             end
        elseif macro.dehom == 4 % Machine Learning Generated Topology
            
        end
        save_results(foutput,dehom_history,[],[],[]);
        final_history{end+1,1} = 'multiscale_ini'; final_history{end,2} = ones(size(xMulti))*macro.volfrac;
        final_history{end+1,1} = 'multiscale_final'; final_history{end,2} = xMulti;
        
        % Post-process final de-homogenized topology
        if (macro.final_comp || macro.rm_low_SE > 0 || macro.rm_skel)
            if macro.final_comp
                subtime = tic;
                [c0, ~] = SS_FEA(final_history{end-1,2},multi,mat); % add command line output eventually for final_comp operations
                final_history{end-1,4} = c0; final_history{end-1,5} = toc(subtime); final_history{end-1,6} = c0/c0;
            end
            
            % Skeletonization post processing
            if macro.rm_skel
                subtime = tic;
                xMulti = PostProcess_Skel(xMulti>macro.den_threshold);
                dehom_history{end+1,1} = 'DH_PP0_Skel_Multi'; dehom_history{end,2} = xMacro; dehom_history{end,3} = xMicro; dehom_history{end,4} = xMulti; dehom_history{end,5} = xWeight; dehom_history{end,6} = toc(subtime);
                if macro.displayflag
                    figure(10); clf;
                    display_top(xMulti); title('DH PP0 Skel Multi'); drawnow;
                    save_results(foutput,[],'DH_PP0_Skel_Multi',[],[]);
                end
            end
            
            % xMulti FEA analysis after PPskel, needed for rm_low_SE or final_comp
            if (macro.final_comp || macro.rm_low_SE > 0)
                subtime = tic;
                [c, c_vec] = SS_FEA(xMulti,multi,mat); % initial FEA analysis needed for rm_low_SE or final_comp
            else
                c = [];
            end

            % Remove low strain energy solids
            if macro.rm_low_SE > 0
                filter_neighborhood = [0 1 0; 1 1 1; 0 1 0];
                xMulti(reshape(c_vec(:,1),size(xMulti)) < mean(c_vec(:,1))*0.001) = 0;
                xMulti = imopen(xMulti,filter_neighborhood); % open
                xMulti = imclose(xMulti,filter_neighborhood); % close
                dehom_history{end+1,1} = 'DH_PP1_Post_Multi'; dehom_history{end,2} = xMacro; dehom_history{end,3} = xMicro; dehom_history{end,4} = xMulti; dehom_history{end,5} = xWeight; dehom_history{end,6} = toc(subtime); dehom_history{end,7} = c_vec;
                if macro.displayflag
                    figure(10); clf;
                    display_top(xMulti); title('DH PP1 Post Multi'); drawnow;
                    save_results(foutput,[],'DH_PP1_Post_Multi',[],[]);
                end

                for f = 2:macro.rm_low_SE % repeat rm_low_SE iterations
                    subtime = tic;
                    [c, c_vec] = SS_FEA(xMulti,multi,mat);
                    xMulti(reshape(c_vec(:,1),size(xMulti)) < mean(c_vec(:,1))*0.001) = 0;
                    xMulti = imopen(xMulti,filter_neighborhood); % open
                    xMulti = imclose(xMulti,filter_neighborhood); % close
                    dehom_history{end+1,1} = ['DH_PP' num2str(f) '_Post_Multi']; dehom_history{end,2} = xMacro; dehom_history{end,3} = xMicro; dehom_history{end,4} = xMulti; dehom_history{end,5} = xWeight; dehom_history{end,6} = toc(subtime); dehom_history{end,7} = c_vec;
                    if macro.displayflag
                        figure(10); clf;
                        display_top(xMulti); title(['DH PP' num2str(f) ' Post Multi']); drawnow;
                        save_results(foutput,[],['DH_PP' num2str(f) '_Post_Multi'],[],[]);
                    end
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