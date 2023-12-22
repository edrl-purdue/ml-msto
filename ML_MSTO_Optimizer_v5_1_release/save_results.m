% SAVE ML-MSTO RESULTS
function foutput = save_results(foutput,var,fig_name,macro,micro)
    % Initialize Call: foutput = save_results([],[],[],macro,micro) % Creates output folder and initial save file
    % Save Variable Call: save_results(foutput,var,[],[],[])
    % Save Figure Call: save_results(foutput,[],fig_name,[],[])
    % Save Variable and Figure Call: save_results(foutput,var,fig_name,[],[])

    %% CREATE OUTPUT FOLDER
    if isempty(foutput) && isempty(micro) % ssto only (no micro)
        if all(macro.phys == 1)
            physshort = 'Mecha';
        elseif all(macro.phys == 2)
            physshort = 'Therm';
        elseif all(macro.phys == 3)
            physshort = 'Fluid';
        elseif all(macro.phys == [1 2])
            physshort = 'ThMec';
        elseif all(macro.phys == [2 3])
            physshort = 'ThFlu';
        end
        % Output Folder Name
        existing_folder = dir(['Results/' 'SSTO_' physshort '_' num2str(macro.dim) 'D_*']);
        if numel(existing_folder) == 0
            fnum = 1; zeros_num = '00';
        else
            folder_num = cell(numel(existing_folder),1);
            for i = 1:numel(existing_folder)
                folder_num{i} = existing_folder(i).name;
                folder_num{i} = str2double(folder_num{i}(15:17));
            end
            fnum = max(cell2mat(folder_num))+1;
            if length(int2str(fnum)) == 1
                zeros_num = '00';
            elseif length(int2str(fnum)) == 2
                zeros_num = '0';
            elseif length(int2str(fnum)) == 3
                zeros_num = [];
            end
        end
        foutput = ['SSTO_' physshort '_' num2str(macro.dim) 'D_' zeros_num num2str(fnum) '_(nelx' num2str(macro.nelx) '_nely' num2str(macro.nely) '_nelz' num2str(macro.nelz) '_VF' num2str(macro.volfrac) '_denmin' num2str(macro.flt_den_min) '_senmin' num2str(macro.flt_sen_min) '_del' num2str(macro.delta)];
        if ~exist([cd '/Results/' foutput], 'dir')
            mkdir([cd '/Results/' foutput]);
        end
        % Save Initial Output File
        save(['Results/' foutput '/MSTO_Output.mat'],'macro','micro','-v7.3');
    elseif isempty(foutput)
        if all(micro.phys == 1)
            physshort = 'Mecha';
        elseif all(micro.phys == 2)
            physshort = 'Therm';
        elseif all(micro.phys == 3)
            physshort = 'Fluid';
        elseif all(micro.phys == [1 2])
            physshort = 'ThMec';
        elseif all(micro.phys == [2 3])
            physshort = 'ThFlu';
        end
        % Output Folder Name
        existing_folder = dir(['Results/' 'MSTO_' physshort '_' num2str(micro.dim) 'D_*']);
        if numel(existing_folder) == 0
            fnum = 1; zeros_num = '00';
        else
            folder_num = cell(numel(existing_folder),1);
            for i = 1:numel(existing_folder)
                folder_num{i} = existing_folder(i).name;
                folder_num{i} = str2double(folder_num{i}(15:17));
            end
            fnum = max(cell2mat(folder_num))+1;
            if length(int2str(fnum)) == 1
                zeros_num = '00';
            elseif length(int2str(fnum)) == 2
                zeros_num = '0';
            elseif length(int2str(fnum)) == 3
                zeros_num = [];
            end
        end
        foutput = ['MSTO_' physshort '_' num2str(macro.dim) 'D_' zeros_num num2str(fnum) '_(nelx' num2str(macro.nelx) '_nely' num2str(macro.nely) '_nelz' num2str(macro.nelz) '_VF' num2str(macro.volfrac) '_denmin' num2str(macro.flt_den_min) '_senmin' num2str(macro.flt_sen_min) '_del' num2str(macro.delta)];
        if ~exist([cd '/Results/' foutput], 'dir')
            mkdir([cd '/Results/' foutput]);
        end
        % Save Initial Output File
        save(['Results/' foutput '/MSTO_Output.mat'],'macro','micro','-v7.3');
    end

    %% SAVE MSTO VARIABLES
    if ~isempty(var)
        S.(inputname(2)) = var;
        save(['Results/' foutput '/MSTO_Output.mat'],'-append','-struct','S');
    end
    
    %% SAVE MATLAB FIGURE
    if ~isempty(fig_name)
%         print(['Results/' foutput '/' fig_name '.png'], '-dpng', '-r900');
    end
end