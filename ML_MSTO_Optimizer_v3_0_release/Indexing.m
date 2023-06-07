% GENERATE INDEXING NUMBERS FOR CH TENSORS AND PHYSICS MODELS
function [ml_index, ch_index, phys_index] = Indexing(indtype,scale,ortho)
    
    %% indtype = 'part'
    % 2D Ortho Mecha: C11, C21, C22, C33
    % 2D Ortho Therm: k11, k22
    % 2D Ortho ThMec: A11, A21
    % 2D Aniso Mecha: C11, C21, C22, C31, C32, C33
    % 2D Aniso Therm: k11, k21, k22
    % 2D Aniso ThMec: A11, A21, A31

    % 3D Ortho Mecha: C11, C21, C22, C31, C32, C33, C44, C55, C66
    % 3D Ortho Therm: k11, k22, k33
    % 3D Ortho ThMec: A11, A21, A31
    % 3D Aniso Mecha: C11, C21, C22, C31, C32, C33, C41, C42, C43, C44, C51, C52, C53, C54, C55, C61, C62, C63, C64, C65, C66
    % 3D Aniso Therm: k11, k21, k22, k31, k32, k33
    % 3D Aniso ThMec: A11, A21, A31, A41, A51, A61

    %% indtype = 'full'
    % 2D Ortho Mecha: C11, C12, C21, C22, C33
    % 2D Ortho Therm: k11, k22
    % 2D Ortho ThMec: A11, A21
    % 2D Aniso Mecha: C11, C12, C13, C21, C22, C23, C31, C32, C33
    % 2D Aniso Therm: k11, k21, k21, k22
    % 2D Aniso ThMec: A11, A21, A31

    % 3D Ortho Mecha: C11, C12, C13, C21, C22, C33, C31, C32, C33, C44, C55, C66
    % 3D Ortho Therm: k11, k22, k33
    % 3D Ortho ThMec: A11, A21, A31
    % 3D Aniso Mecha: C11, C12, C13, C14, C15, C16, C21, C22, C23, C24, C25, C26, C31, C32, C33, C34, C35, C36, C41, C42, C43, C44, C45, C46, C51, C52, C53, C54, C55, C56, C61, C62, C63, C64, C65, C66
    % 3D Aniso Therm: k11, k12, k13, k21, k22, k23, k31, k32, k33
    % 3D Aniso ThMec: A11, A21, A31, A41, A51, A61

    %% CH Tensor Indexing
    ml_index_part = cell(3,3); % for ML model indexing
    ml_index_full = cell(3,3); % for CH indexing
    ch_index_full = cell(3,3); % for CH indexing
    
    num_LC(:,:,1) = [3, 1; 1, 2]; num_LC(:,:,2) = [6, 1; 1, 3]; % number of load cases for FE. Rows/Cols are physics modes. Pages are scale.dimension, 2D/3D
    
    ml_index_part{1,1} = []; ml_index_full{1,1} = []; ch_index_full{1,1} = [];
    for y = 1:num_LC(1,1,scale.dim-1)
        for x = 1:num_LC(1,1,scale.dim-1)
            if ortho
                if (x <= scale.dim && y <= scale.dim) || y == x
                    ch_index_full{1,1}(size(ch_index_full{1,1},1)+1,1) = y;
                    ch_index_full{1,1}(size(ch_index_full{1,1},1),2) = x;
        
                    if x <= y
                        ml_index_part{1,1}(size(ml_index_part{1,1},1)+1,1) = y;
                        ml_index_part{1,1}(size(ml_index_part{1,1},1),2) = x;
                        ml_index_full{1,1}(size(ml_index_full{1,1},1)+1,1) = y;
                        ml_index_full{1,1}(size(ml_index_full{1,1},1),2) = x;
                    else
                        ml_index_full{1,1}(size(ml_index_full{1,1},1)+1,1) = x;
                        ml_index_full{1,1}(size(ml_index_full{1,1},1),2) = y;
                    end
                end
            else
                ch_index_full{1,1}(size(ch_index_full{1,1},1)+1,1) = y;
                ch_index_full{1,1}(size(ch_index_full{1,1},1),2) = x;
    
                if x <= y
                    ml_index_part{1,1}(size(ml_index_part{1,1},1)+1,1) = y;
                    ml_index_part{1,1}(size(ml_index_part{1,1},1),2) = x;
                    ml_index_full{1,1}(size(ml_index_full{1,1},1)+1,1) = y;
                    ml_index_full{1,1}(size(ml_index_full{1,1},1),2) = x;
                else
                    ml_index_full{1,1}(size(ml_index_full{1,1},1)+1,1) = x;
                    ml_index_full{1,1}(size(ml_index_full{1,1},1),2) = y;
                end
            end
        end
    end

    ml_index_part{2,2} = []; ml_index_full{2,2} = []; ch_index_full{2,2} = [];
    for y = 1:num_LC(2,2,scale.dim-1)
        for x = 1:num_LC(2,2,scale.dim-1)
            if ortho
                if y == x
                    ch_index_full{2,2}(size(ch_index_full{2,2},1)+1,1) = y;
                    ch_index_full{2,2}(size(ch_index_full{2,2},1),2) = x;
        
                    if x <= y
                        ml_index_part{2,2}(size(ml_index_part{2,2},1)+1,1) = y;
                        ml_index_part{2,2}(size(ml_index_part{2,2},1),2) = x;
                        ml_index_full{2,2}(size(ml_index_full{2,2},1)+1,1) = y;
                        ml_index_full{2,2}(size(ml_index_full{2,2},1),2) = x;
                    else
                        ml_index_full{2,2}(size(ml_index_full{2,2},1)+1,1) = x;
                        ml_index_full{2,2}(size(ml_index_full{2,2},1),2) = y;
                    end
                end
            else
                ch_index_full{2,2}(size(ch_index_full{2,2},1)+1,1) = y;
                ch_index_full{2,2}(size(ch_index_full{2,2},1),2) = x;
    
                if x <= y
                    ml_index_part{2,2}(size(ml_index_part{2,2},1)+1,1) = y;
                    ml_index_part{2,2}(size(ml_index_part{2,2},1),2) = x;
                    ml_index_full{2,2}(size(ml_index_full{2,2},1)+1,1) = y;
                    ml_index_full{2,2}(size(ml_index_full{2,2},1),2) = x;
                else
                    ml_index_full{2,2}(size(ml_index_full{2,2},1)+1,1) = x;
                    ml_index_full{2,2}(size(ml_index_full{2,2},1),2) = y;
                end
            end
        end
    end

    ml_index_part{2,1} = []; ml_index_full{2,1} = []; ch_index_full{2,1} = [];
    for y = 1:num_LC(1,1,scale.dim-1)
        if ortho
            if y <= scale.dim
                ml_index_part{2,1}(size(ml_index_part{2,1},1)+1,1) = y;
                ml_index_part{2,1}(size(ml_index_part{2,1},1),2) = 1;
                ml_index_full{2,1}(size(ml_index_full{2,1},1)+1,1) = y;
                ml_index_full{2,1}(size(ml_index_full{2,1},1),2) = 1;
                ch_index_full{2,1}(size(ch_index_full{2,1},1)+1,1) = y;
                ch_index_full{2,1}(size(ch_index_full{2,1},1),2) = 1;
            end
        else
            ml_index_part{2,1}(size(ml_index_part{2,1},1)+1,1) = y;
            ml_index_part{2,1}(size(ml_index_part{2,1},1),2) = 1;
            ml_index_full{2,1}(size(ml_index_full{2,1},1)+1,1) = y;
            ml_index_full{2,1}(size(ml_index_full{2,1},1),2) = 1;
            ch_index_full{2,1}(size(ch_index_full{2,1},1)+1,1) = y;
            ch_index_full{2,1}(size(ch_index_full{2,1},1),2) = 1;
        end
    end
    
    if strcmp(indtype,'part')
        ml_index = ml_index_part;
        ch_index = ml_index_part;
    elseif strcmp(indtype,'full')
        ml_index = ml_index_full;
        ch_index = ch_index_full;
    end

    %% Physics Model Indexing
    phys_index = [];
    if any(scale.phys == 1)
        phys_index(size(phys_index,1)+1,:) = [1, 1];
    end
    if any(scale.phys == 2)
        phys_index(size(phys_index,1)+1,:) = [2, 2];
    end
%     if all(scale.phys == [1 2])
%         phys_index(size(phys_index,1)+1,:) = [2, 1];
%     end
end