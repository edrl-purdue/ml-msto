% CHECK TO MAKE SURE HOMOGENIZED TENSORS ARE VALID
function [valid, validCH] = Check_CH(C,scale,mat,VF)
    tol_ch = eps^0.5; % add/subtract a tolerance to the Cmin and Cmax thresholds
    valid = 0;
    validCH = cell(size(C));
    for C_row = 1:size(C,1)
        for C_col = 1:size(C,2)
            if ~isempty(C{C_row,C_col})
                if C_row == 1 && C_col == 1 % mechanical CH checking
                    Cmin = mat.Emin * Prepare_C(mat.nu,1,scale.dim) - tol_ch;
                    Cmax = mat.Emax * Prepare_C(mat.nu,1,scale.dim) + tol_ch;
                    for i = 1:scale.dim
                        for j = 1:scale.dim
                            if Cmin(i,j) <= C{C_row,C_col}(i,j) && C{C_row,C_col}(i,j) <= Cmax(i,j) % check top left [dim x dim] values of CH
                                valid = 1;
                                validCH = C;
                            else
                                valid = 0;
                                if VF > 0.5
                                    validCH{C_row,C_col}(i,j) = Cmax(i,j);
                                else
                                    validCH{C_row,C_col}(i,j) = Cmin(i,j);
                                end
                                return
                            end
                        end
                    end
                    for k = (scale.dim+1):size(C{C_row,C_col},2)
                        if Cmin(k,k) <= C{C_row,C_col}(k,k) && C{C_row,C_col}(k,k) <= Cmax(k,k) % check bottom right diagonal values of CH
                            valid = 1;
                            validCH = C;
                        else
                            valid = 0;
                            if VF > 0.5
                                validCH{C_row,C_col}(k,k) = Cmax(k,k);
                            else
                                validCH{C_row,C_col}(k,k) = Cmin(k,k);
                            end
                            return
                        end
                    end
                end
                if C_row == 2 && C_col == 2 % thermal CH checking
                    Cmin = mat.kmin * Prepare_C(mat.nu,2,scale.dim) - tol_ch;
                    Cmax = mat.kmax * Prepare_C(mat.nu,2,scale.dim) + tol_ch;
                    for k = 1:scale.dim
                        if Cmin(k,k) <= C{C_row,C_col}(k,k) && C{C_row,C_col}(k,k) <= Cmax(k,k) % check diagonal values of CH
                            valid = 1;
                            validCH = C;
                        else
                            valid = 0;
                            if VF > 0.5
                                validCH{C_row,C_col}(k,k) = Cmax(k,k);
                            else
                                validCH{C_row,C_col}(k,k) = Cmin(k,k);
                            end
                            return
                        end
                    end
                end
                if C_row == 2 && C_col == 1
                    valid = 1; % do not check homogenized thermal strain for feasibility
                    
%                     alphaH = C{1,1}/C{C_row,C_col}'; % Homogenized thermal strain vector
%                     Amin = 1/(scale.dim)/mat.Amin;
%                     Amax = 1/(scale.dim)/mat.Amax;
%                     for k = 1:scale.dim
%                         if Amax <= alphaH(k) && alphaH(k) <= Amin
%                             valid = 1;
%                         else
%                             valid = 0;
%                             return
%                         end
%                     end
                end
            end
        end
    end
end