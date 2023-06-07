% CAP ANN PREDICTIONS WITH H-S BOUNDS
function pred_out = Cap_ANN(pred_in,Cminvec,Cmaxvec)
    % pred_in [rows = number of predictions, cols = number of C properties]
    % for 1 prediction, pred_in is a row vector
    if size(Cminvec,1)==1 && size(Cmaxvec,1)==1
        pred_out = arrayfun(@min,arrayfun(@max,pred_in,repmat(Cminvec,[size(pred_in,1),1])),repmat(Cmaxvec,[size(pred_in,1),1]));
    else
        pred_out = arrayfun(@min,arrayfun(@max,pred_in,Cminvec),Cmaxvec);
    end
    
end