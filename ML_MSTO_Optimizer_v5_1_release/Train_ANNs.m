% TRAIN ANN MODELS
function [ML_model,TR] = Train_ANNs(target_path,SP_para,SP_C,micro,n_uni_train,NN,NL,reps)
    %% PREPARE FEATURES AND TARGETS
    % Features and Targets
    X_0 = SP_para;
    Y_0 = SP_C;
    % Augmentations
    X = cell(micro.N,1);
    Y = cell(micro.N,1); %Homogenized Tensor Output
    WaitMessage = parfor_waitbar(micro.N); %Add ('Waitbar', true) as arguments to turn on a waitbar popup
    fprintf('     |                    Augment TOMs                   |\n'); tic
    parfor i = 1:micro.N
        [X{i},Y{i}] = Augment_TOMs(micro,X_0(i,:),Y_0{i});
        WaitMessage.Send;
    end
    WaitMessage.Destroy;
    fprintf('     |         Done. Elapsed Time: %4.2f seconds.         |\n\n',toc); augtime = toc;
    X = cell2mat(X);
    Y = vertcat(Y{:});
    if micro.dim == 2; aug_total = 8; elseif micro.dim == 3; aug_total = 48; end % number of total augmentations

    %% NUMBER INDEX TESTING AND TRAINING DATA SETS
    n_uni_total = size(X,1)/aug_total; % Total number of samples without augmentations
    if n_uni_train > n_uni_total*0.7 % if n_uni_train is larger than 70% of n_uni_total
        n_uni_train = floor(n_uni_total*0.7); % n_uni_train should not be bigger than 70% of n_uni_total if testing IS planned
%         fprintf('WARNING: n_uni_train IS TOO LARGE. SETTING n_uni_train TO: %d\n',n_uni_train)
    end
    n_uni_test = round(n_uni_train*0.3/0.7); % Total data is 70% training data and 30% testing/validation data
    n_uni_model = n_uni_train + n_uni_test;
    total_indexes = randperm(n_uni_total,n_uni_model);
    n_ind_block_training = total_indexes(1:n_uni_train);
    n_ind_block_testing = total_indexes(n_uni_train+1:end);

    %% INDEX TRAINING
    index_training = [];
    ind1 = 1;
    for i1 = 1:n_uni_train
        if micro.n_aug == 1
            rand_samp = 1;
        elseif micro.n_aug > 1
            rand_samp = sort([1 randperm(aug_total-1,micro.n_aug-1)+1])'; % always include the first augmentation of a sample's augmentation set
        end
        index_training(ind1:ind1 + micro.n_aug-1,1) = aug_total*(n_ind_block_training(i1)-1) + rand_samp;
        ind1 = ind1 + micro.n_aug;    
    end

    %% INDEX TESTING
    index_testing = [];
    ind1 = 1;
    for i1 = 1:n_uni_test
        if micro.n_aug == 1
            rand_samp = 1;
        elseif micro.n_aug > 1
            rand_samp = sort([1 randperm(aug_total-1,micro.n_aug-1)+1])'; % always include the first augmentation of a sample's augmentation set
        end
        index_testing(ind1:(ind1 + micro.n_aug-1),1) = aug_total*(n_ind_block_testing(i1)-1) + rand_samp;
        ind1 = ind1 + micro.n_aug;
    end
    
    %% TRAIN SURROGATE MODELS
    [ML_model,TR] = ANN_model(X,Y,index_training,index_testing,micro,NN,NL,reps);
    
    %% SAVE SURROGATE MODELS
    save(target_path,'-append','ML_model','TR','augtime'); % does not appended augmented X and Y; to save space
end

%% GENERATE ANN MODELS
function [ML_model,TR] = ANN_model(X,Y,index_training,index_testing,micro,NN,NL,reps)
    % Calculate min and max constitutive tensor values
    Cmin = cell(3,3); Cmax = cell(3,3);
    if any(micro.phys == 1)
        Cmin{1,1} = micro.mat.Emin * Prepare_C(micro.mat.nu,1,micro.dim);
        Cmax{1,1} = micro.mat.Emax * Prepare_C(micro.mat.nu,1,micro.dim);
    end
    if any(micro.phys == 2)
        Cmin{2,2} = micro.mat.kmin * Prepare_C(micro.mat.nu,2,micro.dim);
        Cmax{2,2} = micro.mat.kmax * Prepare_C(micro.mat.nu,2,micro.dim);
    end
    if 2 == sum(micro.phys == [1 2])
        Cmin{2,1} = micro.mat.Amin * Prepare_C(micro.mat.nu,1,micro.dim); % might not be finished yet
        Cmax{2,1} = micro.mat.Amax * Prepare_C(micro.mat.nu,1,micro.dim); % might not be finished yet
    end
    [C_lower_HS, C_upper_HS] = Prepare_HS_bounds(micro.mat,micro.dim,micro.ortho);

    % Initalize ANN
    % https://stats.stackexchange.com/questions/181/how-to-choose-the-number-of-hidden-layers-and-nodes-in-a-feedforward-neural-netw
    % Rules of thumb for NN selection
    % NN is the average between the sizes of the input and output neurons
    % NN is between the size of input neurons and size of output neurons
    % NN is 2/3 the size of the input neurons plus the size of the output neurons
    % NN should be less than twice the size of the input neurons
    % NN should be sqrt(input neurons * output neurons) for a 3 layer NN
    hiddensizes = ones(1,NL)*NN;
    net = fitnet(hiddensizes,'trainscg'); %'trainlm' (faster for a smaller number of samples),'trainscg' (faster for a larger number of samples)
%     net.trainParam.epochs = 1e6; % max training iterations
    % https://stats.stackexchange.com/questions/181/how-to-choose-the-number-of-hidden-layers-and-nodes-in-a-feedforward-neural-netw
    % https://www.mathworks.com/help/deeplearning/ug/choose-a-multilayer-neural-network-training-function.html
    % http://www.faqs.org/faqs/ai-faq/neural-nets/part1/preamble.html
    net.trainParam.showWindow = 0;
    net.divideFcn = 'divideind';
    net.divideParam.trainInd = index_training; %70
    net.divideParam.valInd = index_testing(1:ceil(numel(index_testing)/2)); %15
    net.divideParam.testInd = index_testing(ceil(numel(index_testing)/2)+1:end); %15
    % https://www.mathworks.com/help/deeplearning/ug/divide-data-for-optimal-neural-network-training.html
    net.layers{1:end}.transferFcn = 'tansig';
    net.layers{end}.transferFcn = 'purelin'; % https://media.licdn.com/dms/image/C4D12AQGTOQOnaOcXgA/article-inline_image-shrink_400_744/0/1520042651987?e=1683158400&v=beta&t=4_Eh2uBhYFrDLR__esDKCbQ5ji6DETS-ZFKZTYcKAO4

    % Train a ANN for every output
    fprintf('     |                     Train ANNs                    |\n'); mlstart = tic;
    reps_ML_model = cell(reps,1);
    reps_TR = cell(reps,1);
    repsCorr = zeros(reps,1);
    repsMSE = zeros(reps,1);
    WaitMessage = parfor_waitbar(reps); %Add ('Waitbar', true) as arguments to turn on a waitbar popup
    [ml_index, ~, phys_index] = Indexing('part',micro,micro.ortho);
    outnum = 0; for p = 1:size(phys_index,1); for ind = 1:size(ml_index{phys_index(p,1),phys_index(p,2)},1); outnum = outnum + 1; end; end
    Yind = zeros(size(Y,1),outnum); an = 1;
    Cminvec = zeros(1,outnum); Cmaxvec = zeros(1,outnum);
    for p = 1:size(phys_index,1)
        for ind = 1:size(ml_index{phys_index(p,1),phys_index(p,2)},1)
            % Parse SP_C (Y) Data
            for y = 1:size(Y,1)
                Yind(y,an) = Y{y}{phys_index(p,1),phys_index(p,2)}(ml_index{phys_index(p,1),phys_index(p,2)}(ind,1),ml_index{phys_index(p,1),phys_index(p,2)}(ind,2));
            end
            Cminvec(1,an) = Cmin{phys_index(p,1),phys_index(p,2)}(ml_index{phys_index(p,1),phys_index(p,2)}(ind,1),ml_index{phys_index(p,1),phys_index(p,2)}(ind,2));
            Cmaxvec(1,an) = Cmax{phys_index(p,1),phys_index(p,2)}(ml_index{phys_index(p,1),phys_index(p,2)}(ind,1),ml_index{phys_index(p,1),phys_index(p,2)}(ind,2));
            an = an + 1;
        end
    end
    for n = 1:reps
        rng shuffle
        % Fit ANN model
        tStart = tic;
        [ML_model, TR] = train(net,X',Yind','useParallel','yes');
        TR.total_time = toc(tStart);
        TR.Cminvec = Cminvec;
        TR.Cmaxvec = Cmaxvec;
        
        % Correlation
        TR.train_corr = corrcoef(ML_model(X(TR.trainInd,:)')',Yind(TR.trainInd,:));
        TR.val_corr = corrcoef(ML_model(X(TR.valInd,:)')',Yind(TR.valInd,:));
        TR.test_corr = corrcoef(Cap_ANN(ML_model(X(TR.testInd,:)')',C_lower_HS(X(TR.testInd,1)'),C_upper_HS(X(TR.testInd,1)')),Yind(TR.testInd,:));
%             TR.test_corr = corrcoef(Cap_ANN(ML_model(X(TR.testInd,:)')',Cminvec,Cmaxvec),Yind(TR.testInd,:));
        repsCorr(n) = TR.test_corr(2,1);
    
        % MSE
        TR.train_MSE = TR.best_perf;
        TR.val_MSE = TR.best_vperf;
        TR.test_MSE = immse(Cap_ANN(ML_model(X(TR.testInd,:)')',C_lower_HS(X(TR.testInd,1)'),C_upper_HS(X(TR.testInd,1)')),Yind(TR.testInd,:));
%             TR.test_MSE = immse(Cap_ANN(ML_model(X(TR.testInd,:)')',Cminvec,Cmaxvec),Yind(TR.testInd,:));
        repsMSE(n) = TR.test_MSE;

        reps_ML_model{n} = ML_model;
        reps_TR{n} = TR;
        WaitMessage.Send;
    end
    corr_vecs = [1-repsCorr, transpose(1:numel(repsCorr))];
    corr_vecs = [sortrows(corr_vecs,1), transpose(1:numel(repsCorr))]; % [corr, model #, corr rank]
    corr_vecs = sortrows(corr_vecs,2);

    MSE_vecs = [repsMSE, transpose(1:numel(repsMSE))];
    MSE_vecs = [sortrows(MSE_vecs,1), transpose(1:numel(repsMSE))]; % [MSE, model #, MSE rank]
    MSE_vecs = sortrows(MSE_vecs,2);

    rank_vecs = [corr_vecs(:,3) + MSE_vecs(:,3), repsMSE, 1-repsCorr, transpose(1:numel(repsMSE))]; % [total rank, MSE, corr, model #]
    bestmodel = rank_vecs(rank_vecs(:,1)==min(rank_vecs(:,1)),:);
    bestmodel = sortrows(bestmodel,[2,3]);

    ML_model = reps_ML_model{bestmodel(1,4)};
    TR = reps_TR{bestmodel(1,4)};
    WaitMessage.Destroy;
    fprintf('     |         Done. Elapsed Time: %4.2f seconds.         |\n\n',toc(mlstart));
end