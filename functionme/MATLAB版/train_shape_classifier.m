function model = train_shape_classifier(featureCsv, modelOutput)
% TRAIN_SHAPE_CLASSIFIER  Train an ensemble of neural networks to classify
% circular vs irregular particles.
%
% Uses patternnet (Deep Learning Toolbox).  Top-5 nets ensembled.
% Bad datasets (match rate < 60%) are excluded from training.
%
% Usage:
%   model = train_shape_classifier;
%   model = train_shape_classifier('features_table.csv', 'shape_classifier_model.mat');

    if nargin < 1 || isempty(featureCsv)
        scriptDir = fileparts(mfilename('fullpath'));
        if isempty(scriptDir), scriptDir = pwd; end
        featureCsv = fullfile(scriptDir, 'features_table.csv');
    end
    if nargin < 2 || isempty(modelOutput)
        scriptDir = fileparts(mfilename('fullpath'));
        if isempty(scriptDir), scriptDir = pwd; end
        modelOutput = fullfile(scriptDir, 'shape_classifier_model.mat');
    end

    % ---------- load + clean ----------
    if ~isfile(featureCsv)
        error('Feature CSV not found: %s', featureCsv);
    end
    T = readtable(featureCsv);
    fprintf('Loaded %d rows\n', height(T));

    % Keep only matched particles
    T = T(T.match_valid == 1, :);

    % Remove bad datasets (< 12 matched = < 60% match rate)
    [G, dsIds] = findgroups(T.dataset_id);
    keepDs = true(height(T), 1);
    for g = 1:numel(dsIds)
        nDs = sum(G == g);
        if nDs < 12
            keepDs(G == g) = false;
            fprintf('Excluding %s: only %d matched\n', char(dsIds(g)), nDs);
        end
    end
    nRemoved = sum(~keepDs);
    if nRemoved > 0
        fprintf('Removed %d particles from %d bad datasets\n', nRemoved, ...
            numel(unique(T.dataset_id(~keepDs))));
    end
    T = T(keepDs, :);
    fprintf('Using %d matched particles (%d datasets)\n', height(T), ...
        numel(unique(T.dataset_id)));

    % ---------- feature engineering ----------
    areaPx       = T.area_px;
    perimeterPx  = T.perimeter_px;
    equivDiamPx  = T.equiv_diam_px;
    solidity     = T.solidity;
    eccentricity = T.eccentricity;
    extent       = T.extent;
    majorAxis    = T.major_axis_px;
    minorAxis    = T.minor_axis_px;
    circularity  = T.circularity;
    aspectRatio  = T.aspect_ratio;

    roundness    = 4 * areaPx ./ (pi * max(majorAxis.^2, eps));
    normPerimeter = perimeterPx ./ (pi * max(equivDiamPx, eps));
    convexDeficiency = 1 - solidity;
    eccentricity2 = eccentricity .^ 2;
    circ_solid_interact = circularity .* solidity;

    % MIP features (use mean imputation if NaN)
    mipFields = {'mip_mean','mip_std','mip_max','mip_median','mip_contrast','mip_energy', ...
                 'mip_skewness','mip_entropy','mip_gradient','mip_iqr'};

    shapeFeat = [circularity, solidity, eccentricity, extent, aspectRatio, ...
                 roundness, normPerimeter, convexDeficiency, eccentricity2, ...
                 circ_solid_interact];

    nShape = size(shapeFeat, 2);
    nMip = numel(mipFields);
    mipFeat = zeros(height(T), nMip);
    for j = 1:nMip
        if ismember(mipFields{j}, T.Properties.VariableNames)
            mipFeat(:, j) = T.(mipFields{j});
        end
    end
    % Impute NaN MIP features with column mean
    for j = 1:nMip
        bad = ~isfinite(mipFeat(:, j));
        if any(bad), mipFeat(bad, j) = nanmean(mipFeat(:, j)); end
    end

    X = [shapeFeat, mipFeat];

    featureNames = [{'circularity','solidity','eccentricity','extent','aspect_ratio', ...
        'roundness','norm_perimeter','convex_deficiency','eccentricity2','circ_solid_interact'}, ...
        mipFields];

    fprintf('Features (%d): %s\n', numel(featureNames), strjoin(featureNames, ', '));

    % Labels: circular = 1, irregular = 0
    Y = double(strcmp(T.shape_label, 'circular'));
    fprintf('Samples: %d circular, %d irregular\n', sum(Y==1), sum(Y==0));

    % Sample weights: closer match → higher confidence → higher weight
    weights = 1 ./ max(T.match_distance_px, 1);
    weights = weights / mean(weights);

    % ---------- clean ----------
    badRows = any(~isfinite(X), 2);
    if any(badRows)
        fprintf('Removing %d NaN rows\n', sum(badRows));
        X(badRows, :) = []; Y(badRows) = []; weights(badRows) = [];
    end
    N = numel(Y);

    % ---------- standardize ----------
    featureMean = mean(X, 1);
    featureStd  = std(X, 0, 1);
    featureStd(featureStd < eps) = 1;
    Xs = (X - featureMean) ./ featureStd;

    % ================================================================
    %  Architecture selection
    % ================================================================
    fprintf('\n--- Architecture selection (5-fold CV) ---\n');
    rng(42);
    foldId = mod(randperm(N), 5) + 1;

    hiddenSizes = [10, 14, 18, 22, 26];
    bestH = 14; bestCV = 0;

    for h = hiddenSizes
        foldAcc = zeros(5, 1);
        for fold = 1:5
            teIdx = (foldId == fold); trIdx = ~teIdx;
            Xtr = Xs(trIdx, :)'; Ytr = full(ind2vec(Y(trIdx)' + 1));
            Xte = Xs(teIdx, :)';

            net = patternnet(h);
            net.trainParam.showWindow = false;
            net.trainParam.showCommandLine = false;
            net = train(net, Xtr, Ytr);
            pred = net(Xte);
            [~, cls] = max(pred, [], 1);
            foldAcc(fold) = mean((cls - 1)' == Y(teIdx));
        end
        fprintf('  Hidden=%2d  CV=%.1f%%\n', h, mean(foldAcc)*100);
        if mean(foldAcc) > bestCV, bestCV = mean(foldAcc); bestH = h; end
    end
    fprintf('Best: %d hidden (CV=%.1f%%)\n', bestH, bestCV*100);

    % ================================================================
    %  Train ensemble of 20 nets, keep top 5
    % ================================================================
    fprintf('\n--- Training ensemble (20 nets, keeping top 5) ---\n');
    nTrain = 20; nKeep = 5;
    nets = cell(nTrain, 1);
    valAcc = zeros(nTrain, 1);

    for i = 1:nTrain
        holdIdx = randperm(N, round(N * 0.8));
        valIdx  = setdiff(1:N, holdIdx);
        Xtr = Xs(holdIdx, :)';
        Ytr = full(ind2vec(Y(holdIdx)' + 1));
        Xval = Xs(valIdx, :)';
        Yval = Y(valIdx);

        net = patternnet(bestH);
        net.trainParam.showWindow = false;
        net.trainParam.showCommandLine = false;
        net.divideFcn = 'divideind';
        net.divideParam.trainInd = 1:size(Xtr, 2);
        net.divideParam.valInd   = [];
        net.divideParam.testInd  = [];
        net = train(net, Xtr, Ytr);

        pred = net(Xval);
        [~, cls] = max(pred, [], 1);
        valAcc(i) = mean((cls - 1)' == Yval);
        nets{i} = net;
    end

    [sortedAcc, sortIdx] = sort(valAcc, 'descend');
    topNets = nets(sortIdx(1:nKeep));
    fprintf('Top-5 val accuracies: %.1f%% %s\n', mean(sortedAcc(1:nKeep))*100, ...
        sprintf('%.1f%% ', sortedAcc(1:nKeep)*100));

    % ================================================================
    %  Ensemble CV evaluation
    % ================================================================
    fprintf('\n--- Ensemble 5-fold CV ---\n');
    cvPred = zeros(N, 1);
    cvConf = zeros(N, 1);
    evalFoldId = mod(randperm(N), 5) + 1;

    for fold = 1:5
        teIdx = (evalFoldId == fold); trIdx = ~teIdx;
        Xtr = Xs(trIdx, :)'; Ytr = full(ind2vec(Y(trIdx)' + 1));
        Xte = Xs(teIdx, :)';

        % Train & ensemble
        foldNets = cell(nKeep, 1);
        for k = 1:nKeep
            net = patternnet(bestH);
            net.trainParam.showWindow = false;
            net.trainParam.showCommandLine = false;
            net.divideFcn = 'divideind';
            net.divideParam.trainInd = 1:size(Xtr, 2);
            net.divideParam.valInd   = [];
            net.divideParam.testInd  = [];
            foldNets{k} = train(net, Xtr, Ytr);
        end

        % Ensemble prediction: average softmax outputs
        sumOut = zeros(2, sum(teIdx));
        for k = 1:nKeep
            sumOut = sumOut + foldNets{k}(Xte);
        end
        avgOut = sumOut / nKeep;
        [maxVal, cls] = max(avgOut, [], 1);
        cvPred(teIdx) = (cls - 1)';
        cvConf(teIdx) = maxVal';
    end

    cvAccuracy = mean(cvPred == Y);
    confMat = confusionmat(Y, cvPred);
    tp = confMat(2,2); tn = confMat(1,1);
    fp = confMat(1,2); fn = confMat(2,1);
    precision = tp / max(tp + fp, 1);
    recall    = tp / max(tp + fn, 1);
    f1 = 2 * precision * recall / max(precision + recall, 1e-9);

    fprintf('\n=== Final Results (Top-%d Ensemble) ===\n', nKeep);
    fprintf('CV Accuracy:  %.1f%%\n', cvAccuracy*100);
    fprintf('              PredCirc  PredIrr\n');
    fprintf('TrueCirc         %3d       %3d\n', tp, fp);
    fprintf('TrueIrr          %3d       %3d\n', fn, tn);
    fprintf('Precision: %.1f%%, Recall: %.1f%%, F1: %.3f\n', precision*100, recall*100, f1);

    % ---------- feature importance ----------
    imp = zeros(1, numel(featureNames));
    for j = 1:numel(featureNames)
        Xperm = Xs;
        Xperm(:, j) = Xperm(randperm(N), j);
        sumOut = zeros(2, N);
        for k = 1:nKeep
            sumOut = sumOut + topNets{k}(Xperm');
        end
        [~, cls] = max(sumOut, [], 1);
        imp(j) = cvAccuracy - mean((cls - 1)' == Y);
    end
    impZ = (imp - mean(imp)) / max(std(imp), 1e-9);
    fprintf('\nFeature Importance (Z-score):\n');
    for j = 1:numel(featureNames)
        fprintf('  %-22s  %+.2f\n', featureNames{j}, impZ(j));
    end

    % ---------- save model ----------
    model = struct();
    model.nets               = topNets;
    model.nKeep              = nKeep;
    model.nHidden            = bestH;
    model.featureNames       = featureNames;
    model.featureMean        = featureMean;
    model.featureStd         = featureStd;
    model.cvAccuracy         = cvAccuracy;
    model.confusionMat       = confMat;
    model.oobPrecision       = precision;
    model.oobRecall          = recall;
    model.oobF1              = f1;
    model.featureImportanceZ = impZ;
    model.trainingDate       = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    model.nSamples           = N;
    model.nCircular          = sum(Y == 1);
    model.nIrregular         = sum(Y == 0);

    save(modelOutput, 'model');
    fprintf('\nModel saved to: %s\n', modelOutput);

    % ---------- visualization ----------
    figure('Name', 'Shape Classifier (Ensemble)', 'NumberTitle', 'off', ...
        'Position', [100 100 1000 360]);

    subplot(1, 3, 1);
    barh(flip(impZ));
    set(gca, 'YTickLabel', flip(featureNames));
    xlabel('Importance (Z)');
    title(sprintf('Feature Importance  (CV=%.1f%%)', cvAccuracy*100));
    grid on;

    subplot(1, 3, 2);
    imagesc(confMat); colormap(flipud(bone)); colorbar;
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Pred Irr', 'Pred Circ'}, ...
             'YTick', [1 2], 'YTickLabel', {'True Irr', 'True Circ'});
    title(sprintf('P=%.1f%% R=%.1f%% F1=%.3f', precision*100, recall*100, f1));
    for r = 1:2, for cl = 1:2
        text(cl, r, num2str(confMat(r,cl)), 'HorizontalAlign', 'center', ...
            'FontSize', 16, 'FontWeight', 'bold', 'Color', [1 0.4 0]);
    end, end

    subplot(1, 3, 3);
    hold on;
    plot(X(Y==0, 2), X(Y==0, 1), 'r.', 'MarkerSize', 8);
    plot(X(Y==1, 2), X(Y==1, 1), 'b.', 'MarkerSize', 8);
    xlabel('Solidity'); ylabel('Circularity');
    title(sprintf('Feature Space (%d samples)', N));
    legend('Irregular', 'Circular', 'Location', 'best'); grid on;
    drawnow;
end
