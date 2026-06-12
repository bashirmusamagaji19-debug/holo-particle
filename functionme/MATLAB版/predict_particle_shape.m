function [predLabels, confidence] = predict_particle_shape(candidateStats, modelPath, ...
    mipImage, candRoiBoxes)
% PREDICT_PARTICLE_SHAPE  Classify particles using a Top-N ensemble of
% neural networks.
%
% Inputs:
%   candidateStats  – struct array from regionprops
%   modelPath        – path to shape_classifier_model.mat (auto-detected)
%   mipImage         – MIP image for patch intensity features (optional)
%   candRoiBoxes     – N×4 ROI boxes for patch extraction (optional)
%
% Outputs:
%   predLabels  – cell array of 'circular' / 'irregular'
%   confidence  – [0,1] ensemble-averaged softmax probability

    % ---------- load model ----------
    model = [];
    pathsToTry = {};
    if nargin >= 2 && ~isempty(modelPath), pathsToTry{end+1} = modelPath; end
    scriptDir = fileparts(mfilename('fullpath'));
    if isempty(scriptDir), scriptDir = pwd; end
    pathsToTry{end+1} = fullfile(scriptDir, 'shape_classifier_model.mat');

    for p = 1:numel(pathsToTry)
        if isfile(pathsToTry{p})
            try
                loaded = load(pathsToTry{p}, 'model');
                if isfield(loaded, 'model') && isfield(loaded.model, 'nets')
                    model = loaded.model; break;
                end
            catch
            end
        end
    end

    nParticles = numel(candidateStats);
    predLabels = repmat({'irregular'}, nParticles, 1);
    confidence = zeros(nParticles, 1);
    if nParticles == 0, return; end

    if isempty(model)
        % ---- fallback: circularity threshold ----
        for k = 1:nParticles
            s = candidateStats(k);
            circ = 4 * pi * getField(s, 'Area', 0) / max(getField(s, 'Perimeter', 1)^2, eps);
            if circ >= 0.75, predLabels{k} = 'circular'; end
            confidence(k) = min(1, max(0, (circ - 0.4) / 0.5));
        end
        return;
    end

    % ---------- extract raw features ----------
    n = nParticles;
    areaPx = zeros(n,1); perimeterPx = zeros(n,1); equivDiamPx = zeros(n,1);
    solidityV = zeros(n,1); eccentricityV = zeros(n,1);
    extentV = zeros(n,1); majorAxis = zeros(n,1); minorAxis = zeros(n,1);

    for k = 1:n
        s = candidateStats(k);
        areaPx(k)       = getField(s, 'Area', 0);
        perimeterPx(k)  = getField(s, 'Perimeter', 1);
        equivDiamPx(k)  = getField(s, 'EquivDiameter', 1);
        solidityV(k)    = getField(s, 'Solidity', 0.9);
        eccentricityV(k)= getField(s, 'Eccentricity', 0);
        extentV(k)      = getField(s, 'Extent', 0.7);
        majorAxis(k)    = getField(s, 'MajorAxisLength', 1);
        minorAxis(k)    = getField(s, 'MinorAxisLength', 1);
    end

    % --- derived shape features ---
    circularity = 4 * pi * areaPx ./ max(perimeterPx.^2, eps);
    aspectRatio = majorAxis ./ max(minorAxis, eps);
    roundness   = 4 * areaPx ./ (pi * max(majorAxis.^2, eps));
    normPerimeter = perimeterPx ./ (pi * max(equivDiamPx, eps));
    convexDef = 1 - solidityV;
    ecc2 = eccentricityV .^ 2;
    csInteract = circularity .* solidityV;

    shapeFeat = [circularity, solidityV, eccentricityV, extentV, aspectRatio, ...
                 roundness, normPerimeter, convexDef, ecc2, csInteract];

    % --- MIP patch features ---
    nMip = numel(model.featureNames) - size(shapeFeat, 2);
    mipFeat = zeros(n, nMip);
    if nMip > 0 && nargin >= 4 && ~isempty(mipImage) && ~isempty(candRoiBoxes)
        [nyM, nxM] = size(mipImage);
        for k = 1:n
            roi = candRoiBoxes(k, :);
            rx = max(1, round(roi(1))); ry = max(1, round(roi(2)));
            rw = min(nxM - rx + 1, round(roi(3)));
            rh = min(nyM - ry + 1, round(roi(4)));
            if rw > 3 && rh > 3
                p = double(mipImage(ry:ry+rh-1, rx:rx+rw-1));
                pv = p(:);
                mn = min(pv); mx = max(pv);
                [gx, gy] = gradient(p);
                gm = sqrt(gx.^2 + gy.^2);
                counts = histcounts(pv, 32);
                pn = counts / sum(counts); pn(pn==0) = [];
                mipFeat(k, :) = [mean(pv), std(pv), mx, median(pv), ...
                    (mx - mn)/max(mx + mn, eps), sum(pv.^2)/numel(pv), ...
                    mSkewness(pv), -sum(pn.*log2(pn)), mean(gm(:)), mIqr(pv)];
            else
                mipFeat(k, :) = model.featureMean(end-nMip+1:end);
            end
        end
    else
        mipFeat = repmat(model.featureMean(end-nMip+1:end), n, 1);
    end

    X = [shapeFeat, mipFeat];

    % Handle NaN
    for j = 1:size(X, 2)
        bad = ~isfinite(X(:, j));
        if any(bad), X(bad, j) = model.featureMean(j); end
    end

    % Standardize
    Xs = ((X - model.featureMean) ./ max(model.featureStd, eps))';

    % ---------- ensemble predict ----------
    nKeep = model.nKeep;
    sumOut = zeros(2, n);
    for k = 1:nKeep
        sumOut = sumOut + model.nets{k}(Xs);
    end
    avgOut = sumOut / nKeep;
    [maxVal, cls] = max(avgOut, [], 1);

    for k = 1:n
        if cls(k) == 2, predLabels{k} = 'circular'; end
        confidence(k) = maxVal(k);
    end
end

function v = getField(s, fieldName, default)
    if isfield(s, fieldName) && isfinite(s.(fieldName))
        v = double(s.(fieldName));
    else
        v = default;
    end
end

function s = mSkewness(x)
    x = x(:); mu = mean(x);
    s = mean((x - mu).^3) / max(mean((x - mu).^2)^1.5, eps);
end

function q = mIqr(x)
    x = sort(x(:)); n = numel(x);
    q = x(min(n, round(n*0.75))) - x(max(1, round(n*0.25)));
end
