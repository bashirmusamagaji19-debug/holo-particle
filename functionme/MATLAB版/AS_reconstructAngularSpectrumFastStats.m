function summary = AS_reconstructAngularSpectrumFastStats(holoPre, zVec, lambda, pixel, useGPU, params)
% Reference implementation only. The active method-1 baseline for Python
% alignment is Fresnel_reconstructFastStats.m.
%AS_RECONSTRUCTANGULARSPECTRUMFASTSTATS Fast summary reconstruction for method-1.

    if nargin < 5 || isempty(useGPU)
        useGPU = false;
    end
    if nargin < 6 || isempty(params)
        params = struct();
    end
    params = iFillDefaultParams(params);

    if isempty(holoPre) || ~ismatrix(holoPre)
        error('AS_reconstructAngularSpectrumFastStats:InvalidInput', ...
            'holoPre must be a non-empty 2D array.');
    end
    if isempty(zVec)
        error('AS_reconstructAngularSpectrumFastStats:InvalidInput', ...
            'zVec must not be empty.');
    end

    holo = single(holoPre);
    zVec = single(zVec(:).');
    lambda = single(lambda);
    pixel = single(pixel);

    [ny, nx] = size(holo);
    nz = numel(zVec);
    requestedGPU = logical(useGPU);
    canUseGPU = requestedGPU && iHasUsableGpu();

    tAll = tic;
    if canUseGPU
        try
            [summary, stat] = iFastStatsCore(holo, zVec, lambda, pixel, params, true);
            modeMsg = 'GPU';
        catch ME
            warning('AS_reconstructAngularSpectrumFastStats:GPUFallback', ...
                'Fast stats GPU path failed, falling back to CPU. Reason: %s', ME.message);
            [summary, stat] = iFastStatsCore(holo, zVec, lambda, pixel, params, false);
            modeMsg = 'CPU (GPU fallback)';
        end
    else
        if requestedGPU
            fprintf('[FastStats] GPU unavailable, using CPU.\n');
        end
        [summary, stat] = iFastStatsCore(holo, zVec, lambda, pixel, params, false);
        modeMsg = 'CPU';
    end

    summary.mode = modeMsg;
    summary.volumeSize = [ny, nx, nz];
    summary.totalSec = toc(tAll);

    fprintf(['[FastStats][%s] upload %.3fs, grid %.3fs, fft %.3fs, mip %.3fs, ' ...
             'candidate %.3fs, curve %.3fs, gather %.3fs, batch %d\n'], ...
        modeMsg, stat.uploadSec, stat.gridSec, stat.fftSec, stat.mipSec, ...
        stat.candidateSec, stat.curveSec, stat.gatherSec, stat.batchSize);
    fprintf('[FastStats] summary complete. candidates %d, valid %d, total %.3fs\n', ...
        stat.candidateCount, size(summary.coords3D, 1), summary.totalSec);
end

function [summary, stat] = iFastStatsCore(holo, zVec, lambda, pixel, params, useGPU)
    [ny, nx] = size(holo);
    nz = numel(zVec);
    stat = iInitStat();

    tGrid = tic;
    [FX, FY] = iFrequencyGrid(nx, ny, pixel, 'single');
    lambda2 = lambda * lambda;
    if useGPU
        holoW = gpuArray(holo);
        zVecW = gpuArray(zVec);
        FX = gpuArray(FX);
        FY = gpuArray(FY);
        stat.uploadSec = toc(tGrid);
    else
        holoW = holo;
        zVecW = zVec;
    end
    radicand = single(1) - lambda2 .* (FX.^2 + FY.^2);
    passband = single(radicand >= 0);
    kz = sqrt(max(radicand, single(0)));
    phaseBase = complex(single(0), single(2 * pi / lambda)) .* kz;
    stat.gridSec = toc(tGrid);

    tFft = tic;
    U0f = fft2(holoW);
    stat.fftSec = toc(tFft);

    tMip = tic;
    if useGPU
        mipW = gpuArray.zeros(ny, nx, 'single');
        argmaxZW = gpuArray.ones(ny, nx, 'uint16');
        sumW = gpuArray.zeros(ny, nx, 'single');
    else
        mipW = zeros(ny, nx, 'single');
        argmaxZW = ones(ny, nx, 'uint16');
        sumW = zeros(ny, nx, 'single');
    end

    for k = 1:nz
        H = exp(phaseBase .* zVecW(k)) .* passband;
        I = abs(ifft2(U0f .* H)).^2;
        fig = figure;
        for i=1:size(I,3)

            imshow(max(max(I(:,:,i)))-I(:,:,i),[])
            % imshow(I(:,:,i),[])
            pause(0.01)
        end
        close(fig)

        sumW = sumW + I;
        updateMask = I > mipW;
        mipW(updateMask) = I(updateMask);
        argmaxZW(updateMask) = uint16(k);
    end
    stat.mipSec = toc(tMip);
    stat.batchSize = 1;

    tg = tic;
    if useGPU
        mip2D = gather(mipW);
        argmaxZMap = gather(argmaxZW);
        mean2D = gather(sumW ./ max(single(nz), 1));
    else
        mip2D = mipW;
        argmaxZMap = argmaxZW;
        mean2D = sumW ./ max(single(nz), 1);
    end
    stat.gatherSec = toc(tg);

    roiParams = iBuildAdaptiveRoiParams(params, pixel, zVec);
    tCand = tic;
    candidateData = iDetectCandidates(mip2D, roiParams);
    stat.candidateSec = toc(tCand);
    stat.candidateCount = numel(candidateData.candidates);

    if isempty(candidateData.candidates)
        summary = iEmptySummary(candidateData, zVec, [ny, nx, nz]);
        return;
    end

    tCurve = tic;
    candidateData.maxMap = mip2D;
    candidateData.meanMap = mean2D;
    candidateData.argmaxZMap = argmaxZMap;
    [focusIdx, peakRatioApprox] = iEstimateFocusFromMaps(candidateData);
    bestPatches = iCollectBestPatches(U0f, phaseBase, passband, zVecW, ...
        candidateData.candidates, focusIdx, useGPU);
    stat.curveSec = toc(tCurve);

    summary = iFinalizeSummaryFromMaps(candidateData, focusIdx, peakRatioApprox, ...
        bestPatches, zVec, params, [ny, nx, nz]);
end

function bestPatches = iCollectBestPatches(U0f, phaseBase, passband, zVecW, candidates, focusIdx, useGPU)
    numCand = numel(candidates);
    bestPatches = cell(numCand, 1);
    uniqueIdx = unique(double(focusIdx(:)).');

    for zIdx = uniqueIdx
        if zIdx < 1
            continue;
        end

        H = exp(phaseBase .* zVecW(zIdx)) .* passband;
        I = abs(ifft2(U0f .* H)).^2;
        if useGPU
            I = gather(I);
        end

        pickIdx = find(double(focusIdx(:)) == zIdx);
        for k = 1:numel(pickIdx)
            n = pickIdx(k);
            c = candidates(n);
            bestPatches{n} = I(c.y1:c.y2, c.x1:c.x2);
        end
    end
end

function [focusIdx, peakRatioApprox] = iEstimateFocusFromMaps(candidateData)
    numCand = numel(candidateData.candidates);
    focusIdx = ones(numCand, 1, 'uint32');
    peakRatioApprox = ones(numCand, 1);

    for n = 1:numCand
        c = candidateData.candidates(n);
        zPatch = double(candidateData.argmaxZMap(c.y1:c.y2, c.x1:c.x2));
        wPatch = double(candidateData.maxMap(c.y1:c.y2, c.x1:c.x2));
        meanPatch = double(candidateData.meanMap(c.y1:c.y2, c.x1:c.x2));

        zVals = zPatch(:);
        wVals = max(wPatch(:), 0);
        valid = isfinite(zVals) & (zVals >= 1) & isfinite(wVals);
        zVals = zVals(valid);
        wVals = wVals(valid);

        if isempty(zVals)
            focusIdx(n) = uint32(1);
        else
            uniqueZ = unique(zVals);
            weightPerZ = zeros(size(uniqueZ));
            for k = 1:numel(uniqueZ)
                weightPerZ(k) = sum(wVals(zVals == uniqueZ(k)));
            end
            [~, bestIdx] = max(weightPerZ);
            focusIdx(n) = uint32(uniqueZ(bestIdx));
        end

        peakVal = mean(wPatch(:));
        meanVal = mean(meanPatch(:));
        peakRatioApprox(n) = peakVal / max(meanVal, eps);
    end
end

function summary = iFinalizeSummaryFromMaps(candidateData, focusIdx, peakRatioApprox, bestPatches, zVec, params, volumeSize)
    zVecUm = double(zVec(:)) * 1e6;
    numCand = numel(candidateData.candidates);
    coords3D = zeros(numCand, 3);
    axialCurves = cell(numCand, 1);
    validMask = false(numCand, 1);
    stats = candidateData.stats;

    for n = 1:numCand
        c = candidateData.candidates(n);
        maxIdx = min(numel(zVecUm), max(1, double(focusIdx(n))));
        peakRatio = peakRatioApprox(n);
        isAwayFromEdge = (maxIdx > params.edgeMargin) && (maxIdx <= numel(zVecUm) - params.edgeMargin);

        measureInfo = AS_refineMeasurementROI(bestPatches{n}, [c.cxLocal, c.cyLocal], ...
            [1, 1, size(bestPatches{n}, 2), size(bestPatches{n}, 1)], params);
        measureBox = measureInfo.measureBoxLocal;
        measurePatch = bestPatches{n}(measureBox(2):measureBox(2)+measureBox(4)-1, ...
            measureBox(1):measureBox(1)+measureBox(3)-1);
        measureCenter = measureInfo.refinedCenterLocal - [measureBox(1), measureBox(2)] + 1;

        focusDiamPx = iEstimateFocusDiameterFromPatch(measurePatch, ...
            round(measureCenter(1)), round(measureCenter(2)), ...
            params.xyThreshold, stats(n).EquivDiameter);

        stats(n).EquivDiameter = focusDiamPx;
        coords3D(n, :) = [c.centerXY(1), c.centerXY(2), zVecUm(maxIdx)];
        axialCurves{n} = zVecUm(maxIdx);
        validMask(n) = isAwayFromEdge && (peakRatio >= params.minPeakRatio);
    end

    summary = struct();
    summary.img2D = candidateData.mip2D;
    summary.mip2D = candidateData.mip2D;
    summary.bw = candidateData.bw;
    summary.zVecUm = zVecUm;
    summary.volumeSize = volumeSize;
    summary.hasVolume = false;
    summary.validMask = validMask;
    summary.candidateCoords3D = coords3D;
    summary.candidateRoiBoxes = candidateData.roiBoxes;
    summary.candidateStats = stats;
    summary.coords3D = coords3D(validMask, :);
    summary.axialCurves = axialCurves(validMask);
    summary.roiBoxes = candidateData.roiBoxes(validMask, :);
    summary.stats = stats(validMask);
    if isempty(summary.stats)
        summary.diamsPx = zeros(0, 1);
    else
        summary.diamsPx = reshape([summary.stats.EquivDiameter], [], 1);
    end
end

function out = iEmptySummary(candidateData, zVec, volumeSize)
    out = struct();
    out.img2D = candidateData.mip2D;
    out.mip2D = candidateData.mip2D;
    out.bw = candidateData.bw;
    out.zVecUm = double(zVec(:)) * 1e6;
    out.volumeSize = volumeSize;
    out.hasVolume = false;
    out.validMask = false(0, 1);
    out.candidateCoords3D = zeros(0, 3);
    out.candidateRoiBoxes = zeros(0, 4);
    out.candidateStats = struct([]);
    out.coords3D = zeros(0, 3);
    out.axialCurves = {};
    out.roiBoxes = zeros(0, 4);
    out.stats = struct([]);
    out.diamsPx = zeros(0, 1);
end

function out = iDetectCandidates(mip2D, params)
    levelRatio = min(max(params.xyThreshold, 0), 1);
    mipPeak = max(mip2D(:));
    levelAbs = levelRatio * mipPeak;
    bw = mip2D > levelAbs;
    bw = bwareaopen(bw, 1);
    bw = imfill(bw, 'holes');

    stats = regionprops(bw, mip2D, 'Centroid', 'WeightedCentroid', ...
        'EquivDiameter', 'BoundingBox');

    out = struct();
    out.img2D = mip2D;
    out.mip2D = mip2D;
    out.bw = bw;
    out.stats = stats;
    out.roiBoxes = zeros(numel(stats), 4);
    out.candidates = repmat(struct( ...
        'centerXY', [], ...
        'x1', 1, ...
        'x2', 1, ...
        'y1', 1, ...
        'y2', 1, ...
        'roiBox', [1 1 1 1], ...
        'cxLocal', 1, ...
        'cyLocal', 1), numel(stats), 1);

    for n = 1:numel(stats)
        centerXY = stats(n).WeightedCentroid;
        if any(~isfinite(centerXY))
            centerXY = stats(n).Centroid;
        end

        cx = round(centerXY(1));
        cy = round(centerXY(2));
        candidate = struct( ...
            'centerXY', centerXY, ...
            'equivDiameterPx', stats(n).EquivDiameter, ...
            'bboxSizePx', [stats(n).BoundingBox(3), stats(n).BoundingBox(4)]);
        roiInfo = AS_buildAdaptiveROI(mip2D, candidate, params);

        x1 = roiInfo.searchBox(1);
        y1 = roiInfo.searchBox(2);
        x2 = x1 + roiInfo.searchBox(3) - 1;
        y2 = y1 + roiInfo.searchBox(4) - 1;

        out.roiBoxes(n, :) = roiInfo.searchBox;
        out.candidates(n).centerXY = centerXY;
        out.candidates(n).x1 = x1;
        out.candidates(n).x2 = x2;
        out.candidates(n).y1 = y1;
        out.candidates(n).y2 = y2;
        out.candidates(n).roiBox = out.roiBoxes(n, :);
        out.candidates(n).cxLocal = min(x2 - x1 + 1, max(1, cx - x1 + 1));
        out.candidates(n).cyLocal = min(y2 - y1 + 1, max(1, cy - y1 + 1));
    end
end

function roiParams = iBuildAdaptiveRoiParams(params, pixel, zVec)
    roiParams = params;
    roiParams.pix_um = double(pixel) * 1e6;

    if numel(zVec) > 1
        roiParams.dz = double(median(abs(diff(zVec))));
    else
        roiParams.dz = 0;
    end

    if ~isfield(roiParams, 'measureRoiEnergyFraction') || ~isscalar(roiParams.measureRoiEnergyFraction) || ~isfinite(roiParams.measureRoiEnergyFraction)
        roiParams.measureRoiEnergyFraction = 0.92;
    end
    if ~isfield(roiParams, 'searchRoiGrowFactor') || ~isscalar(roiParams.searchRoiGrowFactor) || ~isfinite(roiParams.searchRoiGrowFactor)
        roiParams.searchRoiGrowFactor = 1.25;
    end
    if ~isfield(roiParams, 'edgeEnergyThreshold') || ~isscalar(roiParams.edgeEnergyThreshold) || ~isfinite(roiParams.edgeEnergyThreshold)
        roiParams.edgeEnergyThreshold = 0.20;
    end
    if ~isfield(roiParams, 'searchMinRadiusUm') || ~isscalar(roiParams.searchMinRadiusUm) || ~isfinite(roiParams.searchMinRadiusUm)
        roiParams.searchMinRadiusUm = 30;
    end
    if ~isfield(roiParams, 'roiMaxRadiusPx') || ~isscalar(roiParams.roiMaxRadiusPx) || ~isfinite(roiParams.roiMaxRadiusPx)
        roiParams.roiMaxRadiusPx = max(roiParams.roiMinRadius + 4, ceil(250 / max(roiParams.pix_um, eps)));
    end
end

function focusDiamPx = iEstimateFocusDiameterFromPatch(roiFocus, cxLocal, cyLocal, xyThreshold, fallbackDiamPx)
    if isempty(roiFocus)
        focusDiamPx = fallbackDiamPx;
        return;
    end

    roiFocus = AS_normalizeImage(roiFocus);
    if ~any(roiFocus(:) > 0)
        focusDiamPx = fallbackDiamPx;
        return;
    end

    if exist('graythresh', 'file') == 2
        localLevel = graythresh(roiFocus);
    else
        localLevel = xyThreshold;
    end
    localLevel = max(0.35, min(0.85, max(localLevel, 0.6 * xyThreshold)));

    bwFocus = imbinarize(roiFocus, localLevel);
    bwFocus = bwareaopen(bwFocus, 1);
    bwFocus = imfill(bwFocus, 'holes');
    if ~any(bwFocus(:))
        focusDiamPx = fallbackDiamPx;
        return;
    end

    cc = bwconncomp(bwFocus, 8);
    localStats = regionprops(cc, roiFocus, 'EquivDiameter', 'WeightedCentroid', 'Centroid');
    cxLocal = min(size(bwFocus, 2), max(1, round(cxLocal)));
    cyLocal = min(size(bwFocus, 1), max(1, round(cyLocal)));
    centerIdx = sub2ind(size(bwFocus), cyLocal, cxLocal);

    pickIdx = [];
    for k = 1:cc.NumObjects
        if any(cc.PixelIdxList{k} == centerIdx)
            pickIdx = k;
            break;
        end
    end

    if isempty(pickIdx)
        centers = reshape([localStats.WeightedCentroid], 2, []).';
        if isempty(centers) || any(~isfinite(centers(:)))
            centers = reshape([localStats.Centroid], 2, []).';
        end
        dist2 = (centers(:, 1) - cxLocal).^2 + (centers(:, 2) - cyLocal).^2;
        [~, pickIdx] = min(dist2);
    end

    focusDiamPx = localStats(pickIdx).EquivDiameter;
    if ~(isfinite(focusDiamPx) && focusDiamPx > 0)
        focusDiamPx = fallbackDiamPx;
    end
end

function [FX, FY] = iFrequencyGrid(nx, ny, pixel, dtype)
    fx = ifftshift((-floor(nx/2):ceil(nx/2)-1) ./ (double(nx) * double(pixel)));
    fy = ifftshift((-floor(ny/2):ceil(ny/2)-1) ./ (double(ny) * double(pixel)));
    [FX, FY] = meshgrid(cast(fx, dtype), cast(fy, dtype));
end

function tf = iHasUsableGpu()
    tf = false;
    try
        if exist('gpuDeviceCount', 'file') == 2
            tf = gpuDeviceCount > 0;
        end
    catch
        tf = false;
    end
end

function out = iInitStat()
    out = struct('uploadSec', 0, 'gridSec', 0, 'fftSec', 0, ...
        'mipSec', 0, 'candidateSec', 0, 'curveSec', 0, ...
        'gatherSec', 0, 'batchSize', 1, 'candidateCount', 0);
end

function params = iFillDefaultParams(params)
    if ~isfield(params, 'minPeakRatio') || ~isfinite(params.minPeakRatio)
        params.minPeakRatio = 1.15;
    end
    if ~isfield(params, 'edgeMargin') || ~isfinite(params.edgeMargin)
        params.edgeMargin = 2;
    end
    if ~isfield(params, 'xyThreshold') || ~isfinite(params.xyThreshold)
        params.xyThreshold = 0.85;
    end
    if ~isfield(params, 'roiScale') || ~isfinite(params.roiScale)
        params.roiScale = 0.50;
    end
    if ~isfield(params, 'roiMinRadius') || ~isfinite(params.roiMinRadius)
        params.roiMinRadius = 15;
    end

    params.minPeakRatio = max(params.minPeakRatio, 1);
    params.edgeMargin = max(0, round(params.edgeMargin));
    params.xyThreshold = min(max(params.xyThreshold, 0), 1);
    params.roiScale = max(0.1, params.roiScale);
    params.roiMinRadius = max(1, round(params.roiMinRadius));
end
