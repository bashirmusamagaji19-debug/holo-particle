% Active method-1 baseline kept for Python alignment.
function summary = Fresnel_reconstructFastStats(holoPre, zVec, lambda, pixel, useGPU, params)
%FRESNEL_RECONSTRUCTFASTSTATS Fast particle-statistics summary (Fresnel Method).

    if nargin < 5 || isempty(useGPU)
        useGPU = false;
    end
    if nargin < 6 || isempty(params)
        params = struct();
    end
    params = iFillDefaultParams(params);

    if isempty(holoPre) || ~ismatrix(holoPre)
        error('Fresnel_reconstructFastStats:InvalidInput', ...
            'holoPre must be a non-empty 2D array.');
    end
    if isempty(zVec)
        error('Fresnel_reconstructFastStats:InvalidInput', ...
            'zVec must not be empty.');
    end

    zVec = double(zVec(:).');
    lambda = double(lambda);
    pixel = double(pixel);
    holo = single(holoPre);

    [ny, nx] = size(holo);
    nz = numel(zVec);
    requestedGPU = logical(useGPU);
    canUseGPU = requestedGPU && iHasUsableGpu();

    tAll = tic;
    if canUseGPU
        try
            [summary, stat] = iFastStatsGpuCore(holo, zVec, lambda, pixel, params);
            modeMsg = 'GPU（混合精度）';
        catch ME
            warning('Fresnel_reconstructFastStats:GPUFallback', ...
                'Fresnel fast stats GPU path failed, falling back to CPU. Reason: %s', ME.message);
            [summary, stat] = iFastStatsCpuCore(holo, zVec, lambda, pixel, params);
            modeMsg = 'CPU（回退）';
        end
    else
        if requestedGPU
            fprintf('[菲涅尔快速统计] GPU 不可用，切换为 CPU。\n');
        end
        [summary, stat] = iFastStatsCpuCore(holo, zVec, lambda, pixel, params);
        modeMsg = 'CPU（混合精度）';
    end

    summary.mode = modeMsg;
    summary.volumeSize = [ny, nx, nz];
    summary.totalSec = toc(tAll);

    fprintf(['[菲涅尔快速统计][%s] 频率网格 %.3f 秒，FFT %.3f 秒，MIP %.3f 秒，' ...
             '候选检测 %.3f 秒，轴向曲线 %.3f 秒，结果回传 %.3f 秒，批大小 %d\n'], ...
        modeMsg, stat.gridSec, stat.fftSec, stat.mipSec, ...
        stat.candidateSec, stat.curveSec, stat.gatherSec, stat.batchSize);
    fprintf('[菲涅尔快速统计] 汇总完成。候选 %d 个，有效 %d 个，总耗时 %.3f 秒\n', ...
        stat.candidateCount, size(summary.coords3D, 1), summary.totalSec);
end

function [summary, stat] = iFastStatsGpuCore(holo, zVec, lambda, pixel, params)
    [ny, nx] = size(holo);
    nz = numel(zVec);
    stat = iInitStat();
    g = gpuDevice;

    tGrid = tic;
    [FX, FY] = iFrequencyGrid(nx, ny, pixel, 'double');
    FXG = gpuArray(FX);
    FYG = gpuArray(FY);
    phaseBase = (-1i * pi * lambda) .* (FXG.^2 + FYG.^2);
    k_offset = (2 * pi / lambda);
    stat.gridSec = toc(tGrid);

    tFft = tic;
    U0f = fft2(gpuArray(holo));
    stat.fftSec = toc(tFft);

    batchSize = iChooseBatchSizeGPU(nx, ny, nz, g.AvailableMemory);
    stat.batchSize = batchSize;

    mipG = gpuArray.zeros(ny, nx, 'single');
    argmaxZG = gpuArray.ones(ny, nx, 'uint16');
    sumG = gpuArray.zeros(ny, nx, 'single');

    tMip = tic;
    for s = 1:batchSize:nz
        idx = s:min(s + batchSize - 1, nz);
        zBatch = reshape(zVec(idx), 1, 1, []);
        H = exp(1i * k_offset * zBatch + phaseBase .* zBatch);
        Amp = abs(ifft2(U0f .* single(H)));
        sliceMax = max(max(Amp, [], 1), [], 2);
        if params.invertContrast
            I = bsxfun(@minus, sliceMax, Amp);
        else
            I = Amp;
        end

        sumG = sumG + sum(I, 3);
        [batchMax, localArgmax] = max(I, [], 3);
        updateMask = batchMax > mipG;
        mipG(updateMask) = batchMax(updateMask);
        argLocal = uint16(idx(1) - 1) + uint16(localArgmax);
        argmaxZG(updateMask) = argLocal(updateMask);
    end
    stat.mipSec = toc(tMip);

    tg = tic;
    mip2D = gather(mipG);
    argmaxZMap = gather(argmaxZG);
    mean2D = gather(sumG ./ max(single(nz), 1));
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

    candidateData.maxMap = mip2D;
    candidateData.meanMap = mean2D;
    candidateData.argmaxZMap = argmaxZMap;

    tCurve = tic;
    [focusIdx, peakRatioApprox] = iEstimateFocusFromMaps(candidateData);
    bestPatches = iCollectBestPatchesGpu(U0f, phaseBase, k_offset, zVec, ...
        candidateData.candidates, focusIdx, params.invertContrast);
    stat.curveSec = toc(tCurve);

    summary = iFinalizeSummaryFromMaps(candidateData, focusIdx, peakRatioApprox, ...
        bestPatches, zVec, roiParams, [ny, nx, nz]);
end

function [summary, stat] = iFastStatsCpuCore(holo, zVec, lambda, pixel, params)
    [ny, nx] = size(holo);
    nz = numel(zVec);
    stat = iInitStat();

    tGrid = tic;
    [FX, FY] = iFrequencyGrid(nx, ny, pixel, 'double');
    phaseBase = (-1i * pi * lambda) .* (FX.^2 + FY.^2);
    k_offset = (2 * pi / lambda);
    stat.gridSec = toc(tGrid);

    tFft = tic;
    U0f = fft2(holo);
    stat.fftSec = toc(tFft);

    batchSize = iChooseBatchSizeCPU(nx, ny, nz);
    stat.batchSize = batchSize;
    mip2D = zeros(ny, nx, 'single');
    argmaxZMap = ones(ny, nx, 'uint16');
    sum2D = zeros(ny, nx, 'single');

    tMip = tic;
    for s = 1:batchSize:nz
        idx = s:min(s + batchSize - 1, nz);
        zBatch = reshape(zVec(idx), 1, 1, []);
        H = exp(1i * k_offset * zBatch + phaseBase .* zBatch);
        Amp = abs(ifft2(U0f .* single(H)));
        sliceMax = max(max(Amp, [], 1), [], 2);
        if params.invertContrast
            I = bsxfun(@minus, sliceMax, Amp);
        else
            I = Amp;
        end

        sum2D = sum2D + sum(I, 3);
        [batchMax, localArgmax] = max(I, [], 3);
        updateMask = batchMax > mip2D;
        mip2D(updateMask) = batchMax(updateMask);
        batchArgmax = uint16(idx(1) - 1) + uint16(localArgmax);
        argmaxZMap(updateMask) = batchArgmax(updateMask);
    end
    stat.mipSec = toc(tMip);

    roiParams = iBuildAdaptiveRoiParams(params, pixel, zVec);
    tCand = tic;
    candidateData = iDetectCandidates(mip2D, roiParams);
    stat.candidateSec = toc(tCand);
    stat.candidateCount = numel(candidateData.candidates);

    if isempty(candidateData.candidates)
        summary = iEmptySummary(candidateData, zVec, [ny, nx, nz]);
        return;
    end

    candidateData.maxMap = mip2D;
    candidateData.meanMap = sum2D ./ single(max(1, nz));
    candidateData.argmaxZMap = argmaxZMap;

    tCurve = tic;
    [focusIdx, peakRatioApprox] = iEstimateFocusFromMaps(candidateData);
    bestPatches = iCollectBestPatchesCpu(U0f, phaseBase, k_offset, zVec, ...
        candidateData.candidates, focusIdx, params.invertContrast);
    stat.curveSec = toc(tCurve);

    summary = iFinalizeSummaryFromMaps(candidateData, focusIdx, peakRatioApprox, ...
        bestPatches, zVec, roiParams, [ny, nx, nz]);
end

function bestPatches = iCollectBestPatchesGpu(U0f, phaseBase, k_offset, zVec, candidates, bestPatchIdx, invertContrast)
    numCand = numel(candidates);
    bestPatches = cell(numCand, 1);
    uniqueIdx = unique(double(bestPatchIdx(:)).');

    for zIdx = uniqueIdx
        if zIdx < 1
            continue;
        end
        zNow = zVec(zIdx);
        H = exp(1i * k_offset * zNow + phaseBase .* zNow);
        Amp = abs(ifft2(U0f .* single(H)));
        if invertContrast
            sliceMax = max(max(Amp, [], 1), [], 2);
            I = bsxfun(@minus, sliceMax, Amp);
        else
            I = Amp;
        end

        pickIdx = find(double(bestPatchIdx(:)) == zIdx);
        for k = 1:numel(pickIdx)
            n = pickIdx(k);
            c = candidates(n);
            bestPatches{n} = gather(I(c.y1:c.y2, c.x1:c.x2));
        end
    end
end

function bestPatches = iCollectBestPatchesCpu(U0f, phaseBase, k_offset, zVec, candidates, bestPatchIdx, invertContrast)
    numCand = numel(candidates);
    bestPatches = cell(numCand, 1);
    uniqueIdx = unique(double(bestPatchIdx(:)).');

    for zIdx = uniqueIdx
        if zIdx < 1
            continue;
        end
        zNow = zVec(zIdx);
        H = exp(1i * k_offset * zNow + phaseBase .* zNow);
        Amp = abs(ifft2(U0f .* single(H)));
        if invertContrast
            sliceMax = max(max(Amp, [], 1), [], 2);
            I = bsxfun(@minus, sliceMax, Amp);
        else
            I = Amp;
        end

        pickIdx = find(double(bestPatchIdx(:)) == zIdx);
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

        peakRatioApprox(n) = mean(wPatch(:)) / max(mean(meanPatch(:)), eps);
    end
end

function summary = iFinalizeSummaryFromMaps(candidateData, focusIdx, ~, bestPatches, zVec, params, volumeSize)
    zVecUm = double(zVec(:)) * 1e6;
    numCand = numel(candidateData.candidates);
    coords3D = zeros(numCand, 3);
    axialCurves = cell(numCand, 1);
    validMask = false(numCand, 1);
    stats = candidateData.stats;

    isSim = isfield(params, 'simulationMode') && params.simulationMode;

    for n = 1:numCand
        c = candidateData.candidates(n);
        maxIdx = min(numel(zVecUm), max(1, double(focusIdx(n))));

        if isSim
            % ---- Simulation mode: Heywood equivalent diameter ----
            % D_eq = 2 * sqrt(Area / pi) based on MIP binary mask area.
            % This is the unified metric for all particle shapes (circle +
            % irregular) — same standard used in HACPI cloud-particle
            % characterization.  Irregular particles get a physically
            % meaningful "equivalent circle" diameter from their projected
            % area, not from the focus-patch Otsu threshold (which assumes
            % compact circular focus spots).
            areaPx = stats(n).Area;
            heywoodDiamPx = 2 * sqrt(areaPx / pi);
            stats(n).EquivDiameter = heywoodDiamPx;
            stats(n).Area = areaPx;  % preserve for downstream volume calc
        else
            % ---- Experimental mode: focus-patch Otsu diameter ----
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
        end

        coords3D(n, :) = [c.centerXY(1), c.centerXY(2), zVecUm(maxIdx)];
        axialCurves{n} = zVecUm(maxIdx);
        diamUm = stats(n).EquivDiameter * iGetFieldOr(params, 'pix_um', 1);
        minValidDiamUm = iGetFieldOr(params, 'minValidDiamUm', 0);
        isLargeEnough = diamUm >= minValidDiamUm;
        validMask(n) = isLargeEnough;
    end

    summary = struct();
    summary.img2D = candidateData.img2D;
    summary.mip2D = candidateData.mip2D;
    summary.bw = candidateData.bw;
    summary.zVecUm = zVecUm;
    summary.volumeSize = volumeSize;
    summary.hasVolume = false;
    summary.validMask = validMask;
    % Diagnostic: pass through dual-threshold masks if present
    if isfield(candidateData, 'bwLow')
        summary.bwLow  = candidateData.bwLow;
        summary.bwHigh = candidateData.bwHigh;
    end
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
    out.img2D = candidateData.img2D;
    out.mip2D = candidateData.mip2D;
    out.bw = candidateData.bw;
    out.zVecUm = double(zVec(:)) * 1e6;
    out.volumeSize = volumeSize;
    out.hasVolume = false;
    out.validMask = false(0, 1);
    if isfield(candidateData, 'bwLow')
        out.bwLow  = candidateData.bwLow;
        out.bwHigh = candidateData.bwHigh;
    end
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
    img2D = AS_normalizeImage(mip2D);
    level = min(max(params.xyThreshold, 0), 1);

    isSim = isfield(params, 'simulationMode') && params.simulationMode;

    % ---- Binarization (common) ----
    bw = imfill(bwareaopen(imbinarize(img2D, level), 3), 'holes');

    if isSim
        % ============================================================
        %  Simulation: dual-threshold + conditional marker watershed.
        %
        %  lowLevel  — captures full particle extent (weak bridges +
        %              extended protrusions of irregular particles).
        %  highLevel — bright nuclei used as markers.
        %
        %  For each low-threshold connected component:
        %    1 marker  → keep whole  (single irregular particle, no split)
        %    N markers → marker-controlled watershed  (N nearby particles
        %                 that merged in the low mask)
        %  Guards: markers must be ≥ minSepPx apart; component must be
        %          larger than a single expected particle area.
        % ============================================================

        highLevel = params.xyThreshold;
        lowLevel  = max(0.12, 0.40 * highLevel);

        bwHigh = imbinarize(img2D, highLevel);
        bwLow  = imbinarize(img2D, lowLevel);

        bwLow = imclose(bwLow, strel('disk', 4));
        bwLow = imfill(bwLow, 'holes');
        bwLow = bwareaopen(bwLow, 3);

        % ---- Conditional marker watershed per component ----
        Llow = bwlabel(bwLow, 8);
        nComp = max(Llow(:));
        bwFinal = false(size(bwLow));

        % Guard thresholds
        %   minSepPx: markers must be ≥ this far apart to be considered
        %             separate particles.  For 50 um particles at 3.45 um/px
        %             this is ~22 px (~75 um).  An irregular single particle
        %             would rarely have internal bright spots this far apart.
        minSepPx = max(15, params.minDiamPx * 1.2);
        %   minAreaForSplit: component must be large enough to plausibly
        %                    contain two minimum-diameter particles.
        minAreaForSplit = 2.2 * pi * (0.5 * params.minDiamPx)^2;

        for id = 1:nComp
            comp = (Llow == id);
            compArea = sum(comp(:));

            % Markers inside this component
            compMarkers = bwHigh & comp;
            compMarkers = bwareaopen(compMarkers, 3);
            Lm = bwlabel(compMarkers, 8);
            nMarks = max(Lm(:));

            if nMarks <= 1
                bwFinal = bwFinal | comp;
                continue;
            end

            % Multiple markers — check if truly distinct particles
            markerProps = regionprops(Lm, 'Centroid');
            centroids = cat(1, markerProps.Centroid);
            dx = centroids(:,1) - centroids(:,1)';
            dy = centroids(:,2) - centroids(:,2)';
            dists = sqrt(dx.^2 + dy.^2);
            dists = dists(triu(true(size(dists)), 1));

            if isempty(dists) || min(dists) < minSepPx || compArea < minAreaForSplit
                bwFinal = bwFinal | comp;
                continue;
            end

            % Attempt marker-controlled watershed within this component
            D = -bwdist(~comp);
            D(~comp) = Inf;
            D2 = imimposemin(D, Lm > 0);
            Lw = watershed(D2);
            splitMask = comp;
            splitMask(Lw == 0) = 0;

            % ---- Fragment check: undo split if any piece is too small ----
            fragProps = regionprops(splitMask, 'EquivDiameter');
            fragDiams = [fragProps.EquivDiameter];
            if any(fragDiams < params.minDiamPx)
                % At least one fragment is below the minimum particle size
                % → this was an irregular single particle, not a cluster
                bwFinal = bwFinal | comp;
            else
                bwFinal = bwFinal | splitMask;
            end
        end

        bw = bwFinal;

        stats_raw = regionprops(bw, img2D, 'Centroid', 'WeightedCentroid', ...
            'EquivDiameter', 'BoundingBox', 'Area', 'Perimeter', ...
            'Solidity', 'Eccentricity', 'Extent', 'MajorAxisLength', 'MinorAxisLength', 'ConvexArea');

        out = struct();
        out.img2D = img2D;
        out.mip2D = mip2D;
        out.bw = bw;
        out.bwLow  = bwLow;   % diagnostic: low-threshold mask
        out.bwHigh = bwHigh;  % diagnostic: high-threshold mask

        if isempty(stats_raw)
            out.stats = struct([]);
            out.roiBoxes = zeros(0, 4);
            out.candidates = struct('centerXY', {}, 'x1', {}, 'x2', {}, 'y1', {}, 'y2', {}, ...
                'roiBox', {}, 'cxLocal', {}, 'cyLocal', {});
            return;
        end

        % Filter by diameter only (no circularity)
        all_diams = [stats_raw.EquivDiameter];
        keep_idx = (all_diams >= params.minDiamPx) & (all_diams <= params.maxDiamPx);
        stats = stats_raw(keep_idx);

        L_conn = labelmatrix(bwconncomp(bw));
        bw = ismember(L_conn, find(keep_idx));

        out.stats = stats;
        nCand = numel(stats);
        out.roiBoxes = zeros(nCand, 4);
        out.candidates = repmat(struct('centerXY', [], 'x1', 1, 'x2', 1, 'y1', 1, 'y2', 1, ...
            'roiBox', [1 1 1 1], 'cxLocal', 1, 'cyLocal', 1), nCand, 1);

        [ny, nx] = size(img2D);
        for n = 1:nCand
            centerXY = stats(n).WeightedCentroid;
            if any(~isfinite(centerXY))
                centerXY = stats(n).Centroid;
            end
            cx = round(centerXY(1));
            cy = round(centerXY(2));

            % ---- BoundingBox + fixed margin (no edge-energy expansion) ----
            bbox = stats(n).BoundingBox;
            marginPx = max(8, round(0.35 * max(bbox(3), bbox(4))));
            if isfield(params, 'searchMinRadiusUm') && isfield(params, 'pix_um')
                marginPx = max(marginPx, round(params.searchMinRadiusUm / params.pix_um));
            end

            x1 = max(1, floor(bbox(1)) - marginPx);
            y1 = max(1, floor(bbox(2)) - marginPx);
            x2 = min(nx, ceil(bbox(1) + bbox(3)) + marginPx);
            y2 = min(ny, ceil(bbox(2) + bbox(4)) + marginPx);

            out.roiBoxes(n, :) = [x1, y1, x2 - x1 + 1, y2 - y1 + 1];
            out.candidates(n).centerXY = centerXY;
            out.candidates(n).x1 = x1;
            out.candidates(n).y1 = y1;
            out.candidates(n).x2 = x2;
            out.candidates(n).y2 = y2;
            out.candidates(n).roiBox = out.roiBoxes(n, :);
            out.candidates(n).cxLocal = min(x2 - x1 + 1, max(1, cx - x1 + 1));
            out.candidates(n).cyLocal = min(y2 - y1 + 1, max(1, cy - y1 + 1));
        end

    else
        % ============================================================
        %  Experimental mode: watershed + adaptive ROI (unchanged)
        % ============================================================

        D = -bwdist(~bw);
        mask = imextendedmin(D, params.watershedHmin);
        D2 = imimposemin(D, mask);
        L_watershed = watershed(D2);
        bw(L_watershed == 0) = 0;

        stats_raw = regionprops(bw, img2D, 'Centroid', 'WeightedCentroid', ...
            'EquivDiameter', 'BoundingBox', 'Area', 'Perimeter', ...
            'Solidity', 'Eccentricity', 'Extent', 'MajorAxisLength', 'MinorAxisLength', 'ConvexArea');

        out = struct();
        out.img2D = img2D;
        out.mip2D = mip2D;
        out.bw = bw;

        if isempty(stats_raw)
            out.stats = struct([]);
            out.roiBoxes = zeros(0, 4);
            out.candidates = struct('centerXY', {}, 'x1', {}, 'x2', {}, 'y1', {}, 'y2', {}, ...
                'roiBox', {}, 'cxLocal', {}, 'cyLocal', {});
            return;
        end

        all_diams = [stats_raw.EquivDiameter];
        all_areas = [stats_raw.Area];
        all_perimeters = [stats_raw.Perimeter];
        circularity = (4 * pi * all_areas) ./ (all_perimeters.^2 + eps);
        keep_idx = (all_diams >= params.minDiamPx) & (all_diams <= params.maxDiamPx) ...
                 & (circularity > params.minCircularity);
        stats = stats_raw(keep_idx);

        L_conn = labelmatrix(bwconncomp(bw));
        bw = ismember(L_conn, find(keep_idx));

        out.stats = stats;
        out.roiBoxes = zeros(numel(stats), 4);
        out.candidates = repmat(struct('centerXY', [], 'x1', 1, 'x2', 1, 'y1', 1, 'y2', 1, ...
            'roiBox', [1 1 1 1], 'cxLocal', 1, 'cyLocal', 1), numel(stats), 1);

        for n = 1:numel(stats)
            centerXY = stats(n).WeightedCentroid;
            if any(~isfinite(centerXY))
                centerXY = stats(n).Centroid;
            end

            cx = round(centerXY(1));
            cy = round(centerXY(2));
            candidate = struct('centerXY', centerXY, 'equivDiameterPx', stats(n).EquivDiameter, ...
                'bboxSizePx', [stats(n).BoundingBox(3), stats(n).BoundingBox(4)]);
            roiInfo = AS_buildAdaptiveROI(img2D, candidate, params);

            out.roiBoxes(n, :) = roiInfo.searchBox;
            out.candidates(n).centerXY = centerXY;
            out.candidates(n).x1 = roiInfo.searchBox(1);
            out.candidates(n).y1 = roiInfo.searchBox(2);
            out.candidates(n).x2 = out.candidates(n).x1 + roiInfo.searchBox(3) - 1;
            out.candidates(n).y2 = out.candidates(n).y1 + roiInfo.searchBox(4) - 1;
            out.candidates(n).roiBox = out.roiBoxes(n, :);
            out.candidates(n).cxLocal = min(out.candidates(n).x2 - out.candidates(n).x1 + 1, max(1, cx - out.candidates(n).x1 + 1));
            out.candidates(n).cyLocal = min(out.candidates(n).y2 - out.candidates(n).y1 + 1, max(1, cy - out.candidates(n).y1 + 1));
        end
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
    if ~isfield(roiParams, 'measureRoiEnergyFraction')
        roiParams.measureRoiEnergyFraction = 0.92;
    end
    if ~isfield(roiParams, 'searchRoiGrowFactor')
        roiParams.searchRoiGrowFactor = 1.25;
    end
    if ~isfield(roiParams, 'edgeEnergyThreshold')
        roiParams.edgeEnergyThreshold = 0.20;
    end
    if ~isfield(roiParams, 'searchMinRadiusUm')
        roiParams.searchMinRadiusUm = 30;
    end
    if ~isfield(roiParams, 'roiMaxRadiusPx')
        roiParams.roiMaxRadiusPx = max(roiParams.roiMinRadius + 4, ceil(250 / max(roiParams.pix_um, eps)));
    end
end

function focusDiamPx = iEstimateFocusDiameterFromPatch(roiFocus, cxLocal, cyLocal, xyThreshold, fallbackDiamPx)
    if isempty(roiFocus) || ~any(roiFocus(:) > 0)
        focusDiamPx = fallbackDiamPx;
        return;
    end
    roiFocus = AS_normalizeImage(roiFocus);
    baseLevel = xyThreshold;
    localLevel = max(0.35, min(0.85, max(baseLevel, 0.6 * xyThreshold)));
    bwFocus = imfill(bwareaopen(imbinarize(roiFocus, localLevel), 1), 'holes');
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
        [~, pickIdx] = min((centers(:, 1) - cxLocal).^2 + (centers(:, 2) - cyLocal).^2);
    end

    focusDiamPx = localStats(pickIdx).EquivDiameter;
    if ~(isfinite(focusDiamPx) && focusDiamPx > 0)
        focusDiamPx = fallbackDiamPx;
    end
end

function batchSize = iChooseBatchSizeGPU(nx, ny, nz, availMem)
    perPageBytes = double(nx) * double(ny) * 36;
    usableMem = max(0, 0.65 * double(availMem) - 0.4e9);
    if usableMem <= 0
        batchSize = 1;
        return;
    end
    batchSize = floor(usableMem / perPageBytes);
    batchSize = max(1, min([nz, batchSize, 64]));
end

function batchSize = iChooseBatchSizeCPU(nx, ny, nz)
    perPageBytes = double(nx) * double(ny) * 36;
    batchSize = max(1, min(nz, floor(1.8e9 / perPageBytes)));
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
    end
end

function out = iInitStat()
    out = struct('uploadSec', 0, 'gridSec', 0, 'fftSec', 0, 'mipSec', 0, ...
        'candidateSec', 0, 'curveSec', 0, 'gatherSec', 0, 'batchSize', 1, ...
        'candidateCount', 0);
end

function params = iFillDefaultParams(params)
    if ~isfield(params, 'minPeakRatio') || ~isfinite(params.minPeakRatio)
        params.minPeakRatio = 1.15;
    end
    if ~isfield(params, 'edgeMargin') || ~isfinite(params.edgeMargin)
        params.edgeMargin = 2;
    end
    if ~isfield(params, 'xyThreshold') || ~isfinite(params.xyThreshold)
        params.xyThreshold = 0.70;
    end
    if ~isfield(params, 'roiScale') || ~isfinite(params.roiScale)
        params.roiScale = 0.50;
    end
    if ~isfield(params, 'roiMinRadius') || ~isfinite(params.roiMinRadius)
        params.roiMinRadius = 15;
    end
    if ~isfield(params, 'minValidDiamUm') || ~isfinite(params.minValidDiamUm)
        params.minValidDiamUm = 0;
    end

    % --- 仿真模式标志 ---
    if ~isfield(params, 'simulationMode') || ~isscalar(params.simulationMode) || ~islogical(params.simulationMode)
        if isfield(params, 'simulationMode') && ~islogical(params.simulationMode)
            params.simulationMode = logical(params.simulationMode);
        else
            params.simulationMode = false;
        end
    end

    % --- 候选粒子筛选参数：实验/仿真两套默认值 ---
    if params.simulationMode
        % 仿真模式：放宽筛选，由物理参数和 MIP 对比度决定
        if ~isfield(params, 'minDiamPx') || ~isfinite(params.minDiamPx)
            params.minDiamPx = 3;
        end
        if ~isfield(params, 'maxDiamPx') || ~isfinite(params.maxDiamPx)
            params.maxDiamPx = 200;
        end
        if ~isfield(params, 'minCircularity') || ~isfinite(params.minCircularity)
            params.minCircularity = 0.3;
        end
        % 仿真模式：默认不做对比度反转（仿真粒子在重建振幅中可能是亮区）
        if ~isfield(params, 'invertContrast') || ~isfinite(params.invertContrast)
            params.invertContrast = false;
        end
        % 仿真模式：降低分水岭 H-minima 深度（仿真 MIP 噪声更少、对比度不同）
        if ~isfield(params, 'watershedHmin')
            params.watershedHmin = 0.3;
        end
    else
        % 实验模式：保持与原有硬编码一致，确保向后兼容
        if ~isfield(params, 'minDiamPx') || ~isfinite(params.minDiamPx)
            params.minDiamPx = 20;
        end
        if ~isfield(params, 'maxDiamPx') || ~isfinite(params.maxDiamPx)
            params.maxDiamPx = 35;
        end
        if ~isfield(params, 'minCircularity') || ~isfinite(params.minCircularity)
            params.minCircularity = 0.6;
        end
        if ~isfield(params, 'invertContrast') || ~isfinite(params.invertContrast)
            params.invertContrast = true;
        end
        if ~isfield(params, 'watershedHmin')
            params.watershedHmin = 0.8;
        end
    end

    % 参数范围校验
    params.minPeakRatio = max(params.minPeakRatio, 1);
    params.edgeMargin = max(0, round(params.edgeMargin));
    params.xyThreshold = min(max(params.xyThreshold, 0), 1);
    params.roiScale = max(0.1, params.roiScale);
    params.roiMinRadius = max(1, round(params.roiMinRadius));
    params.minValidDiamUm = max(0, params.minValidDiamUm);
    params.minDiamPx = max(1, round(params.minDiamPx));
    params.maxDiamPx = max(params.minDiamPx, round(params.maxDiamPx));
    params.minCircularity = min(max(params.minCircularity, 0), 1);
    params.watershedHmin = max(0.01, params.watershedHmin);
end

function value = iGetFieldOr(s, name, fallback)
    if isfield(s, name) && ~isempty(s.(name)) && all(isfinite(s.(name)(:)))
        value = s.(name);
    else
        value = fallback;
    end
end
