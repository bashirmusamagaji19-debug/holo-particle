function [coords3D, axialCurves, roiBoxes, bw, img2D, stats, validMask] = ...
% Reference implementation only. The active method-1 baseline for Python
% alignment is Fresnel_reconstructFastStats.m.
    AS_localizeParticles3D(volData, zVecUm, minPeakRatio, edgeMargin, xyThreshold, roiScale, roiMinRadius, roiOptions)
%AS_LOCALIZEPARTICLES3D Localize particles from a reconstructed 3D volume.

    if nargin < 8 || isempty(roiOptions)
        roiOptions = struct();
    end
    if nargin < 7 || isempty(roiMinRadius)
        roiMinRadius = 15;
    end
    if nargin < 6 || isempty(roiScale)
        roiScale = 0.50;
    end
    if nargin < 5 || isempty(xyThreshold)
        xyThreshold = 0.85;
    end
    if nargin < 4 || isempty(edgeMargin)
        edgeMargin = 2;
    end
    if nargin < 3 || isempty(minPeakRatio)
        minPeakRatio = 1.15;
    end

    minPeakRatio = max(minPeakRatio, 1);
    edgeMargin = max(0, round(edgeMargin));
    xyThreshold = min(max(xyThreshold, 0), 1);
    roiScale = max(0.1, roiScale);
    roiMinRadius = max(1, round(roiMinRadius));

    mip2D = max(volData, [], 3);
    img2D = mip2D;
    mipPeak = max(mip2D(:));
    levelAbs = xyThreshold * mipPeak;
    bw = mip2D > levelAbs;
    bw = bwareaopen(bw, 1);
    bw = imfill(bw, 'holes');
    figure,imshow(mip2D,[[]])
    stats = regionprops(bw, mip2D, 'Centroid', 'WeightedCentroid', ...
        'EquivDiameter', 'BoundingBox');

    if isempty(stats)
        coords3D = zeros(0, 3);
        axialCurves = {};
        roiBoxes = zeros(0, 4);
        validMask = false(0, 1);
        return;
    end

    [~, ~, nZ] = size(volData);
    numCand = numel(stats);
    coords3D = zeros(numCand, 3);
    axialCurves = cell(numCand, 1);
    roiBoxes = zeros(numCand, 4);
    validMask = false(numCand, 1);

    roiParams = iBuildAdaptiveRoiParams(roiScale, roiMinRadius, zVecUm, roiOptions);

    for n = 1:numCand
        centerXY = stats(n).WeightedCentroid;
        if any(~isfinite(centerXY))
            centerXY = stats(n).Centroid;
        end

        candidate = struct( ...
            'centerXY', centerXY, ...
            'equivDiameterPx', stats(n).EquivDiameter, ...
            'bboxSizePx', [stats(n).BoundingBox(3), stats(n).BoundingBox(4)]);
        roiInfo = AS_buildAdaptiveROI(mip2D, candidate, roiParams);
        searchBox = roiInfo.searchBox;
        roiBoxes(n, :) = searchBox;

        x1 = searchBox(1);
        y1 = searchBox(2);
        x2 = x1 + searchBox(3) - 1;
        y2 = y1 + searchBox(4) - 1;

        axialCurve = zeros(nZ, 1);
        for k = 1:nZ
            roiSlice = volData(y1:y2, x1:x2, k);
            axialCurve(k) = mean(roiSlice(:));
        end
        axialCurve = smoothdata(axialCurve, 'movmean', 5);
        [peakVal, maxIdx] = max(axialCurve);
        meanVal = mean(axialCurve);
        peakRatio = peakVal / max(meanVal, eps);
        isAwayFromEdge = (maxIdx > edgeMargin) && (maxIdx <= nZ - edgeMargin);

        focusSlice = volData(:, :, maxIdx);
        measureInfo = AS_refineMeasurementROI(focusSlice, centerXY, searchBox, roiParams);
        measureBox = measureInfo.measureBox;
        mx1 = measureBox(1);
        my1 = measureBox(2);
        mx2 = mx1 + measureBox(3) - 1;
        my2 = my1 + measureBox(4) - 1;
        refinedXY = measureInfo.refinedCenterXY;

        stats(n).EquivDiameter = iEstimateFocusDiameter( ...
            focusSlice, mx1, mx2, my1, my2, round(refinedXY(1)), round(refinedXY(2)), ...
            xyThreshold, stats(n).EquivDiameter);

        validMask(n) = isAwayFromEdge && (peakRatio >= minPeakRatio);
        coords3D(n, :) = [refinedXY(1), refinedXY(2), zVecUm(maxIdx)];
        axialCurves{n} = axialCurve;
        % fprintf('粒子 #%d: 初步直径=%.2f px, 精测直径=%.2f px, 所在层=%d \n', ...
        %     n, initialDiam, finalDiam, maxIdx);
        % t = 1
    end
end

function roiParams = iBuildAdaptiveRoiParams(roiScale, roiMinRadius, zVecUm, roiOptions)
    roiParams = roiOptions;
    roiParams.roiScale = roiScale;
    roiParams.roiMinRadius = roiMinRadius;

    if ~isfield(roiParams, 'measureRoiEnergyFraction') || ~isscalar(roiParams.measureRoiEnergyFraction) || ~isfinite(roiParams.measureRoiEnergyFraction)
        roiParams.measureRoiEnergyFraction = 0.92;
    end
    if ~isfield(roiParams, 'searchRoiGrowFactor') || ~isscalar(roiParams.searchRoiGrowFactor) || ~isfinite(roiParams.searchRoiGrowFactor)
        roiParams.searchRoiGrowFactor = 1.25;
    end
    if ~isfield(roiParams, 'edgeEnergyThreshold') || ~isscalar(roiParams.edgeEnergyThreshold) || ~isfinite(roiParams.edgeEnergyThreshold)
        roiParams.edgeEnergyThreshold = 0.20;
    end
    if ~isfield(roiParams, 'pix_um') || ~isscalar(roiParams.pix_um) || ~isfinite(roiParams.pix_um)
        roiParams.pix_um = 1;
    end
    if ~isfield(roiParams, 'searchMinRadiusUm') || ~isscalar(roiParams.searchMinRadiusUm) || ~isfinite(roiParams.searchMinRadiusUm)
        roiParams.searchMinRadiusUm = 30;
    end
    if ~isfield(roiParams, 'roiMaxRadiusPx') || ~isscalar(roiParams.roiMaxRadiusPx) || ~isfinite(roiParams.roiMaxRadiusPx)
        roiParams.roiMaxRadiusPx = max(roiMinRadius + 4, ceil(250 / max(roiParams.pix_um, eps)));
    end
    if ~isfield(roiParams, 'dz') || ~isscalar(roiParams.dz) || ~isfinite(roiParams.dz)
        if numel(zVecUm) > 1
            roiParams.dz = median(abs(diff(zVecUm))) * 1e-6;
        else
            roiParams.dz = 0;
        end
    end
end

function focusDiamPx = iEstimateFocusDiameter(focusSlice, x1, x2, y1, y2, cx, cy, xyThreshold, fallbackDiamPx)
    roiFocus = focusSlice(y1:y2, x1:x2);
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
    % --- 临时增加：看到这个粒子的测量形状 ---
    figure; 
    subplot(1,2,1); imshow(roiFocus, []); title('局部聚焦灰�?);
    subplot(1,2,2); imshow(bwFocus); title(['局部二值化，直径：', num2str(fallbackDiamPx)]);
    % ---------------------------------------
    bwFocus = bwareaopen(bwFocus, 1);
    bwFocus = imfill(bwFocus, 'holes');
    if ~any(bwFocus(:))
        focusDiamPx = fallbackDiamPx;
        return;
    end

    cc = bwconncomp(bwFocus, 8);
    localStats = regionprops(cc, roiFocus, 'EquivDiameter', 'WeightedCentroid', 'Centroid');

    cxLocal = min(size(bwFocus, 2), max(1, round(cx - x1 + 1)));
    cyLocal = min(size(bwFocus, 1), max(1, round(cy - y1 + 1)));
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
