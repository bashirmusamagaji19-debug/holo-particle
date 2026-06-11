function roi = AS_buildAdaptiveROI(img2D, candidate, params)
%AS_BUILDADAPTIVEROI Build an adaptive search ROI around one candidate.

    if nargin < 3 || isempty(params)
        params = struct();
    end
    params = iFillDefaults(params);

    if isempty(img2D) || ~ismatrix(img2D)
        error('AS_buildAdaptiveROI:InvalidImage', 'img2D must be a non-empty 2D image.');
    end
    if ~isfield(candidate, 'centerXY') || numel(candidate.centerXY) ~= 2
        error('AS_buildAdaptiveROI:InvalidCandidate', 'candidate.centerXY is required.');
    end

    [nRows, nCols] = size(img2D);
    cx = min(nCols, max(1, round(candidate.centerXY(1))));
    cy = min(nRows, max(1, round(candidate.centerXY(2))));

    equivDiameterPx = iGetFieldOr(candidate, 'equivDiameterPx', 2 * params.roiMinRadius);
    bboxSizePx = iGetFieldOr(candidate, 'bboxSizePx', [equivDiameterPx, equivDiameterPx]);
    bboxRadiusPx = ceil(max(bboxSizePx(:)) / 2);
    hintedRadiusPx = ceil(equivDiameterPx * max(params.roiScale, 0.75));
    physicalFloorPx = ceil(params.searchMinRadiusUm / max(params.pix_um, eps));

    radius = max([params.roiMinRadius, bboxRadiusPx, hintedRadiusPx, physicalFloorPx]);
    radius = min(radius, params.roiMaxRadiusPx);

    expansionCount = 0;
    while true
        box = iCenteredBox(size(img2D), cx, cy, radius);
        patch = img2D(box(2):box(2)+box(4)-1, box(1):box(1)+box(3)-1);
        edgeEnergy = iEdgeEnergy(patch);

        if edgeEnergy <= params.edgeEnergyThreshold || radius >= params.roiMaxRadiusPx
            break;
        end

        nextRadius = min(params.roiMaxRadiusPx, ceil(radius * params.searchRoiGrowFactor));
        if nextRadius <= radius
            break;
        end
        radius = nextRadius;
        expansionCount = expansionCount + 1;
    end

    roi = struct();
    roi.searchRadiusPx = radius;
    roi.searchBox = box;
    roi.edgeEnergy = edgeEnergy;
    roi.expansionCount = expansionCount;
    roi.cxLocal = cx - box(1) + 1;
    roi.cyLocal = cy - box(2) + 1;
    roi.usedFallback = (radius == params.roiMinRadius);
end

function params = iFillDefaults(params)
    params.roiScale = iScalarOr(params, 'roiScale', 0.50);
    params.roiMinRadius = max(1, round(iScalarOr(params, 'roiMinRadius', 15)));
    params.pix_um = iScalarOr(params, 'pix_um', 1);
    params.searchMinRadiusUm = iScalarOr(params, 'searchMinRadiusUm', 30);
    params.searchRoiGrowFactor = max(1.05, iScalarOr(params, 'searchRoiGrowFactor', 1.25));
    params.edgeEnergyThreshold = min(max(iScalarOr(params, 'edgeEnergyThreshold', 0.20), 0), 1);

    defaultMaxPx = max(params.roiMinRadius + 4, ...
        ceil(iScalarOr(params, 'searchMaxRadiusUm', 250) / max(params.pix_um, eps)));
    params.roiMaxRadiusPx = max(params.roiMinRadius, round(iScalarOr(params, 'roiMaxRadiusPx', defaultMaxPx)));
end

function value = iScalarOr(s, name, fallback)
    if isfield(s, name) && isscalar(s.(name)) && isfinite(s.(name))
        value = double(s.(name));
    else
        value = fallback;
    end
end

function value = iGetFieldOr(s, name, fallback)
    if isfield(s, name) && ~isempty(s.(name))
        value = double(s.(name));
    else
        value = fallback;
    end
end

function box = iCenteredBox(imSize, cx, cy, radius)
    nRows = imSize(1);
    nCols = imSize(2);
    x1 = max(1, cx - radius);
    x2 = min(nCols, cx + radius);
    y1 = max(1, cy - radius);
    y2 = min(nRows, cy + radius);
    box = [x1, y1, x2 - x1 + 1, y2 - y1 + 1];
end

function score = iEdgeEnergy(patch)
    if isempty(patch)
        score = 0;
        return;
    end

    patch = double(patch);
    if numel(patch) == 1 || size(patch, 1) < 2 || size(patch, 2) < 2
        score = 0;
        return;
    end

    border = [patch(1, :), patch(end, :), patch(2:end-1, 1).', patch(2:end-1, end).'];
    bg = median(patch(:));
    peak = max(patch(:));
    denom = max(peak - bg, eps);
    score = max(0, (max(border) - bg) / denom);
    score = min(score, 1);
end
