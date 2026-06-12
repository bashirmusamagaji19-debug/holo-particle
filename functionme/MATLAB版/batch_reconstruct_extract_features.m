function batch_reconstruct_extract_features(dataRoot, outputCsv)
% BATCH_RECONSTRUCT_EXTRACT_FEATURES  Batch reconstruct all datasets,
% match particles to ground truth, and extract shape features for
% classifier training.
%
% Usage:
%   batch_reconstruct_extract_features;
%   batch_reconstruct_extract_features('..\仿真图像生成', 'features_table.csv');
%
% Output columns:
%   dataset_id           – folder name (e.g. '0017')
%   particle_id          – ground truth particle id
%   shape_label          – 'circular' | 'irregular'
%   area_px              – MIP binary mask Area (pixels)
%   equiv_diam_px        – Equivalent diameter (pixels)
%   circularity          – 4*pi*Area / Perimeter^2
%   solidity             – Area / ConvexArea
%   eccentricity         – ellipse eccentricity (0=circle, ~1=line)
%   extent               – Area / BoundingBoxArea
%   aspect_ratio         – MajorAxisLength / MinorAxisLength
%   major_axis_px        – Major axis length (pixels)
%   minor_axis_px        – Minor axis length (pixels)
%   perimeter_px         – Perimeter (pixels)
%   match_distance_px    – distance to matched ground truth (pixels)
%   match_valid          – 1 if matched, 0 if unmatched

    % ---------- defaults ----------
    if nargin < 1 || isempty(dataRoot)
        % Assume this script lives in MATLAB版/, data is in ../仿真图像生成
        scriptDir = fileparts(mfilename('fullpath'));
        if isempty(scriptDir)
            scriptDir = pwd;
        end
        dataRoot = fullfile(scriptDir, '..', '仿真图像生成');
    end
    if nargin < 2 || isempty(outputCsv)
        outputDir = fileparts(mfilename('fullpath'));
        if isempty(outputDir), outputDir = pwd; end
        outputCsv = fullfile(outputDir, 'features_table.csv');
    end

    % ---------- discover datasets ----------
    % Scan both the flat 00xx folders AND training_set/ sub-folders.
    datasetDirs = {};

    % (a) Flat numbered folders in dataRoot (e.g. 0013, 0015-0027)
    flatDirs = dir(fullfile(dataRoot, '00*'));
    for f = 1:numel(flatDirs)
        if flatDirs(f).isdir
            dsId = flatDirs(f).name;
            matFile = fullfile(flatDirs(f).folder, dsId, [dsId '_hologram.mat']);
            csvFile = fullfile(flatDirs(f).folder, dsId, [dsId '_ground_truth.csv']);
            if isfile(matFile) && isfile(csvFile)
                % Check if CSV has shape_type (skip non-classifier datasets)
                fid = fopen(csvFile, 'r');
                header = fgetl(fid);
                fclose(fid);
                if contains(header, 'shape_type')
                    datasetDirs{end+1} = fullfile(flatDirs(f).folder, dsId); %#ok<AGROW>
                end
            end
        end
    end

    % (b) training_set/ sub-folders
    trainSetRoot = fullfile(dataRoot, 'training_set');
    if isfolder(trainSetRoot)
        trainDirs = dir(fullfile(trainSetRoot, '0*'));
        for f = 1:numel(trainDirs)
            if trainDirs(f).isdir
                dsId = trainDirs(f).name;
                matFile = fullfile(trainDirs(f).folder, dsId, [dsId '_hologram.mat']);
                csvFile = fullfile(trainDirs(f).folder, dsId, [dsId '_ground_truth.csv']);
                if isfile(matFile) && isfile(csvFile)
                    datasetDirs{end+1} = fullfile(trainDirs(f).folder, dsId); %#ok<AGROW>
                end
            end
        end
    end

    nDatasets = numel(datasetDirs);
    fprintf('Found %d datasets with shape_type labels\n', nDatasets);
    if nDatasets == 0
        error('No valid datasets found in %s or %s', dataRoot, trainSetRoot);
    end

    % ---------- reconstruction parameters (matching all datasets) ----------
    lambda    = 638e-9;               % wavelength (m)
    pixel     = 3.45e-6;              % pixel size (m)
    zMin      = 20e-3;                % z min (m)
    zMax      = 35e-3;                % z max (m)
    zStep     = 20e-6;                % z step (m) — 750 layers for fine MIP
    zVec      = zMin : zStep : zMax;
    % Try GPU, fall back to CPU automatically
    useGPU    = true;

    % Simulation-mode params (passed to Fresnel_reconstructFastStats)
    simParams = struct(...
        'simulationMode', true, ...
        'minDiamPx',       3, ...
        'maxDiamPx',     200, ...
        'minCircularity', 0.3, ...
        'invertContrast', false, ...
        'watershedHmin',  0.3, ...
        'xyThreshold',    0.70, ...
        'minPeakRatio',   1.05, ...
        'edgeMargin',     2, ...
        'roiScale',       0.50, ...
        'roiMinRadius',   15, ...
        'pix_um',         3.45);

    % Matching threshold: maximum allowed distance (pixels) between
    % detected particle and ground truth for a valid match.
    matchThresholdPx = 15;  % tighter: reduce label noise

    % ---------- collect all features ----------
    allRows = {};  % cell array of structs

    for d = 1:nDatasets
        dsDir = datasetDirs{d};
        [~, dsId] = fileparts(dsDir);
        matFile   = fullfile(dsDir, [dsId '_hologram.mat']);
        csvFile   = fullfile(dsDir, [dsId '_ground_truth.csv']);

        if ~isfile(matFile)
            fprintf('[%s] .mat file not found: %s — skipping\n', dsId, matFile);
            continue;
        end
        if ~isfile(csvFile)
            fprintf('[%s] ground truth CSV not found: %s — skipping\n', dsId, csvFile);
            continue;
        end

        fprintf('=== [%s] Reconstructing... ===\n', dsId);

        % ---- load hologram (.mat has variable 'hologram') ----
        loaded = load(matFile);
        if isfield(loaded, 'hologram')
            hologram = loaded.hologram;
        else
            % try first numeric variable
            fns = fieldnames(loaded);
            hologram = loaded.(fns{1});
        end

        % Ensure double [0,1] — .mat stores float64 already
        if isa(hologram, 'uint8')
            hologram = double(hologram) / 255;
        end
        hologram = double(hologram);

        % ---- reconstruction ----
        summary = Fresnel_reconstructFastStats(hologram, zVec, lambda, pixel, useGPU, simParams);

        % ---- ground truth ----
        gt = readtable(csvFile);

        % ---- extract candidate stats ----
        if isfield(summary, 'candidateStats') && ~isempty(summary.candidateStats)
            candidateStats = summary.candidateStats;
        else
            candidateStats = summary.stats;
        end
        if isfield(summary, 'candidateCoords3D') && ~isempty(summary.candidateCoords3D)
            candCoords = summary.candidateCoords3D(:, 1:2);  % (x, y) in pixels
        else
            candCoords = summary.coords3D(:, 1:2);
        end
        nCand = size(candCoords, 1);

        % ---- ground truth positions (MATLAB 1-indexed) ----
        gt_x = gt.x_pixel_matlab;
        gt_y = gt.y_pixel_matlab;
        gt_shape = gt.shape_type;  % cell array: 'sphere' / 'aggregate' / 'polyhedron'

        % ---- match each detected particle to nearest ground truth ----
        for c = 1:nCand
            s = candidateStats(c);

            % build feature row
            row = struct();
            row.dataset_id     = dsId;
            row.particle_id    = NaN;
            row.shape_label    = 'unmatched';
            row.match_valid    = 0;
            row.match_distance_px = NaN;

            % ---- features from regionprops ----
            row.area_px        = getField(s, 'Area', NaN);
            row.perimeter_px   = getField(s, 'Perimeter', NaN);
            row.equiv_diam_px  = getField(s, 'EquivDiameter', NaN);
            row.solidity       = getField(s, 'Solidity', NaN);
            row.eccentricity   = getField(s, 'Eccentricity', NaN);
            row.extent         = getField(s, 'Extent', NaN);
            row.major_axis_px  = getField(s, 'MajorAxisLength', NaN);
            row.minor_axis_px  = getField(s, 'MinorAxisLength', NaN);

            % derived shape features
            if ~isnan(row.area_px) && ~isnan(row.perimeter_px) && row.perimeter_px > 0
                row.circularity = 4 * pi * row.area_px / (row.perimeter_px ^ 2);
            else
                row.circularity = NaN;
            end
            if ~isnan(row.major_axis_px) && ~isnan(row.minor_axis_px) && row.minor_axis_px > 0
                row.aspect_ratio = row.major_axis_px / row.minor_axis_px;
            else
                row.aspect_ratio = NaN;
            end

            % ---- MIP patch intensity features ----
            % Extract a tight patch around the candidate from the raw MIP.
            % Irregular particles tend to have dimmer, less uniform MIP spots.
            if isfield(summary, 'mip2D') && ~isempty(summary.mip2D)
                mip = double(summary.mip2D);
                [nyMip, nxMip] = size(mip);
                % Use candidate ROI box for the patch
                if isfield(summary, 'candidateRoiBoxes') && ~isempty(summary.candidateRoiBoxes)
                    roi = summary.candidateRoiBoxes(c, :);
                elseif isfield(summary, 'roiBoxes') && c <= size(summary.roiBoxes, 1)
                    roi = summary.roiBoxes(c, :);
                else
                    roi = [];
                end
                if ~isempty(roi)
                    rx = max(1, round(roi(1))); ry = max(1, round(roi(2)));
                    rw = min(nxMip - rx + 1, round(roi(3)));
                    rh = min(nyMip - ry + 1, round(roi(4)));
                    if rw > 3 && rh > 3
                        patch = double(mip(ry:ry+rh-1, rx:rx+rw-1));
                        patchV = patch(:);
                        row.mip_mean    = mean(patchV);
                        row.mip_std     = std(patchV);
                        row.mip_max     = max(patchV);
                        row.mip_median  = median(patchV);
                        mn = min(patchV);
                        row.mip_contrast = (row.mip_max - mn) / max(row.mip_max + mn, eps);
                        row.mip_energy  = sum(patchV .^ 2) / numel(patchV);
                        % Texture features
                        row.mip_skewness = iSkewness(patchV);
                        row.mip_entropy  = iImageEntropy(patch);
                        % Gradient magnitude (edge info)
                        [gx, gy] = gradient(patch);
                        gm = sqrt(gx.^2 + gy.^2);
                        row.mip_gradient = mean(gm(:));
                        row.mip_iqr       = iManualIqr(patchV);
                    else
                        row.mip_mean = NaN; row.mip_std = NaN; row.mip_max = NaN;
                        row.mip_median = NaN; row.mip_contrast = NaN; row.mip_energy = NaN;
                        row.mip_skewness = NaN; row.mip_entropy = NaN;
                        row.mip_gradient = NaN; row.mip_iqr = NaN;
                    end
                else
                    row.mip_mean = NaN; row.mip_std = NaN; row.mip_max = NaN;
                    row.mip_median = NaN; row.mip_contrast = NaN; row.mip_energy = NaN;
                    row.mip_skewness = NaN; row.mip_entropy = NaN;
                    row.mip_gradient = NaN; row.mip_iqr = NaN;
                end
            else
                row.mip_mean = NaN; row.mip_std = NaN; row.mip_max = NaN;
                row.mip_median = NaN; row.mip_contrast = NaN; row.mip_energy = NaN;
                row.mip_skewness = NaN; row.mip_entropy = NaN;
                row.mip_gradient = NaN; row.mip_iqr = NaN;
            end

            % ---- match to ground truth ----
            cx = candCoords(c, 1);
            cy = candCoords(c, 2);
            dists = sqrt((gt_x - cx).^2 + (gt_y - cy).^2);
            [minDist, minIdx] = min(dists);

            if minDist <= matchThresholdPx
                row.match_valid       = 1;
                row.match_distance_px = minDist;
                row.particle_id       = gt.particle_id(minIdx);
                gtVal = gt_shape{minIdx};
                if iscell(gtVal), gtVal = gtVal{1}; end
                gtType = char(string(gtVal));
                if strcmpi(gtType, 'sphere')
                    row.shape_label = 'circular';
                else
                    row.shape_label = 'irregular';
                end
            end

            allRows{end+1} = row; %#ok<AGROW>
        end

        fprintf('[%s] done: %d candidates, %d matched\n', dsId, nCand, ...
            sum(cellfun(@(r) r.match_valid, allRows(end-nCand+1:end))));
    end

    % ---------- assemble table and write CSV ----------
    if isempty(allRows)
        fprintf('No features extracted. Check that datasets exist.\n');
        return;
    end

    % Convert cell array of scalar structs to struct array, then to table
    if numel(allRows) == 1
        featureTable = struct2table(allRows{1});
    else
        featureTable = struct2table(cat(1, allRows{:}));
    end

    % Write CSV
    writetable(featureTable, outputCsv);
    fprintf('\nFeature table written to: %s\n', outputCsv);
    fprintf('Total rows: %d\n', height(featureTable));
    fprintf('Matched:   %d\n', sum(featureTable.match_valid));
    fprintf('Unmatched: %d\n', sum(~featureTable.match_valid));

    % Quick summary by class
    matched = featureTable(featureTable.match_valid == 1, :);
    nCirc = sum(strcmp(matched.shape_label, 'circular'));
    nIrr  = sum(strcmp(matched.shape_label, 'irregular'));
    fprintf('Circular: %d, Irregular: %d\n', nCirc, nIrr);
end

function s = iSkewness(x)
    % Manual skewness (avoid Statistics Toolbox dependency)
    x = x(:);
    mu = mean(x);
    m3 = mean((x - mu).^3);
    m2 = mean((x - mu).^2);
    s = m3 / max(m2^1.5, eps);
end

function q = iManualIqr(x)
    % Manual IQR (avoid Statistics Toolbox dependency)
    x = sort(x(:));
    n = numel(x);
    q1 = x(max(1, round(n * 0.25)));
    q3 = x(min(n, round(n * 0.75)));
    q = q3 - q1;
end

function e = iImageEntropy(patch)
    % Histogram entropy of a grayscale image patch
    counts = histcounts(patch(:), 32);
    p = counts / sum(counts);
    p(p == 0) = [];
    e = -sum(p .* log2(p));
end

function v = getField(s, fieldName, default)
    if isfield(s, fieldName) && isfinite(s.(fieldName))
        v = double(s.(fieldName));
    else
        v = default;
    end
end
