function hParticleFig = AS_showParticleSpheres(coords3D, diams_um, pix_um, parentAx)
%AS_SHOWPARTICLESPHERES Display localized particles as enlarged spheres.
%   Converts x/y from pixels to micrometers and renders each particle as a
%   visible sphere centered on the recovered 3D coordinate.

    if isempty(coords3D)
        hParticleFig = [];
        return;
    end
    if nargin < 4
        parentAx = [];
    end

    coordsUm = [coords3D(:, 1) * pix_um, coords3D(:, 2) * pix_um, coords3D(:, 3)];
    if isempty(diams_um)
        baseRadius = max(40, 6 * pix_um);
    else
        baseRadius = max([40, 6 * pix_um, median(diams_um) * 1.5]);
    end

    if isempty(parentAx) || ~isgraphics(parentAx, 'axes')
        hParticleFig = figure( ...
            'Name', 'Sparse Particle Field', ...
            'NumberTitle', 'off', ...
            'Color', [0.05 0.07 0.12], ...
            'Position', [1150 120 820 640]);
        ax = axes('Parent', hParticleFig);
    else
        ax = parentAx;
        hParticleFig = ancestor(ax, 'figure');
        cla(ax);
    end
    hold(ax, 'on');

    [sx, sy, sz] = sphere(20);
    colors = lines(max(size(coordsUm, 1), 7));
    for n = 1:size(coordsUm, 1)
        radiusUm = max(baseRadius, diams_um(min(n, numel(diams_um))) * 1.2 / 2);
        surf(ax, ...
            sx * radiusUm + coordsUm(n, 1), ...
            sy * radiusUm + coordsUm(n, 2), ...
            sz * radiusUm + coordsUm(n, 3), ...
            'FaceColor', colors(mod(n-1, size(colors, 1)) + 1, :), ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.95);
    end

    grid(ax, 'on');
    axis(ax, 'equal');
    view(ax, 32, 24);
    xlabel(ax, ['x (' char(181) 'm)']);
    ylabel(ax, ['y (' char(181) 'm)']);
    zlabel(ax, ['z (' char(181) 'm)']);
    title(ax, 'Localized Particle Field');
    set(ax, 'Color', [0.10 0.11 0.15], ...
        'XColor', [0.85 0.90 1.00], ...
        'YColor', [0.85 0.90 1.00], ...
        'ZColor', [0.85 0.90 1.00], ...
        'FontName', 'Times New Roman');
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
    drawnow;
end
