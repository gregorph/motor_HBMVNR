function groupedViolinWithBoxesFlexible(A, B, C, colors, subNames, arrayNames)
% A,B,C: 3×n numeric arrays (same n)
% colors: 3×3 RGB for G1..G3 (optional)
% subNames: 1×3 strings/cellstr (optional)
% arrayNames: 1×3 strings/cellstr (optional)

    % --- defaults & checks ---
    if nargin < 4 || isempty(colors),     colors = lines(3); end
    if nargin < 5 || isempty(subNames),   subNames = ["G1","G2","G3"]; end
    if nargin < 6 || isempty(arrayNames), arrayNames = ["Array A","Array B","Array C"]; end

    assert(size(A,1)==3 && size(B,1)==3 && size(C,1)==3, 'A,B,C must each have 3 rows.');
    assert(size(A,2)==size(B,2) && size(B,2)==size(C,2), 'A,B,C must have same number of columns.');
    assert(isnumeric(A) && isnumeric(B) && isnumeric(C), 'A,B,C must be numeric.');
    assert(size(colors,1)==3 && size(colors,2)==3, 'colors must be 3×3 RGB.');

    subNames   = string(subNames(:)).';
    arrayNames = string(arrayNames(:)).';

    % --- y: numeric data vector (force numeric) ---
    y = [A(:); B(:); C(:)];
    y = double(y); % in case anything weird happened upstream

    % --- Build two-line tick labels for each (Array,Subgroup) position ---
    pairLevels = strings(1,9); % A-G1..G3, B-G1..G3, C-G1..G3
    idx = 1;
    for m = 1:3
        for k = 1:3
            pairLevels(idx) = arrayNames(m) + "\newline" + subNames(k);
            idx = idx + 1;
        end
    end

    % Assign each observation to a position and subgroup
    blocks = {A,B,C};
    pairLab = strings(size(y));  % x labels
    subLab  = strings(size(y));  % color group
    offset = 0;
    for m = 1:3
        for k = 1:3
            v = blocks{m}(k,:); v = v(:);
            n = numel(v);
            pairLab(offset+(1:n)) = arrayNames(m) + "\newline" + subNames(k);
            subLab(offset+(1:n))  = subNames(k);
            offset = offset + n;
        end
    end

    % --- x: categorical positions (keep as categorical) ---
    x      = categorical(pairLab, pairLevels, 'Ordinal', true);
    subCat = categorical(subLab,  subNames,   'Ordinal', true);

    % --- Plot ---
    figure; ax = gca; cla(ax); hold(ax,'on');

    % 1) Violin (NUMERIC first arg, CATEGORICAL second)
    vp = violinplot(y, x);
    for i = 1:numel(vp)
        kIdx = mod(i-1,3)+1;
        vp(i).ViolinColor = colors(kIdx,:);
        if isprop(vp(i),'ViolinAlpha'), vp(i).ViolinAlpha = 0.25; end
    end

    % 2) Box inside each violin
    ax.ColorOrder = colors; % ensure subgroup color order
    boxchart(x, y, 'GroupByColor', subCat, ...
        'BoxFaceColor','flat', 'BoxFaceAlpha',0.12, ...
        'MarkerStyle','none', 'BoxWidth',0.35, ...
        'WhiskerLineStyle','-', 'MedianLineColor','k');

    ylabel('Values');
    title('Violin with Box Overlay (3×n × 3 arrays; two-layer labels)');
    legend(subNames, 'Location','bestoutside', 'Title','Subgroup');
    ax.XTickLabelRotation = 0;
    hold(ax,'off');
end
