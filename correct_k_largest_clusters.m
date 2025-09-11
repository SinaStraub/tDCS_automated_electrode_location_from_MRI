function centroids_voxel = correct_k_largest_clusters(binaryMask, k, imbalance_thresh)
% % Author: Sina Straub, sina.straub@gmail.com, sina.straub@unibe.ch
% % Copyright (c) 2025 Sina Straub. Licensed under the GPL v3.


%%% Initial CC and stats
CC = bwconncomp(binaryMask, 26);
stats = regionprops3(CC, 'Volume', 'PrincipalAxisLength', 'Centroid');
elongation_metric = (stats.PrincipalAxisLength(:,1) + stats.PrincipalAxisLength(:,2)) ./ stats.PrincipalAxisLength(:,3);
too_elongated = elongation_metric > 10;
size_thres_h=sort(stats.Volume, 'descend');

size_thres=mean(size_thres_h(1:min(k,CC.NumObjects)))*2;
too_large = stats.Volume > size_thres;
suspiciousIdx = find(too_elongated | too_large);
preSplitMasks = {};

for i = 1:numel(suspiciousIdx)
    idx = suspiciousIdx(i);
    pix = CC.PixelIdxList{idx};
    [x, y, z] = ind2sub(CC.ImageSize, pix);
    coords = [x(:), y(:), z(:)];

    %%% skip degenerate cases
    if size(coords,1) < 10
        continue;
    end

    try
        clusterLabels = kmeans(coords, 2, 'Replicates', 3);
    catch
        warning('k-means failed on suspicious component %d', idx);
        continue;
    end

    for c = 1:2
        newMask = false(CC.ImageSize);
        sel = coords(clusterLabels == c, :);
        newMask(sub2ind(CC.ImageSize, sel(:,1), sel(:,2), sel(:,3))) = true;
        CCn = bwconncomp(newMask, 26);

        stats_n = regionprops3(CCn, 'Volume', 'PrincipalAxisLength', 'Centroid');
        elongation_new = (stats_n.PrincipalAxisLength(:,1) + stats_n.PrincipalAxisLength(:,2)) ./ stats_n.PrincipalAxisLength(:,3);
        if nnz(newMask) > 0 & elongation_new<8
            preSplitMasks{end+1} = newMask;
        end
    end
end
validIdx = true(CC.NumObjects, 1);
validIdx(suspiciousIdx) = false;
validPixelIdxList = CC.PixelIdxList(validIdx);
validMasks = cellfun(@(p) buildMaskFromIdx(p, CC.ImageSize), validPixelIdxList, 'UniformOutput', false);

allMasks = [validMasks, preSplitMasks];

binaryMask_cleaned = false(CC.ImageSize);
for i = 1:numel(allMasks)
    binaryMask_cleaned = binaryMask_cleaned | allMasks{i};
end


%%% Get initial connected components of clean mask - may be redundant
CC = bwconncomp(binaryMask_cleaned, 26);
numPixels = cellfun(@numel, CC.PixelIdxList);
[sortedSizes, sortedIdx] = sort(numPixels, 'descend');
%%%  selectedIdx = sortedIdx(1:min(k, numel(sortedIdx)));
selectedIdx = sortedIdx(1:CC.NumObjects);

if numel(selectedIdx) < k-2
    error('Not enough clusters in the mask to extract %d components.', k);
end

%  % Create initial cluster masks
clusterMasks = cell(1, CC.NumObjects);
clusterSizes = zeros(1, CC.NumObjects);
for i = 1:CC.NumObjects
    mask = false(size(binaryMask_cleaned));
    mask(CC.PixelIdxList{selectedIdx(i)}) = true;
    clusterMasks{i} = mask;
    clusterSizes(i) = sortedSizes(i);
end

%%% Remove small and mark large clusters
medianSize = median(clusterSizes);
tooSmall = clusterSizes < (1 - imbalance_thresh) * medianSize;
clusterMasks = clusterMasks(~tooSmall);
clusterSizes = clusterSizes(~tooSmall);

%%% Subdivide large clusters using kmeans until we have k total
while numel(clusterMasks) < k
    [~, idx] = max(clusterSizes);
    maskToSplit = clusterMasks{idx};

    [x, y, z] = ind2sub(size(binaryMask_cleaned), find(maskToSplit));
    coords = [x(:), y(:), z(:)];

    try
        clusterLabels = kmeans(coords, 2, 'MaxIter', 100, 'Replicates', 3);
    catch
        warning('k-means failed on large cluster.');
        break;
    end

    %%% Build two masks from the split
    newMasks = {false(size(binaryMask_cleaned)), false(size(binaryMask_cleaned))};
    for ii = 1:2
        vox = coords(clusterLabels == ii, :);
        newMasks{ii}(sub2ind(size(binaryMask_cleaned), vox(:,1), vox(:,2), vox(:,3))) = true;
    end

    %%% Replace the large mask with two smaller ones
    clusterMasks(idx) = [];
    clusterSizes(idx) = [];
    for m = 1:2
        clusterMasks{end+1} = newMasks{m};
        clusterSizes(end+1) = nnz(newMasks{m});
    end
end

%%% Compute centroids from the final cluster masks
centroids_voxel = zeros(k, 3);
for i = 1:k
    CC1 = bwconncomp(clusterMasks{i}, 26);
    stats = regionprops3(CC1, 'Centroid');
    if height(stats) ~= 1
        warning('Multiple components found in one mask after kmeans. Using largest one.');
        sizes = cellfun(@numel, CC1.PixelIdxList);
        [~, maxIdx] = max(sizes);
        tmpMask = false(size(binaryMask_cleaned));
        tmpMask(CC1.PixelIdxList{maxIdx}) = true;
        stats = regionprops3(tmpMask, 'Centroid');
    end
    %%% Swap X and Y here if necessary
    centroids_voxel(i,1) = stats.Centroid(2); % X (column 2)
    centroids_voxel(i,2) = stats.Centroid(1); % Y (column 1)
    centroids_voxel(i,3) = stats.Centroid(3); % Z
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% less redundant, works well for good data:
% function centroids_voxel = correct_k_largest_clusters(binaryMask, k, imbalance_thresh)
%     % Get initial connected components
%     CC = bwconncomp(binaryMask, 26);
%     numPixels = cellfun(@numel, CC.PixelIdxList);
%     [sortedSizes, sortedIdx] = sort(numPixels, 'descend');
%     selectedIdx = sortedIdx(1:min(k, numel(sortedIdx)));
%
%     if numel(selectedIdx) < k-2
%         error('Not enough clusters in the mask to extract %d components.', k);
%     end
%
%     % Create initial cluster masks
%     clusterMasks = cell(1, k);
%     clusterSizes = zeros(1, k);
%     for i = 1:k
%         mask = false(size(binaryMask));
%         mask(CC.PixelIdxList{selectedIdx(i)}) = true;
%         clusterMasks{i} = mask;
%         clusterSizes(i) = sortedSizes(i);
%     end
%
%     % Remove small and mark large clusters
%     medianSize = median(clusterSizes);
%     tooSmall = clusterSizes < (1 - imbalance_thresh) * medianSize;
%     clusterMasks = clusterMasks(~tooSmall);
%     clusterSizes = clusterSizes(~tooSmall);
%
%     % Subdivide large clusters using kmeans until we have k total
%     while numel(clusterMasks) < k
%         [~, idx] = max(clusterSizes);
%         maskToSplit = clusterMasks{idx};
%
%         [x, y, z] = ind2sub(size(binaryMask), find(maskToSplit));
%         coords = [x(:), y(:), z(:)];
%
%         try
%             clusterLabels = kmeans(coords, 2, 'MaxIter', 100, 'Replicates', 3);
%         catch
%             warning('k-means failed on large cluster.');
%             break;
%         end
%
%         % Build two masks from the split
%         newMasks = {false(size(binaryMask)), false(size(binaryMask))};
%         for ii = 1:2
%             vox = coords(clusterLabels == ii, :);
%             newMasks{ii}(sub2ind(size(binaryMask), vox(:,1), vox(:,2), vox(:,3))) = true;
%         end
%
%         % Replace the large mask with two smaller ones
%         clusterMasks(idx) = [];
%         clusterSizes(idx) = [];
%         for m = 1:2
%             clusterMasks{end+1} = newMasks{m};
%             clusterSizes(end+1) = nnz(newMasks{m});
%         end
%     end
%
%     % Compute centroids from the final cluster masks
%     centroids_voxel = zeros(k, 3);
%     for i = 1:k
%         CC1 = bwconncomp(clusterMasks{i}, 26);
%         stats = regionprops3(CC1, 'Centroid');
%         if height(stats) ~= 1
%             warning('Multiple components found in one mask after kmeans. Using largest one.');
%             sizes = cellfun(@numel, CC1.PixelIdxList);
%             [~, maxIdx] = max(sizes);
%             tmpMask = false(size(binaryMask));
%             tmpMask(CC1.PixelIdxList{maxIdx}) = true;
%             stats = regionprops3(tmpMask, 'Centroid');
%         end
%         % Swap X and Y here if necessary
%         centroids_voxel(i,1) = stats.Centroid(2); % X (column 2)
%         centroids_voxel(i,2) = stats.Centroid(1); % Y (column 1)
%         centroids_voxel(i,3) = stats.Centroid(3); % Z
%     end
% end
% %



function mask = buildMaskFromIdx(pixIdx, imSize)
mask = false(imSize);
mask(pixIdx) = true;
end
