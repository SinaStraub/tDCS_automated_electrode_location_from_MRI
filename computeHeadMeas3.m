function computeHeadMeas3(mesh,path_to_fiducials, subject)
% % Author: Sina Straub, sina.straub@gmail.com, sina.straub@unibe.ch
% % Copyright (c) 2025 Sina Straub. Licensed under the GPL v3.

%%% Load mesh 
%%%mesh = mesh_load_gmsh4('sub-h3462.msh');  %%% adjust filename
%%%path_to_fiducials='eeg_positions/Fiducials.csv';
%%%subject='sub-1234';

all_tri = mesh.triangles;            % Nx3
tri_regions = mesh.triangle_regions; % Nx1
vertices = mesh.nodes;               % Mx3

%%% Extract scalp region (e.g., 1005)
scalp_id = 1005;
scalp_mask = tri_regions == scalp_id;
scalp_tri = all_tri(scalp_mask, :);
scalp_vertices = vertices;

%%%Load fiducials 
T = readtable(path_to_fiducials, 'VariableNamingRule','preserve');
get_pos = @(label) [T{strcmpi(T.Var5, label), 2:4}];

nasion = get_pos('Nz');
inion  = get_pos('Iz');
lpa    = get_pos('LPA');
rpa    = get_pos('RPA');

%%% Nearest scalp vertices 
[~, lpa_idx]    = min(vecnorm(scalp_vertices - lpa, 2, 2));
[~, rpa_idx]    = min(vecnorm(scalp_vertices - rpa, 2, 2));

%%% Midline filtering for Nz-Iz 
mid_tol = 5;
midline_idx = find(abs(scalp_vertices(:,1)) < mid_tol);
keep_mid_tri = all(ismember(scalp_tri, midline_idx), 2);
mid_tri = scalp_tri(keep_mid_tri, :);

[~, ~, ~] = unique([midline_idx; mid_tri(:)]);
vmap_mid = containers.Map(midline_idx, 1:numel(midline_idx));
mid_vertices = scalp_vertices(midline_idx, :);
mid_tri_mapped = zeros(size(mid_tri));
for i = 1:3
    mid_tri_mapped(:,i) = cell2mat(values(vmap_mid, num2cell(mid_tri(:,i))));
end

[~, nasion_idx_m] = min(vecnorm(mid_vertices - nasion, 2, 2));
[~, inion_idx_m]  = min(vecnorm(mid_vertices - inion,  2, 2));

Gm = triangulation2graph(mid_tri_mapped, mid_vertices);
path_nodes_m = shortestpath(Gm, nasion_idx_m, inion_idx_m);
path_coords_m = mid_vertices(path_nodes_m, :);
try
path_coords_m = smooth_geodesic_path(path_coords_m, 100);
catch
subject
end
mid_len = sum(vecnorm(diff(path_coords_m),2,2));
fprintf('Midline Nz–Iz arc length: %.2f mm\n', mid_len);

%%% Z-based filtering for crown arc LPA–RPA 
z_margin = 30;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% might require adjustment
z_top = max(scalp_vertices(:,3)) - z_margin;
crown_idx = find(scalp_vertices(:,3) > z_top);
keep_crown_tri = all(ismember(scalp_tri, crown_idx), 2);
crown_tri = scalp_tri(keep_crown_tri, :);

[~, ~, ~] = unique([crown_idx; crown_tri(:)]);
vmap_cr = containers.Map(crown_idx, 1:numel(crown_idx));
crown_vertices = scalp_vertices(crown_idx, :);
crown_tri_mapped = zeros(size(crown_tri));
for i = 1:3
    crown_tri_mapped(:,i) = cell2mat(values(vmap_cr, num2cell(crown_tri(:,i))));
end

[~, lpa_idx_cr] = min(vecnorm(crown_vertices - lpa, 2, 2));
[~, rpa_idx_cr] = min(vecnorm(crown_vertices - rpa, 2, 2));

Gcrown = triangulation2graph(crown_tri_mapped, crown_vertices);
path_nodes_cr = shortestpath(Gcrown, lpa_idx_cr, rpa_idx_cr);
path_coords_cr = crown_vertices(path_nodes_cr, :);
try
path_coords_cr = smooth_geodesic_path(path_coords_cr, 100);
catch
subject
end
geodist_cr = sum(vecnorm(diff(path_coords_cr), 2, 2));

%%% Connect anatomical LPA to crown LPA geodesically on scalp mesh 
Gfull = triangulation2graph(scalp_tri, scalp_vertices);
crown_lpa_global_idx = crown_idx(lpa_idx_cr);
crown_rpa_global_idx = crown_idx(rpa_idx_cr);

%%% LPA to crown-LPA
path_lpa_to_crown = shortestpath(Gfull, lpa_idx, crown_lpa_global_idx);
coords_lpa_to_crown = scalp_vertices(path_lpa_to_crown, :);
try
coords_lpa_to_crown = smooth_geodesic_path(coords_lpa_to_crown, 100);
catch
subject
end
geodist_lpa_seg = sum(vecnorm(diff(coords_lpa_to_crown), 2, 2));

%%% RPA to crown-RPA
path_rpa_to_crown = shortestpath(Gfull, rpa_idx, crown_rpa_global_idx);
coords_rpa_to_crown = scalp_vertices(path_rpa_to_crown, :);
try
coords_rpa_to_crown = smooth_geodesic_path(coords_rpa_to_crown, 100);
catch
subject
end
geodist_rpa_seg = sum(vecnorm(diff(coords_rpa_to_crown), 2, 2));

%%% Total arc from LPA to RPA via crown 
total_arc_dist = geodist_lpa_seg + geodist_cr + geodist_rpa_seg;

%%% Plot 
hold on; axis equal; view(3);
trisurf(scalp_tri, scalp_vertices(:,1), scalp_vertices(:,2), scalp_vertices(:,3), ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', 'cyan');

%%% Fiducials
plot3(lpa(1), lpa(2), lpa(3), 'mo', 'MarkerSize', 10, 'DisplayName', 'LPA');
plot3(rpa(1), rpa(2), rpa(3), 'ko', 'MarkerSize', 10, 'DisplayName', 'RPA');

%%% Crown arc
plot3(path_coords_cr(:,1), path_coords_cr(:,2), path_coords_cr(:,3), 'b-', 'LineWidth', 2, 'DisplayName', 'Crown arc');

%%% Connecting segments
plot3(coords_lpa_to_crown(:,1), coords_lpa_to_crown(:,2), coords_lpa_to_crown(:,3), 'm--', 'LineWidth', 2, 'DisplayName', 'LPA to crown');
plot3(coords_rpa_to_crown(:,1), coords_rpa_to_crown(:,2), coords_rpa_to_crown(:,3), 'k--', 'LineWidth', 2, 'DisplayName', 'RPA to crown');
xlabel('X'); ylabel('Y'); zlabel('Z');

%%% Output distances 
fprintf('Total LPA–RPA arc: %.2f mm\n', total_arc_dist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Total circumference 
%%% Find fractional points along midline and crown arcs 
%%% Nz–Iz points
p_Nz = interp_geodesic_point(path_coords_m, 0.10);  % 10%
p_Iz = interp_geodesic_point(path_coords_m, 0.90);  % 90%

%%% LPA–RPA points
full_coords_lpa_to_rpa = [
    coords_lpa_to_crown;
    path_coords_cr(2:end-1,:);  %%% skip duplicates
    coords_rpa_to_crown(end:-1:1,:)  %%% reversed to match LPA→RPA direction
];
try
full_coords_lpa_to_rpa = smooth_geodesic_path(full_coords_lpa_to_rpa, 100);
catch
subject
end
p_LPA = interp_geodesic_point(full_coords_lpa_to_rpa, 0.10);
p_RPA = interp_geodesic_point(full_coords_lpa_to_rpa, 0.90);

%%% Find nearest scalp vertices to these points
[~, idx_p_Nz] = min(vecnorm(scalp_vertices - p_Nz, 2, 2));
[~, idx_p_Iz] = min(vecnorm(scalp_vertices - p_Iz, 2, 2));
[~, idx_p_LPA] = min(vecnorm(scalp_vertices - p_LPA, 2, 2));
[~, idx_p_RPA] = min(vecnorm(scalp_vertices - p_RPA, 2, 2));

%%% Build full geodesic loop: Iz -> LPA -> Nz -> RPA -> Iz
lambda=5;
Gfull= triangulation2graph_zpenalized(scalp_tri, scalp_vertices, lambda);

seg1 = shortestpath(Gfull, idx_p_Iz, idx_p_LPA);
seg2 = shortestpath(Gfull, idx_p_LPA, idx_p_Nz);
seg3 = shortestpath(Gfull, idx_p_Nz, idx_p_RPA);
seg4 = shortestpath(Gfull, idx_p_RPA, idx_p_Iz);

loop_path = [seg1, seg2(2:end), seg3(2:end), seg4(2:end)];

loop_coords = scalp_vertices(loop_path, :);
try
loop_coords = smooth_geodesic_path(loop_coords, 100);
catch
subject
end
circumference_custom = sum(vecnorm(diff(loop_coords), 2, 2));
fprintf('Custom circumference (10%% offset points): %.2f mm\n', circumference_custom);

%%% Plot the loop
plot3(loop_coords(:,1), loop_coords(:,2), loop_coords(:,3), 'g-', 'LineWidth', 2, 'DisplayName', 'Custom Circumference');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visualization 
trisurf(scalp_tri, scalp_vertices(:,1), scalp_vertices(:,2), scalp_vertices(:,3), ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'cyan');

plot3(nasion(1), nasion(2), nasion(3), 'ro', 'MarkerSize', 10, 'DisplayName', 'Nasion');
plot3(inion(1), inion(2), inion(3), 'go', 'MarkerSize', 10, 'DisplayName', 'Inion');
plot3(lpa(1), lpa(2), lpa(3), 'mo', 'MarkerSize', 10, 'DisplayName', 'LPA');
plot3(rpa(1), rpa(2), rpa(3), 'ko', 'MarkerSize', 10, 'DisplayName', 'RPA');

 plot3(path_coords_m(:,1), path_coords_m(:,2), path_coords_m(:,3), 'r-', 'LineWidth', 2, 'DisplayName', 'Midline Nz–Iz');
 plot3(path_coords_cr(:,1), path_coords_cr(:,2), path_coords_cr(:,3), 'b-', 'LineWidth', 2, 'DisplayName', 'Crown arc LPA–RPA');

xlabel('X'); ylabel('Y'); zlabel('Z');
data=table(string(subject),round(mid_len,2),round(total_arc_dist,2),round(circumference_custom,2), 'VariableNames', { 'Subject','Nz-Iz in mm','LPA-RPA','Circumference in mm'});
 write_or_append_tsv('headMeas_allSub.tsv', data)

%%% Helper function 
function G = triangulation2graph(tri, vertices)
    edges = [tri(:, [1 2]); tri(:, [2 3]); tri(:, [3 1])];
    edges = unique(sort(edges, 2), 'rows');
    dists = vecnorm(vertices(edges(:,1), :) - vertices(edges(:,2), :), 2, 2);
    A = sparse(edges(:,1), edges(:,2), dists, size(vertices,1), size(vertices,1));
    A = max(A, A');
    G = graph(A);
end
function G = triangulation2graph_zpenalized(tri, vertices, lambda)
%%% TRIANGULATION2GRAPH_ZPENALIZED Create a graph with optional z-penalty
%
%%% Inputs:
%%%   tri      - Nx3 triangulation matrix
%%%   vertices - Mx3 vertex coordinates
%%%   lambda   - scalar penalty multiplier for z-differences (e.g., 5)
%
%%% Output:
%%%   G        - MATLAB graph object

    %%% Build undirected edges from triangulation
    edges = [tri(:, [1 2]); tri(:, [2 3]); tri(:, [3 1])];
    edges = unique(sort(edges, 2), 'rows');

    %%% Coordinates of edge endpoints
    v1 = vertices(edges(:,1), :);
    v2 = vertices(edges(:,2), :);

    %%% Compute edge weights: Euclidean distance + lambda * |dz|
    eucl_dist = vecnorm(v1 - v2, 2, 2);
    dz_penalty = abs(v1(:,3) - v2(:,3));
    weights = eucl_dist + lambda * dz_penalty;

    %%% Build symmetric sparse adjacency matrix
    A = sparse(edges(:,1), edges(:,2), weights, size(vertices,1), size(vertices,1));
    A = max(A, A');  %%% ensure undirected

    %%% Create graph object
    G = graph(A);
end

function point = interp_geodesic_point(path_coords, rel_pos)
    %%% Compute cumulative geodesic arc length
    dists = [0; cumsum(vecnorm(diff(path_coords), 2, 2))];
    total_len = dists(end);
    target_len = rel_pos * total_len;

    %%% Interpolate x, y, z separately
    point = interp1(dists, path_coords, target_len);  %%% Should return 1x3
end

function smooth_coords = smooth_geodesic_path(path_coords, n_points)
    %%% Compute cumulative arc length
    dists = [0; cumsum(vecnorm(diff(path_coords), 2, 2))];
    
    %%% Parametric spline for each dimension
    pp_x = spline(dists, path_coords(:,1)');
    pp_y = spline(dists, path_coords(:,2)');
    pp_z = spline(dists, path_coords(:,3)');
    
    %%% Sample points uniformly along arc length
    sample_locs = linspace(0, dists(end), n_points);
    
    %%% Evaluate splines
    smooth_x = ppval(pp_x, sample_locs);
    smooth_y = ppval(pp_y, sample_locs);
    smooth_z = ppval(pp_z, sample_locs);
    
    smooth_coords = [smooth_x', smooth_y', smooth_z'];
end

end