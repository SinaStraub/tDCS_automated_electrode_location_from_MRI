function coords_ordered=assign_electrodes(fullpath, subdir, numElec, elec_list, freezeFirst,coords)
% % Author: Sina Straub, sina.straub@gmail.com, sina.straub@unibe.ch
% % Copyright (c) 2025 Sina Straub. Licensed under the GPL v3.
% % Assigns detected electrodes to the intended EEG list.
% % freezeFirst = 1 : electrode #1, i.e.,anode, is not rematched.
% % coords - either txt '/actual_electrodes.txt' file with coordinates or matrix with columns x,y,z

    %%% Load intended EEG positions
    [epos, ~] = get_standard_eeg_pos([fullpath, subdir], elec_list);
    all_list = epos';    
    %%% Load detected actual electrode coordinates
if isnumeric(coords)
    all_actual=coords;
else
    coords = readmatrix([fullpath, subdir, coords]);
    all_actual = coords(:, 3:5);
end
    
    %%% Prepare distance matrix (XY-plane distances)
    D = pdist2(all_list(:,1:2), all_actual(:,1:2));

    %%% If needed, freeze the first electrode (#1)
    fixed_map = zeros(numElec,1);   % store assignment (actual indices)

    if freezeFirst == 1
        %%% Force electrode 1, actual 1
        fixed_map(1) = 1;

        %%% Remove row/column 1 from Hungarian assignment
        Dsub = D(2:end, 2:end);

        %%% Hungarian matching for electrodes 2..N
        pairs = matchpairs(Dsub, 1e6);

        %%% Record matches (offsetting indices by +1 because of submatrix)
        for k = 1:size(pairs,1)
            intended_idx = pairs(k,1) + 1;   % add 1: original index
            actual_idx   = pairs(k,2) + 1;
            fixed_map(intended_idx) = actual_idx;
        end

    else
        %%% No freeze: match everything normally
        pairs = matchpairs(D, 1e6);
        for k = 1:size(pairs,1)
            fixed_map(pairs(k,1)) = pairs(k,2);
        end
    end
    if any(fixed_map == 0)
        error('%s: no actual electrode assigned to %s', subdir, strjoin(elec_list(fixed_map==0), ', '));
    end

    %%% Build ordered coordinate array
    coords_ordered = zeros(numElec, 3);
    for i = 1:numElec
        coords_ordered(i,:) = all_actual(fixed_map(i), :);
    end
    %%% Write output file
    fid = fopen([fullpath, subdir, '/actual_electrodes_ordered.txt'], 'w');

    for i = 1:numElec
        fprintf(fid, 'electrode %s  %.2f  %.2f  %.2f  5  1\n', ...
            elec_list{i}, ...
            coords_ordered(i,1), coords_ordered(i,2), coords_ordered(i,3));
    end

    fclose(fid);
    for ll=1:numElec
plot3( coords_ordered(ll,1),  coords_ordered(ll,2),  coords_ordered(ll,3), 'mo', 'MarkerSize', ll*2, 'DisplayName', elec_list{ll}, 'Color','b');
    end
  view(190,60)
    title(subdir)
end


