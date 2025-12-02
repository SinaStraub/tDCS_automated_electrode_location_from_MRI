function assign_cathodes(fullpath, subdir, numElec, elec_list, freezeFirst)
% % Author: Sina Straub, sina.straub@gmail.com
% % Assigns detected electrodes to the intended EEG list.
% % freezeFirst = 1 → electrode #1 is fixed and not rematched, i.e.,anode.

    %%% Load intended EEG positions
    [epos, ~] = get_standard_eeg_pos([fullpath, subdir], elec_list);
    all_list = epos';    
    %%% Load detected actual electrode coordinates
    coords = readmatrix([fullpath, subdir, '/actual_electrodes.txt']);
    all_actual = coords(:, 3:5);
    %%% Prepare distance matrix (XY-plane distances)
    D = pdist2(all_list(:,1:2), all_actual(:,1:2));

    %%% If needed, freeze the first electrode (#1)
    fixed_map = zeros(numElec,1);   % store assignment (actual indices)

    if freezeFirst == 1
        %%% Force electrode 1 → actual 1
        fixed_map(1) = 1;

        %%% Remove row/column 1 from Hungarian assignment
        Dsub = D(2:end, 2:end);

        %%% Hungarian matching for electrodes 2..N
        pairs = matchpairs(Dsub, 1e6);

        %%% Record matches (offsetting indices by +1 because of submatrix)
        for k = 1:size(pairs,1)
            intended_idx = pairs(k,1) + 1;   % add 1 → original index
            actual_idx   = pairs(k,2) + 1;
            fixed_map(intended_idx) = actual_idx;
        end

    else
        %%% No freeze → match everything normally
        pairs = matchpairs(D, 1e6);
        for k = 1:size(pairs,1)
            fixed_map(pairs(k,1)) = pairs(k,2);
        end
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
end


