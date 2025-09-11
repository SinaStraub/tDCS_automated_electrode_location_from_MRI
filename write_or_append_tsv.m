function write_or_append_tsv(filename, newData)
% % Author: Sina Straub, sina.straub@gmail.com, sina.straub@unibe.ch
% % Copyright (c) 2025 Sina Straub. Licensed under the GPL v3.

    % write_or_append_tsv: Appends a table to a .tsv file. Writes a new file if it doesn't exist.
    %
    % Usage:
    %   write_or_append_tsv('data.tsv', newTable)
    %
    % Inputs:
    %   filename - String, name of the .tsv file
    %   newData  - MATLAB table, new data to write or append

    if isfile(filename)
        % Read only header of existing file
        fid = fopen(filename, 'r');
        headerLine = fgetl(fid);
        fclose(fid);
        existingVars = strsplit(strtrim(headerLine), '\t');

        newVars = newData.Properties.VariableNames;

        % Check if variable names match (case sensitive)
        if ~isequal(existingVars, newVars)
            warning('Variable names do not match. Overwriting existing file: %s', filename);
            writetable(newData, filename, 'FileType', 'text', 'Delimiter', '\t');
        else
            % Append new data without header
            writetable(newData, filename, 'FileType', 'text', 'Delimiter', '\t', 'WriteVariableNames', false, 'WriteMode', 'append');
        end
    else
        % Write new file with header
        writetable(newData, filename, 'FileType', 'text', 'Delimiter', '\t');
    end
end
