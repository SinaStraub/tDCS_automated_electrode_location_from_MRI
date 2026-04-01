function [elec_pos,elec_list]=get_standard_eeg_pos(path_to_eeg_pos,elec_list)
% % Author: Sina Straub, sina.straub@gmail.com, sina.straub@unibe.ch
% % Copyright (c) 2025 Sina Straub. Licensed under the GPL v3.

%%% reads coordinates from EEG10-10 and plots them
T = readtable([path_to_eeg_pos,'/eeg_positions/EEG10-10_UI_Jurak_2007.csv'], 'VariableNamingRule','preserve');
get_pos = @(label) [T{strcmpi(T.Var5, label), 2:4}];
for ll=1:length(elec_list)
elec_pos(:,ll) = get_pos(elec_list{ll});
end

end