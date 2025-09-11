function S=simNibs_template(fullpath,subject,elist,outputfolder,electype)
%electype='el', 'rec','ring'

% %this function includes code from SimNIBS-4.5 examples "tDCS_advanced": % Copyright (c) 2021 SimNIBS developers. Licensed under the GPL v3. 
% % code modified by Sina Straub, sina.straub@gmail.com, sina.straub@unibe.ch
% % Copyright (c) 2025 Sina Straub. Licensed under the GPL v3.

%%% Create a new tDCS simulation struct
S = sim_struct('SESSION');
if size(elist,1)<2
    enum=length(elist);
else
    enum=size(elist,1);
end
%%% Define tDCS
S.subpath = [fullpath,subject];
S.pathfem = [fullpath,subject,outputfolder];
S.map_to_fsavg = true;
S.map_to_MNI = true;
S.map_to_surf = true;   %%%  Map to subject's middle gray matter surface
S.map_to_vol = true;    %%%  Save as nifti volume
S.tissues_in_niftis = [1,2,3]; % Results in the niftis will be masked
% to only show WM (1), GM (2), CSF(3)
% (standard: only GM)
% To get fields everywhere:
%    S.tissues_in_niftis = 'all'

%S.open_in_gmsh = true; % show results in gmsh (not for the the niftis)
S.fields = 'eEjJ'; % save e-field and current density
%%% elliptic electrodes with holes (defines below):
if electype=="ring" || electype=="el"
    etype='ellipse';
    edim=[12,12,3];
    %%%Option2:rectangle electrodes:
elseif electype=="rec"
    etype='rectangle';
    edim=[14,17.5,1];
end

% add a TDCS simulation
S.poslist{1} = sim_struct('TDCSLIST');
S.poslist{1}.currents = [1e-3, -1e-3]; %%% Current going through each channel, in Ampere

%%%choose different electrode configurations
if  enum==5
    S.poslist{1}.currents = [1e-3, -0.25e-3,-0.25e-3,-0.25e-3,-0.25e-3]; % Current going through each channel, in Ampere
    for ll=1:enum
        % define electrodes - first is anode
        S.poslist{1}.electrode(ll).channelnr = ll;
        if electype=="ring"
            S.poslist{1}.electrode(ll).holes = sim_struct('ELECTRODE');
        end% Define the hole
        if ~isnumeric(elist)
            S.poslist{1}.electrode(ll).centre = elist{ll};
            if electype=="ring"
                S.poslist{1}.electrode(ll).holes.centre = elist{ll}; % Hole is also centered
            end
        else
            S.poslist{1}.electrode(ll).centre = elist(ll,:);
            if electype=="ring"
                S.poslist{1}.electrode(ll).holes.centre = elist(ll,:); % Hole is also centered
            end
        end
        S.poslist{1}.electrode(ll).shape = etype;
        S.poslist{1}.electrode(ll).dimensions = edim(1:2);
        S.poslist{1}.electrode(ll).thickness = edim(3);
        % Hole
        if electype=="ring"
            S.poslist{1}.electrode(ll).holes.shape = 'ellipse'; % Shape of the hole
            S.poslist{1}.electrode(ll).holes.dimensions = [6, 6]; % Diameter of 5mm
        end
    end
elseif enum==2
    S.poslist{1}.currents = [1e-3, -1e-3]; % Current going through each channel, in Ampere
    for ll=1:enum
        % define electrodes - first is anode
        S.poslist{1}.electrode(ll).channelnr = ll;
        if electype=="ring"
            S.poslist{1}.electrode(ll).holes = sim_struct('ELECTRODE'); % Define the hole
        end
        if ~isnumeric(elist)
            S.poslist{1}.electrode(ll).centre = elist{ll};
            if electype=="ring"
                S.poslist{1}.electrode(ll).holes.centre = elist{ll}; % Hole is also centered
            end
        else
            S.poslist{1}.electrode(ll).centre = elist(ll,:);
            if electype=="ring"
                S.poslist{1}.electrode(ll).holes.centre = elist(ll,:); % Hole is also centered
            end
        end
        S.poslist{1}.electrode(ll).shape = etype;
        S.poslist{1}.electrode(ll).dimensions = edim(1:2);
        S.poslist{1}.electrode(ll).thickness = edim(3);
        if electype=="ring"
            S.poslist{1}.electrode(ll).holes.shape = 'ellipse'; % Shape of the hole
            S.poslist{1}.electrode(ll).holes.dimensions = [6, 6]; % Diameter of 5mm
        end
    end

end

end