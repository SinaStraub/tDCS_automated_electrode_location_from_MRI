%%% set environment and add path - do this in a script that calls this
%%% function
% % setenv('SIMNIBSDIR', '/path/to/SimNIBS-4.5/simnibs')
% % setenv('SIMNIBSPYTHON', '/path/to/SimNIBS-4.5/simnibs_env/bin')
% % addpath('/path/to/SimNIBS-4.5/matlab_tools')
% % addpath('/path/to/roast-master')


function tDCS_locate_elec_from_MRI_func(fullpath,subdir,numElec, elec_list,smallblobsize,cutoffpositionz,prior_off,thres_factor,imbalance_thresh)
% % Author: Sina Straub, sina.straub@gmail.com, sina.straub@unibe.ch
% % Copyright (c) 2025 Sina Straub. Licensed under the GPL v3.

%might need to remove "/.../.../SimNIBS-4.5/simnibs_env/lib/python3.11/site-packages/simnibs/examples"
%from path (because it contains repelem.m which is also a Matlab function)

%%%fullpath - location of subject dirs
%%%subdir - list of subjects: m2m_sub-xyz
%%%numElec - number of electrodes
%%% elec_list - list of planned electrode positions such as {'F3','AF3','F1','FC3','F5'}
%%% smallblobsize - 30 (use 30, only adjust for low quality data, e.g. to 50)
%%% cutoffpositionz - no electrodes are expected below this slice index
%%% prior_off - 0 = use prior, 1 = do not use prior
%%% thres_factor - factor that is multiplied with the mean signal intensity
%%% from all voxels that have not been segmented to be tissue. 2 is default and should
%%% work when electrodes/gel are ok visible, otherwise try 0.5, 0.8, 0.9
if nargin<8
    thres_factor=2;
end

if nargin<9
    imbalance_thresh=0.6;
end

if nargin<7
    prior_off=0;
end

figure
h_mesh = mesh_load_gmsh4([fullpath,subdir,'/',subdir(5:end),'.msh']);
try    %%%Get head measures - this is optional
    computeHeadMeas3(h_mesh, [fullpath,subdir,'/','eeg_positions/Fiducials.csv'], subdir(5:end))
catch
    all_tri = h_mesh.triangles;            % Nx3
    tri_regions = h_mesh.triangle_regions; % Nx1
    vertices = h_mesh.nodes;               % Mx3
    %%% Extract scalp region (e.g., 1005)
    scalp_id = 1005;
    scalp_mask = tri_regions == scalp_id;
    scalp_tri = all_tri(scalp_mask, :);
    scalp_vertices = vertices;
    trisurf(scalp_tri, scalp_vertices(:,1), scalp_vertices(:,2), scalp_vertices(:,3), ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'cyan');
end
%%%load tissues, t1 and generate electrode mask and get electrode coord
[epos,elist]=get_standard_eeg_pos([fullpath,subdir],elec_list);
for lll=1:length(elec_list)
    plot3(epos(1,lll), epos(2,lll), epos(3,lll), 'mo', 'MarkerSize', 8, 'DisplayName', elist{lll});
end
%%%extract electrode locations
%%%%tissue seg
h_tissues=load_untouch_nii([fullpath,subdir,'/final_tissues.nii.gz']);
%%%%T1w
h_t1=load_untouch_nii([fullpath,subdir,'/T1.nii.gz']);
%%%initialze mask
mask_elec=ones(size(h_t1.img));
mask_elec(:,:,1:cutoffpositionz)=zeros(size(mask_elec,1),size(mask_elec,2),cutoffpositionz);
mask_elec(h_tissues.img>0)=0;

%%%make sure artifacts are not mistaken for electrodes - this needs to
%%%be adjusted for different electrode locations
vox=h_tissues.hdr.dime.pixdim(2:4);
if prior_off==0
    mask_elec=electrode_loc_prior(mask_elec, h_tissues, vox,elec_list,[fullpath,subdir]);
end
t1w=double(h_t1.img);
t1w=t1w./max(t1w(:));
thres=mean(t1w(mask_elec~=0))*thres_factor;
%%%find image signal from electrode gel:
mask_elec(t1w<thres)=0;    %%% Clean up:
mask_elec=imerode(mask_elec,strel('sphere',1));
mask_elec = bwareaopen(mask_elec, floor(smallblobsize/mean(vox)));  %%% remove small blobs
%%%save electrode segmentation:
h_help=h_tissues;
h_help.img=mask_elec;
save_untouch_nii(h_help,[fullpath,subdir,'/mask_elec.nii']);
%%%connected comps
centroids_voxel = correct_k_largest_clusters(mask_elec, numElec,imbalance_thresh);

%%%Convert coordinates
try
    offset=[h_help.hdr.hist.srow_x(4),h_help.hdr.hist.srow_y(4),h_help.hdr.hist.srow_z(4)];
catch
    offset=[h_help.hdr.hist.qoffset_x,h_help.hdr.hist.qoffset_y,h_help.hdr.hist.qoffset_z];
end
clear coords
clear coords_ordered
A=[h_help.hdr.hist.srow_x(1:3);h_help.hdr.hist.srow_y(1:3);h_help.hdr.hist.srow_z(1:3)];
for ll=1:size(centroids_voxel,1)
    coords(ll,:)=A*((centroids_voxel(ll,:)-1)')+offset';
end

% % %save actual electrode locations - this needs to be adjusdted depending
% % %on the electrode setup, now it works for two electrodes on the
% % %left/right side of the head (e.g. F3 (anode) and F4) or five eletrodes
% % %(middel one is anode)
try
    coords_ordered= get_anode(coords,numElec,elec_list,subdir);
    coords_ordered=assign_electrodes(fullpath, subdir, numElec, elec_list, 1,coords_ordered);
catch
    coords_ordered=assign_electrodes(fullpath, subdir, numElec, elec_list, 0,coords);

end


landmarks=zeros(size(h_help.img));
for cc=1:size( centroids_voxel,1)
    h=round( centroids_voxel(cc,:),0);
    i1=max(1,h(1)-2):min(size(landmarks,1),h(1)+2);i2=max(1,h(2)-2):min(size(landmarks,2),h(2)+2);i3=max(1,h(3)-2):min(size(landmarks,3),h(3)+2);
    landmarks(i1,i2,i3)=cc;
end
h_help.img=landmarks;
save_untouch_nii(h_help,[fullpath,subdir,'/landmarks.nii'])%%%visualize centroids in coordinate system
% % epos list landmarks
for ll=1:size(epos,2)
    epos_nii(:,ll)=inv(A)*(epos(:,ll)-offset')+1;
end

landmarks_list=zeros(size(h_help.img));
for cc=1:size( epos_nii,2)
    h=round( epos_nii(:,cc),0);
    i1=max(1,h(1)-2):min(size(landmarks_list,1),h(1)+2);i2=max(1,h(2)-2):min(size(landmarks_list,2),h(2)+2);i3=max(1,h(3)-2):min(size(landmarks_list,3),h(3)+2);
    landmarks_list(i1,i2,i3)=cc;
end
h_help.img=landmarks_list;
save_untouch_nii(h_help,[fullpath,subdir,'/landmarks_list.nii'])

% % % setup SimNibs simulation and run
% S1=simNibs_template(fullpath,subdir,elec_list,'/elec_pos_list_rec','rec');
% % % run SimNibs
% run_simnibs(S1)
% S2=simNibs_template(fullpath,subdir,coords_ordered,'/elec_pos_actual_rec','rec');
% % % run SimNibs
% run_simnibs(S2)

end
