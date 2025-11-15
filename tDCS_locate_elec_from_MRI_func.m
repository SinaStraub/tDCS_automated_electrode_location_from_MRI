%%% set environment and add path - do this in a scipt that calls this
%%% function
% % setenv('SIMNIBSDIR', '/path/to/SimNIBS-4.5/simnibs')
% % setenv('SIMNIBSPYTHON', '/path/to/SimNIBS-4.5/simnibs_env/bin')
% % addpath('/path/to/SimNIBS-4.5/matlab_tools')
% % addpath('/path/to/roast-master')


function tDCS_locate_elec_from_MRI_func(fullpath,subdir,numElec, elec_list,smallblobsize,cutoffpositionz,thres_factor)
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
%%% thres_factor - factor that is multiplied with the mean signal intensity
%%% from all voxels that have not been segmented to be tissue. 2 is default and should
%%% work when electrodes/gel are ok visible, otherwise try 0.5, 0.8, 0.9
if nargin<7
thres_factor=2;
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
    %%%extract electrode locations
    %%%save mask in RAS
    %%%tissue seg
    h_tissues=load_untouch_nii([fullpath,subdir,'/final_tissues.nii.gz']);
    save_untouch_nii(h_tissues,[fullpath,subdir,'/final_tissues.nii']);
    [mriRAS,isNonRAS] = convertToRAS([fullpath,subdir,'/final_tissues.nii']);
    h_tissues_ras=load_untouch_nii([fullpath,subdir,'/final_tissues_ras.nii']);
    %%%T1w
    h_t1=load_untouch_nii([fullpath,subdir,'/T1.nii.gz']);
    save_untouch_nii(h_t1,[fullpath,subdir,'/T1.nii']);
    [mriRAS,isNonRAS] = convertToRAS([fullpath,subdir,'/T1.nii']);
    h_t1_ras=load_untouch_nii([fullpath,subdir,'/T1_ras.nii']);
    %%%initialze mask
    mask_elec=ones(size(h_t1_ras.img));

    mask_elec(:,:,1:cutoffpositionz)=zeros(size(mask_elec,1),size(mask_elec,2),cutoffpositionz);

    mask_elec(h_tissues_ras.img>0)=0;

    %%%make sure artifacts are not mistaken for electrodes - this needs to
    %%%be adjusted for different electrode locations
    vox=h_tissues_ras.hdr.dime.pixdim(2:4);%%zxy - ras
    mask_elec=electrode_loc_prior(mask_elec, h_tissues_ras, vox,elec_list,[fullpath,subdir]);
    t1w=double(h_t1_ras.img);
    t1w=t1w./max(t1w(:));
    thres=mean(t1w(mask_elec~=0))*thres_factor;
    %%%find image signal from electrode gel:
    mask_elec(t1w<thres)=0;    %%% Clean up:
    mask_elec=imerode(mask_elec,strel('sphere',1));
    mask_elec = bwareaopen(mask_elec, floor(smallblobsize/mean(vox)));  %%% remove small blobs
    %%%save electrode segmentation in ras coordinate system:
    h_help=h_tissues_ras;
    h_help.img=mask_elec;
    save_untouch_nii(h_help,[fullpath,subdir,'/mask_elec_ras.nii']);
    %%%connected comps
    centroids_voxel = correct_k_largest_clusters(mask_elec, numElec,0.6);

    %%%Convert coordinates
    offset=[h_help.hdr.hist.qoffset_x,h_help.hdr.hist.qoffset_y,h_help.hdr.hist.qoffset_z];
    clear coords
    clear coords_ordered
    A=[h_help.hdr.hist.srow_x(1:3);h_help.hdr.hist.srow_y(1:3);h_help.hdr.hist.srow_z(1:3)];
    for ll=1:size(centroids_voxel,1)
        coords(ll,:)=A*((centroids_voxel(ll,:))')+offset';
    end

    hold on;plot3(coords(:,1), coords(:,2), coords(:,3), 'ro','MarkerSize', 6,'Color','b');
    view(190,60)
    title(subdir)

    % % %save actual electrode locations - this needs to be adjusdted depending
    % % %on the electrode setup, now it works for two electrodes on the
    % % %left/right side of the head (e.g. F3 (anode) and F4) or five eletrodes
    % % %(middel on is anode)
    fid = fopen([fullpath,subdir,'/actual_electrodes.txt'], 'w');
    coords_ordered= get_anode(coords,numElec,elec_list,subdir);
    % if numElec==2
    %     if coords(1,1)>coords(2,1)
    %         coords_ordered(1,:)=coords(2,:);%%anode ist first
    %         coords_ordered(2,:)=coords(1,:);
    %     else
    %         coords_ordered(1,:)=coords(1,:);%%anode ist first
    %         coords_ordered(2,:)=coords(2,:);
    %     end
    % elseif numElec==5
    %     if size(coords,1)==5
    %         [coor1ds_sorted,ind1]=sort(coords(:,1));
    %         [coor2ds_sorted,ind2]=sort(coords(:,2));
    %         if ind1(3)==ind2(3)%%%electrode which is in the middle in xy plane
    %             ind=ind1(3);
    %         else
    %             coords
    %             x_val = input('Enter x-coordinate rounded to one decimal place (float): ');%calls for user input
    %             ind_h=find(round(coords(:,1),1)==x_val);
    %             ind=ind_h(1);
    %         end
    %         coords_ordered(1,:)=coords(ind,:);%%%anode ist first
    %         r_ind=2;
    %         for i=1:size(coords,1)
    %             if i~=ind
    %                 coords_ordered(r_ind,:)=coords(i,:);
    %                 r_ind=r_ind+1;
    %             end
    %         end
    %     else
    %         fprintf('Electrodes are missing for subject "%s"', subdir);
    %     end
    % end

    for i = 1:size(coords_ordered,1)
        if i==1
            fprintf(fid, 'electrode A%d  %.2f  %2f  %.2f  5  1\n', i, coords_ordered(i,1), coords_ordered(i,2), coords_ordered(i,3));
        else
            fprintf(fid, 'electrode K%d  %.2f  %.2f  %.2f  5  1\n', i, coords_ordered(i,1), coords_ordered(i,2), coords_ordered(i,3));
        end
    end
    fclose(fid);
    landmarks=zeros(size(h_help.img));
    for cc=1:size( centroids_voxel,1)
        h=round( centroids_voxel(cc,:),0);
        landmarks(h(1)-2:h(1)+2,h(2)-2:h(2)+2,h(3)-2:h(3)+2)=ones(5,5,5)*cc;
    end
    h_help.img=landmarks;
    save_untouch_nii(h_help,[fullpath,subdir,'/landmarks.nii'])%%%visualize centroids in ras coordinate system

    %%setup SimNibs simulation and run
    S1=simNibs_template(fullpath,subdir,elec_list,'/elec_pos_list','ring');
    %%%run SimNibs
    run_simnibs(S1)
    S2=simNibs_template(fullpath,subdir,coords_ordered,'/elec_pos_actual','ring');
    %%%run SimNibs
    run_simnibs(S2)

end
%end