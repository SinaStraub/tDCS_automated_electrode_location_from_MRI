function mask_elec=electrode_loc_prior(mask_elec, h_help, vox,elec_config,datapath)
% % Author: Sina Straub, sina.straub@gmail.com, sina.straub@unibe.ch
% % Copyright (c) 2025 Sina Straub. Licensed under the GPL v3.

help=double(h_help.img);
if  isequal(sort(elec_config), sort({'F3','AF3','F1','FC3','F5'}))
    z_ind(1)=find(squeeze(help(size(help,1)/2,size(help,2)/2,:)),1,'last');%%z-dir
    z_ind(2)=find(squeeze(help(size(help,1)/2,size(help,2)/2+floor(7/vox(1)),:)),1,'last');
    z_ind(3)=find(squeeze(help(size(help,1)/2,size(help,2)/2-floor(7/vox(1)),:)),1,'last');
    z_ind_final=max(z_ind)-floor(63/vox(1));
    yb_ind(1)=find(squeeze(help(size(help,1)/2,:,z_ind_final)),1);%%y-dir
    yb_ind_final=(yb_ind)+floor(44/vox(3));
    yf_ind(1)=find(squeeze(help(size(help,1)/2,:,z_ind_final)),1,'last');
    yf_ind_final=(yf_ind);
    mask_elec(:,1:yb_ind_final,:)=zeros(size(mask_elec,1),yb_ind_final,size(mask_elec,3));
    mask_elec(:,yf_ind_final-floor(0/vox(3)):size(mask_elec,2),:)=zeros(size(mask_elec,1),1+size(mask_elec,2)-(yf_ind_final-floor(0/vox(3))),size(mask_elec,3));%19
    h_help.img=mask_elec;
    save_untouch_nii(h_help,[datapath,'/h1.nii']);%%%this saves a nifty of the preliminary location prior
    %%%this code is for the electrode configuration {'F3','AF3','F1','FC3','F5'}, it removes voxels from the right side of the head
    xr_ind(1)=find(squeeze(help(:,yf_ind_final,z_ind_final)),1,'last');%for left side remove 'last'
    xr_ind_final=(xr_ind)+floor(32/vox(2));%%+x means move to right
    mask_elec(xr_ind_final-floor(48/vox(2)):size(mask_elec,1),:,:)=zeros(1+size(mask_elec,1)-(xr_ind_final-floor(48/vox(2))),size(mask_elec,2),size(mask_elec,3));
    h_help.img=mask_elec;
    save_untouch_nii(h_help,[datapath,'/h2.nii']);%%%this saves a nifty of the location prior

elseif isequal(sort(elec_config), sort({'F3','F4'}))
    z_ind(1)=find(squeeze(help(size(help,1)/2,size(help,2)/2,:)),1,'last');%%z-dir
    z_ind(2)=find(squeeze(help(size(help,1)/2,size(help,2)/2+floor(7/vox(1)),:)),1,'last');
    z_ind(3)=find(squeeze(help(size(help,1)/2,size(help,2)/2-floor(7/vox(1)),:)),1,'last');
    z_ind_final=max(z_ind)-floor(63/vox(1));
    yb_ind(1)=find(squeeze(help(size(help,1)/2,:,z_ind_final)),1);%%y-dir
    yb_ind_final=(yb_ind)+floor(44/vox(3));
    yf_ind(1)=find(squeeze(help(size(help,1)/2,:,z_ind_final)),1,'last');
    yf_ind_final=(yf_ind);
    mask_elec(:,1:yb_ind_final,:)=zeros(size(mask_elec,1),yb_ind_final,size(mask_elec,3));
    mask_elec(:,yf_ind_final-floor(19/vox(3)):size(mask_elec,2),:)=zeros(size(mask_elec,1),1+size(mask_elec,2)-(yf_ind_final-floor(19/vox(3))),size(mask_elec,3));
    h_help.img=mask_elec;
    save_untouch_nii(h_help,[datapath,'/h1.nii']);%%%this saves a nifty of the location prior

%%% more electrode configurations can be added here: elseif ...    
end
end