function coords_ordered= get_anode(coords,numElec,elec_list,subdir)
% % Author: Sina Straub, sina.straub@gmail.com, sina.straub@unibe.ch
% % Copyright (c) 2025 Sina Straub. Licensed under the GPL v3.

if numElec==2 && isequal((elec_list), ({'F3','F4'}))
    if coords(1,1)>coords(2,1)
        coords_ordered(1,:)=coords(2,:);%%anode is first
        coords_ordered(2,:)=coords(1,:);
    else
        coords_ordered(1,:)=coords(1,:);%%anode is first
        coords_ordered(2,:)=coords(2,:);
    end
elseif numElec==5 && isequal((elec_list), ({'F3','AF3','F1','FC3','F5'})) %the order has to equal the order in elec_list, otherwise:
    if size(coords,1)==5
        [coor1ds_sorted,ind1]=sort(coords(:,1));
        [coor2ds_sorted,ind2]=sort(coords(:,2));
        if ind1(3)==ind2(3)%%%electrode which is in the middle in xy plane
            ind=ind1(3);
        else
             error('Skipping get_anode.');
            
        end
        coords_ordered(1,:)=coords(ind,:);%%%anode ist first
        r_ind=2;
        for i=1:size(coords,1)
            if i~=ind
                coords_ordered(r_ind,:)=coords(i,:);
                r_ind=r_ind+1;
            end
        end
    else
        error('Electrodes are missing for subject "%s"', subdir);
    end
    elseif numElec==5 && isequal((elec_list), ({'F4','AF4','F2','FC4','F6'})) %the order has to equal the order in elec_list, otherwise:
    if size(coords,1)==5
        [coor1ds_sorted,ind1]=sort(coords(:,1));
        [coor2ds_sorted,ind2]=sort(coords(:,2));
        if ind1(3)==ind2(3)%%%electrode which is in the middle in xy plane
            ind=ind1(3);
        else
            error('Skipping get_anode.');
            
        end
        coords_ordered(1,:)=coords(ind,:);%%%anode ist first
        r_ind=2;
        for i=1:size(coords,1)
            if i~=ind
                coords_ordered(r_ind,:)=coords(i,:);
                r_ind=r_ind+1;
            end
        end
    else
       error('Electrodes are missing for subject "%s" or have a different layout - skipping get anode', subdir);
    end
%elseif %put code for other setups here
end

end