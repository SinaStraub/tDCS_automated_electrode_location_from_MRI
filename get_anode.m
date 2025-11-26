function coords_ordered= get_anode(coords,numElec,elec_list,subdir)
if numElec==2 && isequal((elec_list), ({'F3','F4'}))
    if coords(1,1)>coords(2,1)
        coords_ordered(1,:)=coords(2,:);%%anode ist first
        coords_ordered(2,:)=coords(1,:);
    else
        coords_ordered(1,:)=coords(1,:);%%anode ist first
        coords_ordered(2,:)=coords(2,:);
    end
elseif numElec==5 && isequal((elec_list), ({'F3','AF3','F1','FC3','F5'})) %the order has to equal the order in elec_list, otherwise:
    %isequal(sort(elec_config), sort({'F3','AF3','F1','FC3','F5'}))
    if size(coords,1)==5
        [coor1ds_sorted,ind1]=sort(coords(:,1));
        [coor2ds_sorted,ind2]=sort(coords(:,2));
        if ind1(3)==ind2(3)%%%electrode which is in the middle in xy plane
            ind=ind1(3);
        else
            coords
            x_val = input('Enter x-coordinate rounded to one decimal place (float): ');%calls for user input
            ind_h=find(round(coords(:,1),1)==x_val);
            ind=ind_h(1);
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
        fprintf('Electrodes are missing for subject "%s"', subdir);
    end
%elseif %put code for other setups here
end

end