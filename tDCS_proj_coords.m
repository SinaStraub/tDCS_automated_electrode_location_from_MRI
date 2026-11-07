function [coords_proj,info]=tDCS_proj_coords(fullpath,subdir,coords,nii_name,ordname)
% % Author: Sina Straub, sina.straub@gmail.com, sina.straub@unibe.ch
% % Copyright (c) 2025 Sina Straub. Licensed under the GPL v3.

%%%Projects coordinates onto the scalp surface. Cleaned up version of the
%%%tDCS_proj_coords script: a function, all electrodes at once, no gradient.
%%%
%%%coords     - N x 3 world coordinates (mm). If omitted or empty they are read
%%%             from actual_electrodes_ordered.txt, which holds the UNPROJECTED
%%%             gel centroids.
%%%nii_name   - optional, e.g. 'landmarks_proj.nii'. Writes 5x5x5 cubes at the
%%%             projected positions, same idea as landmarks.nii, so the result
%%%             can be checked on the T1. Omit to skip writing.
%%%
%%%coords_proj- N x 3 world coordinates on the scalp (mm)
%%%info       - .moved_mm per point, .labels, .missed
%%%
%%%The old script built a surface normal from [Dx,Dy,Dz]=gradient(D). That is
%%%wrong: gradient returns (d/dim2, d/dim1, d/dim3), so the first two components
%%%were swapped relative to p, which is in ind2sub order, and the ray marched
%%%off-axis. No normal is needed - marching from p towards the nearest surface
%%%point does the same job with nothing to get wrong.

if nargin<4
    nii_name='';
end
if nargin<5
    ordname='';
end

%%%tissue file: prefer the standard one, fall back to the RAS copy
tisfile=[fullpath,subdir,'/final_tissues.nii.gz'];
if ~isfile(tisfile)
    tisfile=[fullpath,subdir,'/final_tissues_ras.nii'];
end
tis=load_untouch_nii(tisfile);
mask=tis.img==5;%%%5 = scalp
if ~any(mask(:))
    error('%s: no scalp voxels (label 5) in %s',subdir,tisfile);
end

offset=[tis.hdr.hist.srow_x(4),tis.hdr.hist.srow_y(4),tis.hdr.hist.srow_z(4)];
A=[tis.hdr.hist.srow_x(1:3);tis.hdr.hist.srow_y(1:3);tis.hdr.hist.srow_z(1:3)];

labels={};
if nargin<3 || isempty(coords)
    ordfile=get_ordered_file(fullpath,subdir,ordname);
    if isempty(ordfile)
        error('%s: no actual_electrodes_ordered*.txt and no coords given',subdir);
    end
    fid=fopen(ordfile,'r');
    if fid==-1
        error('%s: cannot open %s',subdir,ordfile);
    end
    C=textscan(fid,'electrode %s %f %f %f %f %f');
    fclose(fid);
    labels=C{1}(:)';
    coords=[C{2} C{3} C{4}];
end
n_pts=size(coords,1);

%%%outer shell of the scalp, once for all points
surface=mask & ~imerode(mask,ones(3,3,3));
[xs,ys,zs]=ind2sub(size(mask),find(surface));
surf_pts=[xs ys zs];

coords_proj=zeros(n_pts,3);
moved_mm=zeros(n_pts,1);
missed=false(n_pts,1);

for ll=1:n_pts

    %%%world -> voxel, 1-based, NOT rounded
    p=(A\(coords(ll,:)'-offset'))'+1;

    %%%nearest surface point
    [~,idx]=min(sum((surf_pts-p).^2,2));
    p0=surf_pts(idx,:);

    d=p0-p;
    if norm(d)>1e-8
        d=d/norm(d);
        step=0.25;%%%voxels
        x=p;
        hit=false;
        for s=1:ceil(4*norm(p0-p)/step)+100
            x=x+step*d;
            xi=round(x);
            if any(xi<1) || xi(1)>size(mask,1) || xi(2)>size(mask,2) || xi(3)>size(mask,3)
                break%%%left the volume
            end
            if mask(xi(1),xi(2),xi(3))
                hit=true;
                break
            end
        end
        if hit
            p_proj=x;
        else
            p_proj=p0;%%%fall back to the nearest surface voxel
            missed(ll)=true;
        end
    else
        p_proj=p;
    end

    %%%voxel -> world
    coords_proj(ll,:)=(A*((p_proj-1)')+offset')';
    moved_mm(ll)=norm(coords_proj(ll,:)-coords(ll,:));
end

if any(missed)
    warning('%s: ray missed the scalp for %d point(s), used the nearest surface voxel',subdir,sum(missed));
end

info.moved_mm=moved_mm;
info.labels=labels;
info.missed=missed;

%%%optional nifti for checking, same 5x5x5 cubes as landmarks.nii
if ~isempty(nii_name)
    sz=size(tis.img);
    vol=zeros(sz);
    for ll=1:n_pts
        v=(A\(coords_proj(ll,:)'-offset'))'+1;
        h=round(v);
        i1=max(1,h(1)-2):min(sz(1),h(1)+2);
        i2=max(1,h(2)-2):min(sz(2),h(2)+2);
        i3=max(1,h(3)-2):min(sz(3),h(3)+2);
        if isempty(i1) || isempty(i2) || isempty(i3)
            warning('%s: projected point %d is outside the volume',subdir,ll);
            continue
        end
        vol(i1,i2,i3)=ll;
    end
    tis.img=vol;
    save_untouch_nii(tis,[fullpath,subdir,'/',nii_name]);
end

end
