function T = check_mask_geometry(fullpath, subdir, numElec, imbalance_thresh, do_plot, write_mesh)
% % Author: Sina Straub
% % Copyright (c) 2026 Sina Straub. Licensed under the GPL v3.

%%%Geometry of the segmented electrode regions of one participant. This is an
%%%analysis of the electrode segmentation (mask_elec.nii) and is not needed for
%%%the extraction of electrode coordinates.
%%%
%%%Run tDCS_locate_elec_from_MRI_func first: this function reads its output
%%%mask_elec.nii and the CHARM files in the participant folder.
%%%
%%%fullpath         - path to the location of the participant directory
%%%subdir           - m2m_ParticipantID
%%%numElec          - number of electrodes in the montage
%%%imbalance_thresh - as used for the coordinate extraction (default 0.6)
%%%do_plot          - 1 (default): surface view of the segmented regions with the
%%%                   scalp, the in-plane axes (blue: long, green: short), the
%%%                   shortest principal axis (red) and the scalp normal (black)
%%%write_mesh       - 1: append the segmented regions to the head mesh and save
%%%                   it as ParticipantID_with_elec.msh for display in SimNIBS
%%%                   (default 0; needs the SimNIBS Matlab functions). The
%%%                   regions are surfaces for visual inspection only and are
%%%                   not conductors in any simulation.
%%%
%%%T - table with one row per electrode:
%%%   vol_mm3      volume of the segmented region
%%%   len_mm       longer in-plane extent
%%%   wid_mm       shorter in-plane extent
%%%   aspect       len_mm / wid_mm
%%%   thk_mm       extent along the scalp normal
%%%   flat         thk_mm / wid_mm
%%%   tilt_deg     angle between the shortest principal axis and the scalp normal
%%%   standoff_mm  height of the centroid above the scalp
%%%   gap_mm       distance from the segmented region to the scalp
%%%                (5th percentile of the voxel-to-scalp distances)
%%%   cx,cy,cz     centroid (world mm)
%%%   ax,ay,az     direction of the longer in-plane axis
%%%
%%%Extents are twice the standard deviation of the voxel coordinates along the
%%%respective direction; for a uniformly filled region the full extent is about
%%%1.73 times that value. The scalp normal is the normal of a plane fitted to
%%%the scalp surface within 10 mm of the scalp point nearest the electrode. The
%%%thickness is measured along that normal and the in-plane extents within the
%%%tangent plane, so that the two cannot be interchanged.
%%%
%%%Needs: NIfTI toolbox (load_untouch_nii), Image Processing Toolbox,
%%%Statistics and Machine Learning Toolbox, correct_k_largest_clusters with its
%%%second output (label volume).
%%%
%%%Example:
%%%T = check_mask_geometry('/data/study/', 'm2m_sub-01', 5);

if nargin<4 || isempty(imbalance_thresh)
    imbalance_thresh = 0.6;
end
if nargin<5 || isempty(do_plot)
    do_plot = 1;
end
if nargin<6 || isempty(write_mesh)
    write_mesh = 0;
end
if fullpath(end)~='/' && fullpath(end)~=filesep
    fullpath = [fullpath,'/'];
end
K = numElec;

tis = load_untouch_nii([fullpath,subdir,'/final_tissues.nii.gz']);
msk = load_untouch_nii([fullpath,subdir,'/mask_elec.nii']);
A   = [tis.hdr.hist.srow_x(1:3);tis.hdr.hist.srow_y(1:3);tis.hdr.hist.srow_z(1:3)];
off = [tis.hdr.hist.srow_x(4);tis.hdr.hist.srow_y(4);tis.hdr.hist.srow_z(4)];
vox = abs(det(A));
w   = @(IJK) (A*(IJK'-1) + off)';                 %%% voxel subs -> world mm

%%% same subdivision into individual electrodes as in the coordinate extraction
[~, L]     = correct_k_largest_clusters(msk.img>0, K, imbalance_thresh);
[mx,my,mz] = ind2sub(size(L),find(L>0));
P   = w([mx my mz]);
lab = double(L(sub2ind(size(L),mx,my,mz)));

scalp = tis.img==5;                               %%% 5 = scalp
srf   = scalp & ~imerode(scalp,ones(3,3,3));
[sx,sy,sz] = ind2sub(size(scalp),find(srf));
S = w([sx sy sz]);


%% ===================================================== measurements
T  = table;
G  = cell(K,1);                                   %%% kept for the figure
for e = 1:K
    Pe = P(lab==e,:);
    if size(Pe,1) < 10, continue, end
    c = mean(Pe,1);

    [Vp,Dp] = eig(cov(Pe));  [~,op] = sort(diag(Dp));  Vp = Vp(:,op);

    %%% scalp patch around the NEAREST scalp point, normal pointing outwards
    [~,i0] = min(sum((S-c).^2,2));  p0 = S(i0,:);
    nb = S(sum((S-p0).^2,2) < 10^2, :);
    [Vs,Ds] = eig(cov(nb));  [~,os] = sort(diag(Ds));  nrm = Vs(:,os(1));
    if dot(nrm,(c-mean(nb,1))') < 0, nrm = -nrm; end

    %%% frame from the scalp, so through-plane and in-plane cannot swap
    tmp = zeros(3,1); [~,kmin] = min(abs(nrm)); tmp(kmin) = 1;
    e1  = tmp - (tmp'*nrm)*nrm;  e1 = e1/norm(e1);
    e2  = cross(nrm,e1);
    Q   = Pe - c;
    thk = 2*std(Q*nrm);                            %%% through-plane
    uv  = [Q*e1, Q*e2];                            %%% in the tangent plane
    [V2,D2] = eig(cov(uv));
    [d2,o2] = sort(diag(D2),'descend');  V2 = V2(:,o2);
    len = 2*sqrt(d2(1));  wid = 2*sqrt(d2(2));
    aL  = e1*V2(1,1) + e2*V2(2,1);                 %%% in-plane long
    aS  = e1*V2(1,2) + e2*V2(2,2);                 %%% in-plane short

    tilt     = acosd(abs(dot(Vp(:,1),nrm)));
    standoff = dot(c-p0, nrm');

    near = S(sum((S-c).^2,2) < 40^2, :);
    gap  = prctile(min(pdist2(Pe,near),[],2), 5);

    T = [T; table(e, size(Pe,1)*vox, len, wid, len/wid, thk, thk/wid, ...
                  tilt, standoff, gap, c(1), c(2), c(3), aL(1), aL(2), aL(3), ...
         'VariableNames', ...
         {'elec','vol_mm3','len_mm','wid_mm','aspect','thk_mm','flat', ...
          'tilt_deg','standoff_mm','gap_mm','cx','cy','cz','ax','ay','az'})]; %#ok<AGROW>
    G{e} = struct('c',c,'aL',aL,'aS',aS,'pmin',Vp(:,1),'nrm',nrm, ...
                  'len',len,'wid',wid,'thk',thk);
end
disp(T)

if ~do_plot && ~write_mesh
    return
end

%%% crop around the electrodes; isosurface returns x = dim2, hence [2 1 3]
pad = 30;
lo  = max(1,[min(mx) min(my) min(mz)]-pad);
hi  = min(size(L),[max(mx) max(my) max(mz)]+pad);
cr  = @(B) B(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3));
v2w = @(Vx) (A*((Vx(:,[2 1 3])+lo-2)') + off)';


%% ===================================================== surface view
if do_plot
    figure; hold on; axis equal; grid on; view(3)
    xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)')
    title(subdir,'Interpreter','none')

    fs = reducepatch(isosurface(cr(double(scalp)),0.5),0.2);
    patch('Vertices',v2w(fs.vertices),'Faces',fs.faces, ...
          'FaceColor',[.85 .80 .75],'EdgeColor','none','FaceAlpha',0.2);

    cols = lines(K);
    h    = gobjects(1,4);
    for e = 1:K
        if isempty(G{e}), continue, end
        fe = isosurface(smooth3(cr(double(L==e)),'box',3),0.5);
        patch('Vertices',v2w(fe.vertices),'Faces',fe.faces, ...
              'FaceColor',cols(e,:),'EdgeColor','none','FaceAlpha',0.6);
        g = G{e};  c = g.c;
        q = @(v,ln,col) quiver3(c(1),c(2),c(3),v(1)*ln,v(2)*ln,v(3)*ln,0, ...
                                'Color',col,'LineWidth',2,'MaxHeadSize',0.6);
        h(1) = q(g.aL,   g.len, 'b');
        h(2) = q(g.aS,   g.wid, 'g');
        h(3) = q(g.pmin, g.thk, 'r');
        h(4) = q(g.nrm,  g.thk, 'k');
    end
    legend(h,{'in-plane long','in-plane short','shortest principal axis', ...
              'scalp normal'},'Location','best')
    camlight; lighting gouraud
end


%% ============================== regions appended to the head mesh
if write_mesh
    subid = subdir(5:end);                        %%% m2m_<id> -> <id>
    m = mesh_load_gmsh4([fullpath,subdir,'/',subid,'.msh']);

    Vnew = []; Fnew = []; Rnew = [];
    for e = 1:K
        fe   = isosurface(smooth3(cr(double(L==e)),'box',3),0.5);
        Fnew = [Fnew; fe.faces + size(m.nodes,1) + size(Vnew,1)]; %#ok<AGROW>
        Vnew = [Vnew; v2w(fe.vertices)];                          %#ok<AGROW>
        Rnew = [Rnew; (2000+e)*ones(size(fe.faces,1),1)];         %#ok<AGROW>
    end

    m.nodes            = [m.nodes;            Vnew];
    m.triangles        = [m.triangles;        Fnew];
    m.triangle_regions = [m.triangle_regions; Rnew];
    m.node_data        = {};            %%% lengths no longer match - drop them
    m.element_data     = {};

    mesh_save_gmsh4(m,[fullpath,subdir,'/',subid,'_with_elec.msh']);
end

end
