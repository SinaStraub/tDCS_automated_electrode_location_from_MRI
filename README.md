The core algorithm is discribed in https://www.biorxiv.org/content/10.64898/2025.12.08.693046v1.full.pdf, please cite, if you use this algorithm. This package provides a few additional functions which have not formally been evaluated and the code available here is an update from what is described in the preprint. One estimates very roughly distances between head landmarks (inion, nasion, left and right tragus) as well as head circumference. Additionally, the assign_electrodes function automatically assigns standard electrode labels (such as F3, AF1, AF3) to electrode coordinates that were previously unlabeled, and it operates independently of the specific electrode montage.
Dependencies: ROAST (https://www.parralab.org/roast/), SimNIBS-4.5 (https://simnibs.github.io/simnibs/build/html/installation/simnibs_installer.html), NIfTI (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)

Example data can be found here: https://doi.org/10.5281/zenodo.17860288

### Evaluation of the electrode segmentation (optional)
`check_mask_geometry.m` measures the geometry of the segmented electrode regions (`mask_elec.nii`) of one participant: in-plane length and width, thickness, aspect ratio, distance to the scalp, tilt, centroid and long in-plane axis. It is not needed for the coordinate extraction and is run after `tDCS_locate_elec_from_MRI_func`. Optionally it shows the regions as surfaces with their axes and appends them to the head mesh (`<ParticipantID>_with_elec.msh`) for visual inspection in SimNIBS.
    T = check_mask_geometry(fullpath, subdir, numElec);
Requires the NIfTI toolbox, the Image Processing Toolbox and the Statistics and Machine Learning Toolbox; the mesh export additionally requires the SimNIBS Matlab functions.
