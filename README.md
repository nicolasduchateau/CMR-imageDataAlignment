Code for the alignment of cardiac MRI data to a common reference, available under the license [CeCILL-B](http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html)

# post-processing & spatial alignment on cardiac MRI data (delayed enhancement, segmentations from CVI42 software).

Release 1.0 = MATLAB scripts

Author: Nicolas Duchateau (nicolas.duchateau@creatis.insa-lyon.fr)

Date: June, 2023

Links to the corresponding publications at: <br/> https://www.creatis.insa-lyon.fr/~duchateau/#publications

------------------------------------------------------------------------------------------------------------------------
**NOTICE:**

This code is made open-access. Comments and bug reports are welcome, as well as feedback on its possible improvements.

Published reports of research using this code (or a modified version) may cite the following article that describes the method: 

*Duchateau N, Viallon M, Petrusca L, Clarysse P, Mewton N, Belle L, Croisille P. Frontiers in Cardiovascular Medicine 2023;10.* 
https://doi.org/10.3389/fcvm.2023.1136760

------------------------------------------------------------------------------------------------------------------------
**ARCHIVE CONTENT:**

### CODE:

**a1_readCVI42Contours.m** = extract data from Dicom files and/or post-process the image and contours to be stored as .mat files

**a2_normalizeContours.m** = manually identify the LV-RV junction and the valve openings on the corresponding slices, and compute the radial / circumferential / long-axis coordinates on each case

**a3_alignData.m** = transport the imaging and segmentation data from each subject's anatomy to a common reference geometry

### DATA: 
https://www.creatis.insa-lyon.fr/~duchateau/DOWNLOADS/DEMOS/Duchateau_FrontCV_2023/

**Excel file** = pre-identified subfolders to the data, and IDs of the basal and apical slices

**DICOM** = where to store your original DICOM files [not provided]

**EXPORT** = where to store your export from the CVI42 software (.cvi42wsx files = XML encoding) [not provided]

**PROCESSED** = where the .mat files are stored after the first script

**P_DATA** = where the other files are stored, including snapshots of the outputs

