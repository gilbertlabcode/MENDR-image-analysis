# MENDR image analysis

This repository contains MATLAB scripts and functions used to perform automatic image analysis of the datasets presented in the article **"An in vitro functional assay to predict in vivo muscle stem cell mediated repair"**. If you use or modify these scripts, we kindly ask you to cite the aforementioned publication.

## Requirements

We have tested these scripts using MATLAB 2019a. You will also need the MATLAB image analysis toolbox, and need a running installation of FIJI to use some of the pre-processing scripts.

Add the folders in the repository to your search path before running any pipeline.

## Contents

* analysis_4x_data: Scripts for analyzing cell coverage in low magnification (4X) images of engineered muscle tissue based on histogram decomposition.

* EdU: Scripts for getting EdU+ (%) cells in 3D images of engineered muscle tissue with sparse content of GFP+ cells.

* nuclei_segmentation: Scripts for performing nuclei segmentation in 3D images of engineered muscle tissue to estimate fusion indices in SAA+ and GFP+ structures.

* aux_functions: Auxilliary functions used by other scripts

  * triangle_th: MATLAB implementation of the triangle method for thresholding images. Code copied from https://www.mathworks.com/matlabcentral/fileexchange/28047-gray-image-thresholding-using-the-triangle-method

  * readTiffStack: Function for loading multi-plane (i.e. 3D) images into MATLAB. Code copied from http://www.matlabtips.com/how-to-load-tiff-stacks-fast-really-fast/

  * hyperellipsoidfit: Function for fitting ellipsoids of any dimension to point clouds. The method is described in *"Direct Least Square Fitting of Hyperellipsoids," IEEE Transactions on Pattern Analysis and Machine Intelligence, DOI: 10.1109/TPAMI.2017.2658574*. Code taken from https://www.mathworks.com/matlabcentral/fileexchange/59186-hyperellipsoidfit and distributed as per the included license agreement.
