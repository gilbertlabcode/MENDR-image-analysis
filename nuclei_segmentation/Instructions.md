# Instructions for analysis of fusion index and nuclei segmentation

## Description

Note: This algorithm is slow, and takes about 5-10 minutes per image. The kmeans algorithm uses a random initialization, so it is possible that multiple runs of the pipeline for the same image result in slightly different results (less than 1% difference).

Steps to use the scripts in this folder to analyze images (40x magnification) of engineered muscle templates to calculate fusion indices. The images must include a nuclear channel and at least a GFP or SAA channel. The scripts assume that the master folder to be analyzed has the following structure:

Master -
  Treatments (e.g. drugs tested) -
    Reps -
      Image files

While your folders might be named differently, we assume this nested structure. Image files to be processed are typically .oib or .oif files saved by a confocal microscope.

The scripts will first create anisotropic 3D single-channel images using ImageJ by processing a batch of images in your folder. Then, a MATLAB function can be run on those images to segment nuclei and classify them as GFP+, SAA+ or double-positive.

## Steps

1. Resample your images to get anisotropic 3D images:

  Open ImageJ (FIJI) and run either the "Reslicer_fusionindex" script. This script will batch-process a nested folder with images (see above) and create a new folder with anisotropic 3D images (single-channel) for the SAA, GFP, and nuclear channel. It also creates a text file with information about the images, including the pixel size in microns (necessary for the following step). See the script for more detailed instructions.

2. Analyze images with MATLAB:

  Open MATLAB and make sure the directory with all the scripts is in your search path, or that it is your current directory. Run the script "FusionLooper"; it will prompt you to select the folder with the resliced images produced by ImageJ, as well as the text file exported by ImageJ. The script will export a new text file with the information of the images, along with total nuclei count, GFP nuclei, SAA nuclei, and double-positive nuclei, for each image .

  The script "FusionLooper" is a wrapper to automate batch analysis of multiple images. Open this script for more details on how each image is analyzed.

  Note: All of the other scripts in this folder all called internally by the segmentation routine, and the user does not have to modify them. The exception is the "segmentationParametersNuclei", which stores the parameters used for the algorithm.
