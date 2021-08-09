# Instructions for analysis of EdU data

## Description

Steps to use the scripts in this folder to analyze images (20x magnification) of engineered muscle templates that have been stained for EdU+ nuclei. The scripts assume that the master folder to be analyzed has the following structure:

Master -
  Treatments (e.g. drugs tested) -
    Reps -
      Image files

While your folders might be named differently, we assume this nested structure. Image files to be processed are typically .oib or .oif files saved by a confocal microscope.

The scripts will first create anisotropic 3D single-channel images using ImageJ by processing a batch of images in your folder. Then, a MATLAB function can be run on those images to count the number of GFP+ objects in them, as well as the number of GFP+ objects that contain EdU+ nuclei.

## Steps

1. Resample your images to get anisotropic 3D images:

  Open ImageJ (FIJI) and run either the "Reslicer_EdU" script. This script will batch-process a nested folder with images (see above) and create a new folder with anisotropic 3D images (single-channel) for the EdU, GFP, and nuclear channel. It also creates a text file with information about the images, including the pixel size in microns (necessary for the following step). See the script for more detailed instructions.

2. Analyze images with MATLAB:

  Open MATLAB and make sure the directory with all the scripts is in your search path, or that it is your current directory. Run the script "EdUGFPLooper"; it will prompt you to select the folder with the resliced images produced by ImageJ, as well as the text file exported by ImageJ. The script will export a new text file with the GFP+EdU+ objects (%).

  The script "EdUGFPLooper" is a wrapper to automate batch analysis of multiple images. Open this script for more details on how each image is analyzed.
