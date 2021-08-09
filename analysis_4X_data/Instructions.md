# Instructions for analysis of fiber coverage

## Description

Steps to use the scripts in this folder to analyze images (4x magnification) of engineered muscle templates. The scripts assume that the master folder to be analyzed has the following structure:

Master -
  Treatments (e.g. drugs tested) -
    Reps -
      Image files

While your folders might be named differently, we assume this nested structure. Image files to be processed are typically .oib or .oif files saved by a confocal microscope.

The scripts will first create maximum projections of your images using ImageJ, and will threshold them and store the percentage of pixels (i.e. coverage) above threshold using the triangle method. The results and max projected images will be stored.

Then, the MATLAB functions will loop through the projected images and threshold them using a histogram decomposition method and also produce a text file with the results. For our manuscript, we found that a combination of both thresholding methods worked well better than a single one. The triangle method worked better with images that had a low fiber coverage, while the MATLAB method worked better for images with more fiber coverage.

## Steps

1. Create max projections of the channel to be analyzed:

  Open ImageJ (FIJI) and run either the "Max_project threshold triangle GFP" or "Max_project threshold triangle SAA" scripts depending on which channel of the images you want to analyze. These scripts will save a text file with results from the triangle method, and will store the projected images in a directory called "projected" under the main directory you processed.

2. Analyze max projections with MATLAB:

  Open MATLAB and make sure the directory with the scripts is in your search path, or that it is your current directory. Run the "ThresholdFolderLooper" script by either opening it and clicking "Run", or by typing ThresholdFolderLooper in the command line. You will be prompted to select a folder where the images are contained; choose the folder where the max projections from step 1 were saved. The MATLAB script will save a text file with the results in the chosen directory.

3. Put together data:

  You can import both text files from steps 1 and 2 into an EXCEL sheet or similar. To select the most appropriate results for each image, use the "flag" column of the MATLAB text file. If flag = TRUE, the MATLAB thresholding was not reliable and you should choose the results from the triangle method (FIJI). If flag = 0. choose the area coverage from the MATLAB functions.
