//Close all open images
run("Close All");
//Select the directory (each folder for each sample will have three images- with multiple planes- called GFP, red and nuclei; we do not care about red for now)
dir = getDirectory('Choose Directory of images');

//Start batch mode
setBatchMode(true);

//////////////////////////////////////////////////////////////////////////////////////////
//Part 1: Open images, blur them and stack each channel for ease of handling
//////////////////////////////////////////////////////////////////////////////////////////

//Open cytoplasm (GFP) multi-plane image

open(dir+"GFP.tif");
rename("Cytoplasm");

//Open images in the multi-plane image
open(dir+"nuclei.tif");
rename("Nuclei");

//Get max projections of each channel (used for display in the end)
selectWindow("Cytoplasm");
run("Z Project...", "projection=[Max Intensity]");
rename("MAX_cytoplasm_original");
run("8-bit");

selectWindow("Nuclei");
run("Z Project...", "projection=[Max Intensity]");
rename("MAX_nuclei_original");
run("8-bit");

//Merge the original nuclei and cytoplasm projections for visualization purposes in the end
run("Merge Channels...", "c2=MAX_cytoplasm_original c1=MAX_nuclei_original create keep ignore");
run("RGB Color"); //Convert to rgb
rename("Original");
selectWindow("MAX_nuclei_original");//Close this channel, not needed anymore
close();

//Blur stack: Median filter of radius 1 for nuclei (determined to provide enough smoothing without deforming the nuclei too much - keeps circular shape)
selectWindow("Nuclei");
run("Median...", "radius=1 stack");

//Blur stack: Median filter of radius 5 for cytoplasm (determined to provide enough smoothing without blurring out thin fibers)
selectWindow("Cytoplasm");
run("Median...", "radius=5 stack");
//Find edges of cytoplasms (for ease of segmentation)
run("Duplicate...","title=Cytoplasm_edges duplicate");
run("Find Edges", "stack");
run("Z Project...", "projection=[Max Intensity]");
run("Enhance Contrast...", "saturated=0.3 normalize equalize");
run("8-bit");
selectWindow("Cytoplasm_edges");
close();

//////////////////////////////////////////////////////////////////////////////////////////
//Part 2: Binarize cytoplasm 
//////////////////////////////////////////////////////////////////////////////////////////
selectWindow("Cytoplasm");
//setAutoThreshold("Otsu dark stack");
setOption("BlackBackground", true);
setThreshold(350,65535);
run("Convert to Mask","stack");
run("Green");

//Fill holes smaller than 1000 pixels (approx 1.5 times the average nucleus of 680 pixels)
run("Duplicate...", "title=Threshold1 duplicate"); //Duplicate to get the mask
run("Invert", "stack");//Invert so the holes become white
run("Analyze Particles...", "size=0-1000 pixel show=Masks clear stack"); //Get particles of desired size
imageCalculator("OR create stack", "Cytoplasm","Mask of Threshold1");
rename("Cytoplasm_filled");
selectWindow("Threshold1");
close();
selectWindow("Mask of Threshold1");
close();
selectWindow("Cytoplasm");
close(); //Close unfilled cytoplasm image

//////////////////////////////////////////////////////////////////////////////////////////
//Part 3: Eliminate nuclei outside of selected cytoplasm regions and threshold
//////////////////////////////////////////////////////////////////////////////////////////

//Threshold nuclei
selectWindow("Nuclei");
setOption("BlackBackground", true);
setThreshold(300,65535);
run("Convert to Mask","stack");
run("Blue");

//Perform "AND" operation with cytoplasm channel (remove nuclei outside of the GFP cells)
imageCalculator("AND create stack", "Nuclei","Cytoplasm_filled");
rename("Nuclei_excluded");

//////////////////////////////////////////////////////////////////////////////////////////
//Part 4: Stack images and produce output
//////////////////////////////////////////////////////////////////////////////////////////

//Z- project cytoplasm
selectWindow("Cytoplasm_filled"); //Cytoplasm
run("Z Project...", "projection=[Max Intensity]");
rename("MAX_cytoplasm");
selectWindow("Cytoplasm_filled");
close(); //Close cytoplasm stack 

//Z-project nuclei and remove artifacts (watershed to try to separate touching nuclei, remove particles)
selectWindow("Nuclei_excluded"); //Nuclei
run("Z Project...", "projection=[Max Intensity]");
run("Watershed");
//Eliminate small artifacts in nuclei channel (smaller than 170 pixels = 0.25*Mean size of nuclei); this includes nuclei objects that are only partially in a fiber
run("Analyze Particles...", "size=0-170 pixel show=Masks clear stack"); //Remove objects smaller than 340 pixels, can be changed
imageCalculator("XOR stack", "MAX_Nuclei_excluded","Mask of MAX_Nuclei_excluded");
rename("MAX_nuclei");
selectWindow("Mask of MAX_Nuclei_excluded");
close();
selectWindow("Nuclei_excluded");
close(); //Close cytoplasm stack 

//Merge the original projection with the edges detected (for easy segmentation)
run("Merge Channels...", "c2=MAX_Cytoplasm_edges c3=MAX_cytoplasm_original create ignore");
run("RGB Color"); //Convert to rgb
rename("Edges");
selectWindow("Composite"); //Close the composite
run("Close");

//Close remaining image
selectWindow("Nuclei");
close(); //Close original nuclei image (the one before eliminating nuclei outside of the fibers)

//Stack the rgb images
run("Images to Stack", "name=Stack_max title=[] use"); //Stack!
selectWindow("Stack_max");

//Exit batch mode
setBatchMode(false);

//Save the produced thresholded image so that it can be opened without running this macro again. Save in the parent directory
f_name=File.getName(dir);
parent=File.getParent(dir);
parent_name=File.getName(parent); 
//saveAs("Tiff", dir+f_name+".tif"); save in directory
saveAs("Tiff", parent+File.separator+parent_name+f_name+".tif"); // save in parent directory
run("Close All");
