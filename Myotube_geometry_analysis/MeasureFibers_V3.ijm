//Get directory to save files (text file)
dir = getDirectory("Choose directory of Input Images");
parentDir = File.getParent(dir);	
//Get list of files in each folder
list = getFileList(dir);
separator = File.separator();	
fileTitle = File.getName(dir);	
dir_out = getDirectory("Choose directory of Outputs");
run("ROI Manager...");
roiManager("Show All with labels");
setTool("polygon");
//Create table where results are stored
title = "results";
f = "["+title+"]"; 
run("New... ", "name="+f+" type=Table"); 
print(f,"\\Headings: File Name\tLength (micro m)\tMean fiber width (micro m)\tStd dev of width (micro m)\tNuclei count");

for(j=0; j<list.length; j++) {
	f_name=File.getName(dir+separator+list[j]);
	//open the file
	open(dir+separator+list[j]);
	//wait for user to select fibers
	waitForUser("Please mark 10 fibers from all sides of the image, each fiber needs to include 3 nuclei.");
	roiManager("Save", dir_out + separator + f_name+ ".zip");
	n_roi=roiManager("count");
	//Rename: When the stack_max is saved, it is saved as Stack_max.tif (or any other name), so we have to change the name in ImageJ for consistency
	rename("Stack_max"); 
	//Start batch mode
	setBatchMode(true);

	////Loop for all ROIs in the ROI manager

	for (i=0;i<n_roi;i++){

		//////////////////////////////////////////////////////////////////////////////////////////
		//Part 1: Select ROI and apply it to get section (fiber) of interest (manual segmentation)
		//////////////////////////////////////////////////////////////////////////////////////////

		//Go to main image
		selectWindow("Stack_max");
		//Select ROI
		roiManager("Select",i);
		//Duplicate subimage and make mask
		run("Duplicate...", "title=Stack_max-1 duplicate");
		run("Create Mask");
		//Apply mask
		imageCalculator("AND stack", "Stack_max-1","Mask");
		selectWindow("Mask");
		close();
		//Break into images again:  Break the duplicated stack and unnecesary channels (edges and original)
		selectWindow("Stack_max-1");
		run("Stack to Images");
		selectWindow("Edges");
		close();
		selectWindow("Original");
		close();

		//////////////////////////////////////////////////////////////////////////////////////////
		//Part 2: Analysis of the cytoplasm channel
		//////////////////////////////////////////////////////////////////////////////////////////

		//// Pre-processing of the segmented fibre

		//Skeletonize cytoplasm, prune small branches, get distance transform and apply selection of skeleton as ROI in distance transform. 
		//Measure mean value (mean distance/2) and std value (std/2), and length (area)
		selectWindow("MAX_cytoplasm");
		//Upon destacking the images are RGB, so we convert them back to 8-bit
		run("8-bit");
		//Since the image was green for convenience, converting it to RGB results in a grey tone that must me made
		//white for processing (fully binarized image), so we threshold it with a threshold of 1.
		setOption("BlackBackground", true);
		setThreshold(1,255);
		run("Convert to Mask");
		//Fill holes
		run("Fill Holes");
		//Eliminate spurious objects resulting from imperfect ROI selection (objects below 300 pixels - adjustable)
		run("Analyze Particles...", "size=300-Infinity pixel show=Masks"); //We only select objects larjer than 300 pixels - the main fibre
		//Apply this mask to eliminate the spurious objects, and then close it
		imageCalculator("AND", "MAX_cytoplasm","Mask of MAX_cytoplasm");
		selectWindow("Mask of MAX_cytoplasm");
		close();

		//// Euclidean distance transform: For every white pixel, get the euclidean distance (in pixels) to the closes black pixels. This equals half width
		//// get length and width of the fiber

		run("Duplicate...", "title=Dist");
		run("Distance Map");
		//Skeletonize
		selectWindow("MAX_cytoplasm");
		run("Skeletonize (2D/3D)");
		//Prune small side branches. The parameter (20 micro m) is the length of branches to be pruned
		run("Pruning ", "threshold=20.0");

		//Get selection of skeleton
		run("Invert"); //Don't know why, but the selection works with the inverted image
		run("Create Selection");
		//Measure in the distance map: Apply the selection corresponding to the skeleton
		selectWindow("Dist");
		run("Restore Selection");
		//Get measurements (does not open the results window, but it is equivalent to run("measure")
		List.setMeasurements;
		//Close distance map
		close();

		//Retrieve results and process them to convert to micro meters according to the conversion factor for 40x images given by fiji
		Length=List.getValue("Area")/0.497;//0.497 is the conversion factor (1 pixel = 0.497 micro meters) for 40 x images, given by Fiji. The area of the selection equals the length of the skeleton times 1 pixel, so this is converted accordingly
		Width=List.getValue("Mean")*0.497*2;//Again using the conversion factor, the distance transform yiels results in pixel distance, so we convert it to micro m and multiply it by two (the distance transform gives half width)
		S_Width=List.getValue("StdDev")*0.497*2;//Same as above for the standard deviation of the width

		//////////////////////////////////////////////////////////////////////////////////////////
		//Part 3: Analysis of the nuclei channel
		//////////////////////////////////////////////////////////////////////////////////////////

		//Count nuclei in Cytoplasm ROI: Distance transform and find maxima
		selectWindow("MAX_nuclei");
		//Upon destacking the images are RGB, so we convert them back to 8-bit
		run("8-bit");
		//Since the image was blue for convenience, converting it to RGB results in a grey tone that must me made
		//white for processing (fully binarized image), so we threshold it with a threshold of 1.
		setOption("BlackBackground", true);
		setThreshold(1,255);
		run("Convert to Mask");
		//Remove speckles or nuclei that are barely in the fibre: The mean fibre size is 680, so the blob must be at least 170 pixels to be considered in the fibre
		//this is done because if the segmentation is not perfect, parts of random nuclei will be included in the mask, which we then eliminate here
		run("Analyze Particles...", "size=170-Infinity pixel show=Masks");
		//Apply this mask to eliminate the spurious objects, and then close it
		imageCalculator("AND", "MAX_nuclei","Mask of MAX_nuclei");
		selectWindow("Mask of MAX_nuclei");
		close();

		//The watershed segmentation can separate some overlapping nuclei better than the distance map, depending on the shape of the blob
		run("Watershed"); 
		//Euclidean distance transform (used to find maxima (count objects and their positions)
		run("Distance Map");
		//Find maxima (noise tolerance 1, since the distance map is not noisy)
		run("Find Maxima...", "noise=1 output=Count");
		//Close distance map
		close();

		//Retrieve results: Number of nuclei in the count
		N = getResult("Count"); 
		selectWindow("Results");
		run("Close");

		//Close all secondary images (other images, leave the main one open)
		selectWindow("Stack_max");
		close("\\Others");
		//Print results in the current txt file
		print(f, f_name + "\t" + Length + "\t"+ Width + "\t" + S_Width + "\t" + N );
	}
	i = 0;
	while (i < n_roi){
		roiManager("Select",i);
		roiManager("delete");
		n_roi = n_roi - 1;
	}
	setBatchMode(false);
	selectWindow("ROI Manager");
	run("Close All");
	//Loop finished, exit batch mode!
}
selectWindow(title);
saveAs("results", dir_out + separator + File.getName(dir) + ".txt");
selectWindow(title);
run("Close");
close("Roi Manager");
