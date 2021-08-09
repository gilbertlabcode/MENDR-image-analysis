/* Reslicer
* Author: Amir Meysami Fard. 
Updated feb 25 2020
*/

/*To fully automate this macro, do the following:
 * 1- Open an oib file, make sure split channels is NOT checked, then click open
 * 2- Close the file
 * 3- Run the following: plugin -> Bio-Formats -> Bio-Formats Windowless Import, choose an oif file 
*/

run("Set Measurements...", "area mean standard min shape area_fraction display redirect=None decimal=3");

setBatchMode(true);															
//Promts user to select an input folder containing the files to be processed
path_parent = getDirectory("Choose Input Directory containing replicates");		
//Extracts the folder name						
parent_name = File.getName(path_parent);	
//Extracts the name of the parent folder (folder containing the input folder)										
path_grandparent = File.getParent(path_parent);		
//Identifies the type of slash used in the file path. Essentially auto determine PC or MAC									
separator = File.separator();																											
path_parent_array = getFileList(path_parent);	
fileTitle = File.getName(path_parent);	
//Generates the file path for the output folder based on input folder name
outputPath = path_grandparent + separator + parent_name + "_resliced";		
//Creates output folder in the Parent folder		
File.makeDirectory(outputPath);												

//Create file where results are stored
f=File.open(outputPath+"-rStacks.txt");
//Print headers
print(f,"Replicate \t Time point \t Image \t NewID \t \width(um) \t Height(um) \t Depth(um)");

//Image counter
ii=0;

for(j=0; j<path_parent_array.length; j++) {	 // loop through replicates
	drug_name=File.getName(path_parent_array[j]);
	path_drug = path_parent + separator + File.getName(path_parent_array[j]);
	path_drug_array = getFileList(path_drug);
	output_drug = outputPath + separator +  File.getName(path_parent_array[j]);
	
	for (p = 0; p<path_drug_array.length; p++){ //loop through timepoints
		rep_name=File.getName(path_drug_array[p]);
		path_rep = path_drug + separator + File.getName(path_drug_array[p]);
		path_rep_array = getFileList(path_rep);
		output_rep = output_drug + separator + File.getName(path_drug_array[p]);
		
		for (q = 0; q<path_rep_array.length; q++){ //loop through images
			fileTitle = File.getName(path_rep_array[q]);							//Gets title of oif file
			
			if (endsWith(fileTitle,".oif") == 1 || endsWith(fileTitle,".oib") == 1){
				//Open image
				open(path_rep + separator + fileTitle); 								
				//Do not split channels; reslice together and then split
				im_name = replace(fileTitle,".oif","");								// deletes .oif for future use
				im_name = replace(fileTitle,".oib","");								// deletes .oib for future use
				ii=ii+1;
				//Convert to 8 bit
				run("8-bit");
				getVoxelSize(width, height, depth, unit);

				run("Reslice [/]...", "output=&width start=Top");
				selectImage(2);
				run("Reslice [/]...", "output=&width start=Top");
				selectImage(3);
				
				//Split channels
				run("Split Channels");

				if (nImages==4){
					//If there are four open images, there was only SAA and Nuclei, so SAA is 3 and Nuclei is 4
					
					//Im3 is SAA
					selectImage(3);
					saveAs("tif", outputPath + separator + "SAA_reslice"+ii+".tif");
					//Im4 is nuc
					selectImage(4);
					saveAs("tif", outputPath + separator + "Nuc_reslice"+ii+".tif");
				
				} else {
					//There are 5 images, so images 3 is GFP, 4 is SAA, 5 is nucle
					//Im3 is GFP
					selectImage(3);
					saveAs("tif", outputPath + separator + "GFP_reslice"+ii+".tif");
					//Im4 is SAA
					selectImage(4);
					saveAs("tif", outputPath + separator + "SAA_reslice"+ii+".tif");
					//Im5 is nuc
					selectImage(5);
					saveAs("tif", outputPath + separator + "Nuc_reslice"+ii+".tif");
					
				}

				print(f,drug_name+"\t"+rep_name+"\t"+im_name+"\t"+d2s(ii,0)+"\t"+width+"\t"+height+"\t"+depth);
				//Close all images
				close("*");

				
			}
		
		}
	
	}

}

