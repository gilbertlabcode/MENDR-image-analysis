/* Macro Name: Max_project Batch Converter
 *  Author: Amir Reza Meysami Fard. Modified by Jose Cadavid to threshold images right after conversion. Keeps projected images regardless of thresholding method
 *  Date of Last Revision: August 11th, 2019
 *  Input: Folder containing OIF files and their supporting folders, and (or) OIB files
 *  Output: Max_projects of images
 */


/*To fully automate this macro, do the following:
 * 1- Open an oib file, make sure split channels is checked, then click open
 * 2- Close the file
 * 3- Run the following: plugin -> Bio-Formats -> Bio-Formats Windowless Import, choose an oif file
*/

run("Set Measurements...", "area mean standard min shape area_fraction display redirect=None decimal=3");

setBatchMode(true);
//Prompts user to select an input folder containing the files to be processed
path_parent = getDirectory("Choose Input Directory containing images from drugs");
//Extracts the folder name
parent_name = File.getName(path_parent);
//Extracts the name of the parent folder (folder containing the input folder)
path_grandparent = File.getParent(path_parent);
//Identifies the type of slash used in the file path. Essentially auto determine PC or MAC
separator = File.separator();
path_parent_array = getFileList(path_parent);
fileTitle = File.getName(path_parent);
//Generates the file path for the output folder based on input folder name
outputPath = path_grandparent + separator + parent_name + "_projected";
//Creates output folder in the Parent folder
File.makeDirectory(outputPath);

//Create file where results are stored
f=File.open(outputPath+"-results Max Projection.txt");
//Print the first line as a dummy identifier (done so the file sorter can know where to start and stop when converting the txt file to an excel table)
print(f,"Threshold results");
//Print headers
print(f,"Drug" + "\t" + "Replicate" + "\t" + "Image name" + "\t" + "% Coverage");

for(j=0; j<path_parent_array.length; j++) {
	drug_name=File.getName(path_parent_array[j]);
	path_drug = path_parent + separator + File.getName(path_parent_array[j]);
	path_drug_array = getFileList(path_drug);
	output_drug = outputPath + separator +  File.getName(path_parent_array[j]);
	File.makeDirectory(output_drug);

	for (p = 0; p<path_drug_array.length; p++){
		rep_name=File.getName(path_drug_array[p]);
		path_rep = path_drug + separator + File.getName(path_drug_array[p]);
		path_rep_array = getFileList(path_rep);
		output_rep = output_drug + separator + File.getName(path_drug_array[p]);
		File.makeDirectory(output_rep);

		for (q = 0; q<path_rep_array.length; q++){
			fileTitle = File.getName(path_rep_array[q]);							//Gets title of oif file

			if (endsWith(fileTitle,".oif") == 1 || endsWith(fileTitle,".oib") == 1){
				open(path_rep + separator + fileTitle); 									//Generally opens 3 stacks, last one being the nuclei (if the order is GFP, SAA, and DRAQ5)
				//Operate on channel 1
				selectImage(1);
				run("8-bit");
				run("Z Project...", "projection=[Max Intensity]");
				fileTitle = replace(fileTitle,".oif","");								// deletes .oif for future use
				fileTitle = replace(fileTitle,".oib","");								// deletes .oib for future use
				saveAs("tif", output_rep + separator + fileTitle + ".tif");
				//Now that the image is saved, we can threshold it with the triangle method and get the %coverage

				setOption("BlackBackground", true);
				setAutoThreshold("Triangle dark no-reset");
				run("Convert to Mask");
				//Measure and set area fraction as output
				List.setMeasurements;
				A_frac=List.getValue("%Area");
    			print(f,drug_name + "\t" + rep_name + "\t" + fileTitle + "\t" + A_frac);
				//Close all images
				close("*");
			}
		}
	}
}

print("  Done!!!");			//Let's user know it's done
