setBatchMode(true);															//

path = getDirectory("Choose Input Directory");								//Prompts user to select an input folder containing the files to be processed
inputFolder = File.getName(path);											//Extracts the folder name
parentPath = File.getParent(path);											//Extracts the name of the parent folder (folder containing the input folder)
separator = File.separator();												//Identifies the type of slash used in the file path. Essentially auto determine PC or MAC															
folderArray = getFileList(path);	
fileTitle = File.getName(path);	
outputPath = parentPath + separator + fileTitle + "_separated";				//Generates the file path for the output folder based on input folder name
File.makeDirectory(outputPath);												//Creates output folder in the Parent folder

for(j=0; j<folderArray.length; j++) {	
	fileTitle = File.getName(folderArray[j]);								//Gets title of oif file
	if (endsWith(fileTitle,".oif") == 1 || endsWith(fileTitle,".oib") == 1){
	open(path + separator + fileTitle); 									//Generally opens 3 stacks, last one being the nuclei

	fileTitle = replace(fileTitle,".oif","");								// deletes .oif for future use
	fileTitle = replace(fileTitle,".oib","");								// deletes .oib for future use
	
	outputPath2 = outputPath + separator + fileTitle;
	File.makeDirectory(outputPath2);
	
	outputPathNuclei = outputPath2 + separator + "nuclei.tif";
	saveAs("Tiff", outputPathNuclei);
	close();
	outputPathGFP = outputPath2 + separator + "GFP.tif";
	saveAs("Tiff", outputPathGFP);
	close();
	}
}
