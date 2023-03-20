// @Boolean(label="Preview/adjust threshold?", value=1) adjustThrsh
// @File(label="Path to image for analysis:", style="file") sFile
// @Integer(label="Organelle membrane/matrix marker channel for segmentation:", value=1) channelMito
// @Integer(label="Channel to quantify:", value=2) channelMeasure
// @Integer(label="Adaptive threshold size", value=40) adptvThrshSize
// @Integer(label="Adaptive threshold subtract", value=-8) adptvThrshSubtract
// @Integer(label="3D object counter threshold", value=2) objCntr3DThrsh

// Usage: 
// *** requires adaptive threshold plugin:
// https://sites.google.com/site/qingzongtseng/adaptivethreshold
// 1. Run the macro
// * Tick box to test adaptive threshold - 3D object threshold seems fine at 2 with 8-bit thresholded image
// 2. At prompt, draw roi around cells where mitos should be detected and add them to the ROI manager (press key [t])

// Michael Wozny 20230320

// clear log, ROI manager, results and set background black
run("Close All");
print("\\Clear");
roiManager("reset");
run("Clear Results");
setOption("BlackBackground", true);

// open file 'sFile'
open(sFile)
img = getImageID;
run("Set Scale...", "distance=0.00 known=0.00 pixel=1 unit=pixel global")

// if adjustThrsh prompt user to preview thresholds and apply
run("Set Scale...", "global");
if (adjustThrsh==true) {
	previewThrsh();
	// open file 'sFile'
	open(sFile);
	exit
}

// prompt user to draw an ROI around cell to test 3D object counter threshold  
waitForUser("Draw ROIs around cells to analyse, save ROI [t] to manager, OK when finished");
run("Select None");

// select sFile, run adpative threshold 
img = getImageID;
selectImage(img);
setBatchMode(true);
thresholdImage(channelMito);
countMitosPerROI();
print("Done processing... check logs and images after thresholding...");
setBatchMode(false);
close();

// apply adaptive threshold to z-slices of 'channelMito' in open image
function thresholdImage(channelMito) {
	// log the threshold values used
	print("Adaptive threshold size = "+adptvThrshSize);
	print("Adaptive threshold subtract = "+adptvThrshSubtract);
	print("3D object counter threshold = "+objCntr3DThrsh);
    
    // image 'channelMito' -> new stack
    run("Duplicate...", "title=id1 duplicate channels="+channelMito);   
    id1 = getImageID;
	selectImage(id1);
	run("8-bit");
	Stack.getDimensions(width, height, channels, slices, frames);
	
	// run adaptive threshold over stack
	for(z=1; z<slices+1; z++){
  		Stack.setSlice(z);
  		run("adaptiveThr ", "using=[Weighted mean] from="+adptvThrshSize+" then="+adptvThrshSubtract+" stack");
		}
}

function countMitosPerROI() {
	Table.create("Results");
	
	// get path and file name from sFile user input
	sFilePath = substring(sFile,0,lastIndexOf(sFile,File.separator));
	sFileName = substring(sFile,lastIndexOf(sFile,File.separator)+1,lastIndexOf(sFile,"."));
	
	// make a directory to save ROI surface stacks
	ROIsurfDir = sFilePath+File.separator+"ROIs_"+sFileName+File.separator;
  	File.makeDirectory(ROIsurfDir);
  	if (!File.exists(ROIsurfDir)){
      	exit("Unable to create directory");
  		}
  	print("");
  	print(ROIsurfDir);
  	
  	// make a directory to save mito measurements
	mitoMeasDir = sFilePath+File.separator+"mitoMeasurements"+File.separator;
  	File.makeDirectory(mitoMeasDir);
  	if (!File.exists(mitoMeasDir)){
      	exit("Unable to create directory");
  		}
  	print("");
  	print(mitoMeasDir);
	
    // name the original stack img and thresholded 'channelMito' stack id1    
    imgWindow = substring(sFile,lastIndexOf(sFile,File.separator)+1,lastIndexOf(sFile,""));
    selectWindow(imgWindow);
	img = getImageID;
    selectWindow("id1");
    id1 = getImageID;
	
	// save ROI positions
    roiManager("Save",ROIsurfDir+sFileName+"_RoiSet.zip");
    
    // iterate through ROIs in manager
    ROIstartingLength = roiManager("count");
    for (i = 0; i < ROIstartingLength; i++) { 
	
		// select thresholded 'channelMito' stack
		selectImage(id1);
	
		// name the merged surface stack after the ROI id
		roiManager("select", i); 
		ROIname = Roi.getName();
		ROIsurfFileName = sFileName+"_Roi_"+ROIname;
		
		// duplicate ROI as stack
		run("Duplicate...", "title=id2 duplicate channels="+channelMito);   
    	id2 = getImageID;
    	selectImage(id2);
    
    	// clear the area outside the ROI (elimate mitos outside of cell)
    	run("Clear Outside", "stack");
		
		// measure max for each slice to test if mitos are present (bug in 3D object counter breaks macro if not mito detected)
		maxAboveZero = 0;
		print("");
		selectImage(id2);
		Stack.getDimensions(width, height, channels, slices, frames);
		for(z=1; z<slices+1; z++){
  			Stack.setSlice(z);
			getStatistics(area, mean, min, max, std, histogram);
			if (max>0){
				
				print("First max = "+max+" detected at z="+z);
				maxAboveZero = 1;
				break
			}
			else{
				print("max = "+max+" at z="+z);
			}
		}
		
		// if max > 0, detect and track mitos with 3D object counter
		cellsAnalysed = 0;
		mitoResultRow=0;
		if (maxAboveZero==1){
			print(ROIname);
			run("3D Objects Counter", "threshold="+objCntr3DThrsh+" slice=1 min.=25 max.=29030400 surfaces statistics summary");
			
			// merge the 3D object surface, id2 stack, duplicate of original 'channelMito' stack and duplicate of 'channelMeasure'
			selectImage(img);
			roiManager("select", i);
			run("Duplicate...", "title=orig duplicate channels="+channelMito);
			run("8-bit");
			
			selectImage(img);
			roiManager("select", i);
			run("Duplicate...", "title=orig2 duplicate channels="+channelMeasure);
			run("8-bit");
			
			run("Merge Channels...", "c1=orig  c2=orig2 c4=[Surface map of id2] c7=id2 create");
			id3 = getImageID;
			
			// save the merged surface stack as TIFF
			selectImage(id3);
			saveAs("TIFF",ROIsurfDir+ROIsurfFileName);	

    		// add results row for ROI
    		// count the number of rows in statistics from 3D obj counter, corresponds to number of objects/mitos
    		selectWindow("Statistics for id2");
    		n=Table.size("Statistics for id2");
    		print("mitos detected = "+n);
    		if (cellsAnalysed==0) {
    			Table.create("mitoResults");
    			cellsAnalysed = 1;
    		}
    		selectWindow("mitoResults");
    		Table.set("ROI", mitoResultRow, ROIname);
			Table.set("MitoCount", mitoResultRow, n);
			mitoResultRow= mitoResultRow + 1;
			
			// add stats from 3D object counter output ROI 
			selectWindow("Statistics for id2");
			statsFilename = "mitoMeasure_"+sFileName+"_Roi_"+ROIname+".csv";
			statsPath = mitoMeasDir+statsFilename;
			saveAs("mitoResults", statsPath);
			close(statsFilename);
			
			// measure fluorescence intensity of channelMeasure within masks produced from adaptive threshold of channelMito
			selectImage(id3);
			run("Duplicate...", "duplicate");
			id4 = getImageID();
			selectImage(id3);
			Stack.getDimensions(width, height, channels, slices, frames);	

			// use channel 4, the output of adaptive threshold, as a selection area and
			// apply selection to channel 2 of id3 (duplicated from channelMeasure of img) for measure
			maxAboveZeroMeasure = 0;
			print("");
			for(z=1; z<slices+1; z++){
  				selectImage(id3);
				Stack.setSlice(z);
				Stack.setChannel(4)
				// use max to test whether a selection can be made. 
				getStatistics(area, mean, min, max, std, histogram);
				print("max = "+max+" for measure");
				run("Select None");
				if (max>0){
						print("Making selection for measurement at z="+z);
						selectWindow("Results");
						Stack.setSlice(z);
						Stack.setChannel(4);
						run("Create Selection");
						Stack.setChannel(2);
						run("Measure");
						run("Select None");
						continue
				}
				if (max<=0){
					print("max<0 for all z - check if there really is no mito!");
					run("Measure");
					setResult("Area", z-1, 0);
					setResult("Mean", z-1, 0);
					setResult("Perim.", z-1, 0);
					setResult("IntDen", z-1, 0);
					setResult("RawIntDen", z-1, 0);
					setResult("MinThr", z-1, 0);
					setResult("MaxThr", z-1, 0);
					print("max = "+max+" at z="+z+", setting result for measure to 0");
					continue
				}
				else{
					print("on z="+z+", something went wrong!");
				}
			}
			resultsMeasureFilename = "resultsMeasure"+sFileName+"_Roi_"+ROIname+".csv";
			resultsMeasurePath = mitoMeasDir+resultsMeasureFilename;
			selectWindow("Results");
			saveAs("Results", resultsMeasurePath);

			close(resultsMeasureFilename);
			run("Clear Results");
			
			// make a duplicate to subtract the background from
			run("Select None");
			selectImage(id3);
			run("Duplicate...", "title=c2Intact duplicate channels=2"); 
			run("Select None");

			// measure within ROI/outside of mitochondria (cell background) 
			for(z=1; z<slices+1; z++){
				selectImage(id3);
				Stack.setSlice(z);
				Stack.setChannel(4)
				// use max to test whether a selection can be made. 
				getStatistics(area, mean, min, max, std, histogram);
				print("max = "+max+" for background");
				run("Select None");
				if (max>0){
					// copy ROI to id3 for measurement
					selectImage(img);
					roiManager("select", i);
					run("Copy");
					selectImage(id3);
					run("Restore Selection");
					roiManager("Add");
					run("Select None");
					// keep the index of the current ROI as currIdxROI
					roiManager("Select", roiManager("count")-1);
					currIdxROI = roiManager("index");
					Roi.setPosition(channelMeasure, z, 1);
					selectImage(id3);
					Stack.setSlice(z);
					Stack.setChannel(channelMeasure);
					run("Select None");
					
					// use channel 4, the output of adaptive threshold, as a selection area
					selectImage(id3);
					Stack.setSlice(z);
					Stack.setChannel(4);
					run("Create Selection");
					Stack.setChannel(channelMeasure);
					roiManager("Add");
					roiManager("Select", roiManager("count")-1);
					mitoIdxROI = roiManager("index");
					Roi.setPosition(channelMeasure, z, 1);

					// combine currIdxROI & mitoIdxROI
					run("Select None");
					roiManager("Select", newArray(currIdxROI,mitoIdxROI));
					selectImage(id3);
					Stack.setSlice(z);
					Stack.setChannel(channelMeasure);
					roiManager("XOR");
					selectImage(id3);
					Stack.setSlice(z);
					Stack.setChannel(channelMeasure);
					roiManager("Add");
					roiManager("Select", roiManager("count")-1); 
					bckgrndIdxROI = roiManager("index");
					Roi.setPosition(channelMeasure, z, 1);
					selectImage(id3);
					Stack.setSlice(z);
					Stack.setChannel(channelMeasure);
					run("Clear Outside", "slice");
					run("Measure");
					run("Select None");
					continue
				}
				if (max<=0){
					print("max<0 for all z - check if there really is no mito!");
					selectWindow("Results");
					run("Measure");
					setResult("Area", z-1, 0);
					setResult("Mean", z-1, 0);
					setResult("Perim.", z-1, 0);
					setResult("IntDen", z-1, 0);
					setResult("RawIntDen", z-1, 0);
					setResult("MinThr", z-1, 0);
					setResult("MaxThr", z-1, 0);
					print("max = "+max+" at z="+z+", setting result for background to 0");
					continue
				}
				else{
					print("on z="+z+", no thresholded image max!");
				}
			}
			
			selectWindow("Results");
			resultsBackgroundFilename = "resultsBackground"+sFileName+"_Roi_"+ROIname+".csv";
			resultsBackgroundPath = mitoMeasDir+resultsBackgroundFilename;
			saveAs("Results", resultsBackgroundPath);
			close(resultsBackgroundFilename);
			run("Clear Results");
			run("Select None");
			
			// duplicate channels for merge
			selectImage(id3);
			run("Duplicate...", "title=c1 duplicate channels=1");
			selectImage(id3);
			run("Duplicate...", "title=c2 duplicate channels=2");
			selectImage(id3);
			run("Duplicate...", "title=c3 duplicate channels=3");
			selectImage(id3);
			run("Duplicate...", "title=c4 duplicate channels=4");

			run("8-bit");
			run("Merge Channels...", "c1=c1 c2=c2Intact c7=c3 c4=c4 c5=c2 create"); 
			saveAs("TIFF",ROIsurfDir+ROIsurfFileName+"_background");	
			close();
			
			// close open images
			
			selectImage(id3);
			close();
			selectImage(id4);
			close();
		}

		// if max is not > 0, report, save thresholded stack without surface stack and record ROI mito count = 0   
		else{
			print(ROIname);
			print("max<0 for all z - check if there really is no mito!");
			
			// merge the thresholded stack with original and save as TIFF
			selectImage(img);
			roiManager("select", i);
			run("Duplicate...", "title=orig duplicate channels="+channelMito);
			run("8-bit");
			run("Merge Channels...", "c1=orig c7=id2 create");
			id3 = getImageID;
			selectImage(id3);
			saveAs("TIFF",ROIsurfDir+ROIsurfFileName);

			// add results row for ROI
    		// no mito detected - add a 0 to mito count for ROI
    		n = 0;
    		print("mitos detected = "+n);
    		if (cellsAnalysed==0) {
    			Table.create("mitoResults");
    			cellsAnalysed = 1;
    		}
    		selectWindow("Results");
    		setResult("ROI", mitoResultRow, ROIname);
			setResult("MitoCount", mitoResultRow, n);
			mitoResultRow= mitoResultRow + 1;
			selectImage(id3);
			close();
		}
    }
    csvDir = sFilePath+File.separator+"mitoCounts"+File.separator;
    File.makeDirectory(csvDir);
    if (!File.exists(csvDir)){
      	exit("Unable to create directory");
    }
  	print("");
  	print(csvDir);

  	// save mito count for ROIs
  	selectWindow("mitoResults");
    saveAs("Results",csvDir+"mitoCount_"+sFileName+".csv");
    close("mitoResults");

    // save background ROI positions
    roiManager("Save",ROIsurfDir+sFileName+"_RoiSetBackgroundROIs.zip");

	// save log
	selectWindow("Log"); 
	saveAs("Text",ROIsurfDir+sFileName+"_log.txt");
}

function previewThrsh() {
	// image 'channelMito' -> new stack
    run("Duplicate...", "title=id1 duplicate channels="+channelMito);   
    id1 = getImageID;
	selectImage(id1);
	run("8-bit");
	
	// prompt the user to open adaptive threshold plugin and test preview on 8 bit 'channelMito' stack
	waitForUser("Open Plugins>adaptiveThreshold plugin and test threshold with preview, remember size & subtract");
	close();
	Dialog.create("Best adaptive threshold values?");
	Dialog.addNumber("Best adaptive threshold size:", adptvThrshSize);
	Dialog.addNumber("Best adaptive threshold subtract:", adptvThrshSubtract);
	Dialog.show();

	// prompt user to draw an ROI around cell to test 3D object counter threshold  
	waitForUser("Draw an ROI around a cell, save ROI [t] to manager");
	run("Select None");
	
	// run adaptive threshold
	img = getImageID;
	selectImage(img);
	run("Duplicate...", "title=id1 duplicate channels="+channelMito); 
	id1 = getImageID;
	selectImage(id1);
	run("8-bit");
	Stack.getDimensions(width, height, channels, slices, frames);
	for(z=1; z<slices+1; z++){
  		Stack.setSlice(z);
  		run("adaptiveThr ", "using=[Weighted mean] from="+adptvThrshSize+" then="+adptvThrshSubtract+" stack");
		}
	
	// run 3D object counter on ROI
	selectImage(id1);
	
	// select ROI and duplicate as stack
	roiManager("select", 0); 
	run("Duplicate...", "title=id2 duplicate channels="+channelMito);
	id2 = getImageID;
	selectImage(id2);
	// clear the area outside the ROI (elimate mitos outside of cell)
    run("Clear Outside", "stack");

	// prompt user to open 3D Object Counter and preview threshold
	waitForUser("Open Analyze>3D Object Counter, remember threshold");
	Dialog.create("Best 3D Object Counter threshold?");
	Dialog.addNumber("Enter threshold:", objCntr3DThrsh);
	Dialog.show();
	close("Statistics for id2");
	run("Close All");
}