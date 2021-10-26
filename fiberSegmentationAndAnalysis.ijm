// ===============================
// AUTHOR : Gisele Miranda
// IVESTIGATORS : Fredrik Bäcklund, Anna Rising
// CREATE DATE : 2020 - 08 - 12
// PURPOSE : An image based analysis of recombinant silk fiber structures
// NOTES : Fiji required - https://fiji.sc
// ===============================

/************************** parameters **************************/
minDist = 2; // minimum diameter of the fiber - used to calculate histogram bins
maxDist = 50; // maximum diameter of the fiber - used to calculate histogram bins
pixPerMic = 2.04081633; // equivalent in pixels to 1 micron 
scale = true; // whether to export the measures in microns or pixels

// thresholding settings for the segmentation of the BF image
thresholdBF_automatic = "Yen"; // automatic threshold
thresholdBF_manual = 103; // manual threshold
useThresholdBF_manual = false; // whether to use the manual threshold of not

// thresholding settings for the segmentation of the POM image
thresholdPOM_automatic = "Yen"; // the threshold method can be different from thresholdBF_automatic
thresholdPOM_manual = 20;
useThresholdPOM_manual = false;

// paramenters of the histograms of area and intensity values of the segmented particles in the POM image 
init_area = 0;
end_area = 20000;
bin_area = 50;
init_intensity = 20;
end_intensity = 255;
bin_intensity = 5;

// settings for diagonal fibers
diagonalFiber = true; // whether fibers are diagonal or not
nPixelsCorner = 10; // number of pixels that should be eliminated from the fiber if endings are in the corners and oriented 45 degrees
nPixelsSide = 20; // number of pixels that should be eliminated from the fiber if endings are not in the corners and oriented 45 degrees
/****************************************************************/

// read input directory
path = getDirectory("Choose a Directory");
dir = getFileList(path);

// create output directory
outPath = path + "segmentation/";
File.makeDirectory(outPath);

// initialize buffers to save measures and histograms
bufferMeasuresBF = "File;Area;Mean;StdDev;Max;Min;Perim;avgDiam1;stdDiam1;minDiam1;maxDiam1;avgDiam2;stdDiam2;minDiam2;maxDiam2;avgDiam;stdDiam;thresholdValue\n";
bufferMeasuresPOM = "File;Area;Mean;StdDev;Max;Min;Perim;avgDiam1;stdDiam1;minDiam1;maxDiam1;avgDiam2;stdDiam2;minDiam2;maxDiam2;avgDiam;stdDiam;thresholdValue\n";
bufferHistBF = "bin;";
bufferHistPOM = "bin;";
bufferHistDist = "bin;";
bufferArea = "";
bufferIntensity = "";
bufferSummaryPOM = "File;Mean Area;Std Area;Min Area;Max Area;Mean intensity;Std Intensity;Min Intensity;Max Intensity\n";
bufferAreaHistPOM = "bin;";
bufferIntensityHistPOM = "bin;";
nBins = 256;
increment = 256/nBins;

// initialize histogram buffer for BF and POM image
startBin = 0;
for(binCount=0; binCount<nBins; binCount++) {
	if(binCount == (nBins-1)) {
		bufferHistBF = bufferHistBF + startBin + "\n";
		bufferHistPOM = bufferHistPOM + startBin + "\n";
	} else {
		bufferHistBF = bufferHistBF + startBin + ";";
		bufferHistPOM = bufferHistPOM + startBin + ";";
	}
	startBin = startBin + increment;
}

// initialize histogram buffer for diameter measure
nBinsDist = maxDist - minDist + 1;
for(binCount=minDist; binCount<=maxDist; binCount++) {
	binVal = binCount;
	if(scale) binVal = binVal/pixPerMic;
	
	if(binCount == maxDist) {
		bufferHistDist = bufferHistDist + binVal + "\n";
	} else {
		bufferHistDist = bufferHistDist + binVal + ";";
	}
}

// initialize histogram buffer for area of segmented objects located inside POM fiber
for(j=init_area; j<end_area; j=j+bin_area) {
	val = j;
	if(scale) val = val * 1/pixPerMic * 1/pixPerMic;
	if((j+bin_area) >= end_area) bufferAreaHistPOM = bufferAreaHistPOM + val + "\n";
	else bufferAreaHistPOM = bufferAreaHistPOM + val + ";";
}

// initialize histogram buffer for intensity of segmented objects inside POM fiber
for(j=init_intensity; j<end_intensity; j=j+bin_intensity) {
	if((j+bin_intensity) >= end_intensity) bufferIntensityHistPOM = bufferIntensityHistPOM + j + "\n";
	else bufferIntensityHistPOM = bufferIntensityHistPOM + j + ";";
}

// iterate over all images of the input directory
for(cont=0; cont<dir.length; cont++) {
	file = dir[cont];
	
	if(matches(file, ".*BF.*") && !endsWith(file, "segmentation/")) { // read BF images and search for the corresponding POM image, with the exact same name
		print(path+file);
		aux = split(file,".");
		
		open(path+file);
		rename("image");
		if(scale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
		
		run("Duplicate...", "RGB");
		rename("RGB");
		run("Split Channels");

		// open corresponding POM image
		POMimage = replace(file, "BF", "PL");
		print(path+POMimage);
		open(path+POMimage);
		rename("POMimage");
		if(scale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");

		run("Duplicate...", "POM-RGB");
		rename("POM-RGB");
		run("Split Channels");

		selectWindow("image");
		run("8-bit");
		run("Duplicate...", "mask");
		rename("mask");

		selectWindow("POMimage");
		run("8-bit");

		selectWindow("mask");
		thres = segmentation(); // BF segmentation
		
		run("Analyze Particles...", "size=0-Infinity clear add");
		roiManager("Show All without labels");

		// find the biggest particle in the binary mask corresponding to the fiber
		biggestROI = 0; // area of the biggest ROI
		biggestInd = 0; // index of the biggest ROI
	
		for (i=0; i<nResults; i++) {
			area_i = getResult("Area", i);
			if(area_i > biggestROI) {
				biggestROI = area_i;
				biggestInd = i;
			}
		}
		
		roiManager("Select", biggestInd);
		run("Create Mask");
		rename("Fiber Mask");

		selectWindow("ROI Manager");
	    run("Close");

		if(saveFiberXYcoordinates()) {
			if(diagonalFiber) {
				cropXYCoordinates(); // crop the borders of the fibers accordingly
			}
			else {
				selectWindow("contour1");
				run("Select None");
				selectWindow("contour2");
				run("Select None");
			}
			close("outerParts");

			selectWindow("Fiber Mask");
			run("Create Selection");
	
			measures = measure("image", file); // measures for the 8-bit image
			bufferMeasuresBF = bufferMeasuresBF + measures;
			saveOverlay("image", outPath, file);
			bufferHistBF = bufferHistBF + histogram("image", file, nBins, true, 0, 256);
			selectWindow("Results"); 
		    run("Close");
		    
			measures = measure("POMimage", POMimage); // measures for the 8-bit POM image
			bufferMeasuresPOM = bufferMeasuresPOM + measures;
			saveOverlay("POMimage", outPath, POMimage);
			bufferHistPOM = bufferHistPOM + histogram("POMimage", POMimage, nBins, true, 0, 256);
	
			selectWindow("Results"); 
		    run("Close");
		    
			// POM segmentation
			selectWindow("POMimage");
			run("Clear Outside");
			run("Select None");
			thresPOM = segmentationPOM();
			run("Analyze Particles...", "clear add");
			selectWindow("POMimage");
			roiManager("Measure");
	
			numROIs = roiManager("count");
			
			bufferArea = bufferArea + POMimage + ";";
			bufferIntensity = bufferIntensity + POMimage + ";";
			bufferAreaHistPOM = bufferAreaHistPOM + POMimage + ";";
			bufferIntensityHistPOM = bufferIntensityHistPOM + POMimage + ";";
			
			areaPOM = newArray(numROIs);
			intensityPOM = newArray(numROIs);
			
			for(row=0; row<(numROIs-1); row++) {
				bufferArea = bufferArea + getResult("Area", row) + ";";
				bufferIntensity = bufferIntensity + getResult("Mean", row) + ";";

				areaPOM[row] = getResult("Area", row);
				intensityPOM[row] = getResult("Mean", row);
			}
			bufferArea = bufferArea + getResult("Area", row) + "\n";
			bufferIntensity = bufferIntensity + getResult("Mean", row) + "\n";
			areaPOM[row] = getResult("Area", row);
			intensityPOM[row] = getResult("Mean", row);
			
			countsArea = getHistogramPOMObjects(areaPOM,init_area,end_area,bin_area);
			countsIntensity = getHistogramPOMObjects(intensityPOM,init_intensity,end_intensity,bin_intensity);

			for(j=0; j<countsArea.length; j++) {
				if(j == (countsArea.length-1)) bufferAreaHistPOM = bufferAreaHistPOM + countsArea[j] + "\n";
				else bufferAreaHistPOM = bufferAreaHistPOM + countsArea[j] + ";";
			}

			for(j=0; j<countsIntensity.length; j++) {
				if(j == (countsIntensity.length-1)) bufferIntensityHistPOM = bufferIntensityHistPOM + countsIntensity[j] + "\n";
				else bufferIntensityHistPOM = bufferIntensityHistPOM + countsIntensity[j] + ";";
			}
		
			selectWindow("POMimage");
			run("From ROI Manager");
			roiManager("Show All without labels");
			saveOverlay("POMimage", outPath, "segmented_"+POMimage);
			roiManager("Save", outPath+POMimage+"RoiSet.zip");
			close("POMmask");

			if(numROIs > 1) {
				run("Summarize");
				avgArea = getResult("Area", numROIs);
				avgIntensity = getResult("Mean", numROIs);
				bufferSummaryPOM = bufferSummaryPOM + POMimage + ";" + avgArea + ";" + getResult("Area", (numROIs+1)) + ";" + getResult("Area", (numROIs+2)) + ";" + getResult("Area", (numROIs+3)) + ";" + avgIntensity + ";" + getResult("Mean", (numROIs+1)) + ";" + getResult("Mean", (numROIs+2)) + ";" + getResult("Mean", (numROIs+3)) + "\n";
			} else {
				avgArea = getResult("Area", 0);
				avgIntensity = getResult("Mean", 0);
				bufferSummaryPOM = bufferSummaryPOM + POMimage + ";" + avgArea + ";0;" + avgArea + ";" + avgArea + ";" + avgIntensity + ";0;" + avgIntensity + ";" + avgIntensity + "\n";
			}

			selectWindow("ROI Manager"); 
		    run("Close");
		    selectWindow("Results"); 
		    run("Close");
	
			// measure diameter
			selectWindow("Fiber Mask");
			run("Select None");
			
			diameter = measureDiameter("contour1", "contour2");
			diameter = diameter + ";" + measureDiameter("contour2", "contour1");
			
			diamVec = split(diameter, ";");
			avgDiam = (parseFloat(diamVec[0]) + parseFloat(diamVec[4])) / 2;
			avgStd = (parseFloat(diamVec[1]) + parseFloat(diamVec[5])) / 2;
			
			// calculate average diameter (from both directions), std and threshold method to the buffer
			bufferMeasuresBF = bufferMeasuresBF + ";" + diameter + ";" + avgDiam + ";" + avgStd + ";" + thres + "\n";
			bufferMeasuresPOM = bufferMeasuresPOM + ";" + diameter + ";" + avgDiam + ";" + avgStd + ";" + thresPOM + "\n";
	
			// create the segmentation mask again to run "measures" on the RGB individual channels
		    selectWindow("image");
			run("Create Mask");
			run("Analyze Particles...", "size=0-Infinity clear add");
			roiManager("Show All without labels");
			roiManager("Set Line Width", 2);
	
		    ids = newArray("RGB (red)", "RGB (green)", "RGB (blue)");
			names = newArray(file+"_red", file+"_green", file+"_blue");
			namesPOM = newArray(POMimage+"_red", POMimage+"_green", POMimage+"_blue");
	
			for(i=0; i<ids.length; i++) {
				measures = measure(ids[i], names[i]);
				bufferMeasuresBF = bufferMeasuresBF + measures + ";" + diameter + ";" + avgDiam + ";" + avgStd + ";" + thres + "\n";
				selectWindow("Results"); 
		    	run("Close");
				
				measures = measure("POM-"+ids[i], namesPOM[i]);
				bufferMeasuresPOM = bufferMeasuresPOM + measures + ";" + diameter + ";" + avgDiam + ";" + avgStd + ";" + thresPOM + "\n";
				
				bufferHistBF = bufferHistBF + histogram(ids[i], names[i], nBins, true, 0, 256);
				bufferHistPOM = bufferHistPOM + histogram("POM-"+ids[i], namesPOM[i], nBins, true, 0, 256);
	
				selectWindow("Results");
	    		run("Close");
			}
			
			selectWindow("ROI Manager"); 
		    run("Close");
	
		    // histogram of the distance transform
		    bin = maxDist-minDist;
			bufferHistDist = bufferHistDist + histogram("EDT-contour1", file, nBinsDist, false, minDist, maxDist);
			bufferHistDist = bufferHistDist + histogram("EDT-contour2", file, nBinsDist, false, minDist, maxDist);
			
			run("Close All");
		}
		run("Close All");
	}
}

// create files and save buffers
summaryPathBF = outPath+"summaryBF.csv";
histPathBF = outPath+"histBF.csv";
summaryPathPOM = outPath+"summaryPOM.csv";
histPathPOM = outPath+"histPOM.csv";
histPathDist = outPath+"histDiameter.csv";
summaryAreaFiberPOM_perObj = outPath+"areaPerFiberPOM_perObj.csv";
summaryIntensityFiberPOM_perObj = outPath+"intensityPerFiberPOM_perObj.csv";
summaryAreaFiberPOM = outPath+"areaIntensityPerFiberPOM.csv";
histAreaFiberPOM = outPath+"histogramAreaPerFiberPOM.csv";
histIntensityFiberPOM = outPath+"histogramIntensityPerFiberPOM.csv";

if(File.exists(summaryPathBF))
	File.delete(summaryPathBF);
if(File.exists(histPathBF))
	File.delete(histPathBF);
if(File.exists(summaryPathPOM))
	File.delete(summaryPathPOM);
if(File.exists(histPathPOM))
	File.delete(histPathPOM);
if(File.exists(histPathDist))
	File.delete(histPathDist);
if(File.exists(summaryAreaFiberPOM_perObj))
	File.delete(summaryAreaFiberPOM_perObj);
if(File.exists(summaryIntensityFiberPOM_perObj))
	File.delete(summaryIntensityFiberPOM_perObj);
if(File.exists(summaryAreaFiberPOM))
	File.delete(summaryAreaFiberPOM);
if(File.exists(histAreaFiberPOM))
	File.delete(histAreaFiberPOM);
if(File.exists(histIntensityFiberPOM))
	File.delete(histIntensityFiberPOM);

summaryFileBF = File.open(summaryPathBF);
print(summaryFileBF, bufferMeasuresBF);
File.close(summaryFileBF);

histFileBF = File.open(histPathBF);
print(histFileBF, bufferHistBF);
File.close(histFileBF);

summaryFilePOM = File.open(summaryPathPOM);
print(summaryFilePOM, bufferMeasuresPOM);
File.close(summaryFilePOM);

histFilePOM = File.open(histPathPOM);
print(histFilePOM, bufferHistPOM);
File.close(histFilePOM);

histFileDist = File.open(histPathDist);
print(histFileDist, bufferHistDist);
File.close(histFileDist);

areaFiberFilePOM_perObj = File.open(summaryAreaFiberPOM_perObj);
print(areaFiberFilePOM_perObj, bufferArea);
File.close(areaFiberFilePOM_perObj);

intesityFiberFilePOM_perObj = File.open(summaryIntensityFiberPOM_perObj);
print(intesityFiberFilePOM_perObj, bufferIntensity);
File.close(intesityFiberFilePOM_perObj);

areaFiberFilePOM = File.open(summaryAreaFiberPOM);
print(areaFiberFilePOM, bufferSummaryPOM);
File.close(areaFiberFilePOM);

histAreaFiberPOMFile = File.open(histAreaFiberPOM);
print(histAreaFiberPOMFile, bufferAreaHistPOM);
File.close(histAreaFiberPOMFile);

histIntensityFiberPOMFile = File.open(histIntensityFiberPOM);
print(histIntensityFiberPOMFile, bufferIntensityHistPOM);
File.close(histIntensityFiberPOMFile);

if(File.exists(path+"contour1.txt"))
	File.delete(path+"contour1.txt");

if(File.exists(path+"contour2.txt"))
	File.delete(path+"contour2.txt");

print("finished");

/************************** functions **************************/

function getHistogramPOMObjects(values, init, end, bin) {
	sizeHist = Math.ceil((end-init)/bin);
	counts = newArray(sizeHist);

	ind = 0;
	for(j=init; j<end; j=j+bin) {
		counts[ind] = 0;
		ind++;
	}

	for(i=0; i<values.length; i++) {
		ind = 0;
		for(j=init; j<=end; j=j+bin) {
			val = j;
			if(scale) val = val * 1/pixPerMic * 1/pixPerMic;
			if(values[i] >= j && values[i] < (j+bin)) {
				counts[ind] = counts[ind] + 1;
			}
			ind++;
		}
	}
	return counts;
}

function saveFiberXYcoordinates() {
	saved = false;
	
	run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
	selectWindow("Fiber Mask");	
	run("Duplicate...", "outerParts");
	rename("outerParts");
	run("Invert");
	
	run("Analyze Particles...", "size="+1000+"-Infinity clear add"); // filter small particles 

	if(roiManager("count") == 2) {
		exportXYCoordinates("outerParts", 0, "contour1");
		exportXYCoordinates("outerParts", 1, "contour2");
		saved = true;
	}

	selectWindow("ROI Manager"); 
    run("Close");
		    
	return saved;
}

function exportXYCoordinates(imageName, roiID, exportName) {
	selectWindow(imageName);
	roiManager("Select", roiID);
	run("Create Mask");
	setOption("BlackBackground", true);
	run("Dilate");
	imageCalculator("AND create", "Fiber Mask","Mask");
	selectWindow("Result of Fiber Mask");
	rename(exportName);
	close("Mask");

	selectWindow(exportName);
	run("Create Selection");
	run("Save XY Coordinates...", "save=["+path+exportName+".txt]");
}

function cropXYCoordinates() {
	
	contour1 = File.openAsString(path+"contour1.txt");
	contour2 = File.openAsString(path+"contour2.txt");

	lines1 = split(contour1,"\n");
	lines2 = split(contour2,"\n");

	lines1 = Array.slice(lines1,1,lines1.length);
	lines2 = Array.slice(lines2,1,lines2.length);

	first_line1 = split(lines1[0],"\t");
	first_line2 = split(lines2[0],"\t");

	last_line1 = split(lines1[(lines1.length-1)],"\t");
	last_line2 = split(lines2[(lines2.length-1)],"\t");

	// get coordinates of first pixels from contour1 and contour2
	x1_i = parseInt(first_line1[0]);
	y1_i = parseInt(first_line1[1]);

	x2_i = parseInt(first_line2[0]);
	y2_i = parseInt(first_line2[1]);

	x1_f = parseInt(last_line1[0]);
	y1_f = parseInt(last_line1[1]);

	x2_f = parseInt(last_line2[0]);
	y2_f = parseInt(last_line2[1]);

	/*************************************************************************/
	// create variables to store which indices should be excluded - up to bottom
	ind1_i = 0;
	ind2_i = 0;

	if(y1_i < y2_i) { // then contour 1 is the upper boundary of the fiber
		if(x1_i != x2_i) {
			if(x1_i > y2_i) {
				ind2_i = nPixelsCorner;
			} else {
				ind1_i = nPixelsCorner;
			}
		} else {
		}
	} else if(y1_i == y2_i) {
		if(x1_i > x2_i) { // then contour 1 is the upper boundary of the fiber
			ind2_i = nPixelsSide;
		} else {
			ind1_i = nPixelsSide;
		}
	} else  { // then contour 2 is the upper boundary of the fiber
		if(x2_i != x1_i) {
			if(x2_i > y1_i) {
				ind1_i = nPixelsCorner;
			} else {
				ind2_i = nPixelsCorner;
			}
		} else {
			ind2_i = nPixelsSide;
		}
	}

	/*************************************************************************/
	// create variables to store which indices should be excluded - bottom up
	ind1_f = lines1.length-1;
	ind2_f = lines2.length-1;

	// get widht and height of contour image
	selectWindow("contour1");
	width = getWidth();
	height = getHeight();

	if(y1_f < y2_f) { // then contour 1 is the upper boundary of the fiber
		if((width-x1_f) != (width-x2_f)) {
			if( (height - y1_f) > (width - x2_f) ) {
				ind2_f = lines2.length-1-nPixelsCorner;
			} else {
				ind1_f = lines1.length-1-nPixelsCorner;
			}
		} else {
			ind2_f = lines2.length-1-nPixelsSide;
		}
	} else if(y1_f == y2_f) {
		if(x1_f > x2_f) { // then contour 1 is the upper boundary of the fiber
			ind1_f = lines1.length-1-nPixelsSide;
		} else {
			ind2_f = lines2.length-1-nPixelsSide;
		}
	} else { // then contour 2 is the upper boundary of the fiber
		if((width-x1_f) != (width-x2_f)) {
			if((height - y2_f) > (width - x1_f)) {
				ind1_f = lines1.length-1-nPixelsCorner;
			} else {
				ind2_f = lines2.length-1-nPixelsCorner;
			}
		} else {
			ind1_f = lines1.length-1-nPixelsSide;
		}
	}

	selectWindow("contour1");
	run("Select None");
	for(a=(lines1.length-1); a>=ind1_f; a--) {
		line = split(lines1[a],"\t");
		setPixel(line[0], line[1], 0);
	}
	for(a=ind1_i; a>=0; a--) {
		line = split(lines1[a],"\t");
		setPixel(line[0], line[1], 0);
	}

	selectWindow("contour2");
	run("Select None");
	for(a=(lines2.length-1); a>=ind2_f; a--) {
		line = split(lines2[a],"\t");
		setPixel(line[0], line[1], 0);
	}
	for(a=ind2_i; a>=0; a--) {
		line = split(lines2[a],"\t");
		setPixel(line[0], line[1], 0);
	}

	close("Fiber Mask");
	imageCalculator("Add create", "contour1","contour2");
	selectWindow("Result of contour1");
	rename("Fiber Mask");
	
	li_1 = split(lines1[ind1_i],"\t");
	lf_1 = split(lines1[ind1_f],"\t");
	li_1_xi = parseInt(li_1[0]);
	li_1_yi = parseInt(li_1[1]);
	lf_1_xf = parseInt(lf_1[0]);
	lf_1_yf = parseInt(lf_1[1]);
	
	li_2 = split(lines2[ind2_i],"\t");
	lf_2 = split(lines2[ind2_f],"\t");
	li_2_xi = parseInt(li_2[0]);
	li_2_yi = parseInt(li_2[1]);
	lf_2_xf = parseInt(lf_2[0]);
	lf_2_yf = parseInt(lf_2[1]);

	makeLine(li_1_xi, li_1_yi, li_2_xi, li_2_yi);
	run("Create Mask");
	imageCalculator("Add", "Fiber Mask","Mask");
	close("Mask");
	makeLine(lf_1_xf, lf_1_yf, lf_2_xf, lf_2_yf);
	run("Create Mask");
	imageCalculator("Add", "Fiber Mask","Mask");
	close("Mask");

	selectWindow("Fiber Mask");
	run("Fill Holes");
}

function saveOverlay(imageName, outPath, fileName) {
	run("Flatten");
	saveAs("png", outPath + fileName + ".png");
	close(fileName + ".png");
}

function measure(imageName, fileName) {
	selectWindow(imageName);
	run("Restore Selection");
		
	run("Set Measurements...", "area mean standard min perimeter redirect=None decimal=3");
	if(!scale) run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
	run("Measure");

	area = getResultString("Area", 0);
	mean = getResultString("Mean", 0);
	std = getResultString("StdDev", 0);
	max = getResultString("Max", 0);
	min = getResultString("Min", 0);
	perim = getResultString("Perim.", 0);
	
	return fileName + ";" + area + ";" + mean + ";" + std + ";" + max + ";" + min + ";" + perim;
}

function measureDiameter(img1, img2) { // diamenter from img1 to img2
	diam = "";
	selectWindow(img1);
	run("Invert");
	run("Exact Euclidean Distance Transform (3D)");
	rename("EDT-"+img1);
	selectWindow(img1);
	run("Invert");
	
	selectWindow(img2);
	run("Set Measurements...", "area mean standard min perimeter redirect=None decimal=3");
	run("Analyze Particles...", "clear add");
	selectWindow("EDT-"+img1);
	run("Conversions...", "weighted");
	setOption("ScaleConversions", false);
	run("16-bit");
	run("From ROI Manager");
	roiManager("Measure");
	meanContour = getResult("Mean", 0);
	stdContour = getResult("StdDev", 0);
	minContour = getResult("Min", 0);
	maxContour = getResult("Max", 0);
	
	if(scale) diam = diam + (meanContour/pixPerMic) + ";" + (stdContour/pixPerMic) + ";" + (minContour/pixPerMic) + ";" + (maxContour/pixPerMic);
	else diam = diam + meanContour + ";" + stdContour + ";" + minContour + ";" + maxContour;

	selectWindow("Results"); 
	run("Close");
	selectWindow("ROI Manager"); 
	run("Close");
	
	return diam;
}

function histogram(imageName, fileName, n_bins, restore, minValue, maxValue) {
	selectWindow(imageName);
	if(restore) run("Restore Selection");

	if(maxValue == 256) getHistogram(values, counts, n_bins);
	else getHistogram(values, counts, n_bins, minValue, maxValue);
	line = fileName + ";";
	for(i=0; i<n_bins; i++) {
		if(i == (n_bins-1)) {
			line = line + counts[i] + "\n";
		}
		else line = line + counts[i] + ";";
	}
	
	return line;
}

function segmentation() {
	
	if(useThresholdBF_manual) {
		setThreshold(0, thresholdBF_manual);	
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Fill Holes");
		return thresholdBF_manual;
	} else {
		setAutoThreshold(thresholdBF_automatic+" dark");
		setOption("BlackBackground", true);
		getThreshold(lower, upper);
		run("Convert to Mask");
		run("Invert");
		run("Fill Holes");
		return lower;
	}
}

function segmentationPOM() {
	// save the area per fiber components
	selectWindow("POMimage");

	run("Duplicate...", "POMmask");
	rename("POMmask");
	run("Restore Selection");		

	if(useThresholdPOM_manual)
		setThreshold(thresholdPOM_manual, 255);	
	else
		setAutoThreshold(thresholdPOM_automatic+" dark");
		
	setOption("BlackBackground", true);
	getThreshold(lower, upper);
	run("Convert to Mask");
		
    return lower;
}