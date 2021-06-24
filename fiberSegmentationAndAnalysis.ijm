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
scale = true;

// segmentation of the BF image
thresholdMethod1 = "Yen"; // automatic threshold
thredholdManual1 = 103; // manual threshold
useThresholdManual1 = false; // true or false - whether to use the manual threshold of not

// segmentation of the PL image
thresholdMethod2 = "Yen"; // can be different from thresholdMethod1
thredholdManual2 = 20;
useThresholdManual2 = false; // true or false

diagonalFiber = true; // true or false - whether set of fibers are diagonal or not
nPixelsCorner = 10; // how many pixels should be eliminated from the fiber if endings are in the corners 
nPixelsSide = 20; // how many pixels should be eliminated from the fiber if endings are not in the corners

// paramenters for the histogram of area and average intensity of the particles segmented in the PL - ending points and bin size
init_area = 0;
end_area = 20000;
bin_area = 50;

init_intensity = 20;
end_intensity = 255;
bin_intensity = 5;
/************************** parameters **************************/

// read input directory
path = getDirectory("Choose a Directory");
dir = getFileList(path);

// create output directory
outPath = path + "segmentation/";
File.makeDirectory(outPath);

// initialize buffers to save measures and histograms
bufferMeasuresBF = "File;Area;Mean;StdDev;Max;Min;Perim;avgDiam1;stdDiam1;minDiam1;maxDiam1;avgDiam2;stdDiam2;minDiam2;maxDiam2;avgDiam;stdDiam;thresholdValue\n";
bufferMeasuresPL = "File;Area;Mean;StdDev;Max;Min;Perim;avgDiam1;stdDiam1;minDiam1;maxDiam1;avgDiam2;stdDiam2;minDiam2;maxDiam2;avgDiam;stdDiam;thresholdValue\n";
bufferHistBF = "bin;";
bufferHistPL = "bin;";
bufferHistDist = "bin;";
bufferArea = "";
bufferIntensity = "";
bufferSummaryPL = "File;Mean Area;Std Area;Min Area;Max Area;Mean intensity;Std Intensity;Min Intensity;Max Intensity\n";
bufferAreaHistPL = "bin;";
bufferIntensityHistPL = "bin;";
nBins = 256;
increment = 256/nBins;

// initialize histogram files for BF and PL image
startBin = 0;
for(binCount=0; binCount<nBins; binCount++) {
	if(binCount == (nBins-1)) {
		bufferHistBF = bufferHistBF + startBin + "\n";
		bufferHistPL = bufferHistPL + startBin + "\n";
	} else {
		bufferHistBF = bufferHistBF + startBin + ";";
		bufferHistPL = bufferHistPL + startBin + ";";
	}
	startBin = startBin + increment;
}

// initialize histogram file for diameter measure
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

// initialize histogram file for area of segmented objects inside PL fiber
for(j=init_area; j<end_area; j=j+bin_area) {
	val = j;
	if(scale) val = val * 1/pixPerMic * 1/pixPerMic;
	if((j+bin_area) >= end_area) bufferAreaHistPL = bufferAreaHistPL + val + "\n";
	else bufferAreaHistPL = bufferAreaHistPL + val + ";";
}

// initialize histogram file for intensity of segmented objects inside PL fiber
for(j=init_intensity; j<end_intensity; j=j+bin_intensity) {
	if((j+bin_intensity) >= end_intensity) bufferIntensityHistPL = bufferIntensityHistPL + j + "\n";
	else bufferIntensityHistPL = bufferIntensityHistPL + j + ";";
}

// iterate over all images of the input directory
for(cont=0; cont<dir.length; cont++) {
	file = dir[cont];
	
	if(matches(file, ".*BF.*") && !endsWith(file, "segmentation/")) { // read BF images and search for the corresponding PL image, with the exact same name
		print(path+file);
		aux = split(file,".");
		
		open(path+file);
		rename("image");
		if(scale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");
		
		run("Duplicate...", "RGB");
		rename("RGB");
		run("Split Channels");

		// open corresponding PL image
		PLimage = replace(file, "BF", "PL");
		print(path+PLimage);
		open(path+PLimage);
		rename("PLimage");
		if(scale) run("Set Scale...", "distance=" + pixPerMic + " known=1 unit=micron");

		run("Duplicate...", "PL-RGB");
		rename("PL-RGB");
		run("Split Channels");

		selectWindow("image");
		run("8-bit");
		run("Duplicate...", "mask");
		rename("mask");

		selectWindow("PLimage");
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
			//run("Interpolate", "interval=1");
	
			measures = measure("image", file); // measures for the 8-bit image
			bufferMeasuresBF = bufferMeasuresBF + measures;
			saveOverlay("image", outPath, file);
			bufferHistBF = bufferHistBF + histogram("image", file, nBins, true, 0, 256);
			selectWindow("Results"); 
		    run("Close");
		    
			measures = measure("PLimage", PLimage); // measures for the 8-bit PL image
			bufferMeasuresPL = bufferMeasuresPL + measures;
			saveOverlay("PLimage", outPath, PLimage);
			bufferHistPL = bufferHistPL + histogram("PLimage", PLimage, nBins, true, 0, 256);
	
			selectWindow("Results"); 
		    run("Close");
		    
			// PL segmentation
			selectWindow("PLimage");
			run("Clear Outside");
			run("Select None");
			thresPL = segmentationPL();
			run("Analyze Particles...", "clear add");
			selectWindow("PLimage");
			roiManager("Measure");
	
			numROIs = roiManager("count");
			
			bufferArea = bufferArea + PLimage + ";";
			bufferIntensity = bufferIntensity + PLimage + ";";
			bufferAreaHistPL = bufferAreaHistPL + PLimage + ";";
			bufferIntensityHistPL = bufferIntensityHistPL + PLimage + ";";
			
			areaPL = newArray(numROIs);
			intensityPL = newArray(numROIs);
			
			for(row=0; row<(numROIs-1); row++) {
				bufferArea = bufferArea + getResult("Area", row) + ";";
				bufferIntensity = bufferIntensity + getResult("Mean", row) + ";";

				areaPL[row] = getResult("Area", row);
				intensityPL[row] = getResult("Mean", row);
			}
			bufferArea = bufferArea + getResult("Area", row) + "\n";
			bufferIntensity = bufferIntensity + getResult("Mean", row) + "\n";
			areaPL[row] = getResult("Area", row);
			intensityPL[row] = getResult("Mean", row);
			
			countsArea = getHistogramPLObjects(areaPL,init_area,end_area,bin_area);
			countsIntensity = getHistogramPLObjects(intensityPL,init_intensity,end_intensity,bin_intensity);

			for(j=0; j<countsArea.length; j++) {
				if(j == (countsArea.length-1)) bufferAreaHistPL = bufferAreaHistPL + countsArea[j] + "\n";
				else bufferAreaHistPL = bufferAreaHistPL + countsArea[j] + ";";
			}

			for(j=0; j<countsIntensity.length; j++) {
				if(j == (countsIntensity.length-1)) bufferIntensityHistPL = bufferIntensityHistPL + countsIntensity[j] + "\n";
				else bufferIntensityHistPL = bufferIntensityHistPL + countsIntensity[j] + ";";
			}
		
			selectWindow("PLimage");
			run("From ROI Manager");
			roiManager("Show All without labels");
			saveOverlay("PLimage", outPath, "segmented_"+PLimage);
			roiManager("Save", outPath+PLimage+"RoiSet.zip");
			close("PLmask");

			if(numROIs > 1) {
				run("Summarize");
				avgArea = getResult("Area", numROIs);
				avgIntensity = getResult("Mean", numROIs);
				bufferSummaryPL = bufferSummaryPL + PLimage + ";" + avgArea + ";" + getResult("Area", (numROIs+1)) + ";" + getResult("Area", (numROIs+2)) + ";" + getResult("Area", (numROIs+3)) + ";" + avgIntensity + ";" + getResult("Mean", (numROIs+1)) + ";" + getResult("Mean", (numROIs+2)) + ";" + getResult("Mean", (numROIs+3)) + "\n";
			} else {
				avgArea = getResult("Area", 0);
				avgIntensity = getResult("Mean", 0);
				bufferSummaryPL = bufferSummaryPL + PLimage + ";" + avgArea + ";0;" + avgArea + ";" + avgArea + ";" + avgIntensity + ";0;" + avgIntensity + ";" + avgIntensity + "\n";
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
			bufferMeasuresPL = bufferMeasuresPL + ";" + diameter + ";" + avgDiam + ";" + avgStd + ";" + thresPL + "\n";
	
			// create the segmentation mask again to run "measures" on the RGB channels, individually
		    selectWindow("image");
			run("Create Mask");
			run("Analyze Particles...", "size=0-Infinity clear add");
			roiManager("Show All without labels");
			roiManager("Set Line Width", 2);
	
		    ids = newArray("RGB (red)", "RGB (green)", "RGB (blue)");
			names = newArray(file+"_red", file+"_green", file+"_blue");
			namesPL = newArray(PLimage+"_red", PLimage+"_green", PLimage+"_blue");
	
			for(i=0; i<ids.length; i++) {
				measures = measure(ids[i], names[i]);
				bufferMeasuresBF = bufferMeasuresBF + measures + ";" + diameter + ";" + avgDiam + ";" + avgStd + ";" + thres + "\n";
				selectWindow("Results"); 
		    	run("Close");
				
				measures = measure("PL-"+ids[i], namesPL[i]);
				bufferMeasuresPL = bufferMeasuresPL + measures + ";" + diameter + ";" + avgDiam + ";" + avgStd + ";" + thresPL + "\n";
				
				bufferHistBF = bufferHistBF + histogram(ids[i], names[i], nBins, true, 0, 256);
				bufferHistPL = bufferHistPL + histogram("PL-"+ids[i], namesPL[i], nBins, true, 0, 256);
	
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
summaryPathPL = outPath+"summaryPL.csv";
histPathPL = outPath+"histPL.csv";
histPathDist = outPath+"histDiameter.csv";
summaryAreaFiberPL_perObj = outPath+"areaPerFiberPL_perObj.csv";
summaryIntensityFiberPL_perObj = outPath+"intensityPerFiberPL_perObj.csv";
summaryAreaFiberPL = outPath+"areaIntensityPerFiberPL.csv";
histAreaFiberPL = outPath+"histogramAreaPerFiberPL.csv";
histIntensityFiberPL = outPath+"histogramIntensityPerFiberPL.csv";

if(File.exists(summaryPathBF))
	File.delete(summaryPathBF);
if(File.exists(histPathBF))
	File.delete(histPathBF);
if(File.exists(summaryPathPL))
	File.delete(summaryPathPL);
if(File.exists(histPathPL))
	File.delete(histPathPL);
if(File.exists(histPathDist))
	File.delete(histPathDist);
if(File.exists(summaryAreaFiberPL_perObj))
	File.delete(summaryAreaFiberPL_perObj);
if(File.exists(summaryIntensityFiberPL_perObj))
	File.delete(summaryIntensityFiberPL_perObj);
if(File.exists(summaryAreaFiberPL))
	File.delete(summaryAreaFiberPL);
if(File.exists(histAreaFiberPL))
	File.delete(histAreaFiberPL);
if(File.exists(histIntensityFiberPL))
	File.delete(histIntensityFiberPL);

summaryFileBF = File.open(summaryPathBF);
print(summaryFileBF, bufferMeasuresBF);
File.close(summaryFileBF);

histFileBF = File.open(histPathBF);
print(histFileBF, bufferHistBF);
File.close(histFileBF);

summaryFilePL = File.open(summaryPathPL);
print(summaryFilePL, bufferMeasuresPL);
File.close(summaryFilePL);

histFilePL = File.open(histPathPL);
print(histFilePL, bufferHistPL);
File.close(histFilePL);

histFileDist = File.open(histPathDist);
print(histFileDist, bufferHistDist);
File.close(histFileDist);

areaFiberFilePL_perObj = File.open(summaryAreaFiberPL_perObj);
print(areaFiberFilePL_perObj, bufferArea);
File.close(areaFiberFilePL_perObj);

intesityFiberFilePL_perObj = File.open(summaryIntensityFiberPL_perObj);
print(intesityFiberFilePL_perObj, bufferIntensity);
File.close(intesityFiberFilePL_perObj);

areaFiberFilePL = File.open(summaryAreaFiberPL);
print(areaFiberFilePL, bufferSummaryPL);
File.close(areaFiberFilePL);

histAreaFiberPLFile = File.open(histAreaFiberPL);
print(histAreaFiberPLFile, bufferAreaHistPL);
File.close(histAreaFiberPLFile);

histIntensityFiberPLFile = File.open(histIntensityFiberPL);
print(histIntensityFiberPLFile, bufferIntensityHistPL);
File.close(histIntensityFiberPLFile);

if(File.exists(path+"contour1.txt"))
	File.delete(path+"contour1.txt");

if(File.exists(path+"contour2.txt"))
	File.delete(path+"contour2.txt");

print("finished");

/************************** functions **************************/

function getHistogramPLObjects(values, init, end, bin) {
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
	//run("Clear Outside");

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
	
	if(useThresholdManual1) {
		setThreshold(0, thredholdManual1);	
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("Fill Holes");
		return thredholdManual1;
	} else {
		setAutoThreshold(thresholdMethod1+" dark");
		setOption("BlackBackground", true);
		getThreshold(lower, upper);
		run("Convert to Mask");
		run("Invert");
		run("Fill Holes");
		return lower;
	}
}

function segmentationPL() {
	// save the area per fiber components
	selectWindow("PLimage");

	run("Duplicate...", "PLmask");
	rename("PLmask");
	run("Restore Selection");		

	if(useThresholdManual2)
		setThreshold(thredholdManual2, 255);	
	else
		setAutoThreshold(thresholdMethod2+" dark");
		
	setOption("BlackBackground", true);
	getThreshold(lower, upper);
	run("Convert to Mask");
		
    return lower;
}