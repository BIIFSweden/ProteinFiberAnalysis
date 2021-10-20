# An image-based analysis of recombinant silk fiber structures

The pipeline provided in this project performs segmentation and quantification of recombinant silk fibers using Fiji and it was implemented in ImageJ Macro Language. The input of the pipeline is a pair of images depicting a sample fiber from two microscopy techniques, brightfield (BF) and crossed polarizers (POM). The brightfield counterpart will be used for segmentation using Yen's threshold method, which yields a binary mask separating the fiber from the background. Next, a second segmentation is applied to the POM image restricted only to the fiber region. This is performed in order to characterize the different regions in the fiber with varying intensity patterns. The whole process is illustrated below

![Pipeline](img/scheme_image.png "Title Text")

Different measures are exctracted from the segmented fiber and fiber regions: area, average intensity, average diameter and the corresponding standard deviation, minimum diameter, maximum diameter and the threshold value used. The fiber diameter was estimated using the Euclidean distance transform (EDT) computed from the upper boundary of the fiber towards the lower boundary and vice-versa. Both values were averaged in order to obtain the final diameter. 

### 1.	Input directory

Figure 1 shows an example of an input directory containing the images to be analyzed. Images should be grouped by pairs of BF (brightfield) and POM (polarized microscopy) images.

### 2.	Software requirements

The software listed below should be installed before running the Fiji script. 

* [Fiji](https://fiji.sc): follow the instructions in the link to download Fiji.

### 3.	Running the pipeline

To run the pipeline, open Fiji and go to Plugins – Macros – Edit... and browse the *fiberSegmentationAndAnalysis.ijm* file to load the script. Then, the interface showed in Fig. 2 will appear. Between lines 10 and 36, different parameters can be modified, as described next:

#### 3.1	Parameters related to size and scale and additional settings

* minDist: the minimum distance between any two pixels of the fiber (minimum diameter)
* maxDist: the minimum distance between any two pixels of the fiber (maximum diameter)
* pixPerMic: number of pixels per micron 
* scale: true of false, whether to use or not the given scale.
* diagonalFiber: true or false - whether set of fibers are diagonal or not
* nPixelsCorner: how many pixels should be eliminated from the fiber if endings are in the corners (check Fig. 3.)
* nPixelsSide: how many pixels should be eliminated from the fiber if endings are not in the corners (check Fig. 3.)

#### 3.2	Parameters related to thresholding

* thresholdMethod1: automatic threshold method for the segmentation of the BF image
* thresholdManual1: user-defined threshold value for the segmentation of the BF image
* useThresholdManual1: true or false - whether or not to use the manual threshold for the BF image, if false then the corresponding automatic threshold is used.
* thresholdMethod2: automatic threshold method for the segmentation of the PL image
* thresholdManual2: user-defined threshold value for the segmentation of the PL image
* useThresholdManual2: true or false - whether or not to use the manual threshold for the PL image, if false then the corresponding automatic threshold is used.

#### 3.3	Parameters related to histogram generation

* init_area:
* end_area:
* bin_area:
* init_intensity:
* end_intensity:
* bin_intensity:


### 4.	Ouput files
