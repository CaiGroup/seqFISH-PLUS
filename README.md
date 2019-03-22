# seqFISH-PLUS
Give a one line description

## Getting Started
* Download all the contents of the folder and add it to your Matlab path.

### Prerequisites
* MATLAB

### Dependencies
* radialcenter.m by Dr. Parthasarathy, which can be downloaded [here](https://media.nature.com/original/nature-assets/nmeth/journal/v9/n7/extref/nmeth.2071-S2.zip).

## Running the Code
* Overview of the steps for the seqFISH+ pipeline *
### Prerequirements
* points spread function of a bead
* thresholding values after Preprocessing Step

### Preprocessing
* BackCorr_After_Analyze_ALL_channels.m
    * a. grab_ims_NO_BEADS.m - grabs all the hybridization images in the 3 channels
    * b. grab_ims_NO_BEADS.m - grab the dapi channel in channel 4 
    * c. imregcorr.m - calculate the dapi transformation to the first reference dapi image
    * d. deconvlucy.m - deconvolve the image with the point spread function
    * e. grag_ims_Blanks.m - grab the blank images
    * f. back_correct_multi_pos.m - get the background corrections
    * g. findDotsBarcodeV2.m - find the dots using a laplacian of gaussians and a regional maxima
    * h. radialcenter.m - super resolve dots with radial center
    * i. transformPointsForward - apply transformation to dots
### Postprocessing
* Callbarcodes_ALL_channels_cell_by_cell.m
    * a. gene_calling_sandbox_3_R.m - colocolaize all the points and recieve all the found barcodes
    * b. handle_ambiguous_2_segged.m - match points to genes and remove ambiguous barcodes
    * c. PointLocations_NEW_no_intensity.m – create the point locations and lists for each gene
    * d. Numseeds_NEW.m – calculate the number of seeds for each point
    * e. filterseedsLinus_marked_up.m – filter the points that have seeds 3 or 4 (not in function)

## License
Free for non-commercial and academic research. The software is bound by the licensing rules of California Institute of Technology (Caltech)

## Acknowledgements
Mike Lawson - Writing and implementing code for the project
Sheel Shah - Developing the algorithm to find the barcodes, finding dots, and much more

## Contact
Linus 
Nico Pierson



