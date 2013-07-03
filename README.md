Tang_et_al_LINCS_cell_scoring
=============================

Code for scoring cells in images from LTR-NucView assay published in "Tang, Y. et al. Differential Determinants of Cancer Cell Insensitivity to Antimitotic Drugs Discriminated by a One-Step Cell Imaging Assay. J Biomol Screen (2013). doi:10.1177/1087057113493804"

There is only one input required to run this code, which is the path to the folder where the input images are located. Within this folder, all the input images are expected to be stored in the same folder named "input_images". Three input channel are expected, and they are DAPI/Hoechst, Lysotracker red (LTR), and NucView 488 respectively. The file names of the three channels need to be listed in three separate text files named "blue_files.txt", "red_files.txt" and "green_files.txt" respectively.

To configure the parameters for this code, the 1st section of the code needs to be edited.
nsite: # of images acquired per well
max_nuc_size: diameter of largest nucleus expected in pixels
mean_nuc_size: average diameter of nuclei in pixels
initial_nuc_int_cutoff: initial threshold adjustment factor
min_nuc_int_cutoff: minimum percentage of pixels removed from threshold adjustment
nuc_int_cutoff_scaler: the scale factor for threshold adjustment as a function of nuclear intensity profile
max_nuc_int_cutoff: maximum percentage of pixels removed from threshold adjustment
min_mito_ffactor: minimum form factor for a mitotic LTR blob
min_mito_mal: minimum major axis length (diameter) for a mitotic LTR blob
max_mito_mal: maximum major axis length (diameter) for a mitotic LTR blob
min_mito_area: minimum area for a mitotic LTR blob
mito_LTR_int_above_bg: minimum LTR channel pixel grayscale intensity above background in 16-bit for mitotic cells
min_mito_LTR_coverage: minimum percent of nuclear region covered by LTR blobs for mitotic cells
fitc_above_bg: minimum NucView channel pixel grayscale intensity above background in 16-bit
fitc_join_radius: maximum distance in pixels between NucView spots getting joined together
min_fitc_area: minimum area of NucView spots for apoptotic cells
min_LTR_int_above_bg: minimum LTR grayscale intensity above background for detecting cytosol regions
max_dead_nuc_LTR_coverage: maximum percent of nuclear area coverage by LTR labeled cytosol regions for a late stage dead cell
min_dead_nuc_area: minimum area for a late stage dead cell

For output, three output folders need to be created with names "output1_images", "output2_images" and "output3_images". Output1_images folder would contain images based on the DAPI/Hoechst channel, with color-coded outline labeling each nucleus. White outline for interphase cell, red for mitotic cell, green for apoptotic cell, and yellow for late stage dead cell. Output2_images folder contains output images based on the NucView channel, with apoptotic objects outlined. Output3_images folder contains output images based on the LTR channel, with mitotic objects outlined.

Numeric data extracted from the image analysis is stored in two csv files. "Summary_site_data.csv" file contains data compiled from all the cells in each image. "Summary_well_data.csv" file contains data compiled from all the images acquired from the same well. Output data includes total cell count, interphase cell count, apoptotic cell count, late stage dead cell count, mitotic cell count, percentage of each phenotype, median value of interphase cell size, and median value of interphase cell form factor (roundness measure).
