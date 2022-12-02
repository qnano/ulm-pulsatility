General information
---
This software belongs to the article:\
Retrieving pulsatility in ultrasound localization microscopy\
DOI: https://doi.org/10.1109/OJUFFC.2022.3221354

Data availability
---
The simulation data underlying Figures 3, 4 and 6b in the article can be found at:\
DOI: https://doi.org/10.4121/21517878

Windows installation
---
The code was developed and tested in MATLAB 2021b, running on Windows 10.

Steps:

1. Clone repository or download software from https://github.com/qnano/ulm-pulsatility. Verify that the code (with the structure defined below) is contained in the directory "/ULM_pulsatility/".

2. Download MATLAB 2021b from: https://nl.mathworks.com/products/new_products/release2021b.html

3. Install the following MATLAB toolboxes:
    - Communications
    - Bioinformatics
    - Image Processing
    - Curve Fitting
    - Signal Processing
    - Statistics and Machine Learning
    - Parallel Computing
    - Computer Vision Toolbox


4. Download the data
    - Download the desired .zip folders from https://doi.org/10.4121/21517878. 
    - Create an empty folder in the directory "/ULM_pulsatility/" named "DataSets"
    - Extract the content (which is a folder with the same name as the .zip) to "/ULM_pulsatility/DataSets/".
    - Verify that the code (with the structure defined below) is contained in the directory "/ULM_pulsatility/" and the data sets in the directory "/ULM_pulsatility/DataSets/".

Data processing
---
We describe the procedure that is needed to reproduce Figures 3, 4 and 6b from the article from the available simulation data.

### The simulation data

The simulation data that is available on https://doi.org/10.4121/21517878 consists of the following files and .zip folders:

1. DataSetFig3. This dataset contains the simulation data needed to reproduce Fig. 3 of the article. 

2. DataSetR1. This dataset contains the simulation data corresponding to the lateral vessel location R1. It is needed to reproduce Fig 4 and Fig. 6b. 

3. DataSetR2. This dataset contains the simulation data corresponding to the lateral vessel location R2. It is needed to reproduce Fig 4 and Fig. 6b. 

4. DataSetR3. This dataset contains the simulation data corresponding to the lateral vessel location R3. It is needed to reproduce Fig 4 and Fig. 6b. 

### The repository

The repository consists of 3 main folders:

1. scripts. This folder contains the MATLAB code needed for reproduction of the mentioned figures of the article.

2. functions. This folder contains two MATLAB functions that are called in the MATLAB code. These functions are solely needed for visualization purposes.

3. PALA-master. This folder contains the MATLAB code of the LOTUS TOOLBOX, corresponding to the article: "Performance benchmarking of microbubble-localization algorithms for ultrasound localization microscopy", by Baptiste Heiles, Arthur Chavignon, Vincent Hingot, Pauline Lopez, Eliott Teston and Olivier Couture. Code by: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. For documentation see: https://github.com/AChavignon/PALA . Small alterations were made to accomodate for our application. Instead of using the function ULM_Track2MatOut.m the new function ULM_Track2MatOut_new.m is used. A new functionality was created that collects all velocity estimations in a single pixel and stores them in an array. The original function solely outputs the average of these estimations. This extention is necessary for the application of the Eulerian Method.

After installation of the simulation data, the "ULM_pulsatility" folder also contains the folder "DataSets" in which the data sets from https://doi.org/10.4121/21517878 are placed.
### Reproducing results

1. To produce the desired figures from the article, open the appropriate scripts:
   - ULM_processing.mat: Figure 3. This script contains the full ULM processing pipeline from input ultrasound images to output ULM reconstructions.
   - LagrangianApproach.mat: Figure 4. This script processes the input ultrasound images into ULM tracks similar as in ULM_processing.mat. Instead of processing these tracks into an image, this script applies the Lagrangian Approach as descriped in the article to obtain instantaneous velocity fluctuations and retrieve pulsatility fraction. 
   - EulerianApproach.mat: Figure 6b. This script processes the input ultrasound images into ULM tracks similar as in ULM_processing.mat. Instead of processing these tracks into an image, this script applies the Eulerian Approach as descriped in the article to obtain velocity distributions (i.e. histograms) and retrieve pulsatility fraction. 

2. Before running the desired MATLAB file verify that 
   - the DATA_folder is set to the path of the folder containing the downloaded data sets. Then select in the line below, the data set you want to use.
   - the REP_folder is set to the path of the repository.

3. Run the appropriate program(s) to process and plot the data.