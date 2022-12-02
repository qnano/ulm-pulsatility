%% Script used for results in Figure 3
% By: Myrthe Wiersma
% corresponding paper:  
% "Retrieving pulsatility in ultrasound localization microscopy"
% IEEE Open Journal of Ultrasonics, Ferroelectrics, and Frequency Control
% Latest update by Myrthe Wiersma on 2022/11/30


% This script takes an input data set containing:
% - PAR:                    A data struct that contains the simulation 
%                           parameters organized per vessel
%       * MB_nr_av:         Nr of MB simulated in the vessel
%       * start:            Start location vessel [m]      
%       * theta:            Orientation vessel CCW direction w.r.t. x-axis
%       * d:                Diameter vessel [m]
%       * l:                Length vessel [m]
%       * type:             Type of vessel 'V' or 'A' for vein/artery resp.
%       * dtheta:           Only for second vessel (in PAR{1,1}{2,1})
%                           Gives angle between the two vessels
%       * f_s:              Sampling frequency (frame rate) [Hz]
%       * f_p:              Ground truth (GT) pulsatility frequency [Hz]
%       * u_t_cp:           GT velocity over time of vessel centerline
%                           [m/s]
%       * u_analytic:       GT analytical expression of velocity profile
%                           over cross-section of vessel. (also _max and
%                           _min variant for max and min velocity profile
%                           during puls. cycle) [m/s]
%       * u_max & u_min:    Resp. GT maximum and minimum velocity of
%                           vessel centerline [m/s]
%       * u_cpix:           GT average velocity of centerline over full
%                           puls. cycle [m/s]
%       * frac_p:           GT pulsatility fraction [-]
%                           
% - Largest_vel:            Largest velocity present during the simulation
%                           (Mostly used for scaling the axes) [m/s]
%
% - Res:                    Super-resolution factor. Dictates ULM pixel size                    
%
% - LOTUS_INPUT:            A folder containing the input B-mode chunks
%       * chunk***.mat:     Data set containing IQ variable of size
%                           180x240x1000, containing 1000 B-mode images of 
%                           dimension 9mmx12mm
% - TOTAL_MB.mat:           Contains MB_loc_conc --> list of all GT MB
%                           locations

% This script outputs:
% - The subfigures used to construct Fig. 6b of the paper. Some additional
%   formatting is necessary.

% Requirements:
% The following MATLAB toolboxes should be installed:
% - Communications
% - Bioinformatics
% - Image Processing
% - Curve Fitting
% - Signal Processing
% - Statistics and Machine Learning
% - Parallel Computing
% - Computer Vision Toolbox

% For ULM processing the LOTUS TOOLBOX is used.
% # OPEN PLATFORM FOR ULTRASOUND LOCALIZATION MICROSCOPY: PERFORMANCE ASSESSMENT OF LOCALIZATION ALGORITHMS
% ![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)
% **AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot.
%
% You can alter the ULM processing parameters by adjusting the ULM struct
% Smooth factor can be adjusted in
% ...\ULM_pulsatility\PALA-master\PALA\PALA_addons\ULM_toolbox\ULM_tracking2D.m

% First part of the code is equivalent to ULM_processing.m


%% 
clear all
close('all')

% Choose folder of data set to perform ULM reconstruction on
% Choose folder of data set to perform ULM reconstruction on
Data_folder = 'C:\Users\Lab3\Documents\DataSets';                           % Change this to the location of the folder containing the datasets on your pc
main_folder = [Data_folder filesep 'DataSetR1'];                            % Select here which dataset you want to proces (DataSetR1 or DataSetR2 or DatasetR3)
REP_folder = 'C:\Users\Lab3\Documents\ULM_pulsatility';                     % Change this to the location of the repository on your pc
PATCHES = 1;                                                                % Is the simulation structured in patches? Keep to 1
                                                                            % The accompanied dataset: 
                                                                            % Only contains MBs in one 3mmx3mm patch
cross_section_PLOT = 1;                                                     % Do you want to draw the cross section? (0/1)
POWER_DOPPLER = 1;                                                          % Do you want to draw the Power Doppler equivalent? (0/1)

% Add folders
addpath(genpath(REP_folder))
addpath(genpath(main_folder))

%% Checking licenses and features
ToolBoxRequires = {'Communications','Bioinformatics','Image Processing','Curve Fitting','Signal Processing','Statistics and Machine Learning','Parallel Computing','Computer Vision Toolbox'};
err = 0;
for featureName=ToolBoxRequires
   IsInstalledToolbox = contains(struct2array(ver), featureName{1});
   if ~IsInstalledToolbox, warning([featureName{1} ' is missing']),err=1;end
end
if err,error('Toolbox are missing.');end;clear ToolBoxRequires featureName IsInstalledToolbox err

%% Start Run of ULM pipeline

% ___________________Load data & create subfolders_________________________
fprintf('--- Loading stuff --- \n\n')

% ------Load MB file ------
cd(main_folder)
find_MB = dir('TOTAL_MB*.*');
load(find_MB.name, 'MB_loc_conc')

% ------Load PAR file ------                                                % includes PAR struct, Largest_vel, Res, US_pix_size
find_PAR = dir('PAR*.*');
if length(find_PAR)>1
    error('Multiple PAR structs found')
end
load(find_PAR.name);
nr_patches = size(PAR,1);
nr_sim = size(PAR,2);
nr_ves = size(PAR{1,1},1);

% ------Load UserParam ------
find_US_par = dir('US_par*.*');
if isempty(find_US_par)
    US_pix_size = [50e-6, 50e-6];
else
load(find_US_par(1).name);
US_pix_size = UserParam.US_pix_size;
end

% ------Find nr of MBs per frame ------
nr_MB = 0;
for ves = 1:size(PAR{1,1},1)
    if isfield(PAR{1,1}{ves,1}, 'MB_nr_av')
        nr_MB = nr_MB + PAR{1,1}{ves,1}.MB_nr_av;
    end
end

% load GT average vel map
find_velAv = dir('velAv*.*');
load([find_velAv.name])

% load GT orientation
find_GT_orientation = dir('GT_theta*.*');
load(find_GT_orientation(1).name)

Results_folder = [main_folder filesep 'ULM_results'];
if not(isfolder(Results_folder))    
    mkdir(Results_folder)
end

%% PARAMETERS

% ------ Load first chunk of US frames ------
LOTUS_input_folder = [main_folder filesep  'LOTUS_INPUT'];
cd(LOTUS_input_folder)
IQfiles = dir('chunk*.*');
load([IQfiles(1).folder filesep IQfiles(1).name]);

% ------ general settings ------
if numel(US_pix_size)==1
    US_pix_size = [US_pix_size US_pix_size];                                % US pixel size equal in both directions
end
ULM_pix_size = US_pix_size./Res;                                            % pixel size of ULM reconstruction in z,x direction [m]
nr_frames = size(MB_loc_conc,2);
SizeOfBloc = size(IQ);                                                      % [nr pixels in z, nr pixels in x, chunk_size]
chunk_size = size(IQ,3);

% ------ patch settings ------                                              % (equal for all simulations)
% Process a patch of 3mmx3mm of the input IQ data.
patch_size = [3e-3, 3e-3];                                                  % size of patch [m]
US_patch_size = patch_size./US_pix_size;                                    % size of patch expressed in nr of US pixels
ULM_patch_size = patch_size./ULM_pix_size;                                  % size of patch expressed in nr of ULM pixels
patch_extremes = [3.7 0.4; 12.7 12.4]*1e-3;                                 % Location of patch in the IQ data (don't change)

% ------ Ground truth vessel boundaries ------
% An outline of the simulated vessels is constructed here from geometrical
% data in PAR
for ves = 1:nr_ves
points_CW{ves}(1,1:2) = [PAR{1,1}{ves,1}.start(1)-cos(PAR{1,1}{ves,1}.theta)*PAR{1,1}{ves,1}.d/2, PAR{1,1}{ves,1}.start(2)-sin(PAR{1,1}{ves,1}.theta)*PAR{1,1}{ves,1}.d/2 ];
points_CW{ves}(2,1:2) = [points_CW{ves}(1,1)-sin(PAR{1,1}{ves,1}.theta)*PAR{1,1}{ves,1}.l, points_CW{ves}(1,2)+cos(PAR{1,1}{ves,1}.theta)*PAR{1,1}{ves,1}.l ];
points_CW{ves}(3,1:2) = [points_CW{ves}(2,1)+cos(PAR{1,1}{ves,1}.theta)*PAR{1,1}{ves,1}.d, points_CW{ves}(2,2)+sin(PAR{1,1}{ves,1}.theta)*PAR{1,1}{ves,1}.d ];
points_CW{ves}(4,1:2) = [points_CW{ves}(3,1)+sin(PAR{1,1}{ves,1}.theta)*PAR{1,1}{ves,1}.l, points_CW{ves}(3,2)-cos(PAR{1,1}{ves,1}.theta)*PAR{1,1}{ves,1}.l ];
points_CW{ves}(5,1:2) = [points_CW{ves}(4,1)-cos(PAR{1,1}{ves,1}.theta)*PAR{1,1}{ves,1}.d, points_CW{ves}(4,2)-sin(PAR{1,1}{ves,1}.theta)*PAR{1,1}{ves,1}.d ];
points_CW_patch{ves} = points_CW{ves} - [patch_extremes(1,1), patch_extremes(1,2)];
points_CW_ULM_pix{ves} = points_CW_patch{ves}./ULM_pix_size +1/2;
end

% ------ ULM reconstruction parameters ------
% For documentation see LOTUS TOOLBOX
% # OPEN PLATFORM FOR ULTRASOUND LOCALIZATION MICROSCOPY: PERFORMANCE ASSESSMENT OF LOCALIZATION ALGORITHMS
% ![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)
% **AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot.

lambda = 1540/(17.8571*1E6);                                                % wavelength = speed of sound/frequency(17.85 MHz)
ScaleOfPixel = US_pix_size/lambda;                                          % pixel size in unit [wavelengths]
framerate = PAR{1,1}{1,1}.f_s;                                              % imaging framerate in [Hz]
Nbuffers = numel(IQfiles);                                                  % number of chunks to process
Max_linking_distance = 1;                                                   % In [pixels]
Max_linking_distance = Max_linking_distance*ScaleOfPixel(1);                % Convert to wavelengths
dt = 1/framerate;

if US_pix_size(1)<lambda                                                    % Rule of Thumb of setting the Detection parameters 
    FWHM = 3;
    NlocalMax = 3;
else
    FWHM = 3;
    NlocalMax = 5;
end
  
ULM = struct('numberOfParticles', nr_MB,...                                 % Number of particles per frame. (30-100)
    'size',[SizeOfBloc(1) SizeOfBloc(2) SizeOfBloc(3)],...                  % size of input data [nb_pixel_z nb_pixel_x nb_frame_per_bloc]
    'scale',[ScaleOfPixel 1/framerate],...                                  % Scale [z x dt], size of pixel in the scaling unit. (here, pixsize = 1*lambda)
    'res',Res,...                                                           % Resolution factor. Typically 10 for final image rendering at lambda/10.
    'SVD_cutoff',[5 SizeOfBloc(3)/4],...                                    % SVD filtering, to be adapted to your clutter/SNR levels
    'max_linking_distance',Max_linking_distance,...                         % Maximum linking distance between two frames to reject pairing, in pixels units (UF.scale(1)). (2-4 pixel).
    'min_length',5,...                                                      % Minimum allowed length of the tracks in time. (5-20 frames)
    'fwhm',[1 1]*FWHM,...                                                   % Size [pixel] of the mask for localization. (3x3 for pixel at lambda, 5x5 at lambda/2). [fmwhz fmwhx]
    'max_gap_closing', 0,...                                                % Allowed gap in microbubbles' pairing. (if you want to skip frames 0)
    'interp_factor',1/Res,...                                               % Interpfactor (decimation of tracks)
    'LocMethod','Radial'...                                                 % Select localization algorithm (WA,Interp,Radial,CurveFitting,NoLocalization)
    );
ULM.parameters.NLocalMax = NlocalMax;                                       % Safeguard on the number of maxLocal in the fwhm*fwhm grid (3 for fwhm=3, 7 for fwhm=5)
ULM.lambda = lambda;
ULM.SRscale = ULM.scale(1)/ULM.res;
ULM.SRsize = round(ULM.size(1:2).*ULM.scale(1:2)/ULM.SRscale);

%% ULM processing
% For documentation see LOTUS TOOLBOX
% # OPEN PLATFORM FOR ULTRASOUND LOCALIZATION MICROSCOPY: PERFORMANCE ASSESSMENT OF LOCALIZATION ALGORITHMS
% ![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)
% **AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot.
tic
fprintf('--- FILTERING, LOCALIZATION AND TRACKING --- \n')

parfor hhh = 1:min(999,Nbuffers) % can be used with parallel pool
    % --------- Loading and Filtering ----------
    tmp = load([IQfiles(hhh).folder filesep IQfiles(hhh).name],'IQ');
    IQ_filt = SVDfilter(tmp.IQ,ULM.SVD_cutoff);tmp = [];
    IQ_filtered{hhh}=IQ_filt;

    % --------- Localization ----------
    [MatTracking] = ULM_localization2D(abs(IQ_filt),ULM); IQ_filt=[];       % MatTracking 1st column is the intensity, 2nd 3rd are the coordinates of the bubble in pixel coordinates, 4th is the frame number
    MatTracking(:,2:3) = (MatTracking(:,2:3)).*ULM.scale(1:2);              % Now MatTracking is in [wavelengths]
    Localizations{hhh} = MatTracking;
    
    % --------- Tracking ----------
    Track_tot_i{hhh} = ULM_tracking2D(MatTracking,ULM);                     % Tracking algorithm (list of tracks)
    Track_tot{hhh} = Track_tot_i{hhh};
    Track_tot_raw_i{hhh} = ULM_tracking2D(MatTracking,ULM, 'nointerp');     % Tracking algorithm (list of tracks)
    Track_tot_raw{hhh} = Track_tot_raw_i{hhh};
    disp(['Finished chunk nr ', num2str(hhh)])
end
fprintf('\n')
Track_tot = cat(1,Track_tot{:});
Track_tot_raw = cat(1,Track_tot_raw{:});

%% Mapping
% -------- Find nr of tracked MB ----------
% % % total_nr_of_tracked_loc = 0;
% % % interp_factor_tracking = 1/ULM.max_linking_distance/ULM.res*.8;
% % % for track = 1:size(Track_tot,1)
% % %     total_nr_of_tracked_loc = total_nr_of_tracked_loc + size(Track_tot{track,1},1)*interp_factor_tracking;
% % % end

% -------- Find all GT MB locations in local patch coordinates --------
% % % MB_loc_conc_tot = cat(1,MB_loc_conc{:})-[3.7e-3, 0.4e-3];
% % % MB_loc_conc_tot_pix = MB_loc_conc_tot./US_pix_size.*ULM.scale(1:2)/ULM.SRscale+[1/2 1/2];

% ------ Transform MB_loc to same format as Localizations -----
% % % for frame = 1:nr_frames
% % %     chunk_nr = ceil(frame/chunk_size);
% % %     Locs_GT = (MB_loc_conc{1,frame}./US_pix_size).*ULM.scale(1:2);
% % %     Localizations_GT{1,chunk_nr}(((frame-(chunk_nr-1)*chunk_size)-1)*nr_MB+1:(frame-(chunk_nr-1)*chunk_size)*nr_MB,2:4) = [Locs_GT, ones(nr_MB,1)*(frame-(chunk_nr-1)*chunk_size)];
% % % end

% -------- Convert tracks into SRpixel ---------------
Track_matout = cellfun(@(x) (x(:,[1 2 3 4]))./[ULM.SRscale ULM.SRscale 1 1] ,Track_tot,'UniformOutput',0);
% % % Loc_GT_matout = cellfun(@(x) (x(:,[2 3 4]))./[ULM.SRscale ULM.SRscale 1] ,Localizations_GT,'UniformOutput',0);

% ---------- Accumulate tracks on the final MatOut grid ---------
fprintf('--- Creating images --- \n\n')
%%
% % % MatOut = ULM_Track2MatOut_new(Track_matout,ULM.SRsize);                     % tracks accumulated on supergrid [z x]
% % % MatOut_loc_GT = ULM_Track2MatOut_new(Loc_GT_matout,ULM.SRsize);             % GT tracks accumulated on in supergrid [z x]
[MatOut_vel,MatOut_vel_distr] = ULM_Track2MatOut_new(Track_matout,ULM.SRsize,'mode','2D_velnorm'); 
% MatOut_vel -->        velocity map, where vel is averaged if multiple tracks that passed the pixel
% MatOut_vel_distr -->  Each pixel collects all velocity measurements of all tracks that passed the pixel

MatOut_vel = MatOut_vel*ULM.lambda;                                         % Convert into [m/s]
Velocity_hist = [];
for i = 1:ULM.SRsize(1)
    for j = 1:ULM.SRsize(2)
        MatOut_vel_distr{i,j} = MatOut_vel_distr{i,j}*ULM.lambda;           % Convert into [m/s]
        Velocity_hist = [Velocity_hist,MatOut_vel_distr{i,j}];
    end
end
%%
figure
if main_folder(end-1:end)== 'R1'
    histogram(Velocity_hist*1000, 'BinWidth', 0.1, 'FaceColor', [0.36,0.27,0.88],   'Normalization',  'pdf')
elseif main_folder(end-1:end)== 'R2'
    histogram(Velocity_hist*1000, 'BinWidth', 0.1, 'FaceColor', [0.07,0.62,1.00],   'Normalization',  'pdf')
elseif main_folder(end-1:end)=='R3'
    histogram(Velocity_hist*1000, 'BinWidth', 0.1, 'FaceColor', [1,0,0],   'Normalization',  'pdf')
end



