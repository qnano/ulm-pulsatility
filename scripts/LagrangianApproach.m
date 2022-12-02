%% Script Pulsatility Retrieval Method 1: The Lagrangian approach
% By: Myrthe Wiersma
% corresponding paper:  
% "Retrieving pulsatility in ultrasound localization microscopy"
% IEEE Open Journal of Ultrasonics, Ferroelectrics, and Frequency Control
% Latest update by Myrthe Wiersma on 2022/11/30

% This script takes 3 input data sets corresponding to lateral locations 
% R1, R2, and R3 (see paper). Each of the three data sets contains the 
% the following variables that will be used in this script:

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
% - TOTAL_MB:               Contains MB_loc_conc --> list of all GT MB
%                           locations

% This script outputs:
% - All subfigures used to construct Fig. 4 of the paper 
%   (saved in SIM01/ULM_results) Only some formatting is needed to
%   reconstruct Fig. 4

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
Data_folder = 'C:\Users\Lab3\Documents\DataSets';                           % Change this to the location of the folder containing the datasets on your pc
Folders = cell(3,1);
Folders{1,1} =  [ Data_folder filesep 'DataSetR1'];              
Folders{2,1} =  [ Data_folder filesep 'DataSetR2'];               
Folders{3,1} =  [ Data_folder filesep 'DataSetR3'];               

REP_folder = 'C:\Users\Lab3\Documents\ULM_pulsatility';                     % Change this to the location of the repository on your pc
DRAW_VEL = 0;

PATCHES = 1;                                                                % Is the simulation structured in patches? Keep to 1
                                                                            % The accompanied dataset: 
                                                                            % Only contains MBs in one 3mmx3mm patch
P_f_estimates = cell(3,1);                                                  % Saves Pf estimates of the 3 data sets
Trajectory = cell(3,2);                                                     % Saves the trajectory to plot 
for R = 1:3
% if main_folder(end-1:end) == 'R3'
%     SECOND_VESSEL = 1;
% else
%     SECOND_VESSEL =0;
% end
if R ==3
    SECOND_VESSEL = 1;
else
    SECOND_VESSEL = 0;
end

main_folder = Folders{R,1};
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
Res = 10;

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
Max_linking_distance = 0.4;                                                 % In [wavelengths]

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

%% Extract patches
smooth_factor =51;
nr_patches
for PATCH = 1:nr_patches
    if PATCH < 10
        patch_name = 'patch0';
    else
        patch_name = 'patch';
    end

    % Extract the patch
    row = ceil(PATCH/4);                                                    % the patches are ordered as [ 1 2 3 4; 5 6 7 8; 9 10 11 12]
    column = PATCH-(row-1)*4;

    % -----------------
    if size(PAR{1,1},1)>1
        if PAR{1,1}{2,1}.type == 'V'
            vel_t_cp = PAR{1,1}{2,1}.u_t_cp;
            par_ves = 2;
            vel_t_cp_sec_ves = PAR{1,1}{1,1}.u_t_cp;
            sec_ves = 1;
            if isfield(PAR{1,1}{2,1}, 'd_theta')
                parallel = 0;
            else
                parallel = 1;
            end
        else
            vel_t_cp = PAR{1,1}{1,1}.u_t_cp;
            par_ves = 1;
            vel_t_cp_sec_ves = PAR{1,1}{2,1}.u_t_cp;
            sec_ves = 2;
            if isfield(PAR{1,1}{2,1}, 'd_theta')
                parallel = 0;
            else
                parallel = 1;
            end
        end
    else
        vel_t_cp = PAR{1,1}{1,1}.u_t_cp;
        par_ves = 1;
    end

    % ------- plotting index ----

    ax_limits_ULM_x = [130.5 170.5];
    ax_limits_ULM_z = [45.5 255.5];
    ax_limits_PD_x = ((ax_limits_ULM_x-0.5)*ULM_pix_size(2))/US_pix_size(2)+0.5;
    ax_limits_PD_z = ((ax_limits_ULM_z-0.5)*ULM_pix_size(1))/US_pix_size(1)+0.5;

    % ------------- Retrieve velocity of time measurements ----------------
 
    for trackk = 1:size(Track_tot,1)
        Track_center = Track_tot_raw{trackk,1};
    
        for i = 1:size(Track_center,1)-1
            Vel_z_Track_center{1,trackk}(i) = (Track_center(i+1,1)-Track_center(i,1))*framerate*lambda*1E3;

        end

        Vel_z_center_smooth{1,trackk} = abs(smooth(Vel_z_Track_center{1,trackk},smooth_factor));
        Track_center_frames{1,trackk} = Track_center(1:end-1,3)'; 
    end
        
    % ------------------ Find GT vel ---------------------
    center_location = 0;% mean(Track_center(:,2))*lambda-(PAR{1,1}{par_ves,1}.start(2)-patch_extremes(1,2));
    if main_folder(end-1:end) == 'R2'
    center_location = 4e-5;
    end
    for ff = 1:length(vel_t_cp)
        if SECOND_VESSEL
        vel_t_GT(ff) = -vel_t_cp_sec_ves(ff)/(PAR{1,1}{sec_ves,1}.d/2)^2 *(center_location^2 - (PAR{1,1}{sec_ves,1}.d/2)^2);
        else
        vel_t_GT(ff) = -vel_t_cp(ff)/(PAR{1,1}{par_ves,1}.d/2)^2 *(center_location^2 - (PAR{1,1}{par_ves,1}.d/2)^2);
        end
    end
    %---------------------
    Single_track_raw = Track_tot_raw;
end



%% Filtering of velocity meurements

v_max = max(vel_t_GT);
v_min = min(vel_t_GT);
v_av = mean(vel_t_GT);

Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
chunk_nr=1;
chunk_tracks_summed = 0;
Max_vel_values = [];
Max_vel_frames = [];
Min_vel_values = [];
Min_vel_frames = [];

for TRACK = 1:size(Single_track_raw,1)
Vel_raw = Vel_z_Track_center{1,TRACK};
L = length(Vel_raw);
if TRACK<=size(Track_tot_raw_i{1,chunk_nr},1)+chunk_tracks_summed

else
    chunk_tracks_summed = chunk_tracks_summed+size(Track_tot_raw_i{1,chunk_nr},1);
    chunk_nr = chunk_nr +1;
    
end
frames = Single_track_raw{TRACK,1}(:,3) + (chunk_nr-1)*chunk_size;
frames = frames(1:end-1);
t = (0:L-1)*T;

% Calculate spectrum of raw velocity estimates
Spectrum = fft(Vel_raw);
P2 = abs(Spectrum/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

% Find f_dx area
f_max = 1/(US_pix_size(1)/v_max);
f_min = 1/(US_pix_size(1)/v_min);
f_av = 1/(US_pix_size(1)/v_av);
f_bw = f_max-f_min;

% SMOOTHING BY MOVING AVERAGE --------------------------------------------

Vel_smooth = abs(smooth(Vel_raw,smooth_factor));

v_smooth{TRACK}=Vel_smooth(ceil(smooth_factor*2/5+0.5):end-ceil(smooth_factor*2/5));
frames_smooth{TRACK}=frames(ceil(smooth_factor*2/5+0.5):end-ceil(smooth_factor*2/5));

% Band stop filter -------------------------------------------------------
vel_filtered = bandstop(Vel_raw,[f_min-0.2*f_bw,f_max+0.2*f_bw],Fs);
if f_bw < 15
    Wo = f_av/(Fs/2);  BW = Wo/2;
    [b,a] = iirnotch(Wo,BW); 
    vel_filtered = filter(b,a,Vel_raw);
end

Spectrum_filt = fft(vel_filtered);
P2_filt = abs(Spectrum_filt/L);
P1_filt = P2_filt(1:L/2+1);
P1_filt(2:end-1) = 2*P1_filt(2:end-1);
f = Fs*(0:(L/2))/L;

% FILTER + SMOOTHING ------------------------------------------------------
Vel_filt_smooth = abs(smooth(vel_filtered,smooth_factor));

v_filt_smooth{TRACK}=Vel_filt_smooth(ceil(smooth_factor*3/5+0.5):end-ceil(smooth_factor*3/5));
frames_filt_smooth{TRACK}=frames(ceil(smooth_factor*3/5+0.5):end-ceil(smooth_factor*3/5));

% Plotting full data set
if DRAW_VEL ==1
    if R == 1
        figure(1)
        hold on
        plot(frames,-Vel_raw)
    elseif R == 2
        figure(11)
        hold on
        plot(frames,-Vel_raw)
    else 
        figure(21)
        hold on
        plot(frames,Vel_raw)
    end
    
    if TRACK ==1
    plot([0:1:length(vel_t_GT)-1],vel_t_GT*1000)
    end
    % plot(frames,vel_t_GT(frames)*1000)
    title('Velocity measurements')
    xlabel('frames')
    ylabel('v [mm/s]')

    if R == 1
        figure(2)
    elseif R == 2
        figure(12)
    else 
        figure(22)
    end
    hold on
    if length(frames)>300       % Only draw spectrum of long tracks
        fill([f_min, f_min,f_max,f_max],[0, max(P1)*1.5, max(P1*1.5), 0], [0.80,0.80,0.80], 'EdgeColor','none')
        hold on
        plot(f,P1) 
    end
    title('Spectrum of raw velocity')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')

    if R == 1
        figure(3)
    elseif R == 2
        figure(13)
    else 
        figure(23)
    end
    hold on
    plot(frames_smooth{TRACK},v_smooth{TRACK})
    hold on
    if TRACK ==1
    plot([0:1:length(vel_t_GT)-1],vel_t_GT*1000)
    end
    title('Smoothed velocity measurements')
    xlabel('frames')
    ylabel('v [mm/s]')

    if R == 1
        figure(4)
    elseif R == 2
        figure(14)
    else 
        figure(24)
    end
    hold on
    plot(frames_filt_smooth{TRACK},v_filt_smooth{TRACK})
    hold on
    if TRACK ==1
    plot([0:1:length(vel_t_GT)-1],vel_t_GT*1000)
    end
    title('Filtered + smoothed velocity measurements')
    xlabel('frames')
    ylabel('v [mm/s]')

end

% CREATE LARGE RESULTS FIGURE
figure(200)
hold on
if R==1
    subplot(3,4,3)
    hold on
    if TRACK ==1
        plot([0:1:length(vel_t_GT)-1],vel_t_GT*1000, 'k')
    end
    if (frames(end)>1000) && (frames(1)<2000)
        subplot(3,4,3)
        hold on
        plot(frames_filt_smooth{TRACK},v_filt_smooth{TRACK},'Color', [0.36,0.27,0.88])
        hold on
        
        title('Filtered + smoothed velocity measurements')
        xlabel('frames')
        ylabel('v [mm/s]')
        axis([1000 2000 7.8 10.2])
    end
end
if R==2
    if frames(1)<600
        subplot(3,4,6)
        hold on
        plot(frames_smooth{TRACK},v_smooth{TRACK},  'Color', [0.07,0.62,1.00])
        hold on
        if TRACK ==1
            plot([0:1:length(vel_t_GT)-1],vel_t_GT*1000,  'k')
        end
        title('Smoothed velocity measurements')
        xlabel('frames')
        ylabel('v [mm/s]')
        axis([0 600 2 5.2])
        
        if length(frames)>300   % Only show spectrum of long tracks
            subplot(3,4,10)
            hold on
            fill([f_min, f_min,f_max,f_max],[0, max(P1)*1.5, max(P1*1.5), 0], [0.80,0.80,0.80], 'EdgeColor','none')
            hold on
            plot(f,P1, 'Color', [0.07,0.62,1.00]) 
            title('Spectrum of raw velocity')
            xlabel('f (Hz)')
            ylabel('|P1(f)|')
            axis([0 150 0 4])
        end
    end
    subplot(3,4,7)
    hold on
    if TRACK ==1
        plot([0:1:length(vel_t_GT)-1],vel_t_GT*1000, 'k')
    end
    if frames(1)<1000
        subplot(3,4,7)
        hold on
        plot(frames_filt_smooth{TRACK},v_filt_smooth{TRACK},  'Color', [0.07,0.62,1.00])
        hold on
        
        title('Filtered + smoothed velocity measurements')
        xlabel('frames')
        ylabel('v [mm/s]')
        axis([0 1000 2 4.5])
    end
    if TRACK == 4
        Trajectory{2,1} = Track_tot_raw{TRACK,1}(:,2)*lambda*1000;
        Trajectory{2,2} = Track_tot_raw{TRACK,1}(:,1)*lambda*1000;
    end
end
if R==3
    subplot(3,4,11)
    hold on
    if TRACK ==1
       plot([0:1:length(vel_t_GT)-1],vel_t_GT*1000, 'k')
    end
    if (frames(end)>4200) && ( frames(1)<5200)
        subplot(3,4,11)
        hold on
        plot(frames_filt_smooth{TRACK},v_filt_smooth{TRACK}, 'Color', [1.00,0.00,0.00])
        hold on
        
        title('Filtered + smoothed velocity measurements')
        xlabel('frames')
        ylabel('v [mm/s]')
        axis([4200 5200 6.5 11])
    end
end




if( R ==1) && (TRACK==1)
load([Data_folder filesep   'DataSetFig3\RAW_tracks.mat'])

subplot(3,4,1)
hold on
for i =1:size(Tracks_full_acq,1)
    if i ==1
        plot(Tracks_full_acq{i,1}(:,2)*ULM.lambda*1000,Tracks_full_acq{i,1}(:,1)*ULM.lambda*1000, 'Color', [0.7,0.7,0.7], 'DisplayName', 'ULM tracks')
    else
        plot(Tracks_full_acq{i,1}(:,2)*ULM.lambda*1000,Tracks_full_acq{i,1}(:,1)*ULM.lambda*1000, 'Color', [0.7,0.7,0.7],  'HandleVisibility', 'Off')

    end
    hold on
end

hold on
axis('equal')
set(gca, 'YDir', 'reverse' )
for vess = 1:nr_ves
   if vess ==1
        plot(points_CW_patch{1,vess}(:,2)*1000,points_CW_patch{1,vess}(:,1)*1000, 'Color','k', 'DisplayName', 'GT vessel boundary')
   else
        plot(points_CW_patch{1,vess}(:,2)*1000,points_CW_patch{1,vess}(:,1)*1000, 'Color','k', 'HandleVisibility','Off')

   end
end


subplot(3,4,2)
hold on

for i =1:size(Tracks_full_acq,1)
    if i ==1
        plot(Tracks_full_acq{i,1}(:,2)*ULM.lambda*1000,Tracks_full_acq{i,1}(:,1)*ULM.lambda*1000, 'Color', [0.7,0.7,0.7], 'DisplayName', 'ULM tracks')
    else
        plot(Tracks_full_acq{i,1}(:,2)*ULM.lambda*1000,Tracks_full_acq{i,1}(:,1)*ULM.lambda*1000, 'Color', [0.7,0.7,0.7],  'HandleVisibility', 'Off')

    end
    hold on
end

hold on
axis('equal')
axis([1.435 1.64  0.48 0.58])
set(gca, 'YDir', 'reverse' )
for vess = 1:nr_ves
   if vess ==1
        plot(points_CW_patch{1,vess}(:,2)*1000,points_CW_patch{1,vess}(:,1)*1000, 'Color','k', 'DisplayName', 'GT vessel boundary')
   else
        plot(points_CW_patch{1,vess}(:,2)*1000,points_CW_patch{1,vess}(:,1)*1000, 'Color','k', 'HandleVisibility','Off')

   end
end
end

if (R==1)&&(TRACK==1)
    subplot(3,4,1)
    hold on
    plot(Track_tot_raw{TRACK,1}(:,2)*lambda*1000,Track_tot_raw{TRACK,1}(:,1)*lambda*1000,  'Color', [0.36,0.27,0.88], 'DisplayName', 'R1')
    subplot(3,4,2)
    hold on
    plot(Track_tot_raw{TRACK,1}(:,2)*lambda*1000,Track_tot_raw{TRACK,1}(:,1)*lambda*1000,  'Color', [0.36,0.27,0.88], 'DisplayName', 'R1')
    legend()
end
if (R==2)&&(TRACK==4)
    subplot(3,4,1)
    hold on
    plot(Track_tot_raw{TRACK,1}(:,2)*lambda*1000,Track_tot_raw{TRACK,1}(:,1)*lambda*1000,  'Color', [0.07,0.62,1.00], 'DisplayName', 'R2')
    subplot(3,4,2)
    hold on
    plot(Track_tot_raw{TRACK,1}(:,2)*lambda*1000,Track_tot_raw{TRACK,1}(:,1)*lambda*1000,  'Color', [0.07,0.62,1.00], 'DisplayName', 'R2')
    legend()
end
if (R==3)&&(TRACK==10)
    subplot(3,4,1)
    hold on
    plot(Track_tot_raw{TRACK,1}(:,2)*lambda*1000,Track_tot_raw{TRACK,1}(:,1)*lambda*1000,  'Color', [1,0,0], 'DisplayName', 'R3')
    subplot(3,4,2)
    hold on
    plot(Track_tot_raw{TRACK,1}(:,2)*lambda*1000,Track_tot_raw{TRACK,1}(:,1)*lambda*1000,  'Color', [1,0,0], 'DisplayName', 'R3')
    legend()
end

end

%%

Max_vel_values = [];
Max_vel_frames = [];
Min_vel_values = [];
Min_vel_frames = [];
for TRACK = 1:size(Single_track_raw,1)
    nr_periods = ceil(size(v_filt_smooth{TRACK},1)/200);   % 200 frames in one cycle 
    for i = 1:nr_periods
        if i == nr_periods %then don't do 200 frames
            v_filtered_cycle = v_filt_smooth{TRACK}((i-1)*200+1:end,1);
            frames_filtered_cycle = frames_filt_smooth{TRACK}((i-1)*200+1:end,1);
        else
            v_filtered_cycle = v_filt_smooth{TRACK}((i-1)*200+1:i*200,1);
            frames_filtered_cycle = frames_filt_smooth{TRACK}((i-1)*200+1:i*200,1);
        end
        if isempty(v_filtered_cycle)==0
        [MAX, INDEX] = max(v_filtered_cycle);
        if (INDEX>10) && (INDEX<length(v_filtered_cycle)-10)
            Max_3_av = (v_filtered_cycle(INDEX-1)+ v_filtered_cycle(INDEX) + v_filtered_cycle(INDEX+1))/3;
            Max_vel_values = [Max_vel_values,Max_3_av];
            Max_vel_frames = [Max_vel_frames,frames_filtered_cycle(INDEX)];
        end

        [MIN, INDEX_MIN] = min(v_filtered_cycle);
        if (INDEX_MIN>10) && (INDEX_MIN<length(v_filtered_cycle)-10)
            Min_3_av = (v_filtered_cycle(INDEX_MIN-1)+ v_filtered_cycle(INDEX_MIN) + v_filtered_cycle(INDEX_MIN+1))/3;
            Min_vel_values = [Min_vel_values,Min_3_av];
            Min_vel_frames = [Min_vel_frames,frames_filtered_cycle(INDEX_MIN)];
        end
        end
    end

end

Min_frames = Min_vel_frames;
Max_frames = Max_vel_frames;
P_f_values = [];
if length(Min_vel_frames)<=length(Max_vel_frames)
    for i = 1:length(Min_vel_frames)
        min_f = Min_frames(i);
        count =0;
        for j = 1:length(Max_vel_frames)
            if (Max_frames(j)>min_f) 
                max_index = j;
                break
            end
        end
        min_f_value = Min_vel_values(i);
        max_f_value = Max_vel_values(max_index);
        P_f = (max_f_value-min_f_value)/((max_f_value+min_f_value)/2);
        if P_f>0
            P_f_values = [P_f_values,P_f];
        end
    end
else
    for i = 1:length(Max_vel_frames)
        max_f = Max_frames(i);
        for j = 1:length(Min_vel_frames)
            if (Min_frames(j)>max_f) 
                min_index = j;
                break
            end
        end
        max_f_value = Max_vel_values(i);
        min_f_value = Min_vel_values(min_index);
        P_f = (max_f_value-min_f_value)/((max_f_value+min_f_value)/2);
        if P_f>0
            P_f_values = [P_f_values,P_f];
        end
    end
    
end

figure
boxplot(P_f_values)
P_f_estimates{R,1}=P_f_values;


clearvars -except P_f_estimates Folders REP_folder PATCHES DRAW_VEL Data_folder ULM
end

%% Plot the boxplots

Pf_R1 = P_f_estimates{1,1}';
Pf_R2 = P_f_estimates{2,1}';
Pf_R3 = P_f_estimates{3,1}';

group = [ ones(size(Pf_R1)); 2*ones(size(Pf_R2)); 3*ones(size(Pf_R3))];
figure(200)
subplot(3,4,4)
boxplot([Pf_R1;Pf_R2;Pf_R3],group)
set(gca,  'XTickLabel', {'R1', 'R2', 'R3'})
