%% Script used for results in Figure 3
% By: Myrthe Wiersma
% corresponding paper:  
% "Retrieving pulsatility in ultrasound localization microscopy"
% IEEE Open Journal of Ultrasonics, Ferroelectrics, and Frequency Control
% Latest update by Myrthe Wiersma on 2022/11/20
%
%
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
% - LOTUS_INPUT_SIM**:      A folder containing the input B-mode chunks
%       * chunk***.mat:     Data set containing IQ variable of size
%                           180x240x1000, containing 1000 B-mode images of 
%                           dimension 9mmx12mm
% - TOTAL_MB_*MB_<date>:    Contains MB_loc_conc --> list of all GT MB
%                           locations
%
% This script outputs:
% - All subfigures used to construct Fig. 3 of the paper 
%   (saved in SIM01/ULM_results) Only some formatting is needed to
%   reconstruct Fig. 3
% - Stores all important results in ULM_output_SIM**_<date>
%       * MatOut            ULM density reconstruction               
%       * MatOut_vel        ULM velocity reconstruction
%       * MatOut_theta      ULM orientation reconstruction
%       * MatOut_loc_GT     Accumulation of all GT MB locations
%       * ULM               Contains all ULM processing parameters
%       * Track_tot         Contains all found MB tracks
%       * Localizations     Contains all found MB localizations
%       * PowDop            Power Doppler image
%       * RMSE_localization RMSE of MB localization [m]
%       * errorx, errorz    All MB localization errors in x & z resp. [m]
%       * SE_vel_rec        Squared errors velocity reconstruction
%       * SE_theta_rec      Squared errors orientation reconstruction

%% 
clear all
close('all')

% Choose folder of data set to perform ULM reconstruction on
main_folder =  'C:\Users\Lab3\Documents\DataSets\DataSetFig3';              % Change this to the location of DataSetFig3 on your pc
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
Max_linking_distance = 0.6136;                                              % In [pixels]
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
total_nr_of_tracked_loc = 0;
interp_factor_tracking = 1/ULM.max_linking_distance/ULM.res*.8;
for track = 1:size(Track_tot,1)
    total_nr_of_tracked_loc = total_nr_of_tracked_loc + size(Track_tot{track,1},1)*interp_factor_tracking;
end

% -------- Find all GT MB locations in local patch coordinates --------
MB_loc_conc_tot = cat(1,MB_loc_conc{:})-[3.7e-3, 0.4e-3];
MB_loc_conc_tot_pix = MB_loc_conc_tot./US_pix_size.*ULM.scale(1:2)/ULM.SRscale+[1/2 1/2];

% ------ Transform MB_loc to same format as Localizations -----
for frame = 1:nr_frames
    
    if frame <= 30000
        nr_MB = 2;
    elseif frame > 30000
        nr_MB = 3;
    end
    
    chunk_nr = ceil(frame/chunk_size);
    Locs_GT = (MB_loc_conc{1,frame}./US_pix_size).*ULM.scale(1:2);
    Localizations_GT{1,chunk_nr}(((frame-(chunk_nr-1)*chunk_size)-1)*nr_MB+1:(frame-(chunk_nr-1)*chunk_size)*nr_MB,2:4) = [Locs_GT, ones(nr_MB,1)*(frame-(chunk_nr-1)*chunk_size)];
end

% -------- Convert tracks into SRpixel ---------------
Track_matout = cellfun(@(x) (x(:,[1 2 3 4]))./[ULM.SRscale ULM.SRscale 1 1] ,Track_tot,'UniformOutput',0);
Loc_GT_matout = cellfun(@(x) (x(:,[2 3 4]))./[ULM.SRscale ULM.SRscale 1] ,Localizations_GT,'UniformOutput',0);

% ---------- Accumulate tracks on the final MatOut grid ---------
fprintf('--- Creating images --- \n\n')

MatOut = ULM_Track2MatOut_new(Track_matout,ULM.SRsize);                     % tracks accumulated on supergrid [z x]
MatOut_loc_GT = ULM_Track2MatOut_new(Loc_GT_matout,ULM.SRsize);             % GT tracks accumulated on in supergrid [z x]
[MatOut_vel,MatOut_vel_distr] = ULM_Track2MatOut_new(Track_matout,ULM.SRsize,'mode','2D_velnorm'); 
% MatOut_vel -->        velocity map, where vel is averaged if multiple tracks that passed the pixel
% MatOut_vel_distr -->  Each pixel collects all velocity measurements of all tracks that passed the pixel

MatOut_vel = MatOut_vel*ULM.lambda;                                         % Convert into [m/s]
for i = 1:ULM.SRsize(1)
    for j = 1:ULM.SRsize(2)
        MatOut_vel_distr{i,j} = MatOut_vel_distr{i,j}*ULM.lambda;           % Convert into [m/s]
    end
end
MatOut_theta = ULM_Track2MatOut_new(Track_matout,ULM.SRsize, 'mode', '2D_theta'); % Orientation map on supergrid
MatOut_theta_or = MatOut_theta;
MatOut_vel_z = ULM_Track2MatOut_new(Track_matout,ULM.SRsize,'mode','2D_vel_z'); % Velocity map showing velocity in z direction supergrid [z x]
MatOut_vel_z = MatOut_vel_z*ULM.lambda;                                         % Convert into [m/s]

% -----------Wrap MatOut_theta between [0,2*pi] ---------------
% Additionally create MatOut_theta_sim that has NaN for empty pixels
% For plotting purpose later on
MatOut_theta_im = MatOut_theta-pi;  
for ii = 1:size(MatOut_theta,1)
    for jj = 1:size(MatOut_theta,2)
        if MatOut(ii,jj)==0
            MatOut_theta_im(ii,jj) = NaN;
        else
            MatOut_theta(ii,jj) = MatOut_theta(ii,jj)-pi;   
        end
        if MatOut_theta(ii,jj)>=2*pi
            MatOut_theta(ii,jj)=MatOut_theta(ii,jj)-2*pi;
            MatOut_theta_im(ii,jj)=MatOut_theta_im(ii,jj)-2*pi;
        elseif MatOut_theta(ii,jj)<0
            MatOut_theta(ii,jj)=MatOut_theta(ii,jj)+2*pi;
            MatOut_theta_im(ii,jj)=MatOut_theta_im(ii,jj)+2*pi;
        end
    end
end
toc

%% PowerDoppler rendering
if POWER_DOPPLER ==1
    tic
    fprintf('--- Generating Power Doppler rendering --- \n\n')
    PowDop = [];
    for hhh=1:min(20,Nbuffers)
        tmp = load([IQfiles(end-hhh).folder filesep IQfiles(end-hhh).name],'IQ');
        IQ_filt = SVDfilter(tmp.IQ,ULM.SVD_cutoff);tmp = [];
        PowDop(:,:,hhh) = sqrt(sum(abs(IQ_filt).^2,3));
    end
    PowDop_Image = mean(PowDop,3).^(1/3);
    fprintf('--- Done with Power Doppler rendering --- \n\n')
    toc
end

%% Metrics --> compute metrics from ULM results
% ------- localization precision --------
total_nr_of_loc = 0;
count = 0;
error_threshold = 100e-6;   
for i = 1: size(Localizations,2)
    total_nr_of_loc = total_nr_of_loc + size(Localizations{1,i},1);
    for f = 1:SizeOfBloc(3)
        frame_nr = (i-1)*SizeOfBloc(3)+f;
        idx = find(Localizations{1,i}(:,4)==f);
        Loc_frame_US = Localizations{1,i}(idx,2:3)./ULM.scale(1:2) + 1/2; 
        Loc_frame_m = (Loc_frame_US-1/2).*US_pix_size;
        if frame_nr == size(MB_loc_conc,2)
            continue
        end
        for loc = 1:size(Loc_frame_m,1)
            count = count+1;
            [errorxz(count), index] = min(sqrt(sum((Loc_frame_m(loc,:)-(MB_loc_conc{1,frame_nr+1}-patch_extremes(1,:))).^2,2)));
            errorz(count) = Loc_frame_m(loc,1)-(MB_loc_conc{1,frame_nr+1}(index,1)-patch_extremes(1,1));
            errorx(count) = Loc_frame_m(loc,2)-(MB_loc_conc{1,frame_nr+1}(index,2)-patch_extremes(1,2));

            if errorxz(count)>error_threshold                               % Then we don't want to take this localization into account for loc prec. calculation --> it will be regarded as a detection loss
                count = count-1;
            end
        end
    end
end

mean_errorz = mean(errorz);
std_errorz = std(errorz);
mean_errorx = mean(errorx);
std_errorx = std(errorx);
mean_error = mean(errorxz);
std_error = std(errorxz);

RMSE_localization = sqrt(mean(errorz.^2+errorx.^2));
nr_locs_fraction = total_nr_of_tracked_loc/total_nr_of_loc;

%% Extract patches & Plot the results figure


for PATCH = 1:nr_patches
    if PATCH < 10
        patch_name = 'patch0';
    else
        patch_name = 'patch';
    end

    % Extract the patch
    row = ceil(PATCH/4);                                                    % the patches are ordered as [ 1 2 3 4; 5 6 7 8; 9 10 11 12], here we only consider patch 1
    column = PATCH-(row-1)*4;   
    Vel_rec = MatOut_vel((row-1)*ULM_patch_size(1)+1:row*ULM_patch_size(1),(column-1)*ULM_patch_size(1)+1:column*ULM_patch_size(1));
    Theta_rec = MatOut_theta((row-1)*ULM_patch_size(1)+1:row*ULM_patch_size(1),(column-1)*ULM_patch_size(1)+1:column*ULM_patch_size(1));
    Theta_rec_im = MatOut_theta_im((row-1)*ULM_patch_size(1)+1:row*ULM_patch_size(1),(column-1)*ULM_patch_size(1)+1:column*ULM_patch_size(1));
    Dens_rec = MatOut((row-1)*ULM_patch_size(1)+1:row*ULM_patch_size(1),(column-1)*ULM_patch_size(1)+1:column*ULM_patch_size(1));
    Loc_rec_GT = MatOut_loc_GT((row-1)*ULM_patch_size(1)+1:row*ULM_patch_size(1),(column-1)*ULM_patch_size(1)+1:column*ULM_patch_size(1));

    % ------------- Metrics per patch ------------------------------
    % load GT average velocity and orientation map
    GT_theta_im = GT_theta;
    for i = 1:size(GT_theta,1)
        for j = 1:size(GT_theta,2)
            if GT_theta(i,j)==0
                GT_theta_im(i,j)=NaN;
            end
        end
    end

    % Compute vel & theta RMSE
    index_vel = 0;
    index_theta = 0;
    index_vel2 = 0;
    index_theta2 = 0;
    for index_z = 1:size(Vel_rec,1)
    for index_x = 1:size(Vel_rec,2)
        if Vel_rec(index_z,index_x)>0
            index_vel = index_vel+1;
            SE_vel_rec(index_vel) = (Vel_rec(index_z,index_x)-velAv(index_z,index_x)).^2;
        end
        if isnan(Theta_rec_im(index_z,index_x)) == 0
            index_theta = index_theta+1;
            abs_error = abs(Theta_rec(index_z,index_x)-GT_theta(index_z,index_x));
            if (abs_error >= 0) && (abs_error <= pi)
                SE_theta_rec(index_theta) = (abs_error).^2;
            elseif (abs_error > pi) && (abs_error <= 3*pi)
                SE_theta_rec(index_theta) = (2*pi-abs_error).^2;
            elseif (abs_error > 3*pi) 
                disp('caution abs error > 3*pi ')
                SE_theta_rec(index_theta) = (4*pi-abs_error).^2;
            else
                error('Theta error falls outside expected interval')
            end
        end
        
        if (Vel_rec(index_z,index_x)>0)&&(velAv(index_z,index_x)>0)         % ONLY consider true vessel locations
            index_vel2 = index_vel2+1;
            SE_vel_rec2(index_vel) = (Vel_rec(index_z,index_x)-velAv(index_z,index_x)).^2;
        end
        if (isnan(Theta_rec_im(index_z,index_x)) == 0) &&(isnan(GT_theta_im(index_z,index_x)) == 0)
            index_theta2 = index_theta2+1;
            abs_error = abs(Theta_rec(index_z,index_x)-GT_theta(index_z,index_x));
            if (abs_error >= 0) && (abs_error <= pi)
                SE_theta_rec2(index_theta) = (abs_error).^2;
            elseif (abs_error > pi) && (abs_error <= 3*pi)
                SE_theta_rec2(index_theta) = (2*pi-abs_error).^2;
            elseif (abs_error > 3*pi) 
                disp('caution abs error > 3*pi ')
                SE_theta_rec2(index_theta) = (4*pi-abs_error).^2;
            else
                error('Theta error falls outside expected interval')
            end
        end
    end
    end
    if exist('SE_vel_rec')
        RMSE_vel_rec = sqrt(mean(SE_vel_rec, 'all'));
    else
        RMSE_vel_rec = NaN;
    end

    if exist('SE_theta_rec')
        RMSE_theta_rec = sqrt(mean(SE_theta_rec, 'all'));
    else
        RMSE_theta_rec = NaN;
    end

    if exist('SE_vel_rec2')
        RMSE_vel_rec2 = sqrt(mean(SE_vel_rec2, 'all'));
    else
        RMSE_vel_rec2 = NaN;
    end

    if exist('SE_theta_rec2')
        RMSE_theta_rec2 = sqrt(mean(SE_theta_rec2, 'all'));
    else
        RMSE_theta_rec2 = NaN;
    end

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

    % ------------- Cross sections -----------------------
    if cross_section_PLOT ==1
    % Define pixels where cross sections need to be taken
    % Note the vessels are vertical!
    z_values = [51 250];
    distances = [45, 30, 15]*1E-6-10E-6;
    side = sqrt((1.5E-3)^2-(40E-6)^2);
    location_cross_section_m = 0.5E-3 + side - distances/(40E-6)*side;         % location along z axis where cross section should be taken
    location_cross_section_ULM = location_cross_section_m/ULM_pix_size(1)+0.5;
    location_cross_section_pix = location_cross_section_ULM+0.5;                                                          
    zz = round(location_cross_section_pix);                                 % z index where to take local cross-section
    x_values = [102 199];
    x_cross = (x_values(1)-1/2)*ULM_pix_size(2):ULM_pix_size(2):(x_values(2)-1/2)*ULM_pix_size(2);      % locations of cross section pixels along z [m]

  
    % extract Density map cross section
    cross_rec_dens = Dens_rec(z_values(1):z_values(2),x_values(1):x_values(2));
    cross_rec_dens_av = mean(cross_rec_dens,1);
    cross_rec_dens_av_norm = cross_rec_dens_av/max(cross_rec_dens_av);

    z_values_pd = [floor((z_values(1)-1)/Res+1), ceil(z_values(2)/Res)];
    x_values_pd = [floor((x_values(1)-1)/Res+1), ceil(x_values(2)/Res)];
    x_cross_pd = (x_values_pd(1))*US_pix_size(2):US_pix_size(2):(x_values_pd(2))*US_pix_size(2);

    PowDop_patch = PowDop_Image((row-1)*US_patch_size(1)+1:row*US_patch_size(1),(column-1)*US_patch_size(1)+1:column*US_patch_size(1));
    cross_rec_powdop = PowDop_patch(z_values_pd(1):z_values_pd(2),x_values_pd(1):x_values_pd(2));
    
    cross_rec_powdop_av = mean(cross_rec_powdop,1);
    cross_rec_powdop_av_norm = cross_rec_powdop_av/max(cross_rec_powdop_av);
    
    for i = 1:numel(zz) %size(zz,1) % Extract cross sections at specific location
        cross_rec_zz_dens(i,:) = Dens_rec(zz(i),x_values(1):x_values(2));
        cross_rec_zz_norm(i,:) = cross_rec_zz_dens(i,:)/max(cross_rec_zz_dens(i,:));
        zz_pd(i) = floor((zz(i)-1)/Res+1);
        cross_rec_powdop_zz(i,:) = PowDop_patch(zz_pd(i),x_values_pd(1):x_values_pd(2));
        cross_rec_powdop_zz_norm(i,:) = (cross_rec_powdop_zz(i,:)-min(cross_rec_powdop_zz(i,:)))/(max(cross_rec_powdop_zz(i,:))-min(cross_rec_powdop_zz(i,:)));
    end
    
    % ---------------- density cross section plot  ----------------------
    % IF YOU WANT SEPERATE VELOCITY CROSS SECTION PLOTS, UNCOMMENT THESE
%     figure
% 
%     for ves = 1:nr_ves
%         bounds = [points_CW{ves}(end-1,2) points_CW{ves}(end,2)]-patch_extremes(1,2);
%         if ves ==1
%             plot(bounds(1)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'DisplayName', 'GT vessel boundaries')
%         else
%             plot(bounds(1)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'HandleVisibility', 'Off')
%         end
%         hold on
%         plot(bounds(2)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'HandleVisibility', 'Off')
%     end
%     for i = 1:numel(zz) %size(zz,1)
%         plot(x_cross,cross_rec_zz_norm(i,:), 'Color', [0, 0.45, 0.64], 'LineStyle', '-', 'DisplayName', ['ULM'])
%         hold on
%     end
%     for j = 1:numel(zz) %size(zz,1)
%         plot(x_cross_pd,cross_rec_powdop_zz_norm(j,:), 'r', 'DisplayName', ['Power Doppler'])
%     end
%     legend()
%     title({['Normalized cross sections of density map']},'interpreter','latex')
%     xticks([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8]*1e-3)
%     xticklabels([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8])
%     axis([1.3*1E-3 1.7*1E-3 0 1.1])
%     xlabel('z [mm]')
%     ylabel('Normalized intensity [-]')
% 
%     savefig([Results_folder, filesep, 'Dens_cross_section_', patch_name, num2str(PATCH), '.fig'])

    % ------ Vel cross section ----------
    % extract Velocity map cross section
    cross_rec = Vel_rec(z_values(1):z_values(2),x_values(1):x_values(2));
   
    cross_rec_av = mean(cross_rec,1);
    cross_rec_zz = Vel_rec(zz(i),x_values(1):x_values(2));
    
    for i =1:length(zz)
         cross_rec_zz(i,:) = Vel_rec(zz(i),x_values(1):x_values(2));
    end
    % finding analytic ground truth velocity profile
    m_GT = [x_cross(1):1e-6:x_cross(end)]+patch_extremes(1,2);                  % variable along cross section to use for analytical solution [m]
    z_GT = [x_cross(1):1e-6:x_cross(end)];                                      % variable along cross section to use for analytical solution local patch coordinates [m]
    
    if nr_ves ==1
        for z_location = 1:length(zz)
            cross_GT{1} = PAR{1,1}{1,1}.u_analytic(m_GT);
            if ~(PAR{1,1}{1,1}.u_max - PAR{1,1}{ves,1}.u_min<1e-10)
                cross_GT_min{1}= PAR{1,1}{1,1}.u_analytic_min(m_GT);
                cross_GT_max{1} = PAR{1,1}{1,1}.u_analytic_max(m_GT);
            end
        end
    
    elseif   (nr_ves==2) && (parallel ==1) % two parallel vessels
        
        for z_location = 1:length(zz)
            % main vessel
            cross_GT{z_location,par_ves} = PAR{1,1}{par_ves,1}.u_analytic(m_GT);
            if ~(PAR{1,1}{par_ves,1}.u_max - PAR{1,1}{par_ves,1}.u_min<1e-10)
                cross_GT_min{z_location,par_ves}= PAR{1,1}{par_ves,1}.u_analytic_min(m_GT);
                cross_GT_max{z_location,par_ves}= PAR{1,1}{par_ves,1}.u_analytic_max(m_GT);
            end
            % second vessel
            cross_GT{z_location,sec_ves} = PAR{1,1}{sec_ves,1}.u_analytic(m_GT);
            if ~(PAR{1,1}{sec_ves,1}.u_max - PAR{1,1}{sec_ves,1}.u_min<1e-10)
                cross_GT_min{z_location,sec_ves}= PAR{1,1}{sec_ves,1}.u_analytic_min(m_GT);
                cross_GT_max{z_location,sec_ves}= PAR{1,1}{sec_ves,1}.u_analytic_max(m_GT);
            end
        end
        
    else            % Thus we have non parallel vessels --> A bit more effort to correctly get vessel boundaries.
        
        for z_location = 1:length(zz)
            % main vessel
            cross_GT{z_location,par_ves} = PAR{1,1}{par_ves,1}.u_analytic(m_GT);
            if ~(PAR{1,1}{par_ves,1}.u_max - PAR{1,1}{par_ves,1}.u_min<1e-10)
                cross_GT_min{z_location,par_ves}= PAR{1,1}{par_ves,1}.u_analytic_min(m_GT);
                cross_GT_max{z_location,par_ves}= PAR{1,1}{par_ves,1}.u_analytic_max(m_GT);
            end
            % second (tilted vessel)
            center_ves = PAR{1,1}{sec_ves,1}.start(2)-sin(PAR{1,1}{2,1}.d_theta)*(location_cross_section_m(z_location)-0.5E-3);
            tilted_m_GT = cos(PAR{1,1}{2,1}.d_theta)*(m_GT-center_ves);

            cross_GT{z_location,sec_ves} = PAR{1,1}{sec_ves,1}.u_analytic_local(tilted_m_GT);               % In this case the second vessel is tilted           
            if ~(PAR{1,1}{sec_ves,1}.u_max - PAR{1,1}{sec_ves,1}.u_min<1e-10)
                cross_GT_min{z_location,sec_ves}= PAR{1,1}{sec_ves,1}.u_analytic_min_local(tilted_m_GT);
                cross_GT_max{z_location,sec_ves}= PAR{1,1}{sec_ves,1}.u_analytic_max_local(tilted_m_GT);
            end
        end
    end

    % ---------------- velocity cross section plot  ----------------------
    % IF YOU WANT SEPERATE VELOCITY CROSS SECTION PLOTS, UNCOMMENT THESE
%     for z_location = 1:size(cross_GT,1) %for multiple cross-sections
%     figure
%     
%     plot(x_cross,cross_rec_zz(z_location,:), 'Color', [0, 0.45, 0.64], 'DisplayName', 'ULM')
%     hold on
%     for ves =1:nr_ves
%         if ves == par_ves
%             plot(z_GT,cross_GT{1,ves}, 'k', 'DisplayName', 'Analytic GT')
%             if ~(PAR{1,1}{ves,1}.u_max - PAR{1,1}{ves,1}.u_min<1e-10)
%                 plot(z_GT,cross_GT_min{1,ves}, 'k:', 'DisplayName', 'Analytic GT min/max')
%                 plot(z_GT,cross_GT_max{1,ves}, 'k:','HandleVisibility','off')
%             end
%         else
%             plot(z_GT,cross_GT{z_location,ves}, 'k', 'DisplayName', 'Analytic GT')
%             if ~(PAR{1,1}{ves,1}.u_max - PAR{1,1}{ves,1}.u_min<1e-10)
%                 plot(z_GT,cross_GT_min{z_location,ves}, 'k:', 'DisplayName', 'Analytic GT min/max')
%                 plot(z_GT,cross_GT_max{z_location,ves}, 'k:','HandleVisibility','off')
%             end
%         end
%     end
%     title({['Velocity along cross section ', num2str(z_location)]},'interpreter','latex')
%     xticks([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8]*1e-3)
%     xticklabels([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8])
%     axis([1.3*1E-3 1.7*1E-3 0 1.1*Largest_vel])
%     xlabel('z [mm]')
%     ylabel('velocity [m/s]')
% 
%     end
    
    % ------------- Reconstructions ----------------------

    ax_limits_ULM_x = [130.5 170.5];
    ax_limits_ULM_z = [45.5 255.5];
    ax_limits_PD_x = ((ax_limits_ULM_x-0.5)*ULM_pix_size(2))/US_pix_size(2)+0.5;
    ax_limits_PD_z = ((ax_limits_ULM_z-0.5)*ULM_pix_size(1))/US_pix_size(1)+0.5;

    % START DRAWING THE LARGE RESULTS FIGURE
    % GT vel------------------------------------------
    figure
    ax(1) = subplot(4,6,1)  ;
    imagesc(velAv*1E3)      % *1E3 for displaying mm/s
    hold on
    axis('equal')
    axis([ax_limits_ULM_x(1) ax_limits_ULM_x(2) ax_limits_ULM_z(1) ax_limits_ULM_z(2)])
    xticks([0.5 50.5 100.5 150.5 200.5 250.5 300.5])
    xticklabels([0 0.5 1 1.5 2 2.5 3])
    yticks([0.5 50.5 100.5 150.5 200.5 250.5 300.5])
    yticklabels([0 0.5 1 1.5 2 2.5 3])
    xlabel('x [mm]')
    ylabel('z [mm]')
    title({['ULM Velocity map']},'interpreter','latex')
    cmap = parula;
    cmap(1,:) = zeros(1,3);
    colormap(ax(1), cmap)
    c=colorbar;
    c.Label.String = 'velocity [mm/s]';

    % Snapshot of an IQ image with ROIs & Localizations
    ROI_size = 5;
    % Plot some localizations!
    IQ_FILE_NUMBER = 31;
    if PATCH == 1
        load([ IQfiles(IQ_FILE_NUMBER).folder filesep IQfiles(IQ_FILE_NUMBER).name])
        nr_fr_to_draw = 1;
        IQ_ten_frames = IQ(:,:,1:nr_fr_to_draw);
        clearvars IQ
        for frame=1:nr_fr_to_draw
            IQ_frame = IQ_ten_frames((row-1)*US_patch_size(1)+1:row*US_patch_size(1),(column-1)*US_patch_size(1)+1:column*US_patch_size(1),frame);
            idx_loc = find(Localizations{1,IQ_FILE_NUMBER}(:,4)==frame);
            Loc_frame_US = Localizations{1,IQ_FILE_NUMBER}(idx_loc,2:3)./ULM.scale(1:2) ;
            Loc_frame_m = (Loc_frame_US-1/2)./US_pix_size;
            nr_loc = size(Loc_frame_US,1);
    
            MB_loc_frame = MB_loc_conc{1,(IQ_FILE_NUMBER-1)*chunk_size+frame+1}-[patch_extremes(1,1) patch_extremes(1,2)];
            MB_loc_frame_US = MB_loc_frame./US_pix_size;
            ROI_center_pixels = ceil(MB_loc_frame_US - 0.5);
            ROI_indices = [];
            Noise_index_min = [];       % Area used to compute SNR
            Noise_index_max = [];       % Area used to compute SNR
            Noise_index = [];
            ROI_index = cell(1,nr_MB);
            for ROI=1:nr_MB
                if ROI > length(ROI_center_pixels)
                    break
                end
                max_idx = US_patch_size(1)*(ROI_center_pixels(ROI,2)-1)+ROI_center_pixels(ROI,1); 
                rows_im = ROI_center_pixels(ROI,1);
                col_im = ROI_center_pixels(ROI,2);
                if ROI_size ==3
                    ROI_index{ROI} = [ max_idx, max_idx+1, max_idx-1, max_idx-rows_im, max_idx-rows_im-1, max_idx-rows_im+1, max_idx+rows_im, max_idx+rows_im-1, max_idx+rows_im+1];
                    ROI_extremes{ROI} = [rows_im-1.5, col_im - 1.5; rows_im-1.5, col_im + 1.5; rows_im + 1.5 , col_im + 1.5; rows_im + 1.5, col_im - 1.5; rows_im-1.5, col_im - 1.5];
                elseif ROI_size ==5
                    for coll = col_im-2:1:col_im+2
                        for roww = rows_im-2:1:rows_im+2
                            ROI_index{ROI} = [ROI_index{ROI}, coll*US_patch_size(1)+roww];
                        end
                    end
                    ROI_extremes{ROI} = [rows_im-2.5, col_im - 2.5; rows_im-2.5, col_im + 2.5; rows_im + 2.5 , col_im + 2.5; rows_im + 2.5, col_im - 2.5; rows_im-2.5, col_im - 2.5];

                else
                    error('ROI not 3x3 or 5x5')
                end
                ROI_index{ROI} = sort(ROI_index{ROI});
                
                ROImin = ROI_index{ROI}(1:ROI_size);
                ROImax = ROI_index{ROI}(end-ROI_size+1:end);
                lower_col = floor(ROImin(1)/US_patch_size(1));
                upper_col = US_patch_size(1)-ceil(ROImax(1)/US_patch_size(1))-1;
                for colll = 2:lower_col
                    Noise_index_min = [Noise_index_min,ROImin-colll*US_patch_size(1)];
                end
                Noise_min_extremes{ROI} = [0.5, rem(ROImin(1),US_patch_size(1))-0.5; -1.5+lower_col, rem(ROImin(1),US_patch_size(1))-0.5; -1.5+lower_col,rem(ROImin(1),US_patch_size(1))+ROI_size-0.5; 0.5, rem(ROImin(1),US_patch_size(1))+ROI_size-0.5;0.5, rem(ROImin(1),US_patch_size(1))-0.5];

                for colll = 2:upper_col
                    Noise_index_max = [Noise_index_max,ROImin+colll*US_patch_size(1)];
                end
                Noise_max_extremes{ROI} = [ceil(ROImax(1)/US_patch_size(1))+0.5, rem(ROImin(1),US_patch_size(1))-0.5; US_patch_size(2)+0.5, rem(ROImin(1),US_patch_size(1))-0.5; US_patch_size(2)+0.5, rem(ROImin(1),US_patch_size(1))+ROI_size-0.5; ceil(ROImax(1)/US_patch_size(1))+0.5, rem(ROImin(1),US_patch_size(1))+ROI_size-0.5 ; ceil(ROImax(1)/US_patch_size(1))+0.5, rem(ROImin(1),US_patch_size(1))-0.5];
                Noise_index= [Noise_index, Noise_index_min, Noise_index_max];
                Noise_index = unique(Noise_index);
                Noise_index_min =[];
                Noise_index_max =[];

                ROI_indices = [ROI_indices,ROI_index{ROI}];
            end

            Signal_mean = mean(IQ_frame(ROI_indices));
            Noise_mean = mean(IQ_frame(Noise_index));
            
            SNR(frame) = mag2db(Signal_mean/Noise_mean);
            disp(['SNR of frame ', num2str(frame), ' isss ', num2str(SNR(frame)),' dB']);
            
            ax(2)=subplot(4,6,2);
            imagesc(mag2db(IQ_ten_frames(:,:,frame))- max(mag2db(IQ_ten_frames(:,:,frame)),[],'all'), [-40 0])
            
            colormap(ax(2),'gray')
            hold on
            plot(Loc_frame_US(:,2), Loc_frame_US(:,1), 'rx', 'DisplayName','Estimated MB location')
            plot(MB_loc_frame_US(:,2), MB_loc_frame_US(:,1), 'gx', 'DisplayName','Ground truth MB location')
            for roi=1:nr_MB
                if roi > length(ROI_center_pixels)
                    break
                end
                if roi ==1
                    plot(ROI_extremes{roi}(:,2),ROI_extremes{roi}(:,1), 'r', 'Displayname', 'Signal ROI')
                    %plot(Noise_min_extremes{roi}(:,1), Noise_min_extremes{roi}(:,2),'Color',[0, 0.45, 0.64],'Displayname', 'Noise ROI')
                    %plot(Noise_max_extremes{roi}(:,1), Noise_max_extremes{roi}(:,2),'Color',[0, 0.45, 0.64],'HandleVisibility', 'Off')

                else
                    plot(ROI_extremes{roi}(:,2),ROI_extremes{roi}(:,1), 'r', 'HandleVisibility', 'Off')
                    %plot(Noise_min_extremes{roi}(:,1), Noise_min_extremes{roi}(:,2),'Color',[0, 0.45, 0.64],'HandleVisibility', 'Off')
                    %plot(Noise_max_extremes{roi}(:,1), Noise_max_extremes{roi}(:,2),'Color',[0, 0.45, 0.64],'HandleVisibility', 'Off')
                end
            end
            legend()
            title({['US frame', num2str(frame)],['$\sigma$ = ', num2str(std_errorx*1E6), ' $\mu m$,  SNR = ', num2str(SNR(frame))]},'interpreter','latex')
            axis('equal')
            axis([(column-1)*ULM_patch_size(2)/Res+0.5+25, column*ULM_patch_size(2)/Res+.5-24, (row-1)*ULM_patch_size(1)/Res+0.5+1, row*ULM_patch_size(1)/Res+.5-5])
            xlabel('x [mm]')
            ylabel('z [mm]')
            zlabels = ([0:0.5e-3:3e-3])*1000;
            z_ticks = [0:0.5e-3:3e-3]/(US_pix_size(1))+0.5;
            yticks(z_ticks);
            yticklabels({num2str(zlabels(1)),num2str(zlabels(2)),num2str(zlabels(3)),num2str(zlabels(4)),num2str(zlabels(5)),num2str(zlabels(6)),num2str(zlabels(7))})

            xlabels = ([0:0.5e-3:3e-3])*1000;
            x_ticks = [0:0.5e-3:3e-3]/(US_pix_size(1))+0.5;
            xticks(x_ticks);
            xticklabels({num2str(xlabels(1)),num2str(xlabels(2)),num2str(xlabels(3)),num2str(xlabels(4)),num2str(xlabels(5)),num2str(xlabels(6)),num2str(xlabels(7))})
            
            c = colorbar;
            c.Label.String = 'Intensity dB';
        end
    end

    % PowDop dens ---------------------------------
    ax(3) = subplot(4,6,3);        
    imagesc(PowDop_Image);
    hold on
    axis('equal')
    axis([ax_limits_PD_x(1) ax_limits_PD_x(2) ax_limits_PD_z(1) ax_limits_PD_z(2)])
    xticks([0.5 50.5 100.5 150.5 200.5 250.5 300.5])
    xticklabels([0 0.5 1 1.5 2 2.5 3])
    yticks([0.5 50.5 100.5 150.5 200.5 250.5 300.5])
    yticklabels([0 0.5 1 1.5 2 2.5 3])
    ylabel('z [mm]')
    xlabel('x [mm]')
    title({['Power Doppler reconstruction']},'interpreter','latex')
 
    max_powdop = round(max(mean(PowDop,3).^(1/3),[],'all'));
    cmap = create_colormapHOT(max_powdop);
    colormap(ax(3),cmap)
    Colorbarticks = [0, max_powdop];
    ColorbarLabel{1,1}= num2str(0);
    ColorbarLabel{1,2}= num2str(1);
    
    c= colorbar('XTick', Colorbarticks, 'XTickLabel', ColorbarLabel);
    c.Label.String = 'Normalized intensity';

    % Plot the dens reconstruction ----------------------------------
    ax(4)=subplot(4,6,4)  ;       % ULM dens
    imagesc(Dens_rec)
    hold on
    axis('equal')
    axis([ax_limits_ULM_x(1) ax_limits_ULM_x(2) ax_limits_ULM_z(1) ax_limits_ULM_z(2)])
    xticks([0.5 50.5 100.5 150.5 200.5 250.5 300.5])
    xticklabels([0 0.5 1 1.5 2 2.5 3])
    yticks([0.5 50.5 100.5 150.5 200.5 250.5 300.5])
    yticklabels([0 0.5 1 1.5 2 2.5 3])
    xlabel('x [mm]')
    ylabel('z [mm]')
    legend()
    title({['ULM Density map']},'interpreter','latex')
    max_dens = round(max(Dens_rec,[],'all'));
    cmap = create_colormapHOT(max_dens);
    colormap(ax(4),cmap)
    if max_dens < 15
        Colorbarticks = ([0:1:max_dens] + 1/2)*max_dens/(max_dens+1);
    elseif max_dens<30
        Colorbarticks = ([0:2:max_dens] + 1/2)*max_dens/(max_dens+1);
    elseif max_dens<60
        Colorbarticks = ([0:4:max_dens] + 1/2)*max_dens/(max_dens+1);
    else 
        Colorbarticks = ([0:10:max_dens]+ 1/2)*max_dens/(max_dens+1);
    end
    for ticks = 1:length(Colorbarticks)
        ColorbarLabel{1,ticks}=num2str(Colorbarticks(ticks)*(max_dens+1)/max_dens-1/2);
    end
    c= colorbar('XTick', Colorbarticks, 'XTickLabel', ColorbarLabel);
    c.Label.String = 'count [-]';

    % ULM vel--------------------------------------
    ax(5) = subplot(4,6,5)  ;      
    
    imagesc(Vel_rec*1E3)
    hold on
    axis('equal')
    axis([ax_limits_ULM_x(1) ax_limits_ULM_x(2) ax_limits_ULM_z(1) ax_limits_ULM_z(2)])
    xticks([0.5 50.5 100.5 150.5 200.5 250.5 300.5])
    xticklabels([0 0.5 1 1.5 2 2.5 3])
    yticks([0.5 50.5 100.5 150.5 200.5 250.5 300.5])
    yticklabels([0 0.5 1 1.5 2 2.5 3])
    xlabel('x [mm]')
    ylabel('z [mm]')
    title({['ULM Velocity map']},'interpreter','latex')
    cmap = parula;
    cmap(1,:) = zeros(1,3);
    colormap(ax(5), cmap)
    c=colorbar;
    c.Label.String = 'velocity [mm/s]';

    % ULM theta-------------------------------------
    ax(6) = subplot(4,6,6)  ;      
    imagesc_theta(Theta_rec_im)
    hold on
    load([REP_folder filesep  'functions\ThetaColorMap.mat'])
    colormap(ax(6),myColorMap)
    axis on
    set(gca, 'Color', [0, 0, 0])
    axis([ax_limits_ULM_x(1) ax_limits_ULM_x(2) ax_limits_ULM_z(1) ax_limits_ULM_z(2)])
    xticks([0.5 50.5 100.5 150.5 200.5 250.5 300.5])
    xticklabels([0 0.5 1 1.5 2 2.5 3])
    yticks([0.5 50.5 100.5 150.5 200.5 250.5 300.5])
    yticklabels([0 0.5 1 1.5 2 2.5 3])
    xlabel('x [mm]')
    ylabel('z [mm]')
    title({['ULM orientation map ']},'interpreter','latex')
    c = colorbar;
    c.Label.String = 'CCW angle w.r.t x-axis [rad]';


    % cross section dens 1-------------------------
    z_location =1;
    ax(7) = subplot(4,6,7)   ;      
    for ves = 1:nr_ves
        if (ves == par_ves) || ((ves == sec_ves) && (parallel ==1))
            bounds = [points_CW{ves}(end-1,2) points_CW{ves}(end,2)]-patch_extremes(1,2);
            if ves ==1
                plot(bounds(1)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'DisplayName', 'GT vessel boundaries')
            else
                plot(bounds(1)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'HandleVisibility', 'Off')
            end
            hold on
            plot(bounds(2)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'HandleVisibility', 'Off')
        elseif (ves == sec_ves) && (parallel ==0)
            center_ves = (PAR{1,1}{sec_ves,1}.start(2)-patch_extremes(1,2))-sin(PAR{1,1}{2,1}.d_theta)*(location_cross_section_m(z_location)-0.5E-3);
            tilted_d = (PAR{1,1}{sec_ves,1}.d/2)/cos(PAR{1,1}{2,1}.d_theta);

            bounds = [center_ves - tilted_d , center_ves + tilted_d];
            if ves ==1
                plot(bounds(1)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'DisplayName', 'GT vessel boundaries')
            else
                plot(bounds(1)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'HandleVisibility', 'Off')
            end
            hold on
            plot(bounds(2)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'HandleVisibility', 'Off')
        end
    end
    histogram('BinEdges', [x_cross-ULM_pix_size(1)/2, x_cross(end)+ULM_pix_size(1)/2],'BinCounts', cross_rec_zz_norm(1,:), 'DisplayName', ['ULM'])
    plot(x_cross_pd,cross_rec_powdop_zz_norm(1,:), 'r', 'DisplayName', ['Power Doppler'])
    
    legend()
    title({['Normalized cross sections of density map']},'interpreter','latex')
    xticks([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8]*1e-3)
    xticklabels([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8])
    axis([1.3*1E-3 1.7*1E-3 0 1.1])
    xlabel('z [mm]')
    ylabel('Normalized intensity [-]')

    % cross section dens 2-------------------------
    z_location = 2;
    ax(8) = subplot(4,6,13)  ;
    for ves = 1:nr_ves
        if (ves == par_ves) || ((ves == sec_ves) && (parallel ==1))
            bounds = [points_CW{ves}(end-1,2) points_CW{ves}(end,2)]-patch_extremes(1,2);
            if ves ==1
                plot(bounds(1)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'DisplayName', 'GT vessel boundaries')
            else
                plot(bounds(1)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'HandleVisibility', 'Off')
            end
            hold on
            plot(bounds(2)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'HandleVisibility', 'Off')
        elseif (ves == sec_ves) && (parallel ==0)
            center_ves = (PAR{1,1}{sec_ves,1}.start(2)-patch_extremes(1,2))-sin(PAR{1,1}{2,1}.d_theta)*(location_cross_section_m(z_location)-0.5E-3);
            tilted_d = (PAR{1,1}{sec_ves,1}.d/2)/cos(PAR{1,1}{2,1}.d_theta);

            bounds = [center_ves - tilted_d , center_ves + tilted_d];
            if ves ==1
                plot(bounds(1)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'DisplayName', 'GT vessel boundaries')
            else
                plot(bounds(1)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'HandleVisibility', 'Off')
            end
            hold on
            plot(bounds(2)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'HandleVisibility', 'Off')
        end
    end
    histogram('BinEdges', [x_cross-ULM_pix_size(1)/2, x_cross(end)+ULM_pix_size(1)/2],'BinCounts', cross_rec_zz_norm(2,:), 'DisplayName', ['ULM'])
    plot(x_cross_pd,cross_rec_powdop_zz_norm(2,:), 'r', 'DisplayName', ['Power Doppler'])
    legend()
    title({['Normalized cross sections of density map']},'interpreter','latex')
    xticks([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8]*1e-3)
    xticklabels([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8])
    axis([1.3*1E-3 1.7*1E-3 0 1.1])
    xlabel('z [mm]')
    ylabel('Normalized intensity [-]')

    % cross section dens 3-------------------------
    z_location = 3;
    ax(9) = subplot(4,6,19)  ;
    for ves = 1:nr_ves
        if (ves == par_ves) || ((ves == sec_ves) && (parallel ==1))
            bounds = [points_CW{ves}(end-1,2) points_CW{ves}(end,2)]-patch_extremes(1,2);
            if ves ==1
                plot(bounds(1)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'DisplayName', 'GT vessel boundaries')
            else
                plot(bounds(1)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'HandleVisibility', 'Off')
            end
            hold on
            plot(bounds(2)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'HandleVisibility', 'Off')
        elseif (ves == sec_ves) && (parallel ==0)
            center_ves = (PAR{1,1}{sec_ves,1}.start(2)-patch_extremes(1,2))-sin(PAR{1,1}{2,1}.d_theta)*(location_cross_section_m(z_location)-0.5E-3);
            tilted_d = (PAR{1,1}{sec_ves,1}.d/2)/cos(PAR{1,1}{2,1}.d_theta);

            bounds = [center_ves - tilted_d , center_ves + tilted_d];
            if ves ==1
                plot(bounds(1)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'DisplayName', 'GT vessel boundaries')
            else
                plot(bounds(1)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'HandleVisibility', 'Off')
            end
            hold on
            plot(bounds(2)*[1 1],[0,1.1], 'Color', [0.14, 0.14, 0.14], 'LineStyle', '--', 'HandleVisibility', 'Off')
        end
    end
    histogram('BinEdges', [x_cross-ULM_pix_size(1)/2, x_cross(end)+ULM_pix_size(1)/2],'BinCounts', cross_rec_zz_norm(3,:), 'DisplayName', ['ULM'])
    plot(x_cross_pd,cross_rec_powdop_zz_norm(3,:), 'r', 'DisplayName', ['Power Doppler'])

    legend()
    title({['Normalized cross sections of density map']},'interpreter','latex')
    xticks([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8]*1e-3)
    xticklabels([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8])
    axis([1.3*1E-3 1.7*1E-3 0 1.1])
    xlabel('z [mm]')
    ylabel('Normalized intensity [-]')

   
    % cross section vel 1---------------------------
    z_location = 1;
    ax(10) = subplot(4,6,9) ;   
    plot(x_cross,cross_rec_zz(1,:), 'Color', [0, 0.45, 0.64], 'DisplayName', 'ULM')
    hold on
    for ves =1:nr_ves
        if ves == par_ves
            plot(z_GT,cross_GT{1,ves}, 'k', 'DisplayName', 'Analytic GT')
            if ~(PAR{1,1}{ves,1}.u_max - PAR{1,1}{ves,1}.u_min<1e-10)
                plot(z_GT,cross_GT_min{1,ves}, 'k:', 'DisplayName', 'Analytic GT min/max')
                plot(z_GT,cross_GT_max{1,ves}, 'k:','HandleVisibility','off')
            end
        else
            plot(z_GT,cross_GT{z_location,ves}, 'k', 'DisplayName', 'Analytic GT')
            if ~(PAR{1,1}{ves,1}.u_max - PAR{1,1}{ves,1}.u_min<1e-10)
                plot(z_GT,cross_GT_min{z_location,ves}, 'k:', 'DisplayName', 'Analytic GT min/max')
                plot(z_GT,cross_GT_max{z_location,ves}, 'k:','HandleVisibility','off')
            end
        end
    end
    title({['Velocity along cross section 1']},'interpreter','latex')
    xticks([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8]*1e-3)
    xticklabels([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8])
    axis([1.3*1E-3 1.7*1E-3 0 1.1*Largest_vel])

    xlabel('z [mm]')
    ylabel('velocity [m/s]')

    % cross section vel 2  ---------------------------- 
    z_location = 2;
    ax(11) = subplot(4,6,15) ;      
    plot(x_cross,cross_rec_zz(2,:), 'Color', [0, 0.45, 0.64], 'DisplayName', 'ULM')
    hold on
    for ves =1:nr_ves
        if ves == par_ves
            plot(z_GT,cross_GT{1,ves}, 'k', 'DisplayName', 'Analytic GT')
            if ~(PAR{1,1}{ves,1}.u_max - PAR{1,1}{ves,1}.u_min<1e-10)
                plot(z_GT,cross_GT_min{1,ves}, 'k:', 'DisplayName', 'Analytic GT min/max')
                plot(z_GT,cross_GT_max{1,ves}, 'k:','HandleVisibility','off')
            end
        else
            plot(z_GT,cross_GT{z_location,ves}, 'k', 'DisplayName', 'Analytic GT')
            if ~(PAR{1,1}{ves,1}.u_max - PAR{1,1}{ves,1}.u_min<1e-10)
                plot(z_GT,cross_GT_min{z_location,ves}, 'k:', 'DisplayName', 'Analytic GT min/max')
                plot(z_GT,cross_GT_max{z_location,ves}, 'k:','HandleVisibility','off')
            end
        end
    end
    title({['Velocity along cross section']},'interpreter','latex')
    xticks([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8]*1e-3)
    xticklabels([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8])
    axis([1.3*1E-3 1.7*1E-3 0 1.1*Largest_vel])

    xlabel('z [mm]')
    ylabel('velocity [m/s]')

    % cross section vel 3  ---------------------------- 
    z_location = 3;
    ax(12) = subplot(4,6,21) ;      
    plot(x_cross,cross_rec_zz(3,:), 'Color', [0, 0.45, 0.64], 'DisplayName', 'ULM')
    hold on
    for ves =1:nr_ves
        if ves == par_ves
            plot(z_GT,cross_GT{1,ves}, 'k', 'DisplayName', 'Analytic GT')
            if ~(PAR{1,1}{ves,1}.u_max - PAR{1,1}{ves,1}.u_min<1e-10)
                plot(z_GT,cross_GT_min{1,ves}, 'k:', 'DisplayName', 'Analytic GT min/max')
                plot(z_GT,cross_GT_max{1,ves}, 'k:','HandleVisibility','off')
            end
        else
            plot(z_GT,cross_GT{z_location,ves}, 'k', 'DisplayName', 'Analytic GT')
            if ~(PAR{1,1}{ves,1}.u_max - PAR{1,1}{ves,1}.u_min<1e-10)
                plot(z_GT,cross_GT_min{z_location,ves}, 'k:', 'DisplayName', 'Analytic GT min/max')
                plot(z_GT,cross_GT_max{z_location,ves}, 'k:','HandleVisibility','off')
            end
        end
    end
    title({['Velocity along cross section']},'interpreter','latex')
    xticks([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8]*1e-3)
    xticklabels([1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8])
    axis([1.3*1E-3 1.7*1E-3 0 1.1*Largest_vel])

    xlabel('z [mm]')
    ylabel('velocity [m/s]')

    
    % Histograms -------------------------------------
    ax(13)=subplot(4,6,11);
    histogram(errorx,'BinWidth',1E-6)
    title({['Histogram of error in x direction']}, 'interpreter', 'latex')
    xlabel('error x [micron]')
    ylabel('Count [-]')

    ax(14)=subplot(4,6,17);
    histogram(errorz,'BinWidth',1E-6)
    xlabel('error z [micron]')
    ylabel('Count [-]')
    title({['Histogram of error in z direction']}, 'interpreter', 'latex')

    % Save the results figure
    savefig([Results_folder, filesep, 'Large_results_figure_', patch_name, num2str(PATCH), '.fig'])

    end
end

save([Results_folder filesep 'ULM_output_' datestr(now,'yyyy_mm_dd_hh') '.mat'],'MatOut','MatOut_vel', 'MatOut_theta','MatOut_loc_GT','ULM','Track_tot','Localizations', 'PowDop', 'RMSE_localization', 'errorx', 'errorz', 'SE_vel_rec', 'SE_theta_rec')


%% Colormap function

function cmap = create_colormapHOT(max_value)

if ceil(max_value+1)==1
    cmap = [0 0 0];
elseif ceil(max_value+1)==2
    cmap = [ 0 0 0 ; 1 0 0 ];
elseif ceil(max_value+1)==3
    cmap = [ 0 0 0 ; 1 0 0 ; 1 1 0];
elseif ceil(max_value+1)==4
    cmap = [ 0 0 0 ; 1 0 0 ; 1 1 0; 1 1 1];
elseif ceil(max_value+1)==5
    cmap = [ 0 0 0 ; 0.5 0 0 ; 1 0 0 ; 1 1 0; 1 1 1];
elseif ceil(max_value+1)==6
    cmap = [ 0 0 0 ; 0.5 0 0 ; 1 0 0 ; 1 0.5 0; 1 1 0; 1 1 1];
elseif ceil(max_value+1)==7
    cmap = [ 0 0 0 ; 0.5 0 0 ; 1 0 0 ; 1 0.5 0; 1 1 0; 1 1 0.5; 1 1 1];
elseif ceil(max_value+1)==8
    cmap = [ 0 0 0 ; 0.3 0 0 ; 0.7 0 0; 1 0 0 ; 1 0.5 0; 1 1 0; 1 1 0.5; 1 1 1];
elseif ceil(max_value+1)>8
    cmaphot = hot(ceil(max_value));
    cmap = [0 0 0;cmaphot];
end

end

