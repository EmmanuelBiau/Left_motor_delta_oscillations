%% I- CREATING HEADMODEL with an mri template FOR EACH PARTICIPANT WITH ITS POLHEMUS, FOR FURTHER SOURCE LOCALIZATION ANALYSIS;
clearvars;clc; addpath C:\toolbox\fieldtrip\; ft_defaults;
close all

% Load variables for the head model; 
cd C:\toolbox\fieldtrip\template\headmodel;
load standard_bem;   %load variable 'vol'
load standard_mri;   %load variable 'bem'

subjects = 1;

cd D:\Maastricht_2016-2018\EEG_AVspeech\ANALYSES_EEG_2019\Source_Analyses\headmodel;
% Prepare layout;
cfg = [];
cfg.layout='easycapM15_source';
lay = ft_prepare_layout(cfg,[]);

% Load the polhemus of the participant and add elect.labels;
s = subjects;
pos_file = ('test_cap.pos');
elec = ft_read_sens(pos_file); 
elec.label = elec.label(1:128);

% Get the fidicual locations marked in the MRI and warp them onto the BEM;
vox_nas   = [92 210 35];     % approximate nasion position;
vox_lpa   = [9 112 26];      % approximate lpa position;
vox_rpa   = [175 112 26];    % approximate rpa position;
vox_trans = mri.transform;   % mri transformation matrix;

head_nas  = ft_warp_apply(vox_trans, vox_nas, 'homogenous'); % Nasion warped onto BEM;
head_lpa  = ft_warp_apply(vox_trans, vox_lpa, 'homogenous'); % LPA warped onto BEM;
head_rpa  = ft_warp_apply(vox_trans, vox_rpa, 'homogenous'); % RPA warped onto BEM;

% Align electrode structure to the BEM using the warped MRI fidicual positions;
bem_fid.elecpos = [head_nas; head_lpa; head_rpa];
bem_fid.label   = {'Nasion', 'LPA', 'RPA'};
bem_fid.unit    = 'mm';
elec = ft_convert_units(elec, 'mm');

cfg = [];
cfg.method = 'template';
cfg.target = bem_fid;
cfg.elec = elec;
cfg.fiducial = {'Nasion', 'LPA', 'RPA'};
elec = ft_electroderealign(cfg);
% convert measures in mm;
elec_ft = ft_convert_units(elec, 'mm');

% Create subject folder and save the re-aligned electrodes;
cd cd D:\Maastricht_2016-2018\EEG_AVspeech\ANALYSES_EEG_2019\Source_Analyses\headmodel\;
save(['subj',num2str(s),'_elec_aligned'], 'elec');

%% Create the head model with the re-aligned electrodes;

% Plot the aligned MRI fiducials;
figure;
ft_plot_mesh(bem_fid.elecpos,'vertexcolor',[1 0 0],'vertexmarker','o') % plot fiducials as red circles
ft_plot_mesh(vol.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
% Plot the re-aligned electrodes over the mesh;
ft_plot_mesh(elec_ft.chanpos,'vertexcolor',[1 0 0],'vertexmarker','o') % plot all electrodes as red circles
ft_plot_mesh(vol.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);

% Interactive align the electrodes to the BEM;
cfg = [];
cfg.method = 'interactive';
cfg.elec = elec_ft;
cfg.headshape = vol.bnd(1);
elec_ft = ft_electroderealign(cfg);

% Plot the re-aligned electrodes;
figure;
ft_plot_mesh(elec_ft.chanpos,'vertexcolor',[1 0 0],'vertexmarker','o') % plot all electrodes as red circles
ft_plot_mesh(vol.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);

% Prepare source model using pre-specifed [x,y,z] grid; 
vol = ft_convert_units(vol,'cm');
elec_ft = ft_convert_units(elec_ft,'cm');
cfg = [];
cfg.grid.xgrid = -7 : 1 : 7;
cfg.grid.ygrid = -10 : 1 : 7;
cfg.grid.zgrid = -6 : 1 : 8;    
cfg.headmodel = vol;
cfg.elec = elec_ft;
grid = ft_prepare_sourcemodel(cfg);

% Save all the headmodel variables of the subject;
cd D:\Maastricht_2016-2018\EEG_AVspeech\ANALYSES_EEG_2019\Source_Analyses\headmodel\;
save elec_ft elec_ft;
save grid grid;
save vol vol;

%% Plot the head model of the participant for a very last check;

% Plot only the grid positions within the brain;
figure;
ft_plot_mesh(grid.pos(grid.inside,:))

% Plot the BEM;
ft_plot_mesh(vol.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(vol.bnd(2),'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(vol.bnd(3),'edgecolor','none','facecolor',[0.4 0.6 0.4],'facealpha',0.6);

% Plot aligned electrodes without modifying;

s = 19;

cd D:\maastricht_AV_speech_2020\source_analyses\headmodel;
load(['subj',num2str(s),'_elec_ft']);
load(['subj',num2str(s),'_grid']);
load(['subj',num2str(s),'_vol']);

cfg = [];
cfg.method = 'interactive';
cfg.elec = elec_ft;
cfg.headshape = vol.bnd(1);
ft_electroderealign(cfg);




