%% source interpolation across participants;
clearvars;clc; close all; addpath C:\toolbox\fieldtrip\; ft_defaults;

%Load mri template from FT and the grid template from one participant (e.g. grid of headmodel);
cd xxx\headmodel\; load standard_mri;
cd xxx\source_analyses\headmodel\; load('grid.mat');

subjects = [1 2 4:11 13:19 21:26];

%Create a structure with all the MI_source data together;
clear allsubj_source

for s = subjects

    %load participant's source data;
    cd(['xxx\subj_',num2str(s)]); load(['source_data_subj',num2str(s)]);
    
    %Select the time-window and frequency band of interest;
    cfg = [];
    cfg.latency     = [3 9];
    cfg.avgovertime = 'yes';
    cfg.frequency   = [2 3];
    cfg.avgoverfreq = 'yes';
    NMAC_source = ft_selectdata(cfg,source_NMAC);
    NMSC_source = ft_selectdata(cfg,source_NMSC);
    HMAC_source = ft_selectdata(cfg,source_HMAC);
    HMSC_source = ft_selectdata(cfg,source_HMSC);
    NMSI_source = ft_selectdata(cfg,source_NMSI);
    HMAI_source = ft_selectdata(cfg,source_HMAI);
    HMSI_source = ft_selectdata(cfg,source_HMSI);
    
    if ~ismember(s,[18 26]); NMAI_source = ft_selectdata(cfg,source_NMAI); end
    
    %Put all participants in a common structure;
    allsubj_source.NMAC{1,s} = NMAC_source;
    allsubj_source.NMSC{1,s} = NMSC_source;
    allsubj_source.HMAC{1,s} = HMAC_source;
    allsubj_source.HMSC{1,s} = HMSC_source;
    allsubj_source.NMSI{1,s} = NMSI_source;
    allsubj_source.HMAI{1,s} = HMAI_source;
    allsubj_source.HMSI{1,s} = HMSI_source;
    
    if ~ismember(s,[18 26]); allsubj_source.NMAI{1,s} = NMAI_source; end

end

msgbox('**** DONE ****');

%-------------------------------------------------------------------------%
allsource_NMCs = cell(1, size(allsubj_source.NMAC,2));
allsource_HMCs = cell(1, size(allsubj_source.HMAC,2));

for s = subjects

    diff_NMCs_temp = allsubj_source.NMAC{s};
    diff_NMCs_temp.powspctrm = [];
    diff_NMCs_temp.powspctrm = (allsubj_source.NMAC{s}.powspctrm - allsubj_source.NMSC{s}.powspctrm)./allsubj_source.NMSC{s}.powspctrm;  
    allsource_NMCs{s} = diff_NMCs_temp;
    
    diff_HMCs_temp = allsubj_source.HMAC{s};
    diff_HMCs_temp.powspctrm = [];
    diff_HMCs_temp.powspctrm = (allsubj_source.HMAC{s}.powspctrm - allsubj_source.HMSC{s}.powspctrm)./allsubj_source.HMSC{s}.powspctrm;   
    allsource_HMCs{s} = diff_HMCs_temp;
    
    clear diff_NMCs_temp diff_HMCs_temp
    
end

%Grand average allsource_NMCs allsource_HMCs;
bad_subj = [3 12 20];
allsource_NMCs(bad_subj) = [];
allsource_HMCs(bad_subj) = [];

cfg = [];
cfg.parameter = 'powspctrm';
GA_diff_NMCs = ft_freqgrandaverage(cfg,allsource_NMCs{:});
GA_diff_HMCs = ft_freqgrandaverage(cfg,allsource_HMCs{:});

% Using template grid;
sourceTmpl = [];
sourceTmpl.inside = grid.inside;
sourceTmpl.dim = grid.dim;
sourceTmpl.pos = grid.pos;
sourceTmpl.unit = grid.unit;

% Apply this source 'template' to sound and movie data and mask the data;
GA_source_NMCs_temp = sourceTmpl;
GA_source_NMCs_temp.powspctrm = nan(size(grid.pos,1),1); 
GA_source_NMCs_temp.powspctrm(sourceTmpl.inside)= GA_diff_NMCs.powspctrm;
GA_source_HMCs_temp = sourceTmpl;
GA_source_HMCs_temp.powspctrm = nan(size(grid.pos,1),1); 
GA_source_HMCs_temp.powspctrm(sourceTmpl.inside)= GA_diff_HMCs.powspctrm;

%Interpolate the parameter 'powspctrm';
cfg = [];
cfg.downsample = 2;
cfg.parameter = 'powspctrm';
GA_source_NMCs = ft_sourceinterpolate(cfg, GA_source_NMCs_temp, mri); 
GA_source_HMCs = ft_sourceinterpolate(cfg, GA_source_HMCs_temp, mri);

%% Plot source localizations in movies and sounds conditions;
close all;
atlas = ft_read_atlas('xxx\ROI_MNI_V4.nii');
cfg = [];
cfg.method = 'surface';                 
cfg.funparameter = 'powspctrm';        
cfg.funcolorlim = 'maxabs';             
cfg.maskparameter = 'powspctrm';
cfg.opacitylim =  [0.05 0.1];
cfg.funcolormap = 'jet';
cfg.atlas = atlas;
cfg.location = 'max';
cfg.crosshair = 'yes';
cfg.markercolor   = [1 1 1];
ft_sourceplot(cfg, GA_source_NMCs);
set(gcf, 'name','NMCs_source','color','w')
ft_sourceplot(cfg, GA_source_HMCs);
set(gcf, 'name','HMCs_source','color','w')

%% export interpolated source data to nifti;

cd xxx\;
cfg = [];
cfg.parameter = 'powspctrm';
cfg.filename = 'NMCs_contrast';
cfg.datatype = 'float';
cfg.filetype = 'nifti';
cfg.vmpversion = 2;
cfg.coordsys = 'spm';
ft_volumewrite(cfg,GA_source_NMCs);

%%