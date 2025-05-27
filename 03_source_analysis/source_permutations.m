%% source cluster-based permutation tests;
clearvars;clc; close all; addpath C:\toolbox\fieldtrip\; ft_defaults;

%Load mri template from FT and the grid template from one participant (e.g. grid of headmodel);
cd xxx\headmodel\; load standard_mri;
cd xxx\headmodel\; load('grid.mat');

subjects = [1 2 4:11 13:19 21:26];

%Create a structure with all the MI_source data together;
clear allsubj_source

for s = subjects

    %load participant's source data;
    cd(['xxx\individual_source_data\subj_',num2str(s)]);
    load(['Source_data_subj',num2str(s)]);
    
    %Select the time-window and frequency band of interest;
    cfg.latency     = [3 9];
    cfg.avgovertime = 'yes';
    cfg.frequency   = [2 3];
    cfg.avgoverfreq = 'yes';
    NMAC_source = ft_selectdata(cfg,source_NMAC);
    NMSC_source = ft_selectdata(cfg,source_NMSC);
    HMAC_source = ft_selectdata(cfg,source_HMAC);
    HMSC_source = ft_selectdata(cfg,source_HMSC);
    
    %Using template grid;
    sourceTmpl = [];
    sourceTmpl.inside = grid.inside;
    sourceTmpl.dim = grid.dim;
    sourceTmpl.pos = grid.pos;
    sourceTmpl.unit = grid.unit;
    
    %Reconstruct pow at each grid source NMCs;
    sourceData_NMAC = sourceTmpl;
    sourceData_NMAC.pow = nan(size(grid.pos,1),1);        
    sourceData_NMAC.pow(sourceTmpl.inside)= NMAC_source.powspctrm;
    sourceData_NMSC = sourceTmpl;
    sourceData_NMSC.pow = nan(size(grid.pos,1),1);        
    sourceData_NMSC.pow(sourceTmpl.inside)= NMSC_source.powspctrm;
    
    %Reconstruct pow at each grid source HMCs;
    sourceData_HMAC = sourceTmpl;
    sourceData_HMAC.pow = nan(size(grid.pos,1),1);        
    sourceData_HMAC.pow(sourceTmpl.inside)= HMAC_source.powspctrm;
    sourceData_HMSC = sourceTmpl;
    sourceData_HMSC.pow = nan(size(grid.pos,1),1);        
    sourceData_HMSC.pow(sourceTmpl.inside)= HMSC_source.powspctrm;
    
    %Put all participants in a common structure;
    allsubj_source.NMAC{1,s} = sourceData_NMAC;
    allsubj_source.NMSC{1,s} = sourceData_NMSC;
    allsubj_source.HMAC{1,s} = sourceData_HMAC;
    allsubj_source.HMSC{1,s} = sourceData_HMSC;

end

msg = msgbox('**** DONE ****');

%------------------------------------------------------------------------------------------------%
clear sourceData_int_NMCs sourceData_int_HMCs stat_source_NMCs stat_source_HMCs

%final subjects to include in the source analyses;
bad_subj = [3 12 20];
allsubj_source.NMAC(bad_subj) = [];
allsubj_source.NMSC(bad_subj) = [];
allsubj_source.HMAC(bad_subj) = [];
allsubj_source.HMSC(bad_subj) = [];

%stats for NMCs;
cfg = [];
cfg.avgoverchannel ='yes';
cfg.dim = allsubj_source.NMAC{1}.dim;
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.parameter = 'pow';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.tail = 0;
cfg.numrandomization = 5000;

nsubj = numel(allsubj_source.NMAC);
cfg.design(1,:) = [1:nsubj 1:nsubj];
cfg.design(2,:) = [ones(1,nsubj)*1 ones(1,nsubj)*2];
cfg.uvar = 1; 
cfg.ivar = 2; 

stat_source_NMCs = ft_sourcestatistics(cfg,allsubj_source.NMAC{:},allsubj_source.NMSC{:});
stat_source_HMCs = ft_sourcestatistics(cfg,allsubj_source.HMAC{:},allsubj_source.HMSC{:});

%Load mri template from FT;
cd xxx\headmodel\; load standard_mri;

%interpolate the parameter 'stat';
cfg = [];
cfg.downsample = 2;      
cfg.parameter = 'stat'; 
sourceData_int_NMCs = ft_sourceinterpolate(cfg,stat_source_NMCs,mri);
sourceData_int_HMCs = ft_sourceinterpolate(cfg,stat_source_HMCs,mri);

%% Plot source localizations in NMCs and HMCs contrasts;
close all;
atlas = ft_read_atlas('C:\toolbox\fieldtrip\template\atlas\aal\ROI_MNI_V4.nii');
%select tvalue range NMCs contrast;
pvalue_NMCs = vertcat(stat_source_NMCs.posclusters(1).prob);
tvalue_threshold_NMCs = min(stat_source_NMCs.stat(stat_source_NMCs.prob == pvalue_NMCs)); 
tvalue_max_NMCs = max(stat_source_NMCs.stat(stat_source_NMCs.prob == pvalue_NMCs));
%select tvalue range HMCs contrast;
pvalue_HMCs = vertcat(stat_source_HMCs.posclusters(1).prob);
tvalue_threshold_HMCs = min(stat_source_HMCs.stat(stat_source_HMCs.prob == pvalue_HMCs)); 
tvalue_max_HMCs = max(stat_source_HMCs.stat(stat_source_HMCs.prob == pvalue_HMCs));

%chose representation to be plotted (surface; ortho; slice);

method = 'surface';

%Plot NMCs Contrast;
cfg = [];
cfg.method = method;       
cfg.funparameter = 'stat';        
cfg.funcolorlim = 'maxabs';             
cfg.maskparameter = 'stat';
cfg.opacitymap =  'rampup';
cfg.funcolormap = 'jet';
cfg.opacitylim = [tvalue_threshold_NMCs round(tvalue_max_NMCs,1)];
cfg.location = 'max';
cfg.crosshair = 'yes';
cfg.atlas = atlas;
ft_sourceplot(cfg,sourceData_int_NMCs);
set(gcf, 'name', ['NMCs_diff', num2str(numel(allsubj_source.NMAC))]);
%Plot HMCs Contrast;
cfg = [];
cfg.method = method;        
cfg.funparameter = 'stat';        
cfg.funcolorlim = 'maxabs';             
cfg.maskparameter = 'stat';
cfg.opacitymap =  'rampup';
cfg.funcolormap = 'jet';
cfg.opacitylim =  [tvalue_threshold_HMCs tvalue_max_HMCs];
cfg.location = 'max';
cfg.crosshair = 'yes';
cfg.atlas = atlas; 
ft_sourceplot(cfg,sourceData_int_HMCs);
set(gcf, 'name', ['HMCs_diff', num2str(numel(allsubj_source.NMAC))]);

%% export interpolated source data to nifti; 

cfg = [];
cfg.parameter = 'stat';
cfg.filename = 'NMCs_contrast';
cfg.datatype = 'float';
cfg.filetype = 'nifti';
cfg.vmpversion = 2;
cfg.coordsys = 'spm';
ft_volumewrite(cfg, sourceData_int_NMCs);

cfg = [];
cfg.parameter = 'stat';
cfg.filename = 'HMCs_contrast';
cfg.datatype = 'float';
cfg.filetype = 'nifti';
cfg.vmpversion = 2;
cfg.coordsys = 'spm';
ft_volumewrite(cfg, sourceData_int_HMCs);

%%