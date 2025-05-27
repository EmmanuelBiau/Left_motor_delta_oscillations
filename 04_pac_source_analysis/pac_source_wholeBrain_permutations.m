%% source grandaverage across participants
clearvars;clc; close all; addpath C:\toolbox\fieldtrip\; ft_defaults;

%Load mri template from FT and the grid template from one participant (e.g. grid of headmodel);
addpath xxx\source_analyses\headmodel; load('grid.mat');
load standard_mri;

%Format Whole brain delta-beat pac data within grid template for each participant;
allsubj_pac_source = [];

subjects = [1 2 4:26];
for s = subjects

    %load participant's source data;
    cd(['xxx\subj_',num2str(s)]);
    load(['subj_',num2str(s),'_pac_source_wb'],'pac_wb_beta');
    
    pac_temp = pac_wb_beta;
    pac_temp = rmfield(pac_temp,{'phase','bin_val','phapow','powpow'});
    pac_temp.powspctrm = permute(pac_temp.powspctrm,[2 3 1 4]);
    pac_temp.powspctrm = pac_temp.powspctrm(:,:,1,1);
    pac_temp.dimord = 'rpt_chan_freq_time';
    pac_source = pac_temp; clear pac_temp clear pac_wb_beta;

    for cond = 1:4
        
        cfg = [];
        cfg.keeptrials = 'no';
        cfg.avgovertime = 'yes';
        cfg.avgoverfreq = 'yes';
       
        if cond == 1 
            pac_nma_temp = ft_freqdescriptives(cfg,pac_source);
            elseif cond ==2 
            pac_nms_temp = ft_freqdescriptives(cfg,pac_source);
            elseif cond ==3 
            pac_hma_temp = ft_freqdescriptives(cfg,pac_source);
            elseif cond ==4 
            pac_hms_temp = ft_freqdescriptives(cfg,pac_source);
        end

    end

    cfg = [];
    cfg.keeptrials = 'no';
    cfg.avgovertime = 'yes';
    cfg.avgoverfreq = 'yes';
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,1) | isequal(x.condition,2), pac_source.trialinfo, 'UniformOutput', false)));
    pac_NMCs_temp = ft_freqdescriptives(cfg,pac_source);
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,3)| isequal(x.condition,4), pac_source.trialinfo, 'UniformOutput', false)));
    pac_HMCs_temp = ft_freqdescriptives(cfg,pac_source);
    
    % Using template grid;
    sourceTmpl = [];
    sourceTmpl.inside = grid.inside;
    sourceTmpl.dim = grid.dim;
    sourceTmpl.pos = grid.pos;
    sourceTmpl.unit = grid.unit;
    
    %Reconstruct pow at each grid source NMCs;
    pac_nma = sourceTmpl;
    pac_nma.pow = nan(size(grid.pos,1),1);        
    pac_nma.pow(sourceTmpl.inside)= pac_nma_temp.powspctrm;
    pac_nms = sourceTmpl;
    pac_nms.pow = nan(size(grid.pos,1),1);        
    pac_nms.pow(sourceTmpl.inside)= pac_nms_temp.powspctrm;
    pac_hma = sourceTmpl;
    pac_hma.pow = nan(size(grid.pos,1),1);        
    pac_hma.pow(sourceTmpl.inside)= pac_hma_temp.powspctrm;
    pac_hms = sourceTmpl;
    pac_hms.pow = nan(size(grid.pos,1),1);        
    pac_hms.pow(sourceTmpl.inside)= pac_hms_temp.powspctrm;
    pac_NMCs = sourceTmpl;
    pac_NMCs.pow = nan(size(grid.pos,1),1);        
    pac_NMCs.pow(sourceTmpl.inside)= pac_NMCs_temp.powspctrm;
    pac_HMCs = sourceTmpl;
    pac_HMCs.pow = nan(size(grid.pos,1),1);        
    pac_HMCs.pow(sourceTmpl.inside)= pac_HMCs_temp.powspctrm;
    
    %Put all participants in a common structure;
    allsubj_pac_source.NMA{1,s} = pac_nma;
    allsubj_pac_source.NMS{1,s} = pac_nms;
    allsubj_pac_source.HMA{1,s} = pac_hma;
    allsubj_pac_source.HMS{1,s} = pac_hms;
    allsubj_pac_source.NMCs{1,s} = pac_NMCs;
    allsubj_pac_source.HMCs{1,s} = pac_HMCs;
    
    clear pac_nma pac_nms pac_hma pac_hms pac_NMCs pac_HMCs

end

%------------------------------------------------------------------------------------------------%
clear sourceData_int_NMCs sourceData_int_HMCs stat_source_NMCs stat_source_HMCs

%final subjects to include in the source analyses;
if size(allsubj_pac_source.NMA,2)~= 23
    bad_subj = [3 12 20];
    allsubj_pac_source.NMA(bad_subj) = [];
    allsubj_pac_source.NMS(bad_subj) = [];
    allsubj_pac_source.HMA(bad_subj) = [];
    allsubj_pac_source.HMS(bad_subj) = [];
    allsubj_pac_source.NMCs(bad_subj) = [];
    allsubj_pac_source.HMCs(bad_subj) = [];
end

%stats for NMCs versus HMCs;
cfg = [];
nsubj = numel(allsubj_pac_source.NMA);
cfg.design(1,:) = [1:nsubj 1:nsubj];
cfg.design(2,:) = [ones(1,nsubj)*1 ones(1,nsubj)*2];
cfg.uvar = 1; 
cfg.ivar = 2; 

cfg.dim = allsubj_pac_source.NMA{1}.dim;
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.parameter = 'pow';
cfg.correctm = 'cluster';
cfg.correcttail = 'alpha';
cfg.clusteralpha = 0.05;
cfg.clustertail = cfg.clusteralpha;
cfg.numrandomization = 5000;
cfg.alpha = 0.05;
cfg.tail = tail;
stat_source_NMCvsHMC = ft_sourcestatistics(cfg,allsubj_pac_source.NMCs{:},allsubj_pac_source.HMCs{:});
stat_source_NMCs = ft_sourcestatistics(cfg,allsubj_pac_source.NMS{:},allsubj_pac_source.NMA{:});
stat_source_HMCs = ft_sourcestatistics(cfg,allsubj_pac_source.HMS{:},allsubj_pac_source.HMA{:});

%Load mri template from FT and the grid template from one participant (e.g. grid of headmodel);
cd C:\toolbox\fieldtrip\template\headmodel\; load standard_mri;

%interpolate the parameter 'stat';
cfg = [];
cfg.downsample = 2;      
cfg.parameter = 'stat'; 
sourceData_int_stat_source_NMCvsHMC = ft_sourceinterpolate(cfg,stat_source_NMCvsHMC ,mri);
sourceData_int_stat_source_NMCs = ft_sourceinterpolate(cfg,stat_source_NMCs,mri);
sourceData_int_stat_source_HMCs = ft_sourceinterpolate(cfg,stat_source_HMCs,mri);

%Plot source localizations in NMCs and HMCs contrasts;
close all;
atlas = ft_read_atlas('C:\toolbox\fieldtrip\template\atlas\aal\ROI_MNI_V4.nii');

%Plot source localizations in NMCs and HMCs contrasts;
close all;
atlas = ft_read_atlas('C:\toolbox\fieldtrip\template\atlas\aal\ROI_MNI_V4.nii');
%select tvalue range NMCs contrast;
pvalue_NMCvsHMC = vertcat(stat_source_NMCvsHMC.posclusters(1).prob);
tvalue_threshold_NMCvsHMC = min(stat_source_NMCvsHMC.stat(stat_source_NMCvsHMC.prob == pvalue_NMCvsHMC)); 
tvalue_max_NMCvsHMC = max(stat_source_NMCvsHMC.stat(stat_source_NMCvsHMC.prob == pvalue_NMCvsHMC));

%chose representation to be plotted (surface; ortho; slice);
method = 'ortho';

%Plot NMCs vs HMCs delta-beta PAC Contrast;
close all;
cfg = [];
cfg.method = method;                 
cfg.funparameter = 'stat';        
cfg.funcolorlim = 'maxabs';             
cfg.maskparameter = 'stat';
cfg.opacitylim = [tvalue_threshold_NMCvsHMC-0.5 round(tvalue_max_NMCvsHMC,4)+0.5];
cfg.funcolormap = 'jet';
cfg.atlas = atlas;
cfg.location = 'max';
cfg.crosshair = 'yes';
cfg.markercolor   = [1 1 1];
ft_sourceplot(cfg, sourceData_int_stat_source_NMCvsHMC);
set(gcf, 'name','delta-beta PAC: NMCs vs. HMCs','color','w')

%Determine significant sources from CBP tests;
%Load grid source headmodel and atlas;
cd xxx\source_analyses\headmodel\; 
load('grid.mat');
atlas = ft_read_atlas('C:\toolbox\fieldtrip\template\atlas\aal\ROI_MNI_V4.nii');

%Call ft_sourceinterpolate to interpolate atlas with your grid headmodel;
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
atlas_sourcemodel = ft_sourceinterpolate(cfg, atlas, grid);

%replace NAN with zeros;
atlas_sourcemodel.tissue(isnan(atlas_sourcemodel.tissue)) = 0;
ids = find(atlas_sourcemodel.tissue); %all interpolate regions
id = atlas_sourcemodel.tissue(ids); %all interpolate regions index
ROI_labels = atlas.tissuelabel(id);

%find significant grid from stats;
grid_inside = stat_source_NMCvsHMC.prob(stat_source_NMCvsHMC.inside);
[indx_grid,~] = find(round(stat_source_NMCvsHMC.prob,4) == round(stat_source_NMCvsHMC.negclusters(1).prob,4));
%create Table to store siggnificant grids information;
significant_grids = cell(length(indx_grid(:,1)),3);
significant_grids = array2table(significant_grids, 'VariableNames',{'grid_index','grid_t_value','grid_label'});

%find corresponding grid index in the interpolated atlas with your source headmodel; 
count =0; clear i;
for i = indx_grid'
    count = count+1;
    if sum(ismember(ids,i)) == 1
    [significant_grids.grid_index{count},~] = find(ismember(ids,i));    
    elseif sum(ismember(ids,i)) == 0
    [~,significant_grids.grid_index{count}] = min(round(abs(ids - i)));
    end
end

%find the ROI_label corresponding to the grid index;
clear i;
for i = 1:length(significant_grids.grid_index)
    significant_grids.grid_t_value{i} = stat_source_NMCvsHMC.stat(indx_grid(i));
    significant_grids.grid_label{i} = ROI_labels(significant_grids.grid_index{i});
end

significant_grids = sortrows(significant_grids,2,'descend');

%export interpolated source data to nifti; 
cfg = [];
cfg.parameter = 'stat';
cfg.filename = 'PAC_NMCvsHMCs';
cfg.datatype = 'float';
cfg.filetype = 'nifti';
cfg.vmpversion = 2;
cfg.coordsys = 'spm';
ft_volumewrite(cfg,sourceData_int_stat_source_NMCvsHMC);

%%