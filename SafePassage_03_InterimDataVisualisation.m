%% Safe Passage; 3) Visualisation of data

% This script plots the data to give an visualistation of the EEG data for:
% all trls (conditions collapsed)
% 1) N clean trials per channel
% Topoplots for:
% a) topoplot for theta power; 4-7Hz
% b) topoplot for alpha power; 8-12Hz
% Power spectra for:
% different ROIs for collapsed and separate conditions, and for different
% age bins

% Created by Rianne Haartsen, November 2021

%%
% Set up local paths to scripts
clear variables

%add eeglab path
addpath(genpath('XXX/eeglab14_1_2b')); 
%add fieldtrip path and set to defaults
addpath('XXX/MATLAB/fieldtrip-20180925'); 
ft_defaults
%add LM code TaskEngine2
addpath(genpath('XXX/TaskEngine2'))
addpath(genpath('XXX/lm_tools'));
addpath(genpath('XXX/eegtools'));
%add RH code
addpath('XXX/Safe_Passage');
% add SP path

cd('XXX/SafePassage')
load SP_preprocReREFavg_Power.mat

% get paths
fullpath_data = cell(1,1);
for ss = 1:height(SP_preprocReREFavg_Power)
    if ~isempty(SP_preprocReREFavg_Power.PreprocEEG{ss})
        if isempty(fullpath_data{1,1})
            fullpath_data{1,1} = SP_preprocReREFavg_Power.PreprocEEG{ss};
        else
            fullpath_data{1,end+1} = SP_preprocReREFavg_Power.PreprocEEG{ss};
        end
    end
end

% load data 
pow_data = cellfun(@load, fullpath_data, 'uniform', false);

%% Number of clean trials per channel
NCT_powvals = cellfun(@(x) x.DataPowb_clean.powspctrm(:,:,1), pow_data, 'uniform', false);
NCT_perch = nan(size(NCT_powvals,2),20);
for ss = 1:size(NCT_powvals,2)
    xxx = ~isnan(NCT_powvals{1,ss});
    NCT_chs = sum(xxx,1);
    Ntot = size(NCT_powvals{1,ss},1);
    NCT_perch(ss,:) = [Ntot NCT_chs];
end
clear xxx NCT_chs Ntot

% visualise the data 
NCT_plots = figure;
subplot(1,2,1)
imagesc(NCT_perch); 
a = colorbar;
a.Label.String = 'Number of clean trials';
xticks([1:1:20])
xticklabels(cat(2,{'Total'},pow_data{1,1}.DataPowb_clean.label'))
ylabel('Subject nr'); xlabel('Channel')
title('Number of clean trials per channel')

subplot(1,2,2)
percNCT_perch = NCT_perch(:,2:end)./NCT_perch(:,1)*100;
imagesc(percNCT_perch); 
b = colorbar;
b.Label.String = 'Percentage of clean trials';
xticks([1:1:20])
xticklabels(pow_data{1,1}.DataPowb_clean.label')
ylabel('Subject nr'); xlabel('Channel')
title('Percentage of clean trials per channel')

% topoplots for NCT percentage
NCT_perch_avg = mean(NCT_perch(:,2:end),1);
percNCT_perch_avg = mean(percNCT_perch,1);

Template_ftdata = pow_data{1,1}.DataPowb_clean;
xxx = repmat(percNCT_perch_avg,[size(Template_ftdata.powspctrm,1),1]);
Template_ftdata.percNTC_avg = zeros(size(Template_ftdata.powspctrm));
for ii = 1:size(Template_ftdata.percNTC_avg,3)
    Template_ftdata.percNTC_avg(:,:,ii) = xxx;
end

load SP_20ch_layout_labels.mat

        cfg = [];
        cfg.layout = 'eeg1010.lay';
        cfg.channel = Template_ftdata.label;
        SP_layout1010Enobio = ft_prepare_layout(cfg);
  
figure
    cfg = [];
    cfg.parameter           = 'percNTC_avg'; 
    cfg.highlight           = 'labels';
    cfg.highlightcolor      = [0 0 0];
    cfg.highlightsize       = 14;
    cfg.highlightfontsize   = 14;
    cfg.colorbar            = 'SouthOutside';
    cfg.zlim                = [90 100];
    cfg.comment             = 'zlim'; 
    cfg.layout              = SP_layout1010Enobio; 
    ft_topoplotER(cfg, Template_ftdata)

%% Topoplots for theta and alpha power (log and absolute)

% get log power data for all trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PowLog_alltrls = cellfun(@(x) x.Power_alltrls.log, pow_data, 'uniform', false);
PowLog_alltrls_mat = nan(19,95,size(PowLog_alltrls,2));
for ii = 1:size(PowLog_alltrls,2)
    PowLog_alltrls_mat(:,:,ii) = cell2mat(PowLog_alltrls(1,ii));
end
   
% Calculate grand averages across all participants
Gavg_freq_alltrls = mean(PowLog_alltrls_mat,3);

Template_ftdata = pow_data{1,1}.DataPow_clean;
Template_ftdata.pow_avg = zeros(size(Template_ftdata.powspctrm));
for tt = 1:size(Template_ftdata.pow_avg,3)
    Template_ftdata.pow_avg(tt,:,:) = Gavg_freq_alltrls;
end

load SP_20ch_layout_labels.mat

        cfg = [];
        cfg.layout = 'eeg1010.lay';
        cfg.channel = Template_ftdata.label;
        SP_layout1010Enobio = ft_prepare_layout(cfg);
  
figure
subplot(2,1,1) % theta power
    cfg = [];
    cfg.parameter           = 'pow_avg'; 
    cfg.highlight           = 'labels';
    cfg.highlightcolor      = [0 0 0];
    cfg.highlightsize       = 8;
    cfg.highlightfontsize   = 12;
    cfg.colorbar            = 'WestOutside';
    cfg.xlim                = [4 7];
    cfg.layout              = SP_layout1010Enobio; 
    ft_topoplotER(cfg, Template_ftdata)
    title('Theta power (log) - all trials')
subplot(2,1,2) % alpha power
    cfg = [];
    cfg.parameter           = 'pow_avg'; 
    cfg.highlight           = 'labels';
    cfg.highlightcolor      = [0 0 0];
    cfg.highlightsize       = 8;
    cfg.highlightfontsize   = 12;
    cfg.colorbar            = 'WestOutside';
    cfg.xlim                = [8 12];
    cfg.layout              = SP_layout1010Enobio; 
    ft_topoplotER(cfg, Template_ftdata)
    title('Alpha power (log) - all trials')


    
% get abs power data for all trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PowAbs_alltrls = cellfun(@(x) x.Power_alltrls.abs, pow_data, 'uniform', false);
PowAbs_alltrls_mat = nan(19,95,size(PowAbs_alltrls,2));
for ii = 1:size(PowAbs_alltrls,2)
    PowAbs_alltrls_mat(:,:,ii) = cell2mat(PowAbs_alltrls(1,ii));
end
   
% Calculate grand averages across all participants
Gavg_freq_alltrls = mean(PowAbs_alltrls_mat,3);

Template_ftdata = pow_data{1,1}.DataPow_clean;
Template_ftdata.pow_avg = zeros(size(Template_ftdata.powspctrm));
for tt = 1:size(Template_ftdata.pow_avg,3)
    Template_ftdata.pow_avg(tt,:,:) = Gavg_freq_alltrls;
end

load SP_20ch_layout_labels.mat

        cfg = [];
        cfg.layout = 'eeg1010.lay';
        cfg.channel = Template_ftdata.label;
        SP_layout1010Enobio = ft_prepare_layout(cfg);
  
figure
subplot(2,1,1) % theta power
    cfg = [];
    cfg.parameter           = 'pow_avg'; 
    cfg.highlight           = 'labels';
    cfg.highlightcolor      = [0 0 0];
    cfg.highlightsize       = 8;
    cfg.highlightfontsize   = 12;
    cfg.colorbar            = 'West';
    cfg.xlim                = [4 7];
    cfg.layout              = SP_layout1010Enobio; 
    ft_topoplotER(cfg, Template_ftdata)
    title('Theta power (abs) - all trials')
subplot(2,1,2) % alpha power
    cfg = [];
    cfg.parameter           = 'pow_avg'; 
    cfg.highlight           = 'labels';
    cfg.highlightcolor      = [0 0 0];
    cfg.highlightsize       = 8;
    cfg.highlightfontsize   = 12;
    cfg.colorbar            = 'West';
    cfg.xlim                = [8 12];
    cfg.layout              = SP_layout1010Enobio; 
    ft_topoplotER(cfg, Template_ftdata)
    title('Alpha power (abs)  - all trials')
    
    
    
    
    
%% Separate conditions (log power)
% get log power data for social and non-social trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PowLog_Strls = cellfun(@(x) x.Powerb_soctrls.log, pow_data, 'uniform', false);
PowLog_Strls_mat = nan(19,95,size(PowLog_Strls,2));
PowLog_NStrls = cellfun(@(x) x.Powerb_nsoctrls.log, pow_data, 'uniform', false);
PowLog_NStrls_mat = nan(19,95,size(PowLog_NStrls,2));
PowLog_Cdiff_mat = nan(19,95,size(PowLog_NStrls,2));
for ii = 1:size(PowLog_Strls,2)
    PowLog_Strls_mat(:,:,ii) = cell2mat(PowLog_Strls(1,ii));
    PowLog_NStrls_mat(:,:,ii) = cell2mat(PowLog_NStrls(1,ii));
    PowLog_Cdiff_mat(:,:,ii) = PowLog_Strls_mat(:,:,ii)-PowLog_NStrls_mat(:,:,ii);
end
   
% Calculate grand averages across all participants
% for relative power
Gavg_freq_Strls = mean(PowLog_Strls_mat,3);
Gavg_freq_NStrls = mean(PowLog_NStrls_mat,3);
Gavg_freq_Cdiff = mean(PowLog_Cdiff_mat,3);

% social data
Template_ftdata = pow_data{1,1}.DataPowb_clean;
Template_ftdata.pow_avg = zeros(size(Template_ftdata.powspctrm));
for tt = 1:size(Template_ftdata.pow_avg,3)
    Template_ftdata.pow_avg(tt,:,:) = Gavg_freq_Strls;
end
  
figure
subplot(3,2,1) % theta power
    cfg = [];
    cfg.parameter           = 'pow_avg'; 
    cfg.highlight           = 'labels';
    cfg.highlightcolor      = [0 0 0];
    cfg.highlightsize       = 8;
    cfg.highlightfontsize   = 12;
    cfg.colorbar            = 'WestOutside';
    cfg.xlim                = [4 7];
    cfg.zlim                = [.55 1.25];
    cfg.layout              = SP_layout1010Enobio; 
    ft_topoplotER(cfg, Template_ftdata)
    title('Theta power (log) - social trials')
subplot(3,2,2) % alpha power
    cfg = [];
    cfg.parameter           = 'pow_avg'; 
    cfg.highlight           = 'labels';
    cfg.highlightcolor      = [0 0 0];
    cfg.highlightsize       = 8;
    cfg.highlightfontsize   = 12;
    cfg.colorbar            = 'WestOutside';
    cfg.xlim                = [8 12];
    cfg.zlim                = [-.05 .37];
    cfg.layout              = SP_layout1010Enobio; 
    ft_topoplotER(cfg, Template_ftdata)
    title('Alpha power (log) - social trials')
    
% non-social data
Template_ftdata = pow_data{1,1}.DataPowb_clean;
Template_ftdata.pow_avg = zeros(size(Template_ftdata.powspctrm));
for tt = 1:size(Template_ftdata.pow_avg,3)
    Template_ftdata.pow_avg(tt,:,:) = Gavg_freq_NStrls;
end
  
subplot(3,2,3) % theta power
    cfg = [];
    cfg.parameter           = 'pow_avg'; 
    cfg.highlight           = 'labels';
    cfg.highlightcolor      = [0 0 0];
    cfg.highlightsize       = 8;
    cfg.highlightfontsize   = 12;
    cfg.colorbar            = 'WestOutside';
    cfg.xlim                = [4 7];
    cfg.zlim                = [.55 1.25];
    cfg.layout              = SP_layout1010Enobio; 
    ft_topoplotER(cfg, Template_ftdata)
    title('Theta power (log) - nonsocial trials')
subplot(3,2,4) % alpha power
    cfg = [];
    cfg.parameter           = 'pow_avg'; 
    cfg.highlight           = 'labels';
    cfg.highlightcolor      = [0 0 0];
    cfg.highlightsize       = 8;
    cfg.highlightfontsize   = 12;
    cfg.colorbar            = 'WestOutside';
    cfg.xlim                = [8 12];
    cfg.zlim                = [-.05 .37];
    cfg.layout              = SP_layout1010Enobio; 
    ft_topoplotER(cfg, Template_ftdata)
    title('Alpha power (log) - nonsocial trials')

% condition differences
Template_ftdata = pow_data{1,1}.DataPowb_clean;
Template_ftdata.pow_avg = zeros(size(Template_ftdata.powspctrm));
for tt = 1:size(Template_ftdata.pow_avg,3)
    Template_ftdata.pow_avg(tt,:,:) = Gavg_freq_Cdiff;
end
  
subplot(3,2,5) % theta power
    cfg = [];
    cfg.parameter           = 'pow_avg'; 
    cfg.highlight           = 'labels';
    cfg.highlightcolor      = [0 0 0];
    cfg.highlightsize       = 8;
    cfg.highlightfontsize   = 12;
    cfg.colorbar            = 'WestOutside';
    cfg.xlim                = [4 7];
    cfg.layout              = SP_layout1010Enobio; 
    ft_topoplotER(cfg, Template_ftdata)
    title('Theta power (log) - condition differences')
subplot(3,2,6) % alpha power
    cfg = [];
    cfg.parameter           = 'pow_avg'; 
    cfg.highlight           = 'labels';
    cfg.highlightcolor      = [0 0 0];
    cfg.highlightsize       = 8;
    cfg.highlightfontsize   = 12;
    cfg.colorbar            = 'WestOutside';
    cfg.xlim                = [8 12];
    cfg.layout              = SP_layout1010Enobio; 
    ft_topoplotER(cfg, Template_ftdata)
    title('Alpha power (log) - condition differences')
    
    
    
%% Power spectra for different RoIs

% All trials
Allppt_Front_all = nan(height(SP_preprocReREFavg_Power),95);
Allppt_Cent_all = nan(height(SP_preprocReREFavg_Power),95);
Allppt_Pari_all = nan(height(SP_preprocReREFavg_Power),95);
Allppt_Occip_all = nan(height(SP_preprocReREFavg_Power),95);
Allppt_WS_all = nan(height(SP_preprocReREFavg_Power),95);

for pp = 1:height(SP_preprocReREFavg_Power)
    if ~isempty(SP_preprocReREFavg_Power.ID{pp})
        Allppt_Front_all(pp,:) = SP_preprocReREFavg_Power.PowSp_LogAll_Front{pp};
        Allppt_Cent_all(pp,:) = SP_preprocReREFavg_Power.PowSp_LogAll_Cent{pp};
        Allppt_Pari_all(pp,:) = SP_preprocReREFavg_Power.PowSp_LogAll_Pari{pp};
        Allppt_Occip_all(pp,:) = SP_preprocReREFavg_Power.PowSp_LogAll_Occip{pp};
        Allppt_WS_all(pp,:) = SP_preprocReREFavg_Power.PowSp_LogAll_WholeScalp{pp};
    end
end

load(SP_preprocReREFavg_Power.PreprocEEG{1}, 'DataPowb_clean')
Freqs = DataPowb_clean.freq;
clear DataPowb_clean

figure
Freqrangexlim = [0 20];
s1 = subplot(5,1,1);
    mns = nanmean(Allppt_WS_all,1);
    std_dev = std(Allppt_WS_all,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.3010 0.7450 0.9330],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.3010 0.7450 0.9330],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    xlim(Freqrangexlim)
    xline(4,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(7,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(8,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    xline(12,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    ylabel('Log Power'); xlabel('Frequency (Hz)')
    title('All trials: Across the whole scalp')

s2 = subplot(5,1,2);
    mns = nanmean(Allppt_Front_all,1);
    std_dev = std(Allppt_Front_all,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.3010 0.7450 0.9330],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.3010 0.7450 0.9330],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    xlim(Freqrangexlim)
    xline(4,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(7,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(8,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    xline(12,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    ylabel('Log Power'); xlabel('Frequency (Hz)')
    title('All trials: Frontal (Fz, F3, F4)')

s3 = subplot(5,1,3);
    mns = nanmean(Allppt_Cent_all,1);
    std_dev = std(Allppt_Cent_all,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.3010 0.7450 0.9330],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.3010 0.7450 0.9330],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    xlim(Freqrangexlim)
    xline(4,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(7,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(8,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    xline(12,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    ylabel('Log Power'); xlabel('Frequency (Hz)')
    title('All trials: Central (Cz, C3, C4)')
    
s4 = subplot(5,1,4);
    mns = nanmean(Allppt_Pari_all,1);
    std_dev = std(Allppt_Pari_all,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.3010 0.7450 0.9330],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.3010 0.7450 0.9330],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    xlim(Freqrangexlim)
    xline(4,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(7,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(8,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    xline(12,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    ylabel('Log Power'); xlabel('Frequency (Hz)')
    title('All trials: Parietal (Pz, P3, P4)')

s5 = subplot(5,1,5);
    mns = nanmean(Allppt_Occip_all,1);
    std_dev = std(Allppt_Occip_all,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.3010 0.7450 0.9330],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.3010 0.7450 0.9330],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    xlim(Freqrangexlim)
    xline(4,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(7,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(8,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    xline(12,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    ylabel('Log Power'); xlabel('Frequency (Hz)')
    title('All trials: Occipital (Oz, PO7, PO8)')
    
linkaxes([s1 s2 s3 s4 s5],'xy')

%% Power spectra for separate conditions (social and non-social)

% All trials
Allppt_SFront = nan(height(SP_preprocReREFavg_Power),95);
Allppt_SCent = nan(height(SP_preprocReREFavg_Power),95);
Allppt_SPari = nan(height(SP_preprocReREFavg_Power),95);
Allppt_SOccip = nan(height(SP_preprocReREFavg_Power),95);
Allppt_SWS = nan(height(SP_preprocReREFavg_Power),95);

Allppt_NSFront = nan(height(SP_preprocReREFavg_Power),95);
Allppt_NSCent = nan(height(SP_preprocReREFavg_Power),95);
Allppt_NSPari = nan(height(SP_preprocReREFavg_Power),95);
Allppt_NSOccip = nan(height(SP_preprocReREFavg_Power),95);
Allppt_NSWS = nan(height(SP_preprocReREFavg_Power),95);

for pp = 1:height(SP_preprocReREFavg_Power)
    if ~isempty(SP_preprocReREFavg_Power.ID{pp})
        Allppt_SFront(pp,:) = SP_preprocReREFavg_Power.PowSp_LogSoc_Front{pp};
        Allppt_SCent(pp,:) = SP_preprocReREFavg_Power.PowSp_LogSoc_Cent{pp};
        Allppt_SPari(pp,:) = SP_preprocReREFavg_Power.PowSp_LogSoc_Pari{pp};
        Allppt_SOccip(pp,:) = SP_preprocReREFavg_Power.PowSp_LogSoc_Occip{pp};
        Allppt_SWS(pp,:) = SP_preprocReREFavg_Power.PowSp_LogSoc_WholeScalp{pp};
        
        Allppt_NSFront(pp,:) = SP_preprocReREFavg_Power.PowSp_LogNSoc_Front{pp};
        Allppt_NSCent(pp,:) = SP_preprocReREFavg_Power.PowSp_LogNSoc_Cent{pp};
        Allppt_NSPari(pp,:) = SP_preprocReREFavg_Power.PowSp_LogNSoc_Pari{pp};
        Allppt_NSOccip(pp,:) = SP_preprocReREFavg_Power.PowSp_LogNSoc_Occip{pp};
        Allppt_NSWS(pp,:) = SP_preprocReREFavg_Power.PowSp_LogNSoc_WholeScalp{pp};
    end
end

load(SP_preprocReREFavg_Power.PreprocEEG{1}, 'DataPowb_clean')
Freqs = DataPowb_clean.freq;
clear DataPowb_clean

figure
if ~exist('Freqrangexlim', 'var')
    Freqrangexlim = [0 20];
end
s1 = subplot(5,1,1);
    mns = nanmean(Allppt_SWS,1);
    std_dev = std(Allppt_SWS,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.4940 0.1840 0.5560],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.4940 0.1840 0.5560],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    mns = nanmean(Allppt_NSWS,1);
    std_dev = std(Allppt_NSWS,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.9290 0.6940 0.1250],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.9290 0.6940 0.1250],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    
    xlim(Freqrangexlim)
    xline(4,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(7,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(8,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    xline(12,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    ylabel('Log Power'); xlabel('Frequency (Hz)')
    title('Soc-NSoc trials: Across the whole scalp')

s2 = subplot(5,1,2);
    mns = nanmean(Allppt_SFront,1);
    std_dev = std(Allppt_SFront,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.4940 0.1840 0.5560],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.4940 0.1840 0.5560],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    mns = nanmean(Allppt_NSFront,1);
    std_dev = std(Allppt_NSFront,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.9290 0.6940 0.1250],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.9290 0.6940 0.1250],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    
    xlim(Freqrangexlim)
    xline(4,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(7,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(8,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    xline(12,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    ylabel('Log Power'); xlabel('Frequency (Hz)')
    title('Soc-NSoc trials: Frontal (Fz, F3, F4)')

s3 = subplot(5,1,3);
    mns = nanmean(Allppt_SCent,1);
    std_dev = std(Allppt_SCent,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.4940 0.1840 0.5560],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.4940 0.1840 0.5560],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    mns = nanmean(Allppt_NSCent,1);
    std_dev = std(Allppt_NSCent,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.9290 0.6940 0.1250],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.9290 0.6940 0.1250],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    
    xlim(Freqrangexlim)
    xline(4,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(7,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(8,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    xline(12,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    ylabel('Log Power'); xlabel('Frequency (Hz)')
    title('Soc-NSoc trials: Central (Cz, C3, C4)')
    
s4 = subplot(5,1,4);
    mns = nanmean(Allppt_SPari,1);
    std_dev = std(Allppt_SPari,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.4940 0.1840 0.5560],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.4940 0.1840 0.5560],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    mns = nanmean(Allppt_NSPari,1);
    std_dev = std(Allppt_NSPari,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.9290 0.6940 0.1250],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.9290 0.6940 0.1250],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    
    xlim(Freqrangexlim)
    xline(4,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(7,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(8,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    xline(12,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    ylabel('Log Power'); xlabel('Frequency (Hz)')
    title('Soc-NSoc trials: Parietal (Pz, P3, P4)')

s5 = subplot(5,1,5);
    mns = nanmean(Allppt_SOccip,1);
    std_dev = std(Allppt_SOccip,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.4940 0.1840 0.5560],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.4940 0.1840 0.5560],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    mns = nanmean(Allppt_NSOccip,1);
    std_dev = std(Allppt_NSOccip,'omitnan');
    curve1 = mns + std_dev;
    curve2 = mns - std_dev;
    Freqs2 = [Freqs, fliplr(Freqs)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(Freqs2, inBetween, [0.9290 0.6940 0.1250],'FaceAlpha',0.1,'LineStyle', 'none');
    hold on;
    plot(Freqs, mns, 'LineStyle', '-', 'Color',[0.9290 0.6940 0.1250],'LineWidth',1);
    clear mns std_dev curve1 curve2 Thresholds2 inBetween
    
    xlim(Freqrangexlim)
    xline(4,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(7,'LineWidth',.5,'Color',[0.9290 0.6940 0.1250],'LineStyle','--')
    xline(8,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    xline(12,'LineWidth',.5,'Color',[0.4660 0.6740 0.1880],'LineStyle','--')
    ylabel('Log Power'); xlabel('Frequency (Hz)')
    title('Soc-NSoc trials: Occipital (Oz, PO7, PO8)')
    
linkaxes([s1 s2 s3 s4 s5],'xy')


%% Power spectra in age bins
 
Ages = cell2mat(SP_preprocReREFavg_Power.Age_yrs);
figure; histogram(Ages,[3 5:1:12]); xlabel('Age (years)'); ylabel('Count')
title('Safe Passage study')

% whole scalp
figure

Age3 = SP_preprocReREFavg_Power.PowSp_LogAll_WholeScalp(Ages < 4);
Age4 = SP_preprocReREFavg_Power.PowSp_LogAll_WholeScalp(Ages >= 4 & ...
    Ages < 5);
Age5 = SP_preprocReREFavg_Power.PowSp_LogAll_WholeScalp(Ages >= 5 & ...
    Ages < 6);
Age6 = SP_preprocReREFavg_Power.PowSp_LogAll_WholeScalp(Ages >= 6 & ...
    Ages < 7);
Age7 = SP_preprocReREFavg_Power.PowSp_LogAll_WholeScalp(Ages >= 7 & ...
    Ages < 8);
Age8 = SP_preprocReREFavg_Power.PowSp_LogAll_WholeScalp(Ages >= 8 & ...
    Ages < 9);
Age9 = SP_preprocReREFavg_Power.PowSp_LogAll_WholeScalp(Ages >= 9 & ...
    Ages < 10);
Age10 = SP_preprocReREFavg_Power.PowSp_LogAll_WholeScalp(Ages >= 10 & ...
    Ages < 11);
Age11 = SP_preprocReREFavg_Power.PowSp_LogAll_WholeScalp(Ages >= 11 & ...
    Ages < 12);

ageCol = jet(9);
plot(Freqs, mean(cell2mat(Age3),1),'Color',ageCol(1,:))
hold on
plot(Freqs, mean(cell2mat(Age4),1),'Color',ageCol(2,:))
plot(Freqs, mean(cell2mat(Age5),1),'Color',ageCol(3,:))
plot(Freqs, mean(cell2mat(Age6),1),'Color',ageCol(4,:))
plot(Freqs, mean(cell2mat(Age7),1),'Color',ageCol(5,:))
plot(Freqs, mean(cell2mat(Age8),1),'Color',ageCol(6,:))
plot(Freqs, mean(cell2mat(Age9),1),'Color',ageCol(7,:))
plot(Freqs, mean(cell2mat(Age10),1),'Color',ageCol(8,:))
plot(Freqs, mean(cell2mat(Age11),1),'Color',ageCol(9,:))
legend({strcat('3 yo N=',num2str(size(Age3,1))), strcat('4 yo N=',num2str(size(Age4,1))), ...
    strcat('5 yo N=',num2str(size(Age5,1))), strcat('6 yo N=',num2str(size(Age6,1))), ...
    strcat('7 yo N=',num2str(size(Age7,1))), strcat('8 yo N=',num2str(size(Age8,1))), ...
    strcat('9 yo N=',num2str(size(Age9,1))), strcat('10 yo N=',num2str(size(Age10,1))), ...
    strcat('11 yo N=',num2str(size(Age11,1)))})
xlabel('Frequency (Hz)'); ylabel('Log Power')
title('All trials: whole scalp')


% frontal
figure

Age3 = SP_preprocReREFavg_Power.PowSp_LogAll_Front(Ages < 4);
Age4 = SP_preprocReREFavg_Power.PowSp_LogAll_Front(Ages >= 4 & ...
    Ages < 5);
Age5 = SP_preprocReREFavg_Power.PowSp_LogAll_Front(Ages >= 5 & ...
    Ages < 6);
Age6 = SP_preprocReREFavg_Power.PowSp_LogAll_Front(Ages >= 6 & ...
    Ages < 7);
Age7 = SP_preprocReREFavg_Power.PowSp_LogAll_Front(Ages >= 7 & ...
    Ages < 8);
Age8 = SP_preprocReREFavg_Power.PowSp_LogAll_Front(Ages >= 8 & ...
    Ages < 9);
Age9 = SP_preprocReREFavg_Power.PowSp_LogAll_Front(Ages >= 9 & ...
    Ages < 10);
Age10 = SP_preprocReREFavg_Power.PowSp_LogAll_Front(Ages >= 10 & ...
    Ages < 11);
Age11 = SP_preprocReREFavg_Power.PowSp_LogAll_Front(Ages >= 11 & ...
    Ages < 12);

ageCol = jet(9);
plot(Freqs, mean(cell2mat(Age3),1),'Color',ageCol(1,:))
hold on
plot(Freqs, mean(cell2mat(Age4),1),'Color',ageCol(2,:))
plot(Freqs, mean(cell2mat(Age5),1),'Color',ageCol(3,:))
plot(Freqs, mean(cell2mat(Age6),1),'Color',ageCol(4,:))
plot(Freqs, mean(cell2mat(Age7),1),'Color',ageCol(5,:))
plot(Freqs, mean(cell2mat(Age8),1),'Color',ageCol(6,:))
plot(Freqs, mean(cell2mat(Age9),1),'Color',ageCol(7,:))
plot(Freqs, mean(cell2mat(Age10),1),'Color',ageCol(8,:))
plot(Freqs, mean(cell2mat(Age11),1),'Color',ageCol(9,:))
legend({strcat('3 yo N=',num2str(size(Age3,1))), strcat('4 yo N=',num2str(size(Age4,1))), ...
    strcat('5 yo N=',num2str(size(Age5,1))), strcat('6 yo N=',num2str(size(Age6,1))), ...
    strcat('7 yo N=',num2str(size(Age7,1))), strcat('8 yo N=',num2str(size(Age8,1))), ...
    strcat('9 yo N=',num2str(size(Age9,1))), strcat('10 yo N=',num2str(size(Age10,1))), ...
    strcat('11 yo N=',num2str(size(Age11,1)))})
xlabel('Frequency (Hz)'); ylabel('Log Power')
title('All trials: Frontal region')


% central
figure

Age3 = SP_preprocReREFavg_Power.PowSp_LogAll_Cent(Ages < 4);
Age4 = SP_preprocReREFavg_Power.PowSp_LogAll_Cent(Ages >= 4 & ...
    Ages < 5);
Age5 = SP_preprocReREFavg_Power.PowSp_LogAll_Cent(Ages >= 5 & ...
    Ages < 6);
Age6 = SP_preprocReREFavg_Power.PowSp_LogAll_Cent(Ages >= 6 & ...
    Ages < 7);
Age7 = SP_preprocReREFavg_Power.PowSp_LogAll_Cent(Ages >= 7 & ...
    Ages < 8);
Age8 = SP_preprocReREFavg_Power.PowSp_LogAll_Cent(Ages >= 8 & ...
    Ages < 9);
Age9 = SP_preprocReREFavg_Power.PowSp_LogAll_Cent(Ages >= 9 & ...
    Ages < 10);
Age10 = SP_preprocReREFavg_Power.PowSp_LogAll_Cent(Ages >= 10 & ...
    Ages < 11);
Age11 = SP_preprocReREFavg_Power.PowSp_LogAll_Cent(Ages >= 11 & ...
    Ages < 12);

ageCol = jet(9);
plot(Freqs, mean(cell2mat(Age3),1),'Color',ageCol(1,:))
hold on
plot(Freqs, mean(cell2mat(Age4),1),'Color',ageCol(2,:))
plot(Freqs, mean(cell2mat(Age5),1),'Color',ageCol(3,:))
plot(Freqs, mean(cell2mat(Age6),1),'Color',ageCol(4,:))
plot(Freqs, mean(cell2mat(Age7),1),'Color',ageCol(5,:))
plot(Freqs, mean(cell2mat(Age8),1),'Color',ageCol(6,:))
plot(Freqs, mean(cell2mat(Age9),1),'Color',ageCol(7,:))
plot(Freqs, mean(cell2mat(Age10),1),'Color',ageCol(8,:))
plot(Freqs, mean(cell2mat(Age11),1),'Color',ageCol(9,:))
legend({strcat('3 yo N=',num2str(size(Age3,1))), strcat('4 yo N=',num2str(size(Age4,1))), ...
    strcat('5 yo N=',num2str(size(Age5,1))), strcat('6 yo N=',num2str(size(Age6,1))), ...
    strcat('7 yo N=',num2str(size(Age7,1))), strcat('8 yo N=',num2str(size(Age8,1))), ...
    strcat('9 yo N=',num2str(size(Age9,1))), strcat('10 yo N=',num2str(size(Age10,1))), ...
    strcat('11 yo N=',num2str(size(Age11,1)))})
xlabel('Frequency (Hz)'); ylabel('Log Power')
title('All trials: Central region')


% parietal
figure

Age3 = SP_preprocReREFavg_Power.PowSp_LogAll_Pari(Ages < 4);
Age4 = SP_preprocReREFavg_Power.PowSp_LogAll_Pari(Ages >= 4 & ...
    Ages < 5);
Age5 = SP_preprocReREFavg_Power.PowSp_LogAll_Pari(Ages >= 5 & ...
    Ages < 6);
Age6 = SP_preprocReREFavg_Power.PowSp_LogAll_Pari(Ages >= 6 & ...
    Ages < 7);
Age7 = SP_preprocReREFavg_Power.PowSp_LogAll_Pari(Ages >= 7 & ...
    Ages < 8);
Age8 = SP_preprocReREFavg_Power.PowSp_LogAll_Pari(Ages >= 8 & ...
    Ages < 9);
Age9 = SP_preprocReREFavg_Power.PowSp_LogAll_Pari(Ages >= 9 & ...
    Ages < 10);
Age10 = SP_preprocReREFavg_Power.PowSp_LogAll_Pari(Ages >= 10 & ...
    Ages < 11);
Age11 = SP_preprocReREFavg_Power.PowSp_LogAll_Pari(Ages >= 11 & ...
    Ages < 12);

ageCol = jet(9);
plot(Freqs, mean(cell2mat(Age3),1),'Color',ageCol(1,:))
hold on
plot(Freqs, mean(cell2mat(Age4),1),'Color',ageCol(2,:))
plot(Freqs, mean(cell2mat(Age5),1),'Color',ageCol(3,:))
plot(Freqs, mean(cell2mat(Age6),1),'Color',ageCol(4,:))
plot(Freqs, mean(cell2mat(Age7),1),'Color',ageCol(5,:))
plot(Freqs, mean(cell2mat(Age8),1),'Color',ageCol(6,:))
plot(Freqs, mean(cell2mat(Age9),1),'Color',ageCol(7,:))
plot(Freqs, mean(cell2mat(Age10),1),'Color',ageCol(8,:))
plot(Freqs, mean(cell2mat(Age11),1),'Color',ageCol(9,:))
legend({strcat('3 yo N=',num2str(size(Age3,1))), strcat('4 yo N=',num2str(size(Age4,1))), ...
    strcat('5 yo N=',num2str(size(Age5,1))), strcat('6 yo N=',num2str(size(Age6,1))), ...
    strcat('7 yo N=',num2str(size(Age7,1))), strcat('8 yo N=',num2str(size(Age8,1))), ...
    strcat('9 yo N=',num2str(size(Age9,1))), strcat('10 yo N=',num2str(size(Age10,1))), ...
    strcat('11 yo N=',num2str(size(Age11,1)))})
xlabel('Frequency (Hz)'); ylabel('Log Power')
title('All trials: Parietal region')



% occipital
figure

Age3 = SP_preprocReREFavg_Power.PowSp_LogAll_Occip(Ages < 4);
Age4 = SP_preprocReREFavg_Power.PowSp_LogAll_Occip(Ages >= 4 & ...
    Ages < 5);
Age5 = SP_preprocReREFavg_Power.PowSp_LogAll_Occip(Ages >= 5 & ...
    Ages < 6);
Age6 = SP_preprocReREFavg_Power.PowSp_LogAll_Occip(Ages >= 6 & ...
    Ages < 7);
Age7 = SP_preprocReREFavg_Power.PowSp_LogAll_Occip(Ages >= 7 & ...
    Ages < 8);
Age8 = SP_preprocReREFavg_Power.PowSp_LogAll_Occip(Ages >= 8 & ...
    Ages < 9);
Age9 = SP_preprocReREFavg_Power.PowSp_LogAll_Occip(Ages >= 9 & ...
    Ages < 10);
Age10 = SP_preprocReREFavg_Power.PowSp_LogAll_Occip(Ages >= 10 & ...
    Ages < 11);
Age11 = SP_preprocReREFavg_Power.PowSp_LogAll_Occip(Ages >= 11 & ...
    Ages < 12);

ageCol = jet(9);
plot(Freqs, mean(cell2mat(Age3),1),'Color',ageCol(1,:))
hold on
plot(Freqs, mean(cell2mat(Age4),1),'Color',ageCol(2,:))
plot(Freqs, mean(cell2mat(Age5),1),'Color',ageCol(3,:))
plot(Freqs, mean(cell2mat(Age6),1),'Color',ageCol(4,:))
plot(Freqs, mean(cell2mat(Age7),1),'Color',ageCol(5,:))
plot(Freqs, mean(cell2mat(Age8),1),'Color',ageCol(6,:))
plot(Freqs, mean(cell2mat(Age9),1),'Color',ageCol(7,:))
plot(Freqs, mean(cell2mat(Age10),1),'Color',ageCol(8,:))
plot(Freqs, mean(cell2mat(Age11),1),'Color',ageCol(9,:))
legend({strcat('3 yo N=',num2str(size(Age3,1))), strcat('4 yo N=',num2str(size(Age4,1))), ...
    strcat('5 yo N=',num2str(size(Age5,1))), strcat('6 yo N=',num2str(size(Age6,1))), ...
    strcat('7 yo N=',num2str(size(Age7,1))), strcat('8 yo N=',num2str(size(Age8,1))), ...
    strcat('9 yo N=',num2str(size(Age9,1))), strcat('10 yo N=',num2str(size(Age10,1))), ...
    strcat('11 yo N=',num2str(size(Age11,1)))})
xlabel('Frequency (Hz)'); ylabel('Log Power')
title('All trials: Occipital region')