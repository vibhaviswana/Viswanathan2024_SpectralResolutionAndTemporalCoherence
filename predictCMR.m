% Derive comodulation masking release (CMR) predictions from the 
% Viswanathan et al. (2022) temporal-coherence-based source segregation model 
% as a function of CI vocoding and current spread. 
% Copyright 2024 Vibha Viswanathan. All rights reserved.

%% Setup

% Model parameters
f_low = 125;
f_high = 8000;
nfilts = 30;
list_CFs = invcams(linspace(cams(f_low), cams(f_high), nfilts));
nCFs = numel(list_CFs);
cohc  = 1.0;    % normal OHC function
cihc  = 1.0;    % normal IHC function
% Species: specify 1 for cat, 2 for human with Shera et al. tuning,
% 3 for human with Glasberg & Moore tuning
species = 2;
noiseType = 0;  % 1 for variable fGn, 0 for fixed fGn (constant spont rate)
spont = 10; % AN fiber spontaneous rates for the Bruce 2018 model
% implnt: 0 for approximate or 1 for full implementation of the
% power-law functions in the synapse (neural adaptation)
implnt = 0;
tabs = 0.6e-3; % absolute refractory period
trel = 0.6e-3; % relative refractory period
fs_model = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
indCF_model = 22;
nrep1 = 30; % number of stimulus repetitions (e.g., 50);
nrep2 = 5;

% Stimulus parameters
OFC_SPL = 48;
T = 0.5; % duration of stimulus (s)
SNRs = [12:-6:-18,-inf]; % signal-to-noise ratios
SNRs_for_plotting = 12:-6:-24;
nSNRs = numel(SNRs);
% CMR test stimuli: 
% inputComod_intact_voc_currSpread and inputCodev_intact_voc_currSpread 
% contain stimuli for the Comodulated and Codeviant conditions, respectively. 
% Each variable has size: number of vocoding conditions x number of SNRs x number of time samples.
load('CMRtestStimuli.mat'); 

% PSTH parameters
psthbinwidth = 1e-3; % binwidth (s)
nTimePtsPerPSTHbin = round(psthbinwidth*fs_model);  % no. of time points per PSTH bin
numPSTHbins = round((T*fs_model)/nTimePtsPerPSTHbin); % total no. of PSTH bins
fs_PSTH = 1/psthbinwidth;

% Vocoding params
nVocodingConds = 3; % Intact, Vocoded, Current Spread conditions
allTitleStrs = {'Intact','Vocoded',sprintf('%s\\newline%s','Current','Spread')}; %'Current Spread'};%

%% Run model

psth_comod = zeros(nVocodingConds,nSNRs,nCFs,nrep2,numPSTHbins);
psth_codev = zeros(nVocodingConds,nSNRs,nCFs,nrep2,numPSTHbins);
CNoutput_comod = zeros(nVocodingConds,nSNRs,nrep2,numPSTHbins);
CNoutput_codev = zeros(nVocodingConds,nSNRs,nrep2,numPSTHbins);
for whichVocCond = 1:nVocodingConds
    for whichSNR = 1:nSNRs
        for whichCF = 1:nCFs
            CF = list_CFs(whichCF);
            for whichrep = 1:nrep2
                % COMOD
                vihc = model_IHC_BEZ2018(inputComod_intact_voc_currSpread(whichVocCond,whichSNR,:),CF,nrep1,...
                    1/fs_model,T,cohc,cihc,species);
                [out1,~,~] = model_Synapse_BEZ2018(vihc,...
                    CF,nrep1,1/fs_model,noiseType,implnt,spont,tabs,trel);
                psth_comod(whichVocCond,whichSNR,whichCF,whichrep,:) = sum(reshape(out1,nTimePtsPerPSTHbin,numPSTHbins)); % calculate proper PSTH (unit: spikes/s)
                psth_comod(whichVocCond,whichSNR,whichCF,whichrep,:) = psth_comod(whichVocCond,whichSNR,whichCF,whichrep,:)/nrep1/psthbinwidth;
                % CODEV
                vihc = model_IHC_BEZ2018(inputCodev_intact_voc_currSpread(whichVocCond,whichSNR,:),CF,nrep1,...
                    1/fs_model,T,cohc,cihc,species);
                [out1,~,~] = model_Synapse_BEZ2018(vihc,...
                    CF,nrep1,1/fs_model,noiseType,implnt,spont,tabs,trel);
                psth_codev(whichVocCond,whichSNR,whichCF,whichrep,:) = sum(reshape(out1,nTimePtsPerPSTHbin,numPSTHbins)); % calculate proper PSTH (unit: spikes/s)
                psth_codev(whichVocCond,whichSNR,whichCF,whichrep,:) = psth_codev(whichVocCond,whichSNR,whichCF,whichrep,:)/nrep1/psthbinwidth;
            end
        end
    end
end
for whichVocCond = 1:nVocodingConds
    for whichSNR = 1:nSNRs
        for whichrep = 1:nrep2
            temp = CNunit(squeeze(psth_comod(whichVocCond,whichSNR,:,whichrep,:)),list_CFs,indCF_model,psthbinwidth);
            CNoutput_comod(whichVocCond,whichSNR,whichrep,:) = temp(1:numPSTHbins);
            temp = CNunit(squeeze(psth_codev(whichVocCond,whichSNR,:,whichrep,:)),list_CFs,indCF_model,psthbinwidth);
            CNoutput_codev(whichVocCond,whichSNR,whichrep,:) = temp(1:numPSTHbins);
        end
    end
end

%% Derive CMR predictions

null_ind = nSNRs; % no-signal SNR condition
cond = {'comod', 'codev'};
nconds = numel(cond);
varyDprimeThresh = linspace(0.2,0.8,5);
CMR = nan(numel(varyDprimeThresh),nVocodingConds,nrep2);
dprime = nan(nVocodingConds,nrep2,nconds,nSNRs);
for whichVocCond = 1:nVocodingConds
    for whichrep = 1:nrep2
        m = zeros(nconds, nSNRs);
        m0 = zeros(nconds, nSNRs);
        s0 = zeros(nconds, nSNRs);
        for c = 1:nconds
            temp0 = squeeze(eval(strcat('CNoutput_', cond{c},'(whichVocCond,:,whichrep,:)')));
            m(c, :) = mean(temp0, 2);
            m0(c, :) = mean(temp0(null_ind, :), 2);
            s0(c, :) = std(temp0(null_ind, :), [], 2);
        end
        dprime(whichVocCond,whichrep,:,:) = max((m - m0) ./ s0, 0);
        % Calculate CMR
        x = SNRs((end-1):-1:1)';
        % Define growthfit params
        a = 0.005:0.005:0.2; % slope parameter
        b = 0:0.02:2; % y-value at which the growth curve plateus
        c = 0; % baseline parameter (starting y-value)
        tt = -24:1:12; % x-axis location parameter (SNR at which the dprime is halfway between c and b)
        y = squeeze(dprime(whichVocCond, whichrep, 1, (end-1):-1:1));
        [~, params1] = growthfit(x, y, a, b, c, tt);
        y = squeeze(dprime(whichVocCond, whichrep, 2, (end-1):-1:1));
        [~, params2] = growthfit(x, y, a, b, c, tt);
        for kdprimethresh = 1:numel(varyDprimeThresh)
            dprime_thresh = varyDprimeThresh(kdprimethresh);
            thresh_comod = growthinv(dprime_thresh, params1);
            thresh_codev = growthinv(dprime_thresh, params2);
            CMR(kdprimethresh,whichVocCond,whichrep) = thresh_codev - thresh_comod;
            if (dprime(:) < dprime_thresh)
                CMR(kdprimethresh,whichVocCond,whichrep) = 0;
            end
        end
    end
end
dprime_aveReps = squeeze(mean(dprime,2));
dprime_steReps = squeeze(std(dprime,[],2)/sqrt(nrep2));
CMR(CMR<0) = 0;
CMR_avethresh = squeeze(mean(CMR,1)); % average over dprime thresholds

figure;
tl = tiledlayout(1,3,'TileSpacing','Tight','Padding','Compact');
for whichVocCond = 1:nVocodingConds
    nexttile;
    p1 = plot(SNRs_for_plotting,squeeze(dprime_aveReps(whichVocCond,1,:)),'s-','linew', 2,'markersize',10,'color',[239,138,98]/255);
    hold on;
    errorbar(SNRs_for_plotting,squeeze(dprime_aveReps(whichVocCond,1,:)),squeeze(dprime_steReps(whichVocCond,1,:)),'linew', 2,'color',[239,138,98]/255);
    hold on;
    p2 = plot(SNRs_for_plotting,squeeze(dprime_aveReps(whichVocCond,2,:)),'^-','linew', 2,'markersize',10,'color',[103,169,207]/255);
    hold on;
    errorbar(SNRs_for_plotting,squeeze(dprime_aveReps(whichVocCond,2,:)),squeeze(dprime_steReps(whichVocCond,2,:)),'linew', 2,'color',[103,169,207]/255);
    hold on;
    if (whichVocCond==1)
        legend([p1,p2],{'Comodulated','Codeviant'},'location','northwest');
    end
    xticks(fliplr(SNRs_for_plotting));
    xticklabels({'No Signal','-18','-12','-6','0','6','12'});
    title(allTitleStrs{whichVocCond});
    xlim([-24,12]);
    set(gca, 'FontSize', 20);
end
xlabel(tl,'SNR (dB)', 'FontSize', 24);
ylabel(tl,'d''', 'FontSize', 24);

figure;
bar(1:nVocodingConds,squeeze(mean(CMR_avethresh,2)),'facecolor',0.5*[1,1,1],'edgecolor',0.5*[1,1,1],'linewidth',2);
hold on;
errorbar(1:nVocodingConds,squeeze(mean(CMR_avethresh,2)),squeeze(std(CMR_avethresh,[],2)/sqrt(nrep2)),'k','linestyle','none','linewidth',2);
xticks(1:nVocodingConds);
xlim([0.5,nVocodingConds+0.5]);
xticklabels(allTitleStrs);
set(gca,'FontSize',46);
xtickangle(0);
ylabel('CMR (dB)','FontSize',48);

%% Helper functions

function y = growth(x, params)
a = params(1);
b = params(2);
c = params(3);
t = params(4);
y = b./(1 + exp(-1*a*(x - t))) + c;
end

function y = growthinv(x, params)
% Inverse of the growth function
a = params(1);
b = params(2);
c = params(3);
t = params(4);
y = t - (1/a)*log(b./(x-c) - 1);
y(x > b) = 94;
y(y > 88) = 94;
end

function [val, fitted] = growthfit(x, y, a, b, c, t)
% a -- Slope parameter
% b -- Y-value at which the growth curve plateus
% c -- Baseline paramater (starting y-value)
% t -- X-axis location parameter (this is the SNR at which the dprime is halfway between c and b)

% Brute force optimization by least L1 norm
na = numel(a);
nb = numel(b);
nc = numel(c);
nt = numel(t);

err = inf*ones(na, nb, nc, nt);
for ka = 1:na
    for kb = 1:nb
        for kc = 1:nc
            for kt = 1:nt
                params = [a(ka), b(kb), c(kc), t(kt)];
                ypred = growth(x, params);
                err(ka,kb,kc,kt) = nansum(abs(y - ypred));
            end
        end
    end
end

[val, ind] = min(err(:));
[oa, ob, oc, ot] = ind2sub([na, nb, nc, nt], ind);
fitted = [a(oa), b(ob), c(oc), t(ot)];
end
