%Simulate data - The third node gives common input to the other nodes (nodes 1 and 2) at a time delay of 1 and 2 samples

cfg = [];
cfg.method      = 'ar';
cfg.ntrials     = 1000;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 3;
cfg.bpfilter    = 'no';
cfg.blc         = 'yes';
cfg.params(:,:,1)      = [0.55   0       0.25; 
                          0      0.55    0.25;
                          0      0       0.55]; %this simulates 3->1 and 3->2 influence at the first time delay
cfg.params(:,:,2)      = [-0.8  0        -0.1; 
                           0   -0.8      -0.1;
                           0    0        -0.8];  %this simulates 3->1 and 3->2 influence at the second time delay
cfg.noisecov     = [1 0 0; 
                    0 1 0;
                    0 0 1];

data = ft_connectivitysimulation(cfg);

%Remove the third channel 

data_12 = data;
for i = 1:cfg.ntrials
    data_12.trial{i}(3,:) = [];
end
data_12.label(end)= [];

%running this model (that lacks the explicit common input term) gives
%identical parametric granger and PSI results, but the nonparametric
%granger diverges. Increasing numiterations improves but does not entriely
%remove this.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate power, coherence, and Granger causality based on parametric and
%non-parametric estimates

%calculate the fourier coefficients (non-parametric derivation of power)
cfg        = [];
cfg.method = 'mtmfft';
cfg.taper  = 'dpss'; 
cfg.output = 'fourier';
cfg.tapsmofrq = 3;
cfg.foilim = [0 100];
freq       = ft_freqanalysis(cfg, data);

%freqdescriptives calculates the power spectrum
cfg         = [];
cfg.complex = 'complex';
cfg.jackknife     = 'yes';

fd          = ft_freqdescriptives(cfg, freq);

%Parametric (auto-regressive model based) derivation of AR coefficients
%multivariate analysis will compute the auto-regressive coefficients and associated noise covariance matrix
%JM: can this interface to calling mvaranalysis be cleaned up?

%Warning: to run this code (ft_mvaranalysis) requires having installed the
%BioSig toolbox and include function 'mvar' in the path:
%Download it at http://sourceforge.net/projects/biosig/

ntrl          = length(data.trial);
nsmp          = size(data.trial{1},2);
data.cfg.trl  = [1:nsmp:(ntrl-1)*nsmp+1;nsmp:nsmp:ntrl*nsmp]';
data.cfg.trl(:,3) = 0;
cfg           = [];
cfg.t_ftimwin = 1;
cfg.toi       = 0.5;
cfg.order     = 2;                           %model order of 2, this is known a priori (we simulated the data using a model order of 2)
mdata         = ft_mvaranalysis(cfg, data);

%calculate cross-spectral density and transfer functions associated with the auto-regressive model
cfg        = [];
cfg.method = 'mvar';
cfg.foi    = [0:100];
mfreq      = ft_freqanalysis(cfg, mdata);

%Phase-slope index calculation
cfg           = [];
cfg.method    = 'psi';
cfg.bandwidth = 4;
psi1 = ft_connectivityanalysis(cfg, freq);
%psi2 = ft_connectivityanalysis(cfg, mfreq);

%Coherence calculation
cfg           = [];
cfg.method    = 'coh';
cfg.complex   = 'abs';
coh1 = ft_connectivityanalysis(cfg, freq); 
coh2 = ft_connectivityanalysis(cfg, mfreq); 

%Granger causality calculation
cfg           = [];
cfg.method    = 'granger';
cfg.granger.sfmethod = 'multivariate';
g1 = ft_connectivityanalysis(cfg, freq);
g1 = ft_checkdata(g1, 'cmbrepresentation', 'full');
g2 = ft_connectivityanalysis(cfg, mfreq);
% 
% cfg = [];
% cfg.method = 'instantaneous_causality';
% ic1 = connectivityanalysis(cfg,freq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now plot the output for the various connectivity measures:
%out put variables 1 - nonparam 2 - param

figure;plot(fd.freq, fd.powspctrm(1,:)); hold on
plot(fd.freq, fd.powspctrm(2,:),'r'); 
plot(fd.freq, fd.powspctrm(3,:),'k'); 
legend('Power ch 1', 'Power ch 2','Power ch 3');
title('Nonparametric Power');

figure;plot(mfreq.freq, squeeze(abs(mfreq.crsspctrm(1,1,:)))); hold on; 
plot(mfreq.freq, squeeze(abs(mfreq.crsspctrm(2,2,:))),'r'); 
plot(mfreq.freq, squeeze(abs(mfreq.crsspctrm(3,3,:))),'k'); 
legend('Power ch 1', 'Power ch 2','Power ch 3');
title('Parametric Power');


figure;plot(g1.freq,squeeze(coh1.cohspctrm(1,2,:))); hold on;
plot(g1.freq,squeeze(coh1.cohspctrm(1,3,:)),'r');
plot(g1.freq,squeeze(coh1.cohspctrm(2,3,:)),'k');
legend('1-2','1-3','2-3');
title('Nonparametric Coherence spectrum');

figure;plot(g1.freq,squeeze(coh2.cohspctrm(1,2,:))); hold on
plot(g1.freq,squeeze(coh2.cohspctrm(1,3,:)),'r'); 
plot(g1.freq,squeeze(coh2.cohspctrm(2,3,:)),'k'); 
legend('1-2','1-3','2-3');
title('Parametric Coherence spectrum');


figure;plot(g1.freq,squeeze(psi1.psispctrm(1,2,:))); hold on;
plot(g1.freq,squeeze(psi1.psispctrm(1,3,:)),'r');
plot(g1.freq,squeeze(psi1.psispctrm(2,3,:)),'k');
legend('1->2','1->3','3->2');title('PSI nonparametric'); 

figure;plot(g1.freq,squeeze(g1.grangerspctrm(1,2,:)));hold on
plot(g1.freq,squeeze(g1.grangerspctrm(2,1,:)),'r');
plot(g1.freq,squeeze(g1.grangerspctrm(3,1,:)),'k');
plot(g1.freq,squeeze(g1.grangerspctrm(3,2,:)),'g');
plot(g1.freq,squeeze(g1.grangerspctrm(1,3,:)),'c');
plot(g1.freq,squeeze(g1.grangerspctrm(2,3,:)),'y');
title('Granger nonparametric estimates');legend('1->2','2->1','3->1','3->2','1->3','2->3')


figure; plot(g1.freq,squeeze(g2.grangerspctrm(1,2,:)));hold on;
plot(g1.freq,squeeze(g2.grangerspctrm(2,1,:)),'r');
plot(g1.freq,squeeze(g2.grangerspctrm(3,1,:)),'k');
plot(g1.freq,squeeze(g2.grangerspctrm(3,2,:)),'g');
plot(g1.freq,squeeze(g2.grangerspctrm(1,3,:)),'c');
plot(g1.freq,squeeze(g2.grangerspctrm(2,3,:)),'y');
legend('1->2','2->1','3->1','3->2','1->3','2->3');title('Granger parametric estimates ');


% 
% 
% figure; plot(g1.freq, squeeze(abs(ic1.instantspctrm(1,2,:)))); hold on;
% plot(g1.freq, squeeze(abs(ic1.instantspctrm(1,3,:))),'r');
% plot(g1.freq, squeeze(abs(ic1.instantspctrm(2,3,:))),'k');
% legend('1-2','1-3','2-3');
% title('Instantaneous Causality');
% 












