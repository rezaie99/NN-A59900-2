cfg = [];
cfg.method      = 'ar';
cfg.ntrials     = 1000;
cfg.triallength = 0.5;
cfg.fsample     = 200;
cfg.nsignal     = 2;
cfg.bpfilter    = 'no';
cfg.blc         = 'yes';
cfg.params(:,:,1)      = [0.55   0.025; 
                          0.025  0.55];
cfg.params(:,:,2)      = [-0.8  -0.1; 
                          -0.1  -0.8];
cfg.noisecov     = [1 0.3;  
                    0.3 1];
cfg.absnoise     =         [0 0;   %increase the amount of independent noise added on top of the sender signal.
                            0 0];
data = ft_connectivitysimulation(cfg);

%simulate two signals that both influence each other with equal weight at lag = 1 and 2 (Figure 8A - Case 2: extra noise)

cfg.absnoise = [1 0;   %increase the amount of independent noise added on top of signal 1
                0 0];
datan = ft_connectivitysimulation(cfg);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate power, coherence, and Granger causality based on parametric and
%non-parametric estimates for "clean" data without extra noise (Case 1, Figure 8B-D)

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
cfg.t_ftimwin = 0.5;
cfg.toi       = 0.25;
cfg.order     = 2;
mdata         = ft_mvaranalysis(cfg, data);

mdata.offset = 0;
cfg        = [];
cfg.method = 'mvar';
cfg.foi    = [0:120];
mfreq      = ft_freqanalysis(cfg, mdata);

%calculate the fourier coefficients (non-parametric derivation of power)
cfg        = [];
cfg.method = 'mtmfft';
cfg.taper  = 'dpss';
cfg.output = 'fourier';
cfg.tapsmofrq = 4;
cfg.foilim = [0 100];
freq       = ft_freqanalysis(cfg, data); 

cfg           = [];
cfg.method    = 'coh';
cfg.complex   = 'abs';
coh1 = ft_connectivityanalysis(cfg, freq); 

csd = ft_checkdata(freq, 'cmbrepresentation', 'fullfast');

cfg           = [];
cfg.method = 'coh';
cohp = ft_connectivityanalysis(cfg, mfreq);

cfg           = [];
cfg.method    = 'granger';
cfg.granger.sfmethod = 'bivariate';
g1 = ft_connectivityanalysis(cfg, csd);

cfg           = [];
cfg.method    = 'granger';
gp = ft_connectivityanalysis(cfg, mfreq);


figure;
plot(freq.freq,squeeze(csd.crsspctrm(1,1,:))); hold on; plot(freq.freq,squeeze(csd.crsspctrm(2,2,:)),'r'); 
title('Non-parametric power estimates, , case 1 (no extra noise)'); legend('chan 1','chan 2');

figure; plot(coh1.freq,squeeze(coh1.cohspctrm(1,2,:))); title('Nonparametric Coherence spectrum, case 1 (no extra noise)');
figure; plot(cohp.freq,squeeze(cohp.cohspctrm(1,2,:))); title('Parametric Coherence spectrum, case 1 (no extra noise)');

figure;plot(g1.freq,squeeze(g1.grangerspctrm(1,:)));hold on
plot(g1.freq,squeeze(g1.grangerspctrm(2,:)),'r');
title('Granger nonparametric estimates, case 1 (no extra noise)');legend('1->2','2->1');

figure;plot(gp.freq,squeeze(gp.grangerspctrm(1,2,:)));hold on
plot(gp.freq,squeeze(gp.grangerspctrm(2,1,:)),'r');
title('Granger parametric estimates, case 1 (no extra noise)');legend('1->2','2->1'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate power, coherence, and Granger causality based on parametric and
%non-parametric estimates for data with extra noise (Case 2, Figure 8E-G)

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
cfg.t_ftimwin = 0.5;
cfg.toi       = 0.25;
cfg.order     = 2;
mdata         = ft_mvaranalysis(cfg, datan);

mdata.offset = 0;
cfg        = [];
cfg.method = 'mvar';
cfg.foi    = [0:120];
mfreq      = ft_freqanalysis(cfg, mdata);

%calculate the fourier coefficients (non-parametric derivation of power)
cfg        = [];
cfg.method = 'mtmfft';
cfg.taper  = 'dpss';
cfg.output = 'fourier';
cfg.tapsmofrq = 4;
cfg.foilim = [0 100];
freq       = ft_freqanalysis(cfg, datan); 

cfg           = [];
cfg.method    = 'coh';
cfg.complex   = 'abs';
coh1 = ft_connectivityanalysis(cfg, freq); 

csd = ft_checkdata(freq, 'cmbrepresentation', 'fullfast');

cfg           = [];
cfg.method = 'coh';
cohp = ft_connectivityanalysis(cfg, mfreq);

cfg           = [];
cfg.method    = 'granger';
cfg.granger.sfmethod = 'bivariate';
g1 = ft_connectivityanalysis(cfg, csd);

cfg           = [];
cfg.method    = 'granger';
gp = ft_connectivityanalysis(cfg, mfreq);


figure;
plot(freq.freq,squeeze(csd.crsspctrm(1,1,:))); hold on; plot(freq.freq,squeeze(csd.crsspctrm(2,2,:)),'r'); 
title('Non-parametric power estimates, case 2 (extra noise on channel 1)'); legend('chan 1','chan 2');

figure; plot(coh1.freq,squeeze(coh1.cohspctrm(1,2,:))); title('Nonparametric Coherence spectrum, case 2 (extra noise on channel 1)');
figure; plot(cohp.freq,squeeze(cohp.cohspctrm(1,2,:))); title('Parametric Coherence spectrum, case 2 (extra noise on channel 1)');

figure;plot(g1.freq,squeeze(g1.grangerspctrm(1,:)));hold on
plot(g1.freq,squeeze(g1.grangerspctrm(2,:)),'r');
title('Granger nonparametric estimates, case 2 (extra noise on channel 1)');legend('1->2','2->1');

figure;plot(gp.freq,squeeze(gp.grangerspctrm(1,2,:)));hold on
plot(gp.freq,squeeze(gp.grangerspctrm(2,1,:)),'r');
title('Granger parametric estimates, case 2 (extra noise on channel 1)');legend('1->2','2->1'); 














