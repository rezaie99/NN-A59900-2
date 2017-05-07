
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Common reference, no bipolar derivation, no neuronal coupling (Figure 5A, red trace)
cfg = [];
cfg.method      = 'linear_mix';               %apply a linear mixture of various components into data
cfg.ntrials     = 1000;                       %simulate 1000 trials
cfg.triallength = 1;                          %each trial is 1 second
cfg.fsample     = 1000;                       %sampling rate is 1000 Hz
cfg.nsignal     = 2;                          %simulate 2 unipolar channels
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [35 55 ];                   %apply a band-pass filter between 35 and 55 Hz to the white noise
cfg.blc         = 'yes';                      %remove the mean of each trial

%Mixing matrix corresponds to:
%Source 1, Signal 2, Common reference, Coupling A, Coupling B

%Linearly mix a signal (third column, the common reference) which will be common to both channels
%Columns 1 and 2 correspond to the un-shared signals between the two sources
%Columns 3 and 4 are intentionally left blank, denoting that the shared component (true coupling) is zero
cfg.mix         = [0   1   1 0     0;
                   1   0   1 0     0];
%We can include delays between the different components if we wish. Here they are left to zero, because we assume an instantaneous mixture model            
cfg.delay       = [0 0 0 0 0    
                   0 0 0 0 0] ; 

cfg.absnoise    = 1;   %the variance that will be added to each channel after linear mixing
data1nc = ft_connectivitysimulation(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Common reference, no bipolar derivation, with neuronal coupling (Figure 5A, blue trace)
%same base parameters as above
cfg = [];
cfg.method      = 'linear_mix';
cfg.ntrials     = 1000;
cfg.triallength = 1;
cfg.fsample     = 1000;
cfg.nsignal     = 2;
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [35 55 ];
cfg.blc         = 'yes';

%Mixing matrix columns corresponds to:
%Source 1, Signal 2, Common reference, Coupling A, Coupling B
%Each row is a simulated channel

%Linearly mix a signal (third column, the common reference) which will be common to both channels
%Columns 1 and 2 correspond to the un-shared signals between the two sources
%Columns 3 and 4 correspond to the shared components (true neuronal coupling) between the sources
cfg.mix         = [0     1   1 4     1;
                   1     0   1 2     3];

cfg.delay       = [0 0 0 0 0    
                   0 0 0 0 0] ; 

cfg.absnoise    = 1;
data1c = ft_connectivitysimulation(cfg);



%%%Calculate coherence for unipolar case - no coupling and with coupling

cfg        = [];
cfg.method = 'mtmfft';
cfg.taper  = 'dpss';
cfg.output = 'fourier';
cfg.tapsmofrq = 3;
cfg.foilim = [0 100];
freq1       = ft_freqanalysis(cfg, data1nc); 
cfg           = [];
cfg.method    = 'coh';
cfg.complex   = 'abs';
coh1 = ft_connectivityanalysis(cfg, freq1); 

cfg        = [];
cfg.method = 'mtmfft';
cfg.taper  = 'dpss'; %'rectwin'; %
cfg.output = 'fourier';
cfg.tapsmofrq = 3;
cfg.foilim = [0 100];
freq2       = ft_freqanalysis(cfg, data1c); 
cfg           = [];
cfg.method    = 'coh';
cfg.complex   = 'abs';
coh2 = ft_connectivityanalysis(cfg, freq2); 

figure; plot(coh1.freq,squeeze(coh1.cohspctrm(1,2,:)),'r'); hold on; plot(coh1.freq,squeeze(coh2.cohspctrm(1,2,:))); 
legend('No Coupling', 'Coupling');
title('Unipolar with common reference');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BIPOLAR CASE

%Common reference, with bipolar derivation, no true neuronal coupling (Figure 5B, red trace)
cfg = [];
cfg.method      = 'linear_mix';
cfg.ntrials     = 1000;
cfg.triallength = 1;
cfg.fsample     = 1000;
cfg.nsignal     = 4;                    %simulate 4 unipolar channels
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [35 55];
cfg.blc         = 'yes';

%Mixing matrix columns corresponds to:
%Source 1, Signal 2, Common reference, Source 3, Source 4
%Each row is a simulated channel

%Linearly mix a signal (third column, the common reference) which will be common to all 4 unipolar channels
%Columns 1 and 2 correspond to the un-shared signals between the sources 1 and 2
%Columns 3 and 4 correspond to the un-shared signals between the sources 3 and 4
cfg.mix         = [1   0   1 0     0;
                   0   1   1 0     0;
                   0   0   1 1     0;
                   0   0   1 0     1;
                   
                   ];
              
cfg.delay       = [0 0 0 0 0 ;   
                   0 0 0 0 0 ;
                   0 0 0 0 0 ;   
                   0 0 0 0 0 ;
                   ] ; 
                         
cfg.absnoise    = 1;

data4 = ft_connectivitysimulation(cfg);

%Calculate bipolar derivatives (channel 1 - channel 2 and channel 3 - channel 4)

databp1 = data1nc;

for tr = 1:length(databp1.trial)
   databp1.trial{tr}(1,:) = data4.trial{tr}(1,:) -  data4.trial{tr}(2,:);
   databp1.trial{tr}(2,:) = data4.trial{tr}(3,:) -  data4.trial{tr}(4,:);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Common reference, with bipolar derivation, with true neuronal coupling (Figure 5B, blue trace)
cfg = [];
cfg.method      = 'linear_mix';
cfg.ntrials     = 1000;
cfg.triallength = 1;
cfg.fsample     = 1000;
cfg.nsignal     = 4;
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [35 55];
cfg.blc         = 'yes';


%Mixing matrix columns corresponds to:
%Source 1, Signal 2, Common reference, Source 3, Source 4, Coupling between sources
%Each row is a simulated channel

%Linearly mix a signal (third column, the common reference) which will be common to all 4 unipolar channels
%Columns 1 and 2 correspond to the un-shared signals between the sources 1 and 2
%Columns 3 and 4 correspond to the un-shared signals between the sources 3 and 4
%Columns 5 and 6 correspond to the true neuronal coupling between the sources

cfg.mix         = [1 0   1 0    0 4
                   0 1   1 0    0 2
                   0 0   1 1    0 3
                   0 0   1 0    1 1
                   
                   ];
              
cfg.delay       = [0 0 0 0 0 0;   
                   0 0 0 0 0 0;
                   0 0 0 0 0 0;   
                   0 0 0 0 0 0;
                   
                   ] ; 
                         
cfg.absnoise    = 1;

data4 = ft_connectivitysimulation(cfg);

%Calculate bipolar derivatives (channel 1 - channel 2 and channel 3 - channel 4)

databp2 = data2nc;
for tr = 1:length(databp2.trial)
   databp2.trial{tr}(1,:) = data4.trial{tr}(1,:) -  data4.trial{tr}(2,:);
   databp2.trial{tr}(2,:) = data4.trial{tr}(3,:) -  data4.trial{tr}(4,:);
    
end

%%%Calculate coherence for bipolar case - no coupling and with coupling

cfg        = [];
cfg.method = 'mtmfft';
cfg.taper  = 'dpss'; %'rectwin'; %
cfg.output = 'fourier';
cfg.tapsmofrq = 3;
cfg.foilim = [0 100];
freq1       = ft_freqanalysis(cfg, databp1); 
cfg           = [];
cfg.method    = 'coh';
cfg.complex   = 'abs';
coh1 = ft_connectivityanalysis(cfg, freq1); 

cfg        = [];
cfg.method = 'mtmfft';
cfg.taper  = 'dpss'; %'rectwin'; %
cfg.output = 'fourier';
cfg.tapsmofrq = 3;
cfg.foilim = [0 100];
freq2       = ft_freqanalysis(cfg, databp2); 
cfg           = [];
cfg.method    = 'coh';
cfg.complex   = 'abs';
coh2 = ft_connectivityanalysis(cfg, freq2); 

figure; plot(coh1.freq,squeeze(coh1.cohspctrm(1,2,:)),'r'); hold on; plot(coh1.freq,squeeze(coh2.cohspctrm(1,2,:))); 
legend('No Coupling', 'Coupling');
title('Bipolar recordings with common reference');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%UNIPOLAR case, with seperate references
%Unique references, no bipolar derivation, no neuronal coupling (Figure 5C, red trace)
cfg = [];
cfg.method      = 'linear_mix';               %apply a linear mixture of various components into data
cfg.ntrials     = 1000;                       %simulate 1000 trials
cfg.triallength = 1;                          %each trial is 1 second
cfg.fsample     = 1000;                       %sampling rate is 1000 Hz
cfg.nsignal     = 2;                          %simulate 2 unipolar channels
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [35 55 ];                   %apply a band-pass filter between 35 and 55 Hz to the white noise
cfg.blc         = 'yes';                      %remove the mean of each trial

%Mixing matrix corresponds to:
%Signal A, Signal B, Common part of both signals, Coupling A, Coupling B

%Columns 1 and 2 correspond to the un-shared signals between the two sources
%Column 3 - no common reference - third column is zero 
%Columns 4 and 5 are intentionally left blank, denoting that the shared component (true coupling) is zero
cfg.mix         = [0   1   0 0     0;
                   1   0   0 0     0];
%We can include delays between the different components if we wish. Here they are left to zero, because we assume an instantaneous mixture model            
cfg.delay       = [0 0 0 0 0    
                   0 0 0 0 0] ; 

cfg.absnoise    = 1;   %the variance that will be added to each channel after linear mixing
data2nc = ft_connectivitysimulation(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unique references, no bipolar derivation, with neuronal coupling (Figure 5C, blue trace)
%same base parameters as above
cfg = [];
cfg.method      = 'linear_mix';
cfg.ntrials     = 1000;
cfg.triallength = 1;
cfg.fsample     = 1000;
cfg.nsignal     = 2;
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [35 55 ];
cfg.blc         = 'yes';

%Mixing matrix corresponds to:
%Signal A, Signal B, Common part of both signals, Coupling A, Coupling B

%Columns 1 and 2 correspond to the un-shared signals between the two sources
%Column 3 - no common reference - third column is zero 
%Columns 4 and 5 denote that the shared component (true coupling) between sources

cfg.mix         = [0     1   0 4     1;
                   1     0   0 2     3];

cfg.delay       = [0 0 0 0 0    
                   0 0 0 0 0] ; 

cfg.absnoise    = 1;
data2c = ft_connectivitysimulation(cfg);



%%%Calculate coherence for case of unipolar with unique references case - no coupling and with coupling

cfg        = [];
cfg.method = 'mtmfft';
cfg.taper  = 'dpss';
cfg.output = 'fourier';
cfg.tapsmofrq = 3;
cfg.foilim = [0 100];
freq1       = ft_freqanalysis(cfg, data2nc); 
cfg           = [];
cfg.method    = 'coh';
cfg.complex   = 'abs';
coh1 = ft_connectivityanalysis(cfg, freq1); 

cfg        = [];
cfg.method = 'mtmfft';
cfg.taper  = 'dpss'; %'rectwin'; %
cfg.output = 'fourier';
cfg.tapsmofrq = 3;
cfg.foilim = [0 100];
freq2       = ft_freqanalysis(cfg, data2c); 
cfg           = [];
cfg.method    = 'coh';
cfg.complex   = 'abs';
coh2 = ft_connectivityanalysis(cfg, freq2); 

figure; plot(coh1.freq,squeeze(coh1.cohspctrm(1,2,:)),'r'); hold on; plot(coh1.freq,squeeze(coh2.cohspctrm(1,2,:))); 
legend('No Coupling', 'Coupling');
title('Unipolar recordings with seperate reference');











