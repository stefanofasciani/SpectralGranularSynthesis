% Spectral Granular Synthesis (C) 2018 Stefano Fasciani, University of Wollongong in Dubai
% Inquiries: stefanofasciani@stefanofasciani.com
% 
% The Spectral Granular Synthesis software can be obtained at
% http://stefanofasciani.com/sgs.html
% 
% This Spectral Granular Synthesis is free software: 
% you can redistribute it and/or modify it under the terms of the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
% 
% If you use the Spectral Granular Synthesis or any part of it in any system
% or publication, please acknowledge its author by adding a reference 
% to this pubblication:
% 
% S. Fasciani, "Spectral Granular Synthesis" in proceedings of
% International Computer Music Conference 2018, Daegu, Korea.



% when the GUI is used, this script is automatically run, synthesis parameters
% and options are taken from the GUI. Otherwise, synthesis parameters and 
% can be fedined below if the script is run without the GUI.

%synthesis parameters and options
if ~exist('app','var')
    Fs_sys = 48000; %System sampling rate Hz
    fileIn = ''; % input sound file
    fileOut = ''; %'out.wav'; % output sound file, leave empty if for no output file generation
    normalize = 1; %if 1 output is normalized to max amplitude 1
    playback = 0; % if equal to 1 playback the sound at the end of the synthesis loop
    plots = 0; % if equal to 1 display signal plots and window shape
    duration_sec = 5; %length of output signal in seconds (any value above 0s)
    outChannels = 2; %1 is mono, 2 is stereo, etc...
    compTimeGran = 0; %if 1 it computes the time domain granulation equivalent as well
    grainSizeIn_smp = 2048; %size of grain in samples extracted from source (any value above 0 and below the lenght of the sound file in samples), better if even
    grainSizeOut_smp = 2048; %size of grain in samples sequenced for synthesis (any value above 0 and below the lenght of the sound file in samples), must be even, better if power of 2
    grainOverlap_smp = 1024; %overlap between consecutive grains in samples (any value between 0 and grainSizeIn_smp-1), better if even

    windowType = 3; %defines window type, 0 hanning, 1 hamming, 2 blackman, 3 kaiser
    kaiserBeta = 7; %defines the window shape of kaiser window only, with 0 is retangular
    gammaFast = 0; % 0 computes more accurate gamma for RTPGHI slow, 1 computes less accurate gamma for RTPGHI but fast (computation once only, out of synthesis loop)

    inversionAlgo = 1; %algorithm for spectrum inversion, 0 SPSI, 1 RTPGHI
    thresholdRtphi = 6; %Threshold of RTPGHI 10^-(thresholdRtphi), default is 6

    grainSelMode = 1; % 1 = fixed position + random spray; 2 = sequential + random spray; 3 = random;
    filePos_perc = 0.1; % file position - percentage between 0 and 1 - to extract grains (mode 1 and 2 only)
    spray_perc = 0.3; % random deviation in percentage between 0 and 1 from the fixed position to extract grains (mode 1 and 2 only)
    sprayMode = 1; % 1 backward and forward, 2 backward only, 3 forward only (mode 1 and 2 only)
    sequenceRate = 0.5; % grain selection rate (positive or negative) with respect grain size (mode 2 only)
    filePosModAmp = 0;% amplitude of file position modulation - percentage between 0 and 1
    filePosModFreq = 0; % frequency of file position modulation 

    grainSizeInVarMode = 1; % 0 no modulation, 1 random, 2 sinusoidal modulated 
    grainSizeVarAmp = 0; % percentage of random deviation of the grain size (percentage of grainSizeIn_smp)
    grainSizeVarFreq = 0; % frequency of modulation (only for sinusoidal modulated) 

    freqShiftVarMode = 0; % 0 fixed , 1 random, 2 sinusoidal modulated - variation of shift of FFT bins 0 to 1,
    freqShiftVarAmp = 0; % percentage of random deviation of the grain size -1 fo 1 (percentage of grainSizeIn_smp)
    freqShiftVarFreq = 0; % frequency of modulation (only for sinusoidal modulated) 
else
    Fs_sys = app.Fs_sys;
    fileIn = app.fileIn;
    normalize = app.normalize;
    fileOut = app.fileOut;
    playback = app.playback;
    plots = app.plots;
    duration_sec = app.duration_sec;
    outChannels = app.outChannels;
    compTimeGran = app.compTimeGran;
    grainSizeIn_smp = app.grainSizeIn_smp;
    grainSizeOut_smp = app.grainSizeOut_smp;
    grainOverlap_smp = app.grainOverlap_smp;
    
    windowType = app.windowType;
    kaiserBeta = app.kaiserBeta;
    gammaFast = app.gammaFast;
    
    inversionAlgo = app.inversionAlgo;
    thresholdRtphi = app.thresholdRtphi;
    
    grainSelMode = app.grainSelMode;
    filePos_perc = app.filePos_perc;
    spray_perc = app.spray_perc;
    sprayMode = app.sprayMode;
    sequenceRate = app.sequenceRate;
    filePosModAmp = app.filePosModAmp;
    filePosModFreq = app.filePosModFreq;
    
    grainSizeInVarMode = app.grainSizeInVarMode;
    grainSizeVarAmp = app.grainSizeVarAmp;
    grainSizeVarFreq = app.grainSizeVarFreq;
    
    freqShiftVarMode = app.freqShiftVarMode;
    freqShiftVarAmp = app.freqShiftVarAmp;
    freqShiftVarFreq = app.freqShiftVarFreq;
end

%parameter check, conversion and sound loading

%window handlers
if windowType == 0
    windowhandler = @hanning; iskaiser = 0;
elseif windowType == 1
    windowhandler = @hamming; iskaiser = 0;
elseif windowType == 2
    windowhandler = @blackman; iskaiser = 0;
elseif windowType == 3
    windowhandler = @kaiser; iskaiser = 1;
else
    error('wrong window type');
end


%make grain sizes even
if (rem(grainSizeIn_smp,2)~=0)
    grainSizeIn_smp = grainSizeIn_smp + 1;
end
if (rem(grainSizeOut_smp,2)~=0)
    grainSizeOut_smp = grainSizeOut_smp + 1;
end

spectrumSize = floor(grainSizeOut_smp/2) + 1;

%load rource file and resample
[snd,Fs_sig]=audioread(fileIn);

if Fs_sig ~= Fs_sys
    snd = resample(snd, Fs_sys, Fs_sig);
end 

sndLength_smp = length(snd); %length of input signal in samples
inChannels = size(snd,2); %number of channels input file

if duration_sec < 0
    error('negative duration');
end

%check input parameters
if (grainSizeIn_smp <= 0) || (grainSizeIn_smp > sndLength_smp)
    error('grain size out of range');
end

if (sum(grainSelMode == [1,2,3])==0)
    error('not supported grain selection mode');
end

if (filePos_perc < 0) || (filePos_perc > 1)
    error('file position out of range');
end

if (spray_perc < 0) || (spray_perc > 1)
    error('spray percentage out of range');
end

duration_smp = round(duration_sec * Fs_sys); %length of output signal in samples
grainHop_smp = grainSizeOut_smp - grainOverlap_smp;
filePos_smp = round(filePos_perc * (sndLength_smp - grainSizeIn_smp - 1)); %file position in samples for mode 1

if ~((grainHop_smp >= 1) && (grainHop_smp <= grainSizeOut_smp))
    error('invalid overlap size');
end

%determine variables used for modulations and random variations
filePosMod_smp = filePos_smp;
filePosModFreq_f=(grainHop_smp/Fs_sys) * filePosModFreq;

grainSizeVarAmp_smp = round(grainSizeIn_smp * grainSizeVarAmp);
grainSizeVarFreq_f = (grainHop_smp/Fs_sys) * grainSizeVarFreq;
grainSizeMod_smp = grainSizeIn_smp;

freqShiftVarAmp_bin = round(spectrumSize * freqShiftVarAmp);
freqShiftVarFreq_f = (grainHop_smp/Fs_sys) * freqShiftVarFreq ;
 
% computing constants to support the selected grain selection mode
if grainSelMode == 1 % fixed + random spray + modulation
    
    if sprayMode == 1 % 1 backward and forward, 2 backward only, 3 forward only
        if (filePos_smp > ((sndLength_smp - grainSizeIn_smp - grainSizeVarAmp_smp - 1) - filePos_smp))
            spray_smp = 2 * round(((sndLength_smp - grainSizeIn_smp - grainSizeVarAmp_smp - 1) - filePos_smp) * spray_perc);
            modulAmp_smp = round(((sndLength_smp - grainSizeIn_smp - grainSizeVarAmp_smp - 1) - filePos_smp) * filePosModAmp);
            sprayOffset_smp = spray_smp / 2;
            modulAmpOffset_smp = 0;
            modulAmpOffset_p = 0;
        else
            spray_smp = 2 * round(filePos_smp * spray_perc);
            modulAmp_smp = round(filePos_smp * filePosModAmp);
            sprayOffset_smp = spray_smp / 2;
            modulAmpOffset_smp = 0;
            modulAmpOffset_p = 0;
        end
    elseif sprayMode == 2
        spray_smp = - round(filePos_smp * spray_perc);
        modulAmp_smp = -0.5 * round(filePos_smp * filePosModAmp);
        sprayOffset_smp = 0;
        modulAmpOffset_smp = round(modulAmp_smp/2);
        modulAmpOffset_p = pi/6;
    elseif sprayMode == 3
        spray_smp = round(((sndLength_smp - grainSizeIn_smp - grainSizeVarAmp_smp - 1) - filePos_smp) * spray_perc);
        modulAmp_smp = 0.5 * round(((sndLength_smp - grainSizeIn_smp - grainSizeVarAmp_smp - 1) - filePos_smp) * filePosModAmp);
        sprayOffset_smp = 0;
        modulAmpOffset_smp = round(modulAmp_smp/2);
        modulAmpOffset_p = - pi/6;
    else
        error('wrong spray mode');
    end
    
elseif grainSelMode == 2 % sequential + random spray + modulation
        
    if sprayMode == 1 % 1 backward and forward, 2 backward only, 3 forward only
        spray_smp = 2 * round(abs(round(grainSizeIn_smp * (sequenceRate))) * spray_perc);
        modulAmp_smp = (round(grainSizeIn_smp * filePosModAmp));
        sprayOffset_smp = spray_smp / 2;
        modulAmpOffset_smp = 0;
        modulAmpOffset_p = 0;
    elseif sprayMode == 2 || sprayMode == 3 
        spray_smp = abs(round(grainSizeIn_smp * (sequenceRate))) * spray_perc * sign(sequenceRate);
        modulAmp_smp = 0.5 * (round(grainSizeIn_smp * filePosModAmp) * sign(sequenceRate+eps));
        sprayOffset_smp = 0;
        modulAmpOffset_smp = 0;
        modulAmpOffset_p = 0;
    else
        error('wrong spray mode');
    end
    
end

% generating window envelopes
windowIn = computewinfromhandler(windowhandler,kaiserBeta,grainSizeIn_smp,iskaiser);
windowOut = computewinfromhandler(windowhandler,kaiserBeta,grainSizeOut_smp,iskaiser);
windowOutDual = fftshift(gabdualwindow(fftshift(windowOut),grainHop_smp,grainSizeOut_smp));

%plotting
if plots == 1
    dispWin(:,1) = [windowOut ; zeros(grainHop_smp,1)]; 
    dispWin(:,2) = [zeros(grainHop_smp,1) ; windowOut];
    figure(1);
    plot(dispWin);
    title('Window Shape and Overlap');
    display_waveforms(snd,2,'Source Sound');
    display_spectrograms(snd,3,'Source Sound');
end

%compute constants and variables to support the selected inversion algorithm
if inversionAlgo == 0
    startPhase = zeros(spectrumSize,1);
elseif inversionAlgo == 1
    tol = 10^(-thresholdRtphi);
    mRange = (1:spectrumSize-2)';
    tGradPlusConsts = 2*pi*grainHop_smp*mRange/grainSizeOut_smp;
    newPhase = zeros(spectrumSize,2);
    tGrad = zeros(spectrumSize,3);
    sliceIn = zeros(spectrumSize,3);
    logs = zeros(spectrumSize,3);
    idx = 1:2;
    gamma = pghi_findgamma_fast(fftshift(windowOut),gammaFast);
    tGradMulConst = grainHop_smp*grainSizeOut_smp/gamma/2;
    fGradMul = @(fGrad) -gamma/(grainHop_smp*grainSizeOut_smp)*fGrad;
    tGradMul = @(tGrad) tGradMulConst*tGrad + tGradPlusConsts;
else
    error('wrong inversion algorithm');
end
    
%define array to store output signals
out_s = zeros(duration_smp,outChannels); 
if compTimeGran == 1
    out_t = zeros(duration_smp+grainSizeVarAmp_smp+grainSizeIn_smp,outChannels); 
end


% the synthesis loop uses two spectrogram inversion techniques (SPSI and
% RTPGHI). The implementation of such techniques is provided by the
% Phaseret library by Zden?k Pr?�a (which also requires ltfat library by 
% Peter L. S�ndergaard and Zden?k Pr?�a). Both techniques support real-time
% spectrogram inversion computation. The following code snippet is an example
% of Spectrogram inversion using Phaseret library working on an entire signal s.
% 
% % g=fftshift(hanning(M));
% % gd=gabdual(g,a,M);
% % gamma = pghi_findgamma(g,a,M); %needed only by RTPGHI
% % [c,Ls]=dgtreal(s,g,a,M,'timeinv');
% % crec=spsi(abs(c),a,M); % for SPSI
% % OR
% % crec=rtpghi(c,gamma,a,M); % for RTPGHI
% % [srec]=idgtreal(crec,gd,a,M,Ls,'timeinv');
%
% in the synthesis loop the above functions has been modified and simplified
% to work (iteratively) only on signle signal windows, enabling fast computation 
% of the synthesis output one grain at a time (real-time ready).

niter = ceil((duration_smp - grainSizeOut_smp) / grainHop_smp); % number of iterarions

%for each output channel a new synthesis loop is performed using a
%different input channel (if available) and recomputing the randomized 
% components of the synthesis loop
for j=1:outChannels 
    
    k=rem((j-1),inChannels)+1; %index of input channel to use
    
    %synthesis loop - the loop is real-time ready (it extract and process and 
    %output one grain at a time).
    for i=0:(niter - 1)
        
        %compute current file position (excluding sequential progress)
        if filePosModAmp ~= 0
            filePosMod_smp = filePos_smp + round(modulAmp_smp * sin(2*pi*filePosModFreq_f*i + modulAmpOffset_p)) + modulAmpOffset_smp; 
        end
        
        %compute current grain size
        if grainSizeInVarMode == 1 % 0 no variation, 1 random, sinusoidal modulated 
            grainSizeMod_smp = abs(round(grainSizeIn_smp + grainSizeVarAmp_smp * randn));
        elseif grainSizeInVarMode == 2
            grainSizeMod_smp = grainSizeIn_smp + round(grainSizeVarAmp_smp * sin(2*pi*grainSizeVarFreq_f*i)); 
        end
        if grainSizeMod_smp < 2
            grainSizeMod_smp = 2;
        end
        
        %extract current grain from source file depending on mode
        if grainSelMode == 1
            grain = getgrain_fix(snd(:,k), grainSizeMod_smp, filePosMod_smp, spray_smp, sprayOffset_smp); 
            elseif grainSelMode == 2
            grain = getgrain_seq(snd(:,k), grainSizeMod_smp, filePosMod_smp, spray_smp, sprayOffset_smp, sequenceRate, sndLength_smp, i);
            else
            grain = getgrain_rnd(snd(:,k), grainSizeMod_smp);
        end
        
        %recompute window if size has changed
        if (length(windowIn) ~= length(grain))
            windowIn = computewinfromhandler(windowhandler,kaiserBeta,length(grain),iskaiser); 
        end
        
        %apply window, compute the magnitude spectrum, and resize the
        %spectrum to match the output grain size
        magSpect = magspect_resampled(grain,windowIn,grainSizeOut_smp,spectrumSize);
        
        %apply circular shifting of frequency bins
        if freqShiftVarMode == 0 % 0 fixed, 1 random, sinusoidal modulated 
            magSpect = circshift(magSpect,freqShiftVarAmp_bin);
        elseif freqShiftVarMode == 1 
            magSpect = circshift(magSpect,round(freqShiftVarAmp_bin*randn));
        elseif freqShiftVarMode == 2
            magSpect = circshift(magSpect,round(freqShiftVarAmp_bin*sin(2*pi*freqShiftVarFreq_f*i)));
        end
        
        %reconstruct phase by spectrogram inversion algorithm
        if inversionAlgo == 0
            [reconFreq,startPhase] = comp_spsi(magSpect,grainHop_smp,grainSizeOut_smp,startPhase);
        elseif inversionAlgo == 1
            sliceInLog = log(magSpect+eps);
            logs(:,1:end-1) = logs(:,2:end);
            logs(:,end) = sliceInLog;
            fGrad = fGradMul((logs(:,3)-logs(:,1))/2);
            tGrad(:,1:end-1) = tGrad(:,2:end);
            tGrad(2:end-1,end) = tGradMul(conv2(logs(:,end),[1;0;-1],'valid'));
            newPhase = comp_rtpghiupdate(logs(:,idx),tGrad(:,idx),fGrad,newPhase,tol(1),grainSizeOut_smp);
            reconFreq = magSpect.*exp(1i*newPhase);      
        else
            error('wrong inversion algorithm');
        end
        
        %obtain reconstruced grain by inverse transform and apply dual window
        reconTime=single_idgtrealtimeinv(reconFreq,windowOutDual,grainSizeOut_smp);
        
        %overlap and add reconstructed grain
        out_idx = (i * grainHop_smp ) + 1; % index of current grain in output array
        out_s(out_idx:(out_idx + grainSizeOut_smp - 1),j) = out_s(out_idx:(out_idx + grainSizeOut_smp - 1),j) + reconTime; %add and overlap current grain to output

        %equivalent time domain granulation
        if compTimeGran == 1
            out_t(out_idx:(out_idx + length(grain) - 1),j) = out_t(out_idx:(out_idx + length(grain) - 1),j) + (grain.*windowIn); %add and overlap current grain to output
        end

    end

end

%normalize amplitude to 1
if normalize == 1
    out_s = out_s./(max(abs(out_s(:))));
    if compTimeGran == 1
        out_t = out_t./(max(out_t(:)));
    end
end

%plot waveform and spectrogram of synthesis output
if plots == 1
    display_waveforms(out_s,4,'Spectral Granulation');
    display_spectrograms(out_s,5,'Spectral Granulation');
    if compTimeGran == 1
        display_waveforms(out_t,6,'Time Domain Granulation');
        display_spectrograms(out_t,7,'Time Domain Granulation');
    end    
end

%playback (if more than 2 channels, it play only the first two)
if playback == 1 
    if outChannels <= 2
        sound(out_s,Fs_sys);
    else
        sound(out_s(:,1:2),Fs_sys);
    end
    pause(duration_sec+1);
    if compTimeGran == 1
        if outChannels <= 2
            sound(out_t,Fs_sys);
        else
            sound(out_t(:,1:2),Fs_sys);
        end
    end
end

% save output to file
if ~isempty(fileOut)
    audiowrite(fileOut,out_s,Fs_sys);
    if compTimeGran == 1
        [filepath,name,ext] = fileparts(fileOut);
        fileOutT = [filepath name '_time' ext];
        audiowrite(fileOutT,out_t,Fs_sys);
    end
end
        
        