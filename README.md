# Spectral Granular Synthesis

The Spectral Granular Synthesis software can be obtained at http://stefanofasciani.com/sgs.html

Spectral Granular Synthesis (C) 2018 Stefano Fasciani, University of Wollongong in Dubai  
Inquiries: stefano@fasciani.xyz


## Description

The MATLAB software included in this repository implements the Spectral Granular Synthesis algorithm described in  
S. Fasciani, "Spectral Granular Synthesis" in proceedings of International Computer Music Conference 2018, Daegu, Korea.

The software has been developed using MATLAB R2017b.

The reporitory includes also MATLAB app and a standalone executable application for WIN and OSX (it requires to install the free MATLAB Runtime)

## Instructions

### MATLAB Script
a) Include the folder 'lib' and its subfolders to the MATLAB search path.  
b) If the library Phaseret (http://ltfat.github.io/phaseret/) is already included included in the MATLAB search path, the subfolder 'lib/phaseret' is not required.  
c) If the library Ltfat (http://ltfat.github.io/) is already included included in the MATLAB search path, the subfolder 'lib/ltfat' is not required.  
d) Specify synthesis parameters and options from line 35 to 69 of SpectralGranularSinth.m.  
e) Enter 'SpectralGranularSinth' in the MATLAB command window or simply run the script 'SpectralGranularSinth.m'.  

### MATLAB Application
a) to c) as above.  
d) Enter 'SpectralGranularSinthGUI' in the MATLAB command window or simply double click on the file 'SpectralGranularSinthGUI.mlapp'.  
e) Specify synthesis parameters and options in the GUI.  
f) Click on the 'Run' button.  

### MATLAB Packaged Application
a) Open the 'SpectralGranularSinth.mlappinstall' in 'Release\MATLABapp' and install the application (it includes all dependancies).  
b) In te MATLAB window select the APPS tab and select 'SpectralGranularSinth' under 'My Apps'. Then follow e) to f) from above.  

### Standalone
The three methos above (recommended) require an existing MATLAB installation. A standalone application is provided as well for both WIN and OSX.  
The executables are 'Release/WIN/files_only/SpectralGranularSinth.exe' and 'Release/OSX/files_only/SpectralGranularSinth.app'. Both require a prior installation of MATLAB Runtime (free, available at http://www.mathworks.com/products/compiler/mcr/index.html), as described in the readme.txt file. ALternatively, bundled web installers are provided in 'Release/WIN/installer' and 'Release/OSX/installer' (not recommended because it takes quite some time to download and install the MATLAB Runtime).  
The standalone application has an interface identical to the MATLAB application described above.


## Synthesis parameters and options
Fields that are inactive have no effect with the current settings.

**Sampling Rate** sampling frequency (Hz) for the sound synthesis, playback and Output File. The Source File is resampled to match this rate.

**Output Duration** duration (Seconds) of the sound synthesis output.

**Output Channels** number of output channels. The sound synthesis is repeated for each output channel using the same parameters, using cyclicly a different channel of the Source File, and recomputing the random component of the synthesis algorithm.

**Source File** absolute path of the Source File from which grains are extracted (supported .wav and .aiff)

**Output File** absolute path of the Output File in which the synthesis output is stored (supported .wav and .aiff). If left empty, the Output File is not generated.

**Normalize Output** If enabled, the synthesis output is normalized to 0dB (i.e. unitaty maximum amplitude of output signal)

**Play Output** If enabled, the the synthesys output is reproduced after synthesis ends. If the number of output channels is greater than two, the playback is limited to the first two channels.

**Comp Time Vers.** If enabled, it computes also the equivalent Time Domain Version of the granular synthesis. This comparison is not meaningful when Input and Output grain sizes are different or varied throughtout the synthesis, or when spectrum shift is applied. If playback is enables, the Time Version is reproduced after the Spectral Version. If an Output File is specified, the Time Domain Version is saved on a separate file with appendix '_time.wav'. If Plots are enabled, the Time Domain Version is displayed as well.

**Plots** If enabled, it displays the window envelopes, waveform and spectrogram of synthesis input and output files. 

**PGHI Gamma Fast** If selected, it computes a less accurate but faster gamma for the PGHI. Gamma is computed only once outside the synthesis loop.

**Input Grain Size** Size of the input grains (Samples) extracted from the Source File.

**Output Grain Size** Size of the output grains (Samples) sequenced in the synthesis output.

**Grain Overlap** Overlap (Samples) between consecutive output grains. This can not exceed the Output Grain Size

**Window Type** Type of window used as amplitude envelope for the grains. Available: Hanning, Hamming, Blackman, Kaiser

**Kaiser Beta** Parameter determining the shape of the kaiser window (when equal to 0 the window is rectangular).

**Phase Recognition Algorithm** Algorithm used to reconstruct the phase from the Spectra of the grain sequence. Available SPSI, PGHI, or Off (using zero phase).

**PGHI Threshold** Threshold for the PGHI alsorithm, equal to 10^-(value).

**Grain Extraction** Mode of grain extraction from Source File: Fixed position (allows random deviation and modulatio), Sequential (allows random deviation and modulatio), Random. 

**Sequential Speed** If Sequential Grain Extraction mode is selected, this parameter determine the rate at which the extraction position is incremented dring synthesis (if set to 1 and In Out grain sizes matche, the synthesys is close to the normal speed Source File playback). Negative values are allowed and determine a backward movement of the extraction position. When the end of file is reached, the extraction continues from the beginning of the file.

**File Position** If Fixed Grain Extraction mode is selected, this determines the position from which grains are extracted. If Sequential Grain Extraction mode is selected, this determines the starting position for the sequential extraction.

**Position Random Deviation** Random deviation from the fixed or sequential position from which the frains are extracted. If equal to 0, no random component is added to the current position. The maximum allowed value is 1 which represent the greatest distance between the selected/current position and the file eventual limitations (e.g. file beginning or end). Not available for Random extraction mode.

**Direction** Direction of the random component with respect to the current position. Available Backward & Forward, Forward, or Backward. This selection setting also affects the maximum absolute deviation. Not available for Random extraction mode.

**Position Modulation Amplitude** Amplitude of the sinusoidal modulation of the grain extraction position. If equal to 0 no position modulation is applied. The maximum allowed value is 1 which represent the greatest distance between the selected/current position and the file eventual limitations (e.g. file beginning or end). Not available for Random extraction mode.

**Position Modulation Frequency** Frequency (Hz) of sinusoidal modulation of the grain extraction position. Not available for Random extraction mode.

**In Grain Size Variation** Enables the variation the Input Grain Size, choosing between Off, Random variation, or sinusoidally Modulated.

**Size Variation Amount** Deviation of the random Input Grain Size variation, or Amplitude of the sinusoidal modulation of the Input Grain Size variation. The value (0 to 1) represent the percentage of Input Grain Size.

**Size Variation Frequency** Frequency (Hz) of sinusoidal modulation of Input Grain Size. Not available for Random variation.

**Spectrum Shift Mode** Set the mode for the circular shift of the bins in the spectrum of the grains, choosing between Fixed shift, Random shift, or sinusoidally Modulated shift.

**Spectrum Shift Amount** Deviation of the random Spectrum Shift, or IAmplitude of the sinusoidal modulation of the Spectrum Shift variation. The value (0 to 1) represent the percentage of thr Spectrum Size, which is floor(Output Grain Size/2).

**Spectrum Shift Frequency** Frequency (Hz) of sinusoidal modulation of the Spectrum shift. Not available for Random shift.

**Run** Runs the synthesis and generate outputs according to parameters and options. Led is red until done (or in case of error).











