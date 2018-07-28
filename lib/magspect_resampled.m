% This file is part of
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

function out = magspect_resampled(signal, window, grainSizeOut_smp, spectrumSize)

    if length(signal) == grainSizeOut_smp
        out = abs(single_dgtrealtimeinv(signal,window,grainSizeOut_smp,spectrumSize));
    else
        out = (single_dgtrealtimeinv(signal,window,length(signal),floor(length(signal)/2) + 1));
        out = abs(interp1(1:length(out),out,linspace(1,length(out),spectrumSize),'spline'))';
    end  

end