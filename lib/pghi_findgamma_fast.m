% The code in this file has been adapted from
% Copyright (C) 2016 Zdenek Prusa <zdenek.prusa@gmail.com>.
% This file is part of PHASERET version 0.2.1
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%   AUTHORS: Zdenek Prusa
%


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


function [gamma, Cg, g] = pghi_findgamma_fast(g,mode)%mode 1 is fast 0 slow
    atheight = findbestgauss(g,mode);
    w = winwidthatheight(g, atheight);
    
    gl = length(g);

    Cg = -pi/4*(w/(gl-1))^2/log(atheight);
    gamma = Cg*gl^2;
end


function [atheight,minorm] = findbestgauss(gnum , mode)

if mode == 1
    atheightrange = 0.01:0.01:0.8;
else
    atheightrange = 0.01:0.001:0.8;
end
    
w = winwidthatheight(gnum, atheightrange);

L = 10*numel(gnum);

gnum = fir2long(normalize(gnum,'inf'),L);
norms = zeros(size(atheightrange));
tfrs = zeros(size(atheightrange));
for ii=1:numel(atheightrange)
    [gausstmp,tfrs(ii)] = pgauss(L,'inf','width',w(ii),'atheight',atheightrange(ii));
    norms(ii) = norm(gnum-gausstmp);
end

[~,idx]=min(norms);

atheight = atheightrange(idx);
minorm = norms(idx);

end

function width = winwidthatheight(gnum,atheight)

if ~isnumeric(gnum) || isempty(gnum) || ~isvector(gnum)
    error('%s: gnum must be a numeric vector.', upper(mfilename));
end

if isempty(atheight) || any(atheight) > 1 || any(atheight) < 0
    error('%s: h must be in the interval [0-1].', upper(mfilename));
end

width = zeros(size(atheight));
for ii=1:numel(atheight)
    gl = numel(gnum);
    gmax = max(gnum);
    frac=  1/atheight(ii);
    fracofmax = gmax/frac;


    ind =find(gnum(1:floor(gl/2)+1)==fracofmax,1,'first');
    if isempty(ind)
        %There is no sample exactly half of the height
        ind1 = find(gnum(1:floor(gl/2)+1)>fracofmax,1,'last');
        ind2 = find(gnum(1:floor(gl/2)+1)<fracofmax,1,'first');
        if isempty(ind2)
            width(ii) = gl;
        else
            rest = 1-(fracofmax-gnum(ind2))/(gnum(ind1)-gnum(ind2));
            width(ii) = 2*(ind1+rest-1);
        end
    else
        width(ii) = 2*(ind-1);
    end
end

end
