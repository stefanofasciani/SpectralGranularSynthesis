% The code in this file has been adapted from
% Copyright (C) 2005-2016 Peter L. Soendergaard <peter@sonderport.dk>.
% This file is part of LTFAT version 2.3.1
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

%   AUTHOR : Peter L. Søndergaard.

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


function [gd]=gabdualwindow(g,a,M)
    
    L=lcm(M,a);
    N=L/a;
    glong2=abs(fir2long_fast(g,L)).^2;
    temp=repmat(sum(reshape(glong2,a,N),2),N,1)*M;
    gd=g./long2fir_fast(temp,M);
    
end

function gout=fir2long_fast(gin,L)

     Lfir=length(gin);
    
     if rem(Lfir,2)==0
         gout=middlepadhp_fast(gin,L);
     else
         gout=middlepadwp_fast(gin,L);
     end
    
end

function gout=long2fir_fast(gin,L)

    if rem(L,2)==0
        gout=middlepadhp_fast(gin,L);
    else
        gout=middlepadwp_fast(gin,L);
    end
    
end

function f=middlepadwp_fast(f,L)

Ls=length(f);
Lorig=Ls;

if L~=Ls
    % ---------------   WPE case --------------------------------------
    
    if Lorig==1
      % Rather trivial case
      f=[f(1,:);zeros(L-1,1)];
      
    else
      if Lorig>L
        % Cut
        
        if mod(L,2)==0
          
          % L even. Use average of endpoints.
          f=[f(1:L/2,:);(f(L/2+1,:)+f(Lorig-L/2+1,:))/2;f(Lorig-L/2+2:Lorig,:)];
          
        else
          
          % No problem, just cut.
          f=[f(1:(L+1)/2,:);f(Lorig-(L-1)/2+1:Lorig,:)];
          
        end    
        
      else
        
        d=L-Lorig;
        
        % Extend
        if mod(Lorig,2)==0
          
          % Lorig even. We must split a value.
          
          f=[f(1:Lorig/2,:);...
             f(Lorig/2+1,:)/2;...
             zeros(d-1,1);...
             f(Lorig/2+1,:)/2;...
             f(Lorig/2+2:Lorig,:)];
          
        else
          % Lorig is odd, we can just insert zeros.
          f=[f(1:(Lorig+1)/2,:);zeros(d,1);f((Lorig+3)/2:Lorig,:)];
          
        end
        
      end
    end
    
end
      
end

function f=middlepadhp_fast(f,L)

Ls=length(f);
Lorig=Ls;

if L~=Ls
    
    % ------------------ HPE case ------------------------------------
    
    if Lorig==1
        f=[f(1,:);zeros(L-1,1)];
    else
      if Lorig>L
        
        d=Lorig-L;
        % Cut
        
        if mod(L,2)==0
          % L even
          
          % No problem, just cut.
          f=[f(1:L/2,:);...
             f(Lorig-L/2+1:Lorig,:);];
          
        else
          
          % Average of endpoints.
          f=[f(1:(L-1)/2,:);(f((L+1)/2,:)+f(Lorig-(L-1)/2,:))/2;...
             f(Lorig-(L-1)/2+1:Lorig,:);];
          
        end;
        
      else
        
        d=L-Lorig;
        
        % Extend
        if mod(Lorig,2)==0 
          
          % Lorig even. We can just insert zeros in the middle.
          
          f=[f(1:Lorig/2,:);...
             zeros(d,1);...
             f(Lorig/2+1:Lorig,:)];
          
        else
          % Lorig odd. We need to split a value in two
          f=[f(1:(Lorig-1)/2,:);...
             f((Lorig+1)/2,:)/2;...
             zeros(d-1,1);...
             f((Lorig+1)/2,:)/2;...
             f((Lorig-1)/2+2:Lorig,:)];
          
        end
        
      end
      
    end

end

end