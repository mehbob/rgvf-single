function [xi,yi] = snakeinterp(x,y,dmax,dmin)
%SNAKEINTERP  Interpolate the snake adaptively
%   [xi,yi] = snakeinterp(x,y,dmax,dmin)
%
%   dmax: the maximum distance between two snake points
%   dmin: the maximum distance between two snake points
%   d(i,i+1)>dmax, then a new point is added between i and i+1
%   d(i,i+1)<dmin, then either i or i+1 is removed 
%  
%   NOTE: the spacing of original curve must be close to the 
%         range defined by dmax and dmin. For arbitrary spacing,
%         try snakeinterp1.
% 
%   See also SNAKEINTERP1

%    there is a bug in the program for points removal

%   Chenyang Xu and Jerry L. Prince, 4/1/95, 6/17/97
%   Copyright (c) 1995-97 by Chenyang Xu and Jerry L. Prince
%   Image Analysis and Communications Lab, Johns Hopkins University

% convert to column vector
x = x(:); y = y(:);

N = length(x);

if(N <= 2)
    xi = 0;
    yi = 0;
    return
end

d = abs(x([2:N 1])- x(:)) + abs(y([2:N 1])- y(:));

% remove the points which distance to neighbor points is shorter than dmin
IDX = (d<dmin);

idx = find(IDX==0);
x = x(idx);
y = y(idx);

isChanged = 1;
while(isChanged==1)

Nlik = length(x);
lik = ones(1,Nlik);
isChanged = 0;
for slik = 3 : Nlik-2
    dlik = abs(x(slik) - x([(slik + 2): Nlik,1:(slik -2)]))+ abs(y(slik) - y([(slik + 2): Nlik,1:(slik -2)]));
    dlikmin = min(dlik);
    if(dlikmin<=dmin)
        idxlik = find(dlik == dlikmin);
        idxlikNext = slik + 1 + idxlik(1);
        if(idxlikNext > Nlik)
            idxlikNext = idxlikNext - Nlik;
        end
        StartPt = min(slik,idxlikNext);
        StopPt = max(slik,idxlikNext);
        
        if(StopPt - StartPt < Nlik/2)
            lik([StartPt:StopPt]) = 0;
            isChanged = 1;
        end
        %save likfile lik;
        x = x(find(lik ~= 0));
        y = y(find(lik ~= 0));
        break;
    end
end

end
 
%%%%%%%%%%%%%%%%%%%%%

N = length(x);

if(N <= 2)
    xi = 0;
    yi = 0;
    return
end

d = abs(x([2:N 1])- x(:)) + abs(y([2:N 1])- y(:));


IDX = (d>dmax);

z = snakeindex(IDX);

p = 1:N+1;

xi = interp1(p,[x;x(1)],z');
yi = interp1(p,[y;y(1)],z');

N = length(xi);
d = abs(xi([2:N 1])- xi(:)) + abs(yi([2:N 1])- yi(:));

while (max(d)>dmax),

    IDX = (d>dmax);
    z = snakeindex(IDX);

    p = 1:N+1;

    xi = interp1(p,[xi;xi(1)],z');
    yi = interp1(p,[yi;yi(1)],z');

    N = length(xi);
    d = abs(xi([2:N 1])- xi(:)) + abs(yi([2:N 1])- yi(:));
end
