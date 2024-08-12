function delta = angdiffdeg(x, y, varargin)
% calculate difference between angles in degree
% if 'absolute', delta is the absolute difference between x, y
% range from [0, 180]

% Copyright Tsam Kiu Pun, 2024. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------

if nargin > 2 
    isnorm = strcmp(varargin(1), 'absolute');
else
    isnorm = false;
end
if isnorm
    normDeg  = mod(x - y,360); 
    delta = min(360 - normDeg , normDeg);
else
    delta = angle(exp(1i*x * pi/180) ./ exp(1i*y * pi/180)) * 180/pi;
end
end
