function ybar = Differ(y,delta)
%
% This code implements the difference operation
%
% Inputs:
%    y:  a signal want to calculate its difference
%    delta:  sampling time interval of y 
% Outputs:
%    ybar:  the signal after difference operation
%
% Authors: Hao Liang
% Last modified by: 21/10/02
%

% Parameter setting
L = length(y);
ybar = zeros(1,L-2);

% Difference operation
for i = 2 : L-1
    ybar(i-1)=(y(i+1)-y(i-1))/(2*delta);
end
ybar = [(y(2)-y(1))/delta,ybar,(y(end)-y(end-1))/delta];

end