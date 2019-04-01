function [r,theta] = global2head(L1,L2,L4,t1,t2,t4)
%% global2head: transforms reference stimulus into head coordinate frame
%   INPUT:
%       bodyAng     : body angles
%       headAng     : head angles
%       refAng      : reference angles
%       tout        : time vector
%       rootdir     : directory to save movie >>> varargin{1} 
%   OUTPUT:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% r = sqrt(L1^(2) + 2*sind(t1+t2)*L1*L2 - 2*cosd(t1-t4)*L1*L4 + L2^(2) - 2*sind(t2+t4)*L2*L4 + L4^(2));
% ang = atand((L4*cosd(t4) - L1*cosd(t1) - L2*sind(t2))/(L4*sind(t4) - L1*sind(t1) - L2*cosd(t2)));

x = L4*sind(t4) - ( L1*sind(t1) + L2*sind(t2) );
y = L4*cosd(t4) - ( L1*cosd(t1) + L2*cosd(t2) );

r = sqrt(x.^(2) + y.^(2));
theta = atan2d(x,y);

% 8*sind(0) - 3*sind(0) - 2*cosd(90)
% 8*cosd(0) - 3*cosd(0) - 2*sind(90)

end