function fidel = fidelity(chan, trueS)

% h*cosA = b, 
% where h is channel amplitude of a point, 
% and A represents the angle between true and channel point
s_precomp = linspace(0, 2*pi, length(chan))';
fidel = cos(trueS - s_precomp)'*chan; 


