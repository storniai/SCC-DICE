function R = icdfRoeBaker(seeds,fbar,fsig,maxT,cs0)
% Copyright (C) 2011 by Robert E. Kopp; distributed under GNU GPL v3
% 
% Takes seeds (between 0 and 1) and inverts the CDF from a Roe-Baker distribution with
% feedback factor with mean fbar and standard deviation fsig, and truncated at the right
% at maxT.
%
% Defaults (based on USG 2010 SCC report):
%     fbar = 0.61979
%     fsig = 0.18407
%     maxT = 10
%     cs0 = 1.2
%
if isfinite(maxT)
	fmax = 1-(cs0/maxT);
	pmax = 0.5 * erfc(-(fmax-fbar)/fsig/sqrt(2));
	seeds = seeds*pmax;
end
f = -sqrt(2).*erfcinv(2*seeds).*fsig + fbar;
R = cs0./(1-f);

