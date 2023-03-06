function varargout = rho(P, C, omega)
%rho Frequency-wise stability margin.
%   [r,rinf,omega] = rho(P,C) or [r,rinf] = rho(P,C,omega) returns the
%   frequency-wise stability margin for the plant-controller pair [P,C] as
%   defined in Vinnicombe (2000, p.68).  This essentially gives a
%   frequency-wise measure of the plant-controller pair's robustness. rinf
%   is the infimum of r, unless [P,C] is unstable in which case it is set
%   to zero and a warning message is generated.  Note that this function
%   assumes a positive feedback convention, in line with the reference.
%
%   rho(P,C) or rho(P,C,omega) without any output arguments will display
%   the results in graphical form.
%
%   REFERENCE
%   1. Vinnicombe, G. (2000) Uncertainty and Feedback: H-Infinity 
%      Loop-Shaping and the v-Gap Metric. London: Imperial College
%
%   See also: sigma.


% Check we can check out Control System Toolbox.
licenseAvailable = logical(license('checkout','control_toolbox'));
assert(licenseAvailable, 'rho:missingLicense', ...
    ['Missing license: could not check out a license for Control ' ...
    'System Toolbox.']);

% Check (and identify) the input argument syntax.
narginchk(2, 3);
nargoutchk(0, 3);

noFrequencyArgGiven = nargin < 3;
graphicalOutputRequired = nargout == 0;

% Check that plant and controller are of the right form.
assert(goodClass(P), 'rho:badInputClass', ...
    'Bad input class: P must be an object of class ss, tf or zpk.');
assert(goodClass(C), 'rho:badInputClass', ...
    'Bad input class: C must be an object of class ss, tf or zpk.');
[ny,nu] = size(P);
assert(isequal(size(C), [nu ny]), 'rho:inconsistentDimensions', ...
    'Inconsistent dimensions: P and C must work together.');

% See whether the closed-loop system is stable.
pcIsStable = isstable(feedback(P, C, +1));

% Form the system matrix whose norm we are to minimize:
P = tidy(P);
C = tidy(C);
G = formGeneralizedPlant(P, C);
genPlantStable = isstable(G);
assert(isequal(pcIsStable, genPlantStable));


% Compute the singular values (and a frequency vector, if none supplied)
if noFrequencyArgGiven
    [allSingValues, omega] = sigma(G);    
else    
    allSingValues = sigma(G, matToRowVec(omega));
end
assert(logical(exist('omega', 'var')));

% Find the maximum singular values and the frequency-wise stability margin
maxSingValues = max(allSingValues, [], 1);
assert(isequal(numel(maxSingValues), numel(omega)));
maxSingValues = reshape(maxSingValues, size(omega));
r = 1.0 ./ maxSingValues;
if pcIsStable        
    rinf = 1.0 / norm(G, inf);
    rinf = min([min(r), rinf]); % prevent numerical inconsistency
else
    warning('rho:unstablePCPair', ...
        'Unstable [P,C] pair: stability margin is zero.');    
    rinf = 0.0;
end


% Produce output in the form required by the user.
if graphicalOutputRequired   
    buildGraphicalOutput(r, rinf, omega, G.Ts);
    varargout = {};
else
    varargout = {r, rinf, omega};
end
      
end



function buildGraphicalOutput(r, rinf, omega, Ts)
%buildGraphicalOutput Build a frequency-wise stability margin plot

% Form the core plot.
theAxis = gca();
semilogx(theAxis, omega, r)
xlabel(theAxis, 'Frequency (rad/s)');
ylabel('Frequency-Wise Stability Margin')
title(theAxis, 'Frequency-Wise Stability Margin');
ylim(theAxis, [0 1]);

% Add a dotted line at the infimum value.
hLine = line([omega(1) omega(end)], [rinf rinf]);
set(hLine, 'Color', [0 0 0], 'LineStyle', ':');

% If we've got a discrete time system, mark the Nyquist frequency.
plotShowsDiscrete = ~isequal(Ts, 0);
if plotShowsDiscrete
    nyquistFreq = pi/abs(Ts);
    hLine = line([nyquistFreq nyquistFreq], [0 1]);
    set(hLine, 'Color', [0 0 0], 'LineStyle', '-');
end
   
end



function v = matToRowVec(m)
%matToRowVec Convert matrix to row vector.

vcol = m(:);
v = vcol.';

end



function success = goodClass(obj)
%goodClass True for objects of class double, ss, tf or zpk.

objClass = class(obj);
success = ismember(objClass, {'double', 'ss', 'tf', 'zpk'});

end


function G = formGeneralizedPlant(P, C)
%formGeneralizedPlant Form the generalized plant.
%   This function forms the 'generalized plant' transfer function given in
%   Vinnicombe (2000, p.8), i.e.
%
%      [P;I] (I-C.P)^(-1) [-C;I]
%
%   To get a minimal realization, this should be formed using feedback
%   transformations. Our equations are:
%
%      y = P.u
%      z = v2 - y
%      u = v1 - C.z => v1 - C.v2 + C.y
% 
%   Putting these together in matrix form we get
%
%                              [ v2 ]
%      [ y ] = [  0  0  0  P ] [ v1 ]
%      [ u ]   [ -C  I  C  0 ] [ y  ]
%                              [ u  ]
%
%   From this, it is easy to eliminate the (y,u) pair in the input using a
%   linear fractional transformation.


[ny,nu] = size(P);

Iu = ss(eye(nu));
Iy = ss(eye(ny));
Iz = ss(eye(nu + ny));

Ou = zeros(nu);
Oy = zeros(ny);
Oyu = zeros(ny, nu);
Ouy = zeros(nu, ny);

G11 = [Oy Oyu Oy ss(P)];
G21 = [Ouy Iu Ouy Ou] + ss(C)*[-Iy Oyu Iy Oyu];

G = lft([Iz;Iz] * [G11;G21], Iz);

end



function obj1 = tidy(obj)
%tidy Convert a double, ss, zpk or tf object into a minimal ss object.

ssObj = ss(obj);
sminrealObj = sminreal(ssObj);
obj1 = minreal(sminrealObj, [], false);

end