function vn = interpolate_alpha1(x0, y0, v0, xn, yn)
% Perform our particular type of interpolation.
% Make an interpolator for height or roughness. We'll add safe far-away
% points so we'll always be interpolating and not extrapolating.
int = scatteredInterpolant([100*[-1 -1 1 1]'; x0], ...
    [100*[-1 1 -1 1]'; y0], ...
    [zeros(4, 1);      v0], 'linear');

% Perform the actual interpolation.
vn = int(xn, yn);
