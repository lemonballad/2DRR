function [apod] = apodfun(x, x0, gL, gH, wid, invert)
%APODFUN Generate asymmetric bandpass or notch filter
%
%   [APOD] = APODFUN(X, X0, GL, GH, WID, INVERT) generates an apodization
%   window with Gaussian decay on each side and a flat top region.
%
%   Inputs:
%       x      - Domain vector (frequencies or time points)
%       x0     - Center position of the filter
%       gL     - Gaussian decay width on the low side
%       gH     - Gaussian decay width on the high side
%       wid    - Total width of the flat top region
%       invert - If true (1), create notch filter; if false (0), bandpass
%
%   Outputs:
%       apod   - Apodization window array (same size as x)
%
%   Example:
%       x = linspace(0, 100, 1000);
%       window = apodfun(x, 50, 10, 10, 20, false);
%       plot(x, window);
%
%   See also: lagwind, bisp3cum

% Calculate flat top boundaries
locent = x0 - (wid / 2);
hicent = x0 + (wid / 2);

% Vectorized computation
fu = zeros(size(x));

% Region 1: below low center - Gaussian decay
mask_low = x < locent;
fu(mask_low) = exp(-0.5 * (x(mask_low) - locent).^2 / gL^2);

% Region 2: flat top
mask_flat = (x >= locent) & (x <= hicent);
fu(mask_flat) = 1;

% Region 3: above high center - Gaussian decay
mask_high = x > hicent;
fu(mask_high) = exp(-0.5 * (x(mask_high) - hicent).^2 / gH^2);

% Apply inversion if requested
if invert
    apod = 1 - fu;
else
    apod = fu;
end

end
