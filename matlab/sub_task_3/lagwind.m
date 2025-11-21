function wind = lagwind(lag, window)
%LAGWIND Compute lag window function for spectral estimation
%
%   WIND = LAGWIND(LAG, WINDOW) computes a lag window vector of length LAG
%   using the specified window function.
%
%   Inputs:
%       lag    - Number of lag elements (integer >= 1)
%       window - Window type string (default: 'uniform'):
%                'uniform' or 'u'    - Uniform window
%                'sasaki' or 's'     - Sasaki window
%                'priestley' or 'p'  - Priestley window
%                'parzen' or 'pa'    - Parzen window
%                'hamming' or 'h'    - Hamming window
%                'gaussian' or 'g'   - Gaussian distribution window
%                'daniell' or 'd'    - Daniell window
%
%   Outputs:
%       wind   - Lag window vector of length lag
%
%   Example:
%       w = lagwind(100, 'hamming');
%       plot(w);
%
%   References:
%       C. L. Nikias, A. P. Petropulu, Higher-Order Spectra Analysis:
%       A Nonlinear Signal Processing Framework, PTR Prentice Hall, 1993.
%
%   Copyright (c) 2000 Tom McMurray
%
%   See also: apodfun, bisp3cum

% Handle default arguments
if nargin < 2
    window = 'uniform';
end

% Validate lag
lag = round(lag);
if lag < 1
    error('lagwind:InvalidLag', 'lag must be >= 1');
end

% Special case: lag = 1
if lag == 1
    wind = 1;
    return
end

% Resolve window type
window = lower(num2str(window));
windowarr = {'uniform', 'sasaki', 'priestley', 'parzen', 'hamming', 'gaussian', 'daniell'};
windowind = find(strncmp(window, windowarr, length(window)));

if isempty(windowind)
    error('lagwind:InvalidWindow', 'Unknown window type: %s', window);
end

% Handle ambiguous match (e.g., 'p' matches both priestley and parzen)
if length(windowind) > 1
    window = 'priestley';
else
    window = windowarr{windowind};
end

lag1 = lag - 1;

% Compute lag window vector
switch window
    case 'uniform'
        wind = ones(1, lag);

    case 'sasaki'
        windlag = (0:lag1) / lag1;
        wind = sin(pi * windlag) / pi + cos(pi * windlag) .* (1 - windlag);

    case 'priestley'
        windlag = (1:lag1) / lag1;
        wind = [1, (sin(pi * windlag) / pi ./ windlag - cos(pi * windlag)) * 3 / pi^2 ./ windlag.^2];

    case 'parzen'
        fixlag12 = fix(lag1 / 2);
        fixlag121 = fixlag12 + 1;
        windlag0 = (0:fixlag12) / lag1;
        windlag1 = 1 - (fixlag121:lag1) / lag1;
        wind = zeros(1, lag);
        wind(1:fixlag121) = 1 - (1 - windlag0) .* windlag0.^2 * 6;
        wind(fixlag121+1:lag) = windlag1.^3 * 2;

    case 'hamming'
        wind = 0.54 + 0.46 * cos(pi * (0:lag1) / lag1);

    case 'gaussian'
        wind = [1, erfc(((1:lag1-1) / lag1 - 0.5) * 8 / sqrt(2)) / 2, 0];

    case 'daniell'
        windlag = (1:lag1) / lag1;
        wind = [1, sin(pi * windlag) / pi ./ windlag];
end

end
