%TEST_SPECTRAL Tests for core spectral functions
%
%   Run with: runtests('test_spectral')
%
%   See also: apodfun, lagwind

classdef test_spectral < matlab.unittest.TestCase

    methods (Test)

        function test_apodfun_flat_top(testCase)
            % Flat top region should be 1.0
            x = linspace(0, 100, 1001);
            result = apodfun(x, 50, 10, 10, 20, false);

            % Check flat top region (40-60)
            flat_mask = (x >= 40) & (x <= 60);
            testCase.verifyEqual(result(flat_mask), ones(1, sum(flat_mask)), 'AbsTol', 1e-10);
        end

        function test_apodfun_gaussian_decay(testCase)
            % Edges should decay
            x = linspace(0, 100, 1001);
            result = apodfun(x, 50, 10, 10, 20, false);

            testCase.verifyLessThan(result(1), 0.01);
            testCase.verifyLessThan(result(end), 0.01);
        end

        function test_apodfun_invert(testCase)
            % Inverted filter should be 1 - normal
            x = linspace(0, 100, 1001);
            normal = apodfun(x, 50, 10, 10, 20, false);
            inverted = apodfun(x, 50, 10, 10, 20, true);

            testCase.verifyEqual(normal + inverted, ones(size(x)), 'AbsTol', 1e-10);
        end

        function test_lagwind_uniform(testCase)
            % Uniform window should be all ones
            w = lagwind(10, 'uniform');
            testCase.verifyEqual(w, ones(1, 10));
        end

        function test_lagwind_hamming_start(testCase)
            % Hamming window should start at 1.0
            w = lagwind(10, 'hamming');
            testCase.verifyEqual(w(1), 1.0, 'AbsTol', 1e-10);
        end

        function test_lagwind_output_length(testCase)
            % Output length should equal lag
            for lag = [5, 10, 20, 100]
                w = lagwind(lag, 'uniform');
                testCase.verifyEqual(length(w), lag);
            end
        end

        function test_lagwind_single(testCase)
            % Lag of 1 should return 1
            w = lagwind(1, 'parzen');
            testCase.verifyEqual(w, 1);
        end

        function test_lagwind_all_types(testCase)
            % All window types should produce valid output
            windows = {'uniform', 'sasaki', 'priestley', 'parzen', ...
                       'hamming', 'gaussian', 'daniell'};
            for i = 1:length(windows)
                w = lagwind(20, windows{i});
                testCase.verifyEqual(length(w), 20);
                testCase.verifyFalse(any(isnan(w)));
            end
        end

    end

end
