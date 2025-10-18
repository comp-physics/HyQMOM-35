function tests = test_realizability
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    script_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(script_dir, 'src'));
    addpath(fullfile(script_dir, 'src', 'autogen'));
    addpath(fullfile(script_dir, 'tests', 'utils'));
    
    testCase.TestData.tol = 1e-4;
end

function test_realizability_S2_valid(testCase)
    S110 = 0.5;
    S101 = 0.3;
    S011 = 0.4;
    
    [S110r, S101r, S011r, S2r] = realizability('S2', S110, S101, S011);
    
    verifyGreaterThanOrEqual(testCase, S2r, -testCase.TestData.tol, 'S2 should be non-negative');
    verifyEqual(testCase, S110r, S110, 'AbsTol', testCase.TestData.tol, 'S110 should be unchanged');
    
    fprintf('PASS: realizability S2 valid test\n');
end

function test_realizability_S2_boundary(testCase)
    S110 = 0.99;
    S101 = 0.99;
    S011 = 0.99;
    
    [S110r, S101r, S011r, S2r] = realizability('S2', S110, S101, S011);
    
    verifyGreaterThanOrEqual(testCase, S2r, -testCase.TestData.tol, 'S2 should be non-negative');
    verifyLessThanOrEqual(testCase, abs(S110r), 1.0 + testCase.TestData.tol, 'S110 in bounds');
    verifyLessThanOrEqual(testCase, abs(S101r), 1.0 + testCase.TestData.tol, 'S101 in bounds');
    verifyLessThanOrEqual(testCase, abs(S011r), 1.0 + testCase.TestData.tol, 'S011 in bounds');
    
    fprintf('PASS: realizability S2 boundary test\n');
end

function test_realizability_2D_valid(testCase)
    S30 = 0.5; S40 = 3.5; S11 = 0.3;
    S21 = 0.2; S31 = 0.1; S12 = 0.15;
    S22 = 1.2; S03 = 0.4; S13 = 0.1; S04 = 3.2;
    
    [S21r, S12r, S31r, S22r, S13r] = realizability('2D', S30, S40, S11, S21, S31, S12, S22, S03, S13, S04);
    
    verifyTrue(testCase, isfinite(S21r), 'S21r should be finite');
    verifyTrue(testCase, isfinite(S22r), 'S22r should be finite');
    
    fprintf('PASS: realizability 2D valid test\n');
end

function test_realizability_S111(testCase)
    S110 = 0.5; S101 = 0.3; S011 = 0.4;
    S210 = 0.2; S201 = 0.15; S120 = 0.18;
    S021 = 0.12; S102 = 0.14; S012 = 0.16;
    S111 = 0.1;
    
    S111r = realizability('S111', S110, S101, S011, S210, S201, S120, S021, S102, S012, S111);
    
    verifyTrue(testCase, isfinite(S111r), 'S111r should be finite');
    
    fprintf('PASS: realizability S111 test\n');
end

function test_edge_corner_correction_integration(testCase)
    S110r = 0.99; S101r = 0.3; S011r = 0.4;
    R110 = -0.01; R101 = 1 - S101r^2; R011 = 1 - S011r^2;
    
    [S300r1, S030r1, S400r1, S040r1] = deal(0.5, 0.6, 3.0, 3.2);
    [S300r2, S003r2, S400r2, S004r2] = deal(0.5, 0.5, 3.0, 3.1);
    [S030r3, S003r3, S040r3, S004r3] = deal(0.6, 0.5, 3.2, 3.1);
    [S300r, S030r, S003r, S400r, S040r, S004r] = deal(0.5, 0.6, 0.5, 3.0, 3.2, 3.1);
    [S210r, S201r, S120r, S111r, S102r, S021r, S012r] = deal(0.1, 0.1, 0.1, 0.05, 0.1, 0.1, 0.1);
    [S310r, S301r, S220r, S211r, S202r, S130r, S121r, S112r, S103r] = deal(0.05, 0.05, 1.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05);
    [S031r, S022r, S013r] = deal(0.05, 0.05, 0.05);
    
    [S110, S101, S011, S300, S030, S003, S400, S040, S004, ...
     S210, S201, S120, S111, S102, S021, S012, ...
     S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
        edge_corner_correction(R110, R101, R011, S110r, S101r, S011r, ...
                               S300r1, S030r1, S400r1, S040r1, ...
                               S300r2, S003r2, S400r2, S004r2, ...
                               S030r3, S003r3, S040r3, S004r3, ...
                               S300r, S030r, S003r, S400r, S040r, S004r, ...
                               S210r, S201r, S120r, S111r, S102r, S021r, S012r, ...
                               S310r, S301r, S220r, S211r, S202r, S130r, S121r, S112r, S103r, ...
                               S031r, S022r, S013r);
    
    verifyEqual(testCase, abs(S110), 1.0, 'AbsTol', 0.01, 'S110 should be at edge');
    verifyTrue(testCase, isfinite(S300), 'S300 should be finite');
    
    S2 = 1 + 2*S110*S101*S011 - (S110^2 + S101^2 + S011^2);
    verifyGreaterThanOrEqual(testCase, S2, -testCase.TestData.tol, 'S2 realizability');
    
    fprintf('PASS: edge_corner_correction integration test\n');
end

function test_realizability_preserves_mass(testCase)
    rho = 1.5;
    u = 0.1; v = 0.2; w = 0.3;
    T = 1.0;
    
    M = InitializeM4_35(rho, u, v, w, T, 0, 0, T, 0, T);
    
    [~, ~, ~, M_real] = Flux_closure35_and_realizable_3D(M, 1, 0.5);
    
    verifyEqual(testCase, M_real(1), M(1), 'AbsTol', testCase.TestData.tol, ...
                'Realizability should preserve mass');
    
    fprintf('PASS: realizability preserves mass test\n');
end
