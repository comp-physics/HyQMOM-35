function tests = test_diagnostics
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    script_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(script_dir, 'src'));
    addpath(fullfile(script_dir, 'src', 'autogen'));
    addpath(fullfile(script_dir, 'tests', 'utils'));
    
    testCase.TestData.tol = 1e-10;
end

function test_diagnostics_check2D_valid(testCase)
    S30 = 0.5; S40 = 3.5; S11 = 0.3;
    S21 = 0.2; S31 = 0.1; S12 = 0.15;
    S22 = 1.2; S03 = 0.4; S13 = 0.1; S04 = 3.2;
    
    [S30_out, S40_out, S11_out, S21_out, S31_out, S12_out, S22_out, S03_out, S13_out, S04_out] = ...
        diagnostics('check2D', S30, S40, S11, S21, S31, S12, S22, S03, S13, S04);
    
    verifyTrue(testCase, isfinite(S30_out), 'S30 should be finite');
    verifyTrue(testCase, isfinite(S40_out), 'S40 should be finite');
    verifyLessThanOrEqual(testCase, abs(S11_out), 1.0 + testCase.TestData.tol, ...
                          'S11 should be in bounds');
    
    fprintf('PASS: diagnostics check2D valid test\n');
end

function test_diagnostics_check2D_boundary(testCase)
    S30 = 0.5; S40 = 3.5; S11 = 1.5;
    S21 = 0.2; S31 = 0.1; S12 = 0.15;
    S22 = 1.2; S03 = 0.4; S13 = 0.1; S04 = 3.2;
    
    [S30_out, S40_out, S11_out, S21_out, S31_out, S12_out, S22_out, S03_out, S13_out, S04_out] = ...
        diagnostics('check2D', S30, S40, S11, S21, S31, S12, S22, S03, S13, S04);
    
    verifyEqual(testCase, abs(S11_out), 1.0, 'AbsTol', 0.01, ...
                'S11 should be corrected to boundary');
    
    fprintf('PASS: diagnostics check2D boundary test\n');
end

function test_diagnostics_check2D_all_planes(testCase)
    S300 = 0.5; S400 = 3.5; S110 = 0.3; S210 = 0.2; S310 = 0.1;
    S120 = 0.15; S220 = 1.2; S030 = 0.4; S130 = 0.1; S040 = 3.2;
    S101 = 0.25; S201 = 0.18; S301 = 0.08; S102 = 0.12; S202 = 0.09;
    S003 = 0.45; S103 = 0.11; S004 = 3.1; S011 = 0.28; S021 = 0.16;
    S031 = 0.09; S012 = 0.13; S022 = 1.15; S013 = 0.10;
    
    [S300r1, S400r1, S110r, S210r, S310r, S120r, S220r, S030r1, S130r, S040r1, ...
     S300r2, S400r2, S101r, S201r, S301r, S102r, S202r, S003r2, S103r, S004r2, ...
     S030r3, S040r3, S011r, S021r, S031r, S012r, S022r, S003r3, S013r, S004r3] = ...
        diagnostics('check2D_all_planes', S300, S400, S110, S210, S310, S120, S220, S030, S130, S040, ...
                    S101, S201, S301, S102, S202, S003, S103, S004, S011, S021, S031, S012, S022, S013);
    
    verifyTrue(testCase, isfinite(S300r1), 'All outputs should be finite');
    verifyLessThanOrEqual(testCase, abs(S110r), 1.0 + testCase.TestData.tol, 'S110 in bounds');
    verifyLessThanOrEqual(testCase, abs(S101r), 1.0 + testCase.TestData.tol, 'S101 in bounds');
    verifyLessThanOrEqual(testCase, abs(S011r), 1.0 + testCase.TestData.tol, 'S011 in bounds');
    
    fprintf('PASS: diagnostics check2D all planes test\n');
end

