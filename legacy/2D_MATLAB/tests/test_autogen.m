function tests = test_autogen
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    script_dir = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(script_dir, 'src'));
    addpath(fullfile(script_dir, 'src', 'autogen'));
    addpath(fullfile(script_dir, 'tests', 'utils'));
    
    testCase.TestData.tol = 1e-10;
end

function test_M4toC4_3D_consistency(testCase)
    rho = 1.0;
    u = 0.1; v = 0.2; w = 0.3;
    C200 = 1.2; C110 = 0.3; C101 = 0.2;
    C020 = 1.5; C011 = 0.4; C002 = 1.8;
    
    M = InitializeM4_35(rho, u, v, w, C200, C110, C101, C020, C011, C002);
    m = moment_struct('from_vector', M);
    
    C_array = M4toC4_3D(m.M000, m.M100, m.M200, m.M300, m.M400, ...
                        m.M010, m.M110, m.M210, m.M310, m.M020, m.M120, m.M220, ...
                        m.M030, m.M130, m.M040, m.M001, m.M101, m.M201, m.M301, ...
                        m.M002, m.M102, m.M202, m.M003, m.M103, m.M004, ...
                        m.M011, m.M111, m.M211, m.M021, m.M121, m.M031, ...
                        m.M012, m.M112, m.M013, m.M022);
    
    verifyTrue(testCase, all(isfinite(C_array(:))), 'All central moments should be finite');
    
    fprintf('PASS: M4toC4_3D consistency test\n');
end

function test_S4toC4_3D_r_consistency(testCase)
    C200 = 1.2; C110 = 0.3; C101 = 0.2;
    C020 = 1.5; C011 = 0.4; C002 = 1.8;
    
    S300=0; S210=0; S201=0; S120=0; S111=0; S102=0; S030=0; S021=0; S012=0; S003=0;
    S400=3; S310=0; S301=0; S220=1; S211=0; S202=1;
    S130=0; S121=0; S112=0; S103=0; S040=3; S031=0; S022=1; S013=0; S004=3;
    
    C_array = S4toC4_3D_r(C200, C110, C101, C020, C011, C002, ...
                          S300, S210, S201, S120, S111, S102, S030, S021, S012, S003, ...
                          S400, S310, S301, S220, S211, S202, S130, S121, S112, S103, ...
                          S040, S031, S022, S013, S004);
    
    verifyTrue(testCase, all(isfinite(C_array(:))), 'All central moments should be finite');
    
    fprintf('PASS: S4toC4_3D_r consistency test\n');
end
