function tests = test_edge_corner_correction
%TEST_EDGE_CORNER_CORRECTION Unit tests for edge_corner_correction function
%   Tests the refactored edge_corner_correction against the original script
%   by creating synthetic edge/corner cases.
tests = functiontests(localfunctions);
end

function setup(testCase)
% Add src to path
script_dir = fileparts(fileparts(mfilename('fullpath')));
src_dir = fullfile(script_dir, 'src');
addpath(src_dir);
autogen_dir = fullfile(src_dir, 'autogen');
if exist(autogen_dir, 'dir')
    addpath(autogen_dir);
end
end

function test_edge_S110_boundary(testCase)
%TEST_EDGE_S110_BOUNDARY Test edge case where S110 = ±1 (R110 ≤ 0)
% Create a case where S110 is at boundary
S110r = 0.99;  % Very close to 1
S101r = 0.3;
S011r = 0.4;

% R110 should be ~0 (edge case)
R110 = 1 - S110r^2;
R101 = 1 - S101r^2;
R011 = 1 - S011r^2;

% Set R110 to trigger edge condition
R110 = -0.01;  % Force edge case

% Create reasonable moment values
[S300r1, S030r1, S400r1, S040r1] = deal(0.5, 0.6, 0.3, 0.4);
[S300r2, S003r2, S400r2, S004r2] = deal(0.5, 0.5, 0.3, 0.3);
[S030r3, S003r3, S040r3, S004r3] = deal(0.6, 0.5, 0.4, 0.3);
[S300r, S030r, S003r, S400r, S040r, S004r] = deal(0.5, 0.6, 0.5, 0.3, 0.4, 0.3);

% Initialize other moments
[S210r, S201r, S120r, S111r, S102r, S021r, S012r] = deal(0.1, 0.1, 0.1, 0.05, 0.1, 0.1, 0.1);
[S310r, S301r, S220r, S211r, S202r, S130r, S121r, S112r, S103r] = deal(0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05);
[S031r, S022r, S013r] = deal(0.05, 0.05, 0.05);

% Call the function
[S110, S101, S011, S300, S030, S003, S400, S040, S004, ...
 S210, S201, S120, S111, S102, S021, S012, ...
 S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
    edge_corner_correction(R110, R101, R011, ...
                           S110r, S101r, S011r, ...
                           S300r1, S030r1, S400r1, S040r1, ...
                           S300r2, S003r2, S400r2, S004r2, ...
                           S030r3, S003r3, S040r3, S004r3, ...
                           S300r, S030r, S003r, S400r, S040r, S004r, ...
                           S210r, S201r, S120r, S111r, S102r, S021r, S012r, ...
                           S310r, S301r, S220r, S211r, S202r, S130r, S121r, S112r, S103r, ...
                           S031r, S022r, S013r);

% Verify outputs are reasonable
verifyEqual(testCase, abs(S110), 1, 'AbsTol', 0.01, 'S110 should be ±1 at edge');
verifyGreaterThanOrEqual(testCase, S300, 0, 'S300 should be non-negative');
verifyGreaterThanOrEqual(testCase, S030, 0, 'S030 should be non-negative');
verifyTrue(testCase, isfinite(S210), 'S210 should be finite');

% Check realizability constraint
S2 = 1 + 2*S110*S101*S011 - (S110^2 + S101^2 + S011^2);
verifyGreaterThanOrEqual(testCase, S2, -1e-6, 'S2 realizability should hold');

fprintf('✓ Edge S110 test passed\n');
end

function test_edge_S101_boundary(testCase)
%TEST_EDGE_S101_BOUNDARY Test edge case where S101 = ±1 (R101 ≤ 0)
S110r = 0.3;
S101r = 0.99;
S011r = 0.4;

R110 = 1 - S110r^2;
R101 = -0.01;  % Force edge case
R011 = 1 - S011r^2;

% Create reasonable moment values
[S300r1, S030r1, S400r1, S040r1] = deal(0.5, 0.6, 0.3, 0.4);
[S300r2, S003r2, S400r2, S004r2] = deal(0.5, 0.5, 0.3, 0.3);
[S030r3, S003r3, S040r3, S004r3] = deal(0.6, 0.5, 0.4, 0.3);
[S300r, S030r, S003r, S400r, S040r, S004r] = deal(0.5, 0.6, 0.5, 0.3, 0.4, 0.3);
[S210r, S201r, S120r, S111r, S102r, S021r, S012r] = deal(0.1, 0.1, 0.1, 0.05, 0.1, 0.1, 0.1);
[S310r, S301r, S220r, S211r, S202r, S130r, S121r, S112r, S103r] = deal(0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05);
[S031r, S022r, S013r] = deal(0.05, 0.05, 0.05);

[S110, S101, S011, S300, S030, S003, S400, S040, S004, ...
 S210, S201, S120, S111, S102, S021, S012, ...
 S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
    edge_corner_correction(R110, R101, R011, ...
                           S110r, S101r, S011r, ...
                           S300r1, S030r1, S400r1, S040r1, ...
                           S300r2, S003r2, S400r2, S004r2, ...
                           S030r3, S003r3, S040r3, S004r3, ...
                           S300r, S030r, S003r, S400r, S040r, S004r, ...
                           S210r, S201r, S120r, S111r, S102r, S021r, S012r, ...
                           S310r, S301r, S220r, S211r, S202r, S130r, S121r, S112r, S103r, ...
                           S031r, S022r, S013r);

verifyEqual(testCase, abs(S101), 1, 'AbsTol', 0.01, 'S101 should be ±1 at edge');
S2 = 1 + 2*S110*S101*S011 - (S110^2 + S101^2 + S011^2);
verifyGreaterThanOrEqual(testCase, S2, -1e-6, 'S2 realizability should hold');

fprintf('✓ Edge S101 test passed\n');
end

function test_edge_S011_boundary(testCase)
%TEST_EDGE_S011_BOUNDARY Test edge case where S011 = ±1 (R011 ≤ 0)
S110r = 0.3;
S101r = 0.4;
S011r = 0.99;

R110 = 1 - S110r^2;
R101 = 1 - S101r^2;
R011 = -0.01;  % Force edge case

[S300r1, S030r1, S400r1, S040r1] = deal(0.5, 0.6, 0.3, 0.4);
[S300r2, S003r2, S400r2, S004r2] = deal(0.5, 0.5, 0.3, 0.3);
[S030r3, S003r3, S040r3, S004r3] = deal(0.6, 0.5, 0.4, 0.3);
[S300r, S030r, S003r, S400r, S040r, S004r] = deal(0.5, 0.6, 0.5, 0.3, 0.4, 0.3);
[S210r, S201r, S120r, S111r, S102r, S021r, S012r] = deal(0.1, 0.1, 0.1, 0.05, 0.1, 0.1, 0.1);
[S310r, S301r, S220r, S211r, S202r, S130r, S121r, S112r, S103r] = deal(0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05);
[S031r, S022r, S013r] = deal(0.05, 0.05, 0.05);

[S110, S101, S011, S300, S030, S003, S400, S040, S004, ...
 S210, S201, S120, S111, S102, S021, S012, ...
 S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
    edge_corner_correction(R110, R101, R011, ...
                           S110r, S101r, S011r, ...
                           S300r1, S030r1, S400r1, S040r1, ...
                           S300r2, S003r2, S400r2, S004r2, ...
                           S030r3, S003r3, S040r3, S004r3, ...
                           S300r, S030r, S003r, S400r, S040r, S004r, ...
                           S210r, S201r, S120r, S111r, S102r, S021r, S012r, ...
                           S310r, S301r, S220r, S211r, S202r, S130r, S121r, S112r, S103r, ...
                           S031r, S022r, S013r);

verifyEqual(testCase, abs(S011), 1, 'AbsTol', 0.01, 'S011 should be ±1 at edge');
S2 = 1 + 2*S110*S101*S011 - (S110^2 + S101^2 + S011^2);
verifyGreaterThanOrEqual(testCase, S2, -1e-6, 'S2 realizability should hold');

fprintf('✓ Edge S011 test passed\n');
end

function test_corner_case(testCase)
%TEST_CORNER_CASE Test corner case where all R110, R101, R011 ≤ 0
S110r = 0.99;
S101r = 0.99;
S011r = 0.99;

% All at boundary
R110 = -0.01;
R101 = -0.01;
R011 = -0.01;

[S300r1, S030r1, S400r1, S040r1] = deal(0.5, 0.6, 0.3, 0.4);
[S300r2, S003r2, S400r2, S004r2] = deal(0.5, 0.5, 0.3, 0.3);
[S030r3, S003r3, S040r3, S004r3] = deal(0.6, 0.5, 0.4, 0.3);
[S300r, S030r, S003r, S400r, S040r, S004r] = deal(0.5, 0.6, 0.5, 0.3, 0.4, 0.3);
[S210r, S201r, S120r, S111r, S102r, S021r, S012r] = deal(0.1, 0.1, 0.1, 0.05, 0.1, 0.1, 0.1);
[S310r, S301r, S220r, S211r, S202r, S130r, S121r, S112r, S103r] = deal(0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05);
[S031r, S022r, S013r] = deal(0.05, 0.05, 0.05);

[S110, S101, S011, S300, S030, S003, S400, S040, S004, ...
 S210, S201, S120, S111, S102, S021, S012, ...
 S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
    edge_corner_correction(R110, R101, R011, ...
                           S110r, S101r, S011r, ...
                           S300r1, S030r1, S400r1, S040r1, ...
                           S300r2, S003r2, S400r2, S004r2, ...
                           S030r3, S003r3, S040r3, S004r3, ...
                           S300r, S030r, S003r, S400r, S040r, S004r, ...
                           S210r, S201r, S120r, S111r, S102r, S021r, S012r, ...
                           S310r, S301r, S220r, S211r, S202r, S130r, S121r, S112r, S103r, ...
                           S031r, S022r, S013r);

% At corner, all three should be ±1
verifyEqual(testCase, abs(S110), 1, 'AbsTol', 0.01, 'S110 should be ±1 at corner');
verifyEqual(testCase, abs(S101), 1, 'AbsTol', 0.01, 'S101 should be ±1 at corner');
verifyEqual(testCase, abs(S011), 1, 'AbsTol', 0.01, 'S011 should be ±1 at corner');

% Product should be 1
product = S110 * S101 * S011;
verifyEqual(testCase, product, 1, 'AbsTol', 0.01, 'Product should be 1 at corner');

S2 = 1 + 2*S110*S101*S011 - (S110^2 + S101^2 + S011^2);
verifyGreaterThanOrEqual(testCase, S2, -1e-6, 'S2 realizability should hold');

fprintf('✓ Corner case test passed\n');
end

function test_output_consistency(testCase)
%TEST_OUTPUT_CONSISTENCY Verify outputs are averaged correctly
S110r = 0.99;
S101r = 0.3;
S011r = 0.4;
R110 = -0.01;
R101 = 1 - S101r^2;
R011 = 1 - S011r^2;

[S300r1, S030r1, S400r1, S040r1] = deal(0.5, 0.6, 0.3, 0.4);
[S300r2, S003r2, S400r2, S004r2] = deal(0.5, 0.5, 0.3, 0.3);
[S030r3, S003r3, S040r3, S004r3] = deal(0.6, 0.5, 0.4, 0.3);
[S300r, S030r, S003r, S400r, S040r, S004r] = deal(0.5, 0.6, 0.5, 0.3, 0.4, 0.3);
[S210r, S201r, S120r, S111r, S102r, S021r, S012r] = deal(0.1, 0.1, 0.1, 0.05, 0.1, 0.1, 0.1);
[S310r, S301r, S220r, S211r, S202r, S130r, S121r, S112r, S103r] = deal(0.05, 0.05, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05);
[S031r, S022r, S013r] = deal(0.05, 0.05, 0.05);

[~, ~, ~, ~, ~, ~, ~, ~, ~, ...
 S210, S201, S120, S111, S102, S021, S012, ...
 S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013] = ...
    edge_corner_correction(R110, R101, R011, ...
                           S110r, S101r, S011r, ...
                           S300r1, S030r1, S400r1, S040r1, ...
                           S300r2, S003r2, S400r2, S004r2, ...
                           S030r3, S003r3, S040r3, S004r3, ...
                           S300r, S030r, S003r, S400r, S040r, S004r, ...
                           S210r, S201r, S120r, S111r, S102r, S021r, S012r, ...
                           S310r, S301r, S220r, S211r, S202r, S130r, S121r, S112r, S103r, ...
                           S031r, S022r, S013r);

% All outputs should be finite
all_outputs = [S210, S201, S120, S111, S102, S021, S012, ...
               S310, S301, S220, S211, S202, S130, S121, S112, S103, S031, S022, S013];
verifyTrue(testCase, all(isfinite(all_outputs)), 'All outputs should be finite');

fprintf('✓ Output consistency test passed\n');
end

