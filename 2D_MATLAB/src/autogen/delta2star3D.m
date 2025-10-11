function E = delta2star3D(s300,s400,s110,s210,s310,s120,s220,s030,s130,s040,s101,s201,s301,s102,s202,s003,s103,s004,s011,s111,s211,s021,s121,s031,s012,s112,s013,s022)
%delta2star3D - Wrapper that uses MEX implementation if available, falls back to MATLAB
%    E = delta2star3D(S300,S400,S110,S210,S310,S120,S220,S030,S130,S040,S101,S201,S301,S102,S202,S003,S103,S004,S011,S111,S211,S021,S121,S031,S012,S112,S013,S022)
%
%    Performance: MEX version is ~10-100x faster than pure MATLAB
%    Falls back to optimized MATLAB version if MEX is not available

persistent use_mex;

if isempty(use_mex)
    % Check if MEX version exists (only once at first call)
    use_mex = exist('delta2star3D_mex', 'file') == 3;
    if use_mex
        fprintf('[delta2star3D] Using fast MEX implementation\n');
    else
        fprintf('[delta2star3D] MEX not found, using MATLAB fallback\n');
    end
end

if use_mex
    % Use fast MEX implementation
    E = delta2star3D_mex(s300,s400,s110,s210,s310,s120,s220,s030,s130,s040,s101,s201,s301,s102,s202,s003,s103,s004,s011,s111,s211,s021,s121,s031,s012,s112,s013,s022);
else
    % Fall back to optimized MATLAB version
    E = delta2star3D_matlab(s300,s400,s110,s210,s310,s120,s220,s030,s130,s040,s101,s201,s301,s102,s202,s003,s103,s004,s011,s111,s211,s021,s121,s031,s012,s112,s013,s022);
end

end

