function [Px, Py] = choose_process_grid(nl)
% CHOOSE_PROCESS_GRID Choose nearly-square factorization for MPI process grid
%
% Syntax:
%   [Px, Py] = choose_process_grid(nl)
%
% Description:
%   Finds a factorization Px*Py=nl where Px and Py are as close as possible
%   to sqrt(nl). This creates a nearly-square 2D MPI process grid.
%
% Inputs:
%   nl - Total number of MPI processes
%
% Outputs:
%   Px - Number of processes in x-direction
%   Py - Number of processes in y-direction
%
% Example:
%   [Px, Py] = choose_process_grid(4)  % Returns Px=2, Py=2
%   [Px, Py] = choose_process_grid(6)  % Returns Px=2, Py=3

    bestDiff = inf;
    Px = 1; 
    Py = nl;
    
    for p = 1:nl
        if mod(nl, p) == 0
            q = nl / p;
            d = abs(p - q);
            if d < bestDiff
                bestDiff = d;
                Px = p;
                Py = q;
            end
        end
    end
end

