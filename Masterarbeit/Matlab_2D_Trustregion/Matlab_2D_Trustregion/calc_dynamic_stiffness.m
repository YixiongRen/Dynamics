function S = calc_dynamic_stiffness(freq,K,M,D)
% Calculates the dynamic stiffness matrix
%==========================================================================
% Inputs
%--------------------------------------------------------------------------
% freq  = frequency
% K     = Stiffness matrix
% M     = Mass matrix
%==========================================================================
% Outputs
%--------------------------------------------------------------------------
% S     = Dynamic stiffness matrix
%==========================================================================

omega = 2*pi*freq;
S = K - omega^2*M + omega*D;

end

