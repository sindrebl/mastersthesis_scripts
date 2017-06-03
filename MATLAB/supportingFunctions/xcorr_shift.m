function [sample_shift] = xcorr_shift(shifted_signal,reference_signal)
% xcorr_shift returns the shift lag between a shifted signal and reference
% signal using xcorr() and sample_shift interpolates between the discrete 
% time steps.

[CorrValue, Lags] = xcorr(shifted_signal,reference_signal);
[~, C_maxidx] = max(CorrValue);

if C_maxidx == 0
    error('C_maxidx = 0. Are the rf samples filtered and windowed?')
end

% Interpolation to find maximum between time steps:
sample_shift = ( ...
                (...
                    (CorrValue(C_maxidx-1) - CorrValue(C_maxidx+1))...
                    /( 2*(CorrValue(C_maxidx-1) -2*CorrValue(C_maxidx) ...
                      +CorrValue(C_maxidx+1)) )...
                 ) + Lags(C_maxidx)...
                );
end