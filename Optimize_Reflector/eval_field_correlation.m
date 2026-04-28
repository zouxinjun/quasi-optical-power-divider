function [score, corr_coeff] = eval_field_correlation(H3_field)
%==========================================================================
% eval_field_correlation  —— Normalized correlation coefficient evaluation
%                            between the H3 field and the target pattern
%
% Function: Compute the normalized correlation coefficient between a given
%           H3 magnetic field distribution and the target field pattern,
%           and return its real part as the optimization fitness score.
%
% Principle:
%   After applying a reference phase offset to the H3 field, compute the
%   normalized cross-correlation with the weighted target pattern:
%
%     k = sum(H3_shifted * conj(H_target)) / sqrt(sum|H3|^2 * sum|H_target|^2)
%
%   Fitness score = 100 * Re(k), in the range ~[0, 100]; higher is better.
%
% Dependency files (must be in the same directory as this function):
%   Weight.txt      — weighting window function (1x153, 1 in center region, 0 at both ends)
%   Divider_C.txt   — target field pattern (1x153 complex vector)
%
% Input:
%   H3_field    — 1x153 complex vector, H-field distribution on H3 plane
%                 (computed by calc_H_field)
%
% Output:
%   score       — scalar, fitness score = 100 * Re(corr_coeff); higher is better
%   corr_coeff  — complex scalar, normalized cross-correlation coefficient
%==========================================================================

% Reference phase offset (degrees), used to compensate the fixed phase
% difference between the H3 field and the target pattern
phase_offset_deg = 100;

j_unit = sqrt(-1);

% Load the weighting window function (focus on the center region, suppress edge noise)
weight_window = dlmread('Weight.txt');

% Apply the reference phase offset to the H3 field
H3_shifted = H3_field * exp(j_unit * phase_offset_deg * pi / 180);

% Load the target field pattern and apply the weighting window
H_target = dlmread('Divider_C.txt') .* weight_window;

% Compute the normalized cross-correlation coefficient
numerator   = sum(H3_shifted .* conj(H_target));
denominator = sqrt(sum(abs(H3_shifted).^2) * sum(abs(H_target).^2));
corr_coeff  = numerator / denominator;

% Fitness score: take the real part of the correlation coefficient and
% scale by 100 to map it to the ~[0, 100] range
score = 100 * real(corr_coeff);

end
