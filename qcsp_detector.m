function [maxScore, freqOffset, allScores, frameDetected, freqIndex, maxFound, chipFromMax, chipSinceLastDet] ...
		= qcsp_detector(chipBuffer, pn, N, pOmega, threshold, stepDenominator, stepNumerator, normed)
% QCSP_DETECTOR compute the time-sliding QCSP score of chipBuffer related to pn for a frame of N symbols.
%
%   The score is computed using time sliding correlations accumulated for N symbols.
%   Every chip of chipBuffer are processed, thus resulting in a maxScore vector such that:
%
%              length(maxScore) == length(chipBuffer).
%
%   pOmega rotation hypotheses are computed, each separated by stepNumerator / stepDenominator * pi.
%   Hypotheses are centered on 0, thus when pOmega is odd, hypothese {pOmega - 1 / 2} is without rotation.
%
%   Example: stepDenominator = 2, stepNumerator = 1 <=> step = pi / 2
%            if pOmega = 3: hypotheses = [ -pi / 2, 0, pi / 2 ]
%            if pOmega = 2: hypotheses = [ -pi / 4  ,  pi / 4 ]
%
%   If normed is 'true', data are normalized using a length(pn)-L2 sliding norm.
%   No normalization is done otherwise.
%
%   OUTPUTS:
%     - maxScore contains the maximum score for the best frequency hypotheses.
%       If pOmega = 1, it is undefined, and the user should use allScores instead.
%     - freqOffset contains the best frequency hypothese. If pOmega = 1, it is undefined.
%     - allScores is assured to contain all the pOmega * (length(chipBuffer) - length(pn) + 1) scores,
%       in a matrix of size [ pOmega, length(chipBuffer) ].
%     - frameDetected is a logical vector, set at 'true' when a frame is detectd, and at 'false' otherwise.
%     - freqIndex is a uint32 vector. Each element corresponds to the rotation index in the rotation vector
%       [ -int(pOmega / 2) * step, -int(pOmega / 2) * step + step,  ... , int(pOmega / 2) * step ].
%     - maxFound is a logical vector, set at 'true' when a max is found, 'false' otherwise.
%     - chipFromMax is an internal metric that represent the discrete time between detection time and the max.
%     - chipSinceLastDet is an internal metric representing the discrete time since the last detection
%
