function beam_steering_nb( b, theta, target)
%BEAM_STEERING_NB(b, theta, target) Calculate constrained nb beamformer
%
% Use beamsteering to obtain a desired beampattern. The beamformer has to
% fullfill the following linear constraints:
%  w' * A(theta) = target
% where:
%   w are the narrow band beamformer weights, stored as a column vector in
%   the property b.nb_weights.
%   A(theta) is the array response `matrix' containing the arrar response
%     vector for each theta
%   target is a vector containing the desired beampattern values for theta.
f = b.nb_frequency;
A = array_response_vector(b,theta, f);
b.nb_weights = A*(A'*A)^(-1)*target.';
end
