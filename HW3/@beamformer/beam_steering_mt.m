function beam_steering_mt(b, theta, target)
%BEAM_STEERING_MT(b, theta, target) Calculate constrained mt beamformer
%
% Use beamsteering to obtain a desired beampattern for a multi-tone
% beamformer.
%
% The beamformer has to fullfill the following linear constraints for each
% frequency:
%  w' * A(theta) = target
% where:
%   w are the multi-tone beamformer weights that should be stored in the
%   property b.mt_weights with for each frequency a column.
%   A(theta) is the array response `matrix' containing the array response
%     vector for each theta
%   target is a vector containing the desired beampattern values for theta.
f = b.mt_frequency;
n_mt_fre = size(f,2);
for i = 1:n_mt_fre
    A = array_response_vector(b,theta, f(i));
    b.mt_weights(:,i) = A*(A'*A)^(-1)*target.';
end
end

