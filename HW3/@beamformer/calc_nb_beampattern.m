function calc_nb_beampattern(b)
%CALC_NB_BEAMPATTERN(b) Calculate the beamformer coefficients
%
% Calculate the beampattern coefficient for each angle that is stored in
% b.angles.
%
% The narrowband beamformer weights are available through b.nb_weights
%
% If you need array response vectors you can calculate them using the
% function: b.array_response_vector(theta, f)
%
% Store the resulting beampattern in b.nb_beampattern.
J = b.n_sensor;
theta =b.angles;    %  size 1*361
f = b.nb_frequency;
w = b.nb_weights;      % size 4*1
A = array_response_vector(b,theta, f);  % size 4*361
%B_theta = abs(A' * w).^2 ./(J^2);    % size 361*1
B_theta = abs(A' * w).^2;
% B_theta = abs(A.' * conj(w)).^2;
b.nb_beampattern = B_theta;
end

