function calc_mt_beampattern(b)
%CALC_MT_BEAMPATTERN(b) Calculate a multitone beampattern
%
% Calculate the beampattern for each frequency in b.mt_frequency
% The matrix b.mt_weights contains the weigth vector per frequency as
% column vectors, i.e.,
%   b.mt_weights(:,f1)
%   b.mt_weights(:,f2)
%
% The output b.mt_beampattern should contain the beampattern per frequency
% as row vectors, i.e.,
%   b.mt_beampattern(f1,:) = ...
%   b.mt_beampattern(f2,:) = ...
%
J = b.n_sensor;
f = b.mt_frequency;  %size 1*50
n_mt_fre = size(f,2);
%b.mt_weights = ones(J,n_mt_fre);
theta =b.angles;    %  size 1*361
w = b.mt_weights;      % size 4*50
b.mt_beampattern = zeros(n_mt_fre,size(theta,2));   %size 50*361
for i = 1:n_mt_fre
    A = array_response_vector(b,theta, f(i));  % size 4*361
    %B_theta = abs(A.' * conj(w)).^2 ./(J^2);    % size 361*1
    B_theta = abs(A' * w(:,i)).^2 ./(J^2);    % size 361*1
    b.mt_beampattern(i,:) = B_theta.';        
end
end

