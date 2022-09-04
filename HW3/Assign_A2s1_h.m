%% Assignment A2 Scenario 1: Narrowband beamformers
%
%   Implement tools for narrowband beamformer design and evaluation
%   Result plots: sensor locations, beampatterns in plot format.
%   Note:   Each step uses different functions of the beamformermer class.
%           These functions are denoted with each step.
%

clear all;
close all;
clear classes; 
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign A2s1 b)
% 1. Setup the ULA settings as in figure 3.1 of the assignment

J = 4;                  % Number of sensors
dy = 0.034;                 % meters of element spacing in x-direction
dx = 0;                 % meters of element spacing in y-direction
nb_f = 2500;               % narrowband (nb) frequency in Hz
c = 340;                   % the sound speed 

% Setup an ULA array from the settings and plot the array configuration
my_array = arrays.ULA(J,dx,dy);

% Use the plot function that belongs to the array class (in @array
% directory)
figure;
my_array.plot();

% Create a beamformer object and put settings in the beamformer object.
b = beamformer;
set(b, 'array',         my_array);
set(b, 'angles',        -180:1:180);
set(b, 'nb_frequency',  nb_f);
set(b, 'n_sensor',      J);
set(b, 'sound_speed',   c);

% Display all properties of the beamformer b:
%b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign A2s1 c)
% 1. Implement the array_response_vector.m method that is located in the
%    @beamformer folder.
% 2. Verify the result of the matlab function with your answer at a)

% Remove this return to continu with the assignment
%return;

theta_d = 30;
A = b.array_response_vector(theta_d, nb_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign A2s1 d)
% plot my hand calculate result with weights = 1/J
theta = (b.angles*pi/180).';
%my_pattern = 1/J^4 * abs(sin(sin(theta)*J*pi/4)./sin(sin(theta)*pi/4)).^2;  % size 361*1
my_pattern = 1/J^4 * abs(sin(cos(theta)*J*pi/4)./sin(cos(theta)*pi/4)).^2;   % change sensor position, so the projection is cos(theta)
B = 10*log10(my_pattern);
linspec = {'b-','LineWidth',2};
figure
plot(b.angles, B, linspec{:});
axis tight
title('My Beampattern by hand-calculated')
xlabel('Angle in [degrees]');
ylabel('Beamformer gain in [dB]');
%b.nb_beampattern = my_pattern;
%figure;
%b.plot_nb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign A2s1 e)
% 1. Implement the calc_nb_beampattern.m method that is located in the
%    @beamformer folder.
% 2. Verify the result of the matlab function with your answer at d)

% Remove this return to continu with the assignment
%return;
 
% Set the beamformer weights to 1/J
%b.nb_weights = ones(J,1);
b.nb_weights = 1/J*ones(J,1);
b.calc_nb_beampattern;

% Use the plot function that belongs to the narrowband beamformer
figure;
b.plot_nb;  % the weight is 1/J ,so when theta = 0, we don't get the 0dB gain.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign A2s1 f)
% 1. Implement the beam_steering_nb.m method that is located in the
%    @beamformer folder.
% 2. Verify the result of the matlab function by visual inspection.

% Remove this return to continu with the assignment
%return;

theta = [30]; % row vector containing angles for which constraints hold
target = [1]; % row vector containing target values for the beampattern
b.beam_steering_nb(theta, target);
b.calc_nb_beampattern;

% Use the plot function that belongs to the narrowband beamformer
figure;
b.plot_nb([],theta, {'k-.','LineWidth',2});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign A2s1 g)
% 1. Add the undesired source direction and make sure that the beamformer
% has unity response at 30 degrees and a zero response at -50 degrees.

% Remove this return to continu with the assignment
%return;

theta = [30,-50]; % row vector containing angles for which constraints hold
target = [1,0]; % row vector containing target values for the beampattern
b.beam_steering_nb(theta, target);
b.calc_nb_beampattern;

% Use the plot function that belongs to the narrowband beamformer
figure;
b.plot_nb([],theta, {'k-.','LineWidth',2});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













