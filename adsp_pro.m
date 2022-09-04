clear all;
close all;
% clear classes; 
clc;
%% Array Configuration                                
Scenario = 'all';
switch Scenario
    case 'left' 
        dx = 0;
        dy = 0.01; 
        J = 2;
        my_array = arrays.ULA(J,dx,dy);
        disp('The left two microphones are ready !');
        figure(1);
        my_array.plot();  
        flag = 0;
    case 'upper'
        dx = 0.17;
        dy = 0; 
        J = 2;
        my_array = arrays.ULA(J,dx,dy);
        disp('The upper two microphones are ready !');
        figure(1);
        my_array.plot();  
        flag = 2;
    case 'all'
        dx = 0.17;
        dy = 0.01; 
        J = 4;
        my_array = arrays.square_array(J,dx,dy);
        disp('All four microphones are ready !');
        figure(1);
        my_array.plot();  
        flag = 1;
end


%% Beamformer
mt_f = [500 1000 2000 4000];            % multiple tones in Hz
c = 340;                                % the sound speed 
theta = [0,45];                         % row vector containing angles for which constraints hold
target = [1,0];                         % row vector containing target values for the beampattern

b = beamformer;
set(b, 'array',         my_array);
set(b, 'angles',        -180:1:180);
set(b, 'n_sensor',      J);
set(b, 'sound_speed',   c);
set(b, 'mt_frequency',  mt_f);

%% Pattern with Weights 1/J
figure(2)      % use narrow band to plot all frequcy multi-tone
for i = 1:length(mt_f)
nb_f = mt_f(i);
set(b, 'nb_frequency',  nb_f); 
b.nb_weights = 1/J*ones(J,1);
b.calc_nb_beampattern;
Beam = b.nb_beampattern;
plot(b.angles, 10*log10(Beam),'LineWidth',2);
%b.plot_nb; 
hold on;
end
legend('500Hz','1000Hz','2000Hz','4000Hz','Location','southwest')
axis tight
title('Beampattern in different frequency')
xlabel('Angle in [degrees]');
ylabel('Beamformer gain in [dB]');

b.mt_weights = 1/J*ones(J,length(mt_f)); 
b.calc_mt_beampattern;    %multi-tone
figure(3)
b.plot_mt;

for i = 1:length(mt_f)   %each frequcy to plot polar
nb_f = mt_f(i);
set(b, 'nb_frequency',  nb_f); 
b.nb_weights = 1/J*ones(J,1);
b.calc_nb_beampattern;
Beam = b.nb_beampattern;
figure(100+i)
polarplot(b.angles.'*pi/180,Beam,'LineWidth',2);
rlim([0 1])
title(['Polar plot in ',num2str(mt_f(i)),' Hz using ', Scenario, ' sensors']);
end


%% Null-steering pattern
figure(4)        % use narrow band to plot all frequcy multi-tone
for i = 1:length(mt_f)
nb_f = mt_f(i);
set(b, 'nb_frequency',  nb_f); 
b.beam_steering_nb(theta, target);
b.calc_nb_beampattern;
Beam = b.nb_beampattern;
plot(b.angles, 10*log10(Beam),'LineWidth',2);
hold on;
end
plot([theta;theta], [-350 max(10*log(Beam))]'* ones(1,length(theta)),'k--','LineWidth',2);
legend('500Hz','1000Hz','2000Hz','4000Hz','Location','southwest')
axis tight
title('Beampattern in different frequency')
xlabel('Angle in [degrees]');
ylabel('Beamformer gain in [dB]');


b.beam_steering_mt(theta, target);
b.calc_mt_beampattern;     %multi-tone
figure(5);
b.plot_mt(theta, {'k-.','LineWidth',2});

for i = 1:length(mt_f)
nb_f = mt_f(i);        
set(b, 'nb_frequency',  nb_f); 
b.beam_steering_nb(theta, target);
b.calc_nb_beampattern;
Beam = b.nb_beampattern;
figure(200+i)
polarplot(b.angles.'*pi/180,Beam,'LineWidth',2);    %each frequcy to plot polar
rlim([0 1])
title(['Null-steering Polar plot in ',num2str(mt_f(i)),' Hz using ', Scenario, ' sensors under ',num2str(theta(2)),' deg interference']);
end


%% Fractional Delay
if flag ~= 2           %upper sensors donot need delay
Dataset= load('DatasetAssignBs2.mat');
fs = Dataset.fs;
time_align = dy/c * fs;    %desired source delay
tau = time_align;
%tau = [3,5.15];
%tau = [15.551,42.443];
disp(['Time align = ',num2str(tau),' samples !']);

legend_the = cell(1,length(tau));   %ledend can not accumulate in a loop
legend_fra = cell(1,length(tau));   %so to set a cell
for i = 1:length(tau)
[h,my_h] = Fractional_delay(tau(i));
my_h(isnan(my_h)) = 1;
legend_the{i} = ['Theoretical ' num2str(tau(i)) ' samples'];
legend_fra{i} = ['Fractional ' num2str(tau(i)) ' samples'];
figure(50)
plot(my_h,'LineWidth',1)
hold on;
stem(h,'LineWidth',1)
hold on;
figure(51)
plot(my_h(1:end-183),'LineWidth',1)  
hold on;
%stem(circshift(h,length(h)-183),'LineWidth',1)
stem(h(184:end),'LineWidth',1)    %align by cutting first 183 points
hold on;
end
%legend({legend_the});
if length(tau) == 1
    figure(50)
    legend([legend_the, legend_fra]);
    figure(51)
    legend([legend_the, legend_fra]);
else
figure(50)
legend(legend_the{1}, legend_fra{1},legend_the{length(tau)},legend_fra{length(tau)});
figure(51)
legend(legend_the{1}, legend_fra{1},legend_the{length(tau)},legend_fra{length(tau)});
end
%h = circshift(h,length(h)-183);   %using circshift to do align
h = h(184:end);              % using cut to do align

%% Time Align
x = Dataset.x(1:J,:);
u = x;
u(1,:) = filter(h,1,x(1,:).').';                %align sensor 1
%u(1,:) = filter(my_h,1,x(1,:).').';
if flag == 1                                    %when using four sensors 
    u(J,:) = filter(h,1,x(J,:).').';            %align sensor 4
    %u(J,:) = filter(my_h,1,x(J,:).').';   
end
v0 = 1/J * sum(u);                              %spatial delay and sum
d = v0;
d(:,1:end-1) = v0(:,2:end);                    %integer delay for causality
d(:,end) = v0(:,1);    % d(:,end) = 0;

%d(:,1:end-183) = v0(:,184:end);               %integer delay for causality
%d(:,end-182:end) = v0(:,1:183);    


%v0 = (v0 - sum(v0)./length(v0))./std(v0);

if 0
soundsc(v0./max(v0),fs);
disp('Listen,the upper branch audio is playing !');
end

%audiowrite('GSC_upper_branch.wav',v0,fs);
%v = u(1,:) - u(J,:);
%soundsc(v,fs);
%B = [-1 1];
B(:,1) = -1*ones(J-1,1);                      %generate blocking matrix
B(:,2:J) = eye(J-1);
v = B*u;

%% Adaptive
%w_plot = [];
state = Init_A1s2(J);
w_plot = zeros(state.nw,J-1,length(d));       %to plot the weights
for i = 1:length(d)                           %the number of updates         
    state = Update_A1s2(state,d(i),v(:,i).');
    %w_plot = [w_plot,state.w]; 
    w_plot(:,:,i) = state.w;
end 
w_a = state.w;
g_temp = zeros(J-1,length(d));
for i = 1:J-1
    g_temp(i,:) = filter(w_a(:,i),1,v(i,:));    
end
%g = filter(w_a,1,v);
g = 1/(J-1) * sum(g_temp,1);
y = d - g;

soundsc(y./max(y),fs);
disp('Listen,the output audio is playing !');

figure
%for i = 1:state.nw
%    plot(w_plot(i,:));
%    hold on;
%end
for i = 1:J-1
subplot(J-1,1,i);
    for j = 1:state.nw
        plot(squeeze(w_plot(j,i,:)),'LineWidth',1.5);
        hold on;
    end
title(['Sensor ' num2str(i+1)]);
end

else
   disp('Warning!! There is no delay!');
end


%% Trans function
theta = 0;                      % desired source angle
n_fft = 256;
mt_f = linspace(1,fs/2,n_fft/2);

set(b, 'mt_frequency',  mt_f);

Y = fft(y,n_fft);
X = fft(x.',n_fft).';
Trans = Y./X;
b.mt_weights = Trans(:,1:n_fft/2);

%for i = 1:J
%    b.mt_weights(i,:) = fft(y)./fft(x(i,:));
%end
%b.mt_weights = b.mt_weights(:,1:end/2);
%for i = 1:J
%    b.mt_weights(i,:) = ifft(fft(y)./fft(x(i,:)));
%end
%for i =1:J
%b.mt_weights(i,:) = y./x(i,:);
%end

b.calc_mt_beampattern;
figure(1000);
b.plot_mt(theta, {'k-.','LineWidth',2});


%% Function
function [h,my_h] = Fractional_delay(tau)
dmax = ceil(tau);
d_int = fix(tau);
index_dot = find(num2str(tau)== '.');
if isempty(index_dot)
    L = 10;
    fd = 0;
else
    L =10^(length(num2str(tau)) - index_dot);
    fd = fix((tau - d_int)*L);    
end
h = delay(d_int,fd,L,dmax);
%H =10*log10(abs(fft(h)));
%% Analytical expression
len_filter = length(h);
%half = fix((len_filter-1)/2);
%n = linspace(-half,half,len_filter);
n = 1:len_filter;
%n = 1:len_filter-183;
my_h = sin(pi*(n - tau))./(pi*(n - tau));
end

function [ state ] = Init_A1s2(J)
% Initialize all variables that you want to use: 
deta = 0.1;        %0.1
state.J = J;
state.alpha =0.00008;  %0.0008    
state.nw =5;                            %The smaller the better£¬but the running time becames bigger
state.var = 0.5*ones(1,state.J-1); 
state.w = zeros(state.nw,state.J-1);  
state.x = zeros(state.nw,state.J-1);    %size 5*3 which means Time*Spatial
state.R_inv = deta * eye(state.nw);
end

function state = Update_A1s2( state, x1, x2 )
alpha = state.alpha;
e = x1;  
state.x(2:end,:) = state.x(1:end-1,:);     %
state.x(1,:) = x2;                         %each update only transfer a row
%r = e - state.w.' * state.x;
r = e - (1/(state.J-1)) * sum(diag(state.w.' * state.x));

rls;
%nlms;
%% RLS
function rls
lamda = 0.9999999;             %deta = 0.1
%lamda = 1;
gain = (state.R_inv * state.x)./(lamda^2 + diag(state.x.' * state.R_inv * state.x)).';
state.R_inv = (state.R_inv - gain * state.x.' * state.R_inv)./(lamda^2);
state.w = state.w + r * gain; 
end

%% NLMS
function nlms 
beta = 0.9999999;  %0.5
%beta = 1;
state.var = beta*state.var + (1-beta)*((diag(state.x.'*state.x)).'./state.nw); 
state.w = state.w + 2*alpha*r*state.x./state.var;
end
end























