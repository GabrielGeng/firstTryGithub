clear all;
close all;
clear classes; 
%clc;
%% Array Configuration                           
Scenario = 'all';
switch Scenario
    case 'left' 
        dx = 0;
        dy = 0.01; 
        J = 2;
        my_array = arrays.ULA(J,dx,dy);
        disp('The left two microphones are ready !');
        figure;
        my_array.plot();  
        flag = 0;
    case 'upper'
        dx = 0.17;
        dy = 0; 
        J = 2;
        my_array = arrays.ULA(J,dx,dy);
        disp('The upper two microphones are ready !');
        figure;
        my_array.plot();  
        flag = 0;
    case 'all'
        dx = 0.17;
        dy = 0.01; 
        J = 4;
        my_array = arrays.square_array(J,dx,dy);
        disp('All four microphones are ready !');
        figure;
        my_array.plot();  
        flag = 1;
end

mt_f = [500 1000 2000 4000];            % multiple tones in Hz
c = 340;                                % the sound speed 
%% Fractional 
Dataset= load('DatasetAssignBs2.mat');
fs = Dataset.fs;
time_align = dy/c * fs;
tau = time_align;
%tau = 5.15
dmax = ceil(tau);
d = fix(tau);
fd = fix((tau - d)*10^(length(num2str(tau)) - find( num2str(tau)== '.')));
L =10^(length(num2str(tau)) - find( num2str(tau)== '.'));
h = delay(d,fd,L,dmax);
%H =10*log10(abs(fft(h)));

%% Scenario 2:
x = Dataset.x(1:J,:);
u = x;
u(1,:) = filter(h,1,x(1,:).').';
if flag
    u(J,:) = filter(h,1,x(J,:).').';
end
v0 = 1/J * sum(u);
d = v0;
d(:,1:end-1) = v0(:,2:end);
d(:,end) = v0(:,1);    % d(:,end) = 0;
%v0 = (v0 - sum(v0)./length(v0))./std(v0);
soundsc(v0,fs);
%audiowrite('GSC_upper_branch.wav',v0,fs);
disp('Listen,the upper branch audio is playing !');
%v = u(1,:) - u(J,:);
%soundsc(v,fs);
%B = [-1 1];
B(:,1) = -1*ones(J-1,1);
B(:,2:J) = eye(J-1);
v = B*u;

%% Adaptive
%w_plot = [];
state = Init_A1s2(J);
w_plot = zeros(state.nw,J-1,length(d));
for i = 1:length(d)
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
soundsc(y,fs);
disp('Listen,the output audio is playing !');
figure
%for i = 1:state.nw
%    plot(w_plot(i,:));
%    hold on;
%end
for i = 1:J-1
subplot(J-1,1,i);
    for j = 1:state.nw
        plot(squeeze(w_plot(j,i,:)));
        hold on;
    end
title(['Sensor ' num2str(i+1)]);
end

function [ state ] = Init_A1s2(J)
% Initialize all variables that you want to use: 
deta = 0.1;        %0.1
state.J = J;
state.alpha =0.0008;  %0.0008                                                             
state.nw =5;    %5           % The smaller the better��but the running time becames bigger and the plot only display 5sec
state.w = zeros(state.nw,state.J-1);  
state.x = zeros(state.nw,state.J-1);
state.R_inv = deta * eye(state.nw);
end

function state = Update_A1s2( state, x1, x2 )
alpha = state.alpha;
e = x1;  
state.x(2:end,:) = state.x(1:end-1,:);     % this is a vector
state.x(1,:) = x2; 
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
beta = 0.5;  %0.5
state.var = beta*state.var + (1-beta)*((state.x.'*state.x)./state.nw); 
state.w = state.w + (2*alpha./state.var)*r*state.x;
end
end