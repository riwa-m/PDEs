clc;
clear all;  
close all; 

%% Example 

h=0.1;
k=0.0005;
%k=0.01; %If you TRY, the solution booms
L=1;
T=1;

m = L/h; 
n = T/k;

alpha=1;


%% Calling The User-Defined Code for the PDEs

[XX,TT,w,  lambda, ExactSolution, AbsoluteError] = Parabolic_Forward(L,T,m,n,alpha);

%% Plotting the Exact Solution
figure;
mesh(XX,TT, ExactSolution);
title('Exact Solution');
xlabel('x');
ylabel('t');
zlabel('u(x, t)');
colorbar;

% Plotting the Numerical Solution
figure;
mesh(XX,TT, w);
title('Numerical Solution');
xlabel('x');
ylabel('t');
zlabel('u(x, t)');
colorbar;

% Plotting the Absolute Error
figure;
mesh(XX,TT,AbsoluteError);
title('Absolute Error');
xlabel('x');
ylabel('t');
zlabel('u(x, t)');
colorbar;
 