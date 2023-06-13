clc;
clear all;  
close all; 

%% Example  INPUT 

h=0.01;
k=0.04; 
L=1;
T=1;

m = L/h; 
n = T/k;

alpha=1/pi; 


%% Calling The User-Defined Code for the PDEs

[XX,TT,w,  lambda, ExactSolution, AbsoluteError] = Parabolic_Backward(L,T,m,n,alpha);
%[XX,TT,w,  lambda] = Parabolic_Backward(L,T,m,n,alpha);
% %% Plotting the Exact Solution
% figure;
% mesh(XX,TT, ExactSolution);
% title('Exact Solution');
% xlabel('x');
% ylabel('t');
% zlabel('u(x, t)');
% colorbar;
% 
% % Plotting the Numerical Solution
% figure;
% mesh(XX,TT, w);
% title('Numerical Solution');
% xlabel('x');
% ylabel('t');
% zlabel('u(x, t)');
% colorbar;
% 
% % Plotting the Absolute Error
% figure;
% mesh(XX,TT,AbsoluteError);
% title('Absolute Error');
% xlabel('x');
% ylabel('t');
% zlabel('u(x, t)');
% colorbar;
% 