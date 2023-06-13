%% PARABOLIC PDE- POISSON EQUATION

% To find an approximate the solution of the considered Equation
% u_t - alpha^2 u_xx = 0
%
% Conditions at x values
% u(a,y)=u_x_LEFT(x,t),     that is, the value of the u(x,t) at the x=a (LEFT endpoint)
% u(b,y)=u_x_RIGHT(x,t)     that is, the value of the u(x,t) at the x=b (RIGHT endpoint)
%
% Initial Condition 
% u_t_0,       it is equal to u(x,0)  , that is, the value of the u(x,t) at the t=0 
%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: 
%        Endpoints a,b on x-axis; 
%        Endpoints c,d on t-axis;
%        N: the number of the intervals on x-axis (h = (b-a)/N);
%        M: the number of the intervals on t-axis (k = (d-c)/M); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT: 
%         Approximation w (numerical solution) to y(x(i)) (Exact solution);  
%         for each i=0,1,...,N.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% function   [OUTPUTS]  = FUNCTIONNAME(INPUTS)  
    
function [XX,TT,w,  lambda, ExactSolution, AbsoluteError] = Parabolic_Backward(L,T,m,n,alpha)
%function [XX,TT,w,  lambda] = Parabolic_Backward(L,T,m,n,alpha)
 %% EmBEDDED INPUTS 
%       u_x_LEFT = @(x,t) .....;  % PLEASE ENTER the value of left endpoint, that is, u(a,y)=u_x_LEFT(x,t);
%       u_x_RIGHT = @(x,t) .....;  % PLEASE ENTER the value of right endpoint, that is, u(b,y)=u_x_RIGHT(x,t)
%       u_t_0= @(x,t) .....;  % PLEASE ENTER u(x,0)  

%       You should also enter the Exact Solution. Please see at the end of this Script.

%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EX1- c
u_x_LEFT = @(x,t) 0;
u_x_RIGHT = @(x,t) 0; 
u_t_0 = @(x,t) cos(pi*(x-0.5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% MAIN CODE
h = L/m; 
k = T/n; 

x = 0:h:L;
t = 0:k:T;

lambda = (alpha^2)*k/(h^2);
mu = (1+2*lambda); 
% Bundary conditions

for j = 1:n+1
    w(1,j) = u_x_LEFT(x(1),t(j));
    w(m+1,j) = u_x_RIGHT(x(m+1),t(j));
end


% Initial condition
for i = 1:m+1
    w(i,1) = u_t_0(x(i),t(1));
end 


% Coefficient Matrix

for i = 1:m-1
    K(i,i) = mu;
end

for i = 1:m-2
    K(i,i+1) = -lambda;
    K(i+1,i) = -lambda;
end

% Algebraic equation system  

for j = 2:n+1
    
    F(1,1) = w(2,j-1) + lambda*w(1,j);
    F(m-1,1) = w(m,j-1) + lambda*w(m+1,j);
    
    for i = 2:m-2
        F(i,1) = w(i+1,j-1);
    end
    
    ww = K\F; 
    
    for i = 2:m
        w(i,j) = ww(i-1,1);
    end
    
end


[TT,XX] = meshgrid(t,x); 

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% ExactSolution = ......; % Please enter Exact Solution 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXAMPLE   
ExactSolution = exp(-TT).*cos(pi*(XX-0.5)); 
 
%%    
AbsoluteError = abs(ExactSolution-w);

%%
% Creating the table
tableData = zeros( m+1,n+1, 7);  
%tableData = zeros( m+1,n+1, 5); 

for i = 1:m+1
    for j = 1:n+1
        tableData(i, j, 1) = i;                       % i
        tableData(i, j, 2) = j;                       % j
        tableData(i, j, 3) = XX(i, j);                   % X(i)
        tableData(i, j, 4) = TT(i, j);                   % Y(j)
        tableData(i, j, 5) = w(i, j);                   % w(i, j)
        tableData(i, j, 6) = ExactSolution(i, j);       % ExactSolution(i, j)
        tableData(i, j, 7) = AbsoluteError(i, j);       % AbsoluteError(i, j)
    end
end

% Creating the table format
tableFormat = {'i', 'j', 'X(i)', 'T(j)', 'w(i, j)', 'ExactSolution(i, j)', 'AbsoluteError(i, j)'};
%tableFormat = {'i', 'j', 'X(i)', 'T(j)', 'w(i, j)'};

% ReShaping the table
dataTable = array2table(reshape(tableData, [], 7), 'VariableNames', tableFormat);
%dataTable = array2table(reshape(tableData, [], 5), 'VariableNames', tableFormat);


% Displaying the table
disp(dataTable);