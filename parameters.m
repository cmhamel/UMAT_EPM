time_end = 1000

n_steps = 100;

delta_t = time_end / n_steps

q = 1
plot_var = 'k-';

Gshear_eq = 1.0
Kbulk_eq = 1e3*Gshear_eq
%
nVisco = 1

%Gshear_neq = [5.0, 10.0, 15.0, 20.0];
%Gshear_neq = [15.0 10.0 5.0]
Gshear_neq = [1.0];
Kbulk_neq = 1e3*Gshear_neq;
%tau = [500.0, 500.0, 500.0, 500.0];
%tau = [500.0e-12 1500.0e-12 3000.0e-12]
tau = [10];
% for uniaxial tension case
%
lambda_max = 4.0;

gamma_max = 1.0;

delta_lambda = lambda_max^(1.0/(n_steps));

delta_gamma = 1.0 / n_steps;

F0 = eye(3);
deltaF = [delta_lambda 0 0; ...
          0 1/sqrt(delta_lambda) 0; ...
          0 0 1/sqrt(delta_lambda)];
%       
% deltaF = [1 delta_gamma 0;...
%           0 1 0;...
%           0 0 1]
      
% deltaF = [1 0 0;...
%           0 1 0;...
%           0 0 1];
% 
% F0 = [4 0 0;...
%       0 0.5 0;...
%       0 0 .5];