%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RBC Model from BMR.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define some parameters 
beta = 0.99;
alpha = 0.33;
delta = 0.015;
eta = 1;
rho = 0.95;
sigma = 1;


%% Steady State values
% standard RBC Model has 5 steady state values
% SS ratios useful for undetermined coefficients solution
Rstar = 1/beta;
IYstar = (alpha*delta)/(beta^-1 -1 + delta);
IKstar = delta;
CYstar = 1- (alpha*delta)/((beta^-1) - 1+delta);
YKstar = (beta^-1 -1 + delta)/alpha;
%%
% we want are log linarised equations in the following form

%1: 0  = AXt + BXt-1 + CYt + Dzt
%2: 0 = Et[FXt+1 + GXt + HXt-1 + JYt+1 + KYt + Lzt+1 + Mzt]
%3: zt = Nzt-1 + et

% Where Xt = [kt]; Yt = [ct yt nt rt it]' ; zt = [at]

A = [0 0 1 0 0]';
B = [0 -alpha -(1-delta) 0 (alpha/Rstar)*(YKstar)]';
C = [-CYstar 0 0 eta 0;
    1 1 0 -1 (alpha/Rstar)*(YKstar);
     0 -(1-alpha) 0 1 0;
     0 0 0 0 1;
     -IYstar 0 -IKstar 0 0]';
 
 D = [0 -1 0 0 0]';
 
 
 % second block of equations
 % on 0 = ct  -Etct+1 + 1/eta Ert+1
 F = [0];
 G = [0];
 H = [0];
 L = [0];
 M = [0];
 
 J = [-1 0 0 1/eta 0];
 K = [1 0 0 0 0];
 N = [rho];
 
 % number of endogenous states
 % uhlig code m is the number of cols in A matrix
 m = size(A,2);
[l_equ ,n_endog]= size(C);
k_exog  = min(size(N));      % Technology variable
 %% Solve using Uhlig's method
 % different algotrithm depending on number of variables
 if l_equ == 0
     phi = F;
     gamma = -G;
     theta = -H; 
 Xi = [gamma, theta;
       eye(m), zeros(m)];
 Delta = [phi, zeros(m)
           zeros(m), eye(m)];
 else
 C_0 = (null(C'))';
 
 phi = [zeros(l_equ - n_endog, m)
        F - J*pinv(C) * A];
 gamma = [ C_0 * A
        J*pinv(C)*B - G + K*pinv(C)*A];
    
 theta = [C_0 * B
     K*pinv(C)*B-H];
 
 Xi = [gamma , theta
        eye(m), zeros(m)];
    
 Delta = [phi, zeros(m)
          zeros(m), eye(m)];
 end;

%% decomposition
% Omega is eigenvector, Lamdba is eigenvalues
[Omega, Lambda] = eig(Xi,Delta, 'qz');
Lambda2 = diag(Lambda);
bignum = 1e+07;

%% need to deal with infinite values that eig returns
eigvalues = abs(Lambda);
eigreal = real(Lambda);
%if isinf(eigvalues)
    [row,col] = find(isnan(eigvalues));
    Lambda(row,col) = bignum;
    
%% sort eigenvalues // Add in if statements for imaginary numbers
[Xi_abs, Xi_index] = sort(abs(diag(Lambda)));
% puts in diagnol
Xi_sortvec = Omega(1:2*m, Xi_index);    
Xi_sortval = diag(Lambda(Xi_index, Xi_index));
Xi_select = 1:m;
%if imag(Xi_sortval(m) - conj(Xi_sortval(m+1))) < Tol

Lambda_mat = diag(Xi_sortval(Xi_select));
Omega_mat = [Xi_sortvec((m+1):(2*m),Xi_select)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is the matrix we wanted to find.
P  = Omega_mat*Lambda_mat*pinv(Omega_mat);

% We can now solve for Q,R,S using P
R = -pinv(C) * (A*P+B);
V = [kron(eye(k_exog), A), kron(eye(k_exog),C)
     kron(N',F) + kron(eye(k_exog),... % on the same row dimension
     F*P+J*R+G),kron(N',J)+kron(eye(k_exog), K)];

LN_plus_MM = L*N+M; 
QS_vec = -V \ [D(:)
               LN_plus_MM(:)];
           
% same as inv(-V)* [D(:) ; LN_plus_MM(:)];, just computed differently

% Q is part of the control variable capital
Q = reshape(QS_vec(1:m*k_exog), m, k_exog);

% S is part of the endogenous variables
S= reshape(QS_vec((m*k_exog+1):((m+n_endog * k_exog))), n_endog, k_exog);

W = [eye(m), zeros(m,k_exog)
    R*pinv(P), (S-R*pinv(P)*Q)
    zeros(k_exog,m), eye(k_exog)];


% slightly different results to BMR package

%% Build State Space structure Xt = FFXt-1 + JJet

FF = [P, zeros(m,n_endog), Q*N;
      R, zeros(n_endog, n_endog), S*N
      zeros(k_exog, m), zeros(k_exog, n_endog), eye(k_exog)];

JJ = [zeros(m,m+n_endog), Q;
       zeros(n_endog, m+n_endog), S;
       zeros(k_exog, m+n_endog), eye(k_exog)];
   
%% can produce IRF with this form.
FF_Lag = [P, zeros(m,n_endog), zeros(m,k_exog)
          R, zeros(n_endog,n_endog), zeros(n_endog,k_exog)
          zeros(k_exog, m+n_endog), N];



JJ_contemp = eye(m+n_endog+k_exog) + ...
    [zeros(m,m+n_endog), Q;
       zeros(n_endog, m+n_endog), S;
       zeros(k_exog, m+n_endog), zeros(k_exog,k_exog)];

Init_date = 0;
Horizon = 60;
Period = 12; 
response = zeros(m + n_endog+k_exog, Horizon);
response(m+n_endog+k_exog,1) = 1;           % Tech shock
Time_axis = (-Init_date: (Horizon-1))/Period;
response_mat = [];

response(:,1) = JJ_contemp*response(:,1);

%% IRF Loop to tech shock
for time_counter = 2:Horizon
    % shock to tech is first col last row
    response(:,time_counter) = JJ_contemp * FF_Lag * response(:,time_counter-1);
end;
response_mat = [response_mat;
                response];

%% Plot IRF    
plot(response_mat(1,:));
    
plot(response_mat(2,:));
    
plot(response_mat(3,:));
    
plot(response_mat(4,:));
    
plot(response_mat(5,:));

plot(response_mat(6,:));

plot(response_mat(7,:));
  
  


   
