%
% Accelerated iterative Jacobi matrix inversion for inter-carrier
% interference (ICI)
%
% Reference:
%   Molisch, A.F.; Toeltsch, M.; Vermani, S.., "Iterative Methods for Cancellation
%   of Intercarrier Interference in OFDM Systems," IEEE Transactions on Vehicular 
%   Technology, vol.56, no.4, pp.2158,2167, July 2007
%
% Author: Mo Alloulah
% Date: 220514
%
%%

clear all
clc
close all

%% Create diagonal matrix

M = 64;

diagBoostFactor = M/2;

ITER = 150;
epsilon = 1e-5;
normType = Inf;

% Create a random M-by-M matrix.
zMag = rand();
zPhase = 2*pi*rand();
H_mtrx = (rand(M, M) + 1i*rand(M, M))*zMag/diagBoostFactor + eye(M)*zMag*exp(1i*zPhase);
Y_vect = rand(M, 1) + 1i*rand(M, 1);

% Matlab native reference
X_ref = H_mtrx \ Y_vect;

%% Visualise & Inspect

Z = 20*log10( abs(H_mtrx) + eps );

surfView = [90 90];

figure();
surf( Z )
axis([1 M 1 M]);
view(surfView);
grid();

%% Without accelation

X = jacobi( H_mtrx, Y_vect, ITER, epsilon, 0, normType );

%% With acceleration

X_acc = jacobi( H_mtrx, Y_vect, ITER, epsilon, 1, normType );

%% 

err = norm(X_ref - X_acc, Inf);
fprintf('Jacobi method error = %e\n', err);
