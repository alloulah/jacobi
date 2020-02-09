%
% Accelerated iterative Jacobi matrix inversion subroutine
%
% Reference:
%   Molisch, A.F.; Toeltsch, M.; Vermani, S.., "Iterative Methods for Cancellation
%   of Intercarrier Interference in OFDM Systems," IEEE Transactions on Vehicular 
%   Technology, vol.56, no.4, pp.2158,2167, July 2007
%
% Author: Mo Alloulah
% Date: 220514
%
function X_vect = jacobi( H_mtrx, Y_vect, N, err, isAccelConvrg, errNormType )
    % Solve: 
    %   Y_vect = H_mtrx * X_vect
    
    % The diagonal entries of H_mtrx and their inverses
    D_vect = diag( H_mtrx );
    
    if ~all( D_vect )
        error 'at least one diagonal entry is zero';
    end
    
    Dinv_vect = D_vect.^-1;
    % Matrix of off-diagonal entires of H
    D_mtrx = diag( D_vect );
    Hoff_mtrx = H_mtrx - D_mtrx;
    
    % Use
    %   D_vect.^-1 * Y_vect 
    % as the first approximation to X_vect
    Yinvd_vect = Dinv_vect .* Y_vect;
    X_vect = Yinvd_vect;
    
    %                   -1
    % Iterate X_vect = D_vect (Y_vect - H_mtrx * X_vect)
    %                                    off
    for k = 1:N
        Xprev_vect = X_vect;
        X_vect = Yinvd_vect - Dinv_vect .* (Hoff_mtrx * X_vect);
        
        if ( isAccelConvrg )
            % delayline
            switch ( mod(k, 4) )
                case 1
                    Xim3_vect = X_vect;
                case 2
                    Xim2_vect = X_vect;
                case 3
                    Xim1_vect = X_vect;
                case 0
                    Xi_vect = X_vect;
                    % Accelerate convergence
                    fprintf('Acceleration %d!\n', floor(k/4));
                    Xprev_vect = Xi_vect;
                    [X_vect, alpha_vect] = ...
                        accelConvrg_2nd( Xi_vect, Xim1_vect, Xim2_vect, Xim3_vect );     
            end   
        end
        
        if ( norm( X_vect - Xprev_vect, errNormType ) < err )
            fprintf('Finished in %d interations.\n', k);
            return;
        end
    end
    
    fprintf('The method did not converge.\n');
    X_vect = -1;
end