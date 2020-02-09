%
% Acceleration of convergence for iterative Jacobi matrix inversion
%
% Reference:
%   Molisch, A.F.; Toeltsch, M.; Vermani, S.., "Iterative Methods for Cancellation
%   of Intercarrier Interference in OFDM Systems," IEEE Transactions on Vehicular 
%   Technology, vol.56, no.4, pp.2158,2167, July 2007
%
% Author: Mo Alloulah
% Date: 220514
%
function [ X_vect, alpha_vect ] = accelConvrg_2nd( Xi_vect, Xim1_vect, Xim2_vect, Xim3_vect )
    
    % compute deltas
    deltaXi_vect = Xi_vect - Xim1_vect;
    deltaXim1_vect = Xim1_vect - Xim2_vect;
    deltaXim2_vect = Xim2_vect - Xim3_vect;
    
    % compute Mij
    M11_acc = sum( (1./deltaXim2_vect) .* ...
        ( (deltaXi_vect - deltaXim1_vect) .* (deltaXi_vect - deltaXim1_vect) ) );
    M12_acc = sum( (1./deltaXim2_vect) .* ...
        ( (deltaXi_vect - deltaXim1_vect) .* (deltaXi_vect - deltaXim2_vect) ) );
    M21_acc = sum( (1./deltaXim2_vect) .* ...
        ( (deltaXi_vect - deltaXim2_vect) .* (deltaXi_vect - deltaXim1_vect) ) );
    M22_acc = sum( (1./deltaXim2_vect) .* ...
        ( (deltaXi_vect - deltaXim2_vect) .* (deltaXi_vect - deltaXim2_vect) ) );
    % pack in matrix
    Macc_mtrx = [
        M11_acc M12_acc;
        M21_acc M22_acc;
        ];
    % invert
    Macc_det = M11_acc*M22_acc - M21_acc*M12_acc;
    MaccInv_mtrx = (1/Macc_det) * [  M22_acc  -M12_acc;
                                    -M21_acc   M11_acc; ];
    
    % compute bi's
    b1_acc = sum( (1./deltaXim2_vect) .* (deltaXi_vect .* (deltaXi_vect - deltaXim1_vect)) );
    b2_acc = sum( (1./deltaXim2_vect) .* (deltaXi_vect .* (deltaXi_vect - deltaXim2_vect)) );
    % pack in vector
    bacc_vect = [b1_acc; b2_acc];
    
    % compute alpha
    alphaPrime_vect = MaccInv_mtrx * bacc_vect;
    alpha0 = 1 - sum(alphaPrime_vect);
    alpha_vect = [alpha0; alphaPrime_vect];
    
    % linearly combine to accelerate convergence
    X_vect = ...
        alpha_vect(1) * Xi_vect + ...
        alpha_vect(2) * Xim1_vect + ...
        alpha_vect(3) * Xim2_vect;
    
end

