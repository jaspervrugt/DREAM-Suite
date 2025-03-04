function [R,E] = gchol(A)
% Computes the Gill-Murray generalized Cholesky decomposition of 
% M = A + E = dot(L,L') where
% M is safely symmetric positive definite (SPD) and well conditioned
% E is small (zero if A is already SPD and not much larger than the most
% negative eigenvalues of A)
% In some implementations, e is the diagonal of the correction matrix,
% thus, E = diag(e); 

% Translated from the following Gauss code of J. Gill::</span>
%   ****************************************************************
%   * This is the Gill/Murray cholesky routine.  Reference:        *
%   *							                                   *
%   * Gill, Jeff and Gary King. ``What to do When Your Hessian     *
%   * is Not Invertible: Alternatives to Model Respecification     *
%   * in Nonlinear Estimation,&#39;&#39; Sociological Methods and  *
%   * Research, Vol. 32, No. 1 (2004): Pp. 54--87.                 *
%   *							                                   *
%   * Abstract						                               *
%   * 							                                   *
%   * What should a researcher do when statistical analysis        *
%   * software terminates before completion with a message that    *
%   * the Hessian is not invertable? The standard textbook advice  *
%   * is to respecify the model, but this is another way of saying *
%   * that the researcher should change the question being asked.  *
%   * Obviously, however, computer programs should not be in the   *
%   * business of deciding what questions are worthy of            *
%   * study. Although noninvertable Hessians are sometimes signals *
%   * of poorly posed questions, nonsensical models, or            *
%   * inappropriate estimators, they also frequently occur when    *
%   * information about the quantities of interest exists in the   *
%   * data, through the likelihood function. We explain the        *
%   * problem in some detail and lay out two preliminary proposals *
%   * for ways of dealing with noninvertable Hessians without      *
%   * changing the question asked.                                 *
%   ****************************************************************

% Written by Jasper A. Vrugt, May 6, 2023

n = size(A,1);
R = eye(n);
E = zeros(n,n);
%norm_A = max(sum(abs(A)));
gamm = max(abs(diag(A)));
%delta = eps*max(1,norm_A); % R-code? %max(max(eps*norm_A,eps));
delta = sqrt(eps)*norm(A,'fro');
for j = 1:n
    theta_j = 0;
    for i = 1:n
        sumR = 0;
        for k = 1:i-1 %(i-1) in Python but not R code!!
            sumR = sumR + R(i,k)*R(j,k);
        end
        if i <= j
            R(j,i) = (A(i,j) - sumR)/R(i,i);
        end
        theta_j = max(theta_j,A(i,j) - sumR);
% %         %if (A(i,j)- sumR) > theta_j
% %         %    theta_j = A(i,j) - sumR;
% %         %end
%         if i > j
%             R(i,j) = 0;
%         end
    end
    theta_j2 = theta_j^2;
    sumR = 0;
    for k = 1:j-1 % j-1
        sumR = sumR + R(j,k)^2;
    end
    phi_j = A(j,j) - sumR;
%     if j + 1 < n
%         xi_j = 0;
%         for k = j+1:n-1
%             xi_j = 
    if j+1 <= n % <= n?
        xi_j = max(abs(A(j+1:n,j)));
    else
        xi_j = abs(A(n,j)); % -1 means last element in sequence
    end
%    beta_j = sqrt(max(max(gamm*(xi_j/n)*eps)));
    beta_j2 = max(max(gamm,xi_j/n),eps);
    if delta >= max(abs(phi_j),theta_j2/beta_j2)
        E(j,j) = delta - phi_j;
    elseif abs(phi_j) >= max( delta^2/beta_j2 , delta)
        E(j,j) = abs(phi_j) - phi_j;
    elseif max(delta,abs(phi_j)) < (theta_j2/beta_j2) 
        E(j,j) = (theta_j2/beta_j2) - phi_j;
    end
    R(j,j) = sqrt(A(j,j) - sumR + E(j,j));
end
R = R';
% R is the Cholesky decomposition