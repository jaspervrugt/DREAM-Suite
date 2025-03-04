% Going to zero
syms S A B Ks Ki I t
limit([S^2*(1-A)/(2*B) + A*(Ks-Ki).*(I-Ki*t) - t.*(A+B-1)*(Ks-Ki)^2 ] ./ ...
            ((Ks-Ki)*(A+B-A*B-1)),B,0,'right')
limit(S^2*(1-exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2))/(2*B),B,2,'right')
Answer: (Ki - Ks)*(I - Ki*t)
limit(exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2)*(Ks-Ki).*(I-Ki*t),B,0,'right')
Answer: -(Ki - Ks)*(I - Ki*t)
limit(t.*(exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2)+B-1)*(Ks-Ki)^2,B,0,'right')
Answer: 0
limit((Ks-Ki)*(exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2)+B-exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2)*B-1),B,0,'right')
Answer: 0
% We get 0/0 --> L'hospital's rule
% Do L'hospital's rule
NUM: S^2*(1-A)/(2*B) + A*(Ks-Ki).*(I-Ki*t) - t.*(A+B-1)*(Ks-Ki)^2
diff(S^2*(1-exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2))/(2*B) + exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2)*(Ks-Ki).*(I-Ki*t) - t.*(exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2)+B-1)*(Ks-Ki)^2,B)
ANS: %(S^2*(exp(-(2*B*(Ki - Ks)*(I - Ki*t))/S^2) - 1))/(2*B^2) + t*(Ki - Ks)^2*((2*exp(-(2*B*(Ki - Ks)*(I - Ki*t))/S^2)*(Ki - Ks)*(I - Ki*t))/S^2 - 1) + (2*exp(-(2*B*(Ki - Ks)*(I - Ki*t))/S^2)*(Ki - Ks)^2*(I - Ki*t)^2)/S^2 + (exp(-(2*B*(Ki - Ks)*(I - Ki*t))/S^2)*(Ki - Ks)*(I - Ki*t))/B
limit((S^2*(exp(-(2*B*(Ki - Ks)*(I - Ki*t))/S^2) - 1))/(2*B^2) + t*(Ki - Ks)^2*((2*exp(-(2*B*(Ki - Ks)*(I - Ki*t))/S^2)*(Ki - Ks)*(I - Ki*t))/S^2 - 1) + (2*exp(-(2*B*(Ki - Ks)*(I - Ki*t))/S^2)*(Ki - Ks)^2*(I - Ki*t)^2)/S^2 + (exp(-(2*B*(Ki - Ks)*(I - Ki*t))/S^2)*(Ki - Ks)*(I - Ki*t))/B,B,0,'right')
ANS: t*(Ki - Ks)^2*((2*(Ki - Ks)*(I - Ki*t))/S^2 - 1) + ((Ki - Ks)^2*(I - Ki*t)^2)/S^2
DEN: ((Ks-Ki)*(A+B-A*B-1))
diff((Ks-Ki)*(exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2)+B-exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2)*B-1),B)
ANS: % (Ki - Ks)*(exp(-(2*B*(Ki - Ks)*(I - Ki*t))/S^2) + (2*exp(-(2*B*(Ki - Ks)*(I - Ki*t))/S^2)*(Ki - Ks)*(I - Ki*t))/S^2 - (2*B*exp(-(2*B*(Ki - Ks)*(I - Ki*t))/S^2)*(Ki - Ks)*(I - Ki*t))/S^2 - 1)
limit((Ki - Ks)*(exp(-(2*B*(Ki - Ks)*(I - Ki*t))/S^2) + (2*exp(-(2*B*(Ki - Ks)*(I - Ki*t))/S^2)*(Ki - Ks)*(I - Ki*t))/S^2 - (2*B*exp(-(2*B*(Ki - Ks)*(I - Ki*t))/S^2)*(Ki - Ks)*(I - Ki*t))/S^2 - 1),B,0,'right')
ANS: (2*(Ki - Ks)^2*(I - Ki*t))/S^2
% ANSWER: lim_B->0 : t*(Ki - Ks)^2*((2*(Ki - Ks)*(I - Ki*t))/S^2 - 1) +
%                       ((Ki - Ks)^2*(I - Ki*t)^2)/S^2 / (2*(Ki - Ks)^2*(I - Ki*t))/S^2
% CORRECT - MATCHES JACOBIAN VALUES OF THIRD COLUMN IF BETA --> 0

% NOW DO LIMITS FOR BETA --> 2
syms S A B Ks Ki I t
limit([S^2*(1-A)/(2*B) + A*(Ks-Ki).*(I-Ki*t) - t.*(A+B-1)*(Ks-Ki)^2 ] ./ ...
            ((Ks-Ki)*(A+B-A*B-1)),B,2,'right')
limit(S^2*(1-exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2))/(2*B),B,2,'left')
Answer: -(S^2*(exp(-(4*(Ki - Ks)*(I - Ki*t))/S^2)/2 - 1/2))/2
limit(exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2)*(Ks-Ki).*(I-Ki*t),B,2,'left')
Answer: -exp(-(4*(Ki - Ks)*(I - Ki*t))/S^2)*(Ki - Ks)*(I - Ki*t)
limit(-t.*(exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2)+B-1)*(Ks-Ki)^2,B,2,'left')
Answer: -t*(Ki - Ks)^2*(exp(-(4*(Ki - Ks)*(I - Ki*t))/S^2) + 1)
limit((Ks-Ki)*(exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2)+B-exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2)*B-1),B,2,'left')
Answer: (Ki - Ks)*(exp(-(4*(Ki - Ks)*(I - Ki*t))/S^2) - 1)
% Thus, we get: 
LIMB2 = [ -(S^2*(exp(-(4*(Ki - Ks)*(I - Ki*t))/S^2)/2 - 1/2))/2 - exp(-(4*(Ki - Ks)*(I - Ki*t))/S^2)*(Ki - Ks)*(I - Ki*t) ...
    -t*(Ki - Ks)^2*(exp(-(4*(Ki - Ks)*(I - Ki*t))/S^2) + 1) ] ./ (Ki - Ks)*(exp(-(4*(Ki - Ks)*(I - Ki*t))/S^2) - 1)

limit([0 + Ks*I - 0 ] ./ ...
            ((Ks-0)*(1+B-1*B-1)),B,0,'right')
limit((1-exp(2*B*(Ks-Ki)*(I-Ki*t)/S^2))/(2*B),B,0)       
        
% pIpB = [ S^2*(1-A)/(2*B) + A*(Ks-Ki).*(I-Ki*t) - t.*(A+B-1)*(Ks-Ki)^2 ] ./ ...
%             ((Ks-Ki)*(A+B-A*B-1));