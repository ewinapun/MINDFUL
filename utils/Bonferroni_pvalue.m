function [final_pvalue, pvalues_matrix, F0_vec, F1_vec, betahat_2] = Bonferroni_pvalue(X1,Y1,X2,Y2)
% Bonferroni_pvalue uses Bonferroni method for the multiple hypothesis test
% with the null hypothesis of whether the coefficients of the linear
% regression model on two sets of data (X1,Y1) and (X2,Y2) are equal in
% each dimension
% 
% inputs:
% 
% Xi - 
%     (t*d) matrix of associated neural data where d is the neural data dimension.
% 
% Yi -
%     (t*k) matrix of movement intentions (represent the two segments needed to be compared).

% outputs:
% 
% final_pvalue -
%     p-value of the multiple hypothesis test.
% 
% pvalues_matrix -
%     (d by 2) matrix of p-values whose values in i'th row are
%     p-values of sigma-test (in first column) and beta-test (in second column)
%     done on (Y1,X1(:,i)) and (Y2,X2(:,i)).
% 
% F0_vec & F1_vec 
%     (1 by d) vector of test statistics of sigma (beta) test at each dimension.
%
%--------------------------------------------------------------------------
% History:
%   2022   Copyright Mona Khoshnevis, Tsam Kiu Pun, Brown University.
%--------------------------------------------------------------------------

% augment Y1 and Y2 to include a bias term
Y1 = horzcat(ones(size(Y1,1),1),Y1);
Y2 = horzcat(ones(size(Y2,1),1),Y2);

d = size(X1,2);
T_1 = size(Y1,1);
T_2 = size(Y2,1);
k  = size(Y1,2);

pvalues_matrix = zeros(d,2);

S_1 = Y1' * Y1;
S_2 = Y2' * Y2;

inv_S_1 = inv(S_1);
inv_S_2 = inv(S_2);

betahat_1 = inv_S_1 * Y1' * X1; 
betahat_2 = inv_S_2 * Y2' * X2;

e_1 = X1 - Y1 * betahat_1;
e_2 = X2 - Y2 * betahat_2;

s_1_sq = sum(e_1.^2,1);
s_2_sq = sum(e_2.^2,1);

% sigma-test
F0_vec = s_2_sq/(T_2 - k)./(s_1_sq/(T_1 - k));

a = fcdf(F0_vec,T_2 - k,T_1 - k);

pvalues_matrix(:,1) = 2*min(a, 1-a)';

% beta-test
F1_vec = ((T_1 + T_2-2*k)/k)*(diag((betahat_1 - betahat_2)'*(inv(inv_S_1 + inv_S_2))*(betahat_1-betahat_2))./(s_1_sq + s_2_sq)')';

pvalues_matrix(:,2) = 1-fcdf(F1_vec, k, T_1 + T_2 - 2*k);

final_pvalue = min([min(2*d*pvalues_matrix,[],'all') ,1]);

end
