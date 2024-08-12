function Jdiv = CalcJKL(pd1, pd2)
% calculate KL and symmetric KL for multivariate gaussian
% Copyright Tsam Kiu Pun, 2024. Brown University
% tsam_kiu_pun@brown.edu
% -------------------------------------------------------------------------
KLdiv1 = CalcKL(pd1, pd2);
KLdiv2 = CalcKL(pd2, pd1);
% Jeffrey's (symmetric)
Jdiv = (KLdiv1 + KLdiv2)/2;
end