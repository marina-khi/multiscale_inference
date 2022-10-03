function out = fn_MG(v_psi)

N = length(v_psi);
psiMG = mean(v_psi);
v_psi_exc = v_psi - psiMG; % N x 1
v_psi_excsq = v_psi_exc .* v_psi_exc;
somma = sum(v_psi_excsq);
varianza = somma / N / (N - 1);
psiMG_se = realsqrt(varianza);
psiMG_atval = abs(psiMG / psiMG_se);
psiMG_stars = (psiMG_atval > 1.645) + (psiMG_atval > 1.96) + (psiMG_atval > 2.33);

% save these three scalars into a struc
out.MG = psiMG;
out.MG_se = psiMG_se;
out.MG_stars = psiMG_stars;

end
