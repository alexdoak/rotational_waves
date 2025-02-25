function [Hi,yGP,dydt,jac,psiGP,dpsidt] = aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa) 
dte=t(end)-t(end-1); MN=M*N;

yGP = 2*y(Ey)-y(Ey-M) - dte^2*(dydss(Ey));
dydt(Ey) = (yGP-y(Ey-M))/(2*dte);
jac(Ey) = (dyds(Ey).^2+dydt(Ey).^2); jac=jac';
psiGP = 2*psi(Ey)-psi(Ey-M) - dte^2*(jac(Ey).*pa.vort_fun(psi(Ey),Q));
dpsidt(Ey) =  (psiGP-psi(Ey-M))/(2*dte);
Hi(MN+Ey) = .5*(dpsidt(Ey).^2)./(dyds(Ey).^(2) + dydt(Ey).^(2)) + g*y(Ey) - B;
