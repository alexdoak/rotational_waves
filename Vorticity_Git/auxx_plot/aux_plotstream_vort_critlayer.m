function[] = aux_plotstream_vort_critlayer(X0_list,Y0_list,R,pa,un,de)


if numel(X0_list)>0
X=un.X; Y=un.Y; Psi=un.Psi; psi=un.psi; 
Dpsids = reshape(de.dpsids,pa.M,pa.N);
Dpsidt = reshape(de.dpsidt,pa.M,pa.N);

for p = 1:numel(X0_list)
    X0 = X0_list(p);
    Y0 = Y0_list(p);
   
    I = find((X-X0).^2 + (Y-Y0).^2 < R^2);
    [~,I2] = min(Dpsids(I).^2 + Dpsidt(I).^2);
    I=I(I2);

    psival = un.psi(I);
    contour(X,Y,Psi,[psival-10^(-11),psival+10^(-11)],'b','Linewidth',1.3); 
    hold on;
    contour(-X,Y,Psi,[psival-10^(-11),psival+10^(-11)],'b','Linewidth',1.3); 



  %  contour(pa.wavelength+X,Y,Psi,[psival-10^(-11),psival+10^(-11)],'b','Linewidth',1.3); 
  %  hold on;
  %  contour(pa.wavelength+-X,Y,Psi,[psival-10^(-11),psival+10^(-11)],'b','Linewidth',1.3); 
end
end
