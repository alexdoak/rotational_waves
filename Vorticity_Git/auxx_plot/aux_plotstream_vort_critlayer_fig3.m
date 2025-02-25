function[] = aux_plotstream_vort_critlayer_fig3(X0_list,Y0_list,Yrange_list,R,pa,un,de)


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

    contour(X(:,Yrange_list(p,1):Yrange_list(p,2)),Y(:,Yrange_list(p,1):Yrange_list(p,2)),Psi(:,Yrange_list(p,1):Yrange_list(p,2)),[psival-10^(-11),psival+10^(-11)],'b','Linewidth',1.3); 
    hold on;
    contour(-X(:,Yrange_list(p,1):Yrange_list(p,2)),Y(:,Yrange_list(p,1):Yrange_list(p,2)),Psi(:,Yrange_list(p,1):Yrange_list(p,2)),[psival-10^(-11),psival+10^(-11)],'b','Linewidth',1.3); 
end
end
