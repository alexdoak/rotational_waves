function [] = aux_plotstream_vort_fig3(un,pa,de,fe,Nstream,X0_list,Y0_list,Yrange_list,R)


%k=me.k; psi=un.psi; M=pa.M;
%dpsidst = (psi(k+1+M) + psi(k-1-M) - psi(k-1+M) - psi(k+1-M))/(4*me.ds*me.dt);
%gradpsi = [de.dpsidt,de.dpsids];



M=pa.M;
X=un.X;Y=un.Y;Psi=un.Psi;vort_fun=pa.vort_fun; Q=un.Q; 
wl=pa.wavelength_true;

%set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
%set(0,'defaulttextInterpreter','latex')
% Create a figure
%fig = figure();
%set(fig,'Units','centimeters','Position',[2 2 22 7]);





contour(X,Y,Psi, Nstream,'k'); hold on;
%contour(X,Y,Psi,[-10^(-10),10^(-10)],'k','LineWidth',2);
%contour(X,Y,Psi,[Q-10^(-10),Q+10^(-10)],'k');
contour(-X,Y,Psi,Nstream,'k');
%contour(-X,Y,Psi,[-10^(-10),10^(-10)],'k','LineWidth',2);
%contour(-X,Y,Psi,[Q-10^(-10),Q+10^(-10)],'k');
plot(un.X(:,end),un.Y(:,end),'k','LineWidth',1);
plot(-un.X(:,end),un.Y(:,end),'k','LineWidth',1);
plot(un.X(:,1),un.Y(:,1),'k','LineWidth',2);
plot(-un.X(:,1),un.Y(:,1),'k','LineWidth',2);

%%% Stagnation streamlines
aux_plotstream_vort_critlayer_fig3(X0_list,Y0_list,Yrange_list,R,pa,un,de)

% X=un.X; Y=un.Y; U=un.U; V=un.V; q=un.U.^2+un.V.^2;
% minval = min(q,[],'all')
% [row, col] = find(q == minval);
% row=row(1); col=col(1);
% psival = un.Psi(row,col);
% contour(X,Y,Psi,[psival-10^(-11),psival+10^(-11)],'b','Linewidth',2);  %linspace(psival-.5*10^(-3),psival+.5*10^(-3),17),'b');
% contour(-X,Y,Psi,[psival-10^(-11),psival+10^(-11)],'b','Linewidth',2);  %linspace(psival-.5*10^(-3),psival+.5*10^(-3),17),'b');



%{
I = find(de.dpsids.^2 + de.dpsidt.^2 < 10^(-4));
psival = un.psi(I);
psival2=[psival(1)];
for p=1:numel(psival)-1
    if abs(psival(p+1)-psival(p))<10^(-4)
    else
        psival2=[psival2,psival(p+1)];
    end
end

contour(X,Y,Psi,[psival2],'b','Linewidth',2); 
contour(-X,Y,Psi,[psival2],'b','Linewidth',2); 
%}



%{
contour(wl-X,Y,Psi,Nstream,'k'); hold on;
contour(wl-X,Y,Psi,[-10^(-10),10^(-10)],'k','LineWidth',2);
contour(wl-X,Y,Psi,[Q-10^(-10),Q+10^(-10)],'k');

contour(wl+X,Y,Psi,Nstream,'k'); hold on;
contour(wl+X,Y,Psi,[-10^(-10),10^(-10)],'k','LineWidth',2);
contour(wl+X,Y,Psi,[Q-10^(-10),Q+10^(-10)],'k');
%}

if fe.Stratified==1
    subplot(1,3,3);
    plot(pa.strat(un.psi(M:M:end),un.c),un.y(M:M:end));
    xlabel('rho')
    ylabel('y')
end

%{
contour(wl-X,Y,Psi,Nstream,'k');
%contour(wl-X,Y,Psi,[-10^(-10),10^(-10)],'k','LineWidth',2);
contour(wl-X,Y,Psi,[Q-10^(-10),Q+10^(-10)],'k');

contour(wl+X,Y,Psi,Nstream,'k');
%contour(wl+X,Y,Psi,[-10^(-10),10^(-10)],'k','LineWidth',2);
contour(wl+X,Y,Psi,[Q-10^(-10),Q+10^(-10)],'k');


%contour(2*wl-X,Y,Psi,Nstream,'k');
%contour(2*wl-X,Y,Psi,[-10^(-10),10^(-10)],'k','LineWidth',2);
%contour(2*wl-X,Y,Psi,[Q-10^(-10),Q+10^(-10)],'k');
%}
%{
xlabel('x')
ylabel('y','Rotation',0)
ylim([min(min(Y))-0.04,max(max(Y))+0.04]);
xlim([-max(max(X)),max(max(X))]);
%}

%{
if fe.Stratified==1
rho=un.rho2;
pcolor(X,Y,rho); hold on;
pcolor(-X,Y,rho);
%pcolor(wl-X,Y,rho);
%pcolor(wl+X,Y,rho);
%pcolor(2*wl-X,Y,rho); 
clim([min(min(rho)),max(max(rho))]); alpha .8; shading interp; colorbar;
title('Density and streamlines')
ylim([min(min(Y))-0.04,max(max(Y))+0.04]);
xlim([-max(max(X)),max(max(X))]);

subplot(1,3,2)
pcolor(X,Y,un.vort); hold on;
pcolor(-X,Y,un.vort);
shading interp; colorbar; title('Vorticity');
ylim([min(min(Y))-0.04,max(max(Y))+0.04]);
xlim([-max(max(X)),max(max(X))]);

subplot(1,3,3)
pcolor(un.X,un.Y,(un.U.^2 + un.V.^2).^.5/un.c); hold on;
pcolor(-un.X,un.Y,(un.U.^2 + un.V.^2).^.5/un.c); 
ylim([min(min(Y))-0.04,max(max(Y))+0.04]);
xlim([-max(max(X)),max(max(X))]);
shading interp; colorbar; title('$|U|/c$');

else
%}



%{
subplot(1,2,2)
plot(un.U(end,:),un.Y(end,:),'r','Linewidth',2);
hold on;
plot(un.vort(end,:),un.Y(end,:),'b','Linewidth',2);
%}



%{
subplot(1,2,2)
plot(un.Y(1,:),-un.U(1,:),'b','Linewidth',1);
xlabel('$y$')
ylabel('$-u$','Rotation',0)
xlim([0,un.Y(1,end)])
ylim([0,4.5])
%}


end


