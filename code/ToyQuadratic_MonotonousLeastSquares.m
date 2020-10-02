clear all
close all
% We consider a synthetic dataset of $N=30$ 
% points from the graph of a quadratic function where the  values were 
% generated uniformly on $[-2,2]$. The corresponding $y$-coordinates of the
% graph were  perturbed by additive Gaussian noise. We impose a monotonically
% increasing shape constraint on the interval $[0,2]$, and study the effect
% of the level of the added noise $\xi$ on the desired increasing property of the
% estimate without (KRR) and with (SOC) monotonic shape constraint. 
% Here $\sigma=0.5$ and $\lambda_f=10^{-4}$, while $\xi$ varies in the interval $[0,4]$.

% Computation should take 10s.
%% Definition of the dataset and of the kernel parameters for \xi:=noiseLvl=1
tmin=-2; tmax=2;
load('rng_KRR-SOC.mat'); rng(scurr);
f=@(X) abs(X).^(2);
n=30; X=sort((tmax-tmin)*rand(n,1)+tmin); Yref=f(X); %(1.2-X/tmax).*
ntest=200; xgrid=linspace(-0.2+tmin,tmax+0.1,ntest); Ygrid=f(xgrid);
noiseLvl=[1]; Y=Yref+noiseLvl*randn(n,1);
hgauss = @(u,sig) exp(-u.^2/(2*sig^2));
dhgauss = @(u,sig) u/sig^2.*exp(-u.^2/(2*sig^2));
d2hgauss = @(u,sig) (u.^2-sig^2)/sig^4.*exp(-u.^2/(2*sig^2));
h=hgauss; dh=dhgauss; d2h=d2hgauss;

sigX=0.5; lbdaR=1E-3;
%KRR solution
GXnaive=h(repmat(X,1,n)-repmat(X',n,1),sigX);
GXTestnaive=h(repmat(xgrid',1,n)-repmat(X',ntest,1),sigX);
Anaive=(GXnaive+n*lbdaR*eye(n))\Y;
normBound=1.01*sqrt(Anaive'*GXnaive*Anaive);
%% Definition of the virtual points and solution of SOC through CVXGEN
listNgap=100;%[30:10:50,100];
for ngap=listNgap
xgDelta=linspace(0,10*sigX,1000); deltaMax=xgDelta(find(islocalmin(-dh(xgDelta,sigX)),1));
TminInt=0; TmaxInt=2;
Xgap=linspace(TminInt,TmaxInt,ngap);Xgap=Xgap(2:end-1);
delta=(Xgap(2)-Xgap(1))/2;
if Xgap(2)-Xgap(1)>deltaMax
    Xgap=TminInt:deltaMax:TmaxInt;
    delta=deltaMax;
end
Xgap=Xgap'; Xaug=[X;Xgap]; ngap=size(Xgap,1);
eta=sqrt(-2*d2h(0,sigX)+2*d2h(delta,sigX));
XX_sdist_mat=repmat(Xaug,1,n+ngap)-repmat(Xaug',n+ngap,1);
XXtest_sdist_mat=repmat(xgrid',1,n+ngap)-repmat(Xaug',ntest,1);

GX=h(XX_sdist_mat,sigX); DGX=-dh(XX_sdist_mat,sigX); D2GX=-d2h(XX_sdist_mat,sigX); 
GXtest=h(XXtest_sdist_mat,sigX); DGXtest=-dh(XXtest_sdist_mat,sigX); 
D2GXtest=-d2h(XXtest_sdist_mat,sigX); 

Gtot=[GX(1:n,1:n),-DGX(1:n,n+1:end);-DGX(n+1:end,1:n),D2GX(n+1:end,n+1:end)];
Gtotsym=[GX(1:n,1:n),DGX(1:n,n+1:end);-DGX(n+1:end,1:n),D2GX(n+1:end,n+1:end)];
tol=1E-8; Gcons=sqrtm(Gtotsym+tol*eye(n+ngap));
cvx_begin
    cvx_precision best
    variables A(n+ngap,1)
    minimize( norm(Y-Gtot(1:n,:)*A))
    subject to %
        norm(Gcons*A)<= normBound;
       eta*norm(Gcons*A) <= min(Gtot(n+1:end,:)*[-A(1:n);A(n+1:end)]);
cvx_end
AbestT=A; 
end
%%
figure; hold on;
ofsetInt=0.1; VofsetInt=-1.8; VofsetBar=0.1;
plot([TminInt TmaxInt],[VofsetInt VofsetInt],'g','LineWidth',2,'HandleVisibility','off')
plot([TminInt TminInt],[VofsetInt-VofsetBar VofsetInt+VofsetBar],'g','LineWidth',2,'HandleVisibility','off')
plot([TmaxInt TmaxInt],[VofsetInt-VofsetBar VofsetInt+VofsetBar],'g','LineWidth',2,'HandleVisibility','off')
plot(xgrid,Ygrid,'LineWidth',2)
plot(xgrid,GXTestnaive*Anaive,'k','LineWidth',2)
plot(xgrid,GXtest(:,1:n)*AbestT(1:n)-DGXtest(:,n+1:end)*AbestT(n+1:end),'LineWidth',2)%
plot(X,Y,'+b')

axis([-2.5 2.5 -2 5.5])
lgd=legend([{'noiseless curve'};{'KRR solution'};{'SOC solution'};{'noisy data'}]);
lgd.Location='southwest';
% fig = gcf; fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'KRR_SOC_plot_Gbars','-dpdf')
%% Computation of the proportion and amount of violation of the monotonicity requirement.
noiselvl_arr=[0:0.1:4]; nIter=1000; rng shuffle
ntest=200; Xtest=linspace(0,tmax,ntest); Ygrid=f(Xtest);

loss_table=zeros(length(noiselvl_arr),3,nIter);
for iIter=1:nIter
for iNoise=1:length(noiselvl_arr)
X=sort((tmax-tmin)*rand(n,1)+tmin); Yref=f(X);
Y=Yref+noiselvl_arr(iNoise)*randn(n,1);

GXnaive=h(repmat(X,1,n)-repmat(X',n,1),sigX);
Anaive=(GXnaive+n*lbdaR*eye(n))\Y;

XXtest_sdist_mat=repmat(Xtest',1,n)-repmat(X',ntest,1);
GXTestnaive=h(XXtest_sdist_mat,sigX);
Ytest=GXTestnaive*Anaive;
Xgrid_bool=(Xtest>=0); subYtest=Ytest(Xgrid_bool);
Yincr=cummax(subYtest);
errY=Yincr-subYtest;

loss_table(iNoise,1,iIter)=sum(errY(errY>0))/sum(Xgrid_bool);%L1err
loss_table(iNoise,2,iIter)=sum(errY>0)/sum(Xgrid_bool);%L0err_abs
loss_table(iNoise,3,iIter)=sum(diff(subYtest)<0)/sum(Xgrid_bool);%L0err
end
end

median_table=squeeze(nanmedian(loss_table,3))';
fp_table=squeeze(prctile(loss_table,25,3))';
sp_table=squeeze(prctile(loss_table,75,3))';

figure; hold on
xlabel({'Amplitude of Gaussian noise ($\xi$)'},'Interpreter','latex')
plot(noiselvl_arr,median_table(3,:),'k','LineWidth',2)  
plot(noiselvl_arr,fp_table(3,:),'--k','HandleVisibility','off')
plot(noiselvl_arr,sp_table(3,:),':k','HandleVisibility','off')
ylabel({'Proportion of violation'},'Interpreter','latex')
yyaxis right
plot(noiselvl_arr,median_table(1,:),'b','LineWidth',2)
plot(noiselvl_arr,fp_table(1,:),'--b','HandleVisibility','off')
plot(noiselvl_arr,sp_table(1,:),':b','HandleVisibility','off')
ylabel({'Amount of violation'},'Interpreter','latex')
lgd=legend({'Median proportion of violation',...
    'Median amount of violation'},'Interpreter','latex');
lgd.Location='northwest';
% fig = gcf; fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'KRR_break','-dpdf')