clear all
close all
addpath(genpath('JQR_datasets'))

% In our second example with multiple shape constraints, we focused on the 
% analysis of more than $300$ \tb{aircraft trajectories} which describe the 
% radar-measured altitude ($y$) of aircrafts flying between two cities 
% (Paris and Toulouse) as a function of time ($x$). These trajectories were 
% restricted to their takeoff phase (where the monotone increasing property 
% should hold), giving rise to a total number of samples $N=15657$. We 
% imposed non-crossing and monotonicity property. 

% Computation should take about 10 min for the whole dataset (300 trajectories),
% and 20s for 10 trajectories.

load('ENAC_data_small_crop.mat')
[X,reInd]=sort(X);y=y(reInd);ID=ID(reInd);

minID=1000; maxID=1010;
minID=1000; maxID=2000; %UNCOMMENT TO CONSIDER THE WHOLE DATASET

reInd=(ID<=maxID)&(ID>=minID); X=X(reInd); y=y(reInd); ID=ID(reInd);

[uX,idxX_U,idxU_X]=unique(X);
sigX=3*4; ntest=100; xgrid=linspace(-0.1+min(uX),max(uX)+0.1,ntest); 

nX=length(X);

hgauss = @(u,sig) exp(-u.^2/(2*sig^2)); dhgauss = @(u,sig) u/sig^2.*exp(-u.^2/(2*sig^2));
d2hgauss = @(u,sig) (u.^2-sig^2)/sig^4.*exp(-u.^2/(2*sig^2));
h=hgauss; dh=dhgauss; d2h=d2hgauss;

Xnew=uX;
eta=sqrt(2*h(0,sigX)-2*h((uX(2)-uX(1))/2,sigX));
etaD=sqrt(-2*d2h(0,sigX)+2*d2h((uX(2)-uX(1))/2,sigX));
n=size(Xnew,1);

Q=(0.1:0.2:0.9)'; nq=length(Q); qmat=repmat(Q',nX,1);
%%
Y=y/max(y);
XX_sdist_mat=repmat(Xnew,1,n)-repmat(Xnew',n,1);
XXtest_sdist_mat=repmat(xgrid',1,n)-repmat(Xnew',ntest,1);

GX=h(XX_sdist_mat,sigX); DGX=-dh(XX_sdist_mat,sigX); D2GX=-d2h(XX_sdist_mat,sigX); 
GXtest=h(XXtest_sdist_mat,sigX); DGXtest=-dh(XXtest_sdist_mat,sigX); 
D2GXtest=-d2h(XXtest_sdist_mat,sigX);

tol=1E-4; 
Gtot=[GX,-DGX;-DGX,D2GX]; Gtotsym=[GX,DGX;-DGX,D2GX];
diffmat=(eye(nq,nq-1)-[zeros(nq-1,1),eye(nq-1)]');
tolD=1E-4; Gcons=sqrtm(Gtotsym+tolD*eye(2*n));
matU_X=zeros(nX,n); matU_X((1:nX)'+(idxU_X-1)*nX)=1;
%%
boundNorm=100; boundB=(1E1)*max(abs(Y));
cvx_begin
    variables A(2*n,nq) errorY(nX,nq) b(nq)
    minimize(sum(sum(max(qmat.*errorY,-(1-qmat).*errorY)))/nX)
    subject to
       errorY == matU_X*Gtot(1:n,:)*A +repmat(b',nX,1) - repmat(Y,1,nq);

       norms(Gcons*A,2,1)<= boundNorm;
       norm(b)<= boundB;
       eta*norms(Gcons*A*diffmat,2,1)<= min...
           (Gtot(1:n,:)*A*diffmat+repmat(b'*diffmat,n,1));
       etaD*norms(Gcons*A,2,1)/5 <= min(Gtot(n+1:end,:)*[-A(1:n,:);A(n+1:end,:)]);
cvx_end
AbestT=A; bBestT=b;

Yvals=repmat(bBestT',ntest,1)+GXtest(:,1:n)*AbestT(1:n,:)-DGXtest*AbestT(n+1:end,:);
%%

% ymax=max(y);
% save('QR_enac_data_monotone_new.mat','X','Y','ID','Yvals','xgrid','nX','n','Xnew','ymax')
load('QR_enac_data_monotone_new.mat')
y=Y;ymax=max(y);
figure
hold on;
plot(xgrid,Yvals*ymax,'LineWidth',2,'HandleVisibility','off')%
for i=min(ID):max(ID)
idx_temp=ID==i;
plotBlob=plot(X(idx_temp),y(idx_temp),'--');
plotBlob.Color(4) = 0.3;
end

axis tight
xlabel('Time (s)')
ylabel('Altitude (*100 ft)')
fig = gcf; fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'QR_enac_data_monotone','-dpdf')