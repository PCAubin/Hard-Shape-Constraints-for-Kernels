clear all
close all
addpath(genpath('JQR_datasets'))

% In our \tb{third set of experiments}, we demonstrate the efficiency of our
% proposed SOC estimator on tasks with multiple ($I>1$) hard shape constraints.
% Our first example is drawn from \tb{economics}; we focused on JQR for the
% engel dataset, applying the same parameter optimization as in the second
%experiment. In this benchmark, the $\{(x_n,y_n)\}_{n\in[N]} \subset \R^2$
% pairs correspond to annual household income ($x_n$) and food expenditure
% ($y_n$), preprocessed to have zero mean and unit variance. Engel's law
% postulates a monotone increasing property of $y$ w.r.t.\ $x$, as well as
% concavity. We therefore constrained the quantile functions to be non-crossing,
% and monotonically increasing . 

% Computation should take about 10 min for the whole dataset (220 points),
% and 30s for 50 points.
load("engel.mat")
[X,reInd]=sort(X);y=y(reInd); X(end-10:end)=[];y(end-10:end)=[];
smIndx=unique([randperm(size(X,1),50)]); 
% smIndx=1:size(X,1); %UNCOMMENT TO CONSIDER THE WHOLE DATASET
Y=y(smIndx); nX=length(smIndx);
[X,reInd]=sort(X(smIndx,:));Y=Y(reInd);

Q=(0.1:0.2:0.9)'; nq=length(Q); qmat=repmat(Q',nX,1);
%%
Prct=50; ntest=100; xgrid=linspace(-0.1+min(X),max(X)+0.1,ntest); 
sigX=sqrt(prctile(pdist(X,'squaredeuclidean'), Prct));
hgauss = @(u,sig) exp(-u.^2/(2*sig^2));dhgauss = @(u,sig) u/sig^2.*exp(-u.^2/(2*sig^2));
d2hgauss = @(u,sig) (u.^2-sig^2)/sig^4.*exp(-u.^2/(2*sig^2));
d3hgauss = @(u,sig) (u.^3-3*sig^2*u)/sig^6.*exp(-u.^2/(2*sig^2));
d4hgauss = @(u,sig) (u.^4-6*u.^2+3*sig^2)/sig^8.*exp(-u.^2/(2*sig^2));
h=hgauss; dh=dhgauss; d2h=d2hgauss; d3h=d3hgauss; d4h=d4hgauss;

[Xnew,eta,etaD,etaD2] = CompXgapEta1D_convexity(X,sigX,h,d2h,d4h); n=size(Xnew,1);

%%
XX_sdist_mat=repmat(Xnew,1,n)-repmat(Xnew',n,1);
XXtest_sdist_mat=repmat(xgrid',1,n)-repmat(Xnew',ntest,1);

GX=h(XX_sdist_mat,sigX); DGX=-dh(XX_sdist_mat,sigX); D2GX=-d2h(XX_sdist_mat,sigX); 
GXtest=h(XXtest_sdist_mat,sigX); DGXtest=-dh(XXtest_sdist_mat,sigX); 
D2GXtest=-d2h(XXtest_sdist_mat,sigX);

tol=1E-4; 
Gchol=chol(GX+tol*eye(n)); diffmat=(eye(nq,nq-1)-[zeros(nq-1,1),eye(nq-1)]');
Gtot=[GX,-DGX;-DGX,D2GX]; Gtotsym=[GX,DGX;-DGX,D2GX];
tolD=1E-4; Gcons=sqrtm(Gtotsym+tolD*eye(2*n));
%%
boundNorm=10; boundB=10*max(abs(Y));
cvx_begin
%     cvx_precision low
    variables A(2*n,nq) Apos(nX,nq) Aneg(nX,nq) b(nq)
    minimize(sum(sum(qmat.*Apos + (1-qmat).*Aneg))/nX)
    subject to
       Apos - Aneg == Gtot(1:nX,:)*A +repmat(b',nX,1) - repmat(Y,1,nq);
       Apos >= 0;
       Aneg >= 0;
       norms(Gcons*A,2,1)<= boundNorm;
       norm(b)<= boundB;
       norms(Gcons*A*diffmat,2,1)<= min(repmat(1./eta,1,nq-1).*...
           (Gtot(1:n,:)*A*diffmat+repmat(b'*diffmat,n,1)));
       norms(Gcons*A,2,1) <= min(repmat(1./etaD,1,nq).*(Gtot(n+1:end,:)*[-A(1:n,:);A(n+1:end,:)]));% max(eta)*
cvx_end
AbestT=A; bBestT=b;
Yvals=repmat(b',ntest,1)+GXtest(:,1:n)*AbestT(1:n,:)-DGXtest*AbestT(n+1:end,:);
%%
% save('QR_engel_plot_monotone_crop_new.mat','X','Y','Yvals','xgrid','nX','n','Xnew')
load('QR_engel_plot_monotone_crop_new.mat')
ybottom=-2; off_bottom=0.2;
figure
hold on
ofsetInt=0.1; VofsetInt=ybottom+off_bottom;
plot([min(X) max(X)],[VofsetInt VofsetInt],'g','LineWidth',2,'HandleVisibility','off')
plot([min(X) min(X)],[ybottom max(Y)],'g','LineWidth',2,'HandleVisibility','off')
plot([max(X) max(X)],[ybottom max(Y)],'g','LineWidth',2,'HandleVisibility','off')
plot(X,Y,'+b')
plot(xgrid,Yvals,'HandleVisibility','off','LineWidth',2)
plot(X,repmat(ybottom+off_bottom,nX,1),'k+','HandleVisibility','off')
plot(Xnew(nX+1:end),repmat(ybottom+off_bottom,n-nX,1),'k+')
axis([min(X)-0.2 max(X)+0.2 ybottom max(Y)])

lgd=legend([{'data points $ \{(x_n, y_n)\}_{n\in [N]}$'};...%{'increasing JQR-SOC quantile functions $\{f_j\}_j$'};...
    {'virtual points $\{\tilde{x}_m\}_{m\in [M]}$'}],'Interpreter','latex');
lgd.Location='northwest';
xlabel({'Total income'},'Interpreter','latex')
ylabel({'Income spent on food'},'Interpreter','latex')
fig = gcf; fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'QR_engel_plot_concave_Gbars','-dpdf')