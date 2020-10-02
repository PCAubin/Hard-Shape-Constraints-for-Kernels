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
% monotonically increasing \emph{and} concave. 

% We do a pre-processing of the dataset taking out the linear trend before
% solving the optimization problem (and we adapt the monotonicity
% constraint to take into account this pre-processing). We also take advantage
% that for univariate concave functions increasingness can be imposed by a
% single boundary constraint (rather than adding a large number of SOC
% constraints).

% Computation should take about 10 min for the whole dataset (220 points),
% and 30s for 50 points.

load(strcat("engel",".mat"))

[X,reInd]=sort(X);y=y(reInd); X(end-10:end)=[];y(end-10:end)=[];
f_0=[ones(size(X,1),1),X]\y;
linCoeff=f_0(2); y=y-linCoeff*X;
smIndx=unique([randperm(size(X,1),30)]); 
% smIndx=1:size(X,1); %UNCOMMENT TO CONSIDER THE WHOLE DATASET
Y=y(smIndx); nX=length(smIndx);
[X,reInd]=sort(X(smIndx,:));Y=Y(reInd);
%%
Prct=50; ntest=100; xgrid=linspace(-0.1+min(X),max(X)+0.1,ntest); 
sigX=sqrt(prctile(pdist(X,'squaredeuclidean'), Prct));
hgauss = @(u,sig) exp(-u.^2/(2*sig^2));dhgauss = @(u,sig) u/sig^2.*exp(-u.^2/(2*sig^2));
d2hgauss = @(u,sig) (u.^2-sig^2)/sig^4.*exp(-u.^2/(2*sig^2));
d3hgauss = @(u,sig) (u.^3-3*sig^2*u)/sig^6.*exp(-u.^2/(2*sig^2));
d4hgauss = @(u,sig) -(u.^4-6*sig^2*u.^2+3*sig^4)/sig^8.*exp(-u.^2/(2*sig^2));
h=hgauss; d2h=d2hgauss; d4h=d4hgauss;

[Xnew,eta,etaD,etaD2] = CompXgapEta1D_convexity(X,sigX,h,d2h,d4h); n=size(Xnew,1);

Q=(0.1:0.2:0.9)'; nq=length(Q); qmat=repmat(Q',nX,1);
XX_sdist_mat=repmat(Xnew,1,n)-repmat(Xnew',n,1);
XXtest_sdist_mat=repmat(xgrid',1,n)-repmat(Xnew',ntest,1);

GX=h(XX_sdist_mat,sigX); D2GX=-d2h(XX_sdist_mat,sigX); D4GX=-d4h(XX_sdist_mat,sigX); 
GXtest=h(XXtest_sdist_mat,sigX); D2GXtest=-d2h(XXtest_sdist_mat,sigX); 

dGX_tmax=dhgauss(repmat(max(Xnew),1,n)-repmat(Xnew',1,1),sigX);
d3GX_tmax=d3hgauss(repmat(max(Xnew),1,n)-repmat(Xnew',1,1),sigX);


tol=1E-4; 
Gchol=chol(GX+tol*eye(n)); diffmat=(eye(nq,nq-1)-[zeros(nq-1,1),eye(nq-1)]');
Gtot=[GX,-D2GX;-D2GX,D4GX];
tolD=1E-4; Gcons=sqrtm(Gtot+tolD*eye(2*n));

%%
boundNorm=10; boundB=2*max(abs(Y));
cvx_begin
    cvx_precision low
    variables A(2*n,nq) Apos(nX,nq) Aneg(nX,nq) b(nq)
    minimize(sum(sum(qmat.*Apos + (1-qmat).*Aneg))/nX)
    subject to
       Apos - Aneg == Gtot(1:nX,:)*A +repmat(b',nX,1) - repmat(Y,1,nq);
       Apos >= 0;
       Aneg >= 0;
       -dGX_tmax*A(1:n,:) - d3GX_tmax*A(n+1:end,:) >=-linCoeff;
       norms(Gcons*A,2,1)<= boundNorm;
       norm(b)<= boundB;
       norms(Gcons*A*diffmat,2,1)<= min(repmat(1./eta,1,nq-1).*...
           (Gtot(1:n,:)*A*diffmat+repmat(b'*diffmat,n,1)));
       norms(Gcons*A,2,1) <= min(-repmat(1./etaD2,1,nq).*(Gtot(n+1:end,:)*[A(1:n,:);A(n+1:end,:)]));% max(eta)*
cvx_end
AbestT=A; bBestT=b;

Yvals=linCoeff*xgrid'+repmat(b',ntest,1)+GXtest(:,1:n)*AbestT(1:n,:)-D2GXtest*AbestT(n+1:end,:);
%%
ybottom=-2; off_bottom=0.2; VofsetInt=ybottom+off_bottom;
% save('QR_engel_plot_concave_new.mat','X','Y','Yvals','xgrid','nX','n','Xnew','linCoeff')
% load('QR_engel_plot_concave_crop_new.mat')
figure
hold on
plot([min(X) max(X)],[VofsetInt VofsetInt],'g','LineWidth',2,'HandleVisibility','off')
plot([min(X) min(X)],[ybottom max(Y+linCoeff*X)],'g','LineWidth',2,'HandleVisibility','off')
plot([max(X) max(X)],[ybottom max(Y+linCoeff*X)],'g','LineWidth',2,'HandleVisibility','off')
plot(X,Y+linCoeff*X,'+b')
plot(xgrid,Yvals,'HandleVisibility','off','LineWidth',2)
xlabel({'Total income'},'Interpreter','latex')
ylabel({'Income spent on food'},'Interpreter','latex')
axis([min(X)-0.2 max(X)+0.2 ybottom max(Y+linCoeff*X)])
