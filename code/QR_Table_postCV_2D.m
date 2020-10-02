clear all
close all
addpath(genpath('JQR_datasets'))
% To reproduce Table 1 for the 2D UCI datasets.

% In our \tb{second set of experiments} we considered the problem of JQR 
% where the conditional quantiles are encoded by the pinball loss and the
% shape requirement to fulfill is the non-crossing property. We compared 
% the efficiency of the proposed SOC approach with the PDCD technique 
% \cite{sangnier16joint} which minimizes the same loss \eqref{def_QR_loss} 
% but with a \emph{soft} non-crossing encouraging regularizer. We considered 
% $9$ UCI benchmarks for the JQR task. While the typical dimension in the 
% shape-constrained literature is $d=1$, our datasets were selected with 
% $d\in\{1,2,3\}$. Each dataset was split  into training $(70\%)$ and test 
% $(30\%)$ sets; the split and the experiment were repeated twenty times. 
% For each split, we optimized the hyperparameters 
% $(\sigma,\tilde{\lambda}_f,\tilde{\lambda}_{\bm{b}})$ of SOC, searching 
% over a grid to  minimize the pinball loss through a 5-fold cross validation 
% on the training set. Particularly, the kernel bandwith $\sigma$ was 
% searched over the square root of the deciles of the squared pairwise 
% distance between the points $\{\b x_n\}_{n\in[N]}$. The upper bound 
% $\tilde{\lambda}_f$ on $\sum_{q\in [Q]}\|f_q\|^2_k$ was scanned in the 
% log-scale interval $[-1,2]$. The upper bound $\tilde{\lambda}_{\bm{b}}$ 
% on $\|\bm{b}\|_2$ was kept fixed:  
% $\tilde{\lambda}_{\bm{b}} = 10 \max_{n\in [N]}|y_n|$. We then learned a 
% model on the whole training set and evaluated it on the test set. The 
% covering of $K=\times_{r \in [d]}\left[\min\{(\b x_n)_r\}_{n\in [N]},
% \max\{(\b x_n)_r\}_{n\in[N]}\right]$ was carried out with 
% $\|\cdot\|_2$-balls of radius $\delta$ chosen such that the number $M$ of 
% added points was less than $100$.

% Computation should take 30s with the given parameters for the dataset
% "CobarOre". As it depends on the total number of points, it varies a
% lot between datasets. The total gridsearch could last days.

NAMES_2D = ["CobarOre","topo","caution","snowgeese"];

GridsPrct=80;
GridboundN=1;
maxIter=1;
listFiles=1;
savebool=false;
% GridsPrct=20:10:80; %UNCOMMENT TO PERFORM THE WHOLE GRIDSEARCH
% GridboundN=[0.1,0.5,1,5,10,15,20,40,60,80]; %UNCOMMENT TO PERFORM THE WHOLE GRIDSEARCH
% maxIter=20; %UNCOMMENT TO REPEAT 20 TIMES THE TRAIN/TEST SPLIT 
% listFiles=1:length(NAMES_2D); %UNCOMMENT TO CONSIDER ALL THE 2D DATASETS
% savebool=true; %UNCOMMENT TO SAVE THE RESULTS

hgauss = @(u,sig) exp(-u.^2/(2*sig^2)); h=hgauss;
Q=(0.1:0.2:0.9)'; nq=length(Q);

for nIter=1:maxIter
    tic
for fileNb=listFiles%4:5 
plossCV_hyper=zeros(length(GridsPrct),length(GridboundN));
for sPrct=GridsPrct
for boundN=GridboundN
    load(strcat(NAMES_2D(fileNb),".mat")); Y=y;
    nX=size(X,1); %sPrct=GridsPrct(iGrid); boundN=GridboundN(iGrid);
    sigX=sqrt(prctile(pdist(unique(X),'squaredeuclidean'), sPrct));
    [Xnew,eta] = CompXgapEta2D(X,sigX,h); n=size(Xnew,1);
    PropTrain=0.7; randIdx=randperm(nX); 
    trainIdx=randIdx(1:floor(PropTrain*nX)); testIdx=randIdx(floor(PropTrain*nX)+1:end);
    Xtrain=X(trainIdx,:); ntrain=length(trainIdx);
    Xtest=X(testIdx,:); ntest=length(testIdx);

    tol=1E-4; boundB=10*max(abs(Y));

    GX=h(squareform(pdist(Xnew,'euclidean')),sigX); Gchol=chol(GX+tol*eye(n)); 
    CVX_ploss=@(Xtrain,Xtest)(CVX_pinball_loss_eval_wrapper(Xtrain,Xtest,X,Q,eta,Y,GX,Gchol,boundN,boundB));  
    vals_cv = crossval(CVX_ploss,Xtrain,'k',5);
    plossCV_hyper(find(GridsPrct==sPrct),find(GridboundN==boundN))=mean(vals_cv(:,1));
end
end
minploss=min(min(plossCV_hyper));
[row,col]=find(plossCV_hyper==minploss);
best_sPrct=GridsPrct(row); best_boundN=GridboundN(col);
sigX=sqrt(prctile(pdist(unique(X),'squaredeuclidean'), best_sPrct));
GX=h(squareform(pdist(Xnew,'euclidean')),sigX); Gchol=chol(GX+tol*eye(n));
test_ploss=(CVX_pinball_loss_eval_wrapper(Xtrain,Xtest,X,Q,eta,Y,GX,Gchol,best_boundN,boundB));

disp(strcat("Finished ",NAMES_2D(fileNb)))

if savebool
    saveFileName=['QR_Table_Results_2D.txt'];
    fileID = fopen(saveFileName,'a');
    fprintf(fileID,'%d \t %.2f \t %.2f \t %.2f \t %d \t %.2e \t %d \t %.2e \t %.2e \t %.2e\n',...
            fileNb,test_ploss, best_sPrct, sigX, best_boundN, boundB);
    fclose(fileID);
end
end
elapsedTime=toc
end