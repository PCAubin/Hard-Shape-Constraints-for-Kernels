function result = CVX_pinball_loss_eval_wrapper(Xtrain,Xtest,X,Q,eta,Y,GX,Gchol,boundNorm,boundB)
%CVX_PINBALL_LOSS_EVAL_WRAPPER calls the proper CVX_PINBALL_LOSS_EVAL
%function depending on the dimension of the data (here d<=3)
trainIdx = ismember(X,Xtrain,'rows');
testIdx = ismember(X,Xtest,'rows');
nX=size(X,1);
if size(X,2)==1
    result=CVX_pinball_loss_eval(Q,eta,Y,GX,Gchol,boundNorm,boundB,trainIdx,testIdx);
elseif size(X,2)==2
    result=CVX_pinball_loss_eval_2D(Q,eta,Y,GX,Gchol,boundNorm,boundB,trainIdx,testIdx,nX);
elseif size(X,2)==3
result=CVX_pinball_loss_eval_2D(Q,eta,Y,GX,Gchol,boundNorm,boundB,trainIdx,testIdx,nX);
else
    disp("Dimension too large")
end
end