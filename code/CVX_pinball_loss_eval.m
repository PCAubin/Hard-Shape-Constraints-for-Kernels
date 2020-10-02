function [result, A, b]  =...
    CVX_pinball_loss_eval(Q,eta,Y,GX,Gchol,boundNorm,boundB,trainIdx,testIdx)
%CVX_PINBALL_LOSS_EVAL_1D solves the joint quantile regression problem with
%uniform covering (eta is a constant) for 1D data, with non-crossing SOC 
%constraint and constant biases.
%Computation is done for the training dataset while evaluation of 100* pinball 
%loss is done for the test dataset. A few relevant outputs of cvx are also
%stored.
nq=length(Q); ntrain=sum(trainIdx); ntest=sum(testIdx); n=size(GX,1);
diffmat=(eye(nq,nq-1)-[zeros(nq-1,1),eye(nq-1)]'); qmat=repmat(Q',ntrain,1);
cvx_begin quiet
%     cvx_precision best
    variables A(n,nq) Apos(ntrain,nq) Aneg(ntrain,nq) b(nq)
    minimize(100*sum(sum(qmat.*Apos + (1-qmat).*Aneg))/ntrain)
    subject to
       Apos - Aneg == GX(trainIdx,:)*A +repmat(b',ntrain,1) - repmat(Y(trainIdx),1,nq);
       Apos >= 0;
       Aneg >= 0;
       norms(Gchol*A,2,1)<= boundNorm;
%        sum(sum_square(Gchol*A,2))<= boundNorm;
%        norms(reshape(Gchol*A,[],1))<= boundNorm; 
       norm(b)<= boundB;
       norms(Gchol*A*diffmat,2,1)<= min(repmat(1./eta,1,nq-1).*(GX*A*diffmat+repmat(b'*diffmat,n,1)));% max(eta)*
cvx_end

errorYTest=GX(testIdx,:)*A +repmat(b',ntest,1) - repmat(Y(testIdx),1,nq);
testval=100*sum(sum(max(repmat(Q',ntest,1).*errorYTest,-(1-repmat(Q',ntest,1)).*errorYTest)))/ntest;
result=[testval, cvx_optval, cvx_cputime, cvx_slvitr, cvx_slvtol];
end

