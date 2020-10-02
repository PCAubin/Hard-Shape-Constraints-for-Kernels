function [Xnew,eta] = CompXgapEta2D(Xorig,sigX,h,nbPts_per_dim)
%This function computes the eta parameters for a uniform grid in 2D as
%well as the desired virtual points to have a smaller eta. The eta parameter
%corresponds to 0th derivatives to deal with non-negativity constraints.
%This grid is much cruder than the CompXgapEta1D one as it not adaptive w.r.t.
%the original points.
X=Xorig;
m1=min(Xorig(:,1)); m2=min(Xorig(:,2)); M1=max(Xorig(:,1)); M2=max(Xorig(:,2));
if nargin<4
    nbPts_per_dim=10;
end
deltaDesi=(max(X)-min(X))/nbPts_per_dim;
xCoord=m1+deltaDesi(1)*(1/2+(0:nbPts_per_dim-1));
yCoord=m2+deltaDesi(2)*(1/2+(0:nbPts_per_dim-1));

[Xtemp,Ytemp] = meshgrid(xCoord,yCoord);
Xgap=[Xtemp(:),Ytemp(:)];

if any(max(diff(sort(X)))/2>sigX)
    disp("Some points were too far from the others, new points were added to fill the gap.")
end
if any(deltaDesi>sigX)
    disp(["The new points were still too far for", num2str(sigX), "the property will not be valid on the whole interval."])
%     deltaDesi(deltaDesi>sigX)=sigX/5;
end
Xnew=[X;Xgap]; 
eta=sqrt(2*h(0,sigX)-2*h(norm(deltaDesi/2),sigX));
if not(all(eta))
    disp("Warning: some eta_m=0")
end
end

