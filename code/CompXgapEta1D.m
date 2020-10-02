function [Xnew,eta] = CompXgapEta1D(Xorig,sigX,h,nbPts_per_dim)
%This function computes the eta parameters for a non-uniform grid in 1D as
%well as the desired virtual points to have a smaller eta. The eta parameter
%corresponds to 0th derivatives to deal with non-negativity constraints.
X=Xorig(:,1);
if nargin<4
    nbPts_per_dim=50;
end
deltaDesi=(max(X)-min(X))/nbPts_per_dim;
if (max(diff(sort(X)))/2>sigX)
    disp("Some points were too far from the others, new points were added to fill the gap.")
end
if (deltaDesi>sigX)
    disp("The new points were still too far, the property will not be valid on the whole interval.")
%     deltaDesi(deltaDesi>sigX)=sigX/5;
end
nb_point_to_add=floor(diff(X)/deltaDesi); ngap=sum(nb_point_to_add);
Xgap=zeros(ngap,size(X,2)); count=1;
for i=1:length(nb_point_to_add)
    for j=1:nb_point_to_add(i)
        Xgap(count,:)=X(i,:)+(X(i+1,:)-X(i,:))*j/(nb_point_to_add(i)+1);
        count=count+1;
    end
end
Xnew=[X;Xgap]; [Xnew_sort,Inew]=sort(Xnew); dX=max([diff(Xnew_sort);0],[0;diff(Xnew_sort)])/2;
[~,Inew]=sort(Inew); dX=dX(Inew); dX(dX==0)=min(dX(dX~=0));
eta=sqrt(2*h(0,sigX)-2*h(dX,sigX));
if not(all(eta))
    disp("Warning: some eta_m=0")
end
end

