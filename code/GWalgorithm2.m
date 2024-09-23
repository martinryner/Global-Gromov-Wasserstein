%This function takes the point clouds X and Y and finds the optimal GW
%distance and correspondance. This is just a wrapper to algorithm2

function [OptGamma,GWdist,GWL,GWU,iterations,extremepoints,Gammastars, Gammastarsobj] =GWalgorithm2(X,Y,relativeerror)

%Scaling the data a bit
Xnew = X;
Ynew = Y;
Xnew = Xnew- mean(Xnew')';%This doesn't affect the GW distance, but simplifies 1-norm calculations.
Ynew = Ynew- mean(Ynew')';
scale = max(sum(abs(Xnew),"all"),sum(abs(Ynew),"all"))/size(X,2);
Xnew = Xnew/scale;
Ynew = Ynew/scale;
tt=tic;
[OptGamma,Betsobjsofar,upper,GWU,GWL,Gammastars, Gammastarsobj,UpperXs,retmax,retmin,qptree,iterations] = algorithm2(Xnew,Ynew,relativeerror, scale);


if ~isempty(qptree)
extremepoints = length(qptree.objs);
else
extremepoints = 0;
end
time = toc(tt)
GWdist = GWU(end);









