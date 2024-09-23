%This script runs the algorithm on different test scenarios. Other
%comparative algorithms has been commented out and are not provided, e.g.
%Friesens low rank QAP and Peyrés Entropy regularized interative scheme.


clear all;
s = rng;
rng(1)


problems{1}.d = 2;
problems{1}.d2 = 2;
problems{1}.UN = 1;
problems{1}.tol = 1e-8;
problems{1}.ns = [10 100 500 1000];

problems{2}.d = 2;
problems{2}.d2 = 3;
problems{2}.UN = 1;
problems{2}.tol = 1e-8;
problems{2}.ns = [500];

problems{3}.d = 2;
problems{3}.d2 = 3;
problems{3}.UN = 2;
problems{3}.UNd = [1 ;1; 1];
problems{3}.tol = 1e-8;
problems{3}.ns = [100];

problems{4}.d = 3;
problems{4}.d2 = 3;
problems{4}.UN = 2;
problems{4}.UNd = [1; 1; 0.1];
problems{4}.tol = 1e-2;
problems{4}.ns = [500];%1000];% 100 500];


problems{5}.d = 3;
problems{5}.d2 = 3;
problems{5}.UN = 2;
problems{5}.UNd = [1; 0.5; 0.1];
problems{5}.tol = 1e-2;
problems{5}.ns = [2000];% 100 500 1000];

problems{6}.d = 3;
problems{6}.d2 = 3;
problems{6}.UN = 2;
problems{6}.UNd = [1; 1; 0.5];
problems{6}.tol = 1e-2;
problems{6}.ns = [20 40 ];% 100 500 1000];

problems{7}.d = 2;
problems{7}.d2 =3;
problems{7}.UN = 1;
problems{7}.UNd = [1; 1;1];
problems{7}.tol = 1e-6;
problems{7}.ns = [200 300 400];% 100 500 1000];

for pk =1:7
    disp(['running problem type number' num2str(pk)]);

for nk = 1:length(problems{pk}.ns)
n = problems{pk}.ns(nk);
 disp(['running problem size' num2str(n)])
for tk = 1:5
 disp(['running randomized trial ' num2str(tk)])
%problems{pk}
%n



d =problems{pk}.d; %If d << n it seems stable
if problems{pk}.UN ==2
    X = randn(d,n);
    standarddev = problems{pk}.UNd;
    X = X.*standarddev(1:d);
end
if problems{pk}.UN == 1
 X = 2*(rand(d,100*n)-0.5);
 nrm = sqrt(sum(X.^2,1)) <=1;
 X = X(:,nrm);
 X = X(:,1:n);
end

dd2 = problems{pk}.d2;


if problems{pk}.UN ==2
    Y = randn(dd2,n);
    standarddev = problems{pk}.UNd;
    Y = Y.*standarddev(1:dd2);
end
if problems{pk}.UN == 1
 Y = 2*(rand(dd2,100*n)-0.5);
 nrm = sqrt(sum(Y.^2,1)) <=1;
 Y = Y(:,nrm);
 Y = Y(:,1:n);
end

%X = X- min(X')';%The one norm works better then. This doesn't affect the GW distance
%Y = Y- min(Y')';

ett = ones(n,1);
if size(X,1)> 1
dx = vecnorm(X,2).^2';
else
dx = X.^2';
end
if size(Y,1) > 1
dy = vecnorm(Y,2).^2';
else
dy = Y.^2';
end

C = 2*n*dx*dy'-4*dx*((ett'*Y')*Y)-4*(X'*(X*ett))*dy';
%nCtrue = 4*X'*X*reshape(randP,n,n)*Y'*Y;
%trueobjective =(C(:)+nCtrue(:))'*randP(:) +2*sum(dx)*sum(dy);

d1 = pdist2(X',X','squaredeuclidean');
d2 = pdist2(Y',Y','squaredeuclidean');

t00 = tic;
%[G,opt] = friesenMILP1(X,Y,problems{pk}.tol);
G = [];
if ~isempty(G)
friesendistance = sum(d1.^2,"all")+sum(d2.^2,"all")-2*trace(d1*G*d2*G')
friesentime = toc(t00)
else
friesendistance = NaN;
friesentime = NaN; 
end

t0 = tic;
 %[gammap,a,b,GWp] = perform_gw_sinkhorn(d1,d2,ones(n,1),ones(n,1),00);
 peyretime = toc(t0);
peyredistance  = inf;
% peyredistance = sum(d1.^2,"all")+sum(d2.^2,"all")-2*trace(d1*gammap*d2*gammap')

 
 % mu1 = ones(n,1);
% mu2 = ones(n,1);
% for jj =1:0
% 
% gamma = rand(length(mu1),length(mu2));
% a = ones(length(mu1),1);
% b = ones(length(mu2),1);
% for k = 1:100
% b = mu2./(gamma'*a);
% a = mu1./(gamma*b);
% gg = a.*gamma.*b';
% %err = norm(sum(gg)-mu2)+norm(sum(gg')-mu1)
% end
% gamma = a.*gamma.*b';
% options = [];
% options.gamma_init = gamma;
%     [gammap1,a,b,GWp] = perform_gw_sinkhorn(d1,d2,ones(n,1),ones(n,1),00,options);
% peyredistance1 = sum(d1.^2,"all")+sum(d2.^2,"all")-2*trace(d1*gammap1*d2*gammap1');
% if peyredistance1 <peyredistance
% gammap = gammap1;
% peyredistance = peyredistance1
% end
% end
t0 = tic;
%[run,gd, obj,pool] =QPdirect(X,Y);
directt = toc(t0);
%directdistance = sum(d1.^2,"all")+sum(d2.^2,"all")-2*trace(d1*gd*d2*gd')
directt = NaN;
directdistance = NaN;

global gr;
gr = false;
t0= tic;
[gamma2,opt2,GWL,GWU,iterations,extemepoints] = GWalgorithm2(X,Y,problems{pk}.tol);
%gamma2 = eye(n);
GWU = 1;GWL=1;
%iterations=NaN;
%extemepoints=NaN
timealg2 = toc(t0)
%%%nC = 4*X'*X*reshape(gamma2,n,n)*Y'*Y;
algorithm2distance = sum(d1.^2,"all")+sum(d2.^2,"all")-2*trace(d1*reshape(gamma2,n,n)*d2*reshape(gamma2,n,n)')
alg2relgap = ((GWU(end)-GWL(end))/GWU(end));
alg2iters = length(GWL);


%test to run with Gurobi concave optimization. Commented out here
global gr;
 gr = true;
 t0= tic;
 %[gamma2,opt2] = GWalgorithm2(X,Y,problems{pk}.tol);%,gammap);
 timealg2gr = toc(t0)
 %%%%nC = 4*X'*X*reshape(gamma2,n,n)*Y'*Y;
 algorithm2distancegr = sum(d1.^2,"all")+sum(d2.^2,"all")-2*trace(d1*reshape(gamma2,n,n)*d2*reshape(gamma2,n,n)')
 algorithm2distancegr = NaN;
 timealg2gr  = NaN;



 %This is a test on Peyrés entropy regularized iterative scheme with
 %different starting points
t0 = tic;
xtilde = d1*ones(n,1);
ytilde = d2*ones(n,1);
Ctilde = pdist2(sqrt(xtilde),sqrt(ytilde),"squaredeuclidean");
[~,gamma] = mexEMD(ones(n,1),ones(n,1),Ctilde);
tottests = 1000;
if pk < 7
tottests = 1;
end
pey2iterations = 0;
best = inf;
for k = 1:tottests 
options = [];
options.gamma_init = gamma;
%[gammap1,a,b,GWp] = perform_gw_sinkhorn(d1,d2,ones(n,1),ones(n,1),00,options);
pey2iterations = NaN;
%pey2iterations = pey2iterations+1;%+ length(GWp);
 %peyredistance2 = sum(d1.^2,"all")+sum(d2.^2,"all")-2*trace(d1*gammap1*d2*gammap1');
 peyredistance2 = nan;
 [ (best -algorithm2distance)/algorithm2distance]
if peyredistance2 < best
    best = peyredistance2
end


    if (peyredistance2 -algorithm2distance)/algorithm2distance < problems{pk}.tol;
        break;
    end
    II = eye(n);
    gamma= II(randperm(n),:);
end



%peyredistance2 = inf;
peyretime2= toc(t0);
peyreiters2 = pey2iterations;

problems{pk}.friesentimes(tk,nk) = friesentime;
problems{pk}.friesendistances (tk,nk) = friesendistance ;

problems{pk}.alg2(tk,nk) = timealg2;
problems{pk}.alg2dist(tk,nk) = algorithm2distance;

problems{pk}.alg2extremepoints(tk,nk) = extemepoints;
problems{pk}.alg2gap(tk,nk) = alg2relgap;
problems{pk}.alg2gr(tk,nk) = timealg2gr;
problems{pk}.alg2distgr(tk,nk) = algorithm2distancegr;

problems{pk}.peytime(tk,nk) = peyretime;
problems{pk}.peydist(tk,nk) = peyredistance;
problems{pk}.peydist2(tk,nk) = peyredistance2;
problems{pk}.peytime2(tk,nk) = peyretime2;
problems{pk}.peyiters2(tk,nk) = peyreiters2;

problems{pk}.direct(tk,nk) = directdistance;
problems{pk}.directtime(tk,nk) = directt;
save problemstime.mat problems;

%save algo2.mat alg2;
end
end

end

return;

figure(3)
boxplot(100*((peytime)./alg2),labs,'PlotStyle','compact')

xlabel('problem size n')
ylabel('computational efficiency %')

figure(4)
boxplot(100*((peydist-alg2dist)./alg2dist),labs,'PlotStyle','compact')
xlabel('problem size n')
ylabel('relative error %')

figure(7)
boxplot(100*((friesentimes)./alg2),labs,'PlotStyle','compact')

xlabel('problem size n')
ylabel('computational efficiency')



return;



%Testing an image

coarsegrid = 3;

man=cast(imread("man.png"),'double');

man = squeeze(man(:,:,1));
man = filter2(ones(coarsegrid,coarsegrid),man);
man = man(1:coarsegrid:end,1:coarsegrid:end);

[x1,x2] = meshgrid(1:size(man,2),1:size(man,1));
X = [x1(man<255*coarsegrid*coarsegrid/2)';x2(man<255*coarsegrid*coarsegrid/2)';man(man<255*coarsegrid*coarsegrid/2)'/(255*coarsegrid*coarsegrid)'];

man1=cast(imread("man2.png"),'double');


man1 = squeeze(man1(:,:,1));
man1 = filter2(ones(coarsegrid,coarsegrid),man1);
man1 = man1(1:coarsegrid:end,1:coarsegrid:end);
[x1,x2] = meshgrid(1:size(man1,2),1:size(man1,1));
Y = [x1(man1<255*coarsegrid*coarsegrid/2)';x2(man1<255*coarsegrid*coarsegrid/2)';man1(man1<255*coarsegrid*coarsegrid/2)'/(255*coarsegrid*coarsegrid)'];


X=X(1:2,:);
Y=Y(1:2,:);



maxpoints = 20000;

minsize = min(size(X,2),size(Y,2))
points = min(maxpoints,minsize);
r1 =randperm(size(X,2));

Xnew = X(:,r1(1:points));
r2 =randperm(size(Y,2));
Ynew = Y(:,r2(1:points));

figure(2);
plot(Xnew(1,:),Xnew(2,:),'x')

figure(3);
plot(Ynew(1,:),Ynew(2,:),'x')



[OptGamma0,GWdist,GWL,GWU] =GWalgorithm2(Xnew,Ynew);
GWdistalg2 = GWdist
%[OptGamma0,GWdist] =algorithm4(Xnew/100,Ynew/100);
d1 = pdist2(Xnew',Xnew','squaredeuclidean');
d2 = pdist2(Ynew',Ynew','squaredeuclidean');
[gammap,a,b,GWdist] = perform_gw_sinkhorn(d1,d2,ones(points,1),ones(points,1),0);
peyredistance = sum(d1.^2,"all")+sum(d2.^2,"all")-2*trace(d1*gammap*d2*gammap')
alg2distance = sum(d1.^2,"all")+sum(d2.^2,"all")-2*trace(d1*OptGamma0*d2*OptGamma0')

GWdist
OptGamma = round(OptGamma0);
Xtrans = Xnew*reshape(OptGamma,points,points);

Ynew(1,:)= -Ynew(1,:);% = [-Ynew(1,:);Ynew(2,:)]%Ynew+[-100; 0];%Translate it

figure(5);
plot(Xtrans(1,:),Xtrans(2,:),'x')
figure(6);
plot(Ynew(1,:),Ynew(2,:),'x')
figure(40);
hold off;
iters = 100;

for iteration =1:iters
    Z = Xtrans*(iters-iteration)/iters+Ynew*iteration/iters;
    plot(Z(1,:),-Z(2,:),'x')
    axis tight
    % axis([min([Xtrans(1,:) Ynew(1,:)]) max([Xtrans(1,:) Ynew(1,:)]) min([-Xtrans(2,:) -Ynew(2,:)]) max([-Xtrans(2,:) -Ynew(2,:)])  ])
drawnow;
pause(0.1)

end

figure(21)
hold off;
for k = 1:points

    plot([Xtrans(1,k) Ynew(1,k)],[Xtrans(2,k) Ynew(2,k)]);

    hold on;
end

