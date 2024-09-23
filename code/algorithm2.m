%This is a code implementation of the work in
% Ryner, Martin, Jan Kronqvist, and Johan Karlsson. 
% "Globally solving the Gromov-Wasserstein problem for point clouds in low dimensional Euclidean spaces." 
% Advances in Neural Information Processing Systems 36 (2023): 7930-7946.

% * This software is provided "AS IS" with no warranty of any kind,
% * express or implied, and with no claim as to its suitability for any
% * purpose.


function [highestlowerboundGamma,highestlowerbound,upperbound,GWupperbounds,GWlowerbounds,Gammastars, Gammastarsobj,UpperXs,retmax,retmin,qptree,k] =algorithm2(X,Y, relativeerror, scale,linearterm)
max_iterations = 50000;
dispinterval = 10;%how often should we display progress
highestlowerbound = -inf;
Gap = 100;

X = X-mean(X')';%The one norm works better then. This doesn't affect the GW distance
Y = Y-mean(Y')';


n = size(X,2);
d = size(X,1);
d2=size(Y,1);
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

CX = pdist2(X',X',"squaredeuclidean");
CY = pdist2(Y',Y',"squaredeuclidean");
GW0 = sum(CX.^2,"all")+ sum(CY.^2,"all");

C = 2*n*dx*dy'-4*dx*((ett'*Y')*Y)-4*(X'*(X*ett))*dy';

if exist('linearterm','var')
    C=C-linearterm;
end

%Update the bounding box on W and w0






upperbound = inf;
Gammastars = spalloc(n*n,100,100*n);


disp('Calculating bounds');

Wmax = zeros(d,d2);
Wmin = zeros(d,d2);

for i = 1:d
    for j = 1:d2
        vv = zeros(d,d2);
        vv(i,j)=1;
        cc =2* vv(:)'*kron(Y,X);
        Wmax(i,j) = sum(max(reshape(cc,n,n)));
        [~,Gammastar] = mexEMD(ones(n,1),ones(n,1),-reshape(cc(:),n,n));
        Wmax(i,j) = Gammastar(:)'*cc(:);
        CGammastar = C(:)'*Gammastar(:);
        lowerbound = sum((2*X*Gammastar*Y').*(2*X*Gammastar*Y'),"all")+CGammastar;
        if lowerbound > highestlowerbound
            highestlowerbound = lowerbound;
            highestlowerboundGamma = Gammastar;
        end
    end
end


for i = 1:d
    for j = 1:d2
        vv = zeros(d,d2);
        vv(i,j)=1;
        cc = 2*vv(:)'*kron(Y,X); %index * vec(Z) = 2*index * kron(Y,X)vec(Gamma)
        Wmin(i,j) = sum(min(reshape(cc,n,n)));
        [~,Gammastar] = mexEMD(ones(n,1),ones(n,1),reshape(cc(:),n,n));
        Wmin(i,j) = Gammastar(:)'*cc(:);
        CGammastar = C(:)'*Gammastar(:);
        lowerbound = sum((2*X*Gammastar*Y').*(2*X*Gammastar*Y'),"all")+CGammastar;
        if lowerbound > highestlowerbound
            highestlowerbound = lowerbound;
            highestlowerboundGamma = Gammastar;
        end
    end
end

[~,Gammastar] = mexEMD(ones(n,1),ones(n,1),reshape(C,n,n)/sum(abs(C),"all"));
w0min = Gammastar(:)'*C(:);
CGammastar = C(:)'*Gammastar(:);
lowerbound = sum((2*X*Gammastar*Y').*(2*X*Gammastar*Y'),"all")+CGammastar;
if lowerbound > highestlowerbound
    highestlowerbound = lowerbound;
    highestlowerboundGamma = Gammastar;
end

[~,Gammastar] = mexEMD(ones(n,1),ones(n,1),reshape(-C,n,n)/sum(abs(C),"all"));
w0max = Gammastar(:)'*C(:);
CGammastar = C(:)'*Gammastar(:);
lowerbound = sum((2*X*Gammastar*Y').*(2*X*Gammastar*Y'),"all")+CGammastar;
if lowerbound > highestlowerbound
    highestlowerbound = lowerbound;
    highestlowerboundGamma = Gammastar;
end



retmax = [Wmax(:);w0max];
retmin = [Wmin(:);w0min];

Ac = [];
bc = [];

Q = speye(d*d2+1);
Q(d*d2+1,d*d2+1)=0;


qptree = [];
qptreebak = [];



for k = 1:max_iterations

    qptreebak2 = qptreebak;
    qptreebak = qptree;


    %This is just at test, if the global variable gr is true, we solve the
    %low rank concave problem with Gurobi, which is really slow when the
    %number of constraints increase
    global gr;
    if isempty(gr)
        gr = false;
    end

    if gr

        [opt, x,~] = QP1(Q,[zeros(d*d2,1);1],Ac,bc,[Wmin(:);w0min],[Wmax(:);w0max],'max',size(Ac,1),highestlowerbound,1e-16);

    else



        if ~isempty(Ac)
            

            qptree = bb4(Ac(end,:),bc(end)',[Wmin(:);w0min],[Wmax(:);w0max],qptree);


        else
            qptree = bb4([],[],[Wmin(:);w0min],[Wmax(:);w0max],qptree);
        end

        [~,tin] = max(qptree.objs);
        x = qptree.EPs(:,tin);



        opt = 1;
        if isempty(x)
            opt = 0;
        end

    end

    if opt == 0
        %couldnt find optimum
        break;
    end


    upperbound = x'*Q*x+x(end);
    UpperXs(k,:) = x';
    Wstar = reshape(x(1:d*d2),d,d2);
    w0star = x(end);


    %Test to use a least suqare fitting of the cutting plane, not so
    %effective
    leastsquarefit = false;
    if leastsquarefit
        %Least square fitting
%         cvx_begin quiet
%         variable gamt(n,n)
%         minimize( norm([reshape(2*X*gamt*Y', d*d2,1);C(:)'*gamt(:)]-x))
%         subject to
%         gamt*ones(n,1) == ones(n,1)
%         ones(1,n)*gamt == ones(1,n)
%         gamt>=0
%         cvx_end

        [opt, gamt,~] = QPdirect(X,Y,C(:),x);


        gvv = 2*X*gamt*Y';
        gvv = [gvv(:); C(:)'*gamt(:)];
        riktning = x(:)'-gvv';

        Cnew = 2*X'*reshape(riktning(1:d*d2),d,d2)*Y;
        Cnew = Cnew(:)+riktning(end)*C(:);
        %Cnew = 4*X'*Wstar*Y;
        %Cnew = Cnew(:)+C(:);
        [~,Gammastar] = mexEMD(ones(n,1),ones(n,1),reshape(-Cnew,n,n)/sum(abs(Cnew),"all") );


        newA = riktning;

        newb = 2*X*gamt*Y';
        newb = riktning*[newb(:);C(:)'*gamt(:)];

        nrm = norm(newA);
        newA= newA/nrm;
        newb = newb/nrm;


        Ac(end+1,:) = newA;
        bc(end+1) = newb;
        
        %Make sure the old optuimum is infeasible
        if newA*x< newb
    disp('Something is wrong. ')
        end


        CGammastar = C(:)'*Gammastar(:);
        lowerbound = sum((2*X*Gammastar*Y').*(2*X*Gammastar*Y'),"all")+CGammastar;
        if lowerbound >= highestlowerbound
            highestlowerbound = lowerbound;
            highestlowerboundGamma = Gammastar;

        end


    else



        % Find Gammastar on convex hull with gradient at upper bound

       
        Cnew = 4*X'*Wstar*Y;
        Cnew = Cnew(:)+C(:);
        [~,Gammastar] = mexEMD(ones(n,1),ones(n,1),reshape(-Cnew,n,n)/sum(abs(Cnew),"all") );

        Gammastar = reshape(Gammastar,n,n);

        
        %Add constraint ((2Wstar,1)| (W-2XG*Y^T,w0-<L,G*>)) <= 0
        XGY2 = 2*X*Gammastar*Y';
        CGammastar = C(:)'*Gammastar(:);

        acnew = [2*(Wstar(:))' 1];
        bcnew =(2*Wstar(:)'*XGY2(:)+1*CGammastar);
        nrm = norm(acnew);
        acnew= acnew/nrm;
        bcnew = bcnew/nrm;


        %Add it if it is a new cut , else break,
        if ~isempty(Ac)
            ss = sum(abs(Ac-acnew),2)+abs(bc-bcnew)';
        else
            ss = 1000;
        end


        if sum(ss <sqrt(eps)) == 0
            Ac(end+1,:) = acnew;
            bc(end+1) = bcnew;
         
        else
            break;%It was too close to the existing ones
        end


        lowerbound = sum((2*X*Gammastar*Y').*(2*X*Gammastar*Y'),"all")+CGammastar;
        if lowerbound >= highestlowerbound
            highestlowerbound = lowerbound;
            highestlowerboundGamma = Gammastar;
        end

    end

    GW1 = GW0*scale^4-2*(highestlowerbound+2*sum(dx)*sum(dy))*scale^4;
    GW2 = GW0*scale^4-2*(upperbound+2*sum(dx)*sum(dy))*scale^4;

    Gap = ((GW1-GW2)/GW1);


    if mod(k,dispinterval) == 0
        treesize = 0;
        if ~isempty(qptree)
            treesize = length(qptree.objs);
        end

        disp(['GW Upper bound ' num2str(GW1) ' GW Lower bound ' num2str(GW2) ' Gap ' num2str(Gap) ' tree size ' num2str(treesize)]  );
    end
    GW2 = max(0,GW2);

    GWupperbounds(k) =GW1;
    GWlowerbounds(k) = GW2;
    Gammastars(:,k) = sparse(Gammastar(:));
    Gammastarsobj(k) = GW2;

    plotd =false;
    if plotd
        if d+d2==2
            figure(4040)
            hold off;
            plot([w0min w0max],[Wmin Wmin])
            hold on;
            plot([w0min w0max],[Wmax Wmax])
            plot([w0min w0min],[Wmin Wmax])
            plot([w0max w0max],[Wmin Wmax])
            for kk = 1:size(Ac,1)
                plot([w0min w0max], [bc(kk)/Ac(kk,1)-Ac(kk,2)/Ac(kk,1)*w0min bc(kk)/Ac(kk,1)-Ac(kk,2)/Ac(kk,1)*w0max] );
            end


            plot(w0star,Wstar,'x');
            plot(C(:)'*Gammastar(:),2*X*Gammastar*Y','s')

            axis([w0min w0max Wmin Wmax])
        end
    end


    if( Gap <relativeerror)  && k>3
        break;
    end
end

treesize = 0;
if ~isempty(qptree)
    treesize = length(qptree.objs);
end
disp(['GW Upper bound ' num2str(GW1) ' GW Lower bound ' num2str(GW2) ' Gap ' num2str(Gap) ' tree size ' num2str(treesize)]  );

highestlowerbound=highestlowerbound+ 2*sum(dx)*sum(dy);



end


function [opt, x,obj] = QP1(Q,C,A,b,lb,ub,minmax,lessthansize,objmax,tolerance)
clear model;
%if strcmp(minmax,'max')


params.NonConvex = 2;
model.obj = C(:);
model.Q = Q;
if ~isempty(A)
    model.sense= cast(['<'*ones(lessthansize,1);'='*ones(size(A,1)-lessthansize,1) ],'char');
    model.A = sparse(A);
    model.rhs =b;
else
    model.sense= cast(['='*ones(1,1)],'char');
    model.A = sparse(zeros(1,length(lb)));
    model.rhs = 0;
end
if ~isempty(lb)

    model.lb = lb;
end
if ~isempty(ub)
    model.ub = ub;
end
model.modelsense= minmax;
params.outputflag = 0;

params.MIPGap =1e-16;%tolerance; %1e-16;
params.MIPGapAbs =1e-16;%tolerance;% 1e-20;
params.PoolSearchMode=1;
params.PoolSolutions = 50;

params.MIPGap =1e-20;
params.MIPGapAbs =1e-20;

model.quadcon(1).Qc = -Q;
model.quadcon(1).q = -C;
model.quadcon(1).rhs = -objmax;


result = gurobi(model, params);

if strcmp(result.status,'INFEASIBLE') ==1

    x = [];
    obj = [];


    opt = false;

else
    x = result.x;

    obj = result.objval;
    opt = true;

end
end

function [opt, x,obj,pool] = QPdirect(X,Y,C,x0)
clear model;

%if strcmp(minmax,'max')
n = size(X,2);
d1 = size(X,1);
d2 = size(Y,1);

A0 = [kron(speye(n),ones(1,n))  spalloc(n,d1*d2+1,0) ;kron(ones(1,n),speye(n)) spalloc(n,d1*d2+1,0)];
A1 = [2*kron(Y,X) -eye(d1*d2) zeros(d1*d2,1)];
A2 = [C(:)' zeros(1,d1*d2) -1];
qv = [zeros(n*n,1);ones(d1*d2+1,1)];

model.A = sparse([A0;A2;A1]);
model.rhs = [ones(2*n,1); zeros(d1*d2+1,1)];
model.sense= cast(['='*ones(length(model.rhs),1)],'char');
model.obj = [zeros(n^2,1);-2*x0];
model.Q = sparse(1:(n^2+d1*d2+1),1:(n^2+d1*d2+1),qv);

%params.NonConvex = 2;
model.lb = [zeros(n^2,1);-1000000000*ones(d1*d2+1,1)];
model.ub = [ones(n^2,1);1000000000*ones(d1*d2+1,1)];



model.modelsense= 'min';
params.outputflag = 0;
params.MIPGap =1e-2;
%params.MIPGapAbs =1e-20;
params.PoolSearchMode=1;
params.PoolSolutions = 50;
%params.OutputFlag = 1;


result = gurobi(model, params);



if strcmp(result.status,'INFEASIBLE') ==1
    % disp("empty set")
    x = [];
    obj = [];
    y = [];
    s = [];
    primalslack = [];
    opt = false;
    pool = [];
else
    x =  reshape(result.x(1:n^2),n,n);;

    obj = result.objval;
    opt = true;
    %    pool = result.pool;
    pool = [];

end
end




