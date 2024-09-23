
%This is a code implementation of the work in
% Ryner, Martin, Jan Kronqvist, and Johan Karlsson. 
% "Globally solving the Gromov-Wasserstein problem for point clouds in low dimensional Euclidean spaces." 
% Advances in Neural Information Processing Systems 36 (2023): 7930-7946.

% * This software is provided "AS IS" with no warranty of any kind,
% * express or implied, and with no claim as to its suitability for any
% * purpose.

function [highestlowerboundGamma,highestlowerbound,upperbound,GWupperbounds,GWlowerbounds,Gammastars, Gammastarsobj,UpperXs,retmax,retmin,qptree,k] =algorithm2distr(X,Y, mu,nu,relativeerror, scale,linearterm)
max_iterations = 50000;
dispinterval = 20;%how often should we display progress
highestlowerbound = -inf;
Gap = 100;


X = X-mean(X')';%The one norm works better then. This doesn't affect the GW distance
Y = Y-mean(Y')';


n1 = size(X,2);
n2 = size(Y,2);
d = size(X,1);
d2=size(Y,1);


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
GW0 = mu'*(CX.^2)*mu+ nu'*(CY.^2)*nu;

C = (sum(mu)+sum(nu))*dx*dy'-4*dx*((nu'*Y')*Y)-4*(X'*(X*mu))*dy';

if exist('linearterm','var')
    C=C-linearterm;
end

%Update the bounding box on W and w0

upperbound = inf;
Gammastars = spalloc(n1*n2,100,100*max(n1,n2));


disp('Calculating bounds');

Wmax = zeros(d,d2);
Wmin = zeros(d,d2);

for i = 1:d
    for j = 1:d2
        vv = zeros(d,d2);
        vv(i,j)=1;
        cc =2* vv(:)'*kron(Y,X);
        Wmax(i,j) = sum(max(reshape(cc,n1,n2)));
        [~,Gammastar] = mexEMD(mu,nu,-reshape(cc(:),n1,n2));
        %[~,Gammastar] = SinkhornOMT(mu,nu,-reshape(cc(:),n1,n2));
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
        Wmin(i,j) = sum(min(reshape(cc,n1,n2)));
        [~,Gammastar] = mexEMD(mu,nu,reshape(cc(:),n1,n2));
        %[~,Gammastar] = SinkhornOMT(mu,nu,reshape(cc(:),n1,n2));
        Wmin(i,j) = Gammastar(:)'*cc(:);
        CGammastar = C(:)'*Gammastar(:);
        lowerbound = sum((2*X*Gammastar*Y').*(2*X*Gammastar*Y'),"all")+CGammastar;
        if lowerbound > highestlowerbound
            highestlowerbound = lowerbound;
            highestlowerboundGamma = Gammastar;
        end
    end
end

[~,Gammastar] = mexEMD(mu,nu,reshape(C,n1,n2)/sum(abs(C),"all"));
%[~,Gammastar] = SinkhornOMT(mu,nu,reshape(C,n1,n2)/sum(abs(C),"all"));
w0min = Gammastar(:)'*C(:);
CGammastar = C(:)'*Gammastar(:);
lowerbound = sum((2*X*Gammastar*Y').*(2*X*Gammastar*Y'),"all")+CGammastar;
if lowerbound > highestlowerbound
    highestlowerbound = lowerbound;
    highestlowerboundGamma = Gammastar;
end

[~,Gammastar] = mexEMD(mu,nu,reshape(-C,n1,n2)/sum(abs(C),"all"));
% [~,Gammastar] = SinkhornOMT(mu,nu,reshape(-C,n1,n2)/sum(abs(C),"all"));
w0max = Gammastar(:)'*C(:);
CGammastar = C(:)'*Gammastar(:);
lowerbound = sum((2*X*Gammastar*Y').*(2*X*Gammastar*Y'),"all")+CGammastar;
if lowerbound > highestlowerbound
    highestlowerbound = lowerbound;
    highestlowerboundGamma = Gammastar;
end



retmax = [Wmax(:);w0max];
retmin = [Wmin(:);w0min];

disp('Done calculating bounds')
Ac = [];
bc = [];

Q = speye(d*d2+1);
Q(d*d2+1,d*d2+1)=0;


qptree = [];



x = zeros(d*d2+1,1);%Dummy
Gammastar0 = spalloc(n1,n2,0);

mobmaxqt = 0;







for k = 1:max_iterations






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



        upperbound = x'*Q*x+x(end);


    if opt == 0
        %couldnt find optimum
        break;
    end





    UpperXs(k,:) = x';
    Wstar = reshape(x(1:d*d2),d,d2);
    w0star = x(end);


  


        % Find Gammastar on convex hull with gradient at upper bound


        Cnew = 4*X'*Wstar*Y;
        Cnew = Cnew(:)+C(:);
        [~,Gammastar] = mexEMD(mu,nu,reshape(-Cnew,n1,n2)/sum(abs(Cnew),"all") );
        %[~,Gammastar] = SinkhornOMT(mu,nu,reshape(-Cnew,n1,n2)/sum(abs(Cnew),"all") );

        Gammastar = reshape(Gammastar,n1,n2);

        lowerbound = sum((2*X*Gammastar*Y').*(2*X*Gammastar*Y'),"all")+C(:)'*Gammastar(:);




        %Add constraint ((2Wstar,1)| (W-2XG*Y^T,w0-<L,G*>)) <= 0
        XGY2 = 2*X*Gammastar*Y';
        CGammastar = C(:)'*Gammastar(:);

        acnew = [2*(Wstar(:))' 1];
        bcnew =(2*Wstar(:)'*XGY2(:)+1*CGammastar);
        nrm = norm(acnew);
        acnew= acnew/nrm;
        bcnew = bcnew/nrm;




        %Add it if new, else break,
        if ~isempty(Ac)
            ss = sum(abs(Ac-acnew),2)+abs(bc-bcnew)';
        else
            ss = 1000;
        end


        if sum(ss < sqrt(eps)) > 0
            break;

        end


        Ac(end+1,:) = acnew;
        bc(end+1) = bcnew;



     

        lowerbound = sum((2*X*Gammastar*Y').*(2*X*Gammastar*Y'),"all")+CGammastar;
      
        if lowerbound >= highestlowerbound

            highestlowerbound = lowerbound;
            highestlowerboundGamma = Gammastar;
        end

  
 

    GW1 = GW0*scale^4-2*(highestlowerbound+2*(mu'*dx)*(nu'*dy) )*scale^4;
    GW2 = GW0*scale^4-2*(upperbound+2*(mu'*dx)*(nu'*dy) )*scale^4;

    Gap = ((GW1-GW2)/GW1);


    if mod(k,dispinterval) == 0
        treesize = 0;
        if ~isempty(qptree)
            treesize = length(qptree.objs);
        end


        disp(['GW Upper bound ' num2str(GW1) ' GW Lower bound ' num2str(GW2) ' Gap ' num2str(Gap) ' tree size ' num2str(treesize) ]  );

   
    end
    GW2 = max(0,GW2);

    GWupperbounds(k) =GW1;
    GWlowerbounds(k) = GW2;
    Gammastars(:,k) = sparse(Gammastar(:));
    Gammastarsobj(k) = GW2;

  


    if( Gap <relativeerror)  && k>3
        break;
    end
end

treesize = 0;
if ~isempty(qptree)
    treesize = length(qptree.objs);
end


disp(['GW Upper bound ' num2str(GW1) ' GW Lower bound ' num2str(GW2) ' Gap ' num2str(Gap) ' tree size ' num2str(treesize) ]  );

highestlowerbound=highestlowerbound+ 2*sum(dx)*sum(dy);



end




