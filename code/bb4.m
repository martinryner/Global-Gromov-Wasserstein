%This is the extreme point method, which calculates the new vertices and
%holds the v-space and h-space.
%Input: qptree is the struct that holds the h-space and v-space. Starting
%out it is empty and we create the first structure with the lb (lower
%bound) and ub (upper bound) box. Aend and bend are newly added halfspaces
%Output: qptree is an updated struct when a new half plane was added or a
%new one from the box

function qptree = bb4(Aend,bend,lb,ub,qptree)

if ~exist('qptree','var') || isempty(qptree)
    qptree.EPs = [];
    qptree.lbbases = [];

    qptree.bases = [];
    qptree.objs = [];
    qptree.closematrix = [];
    qptree.used = [];
    qptree.A = [];
    qptree.b = [];

end

EPs = qptree.EPs;
lbbases = qptree.lbbases;

bases = qptree.bases;
objs = qptree.objs;
closematrix = qptree.closematrix;
used = qptree.used;

A = qptree.A;
b = qptree.b;



if isempty(EPs)
    [EPs,lbbases,ubbases,objs] = BoxEps(lb,ub);
    bases = spalloc(1,length(lbbases),0);
    bases = cast(bases,'logical');
    tmpvector = [lbbases;bases];
    closematrix = sparse((tmpvector'*tmpvector==(length(lb)-1)));
    qptree.closematrix = closematrix;
    used= zeros(length(objs),1);
    qptree.used = used;
end



Aew = [A;Aend];
A = Aew;

bew = [b;bend];
b = bew;

%Which extreme points are no longer feasible?
if ~isempty(Aend)
   
    if size(bases,1) < size(A,1)
        bases(size(A,1),1) = 0;
    end

    %Which exremepoints are now infeasible
    infeasiblesbit = A(end,:)*EPs-b(end) >sqrt(eps) ;

    infeasibles= find (infeasiblesbit);
    feasiblebit = ~infeasiblesbit;

    %closematrix contains all points that are adjacent by size(A,2)-1
    %constraints. feasiblebit is a binary vector with feasible points. The
    %hadamardproduct is the feasible adjacent points to every point.
    feasibleadjacents= closematrix(:,infeasiblesbit)&feasiblebit';%for every column k denoting an extreme point, the 1s in column k denotes a feasible adjacent extreme point to point k



    %We reduce the indexing to that we only have the feasible vectors left,
    %the otherones will be removed soon anyway
    reducedfeasibleadjacents =closematrix(feasiblebit,:);


    bend = b(end);
    Aend = A(end,:);

    closematrix1columns = cell(length(infeasibles),1);
    increment = 0;

    newlbbases = cell(1,length(infeasibles));
    newbases = cell(1,length(infeasibles));
    newEPs = cell(1,length(infeasibles));
    newobjs = cell(1,length(infeasibles));


    %run though the infeasible points and create new ones by interpolating
    %to their feasible adjacent ones.
    for tk = 1:length(infeasibles)
        k = infeasibles(tk);



        lookup = feasibleadjacents(:,tk);%was k but with an unreduced feasibleadjacents


        if sum(lookup) > 0
            increment = increment+1;



            lbbit= lbbases(:,k);
            % ubbit = ubbases(:,k);
            Abit = bases(:,k);



            used(lookup) = used(lookup)+1;

            %Do them all at once
            Ls = (bend-Aend*EPs(:,lookup))./(Aend*(EPs(:,k)-EPs(:,lookup)));
            newEP = Ls.*EPs(:,k)+(1-Ls).*EPs(:,lookup);
            newEPs{tk} =newEP;

            newlbbases{tk} = lbbit&lbbases(:,lookup);
           
            newAbits = Abit&bases(:,lookup);
            newAbits(end,:) = 1;
            newbases{tk} =newAbits;


            tmp = newEP.^2;
            newobjs{tk} = sum(tmp(1:end-1,:),1)+newEP(end,:);

            %Test to create closematrix1 below
            closematrix1columns{tk} = find(reducedfeasibleadjacents(:,k))';
           
        end
    end

  
    %Remove the infeasibles from the struct
    EPs(:,infeasibles) = [];
    lbbases(:,infeasibles) = [];
    %  ubbases(:,infeasibles) = [];
    bases(:,infeasibles) = [];
    objs(infeasibles) = [];
    closematrix(:,infeasibles) = [];
    closematrix(infeasibles,:) = [];

    used(infeasibles) = [];

   
    origsz = length(objs);
  

    %collect closematrix1, adjacency from existing points to the new points
    closematrix1columnsv = [closematrix1columns{:}];

   
    %Save them in sequential order. Pre allocation to save time




    lbbases = [lbbases [newlbbases{:}]];
  
    bases = [bases [newbases{:}]];
    EPs = [EPs [newEPs{:}]];
    objs = [objs [newobjs{:}]];
    used(length(objs)) =0;


   
    %pad with zeros if not the correct size
    if size(bases,2) <length(objs)
        bases(1,length(objs)) = false;
    end
   
    if size(lbbases,2) <length(objs)
        lbbases(1,length(objs)) = false;
    end







    %Build closematrix2, can get memory problems, we build it by chunks
    %This one can get large. We can build it by chunks
 
    tmpvector1 = [[newlbbases{:}];[newbases{:}]];

    rbit = (sum(tmpvector1,2)>0); %Reduce it to only work with necessary non-zero elements
    
    tmpvector1 = tmpvector1(rbit,:);


  

    %The original slow one
    % closematrix2o =((tmpvector1'*tmpvector1==(size(A,2)-1)));


    %A better but slow
    SZM = 1000;
    as = cell(ceil(size(tmpvector1,2)/SZM));
    bs = cell(ceil(size(tmpvector1,2)/SZM));

    t1k = 0;


    for t1 = 1:SZM:size(tmpvector1,2)
        t1k =t1k+1;
        len1 = SZM;
        if t1+SZM > size(tmpvector1,2)
            len1 = size(tmpvector1,2)-t1;
        end
        t2k = 0;

        for t2 = 1:SZM:t1
            t2k = t2k+1;
            len2 = SZM;
            if t2+SZM > size(tmpvector1,2)
                len2 = size(tmpvector1,2)-t2;
            end



            [ta, tb, ~] =find( ( tmpvector1(:,t1:t1+len1)' * tmpvector1(:,t2:t2+len2) ==(size(A,2)-1)  ) );
            as{t1k,t2k} = ta(:)'+t1-1;
            bs{t1k,t2k} = tb(:)'+t2-1;

          
        end
    end
    taf = [as{:}];
    tbf = [bs{:}];
   



    %The old one build on the old sparse matrix.
    [closematrixrows,closematrixcols,~] = find(closematrix);
    %Old to new
    cm1b = 1:length(closematrix1columnsv);
    cm1a = closematrix1columnsv;
    %new to new
    cm2a = [taf  tbf];
    cm2b = [tbf  taf];


   %Build the new adjacency matrix
    closematrixrows = [closematrixrows' cm1a (cm1b+origsz ) (cm2a+origsz )];
    closematrixcols = [closematrixcols' (cm1b+origsz) cm1a (cm2b+origsz )];
    %Dont calculate it
    closematrix = sparse(closematrixrows,closematrixcols,true);


  
end




qptree.used = used;
qptree.EPs = EPs;
qptree.lbbases=lbbases;

qptree.bases=bases;
qptree.objs=objs;
qptree.closematrix = closematrix;

qptree.A = A;
qptree.b = b;
end



function [EPs,lbbases,ubbases,objs] = BoxEps(lb,ub)
EPs = zeros(length(lb),2^length(lb));


lbbases = spalloc(length(lb),2^length(lb),100);
ubbases = spalloc(length(lb),2^length(lb),100);


for k = 0:2^length(lb)-1


    bit= bitget(k,1:length(lb))';
    EPs(:,k+1) = (1-bit).*lb+bit.*ub;
    lbbases(:,k+1) = sparse(1-bit);
    ubbases(:,k+1) = sparse(bit);
end

tmp = EPs.^2;
%tmp(end,:) = [];
objs = sum(tmp(1:end-1,:),1)+EPs(end,:);
lbbases = cast(lbbases,'logical');
ubbases = cast(ubbases,'logical');

%Just pack them together, its pointless to keep them separate
lbbases = [lbbases;ubbases];
ubbases =[];%spalloc(1,size(lbbases,2),0);

end

