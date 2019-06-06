
%  function [cid,nr] = kmcluster();
 function [] =km_kmcluster();
 disp('kmeans');
% CSKMEANS K-Means clustering - general method.
% 
% This implements the more general k-means algorithm, where 
% HMEANS is used to find the initial partition and then each
% observation is examined for further improvements in minimizing
% the within-group sum of squares.
%
% [CID,NR,CENTERS] = CSKMEANS(X,K,NC) Performs K-means
% clustering using the data given in X. 
% 
% INPUTS: X is the n x d matrix of data,
% where each row indicates an observation. K indicates
% the number of desired clusters. NC is a k x d matrix for the
% initial cluster centers. If NC is not specified, then the
% centers will be randomly chosen from the observations.
%
% OUTPUTS: CID provides a set of n indexes indicating cluster
% membership for each point. NR is the number of observations
% in each cluster. CENTERS is a matrix, where each row
% corresponds to a cluster center.
%
% See also CSHMEANS
% W. L. and A. R. Martinez, 9/15/01
% Computational Statistics Toolbox 

%[n,d] = size(x);
tic
global result;
global K;   %%%%%%%kmax
global S;
n=max(result(:,2));
nc=[];
% k=input('Input the initial cluster number£º','s');
% k=str2num(k);

for bl=2:K
% reply=input('Whether to specify the output clustering center£ºy/n ','s');
% if reply=='y'	
% 	i=1;
%  while i<=k	
%     	n1=input('','s');
% 	nc(i)=str2num(n1);
% 	i=i+1;
%  end
% else
%     ind = randperm(n);
% % We will add some noise to make it interesting.
% nc = ind(1:k)
% end

ind = randperm(n);
nc = ind(1:bl);  
%if nargin < 2
% Then pick some observations to be the cluster centers.

%end
% set up storage
% integer 1,...,k indicating cluster membership
cid = zeros(1,n); 
% Make this different to get the loop started.
%oldcid = ones(1,n);
% The number in each cluster.
nr = zeros(1,bl); 
% Set up maximum number of iterations.

for i=1:n
for j=1:bl
if i==nc(j)
sim(j)=1;
else
ind=find((result(:,1)==i & result(:,2)==nc(j)) | (result(:,1)==nc(j) & result(:,2)==i));
sim(j)=result(ind,6);
end
end
temp=find(sim==max(sim));
cid(i)=temp(1);
end
cid

maxiter = 10;
iter = 1;
ccid=[];
while ~isequal(cid,ccid) && iter < maxiter
% Implement the hmeans algorithm
% For each point, find the distance to all cluster centers
ccid=cid;
for i = 1:n
for j=1:bl
ind=find(cid==j);
nr(j)=length(ind);
for t=1:length(ind)
if i==ind(t)
sim(t)=1;
else
hang=find((result(:,1)==i & result(:,2)==ind(t))|(result(:,1)==ind(t) & result(:,2)==i));
sim(t)=result(hang,6);
end
end
ave(j)=mean(sim);
end
temp=find(ave==max(ave));
cid(i)=temp(1);
end
iter=iter+1;
end

disp([ ' divided into',num2str(bl), 'clusters']);
cid;
% iter
% nr

BWPsum=0; 
bl;
for kk=1:bl
kkk=find(cid==kk);                  
disp([ 'In',num2str(kk), 'the number of users  ', num2str(kkk)]);   
end
S;

for k1=1:bl
    F1=find(cid==k1);
    L1=length(F1);
    for kk1=1:L1 
        Lavg=1;
        clear avg2;
    for k2=1:bl           
        if k1==k2
%         disp('the same');
            continue;
        else
            sum2=0;    
            F2=find(cid==k2); 
            L2=length(F2);
            for kk2=1:L2
                kk11=double(F1(kk1));
                kk22=double(F2(kk2));
                simlari=S(kk11,kk22);
                sum2=sum2+S(kk11,kk22);
            end
                avg2(Lavg)=sum2/L2;
                Lavg=Lavg+1;
         end            
    end 
    kk1;
    bd(kk11)=min(avg2);
     
    sum1=0;
    for kk1p=1:L1
        if kk1p==kk1
            continue;
        else
            kkp1=double(F1(kk1));
            kkp2=double(F1(kk1p));
            sum1=sum1+S(kkp1,kkp2);      
        end    
            wd(kkp1)=sum1/(L1-1);
    end
        
    end
   
       
    
%     mmm
end
bd;
wd ;
n;
for p=1:n
    bwp(p)=(bd(p)-wd(p))/(bd(p)+wd(p));
end
bwp;
BWP(bl-1)=sum(bwp)/n;
end
BWP
[a b]=min(BWP);
Kopt=b+1;
disp([ 'Optimal clustering number:',num2str(Kopt)]); 
toc

% NK=length(kkk);
% for tt=1:NK
%     sum2=0;
%     KL=1;  
% for jj=1:K
%     sum1=0;
% if jj==kk 
% continue;
% else
% jjj=find(tmpidx==I(jj));
% NJ=length(jjj);
% for uu=1:NJ
% sum1=sum1+S(jjj(uu),kkk(tt));
% end
% end
% lbs(KL)=sum1/NJ;
% KL=KL+1;
% end

% bs=min(lbs);
% for ttt=1:NK
%     if ttt==tt
%         continue;
%     else
%         sum2=sum2+S(kkk(tt),kkk(ttt));
%     end   
% end
% 

% if NK==1 ws=S(kkk(1),kkk(1));
% else ws=sum2/(NK-1);
% end
% 

% BWP=(bs-ws)/(bs+ws);
% BWPsum=BWPsum+BWP;
% end
% 
% 
% BWPavg(bbl)=BWPsum/N;
% bbl=bbl+1;    
% %%%	tmpdpsim=sum(S(sub2ind([N N],notI,tmpidx(notI))));
% %%%	tmpexpref=sum(dS(I));
% %%%	tmpnetsim=tmpdpsim+tmpexpref;
% % else
% %     tmpidx=nan*ones(N,1); tmpnetsim=nan; tmpexpref=nan;
% % end;
% % end;
% BWPavg