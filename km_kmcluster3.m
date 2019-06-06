
function [cid,nr] = km_kmcluster()
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

%[n,d] = size(x);
global result;
n=max(result(:,2));
nc=[];
k=input('Input the initial clustering number:','s');
k=str2num(k);
reply=input('Whether to specify the output clustering center:y/n ','s');
if reply=='y'	
	i=1;
 while i<=k	
    	n1=input('','s');
	nc(i)=str2num(n1);
	i=i+1;
 end
else
    ind = randperm(n);
% We will add some noise to make it interesting.
nc = ind(1:k)
end

    

%if nargin < 2
% Then pick some observations to be the cluster centers.

%end
% set up storage
% integer 1,...,k indicating cluster membership
cid = zeros(1,n); 
% Make this different to get the loop started.
%oldcid = ones(1,n);
% The number in each cluster.
nr = zeros(1,k); 
% Set up maximum number of iterations.

for i=1:n
for j=1:k
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
for j=1:k
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
iter
nr
end