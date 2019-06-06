[fname,dirpath]=uigetfile('*.txt');
openfile=[dirpath fname];
%tic;
A=dlmread(openfile);
totalusers=max(A(:,1));
global result;
result=[];
for n1=1:totalusers
user1=find(A(:,1)==n1);
for n2=n1+1:totalusers
user2=find(A(:,1)==n2);
size1=size(user1,1);
size2=size(user2,1);
i=1;
B=[];
while i<=size1
B1=A(user1(i),:);
B=[B;B1];
i=i+1;
end
j=1;
C=[];
while j<=size2
C1=A(user2(j),:);
C=[C;C1];
j=j+1;
end
%Calculate the similarity between B and C matrices
i=1;
comm=0;
grade=0;
while i<=size1
k1=B(i,2);
k2=find(C(:,2)==k1);
if size(k2,1)==0
i=i+1;
continue;
else
comm=comm+1;
gradeC=C(k2(1),3);
gradeB=B(i,3);
grade=grade+abs(gradeB-gradeC);
i=i+1;
end
end

sims=comm/min(size1,size2);
if comm==0
    GD=0;
else
    GD=grade/(comm*4);
end
simp=sims/exp(GD);

midresult=[n1 n2 comm sims GD simp max(size1,size2) min(size1,size2) max(size1,size2)-min(size1,size2)];
result=[result;midresult];
end
 end
result1=result';
fid=fopen('e:/similarity.txt','wt');
fprintf(fid,'%d       %d       %d       %12.8f       %12.8f       %12.8f       %d       %d       %d\n',result1);
fclose(fid);
%toc;
