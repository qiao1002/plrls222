%%%%%%%%%%%%%%%%%%孔祥真所写的选基因的方法1
clc;
clear;
close all;
load('Data_3_Pollen.mat');
load true_labs;
X=in_X';
load U;
load Genes;
n=249;
selnum=200;%%selnum为所选基因的个数
selemethod=3;%%设置选择基因的方法
rowmax=max(abs(U)')';
if selemethod==1
 %%%%%%%%%%%%%%%%%%孔祥真所写的选基因的方法1

[junk, index]=sort(rowmax,'descend');
selmax=ones(selnum,n);
selectedgenes=cell(selnum,1);
for i=1:selnum
    selmax(i,:)=X(index(i),:);
    selectedgenes{i,:}=Genes(index(i),:);  %%此处记录的是基因名字
end
data=selmax';
save GLU data true_labs;
save  selectedgenes1_200 selectedgenes;
elseif selemethod==2
%%%%%%%%%%%%%%%%%孔祥真所写的选基因的方法2
u = sum(abs(U(:,1:2)),2);
[junk, index]=sort(u,'descend');
 selmax=ones(selnum,249);
selectedgenes=cell(selnum,1);
for i=1:selnum
    selmax(i,:)=X(index(i),:);
    selectedgenes{i,:}=Genes(index(i),:);
end
data=selmax';
save GLU data true_labs;
save  selectedgenes1_200 selectedgenes;
elseif selemethod==3
%%%%%%%%%%%%%%%%%%孔祥真所写的选基因的方法3
u=abs(U);
[junk, index]=sort(u,1,'descend');
index=index(1:selnum,:);
sel=index(:);
s=sortrows(tabulate(sel),-2);
index=s(1:selnum,:);
selmax=ones(selnum,n);
selectedgenes=cell(selnum,1);
for i=1:selnum
    selmax(i,:)=X(index(i),:);
    selectedgenes{i,:}=Genes(index(i),:);
end
data=selmax';
save GLU data true_labs;
save  selectedgenes1_200 selectedgenes;
end



[junk, index]=sort(rowmax,'descend');
selmax=ones(selnum,n);
selectedgenes=cell(selnum,1);
for i=1:selnum
    selmax(i,:)=X(index(i),:);
    selectedgenes{i,:}=Genes(index(i),:);
end
data=selmax';
save GLU data true_labs;
    
%%%%%%%%%%%%%%%将选出的基因写到磁盘文件里
%  filename = ['KselresultCRC_data.txt'];
% fid = fopen(filename, 'wt');
% 
% fprintf(fid,'number = %d\n', selnum);
% for i = 1 :selnum
%    fprintf(fid,'%s\n',selectedgenes{i,:});
% end
% fclose(fid);