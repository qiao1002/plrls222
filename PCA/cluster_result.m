%%计算聚类结果的评价指标ACC，NMI，PURITY
% clc;
% clear;
% load 'V.mat';
% load 'class.mat';
data=V;
gnd=true_labs;
nclass = length(unique(gnd));
iter=20;
acc=zeros(1,20);
nmi=zeros(1,20);
pur=zeros(1,20);
for i=1:iter
  idx=kmeans(data,nclass);
  pre_labels=idx;
  NMI=Cal_NMI_newused(gnd, pre_labels);
%   Purity = purity(max(gnd),result_label, gnd);
  ARI=Contingency_ARI_newused(gnd, pre_labels);
%   clusterresult= ClusteringMeasure(gnd, pre_labels);
  nmi(i)=NMI;
  ari(i)=ARI;
%   acc(i)=clusterresult(1);
%   nmi(i)=clusterresult(2);
%   pur(i)=clusterresult(3);
end
ari_mn=mean(ari);
ari_max=max(ari);
ari_std=std(acc);
nmi_mn=mean(nmi);
nmi_max=max(nmi);
nmi_std=std(nmi);
pur_mn=mean(pur);
pur_std=std(pur);
   