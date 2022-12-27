clear all
% 基因表达数据载入
load(strcat('E:\哈工大交流\投稿论文－1\dataset\CRC',filesep,'CRC_GE','.mat'));
load(strcat('E:\哈工大交流\投稿论文－1\dataset\PAAD',filesep,'PAAD_GE','.mat'));
load(strcat('E:\哈工大交流\投稿论文－1\dataset\ESCA',filesep,'ESCA_GE','.mat'));
load(strcat('E:\哈工大交流\投稿论文－1\dataset\HNSC',filesep,'HNSC_GE','.mat'));
load(strcat('E:\哈工大交流\投稿论文－1\dataset\CHOL',filesep,'CHOL_GE','.mat'));
% 样本标签
load(strcat('E:\哈工大交流\投稿论文－1\dataset\HNSC',filesep,'HNSC_samplecategory','.mat'));
HNSC_sample = SampleCategory;
load(strcat('E:\哈工大交流\投稿论文－1\dataset\PAAD',filesep,'PAAD_samplecategory','.mat'));
PAAD_sample = Samplecategory;
load(strcat('E:\哈工大交流\投稿论文－1\dataset\ESCA',filesep,'ESCA_samplecategory','.mat'));
ESCA_sample = SampleCategory;
load(strcat('E:\哈工大交流\投稿论文－1\dataset\CRC',filesep,'CRCsampleCategory','.mat'));
CRC_sample = COADSampleCategory;
load(strcat('E:\哈工大交流\投稿论文－1\dataset\CHOL',filesep,'CHOL_samplecategory','.mat'));
CHOL_sample = Samplecategory;
% 基因名字,所有数据基因一样
load(strcat('E:\哈工大交流\投稿论文－1\dataset\CHOL',filesep,'CHOL_gelabel','.mat'));
GE_featurename = CHOL_gelabel;
clear Samplecategory CHOL_gelabel;
% 找到癌症样本
index_HNSC_T = find(HNSC_sample==1);
index_PAAD_T = find(PAAD_sample==1);
index_ESCA_T = find(ESCA_sample==1);
index_CRC_T = find(COAD_sample==1);
index_CHOL_T = find(CHOL_sample==1);
% 保留癌症样本
HNSC_GE_T = HNSC_GE(:,index_HNSC_T);
PAAD_GE_T = PAAD_GE(:,index_PAAD_T);
ESCA_GE_T = ESCA_GE(:,index_ESCA_T);
COAD_GE_T = COAD_GE(:,index_COAD_T);
clear ESCA_GE  HNSC_GE PAAD_GE COAD_GE  CHOL_GE HNSC_sample PAAD_sample ESCA_sample CHOL_sample COAD_sample;
% 合并癌症数据（3个数据为例）
H_E_P = [HNSC_GE_T,PAAD_GE_T,ESCA_GE_T];
% 样本重新编号
sample_of_all = [ones(size(index_HNSC_T,1),1);ones(size(index_PAAD_T,1),1)+1;ones(size(index_ESCA_T,1),1)+2];
