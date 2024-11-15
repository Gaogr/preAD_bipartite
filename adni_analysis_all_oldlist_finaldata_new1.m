%%%%%%%%%%%%%%%By Lyuan Xu
function [mat_overall,mat_mean,FCD_W,GE_all,mean_FC,Clus_coef_all] = adni_analysis_all_oldlist_finaldata_new1(ses,nettmp,wmn)

%%%% 
datadir = ['C:\Users\*****\afterbatch_finaldata\'];%%%%%your data folder

load('C:\Users\alxne\Documents\MATLAB\new_proc\net17.mat','net17mark');
mat_overall = [];
%%%%%GMmat
for ii = 1 : length(ses)
    load([datadir ses{ii} '\matr2_afterbatch.mat'],'matrES2new'); %%%%%%%each matrES2new is a Nw*Ng matrix
    %%%%%%%%%%%%%%%%%%%%%%%%Nw is the number of white matter bundles, Ng is the
    %%%%%%%%%%%%%%%%%%%%%%%%number of gray matter ROIs 

    mat_overall = cat(3,mat_overall,matrES2new);
    %     mat_overall = cat(3,mat_overall,matrES2);
end

if nettmp >0
    mat_overall = mat_overall(:,net17mark==nettmp,:);
end

%%%%%%%%%%%%%%%%%%%
if wmn ~=0
    mat_overall(wmn,:,:) = NaN;
end
%%%%%%%%%%%%%%%%%%%

mat_mean = mean(mat_overall,3,'omitnan');
%%%%%FCD
FCD_W = squeeze(mean(mat_overall,2,'omitnan'));
%%%%%%%mean FC
mean_FC = squeeze(mean(mat_overall,[1 2],'omitnan'));
%%%%%GE
thre = 0.1;
GE_all = zeros(size(mat_overall,3),1);
Clus_coef_all = zeros(size(mat_overall,3),1);

for ii = 1 : size(mat_overall,3)
    matA_GG = bip_proj(mat_overall(:,:,ii),thre);
    nn = size(mat_overall(:,:,ii),2);
    [D,~]=distance_wei(1./matA_GG);


    %     Dinv1 = 1./D(:);
    %     Dinv = Dinv1;
    %     Dinv(isinf(Dinv1)) = NaN;
    %     GE_all(ii) = sum(Dinv,'omitnan')/(nn*(nn-1));

    Dinv1 = 1./D;
    Dinv = Dinv1;
    Dinv(isinf(Dinv1)) = NaN;
    Dinv_vec = Dinv(triu(ones(size(Dinv,1)))==0);
    GE_all(ii) = mean(Dinv_vec,'omitnan');
    %%%%%%%%%%%%%%%%%%%Clustering Coef
    matA_GG1 = matA_GG./max(matA_GG(:));
    matA_GG2 = matA_GG1;
    idx = isnan(matA_GG1) | isinf(matA_GG1);
    if any(any(idx))
        matA_GG2(idx)=0;
    end

    n = length(matA_GG2);
    matA_GG2(1:n+1:end)=0;

    Clus_coef_all(ii) = mean(clustering_coef_wd(matA_GG2),'omitnan');

    %     GE_all(ii) = sum(Dinv,'omitnan')/(nn*(nn-1));
end

end