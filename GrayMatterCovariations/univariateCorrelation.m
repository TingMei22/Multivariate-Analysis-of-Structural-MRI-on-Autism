%author: Ting Mei, last edited: 29-08-2023

clear all;clc
%path
modality_names = {'VBM'};
modality_selected = 1;

all_paths_end = {'_100';'_50';'_autodim'};
selected_path = 1;

paths_ends = {'_ICA_asy', '_ICA_asy_withoutID'};
path_choice = 1;

main_path_ica = strcat('..',modality_names{modality_selected},paths_ends{path_choice},'/results_ICA');
path_ica = strcat(main_path_ica,all_paths_end{selected_path});

path_clinical = '/.../clinical_data';

%load data 
sub_courses = load(strcat(path_ica,'/melodic_mix'));
data = load(strcat(path_clinical,'/....txt'));
%column1 subjectID; column2 group; column3 site; column4 sex; column5 fiq;
%column6 age;7 med_use; 8 handedness; 9 ADI_social,10 ADI_communication;11
%ADI_RRB;12 ADOS_total;13 ADOS_SA; 14 ADOS_RRB;15 SRS parent; 16 SRS self; 17 SRS combined;18 RBS;
%19 SSP;
asym_subs = load(strcat(path_clinical,'/...',paths_ends{path_choice},'.txt'));

for sel_d_end = 1:10;
data_end = {'group','med','handedness','ADI','ADOS','SRS_parent','SRS_self','SRS_combined','RBS','SSP'};
midata{sel_d_end}.names = strcat('midata_',data_end{sel_d_end});
midata{sel_d_end}.data = load(strcat(path_clinical,'/missing_',data_end{sel_d_end},'_subls_all.txt'));% midata_group missing of iq
end

Ncomps = size(sub_courses,2);
Nsubs = size(asym_subs,1);

dummy_site = dummyvar(data(:,3));
data = [data(:,1:2) data(:,4:end) dummy_site];

[c,ia,ib] = intersect(data(:,1),asym_subs);
VBM_data = data(ia,:); 

[c,ia] = setdiff(VBM_data(:,1),midata{1}.data);
VBM_data_group = VBM_data(ia,:);% use for group effect,remove sub miss iq
sub_courses_group = sub_courses(ia,:);
VBM_data_group(:,end) = []; %remove the sixth site;only 1 subject

%relabel medication use
indx = find(VBM_data_group(:,6) == 0);
VBM_data_group(indx,6) = 3; %1 use, 2 unknown,3 no use


% separate autism and control group
indx = find(VBM_data_group(:,2) == 1 | VBM_data_group(:,2) == 3);
VBM_data_group(indx,2) = -1;
VBM_data_td = VBM_data_group(indx,:);
sub_courses_td = sub_courses_group(indx,:);
indx = [find(VBM_data_group(:,2) == 2);find(VBM_data_group(:,2) == 4)];
VBM_data_group(indx,2) = 1;
VBM_data_asd = VBM_data_group(indx,:);
sub_courses_asd = sub_courses_group(indx,:);


% sort data for every sub-statistic
sel_d_end = 1;
DATA{sel_d_end}.name = strcat('VBM_data_',data_end{sel_d_end});   
DATA{sel_d_end}.data = [VBM_data_group(:,1:5),VBM_data_group(:,19:23)];
SUBCOURSES{sel_d_end}.name = strcat('sub_courses_',data_end{sel_d_end});
SUBCOURSES{sel_d_end}.data = sub_courses_group;

sel_d_end = 2; % remove missing med + missing iq
[c,ia] = setdiff(VBM_data_group(:,1),midata{sel_d_end}.data,'stable');
DATA{sel_d_end}.name = strcat('VBM_data_',data_end{sel_d_end});   
DATA{sel_d_end}.data = VBM_data_group(ia,:);
SUBCOURSES{sel_d_end}.name = strcat('sub_courses_',data_end{sel_d_end});
SUBCOURSES{sel_d_end}.data = sub_courses_group(ia,:);

sel_d_end = 3; % missing med + missing iq + missing handedness
[c,ia] = setdiff(DATA{2}.data(:,1),midata{sel_d_end}.data,'stable');
DATA{sel_d_end}.name = strcat('VBM_data_',data_end{sel_d_end});   
DATA{sel_d_end}.data = DATA{2}.data(ia,:);
SUBCOURSES{sel_d_end}.name = strcat('sub_courses',data_end{sel_d_end});
SUBCOURSES{sel_d_end}.data = SUBCOURSES{2}.data(ia,:);

DATA{2}.data = [DATA{2}.data(:,1:6),DATA{2}.data(:,19:end)];
DATA{3}.data = [DATA{3}.data(:,1:5),DATA{3}.data(:,7),DATA{3}.data(:,19:end)];

    for sel_d_end = 4:5;
     [c,ia] = setdiff(VBM_data_asd(:,1),midata{sel_d_end}.data,'stable');
     DATA{sel_d_end}.name = strcat('VBM_data_asd_',data_end{sel_d_end});
     DATA{sel_d_end}.data = VBM_data_asd(ia,:);
     SUBCOURSES{sel_d_end}.name = strcat('sub_courses_asd_',data_end{sel_d_end});
     SUBCOURSES{sel_d_end}.data = sub_courses_asd(ia,:);
    end
    
DATA{4}.data = [DATA{4}.data(:,1:5),DATA{4}.data(:,19:23),DATA{4}.data(:,8:10)];
DATA{5}.data = [DATA{5}.data(:,1:5),DATA{5}.data(:,19:23),DATA{5}.data(:,11:13)];
  
for sel_d_end = 6:10;
     [e,ig] = setdiff(VBM_data_asd(:,1),midata{sel_d_end}.data,'stable');
     DATA{sel_d_end}.name = strcat('VBM_data_asd_',data_end{sel_d_end});
     DATA{sel_d_end}.data = VBM_data_asd(ig,:);
     DATA{sel_d_end}.data = [DATA{sel_d_end}.data(:,1:5),DATA{sel_d_end}.data(:,19:23),DATA{sel_d_end}.data(:,sel_d_end+8)];
     SUBCOURSES{sel_d_end}.name = strcat('sub_courses_asd_',data_end{sel_d_end});
     SUBCOURSES{sel_d_end}.data = sub_courses_asd(ig,:);
end


% standardize the data
for i = 1:10;
DATA_nor{i}.name = strcat('VBM_data_asd_',data_end{i});
DATA_nor{i}.data = [DATA{i}.data(:,1:3),zscore(DATA{i}.data(:,4:5)),DATA{i}.data(:,6:end)];
end
DATA_nor{1}.name = strcat('VBM_data_',data_end{1});
DATA_nor{2}.name = strcat('VBM_data_',data_end{2});
DATA_nor{3}.name = strcat('VBM_data_',data_end{3});

for i = 4:10;
DATA_nor{i}.data = [DATA_nor{i}.data(:,1:10),zscore(DATA_nor{i}.data(:,11:end))];
end

% dummy med & handedness
DATA_nor{2}.data = [DATA_nor{2}.data(:,1:5),dummyvar(DATA_nor{2}.data(:,6)),DATA_nor{2}.data(:,7:end)];
DATA_nor{3}.data = [DATA_nor{3}.data(:,1:5),dummyvar(DATA_nor{3}.data(:,6)),DATA_nor{3}.data(:,7:end)];

%name the outcome variables
 for sel_d_pre = 1:3;
    for sel_d_end = 1:10;
         data_pre = {'b_','dev_','stats_'};
         N{sel_d_pre,sel_d_end}.name = strcat(data_pre{sel_d_pre},data_end{sel_d_end}); 
         N{sel_d_pre,sel_d_end}.vals = [];
    end
 end
 

%GLM for group effect & sensitivity analysis
for component = 1:Ncomps;
    for sel_d_end = 1:3;
        X{sel_d_end}.name = strcat('X_',data_end{sel_d_end});
        X{sel_d_end}.matrix = [DATA_nor{sel_d_end}.data(:,3:end),DATA_nor{sel_d_end}.data(:,2)]; 
        [N{1,sel_d_end}.vals,N{2,sel_d_end}.vals,N{3,sel_d_end}.vals] = glmfit(X{sel_d_end}.matrix,SUBCOURSES{sel_d_end}.data(:,component));
        STATISTIC{sel_d_end}.coeffient(component) = N{1,sel_d_end}.vals(end,1);
        STATISTIC{sel_d_end}.pvals(component) = N{3,sel_d_end}.vals.p(end,1);
    end
end

%FDR
for i=1:10
    SIG{1,i}.name = data_end{i};
end 
addpath(genpath('/home/mrstats/tinmei/exercises/statistic/fdr_functions'));
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(STATISTIC{1,1}.pvals,0.05);
im = find(h == 1);
SIG{1,1}.IC = im';
SIG{1,1}.adj_p = adj_p(:,im)';
SIG{1,1}.orig_p = STATISTIC{1}.pvals(:,im)';
for i = 2:3
    SIG{1,i}.IC = im';
    SIG{1,i}.orig_p = STATISTIC{i}.pvals(:,im)';
end


% GLM for group_by_age/sex interaction, sensitivity analysis
age_squa = DATA_nor{1}.data(:,5).^2;
age_inter = DATA_nor{1}.data(:,2) .* DATA_nor{1}.data(:,5);
age_squa_inter = (DATA_nor{1}.data(:,5).^2).* DATA_nor{1}.data(:,2);

%interaction sex*group
D = [DATA_nor{1}.data(:,2),DATA_nor{1}.data(:,3)];
DD = x2fx(D,'interaction',1:2);
DD = [DD(:,3),DD(:,2),DD(:,4)];

for component = 1:Ncomps;
       X_inter = [DATA_nor{1}.data(:,3:end),DATA_nor{1}.data(:,2),age_squa];%check different confounders by replace with the variables
       [B,DEV,STATS] = glmfit(X_inter,SUBCOURSES{1}.data(:,component));
       age2_beta(component) = B(end,1);
       group_beta(component) = B(end-1,1);
       age2_pvals(component) = STATS.p(end,1);   
       group_pvals(component) = STATS.p(end-1,1); 
end
for component = 1:Ncomps;
       X_inter = [DATA_nor{1}.data(:,4:end),DD];
       [B,DEV,STATS] = glmfit(X_inter,SUBCOURSES{1}.data(:,component));
       sex_inter_beta(component) = B(end,1);
       group_inter_beta(component) = B(end-1,1);
       sex_inter_pvals(component) = STATS.p(end,1);   
       group_inter_pvals(component) = STATS.p(end-1,1); 
end

        
%GLM of behavior scales 
%ADI&ADOS
for j = 1:100
    for sel_d_end = 4:5;
        for i = 1:3;
              X{sel_d_end}.name = strcat('X_',data_end{sel_d_end});
              X{sel_d_end}.matrix = [DATA_nor{sel_d_end}.data(:,3:10),SUBCOURSES{sel_d_end}.data(:,j)];
              [N{1,sel_d_end}.vals,N{2,sel_d_end}.vals,N{3,sel_d_end}.vals] = glmfit(X{sel_d_end}.matrix,DATA_nor{sel_d_end}.data(:,i+10));
              STATISTIC{sel_d_end}.name = data_end{sel_d_end};
              STATISTIC{sel_d_end}.coeffient(i,j) = N{1,sel_d_end}.vals(end,1);
              STATISTIC{sel_d_end}.pvals(i,j) = N{3,sel_d_end}.vals.p(end,1);
        end   
    end
end
%SRS&RBS&SSP
for j = 1:100;
    for sel_d_end = 6:10;
           X{sel_d_end}.name = strcat('X_',data_end{sel_d_end});
           X{sel_d_end}.matrix = [DATA_nor{sel_d_end}.data(:,3:10),SUBCOURSES{sel_d_end}.data(:,j)];
           [N{1,sel_d_end}.vals,N{2,sel_d_end}.vals,N{3,sel_d_end}.vals] = glmfit(X{sel_d_end}.matrix,DATA{sel_d_end}.data(:,11));
           STATISTIC{sel_d_end}.name = data_end{sel_d_end};
           STATISTIC{sel_d_end}.coeffient(:,j) = N{1,sel_d_end}.vals(end,1);
           STATISTIC{sel_d_end}.pvals(:,j) = N{3,sel_d_end}.vals.p(end,1); 
        end
    end


%FDR 
for i = 4:10
all_p_vals = [STATISTIC{1,4}.pvals(:);STATISTIC{1,5}.pvals(:);STATISTIC{1,6}.pvals(:);STATISTIC{1,7}.pvals(:);STATISTIC{1,8}.pvals(:);STATISTIC{1,9}.pvals(:);STATISTIC{1,10}.pvals(:);];
[h_clinical, crit_p_clinical, adj_ci_cvrg_clinical, adj_p_clinical]=fdr_bh(STATISTIC{1,i}.pvals(:),0.05);
im = find(h_clinical == 1);
SIG{1,i}.IC = im';
SIG{1,i}.adj_p = adj_p_clinical(im',:);
SIG{1,i}.orig_p = STATISTIC{1,i}.pvals(:,im)';
end
