%author: Ting Mei, last edited: 18-09-2023

clear all; clc; close all

% load subcourses
load('...')); % with subID

% load clinical data
data = load('...txt');
%column1 subjectID; column2 group; column3 site; column4 sex; column5 fiq;
%column6 age;7 ADI_social,8 ADI_communication;9 ADI_RRB;10 ADOS_total;11
%ADOS_SA; 12 ADOS_RRB;13 SRS parent; 14 SRS self; 15 SRS combined;16 RBS;
%17 SSP

%remove subjects without fiq
indx = find(data(:,5) == 777 | data(:,5) == 999);
data(indx,:) = [];

%remove subjects without ADI & ADOS
ADI = [data(:,1),data(:,7:9)];
a = find(ADI(:,2) == 777 | ADI(:,2) == 999);
ADI(a,:) = [];
ADOS = [data(:,1),data(:,11:12)];
b = find(ADOS(:,2) == 777 | ADOS(:,2) == 999);
ADOS(b,:) = [];

%  extract subjects with both ADOS and ADI scores
[C,IA,IB] = intersect(ADOS(:,1),ADI(:,1));
ADOS = ADOS(IA,:);
ADI = ADI(IB,:);
data_AD = [ADI, ADOS(:,2:end)];

% extract the subjects with subcourses, ADI and ADOS 
[C,IA,IB] = intersect(data_AD(:,1),sub_courses(:,1));
sub_courses_AD = sub_courses(IB,:);
sub_courses_AD(:,1) = [];
data_AD = data_AD(IA,:);

% same for SRS&RBS&SSP
data_SR = [data(:,1:2), data(:,13), data(:,16:17)]; %SRS parent
for i = 1:3
    c = find(data_SR(:,i+2) == 777 | data_SR(:,i+2) == 999);
    data_SR(c,:) = [];
end

%analysis of SRS&RBS&SSP was done just in autism group
d = find (data_SR(:,2) == 1 |data_SR(:,2) == 3);
data_SR(d,:) = [];
data_SR(:,2) = [];
[C,IA,IB] = intersect(data_SR(:,1),sub_courses(:,1));
sub_courses_SR = sub_courses(IB,:);
sub_courses_SR(:,1) = [];
data_SR = data_SR(IA,:);

%z-zcore clinical data
data_AD = zscore(data_AD(:,2:end));
data_SR = zscore(data_SR(:,2:end));


%CCA
[A_AD,B_AD,r_AD,U_AD,V_AD,stats_AD] = canoncorr(data_AD, sub_courses_AD);
[A_SR,B_SR,r_SR,U_SR,V_SR,stats_SR] = canoncorr(data_SR, sub_courses_SR);

%hauffe for AD
%brain part
selected_filters=1; %the canonical mode used 
hauffe_B_AD=cov(sub_courses_AD) * B_AD(:,selected_filters) * (cov( (B_AD(:,selected_filters)'*sub_courses_AD')')^-1);
[alb_a2_b alb_b2_b]=sort(abs(hauffe_B_AD(:,selected_filters)),'descend');
alb_b2_b(1:10,1)

figure(1)
gra = bar(1:100,sort(hauffe_B_AD(:,selected_filters),'descend'),'FaceColor',[0.65,0.65,0.65],'EdgeColor',[0.65,0.65,0.65]);
gra.BarWidth = 0.7;
view([-270,90]);

%behaviour part
selected_filters=1; 
hauffe_A_AD=cov(data_AD) * A_AD(:,selected_filters) * (cov((A_AD(:,selected_filters)'*data_AD')')^-1);
[alb_a2_p alb_b2_p]=sort(abs(hauffe_A_AD(:,selected_filters)),'descend');
alb_b2_p(:,1)

figure(2);gcf;
color1 = [0.21875,0.34375,0.34375;0.21875,0.34375,0.34375;0.47917,0.56250,0.56250;0.73958,0.78125,0.78125;0.73958,0.78125,0.78125];
[B,I] = sort(abs(hauffe_A_AD),'descend');
xlim([0,6]);
ylim([-0.8,0.8]);
hold on
for i = 1:5;
    color2(I(i,:),:) = color1(i,:);
    gra = bar(I(i,:),hauffe_A_AD(I(i,:),:),'FaceColor',color2(I(i,:),:),'EdgeColor',color2(I(i,:),:));
    hold on
    gra.BarWidth = 0.5;
end
view([-270,90]);

%hauffe for SR
%brain part
selected_filters=1; 
hauffe_B_SR=cov(sub_courses_SR) * B_SR(:,selected_filters) * (cov((B_SR(:,selected_filters)'*sub_courses_SR')')^-1);

[alb_a3_b alb_b3_b]=sort(abs(hauffe_B_SR(:,selected_filters)),'descend');
alb_b3_b(1:10,1)

figure(3)
gra = bar(1:100,sort(hauffe_B_SR(:,selected_filters),'descend'),'FaceColor',[0.65,0.65,0.65],'EdgeColor',[0.65,0.65,0.65]);
gra.BarWidth = 0.7;
view([-270,90]);

%behaviour part
selected_filters=1; %here use the model order... a bot tricky, in oour case just 1, or the numbr of significant ones
hauffe_A_SR=cov(data_SR) * A_SR(:,selected_filters) * (cov((A_SR(:,selected_filters)'*data_SR')')^-1);

[alb_a3_p alb_b3_p]=sort(abs(hauffe_A_SR(:,selected_filters)),'descend');
alb_b3_p(:,1)

figure(4);gcf;
color3 = [0.21875,0.34375,0.34375;0.47917,0.56250,0.56250;0.73958,0.78125,0.78125];
[B,I] = sort(abs(hauffe_A_SR),'descend');
xlim([0,4]);
ylim([-1,1]);
hold on
for i = 1:3;
    color4(I(i,:),:) = color3(i,:);
    gra = bar(I(i,:),hauffe_A_SR(I(i,:),:),'FaceColor',color4(I(i,:),:),'EdgeColor',color4(I(i,:),:));
    hold on
    gra.BarWidth = 0.5;
end
view([-270,90]);


% prepare for figure
filter_order = 1;
%ADI&ADOS
proj_AD_clinical = data_AD * A_AD(:,filter_order); % CCA variate
proj_AD_brain = sub_courses_AD * B_AD(:,filter_order);
[r_AD,p_AD]=corr(proj_AD_clinical,proj_AD_brain);
%SRS&RBS&SSP
proj_SR_clinical = data_SR * A_SR(:,filter_order); % CCA variate
proj_SR_brain = sub_courses_SR * B_SR(:,filter_order);
[r_SR,p_SR]=corr(proj_SR_clinical,proj_SR_brain);

% make plot of main cca mode scatter plots
color_AD = spring(325);
color_SR = spring(194);
%ADI&ADOS
[sorta,Ia] = sort(data_AD(:,end));
for i = 1:325
    a = Ia(i);
    proj_AD_clinical_sort(i,:) = proj_AD_clinical(a,:);
    proj_AD_brain_sort(i,:) = proj_AD_brain(a,:); 
end
%SRS&RBS&SSP
[sortb,Ib] = sort(data_SR(:,end));
for i = 1:194
    a = Ib(i);
    proj_SR_clinical_sort(i,:) = proj_SR_clinical(a,:);
    proj_SR_brain_sort(i,:) = proj_SR_brain(a,:); 
end

%plot of CCA
%ADI&ADOS
figure;
subplot(2,2,1)
xmin = -4;
xmax = 4;
ymin = -4;
ymax = 4;
scatter(proj_AD_brain_sort,proj_AD_clinical_sort,30,color_AD,'o','filled');
title('Main CCA mode scatter plots of subject courses weights versus ADI & ADOS weights','FontSize',10,'FontWeight','normal');
xlabel('CCA weights: subject courses');
ylabel('CCA weights: ADI & ADOS');
str = sprintf('r = %.2f',r_AD);
text(-3,3,str);
xlim([xmin,xmax]);
ylim([ymin,ymax]);
%SRS&RBS&SSP
subplot(2,2,2);
xmin = -4;
xmax = 4;
ymin = -4;
ymax = 4;
scatter(proj_SR_brain_sort,proj_SR_clinical_sort,30,color_SR,'o','filled');
title('Main CCA mode scatter plots of subject courses weights versus SRS & RBS & SSP weights', 'FontSize',10,'FontWeight','normal');
xlabel('CCA weights: subject courses');
ylabel('CCA weights: SRS & RBS & SSP');
str = sprintf('r = %.2f',r_SR);
text(-3,3,str);
xlim([xmin,xmax]);
ylim([ymin,ymax]);


% permutation test with randperm both clinical and subcourses data
N_AD = size(data_AD,1);
N_SR = size(data_SR,1);

filter_order = 1;
for pt = 1:10000;
    rand_indx1 = randperm(N_AD);
    rand_indx2 = randperm(N_AD);
    rand_beh = data_AD(rand_indx1,:);
    rand_brain = sub_courses_AD(rand_indx2,:);
    [X,Y,r_rand,U,V,STATS] = canoncorr(rand_beh,rand_brain);
    r_perm = corr(rand_beh*X(:,filter_order),rand_brain*Y(:,filter_order));
    R_AD(pt) = r_perm;
end
subplot(2,2,3);
histogram(R_AD,'Normalization','Probability','FaceColor',[0.3010 0.7450 0.9330], 'Edgecolor',[0.3010 0.7450 0.9330]);hold on
scatter(r_AD(:,filter_order),0,'bx','LineWidth',2);
title('Randomized subjects, CCA r-values','FontSize',10,'FontWeight','normal');
p_perm_AD = sum(R_AD > r_AD(:,filter_order))/numel(R_AD);
str = sprintf('p = %.4f',p_perm_AD);
text(0.55,0.08,str);
xlim([0.5,1]);
ylim([0,0.1]);


for pt = 1:10000
    rand_indx1 = randperm(N_SR);
    rand_indx2 = randperm(N_SR);
    rand_beh = data_SR(rand_indx1,:);
    rand_brain = sub_courses_SR(rand_indx2,:);
    [X,Y,r_rand,U,V,STATS] = canoncorr(rand_beh,rand_brain);
    r_perm = corr(rand_beh*X(:,filter_order),rand_brain*Y(:,filter_order));
    R_SR(pt) = r_perm;
end
subplot(2,2,4);
histogram(R_SR,'Normalization','Probability','FaceColor',[0.3010 0.7450 0.9330], 'Edgecolor',[0.3010 0.7450 0.9330]);hold on
scatter(r_SR(:,filter_order),0,'bx','LineWidth',2);
title('Randomized subjects, CCA r-values','FontSize',10,'FontWeight','normal');
p_perm_SR = sum(R_SR > r_SR(:,filter_order))/numel(R_SR);
str = sprintf('p = %.4f',p_perm_SR);
text(0.55,0.08,str);
xlim([0.5,1]);
ylim([0,0.1]);

