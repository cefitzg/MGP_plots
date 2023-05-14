%This script creates and plots the figures for Temporal Dynamics of Faculty Hiring in Mathematics.

clear; clc; close all;

%ANALYSIS

%load in data and compute elite, well-placing and schools for the network analysis
school_list_compute;

%network contruction between 1950-2019 with rolling 10 year window
counter=0; %set counter
for i=1:61
    [hub_ranks,auth_ranks,cit_ind,cmu_ind,com_ind,cor_ind,har_ind,mit_ind,prin_ind,stan_ind,ucb_ind,chi_ind,mich_ind,wash_ind,uwm_ind,yale_ind,~,school_id]=network_gen(fac_data,1950+counter,1959+counter,school_list);
    uni_id{i}=school_id; %record the order of schools relative to the schools_list for each window
    H{i}=hub_ranks; %get the hubs of all schools included in schools_list
    A{i}=auth_ranks; %get the authority of all schools included in schools_list
    %get Cal Tech hubs and authority
    cit_hub(i) = hub_ranks(cit_ind); %cit hub over time
    cit_auth(i) = auth_ranks(cit_ind); %cit auth over time
    %Get CMU hubs and authority by year
    cmu_hub(i) = hub_ranks(cmu_ind); %cmu hub over time
    cmu_auth(i) = auth_ranks(cmu_ind); %cmu auth over time
    %Get Columbia hubs and authority by year
    com_hub(i)=hub_ranks(com_ind); %com hub over time
    com_auth(i)=auth_ranks(com_ind); %com auth over time
    %Get Cornell hubs and authority by year
    cor_hub(i)=hub_ranks(cor_ind); %cor hub over time
    cor_auth(i)=auth_ranks(cor_ind); %cor auth over time
    %Get Harvard hubs and authority by year
    har_hub(i) = hub_ranks(har_ind); %har hub over time
    har_auth(i) = auth_ranks(har_ind); %har auth over time
    %Get MIT hubs and authority by year
    mit_hub(i) = hub_ranks(mit_ind); %mit hub over time
    mit_auth(i) = auth_ranks(mit_ind); %mit authority over time
    %Get Princeton hubs and authority by year
    prin_hub(i) = hub_ranks(prin_ind); %prin hub over time
    prin_auth(i) = auth_ranks(prin_ind); %prin auth over time
    %Get Stanford hubs and authority by year
    stan_hub(i) = hub_ranks(stan_ind); %stan hub over time
    stan_auth(i) = auth_ranks(stan_ind); %stan auth over time
    %Get UCB hubs and authority by year
    ucb_hub(i) = hub_ranks(ucb_ind); %ucb hub over time
    ucb_auth(i) = auth_ranks(ucb_ind); %ucb authority over time.
    %Get UChicago hubs and authority by year
    chi_hub(i) = hub_ranks(chi_ind); %UChicago hub over time
    chi_auth(i) = auth_ranks(chi_ind); %UChicago authority over time.
    %Get Michigan hubs and authority by year
    mich_hub(i) = hub_ranks(mich_ind); %mich hub over time
    mich_auth(i) = auth_ranks(mich_ind); %mich authority over time
    %Get UW hubs and authority by year
    wash_hub(i) = hub_ranks(wash_ind); %wash hub over time
    wash_auth(i) = auth_ranks(wash_ind); %wash authority over time
    %Get UW Madison hubs and authority by year
    uwm_hub(i)=hub_ranks(uwm_ind); %uwm hub over time
    uwm_auth(i)=auth_ranks(uwm_ind); %uwm auth over time
    %Get Yale University hub and authority
    yale_hub(i) = hub_ranks(yale_ind); %yale hub over time
    yale_auth(i) = auth_ranks(yale_ind); %yale authority over time
    start_date(i)=1950+counter; %track date
    counter = counter+1; %update counter
end

close all %close all figures produced by the network constructions

%authority centrality time series for elite departments
auth_ts = vertcat(cit_auth,cmu_auth,com_auth,cor_auth,har_auth,mit_auth,prin_auth,stan_auth, ...
    ucb_auth,chi_auth,mich_auth,wash_auth,uwm_auth,yale_auth)';

%hubs centrality time series for elite departments
hub_ts = vertcat(cit_hub,cmu_hub,com_hub,cor_hub,har_hub,mit_hub,prin_hub,stan_hub, ...
    ucb_hub,chi_hub,mich_hub,wash_hub,uwm_hub,yale_hub)';

%compute the share of hubs and authority held by 14 elite departments by year.
for z =1:length(auth_ts)
    top_14_auth(z) = sum(auth_ts(z,:));
    top_14_hub(z) = sum(hub_ts(z,:));
end

%put hubs and authority scores in order of schools_list
for pp=1:61 %loop through years
    tmp_H = H{pp}; %get hubs for a year for all schools
    tmp_A = A{pp}; %get auth for a year for all schools
    tmp_id_col = uni_id{pp}; %get the order of schools for that year
    for s=1:length(school_list) %loop through all schools in schools_list
        tmp_id = tmp_id_col(s); %get the order relative to schools list for each school
        hubs_in_order(s,pp) = tmp_H(tmp_id); %record the hubs score for each school
        auth_in_order(s,pp) = tmp_A(tmp_id); %record the authority score for each school
    end
end

%preallocate matrices for kendall tau and p values for elite departments
Ken_Mat_Auth = zeros(14,14); %matrix for kendall tau based on kendall tau authority score time series
Ken_Mat_Hub = zeros(14,14); %matrix for kendall tau based on kendall tau hub score time series
Auth_pval = zeros(14,14); %matrix for p values associated with kendall tau based on kendall tau authority score time series
Hub_pval = zeros(14,14); %matrix for p values associated with kendall tau based on kendall tau hub score time series

%compute kendall tau for each combination of elite schools hubs and authority time series
for k = 1:14
    for j = 1:14
        [rho_auth,pval_auth] = corr(auth_ts(:,k),auth_ts(:,j),'Type','Kendall'); %kendall tau for authority scores
        [rho_hub,pval_hub] = corr(hub_ts(:,k),hub_ts(:,j),'Type','Kendall'); %kendall tau for hub scores
        Ken_Mat_Auth(k,j) = rho_auth; %record kendall tau for authority scores pairs
        Auth_pval(k,j) = pval_auth; %record p value associated with kendall tau for authority scores pairs
        Ken_Mat_Hub(k,j)=rho_hub; %record kendall tau for hub score pairs
        Hub_pval(k,j)=pval_hub; %record p value associated with kendall tau for hub score pairs
    end
end

%process Ken_Mat_Auth
Auth_large = Ken_Mat_Auth; %rename Ken_Mat_Auth
Auth_large(Auth_pval >0.05) =0; %if p val is greater than 0.05 set to 0.
Auth_large(abs(Auth_large) <0.5) =0; %if abs(kendall tau score) < 0.5 set to 0
Auth_large(Auth_large == 1)=0; %set diagonal elements (kendall tau of same vector) to zero.
Auth_large = tril(Auth_large); %grab the lower diagonal piece of the authority matrix

%process Ken_Mat_Hub
Hub_large = Ken_Mat_Hub; %rename Ken_Mat_Hub
Hub_large(Hub_pval >0.05) =0; %if p val is greater than 0.05 set to 0.
Hub_large(abs(Hub_large) < 0.5) =0;%if abs(kendall tau score) < 0.5 set to 0
Hub_large(Hub_large == 1)=0; %set diagonal elements (kendall tau of same vector) to zero.
Hub_large = tril(Hub_large); %grab the lower diagonal piece of the hubs matrix.

%compute basic stats about GTF transition rate over time
years = table2array(data(:,2)); %get years
schools = table2array(data(:,1)); %get schools
schools = categorical(schools); %make schools categorical
counter=0; %set counter to 0
for j=1:120 %cycle through 1900-2019
    prof_ind = table2array(data(:,10)); %get flag TRUE or FALSE for professor status
    prof_ind = strcmp(prof_ind,'TRUE'); %convert to 1 and 0.
    prof_ind = prof_ind(years==1900+counter); %filter by year
    fac_prob(j) = sum(prof_ind)/length(prof_ind); %number of profs divided by total phds
    fac_num(j) = sum(prof_ind); %just the numerator
    fac_pos(j) = length(prof_ind); %just the denomator
    start_date_fig1(j) = 1900+counter; %keep track of the data
    counter=counter+1; %increase counter
end

% compute how many schools make up 50% and 90% of the centrality over time.

%preallocate matrices for storage
hub_sum = zeros(2,61); %hub matrix
auth_sum = zeros(2,61); %authority matrix

for ii=1:61 %loop through years
    tmp_col=sort(hubs_in_order(:,ii),'descend'); %order hubs largest to smallest for a given year
    %set counters to 0
    counter=0;
    tmp_v=0;
    c3=0;
    tmp_v3=0;
    for jj=1:length(school_list) %loop through schools
        if tmp_v <=0.5 %go until 0.5
            tmp_v = tmp_v+tmp_col(jj); %add up hub scores until 0.5
            counter = counter+1; %add up number of schools until sum of hub scores = 0.5
        end
        if tmp_v3 <=0.9 %go until 0.9
            tmp_v3 = tmp_v3+tmp_col(jj); %add up hub scores until 0.9
            c3 = c3+1; %add up number of schools until sum of hub scores = 0.9
        end
    end
    hub_sum(1,ii)=counter; %store number of schools
    hub_sum(2,ii)=c3; %store number of schools
end

for iii=1:61 %loop through years
    tmp_col=sort(auth_in_order(:,iii),'descend'); %sort authority from largest to smallest
    %set counters to 0
    counter=0;
    tmp_v=0;
    c3=0;
    tmp_v3=0;
    for jjj=1:length(school_list) %loop through schools
        if tmp_v <=0.5 %until 0.5
            tmp_v = tmp_v+tmp_col(jjj); %add up authority scores until 0.5
            counter = counter+1; %add up number of schools until sum of authority scores = 0.5
        end
        if tmp_v3 <=0.9 %until 0.9
            tmp_v3 = tmp_v3+tmp_col(jjj); %add up authority scores until 0.9
            c3 = c3+1; %add up number of schools until sum of authority scores = 0.9
        end
    end
    auth_sum(1,iii)=counter; %store number of schools
    auth_sum(2,iii)=c3; %store number of schools
end

%GFT by authority score of graduate institution
[auth_ranks,schools]=net_stat(fac_data,1950,2019);
year=data.year; %get year
institution=data.institution; %get phd institution
prof_ind = data.isFacultyTopSchool; %get professor flag
prof_ind = strcmp(prof_ind,'TRUE'); %convert to 1/0

%Probablity of Becoming a Faculty by authority score
faculty_percentage_by_school = zeros(length(schools),1);

for i = 1:length(schools)
    index = strcmpi(schools{i},institution); %get indices for each school
    faculty_percentage_by_school(i) = sum(prof_ind(index))/sum(index); %compute GFT for each school
end

%filter out zeros

nz_ind = find(auth_ranks ~= 0); %find indices of all non-zero authority score
auth_ranks_nz = auth_ranks(nz_ind); %keep non-zero authority scores
faculty_percentage_by_school = faculty_percentage_by_school(nz_ind); %get faculty percentage by school for those with non-zero authority score

%get the R^2 value 
mdl = fitlm(log10(auth_ranks_nz),log10(faculty_percentage_by_school)); 
mdl.Rsquared;

grad_year = fac_data.year; %get year of phd graduation
grad_institution = fac_data.institution; %get grad institution
fac_institution = fac_data.facultyinstitution; %get faculty institution

%build network between 1950 and 2019 to get auth_ranks and schools
[auth_ranks,schools]=net_stat(fac_data,1950,2019);

% PhD Institution Authority Score
grad_institution_auth = zeros(length(grad_institution),1);

for i = 1:length(schools)
    index = strcmpi(schools{i},grad_institution); %find index of each phd school
    grad_institution_auth(index) = auth_ranks(i); %assign phd school correct authority score
end

% Faculty Institution Authority Score
fac_institution_auth = zeros(length(fac_institution),1);

for i = 1:length(schools)
    index = strcmpi(schools{i},fac_institution); %find index of each faculty school
    fac_institution_auth(index) = auth_ranks(i); %assign phd school correct authority score
end

% Graduate Auth vs Faculty Auth over Time (Find people whether their auth moves up)

year_list = unique(grad_year); %get unique phd graduation year
climber_percentage_by_year = zeros(length(year_list),1);

for i = 1:length(year_list)
    index = find(grad_year == year_list(i)); %find index for each grad year
    climber_percentage_by_year(i) = sum(fac_institution_auth(index)> grad_institution_auth(index))/length(index); %compute climber percentage
end

%Averaged graduate production and faculty hires for MIT, Yale, and CMU.

inst = data.institution; %get grad institution 
fac = data.facultyinstitution; %get faculty institution 

isFaculty = data.isFaculty; %use the 0,1 isFaculty flag

%list of schools to plot 
elist_list = {'Carnegie Mellon University', 'Yale University',"Massachusetts Institute of Technology"};

%pre-allocated matrices 
year = zeros(69,3);
num_grads = zeros(69,3);
num_hires = zeros(69,3);

for i = 1:length(elist_list)
    temp_grad_ind = strcmp(inst,elist_list{i}) & strcmp(isFaculty, 'TRUE'); %get indices for people who went to MIT, Yale, or CMU for grad school and became a DG faculty member
    temp_fac_ind = strcmp(fac,elist_list{i}); %get indices for people who became DG faculty that was hired at MIT, Yale, or CMU. 
    
    %filter data 
    temp_grad_data = data(temp_grad_ind,:);
    temp_fac_data = data(temp_fac_ind,:);
    
    %get the years associated with above data 
    temp_years = temp_grad_data.year;
    temp_fac_years = temp_fac_data.year;
    %count up the number of people on a year by year basis
    for j = 1:69
        year(j,i) = 1949+j; %changed to 1949 so it will start in 1950 (j=1)
        num_grads(j,i) = length(find(temp_years == (1949+j))); %changed to 1949 so it will start in 1950 (j=1)
        num_hires(j,i) = length(find(temp_fac_years==(1949+j))); %changed to 1949 so it will start in 1950 (j=1)
    end
    
end

%PLOTTING FIGURES 

%plotting settings
ms = 9; %markersize
fs = 18; %fontsize
fs1 = 6; %fontsize for heat maps
lw = 6; %linewidth
figpos = [100 100 560 420]; %figure position
figpos_sp = [100 100 2*560 420]; %figure position for subplots

%plot raw number of degrees granted and faculty positions recorded
figure(1)
clf
set(gcf,'Position',figpos)
loglog(1900:2019,fac_pos,'LineWidth',lw,'Color','#117A65')
hold on
loglog(1900:2019,fac_num,'LineWidth',lw,'LineStyle',':','Color','#2980B9')
hold on
xlabel('Ph.D. Graduation Year')
ylabel('Count')
set(gca,'fontsize', fs);
axis([1900 2020 0 3200])
legend('Recorded Math Ph.D.s Awarded','Recorded Math Faculty Positions Acquired','Location','southeast')
xticks([1900 1920 1940 1960 1980 2000 2020])

%plot global GTF over time.
figure(2)
clf
set(gcf,'Position',figpos)
plot(1900:2019,fac_prob,'LineWidth',lw,'LineStyle','-','Color','#6C3483 ') %the beginning is dark blue and end is yellow
xlabel('Ph.D. Graduation Year')
ylabel('Graduate to Faculty Transition (GFT) Rate')
axis([1900 2020 0 1])
xticks([1900 1920 1940 1960 1980 2000 2020])
set(gca,'fontsize', fs);

%plot min, med, and max gft transition rate for the well placing schools
figure(3)
clf
set(gcf,'Position',figpos)
plot(start_date_fig3,max_gtf,'LineWidth',lw,'Color','k')
hold on
plot(start_date_fig3,med_gtf,'LineWidth',lw,'Color','r')
plot(start_date_fig3,min_gtf,'LineWidth',lw,'Color','b')
xlabel('Ph.D. Graduation Year')
ylabel('GFT Rate for Well Placing Schools')
set(gca,'fontsize', fs);
axis([start_date_fig3(1) start_date_fig3(end) 0 1])
hold off
legend('Maximum well-placing GFT rate','Median well-placing GFT rate','Minimum well-placing GFT rate')
legend('Location','northeast')
xticks([1950 2015])

%GFT by authority score of grad institution (ALT)
figure(4)
clf
set(gcf,'Position',figpos)
plot(log10(auth_ranks_nz),log10(faculty_percentage_by_school),'o','MarkerSize',ms,'MarkerFaceColor','b')
hold on
plot(log10(linspace(min(auth_ranks_nz),max(auth_ranks_nz),250)),mdl.Coefficients.Estimate(2)*log10(linspace(min(auth_ranks_nz),max(auth_ranks_nz),250))+mdl.Coefficients.Estimate(1),'LineWidth',lw,'LineStyle','-.','Color',[0.5 0 0.5])
xlabel('log_{10}(Authority Score)')
ylabel('log_{10}(GFT Rate)')
set(gca,'FontSize',fs);
axis([min(log10(auth_ranks_nz)) max(log10(auth_ranks_nz)) min(log10(faculty_percentage_by_school)) max(log10(faculty_percentage_by_school))])
% xticks([0.000001 1e-4 1e-2 2e-1])
% xticklabels({'1x10^{-6}','1x10^{-4}','1x10^{-2}','2x10^{-1}'})
% yticks([1e-2 1e-1 5e-1])
% yticklabels({'1x10^{-2}' '1x10^{-1}' '5x10^{-1}'})

%percentage of climbers by year
figure(5)
clf
set(gcf,'Position',figpos)
plot(year_list,climber_percentage_by_year,'o','MarkerSize',ms,'Color','k','MarkerFaceColor','k')
xlim([1950 2019]); ylim([0 0.35]);
xlabel('Ph.D. Graduation Year')
ylabel('Fraction of Graduates Moving Up')
set(gca,'FontSize',fs);

%fraction of authority and hubs held by 50% and 90% of schools over time
figure(6)
clf
set(gcf,'Position',figpos)
plot(1950:2010,hub_sum(1,:),'LineWidth',lw,'LineStyle','-','Color','k')
hold on
plot(1950:2010,hub_sum(2,:),'LineWidth',lw,'LineStyle','-','Color','b')
plot(1950:2010,auth_sum(1,:),'LineWidth',lw,'LineStyle',':','Color','k')
plot(1950:2010,auth_sum(2,:),'LineWidth',lw,'LineStyle',':','Color','b')
xlabel('Year')
ylabel('Number of Departments')
set(gca,'fontsize', fs);
axis([1950 2010 0 length(school_list)])
legend('50% of Hub Centrality','90% of Hub Centrality','50% of Authority Centrality','90% of Authority Centrality','Location','northeast')

%fraction of authority and hubs held by elite schools
figure(7)
clf
set(gcf,'Position',figpos)
plot(1950:2010,top_14_auth,'LineWidth',lw,'LineStyle',':','Color','#FFA500')
hold on
plot(1950:2010,top_14_hub,'LineWidth',lw,'LineStyle','-','Color','#228B22')
xlabel('Year')
ylabel('Fraction of Centrality Held By Elite Depts')
set(gca,'fontsize', fs);
axis([1950 2010 0 1])
legend('Authority Centrality','Hub Centrality')

%plot heat maps of kendall tau based on hub and authority score pairs for elite departments.
figure(8)
clf
set(gcf,'Position',figpos_sp)
%plotting hub
subplot(1,2,1)
myColorMap = jet(256);
myColorMap(65:192,:) = 1; %fix color map so background is white
colormap(myColorMap);
imagesc(Hub_large)
a = colorbar;
ylabel(a,'Kendall \tau','FontSize',fs,'Rotation',270);
ax = gca;
set(gca,'fontsize', fs)
ax.YTick = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
ax.YTickLabel = {'Caltech','CMU','Columbia','Cornell','Harvard','MIT','Princeton','Stanford','Berkeley','Chicago','Michigan','Washington','UW-Madison','Yale'};
ax.XTick = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
ax.XTickLabel = {'Caltech','CMU','Columbia','Cornell','Harvard','MIT','Princeton','Stanford','Berkeley','Chicago','Michigan','Washington','UW-Madison','Yale'};
caxis([-1 1])
%plotting authority
subplot(1,2,2)
myColorMap = jet(256);
myColorMap(65:192,:) = 1; %fix color map so background is white
colormap(myColorMap);
imagesc(Auth_large)
a = colorbar;
ylabel(a,'Kendall \tau','FontSize',fs,'Rotation',270);
ax = gca;
set(gca,'fontsize', fs)
ax.YTick = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
ax.YTickLabel = {'Caltech','CMU','Columbia','Cornell','Harvard','MIT','Princeton','Stanford','Berkeley','Chicago','Michigan','Washington','UW-Madison','Yale'};
ax.XTick = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
ax.XTickLabel = {'Caltech','CMU','Columbia','Cornell','Harvard','MIT','Princeton','Stanford','Berkeley','Chicago','Michigan','Washington','UW-Madison','Yale'};
caxis([-1 1])


%plotting presets for hub and authority time series
slc = cellstr(school_list);
hubs = hubs_in_order;
auths = auth_in_order;
years_f9 = 1950:2010;
color1 = '#D3D3D3';%'#ececec' ; %lighter gray for all schools
color2 = '#999999'; %darker gray for elite schools
highlights = [6,23,106]; % schools to highlight
h2 = [5,10,17,29,32,38,58,50,82,88,98]; %elite schools
%labels = slc(highlights);
elite_lab{1} = 'CMU'; 
elite_lab{2} = 'MIT';
elite_lab{3} = 'Yale';
labels = elite_lab';


%plot hub and authority time series
figure(9);
clf
set(gcf,'Position',figpos_sp)
subplot(1,2,1)
[ii,~,v] = find(hubs');
out = accumarray(ii,v,[],@geomean); % get geometric mean of all schools
po = semilogy(years_f9,hubs',Color=color1,LineWidth=6); %plot all schools
hold all
pe = plot(years_f9,hubs(h2,:),Color=color2,LineWidth=6); %plot just elite schools darker
ax = gca();
ax.ColorOrderIndex = 1;
corder = get(gca,'ColorOrder');
p0 = plot(years_f9,out,'--',LineWidth=6); %plot geometric mean
p1 = plot(years_f9,hubs(highlights,:),LineWidth=6); %plot highlighted schools
legend('-DynamicLegend');
legend([p1;pe(1);po(1);p0],...
    [labels;{'elite schools';'other schools';'geometric mean'}],...
    Location='southwest',FontSize=12);
axis([1950 2010 10^-5 10^0])
xlabel('Date')
ylabel('Hub Scores')
set(gca,'fontsize', fs);

%subplot for authority time series
subplot(1,2,2)
[ii,~,v] = find(auths');
out = accumarray(ii,v,[],@geomean); %get geometric mean of all schools
po = semilogy(years_f9,auths',Color=color1,LineWidth=6); %plot all schools
hold all
pe = plot(years_f9,auths(h2,:),Color=color2,LineWidth=6); %plot just elite schools darker
ax = gca();
ax.ColorOrderIndex = 1;
p0 = plot(years_f9,out,'--',LineWidth=6); %plot geometric mean
p1 = plot(years_f9,auths(highlights,:),LineWidth=6); %plot highlighted schools
legend('-DynamicLegend');
legend([p1;pe(1);po(1);p0],...
    [labels;{'elite schools';'other schools';'geometric mean'}],...
    Location='southwest',FontSize=12);
axis([1950 2010 10^-5 10^0])
xlabel('Date')
ylabel('Authority Scores')
set(gca,'fontsize', fs);

%plot averaged graduate production and averaged faculty hiring for MIT, Yale, and CMU
figure(10)
clf
set(gcf,'Position',figpos_sp)
subplot(1,2,1)
plot(year(:,1),movmean(num_hires(:,1),10),'LineWidth',5,'Color',[0.8500 0.3250 0.0980]);
hold on;
plot(year(:,2),movmean(num_hires(:,2),10),'LineWidth',5,'Color',[0.4940 0.1840 0.5560]);
plot(year(:,3),movmean(num_hires(:,3),10),'LineWidth',5,'Color',[0.9290 0.6940 0.1250]); legend({'CMU','Yale','MIT'},'Location','northwest');
xticks(1950:10:2010); ylim([0 20]);
xlabel('Year'); ylabel('Averaged Faculty Hires');
set(gca,'fontsize', fs);
xlim([1950,2010]);
subplot(1,2,2)
plot(year(:,1),movmean(num_grads(:,1),10),'LineWidth',5,'Color',[0.8500 0.3250 0.0980]);
hold on;
plot(year(:,2),movmean(num_grads(:,2),10),'LineWidth',5,'Color',[0.4940 0.1840 0.5560]);
plot(year(:,3),movmean(num_grads(:,3),10),'LineWidth',5,'Color',[0.9290 0.6940 0.1250]); legend({'CMU','Yale','MIT'},'Location','northwest');
xlabel('Year'); ylabel('Averaged Graduate Production');
xticks(1950:10:2010);
set(gca,'fontsize', fs);
xlim([1950,2010]);


%FUNCTIONS

%network construction function with rolling window of 10 years
function [hub_ranks,auth_ranks,cit_ind,cmu_ind,com_ind,cor_ind,har_ind,mit_ind,prin_ind,stan_ind,ucb_ind,chi_ind,mich_ind,wash_ind,uwm_ind,yale_ind,schools,school_id]=network_gen(data,year1,year2,school_list)
years = table2array(data(:,2)); %get years
doc_inst = table2array(data(:,1)); %get doctoral institution
doc_inst = doc_inst(years <= year2 & years >= year1); %restrict between two years
doc_inst = categorical(doc_inst); %make doc. institute categorical
fac_inst = table2array(data(:,5)); %get faculty institute
fac_inst = fac_inst(years <= year2 & years >= year1); %restrict between two years
fac_inst = categorical(fac_inst); %make faculty institute categorical
d_to_f = horzcat(doc_inst,fac_inst); %column 1 is doc_inst, column 2 is fac_inst
G = digraph(d_to_f(:,2), d_to_f(:,1)); %make network
hub_ranks = centrality(G,'hubs'); %compute hubs score
auth_ranks = centrality(G,'authorities'); %compute authority score
schools=string(G.Nodes{:,1});
%get every school's index relative to school list
for m=1:length(school_list)
    school_id(m) = find(schools == school_list(m));
end

%record elite schools' index
%California Institute of Technology
cit_ind = find(schools == "California Institute of Technology");
%Carnegie Mellon University
cmu_ind = find(schools == "Carnegie Mellon University");
%Columbia University
com_ind = find(schools == "Columbia University");
%Cornell University
cor_ind = find(schools == "Cornell University");
%Harvard University
har_ind = find(schools == "Harvard University");
%Massachusetts Institute of Technology
mit_ind = find(schools == "Massachusetts Institute of Technology");
%Princeton University
prin_ind = find(schools == "Princeton University");
%Stanford University
stan_ind = find(schools == "Stanford University");
%University of California, Berkeley
ucb_ind = find(schools ==  "University of California, Berkeley");
%University of Chicago
chi_ind = find (schools == "The University of Chicago");
%University of Michigan
mich_ind = find(schools == "University of Michigan");
%University of Washington
wash_ind = find(schools == "University of Washington");
%University of Wisconsin-Madison
uwm_ind = find(schools == "University of Wisconsin-Madison");
%Yale University
yale_ind = find(schools == "Yale University");
end

%network construction for two static years
function [auth_ranks,schools]=net_stat(data,year1,year2)
years = table2array(data(:,2)); %get years
doc_inst = table2array(data(:,1)); %get doctoral institute
doc_inst = doc_inst(years <= year2 & years >= year1); %restrict between two years (inclusive)
doc_inst = categorical(doc_inst); %make doc_inst categorical
fac_inst = table2array(data(:,5)); %get faculty institute
fac_inst = fac_inst(years <= year2 & years >= year1); %restrict between two years (inclusive)
fac_inst = categorical(fac_inst); %make fac_inst categorical
d_to_f = horzcat(doc_inst,fac_inst); %column 1 is doc_inst, column 2 is fac_inst
G = digraph(d_to_f(:,2), d_to_f(:,1)); %make network
auth_ranks = centrality(G,'authorities'); %compute authority score
schools=string(G.Nodes{:,1}); %get schools included in network construction
end
