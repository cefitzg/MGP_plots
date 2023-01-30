clear; clc; close all 

%Figure 10 (we cannot release targetlist6.csv publicly)

%Written by Yitong Huang and Cody FitzGerald 

%read in data 
data = readtable('MGPdata.csv'); %publicly released data (anonymized data)
data_cmu_fac = readtable('targetlist6.csv'); %nonanonymized data 

inst = data_cmu_fac.institution; %get grad training department 
fac = data_cmu_fac.facultyinstitution; %get faculty department 

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
    %count up the numbe rof people on a year by year basis
    for j = 1:69
        year(j,i) = 1950+j;
        num_grads(j,i) = length(find(temp_years == (1950+j)));
        num_hires(j,i) = length(find(temp_fac_years==(1950+j)));
    end
    
end

%plotting 

%plotting settings
ms = 9; %markersize
fs = 18; %fontsize
fs1 = 6; %fontsize for heat maps
lw = 6; %linewidth
figpos = [100 100 560 420]; %figure position

figure(101)
clf
set(gcf,'Position',figpos)
plot(year(:,1),movmean(num_grads(:,1),10),'LineWidth',5,'Color',[0.8500 0.3250 0.0980]);
hold on;
plot(year(:,2),movmean(num_grads(:,2),10),'LineWidth',5,'Color',[0.4940 0.1840 0.5560]);
plot(year(:,3),movmean(num_grads(:,3),10),'LineWidth',5,'Color',[0.9290 0.6940 0.1250]); legend({'CMU','Yale','MIT'},'Location','northwest');
xlabel('Year'); ylabel('Averaged Graduate Production');
xticks(1950:10:2010);
set(gca,'fontsize', fs);
xlim([1950,2010]);

figure(102)
clf
set(gcf,'Position',figpos)
plot(year(:,1),movmean(num_hires(:,1),10),'LineWidth',5,'Color',[0.8500 0.3250 0.0980]);
hold on;
plot(year(:,2),movmean(num_hires(:,2),10),'LineWidth',5,'Color',[0.4940 0.1840 0.5560]);
plot(year(:,3),movmean(num_hires(:,3),10),'LineWidth',5,'Color',[0.9290 0.6940 0.1250]); legend({'CMU','Yale','MIT'},'Location','northwest');
xticks(1950:10:2010); ylim([0 20]);
xlabel('Year'); ylabel('Averaged Faculty Hires');
set(gca,'fontsize', fs);
xlim([1950,2010]);