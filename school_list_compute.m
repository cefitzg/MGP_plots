clear; clc; close all;

%This code determines the elite schools, the well-placing schools, and the
%106 schools to use in the network analysis for Temporal Dynamics of
%Faculty Hiring in Mathematics 

%load in data
data = readtable('MGPdata.csv');

%filter the data to records that were trained and became a professor at a top 150 US math department
top_fac_ind = table2array(data(:,12)); %true or false indicating if the entry was a professor at a top 150 US math department.
top_fac_flag = strcmp(top_fac_ind,'TRUE'); %convert true false to 1/0
fac_ind = find(top_fac_flag==1); %find indices such that top_fac_flag==1.
fac_data = data(fac_ind,:); %restrict data set to those indices.

%determine 106 schools to use in the network analysis
counter=0; %set counter
for i=1:61
    [~,schools]=network_gen(fac_data,1950+counter,1959+counter); %construct network for 10 year window starting in 1950 and ending in 2010.
    if i==1
        school_list = schools; %keep schools included in network construction.
    else
        school_list = intersect(school_list,schools); %iteratively intersect school_list with schools included in new network construction.
    end
    counter = counter+1; %update counter
end

%determine elite schools 
dec_move = [0 10 20 30 40 50 60]; %vector to change start and end dates

for k=1:length(dec_move)
    ad_year = dec_move(k); %get number of years to shift 10 year window
    start_year=1950+ad_year; %start year for network construction
    end_year=1959+ad_year; %end year for network construction
    [auth_ranks,schools]=network_gen(fac_data,start_year,end_year); %generate network
    auth_ordered = sort(auth_ranks,'descend'); %order authority largest to smallest
    for j=1:7
        top_7_ind(j) = find(auth_ranks == auth_ordered(j)); %get indicies of the seven largest authority scores
    end
    if k==1
        elite_schools = schools(top_7_ind); %get elite schools
        elites = elite_schools; %rename elite schools
    else
        new_elites = schools(top_7_ind); %get new elite schools
        elites = union(elites,new_elites); %iteratively take the union of elites and new_elites
    end
end
disp("the elite schools include:")
disp(categorical(elites))


%compute the well placing schools
sl = categorical(school_list); %categorical version of school_list
years = table2array(data(:,2)); %get years
schools = table2array(data(:,1)); %get schools from data
schools = categorical(schools); %convert categorical
counter=0; %set counter
for b=1:length(school_list) %loop through schools
    counter=0; %reset counter
    for a=1:66 %loop through years
        start_date_fig3(a) = 1950+counter; %record x axis
        prof_ind = table2array(data(:,12)); %get TRUE/FALSE for professor status
        prof_ind = strcmp(prof_ind,'TRUE'); %convert TRUE/FALSE to 1/0
        prof_ind = prof_ind(years==1950+counter & schools == sl(b)); %get prof_ind flag for every school in school_list every year
        if isempty(prof_ind)==1
            fac_prob_schools(a,b) = -2; %if prof_ind is empty (no students graduated in that year in the MGP database) tag with -2.
            counter=counter+1; %increase counter
        else
            fac_prob_schools(a,b) = sum(prof_ind)/length(prof_ind); %if prof_ind is not empty compute # of eventual faculty members / total number of graduates
            counter=counter+1; %increase counter
        end
    end
end

%compute gft for well placing schools
for zz =1:length(school_list) %loop through the schools for all years
    x_temp = find(fac_prob_schools(:,zz) ==-2 | fac_prob_schools(:,zz) ==0); %find all -2 (did not graduate anyone) or fac_prob_schools = 0 (graduated but no one become a professor)
    if isempty(x_temp) == 1
        keep_ind(zz) = zz; %indices of schools that always graduate someone who will become a prof.
    end
end

keep_ind = nonzeros(keep_ind); %get rid of zeros.
well_placing_schools = sl(keep_ind); %schools that have graduated someone who will become a professor every year between 1950 and 2015.
fac_trans_rate = fac_prob_schools(:,keep_ind); %gft transition rate for well placing schools

for ab=1:66
    min_gtf(ab) = min(fac_trans_rate(ab,:)); %compute minimum gft transition rate for the well placing schools
    med_gtf(ab) = median(fac_trans_rate(ab,:)); %compute the median gft transition rate for the well placing schools
    max_gtf(ab) = max(fac_trans_rate(ab,:)); %compute the maximum gft transition rate for well placing schools
end
disp("The well placing schools include:")
disp(well_placing_schools)

%network generation function.
function [auth_ranks,schools]=network_gen(data,year1,year2)
years = table2array(data(:,2)); %get years
doc_inst = table2array(data(:,1)); %get doctoral institution
doc_inst = doc_inst(years <= year2 & years >= year1); %restrict between two years
doc_inst = categorical(doc_inst); %make doc. institute categorical
fac_inst = table2array(data(:,5)); %get faculty institute
fac_inst = fac_inst(years <= year2 & years >= year1); %restrict between two years
fac_inst = categorical(fac_inst); %make faculty institute categorical
d_to_f = horzcat(doc_inst,fac_inst); %column 1 is doc_inst, column 2 is fac_inst
G = digraph(d_to_f(:,2), d_to_f(:,1)); %make network
schools=string(G.Nodes{:,1}); %get schools represented in the network between year1 and year2
auth_ranks = centrality(G,'authorities'); %compute authority score
end