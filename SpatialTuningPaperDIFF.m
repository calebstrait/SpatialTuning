function SpatialTuningPaperDIFF
close all

%%%%%%LOAD DATA
load /Users/Jessica/Documents/Writing/GradSchool/Data/StagOps/vmPFCdata.mat
load /Users/Jessica/Documents/Writing/GradSchool/Data/StagOps/VSgamblingdata.mat
VSdata = VSgamblingdata;
clear VSgamblingdata
load /Users/Jessica/Documents/Writing/GradSchool/Papers/JNeurophysiology/dACC_data.mat
dACCdata = data;
clear data
load /Users/Jessica/Documents/Writing/GradSchool/Papers/JNeurophysiology/sgACC_data.mat
sgACCdata = data;
clear data
load /Users/Jessica/Documents/Writing/GradSchool/Data/Curiosity/ofc_data.mat
OFCdata = data;
clear data
load /Users/Jessica/Documents/Writing/GradSchool/Data/Curiosity/stim.mat
OFCpsths = stim;
clear data

%%%%%FORMAT DATA
C = {5,1};
for i = 1:5
   C{i}.E1 = {};
   C{i}.E2 = {};
   C{i}.E3 = {};
   C{i}.FirstOnLeft = {};
   C{i}.ChoseLeft = {};
   C{i}.ChosenRewardSize = {};
end
for i = 1:length(vmPFCdata)
    valids = (vmPFCdata{i}.vars(:,11) == 1);
    C{1}.E1{i} = mean(vmPFCdata{i}.psth(valids,150:174),2);
    C{1}.E2{i} = mean(vmPFCdata{i}.psth(valids,200:224),2);
    C{1}.E3{i} = mean(vmPFCdata{i}.psth(valids,265:289),2);
    C{1}.FirstOnLeft{i} = vmPFCdata{i}.vars(valids,7);
    C{1}.ChoseLeft{i} = vmPFCdata{i}.vars(valids,8);
    cho1vs2 = vmPFCdata{i}.vars(valids,9);
    ev1 = vmPFCdata{i}.vars(valids,3);
    C{1}.EV1{i} = ev1;
    ev2 = vmPFCdata{i}.vars(valids,6);
    choRewSize = zeros(size(cho1vs2));
    choRewSize(cho1vs2 == 1) = ev1(cho1vs2 == 1);
    choRewSize(cho1vs2 == 2) = ev2(cho1vs2 == 2);
    C{1}.ChosenRewardSize{i} = choRewSize;
    C{1}.PSTH{i} = vmPFCdata{i}.psth(valids,:);
    C{1}.ABSdiff{i} = abs(ev1-ev2);
end
for i = 1:length(VSdata)
    valids = (VSdata{i}.vars(:,11) == 1);
    C{2}.E1{i} = mean(VSdata{i}.psth(valids,150:174),2);
    C{2}.E2{i} = mean(VSdata{i}.psth(valids,200:224),2);
    C{2}.E3{i} = mean(VSdata{i}.psth(valids,265:289),2);
    C{2}.FirstOnLeft{i} = VSdata{i}.vars(valids,7);
    C{2}.ChoseLeft{i} = VSdata{i}.vars(valids,8);
    cho1vs2 = VSdata{i}.vars(valids,9);
    ev1 = VSdata{i}.vars(valids,3);
    C{2}.EV1{i} = ev1;
    ev2 = VSdata{i}.vars(valids,6);
    choRewSize = zeros(size(cho1vs2));
    choRewSize(cho1vs2 == 1) = ev1(cho1vs2 == 1);
    choRewSize(cho1vs2 == 2) = ev2(cho1vs2 == 2);
    C{2}.ChosenRewardSize{i} = choRewSize;
    C{2}.PSTH{i} = VSdata{i}.psth(valids,:);
    C{2}.ABSdiff{i} = abs(ev1-ev2);
end
for i = 1:length(dACCdata)
    C{3}.E1{i} = mean(dACCdata{i}.psth(:,250:274),2);
    C{3}.E2{i} = mean(dACCdata{i}.psth(:,287:311),2);
    C{3}.E3{i} = mean(dACCdata{i}.psth(:,340:364),2);
    o = dACCdata{i}.vars(:,15);
    f = ones(size(o));
    f(o == 2) = 0;
    C{3}.FirstOnLeft{i} = f;
    chos1 = dACCdata{i}.vars(:,16);
    chosL = zeros(size(chos1));
    chosL(chos1 == 1 & f == 1 | chos1 ~= 1 & f ~= 1) = 1;
    C{3}.ChoseLeft{i} = chosL;
    clear o chos1
    cho1vs2 = dACCdata{i}.vars(:,16);
    evL = dACCdata{i}.vars(:,5) .* dACCdata{i}.vars(:,12) + dACCdata{i}.vars(:,6) .* (1-dACCdata{i}.vars(:,12));
    evR = dACCdata{i}.vars(:,7) .* dACCdata{i}.vars(:,13) + dACCdata{i}.vars(:,8) .* (1-dACCdata{i}.vars(:,13));
    [ev1,ev2] = deal(zeros(size(evL)));
    ev1(f == 1) = evL(f == 1);
    ev1(f ~= 1) = evR(f ~= 1);
    ev2(f == 1) = evR(f == 1);
    ev2(f ~= 1) = evL(f ~= 1);
    choRewSize = zeros(size(cho1vs2));
    choRewSize(cho1vs2 == 1) = ev1(cho1vs2 == 1);
    choRewSize(cho1vs2 == 2) = ev2(cho1vs2 == 2);
    C{3}.EV1{i} = ev1;
    C{3}.ChosenRewardSize{i} = choRewSize;
    C{3}.PSTH{i} = dACCdata{i}.psth;
    C{3}.ABSdiff{i} = abs(ev1-ev2);
end
for i = 1:length(sgACCdata)
    C{4}.E1{i} = mean(sgACCdata{i}.psth(:,250:274),2);
    C{4}.E2{i} = mean(sgACCdata{i}.psth(:,287:311),2);
    C{4}.E3{i} = mean(sgACCdata{i}.psth(:,340:364),2);
    o = sgACCdata{i}.vars(:,15);
    f = ones(size(o));
    f(o == 2) = 0;
    C{4}.FirstOnLeft{i} = f;
    chos1 = sgACCdata{i}.vars(:,16);
    chosL = zeros(size(chos1));
    chosL(chos1 == 1 & f == 1 | chos1 ~= 1 & f ~= 1) = 1;
    C{4}.ChoseLeft{i} = chosL;
    clear o chos1
    cho1vs2 = sgACCdata{i}.vars(:,16);
    evL = sgACCdata{i}.vars(:,5) .* sgACCdata{i}.vars(:,12) + sgACCdata{i}.vars(:,6) .* (1-sgACCdata{i}.vars(:,12));
    evR = sgACCdata{i}.vars(:,7) .* sgACCdata{i}.vars(:,13) + sgACCdata{i}.vars(:,8) .* (1-sgACCdata{i}.vars(:,13));
    [ev1,ev2] = deal(zeros(size(evL)));
    ev1(f == 1) = evL(f == 1);
    ev1(f ~= 1) = evR(f ~= 1);
    ev2(f == 1) = evR(f == 1);
    ev2(f ~= 1) = evL(f ~= 1);
    choRewSize = zeros(size(cho1vs2));
    choRewSize(cho1vs2 == 1) = ev1(cho1vs2 == 1);
    choRewSize(cho1vs2 == 2) = ev2(cho1vs2 == 2);
    C{4}.EV1{i} = ev1;
    C{4}.ChosenRewardSize{i} = choRewSize;
    C{4}.PSTH{i} = sgACCdata{i}.psth;
    C{4}.ABSdiff{i} = abs(ev1-ev2);
end
for i = 1:length(OFCdata)
    C{5}.E1{i} = mean(OFCdata(i).epoch1(:,1:24),2);
    C{5}.E2{i} = mean(OFCdata(i).epoch2(:,1:24),2);
    C{5}.E3{i} = mean(OFCdata(i).epoch4(:,1:24),2);
    C{5}.FirstOnLeft{i} = OFCdata(i).var(9,:)';
    f = OFCdata(i).var(9,:)';
    C{5}.ChoseLeft{i} = OFCdata(i).var(1,:)';
    choL = OFCdata(i).var(1,:)';
    choRewSize = zeros(size(choL));
    evL = OFCdata(i).var(5,:)';
    evR = OFCdata(i).var(4,:)';
    choRewSize(choL == 1) = evL(choL == 1);
    choRewSize(choL == 0) = evR(choL == 0);
    [ev1,ev2] = deal(zeros(size(evL)));
    ev1(f == 1) = evL(f == 1);
    ev1(f ~= 1) = evR(f ~= 1);
    ev2(f == 1) = evR(f == 1);
    ev2(f ~= 1) = evL(f ~= 1);
    C{5}.EV1{i} = ev1;
    C{5}.ChosenRewardSize{i} = choRewSize;
    C{5}.PSTH{i} = OFCpsths(i).psth;
    C{5}.ABSdiff{i} = abs(ev1-ev2);
end
clearvars -except C

%%%%%ANALYSIS
names = {'VM','VS','DA','SA','OF'};
legendX = {'FirstOnLeft','SideChosen'}; %{'EV1','FirstOnLeft','EVchosen','SideChosen'};
count = zeros(5,2,2);
trialSet = {[] []};
for a = 1:5
    for c = 1:length(C{a}.E1)
        trialSet{1} = find(C{a}.ABSdiff{c} > nanmedian(C{a}.ABSdiff{c}));  %EASY
        trialSet{2} = find(C{a}.ABSdiff{c} <= nanmedian(C{a}.ABSdiff{c})); %HARD
        for d = 1:2
            
            spikes = C{a}.E1{c}(trialSet{d});
            
%             [~,~,stats] = glmfit(C{a}.EV1{c}(trialSet{d}),spikes);
%             if stats.p(2) < .05, count(a,1,d) = count(a,1,d) + 1;end
            
            [~,~,stats] = glmfit(C{a}.FirstOnLeft{c}(trialSet{d}),spikes);
            if stats.p(2) < .05, count(a,1,d) = count(a,1,d) + 1;end
            
            spikes = C{a}.E3{c}(trialSet{d});
            
%             [~,~,stats] = glmfit(C{a}.ChosenRewardSize{c}(trialSet{d}),spikes);
%             if stats.p(2) < .05, count(a,3,d) = count(a,3,d) + 1;end
            
            [~,~,stats] = glmfit(C{a}.ChoseLeft{c}(trialSet{d}),spikes);
            if stats.p(2) < .05, count(a,2,d) = count(a,2,d) + 1;end
            
        end
    end
end

for a = 1:5
    for v = 1:2
        [~,pci] = binofit(count(a,v,1),length(C{a}.E1),.05/10);
        if count(a,v,2) / length(C{a}.E1) < pci(1) || count(a,v,2) / length(C{a}.E1) > pci(2)
            fprintf('%s %s\n', names{a}, legendX{v});
        end
    end
end
        
cd /Users/Jessica/Desktop/
keyboard
end