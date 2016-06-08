function AllAreasSlideDims
close all
%Sliding percent of cells encoding EV1, SideOfFirst, EVchosen, SideChosen

%%%%%%LOAD DATA
load /Users/cstrait/Documents/Writing/GradSchool/Data/StagOps/vmPFCdata.mat
load /Users/cstrait/Documents/Writing/GradSchool/Data/StagOps/VSgamblingdata.mat
VSdata = VSgamblingdata;
clear VSgamblingdata
load /Users/cstrait/Desktop/OFC2015/dACC_data.mat
dACCdata = data;
clear data
load /Users/cstrait/Desktop/OFC2015/sgACC_data.mat
sgACCdata = data;
clear data
load /Users/cstrait/Documents/Writing/GradSchool/Data/Curiosity/ofc_data.mat
OFCdata = data;
clear data
load /Users/cstrait/Documents/Writing/GradSchool/Data/Curiosity/stim.mat
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
    C{i}.UnChosenRewardSize = {};
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
    C{1}.EV2{i} = ev2;
    [choRewSize,unchoRewSize] = deal(zeros(size(cho1vs2)));
    choRewSize(cho1vs2 == 1) = ev1(cho1vs2 == 1);
    choRewSize(cho1vs2 == 2) = ev2(cho1vs2 == 2);
    C{1}.ChosenRewardSize{i} = choRewSize;
    C{1}.PSTH{i} = vmPFCdata{i}.psth(valids,:);
    unchoRewSize(cho1vs2 == 2) = ev1(cho1vs2 == 2);
    unchoRewSize(cho1vs2 == 1) = ev2(cho1vs2 == 1);
    C{1}.UnChosenRewardSize{i} = unchoRewSize;
    [evL,evR] = deal(zeros(size(ev1)));
    evL(C{1}.FirstOnLeft{i} == 1) = ev1(C{1}.FirstOnLeft{i} == 1);
    evR(C{1}.FirstOnLeft{i} ~= 1) = ev1(C{1}.FirstOnLeft{i} ~= 1);
    evL(C{1}.FirstOnLeft{i} ~= 1) = ev2(C{1}.FirstOnLeft{i} ~= 1);
    evR(C{1}.FirstOnLeft{i} == 1) = ev2(C{1}.FirstOnLeft{i} == 1);
    C{1}.EVL{i} = evL;
    C{1}.EVR{i} = evR;
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
    C{2}.EV2{i} = ev2;
    [choRewSize,unchoRewSize] = deal(zeros(size(cho1vs2)));
    choRewSize(cho1vs2 == 1) = ev1(cho1vs2 == 1);
    choRewSize(cho1vs2 == 2) = ev2(cho1vs2 == 2);
    C{2}.ChosenRewardSize{i} = choRewSize;
    C{2}.PSTH{i} = VSdata{i}.psth(valids,:);
    unchoRewSize(cho1vs2 == 2) = ev1(cho1vs2 == 2);
    unchoRewSize(cho1vs2 == 1) = ev2(cho1vs2 == 1);
    C{2}.UnChosenRewardSize{i} = unchoRewSize;
    [evL,evR] = deal(zeros(size(ev1)));
    evL(C{2}.FirstOnLeft{i} == 1) = ev1(C{2}.FirstOnLeft{i} == 1);
    evR(C{2}.FirstOnLeft{i} ~= 1) = ev1(C{2}.FirstOnLeft{i} ~= 1);
    evL(C{2}.FirstOnLeft{i} ~= 1) = ev2(C{2}.FirstOnLeft{i} ~= 1);
    evR(C{2}.FirstOnLeft{i} == 1) = ev2(C{2}.FirstOnLeft{i} == 1);
    C{2}.EVL{i} = evL;
    C{2}.EVR{i} = evR;
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
    [choRewSize,unchoRewSize] = deal(zeros(size(cho1vs2)));
    choRewSize(cho1vs2 == 1) = ev1(cho1vs2 == 1);
    choRewSize(cho1vs2 == 2) = ev2(cho1vs2 == 2);
    C{3}.EV1{i} = ev1;
    C{3}.EV2{i} = ev2;
    C{3}.ChosenRewardSize{i} = choRewSize;
    C{3}.PSTH{i} = dACCdata{i}.psth;
    unchoRewSize(cho1vs2 == 2) = ev1(cho1vs2 == 2);
    unchoRewSize(cho1vs2 == 1) = ev2(cho1vs2 == 1);
    C{3}.UnChosenRewardSize{i} = unchoRewSize;
    C{3}.EVL{i} = evL;
    C{3}.EVR{i} = evR;
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
    [choRewSize,unchoRewSize] = deal(zeros(size(cho1vs2)));
    choRewSize(cho1vs2 == 1) = ev1(cho1vs2 == 1);
    choRewSize(cho1vs2 == 2) = ev2(cho1vs2 == 2);
    C{4}.EV1{i} = ev1;
    C{4}.EV2{i} = ev2;
    C{4}.ChosenRewardSize{i} = choRewSize;
    C{4}.PSTH{i} = sgACCdata{i}.psth;
    unchoRewSize(cho1vs2 == 2) = ev1(cho1vs2 == 2);
    unchoRewSize(cho1vs2 == 1) = ev2(cho1vs2 == 1);
    C{4}.UnChosenRewardSize{i} = unchoRewSize;
    C{4}.EVL{i} = evL;
    C{4}.EVR{i} = evR;
end
for i = 1:length(OFCdata)
    C{5}.E1{i} = mean(OFCdata(i).epoch1(:,1:24),2);
    C{5}.E2{i} = mean(OFCdata(i).epoch2(:,1:24),2);
    C{5}.E3{i} = mean(OFCdata(i).epoch4(:,1:24),2);
    C{5}.FirstOnLeft{i} = OFCdata(i).var(9,:)';
    f = OFCdata(i).var(9,:)';
    C{5}.ChoseLeft{i} = OFCdata(i).var(1,:)';
    choL = OFCdata(i).var(1,:)';
    [choRewSize,unchoRewSize] = deal(zeros(size(choL)));
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
    C{5}.EV2{i} = ev2;
    C{5}.ChosenRewardSize{i} = choRewSize;
    C{5}.PSTH{i} = OFCpsths(i).psth;
    unchoRewSize(choL == 0) = evL(choL == 0);
    unchoRewSize(choL == 1) = evR(choL == 1);
    C{5}.UnChosenRewardSize{i} = unchoRewSize;
    C{5}.EVL{i} = evL;
    C{5}.EVR{i} = evR;
end

C{1}.leng = length(vmPFCdata);
C{2}.leng = length(VSdata);
C{3}.leng = length(dACCdata);
C{4}.leng = length(sgACCdata);
C{5}.leng = length(OFCdata);

clearvars -except C

%%%%%ANALYSIS
names = {'VM','VS','DA','SA','OF'};
vlineCords = [250 287];
pthresh = .05 / 20;
for a = [1 2 5] %*%*%*%
    figure
    hold on
    toplotMI = zeros(700,1);
    tostarMI = nan(700,1);
    for b = 1:700-24 % # bins in plot
        [R1s,R2s] = deal(zeros(C{a}.leng,1));
        for c = 1:C{a}.leng % # cells
            
            psth = C{a}.PSTH{c};
            spikes = mean(psth(:,b:b+24),2); % Firing rate/trial in 50ms bins, or 25 bins
            
            var1 = C{a}.EVL{c};
            var2 = C{a}.EVR{c};
            [R1,~] = corrcoef(spikes, var1);
            [R2,~] = corrcoef(spikes, var2);
            
            R1s(c) = R1(2,1);
            R2s(c) = R2(2,1);
            
        end
        
        %Are EV1 coefs correlated with EV2 coefs?
        [R,P] = corrcoef(R1s, R2s);
        toplotMI(b) = R(2,1);
        if P(2,1) < .05, tostarMI(b) = .7;end
    end
    
    plot(toplotMI,'b-');
    plot(tostarMI,'c*');
    ylabel('R coef','FontSize',14);
    xlabel('time','FontSize',14);
    title(sprintf('%s EVL & EVR',names{a}),'FontSize',14);
    axis([0 Inf -1 1]);
    hlineColor(100 * binoinv((1-pthresh),length(C{a}.E1),0.05) / length(C{a}.E1),[.5 .5 .5]);
    vlineColor(vlineCords(1),[.5 .5 .5]);
    vlineColor(vlineCords(2),[.5 .5 .5]);
    hold off
end

cd /Users/cstrait/Desktop/
keyboard
end