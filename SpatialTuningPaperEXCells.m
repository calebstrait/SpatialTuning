function SpatialTuningPaperEXCells
close all

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
end
clearvars -except C

names = {'VM','VS','DA','SA','OF'};
colors = {'b-','g-','r-','m-'};
smo = 20;
vlineCords = [250 287 330] - 188; %offer 1, offer 2, both offers on
for a = 1:5
    for c = 1:length(C{a}.E1)
        f{1} = find(C{a}.FirstOnLeft{c} == 1 & C{a}.EV1{c} >= prctile(C{a}.EV1{c},66));
        f{2} = find(C{a}.FirstOnLeft{c} == 1 & C{a}.EV1{c} < prctile(C{a}.EV1{c},33));
        f{3} = find(C{a}.FirstOnLeft{c} ~= 1 & C{a}.EV1{c} >= prctile(C{a}.EV1{c},66));
        f{4} = find(C{a}.FirstOnLeft{c} ~= 1 & C{a}.EV1{c} < prctile(C{a}.EV1{c},33));
        PSTH = C{a}.PSTH{c};
        if a>2
            PSTH = [PSTH(:,88:end) nan(size(PSTH,1),87)];
        end
        if (a == 5 && c == 31) % max(diffs) > thresh || %% ((a == 1 && c == 59) || (a == 2 && c == 87) || (a == 3 && c == 36) || (a == 4 && c == 102) ||
            figure
            hold on
            h = cell(4,1);
            for x = 1:4
                toplot = mean(PSTH(f{x},100:350))/.02;
                toplot = smooth(toplot,smo);
                tp_SEM = (std(PSTH(f{x},100:350)) / sqrt(length(PSTH(f{x},100:350))))/.02;
                tp_SEM = smooth(tp_SEM,smo);
                h{x} = plot(toplot,colors{x});
                plot(toplot+tp_SEM,colors{x})
                plot(toplot-tp_SEM,colors{x})
            end
            vlineColor(vlineCords(1),[.5 .5 .5]);
            vlineColor(vlineCords(2),[.5 .5 .5]);
            vlineColor(vlineCords(3),[.5 .5 .5]); %1.6s between offer 1 on and both offers on
            ylabel('mean spikes','FontSize',14);
            xlabel('time','FontSize',14);
            title(sprintf('Example cell %s%i',names{a},c),'FontSize',14);
            legend([h{1} h{2} h{3} h{4}],{'Op1 on left, EV1 high','Op1 on left, EV1 low','Op1 on right, EV1 high','Op1 on right, EV1 low'});
            hold off
        end
    end
end

cd /Users/cstrait/Desktop/
keyboard
end