function SpatialTuningPaperLATENCY
close all
%Latency encoding SideOfFirst, SideChosen
areas = 1:5;

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
clear stim

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
    C{1}.PSTH{i} = vmPFCdata{i}.psth(valids,130:end);
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
    C{2}.PSTH{i} = VSdata{i}.psth(valids,100:end);
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
    C{3}.PSTH{i} = dACCdata{i}.psth(:,100:end);
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
    C{4}.PSTH{i} = sgACCdata{i}.psth(:,100:end);
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
    C{5}.PSTH{i} = OFCpsths(i).psth(:,100:end);
end
clearvars -except C areas

%%%%%ANALYSIS
% width = 600;
% count = zeros(5,2,width-24);
% for a = areas
%     for x = 1:(width-24)
%         home
%         fprintf('a=%i x=%i\n',a,x);
%         for c = 1:length(C{a}.E1)
%             
%             PSTH = C{a}.PSTH{c};
%             spikes = sum(PSTH(:,x:x+24),2);
%             
%             [~,~,stats] = glmfit(C{a}.FirstOnLeft{c},spikes);
%             if stats.p(2) < .05, count(a,1,x) = count(a,1,x) + 1;end
%             
%             [~,~,stats] = glmfit(C{a}.ChoseLeft{c},spikes);
%             if stats.p(2) < .05, count(a,2,x) = count(a,2,x) + 1;end
%             
%         end
%     end
% end
% save /Users/Jessica/Documents/Writing/GradSchool/Papers/JNeurophysiology/countLat.mat count
load /Users/Jessica/Documents/Writing/GradSchool/Papers/JNeurophysiology/countLat.mat count

%%%%%FIGURES
names = {'VM','VS','DA','SA','OF'};
colors = {'b-','g-','r-','m-','c-'};
legendX = {'FirstOnLeft','SideChosen'};
vlineCords = [150 187];
firstSig = zeros(5,2);
minX = 125;
pthresh = .05;
[maxes,maxesInteger,denoms] = deal([0 0 0 0 0]);
for i = 1:2
    figure
    hold on
    firstSigsALL = {5,1};
    for a = areas
        
        %Trim pre-bin 124
        count(a,i,1:124) = 0;
        
        %Timing analysis
        for x = minX:576
            if (1 - binocdf(count(a,i,x),length(C{a}.E1),0.05)) <= pthresh && firstSig(a,i) == 0
                firstSig(a,i) = x;
            end
        end
        secs = (firstSig(a,i) - 150)* .02;
        fprintf('%s %s: %5.2f sec\n',names{a},legendX{i},secs);
        
        toplot = 100 * count(a,i,:) / length(C{a}.E1);
        toplot2 = permute(toplot,[3 2 1]);
        plot(toplot2,colors{a});
        maxes(a) = max(toplot2);
        maxesInteger(a) = max(count(a,i,:));
        denoms(a) = length(C{a}.E1);
        
        %BOOTSTRAP
        bootn = 5;
        firstSigs = zeros(1,bootn);
        for b = 1:bootn
            toshuffle = permute(count(a,i,:),[3 2 1]);
            shuffled = toshuffle(randperm(length(toshuffle)));  
            for x = minX:576
                if (1 - binocdf(shuffled(x),length(C{a}.E1),0.05)) <= pthresh && firstSigs(b) == 0
                    firstSigs(b) = x;
                end
            end
        end
        firstSigsALL{a} = firstSigs;
        
    end
    ylabel('% of cells','FontSize',14);
    xlabel('time','FontSize',14);
    title(sprintf('%s',legendX{i}),'FontSize',14);
    legend(names);
    hlineColor(100 * binoinv((1-pthresh),length(C{a}.E1),0.05) / length(C{a}.E1),[.5 .5 .5]);
    vlineColor(vlineCords(1),[.5 .5 .5]);
    vlineColor(vlineCords(2),[.5 .5 .5]);
    hold off
    disp(maxes);
    disp(maxesInteger);
    disp(denoms)
    for z = areas
        for y = areas
            if z < y
                fprintf('%s %s Boot:%5.2f Real:%5.2f\n', names{z}, names{y}, prctile(abs(firstSigsALL{z}-firstSigsALL{y}),95),abs(firstSig(z,i)-firstSig(y,i)));
                [~,pVal,X2] = ChiSquareTest([maxes(z) denoms(z); maxes(y) denoms(y)], 0.05);
                fprintf('      Maxes X2:%5.2f p:%5.2f\n', X2, pVal);
            end
        end
    end
    if i == 1, fprintf('\n-------\n\n');end
end

cd /Users/Jessica/Desktop/
keyboard
end