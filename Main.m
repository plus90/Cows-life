clear;
clc;
DeadCows=[];
num=1;
%  herd=initialize(MAXHERD);
load('herd.mat')

%% simulation cinstants
CountOfGroups=3;
LenghntOfherd=100;%  maximum herd size, including youngstock <= 65,535
CountOfSimulationDays=365;
GroupingPeriodLenghnt=30;
%%
dayInformation=[];
for c=1:LenghntOfherd
    herd(c)=herd(c).NewDay();
    dayInformation=[dayInformation;herd(c).record()]
end
dayOfpriode=1;
period=1;
rationOfCurrentperiod=[];
groupofEachCow=[];
currentStateOfherd=dayInformation(:,20:21);
[groupofEachCow,rationOfCurrentperiod]=kmeans(currentStateOfherd,CountOfGroups);
herdGrouped=herd;
herdUngrouped=herd;
history=[];

for day=1:CountOfSimulationDays
    dayInformation=[];
    for c=1:LenghntOfherd
        if herd(c).number == 0
            herd(c).number=num;
            num=num+1;
        end
        herdGrouped(c)=herdGrouped(c).NewDay();
        
        dayInformation=[dayInformation;[day,period,dayOfpriode,herdGrouped(c).record()]];
        if length(herdGrouped(c).event) == 0
            DeadCows{length(DeadCows)+1}=herdGrouped(c);
            herdGrouped(c)=cow('id',22,[],390+randi(60));
        end
    end

    if dayOfpriode==GroupingPeriodLenghnt
        currentStateOfherd=dayInformation(length(dayInformation)-LenghntOfherd:length(dayInformation),22:23);
        [groupofEachCow,rationOfCurrentperiod]=kmeans(CurrentStateOfherd,CountOfGroups);
        period=period+1;
    else
        dayOfpriode=dayOfpriode+1;
    end
    history=[history;dayInformation];
end

