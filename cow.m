classdef cow
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        state
        parturition
        event
        eventinfo
        ID
        age
        doe%DependenceOfEvents
        currentParturition
        baseage
        Delete
        milkrec
        pa
        tee_pee
        mye
        milch
        weights
        BCS
        heifer
        fcm
        DMI
        number
    end
    
    methods
        function b=getBCS(obj)
            b=1;
        end
        
        function lact=Lact(obj)
            [lact,~]=size(obj.parturition);%s=tedad shkam 
            lact=lact+1;
        end
        function obj = cow(id,S,P,age)
            obj.heifer=1;
            obj.ID=id;
            obj.number=0;
            obj.state=S;
            obj.parturition=P;
            obj.pa= normrnd(100,0.9);
            obj.tee_pee=0;
            obj.mye=0;
            obj.eventinfo=[11,12,1;12,13,2;13,14,3;12,21,4;13,21,4;14,21,4;14,14,3;
                            21,22,5;21,14,7;21,11,8;22,11,8;22,11,9;22,11,6;
                            11,31,10;12,31,10;13,31,10;14,31,10;21,31,10;
                            22,31,10;11,31,11;12,31,11;13,31,11;14,31,11;
                            21,31,11;22,31,11];
            milch=0;
            %DependenceOfEvents
            obj.doe=[1,0,0,0,0,0,0,0,0,0,0;%avalin tokhmakgozari
                     1,1,0,0,0,0,0,0,0,0,0;%dovomin tokhmakgozari
                     1,1,0,0,0,0,0,0,0,0,0;%sevomin tokhmakgozari
                     1,1,1,1,0,0,0,0,0,0,0;%talghih va abestani
                     0,0,0,0,1,0,1,1,0,0,0;%rozhay khoshk
                     0,0,0,0,0,1,1,1,1,0,0;%zayeman
                     0,0,0,0,1,1,1,1,1,0,0;%talghih namovafagh
                     0,0,0,0,0,1,0,1,1,0,0;%seght
                     0,0,0,0,0,1,0,0,1,0,0;%mordeh zaei
                     1,1,1,1,1,1,1,1,1,1,1;%hazf
                     1,1,1,1,1,1,1,1,1,1,1]%marg
            obj.event=[age,6];
            obj.age=age-1;
        end
        function obj = destiny(obj)
            switch obj.state
                case 11
                    %hazf dar asar marg
                    [s,~]=size(obj.parturition);%s=tedad shkam 
                    s=s+1;
                    if length(obj.currentParturition)>0 
                        obj.parturition=[obj.parturition;obj.currentParturition];
                        obj.currentParturition=[];
                    end
                    obj.currentParturition(1)=s;
                    obj.currentParturition(12)=0;
                    obj.currentParturition(8)=s;
                    obj.baseage=obj.age;
                    p=[1,0.039;2,0.056;3,0.085];%ehtemal marg dar shakam hay mokhtalef
                    
                    if s>3
                        pod=0.117;
                    else
                        pod=p(s,2);
                    end
                    if(rand<pod)
                       obj.event=[obj.event;[obj.age+randi(500),11]];
                    end
                    %hazf ejbarei
                    p=[1,0.169;2,0.233;3,0.301];%ehtemal hazf dar shakam hay mokhtalef
                    if s>3
                        pod=0.408;
                    else
                        pod=p(s,2);
                    end
                    if(rand<pod)
                       obj.event=[obj.event;[obj.age+randi(500),10]];
                    end
                    
                    %ejad rokhdad avalin tokhmakgozari
                    m = 19; % mean
                    v = 121; % variance
                    mu = log((m^2)/sqrt(v+m^2));
                    sigma = sqrt(log(v/(m^2)+1));
                    r = lognrnd(mu,sigma);
                    obj.event=[obj.event;[obj.age+fix(r),1]];
                case 12
                    if(rand<0.5)
                        %tashkhis fahli va talghih
                        mu = 21;
                        sigma =4;
                        r = normrnd(mu,sigma);
                        if(rand<0.4)
                            obj.event=[obj.event;[obj.age+fix(r),4]];
                            obj.currentParturition(6)=obj.age+fix(r)-obj.baseage;
                        else
                            if(rand<0.2)
                            mu = 10.5;
                            sigma = 2;
                            r = normrnd(mu,sigma);
                            obj.event=[obj.event;[obj.age+fix(r),2]];
                        else
                            mu = 18.5;
                            sigma = 3;
                            r = normrnd(mu,sigma);
                            obj.event=[obj.event;[obj.age+fix(r),2]];
                        end
                        end
                    else
                        %tokhmakgozari hay dovom
                        if(rand<0.2)
                            mu = 10.5;
                            sigma = 2;
                            r = normrnd(mu,sigma);
                            obj.event=[obj.event;[obj.age+fix(r),2]];
                        else
                            mu = 18.5;
                            sigma = 3;
                            r = normrnd(mu,sigma);
                            obj.event=[obj.event;[obj.age+fix(r),2]];
                        end
                    end
                case 13
                    if(rand<0.5)
                        %tashkhis fahli va talghih
                        mu = 21.5;
                        sigma =4;
                        r = normrnd(mu,sigma);
                        if(rand<0.4)
                            obj.event=[obj.event;[obj.age+fix(r),4]];
                            obj.currentParturition(6)=obj.age+fix(r)-obj.baseage;
                        else
                            %eslah shavad
                            mu = 10.5;
                            sigma = 2.5;
                            r = normrnd(mu,sigma);
                            obj.event=[obj.event;[obj.age+fix(r),3]];
                        end
                    else
                        %tokhmak gozari badi
                        mu = 10.5;
                        sigma = 2.5;
                        r = normrnd(mu,sigma);
                        obj.event=[obj.event;[obj.age+fix(r),3]];
                    end
                case 14
                    if(rand<0.5)
                        if(rand<0.4)
                            mu = 21;
                            sigma =4;
                            r = normrnd(mu,sigma);
                            obj.event(find(obj.event(:,2)==3),1)=obj.age+fix(r);
                            obj.event(find(obj.event(:,2)==3),2)=4;
                            obj.currentParturition(6)=obj.age+fix(r)-obj.baseage;
                        else
                            mu = 26;
                            sigma = 8;
                            r = normrnd(mu,sigma);
                            obj.event(find(obj.event(:,2)==3),1)=obj.age+fix(r);
%                             obj.event=[obj.event;[obj.age+fix(r),3]];
                        end
                    else
                        mu = 21;
                        sigma = 2.5;
                        r = normrnd(mu,sigma);
                        obj.event(find(obj.event(:,2)==3),1)=obj.age+fix(r);
%                         obj.event=[obj.event;[obj.age+fix(r),3]];
                    end
                case 21
                    obj.tee_pee=obj.age;
                    if(rand<0.08)
                        %????? ?? ??? ???? ??????
                      obj.event=[obj.event;[obj.age+randi(12)+30,7]];
                    else
                        %????? ????
                        mu = 278;
                        sigma =6;
                        r = normrnd(mu,sigma);
                        obj.event=[obj.event;[obj.age+fix(r),6]];
                        obj.currentParturition(2)=obj.age+fix(r)-obj.baseage;
                        obj.currentParturition(4)=obj.currentParturition(2)-obj.currentParturition(6);
                        if obj.currentParturition(1)>1
                            obj.currentParturition(8)=obj.parturition(end,8);
                        else
                            obj.currentParturition(8)=obj.age+fix(r);
                        end
                        obj.currentParturition(5)=1;
                        obj.currentParturition(7)=1;
                        dod=find(obj.event(:,2)==10);
                        if length(dod)>0
                            obj.event(dod,1)=obj.age+randi(fix(r));
                        end
                        dod=find(obj.event(:,2)==11);
                        if length(dod)>0
                            obj.event(dod,1)=obj.age+randi(fix(r));
                        end
                        tmpi=randi(10);
                        obj.event=[obj.event;[obj.age+fix(r)-tmpi-45,5]];
                        obj.currentParturition(3)=tmpi+45;
                        
                        p=[1,0.0201;2,0.0074;3,0.0055;4,0.0039];
                        [s,~]=size(obj.parturition);%s=tedad shkam 
                        s=s+1;
                        if s>4
                            pod=0.0023;
                        else
                            pod=p(s,2);
                        end
                        if(rand<pod)
                            obj.event=[obj.event;[obj.age+fix(r)-randi(fix(r)-260),9]];
                        end
                    end 
                case 22
                    obj.tee_pee=0;
                case 31
                    obj.currentParturition(9)=randi(2);
                case 32
                    obj.currentParturition(9)=randi(2);
            end
                
        end
        function obj = NewDay(obj)
            obj.age=obj.age+1;
            eventind=find(obj.event(:,1) == obj.age);
            CurrentDayEvents=obj.event(eventind,:);
            [coe , ~]=size(CurrentDayEvents);
            PossibleActions=obj.eventinfo(find(obj.eventinfo(:,1)==obj.state),:);
            for i=1:coe
                currentAction=PossibleActions(find(PossibleActions(:,3)==CurrentDayEvents(i,2)),:);
                 obj.state=currentAction(1,2);
                if currentAction(1,3)==11
                    obj.Delete=[20,obj.age];
                end
                if currentAction(1,3)==6
                    if obj.heifer ==0
                        obj.milch=1;
                    else
                        obj.milch=0;
                    end
                end
                if currentAction(1,3)==5
                    obj.milch=0;
                end
                if currentAction(1,3)==10
                    p=[16.33,28.83,24.39,13.91,6.45,10.09]/100;
                    sum=p(1);
                    r=rand;
                    for j=1:length(p)
                       if r<=sum
                           break;
                       else
                          sum=sum+p(j+1); 
                       end
                    end
                    obj.Delete=[10+j,obj.age];
                end
                obj=obj.destiny();
                does=find(obj.doe(currentAction(1,3),:)==1);
                for i=1:length(does)
                    obj.event(find(obj.event(:,2)==does(i)),:)=[];
                end
            end
        end
        function weight = body_weight(obj)
            weight=obj.weights(end);
        end
        function obj = UpdateMYE(obj)
            obj.mye = 0.9896 * obj.mye + normrnd(0,0.35);
        end
        function weight = body_weight_heifer(obj)
            gain= 0.799;
            age=obj.age;
            if age > 730
                age=730;
            end
            weight=84+(gain*(age-56));  
        end
        function weight = bodyweight(obj)
            age	= obj.age;				%age in days
            A  = 0.0;
            Yo = 0.0;		%bodyweight formula parameters
            K  = 0.0;
            P1 = 0.0;
            P2 = 0.0;
            P3 = 0.0;
            Tpc = 0.0;
            fage = 0.0;
            flactation = 0.0;
            fpregnancy = 0.0;
            weight	= 0.0;			%kg body weight
            dim=obj.age-obj.baseage+1;
            if obj.tee_pee > 0
                if obj.age - obj.tee_pee > 50
                   Tpc = obj.age - obj.tee_pee - 50;
                end
            end
            lac=obj.Lact()-1;
            if lac == 0
                weight = obj.body_weight_heifer();
            elseif lac == 1
                A  = 600;
                Yo = 42;
                K  = 0.0039;
                P1 = 20;
                P2 = 65;
                P3 = 0.0187;
                fage = A * (1 - (1 - (Yo / A^ 0.33333)) * exp(- K * obj.age)^ 3);
                flactation = P1 * dim / P2 * exp(1 - dim / P2);
                fpregnancy = (P3^ 3) * (Tpc^ 3);
                weight = fage - flactation + fpregnancy;
            else
                A  = 600;
                Yo = 42;
                K  = 0.0060;
                P1 = 40;
                P2 = 70;
                P3 = 0.0187;
                fage = A * (1 - (1 - (Yo / A^ 0.33333)) * exp(- K *obj.age)^3);
                flactation = P1 * dim / P2 * exp(1 -  dim / P2);
                fpregnancy = (P3^3) * (Tpc^3);
                weight = fage - flactation + fpregnancy;
            end
        end
        function totalmilk = milkcal(obj)
            lact=obj.Lact();
            switch(lact)
                case 1
                    a = 16.27183810;
                    b = 0.19003941;
                    c = 0.00142474;
                    d = 509.2815505;
                case 2
                    a = 23.80860333;
                    b = 0.16521530;
                    c = 0.00220808;
                    d = 393.1686101;
                otherwise
                    a = 21.82765557;
                    b = 0.22098927;
                    c = 0.00300247;
                    d = 355.7493150;
            end
            dim=obj.age-obj.baseage+1;
            if obj.tee_pee>0
                tepe=obj.age-obj.tee_pee;
            else
                tepe=0;
            end
            y1 = a * (dim^b) * exp((-c)* dim);
            y2 = 1 + ((tepe/d)^2);
            y4 = obj.pa / 100;
            y5 = (y1/y2) * y4;
            totalmilk= y5 + obj.mye;
        end
        function fcm = fcmcal(obj)
%             fcm = (0.454 * obj.milkcal()) + (17.06 * obj.fatCal() * obj.milkcal());
            if obj.heifer==1
                fcm=0;
            else
                fcm = (0.454 * obj.milkcal()) + (17.06 * (obj.fatCal()/100) * obj.milkcal());
            end
        end
        function pctprotein = proteinCal(obj)
            lact=obj.Lact();
            switch(lact)
                case 1
                    a = 5.474283798;
                    b = -0.556080785;
                    c = -0.055906219;
                    d = 0.498971891;
                case 2
                    a = 6.239086473;
                    b = -0.558322461;
                    c = -0.059773830;
                    d = 0.374620207;
                otherwise
                    a = 6.645518738;
                    b = -0.577692291;
                    c = -0.062138807;
                    d = 0.347863878;
            end
            dim=obj.age-obj.baseage+1;
            if dim <= 1		%protects pow_(log(daysim),d) from being zero
                pctprotein = a *(2^b) * exp(-c * (log(2)^2)) * (log(2)^d);
            else
                pctprotein = a * (dim^b) * exp(-c *(log(dim)^2)) *(log(dim)^d);
            end
        end
        function pctfat = fatCal(obj)
            a = 12.86;
            b = -1.081;
            c = -0.0926;
            d = 1.107;
            dim=obj.age-obj.baseage+1;
            if dim <= 1
                pctfat = a * (2^b) * exp(-c * (log(2)^2)) * (log(2)^d);
            else
                pctfat = a *(dim^b) * exp(-c * (log(dim)^2)) * (log(dim)^d);
            end
        end
        function be= bodyEnergy(obj,bcs)
            EBW=0.817* obj.body_weight();
%           be=(0.037683*bcs*EBW*9.4)+(0.200886-0.006676*bcs*EBW*5.55);
            be=((0.037683*bcs*9.4)+((0.200886-0.006676*bcs)*5.55))*EBW;
        end
        function dmi = DMIcow(obj)
            if obj.state == 22
                dmi=11;
            else
                dim=obj.age-obj.baseage+1;
                 if obj.heifer==1
                     dmi = (0.3720 * obj.fcmcal() + 0.0968 * (obj.body_weight() ^ 0.75)) * 1;
                 else
                    dmi = (0.3720 * obj.fcmcal() + 0.0968 * (obj.body_weight() ^ 0.75)) * (1 - exp(-0.192 * (dim/7 + 3.67)));
                 end
                if dmi==0
                    dmi=0;
                end
            end
        end
        function bcs = BCScal(obj)
            if obj.heifer==1
                td = makedist("Triangular","a",3.3,"b",3.5,"c",3.7);
                bcs=random(td);
            else
                dim=obj.age-obj.baseage+1;
                if dim==1
                    td = makedist("Triangular","a",3,"b",3.5,"c",4);
                    bcs=random(td);
                elseif dim>1 && dim<=70
                    td = makedist("Triangular","a",2,"b",2.5,"c",3);
                    bcs=random(td);
                elseif dim>70 && dim<250
                    td = makedist("Triangular","a",2.5,"b",3,"c",3.5);
                    bcs=random(td);
                else
                    td = makedist("Triangular","a",3,"b",3.5,"c",4);
                    bcs=random(td);
                end
            end
        end
        function nemilk = NEMilk(obj,milk,fat,protein) %completed 
            if obj.milch==1
                
                nemilk = milk * ( 0.0929 * (fat/100) +0.0547*(protein/100)+0.192);
%                 nemilk = milk * ( 0.0929 * (fat / milk * 100 ) +0.0547*(protein / milk * 100 )+0.192);
%                 nemilk = milk * ( 0.36 + 0.0969 * (fat / milk * 100 ) );
            else
                nemilk=0;
            end
        end
        function nemaint = NEMaint(obj,bw,BCS) %completed 
%             if obj.heifer==1
%                 SBW = 0.96*bw;
%                 COMP = 0.8+((BCS-1)*0.05);
%                 nemaint = (0.086*(SBW^0.75)*COMP)+18;
%             else
%                 nemaint = 0.079 * (bw ^ 0.75);
%             end
            nemaint = 0.079 * (bw ^ 0.75);
        end
        function neprog = NEProg(obj,pd)%completed 
%             if pd > 150
%                 A = int32(pd-150);
%                 B = int32(30);
%                 c=0.001905/8;
%                 switch idivide(A, B, 'fix')
%                     case 0
%                         neprog=c*8;
%                     case 1
%                         neprog=c*10;
%                     case 2
%                         neprog=c*15;
%                     otherwise
%                         neprog=c*20;
%                 end
%             else
%                 neprog= 0 ;
%             end
            if pd >= 190
                bw=obj.body_weight();
                MBW=0.55*bw;
                CBW=0.06275*MBW;
                neprog=((0.00318*pd-0.352)*(CBW/45))/0.218;
            else
                neprog= 0 ;
            end
        end
        function negrow = NEGrow(obj) %completed 
            if obj.heifer==1
                bw=obj.body_weight();
                MW=1.8182*bw;                
                ADG=((MW*0.92)-(bw))/300;
                EQEBG=(0.956*(ADG^1.097));
                SBW=0.96*bw;
                EQEBW=(0.891*(SBW^0.75));
                negrow=0.0635*(EQEBW^0.75)*(EQEBG^1.097);
            else
                negrow=0;
            end
        end
        function mplact = MPLact(obj,milk,protein)
            if obj.milch==1 && obj.heifer==0
                MTP=(protein/100);
                mplact=(milk*MTP/0.67)*1000;
            else
                mplact=0;
            end
        end
        function mpmaint = MPMaint(obj,dmi,NEL)
            bw=obj.body_weight();
%             MW=1.8182 * bw; 
%             pd=obj.age-obj.tee_pee;
            TDNp=(((NEL/dmi)*0.92)+0.12)/0.0245;
            CP=dmi*(TDNp/100)*130;
            MPBact=CP*0.64;
            mpmaint=(4.1*(bw^0.50))+(0.3*(bw^0.6))+(dmi*30-0.5*((MPBact/0.8)-MPBact))+((0.4*11.8*dmi)/0.67);
%             if obj.heifer==0
%                 CBW=MW*0.06275;
%                 CW=(18+(pd-190)*0.665)*(CBW/45);
%                 mpmaint=(0.3*((bw-CW)^0.6))+(4.1*((bw-CW)^0.5))+(dmi*30);
%             else
%                 mpmaint=(4.1*(bw^0.50))+(0.3*(bw^0.6))+(dmi*30);
%             end
        end
        function mpgrow = MPGrow(obj)
            bw=obj.body_weight();
            MW=1.8182*bw;  
            ADG=((MW*0.92)-(bw))/300;
            MSBW=MW*0.96;
            SRWtoMSBW=478/MSBW;
            CBW=MW*0.06275;
            SBW=bw*0.96;
            CW=(18+((obj.age-obj.baseage)-190)*0.665)*(CBW/45);
            EQSBW=(SBW-CW)*SRWtoMSBW;
            if EQSBW<=478
                EFFMP_NPg=((83.4-(0.114*EQSBW))/100);
            else
                EFFMP_NPg=0.28908;
            end
            NPg=ADG*(268-(29.4*(obj.NEGrow()/ADG)));
            mpgrow=NPg/EFFMP_NPg;
        end
        function mpprog = MPProg(obj,pd)
           if pd >= 190
                bw=obj.body_weight();
                MBW=0.55*bw;
                CBW=0.06275*MBW;
                mpprog=((((0.69*pd)-69.2)*(CBW/45))/0.33);
            else
                mpprog= 0 ;
            end
        end
        function res= RDP(obj,dmi,NEL)
            TDNp=(((NEL/dmi)*0.92)+0.12)/0.0245;
            res=dmi*TDNp*153;
        end
        function res = RUP(obj,MP,dmi,NEL)
            TDNp=(((NEL/dmi)*0.92)+0.12)/0.0245;
            CP=dmi*(TDNp/100)*130;
            MPBact=CP*0.64;
            res=MP-MPBact;
        end
        function res=BE2Weight(BE)
            a=227.0947;
            b1=1.72E-08;
            b2=-9.53E-05;
            b3=0.307526;
            res=a+b1*BE*b2*BE*b3*BE;
        end
        function res=BE2BCS(BE)
            a=-0.68;
            b1=2.75E-10;
            b2=-1.52E-06;
            b3=4.90E-03;
            res=a+b1*BE*b2*BE*b3*BE;
        end
        function res=record(obj)
            if obj.tee_pee>0
                pd=obj.age-obj.tee_pee;% baresi shavad
            else
                pd=0;
            end
            dim=obj.age-obj.baseage+1;
            weight=obj.body_weight();
            BCS=obj.BCScal();
            energy=obj.bodyEnergy(BCS);
            obj.UpdateMYE();
            if obj.milch==1
                milk=obj.milkcal();
                protein=obj.proteinCal();
                fat=obj.fatCal();
            else
               milk =0;
               protein=0;
               fat=0;
            end
            fcm=obj.fcmcal();
            dmi=obj.DMIcow();
            neMilk = obj.NEMilk(milk,fat,protein);
            neMaint = obj.NEMaint(weight,BCS);
            neProg = obj.NEProg(dim);
            neGrow = obj.NEGrow();
            ne=neMilk+neMaint+neProg+neGrow;
            
            mpLact=obj.MPLact(milk,protein);
            mpMaint=obj.MPMaint(dmi,neMilk);
            mpGrow=obj.MPGrow();
            mpprog=obj.MPProg(pd);
            mp=mpLact+mpMaint+mpGrow+mpprog;
            rup=obj.RUP(mp,dmi,neMilk);
            rdp=obj.RDP(dmi,neMilk);
            res=[obj.number obj.age pd dim weight BCS energy milk protein fat fcm dmi neMilk neMaint neProg neGrow  mpLact mpMaint mpGrow ne/dmi mp rdp rup];
        end
    end
end

