function cattle=initialize(n)
    load('zayesh.mat');
    cows=unique(Zayesh(:,1));
    [NOC , ~]=size(cows);
    for i=1:NOC
        coc(i)=sum(Zayesh{:,1}==cows{i,1});
    end
    for i=1:n
        par = Select([1,2,3,4],[34,25,11,22],rand);
        indofcurent=find(coc>par);
        ind=find(Zayesh{:,1}==cows{indofcurent(randi(length(indofcurent))),1});
        Zayeshs=Zayesh(ind,:);
        f=1;
        for j=2:par
            if isnan(Zayeshs{j,6})==1
                Zayeshs{j,6}=230;
            end
        end
        age=Zayeshs{1,14};
        for j=1:par
            if j>1
                age=age+Zayeshs{j,6};
            end
            tmp(j,1)=Zayeshs(j,3);
            tmp(j,2)=Zayeshs(j,6);
            tmp(j,3)=Zayeshs(j,7);
            tmp(j,4)=Zayeshs(j,8);
            tmp(j,5)=Zayeshs(j,10);
            tmp(j,6)=Zayeshs(j,12);
            tmp(j,7)=Zayeshs(j,13);
            tmp(j,8)=Zayeshs(j,14);
            tmp(j,9)=Zayeshs(j,15);
            tmp(j,10)=Zayeshs(j,16);
            tmp(j,11)=Zayeshs(j,17);
            tmp(j,12)=Zayeshs(j,18);
        end
        cattle(i)=cow(Zayeshs{1,1},22,table2array(tmp),age);
        cattle(i).milch=1;
        cattle(i).heifer=0;
    end
end

