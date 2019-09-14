clear;
clc;
tic;


data=load('AAPL.txt');

originaldata=data;
data=transform2(data);
data(:,1:31)=suoyouliezuidazuixiao(data(:,1:31));

len=size(data,1);

traino=data(1:253,:);
test=data(254:377,:);
originaltestopenprice=originaldata(313:436,2);

trainti=traino(:,1:31);

for i=1:31
 %   seed(i).real=ningjucengcijulei(trainti(:,i),0.005);
  %  seed(i).real=ningjucengcijulei(trainti(:,i),0.004);
    seed(i).real=ningjucengcijulei(trainti(:,i),0.003);
 %   seed(i).real=ningjucengcijulei(trainti(:,i),0.001);
%  seed(i).real=ningjucengcijulei(trainti(:,i),0.002);
 %   seed(i).real=ningjucengcijulei(trainti(:,i),0.006);
 %     seed(i).real=ningjucengcijulei(trainti(:,i),0.007);
%     seed(i).real=ningjucengcijulei(train(:,i),0.008);
%     seed(i).real=ningjucengcijulei(trainti(:,i),0.009);
  %   seed(i).real=ningjucengcijulei(trainti(:,i),0.01);
    seed(i).imag=i;
    toc;
end

seedtemp=seed;

for i=1:31
     temp=seedtemp(i).real; 
     for j=1:size(temp,1)
          if sum(temp(j,:)>0)<6
              seedtemp(i).real(j,:)=zeros(1,size(temp,2));
           end
     end
     seedtemp(i).real(all(seedtemp(i).real==0,2),:) = []; 
end

FR=traino(:,32:51);
for i=1:size(FR,1)
    for j=1:size(FR,2)
       if(FR(i,j)>0.005)
          FR(i,j)=3;
       end
       if(FR(i,j)<-0.005)
           FR(i,j)=1;
       end
       if (FR(i,j)<=0.005&&FR(i,j)>=-0.005)
           FR(i,j)=2;
       end
    end
end

deta=0.91;  %deta<1
i=1;
for m=1:size(seedtemp,2) %31
     for n=1:size(seedtemp(m).real,1) %
          temp2=seedtemp(m).real(n,:); %

         [score, expand]=MES(temp2,deta,seedtemp);
         rule(i).real=expand;
         rule(i).imag=score;

          toc;
          i=i+1
     end
end 
 
n=1;
 for j=1:size(rule,2)
    if ~isempty(rule(j).real)
        rule2(n)=rule(j);
        n=n+1;
    end
 end

 allrule=[];

for i=1:size(rule2,2)
    temp=rule2(i).real;
    
    for j=1:size(temp,2)
       
        if sum(temp(:,j))==0
          finalrule(i,j)=0;
        else
         finalrule(i,j)=mean(trainti(temp(:,j),j));
     
          cindex=temp(:,j);        
        end
end 
    trend=qiurule(FR,cindex);
   for m=1:length(trend)
     if trend(m)>0
         tr=[finalrule(i,:) trend(m)];
         allrule=[allrule; tr];
     end
   end
end



for j=1:size(trainti,1)
    for i=1:size(allrule,1)
     d(j,i)=qiujuli(allrule(i,1:31),trainti(j,:));
    end
end

for j=1:size(test,1)
    for i=1:size(allrule,1)
     dtest(j,i)=qiujuli(allrule(i,1:31),test(j,1:31));
    end
end

class=FR(:,20);
P=d';
T=class';
[p1,minp,maxp,t1,mint,maxt]=premnmx(P,T);

% net=newff(minmax(P),[8,1],{'tansig','purelin'},'trainlm');
net=newff(minmax(P),[round(sqrt(size(allrule,1)+1)+5),1],{'tansig','purelin'},'trainlm');
net.trainParam.epochs = 10000;
net.trainParam.goal=0.01;
[net,tr]=train(net,p1,t1);

%a=dtest(1,:)';
a=dtest';
a=premnmx(a);
b=sim(net,a);
c=postmnmx(b,mint,maxt);


c=round(c);
haspos=-1;
pr=0;
so=[];
bo=[];
for i=1:length(c)-1
    if haspos==1
        	 if(c(i)==1)
                 so=[so i+1];
                 sp=originaltestopenprice(i+1);
                 pr=pr+(sp-bp)/bp;
                 haspos=-1;
             end
    else
             if(c(i)==3)    
                 bo=[bo i+1];
                 bp=originaltestopenprice(i+1);
                 haspos=1;
             end
    end
end










