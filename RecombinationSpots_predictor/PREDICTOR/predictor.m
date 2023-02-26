clc;
clear all;
close all;
TrainFeatureVector=importdata('TrainFeatureVectorRecombinationSpots.mat');
TrainLabel= importdata('LabelRecombinationSpots.mat');
EmpiricalPSTNP=importdata('EmpiricalPSTNP.mat');
Position=importdata('Position.mat');
all_Triple=importdata('allTriple.mat');

%Enter the file that need to be tested
[FileName, PathName, FilterIndex]=uigetfile('.txt','Select file to open');
[head, seq]=fastaread(FileName);

%The number of sequences
SeqNumber=length(seq); 
%--------------------------------------------
n=length(head); 
t=0;
for i=1:SeqNumber
    Se=seq{1,i}; %The i-th sequence
    L=length(Se);%The length of this sequence
    if L>=131
       nn=L-131+1;
     for j=1:nn
        ss=Se(j:j+130);
        Name{1,t+j}=char(head(i));
        testD{1,t+j}=(ss);
     end   
      t=t+nn;
    end
end
n=length(Name);
%%%%%%%%%%%%%%%%%%%%%%%

AA='ACGT';
V=131;%Sample length
PPT1=zeros(n,129);
PPT2=zeros(n,4^3);
%Feature set of the test sample 
PPT11=zeros(n,172);

%PSTNP
for m=1:n
    Se=testD{1,m};
    for j=1:V-2 
        t1=Se(j);
        k1=strfind(AA,t1);
        t2=Se(j+1);
        k2=strfind(AA,t2);
        t3=Se(j+2);
        k3=strfind(AA,t3);
        PPT1(m,j)=EmpiricalPSTNP(16*(k1-1)+4*(k2-1)+k3,j);
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EIIP
EIIP_ACGT=[0.126 0.134 0.0806 0.1335];
[n3,m3]=size(EIIP_ACGT);
for i=1:m3
    for j=1:m3
        for k=1:m3
        EIIP(n3,k+4*(j-1)+16*(i-1))=EIIP_ACGT(n3,i)+EIIP_ACGT(n3,j)+EIIP_ACGT(n3,k);
        end
    end
end
for m=1:n
    S=testD{1,m};
    for j=1:V-2
    a1=S(j);
    a2=S(j+1);
    a3=S(j+2);
    a=strcat(a1,a2,a3);
    g=strmatch(a,all_Triple,'exact');
    PPT2(m,g)=PPT2(m,g)+1;
    end
end
PPT2=PPT2/129;
PPT2=repmat(EIIP,n,1).*PPT2;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %NC
 K1=[];
 for m=1:n
    A=testD{1,m};
    y=zeros(1,4);
    Amat=A
    for j=1:4
        label=AA(j);
        for h=1:V
            if (Amat(:,h)==label(1))
                y(1,j)=y(1,j)+1;
            end
        end
    end
    K1=[K1;y];
end
 PPT3=K1./V;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %DNC
 K2=[];
for m=1:n
    A=testD{1,m};
    y=zeros(4,4);
    Amat=A;
    for i=1:4
        for j=1:4
            label=[AA(i),AA(j)];
            for h=1:V-1
                if (Amat(:,h)==label(1,1))&&(Amat(:,h+1)==label(1,2))
                    y(i,j)=y(i,j)+1;
                end
            end
        end
    end
    q=y(:)';
    K2=[K2;q];
end
 PPT4=K2./(V-1);
 %%%%%
 PPT=[PPT1,PPT2,PPT3,PPT4];%Combine four kinds feature sets
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=1:172
        k=Position(1,p);
        PPT11(:,p)=PPT(:,k);
 end
%---------------------------
xtest=PPT11;
ytest=zeros(n,1);
model=svmtrain(TrainLabel,TrainFeatureVector,'-c 45.2548 -g 4 -w1 1 -w-1 1');
[predict_label, accuracy, dec_values]=svmpredict(ytest,xtest,model);
clc;
uu=0;
for m=1:n
    
    if predict_label(m)==1
        fprintf(2,'>%s \n',Name{1,m});
        fprintf(2,' %s Hotspot \n',testD{1,m});
    else
        fprintf('>%s \n %s Coldspot \n',Name{1,m},testD{1,m});
    end
end


