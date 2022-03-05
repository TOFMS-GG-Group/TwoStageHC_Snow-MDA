%HEADER%

%***********************************%
%Section Below does the 1D Clustering.
%User parameters for cutoff of clusters (1D and 2D) and occurance frequency
%for inclusion of elements in representative mmNP signature of each cluster
CutoffVal=0.6;
CutoffVal_2=0.5;
OccurenceFreq=0.1;

%FileInformation
ClusterMetaDataFile = 'ClusterInfoFile.xlsx';
Filenames=readmatrix(ClusterMetaDataFile,'Sheet','Filename','Range',2,'OutputType','string');
Filenames=Filenames(:,1);
ElementNames=readmatrix(ClusterMetaDataFile,'Sheet','Elements','OutputType','string');
q_Plasma=readmatrix(ClusterMetaDataFile,'Sheet','q_plasma','OutputType','double');
PNCCorrFactors_mL=q_Plasma(:,1).*q_Plasma(:,2);
PNCCorrFactors_mL=PNCCorrFactors_mL/1000;
NbrElements=size(ElementNames,1);

%Section Below does the 1D Clustering.
ElementNames=reshape(ElementNames,1,NbrElements);

ClstrProxyMassesCompile=[];
ClstrProxyNamesCompile=[];
ClstrNbrParticlesCompile=[];
ClstrPNCsCompile=[];
ClstrParticles3DCompile={};
for b=1:size(Filenames,1)

    NPData=xlsread(append(Filenames(b),'.csv'));
        
    Z = linkage(NPData,'average','correlation');
    ClstrAssign=cluster(Z,'Cutoff',CutoffVal,'criterion','distance');
    NbrClusters=max(ClstrAssign);
    figure();
    h=dendrogram(Z,'ColorThreshold',CutoffVal);
    TitleDendrogram = append(Filenames(b),'-Dendrogram');
    title(TitleDendrogram);
    saveas(gcf,TitleDendrogram);
    close(gcf);

    ClstrProxyMasses = [];
    ClstrProxyNames =[];
    ClstrDataNames =[];
    ClstrNbrParticles = [];
    ClstrParticles3D={};

    for i=1:NbrClusters
        index=find(ClstrAssign==i);
        df=NPData(index,:);
        SizeClstr=size(df);
        if SizeClstr(1)>1
            NbrEventsInCluster=sum(df>0);
            SumMassesInCluster=sum(df);
            df(df==0)=NaN;
            OneClassProxyMeanMass=nanmean(df);
            df(isnan(df))=0;
            ClusterOccurance=(NbrEventsInCluster/SizeClstr(1));
            MaskOccurance=ClusterOccurance;
            MaskOccurance(ClusterOccurance<OccurenceFreq)=0;
            MaskOccurance(MaskOccurance>OccurenceFreq)=1;

                for j=1:size(MaskOccurance)
                    k=find(MaskOccurance==1);
                end

                for n=1:size(k,2)
                    ClusterElements=ElementNames(1,k);
                    ElementMasses_Cluster=df(:,k);
                    Nbr_Cluster=NbrEventsInCluster(1,k);
                end

            [temp, order] = sort(Nbr_Cluster(1,:),'descend');
            ClusterElementsSort = ClusterElements(:,order);

            if size(ClusterElements,2)>2
                OneClassProxyName = append(Filenames(b),'-',ClusterElementsSort(1),ClusterElementsSort(2),ClusterElementsSort(3));
                OneClassDataName = append(Filenames(b),'_',ClusterElementsSort(1),ClusterElementsSort(2),ClusterElementsSort(3),'_Data');
            elseif size(ClusterElements,2)>1
                OneClassProxyName = append(Filenames(b),'-',ClusterElementsSort(1),ClusterElementsSort(2));
                OneClassDataName = append(Filenames(b),'_',ClusterElementsSort(1),ClusterElementsSort(2),'_Data');
            else
                OneClassProxyName = append(Filenames(b),'-',ClusterElementsSort(1));
                OneClassDataName = append(Filenames(b),'_',ClusterElementsSort(1),'_Data');
            end

           
           OneClassProxyMeanMass=OneClassProxyMeanMass.*MaskOccurance;
           OneClassProxyMeanMass(isnan(OneClassProxyMeanMass))=0;
           

           OneClassNbrParticles = SizeClstr(1);
           OneClassNbrParticles(isnan(OneClassNbrParticles))=0;
           
            currentFile = sprintf(OneClassDataName,'.mat');
            save(currentFile,'df');

        ClstrProxyMasses = cat(1,ClstrProxyMasses, OneClassProxyMeanMass);
        ClstrProxyNames = cat(1,ClstrProxyNames,OneClassProxyName);
        ClstrDataNames = cat(1,ClstrDataNames,OneClassDataName);
        ClstrNbrParticles = cat(1,ClstrNbrParticles,OneClassNbrParticles);
        ClstrParticles3D = cat(3,ClstrParticles3D,{df});

        end
    end
    
ClstrPNCs = ClstrNbrParticles/PNCCorrFactors_mL(b);

ClstrProxyMassesCompile=cat(1,ClstrProxyMassesCompile,ClstrProxyMasses);
ClstrProxyNamesCompile=cat(1,ClstrProxyNamesCompile,ClstrProxyNames);  
ClstrNbrParticlesCompile=cat(1,ClstrNbrParticlesCompile,ClstrNbrParticles);
ClstrPNCsCompile=cat(1,ClstrPNCsCompile,ClstrPNCs);
ClstrParticles3DCompile=cat(3,ClstrParticles3DCompile,ClstrParticles3D);
end

clear ClstrProxyMasses ClstrProxyNames ClstrDataNames ClstrNbrParticles ClstrParticles3D 
clear OneClassProxyMass ClstrPNCs currentFile OneClassProxyName OneClassDataName NbrClusters 
clear ClusterAssign TitleDendrogram NPData NbrEventsInCluster MaskOccurance SumMassesInCluster
clear ElementMasses_Cluster OccurenceFreq OneClassNbrParticles NbrCluster
clear ClusterElements ClusterElementsSort ClustOccurance SizeClstr df  

%Section Below does the 2D Clustering.
Z = linkage(ClstrProxyMassesCompile,'average','correlation');
D = pdist(ClstrProxyMassesCompile);
LeafOrder = optimalleaforder(Z,D);
ClstrAssign=cluster(Z,'Cutoff',CutoffVal_2,'criterion','distance');

NbrClusters=max(ClstrAssign);
NbrProxyParticles=size(ClstrProxyMassesCompile);
NbrProxyParticles=NbrProxyParticles(1);

ClstrAssignSorted=ClstrAssign(LeafOrder,:);
ClstrProxyMassesSorted = ClstrProxyMassesCompile(LeafOrder,:);
ClstrNbrParticlesSorted=ClstrNbrParticlesCompile(LeafOrder,:);
ClstrProxyNamesSorted=ClstrProxyNamesCompile(LeafOrder,:);
ClstrPNCsSorted=ClstrPNCsCompile(LeafOrder,:);
ClstrParticles3DSorted = ClstrParticles3DCompile(:,:,LeafOrder);

ParticlesInClusters={};
for k = 1:NbrClusters
    index=find(ClstrAssignSorted==k);
    OneClassParticles=[];
    for i=1:size(index)
        ParticleElementMasses=cell2mat(ClstrParticles3DSorted(:,:,index(i)));
        OneClassParticles=vertcat(OneClassParticles,ParticleElementMasses);
    end
   ParticlesInClusters=cat(3,ParticlesInClusters,OneClassParticles);
end
PlotClstrOrder=unique(ClstrAssignSorted,'stable');
ParticlesInClusters=ParticlesInClusters(:,:,PlotClstrOrder);

for n=1:NbrClusters
    Data=cell2mat(ParticlesInClusters(:,:,n));
    DataName=strcat('2DClstr_',num2str(n));
    currentFile = sprintf(DataName,'.mat');
    save(currentFile,'Data');
end

clear ParticleElementMasses OneClassParticles OneClassProxyMeanMass


%Section Below is plotting
mycolormap=[1 1 1
    0.817661881446838 0.803669035434723 0.959092855453491
    0.635323822498322 0.607338070869446 0.918185710906982
    0.45298570394516 0.411007136106491 0.877278566360474
    0.270647615194321 0.214676186442375 0.836371421813965
    0.275114297866821 0.234238088130951 0.870985686779022
    0.278299987316132 0.255871415138245 0.899071455001831
    0.27572038769722 0.278222441673279 0.91279798746109
    0.27314081788063 0.300573468208313 0.926524519920349
    0.270561218261719 0.322924494743347 0.940251052379608
    0.267981618642807 0.345275491476059 0.953977525234222
    0.265402019023895 0.367626518011093 0.967704057693481
    0.262822449207306 0.389977544546127 0.98143059015274
    0.260242849588394 0.412328571081161 0.995157122612
    0.244033336639404 0.435833334922791 0.998833358287811
    0.220642849802971 0.460257142782211 0.997285723686218
    0.196333333849907 0.484719038009644 0.989152371883392
    0.183404758572578 0.507371425628662 0.979795217514038
    0.179921418428421 0.528638124465942 0.965907096862793
    0.176438093185425 0.549904763698578 0.952019035816193
    0.168742850422859 0.570261895656586 0.935871422290802
    0.153999999165535 0.590200006961823 0.921800017356873
    0.146011903882027 0.608914256095886 0.909545242786407
    0.138023808598518 0.627628564834595 0.897290468215942
    0.123752377927303 0.645028591156006 0.884787321090698
    0.109480954706669 0.662428557872772 0.872284114360809
    0.0952095240354538 0.679828584194183 0.859780967235565
    0.0688714310526848 0.694771409034729 0.839357137680054
    0.0296666659414768 0.708166658878326 0.81633335351944
    0.00357142859138548 0.72026664018631 0.791700005531311
    0.00665714265778661 0.731214284896851 0.766014277935028
    0.0433285720646381 0.741095244884491 0.739409506320953
    0.0963952392339706 0.75 0.712038099765778
    0.140771433711052 0.758400022983551 0.684157133102417
    0.171700000762939 0.766961932182312 0.655442833900452
    0.193892866373062 0.775630950927734 0.623871445655823
    0.216085717082024 0.784300029277802 0.592299997806549
    0.246957138180733 0.791795253753662 0.55674284696579
    0.290614277124405 0.797290503978729 0.518828570842743
    0.345095902681351 0.794716358184814 0.476524502038956
    0.399577528238297 0.7921422123909 0.434220403432846
    0.454059153795242 0.789568066596985 0.391916334629059
    0.50854080915451 0.786993861198425 0.349612236022949
    0.563022434711456 0.78441971540451 0.307308167219162
    0.617504060268402 0.781845569610596 0.265004068613052
    0.671985685825348 0.779271423816681 0.222699999809265
    0.724200010299683 0.769842863082886 0.191028565168381
    0.773833334445953 0.759804785251617 0.164609521627426
    0.820314288139343 0.74981427192688 0.153528571128845
    0.863433361053467 0.740599989891052 0.159633338451385
    0.903542876243591 0.733028590679169 0.177414283156395
    0.939257144927979 0.728785693645477 0.209957137703896
    0.972757160663605 0.729771435260773 0.239442855119705
    0.995647609233856 0.743371427059174 0.237147614359856
    0.989874601364136 0.767646014690399 0.223235711455345
    0.984101593494415 0.79192066192627 0.209323808550835
    0.978328585624695 0.816195249557495 0.195411905646324
    0.97255551815033 0.840469837188721 0.181500002741814
    0.966782510280609 0.864744484424591 0.167588099837303
    0.961009502410889 0.889019072055817 0.153676196932793
    0.963711082935333 0.912888884544373 0.137904763221741
    0.966412723064423 0.936758756637573 0.122133336961269
    0.969114303588867 0.960628569126129 0.106361903250217
    0.976899981498718 0.983900010585785 0.0804999992251396
    ];

figure();
h=dendrogram(Z,0,'Reorder',LeafOrder,'Labels',ClstrProxyNamesCompile,'orientation', 'left','ColorThreshold',CutoffVal_2);
saveas(gcf,'2D-Clustering');
close(gcf);

figure();
imagesc(flip(ClstrProxyMassesSorted));
colorbar
m=colorbar;
set(m,'fontsize',14);
m.Label.String = 'Attograms';
set(m,'LineWidth', 1);
set(m,'TickLength',0.02);
set(m,'Location','eastoutside');
XTickArray=1:NbrElements;
YTickArray=1:size(ClstrProxyNamesSorted);
set(gca,'colorscale','log')
set(gca,'XTick',XTickArray);
set(gca,'YTick',YTickArray);
xticklabels(ElementNames)
yticklabels(flip(ClstrProxyNamesSorted))
saveas(gcf,'ProxyMassHeatmap');
close(gcf);

XTickArray=1:NbrElements;
NbrProxies=size(ClstrProxyNamesSorted);
NbrProxies=NbrProxies(1);
YTickArray=1:size(ClstrProxyNamesSorted);

figure();
b=barh(ClstrNbrParticlesSorted);
b.FaceColor = [0.9290 0.6940 0.1250];
set(gca,'YTick',YTickArray);
yticklabels(ClstrProxyNamesSorted);
saveas(gcf,'ClstrNbrParticles');
close(gcf);

MaxCorr=max(Z);
MaxCorr=MaxCorr(3);
XMax=MaxCorr+0.1;
YMax=NbrProxies+0.5;

figure();
subplot(1,20,[1,2])
h=dendrogram(Z,0,'Reorder',LeafOrder,'orientation', 'left','ColorThreshold',CutoffVal_2);
ClstrAssign=cluster(Z,'Cutoff',CutoffVal_2,'criterion','distance');
set(gca,'YTick',YTickArray);
set(gca,'fontsize',14)
set(gca,'LineWidth',1.5)
set(h,'LineWidth',1.5);
axis([0 XMax 0.5 YMax]);
yticklabels([]);
axis off

subplot(1,20,[5,16],'align')
imagesc(flip(ClstrProxyMassesSorted));
CBTicksMin=10;
CBTicksMax=10^max((ceil(log10(max(ClstrProxyMassesSorted)))));
CBTickArray=10.^(linspace(log10(CBTicksMin),log10(CBTicksMax),log10(CBTicksMax)-log10(CBTicksMin)+1));
colorbar
m=colorbar;
set(m,'fontsize',14);
m.Label.String = 'Attograms';
set(m,'LineWidth', 0.5);
set(m,'TickLength',0.02);
set(m,'Ticks',CBTickArray);
set(m,'Location','eastoutside');
set(gca,'colorscale','log')
caxis([CBTicksMin CBTicksMax]);
colormap(mycolormap);
set(gca,'XTick',XTickArray);
set(gca,'YTick',YTickArray);
set(gca,'fontsize',14)
xticklabels(ElementNames)
yticklabels(flip(ClstrProxyNamesSorted))

XMin=10^(floor(log10(min(ClstrPNCsSorted))));
XMax=10^(ceil(log10(max(ClstrPNCsSorted))));
XTickArray=10.^(linspace(log10(XMin),log10(XMax),log10(XMax)-log10(XMin)+1));
subplot(1,20,[18,20], 'align')
b=barh(ClstrPNCsSorted);
b.FaceColor = [0.9290 0.6940 0.1250];
set(gca,'YTick',YTickArray);
axis([XMin XMax 0.5 YMax]);
set(gca,'YTick',[]);
set(gca,'XTick',XTickArray);
set(gca,'xscale','log')
set(gca,'fontsize',14)
set(gca,'LineWidth',1)
ax=gca;
ax.XLabel.String = 'Particles/mL';
set(gca,'fontsize',14)
saveas(gcf,'Cluster-Heatmap');

clear ClstrProxyMassesCompile ClstrProxyNamesCompile ClstrNbrParticlesCompile
clear ClstrPNCsCompile ClstrParticles3DCompile Data
clear PlotClstrOrder NbrProxies MaxCorr
clear ClstrAssign currentFile index
clear i j k h n
clear temp Z order D DataName
clear XMax XMin XTickArray YMax YTickArray
clear CBTickArray CBTicksMax CBTicksMin
clear ax b m



