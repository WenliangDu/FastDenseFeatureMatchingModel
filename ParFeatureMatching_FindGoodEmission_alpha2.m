function [GoodClassidx,tform1,tform2,QualifiedSampling,Sampled_Location1,Sampled_Location2,SearchingRange,RealSearchingRange,k,EpiLines1ab,EpiLines2ab] = ParFeatureMatching_FindGoodEmission_alpha2(Sampled_Location1,Sampled_Location2,SizeI2,Location1)
%{
2018/09/12
ParFeatureMatching_FindGoodEmission_alpha2
1. Remove outliers


2018/09/10
ParFeatureMatching_FindGoodEmission_alpha1
1. For real images
2. Find good sampling with accordant transmation angle by k-means

2018/09/10
ParFeatureMatching_Model_FindGoodEmission_alpha1
1. For simulated datasets
2. Find good sampling with accordant transmation angle by k-means


%}
k = 0;
QualifiedSampling = false;
GoodClassidx = [];
tform1 = [];
tform2 = [];
EpiLines1ab = [];
EpiLines2ab = [];
SearchingRange = [];
RealSearchingRange = [];

Sampled_Location1Old = Sampled_Location1;
Sampled_Location2Old = Sampled_Location2;

while ~QualifiedSampling
    k = k + 1;
%     [idx,~,~,~] = kmeans(Vector_Location,k,'Distance','cosine');
    Vector_Location = Sampled_Location1 - Sampled_Location2;
    idx = kmeans(Vector_Location,k,'Distance','cosine');
    %% Get class number
    ClassN = zeros(1,k);
    Classidx = cell(1,k);
    for i = 1:k
        Classidx{i} = (idx==i);
        ClassN(i) = sum(Classidx{i});
    end
    if max(ClassN) < 30
        x = 1;
        break
    end
    %%
    tformTemp = cell(k,2);
    EpiLinesabTemp = cell(k,2);
    QualifiedSamplingTemp = false(1,k);
    SearchingRangeTemp = Inf(1,k);
    DifferMaxTemp = Inf(1,k);
    MinimumSearchingRangeTemp = zeros(1,k);
    WrongIdxTemp = false(size(Sampled_Location1,1),1);
%     figure,
    for i = 1:k
        CurrentWrong = true;
        
        if ClassN(i) >= 30
            
            while CurrentWrong && ClassN(i) > 30
                Current_Sampled_Location1 = Sampled_Location1(Classidx{i},:);
                Current_Sampled_Location2 = Sampled_Location2(Classidx{i},:);
                %%

                x = 1;
                %%
                Current_Fundemantal = RANSAC_Norm8Point_alpha2(Current_Sampled_Location1, Current_Sampled_Location2);
                [t1, t2] = estimateUncalibratedRectification(Current_Fundemantal, Current_Sampled_Location1, Current_Sampled_Location2,SizeI2);

                tformTemp{i,1} = projective2d(t1);
                tformTemp{i,2} = projective2d(t2);

                Sampled_ProjectedLocation1 = transformPointsForward(tformTemp{i,1}, Sampled_Location1);
                Sampled_ProjectedLocation2 = transformPointsForward(tformTemp{i,2}, Sampled_Location2);

                %%
                Current_Sampled_ProjectedLocation1 = Sampled_ProjectedLocation1(Classidx{i},:);
                Current_Sampled_ProjectedLocation2 = Sampled_ProjectedLocation2(Classidx{i},:);

                Qualified_Sampled_ProjectedLocation1 = Sampled_ProjectedLocation1(~WrongIdxTemp,:);
                Qualified_Sampled_ProjectedLocation2 = Sampled_ProjectedLocation2(~WrongIdxTemp,:);

                [EpiLinesabTemp{i,1},EpiLinesabTemp{i,2}] = CalculateEpipolarLine_alpha2(Current_Sampled_ProjectedLocation1,Current_Sampled_ProjectedLocation2,Qualified_Sampled_ProjectedLocation1,Qualified_Sampled_ProjectedLocation2);



                MayQualified1Temp = (abs(EpiLinesabTemp{i,1}(:,1))<0.01);
                MayQualified2Temp = (abs(EpiLinesabTemp{i,2}(:,1))<0.01);
                MayQualifiedPairsTemp = (MayQualified1Temp&MayQualified2Temp);
                QualifiedSamplingTemp(i) = all(MayQualifiedPairsTemp);

                %% Remove outliers
                ProjectedLocation1 = transformPointsForward(tformTemp{i,1}, Location1);
                SearchingRange = ParFeatureMatching_CalcuSearchingRange_alpha2(ProjectedLocation1,EpiLinesabTemp{i,1});
                MinimumSearchingRange = ((max(ProjectedLocation1(:,2)) - min(ProjectedLocation1(:,2)))/size(ProjectedLocation1,1))*50;
                Diff = abs(Sampled_ProjectedLocation1(:,2) - Sampled_ProjectedLocation2(:,2));
                WrongIdx = Diff>max([SearchingRange MinimumSearchingRange]);

                if any(WrongIdx)              
                    WrongIdxTemp = WrongIdxTemp | WrongIdx;
                    for j = 1:k
                        if i==j && any(WrongIdx(Classidx{j}))
                            x = 1;
%                             CurrentWrong = true;
                        else
                            CurrentWrong = false;
                        end
                        x = 1;
                        ClassN(j) = ClassN(j) - sum(Classidx{j}(WrongIdx));
                        Classidx{j}(WrongIdx) = false;
                        x = 1;
                    end

                    x =1;
                else
                    CurrentWrong = false;
                end
            end
            
            x = 1;
            
            if QualifiedSamplingTemp(i)
                Qualified_Sampled_ProjectedLocation1 = Sampled_ProjectedLocation1(~WrongIdxTemp,:);
                Qualified_Sampled_ProjectedLocation2 = Sampled_ProjectedLocation2(~WrongIdxTemp,:);
                
                DifferMaxTemp(i) = max(abs(Qualified_Sampled_ProjectedLocation1(:,2)-Qualified_Sampled_ProjectedLocation2(:,2)));
%                 ProjectedLocation1 = transformPointsForward(tformTemp{i,1}, Location1);
                SearchingRangeTemp(i) = SearchingRange;
                MinimumSearchingRangeTemp(i) = MinimumSearchingRange;
%                 SearchingRangeTemp(i) = ParFeatureMatching_CalcuSearchingRange_alpha2(ProjectedLocation1,EpiLinesabTemp{i,1});
%                 MinimumSearchingRangeTemp(i) = ((max(ProjectedLocation1(:,2)) - min(ProjectedLocation1(:,2)))/size(ProjectedLocation1,1))*50;
                
%                 Diff = abs(Sampled_ProjectedLocation1(:,2) - Sampled_ProjectedLocation2(:,2));
%                 max(Diff)
                
%                 figure,
%                 plot(Current_Sampled_Location1(:,1),Current_Sampled_Location1(:,2),'r*');hold on
%                 plot(Current_Sampled_Location2(:,1),Current_Sampled_Location2(:,2),'k*');
%                 line([Current_Sampled_Location1(:,1)';Current_Sampled_Location2(:,1)'],[Current_Sampled_Location1(:,2)';Current_Sampled_Location2(:,2)'],'Color','g');
%                 
%                 figure,
%                 plot(Qualified_Sampled_ProjectedLocation1(:,1),Qualified_Sampled_ProjectedLocation1(:,2),'r*');hold on
%                 plot(Qualified_Sampled_ProjectedLocation2(:,1),Qualified_Sampled_ProjectedLocation2(:,2),'k*');
%                 line([Qualified_Sampled_ProjectedLocation1(:,1)';Qualified_Sampled_ProjectedLocation2(:,1)'],[Qualified_Sampled_ProjectedLocation1(:,2)';Qualified_Sampled_ProjectedLocation2(:,2)'],'Color','g');
                x = 1;
            end
            
        end
    end
    if any(WrongIdxTemp)
        Sampled_Location1(WrongIdxTemp,:) = [];
        Sampled_Location2(WrongIdxTemp,:) = [];
        for j = 1:k
            Classidx{j}(WrongIdxTemp) = [];
        end
    end

    if any(QualifiedSamplingTemp)
        QualifiedSampling = true;
        [~,GoodClass] = min(DifferMaxTemp);
        GoodClassidx = Classidx{GoodClass};
        tform1 = tformTemp{GoodClass,1};
        tform2 = tformTemp{GoodClass,2};
        EpiLines1ab = EpiLinesabTemp{GoodClass,1};
        EpiLines2ab = EpiLinesabTemp{GoodClass,2};
        SearchingRange = SearchingRangeTemp(GoodClass);
        RealSearchingRange = max([SearchingRange MinimumSearchingRangeTemp(GoodClass)]);
%         x = 1;
        %% get real searching range
        
        %%
        Sampled_ProjectedLocation1 = transformPointsForward(tformTemp{GoodClass,1}, Sampled_Location1);
        Sampled_ProjectedLocation2 = transformPointsForward(tformTemp{GoodClass,2}, Sampled_Location2);
        
%         figure,
%         plot(Sampled_Location1Old(:,1),Sampled_Location1Old(:,2),'r*');hold on
%         plot(Sampled_Location2Old(:,1),Sampled_Location2Old(:,2),'k*');
%         line([Sampled_Location1Old(:,1)';Sampled_Location2Old(:,1)'],[Sampled_Location1Old(:,2)';Sampled_Location2Old(:,2)'],'Color','g');
%         
%         figure,
%         plot(Sampled_Location1(:,1),Sampled_Location1(:,2),'r*');hold on
%         plot(Sampled_Location2(:,1),Sampled_Location2(:,2),'k*');
%         line([Sampled_Location1(:,1)';Sampled_Location2(:,1)'],[Sampled_Location1(:,2)';Sampled_Location2(:,2)'],'Color','g');
%         
%         figure,
%         plot(Sampled_ProjectedLocation1(:,1),Sampled_ProjectedLocation1(:,2),'r*');hold on
%         plot(Sampled_ProjectedLocation2(:,1),Sampled_ProjectedLocation2(:,2),'k*');
%         line([Sampled_ProjectedLocation1(:,1)';Sampled_ProjectedLocation2(:,1)'],[Sampled_ProjectedLocation1(:,2)';Sampled_ProjectedLocation2(:,2)'],'Color','g');
%         RemainRatio = (max(Sampled_Location1(:,2))-min(Sampled_Location1(:,2)))/(max(Sampled_Location1Old(:,2))-min(Sampled_Location1Old(:,2)))
%         x = 1;
    end
end

%%
