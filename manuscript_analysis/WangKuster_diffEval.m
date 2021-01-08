work_dir='/labs/halperin/AlternativeSplicing/WangKuster2019/';
pept_detect_anno='WangKuster.top.MSDetection.bisbeeDiff.csv'
counts_dir='bisbee/counts/WangKuster.'
spladder_test_dir='spladder2_test/testing*/test_results_extended_C3_'
sample_table='bisbee/WangKuster_testSamples.txt'

cd(work_dir)

%%% calculate t-test p-values and fdr on PSI values and extract mean PSI
%%% and mean read depths per group
peptides=readtable(pept_detect_anno);
ll_idx=startsWith(peptides.Properties.VariableNames,'ll_ratio');
comp_list=strrep(peptides.Properties.VariableNames(ll_idx),'ll_ratio_','');
peptides.ttest_p=nan(height(peptides),length(comp_list));
peptides.mean_psi_diff=nan(height(peptides),length(comp_list));
peptides.mean_psi_diff_conf=nan(height(peptides),length(comp_list));
peptides.mean_depth=nan(height(peptides),length(comp_list));
countsFiles=struct2table(dir([counts_dir '*.csv']))
sample_table=readtable(sample_table,'ReadVariableNames',0,'Delimiter','\t');
for i=1:height(countsFiles)
    opt=detectImportOptions([countsFiles.folder{i} '/' countsFiles.name{i}],'FileType','text');
    idx=find(strcmp(opt.VariableNames,'event_jid'));
    opt.SelectedVariableNames=opt.VariableNames(idx:end);
    countsData=readtable([countsFiles.folder{i} '/' countsFiles.name{i}],opt);
    iso1_idx=endsWith(countsData.Properties.VariableNames,'iso1');
    iso2_idx=endsWith(countsData.Properties.VariableNames,'iso2');
    PSI=countsData{:,iso1_idx}./(countsData{:,iso1_idx}+countsData{:,iso2_idx});
    meanDepth=mean(countsData{:,iso1_idx}+countsData{:,iso2_idx},2);
    PSI_D10=PSI;
    PSI_D10((countsData{:,iso1_idx}+countsData{:,iso2_idx})<10)=NaN;
    sample_names=strrep(countsData.Properties.VariableNames(iso1_idx),'_iso1','');
    for j=1:length(comp_list)
        groups=strsplit(comp_list{j},'vs');
        [~,locb]=ismember(sample_table{:,2},groups);
        set1=ismember(sample_names,strrep(sample_table{locb==1,1},'.','_'));
        set2=ismember(sample_names,strrep(sample_table{locb==2,1},'.','_'));  
        ttest_p=mattest(PSI(:,set1),PSI(:,set2));
        ttest_p_conf=mattest(PSI_D10(:,set1),PSI_D10(:,set2));
        [lia,locb]=ismember(peptides.event_jid,countsData.event_jid);
        peptides.ttest_p(lia,j)=ttest_p(locb(lia));
        peptides.ttest_p_conf(lia,j)=ttest_p_conf(locb(lia));
        peptides.mean_psi_diff(lia,j)=nanmean(PSI(locb(lia),set1),2)-nanmean(PSI(locb(lia),set2),2);
        peptides.mean_psi_diff_conf(lia,j)=nanmean(PSI_D10(locb(lia),set1),2)-nanmean(PSI_D10(locb(lia),set2),2);      
        peptides.mean_depth(lia,j)=meanDepth(locb(lia),:);
    end
    i
end
peptides.ttest_fdr=nan(size(peptides.ttest_p));
peptides.ttest_fdr_conf=nan(size(peptides.ttest_p));
for i=1:size(peptides.ttest_p,2)
    peptides.ttest_fdr(:,i)=mafdr(peptides.ttest_p(:,i));
    peptides.ttest_fdr_conf(:,i)=mafdr(peptides.ttest_p_conf(:,i));
end
peptides.bb_lr=peptides{:,ll_idx};

%%%% identify true tissue specific events by MS
msTissues={'Brain', 'Colon','Duodenum','Ovary','Rectum','SI','Tonsil'};
for i=1:length(comp_list)
    groups=strsplit(comp_list{i},'vs');
    msPairs(i,1)=find(strcmp(msTissues,groups{1}));
    msPairs(i,2)=find(strcmp(msTissues,groups{2}));
end
peptides.msSwitch=zeros(height(peptides),length(msPairs));
iso1_idx=find(startsWith(peptides.Properties.VariableNames,'iso1_maxAb'));
iso2_idx=find(startsWith(peptides.Properties.VariableNames,'iso2_maxAb'));
peptides.iso1_only=peptides{:,iso1_idx}>0 & isnan(peptides{:,iso2_idx});
peptides.iso2_only=peptides{:,iso2_idx}>0 & isnan(peptides{:,iso1_idx});
for i=1:length(msPairs)
    peptides.msSwitch(peptides.iso1_only(:,msPairs(i,1)) & peptides.iso2_only(:,msPairs(i,2)),i)=1;
    peptides.msSwitch(peptides.iso2_only(:,msPairs(i,1)) & peptides.iso1_only(:,msPairs(i,2)),i)=-1;
end

%%% read in spladder test results
spladderFiles=struct2table(dir([spladder_test_dir '*.tsv']));
spladderFiles.comp=extractAfter(spladderFiles.folder,'testing_');
peptides.spladder_p=nan(height(peptides),length(comp_list));
peptides.spladder_p_adj=nan(height(peptides),length(comp_list));
for i=18:height(spladderFiles)
    opt=detectImportOptions([spladderFiles.folder{i} '/' spladderFiles.name{i}],'FileType','text');
    opt=setvartype(opt,3:22,'double');
    spladderRes=readtable([spladderFiles.folder{i} '/' spladderFiles.name{i}],opt);
    cIdx=strcmpi(spladderFiles.comp(i),strrep(comp_list,'vs','_'));
    if sum(cIdx)==0
        curr_tissues=strsplit(spladderFiles.comp{i},'_');
        curr_comp=strcat(curr_tissues{2},'vs',curr_tissues{1});
        cIdx=strcmpi(curr_comp,comp_list);
    end
    [lia,locb]=ismember(peptides.event_id,spladderRes.event_id);
    peptides.spladder_p(lia,cIdx)=spladderRes.p_val(locb(lia));
    peptides.spladder_p_adj(lia,cIdx)=spladderRes.p_val_adj(locb(lia));
end

%%% ROC curve analysis
[ppv_bb,tpr_bb,t_bb,auc_bb]=perfcurve(abs(peptides.msSwitch(:)),peptides.bb_lr(:),1,'XCrit','ppv','TVals',[0:0.1:20]);
[ppv_sp,tpr_sp,t_sp,auc_sp]=perfcurve(abs(peptides.msSwitch(:)),-log(peptides.spladder_p_adj(:)+realmin),1,'XCrit','ppv','TVals',[0:0.1:20]);
[ppv_tt,tpr_tt,t_tt,auc_tt]=perfcurve(abs(peptides.msSwitch(:)),-log(peptides.ttest_p(:)+realmin),1,'XCrit','ppv','TVals',[0:0.1:20]);
[ppv_tt_conf,tpr_tt_conf,t_tt_conf,auc_tt_conf]=perfcurve(abs(peptides.msSwitch(:)),-log(peptides.ttest_p_conf(:)+realmin),1,'XCrit','ppv','TVals',[0:0.1:20]);
[fp_bb,tp_bb]=perfcurve(abs(peptides.msSwitch(:)),peptides.bb_lr(:),1,'XCrit','fp','YCrit','tp','TVals',[0:0.1:20]);
[fp_sp,tp_sp]=perfcurve(abs(peptides.msSwitch(:)),-log(peptides.spladder_p_adj(:)+realmin),1,'XCrit','fp','YCrit','tp','TVals',[0:0.1:20]);
[fp_tt,tp_tt]=perfcurve(abs(peptides.msSwitch(:)),-log(peptides.ttest_p(:)+realmin),1,'XCrit','fp','YCrit','tp','TVals',[0:0.1:20]);
[fp_tt_conf,tp_tt_conf]=perfcurve(abs(peptides.msSwitch(:)),-log(peptides.ttest_p_conf(:)+realmin),1,'XCrit','fp','YCrit','tp','TVals',[0:0.1:20]);

%%%% plot figure 3a
plot(log(fp_bb+tp_bb),log(tp_bb),'k');
hold on;
plot(log(fp_sp+tp_sp),log(tp_sp),'r');
plot(log(fp_tt+tp_tt),log(tp_tt),'b');
plot(log(fp_tt_conf+tp_tt_conf),log(tp_tt_conf),'c');

set(gca,'XTick',log(10.^[1:8]),'XTickLabel',10.^[1:8],'FontSize',6);
set(gca,'YTick',log(2.^[1:8]),'YTickLabel',2.^[1:8]);
title('Differential Splicing Evaluation')
ylabel('Confirmed Events Passing Thresh');
xlabel('Total Events Passing Thresh');
legend({'bbd','sp','tt','tt-d10'},'Location','northwest');

set(gcf,'PaperSize',[3.3 3.3]);
print('Fig3a_TPvsTotal.pdf','-dpdf','-r300','-fillpage');
print('Fig3a_TPvsTotal.svg','-dsvg','-r300');
print('Fig3a_TPvsTotal.png','-dpng','-r300');


%%% plot figure 4
subplot(2,2,1)
peptides.ll_ratio_BrainvsSI(isnan(peptides.ll_ratio_BrainvsSI))=0;
s=scatter(peptides.mean_psi_diff(:,5),peptides.ll_ratio_BrainvsSI,20,log(peptides.mean_depth(:,5)),'filled')
colormap('jet');
s.MarkerFaceAlpha=0.5;
hold on;
tIdx=abs(peptides.msSwitch(:,5))==1;
scatter(peptides.mean_psi_diff(tIdx,5),peptides.ll_ratio_BrainvsSI(tIdx),20,log(peptides.mean_depth(tIdx,5)),'filled')
scatter(peptides.mean_psi_diff(tIdx,5),peptides.ll_ratio_BrainvsSI(tIdx),20,'ko');
ylim([0 20]);
c=colorbar('Ticks',log(10.^[0:5]),'TickLabels',10.^[0:5]);
c.Label.String='mean depth';
c.Label.Position=[2.3429 3.4539 0];
xlabel('mean(PSI Brain) - mean(PSI SI)');
ylabel('BBD LR Brain vs SI');
caxis(log([1 1E3]))
title('Bisbee Diff');
subplot(2,2,2)
peptides.spladder_p_adj(isnan(peptides.spladder_p_adj))=1;
s=scatter(peptides.mean_psi_diff(:,5),-log(peptides.spladder_p_adj(:,5)+1E-15),20,log(peptides.mean_depth(:,5)),'filled')
colormap('jet');
s.MarkerFaceAlpha=0.5;
hold on;
tIdx=abs(peptides.msSwitch(:,5))==1;
scatter(peptides.mean_psi_diff(tIdx,5),-log(peptides.spladder_p_adj(tIdx,5)+1E-15),20,log(peptides.mean_depth(tIdx,5)),'filled')
scatter(peptides.mean_psi_diff(tIdx,5),-log(peptides.spladder_p_adj(tIdx,5)+1E-15),20,'ko');
set(gca,'YTick',fliplr(-log(10.^[-14:2:0])),'YTickLabel',fliplr(10.^[-14:2:0]));
c=colorbar('Ticks',log(10.^[0:5]),'TickLabels',10.^[0:5]);
c.Label.String='mean depth';
c.Label.Position=[2.3429 3.4539 0];
xlabel('mean(PSI Brain) - mean(PSI SI)');
ylabel('SplAdder p-adj');
caxis(log([1 1E3]))
title('SplAdder test');
axis tight;
subplot(2,2,3)
peptides.ttest_p(isnan(peptides.ttest_p))=1;
s=scatter(peptides.mean_psi_diff(:,5),-log(peptides.ttest_p(:,5)),20,log(peptides.mean_depth(:,5)),'filled')
colormap('jet');
s.MarkerFaceAlpha=0.5;
hold on;
tIdx=abs(peptides.msSwitch(:,5))==1;
scatter(peptides.mean_psi_diff(tIdx,5),-log(peptides.ttest_p(tIdx,5)),20,log(peptides.mean_depth(tIdx,5)),'filled')
scatter(peptides.mean_psi_diff(tIdx,5),-log(peptides.ttest_p(tIdx,5)),20,'ko');
set(gca,'YTick',fliplr(-log(10.^[-6:0])),'YTickLabel',fliplr(10.^[-6:0]));
c=colorbar('Ticks',log(10.^[0:5]),'TickLabels',10.^[0:5]);
c.Label.String='mean depth';
c.Label.Position=[2.3429 3.4539 0];
xlabel('mean(PSI Brain) - mean(PSI SI)');
ylabel('t-test pvalue');
caxis(log([1 1E3]))
title('t-test all');
subplot(2,2,4)
peptides.ttest_p_conf(isnan(peptides.ttest_p_conf))=1;
scatter(peptides.mean_psi_diff(:,5),-log(peptides.ttest_p_conf(:,5)),20,log(peptides.mean_depth(:,5)),'filled')
colormap('jet');
s.MarkerFaceAlpha=0.5;
hold on;
tIdx=abs(peptides.msSwitch(:,5))==1;
scatter(peptides.mean_psi_diff(tIdx,5),-log(peptides.ttest_p_conf(tIdx,5)),20,log(peptides.mean_depth(tIdx,5)),'filled')
scatter(peptides.mean_psi_diff(tIdx,5),-log(peptides.ttest_p_conf(tIdx,5)),20,'ko');
set(gca,'YTick',fliplr(-log(10.^[-6:0])),'YTickLabel',fliplr(10.^[-6:0]));
c=colorbar('Ticks',log(10.^[0:5]),'TickLabels',10.^[0:5]);
c.Label.String='mean depth';
c.Label.Position=[2.3429 3.4539 0];
xlabel('mean(PSI Brain) - mean(PSI SI)');
ylabel('t-test D10 pvalue');
title('t-test depth>10');
caxis(log([1 1E3]))

set(gcf,'PaperSize',[8.5 8]);
print('Fig4_diffEval_volcanoPlots.pdf','-dpdf','-r300','-fillpage');
print('Fig4_diffEval_volcanoPlots.svg','-dsvg','-r300');
print('Fig4_diffEval_volcanoPlots.png','-dpng','-r300');

close(gcf);
save('WangKuster_diffEval.mat');

%%% create table S1
summary=table();
summary.testName={'bb';'sp';'tt';'tt-D10'};
summary.thresh1={'lr>8';'p<01';'p<01';'p<01'};
summary.conf_pass1(1)=sum(abs(peptides.bb_lr(abs(peptides.msSwitch)==1))>8);
summary.conf_pass1(2)=sum(peptides.spladder_p(abs(peptides.msSwitch)==1)<0.01);
summary.conf_pass1(3)=sum(peptides.ttest_p(abs(peptides.msSwitch)==1)<0.01);
summary.conf_pass1(4)=sum(peptides.ttest_p_conf(abs(peptides.msSwitch)==1)<0.01);
summary.tot_pass1(1)=sum(abs(peptides.bb_lr(:))>8);
summary.tot_pass1(2)=sum(peptides.spladder_p(:)<0.01);
summary.tot_pass1(3)=sum(peptides.ttest_p(:)<0.01);
summary.tot_pass1(4)=sum(peptides.ttest_p_conf(:)<0.01);

summary.thresh2={'lr>12';'p_adj<05';'fdr<05';'fdr<05'};
summary.conf_pass2(1)=sum(abs(peptides.bb_lr(abs(peptides.msSwitch)==1))>12);
summary.conf_pass2(2)=sum(peptides.spladder_p_adj(abs(peptides.msSwitch)==1)<0.05);
summary.conf_pass2(3)=sum(peptides.ttest_fdr(abs(peptides.msSwitch)==1)<0.05);
summary.conf_pass2(4)=sum(peptides.ttest_fdr_conf(abs(peptides.msSwitch)==1)<0.05);
summary.tot_pass2(1)=sum(abs(peptides.bb_lr(:))>12);
summary.tot_pass2(2)=sum(peptides.spladder_p_adj(:)<0.05);
summary.tot_pass2(3)=sum(peptides.ttest_fdr(:)<0.05);
summary.tot_pass2(4)=sum(peptides.ttest_fdr_conf(:)<0.05);
writetable(summary,'WangKuster_diffEval_threshSummaryTable.csv');

save('WangKuster_diffEval.mat');

