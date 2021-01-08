work_dir='/labs/halperin/AlternativeSplicing/WangKuster2019/';
pept_detect='WangKuster.top.MSdetectionAll.csv'
outscores_dir='bisbee/outlier/scores/WangKuster.'
ref_counts='bisbee/counts/lowerGI'
test_counts='bisbee/counts/WangKuster'

cd(work_dir)

%%% read ms peptide detection data
peptides=readtable(pept_detect);

%%% read outlier scores
bbFiles=struct2table(dir([outscores_dir '*.csv']));
bbRes={};
for i=1:height(bbFiles)
    filename=[bbFiles.folder{i} '/' bbFiles.name{i}];
    opt=detectImportOptions(filename);
    idx=find(strcmp(opt.VariableNames,'event_jid'));
    opt.SelectedVariableNames=opt.VariableNames(idx:end);
    opt=setvartype(opt,opt.SelectedVariableNames(2:end),'double');
    bbRes{i}=readtable(filename,opt);
end
bbRes=vertcat(bbRes{:}); 

%%% define tissue order and gi and other tissues
msTissues={'Brain', 'Colon','Duodenum','Ovary','Rectum','SmallIntestine','Tonsil'};
msGIidx=[2 3 5 6];
msOidx=[1 4 7];

%%% add outlier scores to peptide table
[lia,locb]=ismember(peptides.event_jid,bbRes.event_jid);
peptides.bb_lr=nan(height(peptides),width(bbRes)-1);
peptides.bb_lr(lia,:)=bbRes{locb(lia),2:end};
clear bbRes;

%%% determine true outliers by ms
peptides.msOut=zeros(height(peptides),length(msOidx));
iso1_idx=find(startsWith(peptides.Properties.VariableNames,'iso1_maxAb'));
iso2_idx=find(startsWith(peptides.Properties.VariableNames,'iso2_maxAb'));
peptides.iso1_only=peptides{:,iso1_idx}>0 & isnan(peptides{:,iso2_idx});
peptides.iso2_only=peptides{:,iso2_idx}>0 & isnan(peptides{:,iso1_idx});
for i=1:length(msOidx)
    peptides.msOut(max(peptides.iso2_only(:,msGIidx),[],2)==1 & sum(~isnan(peptides{:,iso1_idx(msGIidx)}),2)==0  & peptides.iso1_only(:,msOidx(i)),i)=1;
    peptides.msOut(max(peptides.iso1_only(:,msGIidx),[],2)==1 & sum(~isnan(peptides{:,iso2_idx(msGIidx)}),2)==0  & peptides.iso2_only(:,msOidx(i)),i)=-1;
end

%%% read reference sample counts
refFiles=struct2table(dir([ref_counts '*.csv']));
refCounts={};
for i=1:height(refFiles)
    opt=detectImportOptions([refFiles.folder{i} '/' refFiles.name{i}],'FileType','text');
    iso1_idx=find(endsWith(opt.VariableNames,'_iso1'));
    iso2_idx=find(endsWith(opt.VariableNames,'_iso2'));
    opt=setvartype(opt,opt.VariableNames([iso1_idx iso2_idx]),'double');
    opt.SelectedVariableNames=['event_jid' opt.VariableNames([iso1_idx iso2_idx])];
    refCounts{i}=readtable([refFiles.folder{i} '/' refFiles.name{i}],opt);
    i
end
refCounts=vertcat(refCounts{:});

%%% calculate reference sample PSI
iso1_idx=find(endsWith(refCounts.Properties.VariableNames,'_iso1'));
iso2_idx=find(endsWith(refCounts.Properties.VariableNames,'_iso2'));
refPSI=refCounts{:,iso1_idx}./(refCounts{:,iso1_idx}+refCounts{:,iso2_idx});
refPSI_D10=refPSI;
refPSI_D10(refCounts{:,iso1_idx}+refCounts{:,iso2_idx}<10)=NaN;

%%% read in test sample counts
testFiles=struct2table(dir([test_counts '*.csv']));
testCounts={};
for i=1:height(testFiles)
    opt=detectImportOptions([testFiles.folder{i} '/' testFiles.name{i}],'FileType','text');
    iso1_idx=find(endsWith(opt.VariableNames,'_iso1'));
    iso2_idx=find(endsWith(opt.VariableNames,'_iso2'));
    opt=setvartype(opt,opt.VariableNames([iso1_idx iso2_idx]),'double');
    opt.SelectedVariableNames=['event_jid' opt.VariableNames([iso1_idx iso2_idx])];
    testCounts{i}=readtable([testFiles.folder{i} '/' testFiles.name{i}],opt);
    i
end
testCounts=vertcat(testCounts{:});

%%% calculate test sample PSI
iso1_idx=find(endsWith(testCounts.Properties.VariableNames,'_iso1'));
iso2_idx=find(endsWith(testCounts.Properties.VariableNames,'_iso2'));
testPSI=testCounts{:,iso1_idx}./(testCounts{:,iso1_idx}+testCounts{:,iso2_idx});
testPSI_D10=testPSI;
testPSI_D10(testCounts{:,iso1_idx}+testCounts{:,iso2_idx}<10)=NaN;

%%% verify event_jid columns are identical between test and ref sets
sum(~strcmp(refCounts.event_jid,testCounts.event_jid))

%%% find median absolute deviation and interquartile outlier scores
mad=(testPSI-nanmedian(refPSI,2))./max(nanmedian(abs(refPSI-nanmedian(refPSI,2)),2),0.01);
q25=quantile(refPSI,0.25,2);
q75=quantile(refPSI,0.75,2);
iqrNorm=zeros(size(testPSI));
for i=1:size(iqrNorm,2)
    lIdx=testPSI(:,i)<q25;
    iqrNorm(lIdx,i)=(testPSI(lIdx,i)-q25(lIdx))./max(q75(lIdx)-q25(lIdx),0.01);
    hIdx=testPSI(:,i)>q75;
    iqrNorm(hIdx,i)=(testPSI(hIdx,i)-q75(hIdx))./max(q75(hIdx)-q25(hIdx),0.01);
end

%%% find median absolute deviation and interquartile outlier scores
%%% excluding PSI values with coverage less than 10
mad_conf=(testPSI_D10-nanmedian(refPSI_D10,2))./max(nanmedian(abs(refPSI_D10-nanmedian(refPSI_D10,2)),2),0.01);
q25_conf=quantile(refPSI_D10,0.25,2);
q75_conf=quantile(refPSI_D10,0.75,2);
iqrNorm_conf=zeros(size(testPSI_D10));
for i=1:size(iqrNorm,2)
    lIdx=testPSI_D10(:,i)<q25_conf;
    iqrNorm_conf(lIdx,i)=(testPSI_D10(lIdx,i)-q25_conf(lIdx))./max(q75_conf(lIdx)-q25_conf(lIdx),0.01);
    hIdx=testPSI_D10(:,i)>q75_conf;
    iqrNorm_conf(hIdx,i)=(testPSI_D10(hIdx,i)-q75_conf(hIdx))./max(q75_conf(hIdx)-q25_conf(hIdx),0.01);
end

%%% add mad and iqr scores to peptide table
[lia,locb]=ismember(peptides.event_jid,testCounts.event_jid);
peptides.mad=mad(locb,:);
peptides.iqrNorm=iqrNorm(locb,:);
peptides.mad_conf=mad_conf(locb,:);
peptides.iqrNorm_conf=iqrNorm_conf(locb,:);
samples=strrep(testCounts.Properties.VariableNames(iso1_idx),'_iso1','');
clear testCounts refCounts;

%%% make table of true outliers
[~,tIdx]=ismember(extractBefore(samples,'_'),lower(msTissues));
peptides.trueOut=zeros(size(peptides.bb_lr));
for i=1:length(msOidx)
    idx=find(msOidx(i)==tIdx);
    peptides.trueOut(:,idx)=repmat(peptides.msOut(:,i),1,length(idx));
end

%%% ROC curve analysis
[ppv_bb,tpr_bb,t_bb,auc_bb]=perfcurve(abs(peptides.trueOut(:)),abs(peptides.bb_lr(:)),1,'XCrit','ppv','TVals',[0:600]);
[ppv_mad,tpr_mad,t_mad,auc_mad]=perfcurve(abs(peptides.trueOut(:)),peptides.mad(:),1,'XCrit','ppv','TVals',[0:0.25:100]);
[ppv_mad_conf,tpr_mad_conf,t_mad_conf,auc_mad_conf]=perfcurve(abs(peptides.trueOut(:)),peptides.mad_conf(:),1,'XCrit','ppv','TVals',[0:0.25:100]);
[ppv_iqr,tpr_iqr,t_iqr,auc_iqr]=perfcurve(abs(peptides.trueOut(:)),peptides.iqrNorm(:),1,'XCrit','ppv','TVals',[0:0.25:100]);
[ppv_iqr_conf,tpr_iqr_conf,t_iqr_conf,auc_iqr_conf]=perfcurve(abs(peptides.trueOut(:)),peptides.iqrNorm_conf(:),1,'XCrit','ppv','TVals',[0:0.25:100]);

[fp_bb,tp_bb]=perfcurve(abs(peptides.trueOut(:)),abs(peptides.bb_lr(:)),1,'XCrit','fp','YCrit','tp','TVals',[0:600]);
[fp_mad,tp_mad]=perfcurve(abs(peptides.trueOut(:)),peptides.mad(:),1,'XCrit','fp','YCrit','tp','TVals',[0:0.25:100]);
[fp_mad_conf,tp_mad_conf]=perfcurve(abs(peptides.trueOut(:)),peptides.mad_conf(:),1,'XCrit','fp','YCrit','tp','TVals',[0:0.25:100]);
[fp_iqr,tp_iqr]=perfcurve(abs(peptides.trueOut(:)),peptides.iqrNorm(:),1,'XCrit','fp','YCrit','tp','TVals',[0:0.25:100]);
[fp_iqr_conf,tp_iqr_conf]=perfcurve(abs(peptides.trueOut(:)),peptides.iqrNorm_conf(:),1,'XCrit','fp','YCrit','tp','TVals',[0:0.25:100]);

%%% plot figure 3b
plot(log(fp_bb+tp_bb),log(tp_bb),'k');
hold on;
plot(log(fp_mad+tp_mad),log(tp_mad),'r');
plot(log(fp_mad_conf+tp_mad_conf),log(tp_mad_conf),'m');
plot(log(fp_iqr+tp_iqr),log(tp_iqr),'b');
plot(log(fp_iqr_conf+tp_iqr_conf),log(tp_iqr_conf),'c');
set(gca,'XTick',log(10.^[1:8]),'XTickLabel',10.^[1:8],'FontSize',6);
set(gca,'YTick',log(2.^[1:8]),'YTickLabel',2.^[1:8]);
title('Splice Outlier Evaluation');
ylabel('Confirmed Events Passing Thresh');
xlabel('Total Events Passing Thresh');
legend({'bbo','mad','mad-d10','iqr','iqr-d10'},'Location','northwest');

set(gcf,'PaperSize',[3.3 3.3]);
print('Fig3b_outlierEval_TPvsTotal.pdf','-dpdf','-r300','-fillpage');
print('Fig3b_outlierEval_TPvsTotal.svg','-dsvg','-r300');
print('Fig3b_outlierEval_TPvsTotal.png','-dpng','-r300');
close(gcf);

%%%% create summary table
summary=table();
summary.testName={'bb';'mad';'mad-D10';'iqr';'iqr-D10'};
summary.thresh1={'10';'20';'20';'10';'8'};
summary.conf_pass1(1)=sum(abs(peptides.bb_lr(abs(peptides.trueOut)==1))>10);
summary.conf_pass1(2)=sum(peptides.mad(abs(peptides.trueOut)==1)>20);
summary.conf_pass1(3)=sum(peptides.mad_conf(abs(peptides.trueOut)==1)>20);
summary.conf_pass1(4)=sum(peptides.iqrNorm(abs(peptides.trueOut)==1)>10);
summary.conf_pass1(5)=sum(peptides.iqrNorm_conf(abs(peptides.trueOut)==1)>8);
summary.tot_pass1(1)=sum(abs(peptides.bb_lr(:))>10);
summary.tot_pass1(2)=sum(peptides.mad(:)>20);
summary.tot_pass1(3)=sum(peptides.mad_conf(:)>20);
summary.tot_pass1(4)=sum(peptides.iqrNorm(:)>10);
summary.tot_pass1(5)=sum(peptides.iqrNorm_conf(:)>8);
summary.thresh2={'14';'40';'25';'20';'10'};
summary.conf_pass2(1)=sum(abs(peptides.bb_lr(abs(peptides.trueOut)==1))>14);
summary.conf_pass2(2)=sum(peptides.mad(abs(peptides.trueOut)==1)>40);
summary.conf_pass2(3)=sum(peptides.mad_conf(abs(peptides.trueOut)==1)>25);
summary.conf_pass2(4)=sum(peptides.iqrNorm(abs(peptides.trueOut)==1)>20);
summary.conf_pass2(5)=sum(peptides.iqrNorm_conf(abs(peptides.trueOut)==1)>10);
summary.tot_pass2(1)=sum(abs(peptides.bb_lr(:))>14);
summary.tot_pass2(2)=sum(peptides.mad(:)>40);
summary.tot_pass2(3)=sum(peptides.mad_conf(:)>25);
summary.tot_pass2(4)=sum(peptides.iqrNorm(:)>20);
summary.tot_pass2(5)=sum(peptides.iqrNorm_conf(:)>10)

writetable(summary,'WangKuster_outlierEval_threshSummary.csv');
save('WangKuster_outlierEval.mat');
