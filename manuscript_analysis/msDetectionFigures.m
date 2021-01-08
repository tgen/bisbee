work_dir='/labs/halperin/AlternativeSplicing/WangKuster2019/'
ms_detect='WangKuster.top.sepNovel.MSdetectionAll.v2.csv'
anno_file='bisbee/WangKuster.comb.bisbeeDiff.threshnan.anno.csv'

cd(work_dir)

%%% read ms peptide detection table
predPeptides=readtable(ms_detect);
predPeptides.event_jid=extractAfter(predPeptides.effectId,'_');

%%% read bisbee annotations and join to ms peptide detection
bbDiff=readtable(anno_file);
predPeptides=join(predPeptides,bbDiff,'Keys',{'event_jid'},'RightVariables',2:19);

%%% define event groups for novel
novel_idx=strcmp(predPeptides.aa_change_type,"Novel");
predPeptides.event_group=predPeptides.event_cat;
predPeptides.event_group(novel_idx & endsWith(predPeptides.event_cat,'ExonInc') & strcmp(predPeptides.refIsoform,'iso1'))={'ExonSkip'};
predPeptides.event_group(novel_idx & endsWith(predPeptides.event_cat,'ExonInc') & strcmp(predPeptides.refIsoform,'iso2'))={'ExonInc'};
predPeptides.event_group(novel_idx & endsWith(predPeptides.event_cat,'IntronRet') & strcmp(predPeptides.refIsoform,'iso1'))={'IntronExc'};
predPeptides.event_group(novel_idx & endsWith(predPeptides.event_cat,'IntronRet') & strcmp(predPeptides.refIsoform,'iso2'))={'IntronRet'};
[tbl,~,~,labels]=crosstab(predPeptides.effect_cat(novel_idx),predPeptides.event_group(novel_idx));
totalNovel=array2table(tbl,'RowNames',labels(1:size(tbl,1),1),'VariableNames',labels(1:size(tbl,2),2));

%%% define event groups for non-novel
canon_idx=strcmp(predPeptides.aa_change_type,"Canonical") | strcmp(predPeptides.aa_change_type,"Other");
predPeptides.event_group(canon_idx & endsWith(predPeptides.event_cat,'ExonInc'))={'ExonIncSkip'};
predPeptides.event_group(canon_idx & endsWith(predPeptides.event_cat,'IntronRet'))={'IntronRetExc'};
predPeptides.effect_group=predPeptides.effect_cat;
predPeptides.effect_group(canon_idx & (strcmp(predPeptides.effect_cat,'Insertion') | strcmp(predPeptides.effect_cat,'Deletion')))={'InDel'};
[tbl,~,~,labels]=crosstab(predPeptides.effect_group(canon_idx),predPeptides.event_group(canon_idx));
totalNotNovel=array2table(tbl,'RowNames',labels(1:size(tbl,1),1),'VariableNames',labels(1:size(tbl,2),2));

%%% find ms isoform detection columns
iso1Ms_idx=startsWith(predPeptides.Properties.VariableNames,'iso1_maxAb');
iso2Ms_idx=startsWith(predPeptides.Properties.VariableNames,'iso2_maxAb');

%%% find rows were both isoforms are detected by MS
detect_both_idx=max(predPeptides{:,iso1Ms_idx},[],2)>0 & max(predPeptides{:,iso2Ms_idx},[],2)>0;

%%% make summary table of events where both isoforms are detected for novel
[tbl,~,~,labels]=crosstab(predPeptides.effect_cat(novel_idx & detect_both_idx),predPeptides.event_group(novel_idx & detect_both_idx));
detectNovel=array2table(tbl,'RowNames',labels(1:size(tbl,1),1),'VariableNames',labels(1:size(tbl,2),2));

%%% count number of events where either isoform is detected
detect_idx=max(predPeptides{:,iso1Ms_idx | iso2Ms_idx},[],2)>0;
sum(detect_idx)

%%% make summary table of events where both isoforms are detected for
%%% non-novel events
[tbl,~,~,labels]=crosstab(predPeptides.effect_group(canon_idx & detect_both_idx),predPeptides.event_group(detect_both_idx & canon_idx));
detectNotNovel=array2table(tbl,'RowNames',labels(1:size(tbl,1),1),'VariableNames',labels(1:size(tbl,2),2));

%%% write output files
writetable(totalNovel,'novelPeptByEventAndEffectCat.csv','WriteRowNames',1);
writetable(totalNotNovel,'canonPeptByEventAndEffectCat.csv','WriteRowNames',1);
writetable(detectNovel,'detectNovelPeptByEventAndEffectCat.csv','WriteRowNames',1);
writetable(detectNotNovel,'detectCanonPeptByEventAndEffectCat.csv','WriteRowNames',1);
writetable(predPeptides,'WangKuster.top.MSDetection.bisbeeDiff.csv');

%%% generate old figure 2
effects=sort(detectNovel.Properties.RowNames);
event_cats_novel=sort(detectNovel.Properties.VariableNames);
x=[1:length(effects)]'*ones(1,length(event_cats_novel));
y=ones(length(effects),1)*[1:length(event_cats_novel)];
subplot(1,2,1)
scatter(x(:),y(:),20*log(reshape(totalNovel{effects,event_cats_novel},[],1)+2),log(reshape(detectNovel{effects,event_cats_novel},[],1)./reshape(totalNovel{effects,event_cats_novel},[],1)+1E-5),'filled');
text(x(:),y(:),strcat('$\frac{',num2str(reshape(detectNovel{effects,event_cats_novel},[],1)),'}{',num2str(reshape(totalNovel{effects,event_cats_novel},[],1),'%-d'),'}$'),'HorizontalAlignment','center','FontSize',8,'Interpreter','latex')
colormap('cool');
c=colorbar('Ticks',log(10.^[-4:0]+1E-5),'TickLabels',num2str(10.^[-4:0]','%.0E'));
c.Label.String='% detected';
c.AxisLocation='in'
c.FontSize=6;
caxis(log([1E-5 1]))
axis([0 length(effects)+2 0 length(event_cats_novel)+1])
set(gca,'XTick',1:length(effects),'XTickLabel',effects,'XTickLabelRotation',90);
set(gca,'YTick',1:length(event_cats_novel),'YTickLabel',event_cats_novel,'FontSize',6);
title('Novel Events');
axis('square');

effects=sort(detectNotNovel.Properties.RowNames);
event_cats=sort(detectNotNovel.Properties.VariableNames);
x=[1:length(effects)]'*ones(1,length(event_cats));
y=ones(length(effects),1)*[1:length(event_cats)];
subplot(1,2,2)
scatter(x(:),y(:),20*log(reshape(totalNotNovel{effects,event_cats},[],1)+2),log(reshape(detectNotNovel{effects,event_cats},[],1)./reshape(totalNotNovel{effects,event_cats},[],1)+1E-5),'filled');
text(x(:),y(:),strcat('$\frac{',num2str(reshape(detectNotNovel{effects,event_cats},[],1)),'}{',num2str(reshape(totalNotNovel{effects,event_cats},[],1),'%-d'),'}$'),'HorizontalAlignment','center','FontSize',8,'Interpreter','latex')
colormap('cool');
c=colorbar('Ticks',log(10.^[-4:0]+1E-5),'TickLabels',num2str(10.^[-4:0]','%.0E'))
c.Label.String='% detected';
c.AxisLocation='in'
c.FontSize=6;
caxis(log([1E-5 1]))
caxis(log([1E-5 1]))
axis([0 length(effects)+2 0 length(event_cats)+1])
set(gca,'XTick',1:length(effects),'XTickLabel',effects,'XTickLabelRotation',90);
set(gca,'YTick',1:length(event_cats),'YTickLabel',event_cats,'FontSize',6);
title('Known Isoform Events');
axis('square');

set(gcf,'papersize',[6.69 3.3])
print('Fig2_eventEffectCounts.pdf','-dpdf','-fillpage','-r300');
print('Fig2_eventEffectCounts.eps','-dsvg','-r300');
print('Fig2_eventEffectCounts.png','-dpng','-r300');
close(gcf);

%%% seperate event categories by ref and alt perspective
predPeptides.event_cat_alt=predPeptides.event_cat;
predPeptides.event_cat_alt(endsWith(predPeptides.event_cat,'ExonInc') & strcmp(predPeptides.refIsoform,'iso1'))={'ExonSkip'};
predPeptides.event_cat_alt(endsWith(predPeptides.event_cat,'ExonInc') & strcmp(predPeptides.refIsoform,'iso2'))={'ExonInc'};
predPeptides.event_cat_alt(endsWith(predPeptides.event_cat,'IntronRet') & strcmp(predPeptides.refIsoform,'iso1'))={'IntronExc'};
predPeptides.event_cat_alt(endsWith(predPeptides.event_cat,'IntronRet') & strcmp(predPeptides.refIsoform,'iso2'))={'IntronRet'};
predPeptides.event_cat_ref=predPeptides.event_cat;
predPeptides.event_cat_ref(endsWith(predPeptides.event_cat,'ExonInc') & strcmp(predPeptides.refIsoform,'iso2'))={'ExonSkip'};
predPeptides.event_cat_ref(endsWith(predPeptides.event_cat,'ExonInc') & strcmp(predPeptides.refIsoform,'iso1'))={'ExonInc'};
predPeptides.event_cat_ref(endsWith(predPeptides.event_cat,'IntronRet') & strcmp(predPeptides.refIsoform,'iso2'))={'IntronExc'};
predPeptides.event_cat_ref(endsWith(predPeptides.event_cat,'IntronRet') & strcmp(predPeptides.refIsoform,'iso1'))={'IntronRet'};

%%% get indexes of event categories
if_idx=strcmp(predPeptides.effect_cat,'Deletion') | strcmp(predPeptides.effect_cat,'Insertion') | strcmp(predPeptides.effect_cat,'Substitution');
fs_idx=strcmp(predPeptides.effect_cat,'FrameDisruption');
wt_idx=strcmp(predPeptides.aa_change_type,'Canonical') | strcmp(predPeptides.aa_change_type,'Other');
nv_idx=strcmp(predPeptides.aa_change_type,'Novel');

%%% determine where ref and alt peptides are detected
iso1Ms_idx=startsWith(predPeptides.Properties.VariableNames,'iso1_maxAb');
iso2Ms_idx=startsWith(predPeptides.Properties.VariableNames,'iso2_maxAb');
ref_detect=(max(predPeptides{:,iso1Ms_idx},[],2)>0 & strcmp(predPeptides.refIsoform,'iso1')) | (max(predPeptides{:,iso2Ms_idx},[],2)>0 & strcmp(predPeptides.refIsoform,'iso2'));
alt_detect=(max(predPeptides{:,iso1Ms_idx},[],2)>0 & strcmp(predPeptides.refIsoform,'iso2')) | (max(predPeptides{:,iso2Ms_idx},[],2)>0 & strcmp(predPeptides.refIsoform,'iso1'));
groups=sort([strcat(event_cats_novel','true'); strcat(event_cats_novel','false')]);

predPeptides.detectCat=repmat({'NonCoding'},height(predPeptides),1);
predPeptides.detectCat(novel_idx & ref_detect)={'Novel-RefOnly'};
predPeptides.detectCat(novel_idx & alt_detect)={'Novel-AltOnly'};
predPeptides.detectCat(novel_idx & detect_both_idx)={'Novel-bothIso'};
predPeptides.detectCat(novel_idx & ~detect_idx)={'Novel-none'};
predPeptides.detectCat(canon_idx & ~detect_idx)={'Canon-none'};
predPeptides.detectCat(canon_idx & detect_idx)={'Canon-OneIso'};
predPeptides.detectCat(canon_idx & detect_both_idx)={'Canon-bothIso'};

%%% generate new Figure 2
detectCounts=cell2table(tabulate(predPeptides.detectCat));
detectCounts=detectCounts([5 3 4 7 2 6 1 8],:); 
subplot(3,2,1)
labels=cell(3,2);
for i=2:4
    labels(i-1,:)={strrep(detectCounts.Var1{i},'Canon-',''),num2str(detectCounts.Var2(i))};
end
p=pie(detectCounts.Var2(2:4),[0 0 1]);
for i=1:length(p)/2
    p(i*2).String=labels(i,:);
    
end
text(-1,2,'Canonical Event Detection','FontSize',12,'HorizontalAlignment','left','FontWeight','bold');
subplot(3,2,2)
labels=cell(4,2);
for i=5:8
    labels(i-4,:)={strrep(detectCounts.Var1{i},'Novel-',''),num2str(detectCounts.Var2(i))};
end
p=pie(detectCounts.Var2(5:8),[0 0 0 1]);
for i=1:length(p)/2
    p(i*2).String=labels(i,:);
end
p(3).CData=3*ones(size(p(3).CData));
p(5).CData=2*ones(size(p(5).CData));
text(-1,2,'Novel Event Detection','FontSize',12,'HorizontalAlignment','left','FontWeight','bold');

colors=jet(6);
nnIdx=[2 4 5];
nIdx=[1 3 4 5 6];

subplot(3,2,3)
[~,ord]=sort(totalNotNovel.Properties.VariableNames);
b=bar(totalNotNovel{:,ord}','stacked');
for i=1:length(b)
    b(i).FaceColor=colors(nnIdx(i),:);
end
set(gca,'XTickLabel',totalNotNovel.Properties.VariableNames(ord),'XTickLabelRotation',90);
legend(totalNotNovel.Properties.RowNames,'Location','northoutside');
title('All Known Isoform Events')

subplot(3,2,4)
[~,ord]=sort(totalNovel.Properties.VariableNames);
b=bar(totalNovel{:,ord}','stacked');
for i=1:length(b)
    b(i).FaceColor=colors(nIdx(i),:);
end
set(gca,'XTickLabel',totalNovel.Properties.VariableNames(ord),'XTickLabelRotation',90);
legend(totalNovel.Properties.RowNames,'Location','northoutside');
title('All Novel Isoform Events')

subplot(3,2,5)
[~,ord]=sort(detectNotNovel.Properties.VariableNames);
b=bar(detectNotNovel{:,ord}','stacked');
for i=1:length(b)
    b(i).FaceColor=colors(nnIdx(i),:);
end
set(gca,'XTickLabel',detectNotNovel.Properties.VariableNames(ord),'XTickLabelRotation',90);
legend(detectNotNovel.Properties.RowNames,'Location','northoutside');
title('MS Detected Known Isoform Events')

subplot(3,2,6)
[~,ord]=sort(detectNovel.Properties.VariableNames);
b=bar(detectNovel{:,ord}','stacked');
for i=1:length(b)
    b(i).FaceColor=colors(nIdx(i),:);
end
set(gca,'XTickLabel',detectNovel.Properties.VariableNames(ord),'XTickLabelRotation',90);
legend(detectNovel.Properties.RowNames,'Location','northoutside');
title('MS Detected Isoform Events')

print('Fig2_eventDetection.pdf','-dpdf','-fillpage','-r300');
print('Fig2_eventDetection.svg','-dsvg','-r300');
print('Fig2_eventDetection.png','-dpng','-r300');
close(gcf);

%%% plot Supplemental Figure 2
subplot(2,2,1);
boxplot(log([predPeptides.insSeqLen(if_idx & wt_idx); predPeptides.delSeqLen(if_idx)]+1),strcat([predPeptides.event_cat_alt(if_idx & wt_idx); predPeptides.event_cat_ref(if_idx)],string([alt_detect(if_idx & wt_idx); ref_detect(if_idx)])),'GroupOrder',groups,'ColorGroup',[alt_detect(if_idx & wt_idx); ref_detect(if_idx)],'PlotStyle','compact')
set(gca,'YTick',log(2.^[0:2:14]+1),'YTickLabel',2.^[0:2:14])
set(gca,'XTick',1.5:2:0.5+2*length(event_cats_novel),'XTickLabel',event_cats_novel,'XTickLabelRotation',90);
ylim(log([0 2.^15]+1));
ylabel('alt aa');
title('InFrame Known');

subplot(2,2,2);
boxplot(log(predPeptides.insSeqLen(if_idx & nv_idx)+1),strcat(predPeptides.event_cat_alt(if_idx & nv_idx),string(alt_detect(if_idx & nv_idx))),'GroupOrder',groups,'ColorGroup',alt_detect(if_idx & nv_idx),'PlotStyle','compact')
set(gca,'YTick',log(2.^[0:2:14]+1),'YTickLabel',2.^[0:2:14])
set(gca,'XTick',1.5:2:0.5+2*length(event_cats_novel),'XTickLabel',event_cats_novel,'XTickLabelRotation',90);
ylim(log([0 2.^15]+1));
ylabel('alt aa');
title('InFrame Novel');

subplot(2,2,3);
boxplot(log([predPeptides.insSeqLen(fs_idx & wt_idx); predPeptides.delSeqLen(fs_idx)]+1),strcat([predPeptides.event_cat_alt(fs_idx & wt_idx); predPeptides.event_cat_ref(fs_idx)],string([alt_detect(fs_idx & wt_idx); ref_detect(fs_idx)])),'GroupOrder',groups,'ColorGroup',[alt_detect(fs_idx & wt_idx); ref_detect(fs_idx)],'PlotStyle','compact')
set(gca,'YTick',log(2.^[0:2:14]+1),'YTickLabel',2.^[0:2:14])
set(gca,'XTick',1.5:2:0.5+2*length(event_cats_novel),'XTickLabel',event_cats_novel,'XTickLabelRotation',90);
ylim(log([0 2.^15]+1));
ylabel('alt aa');
title('FrameDisrupt Known');

subplot(2,2,4);
boxplot(log(predPeptides.insSeqLen(fs_idx & nv_idx)+1),strcat(predPeptides.event_cat_alt(fs_idx & nv_idx),string(alt_detect(fs_idx & nv_idx))),'GroupOrder',groups,'ColorGroup',alt_detect(fs_idx & nv_idx),'PlotStyle','compact')
set(gca,'YTick',log(2.^[0:2:14]+1),'YTickLabel',2.^[0:2:14])
set(gca,'XTick',1.5:2:0.5+2*length(event_cats_novel),'XTickLabel',event_cats_novel,'XTickLabelRotation',90);
ylim(log([0 2.^15]+1));
ylabel('alt aa');
title('FrameDisrupt Novel');

set(gcf,'papersize',[6.69 6.69])
print('FigS2_altSeqLengths.pdf','-dpdf','-fillpage','-r300');
print('FigS2_altSeqLengths.svg','-dsvg','-r300');
print('FigS2_altSeqLengths.png','-dpng','-r300');
close(gcf);

%%% create table by ms detected peptide
msPept={};
for i=1:height(predPeptides)
    if ~isempty(predPeptides.iso1_seq{i})
        pept=strsplit(predPeptides.iso1_seq{i},',');
        if strcmp(predPeptides.refIsoform{i},'iso1')
            effect_cat=replace(predPeptides.effect_cat{i},{'Insertion','Deletion','FrameDisruption','Truncation'},{'Deletion','Insertion','NonFrameDisruption','NonTruncation'});
            info={predPeptides.gene{i}, predPeptides.event_cat_ref{i},effect_cat,replace(predPeptides.aa_change_type{i},'Novel','NotNovel')};
        else
            info={predPeptides.gene{i}, predPeptides.event_cat_alt{i}, predPeptides.effect_cat{i},predPeptides.aa_change_type{i}};
        end
        msPept=[msPept;[pept' repmat(info,length(pept),1)]];
    end
    if ~isempty(predPeptides.iso2_seq{i})
        pept=strsplit(predPeptides.iso2_seq{i},',');
        if strcmp(predPeptides.refIsoform{i},'iso2')
            effect_cat=replace(predPeptides.effect_cat{i},{'Insertion','Deletion','FrameDisruption','Truncation'},{'Deletion','Insertion','NonFrameDisruption','NonTruncation'});
            info={predPeptides.gene{i}, predPeptides.event_cat_ref{i},effect_cat,replace(predPeptides.aa_change_type{i},'Novel','NotNovel')};
        else
            info={predPeptides.gene{i}, predPeptides.event_cat_alt{i}, predPeptides.effect_cat{i},predPeptides.aa_change_type{i}};
        end
        msPept=[msPept;[pept' repmat(info,length(pept),1)]];
    end    
    if ~isempty(predPeptides.both_seq{i})
        pept=strsplit(predPeptides.both_seq{i},',');
        info={predPeptides.gene{i},'None','None','None'};
        msPept=[msPept;[pept' repmat(info,length(pept),1)]];
    end
    i
end
msPept=cell2table(msPept,'VariableNames',{'peptide','gene','event_cat','effect_cat','aa_change_type'});
writetable(msPept,'WangKuster.MSpeptides.csv');

%%% count number of genes that each peptide maps to
[G,u_pept]=findgroups(msPept.peptide);
gene_count=splitapply(@(x) length(unique(x)),msPept.gene,G);
msPeptUnique=table(u_pept,gene_count,'VariableNames',{'peptide','gene_count'});
msPept.aa_change_type=categorical(msPept.aa_change_type,{'Canonical','NotNovel','Other','Novel','None'},'Ordinal',1);
msPept.effect_cat=categorical(msPept.effect_cat,{'Deletion','Insertion','Substitution','Truncation','FrameDisruption','NonTruncation','NonFrameDisruption','None'},'Ordinal',1);

msPept=sortrows(msPept,{'peptide','aa_change_type','effect_cat'});
[lia,locb]=ismember(msPeptUnique.peptide,msPept.peptide);
msPeptUnique.event_cat=msPept.event_cat(locb);
msPeptUnique.effect_cat=msPept.effect_cat(locb);
msPeptUnique.aa_change_type=msPept.aa_change_type(locb);

[tbl,~,~,labels]=crosstab(msPeptUnique.event_cat,msPeptUnique.effect_cat,msPeptUnique.aa_change_type);
msPeptSummary={};
for i=1:size(tbl,3)
    msPeptSummary{i}=array2table(tbl(:,:,i),'VariableNames',labels(:,2));
    msPeptSummary{i}.effect_cat=labels(:,1);
    msPeptSummary{i}.type=repmat(labels(i,3),size(labels,1),1);
end
msPeptSummary=vertcat(msPeptSummary{:});

writetable(msPeptSummary,'msPeptSummary.csv');
writetable(msPept,'msPeptCategories.csv');
writetable(msPeptUnique,'msPeptUnique.csv');

%%% get peptide mapping counts
peptCounts=readtable('WangKusterSpladder2.sepNovel.prot_peptDetection.counts.csv');
canon_pept_idx=msPeptUnique.aa_change_type=='Canonical' | msPeptUnique.aa_change_type=='Other';
peptCounts.canonicalIso=sum(canon_pept_idx);
peptCounts.noIso=sum(msPeptUnique.aa_change_type=='None');
peptCounts.novelAltIso=sum(msPeptUnique.aa_change_type=='Novel');
peptCounts.novelRefIso=sum(msPeptUnique.aa_change_type=='NotNovel');
peptCounts=peptCounts(:,[1 4:7]);


%%% plot figure supplemental figure 1
subplot(1,3,1)
p=pie(peptCounts{:,:});
for i=1:length(p)/2
    p(i*2).String={peptCounts.Properties.VariableNames{i},num2str(peptCounts{:,i})};
end
colormap('jet');
p(3).CData=4*ones(size(p(3).CData));
p(5).CData=2*ones(size(p(5).CData));
p(7).CData=5*ones(size(p(7).CData));
p(9).CData=3*ones(size(p(9).CData));
title('MS detected peptide mapping');
[tbl,~,~,labels]=crosstab(msPeptUnique.event_cat(canon_pept_idx),msPeptUnique.effect_cat(canon_pept_idx));
canonPeptCounts=array2table(tbl,'RowNames',labels(1:size(tbl,1),1),'VariableNames',labels(:,2));
subplot(1,3,2)
canonPeptCounts=canonPeptCounts(:,sum(canonPeptCounts{:,:})>0);
[~,ord]=sort(canonPeptCounts.Properties.RowNames);
b=bar(canonPeptCounts{ord,:},'stacked');
for i=1:length(b)
    b(i).FaceColor=colors(nIdx(i),:);
end
set(gca,'XTickLabel',canonPeptCounts.Properties.RowNames(ord),'XTickLabelRotation',90);
legend(canonPeptCounts.Properties.VariableNames,'Location','northoutside');
title('MS Pept mapping to canon iso')
[tbl,~,~,labels]=crosstab(msPeptUnique.event_cat(msPeptUnique.aa_change_type=='Novel'),msPeptUnique.effect_cat(msPeptUnique.aa_change_type=='Novel'));
novelPeptCounts=array2table(tbl,'RowNames',labels(1:size(tbl,1),1),'VariableNames',labels(:,2));
subplot(1,3,3)
novelPeptCounts=novelPeptCounts(:,sum(novelPeptCounts{:,:})>0);
[~,ord]=sort(novelPeptCounts.Properties.RowNames);
b=bar(novelPeptCounts{ord,:},'stacked');
for i=1:length(b)
    b(i).FaceColor=colors(nIdx(i),:);
end
set(gca,'XTickLabel',novelPeptCounts.Properties.RowNames(ord),'XTickLabelRotation',90);
legend(novelPeptCounts.Properties.VariableNames,'Location','northoutside');
title('MS Pept mapping to novel alt iso')
orient('landscape');
print('FigS1_MsPeptMappingBarPlots.pdf','-dpdf','-r300','-fillpage');
