work_dir='/labs/halperin/AlternativeSplicing/WangKuster2019/';
ms_res_dir='MS_search_results_sepNovel';
top_fasta='bisbee/WangKuster.top.fasta'
top_peptides='bisbee/prot/WangKuster.top.peptides.csv'

cd(work_dir);

%%% read in ms results
msFiles=struct2table(dir([ms_res_dir '/*.txt']));
msPept={};
for i=1:height(msFiles)
    file=[msFiles.folder{i} '/' msFiles.name{i}];
    opt=detectImportOptions(file);
    aIdx=find(startsWith(opt.VariableNames,'Abundance_F'));
    opt.SelectedVariableNames={'AnnotatedSequence',opt.VariableNames{aIdx}};
    opt=setvartype(opt,aIdx,'double');
    msPept{i}=readtable(file,opt);
    msPept{i}.tissue=cellstr(repmat(extractAfter(opt.VariableNames{aIdx},'Sample_'),height(msPept{i}),1));
    msPept{i}.Properties.VariableNames{2}='Abundance';
    msPept{i}.novel=cellstr(repmat(extractBetween(file,'Bisbee','Ms'),height(msPept{i}),1));
    i
end
msPept=vertcat(msPept{:});
msPept.Abundance(isnan(msPept.Abundance))=0;
msPept.Sequence=extractBetween(msPept.AnnotatedSequence,'.','.');
wtPeptMs=unique(msPept.Sequence(strcmp(msPept.novel,'WT')));
lia=ismember(msPept.Sequence,wtPeptMs);
novelPeptMs=unique(msPept.Sequence(~lia));

msPept=unstack(msPept,'Abundance','tissue');
msPept.novel=ismember(msPept.Sequence,novelPeptMs);

%%% read in bisbee predicted protien sequences
topProt=struct2table(fastaread(top_fasta));

%%% find which sequences contain exact matches to ms detected peptides
peptDetect=sparse(height(topProt),height(msPept));
parfor i=1:height(msPept)
    peptDetect(:,i)=contains(topProt.Sequence,msPept.Sequence(i));
    %i
end
save('peptDetect.mat','peptDetect');
[r,c,val]=find(peptDetect);
peptDetectTable=table(topProt.Header(r),msPept.Sequence(c),val);
peptDetectTable.Properties.VariableNames={'protHeader','msPeptSeq','count'};
writetable(peptDetectTable,'WangKusterSpladder2.sepNovel.prot_peptDetection.csv');

%%% remove peptides that map to multiple genes and count novel peptide mapping
peptDetectTable.novel=endsWith(peptDetectTable.protHeader,'NS=True');
peptDetectTable.gene=extractBetween(peptDetectTable.protHeader,'GN=',' ');
[G,u_pept]=findgroups(peptDetectTable.msPeptSeq);
gene_count=splitapply(@(x) length(unique(x)),peptDetectTable.gene,G);
novel_only=splitapply(@(x) min(x),peptDetectTable.novel,G);
gene_specific_pept=u_pept(gene_count==1);
counts.multiGene=sum(gene_count>1);
counts.notNovel=sum(novel_only==0 & gene_count==1);
counts.novel=sum(novel_only==1 & gene_count==1);
writetable(struct2table(counts),'WangKusterSpladder2.sepNovel.prot_peptDetection.counts.csv');
lia=ismember(msPept.Sequence,gene_specific_pept);
msPeptSelect=msPept(lia,:);
[lia,locb]=ismember(peptDetectTable.msPeptSeq,msPeptSelect.Sequence);
[prot_header,ia,ic]=unique(peptDetectTable.protHeader);
peptDetect=sparse(ic(lia),locb(lia),ones(sum(lia),1));

%%% read in bisbee splice event data
predPeptides=readtable(top_peptides,'Delimiter',',');
predPeptides.iso1_count=nan(height(predPeptides),1);
predPeptides.iso1_seq=cell(height(predPeptides),1);
predPeptides.iso1_maxAb=nan(height(predPeptides),height(msFiles));
predPeptides.iso2_count=nan(height(predPeptides),1);
predPeptides.iso2_seq=cell(height(predPeptides),1);
predPeptides.iso2_maxAb=nan(height(predPeptides),height(msFiles));
predPeptides.both_count=nan(height(predPeptides),1);
predPeptides.both_seq=cell(height(predPeptides),1);
predPeptides.both_maxAb=nan(height(predPeptides),height(msFiles));

%%% map ms detected peptides to splice events
for i=1:height(predPeptides)
    iso1_idx=find(startsWith(prot_header,[predPeptides.iso1_header{i} ' ']));
    iso2_idx=find(startsWith(prot_header,[predPeptides.iso2_header{i} ' ']));
    if isempty(iso1_idx) | isempty(iso2_idx)
        continue
    end
    iso1_pept=peptDetect(iso1_idx,:) & ~peptDetect(iso2_idx,:);
    predPeptides.iso1_count(i)=sum(iso1_pept);
    predPeptides.iso1_seq{i}=strjoin(msPeptSelect.Sequence(iso1_pept),',');
    if predPeptides.iso1_count(i)>0
        for j=1:height(msFiles)/2
            predPeptides.iso1_maxAb(i,j)=max(msPeptSelect{iso1_pept,j+3});
        end
    end
    iso2_pept=~peptDetect(iso1_idx,:) & peptDetect(iso2_idx,:);
    predPeptides.iso2_count(i)=sum(iso2_pept);
    predPeptides.iso2_seq{i}=strjoin(msPeptSelect.Sequence(iso2_pept),',');
    if predPeptides.iso2_count(i)>0
       for j=1:height(msFiles)/2
            predPeptides.iso2_maxAb(i,j)=max(msPeptSelect{iso2_pept,j+3});
        end
    end
    both_pept=peptDetect(iso1_idx,:) & peptDetect(iso2_idx,:);
    predPeptides.both_count(i)=sum(both_pept);
    predPeptides.both_seq{i}=strjoin(msPeptSelect.Sequence(both_pept),',');
    if predPeptides.both_count(i)>0
       for j=1:height(msFiles)/2
            predPeptides.both_maxAb(i,j)=max(msPeptSelect{both_pept,j+3});
        end
    end
    i./height(predPeptides)
end

%%% write output file
writetable(predPeptides,'WangKuster.top.sepNovel.MSdetectionAll.v2.csv');
