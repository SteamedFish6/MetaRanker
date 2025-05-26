#MetaRanker v0.1
#Author: Zhenyu Guo

import os
import time
import subprocess
import gzip
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO

## Requirements for shell:
# - blastn
# - cd-hit
# - bwa
# - minimap2 #for long reads
# - samtools
# - pigz #accelerate parsing gzipped fastq file

def _getRankerArgs():
    parser = argparse.ArgumentParser(description='''MetaRanker v0.1''')
    parser.add_argument('-c', "--contigs", type=str, help="[File] name of input contigs file in fasta format", required=True)
    parser.add_argument('-r', "--reads1", type=str, help="[File] name of input reads1 file in fastq(.gz) format", required=True)
    parser.add_argument('-R', "--reads2", type=str, help="[File] name of input reads2 file in fastq(.gz) format, required for pair-ended reads", default='')
    parser.add_argument('-o', "--output", type=str, help="[Dir] directory of output results", required=True)
    parser.add_argument("--nanopore", action='store_true', help="[Flag] process nanopore reads")
    parser.add_argument("--pacbio", action='store_true', help="[Flag] process pacbio reads")
    
    parser.add_argument("--minnum", type=int, help="[Int] stop running pipeline if input contigs fewer than setting threshold (default: 2000)", default=2000)
    parser.add_argument("--minlen", type=int, help="[Int] discard input contigs shorter than setting threshold (default: 500)", default=500)
    parser.add_argument("--no_rename_contigs", action='store_true', help="[Flag] do not rename contigs names with contigs file name")
    
    parser.add_argument("--blast_evalue", type=float, help="[Float] evalue of blastn (default: 0.0001)", default=0.0001)
    parser.add_argument("--blast_identity", type=float, help="[Float] minimum identity of blastn (default: 0.85)", default=0.85)
    parser.add_argument("--blast_cover_len", type=int, help="[Int] minimum cover length of blastn (default: 75)", default=75)
    parser.add_argument("--cdhit_identity", type=float, help="[Float] minimum identity of cd-hit (default: 0.85)", default=0.85)
    
    parser.add_argument("--weight", nargs=3, type=float, help="[Float Float Float] weight of ARG, MGE, VF database (default: 50000 10000 20000)", default=[5e4, 1e4, 2e4])
    
    parser.add_argument('-t', "--threads", type=int, help="[Int] number of threads for BLAST, CD-HIT, BWA, minimap2 and samtools (default: 16)", default=16)
    parser.add_argument("--force", action='store_true', help="[Flag] force to cover existing output files of BLAST, CD-HIT, BWA, minimap2 and samtools")
    
    parser.add_argument('--version', action='version', version="MetaRanker_v0.1")
    
    return parser.parse_args()


def PrepareBlastDB(dbpath: str, dblist: str, is_cover=False) -> None:
    extend_name = "nhr"
    dbtype = "nucl"
    printTime()
    print("Preparing BLAST Database...", end=' ')
    for dbname in dblist:
        dbfile = os.path.join(dbpath, "{}.fasta".format(dbname))
        checkCallCMD("makeblastdb -dbtype {} -parse_seqids -in {}".format(dbtype, dbfile), "{}.{}".format(dbfile, extend_name), is_cover_old=is_cover, print_skipped=False)
    print("BLAST Database is Ready.")

class Sample:
    def __init__(self, contigs_fname: str, dbpath: str, dblist: list, reads_fname1: str, reads_fname2: str='', is_long_reads=False, reads_type="illumina"):
        self.contigs_fname = contigs_fname
        self.dbpath = dbpath
        self.dblist = dblist
        self.reads_fname1 = reads_fname1
        self.reads_fname2 = reads_fname2
        self.is_long_reads = is_long_reads
        self.reads_type = reads_type
        
        self.name_tag = os.path.basename(contigs_fname).rsplit('.', 1)[0]
        self.blastmethod = "blastn"
        self.M8_fdict = {} #key: db, value: M8_fname
        self.REseq_fdict = {} #key: db, value: RiskElements.fa
        self.depth_fdict = {} #key: db, value: depth.tsv
        self.rank_fname = None
        # self.extract_fname = None
        
        self.ncontig = 1
        self.total_reads_num = 1
        self.total_reads_bases = 1
        self.REbases_ddict = {} #key: db, value: dict REseqs base num
    
    def SeqPreprocess(self, outpath: str, minlen=0, is_rename_contigs=False) -> None:
        printTime()
        print("Checking input fastq reads and contigs fasta file...")
        total_reads_num1, total_reads_bases1 = countFastqReads_Bases_pigzAWK(self.reads_fname1, is_gzipped=checkGZipped(self.reads_fname1))
        if self.reads_fname2 != '':
            total_reads_num2, total_reads_bases2 = countFastqReads_Bases_pigzAWK(self.reads_fname2, is_gzipped=checkGZipped(self.reads_fname2))
        else:
            total_reads_num2, total_reads_bases2 = 0, 0
        self.total_reads_num = total_reads_num1 + total_reads_num2
        self.total_reads_bases = total_reads_bases1 + total_reads_bases2
        print("{} reads, {} bases found in {} {}".format(self.total_reads_num, self.total_reads_bases, self.reads_fname1, self.reads_fname2))
        
        record_dict = readFastaFile(self.contigs_fname, to_which='dict')
        
        if minlen <= 0:
            ncontig = len(record_dict)
            print("{} sequences found in {}".format(ncontig, self.contigs_fname))
            if is_rename_contigs:
                outname = os.path.join(outpath, "rename_{}".format(minlen, os.path.basename(self.contigs_fname)))
                record_dict = renameFastaSeqs(record_dict, self.name_tag)
                writeFastaFromDict(record_dict, outname)
                self.contigs_fname = outname
        else:
            new_record_dict = {}
            for seqname in record_dict:
                seq = record_dict[seqname]
                if len(seq) >= minlen:
                    new_record_dict[seqname] = seq
            if is_rename_contigs:
                new_record_dict = renameFastaSeqs(new_record_dict, self.name_tag)
                outname = os.path.join(outpath, "rename_filter{}_{}".format(minlen, os.path.basename(self.contigs_fname)))
            else:
                outname = os.path.join(outpath, "filter{}_{}".format(minlen, os.path.basename(self.contigs_fname)))
            # SeqIO.write(new_record_list, outname, "fasta")
            writeFastaFromDict(new_record_dict, outname)
            print("{} was filtered by >= {} nt/aa, new fasta is {}".format(self.contigs_fname, minlen, outname))
            self.contigs_fname = outname
            ncontig = len(new_record_dict)
            print("{} sequences found in {}".format(ncontig, self.contigs_fname))
        self.ncontig = ncontig
    
    def Blast(self, outpath: str, is_cover_old=False) -> None:
        for dbname in self.dblist:
            printTime()
            db_fname = os.path.join(self.dbpath, "{}.fasta".format(dbname))
            print("Processing {}: {}, in db {}".format(self.blastmethod, self.contigs_fname, dbname))
            tmp_outname = os.path.join(outpath, "tmp_{}_{}_{}.tsv".format(self.blastmethod, self.name_tag, dbname))
            outname = os.path.join(outpath, "out_{}_{}_{}.tsv".format(self.blastmethod, self.name_tag, dbname))
            is_rename = checkCallCMD("{} -task {} -query {} -db {} -out {} -outfmt \"6\" -evalue {} -perc_identity {} -num_threads {}".format(self.blastmethod, self.blastmethod, self.contigs_fname, db_fname, tmp_outname, BLAST_EVALUE, PERC_IDENTITY, NUM_THREAD), outname, is_cover_old=is_cover_old)
            if is_rename:
                if is_cover_old and os.path.exists(outname):
                    os.remove(outname)
                os.rename(tmp_outname, outname)
            self.M8_fdict[dbname] = outname
    
    def M8Preprocess(self, outpath: str, dbmask_dict: dict, is_pos_max=False, is_group_max=False, criteria: str=None, sortby='score', noncross_len=75) -> None:
        printTime()
        print("Handling BLAST result...")
        column_index = ['query_ID', 'target_ID', 'identity', 'cover_len', 'not_match', 'gap', 'start_query', 'end_query', 'start_target', 'end_target', 'evalue', 'score']
        for dbname in self.dblist:
            M8file = self.M8_fdict[dbname]
            
            ##FilterM8
            sheet = pd.read_csv(M8file, sep='\t', header=None, names=column_index)
            outname = os.path.join(outpath, "agm_bp_sfa_"+os.path.basename(M8file))
            if criteria != None:
                sheet = sheet.query(criteria)
            sheet = sheet.sort_values(by=sortby, ascending=False)
            
            if is_group_max == True: #sort, ascending=False, from first row, reserve query_ID first time meets
                outname = os.path.join(outpath, "agm_bp_sfg_"+os.path.basename(M8file))
                groupmax_name_dict = {}
                for i in sheet.index:
                    row = sheet.loc[i]
                    groupmax_name = row['query_ID']
                    if groupmax_name not in groupmax_name_dict:
                        groupmax_name_dict[groupmax_name] = 1
                sheet = sheet.loc[list(groupmax_name_dict.keys())]
            
            elif is_pos_max == True: #remove duplicate query_ID by start_pos & end_pos
                outname = os.path.join(outpath, "agm_bp_sfp_"+os.path.basename(M8file))
                newsheet = pd.DataFrame(columns=column_index)
                all_contigs = sheet['query_ID'].unique()
                for contig_ID in all_contigs:
                    contig_df = sheet[sheet['query_ID']==contig_ID]
                    ##check row dimension by shape
                    if contig_df.shape[0] == 1:
                        newsheet.loc[contig_df.index[0]] = contig_df.iloc[0]
                    else: #record mimimun & maximum position as covered range
                        contig_df = contig_df.sort_values(by=sortby, ascending=False)
                        newsheet.loc[contig_df.index[0]] = contig_df.iloc[0]
                        s1, e1 = sorted((contig_df.iloc[0]['start_query'], contig_df.iloc[0]['end_query']))
                        for i in contig_df.index[1:]:
                            contig_row = contig_df.loc[i]
                            s2, e2 = sorted((contig_row['start_query'], contig_row['end_query']))
                            if s1 - s2 >= noncross_len:
                                newsheet.loc[i] = contig_row
                                s1 = s2
                            elif e2 - e1 >= noncross_len:
                                newsheet.loc[i] = contig_row
                                e1 = e2
                sheet = newsheet
            
            ##Name/LengthBackpaste & AddGeneName
            target_id_list = []
            target_len_list = []
            gene_name_list = []
            old_target_ids = sheet['target_ID']
            for masked_target_ID in old_target_ids:
                target_id_list.append(dbmask_dict[masked_target_ID][0]) #targetID_original
                target_len_list.append(int(dbmask_dict[masked_target_ID][1])) #target_len
                gene_name_list.append(dbmask_dict[masked_target_ID][2])
            sheet['target_ID'] = target_id_list
            del target_id_list
            sheet.insert(2, 'gene_name', gene_name_list)
            del gene_name_list
            sheet.insert(3, 'target_len', target_len_list)
            del target_len_list
            
            sheet.to_csv(outname, index=False, sep='\t')
            self.M8_fdict[dbname] = outname
    
    def PickRESeqs(self, output_dir: str):
        printTime()
        print("Picking Risk Elements from Contigs...")
        raw_fasta_dict = readFastaFile(self.contigs_fname, to_which='dict')
        
        for dbname in self.dblist:
            M8file = self.M8_fdict[dbname]
            sheet = pd.read_csv(M8file, sep='\t')
            new_fasta_dict = {}
            for i in range(len(sheet.index)):
                contig_ID = sheet.loc[sheet.index[i], 'query_ID']
                gene_name = sheet.loc[sheet.index[i], 'gene_name']
                start_pos, end_pos = sorted((sheet.loc[sheet.index[i], 'start_query'], sheet.loc[sheet.index[i], 'end_query']))
                seq = raw_fasta_dict[contig_ID]
                new_seq_name = "{}#{}#{}#{}".format(contig_ID, gene_name, start_pos, end_pos)
                ##cut sequence by M8 start_pos & end_pos
                new_fasta_dict[new_seq_name] = seq[start_pos-1: end_pos]
            out_fasta_fname = os.path.join(output_dir, "{}.{}.fa".format(self.name_tag, dbname))
            writeFastaFromDict(new_fasta_dict, out_fasta_fname)
            self.REseq_fdict[dbname] = out_fasta_fname
    
    def Cdhit(self, output_dir: str, is_cover_old=False):
        printTime()
        print("Processing CD-HIT...")
        for dbname in self.dblist:
            REseq_fname = self.REseq_fdict[dbname]
            tmp_cdhit_REseq_fname = os.path.join(output_dir, "tmp_{}.{}.cdhit.fa".format(self.name_tag, dbname))
            cdhit_REseq_fname = os.path.join(output_dir, "{}.{}.cdhit.fa".format(self.name_tag, dbname))
            is_rename = checkCallCMD("cd-hit -i {} -o {} -c {} -aS {} -n 5 -d 0 -g 1 -M 0 -T {}".format(REseq_fname, tmp_cdhit_REseq_fname, CDHIT_IDENTITY, CDHIT_IDENTITY, NUM_THREAD), cdhit_REseq_fname, is_cover_old=is_cover_old)
            if is_rename:
                if is_cover_old and os.path.exists(cdhit_REseq_fname):
                    os.remove(cdhit_REseq_fname)
                os.rename(tmp_cdhit_REseq_fname, cdhit_REseq_fname)
            self.REseq_fdict[dbname] = cdhit_REseq_fname
    
    def GetDepth(self, output_dir: str, is_cover_old=False):
        for dbname in self.dblist:
            cdhit_REseq_fname = self.REseq_fdict[dbname]
            self.REbases_ddict[dbname] = countFastaBases(cdhit_REseq_fname)
            
            tmp_bam_fname = os.path.join(output_dir, "tmp_{}.{}.bam".format(self.name_tag, dbname))
            tmp_depth_fname = os.path.join(output_dir, "tmp_{}.{}.depth.tsv".format(self.name_tag, dbname))
            bam_fname = os.path.join(output_dir, "{}.{}.bam".format(self.name_tag, dbname))
            depth_fname = os.path.join(output_dir, "{}.{}.depth.tsv".format(self.name_tag, dbname))
            if self.reads_type not in ("nanopore", "pacbio"):
                try:
                    checkCallCMD("bwa index {}".format(cdhit_REseq_fname), depth_fname, is_cover_old=is_cover_old)
                    cmd1 = "bwa mem -t {} -v 0 {} {} {} | samtools view -bS --threads {} | samtools sort --threads {} -o {}".format(NUM_THREAD, cdhit_REseq_fname, self.reads_fname1, self.reads_fname2, NUM_THREAD, NUM_THREAD, tmp_bam_fname)
                    printTime()
                    print("Processing BWA, {} against {} {}".format(cdhit_REseq_fname, self.reads_fname1, self.reads_fname2))
                    is_rename = checkCallCMD(cmd1, depth_fname, is_cover_old=is_cover_old)
                    if is_rename:
                        if is_cover_old and os.path.exists(bam_fname):
                            os.remove(bam_fname)
                        os.rename(tmp_bam_fname, bam_fname)
                    checkCallCMD("samtools index {}".format(bam_fname), depth_fname, is_cover_old=is_cover_old)
                    cmd2 = "samtools depth -a {} | ".format(bam_fname) + '''awk 'BEGIN {OFS="\t"} {sum[$1]+=$3; count[$1]++} END {for (gene in sum) print gene, sum[gene]/count[gene]}' ''' + "> {}".format(tmp_depth_fname)
                    printTime()
                    print("Processing samtools depth, in: {}".format(bam_fname))
                    is_rename = checkCallCMD(cmd2, depth_fname, is_cover_old=is_cover_old)
                    if is_rename:
                        if is_cover_old and os.path.exists(depth_fname):
                            os.remove(depth_fname)
                        os.rename(tmp_depth_fname, depth_fname)
                except:
                    print("Failed to calculate depth.")
                    depth_fname = None
                else:
                    # os.remove(cdhit_REseq_fname+'.clstr')
                    if os.path.exists(cdhit_REseq_fname+'.amb'):
                        os.remove(cdhit_REseq_fname+'.amb')
                        os.remove(cdhit_REseq_fname+'.ann')
                        os.remove(cdhit_REseq_fname+'.bwt')
                        os.remove(cdhit_REseq_fname+'.pac')
                        os.remove(cdhit_REseq_fname+'.sa')
                    if os.path.exists(bam_fname):
                        os.remove(bam_fname)
                        os.remove(bam_fname+'.bai')
            else:
                try:
                    if self.reads_type == 'nanopore':
                        map_option = 'map-ont'
                    elif self.reads_type == 'pacbio':
                        map_option = 'map-pb'
                    cmd1 = "minimap2 -ax {} -t {} {} {} | samtools view -bS --threads {} | samtools sort --threads {} -o {}".format(map_option, NUM_THREAD, cdhit_REseq_fname, self.reads_fname1, NUM_THREAD, NUM_THREAD, tmp_bam_fname)
                    printTime()
                    print("Processing minimap2, {} against {}".format(cdhit_REseq_fname, self.reads_fname1))
                    is_rename = checkCallCMD(cmd1, depth_fname, is_cover_old=is_cover_old)
                    if is_rename:
                        if is_cover_old and os.path.exists(bam_fname):
                            os.remove(bam_fname)
                        os.rename(tmp_bam_fname, bam_fname)
                    checkCallCMD("samtools index {}".format(bam_fname), depth_fname, is_cover_old=is_cover_old)
                    cmd2 = "samtools depth -a {} | ".format(bam_fname) + '''awk 'BEGIN {OFS="\t"} {sum[$1]+=$3; count[$1]++} END {for (gene in sum) print gene, sum[gene]/count[gene]}' ''' + "> {}".format(tmp_depth_fname)
                    printTime()
                    print("Processing samtools depth, in: {}".format(bam_fname))
                    is_rename = checkCallCMD(cmd2, depth_fname, is_cover_old=is_cover_old)
                    if is_rename:
                        if is_cover_old and os.path.exists(depth_fname):
                            os.remove(depth_fname)
                        os.rename(tmp_depth_fname, depth_fname)
                except:
                    print("Failed to calculate depth.")
                    depth_fname = None
                else:
                    if os.path.exists(bam_fname):
                        os.remove(bam_fname)
                        os.remove(bam_fname+'.bai')
            self.depth_fdict[dbname] = depth_fname
    
    def RankRisk(self, outpath: str, weight_dict: dict=None) -> None:
        printTime()
        print("Generating risk rank result...")
        ##Risk_df is built with M8 files, to calculate coocur_score
        Qdb_dict_all = {}
        query_ID_all = {}
        for dbname in self.dblist:
            M8file = self.M8_fdict[dbname]
            sheet = pd.read_csv(M8file, sep='\t')
            Qdb_dict = {}
            for i in range(len(sheet.index)):
                query_ID = sheet.loc[sheet.index[i], 'query_ID']
                query_ID_all[query_ID] = 1
                if query_ID not in Qdb_dict:
                    Qdb_dict[query_ID] = 1
                else:
                    Qdb_dict[query_ID] += 1
            Qdb_dict_all[dbname] = Qdb_dict
        
        Risk_df = pd.DataFrame(columns=self.dblist, index=list(query_ID_all.keys()))
        Risk_df.index.name = 'query_ID'
        for dbname in Qdb_dict_all:
            Qdb_dict = Qdb_dict_all[dbname]
            Risk_df[dbname] = pd.Series(Qdb_dict)
        if Risk_df.empty == False:
            Risk_df = Risk_df.astype(np.float64)
            Risk_df.fillna(0.0, inplace=True)
            coocur_score = np.sum(Risk_df.values) / Risk_df.shape[0] ##
        else:
            coocur_score = 1
        
        outname_R_df = os.path.join(outpath, "RiskMatrix_{}.csv".format(self.name_tag))
        Risk_df.to_csv(outname_R_df)
        
        if weight_dict == None or len(weight_dict) != len(self.dblist):
            weight_dict = {dbname: 1 for dbname in self.dblist}
        
        db_total_weighted_depth_list = []
        for dbname in self.dblist:
            depth_fname = self.depth_fdict[dbname]
            if depth_fname:
                db_depth_df = pd.read_csv(depth_fname, sep='\t', index_col=0, header=None, names=['Element', 'Depth'])
                # db_total_weighted_depth = float(db_depth_df.sum().iloc[0]) * weight_dict[dbname] / self.total_reads_num
                db_basedepth_series = db_depth_df['Depth'] * pd.Series(self.REbases_ddict[dbname])
                db_basedepth_series = db_basedepth_series.dropna()
                db_total_weighted_depth = float(db_basedepth_series.sum()) * weight_dict[dbname] / self.total_reads_bases
            else:
                db_total_weighted_depth = 0.0
            db_total_weighted_depth_list.append(db_total_weighted_depth)
        
        risk_module = np.sqrt(np.sum(np.square(db_total_weighted_depth_list)))
        risk_index = risk_module * coocur_score
        
        outname_stat = os.path.join(outpath, "RiskStat_{}.tsv".format(self.name_tag))
        ws = open(outname_stat, 'w')
        ws.write("SampleName\t{}\tRiskModule\tCoocurScore\tRiskIndex\tReadsNum\tBasesNum\tContigsNum\n".format("\t".join(self.dblist)))
        ws.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.name_tag, "\t".join([str(v) for v in db_total_weighted_depth_list]), 
                                                           risk_module, coocur_score, risk_index, self.total_reads_num, self.total_reads_bases, self.ncontig))
        self.rank_fname = outname_R_df
    
    def ExtractRiskSeqs(self, outpath: str, mark_num_dict: dict) -> None:
        printTime()
        print("Dumping high risk contigs...")
        record_dict = readFastaFile(self.contigs_fname, to_which='dict')
        outname_fasta = os.path.join(outpath, f"RiskSeqs.{self.name_tag}.fasta")
        Risk_norm_df = pd.read_csv(self.rank_fname, index_col='query_ID')
        M8_df_dict = {}
        for dbname in self.dblist:
            M8_fname = self.M8_fdict[dbname]
            M8_df = pd.read_csv(M8_fname, sep='\t', index_col='query_ID')
            M8_df_dict[dbname] = M8_df
        
        new_record_dict = {}
        for query_ID in Risk_norm_df.index:
            # seq = record_dict[query_ID]
            risk_norm_row = Risk_norm_df.loc[query_ID]
            if len(risk_norm_row[risk_norm_row > 0]) > 1: #matched db num >1
                db_info_str_list = []
                for dbname in self.dblist:
                    M8_df = M8_df_dict[dbname]
                    dbmark = mark_num_dict[dbname]
                    
                    try:
                        use_rows = M8_df.loc[query_ID]
                    except:
                        pass
                    else:
                        if len(use_rows.shape) > 1:
                        #     use_rows = use_rows.iloc[0]
                            for idx in range(use_rows.shape[0]): #check row num
                                use_rows_single = use_rows.iloc[idx]
                                gene_name = use_rows_single['gene_name']
                                start_query = use_rows_single['start_query']
                                end_query = use_rows_single['end_query']
                                db_info_str_list.append("{}{},{},{}".format(dbmark, gene_name, start_query, end_query))
                        else:
                            gene_name = use_rows['gene_name']
                            start_query = use_rows['start_query']
                            end_query = use_rows['end_query']
                            db_info_str_list.append("{}{},{},{}".format(dbmark, gene_name, start_query, end_query)) #add db mark char before gene name
                new_seqname = "{}| {}".format(query_ID, "; ".join(db_info_str_list))
                new_record_dict[new_seqname] = record_dict[query_ID]
        writeFastaFromDict(new_record_dict, outname_fasta)
    
    def AddCateName(self, outpath: str, cate_dict: dict, sepline=False) -> None:
        printTime()
        print("Adding category annotation to BLAST result...")
        cate_list_dict = {"CARD": ['AMR Gene Family', 'Drug Class', 'Resistance Mechanism'],
                          "VFDB": ['VF_name', 'VF_bacteria', 'VF_category'],}
        checksep_dict_dict = {"CARD": {'AMR Gene Family': False, 'Drug Class': True, 'Resistance Mechanism': True},
                              "VFDB": {'VF_name': False, 'VF_bacteria': False, 'VF_category': False},}
        
        for dbname in self.dblist:
            if dbname in cate_list_dict:
                cate_names = cate_list_dict[dbname]
                cate_names_num = len(cate_names)
                M8_fname = self.M8_fdict[dbname]
                sheet = pd.read_csv(M8_fname, sep='\t')
                col_index = sheet.columns
                
                for cate_i, cate_type in enumerate(cate_names):
                    outname = os.path.join(outpath, f"{os.path.basename(M8_fname).rsplit('.', 1)[0]}_{cate_type}.tsv")
                    newsheet = pd.DataFrame(columns=col_index)
                    checksep = checksep_dict_dict[dbname][cate_type]
                    for i in range(len(sheet.index)):
                        row = sheet.iloc[i].copy()
                        cate_num = extractCateNum(row['target_ID'], dbname)
                        try:
                            cate_name = renameCateName(cate_dict[cate_num][cate_i], dbname, cate_type)
                        except:
                            cate_name = 'undefined'
                        
                        if sepline == True and checksep == True and ';' in cate_name:
                            ARnameLS = cate_name.split(';')
                        else:
                            ARnameLS = [cate_name]
                        for cate_name in ARnameLS:
                            row['gene_name'] = cate_name
                            newsheet = pd.concat([newsheet, row.to_frame().T], ignore_index=True)
                    
                    # self.M8_fdict[f"{dbname}_{cate_type}"] = outname
                    newsheet.to_csv(outname, index=False, sep='\t')
                
                add_cols = [[] for idx in range(cate_names_num)]
                outname_all = os.path.join(outpath, f"{os.path.basename(M8_fname).rsplit('.', 1)[0]}_alltype.tsv")
                
                for i in range(len(sheet.index)):
                    row = sheet.iloc[i]
                    cate_num = extractCateNum(row['target_ID'], dbname)
                    for idx in range(cate_names_num):
                        try:
                            cate_name = renameCateName(cate_dict[cate_num][idx], dbname, cate_names[idx])
                        except:
                            cate_name = 'undefined'
                        add_cols[idx].append(cate_name)
                    
                for idx in range(cate_names_num):
                    sheet.insert(idx+3, cate_names[idx], add_cols[idx])
                
                sheet.to_csv(outname_all, index=False, sep='\t')
            
            elif dbname == "MGE":
                M8_fname = self.M8_fdict[dbname]
                outname = os.path.join(outpath, f"{os.path.basename(M8_fname).rsplit('.', 1)[0]}_MGE_category.tsv")
                sheet = pd.read_csv(M8_fname, sep='\t')
                col_index = sheet.columns
                gene_names = sheet['gene_name']
                cate_name_list = []
                for genename in gene_names:
                    cate_name_list.append(renameCateName(genename, "MGE"))
                sheet['gene_name'] = cate_name_list
                
                # self.M8_fdict["MGE_category"] = outname
                sheet.to_csv(outname, index=False, sep='\t')
    
    def CalcBPM(self, outpath: str, is_sorted=False) -> None:
        printTime()
        print("Calculating BPM abundance...")
        total_bases_num = self.total_reads_bases
        for dbname in self.depth_fdict:
            depth_fname = self.depth_fdict[dbname]
            if depth_fname:
                if os.path.getsize(depth_fname) > 0:
                    outname = os.path.join(outpath, f"BPM.{self.name_tag}.{dbname}.tsv")
                    
                    sheet = pd.read_csv(depth_fname, sep='\t', header=None, index_col=0)
                    sheet.index.name = "Element"
                    sheet.columns = ["depth"]
                    seq_compo_list = [seqname.split('#') for seqname in sheet.index]
                    sheet['length'] = [abs(int(seq_compo[3]) - int(seq_compo[2])) +1 for seq_compo in seq_compo_list]
                    
                    BPM_sheet = pd.DataFrame(index=sheet.index)
                    BPM_sheet["BPM"] = sheet["depth"] * sheet["length"] * 1e6 / total_bases_num
                    
                    if is_sorted:
                        BPM_sheet = BPM_sheet.sort_values(by=BPM_sheet["BPM"], ascending=False)
                    BPM_sheet.to_csv(outname, sep='\t')


def checkCallCMD(cmd: str, outname: str, is_cover_old=False, print_skipped=True) -> bool:
    '''If outname not exists, then call cmd
    '''
    if os.path.exists(outname) == False or is_cover_old == True:
        subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    else:
        if print_skipped:
            print("[[Skipped]]")
        return False

def printTime() -> None:
    print("[{}]".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())), end=' ')

def makeDir(dirparh: str, subname: str) -> str:
    make_path = os.path.join(dirparh, subname)
    if os.path.isdir(make_path) == False:
        os.mkdir(make_path)
    return make_path

def readSheetFile(fname: str, sep='\t', remove_header=True) -> list:
    tsv_sheet = []
    f1 = open(fname, 'r')
    for line in f1.readlines():
        tsv_line = line[:-1].split(sep)
        tsv_sheet.append(tsv_line)
    f1.close()
    if remove_header == True:
        return tsv_sheet[1:] #Remove titie
    else:
        return tsv_sheet

def loadDbMaskDict(fname: str) -> dict: #pandas is slower here
    len_dict = {}
    len_sheet = readSheetFile(fname)
    for line in len_sheet:
        len_dict[line[0]] = [line[1], line[2], line[3]]
    return len_dict

def loadCateDict(fname_dict: dict) -> dict: #pandas is slower here
    cate_dict = {}
    for dbname in fname_dict:
        cate_sheet = readSheetFile(fname_dict[dbname])
        if dbname == "CARD":
            for line in cate_sheet:
                cate_dict[line[0]] = [line[8], line[9], line[10]]
        elif dbname == "VFDB":
            for line in cate_sheet:
                cate_dict[line[0]] = [line[1], line[3], line[5]]
    return cate_dict

def readFastaFile(fasta_fname: str, to_which='dict'):
    if to_which == 'dict':
        return {record.description: str(record.seq) for record in SeqIO.parse(fasta_fname, 'fasta')}
    elif to_which == 'list':
        return [record for record in SeqIO.parse(fasta_fname, "fasta")]

def writeFastaFromDict(fasta_dict: dict, outname: str) -> None:
    f2 = open(outname, 'w')
    for item in fasta_dict:
        f2.write(">{}\n".format(item))
        seq = (fasta_dict[item]).upper()
        sep_seqLS = []
        row_seqlen_max = 80
        while len(seq) >= row_seqlen_max:
            sep_seqLS.append(seq[:row_seqlen_max])
            seq = seq[row_seqlen_max:]
        if len(seq) > 0:
            sep_seqLS.append(seq)
        for sep_seq in sep_seqLS:
            f2.write(sep_seq+'\n')
    f2.close()

def renameFastaSeqs(raw_fasta_dict:dict, name_prefix: str) -> dict:
    new_fasta_dict = {}
    i = 0
    for seqname in raw_fasta_dict:
        seq = raw_fasta_dict[seqname]
        i += 1
        new_fasta_dict["{}_{}".format(name_prefix, i)] = seq
    return new_fasta_dict

def extractCateNum(target_id: str, dbname: str) -> str:
    cate_num = target_id
    if dbname == "CARD":
        cate_num = "ARO:" + target_id.split(':')[1][:7]
    elif dbname == "VFDB":
        cate_num = 'VF' + target_id.split(' (VF')[1].split(') -')[0]
    return cate_num

def renameCateName(cate_name: str, dbname: str, cate_type: str=None) -> str:
    if dbname == "CARD":
        if cate_type == "AMR Gene Family":
            if 'antibiotic efflux pump' in cate_name:
                sp = cate_name.split('(')[1].split(')')[0]
                cate_name = sp + ' type drug efflux'
            elif 'beta-lactamase' in cate_name:
                cate_name = 'beta-lactam'
            elif 'tetracycline' in cate_name:
                cate_name = 'Tetracycline'
            elif len(cate_name) < 10:
                cate_name = 'Aminoglycoside'
            else:
                cate_name = 'Others'
    
    elif dbname == "MGE":
        prefix = cate_name.split('_', 1)[0]
        # if prefix in ('plasmid', 'vir', 'proph'):
        #     cate_name = prefix #ACLAME
        if prefix[:2] == 'IS':
            cate_name = 'IS'
        elif prefix[:3] in ('Inc', 'rep', 'Col'):
            cate_name = 'plasmid'
        elif prefix[:4] == 'MITE':
            cate_name = 'IS'
        elif prefix[:2] in ('Tn', 'In', 'kappa'):
            cate_name = 'Tn'
        elif prefix[:9] == 'INTEGRALL':
            cate_name = 'INTEGRALL'
        elif prefix[0] == 'p':
            cate_name = 'plasmid'
        # else:
        #     cate_name = 'ICEberg'
    
    return cate_name

def checkGZipped(file_path):
    ''' Returns True if file is gzipped and False otherwise.
         The result is inferred from the first two bits in the file read
         from the input path.
         On unix systems this should be: 1f 8b
         Theoretically there could be exceptions to this test but it is
         unlikely and impossible if the input files are otherwise expected
         to be encoded in utf-8.
    '''
    with open(file_path, mode='rb') as fh:
        bit_start = fh.read(2)
    if(bit_start == b'\x1f\x8b'):
        return True
    else:
        return False

# def countFastqReads_Bases(fname: str, is_gzipped=True) -> tuple[int, int]:
def countFastqReads_Bases(fname: str, is_gzipped=True):
    if is_gzipped:
        handle = gzip.open(fname, "rt")
    else:
        handle = open(fname, 'r')
    total_reads_num = 0
    double_seq_len = 0
    for line in handle:
        if line[0] == '@':
            total_reads_num += 1
        if line[0] not in ('@', '+', ' ', '\n'):
            double_seq_len += (len(line)-1)
    total_seq_len = double_seq_len // 2
    return total_reads_num, total_seq_len

def countFastqReads_Bases_pigzAWK(fname: str, is_gzipped=True):
    '''Try using pigz & awk to accelerate gzipped file reading and reads & bases counting'''
    awk_cmd = r"awk 'BEGIN {reads=0;bases=0} NR%4==1 && /^@/ {reads++} NR%4==2 {bases+=length} END {print reads,bases}'"
    if is_gzipped:
        cmd = "pigz -p {} -dck {} | {}".format(NUM_THREAD, fname, awk_cmd)
    else:
        cmd = "cat {} | {}".format(fname, awk_cmd)
    
    try:
        proc = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        print("Parsing (gzipped) reads file: failed to run shell command, trying default functions.")
        total_reads_num, total_seq_len = countFastqReads_Bases(fname, is_gzipped=is_gzipped)
    else:
        output = proc.stdout.strip()
        if not output:
            raise ValueError("Wrong format of fastq file")
        total_reads_num, total_seq_len = map(int, map(float, output.split())) #e.g., '5.43923e+09' must be converted to float first, then to int
    
    if total_reads_num == 0:
        raise ValueError("Num of total reads is 0, please check input reads file.")
    return total_reads_num, total_seq_len

def countFastaBases(fasta_fname: str) -> dict:
    fasta_dict = readFastaFile(fasta_fname, to_which='dict')
    seq_len_dict = {}
    for seqname in fasta_dict:
        seq_len_dict[seqname] = len(fasta_dict[seqname])
    return seq_len_dict


if __name__ == "__main__":
    params = _getRankerArgs()
    
    t1 = time.time()
    ## Initialization
    if params.nanopore and params.pacbio:
        print("Conflict arguments '--nanopore' and '--pacbio'")
        exit()
    elif params.nanopore:
        Reads_Type = "nanopore"
    elif params.pacbio:
        Reads_Type = "pacbio"
    else:
        Reads_Type = "illumina"
    
    is_Cover_Old_File = True if params.force else False
    is_Rename_Contig = False if params.no_rename_contigs else True
    
    Out_Path = params.output
    Contigs_Fname = params.contigs
    Reads_Fname1 = params.reads1
    Reads_Fname2 = params.reads2
    
    Seq_MinLen = params.minlen
    NUM_THREAD = params.threads
    BLAST_EVALUE =params.blast_evalue
    PERC_IDENTITY = 80
    CDHIT_IDENTITY = params.cdhit_identity
    Filter_Criteria = "identity>={} & cover_len>={}".format(100*params.blast_identity, params.blast_cover_len) #for BLAST M8 result
    Colinear_Min_NonCross_Len = params.blast_cover_len #for BLAST M8 result
    
    Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
    Db_Path = os.path.join(Program_Dir_Path, "ranker_db", "db_fasta")
    Db_List = ["CARD", "MGE", "VFDB"]
    Weight_List = params.weight
    # Weight_List_Str = params.weight
    # Weight_List = [int(Weight_Str) for Weight_Str in Weight_List_Str.split(',')]
    # if len(Weight_List) != 3:
    #     print("Failed to set database weight.")
    #     exit()
    
    Weight_Dict = {Db_List[i]: Weight_List[i] for i in range(len(Db_List))}
    
    Seqname_Lendict_Fname = os.path.join(Program_Dir_Path, "ranker_db", "ranker_blastdb_seqname_length.tsv")
    Category_Fname_Dict = {"CARD": os.path.join(Program_Dir_Path, "ranker_db", "aro_index.tsv"),
                           "VFDB": os.path.join(Program_Dir_Path, "ranker_db", "VFs.tsv"),
                           }
    is_Category_Sepline = True
    
    Dump_Gene_DB_Mark_Dict = {"CARD": '1', "MGE": '2', "VFDB": '3'}
    
    ## Start Run
    PrepareBlastDB(Db_Path, Db_List, is_cover=is_Cover_Old_File)
    
    Tmp_Path = makeDir(Out_Path, "temp")
    M8_Path = makeDir(Out_Path, "output_M8")
    M8pp_Path = makeDir(Out_Path, "preprocessed_M8")
    Risk_Elements_Path = makeDir(Out_Path, "risk_elements")
    Risk_Result_Path = makeDir(Out_Path, "risk_result")
    Extract_Fasta_Path = makeDir(Out_Path, "coocur_structures")
    M8cate_Path = makeDir(M8pp_Path, "categorized_M8")
    BPM_Path = makeDir(Out_Path, "BPM")
    
    DbMask_Dict = loadDbMaskDict(Seqname_Lendict_Fname)
    Cate_Dict = loadCateDict(Category_Fname_Dict)
    
    
    Sample_Obj = Sample(Contigs_Fname, Db_Path, Db_List, Reads_Fname1, Reads_Fname2, reads_type=Reads_Type)
    Sample_Obj.SeqPreprocess(Tmp_Path, minlen=Seq_MinLen, is_rename_contigs=is_Rename_Contig)
    
    if Sample_Obj.ncontig < params.minnum:
        print("No enough contigs!")
    else:
        Sample_Obj.Blast(M8_Path, is_cover_old=is_Cover_Old_File)
        Sample_Obj.M8Preprocess(M8pp_Path, DbMask_Dict, is_pos_max=True, criteria=Filter_Criteria, noncross_len=Colinear_Min_NonCross_Len)
        Sample_Obj.PickRESeqs(Risk_Elements_Path)
        Sample_Obj.Cdhit(Risk_Elements_Path, is_cover_old=is_Cover_Old_File)
        Sample_Obj.GetDepth(Risk_Elements_Path, is_cover_old=is_Cover_Old_File)
        Sample_Obj.RankRisk(Risk_Result_Path, weight_dict=Weight_Dict)
        Sample_Obj.ExtractRiskSeqs(Extract_Fasta_Path, Dump_Gene_DB_Mark_Dict)
        Sample_Obj.AddCateName(M8cate_Path, Cate_Dict, sepline=is_Category_Sepline)
        Sample_Obj.CalcBPM(BPM_Path)
    
    t2 = time.time()
    printTime()
    print("All done. Time elapsed: {:.2f}s.".format(t2-t1))
    
