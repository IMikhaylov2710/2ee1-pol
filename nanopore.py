import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm
import argparse
import multiprocessing
from Bio.SeqRecord import SeqRecord
import sys
from Bio import AlignIO
from Bio.Align import AlignInfo, PairwiseAligner
import numpy as np
from datetime import date, time, datetime

current_directory = str(os.getcwd())
current_analysis_directory = str(os.getcwd())+'/results_Nanopore_'+str(datetime.now()).replace(' ', '_').replace(':', '').split('.')[0]+'/'

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--InPath", help = "path to input file",
                   nargs='?',
                   type = str)
parser.add_argument("-b", "--Barcodes", help = "path to barcodes file",
                   nargs='?',
                   type = str)
parser.add_argument("-r", "--Reference", help = "path to refernce file, absolute path, default is /home/imikhailov/REFs/mirna/mature.fasta",
                   nargs='?',
                   type = str, 
                   default = '/home/imikhailov/REFs/RNA_trans/Homo_sapiens.GRCh38.cdna.all.fa.gz ')
parser.add_argument("-oroot", "--OutPathRoot", help = "path to write all data into, default - this directory/results/Nanopore_CURRENTDATETIME",
                   nargs='?',
                   type = str, 
                   default = current_analysis_directory)
parser.add_argument("-umicons", "--UMIConsensus", help = "UMI consensus sequence for integrity checking",
                   nargs='?',
                   type = str, 
                   default = 'NNNYRNNNYRNNNYRNNN')
parser.add_argument("-fwadapter", "--ForwardAdapterSeq", help = "sequence of forward adapter, default = CAAGCAGAAGACGGCATACGAGAT",
                   nargs='?',
                   type = str, 
                   default = 'CAAGCAGAAGACGGCATACGAGAT')
parser.add_argument("-revadapter", "--ReverseAdapterSeq", help = "sequence of reverse adapter, default = AATGATACGGCGACCACCGAGATC",
                   nargs='?',
                   type = str, 
                   default = 'AATGATACGGCGACCACCGAGATC')
parser.add_argument("-mmcount", "--MisMatch", help = "number of mismatches allowed in Hamming adapter search, default is 5",
                   nargs='?',
                   type = int, 
                   default =5)
parser.add_argument("-mmcountindex", "--MisMatchIndex", help = "number of mismatches allowed in Hamming index search, also with incomplete index sequence, default is 3",
                   nargs='?',
                   type = int, 
                   default = 3)
parser.add_argument("-p", "--PIdent", help = "clustering parameters, default = 90 percent",
                   nargs='?',
                   type = float, 
                   default = 0.9)
parser.add_argument("-m", "--Mapping", help = "mapping parameters, possible arguments - blast, bwa and minimap2, default - blast",
                   nargs='?',
                   type = str, 
                   default = 'blast')
parser.add_argument("-frags", "--MinAsignedFrags", help = "number of assigned frags for salmon quantification, default is 1, do not change it",
                   nargs='?',
                   type = str,
                   default = '1')
parser.add_argument("-ws", "--WordSize", help = "word size for blast quantification",
                   nargs='?',
                   type = str,
                   default = '6')
parser.add_argument("-batch", "--BatchSize", help = "size of batch to process simultaneously, default to 8. Requires fine-tuning",
                   nargs='?',
                   type = str,
                   default = '8')
parser.add_argument("-t", "--Threads", help = "number of threads to use for non-parallel quantifications, default - 16",
                   nargs='?',
                   type = str,
                   default = '16')
parser.add_argument("--trim", action = "store_true", help = "trim poly-A and poly-T")
parser.add_argument("--noquant", action = "store_true", help = "do not preform quantification for analysed data")
parser.add_argument("--muffle", action = "store_true", help = "omit all outputs")
args = parser.parse_args()

#colors for messages
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'  
    
umiis = []    
    
OutPathDemultiplex = os.path.join(args.OutPathRoot, 'demultiplexed/')   
OutPathUMI = os.path.join(args.OutPathRoot, 'UMI/')    
OutPathReduced = os.path.join(args.OutPathRoot, 'reduced/')
OutPathMapped = os.path.join(args.OutPathRoot, 'mapped/')
OutPathTmp = os.path.join(args.OutPathRoot, 'tmp/')    

if not os.path.exists(args.OutPathRoot):
    os.system('mkdir %s' % args.OutPathRoot)
if not os.path.exists(os.path.join(args.OutPathRoot, 'demultiplexed/')):
    os.system('mkdir %s' % os.path.join(args.OutPathRoot, 'demultiplexed/'))
    OutPathDemultiplex = os.path.join(args.OutPathRoot, 'demultiplexed/')
if not os.path.exists(os.path.join(args.OutPathRoot, 'UMI/')):
    os.system('mkdir %s' % os.path.join(args.OutPathRoot, 'UMI/'))
    OutPathUMI = os.path.join(args.OutPathRoot, 'UMI/')
if not os.path.exists(os.path.join(args.OutPathRoot, 'reduced/')):
    os.system('mkdir %s' % os.path.join(args.OutPathRoot, 'reduced/'))
    OutPathReduced = os.path.join(args.OutPathRoot, 'reduced/')
if not os.path.exists(os.path.join(args.OutPathRoot, 'mapped/')):
    os.system('mkdir %s' % os.path.join(args.OutPathRoot, 'mapped/'))
    OutPathMapped = os.path.join(args.OutPathRoot, 'mapped/')
if not os.path.exists(os.path.join(args.OutPathRoot, 'tmp/')):
    os.system('mkdir %s' % os.path.join(args.OutPathRoot, 'tmp/'))
    OutPathTmp = os.path.join(args.OutPathRoot, 'tmp/')    
    
ref_umi_fw = args.UMIConsensus
ref_umi_rv = args.UMIConsensus

nucleotides_degenerate = {}
nucleotides_degenerate['N'] = ['A', 'T', 'G', 'C']
nucleotides_degenerate['Y'] = ['T', 'C']
nucleotides_degenerate['R'] = ['A', 'G']

#beta
def CalculateBatchSize(threads):
    batchSize = math.floor(threads/2)
    return batchSize

def GenerateBatch(fils):
    c = 0
    batch = 1
    batch_dic = {}
    for fil in fils:
        c+=1
        if c%int(args.BatchSize) != 0:
            try:
                batch_dic[batch].append(fil)
            except:
                batch_dic[batch] = [fil]
        else:
            batch_dic[batch].append(fil)
            batch +=1
    return batch_dic

def CallbackMultiprocess(generatedBatchDic, batchNumber, function, *functionArgs):
    print(f'{bcolors.OKGREEN}started working with batch %s/%s' % (batchNumber, len(generatedBatchDic)))
    p = multiprocessing.Pool()
    c=0
    for fil in generatedBatchDic[batchNumber]:
        c+=1
        print(f'{bcolors.OKBLUE}\tStarting subprocess %s/%s' % (c, len(generatedBatchDic[batchNumber])) + bcolors.ENDC)
        p.apply_async(function, [fil, functionArgs])
    p.close()
    p.join()

def LoadBarcodes(path):
    with open(path, 'r') as handle:
        bar_fw = {}
        bar_rv = {}
        for lin in handle:
            lins = lin.strip().split(' ')
            if 'fw' in lins[0]:
                bar_fw[lins[1]] = lins[0]
            elif 'rv' in lins[0]:
                bar_rv[str(Seq(lins[1]).reverse_complement())] = lins[0]
    return bar_fw, bar_rv

def LoadBarcodesRev(path):
    with open(path, 'r') as handle:
        bar_fw = {}
        bar_rv = {}
        for lin in handle:
            lins = lin.strip().split(' ')
            if 'fw' in lins[0]:
                bar_rv[str(Seq(lins[1]).reverse_complement())] = lins[0]
            elif 'rv' in lins[0]:
                bar_fw[lins[1]] = lins[0]
    return bar_fw, bar_rv

def CheckUmiScore(umi, degref):
    umi_score = 0
    for e, let in enumerate(umi):
        if let in nucleotides_degenerate[degref[e]]:
            umi_score+=1
    return umi_score, len(umi)

def Hamming_check(reference, sequence):
    errors = 0
    for i in range(len(reference)):
        if reference[i] != sequence[i]:
            errors += 1
    return errors 

def Hamming_check_incomplete(reference, sequence):
    errors = 0
    for i in range(len(sequence)):
        if reference[i] != sequence[i]:
            errors += 1
    return errors 

def GenerateWindowsList(seq, length):
    windows = []
    for i in range(len(seq)-length):
        candidate = seq[i:i+length]
        windows.append(candidate)
    return windows

def GenerateUMI(sequence, UMIid):
    my_seq = SeqRecord(Seq(sequence), id = str(UMIid))
    return my_seq

def RightTrimT(sequence):
    tcounter = 0
    trimmed_status = 'NA'
    for e, l in enumerate(sequence.seq):
        if l == 'T':
            tcounter +=1
        else:
            if tcounter > 6:
                new_rec = sequence[e:]
                trimmed_status = 'trimmed'
                tcounter = 0
            else:
                tcounter = 0
    if trimmed_status == 'NA':
        TrimmedSeq = sequence
    else:
        TrimmedSeq = new_rec
        trimmed_status = 'NA'
    return TrimmedSeq
    
def LeftTrimA(sequence):
    notacounter = 0
    acounter = 0
    trimmed_status = 'NA'
    for e, l in enumerate(sequence.seq):
        if trimmed_status == 'NA':
            if l != 'A':
                notacounter += 1
                acounter = 0
            else:
                acounter +=1
            if acounter > 6:
                new_rec = sequence[:e-6]
                trimmed_status = 'trimmed'
    if trimmed_status == 'NA':
        TrimmedSeq = sequence
    else:
        if len(new_rec.seq) == 0:
            TrimmedSeq = sequence
        else:
            TrimmedSeq = new_rec
        trimmed_status = 'NA'
    return TrimmedSeq
 
demultiplexMetaInfo = {}
    
def FindFwOrientation(file, bf, br):
    demultiplexed = {}
    demultiplexed_counts = {}
    rv = Seq(args.ReverseAdapterSeq)
    rvrc = rv.reverse_complement()
    recs = [record for record in SeqIO.parse(file, 'fastq')]
    c = 0
    fw = 0
    rv = 0
    for record in tqdm(recs):
        fw_found = 0
        rv_found = 0
        windows = GenerateWindowsList(str(record.seq), len(args.ForwardAdapterSeq))      
        for e, w in enumerate(windows):
            errsfw = Hamming_check(w, args.ForwardAdapterSeq)
            errsrv = Hamming_check(w, str(rvrc))
            if errsfw < args.MisMatch:
                start = e
                start_insert = e+len(args.ForwardAdapterSeq)+len(ref_umi_fw)
                fw_found = 1
                fw+=1
                barcode_fw = str(record.seq[e-12:e])
                umi_fw = str(record.seq[e+len(args.ForwardAdapterSeq):e+len(args.ForwardAdapterSeq)+18])
            if errsrv < args.MisMatch:
                end = e+20
                end_insert = e-len(ref_umi_rv)
                rv_found = 1
                rv+=1
                barcode_rv = str(record.seq[e+len(rvrc):e+len(rvrc)+12])
                umi_rv = str(record.seq[e-18:e])
        if fw_found == 1 and rv_found ==1:
            if len(umi_fw) == len(ref_umi_fw) and len(umi_rv) == len(ref_umi_rv):
                c+=1
                demultiplexMetaInfo[record.id] = 'fw'
                fw_umi_integrity = CheckUmiScore(umi_fw, ref_umi_fw)
                rv_umi_integrity = CheckUmiScore(umi_rv, ref_umi_fw)
                insert = record[start_insert:end_insert]
                if len(insert) > 0:
                    InsertTrimmed = LeftTrimA(insert)
                    if not len(InsertTrimmed) == 0:
                        if args.muffle != True:
                            print('===amplicon===')
                            print(record.seq[start:end])
                            print('===insert===')
                            print(insert.seq+'>'+InsertTrimmed.seq)
                            print('===META INFO===')
                            print('barcodes: fw %s, rv %s, umi_fw %s, umi_rv %s' % (barcode_fw, 
                                                                                   barcode_rv, 
                                                                                   umi_fw, 
                                                                                   umi_rv))
                            print('UMI integrity: fw %s/%s, rv %s/%s' % (fw_umi_integrity[0], 
                                                                         fw_umi_integrity[1],
                                                                         rv_umi_integrity[0],
                                                                         rv_umi_integrity[1]))
                        total_umi = umi_fw+umi_rv
                        umifasta = GenerateUMI(total_umi, record.id)
                        if args.muffle != True:
                            print(umifasta.id, len(umifasta.seq))
                        umiis.append(umifasta)
                        pr1 = 0
                        pr2 = 0
                        if len(barcode_fw)==12 and len(barcode_rv)==12:
                            for b in bf.keys():
                                ers = Hamming_check(b, barcode_fw)
                                if ers < args.MisMatchIndex:
                                    pr1 = bf[b]
                            for b in br.keys():
                                ers = Hamming_check(b, barcode_rv)
                                if ers < args.MisMatchIndex:
                                    pr2 = br[b]
                        elif len(barcode_fw)>8 and len(barcode_rv)>8:
                            for b in bf.keys():
                                ers = Hamming_check_incomplete(b, barcode_fw)
                                if ers < args.MisMatchIndex:
                                    pr1 = bf[b]
                            for b in br.keys():
                                ers = Hamming_check_incomplete(b, barcode_rv)
                                if ers < args.MisMatchIndex:
                                    pr2 = br[b]
                        elif args.muffle != True:
                            print(len(barcode_fw), len(barcode_rv), 'omitted')
                        if pr1 != 0 and pr2 != 0:
                            if args.trim == True:
                                try:
                                    demultiplexed_counts[pr1+','+pr2]+=1
                                    demultiplexed[pr1+','+pr2].append(InsertTrimmed)
                                except:
                                    demultiplexed_counts[pr1+','+pr2]=1
                                    demultiplexed[pr1+','+pr2] = [InsertTrimmed]
                            else:
                                try:
                                    demultiplexed_counts[pr1+','+pr2]+=1
                                    demultiplexed[pr1+','+pr2].append(insert)
                                except:
                                    demultiplexed_counts[pr1+','+pr2]=1
                                    demultiplexed[pr1+','+pr2] = [insert]
                
    print('%s reads processed, %s with fw adapter, %s with rv adapter, %s with both' % (len(recs), fw, rv, c))
    return demultiplexed, demultiplexed_counts
            
def FindRvOrientation(file, bf, br):
    demultiplexed = {}
    demultiplexed_counts = {}
    rv = Seq(args.ForwardAdapterSeq)
    rvrc = rv.reverse_complement()
    recs = [record for record in SeqIO.parse(file, 'fastq')]
    c = 0
    fw = 0
    rv = 0
    for record in tqdm(recs):
        fw_found = 0
        rv_found = 0
        windows = GenerateWindowsList(str(record.seq), len(args.ReverseAdapterSeq))
        for e, w in enumerate(windows):
            errsfw = Hamming_check(w, args.ReverseAdapterSeq)
            if errsfw < args.MisMatch:
                start = e
                start_insert = e+len(args.ReverseAdapterSeq)+len(ref_umi_rv)
                fw_found = 1
                fw+=1
                barcode_fw = str(record.seq[e-12:e])
                umi_fw = str(record.seq[e+len(args.ReverseAdapterSeq):e+len(args.ReverseAdapterSeq)+len(ref_umi_rv)].reverse_complement())
            errsrv = Hamming_check(w, str(rvrc))
            if errsrv < args.MisMatch:
                end = e+20
                end_insert = e-len(ref_umi_rv)
                rv_found = 1
                rv+=1
                barcode_rv = str(record.seq[e+len(rvrc):e+len(rvrc)+12])
                umi_rv = str(record.seq[e-len(ref_umi_fw):e].reverse_complement())
        if fw_found == 1 and rv_found ==1:
            if len(umi_fw) == len(ref_umi_fw) and len(umi_rv) == len(ref_umi_fw):
                demultiplexMetaInfo[record.id] = 'rv'
                c+=1
                fw_umi_integrity = CheckUmiScore(umi_fw, ref_umi_fw)
                rv_umi_integrity = CheckUmiScore(umi_rv, ref_umi_fw)
                insert = record[start_insert:end_insert]
                if len(insert) > 0:
                    InsertTrimmed = RightTrimT(insert)
                    if not len(InsertTrimmed) == 0:
                        if args.muffle != True:
                            print('===amplicon===')
                            print(record.seq[start:end])
                            print('===insert===')
                            print(insert.seq+'>'+InsertTrimmed.seq)
                            print('===META INFO===')
                            print('barcodes: fw %s, rv %s, umi_fw %s, umi_rv %s' % (barcode_fw, 
                                                                                   barcode_rv, 
                                                                                   umi_fw, 
                                                                                   umi_rv))
                            print('UMI integrity: fw %s/%s, rv %s/%s' % (fw_umi_integrity[0], 
                                                                         fw_umi_integrity[1],
                                                                         rv_umi_integrity[0],
                                                                         rv_umi_integrity[1]))
                        total_umi = umi_fw+umi_rv
                        umifasta = GenerateUMI(str(Seq(total_umi).reverse_complement()), record.id)
                        if args.muffle != True:
                            print(umifasta.id, len(umifasta.seq))
                        umiis.append(umifasta)
                        pr1 = 0
                        pr2 = 0
                        if len(barcode_fw)==12 and len(barcode_rv)==12:
                            for b in bf.keys():
                                ers = Hamming_check(b, barcode_fw)
                                if ers < args.MisMatchIndex:
                                    pr2 = bf[b]
                            for b in br.keys():
                                ers = Hamming_check(b, barcode_rv)
                                if ers < args.MisMatchIndex:
                                    pr1 = br[b]
                        elif len(barcode_fw)>8 and len(barcode_rv)>8:
                            for b in bf.keys():
                                ers = Hamming_check_incomplete(b, barcode_fw)
                                if ers < args.MisMatchIndex:
                                    pr2 = bf[b]
                            for b in br.keys():
                                ers = Hamming_check_incomplete(b, barcode_rv)
                                if ers < args.MisMatchIndex:
                                    pr1 = br[b]
                        elif args.muffle != True:
                            print(len(barcode_fw), len(barcode_rv), 'omitted')
                        if pr1 != 0 and pr2 != 0:
                            if args.trim == True:
                                try:
                                    demultiplexed_counts[pr1+','+pr2+'-rv']+=1
                                    demultiplexed[pr1+','+pr2+'-rv'].append(InsertTrimmed)
                                except:
                                    demultiplexed_counts[pr1+','+pr2+'-rv']=1
                                    demultiplexed[pr1+','+pr2+'-rv'] = [InsertTrimmed]
                            else:
                                try:
                                    demultiplexed_counts[pr1+','+pr2+'-rv']+=1
                                    demultiplexed[pr1+','+pr2+'-rv'].append(insert)
                                except:
                                    demultiplexed_counts[pr1+','+pr2+'-rv']=1
                                    demultiplexed[pr1+','+pr2+'-rv'] = [insert]
                                
    print('%s reads processed, %s with fw adapter, %s with rv adapter, %s with both' % (len(recs), fw, rv, c))
    return demultiplexed, demultiplexed_counts

def WriteFastq(DemultiplexedData, path):
    for k in DemultiplexedData:
        SeqIO.write(DemultiplexedData[k], os.path.join(path, k), 'fastq')
        
def ClusterUMIData(file, percent):
    os.system('cd-hit -i %s -o %s -c %s' % (os.path.join(OutPathUMI, file), 
                                            os.path.join(OutPathUMI, file.split('.')[0]), 
                                                        percent))
def GetClusteringInfo(ClusteringData):
    ClusteredData = {}
    with open(os.path.join(OutPathUMI, ClusteringData), 'r') as handle: 
        for lin in handle:
            lins = lin.strip()
            if lins.startswith('>'):
                cluster = lins
            else:
                name = lins.split('\t')[1].split('>')[1].split('...')[0]
                ClusteredData[name] = cluster
    return ClusteredData

def ReduceByClusteringInfo(ClusteredData, DemultiplexedFastqFile, ReducedFastqPath, counter):
    recs = [record for record in SeqIO.parse(os.path.join(OutPathDemultiplex, DemultiplexedFastqFile), 'fastq')]
    ClusteredRecords = {}
    check = False
    for record in recs:
        for k in ClusteredData:
            check = True
            if k in str(record.id):
                try:
                    ClusteredRecords[ClusteredData[k]].append(record)
                except:
                    ClusteredRecords[ClusteredData[k]] = [record]
        if args.muffle != True:
            print(check)
    ClusteredResultingRecords = []
    c = 0 
    for k in ClusteredRecords:
        if len(ClusteredRecords[k]) > 1:
            c+=1
            name = os.path.join(OutPathTmp, 'test_'+str(c)+'.tmp.fasta')
            name_to_be = name.split('.')[0]+'.aln.tmp.fasta'
            SeqIO.write(ClusteredRecords[k], name, 'fasta')
            if not os.path.exists(name_to_be):
                os.system('mafft --auto %s > %s' % (name, name_to_be))
            alignment = AlignIO.read(name_to_be, 'fasta')
            summary_align = AlignInfo.SummaryInfo(alignment)
            consensus = summary_align.dumb_consensus(threshold = 0.6)
            ConsensusScores = {}
            for kk in ClusteredRecords[k]:
                if not len(kk.seq) == 0:
                    aligner = PairwiseAligner()
                    alignments = aligner.align(str(kk.seq), consensus.upper())
                    alignment = alignments[0]
                    try:
                        ConsensusScores[alignment.score].append(kk)
                    except:
                        ConsensusScores[alignment.score] = [kk]
            if not ConsensusScores:
                SortedInstance = ConsensusScores[sorted(ConsensusScores)[-1]]
                QScoreDict = {}
                for record in SortedInstance:
                    QScoreDict[np.mean(record.letter_annotations["phred_quality"])] = record
                ClusteredResultingRecords.append(QScoreDict[sorted(QScoreDict, reverse = True)[0]])
        else:
            ClusteredResultingRecords.append(ClusteredRecords[k][0])
    SeqIO.write(ClusteredResultingRecords, os.path.join(ReducedFastqPath, DemultiplexedFastqFile.split('.')[0])+'.reduced.fastq', 'fastq')

#only for bwa-based mapping quantification
def MapReducedReadsMM(fil, inpath, outpath, ref):
    os.system('minimap2 -ax map-ont %s %s > %s' % (ref, os.path.join(inpath, fil), os.path.join(outpath, fil.split('.')[0]+'.sam')))
    
def MapReducedReadsBwa(fil, inpath, outpath, ref):
    os.system('bwa mem -t 2 -x ont2d %s %s > %s' % (ref, os.path.join(inpath, fil), os.path.join(outpath, fil.split('.')[0]+'.sam')))
    
def QuantifyReads(fil, inpath, outpath, ref, frags):
    os.system('salmon quant -t %s -l A -a %s -o %s --minAssignedFrags %s' % (ref, os.path.join(inpath, fil), os.path.join(outpath, fil.split('.')[0]), frags))
    
#for blast-based quantification

def ConvertFormat(fil, inpath, outpath):
    os.system('seqtk seq -a %s > %s' % (os.path.join(inpath, fil), os.path.join(outpath, fil.split('.')[0]+'.fasta')))

def BlastReads(fil, db, inpath, outpath, word_size, threads):
    os.system('blastn -db %s -query %s -outfmt 6 -word_size %s -out %s -num_threads %s' % (db, os.path.join(inpath, fil), word_size, os.path.join(outpath, fil.split('.')[0]+'.tsv'), args.Threads))

def ParseBlast(fil):
    with open(fil, 'r') as handle:
        identity = {}
        length_of_alignement = {}
        found  = []
        stat = {}
        c = 0
        for lin in handle:
            lins = lin.strip().split('\t')
            if not lins[0] in found:
                c+=1
                try:
                    stat[lins[1]]+=1
                except:
                    stat[lins[1]] = 1
                found.append(lins[0])
                try:
                    identity[lins[2]]+=1
                except:
                    identity[lins[2]]=1
                try:
                    length_of_alignement[lins[3]]+=1
                except:
                    length_of_alignement[lins[3]]=1
    stat_norm = {}
    for k in stat:
        stat_norm[k] = [stat[k]/c*100, stat[k]]
    stat_norm_df = pd.DataFrame.from_dict(stat_norm, orient = 'index', columns = ['percent', 'reads'])
    identity_df = pd.DataFrame.from_dict(identity, orient = 'index', columns = ['reads'])
    length_of_alignement_df = pd.DataFrame.from_dict(length_of_alignement, orient = 'index', columns = ['reads'])
    return stat_norm_df, identity_df, length_of_alignement_df
    
    
starttime = datetime.now()
print('started at %s' % starttime)
print('performing barcodes loading at %s' % datetime.now())
barcodesfw = LoadBarcodes(args.Barcodes)
barcodesrv = LoadBarcodesRev(args.Barcodes)
print('barcodes loaded')
print('starting demultiplexing in fw orientation at %s' % datetime.now()) 
barsfw = FindFwOrientation(args.InPath, barcodesfw[0], barcodesfw[1])
print(barsfw[1])
print('starting demultiplexing in rv orientation at %s' % datetime.now()) 
barsrv = FindRvOrientation(args.InPath, barcodesrv[0], barcodesrv[1])
print(barsrv[1])
SeqIO.write(umiis, os.path.join(OutPathUMI, 'UMIs.fasta'), 'fasta')
print('demultiplection performed')
print('writing resulting trimmed fastq files at %s' % datetime.now())
WriteFastq(barsfw[0], OutPathDemultiplex)
WriteFastq(barsrv[0], OutPathDemultiplex)
print('done')
print('clustering UMI at %s' % datetime.now())
ClusterUMIData('UMIs.fasta', args.PIdent)
print('done')
print('getting clustering info %s' % datetime.now())
clust = GetClusteringInfo('UMIs.clstr')
print('started reducing reads to representative sequences at %s' % datetime.now()) 

fils = [fil for fil in os.listdir(OutPathDemultiplex)]

subfils = GenerateBatch(fils)
gc = 0 
for k in subfils:
    print(f'{bcolors.OKGREEN}started working with batch %s/%s' % (k, len(subfils)))
    p = multiprocessing.Pool()
    c=0
    for fil in subfils[k]:
        c+=1
        gc+=1
        print(f'{bcolors.OKBLUE}\tStarting subprocess %s/%s' % (c, len(subfils[k])) + bcolors.ENDC)
        p.apply_async(ReduceByClusteringInfo, [clust, fil, OutPathReduced, gc])
    p.close()
    p.join()
print(f'{bcolors.OKGREEN}done\n' + bcolors.ENDC)

print('started mapping at %s' % datetime.now())
fils = [fil for fil in os.listdir(OutPathReduced) if fil.endswith('.reduced.fastq')]
subfils = GenerateBatch(fils)
if args.Mapping == 'blast':
    print('converting format at %s' % datetime.now())
    for k in subfils:
        print(f'{bcolors.OKGREEN}started working with batch %s/%s' % (k, len(subfils)))
        p = multiprocessing.Pool()
        c=0
        for fil in subfils[k]:
            c+=1
            print(f'{bcolors.OKBLUE}\tStarting subprocess %s/%s' % (c, len(subfils[k])) + bcolors.ENDC)
            p.apply_async(ConvertFormat, [fil, OutPathReduced, OutPathReduced])
        p.close()
        p.join()
    print(f'{bcolors.OKGREEN}done\n' + bcolors.ENDC)  
    print('performing blast quantification at %s' % datetime.now())
    fils = [fil for fil in os.listdir(OutPathReduced) if fil.endswith('.fasta')]
    subfils = GenerateBatch(fils)      
    for k in subfils:
        print(f'{bcolors.OKGREEN}started working with batch %s/%s' % (k, len(subfils)))
        p = multiprocessing.Pool()
        c=0
        for fil in subfils[k]:
            c+=1
            print(f'{bcolors.OKBLUE}\tStarting subprocess %s/%s' % (c, len(subfils[k])) + bcolors.ENDC)
            p.apply_async(BlastReads, [fil, args.Reference, OutPathReduced, OutPathReduced, args.WordSize, 2])
        p.close()
        p.join()
    print(f'{bcolors.OKGREEN}done\n' + bcolors.ENDC)  
    print('parsing blast results at %s' % datetime.now())
    fils = [fil for fil in os.listdir(OutPathReduced) if fil.endswith('.tsv')]
    for fil in tqdm(fils):
        parsed = ParseBlast(os.path.join(OutPathReduced, fil))
        with pd.ExcelWriter(os.path.join(args.OutPathRoot, fil.split('.')[0]+'.xlsx')) as writer:
            parsed[0].to_excel(writer, sheet_name="GeneralStatistics")
            parsed[1].to_excel(writer, sheet_name="IdentityStatistics")
            parsed[2].to_excel(writer, sheet_name="LengthStatistics")
    print(f'{bcolors.OKGREEN}done\n' + bcolors.ENDC) 
            
elif args.Mapping == 'minimap2' or args.Mapping == 'mm2':
    for k in subfils:
        print(f'{bcolors.OKGREEN}started working with batch %s/%s' % (k, len(subfils)))
        p = multiprocessing.Pool()
        c=0
        for fil in subfils[k]:
            c+=1
            print(f'{bcolors.OKBLUE}\tStarting subprocess %s/%s' % (c, len(subfils[k])) + bcolors.ENDC)
            p.apply_async(MapReducedReadsMM, [fil, OutPathReduced, OutPathMapped, args.Reference])
        p.close()
        p.join()
    print(f'{bcolors.OKGREEN}done\n' + bcolors.ENDC)
elif args.Mapping =='bwa':
    for k in subfils:
        print(f'{bcolors.OKGREEN}started working with batch %s/%s' % (k, len(subfils)))
        p = multiprocessing.Pool()
        c=0
        for fil in subfils[k]:
            c+=1
            print(f'{bcolors.OKBLUE}\tStarting subprocess %s/%s' % (c, len(subfils[k])) + bcolors.ENDC)
            p.apply_async(MapReducedReadsBwa, [fil, OutPathReduced, OutPathMapped, args.Reference])
        p.close()
        p.join()
    print(f'{bcolors.OKGREEN}done\n' + bcolors.ENDC)
            
if not args.noquant and args.Mapping != 'blast':
    print('started quantification at %s' % datetime.now())
    fils = [fil for fil in os.listdir(OutPathMapped) if fil.endswith('.sam')]
    subfils = GenerateBatch(fils)
    for k in subfils:
        print(f'{bcolors.OKGREEN}started working with batch %s/%s' % (k, len(subfils)))
        p = multiprocessing.Pool()
        c=0
        for fil in subfils[k]:
            c+=1
            print(f'{bcolors.OKBLUE}\tStarting subprocess %s/%s' % (c, len(subfils[k])) + bcolors.ENDC)
            p.apply_async(QuantifyReads, [fil, OutPathMapped, OutPathMapped, args.Reference, args.MinAsignedFrags])
        p.close()
        p.join()
    print(f'{bcolors.OKGREEN}done\n' + bcolors.ENDC)