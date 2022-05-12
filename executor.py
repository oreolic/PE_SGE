#%%
import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor


#%%

class general:
    def __init__(self):
        self.codon = {'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu', 'CTT': 'Leu', 'CTC': 'Leu','CTA': 'Leu',
            'CTG': 'Leu', 'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met', 'GTT': 'Val',
            'GTC': 'Val','GTA': 'Val', 'GTG': 'Val', 'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
            'CCT': 'Pro','CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro', 'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr',
            'ACG': 'Thr','GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala', 'TAT': 'Tyr', 'TAC': 'Tyr',
            'TAA': 'STOP','TAG': 'STOP', 'CAT': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln', 'AAT': 'Asn',
            'AAC': 'Asn','AAA': 'Lys', 'AAG': 'Lys', 'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
            'TGT': 'Cys','TGC': 'Cys', 'TGA': 'STOP', 'TGG': 'Trp', 'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg',
            'CGG': 'Arg','AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg', 'GGT': 'Gly', 'GGC': 'Gly',
            'GGA': 'Gly', 'GGG': 'Gly'}
        self.codon_short = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 
                        'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V','STOP':'*'}
        self.codon_abb = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R',
                       'AGT': 'S', 'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I', 'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P',
                       'CCT': 'P', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E',
                       'GAT': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V',
                       'GTT': 'V', 'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TGA': '*', 'TGC': 'C', 'TGG': 'W',
                       'TGT': 'C', 'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}

    def read_fasta(self,acc):
        hg19 = pd.read_csv('EssentialData/hg19_Accession.txt',sep='\t',index_col=1)
        hg38 = pd.read_csv('EssentialData/hg38_Accession.txt',sep='\t',index_col=1)

        if acc in list(hg19.index):
            num = hg19.loc[acc,'Chromosome'].split()[1]        
            fasta = self._fasta19(num)
            return fasta
        elif acc in list(hg38.index):   
            num = hg38.loc[acc,'Chromosome'].split()[1]           
            fasta = self._fasta38(num)
            return fasta
        else:
            return 'FALSE'
        

        
    def _fasta38(self,num):
        total = ''
        with open('EssentialData/hg38/chr{}.fa'.format(num) , 'r') as fasta:
            a = fasta.readlines()

            n = 1
            while n < len(a):
                eachline = a[n].rstrip()
                total += eachline
                n += 1

        return total

    def _fasta19(self,num):
        total = ''
        with open('EssentialData/hg19/hg19_{}.fasta'.format(num) , 'r') as fasta:
            a = fasta.readlines()
            n = 1
            while n < len(a):
                eachline = a[n].rstrip()
                total += eachline
                n += 1
        return total


    def rt_sequence(self,seq):
        rt = ''
        dic = {'A':'T','G':'C','T':'A','C':'G','N':'N','a':'t','g':'c','t':'a','c':'g'}
        for i in seq[::-1]:
            rt += dic[i]
        return rt



    def pam_sequence(self,pamtype):
        if pamtype == 'NGG':
            pam_sequence = ['NGG']
        elif pamtype == 'NRCH':
            fp_pam = []
            for i in ['AA','CA','GA','TA','AC','CC','GC','TC','AG','CG' ,'GG', 'TG', 'AT', 'CT' ,'GT', 'TT']:
                for k in ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']:
                    pam = i+k
                    fp_pam.append(pam)
            rt_pam = [self.rt_sequence(i) for i in fp_pam]
            
        return fp_pam,rt_pam

class subsitution():
    def single_snv_in_one_codon(self,seq):
        lst = []
        for i in ['a','t','g','c']:
            lst.append((0,i+seq[1:]))
            lst.append((1,seq[0]+i+seq[-1]))
            lst.append((2,seq[:2]+i))

        return lst





class extract_cds_from_all_variants:
    
    def _extract_all_cds(self,eachdf):
        final = {}
        for i in eachdf.index:
            nc_accession = eachdf.loc[i,'nc_accession']
            fasta = general().read_fasta(nc_accession)
            
            strand = eachdf.loc[i,'cds_strand']

            cds_loc_lst = eachdf.loc[i,'cds_locations']
            cds_loc_lst = cds_loc_lst[1:-1]
            cds_loc_lst = cds_loc_lst.split(', ')

            dic = {}
                
            for idx,cds in enumerate(cds_loc_lst):
                cds = cds.split('-')
                start = int(cds[0])
                end = int(cds[1])                

                if strand == '-':
                    read = fasta[start:end+1]
                    read = general().rt_sequence(read)
                    idx = len(cds_loc_lst)-idx
                    dic[idx] = [start,end]
                    dic[idx].append(read)

                    extended = fasta[start-100:end+1+100]
                    extended = general().rt_sequence(extended)

                    dic[idx].append(extended)

                else:
                    read = fasta[start:end+1]
                    dic[idx+1] = [start,end]
                    dic[idx+1].append(read)

                    extended = fasta[start-100:end+100]
                    dic[idx+1].append(extended)

            final[i] = dic
        return final

    def _make_each_variant_df(self,eachvar_dic):
        n = 1
        k = 0
        newdic = {}
        while n < len(eachvar_dic)+1:
            cds_info = eachvar_dic[n]
            seq = cds_info[2]
            if len(seq) < 3:
                return pd.DataFrame()
            else:
                if n == 1:
                    if seq[:3] != 'ATG':
                        return pd.DataFrame()
                    else:
                        pass
                else:
                    pass


            l = k
            if l == 0:
                true_cds = seq
            else:
                true_cds = seq[3-l:]

            k = self._check_cds_remnant(seq,k)
            if k == 0:
                true_cds = true_cds
            else:
                true_cds = true_cds[:-k]

            cds_info.append(k)
            cds_info.append(true_cds)
            
            newdic[n] = cds_info
            n += 1

        lst = []
        for i in newdic:
            each_info = eachvar_dic[i]
            lst.append(each_info)
        df = pd.DataFrame(lst,columns = ['start','end','CDS_original','CDS_extended','Remnant_NT','trimmed_CDS'])
        return df

    def _check_cds_remnant(self,seq,prev_rem):
        seq = seq[3-prev_rem:]
        rem = len(seq)%3
        return rem 

    def cds_extractor(self,eachdf):
        p53 = self._extract_all_cds(eachdf)
        dic = {}
        for var in p53:
            eachvar_dic = p53[var]
            eachvardf = self._make_each_variant_df(eachvar_dic)
            dic[var] = eachvardf

        return dic


class editable_pam():
    def codon_start_number_and_search_area(self,vardf):
        n = 0
        codon_num = 0
        while n < (vardf.shape[0]):
            remnant = vardf.loc[n,'Remnant_NT']
            trimmed = vardf.loc[n,'trimmed_CDS']

            no_codon = len(trimmed)/3

            if n == 0:
                vardf.loc[n,'CodonStart'] = 1
            else:
                vardf.loc[n,'CodonStart'] = codon_num+1

            if remnant != 0:
                codon_num += no_codon+1
            else:
                codon_num += no_codon

            n+=1
        return vardf
    def make_codon(self,cds,start,trimmed_pos):
        n = 0
        
        codon_num = start
        codon_idx = trimmed_pos

        lst = []
        while n < len(cds)-3+1:
            codon = cds[n:n+3]
            lst.append([codon_num,codon_idx,codon])  

            n += 3
            codon_num += 1
            codon_idx += 3
        
        return lst
    def search_fp_pam(self,pam,pam_window):
        n = 0
        lst = []
        while n < len(pam_window)-len(pam)+1:
            if pam_window[n:n+len(pam)] == pam:
                lst.append(n)
            else:
                pass
            n+= 1

        return lst
    def find_every_targetable_rp_codon(self,window,pam,cds,codon_start,cds_ext):
        trimmed_pos = cds_ext.find(cds)

        codon_lst = self.make_codon(cds,codon_start,trimmed_pos)

        lst = []
        for eachlst in codon_lst:
            ## eachlst = [codon_number, codon_index_in_cds_ext,codon]
            ## 158:179 ##164
            pam_window = cds_ext[eachlst[1]+window[0]-7:eachlst[1]+window[1]+len(pam)-5]


            pam_window_index = cds_ext.find(pam_window)

            fp_pamlst,rt_pamlst = general().pam_sequence(pam)

            if pam[0] == 'N':
                pam_window = pam_window[:-1]
                rt_pamlst = [i[:-1] for i in rt_pamlst]
                rt_pamlst = list(set(rt_pamlst))
                pam_window_index = pam_window_index
            else:
                pass

            ## PAM = rt(RCH)
            for eachpam in rt_pamlst:
                possible_pam_lst = self.search_fp_pam(eachpam,pam_window)
                if len(possible_pam_lst) != 0:
                    for eachpam_index in possible_pam_lst:
                        eachresult = []
                        eachpam_index = eachpam_index+pam_window_index
                        gx20 = cds_ext[len(pam)+eachpam_index:len(pam)+20+eachpam_index]
                        gx20 = general().rt_sequence(gx20)
                        real_pam = cds_ext[len(pam)+eachpam_index-len(pam):len(pam)+eachpam_index]
                        real_pam = general().rt_sequence(real_pam)
                        eachresult.extend(eachlst)
                        eachresult.extend([real_pam,gx20,'RP'])
                        lst.append(eachresult)
        if len(lst) == 0:
            return pd.DataFrame()
        else:
            return pd.DataFrame(lst,columns = ['codon_number','codon_index','codon','PAM','GX20','GX20_strand'])   
    def find_every_targetable_fp_codon(self,window,pam,cds,codon_start,cds_ext):

        trimmed_pos = cds_ext.find(cds)

        codon_lst = self.make_codon(cds,codon_start,trimmed_pos)

        lst = []
        for eachlst in codon_lst:
            ## eachlst = [codon_number, codon_index_in_cds_ext,codon]
            pam_window = cds_ext[eachlst[1]-window[1]+4:eachlst[1]-window[0]+6+len(pam)]
            pam_window_index = cds_ext.find(pam_window)

            fp_pamlst,rt_pamlst = general().pam_sequence(pam)

            if pam[0] == 'N':
                pam_window = pam_window[1:]
                fp_pamlst = [i[1:] for i in fp_pamlst]
                fp_pamlst = list(set(fp_pamlst))
                pam_window_index = pam_window_index+1
            else:
                pass
            ## PAM = RCH
            for eachpam in fp_pamlst:
                possible_pam_lst = self.search_fp_pam(eachpam,pam_window)
                ## possible_pam_lst = [] / [1] / [1,4]...

                if len(possible_pam_lst) != 0:
                    for eachpam_index in possible_pam_lst:                   
                        eachpam_index = eachpam_index+pam_window_index
                        gx20 = cds_ext[eachpam_index-1-20:eachpam_index-1]
                        real_pam = cds_ext[eachpam_index-1:eachpam_index-1+len(pam)]
                        eachresult = eachlst+[real_pam,gx20,'FP']
                        lst.append(eachresult)
                else: pass
            
        
        if len(lst) == 0:
            return pd.DataFrame()
        else:
            return pd.DataFrame(lst,columns = ['codon_number','codon_index','codon','PAM','GX20','GX20_strand'])
    def every_editable_pam(self,dic,eachdf,window,pam):
        lst = []
        var_num = list(dic.keys())[0]
                    
        ccds_id = eachdf.loc[var_num,'ccds_id']
        nc_acc = eachdf.loc[var_num,'nc_accession']
        chrom = eachdf.loc[var_num,'#chromosome']
        strand = eachdf.loc[var_num,'cds_strand']

        eachvar = dic[var_num]
        if eachvar.shape[0] == 0:
            pass
        else:            
            eachvar = self.codon_start_number_and_search_area(eachvar)
            if strand == '-':
                eachvar = eachvar.rename(columns = {'start':'end','end':'start'})            
            else:
                pass

            eachlst = []

            for idx in eachvar.index:
                cds_num = idx
                start = eachvar.loc[idx,'start']
                end = eachvar.loc[idx,'end']

    

                cds = eachvar.loc[idx,'trimmed_CDS']
                cds_start = eachvar.loc[idx,'CDS_original'].find(cds)

                codon_start = eachvar.loc[idx,'CodonStart']
                cds_ext = eachvar.loc[idx,'CDS_extended']

                targetable_codon_fp = self.find_every_targetable_fp_codon(window,pam,cds,codon_start,cds_ext)
                targetable_codon_rp = self.find_every_targetable_rp_codon(window,pam,cds,codon_start,cds_ext)

                targetable_codon_fprp = pd.concat([targetable_codon_fp,targetable_codon_rp])
                targetable_codon_fprp = targetable_codon_fprp.copy()
                targetable_codon_fprp['CDS_extended'] = cds_ext
                targetable_codon_fprp['nc_accession'] = nc_acc
                targetable_codon_fprp['chr'] = chrom
                targetable_codon_fprp['start'] = start
                targetable_codon_fprp['end'] = end
                targetable_codon_fprp['CDS_start'] = cds_start
                targetable_codon_fprp['Exon'] = cds_num
                eachlst.append(targetable_codon_fprp)


            targetable_codon = pd.concat(eachlst)
            targetable_codon = targetable_codon.copy()
            targetable_codon['ccds_id'] = ccds_id
            targetable_codon['strand'] = strand

            lst.append(targetable_codon) 
        if len(lst) == 0:
            return pd.DataFrame()
        else:
            return pd.concat(lst)

def make_every_single_subsitution(editable_pam_df):
    columns = list(editable_pam_df.columns)
    lst = []
    for idx in range(editable_pam_df.shape[0]):
        eachline = editable_pam_df.iloc[idx]
        codon = eachline['codon']
        snvlst = subsitution().single_snv_in_one_codon(codon)
        codon_index = eachline['codon_index']
        cds_ext = eachline['CDS_extended']


        for tup in snvlst:
            snv_pos = tup[0]
            snv = tup[1]
            
            if snv_pos == 0:
                intended_seq = cds_ext[codon_index-60:codon_index]+snv+cds_ext[codon_index+3:codon_index+3+58]
                wt_seq = cds_ext[codon_index-60:codon_index]+codon+cds_ext[codon_index+3:codon_index+3+58]
            elif snv_pos == 1:
                intended_seq = cds_ext[codon_index-59:codon_index]+snv+cds_ext[codon_index+3:codon_index+3+59]
                wt_seq = cds_ext[codon_index-59:codon_index]+codon+cds_ext[codon_index+3:codon_index+3+59]

            elif snv_pos == 2:
                intended_seq = cds_ext[codon_index-58:codon_index]+snv+cds_ext[codon_index+3:codon_index+3+60]
                wt_seq = cds_ext[codon_index-58:codon_index]+codon+cds_ext[codon_index+3:codon_index+3+60]

            else:
                intended_seq = 'FALSE'

            if intended_seq == 'FALSE':
                pass
            else:
                eachresult = []
                eachresult.extend(eachline)
                eachresult.extend([snv,snv_pos,wt_seq,intended_seq])
                lst.append(eachresult)

    columns.extend(['SNV','SNV_index','WT','Intended_Edit_Seq'])
    result = pd.DataFrame(lst,columns = columns)
    return result    

def extract_editable_info(gene,window,pam):
    df = pd.read_csv('EssentialData/CCDS.current.txt',sep='\t')
    df = df[df['cds_locations']!= '-']
    eachdf = df[df['gene']==gene]
    try:
        dic = extract_cds_from_all_variants().cds_extractor(eachdf)
        editable_pam_df = editable_pam().every_editable_pam(dic,eachdf,window,pam)
        final = make_every_single_subsitution(editable_pam_df)

        final = final.drop(columns = ['CDS_extended'])
        final = final.copy()
        final['gene'] = gene
        final['Cas9'] = pam
    
    except:
        final = pd.DataFrame()

    return final
# %%


def main(window,pam):
    genelst = pd.read_csv('Input/PE/input.txt',sep='\t')
    futs = []
    with ProcessPoolExecutor(max_workers=30) as executor:
        for gene in genelst['gene']:
            fut = executor.submit(extract_editable_info,gene,window,pam)
            futs.append(fut)
    
    merged = []
    for fut in futs:
        merged.append(fut.result())

    return pd.concat(merged)
    

# %%
if __name__ == "__main__":
    window = (1,16) ## input var
    pam = 'NRCH' ## input var
    o = 'PE_NRCH' ## input var
    output = 'Library/{}.txt'.format(o)
    df = main(window,pam)
    df.to_csv(output,sep='\t',index=None)

#%%
