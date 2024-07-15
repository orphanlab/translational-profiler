import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Bio import SeqIO

def readDist(bamfile, verbose=True):
    
    aln = pysam.AlignmentFile(bamfile, 'rb')
    r={}
    for read in aln.fetch():
        try: 
            r[read.query_length] = r[read.query_length] + 1
        except KeyError:
            r[read.query_length] = 1
    if verbose:
        plt.bar(list(r.keys()), list(r.values()))
        print(r)

    return r


def read_prodigal(annfile, gene_type="CDS"):
    a = pd.read_table(annfile, names=['contig','source','type','start','stop','a','direction','b','tags'], sep="\t",
                    comment='#')
    
    # guess if gff
    if len(a.type.unique()) > 1:
        try:
            a = a.loc[a.type.isin(gene_type)]
        except TypeError:
            a = a.loc[a.type == gene_type]
    # indexes start at 0, stop is exclusive
    a.start = a.start - 1
    return a


def featureCounts(bamfile, annfile, paired=False, minlength=0, maxlength=np.inf, rpkm=False):

    if isinstance(annfile, str):
        annfile = read_prodigal(annfile)
    if isinstance(bamfile, str):
        aln = pysam.AlignmentFile(bamfile, 'rb')
    else:
        aln = bamfile

    counts=[]
    lens=[]
    for contig, start, stop, d in annfile[['contig','start','stop','direction']].itertuples(index=False):
        lens=lens + [stop-start]
        if d=='+':
            if paired:
                seen = [r.query_name for r in aln.fetch(contig, start, stop) 
                        if r.is_forward and r.is_mapped and r.is_read1 and r.is_paired
                        and minlength < r.query_alignment_length < maxlength]
                counts = counts + [len(seen)]

                counts[-1] = counts[-1] + len([r.reference_start for r in aln.fetch(contig, start, stop) 
                                               if r.query_name not in seen 
                                               and r.is_reverse and r.is_mapped and r.is_paired  
                                               and r.is_read2 and minlength < r.query_alignment_length < maxlength])
            else:
                counts = counts + [len([r.reference_start for r in aln.fetch(contig, start, stop) 
                                        if r.is_forward and r.is_mapped 
                                        and minlength < r.query_alignment_length < maxlength])]

        if d=='-':
            if paired:
                seen = [r.query_name for r in aln.fetch(contig, start, stop) 
                        if r.is_reverse and r.is_mapped and r.is_read1 and r.is_paired 
                        and minlength < r.query_alignment_length < maxlength]
                counts = counts + [len(seen)]

                counts[-1] = counts[-1] + len([r.reference_start for r in aln.fetch(contig, start, stop) 
                                        if r.query_name not in seen 
                                               and r.is_forward and r.is_mapped and r.is_paired and r.is_read2 
                                               and minlength < r.query_alignment_length < maxmaxlengthread])
            else:
                counts = counts + [len([r.reference_start for r in aln.fetch(contig, start, stop) 
                                        if r.is_reverse and r.is_mapped 
                                        and minlength < r.query_alignment_length < maxlength])]
                
    if rpkm:
        denom = sum([int(x.split('\t')[2]) for x in pysam.idxstats(bamfile).split('\n') if len(x) > 2])
        
        counts = np.array(counts) * (1/(np.array(lens)/1000)) * (1/(denom/1e6))  # rpkm
    
    annfile['counts'] = counts
    
    return annfile


def geneClear(start, stop, d, w, gmat):
    if d == '+':
        sclear = gmat[start-w-1:start,1].sum() == 0
        spclear = gmat[stop+1:stop+1+w,1].sum() == 0
    else:
        sclear = gmat[(start-w-1):start,2].sum() == 0
        spclear = gmat[(stop+1):(stop+1+w),2].sum() == 0

    bclear = sclear & spclear
    
    return bclear
    

def checkGeneSpacing(annfile, window, inverse=False):
    if isinstance(annfile, str):
        annfile = read_prodigal(annfile)
        
    gmat = np.stack((np.arange(annfile.stop.max()), np.zeros(annfile.stop.max()), 
                     np.zeros(annfile.stop.max())), axis=1)
    
    for s,sp,d in annfile[['start','stop','direction']].itertuples(index=False):
        if d=='+':
            gmat[s:sp,1] = 1
        else:
            gmat[s:sp,2] = 1
            
    keeps=[]
    for i,g in enumerate(annfile.itertuples(index=False)):
        keeps = keeps + [geneClear(g.start, g.stop, g.direction, w=window, gmat=gmat)]
        
    if inverse:
        keeps = [~k for k in keeps]
    return annfile[keeps]


def read_fasta(fasta, trim=False):
    #fasta_dict = {}
    seqio = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    fasta_dict = {k:str(v.seq) for k,v in seqio.items()}
    #[str(my_dict[k].seq)[0:5] for k,v in my_dict.items()]
    #     with open(fasta, 'r') as f:
    #         for l in f.readlines():
    #             if l.startswith('>'):
    #                 defline = l.rstrip().replace('>', '')
    #                 if trim:
    #                     defline = defline.split(' ')[0]
    #                 fasta_dict[defline] = ''
    #             else:
    #                 fasta_dict[defline] = fasta_dict[defline] + l.rstrip()
    return fasta_dict


def add_stop_start_to_gff(annfile, fasta_path, trim=False):
    
    if isinstance(fasta_path, str):
        fa = read_fasta(fasta_path, trim=trim)
    else: fa = fasta_path
    
    if isinstance(annfile, str):
        annfile = read_prodigal(annfile)
    
    aa_all = np.empty((annfile.shape[0],2), dtype=str)
    codons_all = []#np.empty((annfile.shape[0],2), dtype=str)
    for i, r in enumerate(annfile[['contig','start','stop','direction']].itertuples(index=False)):
        contig, start, stop, d = r
        if d=='+':
            cdn = np.array([fa[contig][i:i + 3] for i in [start, stop-3]])
            aa = np.array([codontab[c] for c in cdn])
            #print(np.array(cdn)[np.newaxis,:])
        else:
            cdn = np.array([fa[contig][i-3:i][::-1] for i in [stop, start + 3]])
            aa = np.array([anticodontab[c] for c in cdn])
            cdn = [''.join([complement[y] for y in x]) for x in cdn]
        
        aa_all[i,:] = aa #[[0,len(aa)-1]]
        codons_all = codons_all + [cdn]
        #print(np.array(cdn).shape)
        #codons_all[i,:] = (np.array(cdn)) #[[0,len(cdn)-1]]
        #         else:
        #             start_codon = fa[contig][stop-3:stop][::-1]
        #             start_codon = ''.join([complement[x] for x in start_codon])
        #             stop_codon = fa[contig][start:(start + 3)][::-1]
        #             stop_codon = ''.join([complement[x] for x in stop_codon])
    #print(codons_all)
    annfile[['start_aa','stop_aa']] = aa_all
    annfile[['start_cdn','stop_cdn']] = codons_all
    
    return annfile


def whereAreX(annfile, fasta_path, codon='O', trim=False):
    if isinstance(fasta_path, str):
        fa = read_fasta(fasta_path, trim=trim)
    else: fa = fasta_path
    if isinstance(annfile, str):
        annfile = read_prodigal(annfile)
    
    codons_of_interest = []
    for contig, start, stop, d in annfile[['contig','start','stop','direction']].itertuples(index=False):
        if d=='+':
            codons = np.array([codontab[fa[contig][i:i + 3]] for i in range(start,stop,3)])
        else:
            codons = np.array([anticodontab[fa[contig][i-3:i][::-1]] for i in range(stop,start,-3)])
        
        codons_of_interest = codons_of_interest + [np.where(np.isin(codons,np.array(codon)))[0]]
    codons_of_interest = [list(x) for x in codons_of_interest]

    annfile['cdn_pos'] = codons_of_interest
    return annfile


def translate(start,stop, contig, d, fasta_path):
    if isinstance(fasta_path, str):
        fa = read_fasta(fasta_path, trim=trim)
    else: fa = fasta_path

    if d=='+':
        codons = np.array([codontab[fa[contig][i:i + 3]] for i in range(start,stop,3)])
    else:
        codons = np.array([anticodontab[fa[contig][i-3:i][::-1]] for i in range(stop,start,-3)])

    return codons


def geneFromBam(bamfile, contig, start, stop, d, window, offset=0, footprint_size=24, fiveprime=True, 
                codon_resolved=False):
    
    p_out=np.zeros((stop-start)+2*window)

    if isinstance(bamfile, str):
        aln = pysam.AlignmentFile(bamfile, 'rb')
    else:
        aln = bamfile
    poff=1  # so zero ends up 50, not 49

    if not codon_resolved:
        p_out = np.sum(aln.count_coverage(contig, start-window, stop+window), axis=0)
        
    else:
        if d=='+':
            if fiveprime:
                c_f = list(np.unique([r.reference_start + offset for r in aln.fetch(contig, start-window, stop+window) 
                              if np.isin(r.query_length, footprint_size) and r.is_forward], return_counts=True))
            else: 
                c_f = list(np.unique([r.reference_end - offset for r in aln.fetch(contig, start-window, stop+window) 
                              if np.isin(r.query_length, footprint_size) and r.is_forward], return_counts=True)) 
            c_f[0] = c_f[0][np.where(c_f[0] < (stop + window))]
            c_f[1] = c_f[1][np.where(c_f[0] < (stop + window))]
        else:
            if fiveprime:
                c_f = list(np.unique([r.reference_end - offset for r in aln.fetch(contig, start-window+1, stop+window+1) 
                              if np.isin(r.query_length, footprint_size) and r.is_reverse], return_counts=True))
            else:
                c_f = list(np.unique([r.reference_start + offset for r in aln.fetch(contig, start-window+1, stop+window+1) 
                              if np.isin(r.query_length, footprint_size) and r.is_reverse], return_counts=True))
            c_f[0] = c_f[0][np.where(c_f[0] < (stop + window))]
            c_f[1] = c_f[1][np.where(c_f[0] < (stop + window))]
            poff = 0
            
        if len(c_f[0]) > 0:
            p_out[c_f[0] - (start - window) - poff] = c_f[1]

    return p_out

def transformCoverages(cov, d, window, rpm=False, meannorm=False, maxnorm=False):
    
    if d == '-':
        cov = cov[::-1]

    if rpm:
        if isinstance(rpm,int):
            #cov = cov * 1/rpm * 1e6  # rpm  #OLD
            cov = cov / (rpm/1e6)  # rpm  
        else:
            cov = cov * 1
            print("Not transforming, if you use RPM need to pass the number of reads in sample")
        #cov = cov / ((len(cov) - 2*window) / 1000)  # per k
    if meannorm:
        if any(cov[window:-window] > 0):
            cov = cov / np.sum(cov[window:-window], axis=0)
    if maxnorm:
        if any(cov[window:-window] > 0):
            cov = cov / np.max(cov, axis=0)
        
    if window:
        cov1 = cov[0:(3*window)]
        cov2 = cov[(-3*window):]
    else:
        cov1 = cov
        cov2 = cov
    return cov1, cov2


def metagene(bamfile, annfile, window, footprint_size=24, offset=0, num_missing=0, fasta_path=False, minrpm=0, geneposoffset=0,
             fiveprime=True, codon_resolved=False, meannorm=False, maxnorm=False, rpm=False, verbose=False):
    
    covs_start = np.zeros((window*3,1))
    covs_stop = np.zeros((window*3,1))
    
    if isinstance(annfile, str):
        annfile = read_prodigal(annfile)
    if isinstance(bamfile, str):
        aln=pysam.AlignmentFile(bamfile, 'rb')
    else:
        aln = bamfile
        
    if isinstance(fasta_path, str):
        fasta_path = read_fasta(fasta_path, trim=True)
        annfile = add_stop_start_to_gff(annfile=annfile, fasta_path=fasta_path, trim=True)
        #annfile = checkGeneSpacing(annfile, window=window)
        print('only M...* genes')
        annfile = annfile.loc[(annfile.start_cdn == 'ATG') & (annfile.stop_aa == '*')]  # subset first to make fast

    
    if minrpm > 0:
        print("Only showing reads from genes with rpkm > " + str(minrpm))
        annfile = featureCounts(annfile=annfile, bamfile=bamfile, paired=False, rpkm=True)
        annfile = annfile.loc[annfile.counts > minrpm]
    
    if rpm:
        denom = sum([int(x.split('\t')[2]) for x in pysam.idxstats(bamfile).split('\n') if len(x) > 2])
        if codon_resolved:
            rlen_counts = readDist(bamfile=bamfile, verbose=False) 
            denom = int(denom * (rlen_counts[footprint_size]/sum(rlen_counts.values())))

    else:
        denom = False

    for contig, start, stop, d in annfile[['contig','start','stop','direction']].itertuples(index=False):
        if start < window:
            print("Excluding a gene too close to end of contig")
            continue
        covs = geneFromBam(aln, contig=contig, start=start + geneposoffset, stop=stop + geneposoffset, window=window, d=d, fiveprime=fiveprime,
                             codon_resolved=codon_resolved, offset=offset, footprint_size=footprint_size)
        c_start, c_stop = transformCoverages(covs, d=d, window=window, 
                                             meannorm=meannorm, maxnorm=maxnorm, rpm=denom)
        if sum(c_start == 0) < (3*window-num_missing) and len(c_start == 3*window) and len(c_stop == 3*window):
            covs_start = np.concatenate((covs_start, c_start[:,np.newaxis]), axis=1)
            covs_stop = np.concatenate((covs_stop, c_stop[:,np.newaxis]), axis=1)

    if verbose:
        print(covs_start.shape)
    return covs_start[:,1:], covs_stop[:,1:]

def metaheatmap(bamfile, annfile, window, codon_resolved=True, offset=0, num_missing=0, footprint_sizes='auto',
             fiveprime=True, meannorm=False, maxnorm=False, rpm=False):
    
    h_start = np.zeros((window*3,1))
    h_stop = np.zeros((window*3,1))
    
    if isinstance(annfile, str):
        annfile = read_prodigal(annfile)
    if footprint_sizes == 'auto':
        hist = readDist(bamfile = bamfile, verbose=False)
        footprint_sizes = sorted([h for h in hist.keys() if hist[h] > 1000])
    
    for fp in footprint_sizes:
        c1, c2 = metagene(bamfile=bamfile, annfile=annfile, meannorm=meannorm, codon_resolved=codon_resolved, 
                          footprint_size=fp, rpm=rpm, fiveprime=fiveprime, maxnorm=maxnorm, offset=offset, 
                          window=window, num_missing=num_missing)
        
        c1 = np.sum(c1, axis=1)
        c2 = np.sum(c2, axis=1)
        
        h_start = np.concatenate((h_start, c1[:,np.newaxis]), axis=1)
        h_stop = np.concatenate((h_stop, c2[:,np.newaxis]), axis=1)
        
    return h_start[:,1:], h_stop[:,1:], footprint_sizes

            
def plot_metagene(covs_start, covs_stop, bg_start=None, bg_stop=None, window=None, codonlines=True, zero='middle', geneposoffset=0,
                 axis_units=None, colors=['black','#dddddddd'], bg_fill=True, norenorm=False):
    fig, ax = plt.subplots(2,1, figsize=(14,6))
    
    ax[0].plot(np.mean(covs_start, axis=1), c=colors[0])
    ax[1].plot(np.mean(covs_stop, axis=1), c=colors[0])
    
    _ = ax[0].text(x=0.025, y=0.95, s='n: ' + str(covs_start.shape[1]) + ' genes', 
                    ha='left', va='center', c=colors[0], transform=ax[0].transAxes)
    _ = ax[1].text(x=0.025, y=0.95, s='n: ' + str(covs_stop.shape[1]) + ' genes', 
                ha='left', va='center', c=colors[0], transform=ax[1].transAxes) 
    
    foreground = 'Footprints'
    background = 'Coverage'
    if axis_units:
        foreground = foreground + ' (' + axis_units + ')'
        background = background + ' (' + axis_units + ')'
        
    ax[0].set_ylabel(foreground)
    ax[1].set_ylabel(foreground)
    
    ax[0].set_xlabel('nt')
    ax[1].set_xlabel('nt')
    
    if bg_start is not None:        
        y2off = round((bg_start.max()/covs_start.max())/5) * 5
        if y2off == 0:
            y2off = round((covs_start.max()/bg_start.max())/5) * 5
            try:
                y2off = 1/y2off
            except ZeroDivisionError:
                y2off = 1
        if norenorm:
            y2off = 1
        
        if bg_fill:
            ax[0].fill_between(x=np.arange(bg_start.shape[0]),y1=0, y2=np.mean(bg_start, axis=1)/y2off, fc=colors[1])
            ax[1].fill_between(x=np.arange(bg_start.shape[0]),y1=0, y2=np.mean(bg_stop, axis=1)/y2off, fc=colors[1])
        else:
            ax[0].plot(np.mean(bg_start, axis=1)/y2off, c=colors[1])
            ax[1].plot(np.mean(bg_stop, axis=1)/y2off, c=colors[1])
            
        _ = ax[0].text(x=0.025, y=0.88, s='n: ' + str(bg_start.shape[1]) + ' genes', 
                       ha='left', va='center', c=colors[1], transform=ax[0].transAxes)
        _ = ax[1].text(x=0.025, y=0.88, s='n: ' + str(bg_stop.shape[1]) + ' genes', 
                       ha='left', va='center', c=colors[1], transform=ax[1].transAxes) 
            
        
        secax0 = ax[0].secondary_yaxis('right', functions=(lambda x: x*y2off, lambda x: x/y2off))
        secax1 = ax[1].secondary_yaxis('right', functions=(lambda x: x*y2off, lambda x: x/y2off))

        _ = secax0.set_ylabel(background, c=colors[1])
        _ = secax1.set_ylabel(background, c=colors[1])

    if window:
        initsign = '-'
        if window < geneposoffset:
            initsign = '+'
        points_x = [0, window, 2*window, 3*window - 1]
        lab_start = [initsign + str(-1*window + geneposoffset), 
                     "+"+str(0+ geneposoffset), '+'+str(1*window + geneposoffset), '+'+str(2*window-1 + geneposoffset)]
        lab_stop = [str(-2*window + geneposoffset), str(-1*window + geneposoffset), str(0+ geneposoffset),
                    '+'+str(1*window - 1 + geneposoffset)]
        _ = ax[0].set_xticks(points_x, lab_start) 
        _ = ax[1].set_xticks(points_x, lab_stop) 

        if codonlines:
            [ax[0].axvline(y, c='#bbbbbb', linewidth=0.6, zorder=0) for y in 
             np.where(np.arange(0,3*window) - window + 0 + geneposoffset >=0)[0][::3]]#range(3,2*window,3)]
            [ax[1].axvline(y, c='#bbbbbb', linewidth=0.6, zorder=0) for y in range(1,2*window,3)]

        if geneposoffset == 0:
            ax[0].axvline(window, c='red', zorder=0)
            ax[1].axvline(2*window, c='red', zorder=0)


def plot_heatmap(covs, footprints_used, maxnorm=False, window=False, isForward=True):
    fig, ax = plt.subplots(1,1, figsize=(14,6))

    if maxnorm:
        covs = covs/np.max(covs, axis=0)
    hmap = ax.imshow(np.transpose(covs)[::-1], cmap='cividis')
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1%", pad=0.05)
    fig.colorbar(hmap, cax=cax)
    
    _ = ax.set_yticks(range(0,len(footprints_used),2), footprints_used[::-2])
    
    if window:
        points_x = [0, window, 2*window, 3*window - 1]
        lab_start = [str(-1*window), '0', '+'+str(1*window), '+'+str(2*window-1)]
        lab_stop = [str(-2*window), str(-1*window), '0', '+'+str(1*window - 1)]
        
        if isForward:
            _ = ax.set_xticks(points_x, lab_start)
        else:
            _ = ax.set_xticks(points_x, lab_stop)

    
    plt.show()


def fullgene(query, bamfile, annfile, footprint_size=24, window=50, offset=0, rpm=False, codon_resolved=True):
    if isinstance(annfile, str):
        annfile=read_prodigal(annfile)
    if rpm:
        fr = int(pysam.idxstats(bamfile).split('\t')[2])
    else:
        fr = rpm
    
    startstop=[]
    directions=[]
    for contig,start,stop,d,t in annfile.loc[annfile.tags.str.contains(query)][
        ['contig','start','stop','direction','tags']].itertuples(index=False):
        startstop = startstop + [start, stop]
        directions = directions + [d]
        
    minstart = min(startstop)
    maxstart = max(startstop)
    
        
    covs = geneFromBam(bamfile, contig, minstart, maxstart, d, window=window, 
                   footprint_size=footprint_size, offset=offset, codon_resolved=codon_resolved)
    
    if codon_resolved and len(set(directions)) > 1:
        print("Multipe directions in window! \nShowing orientation as " + d + " but COUNTING STRANDLESS!")
        if d == "+":
            d2 = "-"
        else:
            d2 = "+"
            scovs = geneFromBam(bamfile, contig, minstart, maxstart, d2, window=window, 
                                      footprint_size=footprint_size, offset=offset, codon_resolved=codon_resolved)
            covs = covs + scovs

    else:
        print('Direction of region: ' + d)
    
    # flip genes & rpm normalize, fr gives  the number of total mapped reads, w=None means no window: whole thing
    covs, _ = transformCoverages(covs, d, window=None, rpm=fr, meannorm=False, maxnorm=False)
    
    ## flip gene coords
    if d == "-":
        startstop = np.array(startstop[::-1])
        startstop = startstop - startstop[0]
        startstop = startstop * -1
    else:
        startstop = np.array(startstop)
        startstop = startstop - startstop[0]
        
    startstop = startstop + window
        
    return covs, startstop
    

def plot_genes(query, bamfile, annfile, footprint_size, txfile=None, window=50, offset=0, rpm=False, codon_resolved=True,
              cols=['red','#919100','green','blue','black'], axis_units=False, colors=['black','#dddddd'], bg_fill=True,
              norenorm=False, y1lab='Footprints', y2lab='Coverage', genelabs=False, stickmode=True, 
               residues=False, fasta=False, moving_average=0):
    
    if isinstance(annfile, str):
        annfile = read_prodigal(annfile) 

    covs, startstop = fullgene(query=query, bamfile=bamfile, annfile=annfile, footprint_size=footprint_size, 
                               window=window, offset=offset, rpm=rpm, codon_resolved=codon_resolved)
    
    # PLOT
    fig, ax = plt.subplots(1,1, figsize=(14,6))
    
    if len(cols) < (len(startstop)/2):
        print("Need this many colors: " + str(int(len(startstop)/2)))
        cols = np.tile(['#00004a', '#aacdef'], int(len(startstop)/2))
        
    _ = [ax.axvspan(ymin=0.01, ymax=0.04, 
                    xmin=startstop[i*2],# + w, 
                    xmax=startstop[(i*2)+1],# + w,
                   color=cols[i]) 
         for i in range(int(len(startstop)/2))]
    
    if genelabs:
        _ = [ax.text(x = (startstop[i*2] + startstop[(i*2)+1]) / 2,
                     y = 0.025, s=genelabs[ii], c='white', ha='center', va='center', transform=ax.get_xaxis_transform()) 
             for ii,i in enumerate(range(int(len(startstop)/2)))]
    
    foreground = y1lab
    background = y2lab
    if axis_units:
        foreground = foreground + ' (' + axis_units + ')'
        background = background + ' (' + axis_units + ')'
    _ = ax.set_ylabel(foreground)
    _ = ax.set_xlabel('nt')
    
    # if transcriptome for background
    if txfile:
        tx_covs, startstop = fullgene(query=query, bamfile=txfile, annfile=annfile, footprint_size=footprint_size, 
                                      window=window, offset=offset, rpm=rpm, codon_resolved=not bg_fill)
        
        # figure out offset between pileup and footprint-style cov
        y2off = round((tx_covs.max()/covs.max())/5) * 5
        if y2off == 0:
            y2off = round((covs.max()/tx_covs.max())/5) * 5
            try:
                y2off = 1/y2off
            except ZeroDivisionError:
                y2off = 1
        if norenorm:
            y2off=1
        if bg_fill:
            ax.fill_between(x=np.arange(len(tx_covs)),y1=0, y2=tx_covs/y2off, fc=colors[1])
        if stickmode:
            _ = [ax.plot((i,i), (0,x), c=colors[1], linewidth=1) for i,x in enumerate(tx_covs/y2off)]
        else:
            _ = ax.plot(tx_covs/y2off, c=colors[1], linewidth=0.8)
        
        secax = ax.secondary_yaxis('right', functions=(lambda x: x*y2off, lambda x: x/y2off))
        _ = secax.set_ylabel(background, c=colors[1])
    
    if moving_average > 0:
        print(covs.shape)
        ma_covs = np.lib.stride_tricks.sliding_window_view(covs,moving_average).mean(axis=1)
        print(ma_covs.shape)
        ma_covs = np.concatenate((np.repeat(np.nan, int((moving_average-1)/2)), ma_covs, np.repeat(np.nan, int((moving_average-1)/2))))
        ax.plot(ma_covs, c=colors[0], linewidth=1.2)
        
    if stickmode and codon_resolved:
        _ = [ax.plot((i,i), (0,x), c=colors[0], linewidth=0.8) for i,x in enumerate(covs)]
    else:
        _ = ax.plot(covs, c=colors[0], linewidth=0.8)
    
    if residues:
        if not fasta:
            print('missing "fasta=PATH" argument, stopping')
            return
        queryres = whereAreX(annfile.loc[annfile.tags.str.contains(query)], fasta, residues, trim=True)
        if queryres.direction.values[0] == '-':
            cdn_pos = queryres.cdn_pos[::-1]
        else:
            cdn_pos = queryres.cdn_pos
        
        for i,g in enumerate(cdn_pos):
            [ax.axvspan((3*x) + startstop[i*2], (3*x) + startstop[i*2], 0.05, 1, color='#039aDa69') for x in g if len(g) > 0]


def metacodon(bamfile, annfile, fasta_path, window=50, defline_trim=True, codon='M', query='', minrpm=0, rpm=True, fiveprime=True, 
    codon_resolved=False, offset=0, footprint_size=range(18,30), meannorm=False, maxnorm=False, normalgenes=True):

    if isinstance(fasta_path, str):
        fasta_path = read_fasta(fasta_path)
        
    ann = add_stop_start_to_gff(annfile=annfile, fasta_path=fasta_path, trim=defline_trim)
    ann = whereAreX(annfile=ann, fasta_path=fasta_path, codon=codon, trim=defline_trim)
    ann = ann.loc[ann.cdn_pos.str.len() >= 1]
    
    if codon == 'M' or codon == '*':
        ann = checkGeneSpacing(ann, window=window)
    elif not codon == 'O':
        ann = ann.loc[(abs(ann.stop-ann.start) - ann.cdn_pos.values[0][0]*3) >= 1*window]
    
    if normalgenes and not codon == 'O':
        print('only M...* genes')
        ann = ann.loc[(ann.tags.str.contains(query)) & (ann.start_cdn == 'ATG') & (ann.stop_aa == '*')]  # subset first to make fast
    else:
        ann = ann.loc[ann.tags.str.contains(query)]
    
    if minrpm > 0:
        print("Only showing reads with rpm > " + str(minrpm))
        ann = featureCounts(annfile=ann, bamfile=bamfile, paired=False, rpkm=True)
        ann = ann.loc[ann.counts > minrpm]

    if rpm:
        rpm = sum([int(x.split('\t')[2]) for x in pysam.idxstats(bamfile).split('\n') if len(x) > 2])

    covs = np.zeros((window*2,1))
    for contig, start, stop, cdn, d in ann[['contig','start','stop','cdn_pos','direction']].itertuples(index=False):
        if d == '+':
            cdnpos = (cdn[0] + 0)*3 + start
        else:
            cdnpos = stop - (cdn[0]+0)*3
        covs_gene = geneFromBam(bamfile, contig, cdnpos, cdnpos, d, window, 
                                offset=offset, footprint_size=footprint_size, fiveprime=fiveprime, codon_resolved=codon_resolved)

        covs_gene,_ = transformCoverages(covs_gene, d, window, rpm=rpm, meannorm=meannorm, maxnorm=maxnorm)

        covs = np.concatenate((covs, covs_gene[:,np.newaxis]), axis=1)

    return covs[:,1:]


def plot_metacodon(covs_start, bg_start=None, window=None, codonlines=True, autoaxis=True, traces=True,
                 axis_units=None, colors=['black','#dddddddd'], bg_fill=True, indiv_genes=True, reasonablemode=True):
    
    panels = 1
    toomany = covs_start.shape[1]
    
    if (toomany > 500) and reasonablemode and traces:
        print('first 500 only')
        covs_start[0:500,:]
        
    if not indiv_genes:
        toomany=1
    elif toomany > 15:
        toomany = 16
    fig, ax = plt.subplots(toomany,1, figsize=(14, 3*toomany))
    
    if toomany == 1:
        ax = [ax]
    
    ax[0].plot(np.mean(covs_start, axis=1), c=colors[0])
    ax[0].text(0.01, 0.92, 'N=' + str(covs_start.shape[1]) + ' genes', horizontalalignment='left', verticalalignment='center', 
               transform=ax[0].transAxes)
    
    if traces:
        _ = [ax[1].plot(x, c=colors[0], alpha=0.2) for x in covs_start.transpose()]
    
    if indiv_genes:
        _ = [ax[i+2].plot(covs_start[:,i], c=colors[0]) for i in range(toomany)]
    
    foreground = 'Footprints'
    background = 'Footprints'
    if axis_units:
        foreground = foreground + ' (' + axis_units + ')'
        background = background + ' (' + axis_units + ')'
        
    ax[0].set_ylabel(foreground, c=colors[0])
    if traces:
        ax[1].set_ylabel(foreground, c=colors[0])
    
    if bg_start is not None:
        y2off=1
        if autoaxis:
            y2off = round((bg_start.max()/covs_start.max())/5) * 5
            if y2off == 0:
                y2off = round((covs_start.max()/bg_start.max())/5) * 5
                try:
                    y2off = 1/y2off
                except ZeroDivisionError:
                    y2off = 1
        
        if bg_fill:
            ax[0].fill_between(x=np.arange(bg_start.shape[0]),y1=0, y2=np.mean(bg_start, axis=1)/y2off, fc=colors[1])
            if traces:
                ax[1].fill_between(x=np.arange(bg_start.shape[0]),y1=0, y2=np.mean(bg_stop, axis=1)/y2off, fc=colors[1])
        else:
            ax[0].plot(np.mean(bg_start, axis=1)/y2off, c=colors[1])
            if traces:
                _ = [ax[1].plot(x, c=colors[1], alpha=0.2) for x in bg_start.transpose()]
            if indiv_genes:
                _ = [ax[i+2].plot(bg_start[:,i], c=colors[1]) for i in range(toomany)]
        
        if indiv_genes:
            _ = [ax[i+2].text(x=0.05, y=0.9, s=str(i), ha='center', va='center', transform=ax[i+2].transAxes) 
                 for i in range(toomany)]

        
        secax0 = ax[0].secondary_yaxis('right', functions=(lambda x: x*y2off, lambda x: x/y2off))
        _ = secax0.set_ylabel(background, c=colors[1])
        
        if traces:
            secax1 = ax[1].secondary_yaxis('right')
            _ = secax1.set_ylabel(background, c=colors[1])
        if indiv_genes:
            _ = [ax[i+2].secondary_yaxis('right') for i in range(toomany)]


        

    if window:
        points_x = [0, window, 2*window - 1]
        lab_start = [str(-1*window), '0', '+'+str(1*window-1)]
        
        _ = [l.set_xticks(points_x, lab_start) for l in ax]
        _ = [l.axvline(window, c='red', zorder=0) for l in ax]
        
        if codonlines:
            [ax[0].axvline(window+y, c='#bbbbbb', linewidth=0.6, zorder=0) for y in range(3,2*window,3)]
            if traces:
                [ax[1].axvline(y, c='#bbbbbb', linewidth=0.6, zorder=0) for y in range(1,2*window,3)]
            
        
def codonhalves(bamfile, annfile, fasta_path, defline_trim=True, codon='M', query='', minrpm=0, rpm=True, fiveprime=True, 
    codon_resolved=False, offset=0, footprint_size=range(18,30), normalgenes=True):

    if isinstance(fasta_path, str):
        fasta_path = read_fasta(fasta_path)
        
    ann = add_stop_start_to_gff(annfile=annfile, fasta_path=fasta_path, trim=defline_trim)
    ann = whereAreX(annfile=ann, fasta_path=fasta_path, codon=codon, trim=defline_trim)
    ann = ann.loc[ann.cdn_pos.str.len() >= 1]
    
    if normalgenes:# and not codon == 'O':
        print('only M...* genes')
        ann = ann.loc[(ann.tags.str.contains(query)) & (ann.start_cdn == 'ATG')] # subset first to make fast
    else:
        ann = ann.loc[ann.tags.str.contains(query)]
        
    
    if minrpm > 0:
        print("Only showing reads with rpm > " + str(minrpm))
        ann = featureCounts(annfile=ann, bamfile=bamfile, paired=False, rpkm=True)
        ann = ann.loc[ann.counts > minrpm]

    if rpm:
        totalreads = sum([int(x.split('\t')[2]) for x in pysam.idxstats(bamfile).split('\n') if len(x) > 2])
    
    isstop=0
    if codon == '*':
        isstop=1
        ann = checkGeneSpacing(ann, window=350)

    covs_before = []
    covs_after = []
    for contig, start, stop, cdn, d in ann[['contig','start','stop','cdn_pos','direction']].itertuples(index=False):
        if len(cdn) > 0 and cdn[-1] > 83:
            cdnloc = cdn[max(0,min(np.where(np.array(cdn)*3 > 250)[0]))] + (isstop * 1)
        else:
            cdnloc = cdn[0] + (isstop * 1)
        if d == '+':
            cdnpos = (cdnloc + 0)*3 + start
            covs_b = geneFromBam(bamfile, contig, start, cdnpos, d, window=0, 
                        offset=offset, footprint_size=footprint_size, fiveprime=fiveprime, codon_resolved=codon_resolved)
            covs_a = geneFromBam(bamfile, contig, cdnpos, stop+(250*isstop), d, window=0, 
                        offset=offset, footprint_size=footprint_size, fiveprime=fiveprime, codon_resolved=codon_resolved)
        else:
            cdnpos = stop - (cdnloc+0)*3
            covs_b = geneFromBam(bamfile, contig, cdnpos, stop, d, window=0, 
                        offset=offset, footprint_size=footprint_size, fiveprime=fiveprime, codon_resolved=codon_resolved)
            covs_a = geneFromBam(bamfile, contig, start-(250*isstop), cdnpos, d, window=0, 
                        offset=offset, footprint_size=footprint_size, fiveprime=fiveprime, codon_resolved=codon_resolved)
        
        if rpm:
            covs_b = (covs_b * (totalreads/1e6))#/len(covs_b)  # rpkm
            covs_a = (covs_a * (totalreads/1e6))#/len(covs_a)  # rpkm
        
        covs_before = covs_before + [covs_b]
        covs_after = covs_after + [covs_a]

    return covs_before, covs_after

def codonslice(x, y, window=250):
    fullgenes = np.zeros((2*window,1))
    
    for i,b in enumerate(x):    
        a = y[i]
        if np.sum(b) == 0:
            continue
        
        if len(b) < window or len(a) < window:
            continue
        
        fullgene=np.concatenate((b[-window:],a[:window]))
        if np.nansum(fullgene) == 0:
            continue
        
        fullgene=np.cumsum(fullgene)/np.nansum(fullgene)
        before50=fullgene[:window]
        after50=fullgene[window:]
        fullgenes = np.concatenate((fullgenes, np.concatenate((before50, after50))[:, np.newaxis]), axis=1)
        
    fullgenes=fullgenes[:,1:]
    return fullgenes


def plotPrePostAA(residue_dict, window=100, 
                  colors={'K':'green','W':'brown','P':'purple','O':'red', 'L':'gray', 'S':'orange', 'F':'pink','R':'#7a1a00','*':'#303030'},
                 n_best=10000, traces=False):
    fig, ax = plt.subplots(1,1)

    txty = 0.95
    for k,v in residue_dict.items():
        if v.shape[1] > n_best:
            #diffs = np.diff(v)
            nonzeros = np.sum((v[:window,:] < 0.05) + (v[:window,:] > 0.95), axis=0)
            nz, nz_counts = np.unique(nonzeros, return_counts=True)
            nonzeros_min = nz[np.where(np.cumsum(nz_counts) >= n_best)[0][0]]
            nonzeros_id = np.where(nonzeros <= nonzeros_min)[0]
            
            v_best = v[:, nonzeros_id]
        else:
            v_best = v
            
        m=v_best.mean(axis=1)
        sd=v_best.std(axis=1)

        before=m[:window]
        after=m[window:]
        
        slopeb=round(np.polyfit(x=np.arange(window-20),y=before[20:], deg=1)[0],6)
        slopea=round(np.polyfit(x=np.arange(window-20),y=after[:-20], deg=1)[0],6)

        _ = ax.axvline(window, c='black', zorder=0, linewidth=0.7)
        ax.fill_between(y1=m-sd, y2=m+sd, x=np.arange(len(m)), fc=colors[k], alpha=0.1)
        if traces:
            # pick random subset of 20 to make traces
            [ax.plot(x, c=colors[k], alpha=0.1, linewidth=0.5) for x in v_best.transpose()[np.random.randint(0,v_best.shape[1], 20)]]
        ax.plot(m, c=colors[k])

        #_ = ax.axvline(window-25, c='gray', linestyle='--')
        #_ = ax.axvline(window+25, c='gray', linestyle='--')
        _ = ax.text(x=0.025, y=txty, s=k + ' (n: ' + str(v_best.shape[1]) + ')', 
                    ha='left', va='center', c=colors[k], transform=ax.transAxes) 
        txty = txty - 0.05
        
        points_x = [0, int(window/2), window, int(1.5*window), 2*window - 1]
        lab_start = [str(-1*window), '-' + str(int(window/2)), '0', '+'+str(int(.5*window)), '+'+str(1*window-1)]
        
        _ = ax.set_xticks(points_x, lab_start)
        
        _ = ax.set_xlabel('nt')
        _ = ax.set_ylabel('Cumulative coverage (normalized)')
        
        
codontab = {
    'TCA': 'S',
    'TCC': 'S',
    'TCG': 'S',
    'TCT': 'S',
    'TTC': 'F',
    'TTT': 'F',
    'TTA': 'L',
    'TTG': 'L',
    'TAC': 'Y',
    'TAT': 'Y',
    'TAA': '*', # STOP
    'TAG': 'O', # PYL / STOP
    'TGC': 'C',
    'TGT': 'C',
    'TGA': '*', # STOP
    'TGG': 'W',
    'CTA': 'L',
    'CTC': 'L',
    'CTG': 'L',
    'CTT': 'L',
    'CCA': 'P',
    'CCC': 'P',
    'CCG': 'P',
    'CCT': 'P',
    'CAC': 'H',
    'CAT': 'H',
    'CAA': 'Q',
    'CAG': 'Q',
    'CGA': 'R',
    'CGC': 'R',
    'CGG': 'R',
    'CGT': 'R',
    'ATA': 'I',
    'ATC': 'I',
    'ATT': 'I',
    'ATG': 'M',
    'ACA': 'T',
    'ACC': 'T',
    'ACG': 'T',
    'ACT': 'T',
    'AAC': 'N',
    'AAT': 'N',
    'AAA': 'K',
    'AAG': 'K',
    'AGC': 'S',
    'AGT': 'S',
    'AGA': 'R',
    'AGG': 'R',
    'GTA': 'V',
    'GTC': 'V',
    'GTG': 'V',
    'GTT': 'V',
    'GCA': 'A',
    'GCC': 'A',
    'GCG': 'A',
    'GCT': 'A',
    'GAC': 'D',
    'GAT': 'D',
    'GAA': 'E',
    'GAG': 'E',
    'GGA': 'G',
    'GGC': 'G',
    'GGG': 'G',
    'GGT': 'G' 
}

anticodontab={'AGT': 'S', 'AGG': 'S', 'AGC': 'S', 'AGA': 'S', 'AAG': 'F', 'AAA': 'F', 'AAT': 'L', 'AAC': 'L', 'ATG': 'Y', 'ATA': 'Y', 'ATT': '*', 'ATC': 'O', 'ACG': 'C', 'ACA': 'C', 'ACT': '*', 'ACC': 'W', 'GAT': 'L', 'GAG': 'L', 'GAC': 'L', 'GAA': 'L', 'GGT': 'P', 'GGG': 'P', 'GGC': 'P', 'GGA': 'P', 'GTG': 'H', 'GTA': 'H', 'GTT': 'Q', 'GTC': 'Q', 'GCT': 'R', 'GCG': 'R', 'GCC': 'R', 'GCA': 'R', 'TAT': 'I', 'TAG': 'I', 'TAA': 'I', 'TAC': 'M', 'TGT': 'T', 'TGG': 'T', 'TGC': 'T', 'TGA': 'T', 'TTG': 'N', 'TTA': 'N', 'TTT': 'K', 'TTC': 'K', 'TCG': 'S', 'TCA': 'S', 'TCT': 'R', 'TCC': 'R', 'CAT': 'V', 'CAG': 'V', 'CAC': 'V', 'CAA': 'V', 'CGT': 'A', 'CGG': 'A', 'CGC': 'A', 'CGA': 'A', 'CTG': 'D', 'CTA': 'D', 'CTT': 'E', 'CTC': 'E', 'CCT': 'G', 'CCG': 'G', 'CCC': 'G', 'CCA': 'G'}

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

