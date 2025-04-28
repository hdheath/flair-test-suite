import sys
import matplotlib.pyplot as plt

"""
0 isoform
1 chrom
2 strand
3 length
4 exons
5 structural_category
6 associated_gene
7 associated_transcript
8 ref_length
9 ref_exons
10 diff_to_TSS
11 diff_to_TTS
12 diff_to_gene_TSS
13 diff_to_gene_TTS
14 subcategory
15 RTS_stage
16 all_canonical
17 min_sample_cov
18 min_cov
19 min_cov_pos
20 sd_cov
21 FL
22 n_indels
23 n_indels_junc
24 bite
25 iso_exp
26 gene_exp
27 ratio_exp
28 FSM_class
29 coding
30 ORF_length
31 CDS_length
32 CDS_start
33 CDS_end
34 CDS_genomic_start
35 CDS_genomic_end
36 predicted_NMD
37 perc_A_downstream_TTS
38 seq_A_downstream_TTS
39 dist_to_cage_peak
40 within_cage_peak
41 pos_cage_peak
42 dist_to_polya_site
43 within_polya_site
44 polyA_motif
45 polyA_dist
46 ORF_seq
47 TSS_genomic_coord
48 TTS_genomic_coord
49 experiment_id
50 entry_id
51 LRGASP_id
"""



manifestfile = sys.argv[1]
summaryprefix = sys.argv[2]

files = []
for line in open(manifestfile):
    files.append(line.rstrip())

isoclasses = ['FSM', 'ISM', 'NIC', 'NNC', 'Genic_Genomic', 'Genic_Intron', 'Antisense', 'Fusion','Intergenic']

t = set()
alllines = []
for file in files:
    tot, cagesup, polyasup = 0, 0, 0
    srtm, sntm = 0, 0
    isoclasscounts = {x: 0 for x in isoclasses}
    for line in open(file):
        line = line.rstrip().split('\t')
        if line[0] != 'isoform':
            tot += 1
            isoclass = line[5]
            isoclasscounts[isoclass] += 1
            hascagesup = line[40] == 'True' or (line[39] != 'NA' and abs(float(line[39])) <=50)
            haspolyasup = line[43] == 'True' or (line[42] != 'NA' and abs(float(line[42])) <=50)
            annottsssup = (line[10] != 'NA' and abs(float(line[10])) <=50) or (line[11] != 'NA' and abs(float(line[11])) <=50)
            annotttssup = (line[10] != 'NA' and abs(float(line[11])) <=50) or (line[11] != 'NA' and abs(float(line[13])) <=50)
            if isoclass in {'FSM', 'ISM'}:
                if (hascagesup or annottsssup) and (haspolyasup or annotttssup): srtm += 1
            elif (hascagesup or annottsssup) and (haspolyasup or annotttssup): sntm += 1

            if hascagesup: cagesup += 1
            if haspolyasup: polyasup += 1

    #         print(isoclass, line[39:44])#line[40], line[43])
    #         if tot > 100: break
    # break
    print(round(srtm/tot, 3), round(sntm/tot, 3), round((srtm+sntm)/tot, 3))
    outline = [file.split('_classification.txt')[0],
               round(srtm/tot, 3), round(sntm/tot, 3),
               round(cagesup/tot, 3), round(polyasup/tot, 3)] \
              + [isoclasscounts[x] for x in isoclasses] + [tot]
    alllines.append(outline)


out = open(summaryprefix + '.tsv', 'w')
out.write('\t'.join(['filename', 'fracSRTM', 'fracSNTM', "5'endsupport", "3'endsupport"] + isoclasses + ['total_isoforms']) + '\n')
for outline in alllines:
    out.write('\t'.join([str(x) for x in outline]) + '\n')
out.close()


fig, axs = plt.subplots(3)
fig.set_size_inches((6,10))
xvals = list(range(len(alllines)))

allsrtm, allsntm = [x[1] for x in alllines], [x[2] for x in alllines]
axs[0].bar(xvals, allsrtm,  label="SRTM")
axs[0].bar(xvals, allsntm, bottom=allsrtm,  label="SNTM")
axs[0].set_ylim((0,1))
axs[0].legend()
axs[0].set_xticks(xvals, ['' for x in alllines])
box = axs[0].get_position()
axs[0].set_position([box.x0 , box.y0+(box.height*(1/3)), box.width*0.8, box.height * (2/3)])
axs[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))


allcage, allpolya = [x[3] for x in alllines], [x[4] for x in alllines]
axs[1].bar([x-0.2 for x in xvals], allcage, width=0.3, label="5'endsupport")
axs[1].bar([x+0.2 for x in xvals], allpolya, width=0.3, label="3'endsupport")
axs[1].set_ylim((0,1))
axs[1].legend()
axs[1].set_xticks(xvals, ['' for x in alllines])
box = axs[1].get_position()
axs[1].set_position([box.x0 , box.y0+(box.height*(2/3)), box.width*0.8, box.height * (2/3)])
axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))

bottom = [0 for x in range(len(alllines))]
for i in range(len(isoclasses)):
    vals = [x[i+5] for x in alllines]
    axs[2].bar(xvals, vals, bottom=bottom, label=isoclasses[i])
    bottom = [bottom[x] + vals[x] for x in range(len(vals))]
axs[2].legend()
axs[2].set_xticks(xvals, [x[0].split('/')[-1] for x in alllines], rotation=90, ha='center')
box = axs[2].get_position()
axs[2].set_position([box.x0 , box.y0+(box.height*(1)), box.width*0.8, box.height*(2/3)])
axs[2].legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.savefig(summaryprefix + '.png', dpi=600)

