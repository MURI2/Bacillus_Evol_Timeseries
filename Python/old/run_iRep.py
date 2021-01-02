from __future__ import division
import os, csv
import pandas as pd
import  matplotlib.pyplot as plt
from scipy import stats

mydir = os.path.expanduser('~/GitHub/Task2/PoolPopSeq/')

def get_iRep(day = 'D300'):
    bwa_fna_dir = mydir + 'data/bwa_fna/' + day + '/'
    for filename in os.listdir(bwa_fna_dir):
        if filename.endswith('_mapped_sort_NOdup_sort_merged_sort.sam'):
            sam_path = os.path.join(bwa_fna_dir, filename)
            strain = filename[9]
            if strain == 'B' or strain == 'S':
                fna_path = mydir + "data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.fna"
            elif strain == 'C':
                fna_path = mydir + "data/reference_assemblies_task2/Caulobacter_crescentus_NA1000/GCA_000022005.1_ASM2200v1_genomic.fna"
            elif strain == 'D':
                fna_path = mydir + "data/reference_assemblies_task2/Deinococcus_radiodurans_BAA816/GCA_000008565.1_ASM856v1_genomic.fna"
            elif strain == 'F':
                fna_path = mydir + "data/reference_assemblies_task2/Pedobacter_sp_KBS0701/G-Chr1.fna"
            elif strain == 'J':
                fna_path = mydir + "data/reference_assemblies_task2/Janthinobacterium_sp_KBS0711/KBS0711_2015_SoilGenomes_Annotate/G-Chr1.fna"
            elif strain == 'P':
                fna_path = mydir + "data/reference_assemblies_task2/Pseudomonas_sp_KBS0710/G-Chr1.fna"
            else:
                print("strain not recognized")
                continue
            iRep_path = mydir + 'data/iRep/' + day + '/'
            if not os.path.exists(iRep_path):
                os.makedirs(iRep_path)
            sample = "_".join(filename.split("_", 2)[:2])
            iRep_sample_path =  mydir + 'data/iRep/' + day + '/' + sample
            if not os.path.exists(iRep_sample_path):
                os.makedirs(iRep_sample_path)

            # remove regions of very high coverage (> 3 standard deviations with poisson dist)
            #samtools depth *bamfile*  |  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}'
            #cov_mean_cmnd = "samtools depth " + sam_path + " |  awk '{sum+=$3} END {print sum/NR}'"
            #cov_std_cmnd = "samtools depth " + sam_path + " |  awk '{sumsq+=$3*$3} END { print sqrt(sumsq/NR - (sum/NR)**2)}'"
            #cov_mean = float(os.popen(cov_mean_cmnd).read().strip())
            #cov_std = float(os.popen(cov_std_cmnd).read().strip())
            #cutoff = str(cov_mean + (2 * cov_std))
            #cutoff = str(2 * cov_mean)
            command = "iRep -t 6 -f " + fna_path + " -s " + sam_path + " -o " + iRep_sample_path + '/' + sample + '.iRep' + ' -ff'
            print(sample)
            os.system(command)


def clean_iRep():
    days = ['D100', 'D200', 'D300']
    OUT = open(mydir + 'data/iRep/iRep_params.txt', 'w')
    writer = csv.writer(OUT, delimiter= '\t')
    cols = ['Line', 'Strain', 'Treatment', 'Rep', 'Day', 'iRep', \
            'iRep_unfiltered', 'iRep_noGC', 'r2', 'Mean_cov', 'Windows_perc', \
            'Frags_mbp', 'GC_bias', 'GC_r2']
    writer.writerow(cols)
    for day in days:
        iRep_sample_path =  mydir + 'data/iRep/' + day + '/'
        for subdir, dirs, files in os.walk(iRep_sample_path):
            for f in files:
                if f.endswith('.iRep.tsv'):
                    IN = os.path.join(subdir, f)
                    line = subdir.split('/')[-1]
                    data = [line, line[9], line[8], line[10], line[-3:]]
                    for line in open(IN, 'r'):
                        line = line.strip()
                        if len(line) == 1:
                            continue
                        #print(len(line))
                        #print(line)
                        line = line.split('\t')
                        if '.fna' in line[0]:
                            data.append(line[1])
                    writer.writerow(data)
    OUT.close()


def plot_iRep():
    df = pd.read_csv(mydir + 'data/iRep/iRep_params.txt', sep = '\t')
    df = df[(df.Line != 'Sample_L2B3-100')]
    B = df[(df.Strain == 'B') & (df.Day == '100')]
    S = df[(df.Strain == 'S') & (df.Day == '100')]
    B_x = B.Treatment.values
    B_y = B.iRep_unfiltered.values
    S_x = S.Treatment.values
    S_y = S.iRep_unfiltered.values

    def setBoxColors(bp):
        plt.setp(bp['boxes'][0], color='blue')
        plt.setp(bp['caps'][0], color='blue')
        plt.setp(bp['caps'][1], color='blue')
        plt.setp(bp['whiskers'][0], color='blue')
        plt.setp(bp['whiskers'][1], color='blue')
        plt.setp(bp['fliers'][0], color='blue')

        plt.setp(bp['medians'][0], color='blue')

        plt.setp(bp['boxes'][1], color='red')
        plt.setp(bp['caps'][2], color='red')
        plt.setp(bp['caps'][3], color='red')
        plt.setp(bp['whiskers'][2], color='red')
        plt.setp(bp['whiskers'][3], color='red')
        plt.setp(bp['fliers'][1], color='red')
        #plt.setp(bp['fliers'][2], color='red')
        #plt.setp(bp['fliers'][3], color='red')
        plt.setp(bp['medians'][1], color='red')

    # Some fake data to plot
    B_0 = B[(B.Treatment == 0)]
    B_1 = B[(B.Treatment == 1)]
    B_2 = B[(B.Treatment == 2)]
    S_0 = S[(S.Treatment == 0)]
    S_1 = S[(S.Treatment == 1)]
    S_2 = S[(S.Treatment == 2)]

    zero = [B_0.iRep_unfiltered.tolist(), S_0.iRep_unfiltered.tolist()]
    one = [B_1.iRep_unfiltered.values, S_1.iRep_unfiltered.values]
    two = [B_2.iRep_unfiltered.values, S_2.iRep_unfiltered.values]

    fig = plt.figure()
    ax = plt.axes()
    plt.hold(True)

    # first boxplot pair
    bp = plt.boxplot(zero, positions = [1, 2], widths = 0.6)
    setBoxColors(bp)

    # second boxplot pair
    bp = plt.boxplot(one, positions = [4, 5], widths = 0.6)
    setBoxColors(bp)

    # thrid boxplot pair
    bp = plt.boxplot(two, positions = [7, 8], widths = 0.6)
    setBoxColors(bp)

    # set axes limits and labels
    plt.xlim(0,9)
    plt.ylim(1,2.2)
    ax.set_xticklabels(['1-day', '10-day', '100-day'])
    ax.set_xticks([1.5, 4.5, 7.5])
    ax.set_ylabel('Index of replication', fontsize = 20)

    # draw temporary red and blue lines and use them to create a legend
    hB, = plt.plot([1,1],'b-')
    hR, = plt.plot([1,1],'r-')
    plt.legend((hB, hR),('Bacillus', 'knockout'))
    hB.set_visible(False)
    hR.set_visible(False)
    fig_name = mydir + 'figs/S_B_iRep_box.png'
    plt.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)

    plt.close()

def plot_iRep_time():
    df = pd.read_csv(mydir + 'data/iRep/iRep_params.txt', sep = '\t')
    df = df[(df.Line != 'Sample_L2B3-100')]
    B_1 = df[(df.Strain == 'B') & (df.Treatment == 0)]
    B_10 = df[(df.Strain == 'B') & (df.Treatment == 1)]
    B_100 = df[(df.Strain == 'B') & (df.Treatment == 2)]

    fig, ax = plt.subplots()
    treatment_names = ['1-Day', '10-Day', '100-Day']
    treatment_colors = ['#87CEEB', '#FFA500', '#FF6347']
    ax.plot(B_1.Day.values, B_1.iRep_unfiltered.values, marker='o', alpha = 0.8, \
        linestyle='', ms=12, c = '#87CEEB')
    ax.plot(B_10.Day.values, B_10.iRep_unfiltered.values, marker='o', alpha = 0.8, \
        linestyle='', ms=12, c = '#FFA500')
    ax.plot(B_100.Day.values, B_100.iRep_unfiltered.values, marker='o', alpha = 0.8, \
        linestyle='', ms=12, c = '#FF6347')

    fig_name = mydir + 'figs/B_iRep_scatter.png'
    plt.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)

    plt.close()

#get_iRep(day = 'D300')
#clean_iRep()
#plot_iRep()
plot_iRep_time()
