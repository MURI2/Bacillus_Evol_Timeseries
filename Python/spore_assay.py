from __future__ import division
import os
import pandas as pd
import  matplotlib.pyplot as plt
import bacillus_tools as bt
#mydir = os.path.expanduser("~/GitHub/Task2")

def make_fig():
    path_IN = bt.get_path() + '/data/spore_assay/Sporulation_170912_long.txt'
    IN = pd.read_csv(path_IN, sep = '\t')
    #d100
    IN_0B1_100 = IN.loc[(IN['Pop'] == '0B1') & (IN['Day'] == 100)]
    IN_2B1_100 = IN.loc[(IN['Pop'] == '2B1') & (IN['Day'] == 100)]
    IN_mean_0B1_100 = IN_0B1_100['Vegetative_percent'].groupby(IN_0B1_100['Time_hours']).mean().reset_index()
    IN_mean_2B1_100 = IN_2B1_100['Vegetative_percent'].groupby(IN_2B1_100['Time_hours']).mean().reset_index()
    IN_std_0B1_100 = IN_0B1_100['Vegetative_percent'].groupby(IN_0B1_100['Time_hours']).std().reset_index()
    IN_std_2B1_100 = IN_2B1_100['Vegetative_percent'].groupby(IN_2B1_100['Time_hours']).std().reset_index()
    # Day 500
    IN_0B1_500 = IN.loc[(IN['Pop'] == '0B1') & (IN['Day'] == 500)]
    IN_2B1_500 = IN.loc[(IN['Pop'] == '2B1') & (IN['Day'] == 500)]
    IN_mean_0B1_500 = IN_0B1_500['Vegetative_percent'].groupby(IN_0B1_500['Time_hours']).mean().reset_index()
    IN_mean_2B1_500 = IN_2B1_500['Vegetative_percent'].groupby(IN_2B1_500['Time_hours']).mean().reset_index()
    IN_std_0B1_500 = IN_0B1_500['Vegetative_percent'].groupby(IN_0B1_500['Time_hours']).std().reset_index()
    IN_std_2B1_500 = IN_2B1_500['Vegetative_percent'].groupby(IN_2B1_500['Time_hours']).std().reset_index()

    fig = plt.figure(figsize=(8,3))

    #plt.scatter(IN_mean_0B1.Time_hours.values, IN_mean_0B1.Vegetative_percent.values, c='#87CEEB', marker='o', label='_nolegend_', s = 60)
    plt.plot(IN_mean_0B1_100.Time_hours.values, 1.001- IN_mean_0B1_100.Vegetative_percent.values, \
        'b-',  c='#87CEEB')
    plt.plot(IN_mean_2B1_100.Time_hours.values, 1.001- IN_mean_2B1_100.Vegetative_percent.values, \
        'b-',  c = '#FF6347')
    plt.errorbar(IN_mean_0B1_100.Time_hours.values, 1.001- IN_mean_0B1_100.Vegetative_percent.values, \
        IN_std_0B1_100.Vegetative_percent.values,  linestyle='None', marker='o', c='#87CEEB', elinewidth=1.5, label="1-day WT, day 100",)
    plt.errorbar(IN_mean_2B1_100.Time_hours.values, 1.001- IN_mean_2B1_100.Vegetative_percent.values, \
        IN_std_2B1_100.Vegetative_percent.values, linestyle='None', marker='o', c = '#FF6347', elinewidth=1.5, label="100-day WT, day 100",)

    plt.plot(IN_mean_0B1_500.Time_hours.values, 1.001- IN_mean_0B1_500.Vegetative_percent.values, \
        'b-',  c='#87CEEB')
    plt.plot(IN_mean_2B1_500.Time_hours.values, 1.001- IN_mean_2B1_500.Vegetative_percent.values, \
        'b-',  c = '#FF6347')
    plt.errorbar(IN_mean_0B1_500.Time_hours.values, 1.001- IN_mean_0B1_500.Vegetative_percent.values, \
        IN_std_0B1_500.Vegetative_percent.values,  linestyle='None', marker='v', c='#87CEEB', elinewidth=1.5, label="1-day WT, day 500",)
    plt.errorbar(IN_mean_2B1_500.Time_hours.values, 1.001- IN_mean_2B1_500.Vegetative_percent.values, \
        IN_std_2B1_500.Vegetative_percent.values, linestyle='None', marker='v', c = '#FF6347', elinewidth=1.5, label="100-day WT, day 500",)
    #plt.title('Bacillus sporulation', fontsize = 24)
    plt.xlabel('Time (hours)', fontsize = 18)
    plt.ylabel('Percent spores, ' + r'$log_{10}$', fontsize = 18)
    plt.ylim(0.0008, 1.1)
    plt.yscale('log')
    plt.legend(numpoints=1, prop={'size':10},  loc='lower right', frameon=False)
    fig_name = bt.get_path() + '/figs/test.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def make_fig_army():
    #path_IN = mydir + 'data/Sporulation_170912_long.txt'
    path_IN = bt.get_path() + '/data/spore_assay/Sporulation_170912_long.txt'
    IN = pd.read_csv(path_IN, sep = '\t')
    # Day 500
    IN_0B1_500 = IN.loc[(IN['Pop'] == '0B1') & (IN['Day'] == 500)]
    IN_2B1_500 = IN.loc[(IN['Pop'] == '2B1') & (IN['Day'] == 500)]
    IN_mean_0B1_500 = IN_0B1_500['Vegetative_percent'].groupby(IN_0B1_500['Time_hours']).mean().reset_index()
    IN_mean_2B1_500 = IN_2B1_500['Vegetative_percent'].groupby(IN_2B1_500['Time_hours']).mean().reset_index()
    IN_std_0B1_500 = IN_0B1_500['Vegetative_percent'].groupby(IN_0B1_500['Time_hours']).std().reset_index()
    IN_std_2B1_500 = IN_2B1_500['Vegetative_percent'].groupby(IN_2B1_500['Time_hours']).std().reset_index()

    fig = plt.figure()
    plt.plot(IN_mean_0B1_500.Time_hours.values, 1.001- IN_mean_0B1_500.Vegetative_percent.values, \
        'b-',  c='#87CEEB')
    plt.plot(IN_mean_2B1_500.Time_hours.values, 1.001- IN_mean_2B1_500.Vegetative_percent.values, \
        'b-',  c = '#FF6347')

    plt.errorbar(IN_mean_0B1_500.Time_hours.values, 1.001- IN_mean_0B1_500.Vegetative_percent.values, \
        IN_std_0B1_500.Vegetative_percent.values,  linestyle='None', marker='v', c='#87CEEB', elinewidth=1.5, label="1-day transfers",)
    plt.errorbar(IN_mean_2B1_500.Time_hours.values, 1.001- IN_mean_2B1_500.Vegetative_percent.values, \
        IN_std_2B1_500.Vegetative_percent.values, linestyle='None', marker='o', c = '#FF6347', elinewidth=1.5, label="100-day transfers",)
    plt.title('Spore accumulation curves, day 500', fontsize = 20)
    plt.xlabel('Time (hrs)', fontsize = 18)
    plt.ylabel('Percent spores, '+ r'$log_{10}$', fontsize = 18)
    plt.ylim(0.0008, 1.1)
    plt.xlim(0, 400)
    plt.yscale('log', basey = 10)
    plt.legend(numpoints=1, prop={'size':12},  loc='upper left', frameon=False)
    fig_name = mydir + 'figs/army_update_spore.png'
    fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()


#make_fig_army()
make_fig()
