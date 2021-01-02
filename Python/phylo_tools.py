from __future__ import division
import os
from collections import Counter
import numpy as np
import colorsys

import matplotlib.colors as clr
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

import scipy.stats as stats
from statsmodels.base.model import GenericLikelihoodModel

from scipy.linalg import block_diag
from sklearn.metrics.pairwise import euclidean_distances

import parse_file
#import timecourse_util


np.random.seed(123456789)



def get_path():
    return os.path.expanduser("~/GitHub/Bacillus_Evol_Timeseries")

taxa = ['B','S']
treatments = ['0', '1', '2']
replicates = ['1','2','3','4','5']

sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i']

colors_dict = {'0':'#87CEEB', '1': '#FFA500', '2':'#FF6347'}

samples_to_remove = {'1B4':[900],
                        '1B5':[1000],
                        '2S2':[700],
                        '2S3':[700],
                        '2S4':[700],
                        '2S5':[700]}

treatment_label_dict = {'0': '1-day', '1':'10-days', '2':'100-days'}


linestyle_dict = {'B':'--', 'S':':'}

B_S_generation_dict = {'B': { '0':3321, '1':1171, '2': 107},
                        'S': {'0':3321, '1':544, '2':163} }

def get_B_S_generations(strain, treatment, day_cutoff=500):

    return B_S_generation_dict[strain][treatment] * (day_cutoff/1000)




# Michaelis-Menten
def hyperbolic_michaelis_menten(t_star, b0, K, v_max):
    t_star = np.asarray(t_star)

    return b0 + ( (v_max*t_star)/(t_star+K) )


# function to generate confidence intervals based on Fisher Information criteria
def CI_FIC(results):
    # standard errors = square root of the diagnol of a variance-covariance matrix
    ses = np.sqrt(np.absolute(np.diagonal(results.cov_params())))
    cfs = results.params
    lw = cfs - (1.96*ses)
    up = cfs +(1.96*ses)
    return (ses, lw, up)


class fit_hyperbolic_michaelis_menten(GenericLikelihoodModel):
    def __init__(self, endog, exog, **kwds):
        super(fit_hyperbolic_michaelis_menten, self).__init__(endog, exog, **kwds)
        #print len(exog)

    def nloglikeobs(self, params):
        b0 = params[0]
        K = params[1]
        v_max = params[2]
        z = params[3]
        # probability density function (pdf) is the same as dnorm
        exog_pred = hyperbolic_michaelis_menten(self.endog, b0 = b0, K = K, v_max = v_max)
        # need to flatten the exogenous variable
        LL = -stats.norm.logpdf(self.exog.flatten(), loc=exog_pred, scale=np.exp(z))
        return LL

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, method="bfgs", **kwds):

        if start_params is None:
            b0_start = 1
            K_start = 2
            v_max_start = 0.5
            z_start = 0.8

            start_params = np.array([b0_start, K_start, v_max_start, z_start])

        return super(fit_hyperbolic_michaelis_menten, self).fit(start_params=start_params,
                                maxiter=maxiter, method = method, maxfun=maxfun,
                                **kwds)




def fit_hyperbolic_michaelis_menten_best_parameters(t_array, M_array, interceptGuess=0.5):
    # assumes that you've rescaled the x axis so values start at zero
    K_start_list = [0.05,0.1,1]
    v_max_start_list = [-5,-0.5,0.1,5,10,20]
    z_start_list = [-2,-1,-0.5]
    # and while keeping the following initial values constant
    b0_start = interceptGuess

    model = fit_hyperbolic_michaelis_menten(t_array, M_array)


    def loop_params(_model, _b0_start, K_start_list, v_max_start_list, z_start_list):
        _results = []
        for _K_start in K_start_list:
            for _v_max_start in v_max_start_list:
                for _z_start in z_start_list:
                    # b0, A, umax, L, z
                    start_params = [_b0_start, _K_start, _v_max_start, _z_start]
                    _result = model.fit(start_params = start_params, method="lbfgs", \
                        bounds= [(-5,5), (0.01,10000), (0.01,1000), (-20, 20)], \
                        disp = False)

                    #_result = model.fit(start_params = start_params,  method="bfgs", disp=False)

                    _results.append(_result)
        AICs = [_result.aic for _result in _results]
        _best = _results[AICs.index(min(AICs))]
        return _best


    best = loop_params(model, b0_start, K_start_list, v_max_start_list, z_start_list)


    best_CI_FIC = CI_FIC(best)
    best_CI = best.conf_int()
    best_params = best.params


    #print(best_CI)
    #print(best_CI[1,:])

    # order of paramters
    # b0, K, v_max, z
    ses_V_max = best_CI_FIC[0][2]
    ses_K = best_CI_FIC[0][1]

    #best_CI_FIC_lower_V_max = best_CI_FIC[1][2]
    #best_CI_FIC_upper_V_max = best_CI_FIC[2][2]


    #best_CI_lower_V_max = best_CI[2,0]
    #best_CI_upper_V_max = best_CI[2,1]

    # b0, Km V_max, z, CIs
    #return best_params[0], best_params[1], best_params[2], best_params[3], best_CI_lower_V_max, best_CI_upper_V_max

    return best_params[0], best_params[1], best_params[2], best_params[3], ses_K, ses_V_max







def round_sf(number, significant):
    return round(number, significant - len(str(number)))



def hyperbolic_trajectory(ts,log_m):
    result = numpy.log(1+exp(logb)*ts)*exp(loga)
    return result


def fit_hyperbolic_trajectory(ts, log_m, t_min=100):
    ts_adjusted = ts - t_min

    loga0 = 0
    logb0 = log((exp(xs[-1])-1)/ts[-1])
    xmin = fmin(lambda x: numpy.square(xs-powerlaw_trajectory(ts,x[0],x[1])).sum(),numpy.array([loga0,logb0]))
    a = exp(xmin[0])
    b = exp(xmin[1])
    return a,b







latex_formal_dict = {  'B': r'$\mathit{Bacillus\, subtilis} \; \mathrm{NCIB \, 3610 \, \mathrm{WT}}$',
                'S': r'$\mathit{Bacillus\, subtilis} \; \mathrm{NCIB \, 3610} \, \Delta spo0A $'}



latex_dict = {  'B': r'$\mathit{B.\, subtilis} \; \mathrm{WT}$',
                'S': r'$\mathit{B.\, subtilis} \; \Delta spo0A $'}


latex_bold_dict = {  'B': r"$\mathbf{\textit{B.\, subtilis} \; WT}$",
                'S': r'$\mathit{B.\, subtilis} \; \Delta spo0A$'}




#latex_genus_dict = {  'B': r'$\mathit{Bacillus} \, \mathrm{wt} $',
#                'S': r'$\mathit{Bacillus} \, \Delta \mathrm{spo0A} $',
#                'C': r'$\mathit{Caulobacter}$',
#                'D': r'$\mathit{Deinococcus}$',
#                'P': r'$\mathit{Pseudomonas}$',
#                'F': r'$\mathit{Pedobacter}$',
#                'J': r'$\mathit{Janthinobacterium}$'
#                }


#latex_genus_dict = {  'B': r'$\mathit{Bacillus} $'  }


#latex_genus_bold_dict = {'B': r'$\mathbf{\mathit{Bacillus} }$',
#                'S': r'$\mathbf{\mathit{Bacillus} \, \Delta \mathrm{spo0A}} $',
#                'C': r'$\mathbf{\mathit{Caulobacter}}$',
#                'D': r'$\mathbf{\mathit{Deinococcus}}$',
#                'P': r'$\mathbf{\mathit{Pseudomonas}}$',
#                'F': r'$\mathbf{\mathit{Pedobacter}}$',
#                'J': r'$\mathbf{\mathit{Janthinobacterium} }$'}


latex_genus_bold_dict = {'B': r'$\mathbf{\mathit{Bacillus} }$'}




genus_dict = {'B':'Bacillus',
            'C':'Caulobacter',
            'D':'Deinococcus',
            'F':'Pedobacter',
            'J':'Janthinobacterium',
            'P':'Pseudomonas'}


def get_p_value_latex(p_value, alpha=0.05):

    if p_value < alpha:
        return r'$\mathrm{p} < 0.05$'

    else:
        return r'$\mathrm{p} \nless 0.05$'













def run_permanova(PC_space, N_list, iter=10000):

    N_list = np.asarray(N_list)

    F_obs = get_F_2(PC_space, N_list)

    F_permute_list = []

    for i in range(iter):

        PC_space_permute = PC_space[np.random.permutation(PC_space.shape[0]),:]
        F_permute_list.append(get_F_2(PC_space_permute, N_list))

    p = len([j for j in F_permute_list if j > F_obs]) / iter

    return F_obs, p


def get_F_2(PC_space, N_list):
    '''
    Modified F-statistic from Anderson et al., 2017 doi: 10.1111/anzs.12176
    Function assumes that the rows of the count matrix are sorted by group
    i.e., group one is first N1 rows, group two is N2, etc
    '''
    #N = N1+N2
    N = sum(N_list)
    dist_matrix = euclidean_distances(PC_space, PC_space)
    A = -(1/2) * (dist_matrix ** 2)
    I = np.identity(N)
    J_N = np.full((N, N), 1)
    G = (I - ((1/N) * J_N )) @ A @ (I - ((1/N) * J_N ))
    # n matrix list
    n_list = []
    for N_i in N_list:
        n_list.append((1/N_i) * np.full((N_i, N_i), 1))
    #n1 = (1/N1) * np.full((N1, N1), 1)
    #n2 = (1/N2) * np.full((N2, N2), 1)
    #H = block_diag(n1, n2) - ((1/N) * J_N )
    H = block_diag(*n_list) - ((1/N) * J_N )
    # indicator matrices
    # get V matrices
    V_list = []
    for i in range(len(N_list)):
        if i == 0:
            U_i = np.diag( N_list[i]*[1] + sum(N_list[i+1:])*[0])
        elif i == len(N_list) - 1:
            U_i = np.diag( sum(N_list[:i])*[0] + N_list[i]*[1] )
        else:
            U_i = np.diag( sum(N_list[:i])*[0] + N_list[i]*[1] +  sum(N_list[i+1:])*[0])

        V_i = np.trace(((I - H) @ U_i @ (I - H)) @ G ) / (N_list[i]-1)
        V_list.append(V_i)

    F_2 = np.trace(H @ G) / sum( [ (1 - (N_list[i]/N) ) *  V_list[i] for i in range(len(N_list)) ]  )



    return F_2







def get_taxon_ls(taxon):

    if taxon == 'S':
        ls = ':'
    else:
        ls ='--'

    return ls



def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])



#def samples_to_remove(population):


#    if population not in population_dict:
#        return None
#    else:
#        return population_dict[population]



#'0D1':[600],
#'0D2':[900],
#'0D3':[900],
#'0D4':[900],
#'2S2':[700],
#'2S3':[700],
#'2S4':[700],
#'2S5':[700],
#'1P5':[700],
#'2P1':[500],
#'1F5':[500],
#'1F2':[700],
#'2J4':[500]



def get_treatment_name(treatment):
    treatment_dict = {'0': '1-day',
                        '1': '10-day',
                        '2': '100-day'}
    return treatment_dict[str(treatment)]



def plot_species_marker(taxon):

    plot_species_marker_dict = {"B": "o",
                                "S": "o"}

    return plot_species_marker_dict[taxon]


def plot_species_fillstyle(taxon):

    plot_species_fillstyle_dict = {"B": "full",
                                "S": "none"}

    return plot_species_fillstyle_dict[taxon]


def get_colors(treatment):
    get_colors_dict = {'0':'#87CEEB', '1': '#FFA500', '2':'#FF6347'}
    return get_colors_dict[treatment]


#def get_scatter_edge_color(strain, treatment):



def get_scatter_facecolor(taxon, treatment):

    if taxon == 'S':
        return 'white'
    else:
        return get_colors(treatment)




def get_genome_size(taxon):
    genome_size_dict = {"B": 4299822,
                        "S": 4299822}

    return genome_size_dict[taxon]





def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.
    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.
    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.
    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.
    Returns
    -------
    matplotlib.patches.Ellipse
    Other parameters
    ----------------
    kwargs : `~matplotlib.patches.Patch` properties
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor=facecolor,
        **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)




def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """

    try:
        c = clr.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*clr.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])




def mut_freq_colormap():
    #cmap = clr.LinearSegmentedColormap.from_list('Zissou1', ["#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"], N=256)
    cmap = clr.LinearSegmentedColormap.from_list('Darjeeling1', ["#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"], N=256)

    # sample from cmap using uniform dist b/w 0 and 1
    u = np.random.uniform()
    rgb = '#%02x%02x%02x' % tuple([int(x * 100) for x in list(cmap(u))[:-1]])
    #tuple([int(x * 100) for x in list(cmap(u))[:-1]])
    # RGB six digit code
    return rgb





def get_populations(taxon):

    pop_dict = {"B": ['0B1', '0B2', '0B3', '0B4', '0B5',
                        '1B1', '1B2', '1B4', '1B4', '1B5',
                        '2B1', '2B2', '2B4', '2B4', '2B5'],
                "S": ['0S1', '0S2', '0S3', '0S4', '0S5',
                        '1S1', '1S2', '1S4', '1S4', '1S5',
                        '2S1', '2S2', '2S4', '2S4', '2S5']
                }



def get_ref_gbff_dict(taxon):

    ref_dict = {"B": "data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCF_002055965.1_ASM205596v1_genomic.gbff",
                "S": "data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCF_002055965.1_ASM205596v1_genomic.gbff"}

    return ref_dict[taxon]



def get_ref_fna_dict():

    ref_dict = {"B": "data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCF_002055965.1_ASM205596v1_genomic.fna",
                "S": "data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCF_002055965.1_ASM205596v1_genomic.fna"}
    return ref_dict





def common_entries(*dcts):
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)

def get_bacillus_mut_bias():
    return  1.2739

def get_bacillus_mut_rate():
    return 3.35 * (10**-10)

def get_bacillus_indel_rate():
    return 1.20 * (10**-10)



def get_ts_tv_dict():

    ts_tv_dict = {
    ("A", "C"):"V", ("A", "G"):"S", ("A", "T"):"V",
    ("C", "A"):"V", ("C", "G"):"V", ("C", "T"):"S",
    ("G", "A"):"S", ("G", "C"):"V", ("G", "T"):"V",
    ("T", "A"):"V", ("T", "C"):"S", ("T", "G"):"V"}

    return ts_tv_dict

def get_codon_dict():
    # translation table 11
    codon_dict = {
        "TTT":"F", "TCT":"S", "TAT":"Y", "TGT":"C",
        "TTC":"F", "TCC":"S", "TAC":"Y", "TGC":"C",
        "TTA":"L", "TCA":"S", "TAA":"*", "TGA":"*",
        "TTG":"L", "TCG":"S", "TAG":"*", "TGG":"W",

        "CTT":"L", "CCT":"P", "CAT":"H", "CGT":"R",
        "CTC":"L", "CCC":"P", "CAC":"H", "CGC":"R",
        "CTA":"L", "CCA":"P", "CAA":"Q", "CGA":"R",
        "CTG":"L", "CCG":"P", "CAG":"Q", "CGG":"R",

        "ATT":"I", "ACT":"T", "AAT":"N", "AGT":"S",
        "ATC":"I", "ACC":"T", "AAC":"N", "AGC":"S",
        "ATA":"I", "ACA":"T", "AAA":"K", "AGA":"R",
        "ATG":"M", "ACG":"T", "AAG":"K", "AGG":"R",

        "GTT":"V", "GCT":"A", "GAT":"D", "GGT":"G",
        "GTC":"V", "GCC":"A", "GAC":"D", "GGC":"G",
        "GTA":"V", "GCA":"A", "GAA":"E", "GGA":"G",
        "GTG":"V", "GCG":"A", "GAG":"E", "GGG":"G"
        }

    return codon_dict




class classFASTA:

    def __init__(self, fileFASTA):
        self.fileFASTA = fileFASTA

    def readFASTA(self):
        '''Checks for fasta by file extension'''
        file_lower = self.fileFASTA.lower()
        '''Check for three most common fasta file extensions'''
        if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.fna') or \
        file_lower.endswith('.faa') or file_lower.endswith('.ffn'):
            with open(self.fileFASTA, "r") as f:
                return self.ParseFASTA(f)
        else:
            print("Not in FASTA format.")

    def ParseFASTA(self, fileFASTA):
        '''Gets the sequence name and sequence from a FASTA formatted file'''
        fasta_list=[]
        for line in fileFASTA:
            if line[0] == '>':
                try:
                    fasta_list.append(current_dna)
            	#pass if an error comes up
                except UnboundLocalError:
                    #print "Inproper file format."
                    pass
                current_dna = [line.lstrip('>').rstrip('\n'),'']
            else:
                current_dna[1] += "".join(line.split())
        fasta_list.append(current_dna)
        '''Returns fasa as nested list, containing line identifier \
            and sequence'''
        return fasta_list
