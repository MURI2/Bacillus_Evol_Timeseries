
import phylo_tools as pt
import pandas as pd

def friedman(data=None, dv=None, within=None, subject=None):
    """Friedman test for repeated measurements.
    Parameters
    ----------
    data : pandas DataFrame
        DataFrame
    dv : string
        Name of column containing the dependant variable.
    within : string
        Name of column containing the within-subject factor.
    subject : string
        Name of column containing the subject identifier.
    Returns
    -------
    stats : DataFrame
        Test summary ::
        'Q' : The Friedman Q statistic, corrected for ties
        'p-unc' : Uncorrected p-value
        'dof' : degrees of freedom
    Notes
    -----
    The Friedman test is used for one-way repeated measures ANOVA by ranks.
    Data are expected to be in long-format.
    Note that if the dataset contains one or more other within subject
    factors, an automatic collapsing to the mean is applied on the dependant
    variable (same behavior as the ezANOVA R package). As such, results can
    differ from those of JASP. If you can, always double-check the results.
    Due to the assumption that the test statistic has a chi squared
    distribution, the p-value is only reliable for n > 10 and more than 6
    repeated measurements.
    NaN values are automatically removed.
    Examples
    --------
    Compute the Friedman test for repeated measurements.
    >>> from pingouin import friedman, read_dataset
    >>> df = read_dataset('rm_anova')
    >>> friedman(data=df, dv='DesireToKill', within='Disgustingness',
    ...          subject='Subject')
                      Source  ddof1      Q     p-unc
    Friedman  Disgustingness      1  9.228  0.002384
    """

    # Collapse to the mean
    data = data.groupby([subject, within]).mean().reset_index()

    # Remove NaN
    if data[dv].isnull().any():
        data = remove_rm_na(dv=dv, within=within, subject=subject,
                            data=data[[subject, within, dv]])

    # Extract number of groups and total sample size
    grp = data.groupby(within)[dv]
    rm = list(data[within].unique())
    k = len(rm)
    X = np.array([grp.get_group(r).values for r in rm]).T
    n = X.shape[0]

    # Rank per subject
    ranked = np.zeros(X.shape)
    for i in range(n):
        ranked[i] = scipy.stats.rankdata(X[i, :])

    ssbn = (ranked.sum(axis=0)**2).sum()

    # Compute the test statistic
    Q = (12 / (n * k * (k + 1))) * ssbn - 3 * n * (k + 1)

    # Correct for ties
    ties = 0
    for i in range(n):
        replist, repnum = scipy.stats.find_repeats(X[i])
        for t in repnum:
            ties += t * (t * t - 1)

    c = 1 - ties / float(k * (k * k - 1) * n)
    Q /= c

    # Approximate the p-value
    ddof1 = k - 1
    p_unc = scipy.stats.chi2.sf(Q, ddof1)

    # Create output dataframe
    stats = pd.DataFrame({'Source': within,
                          'ddof1': ddof1,
                          'Q': np.round(Q, 3),
                          'p-unc': p_unc,
                          }, index=['Friedman'])

    col_order = ['Source', 'ddof1', 'Q', 'p-unc']

    stats = stats.reindex(columns=col_order)
    stats.dropna(how='all', axis=1, inplace=True)

    return stats



data = pd.read_csv(pt.get_path() +'/data/rm_anova.csv', sep=',' )

print(data)
