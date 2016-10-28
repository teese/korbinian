import numpy as np
import pandas as pd
from multiprocessing import Pool
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import plotly.plotly as py
import plotly.graph_objs as go


def calc_AAIMON_aa_prop_norm_factor(aa_ident):
    """Calculates the normalization factor for each selected homologue of that protein.

    Parameters
    ----------
    aa_ident : float or pd.Series
        observed conservation = FASTA gapped percentage identity of each homologue

    Returns
    -------
    aa_prop_norm_factor : float
        calculated normalization factor of AAIMON ratio for each homologue based on its gapped identity
    """

    # percentage of observed 'invisible' aa substitution within TM / non-TM region
    rand_TM = 0.13
    rand_nonTM = 0.05

    observed_changes = 1 - aa_ident
    real_changes_TM = observed_changes / (1 - rand_TM)
    real_changes_nonTM = observed_changes / (1 - rand_nonTM)
    # calculation of the invisible changes (difference between observed and real changes)
    X_TM = real_changes_TM * rand_TM
    X_nonTM = real_changes_nonTM * rand_nonTM
    # calculation of real conservation (observed conservation - invisible changes)
    real_conserv_TM = aa_ident - X_TM
    real_conserv_nonTM = aa_ident - X_nonTM
    # calculation of normalization factor
    aa_prop_norm_factor = real_conserv_nonTM / real_conserv_TM

    return aa_prop_norm_factor


def create_graph_for_norm_factor():
    arr = np.linspace(0.14, 1.0, 87)
#    print(arr)
    ser = pd.Series(arr)
    norms = calc_AAIMON_aa_prop_norm_factor(ser)

#    fig, ax = plt.subplots()
#    ax.plot(norms)
#    ax.set_xticklabels(arr, fontsize=10)

    trace1 = go.Scatter(
        x=arr,
        y=norms,
        mode='lines+markers',
        name='norm_factor',
        marker=dict(
            size=0.5,
        ),
    )
    data = [trace1]

    layout = go.Layout(
        autosize=True,
        yaxis=dict(
            title='normalization factor',
            titlefont=dict(
                size=20,
                color='black'
            ),
            range=[0,5],
        ),
        xaxis=dict(
            title='% identity',
            titlefont=dict(
                size=20,
                color='black'
            ),
        ),
        paper_bgcolor='rgb(255,255,255)',
        plot_bgcolor='rgb(255,255,255)',
        showlegend=False
    )

    ax = go.Figure(data=data, layout=layout)
    # py.iplot(ax)

    py.image.save_as(ax, filename='norm_factor_ver3.png')

#create_graph_for_norm_factor()


def create_graph_for_normalized_AAIMON(acc, norm_AAIMON, per_ident):
    fig, ax = plt.subplots()
    ax.scatter(per_ident, norm_AAIMON)
#    ax.set_xticklabels(per_ident, rotation=45)
    ax.set_ylabel('normalized AAIMON', fontsize=14)
    pylab.ylim([0.6, 1.5])
    fig.savefig(r'C:\Users\coffee oder tee\Dropbox\Thesis Shenger\my_documents\plot_folder\%s_normalized_AAIMON.png' % acc,
                format="png", dpi=200)
