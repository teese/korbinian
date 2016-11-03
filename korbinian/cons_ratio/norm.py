import pandas as pd
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import os


def calc_AAIMON_aa_prop_norm_factor(aa_ident, rand_TM, rand_nonTM):
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
#    rand_TM = 0.13
#    rand_nonTM = 0.05

    observed_changes = 100 - (aa_ident*100)
    real_changes_TM = observed_changes / (1 - rand_TM)
    real_changes_nonTM = observed_changes / (1 - rand_nonTM)
    # calculation of the invisible changes (difference between observed and real changes)
    invisible_changes_TM = real_changes_TM * rand_TM
    invisible_changes_nonTM = real_changes_nonTM * rand_nonTM
    # calculation of real conservation (observed conservation - invisible changes)
    real_conserv_TM = aa_ident*100 - invisible_changes_TM
    real_conserv_nonTM = aa_ident*100 - invisible_changes_nonTM
    # calculation of normalization factor
    aa_prop_norm_factor = real_conserv_nonTM / real_conserv_TM

    return aa_prop_norm_factor


def save_graph_for_normalized_AAIMON(acc, AAIMON, norm_AAIMON, aa_ident, zipout, protein_name):
    fig, ax = plt.subplots()
    perc_ident = aa_ident*100
    ax.scatter(perc_ident, norm_AAIMON)
#    ax.set_xticklabels(per_ident, rotation=45)
    ax.set_ylabel('normalised AAIMON', fontsize=14)
    pylab.ylim([0.6, 1.5])

    ax2 = ax.twinx()
    ax2.scatter(perc_ident, AAIMON, color='red', marker='^', alpha=0.5)
    pylab.ylim([0.6, 1.5])
    ax.set_xlabel('% identity of homologue', fontsize=14)
    xlim_min = 40
    xlim_max = 100
    ax.set_xlim(xlim_min, xlim_max)

#    fig.savefig(r'C:\Users\coffee oder tee\Dropbox\Undergraduate\Bachelor thesis\my_documents\plot_folder\%s_normalised_AAIMON.png' % acc,
#                format="png", dpi=200)
    fig.savefig(protein_name + '_AAIMON_normalisation.png', format = 'png', dpi = 200)
    plt.close('all')
    zipout.write(protein_name + '_AAIMON_normalisation.png',
                 arcname=os.path.basename(protein_name) + '_AAIMON_normalisation.png')
    os.remove(protein_name + '_AAIMON_normalisation.png')


