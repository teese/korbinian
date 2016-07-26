import os

def filter_and_save_fastA(df, dfs, dfs_filt, acc, TMD, tar_out):
    # start with the same dataframe copy that has filtered for gapped identity, etc (the AIMAN ratio is not necessary for the FastA saving)
    # remove hits lacking sequence and also remove hits with too many gaps in TMD from either query or match
    dfs_filt_FastA = dfs_filt.loc[dfs['%s_SW_match_seq' % TMD].notnull()].loc[
        dfs['%s_SW_query_acceptable_num_gaps' % TMD]].loc[
        dfs['%s_SW_match_acceptable_num_gaps' % TMD]]
    # setup the file names again. Note that the file should already exist, and the query sequence included.
    fasta_file = df.loc[acc, 'fasta_file_BASENAME'] + '%s.fas' % TMD
    fasta_file_path = df.loc[acc, 'fasta_file_BASENAMEPATH'] + '%s.fas' % TMD
    # setup filenames for fastA plus surrounding sequence (interfacial region)
    fasta_file_plus_surr = df.loc[acc, 'fasta_file_plus_surr_path_BASENAME'] + '%s.fas' % TMD
    fasta_file_plus_surr_path = df.loc[acc, 'fasta_file_plus_surr_path_BASENAMEPATH'] + '%s.fas' % TMD

    with open(fasta_file_path, 'w') as f:
        # add original query seq to fasta file. Note that the first SIMAP hit is excluded below.
        f.write('>0000_%s_%s_uniprot_query\n%s\n' % (
            df.loc[acc, 'A2_protein_name'], TMD, df.loc[acc, '%s_seq' % TMD]))
        # for each hit, add the sequence to the fastA file
        for hit in dfs_filt_FastA.loc[dfs['%s_SW_match_seq' % TMD].notnull()].index[1:]:
            # add the original query seq
            f.write('>%04d_%s_%s\n%s\n' % (hit, str(dfs_filt_FastA.loc[hit, 'A2_organism'])[:30],
                                           str(dfs_filt_FastA.loc[hit, 'A4_description'])[:30],
                                           dfs_filt_FastA.loc[hit, '%s_SW_match_seq' % TMD]))
            # logging.info('saved ' + fasta_file_path)
    tar_out.add(fasta_file_path, arcname=fasta_file)
    os.remove(fasta_file_path)

    dfs_filt_FastA_plus_surr = dfs_filt.loc[dfs['%s_SW_match_seq_plus_surr' % TMD].notnull()]
    with open(fasta_file_plus_surr_path, 'w') as f:
        # add the original query seq
        f.write('>00_%s_query_seq\n%s\n' % (df.loc[acc, 'A2_protein_name'],
                                            df.loc[acc, '%s_with_surrounding_seq' % TMD]))
        for hit in dfs_filt_FastA.loc[dfs_filt_FastA['%s_SW_match_seq_plus_surr' % TMD].notnull()].index:
            f.write('>%04d_%s_%s\n%s\n' % (hit, str(dfs_filt_FastA.loc[hit, 'A2_organism'])[:30],
                                           str(dfs_filt_FastA.loc[hit, 'A4_description'])[:30],
                                           dfs_filt_FastA.loc[hit, '%s_SW_match_seq_plus_surr' % TMD]))
            # logging.info('saved ' + fasta_file_plus_surr_path)
    tar_out.add(fasta_file_plus_surr_path, arcname=fasta_file_plus_surr)
    os.remove(fasta_file_plus_surr_path)
    return dfs_filt_FastA