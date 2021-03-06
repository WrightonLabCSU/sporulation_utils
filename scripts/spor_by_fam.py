#!python
"""Contains main entry point to the program and helper functions"""
import os
import click
import pandas as pd
import numpy as np
from pathlib import Path
import logging



def taxa_split(data):
    taxa_level = {
        'd': 'domain',
        'p': 'phyla',
        'c': 'class',
        'o': 'order',
        'f': 'family',
        'g': 'genus',
        's': 'species',
    }

    taxa_data = data.apply(
        (lambda x: {taxa_level[i.split('__', 2)[0]]: i.split('__', 2)[1]
                    for i in x['GTDB_string'].split(';')})
        , axis=1, result_type='expand')
    taxa_data.replace(r'^\s*$', np.nan, regex=True, inplace=True)
    data = pd.concat([
        data,
        taxa_data],
                     axis=1)
    return data

def add_auxiliary_score(data, amg_summary):
    amg_data = pd.read_csv(amg_summary, sep='\t')
    amg_data = amg_data[['scaffold', 'auxiliary_score']]
    amg_data = amg_data.groupby('scaffold').apply(
        lambda x: ", ".join([str(i) for i in x['auxiliary_score'].values]))
    amg_data = pd.DataFrame(amg_data, columns=['auxiliary_score'])
    amg_data.reset_index(inplace=True)
    amg_data.rename(columns={'scaffold': 'Virus ID'}, inplace=True)
    untrimed_names = set(amg_data['Virus ID'])
    trim_len = data['Virus ID'].str.len().max()
    amg_data['Virus ID'] = amg_data['Virus ID'].str.slice(0, trim_len)
    unmached_names_no_trim = set(data['Virus ID']) - untrimed_names
    unmached_names_trim = set(data['Virus ID']) - set(amg_data['Virus ID'])
    logging.info("The number of observations in the merged filtered data set"
                 " that do not match those in the amg_data set before trimming"
                 " strings to the %i character limit: %i",
                 trim_len, len(unmached_names_no_trim))
    logging.info("The number of observations in the merged filtered data set"
                 " that are not in the amg_data set after trimming, and so"
                 " will not get auxiliary score data: %i",
                 len(unmached_names_trim))
    dups = list(amg_data['Virus ID'][amg_data['Virus ID'].duplicated(keep='first')])
    if len(dups) > 0:
        logging.error("There are %i duplicated ids in the"
                      " amg_data after trimming, and so the data can't"
                      " be merged.\nThe duplicates are: %s",
                      len(dups), str(dups))
    if len(unmached_names_trim) > 0:
        logging.warning("These names in the merged filtered data set"
                        " match none in the amg_data set: %s",
                        str(unmached_names_trim))
    data = pd.merge(data, amg_data, on='Virus ID',  how='left')
    return data

def add_sporulation_label(data, sporulator_def):
    is_spor = pd.read_csv(sporulator_def)
    # spilt family out
    data = taxa_split(data)
    is_spor['gtdb_f'] = is_spor['gtdb_f'].str.replace('f__', '', regex=False)
    # unique to family, TODO can family be duplicated with different values
    # for other parts of the taxonomy
    is_spor = is_spor[['f_spor', 'gtdb_f']].drop_duplicates()
    is_spor.columns = ['sporulator', 'family']
    is_spor.fillna('Multiple lifestyle types', inplace=True)
    # is_spor.dropna(inplace=True)
    logging.info("The number of observations in the sporulation by family"
                 " file that have matches in the VirHost+taxonomy file: %i",
                 len(set(is_spor['family']).intersection(set(data['family']))))
    data = data.merge(is_spor, how='left', on='family')
    # data.fillna('unkown', inplace=True)
    return data

# @click.option('-o', '--output_path', type=click.Path(exists=False),
#               multiple=True, default='virhost_sporulators',
#               help="The output path")
# @click.argument('-p', type=click.Path(exists=True))

@click.command()
@click.argument('virhost', type=click.Path(exists=True), required=1)
                #, help='VirHost output for filtering and transformation')
@click.argument('taxonomy', type=click.Path(exists=True), required=1)
                #, help="Taxonomy file for the MAGs in the VirHost output")
@click.argument('output_path', type=click.Path(exists=False), required=1)
                #, help="Path for the output folder. This can't alredy exist.")
@click.option('-s', '--sporulator_def', type=click.Path(exists=True),
              default=None, help="The sporulator definitions.")
@click.option('--amg_summary', type=click.Path(exists=True),
              default=None, help="The sporulator definitions.")
@click.option('-l', '--score_limit', default=0.2,
              help="Maximum score to except.")
@click.option('--taxonomy_type', default=1,
              help="Maximum score to except.")
@click.option('--count_only/--count_and_write', default=False,
              help="Only make counts so you can save a small amount of time.")
def parse_spor_by_fam(virhost:str,
                      taxonomy:str,
                      sporulator_def:str=None,
                      amg_summary:str=None,
                      score_limit:float=0.2,
                      output_path:str='virhost_sporulators',
                      taxonomy_type:int=2,
                      count_only:bool=False
                      ):
    """
    :param virhost_out: Output CSV from virhostmatcher
    :param taxonomy: A taxonomy that matches to the
    :param sporulator_def:
    """
    # Read transformed virhost data
    outpath = Path(output_path)
    outpath.mkdir(parents=True, exist_ok=False)
    logging.basicConfig(filename=str( outpath / 'data_process.log'), filemode='w',
                        format='%(asctime)s %(message)s', level=logging.INFO)
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger().addHandler(console)
    vh_data = pd.read_csv(virhost, index_col=0).T
    vh_data.index = (vh_data.index
                     .str.replace('\.fasta$', '', regex=True)
                     .str.replace('\.fast$', '', regex=True)
                     .str.replace('\.fas$', '', regex=True)
                     .str.replace('\.fa$', '', regex=True)
                     .str.replace('\.fna$', '', regex=True)
                     .str.replace('\.fn$', '', regex=True)
                     .str.replace('\.f$', '', regex=True)
                     .str.replace('\.$', '', regex=True))
    if taxonomy_type == 2:
        tx_data = pd.read_csv(taxonomy, index_col=0)
    if taxonomy_type == 1:
        tx_data = pd.read_csv(taxonomy, sep='\t' )
        tx_data.set_index('user_genome', inplace=True)
        tx_data = tx_data[['classification']]
        tx_data.columns = ['GTDB_string']
    vh_data.index = vh_data.index.astype(str)
    tx_data.index = tx_data.index.astype(str)
    data = pd.merge(vh_data, tx_data, right_index=True, left_index=True)
    data.reset_index(inplace=True)

    # setup columns
    data_value_vars = set(data.columns) - {'GTDB_string'}
    data_id_vars = ['GTDB_string', 'index']
    final_column_order = ['variable', 'value', 'index',  'GTDB_string']
    final_column_names = ['Virus ID', 'Score', 'Host ID',  'Taxonomy_String']

    # add to summary


    logging.info("Number of unique mags in VirHost file: %i",
                 len(vh_data.index.unique()))
    logging.info("Number of unique mags in taxonomy file: %i",
                 len(tx_data.index.unique()))
    if len(data) < 1:
        logging.error("The names of the mags in virHost don't match those"
                      " in the taxonomy. This prevents them from merging!")
    logging.info("Number of unique mags in both VirHost and taxonomy: %i",
                 len(data.index.unique()))
    logging.info("Number of VirHost mags dropped because there is"
                 " no match in taxonomy: %i",
                 len(set(vh_data.index) - set(data.index)))


    # sporulator or non sporulator flag
    if sporulator_def is not None:
        data = add_sporulation_label(data, sporulator_def)
        data_id_vars += ['sporulator']
        final_column_order += ['sporulator']
        final_column_names += ['Sporulator']

    output = data.melt(id_vars=data_id_vars, value_vars=data_value_vars)
    output = output[final_column_order]
    # convert to 4 columns with melt (mag, vmag, tax_str, spor_idic, score)
    output.columns = final_column_names
    pre_filter_len = len(output)
    output = output[output['Score'] <= score_limit]
    output['Virus ID'] = (output['Virus ID']
                          .str.replace(r'\.fasta$', '', regex=True)
                          .str.replace(r'\.fast$', '', regex=True)
                          .str.replace(r'\.fas$', '', regex=True)
                          .str.replace(r'\.fa$', '', regex=True)
                          .str.replace('\.fna$', '', regex=True)
                          .str.replace('\.fn$', '', regex=True)
                          .str.replace(r'\.f$', '', regex=True)
                          .str.replace(r'\.$', '', regex=True))
    logging.info("The number of observations removed by the value filter: %i",
                 (pre_filter_len - len(output)))

    
    logging.info("The number of observations removed by the value filter: %i",
                 (pre_filter_len - len(output)))
    counts = output['Sporulator'].value_counts(dropna=False)
    spor_counts = pd.DataFrame({"Labeled Sporulator": counts[True],
                                "Labeled Non-Sporulator": counts[False], 
                                "Unlabeled Assumed Non-Sporulator": counts[np.NaN], 
                                "All Non-Sporulator": counts[False] + counts[np.NaN]
                                }, index=['Count']).T

    spor_counts['Percent'] =  (spor_counts['Count']/len(output)) * 100
    logging.info(spor_counts.T)
    spor_counts.T.to_csv(outpath / 'spor_count_percent.csv')
    if count_only:
        return
    
    if amg_summary is not None:
        output = add_auxiliary_score(output, amg_summary)
    output.to_csv(outpath / 'mag_vmag_info.csv', index=False)


if __name__ == '__main__':
    parse_spor_by_fam()

# varry small number of matches but that is fine this is a test case
# TODO add option to skip sporulation
# TODO alow the option to take minimum
# TODO filter to value <= 0.2
