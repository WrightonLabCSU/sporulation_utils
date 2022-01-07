#!python
"""Contains main entry point to the program and helper functions"""
import os
import click
import pandas as pd
import numpy as np
from pathlib import Path
import logging

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


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


def add_sporulation_label(data, sporulator_def):
    is_spor = pd.read_csv(sporulator_def)
    # spilt family out
    data = taxa_split(data)
    is_spor['gtdb_f'] = is_spor['gtdb_f'].str.replace('f__', '', regex=False)
    # unique to family, TODO can family be duplicated with different values
    # for other parts of the taxonomy
    is_spor = is_spor[['f_spor', 'gtdb_f']].drop_duplicates()
    is_spor.columns = ['sporulator', 'family']
    is_spor.dropna(inplace=True)
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
@click.option('-l', '--score_limit', default=0.2,
              help="Maximum score to except.")
def parse_spor_by_fam(virhost:str,
                      taxonomy:str,
                      sporulator_def:str=None,
                      score_limit:float=0.2,
                      output_path:str='virhost_sporulators',
                      ):
    """
    :param virhost_out: Output CSV from virhostmatcher
    :param taxonomy: A taxonomy that matches to the
    :param sporulator_def:
    """
    # Read transformed virhost data
    vh_data = pd.read_csv(virhost, index_col=0).T
    vh_data.index = vh_data.index.str.replace('.fa', '', regex=False)
    tx_data = pd.read_csv(taxonomy, index_col=0)
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
    logging.info("Number of unique mags in both VirHost and taxonomy: %i",
                 len(data.index.unique()))
    logging.info("Number of VirHost mags dropped because there is"
                 " no match in taxonomy: %i",
                 len(set(vh_data.index) - set(data.index)))


    # sporulator or non sporulator flage
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
    logging.info("The number of observations removed by the value filter: %i",
                 (pre_filter_len - len(output)))

    outpath = Path(output_path)
    outpath.mkdir(parents=True, exist_ok=False)
    output.to_csv(outpath / 'mag_vmag_info.csv', index=False)


if __name__ == '__main__':
    parse_spor_by_fam()

# varry small number of matches but that is fine this is a test case
# TODO add option to skip sporulation
# TODO alow the option to take minimum
# TODO filter to value <= 0.2
