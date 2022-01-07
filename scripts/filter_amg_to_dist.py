#!python
"""Contains main entry point to the program and helper functions"""
import os
import click
import pandas as pd
from pathlib import Path


@click.command()
@click.argument('amg_summary_in', type=click.Path(exists=True), required=1)
                #, help='VirHost output for filtering and transformation')
@click.argument('dist_modual', type=click.Path(exists=True), required=1)
# , multiple=True
@click.argument('spor_amg_out', type=click.Path(exists=False), required=1)
                #, help="Path for the output folder. This can't alredy exist.")
@click.option('--only_names', is_flag=True, help="Only give gene name, output will be ignored")
@click.option('--dram', is_flag=True, help="if the file is not DRAM-v but dram. Use the summary xlsx and this flag.")
def spor_from_dist(amg_summary_in:str,
                   dist_modual:str,
                      spor_amg_out:str,
                      only_names:bool=False,
                      dram:bool=False
                      ):
    """
    :param amg_summary_in:
    :param spor_amg_out:
    :param only_names:
    """
    if dram:
        data = pd.read_excel(amg_summary_in, sheet_name='MISC')
    else:
        data = pd.read_csv(amg_summary_in, sep='\t')
    headers_df = pd.concat([pd.read_csv(i, sep='\t') for i in [dist_modual]])
    headers = set(headers_df['header'].dropna())
    in_mod = data.apply(lambda x: x['header'] in headers, axis=1)
    udata = data[in_mod]
    if only_names:
        for i in udata['gene'].unique():
            print(i)
        return
    udata.to_csv(spor_amg_out, index=False, sep='\t')


if __name__ == '__main__':
    spor_from_dist()

