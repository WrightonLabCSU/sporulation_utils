# Sporulation Utilities

DRAM-v and DRAM are powerful tools that generate massive amounts of data and can distill that information into summaries of organismal metabolisms.
With much power comes much complexity, great time and effort has gone into the metabolic pathways that are included with DRAM/DRAM-v, but these are only the surface of DRAMs potential.

In this Repository we are collecting a set of tools that allow DRAM-v, and to an extent DRAM, to explore sporulation, throgh auxiliary metabolic genes (AMGs) that are key to it.
The first step of this process is simply to run DRAM-v on the V-MAGs of interest, and then to use DRAM-v's distillation tool to apply a custom distillation module.
However, after that point things get tricky.
 Most of the tools in this repo are to help with those tricky next steps.

## Setting up Environment

You need to somehow install click, numpy, and pandas in your python environment, the how is up to you, conda or pip should both work.
If you are in the lab and working on Zenith you can use the scripts conda environment, which exists for this exact reason.
You should contact Rory Flynn to get access, or whoever is maintaining it if you are from the future.

## Getting started

The tools in this repository are described below with examples.
In order to use the example code, you will need to download the full github repo using the command below.

```
git clone https://github.com/WrightonLabCSU/sporulation_utils.git
```

However, if you want to use one tool or file, and don't need the examples, just use the wget links in the tool sections.


## Sporulation distillate module
To download this file:

```
wget https://raw.githubusercontent.com/WrightonLabCSU/sporulation_utils/main/scripts/filter_amg_to_dist.py
```

In the root of the repository you will find the `sporulation_distillate_module.tsv` that is used to create a custom dram/dram-v distillate with sporulation information included.



And to use this file with `DRAM-v distill` structure your command like falows:

```
DRAM-v.py distill \
          --custom_distillate 'sporulation_distillate_module.tsv' \
          -i 'dram_output/annotations.tsv' \
          -o 'out_distillate_location'
```

[Click here to learn more about how the DRAM-v distillation step works.](https://github.com/WrightonLabCSU/DRAM/wiki/3b.-Running-DRAM-v#dram-v-distill)


## Filter a DRAM AMG/Metabolic Summary File With a Custom Distillate
To download this file:

```
wget https://raw.githubusercontent.com/WrightonLabCSU/sporulation_utils/main/
```

This script is not limited to this sporulation project, but it is still a key part.
When you run DRAM-v with a custom distillate sheet you still get all the AMG information that would normally come with DRAM-v.
This script will take only the AMG information related to the distillate sheet that it is given as an argument.

If you are in the root of the downloaded git repository, you can test this script with this command.

```
python scripts/filter_amg_to_dist.py \
        example_data/filter_amg_to_dist/amg_summary.tsv \ # AMG summary file
        sporulation_distillate_module.tsv \ # distillations module from this repo
        sporulation_amg_summary.tsv # name of the output file

```

If you are only interested in the genes, you can use the --names_only flag which will spit out the gene names related to your distillate module to standard out.
This output can be piped into other commands, here is an example where we count them

```
python scripts/filter_amg_to_dist.py \
        example_data/filter_amg_to_dist/amg_summary.tsv \
        sporulation_distillate_module.tsv \
        NA \
        --only_names | wc -l
```

## Filter and Transform a VirHostMatcher File
To download this file:

```
wget https://raw.githubusercontent.com/WrightonLabCSU/sporulation_utils/main/scripts/spor_by_fam.py
```

When VirHostMatcher is run  sporulation AMG containing viruses, identified with DRAM-v, the output needs to be transformed in order to be useible.
This script performs several transformations.

 * Using a required taxonomy file, with columns 'vMAG_id' and 'GTDB_string' the script will merge in the GTDB_string dropping vMAGs that don't match.
 * Using an optional sporulation_list file containing GTDB strings and sporulation status, the script will add in spoulation status based on the family level of the taxonomy.
 * Finally, the script will filter the data based on d2star measurement, <= 0.2 by default, and give 5 columns of output: 
     * Col1 - Virus ID
     * Col2 - d2* score
     * Col3 - Host ID
     * Col4 - Host GTDB Taxonomy String 
     * Col5 - T/F on Sporulation/non-sporulating

An example of how this script is run would be:
```
# In the case where we do have sporulators by family
rm -r test_output_with_spor_file
python scripts/spor_by_fam.py \
        example_data/spor_by_fam/d2star_k6_1250.csv \
        example_data/spor_by_fam/test_gtdb_strings.csv \
        test_output_with_spor_file \
        -s example_data/spor_by_fam/test_sporlist.csv

# An example where we have no sporulation list.
rm -r test_output_no_spor_file
python scripts/spor_by_fam.py \
        example_data/spor_by_fam/d2star_k6_1250.csv \
        example_data/spor_by_fam/test_gtdb_strings.csv \
        test_output_no_spor_file

```

