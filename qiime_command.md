## This page describe qiime command that used in the manuscript


# de-miltiplex samples
This is done by python script that written by Jin

# merge paired end by Pandaseq
```
pandaseq -f 2014188001.R1.fastq -r 2014188001.R2.fastq -w 2014188001.merged.fasta -g 2014188001.log -B -A simple_bayesian -l 253 -L 253 -t 0.9
```
#combine files
```
add_qiime_labels.py -i merged/ -m jin_mapping_with_description.unix.file.txt -c '#SampleID.merged.fasta' -n 1
```
# open OTU pick
```
pick_open_reference_otus.py -i combined_seqs.fna -o open_ref
```
