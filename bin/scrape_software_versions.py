#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re
from os import path

# TODO Add additional regexes for new tools in process get_software_versions
regexes = {
    'Pipeline': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
	'Bowtie2': ['v_bowtie2.txt', r"version (\S+)"],
	'Cutadapt': ['v_cutadapt.txt', r"(\S+)"],
	'STAR': ['v_star.txt', r"(\S+)"],
	'samtools': ['v_samtools.txt', r"samtools (\S+)"],
	'bedtools': ['v_bedtools.txt', r"bedtools v(\S+)"],
	'deeptools': ['v_deeptools.txt', r"deeptools (\S+)"],
	'Python': ['v_python.txt', r"Python (\S+)"],
    'R': ['v_R.txt', r"R version (\S+)"]
}

results = OrderedDict()
results['Pipeline'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['Bowtie2'] = '<span style="color:#999999;\">N/A</span>'
results['STAR'] = '<span style="color:#999999;\">N/A</span>'
results['samtools'] = '<span style="color:#999999;\">N/A</span>'
results['bedtools'] = '<span style="color:#999999;\">N/A</span>'
results['deeptools'] = '<span style="color:#999999;\">N/A</span>'
results['Python'] = '<span style="color:#999999;\">N/A</span>'
results['R'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    if path.isfile(v[0]):
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'scChip-seq pipeline software versions'
section_name: 'Software Versions'
section_href: 'https://gitlab.curie.fr/data-analysis/chip-seq/t'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")