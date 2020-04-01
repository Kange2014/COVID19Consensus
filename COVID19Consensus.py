#!/usr/bin/python
# Copyright (C) 2019 Ion Torrent Systems, Inc. All Rights Reserved
#
# Plugin: COVID19Consensus
# This plugin is developed for Ion Ampliseq Coronavirus panel sequencing data assembly
#
# Author: Lucius Zheng
# Last modified: 2020/03/14
#

import json
import os
from django.utils.functional import cached_property
from ion.plugin import *
import subprocess
from subprocess import check_output
import shutil

from Bio import SeqIO
from Bio.Seq import Seq

from django.conf import settings
from django.template.loader import render_to_string

def createReport(reportName,reportTemplate,reportData):
    with open(reportName,'w') as bcsum:
        bcsum.write( render_to_string(reportTemplate,reportData) )

class COVID19Consensus(IonPlugin):
    # The version number for this plugin
    version = "1.0.0.0"
    
    # this plugin can run on fullchip runs, thumbnail runs, and composite (merged via project page) runs
    # note that when the plugin is manually launched, only the 'launch' method will be called
    runtypes = [RunType.FULLCHIP, RunType.THUMB, RunType.COMPOSITE]

    # specify when the plugin is called.  For log parsing, stay simple and just get called when the run completes.
    # but can also be called before the run starts, at the block level, or after all other default plugins run
    runlevels = [RunLevel.DEFAULT]
    
    # a simple cached version of the start plugin property
    @cached_property
    def startplugin_json(self):
        return self.startplugin

    @cached_property
    def barcodes_json(self):
        with open('barcodes.json', 'r') as barcodes_handle:
            return json.load(barcodes_handle)
    
    def launch(self, data=None):
        """This is the primary launch method for the plugin."""
    
        # configure django to use the templates folder        
        #settings.configure(TEMPLATE_DIRS=(self.startplugin["runinfo"]["plugin_dir"] + '/templates'),)
        
        if not settings.configured:
            settings.configure( DEBUG=False, TEMPLATE_DEBUG=False,
                                INSTALLED_APPS=('django.contrib.humanize',),
                                TEMPLATE_DIRS=(os.path.join(self.startplugin["runinfo"]["plugin_dir"],'templates'),) 
                            )
        
        
        # define plugin results DIR variable
        results_dir = self.startplugin_json['runinfo']['results_dir']
        
        # start to analyze bam files
        
        Assembly_report = []
        for barcode_name, barcode_values in self.barcodes_json.iteritems():
            # do you work per barcode here!    
            
            # first check to see if the barcode was excluded using the frame work barcodes configuration table        
            selected = True
            barcodeData = self.startplugin_json['pluginconfig'].get('barcodetable',None)
            if barcodeData:
                #print(barcodeData)
                for bc in barcodeData:
                    if  bc.get('barcode_name',"") == barcode_name:
                        selected = bc.get('selected',True)
                        break
            
            if not selected:
                continue
            
            print("Barcode Name: " + barcode_name)
            print("Bam Filepath: " + barcode_values['bam_filepath'])
            print("Read count: " + str(barcode_values['read_count']))
            
            # if no BAM file or file size is 0, then skip the sample
            if not os.path.exists(barcode_values['bam_filepath']):
                print "BAM file does not exist. We will skip the sample in the followed analysis.\n"
                continue
            if os.path.getsize(barcode_values['bam_filepath']) == 0:
                print "BAM file size is zero. We will skip the sample in the followed analysis.\n"
                continue
                        
            if barcode_values['read_count'] < 10000:
                print "BAM file reads number is less than 10000. We will skip the sample in the followed analysis.\n"
                continue
            
            if self.startplugin_json['pluginconfig']['chip'] == "Proton P1":
                parameters = self.startplugin_json['runinfo']['plugin_dir'] + "/reference/ampliseq_germline_cov_p1_parameters.json"
                
            if self.startplugin_json['pluginconfig']['chip'] == "PGM":
                parameters = self.startplugin_json['runinfo']['plugin_dir'] + "/reference/ampliseq_germline_cov_pgm_parameters.json"
            
            if self.startplugin_json['pluginconfig']['chip'] == "510" or self.startplugin_json['pluginconfig']['chip'] == "520" or self.startplugin_json['pluginconfig']['chip'] == "530":
                parameters = self.startplugin_json['runinfo']['plugin_dir'] + "/reference/ampliseq_germline_cov_530_parameters.json"
            
            if self.startplugin_json['pluginconfig']['chip'] == "540":
                parameters = self.startplugin_json['runinfo']['plugin_dir'] + "/reference/ampliseq_germline_cov_540_parameters.json"
            
            if self.startplugin_json['pluginconfig']['chip'] == "550":
                parameters = self.startplugin_json['runinfo']['plugin_dir'] + "/reference/ampliseq_germline_cov_550_parameters.json"
            
            # 1st align-call cycle: TMAP aligns the data to SARS-CoV-2 reference genome, then TVC calls variants
            
            reference = self.startplugin_json['runinfo']['plugin_dir'] + "/reference/hg19_human_trxo_and_cor/hg19_human_trxo_and_cor.fasta"
            target_bed = self.startplugin_json['runinfo']['plugin_dir'] + "/reference/Coronavirus.20200124.bed"
            
            ## tmap
            aligned_filename = barcode_name + '.1.bam'
            align_cmd = "tmap mapall -n 12 -f " + reference + \
                        " -r " + barcode_values['bam_filepath'] + \
                        " -v -Y -u --prefix-exclude 5 -o 1 -J 25 --end-repair 15 --do-repeat-clip --context stage1 map4" + \
                        " | samtools sort -T align.sorted -@ 12 -o " + aligned_filename
            
            align_results = check_output(align_cmd, shell=True, cwd=results_dir)
            
            ## tvc
            call_cmd = "/results/plugins/variantCaller/bin/variant_caller_pipeline.py --input-bam " + aligned_filename + \
                        " --reference-fasta " + reference + \
                        " --output-dir variantCall1 " + \
                        " --parameters-file " + parameters + \
                        " --region-bed " + target_bed + \
                        " --bin-dir /results/plugins/variantCaller/bin " + \
                        " --num-threads 12"
            
            call_results = check_output(call_cmd, shell=True, cwd=results_dir)
            
            ## refine the reference using major alleles
            
            refine_cmd = "cat " + reference + " | vcf-consensus variantCall1/TSVC_variants.vcf.gz > variantCall1/TSVC_variants.consensus.fasta"
            refine_results = check_output(refine_cmd, shell=True, cwd=results_dir)
            
            # 2nd align-call cycle to minimize reference bias in the assembly
            ## first build genome index on the new reference
            build_cmd = "build_genome_index.pl -f variantCall1/TSVC_variants.consensus.fasta -s Coronavirus_with_expr_ctrl1 -l Coronavirus_with_expr_ctrl -v 1.0" 
            build_results = check_output(build_cmd, shell=True, cwd=results_dir)
            
            ## 2nd tmap
            aligned_filename = barcode_name + '.2.bam'
            align_cmd = "tmap mapall -n 12 -f " + "Coronavirus_with_expr_ctrl1/Coronavirus_with_expr_ctrl1.fasta " + \
                        "-r " + barcode_values['bam_filepath'] + \
                        " -v -Y -u --prefix-exclude 5 -o 1 -J 25 --end-repair 15 --do-repeat-clip --context stage1 map4" + \
                        " | samtools sort -T align.sorted -@ 12 -o " + aligned_filename
            
            align_results = check_output(align_cmd, shell=True, cwd=results_dir)
            
            ## 2nd tvc
            call_cmd = "/results/plugins/variantCaller/bin/variant_caller_pipeline.py --input-bam " + aligned_filename + \
                        " --reference-fasta " + " Coronavirus_with_expr_ctrl1/Coronavirus_with_expr_ctrl1.fasta " + \
                        " --output-dir variantCall2 " + \
                        " --parameters-file " + parameters + \
                        " --region-bed " + target_bed + \
                        " --bin-dir /results/plugins/variantCaller/bin " + \
                        " --num-threads 12"
            
            call_results = check_output(call_cmd, shell=True, cwd=results_dir)
            
            ## refine the reference again using major alleles to get a consensus genome
            
            # extract 2019-nCoV consensus sequence in the 1st round
            consensus1 = results_dir + "/variantCall1/TSVC_variants.consensus.2019-nCoV.fasta"
            with open(results_dir + "/variantCall1/TSVC_variants.consensus.fasta", "rU") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.id == "2019-nCoV":
                        SeqIO.write(record, consensus1, "fasta")
                        break
            
            refine_cmd = "cat " + consensus1 + " | vcf-consensus variantCall2/TSVC_variants.vcf.gz > variantCall2/TSVC_variants.consensus.fasta"
            refine_results = check_output(refine_cmd, shell=True, cwd=results_dir)
            
            # Check the depth of each base, ensure every base's depth is >= cutoff
            
            mapQ = self.startplugin_json['pluginconfig']['mapQ']
            minCov = self.startplugin_json['pluginconfig']['minCov']
            
            depth_cmd = "samtools depth -Q " + mapQ + " -b" + target_bed + " " + aligned_filename + " | awk -v FS=\"\t\" -v OFS=\"\t\" -v depth=" + minCov + " '{if($3<depth) print $1,$2-1,$2;}' > " + barcode_name + ".lowDepth.bed"
            depth_results = check_output(depth_cmd, shell=True, cwd=results_dir)
            
            # Mask the consensus genome with "N" at bases where coverage is below the minCov parameter
            if os.path.getsize(results_dir + "/" + barcode_name + ".lowDepth.bed") == 0:
                mask_results = check_output(["cp","variantCall2/TSVC_variants.consensus.fasta", "variantCall2/TSVC_variants.consensus.masked.fasta"
                                            ], cwd=results_dir)
            else:
                mask_results = check_output(["bedtools", "maskfasta", "-fi", "variantCall2/TSVC_variants.consensus.fasta", 
                                        "-fo", "variantCall2/TSVC_variants.consensus.masked.fasta",
                                        "-bed", barcode_name + ".lowDepth.bed"
                                            ], cwd=results_dir)
            
            # Note the panel covers 42-29842 (1-29903), so delete extra bases on the two ends ([1,41], [length-60,length])
            consen_seq = SeqIO.read(results_dir + "/variantCall2/TSVC_variants.consensus.masked.fasta","fasta")
            length = len(consen_seq.seq)
            slice_seq = str(consen_seq.seq)[41:(length-61)]
            consen_seq.seq = Seq(slice_seq)
            consen_seq.id = barcode_name + "-consensus"
            consensus_output = os.path.join(results_dir, barcode_name + "-consensus.fasta")
            SeqIO.write(consen_seq, consensus_output, "fasta")
            
            # delete intermediate files, and more results to the filefolder
            os.remove(os.path.join(results_dir, barcode_name + '.1.bam'))
            os.remove(os.path.join(results_dir, barcode_name + '.1.bam.bai'))
            shutil.rmtree(os.path.join(results_dir,"Coronavirus_with_expr_ctrl1"))
            
            os.mkdir(os.path.join(results_dir,barcode_name))
            shutil.move(os.path.join(results_dir,barcode_name + '.2.bam'),os.path.join(results_dir,barcode_name))
            shutil.move(os.path.join(results_dir,barcode_name + '.2.bam.bai'),os.path.join(results_dir,barcode_name))
            shutil.move(os.path.join(results_dir,"variantCall1"),os.path.join(results_dir,barcode_name))
            shutil.move(os.path.join(results_dir,"variantCall2"),os.path.join(results_dir,barcode_name))
            shutil.move(consensus_output,os.path.join(results_dir,barcode_name))
            shutil.move(os.path.join(results_dir,barcode_name + ".lowDepth.bed"),os.path.join(results_dir,barcode_name))
            
            data_entry = {}
            data_entry['barcode'] = barcode_name
            data_entry['sample'] = barcode_values['sample']
            data_entry['assembly'] = os.path.join(barcode_name,barcode_name + "-consensus.fasta")
            
            Assembly_report.append(data_entry)
            
        ###################################################################################################################
        ## output in HTML
        ###################################################################################################################
        
        render_context = {
            "autorefresh" : False,
            "assemblyData" : Assembly_report
        }
        
        createReport(os.path.join(results_dir,'assembly_report.html'), 'barcode_summary_all.html', render_context )
        
        return True
    
    # Return list of columns you want the plugin table UI to show.
    # Columns will be displayed in the order listed.
    def barcodetable_columns(self):
        return [
            { "field": "selected", "editable": True },
            { "field": "barcode_name", "editable": False },
            { "field": "sample", "editable": False } ]

# Devel use - running directly
if __name__ == "__main__":
    PluginCLI()
