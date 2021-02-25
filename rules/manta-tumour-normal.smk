"""
Tumour normal SV calling pipeline using manta
"""

from pathlib import Path

################################################################################
## Prelim setup
################################################################################
with open(config['donor_list']) as f:
    donors = [x.rstrip() for x in f.readlines()]
manifest_dir = config['manifest_dir']
outdir = config['outdir']



################################################################################
## Rules
################################################################################

rule all:
    input:
        expand(manifest_dir+'/{donor}-tumour-normal.tsv', donor=donors),
        expand(outdir+"/{donor}/results/variants/diploidSV.vcf.gz", donor=donors),
        expand(outdir+"/{donor}/results/variants/somaticSV.vcf.gz", donor=donors),
        expand(outdir+"/{donor}/results/variants/candidateSV.vcf.gz", donor=donors),
        expand(outdir+"/{donor}/results/variants/candidateSmallIndels.vcf.gz", donor=donors)
        expand(outdir+"/{donor}/job-finished.txt", donor=donors)

rule GetReference:
    output:
        fasta = temp(outdir+"/ref/hs37d5.fa"),
        fai = temp(outdir+"/ref/hs37d5.fa.fai"),
    shell:
        """
        aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa {output.fasta}
        aws s3 cp s3://layerlabcu/ref/genomes/hs37d5/hs37d5.fa.fai {output.fai}
        """

rule GetBams:
    """
    Download tumour normal pair.  Rename the files for better
    rule interoperability.
    """
    resources:
        num_downloads = 1
    params: 
        odir = outdir
    threads: 1
    input:
        manifest = manifest_dir+'/{donor}-tumour-normal.tsv'
    output:
        tumour_bam = temp(outdir+"/{donor}/bam/{donor}-tumour.bam"),
        tumour_bai = temp(outdir+"/{donor}/bam/{donor}-tumour.bam.bai"),
        normal_bam = temp(outdir+"/{donor}/bam/{donor}-normal.bam"),
        normal_bai = temp(outdir+"/{donor}/bam/{donor}-normal.bam.bai")
    shell:
        # Download bams
        # Rename bams with specimen type and donor id
        # Normal bam is the first entry, tumour is second
        """
        while [[ $(find {params.odir} -name '*.bam' | wc -l) -ge {config[max_bams]} ]]; do
            sleep 5
        done

        outdir=$(dirname {output.tumour_bam})
        score-client download \
            --validate false \
            --output-dir $outdir \
            --manifest {input.manifest}

        normal=$(sed '2q;d' {input.manifest} | cut -f5)
        tumour=$(sed '3q;d' {input.manifest} | cut -f5)

        mv $outdir/$normal {output.normal_bam}
        mv $outdir/$normal.bai {output.normal_bai}
        mv $outdir/$tumour {output.tumour_bam}
        mv $outdir/$tumour.bai {output.tumour_bai}
        """


rule RunManta:
    """
    For a given donor, run the config and SV calling steps of manta
    in Tumour/Normal mode.
    """
    input:
        manta_install_path = config['manta_install_path'],
        normal_bam = outdir+"/{donor}/bam/{donor}-normal.bam",
        normal_bai = outdir+"/{donor}/bam/{donor}-normal.bam.bai",
        tumour_bam = outdir+"/{donor}/bam/{donor}-tumour.bam",
        tumour_bai = outdir+"/{donor}/bam/{donor}-tumour.bam.bai",
        fasta = rules.GetReference.output.fasta,
        fai = rules.GetReference.output.fai,
    params:
        runDir = outdir+"/{donor}"
    output:
        outdir+"/{donor}/results/variants/diploidSV.vcf.gz",
        outdir+"/{donor}/results/variants/somaticSV.vcf.gz",
        outdir+"/{donor}/results/variants/candidateSV.vcf.gz",
        outdir+"/{donor}/results/variants/candidateSmallIndels.vcf.gz"
    resources:
        manta_running = 1,
    threads:
        workflow.cores - 1
    shell:
        """
        {input.manta_install_path}/bin/configManta.py \
            --normalBam {input.normal_bam} \
            --tumorBam {input.tumour_bam} \
            --referenceFasta {input.fasta} \
            --runDir {params.runDir}
        {params.runDir}/runWorkflow.py -j {threads}
        """

rule UploadResults:
    input:
        diploidSV = outdir+"/{donor}/results/variants/diploidSV.vcf.gz",
        somaticSV = outdir+"/{donor}/results/variants/somaticSV.vcf.gz",
        candidateSV = outdir+"/{donor}/results/variants/candidateSV.vcf.gz",
        candidateSmallIndels = outdir+"/{donor}/results/variants/candidateSmallIndels.vcf.gz"
    output:
        outdir+"/{donor}/job-finished.txt"
    shell:
        """
        aws s3 cp {input.diploidSV} s3://layerlabcu/icgc/manta/{donor}/diploidSV.vcf.gz
        aws s3 cp {input.somaticSV} s3://layerlabcu/icgc/manta/{donor}/somaticSV.vcf.gz
        aws s3 cp {input.candidateSV} s3://layerlabcu/icgc/manta/{donor}/candidateSV.vcf.gz
        aws s3 cp {input.candidateSmallIndels} s3://layerlabcu/icgc/manta/{donor}/candidateSmallIndels.vcf.gz
        
        touch {output}
        """
        
