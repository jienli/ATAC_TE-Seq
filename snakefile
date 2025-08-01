# Snakefile

############################################################
#  CONFIGURATION
############################################################
import os
with open("/oscar/data/jsedivy/jienli/Hayflick/atac_analysis_H1/ATAC_TE-Seq/cwd_debug.txt", "w") as f:
    f.write("Snakemake working directory at parse time: " + os.getcwd() + "\n")
os.chdir("/oscar/data/jsedivy/jienli/Hayflick/atac_analysis_H1/ATAC_TE-Seq")

# Load YAML config:
# import os
# configfile: os.path.join(os.path.dirname(__file__), "config.yaml")
# configfile: "config.yaml"
configfile: "/oscar/data/jsedivy/jienli/Hayflick/atac_analysis_H1/ATAC_TE-Seq/config.yaml"


# for sample in config["samples"]:
#     print("Parsed sample:", sample)

# Samples should be defined in config.yaml, e.g.:
# samples:
#   rep1:
#     fastq: "fastq/rep1.SE.fastq.gz"
#     paired: false
#   rep2:
#     fastq1: "fastq/rep2.R1.fastq.gz"
#     fastq2: "fastq/rep2.R2.fastq.gz"
#     paired: true

SAMPLES = config["samples"]
SAMPLES_flat = [s for group in SAMPLES.values() for s in group]    # flattens the list of samples to 1D list of samples

SAMPLES_FASTQ_DIR = config["sample_fastq_dir"]

BWT2_IDX    = config["bwt2_index"]
GENOME_SIZE = config["genome_size"]
CHRSIZES    = config["chrom_sizes"]
# BLACKLIST   = config["blacklist"]
# REF_FA      = config["reference_fasta"]
# TSS_BED     = config["tss_bed"]
# DNASE_BED   = config["dnase_bed"]
# ENH_BED     = config["enhancer_bed"]
# PROM_BED    = config["promoter_bed"]
SUBSAMPLE   = config.get("subsample_reads", 25000000)
THREADS     = config.get("threads", 8)
MEM_MB      = config.get("mem_mb", 8000)
IDR_THRESH  = config.get("idr_threshold", 0.05)
SMOOTH_WIN  = config.get("smooth_window", 150)

ADAPTOR_ERR_RATE = config.get("adaptor_err_rate", 0.2)
MULTIMAPPING = config.get("multimapping", 100)

MACS2_PVAL_THRESHOLD = config.get("macs2_pval_threshold", 1e-2)
############################################################
#  RULE: aggregate final targets
############################################################

# rule all:
#     input:
#         # flagstat QC for each rep
#         expand("qc/{sample}.flagstat.qc", sample=SAMPLES),
#         # peak calls (narrowPeak) for each rep
#         expand("peaks/{sample}.narrowPeak.gz", sample=SAMPLES),
#         # IDR‐filtered peaks (conservative) for pooled
#         "idr/pooled.conservative.narrowPeak.gz",
#         # bigWig tracks
#         expand("tracks/{sample}.positive.bigwig", sample=SAMPLES),
#         expand("tracks/{sample}.negative.bigwig", sample=SAMPLES),
#         # TSS enrichment plots
#         expand("tss/{sample}.tss_enrich.png", sample=SAMPLES)


# actual all
# rule all:
#     input:
#         # flagstat QC for each rep
#         expand("qc/{sample}.PE.flagstat.qc", sample=SAMPLES_flat),
#         # peak calls (narrowPeak) for each rep
#         expand("peaks/{sample}.narrowPeak.gz", sample=SAMPLES_flat)


# # all except 92 (post PE.bam for now)
# SAMPLES_flat_ = ["SRR14646390", "SRR14646391", "SRR14646392", "SRR14646405", "SRR14646406"]
# rule all:
#     input:
#         # flagstat QC for each rep
#         expand("qc/{sample}.PE.flagstat.qc", sample=SAMPLES_flat_),
#         # peak calls (narrowPeak) for each rep
#         expand("peaks/{sample}.narrowPeak.gz", sample=SAMPLES_flat_),
#         expand("filtered/{sample}.PE.final.bam.bai", sample=SAMPLES_flat_),
#         expand("tracks/{sample}.fc.signal.bw", sample=SAMPLES_flat_),
#         expand("tracks/{sample}.pval.signal.bw", sample=SAMPLES_flat_),
#         "diff/deseq2_results.csv"


# # only 92: still at bowtie2 for now
# SAMPLES_flat__ = ["SRR14646392"]
# rule all_92:
#     input:
#         # flagstat QC for each rep
#         expand("qc/{sample}.PE.flagstat.qc", sample=SAMPLES_flat__),
#         # peak calls (narrowPeak) for each rep
#         expand("peaks/{sample}.narrowPeak.gz", sample=SAMPLES_flat__),
#         expand("filtered/{sample}.PE.final.bam.bai", sample=SAMPLES_flat__),
#         expand("tracks/{sample}.fc.signal.bw", sample=SAMPLES_flat__),
#         expand("tracks/{sample}.pval.signal.bw", sample=SAMPLES_flat__)



# all except 92 (post PE.bam for now)
SAMPLES_flat_ = SAMPLES_flat
rule all:
    input:
        # flagstat QC for each rep
        expand("qc/{sample}.PE.flagstat.qc", sample=SAMPLES_flat_),
        # peak calls (narrowPeak) for each rep
        expand("peaks/{sample}.narrowPeak.gz", sample=SAMPLES_flat_),
        expand("filtered/{sample}.PE.final.bam.bai", sample=SAMPLES_flat_),
        expand("tracks/{sample}.fc.signal.bw", sample=SAMPLES_flat_),
        expand("tracks/{sample}.pval.signal.bw", sample=SAMPLES_flat_),
        "analysis/counts.tsv",
        "analysis/merged_results_count.csv",
        expand("idr/{condition}_qc.txt", condition=SAMPLES)

############################################################
#  0a. Adapter detection                              
############################################################

#if paried end reads, use the first read
# otherwise use the single end read
rule detect_adapter:
    input:
        # fastq = lambda wildcards: (
        #     SAMPLES[wildcards.sample]["fastq"]
        #     if not SAMPLES[wildcards.sample].get("paired", False)
        #     else SAMPLES[wildcards.sample]["fastq1"]
        # )
        fastq = lambda wc: os.path.join(SAMPLES_FASTQ_DIR, f"{wc.sample}_1.fastq.gz")
    output:
        "adapters/{sample}.adapter.txt"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    shell:
        """
        mkdir -p adapters
        python3 {config[script_dir]}/detect_adapter.py {input.fastq} > {output}
        """


############################################################
#  0b. Adapter trimming                                
############################################################

rule trim_adapters_SE:
    # prevent wildcard matching for .R1 and .R2 suffixes (which should go to trim_adapters_PE)
    wildcard_constraints:
        sample=r"(?!.*\.R[12]).+"
    input:
        fastq = lambda wc: os.path.join(SAMPLES_FASTQ_DIR, f"{wc.sample}.fastq.gz"),
        adapter = rules.detect_adapter.output
    params:
        adaptor_err_rate = ADAPTOR_ERR_RATE
    output:
        "trimmed/{sample}.trim.fastq.gz"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    shell:
        """
        mkdir -p trimmed
        cutadapt -m 5 -e {params.adaptor_err_rate} \
            -a $(cat {input.adapter}) --cores {threads} --output {output} {input.fastq}
        """

rule trim_adapters_PE:
    input:
        fastq1 = lambda wc: os.path.join(SAMPLES_FASTQ_DIR, f"{wc.sample}_1.fastq.gz"),
        fastq2 = lambda wc: os.path.join(SAMPLES_FASTQ_DIR, f"{wc.sample}_2.fastq.gz"),
        # adapter = rules.detect_adapter.output
        adapter = lambda wc: f"adapters/{wc.sample}.adapter.txt"
    params:
        adaptor_err_rate = ADAPTOR_ERR_RATE
    output:
        trimmed_r1 = "trimmed/{sample}.R1.trim.fastq.gz",
        trimmed_r2 = "trimmed/{sample}.R2.trim.fastq.gz"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    shell:
        """
        mkdir -p trimmed
        cutadapt -m 5 -e {params.adaptor_err_rate} \
            -a $(cat {input.adapter}) -A $(cat {input.adapter}) \
            -o {output.trimmed_r1} -p {output.trimmed_r2} \
            --cores {threads} \
            {input.fastq1} {input.fastq2}
        """


############################################################
#  1a. Alignment (SE / PE)                        
############################################################

# rule bowtie2_align_SE:
#     input:
#         fastq = "trimmed/{sample}.trim.fastq.gz"
#     output:
#         bam = "aligned/{sample}.SE.bam",
#         flagstat = "qc/{sample}.SE.flagstat.qc",
#         non_mito = "aligned/{sample}.SE.non_mito.bam"
#     params:
#         multimapping = MULTIMAPPING,
#         bwt2_index = BWT2_IDX
#     threads: THREADS
#     resources:
#         mem_mb=MEM_MB
#     conda:
#         "ATAC_Core"
#     shell:
#         """
#         mkdir -p aligned qc
#         bowtie2 -k {params.multimapping} --mm -x {params.bwt2_index} --threads {threads} -U <(zcat -f {input.fastq}) 2> aligned/{wildcards.sample}.align.log \
#           | samtools view -Su - \
#           | samtools sort -@ {threads} -o {output.bam}

#         # flagstat QC
#         samtools sort -n -@ {threads} {output.bam} -O SAM \
#           | SAMstats --sorted_sam_file - --outf {output.flagstat}

#         # extract non‐mito
#         samtools idxstats {output.bam} \
#           | cut -f 1 \
#           | grep -v -P "^chrM$" \
#           | xargs samtools view {output.bam} -@ {threads} -b \
#           > {output.non_mito}
#         """




# New: Bowtie2 mapping PE reads to unsorted BAM
rule bowtie2_map_PE:
    input:
        fastq1 = "trimmed/{sample}.R1.trim.fastq.gz",
        fastq2 = "trimmed/{sample}.R2.trim.fastq.gz"
    output:
        bam_unsorted = temp("aligned/{sample}.unsorted.bam")
    params:
        multimapping = MULTIMAPPING,
        bwt2_index = BWT2_IDX
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        mkdir -p aligned
        
        bowtie2 -q -k {params.multimapping} -X2000 --mm -x {params.bwt2_index} --threads {threads} \
          -1 {input.fastq1} -2 {input.fastq2} 2> aligned/{wildcards.sample}.align.log \
          | samtools view -Su - > {output.bam_unsorted}
        """

# New: Sort BAM
rule sort_bam_PE:
    input:
        bam_unsorted = "aligned/{sample}.unsorted.bam"
    output:
        bam = "aligned/{sample}.PE.bam"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        m_per_thread=$(( ({resources.mem_mb} * 8 / 10) / {threads} ))M

        samtools sort -@ {threads} -m $m_per_thread -o {output.bam} {input.bam_unsorted}
        """


# New: Flagstat QC
rule flagstat_qc_PE:
    input:
        bam = "aligned/{sample}.PE.bam"
    output:
        flagstat = "qc/{sample}.PE.flagstat.qc"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        m_per_thread=$(( ({resources.mem_mb} * 7 / 10) / {threads} ))M
        mkdir -p qc
        samtools sort -n -@ {threads} -m $m_per_thread {input.bam} -O SAM \
          | SAMstats --sorted_sam_file - --outf {output.flagstat}
        """

# New: Extract non-mito reads
rule extract_non_mito_PE:
    input:
        bam = "aligned/{sample}.PE.bam"
    output:
        non_mito = "aligned/{sample}.PE.non_mito.bam"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        samtools view -@ {threads} -b {input.bam} $(samtools idxstats {input.bam} | cut -f 1 | grep -v "^chrM$") > {output.non_mito}
        """




        
############################################################
#  1c. Fraction mitochondrial reads                    
############################################################

rule frac_mito:
    input:
        non_mito = lambda wc: (
            f"aligned/{wc.sample}.PE.non_mito.bam" if SAMPLES[wc.sample]["paired"]
            else f"aligned/{wc.sample}.SE.non_mito.bam"
        ),
        mito = lambda wc: (
            f"aligned/{wc.sample}.PE.bam" if SAMPLES[wc.sample]["paired"]
            else f"aligned/{wc.sample}.SE.bam"
        )
    output:
        "qc/{sample}.frac_mito.txt"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        mkdir -p qc
        samtools sort -n -@ {threads} {input.non_mito} -O SAM \
          | SAMstats --sorted_sam_file - --outf qc/{wildcards.sample}.nonmito.qc

        samtools sort -n -@ {threads} {input.mito} -O SAM \
          | SAMstats --sorted_sam_file - --outf qc/{wildcards.sample}.mito.qc

        Rn=$(awk '$2=="mapped"{sum+=$3}END{{print sum}}' qc/{wildcards.sample}.nonmito.qc)
        Rm=$(awk '$2=="mapped"{sum+=$3}END{{print sum}}' qc/{wildcards.sample}.mito.qc)
        echo -e "$Rm\t$Rn\t$(echo "$Rm/($Rm+$Rn)" | bc -l)" > {output}
        """


############################################################
#  1b. Post‐alignment filtering (dedup + MAPQ)          
############################################################

# rule dedup_filter_PE:
#     input:
#         bam = "aligned/{sample}.PE.bam"
#     output:
#         final_bam = "filtered/{sample}.PE.final.bam",
#         dup_metrics = "qc/{sample}.PE.dup.qc",
#         pbc = "qc/{sample}.PE.pbc.qc"
#     threads: THREADS
#     params:
#         MAPQ_THRESH = 20
#     resources:
#         mem_mb=MEM_MB
#     conda:
#         "ATAC_Core"
#     shell: 
#         """
#         m_per_thread=$(( ({resources.mem_mb} * 8 / 10) / {threads} ))M

#         # Remove unmapped, mate unmapped, QC-fail, keep properly paired (-f 2)
#         samtools view -F 524 -f 2 -q {params.MAPQ_THRESH} -u {input.bam} \
#             | samtools sort -n -@ {threads} -m $m_per_thread -o filtered/{wildcards.sample}.nmsrt.bam


#         # Assign multimappers (remove all that maps to more than 1000 loci and keeps the rests the same) and fixmate
#         samtools view -h filtered/{wildcards.sample}.nmsrt.bam \
#             | python3 {config[script_dir]}/assign_multimappers.py -k {MULTIMAPPING} --paired-end \
#             | samtools fixmate - filtered/{wildcards.sample}.fixmate.bam

#         # ====================================================================    
#         # NOTE: This is a major change from the original script, which used fixmate -r, and removes all secondary reads.
#         # ====================================================================    

#         # Remove orphan pairs and discordant mapping (purpose of the original code from ENCODE)
#         # But here, primarily this just sorts the reads by coordinate (previousely sorted by name)

#         samtools view -F 524 -f 2 -u filtered/{wildcards.sample}.fixmate.bam \
#             | samtools sort -@ {threads} -m $m_per_thread -o filtered/{wildcards.sample}.filt.bam
#         rm filtered/{wildcards.sample}.nmsrt.bam filtered/{wildcards.sample}.fixmate.bam

        
#         # Mark PCR duplicates (only marks the FLAG, does not remove them)
#         java -Xmx4G -jar {config[picard]} MarkDuplicates \
#             INPUT=filtered/{wildcards.sample}.filt.bam \
#             OUTPUT=filtered/{wildcards.sample}.dupmark.bam \
#             METRICS_FILE={output.dup_metrics} \
#             VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true \
#             REMOVE_DUPLICATES=false
#         mv filtered/{wildcards.sample}.dupmark.bam filtered/{wildcards.sample}.filt.bam

#         # Remove duplicates and save final BAM
#         # ====================================================================    
#         # NOTE: This is a major change from the original script, which used -F 1804, and removes all PCR duplicates.
#         # ====================================================================    
#         # In our version, this only compress the bam (-u --> -b)
#         samtools view -F 524 -f 2 -b filtered/{wildcards.sample}.filt.bam \
#             > {output.final_bam}
#         # TODO consolidate this

#         samtools index {output.final_bam}

#         # Flagstat QC
#         samtools sort -n --threads {threads} -m $m_per_thread {output.final_bam} -O SAM \
#             | SAMstats --sorted_sam_file - --outf {output.final_bam}.flagstat.qc

#         # Library complexity (PBC)
#         # mt = Total number of read pairs
#         # m0 = Distinct read pairs (unique fragments)
#         # m1 = Fragments seen once (singleton)
#         # m2 = Fragments seen twice
#         # m0/mt: NRF = Non-redundant fraction
#         # m1/m0: PBC1 = Proportion of unique among distinct
#         # m1/m2: PBC2 = Singleton-to-doubleton ratio
#         samtools sort -n -@ {threads} -m $m_per_thread filtered/{wildcards.sample}.filt.bam -o filtered/{wildcards.sample}.srt.tmp.bam
#         bedtools bamtobed -bedpe -i filtered/{wildcards.sample}.srt.tmp.bam \
#             | awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$4,$6,$9,$10}}' \
#             | grep -v 'chrM' \
#             | sort | uniq -c \
#             | awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0}} \
#             ($1==1){{m1+=1}} ($1==2){{m2+=1}} {{m0+=1;mt+=$1}} \
#             END{{printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}}' \
#             > {output.pbc}
#         rm filtered/{wildcards.sample}.srt.tmp.bam filtered/{wildcards.sample}.filt.bam
#         """

    
# # Step 1: Filter and name-sort paired reads
# rule filter_pairs_PE_dedup1:
#     input:
#         bam = "aligned/{sample}.PE.bam"
#     output:
#         nmsrt = "filtered/{sample}.nmsrt.bam"
#     threads: THREADS
#     resources:
#         mem_mb = MEM_MB
#     params:
#         MAPQ_THRESH = 20
#     conda:
#         "ATAC_Core"
#     shell:
#         """
#         m_per_thread=$(( ({resources.mem_mb} * 7 / 10) / {threads} ))M
#         samtools view -F 524 -f 2 -u {input.bam} \
#             | samtools sort -n -@ {threads} -m $m_per_thread -o {output.nmsrt}
#         """

# # Step 2: Assign multimappers and fixmate
# rule assign_multimappers_PE_dedup2:
#     input:
#         nmsrt = rules.filter_pairs_PE_dedup1.output.nmsrt
#     output:
#         fixmate = "filtered/{sample}.fixmate.bam"
#     threads: THREADS
#     conda:
#         "ATAC_Core"
#     shell:
#         """
#         # samtools view -h {input.nmsrt} \
#         #     | python3 {config[my_script_dir]}/assign_multimappers.py -k 1 --paired-end \
#         #     | samtools fixmate - {output.fixmate}
#         samtools view -h -F 780 {input.nmsrt} \
#             | samtools fixmate - {output.fixmate}
#         """
#         # This assign_multimappers.py script is a custom script that caps multimappers to 1 loci by randomly assigning one of the loci to the read.

# # Step 3: Remove orphans and coordinate-sort
# rule sort_coord_PE_dedup3:
#     input:
#         nmsrt = rules.filter_pairs_PE_dedup1.output.nmsrt,
#         fixmate = rules.assign_multimappers_PE_dedup2.output.fixmate
#     output:
#         filt = "filtered/{sample}.filt.bam"
#     threads: THREADS
#     resources:
#         mem_mb = MEM_MB
#     conda:
#         "ATAC_Core"
#     shell:
#         """
#         m_per_thread=$(( ({resources.mem_mb} * 7 / 10) / {threads} ))M
#         samtools view -F 524 -f 2 -u {input.fixmate} \
#             | samtools sort -@ {threads} -m $m_per_thread -o {output.filt}
#         rm {input.nmsrt} {input.fixmate}
#         """

# # Step 4: Mark PCR duplicates
# rule mark_duplicates_PE_dedup4:
#     input:
#         filt = rules.sort_coord_PE_dedup3.output.filt
#     output:
#         dupmark = "filtered/{sample}.dupmark.bam",
#         dup_metrics = "qc/{sample}.PE.dup.qc"
#     threads: THREADS
#     conda:
#         "ATAC_Core"
#     shell:
#         """
#         java -Xmx4G -jar {config[picard]} MarkDuplicates \
#             INPUT={input.filt} OUTPUT={output.dupmark} METRICS_FILE={output.dup_metrics} \
#             VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
#         """
#         # picard MarkDuplicates requires position sorted BAM, so we need to sort the input BAM by name first

# # Step 5: Filter duplicates to final BAM
# rule filter_dups_PE_dedup5:
#     input:
#         dupmark = rules.mark_duplicates_PE_dedup4.output.dupmark
#     output:
#         final_bam = "filtered/{sample}.PE.final.bam"
#     threads: THREADS
#     resources:
#         mem_mb = MEM_MB
#     conda:
#         "ATAC_Core"
#     shell:
#         """
#         m_per_thread=$(( ({resources.mem_mb} * 7 / 10) / {threads} ))M

#         # samtools view -F 1548 -f 2 -b {input.dupmark} > {output.final_bam}
#         samtools view -F 1548 -f 2 -u {input.dupmark} | samtools sort -@ {threads} -m $m_per_thread -o {output.final_bam}
#         rm {input.dupmark}
#         """


# Step 1: Filter for good read pairs and mark PCR duplicates with picard
rule filter_markDup_PE_dedup1:
    input:
        bam = "aligned/{sample}.PE.bam"
    output:
        dupmark = temp("filtered/{sample}.dupmark.bam"),
        dup_metrics = "qc/{sample}.PE.dup.qc",
        filt_1 = temp("filtered/{sample}.filt_1.bam")
    threads: THREADS
    resources:
        mem_mb = MEM_MB
    params:
        MAPQ_THRESH = 20
    conda:
        "ATAC_Core"
    shell:
        """
        samtools view -F 524 -f 2 -u {input.bam} -o {output.filt_1}

        java -Xmx4G -jar {config[picard]} MarkDuplicates \
            INPUT={output.filt_1} OUTPUT={output.dupmark} METRICS_FILE={output.dup_metrics} \
            VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
        """


# Step 2: (optional filter PCR duplicates) and name-sort paired reads
rule filter_nmsrt_PE_dedup2:
    input:
        dupmark = rules.filter_markDup_PE_dedup1.output.dupmark,
    output:
        nmsrt = temp("filtered/{sample}.nmsrt.bam")
    threads: THREADS
    resources:
        mem_mb = MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        m_per_thread=$(( ({resources.mem_mb} * 7 / 10) / {threads} ))M
        samtools view -F 1548 -f 2 -u {input.dupmark} \
            | samtools sort -n -@ {threads} -m $m_per_thread -o {output.nmsrt}

        # # Use following lines instead if you want to filter out Duplicates after multimapping assignment
        # samtools view -F 524 -f 2 -u {input.dupmark} \
        #     | samtools sort -n -@ {threads} -m $m_per_thread -o {output.nmsrt}
        """


# Step 3: Assign multimappers and fixmate
rule assign_multimappers_PE_dedup3:
    input:
        nmsrt = rules.filter_nmsrt_PE_dedup2.output.nmsrt
    output:
        multi_assigned = "filtered/{sample}.multi_assigned.bam"
    threads: THREADS
    conda:
        "ATAC_Core"
    shell:
        """
        # This assign_multimappers.py script is a custom script that caps multimappers to 1 loci by randomly assigning one of the loci to the read.
        # samtools view -h {input.nmsrt} \
        #     | python3 {config[my_script_dir]}/assign_multimappers.py -k 1 --paired-end \
        #     | samtools fixmate - {output.multi_assigned}

        # NOTE: This keeps only primary alignments and removes secondary alignments
        samtools view -h -F 780 {input.nmsrt} \
            | samtools fixmate - {output.multi_assigned}
        """


# Step 4: Coordinate Sort to final BAM
rule final_bam_PE_dedup4:
    input:
        multi_assigned = rules.assign_multimappers_PE_dedup3.output.multi_assigned
    output:
        final_bam = "filtered/{sample}.PE.final.bam"
    threads: THREADS
    resources:
        mem_mb = MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        m_per_thread=$(( ({resources.mem_mb} * 7 / 10) / {threads} ))M
        samtools view -u {input.multi_assigned} | samtools sort -@ {threads} -m $m_per_thread -o {output.final_bam}
        rm {input.multi_assigned}
        """

# Step 65: Index final BAM
rule index_finalBam_PE:
    input:
        bam = rules.final_bam_PE_dedup4.output.final_bam
    output:
        bai = "filtered/{sample}.PE.final.bam.bai"
    threads: THREADS
    resources:
        mem_mb = MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        samtools index {input.bam}
        """

# Step 6: Flagstat QC on deduplicated BAM
rule flagstat_dedup_PE:
    input:
        bam = rules.final_bam_PE_dedup4.output.final_bam
    output:
        flagstat = "qc/{sample}.PE.flagstat_dedup.qc"
    threads: THREADS
    resources:
        mem_mb = MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        m_per_thread=$(( ({resources.mem_mb} * 7 / 10) / {threads} ))M
        samtools sort -n -@ {threads} -m $m_per_thread {input.bam} -O SAM \
          | SAMstats --sorted_sam_file - --outf {output.flagstat}
        """

# Step 7: Compute library complexity (PBC)
rule pbc_LibComplexity_PE:
    input:
        bam = rules.final_bam_PE_dedup4.output.final_bam
    output:
        pbc = "qc/{sample}.PE.pbc.qc"
    threads: THREADS
    resources:
        mem_mb = MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        m_per_thread=$(( ({resources.mem_mb} * 7 / 10) / {threads} ))M
        samtools sort -n -@ {threads} -m $m_per_thread {input.bam} -o filtered/{wildcards.sample}.PE.srt.tmp.bam
        bedtools bamtobed -bedpe -i filtered/{wildcards.sample}.PE.srt.tmp.bam \
            | awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$4,$6,$9,$10}}' \
            | grep -v 'chrM' \
            | sort | uniq -c \
            | awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0}} \
                ($1==1){{m1+=1}} ($1==2){{m2+=1}} {{m0+=1;mt+=$1}} \
                END{{printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}}' \
            > {output.pbc}
        rm filtered/{wildcards.sample}.PE.srt.tmp.bam
        """










############################################################
#  2a. BAM → tagAlign (SE / PE)                        
############################################################

# rule bam_to_tag:
#     input:
#         bam = rules.dedup_filter.output.final_bam
#     output:
#         tag = "tagAlign/{sample}.SE.tagAlign.gz"
#     threads: 1
#     shell:
#         """
#         bedtools bamtobed -i {input.bam} \
#           | awk 'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print}}' \
#           | gzip -nc \
#           > {output.tag}
#         """


# Note: Change from Original: This rule now always filters out chrM reads
############################################################
#  TagAlign File Format: (separates paried reads)
#  chr#   start    end     name    score   strand
#  chr1    1000    1001    N       1000    +
#  chr1    1002    1003    N       1000    -
############################################################
rule bam_to_tag_PE:
    input:
        bam = rules.final_bam_PE_dedup4.output.final_bam
    output:
        bedpe = "bedpe/{sample}.bedpe.gz",
        tag   = "tagAlign/{sample}.PE.tagAlign.gz",
        tag_subsample = "tagAlign/{sample}.PE.subsample.tagAlign.gz"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        mkdir -p bedpe tagAlign
        # bedtools bamtobed -bedpe -mate1 -i {input.bam} | pigz -p {threads} -nc > {output.bedpe}   # This is wrong, bedtools bamtobed need name sorted bam

        m_per_thread=$(( ({resources.mem_mb} * 7 / 10) / {threads} ))M
        samtools sort -n -@ {threads} -m $m_per_thread -o filtered/{wildcards.sample}.PE.namesorted.bam {input.bam}
        bedtools bamtobed -bedpe -mate1 -i filtered/{wildcards.sample}.PE.namesorted.bam | pigz -p {threads} -nc > {output.bedpe}
        rm filtered/{wildcards.sample}.PE.namesorted.bam

        zcat {output.bedpe} \
        | grep -v “chrM” \
        | awk 'BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n", $1,$2,$3,$9,$4,$5,$6,$10}}' \
        | pigz -p {threads} -nc > {output.tag}

        zcat {output.bedpe} \
        | grep -v “chrM” \
        | awk 'BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n", $1,$2,$3,$9,$4,$5,$6,$10}}' \
        | shuf -n {SUBSAMPLE} --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f {output.bedpe} | wc -c) -nosalt </dev/zero 2>/dev/null) \
        | pigz -p {threads} -nc > {output.tag_subsample}
        """


############################################################
#  2b. Cross‐correlation QC                          
############################################################

# rule xcor:
#     input:
#         tag = rules.bam_to_tag.output.tag
#     output:
#         qc   = "qc/{sample}.xcor.qc",
#         plot = "qc/{sample}.xcor.plot.pdf"
#     threads: THREADS
#     shell:
#         """
#         Rscript run_spp.R \
#           -c={input.tag} -p={threads} -filtchr=chrM \
#           -savp={output.plot} -out={output.qc}
#         sed -r 's/,[^\\t]+//g' {output.qc} > temp && mv temp {output.qc}
#         """


############################################################
#  2c. Self‐pseudoreplicates (SE / PE)             
############################################################

# rule spr_SE:
#     input:
#         tag = rules.bam_to_tag.output.tag
#     output:
#         pr1 = "tagAlign/{sample}.SE.pr1.tagAlign.gz",
#         pr2 = "tagAlign/{sample}.SE.pr2.tagAlign.gz"
#     shell:
#         """
#         nlines=$(zcat {input.tag} | wc -l)
#         nlines=$(( (nlines+1)/2 ))
#         zcat {input.tag} | shuf … | split -d -l $nlines {wildcards.sample}.
#         gzip -nc {wildcards.sample}.00 > {output.pr1}
#         gzip -nc {wildcards.sample}.01 > {output.pr2}
#         """


rule spr_PE:
    input:
        tag = rules.bam_to_tag_PE.output.tag
    output:
        pr1 = "tagAlign/{sample}.PE2SE.pr1.tagAlign.gz",
        pr2 = "tagAlign/{sample}.PE2SE.pr2.tagAlign.gz"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        mkdir -p tagAlign
        PR_PREFIX="tagAlign/{wildcards.sample}.filt.nodup"
        PR1_TA_FILE="tagAlign/{wildcards.sample}.PE2SE.pr1.tagAlign.gz"
        PR2_TA_FILE="tagAlign/{wildcards.sample}.PE2SE.pr2.tagAlign.gz"
        joined="tagAlign/{wildcards.sample}.temp.bedpe"

        # Create pseudo-BEDPE by combining every two lines into one
        zcat {input.tag} | sed 'N;s/\\n/\\t/' | gzip -nc > $joined

        # Count total read pairs
        nlines=$(zcat $joined | wc -l)
        nlines=$(( (nlines + 1) / 2 ))

        # Shuffle and split BEDPE into two parts
        zcat $joined | shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f {input.tag} | wc -c) -nosalt </dev/zero 2>/dev/null) \
          | split -d -l $nlines - $PR_PREFIX
        # this will create two files: {{PR_PREFIX}}00 and {{PR_PREFIX}}01

        # Convert each half to tagAlign format
        awk 'BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' ${{PR_PREFIX}}00 \
          | gzip -nc > $PR1_TA_FILE
        rm ${{PR_PREFIX}}00

        awk 'BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' ${{PR_PREFIX}}01 \
          | gzip -nc > $PR2_TA_FILE
        rm ${{PR_PREFIX}}01
        rm $joined
        """


#######
# NOTE: From now on, there will NOT be any PE vs SE distinction. All reads in TagAlign format will be treated as single-end reads. (Which makes sense for ATAC-seq data)
# This is how ENCODE ATAC-seq pipeline works.
#######

############################################################
#  2d. Pooled replicates + pseudoreps               
############################################################

# Note: This rule updates the original to dynamically pool all replicates
# TODO: this is wrong for now, Need to update smaple sheet and use that info to pool the right samples
rule poll_fullTA:
    input:
        full_TA = lambda wc: expand("tagAlign/{sample}.PE.tagAlign.gz", sample=SAMPLES[wc.condition])
    output:
        pooled = "tagAlign/{condition}.pooled.tagAlign.gz"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        # Pool all replicates into a single tagAlign
        zcat {input.full_TA} | gzip -nc > {output.pooled}
        """



rule poll_PRs:
    input:
        pr1    = lambda wc: expand("tagAlign/{sample}.PE2SE.pr1.tagAlign.gz", sample=SAMPLES[wc.condition]),
        pr2    = lambda wc: expand("tagAlign/{sample}.PE2SE.pr2.tagAlign.gz", sample=SAMPLES[wc.condition])
    output:
        ppr1   = "tagAlign/{condition}.pooled.pr1.tagAlign.gz",
        ppr2   = "tagAlign/{condition}.pooled.pr2.tagAlign.gz"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        # Pool first pseudo-replicates
        zcat {input.pr1} | gzip -nc > {output.ppr1}
        
        # Pool second pseudo-replicates
        zcat {input.pr2} | gzip -nc > {output.ppr2}
        """


############################################################
#  2e. Tn5 shifting                                
############################################################

rule tn5_shift:
    # prevent wildcard matching for pooled samples or pr samples
    wildcard_constraints:
        # sample=r"^(?!.*\.pooled)(?!.*\.pr[12]$).+"
        sample=r"(?!.*\.(?:pooled|pr[12])).+"
        # sample = "[^.]+"
        # sample = "(SRR13601448|SRR13601447)"
    input:
        tag = "tagAlign/{sample}.PE.tagAlign.gz"
    output:
        shifted = "tagAlign/{sample}.tn5.tagAlign.gz"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        zcat {input.tag} \
          | awk -F '\\t' 'BEGIN{{OFS=FS}}{{if($6=="+"){{$2+=4}} else if($6=="-"){{$3-=5}} print}}' \
          | gzip -nc > {output.shifted}
        """

rule tn5_shift_pr:
    # prevent wildcard matching for pooled samples
    wildcard_constraints:
        sample=r"(?!.*\.pooled).+"
    input:
        tag1 = "tagAlign/{sample}.PE2SE.pr1.tagAlign.gz",
        tag2 = "tagAlign/{sample}.PE2SE.pr2.tagAlign.gz"
    output:
        shifted1 = "tagAlign/{sample}.pr1.tn5.tagAlign.gz",
        shifted2 = "tagAlign/{sample}.pr2.tn5.tagAlign.gz"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        zcat {input.tag1} \
          | awk -F '\\t' 'BEGIN{{OFS=FS}}{{if($6=="+"){{$2+=4}} else if($6=="-"){{$3-=5}} print}}' \
          | gzip -nc > {output.shifted1}
        zcat {input.tag2} \
          | awk -F '\\t' 'BEGIN{{OFS=FS}}{{if($6=="+"){{$2+=4}} else if($6=="-"){{$3-=5}} print}}' \
          | gzip -nc > {output.shifted2}
        """


rule tn5_shift_pooled:
    input:
        tag = "tagAlign/{condition}.pooled.tagAlign.gz",
        tag1 = "tagAlign/{condition}.pooled.pr1.tagAlign.gz",
        tag2 = "tagAlign/{condition}.pooled.pr2.tagAlign.gz"
    output:
        shifted = "tagAlign/{condition}.pooled.tn5.tagAlign.gz",
        shifted1 = "tagAlign/{condition}.pooled.pr1.tn5.tagAlign.gz",
        shifted2 = "tagAlign/{condition}.pooled.pr2.tn5.tagAlign.gz"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        zcat {input.tag} \
          | awk -F '\\t' 'BEGIN{{OFS=FS}}{{if($6=="+"){{$2+=4}} else if($6=="-"){{$3-=5}} print}}' \
          | gzip -nc > {output.shifted}

        zcat {input.tag1} \
          | awk -F '\\t' 'BEGIN{{OFS=FS}}{{if($6=="+"){{$2+=4}} else if($6=="-"){{$3-=5}} print}}' \
          | gzip -nc > {output.shifted1}

        zcat {input.tag2} \
          | awk -F '\\t' 'BEGIN{{OFS=FS}}{{if($6=="+"){{$2+=4}} else if($6=="-"){{$3-=5}} print}}' \
          | gzip -nc > {output.shifted2}
        """

############################################################
#  2f. JSD (plotFingerprint)                           
############################################################

# rule jsd:
#     input:
#         bam1 = rules.dedup_filter.output.final_bam,
#         bam2 = lambda wc: rules.dedup_filter.output.final_bam.replace(wc.sample, "rep2")
#     output:
#         plot = "qc/{sample}.jsd.png",
#         log  = "qc/{sample}.jsd.log"
#     threads: THREADS
#     shell:
#         """
#         bedtools intersect -v -abam {input.bam1} -b {BLACKLIST} > tmp1.bam
#         bedtools intersect -v -abam {input.bam2} -b {BLACKLIST} > tmp2.bam
#         plotFingerprint -b tmp1.bam tmp2.bam \
#           --outQualityMetrics {output.log} \
#           --plotFile {output.plot} \
#           --minMappingQuality 30 --numberOfProcessors {threads}
#         """


############################################################
#  2g. GC bias                                         
############################################################

# rule gc_bias:
#     input:
#         bam = rules.dedup_filter.output.final_bam
#     output:
#         plot = "qc/{sample}.gc_plot.png",
#         log  = "qc/{sample}.gc.txt"
#     threads: THREADS
#     shell:
#         """
#         java -Xmx6G -jar {config[picard]}/CollectGcBiasMetrics.jar \
#           R={REF_FA} I={input.bam} O={output.log} CHART=tmp.gc.pdf
#         python plot_gc.py {output.log} {wildcards.sample}
#         """


############################################################
#  2h. Fragment‐length stats (PE only)                 
############################################################

# rule fraglen_stat:
#     input:
#         bam = rules.dedup_filter.output.final_bam
#     output:
#         qc   = "qc/{sample}.nucleosomal.qc",
#         plot = "qc/{sample}.fraglen_dist.png"
#     threads: THREADS
#     shell:
#         """
#         java -Xmx6G -jar {config[picard]}/CollectInsertSizeMetrics.jar \
#           INPUT={input.bam} OUTPUT=tmp.inserts.txt H=tmp.inserts.pdf \
#           VERBOSITY=ERROR QUIET=TRUE W=1000 STOP_AFTER=5000000
#         python fraglen_qc.py tmp.inserts.txt {wildcards.sample}
#         """


############################################################
#  3a. Peak calling (MACS2)                          
############################################################

rule macs2_callpeak:
    input:
        tag = "tagAlign/{sample}.tn5.tagAlign.gz"
    output:
        raw_narrow = "macs2_temp/{sample}_peaks.narrowPeak",
        treat_bdg = "macs2_temp/{sample}_treat_pileup.bdg",
        control_bdg = "macs2_temp/{sample}_control_lambda.bdg",
        summits = "macs2_temp/{sample}_summits.bed",
        xls = "macs2_temp/{sample}_peaks.xls"
    params:
        genome_size = GENOME_SIZE,
        pval_thresh = MACS2_PVAL_THRESHOLD,
        smooth_win = SMOOTH_WIN
    threads: 1
    resources:
        mem_mb = 20
    conda:
        "ATAC_Core"
    shell:
        """
        shiftsize=$((-1 * {params.smooth_win} / 2))
        macs2 callpeak -t {input.tag} -f BED -n "macs2_temp/{wildcards.sample}" \
            -g {params.genome_size} -p {params.pval_thresh} \
            --shift $shiftsize --extsize {params.smooth_win} \
            --nomodel -B --SPMR --keep-dup all --call-summits
        """
# NOTE: --call-summits is used to generate summit peaks, which can produce multiple entries for the same peak region. This will then be consolidated in the consensus peak calling step. (kinda useless, only affect is that it caps the 300k peaks based on #summits, but this is what the encode pipeline does)

rule macs2_cap_peaks:
    input:
        raw = "macs2_temp/{sample}_peaks.narrowPeak"
    output:
        narrow = "peaks/{sample}.narrowPeak.gz"
    params:
        npeaks = 300000
    threads: 1
    shell:
        """
        (sort -k8gr {input.raw} | head -n {params.npeaks} \
            | awk 'BEGIN{{OFS="\\t"}}{{$4="Peak_"NR;print}}' \
            | gzip -nc > {output.narrow}) || true
        """

rule macs2_make_fc_bw:
    input:
        treat = "macs2_temp/{sample}_treat_pileup.bdg",
        control = "macs2_temp/{sample}_control_lambda.bdg"
    output:
        fc_bw = "tracks/{sample}.fc.signal.bw"
    threads: THREADS
    resources:
        mem_mb = MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        macs2 bdgcmp -t {input.treat} -c {input.control} \
            --o-prefix "macs2_temp/{wildcards.sample}" -m FE
        slopBed -i macs2_temp/{wildcards.sample}_FE.bdg -g {CHRSIZES} -b 0 \
            | bedClip stdin {CHRSIZES} macs2_temp/{wildcards.sample}.fc.signal.bedgraph
        sort -k1,1 -k2,2n macs2_temp/{wildcards.sample}.fc.signal.bedgraph \
            | bedtools merge -i - -d -1 -c 4 -o mean \
            > macs2_temp/{wildcards.sample}.fc.signal.srt.bedgraph
        bedGraphToBigWig macs2_temp/{wildcards.sample}.fc.signal.srt.bedgraph \
            {CHRSIZES} {output.fc_bw}
        # rm -f macs2_temp/{wildcards.sample}.fc.signal.bedgraph \
        #        macs2_temp/{wildcards.sample}.fc.signal.srt.bedgraph \
        #        macs2_temp/{wildcards.sample}_FE.bdg
        """

rule macs2_make_pval_bw:
    input:
        tag = "tagAlign/{sample}.tn5.tagAlign.gz",
        treat = "macs2_temp/{sample}_treat_pileup.bdg",
        control = "macs2_temp/{sample}_control_lambda.bdg"
    output:
        pval_bw = "tracks/{sample}.pval.signal.bw"
    threads: THREADS
    resources:
        mem_mb = MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        sval=$(wc -l <(zcat -f {input.tag}) | awk '{{printf "%f", $1/1000000}}')
        macs2 bdgcmp -t {input.treat} -c {input.control} \
            --o-prefix "macs2_temp/{wildcards.sample}" -m ppois -S $sval
        slopBed -i macs2_temp/{wildcards.sample}_ppois.bdg -g {CHRSIZES} -b 0 \
            | bedClip stdin {CHRSIZES} macs2_temp/{wildcards.sample}.pval.signal.bedgraph
        sort -k1,1 -k2,2n macs2_temp/{wildcards.sample}.pval.signal.bedgraph \
            | bedtools merge -i - -d -1 -c 4 -o mean \
            > macs2_temp/{wildcards.sample}.pval.signal.srt.bedgraph
        bedGraphToBigWig macs2_temp/{wildcards.sample}.pval.signal.srt.bedgraph \
            {CHRSIZES} {output.pval_bw}
        # rm -f macs2_temp/{wildcards.sample}.pval.signal.bedgraph \
        #        macs2_temp/{wildcards.sample}.pval.signal.srt.bedgraph \
        #        macs2_temp/{wildcards.sample}_ppois.bdg \
        #        {input.treat} {input.control}
        """


############################################################
#  3b. Blacklist filtering                           
############################################################

# rule filter_peaks:
#     input:
#         narrow = rules.macs2_narrow.output.narrow
#     output:
#         filt = "peaks/{sample}.filt.narrowPeak.gz"
#     shell:
#         """
#         bedtools intersect -v -a {input.narrow} -b {BLACKLIST} \
#           | awk 'BEGIN{{OFS="\\t"}}{{$5=min($5,1000);print}}' \
#           | grep -P 'chr[0-9XY]+' \
#           | gzip -nc > {output.filt}
#         """


############################################################
#  3c. BED → bigBed conversion                        
############################################################

# rule make_bigbed:
#     input:
#         filt = rules.filter_peaks.output.filt
#     output:
#         bb = "peaks/{sample}.narrowPeak.bb"
#     shell:
#         """
#         zcat {input.filt} | sort -k1,1 -k2,2n > tmp.bed
#         bedToBigBed tmp.bed {CHRSIZES} {output.bb}
#         """


############################################################
#  3d. Naïve overlap threshold                        
############################################################

# rule naive_overlap:
#     input:
#         pooled = "peaks/pooled.narrowPeak.gz",
#         rep1   = rules.filter_peaks.output.filt,
#         rep2   = rules.filter_peaks.output.filt.replace("{sample}","rep2")
#     output:
#         overlap = "peaks/pooledInReps.narrowPeak.gz"
#     shell:
#         """
#         intersectBed -wo -a {input.pooled} -b {input.rep1} | \
#           awk '…fraction ≥ 0.5…' | cut -f1-10 \
#           | intersectBed -wo -a - -b {input.rep2} \
#           | awk '…fraction ≥ 0.5…' | cut -f1-10 \
#           | sort | uniq \
#           | bedtools intersect -v -b {BLACKLIST} \
#           | gzip -nc > {output.overlap}
#         """






############################################################
#  4. IDR (true & pseudo)
############################################################

# IDR helper for true replicate pairs
import itertools
import math

def get_true_pairs(condition):
    samples = SAMPLES[condition]
    if len(samples) < 2:
        raise ValueError(f"Condition {condition} has fewer than two replicates for IDR.")
    return list(itertools.combinations(samples, 2))

rule idr_true_pair:
    """
    Run IDR on each pair of true replicates for a condition.
    """
    input:
        peak1=lambda wc: f"peaks/{wc.rep1}.narrowPeak.gz",
        peak2=lambda wc: f"peaks/{wc.rep2}.narrowPeak.gz",
        pooled="peaks/{condition}.pooled.narrowPeak.gz"
    output:
        idr_txt="idr/{condition}.{rep1}-{rep2}.rep.idr.txt"
    params:
        idr_thresh=IDR_THRESH
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "idr_env"
    shell:
        """
        mkdir -p idr
        idr --samples {input.peak1} {input.peak2} \
            --peak-list {input.pooled} \
            --input-file-type narrowPeak \
            --rank p.value \
            --soft-idr-threshold {params.idr_thresh} \
            --output-file {output.idr_txt}
        """

rule idr_pooled_pseudoreps:
    """
    Run IDR on the pooled pseudoreplicates for each condition.
    """
    input:
        ppr1="peaks/{condition}.pooled.pr1.narrowPeak.gz",
        ppr2="peaks/{condition}.pooled.pr2.narrowPeak.gz",
        pooled="peaks/{condition}.pooled.narrowPeak.gz"
    output:
        idr_txt="idr/{condition}.ppr.idr.txt"
    params:
        idr_thresh=IDR_THRESH
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "idr_env"
    shell:
        """
        mkdir -p idr
        idr --samples {input.ppr1} {input.ppr2} \
            --peak-list {input.pooled} \
            --input-file-type narrowPeak \
            --rank p.value \
            --soft-idr-threshold {params.idr_thresh} \
            --output-file {output.idr_txt}
        """

rule idr_self_pseudoreps:
    """
    Run IDR on the pooled pseudoreplicates for each condition.
    """
    input:
        pr1="peaks/{sample}.pr1.narrowPeak.gz",
        pr2="peaks/{sample}.pr2.narrowPeak.gz",
        peak="peaks/{sample}.narrowPeak.gz"
    output:
        idr_txt="idr/{sample}.selfpr.idr.txt"
    params:
        idr_thresh=IDR_THRESH
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "idr_env"
    shell:
        """
        mkdir -p idr
        idr --samples {input.pr1} {input.pr2} \
            --peak-list {input.peak} \
            --input-file-type narrowPeak \
            --rank p.value \
            --soft-idr-threshold {params.idr_thresh} \
            --output-file {output.idr_txt}
        """

rule idr_finalize:
    """
    Generate conservative and optimal peak sets based on IDR results
    from true replicates and pseudoreplicates.
    """
    input:
        pooled="peaks/{condition}.pooled.narrowPeak.gz",
        rep_idrs=lambda wc: [
            f"idr/{wc.condition}.{a}-{b}.rep.idr.txt"
            for a, b in get_true_pairs(wc.condition)
        ],
        ppr_idr="idr/{condition}.ppr.idr.txt"
    output:
        cons="idr/{condition}.pooled.conservative.narrowPeak.gz",
        ppr_set="idr/{condition}.pooled.pprset.narrowPeak.gz",
        opt="idr/{condition}.pooled.optimal.narrowPeak.gz"
    params:
        idr_thresh=IDR_THRESH
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        mkdir -p idr
        IDR_THRESH_TRANSFORMED=$(awk -v p={params.idr_thresh} 'BEGIN{{print -log(p)/log(10)}}')
        # Determine best true-replicate pair by IDR peak count
        best=$(for f in {input.rep_idrs}; do
            c=$(awk 'NR>1 && $12>=IDR_THRESH_TRANSFORMED {{c++}} END {{print c}}' $f)
            echo -e "$c\t$f"
        done | sort -k1,1nr | head -n1 | cut -f2)
        # Convert the selected IDR result to narrowPeak for the conservative set
        awk -v t=$IDR_THRESH_TRANSFORMED 'BEGIN{{OFS="\t"}} $12>=t {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' $best \
            | sort | uniq | gzip -nc > {output.cons}

        # Convert the pseudoreplicate IDR result to narrowPeak for the ppr set
        awk -v t=$IDR_THRESH_TRANSFORMED 'BEGIN{{OFS="\t"}} $12>=t {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {input.ppr_idr} \
            | sort | uniq | gzip -nc > {output.ppr_set}

        # The optimal set is the set with the longest IDR peak count between the conservative set and the pseudoreplicate set
        n_cons=$(zcat {output.cons} | wc -l)
        n_ppr=$(zcat {output.ppr_set} | wc -l)
        if [ "$n_cons" -ge "$n_ppr" ]; then
            cp {output.cons} {output.opt}
        else
            cp {output.ppr_set} {output.opt}
        fi
        """


rule idr_qc:
    """
    Compute IDR QC metrics:
      N1/N2 = number of peaks passing threshold in each replicate’s self-pseudo-replicate IDR
      Nt     = peak count from best true-replicate IDR (conservative set)
      Np     = peak count from pooled pseudo-replicate IDR
    Then:
      RescueRatio = max(Nt,Np)/min(Nt,Np)
      SelfRatio   = max(N1,N2)/min(N1,N2)
      FLAG = -1 if both ratios >2, 0 if either >2, else 1
    """
    input:
        cons    = "idr/{condition}.pooled.conservative.narrowPeak.gz",
        ppr_set = "idr/{condition}.pooled.pprset.narrowPeak.gz",
        self_pr = lambda wc: expand("idr/{sample}.selfpr.idr.txt", sample=SAMPLES[wc.condition])
    output:
        idr_qc = "idr/{condition}_qc.txt"
    params:
        idr_thresh = IDR_THRESH
    threads: THREADS
    resources:
        mem_mb = MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        mkdir -p idr
        # transform p‐value threshold to –log10 scale for column 12
        IDR_THRESH_TRANS=$(awk -v p={params.idr_thresh} 'BEGIN{{print -log(p)/log(10)}}')

        # counts for conservative (Nt) and pooled pseudo (Np) sets
        Nt=$(zcat {input.cons}    | wc -l)
        Np=$(zcat {input.ppr_set} | wc -l)

        # self‐consistency counts for all replicates
        self_prs=({input.self_pr})
        self_counts=()
        for f in "${{self_prs[@]}}"; do
            c=$(awk 'NR>1 && $12>=IDR_THRESH_TRANS{{c++}} END{{print c}}' "$f")
            self_counts+=("$c")
        done
        # join counts with comma
        SC=$(IFS=,; echo "${{self_counts[*]}}")

        # compute ratios
        RescueRatio=$(awk -v n1=$Nt -v n2=$Np 'BEGIN{{print (n2>n1?n2:n1)/(n2<n1?n2:n1)}}')
        SelfRatio=$(awk -v s="$SC" 'BEGIN{{split(s,a,","); min=a[1]; max=a[1]; for(i in a){{if(a[i]<min)min=a[i]; if(a[i]>max)max=a[i]}}; print max/min}}')

        # flag reproducibility: -1=FAIL, 0=Borderline, 1=PASS
        FLAG=$(awk -v rr=$RescueRatio -v sr=$SelfRatio 'BEGIN{{if(rr>2 && sr>2) print -1; else if(rr>2 || sr>2) print 0; else print 1}}')

        # write header and metrics
        echo -e "self_counts\tNt\tNp\tRescueRatio\tSelfConsistencyRatio\tFlag" > {output.idr_qc}
        echo -e "$SC\t$Nt\t$Np\t$RescueRatio\t$SelfRatio\t$FLAG" >> {output.idr_qc}
        """

############################################################
#  Differential Accessibility Analysis
############################################################

# 1. Build consensus peaks across all samples
rule consensus_peaks:
    input:
        # All narrowPeak.gz for all samples in both groups
        # peaks = expand("peaks/{sample}.narrowPeak.gz", sample=SAMPLES_flat)
        peaks = expand("idr/{condition}.pooled.optimal.narrowPeak.gz", condition=SAMPLES)
    output:
        consensus = "analysis/consensus_peaks.bed"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        mkdir -p diff
        pigz -p $(( {threads} / 3 )) -dc {input.peaks} | cut -f1-3 | sort --parallel=$(( {threads} / 3 )) -k1,1 -k2,2n | bedtools merge -i - > {output.consensus}
        """

# 2. Generate count matrix for consensus peaks
rule count_matrix:
    input:
        consensus = rules.consensus_peaks.output.consensus,
        bams = expand("filtered/{sample}.PE.final.bam", sample=SAMPLES_flat)
    output:
        "analysis/counts.tsv"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        # Convert BED to SAF format for featureCounts
        awk 'BEGIN{{OFS="\\t"; print "GeneID","Chr","Start","End","Strand"}} {{print "peak"NR,$1,$2+1,$3,"."}}' {input.consensus} > analysis/consensus_peaks.saf
        featureCounts -p -a analysis/consensus_peaks.saf -o {output} -F SAF -T {threads} {input.bams}
        """

# 3. Perform differential accessibility analysis
rule diff_accessibility:
    input:
        counts = rules.count_matrix.output,
        metadata = "analysis/metadata.csv",
        script = "scripts/differential_accessibility.r"
    output:
        "analysis/deseq2/deseq2_results.csv"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        mkdir -p analysis/deseq2
        Rscript scripts/differential_accessibility.r {input.counts} {input.metadata} {output} > analysis/deseq2/deseq2.log 2>&1
        """




############################################################
#  Downstream TE/TE family analysis | per sample
############################################################

# Ensure directory function is available
from snakemake.io import directory

rule merge_te_annotation:
    input:
        counts=rules.count_matrix.output,
        annot=config["te_annotation"],
        script = "scripts/merge_count_te_annotation.py"
    output:
        merged="analysis/merged_results_count.csv",
        merged_multi="analysis/merged_results_count_multi.csv"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        mkdir -p analysis
        python3 {input.script} \
          --counts {input.counts} \
          --annot {input.annot} \
          --out {output.merged}
        """




############################################################
#  Downstream TE/TE family analysis for deseq2  | grouped
############################################################

# Ensure directory function is available
from snakemake.io import directory

rule merge_te_annotation_deseq2:
    input:
        deseq2= rules.diff_accessibility.output,
        counts= rules.count_matrix.output,
        annot=config["te_annotation"],
        script = "scripts/merge_deseq2_te_annotation.py"
    output:
        merged="analysis/merged_results_deseq2.csv",
        merged_multi="analysis/merged_results_deseq2_multi.csv"
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        mkdir -p analysis
        python3 {input.script} \
          --deseq2 {input.deseq2} \
          --counts {input.counts} \
          --annot {input.annot} \
          --out {output.merged}
        """

rule analyze_peak_enrichment_deseq2:
    input:
        merged=rules.merge_te_annotation_deseq2.output.merged,
    output:
        analysis_dir=directory("analysis/plots")
    threads: THREADS
    resources:
        mem_mb=MEM_MB
    conda:
        "ATAC_Core"
    shell:
        """
        mkdir -p {output.analysis_dir}
        python3 {config[my_script_dir]}/analyze_peak_enrichment.py \
          --input {input.merged} \
          --outdir {output.analysis_dir}
        """




############################################################
#  4g. FRiP                                             
############################################################

# rule frip:
#     input:
#         tag = rules.tn5_shift.output.shifted,
#         idr = rules.idr.output.cons
#     output:
#         frip = "qc/{sample}.frip.txt"
#     shell:
#         """
#         val1=$(bedtools intersect -a <(zcat -f {input.tag}) \
#                -b <(zcat -f {input.idr}) -wa -u | wc -l)
#         val2=$(zcat {input.tag} | wc -l)
#         echo $(awk 'BEGIN{{print {val1}/{val2}}}') > {output.frip}
#         """


############################################################
#  5. Coverage bigWig                              
############################################################

# rule count_signal:
#     input:
#         tag = rules.tn5_shift.output.shifted
#     output:
#         pos = "tracks/{sample}.positive.bigwig",
#         neg = "tracks/{sample}.negative.bigwig"
#     shell:
#         """
#         zcat {input.tag} | sort -k1,1 -k2,2n \
#           | bedtools genomecov -5 -bg -strand + -g {CHRSIZES} -i - \
#           > tmp.pos.bg
#         bedGraphToBigWig tmp.pos.bg {CHRSIZES} {output.pos}

#         zcat {input.tag} | sort -k1,1 -k2,2n \
#           | bedtools genomecov -5 -bg -strand - -g {CHRSIZES} -i - \
#           > tmp.neg.bg
#         bedGraphToBigWig tmp.neg.bg {CHRSIZES} {output.neg}
#         """


############################################################
#  7. TSS enrichment (annotation)                   
############################################################

# rule tss_enrich:
#     input:
#         bam = rules.dedup_filter.output.final_bam,
#         tss = TSS_BED
#     output:
#         plot = "tss/{sample}.tss_enrich.png",
#         log  = "tss/{sample}.tss_enrich.qc"
#     shell:
#         """
#         python make_tss_plot.py {input.bam} {input.tss} \
#           {wildcards.sample} {CHRSIZES} {config[read_len]}
#         """


############################################################
#  7. Fraction in annotated regions                 
############################################################

# rule annot_enrich:
#     input:
#         tag = rules.tn5_shift.output.shifted
#     output:
#         annot = "qc/{sample}.annotation.txt"
#     shell:
#         """
#         for bed in {DNASE_BED} {PROM_BED} {ENH_BED} {BLACKLIST}; do
#           printf "$(basename $bed)\\t"
#           bedtools sort -i $bed | bedtools merge -i - \
#             | bedtools intersect -u -a <(zcat {input.tag}) -b - \
#             | wc -l
#         done > {output.annot}
#         """


############################################################
#  USAGE on Slurm HPC
#
# snakemake --jobs 100 \
#           --cluster "sbatch -c {threads} --mem={resources.mem_mb} \
#                    --time={resources.time or '02:00:00'}" \
#           --latency-wait 60
############################################################