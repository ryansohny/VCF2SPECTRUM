def get_base_sample_name_rule(sample_key):
    """Extract base sample name from sample key (removes signature suffixes)"""
    for suffix in ['_SBS96', '_ID83', '_DBS78']:
        if sample_key.endswith(suffix):
            return sample_key[:-len(suffix)]
    return sample_key

rule vcf_to_sbs96:
    input:
        vcf=lambda wildcards: config["samples"][wildcards.sample_key]["vcf"],
        ref=config["reference_genome"]
    output:
        sbs96="{output_dir}/01.SBS/{sample_key}.SBS96"
    wildcard_constraints:
        sample_key="[^/]+"
    resources:
        mem_mb=get_mem_mb
    run:
        sample_sigs = config["samples"][wildcards.sample_key].get("signatures", [])
        if "SBS96" in sample_sigs:
            shell("""
            mutyper variants \
              --k 3 \
              {input.ref} \
              {input.vcf} \
              | mutyper spectra - --population \
              | sed -n 1,2p - \
              > {output.sbs96}
            """)
        else:
            # Create empty file if SBS96 not requested
            shell("touch {output.sbs96}")

def find_sample_key_for_sbs96(wildcards):
    """Find the sample key that corresponds to this base sample and has SBS96"""
    for sample_key, sample_config in config["samples"].items():
        if (get_base_sample_name_rule(sample_key) == wildcards.sample_base and 
            "SBS96" in sample_config.get("signatures", [])):
            return f"{wildcards.output_dir}/01.SBS/{sample_key}.SBS96"
    return f"{wildcards.output_dir}/01.SBS/dummy.SBS96"

rule process_sbs96:
    input:
        sbs96=find_sample_key_for_sbs96
    output:
        matrix="{output_dir}/{sample_base}_SBS96/{sample_base}.SBS96.all",
        plot="{output_dir}/{sample_base}_SBS96/SBS_96_plots_{sample_base}.pdf",
        plot_pct="{output_dir}/{sample_base}_SBS96/SBS_96_plots_{sample_base}.percentage.pdf"
    params:
        sample="{sample_base}",
        output_dir="{output_dir}/{sample_base}_SBS96"
    resources:
        mem_mb=get_mem_mb
    script:
        "../scripts/process_sbs96.py"


rule vcf_to_id83:
    input:
        vcf=lambda wildcards: config["samples"][wildcards.sample_key]["vcf"],
        ref=config["reference_genome"]
    output:
        id83="{output_dir}/03.ID/{sample_key}.ID83.tsv"
    wildcard_constraints:
        sample_key="[^/]+"
    resources:
        mem_mb=get_mem_mb
    run:
        sample_sigs = config["samples"][wildcards.sample_key].get("signatures", [])
        if "ID83" in sample_sigs:
            shell("""
            python -c "
import pandas as pd
import gzip as gz
import io

def read_vcf(path):
    if path.endswith('.gz'): 
        with gz.open(path, 'rb') as f:
            lines = [l.decode('utf-8') for l in f if not l.startswith(b'##')]
            return pd.read_csv(
                io.StringIO(''.join(lines)),
                dtype={{'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                       'QUAL': str, 'FILTER': str, 'INFO': str}},
                       sep='\t'
                       ).rename(columns={{'#CHROM': 'CHROM'}})
    else:
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
            return pd.read_csv(
                io.StringIO(''.join(lines)),
                dtype={{'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                       'QUAL': str, 'FILTER': str, 'INFO': str}},
                       sep='\t'
                       ).rename(columns={{'#CHROM': 'CHROM'}})

# Read VCF and filter for indels only
vcf = read_vcf('{input.vcf}')
indels = vcf[(vcf['REF'].str.len() != vcf['ALT'].str.len())].copy()
indels.to_csv('{output.id83}', sep='\t', index=False)
"
            """)
        else:
            shell("touch {output.id83}")


def find_sample_key_for_id83(wildcards):
    """Find the sample key that corresponds to this base sample and has ID83"""
    for sample_key, sample_config in config["samples"].items():
        if (get_base_sample_name_rule(sample_key) == wildcards.sample_base and 
            "ID83" in sample_config.get("signatures", [])):
            return f"{wildcards.output_dir}/03.ID/{sample_key}.ID83.tsv"
    return f"{wildcards.output_dir}/03.ID/dummy.ID83.tsv"

rule process_id83:
    input:
        id83=find_sample_key_for_id83,
        ref=config["reference_genome"]
    output:
        matrix="{output_dir}/{sample_base}_ID83/{sample_base}.ID83.all",
        plot="{output_dir}/{sample_base}_ID83/ID_83_plots_{sample_base}.pdf",
        plot_pct="{output_dir}/{sample_base}_ID83/ID_83_plots_{sample_base}.percentage.pdf"
    params:
        sample="{sample_base}",
        output_dir="{output_dir}/{sample_base}_ID83"
    resources:
        mem_mb=get_mem_mb
    script:
        "../scripts/process_id83.py"


def find_sample_key_for_dbs78(wildcards):
    """Find the sample key that corresponds to this base sample and has DBS78"""
    for sample_key, sample_config in config["samples"].items():
        if (get_base_sample_name_rule(sample_key) == wildcards.sample_base and 
            "DBS78" in sample_config.get("signatures", [])):
            return config["samples"][sample_key]["vcf"]
    return "dummy.vcf"

rule process_dbs78:
    input:
        vcf=find_sample_key_for_dbs78
    output:
        matrix="{output_dir}/{sample_base}_DBS78/{sample_base}.DBS78.all",
        plot="{output_dir}/{sample_base}_DBS78/DBS_78_plots_{sample_base}.pdf",
        plot_pct="{output_dir}/{sample_base}_DBS78/DBS_78_plots_{sample_base}.percentage.pdf"
    params:
        sample="{sample_base}",
        output_dir="{output_dir}/{sample_base}_DBS78"
    resources:
        mem_mb=get_mem_mb
    script:
        "../scripts/process_dbs78.py"

rule all_samples:
    input:
        expand("{output_dir}/{sample}/{sample}.SBS96.all",
               output_dir=config["output_dir"],
               sample=config["samples"].keys()),
        expand("{output_dir}/{sample}/SBS_96_plots_{sample}.pdf",
               output_dir=config["output_dir"],
               sample=config["samples"].keys()),
        expand("{output_dir}/{sample}/SBS_96_plots_{sample}.percentage.pdf",
               output_dir=config["output_dir"],
               sample=config["samples"].keys())
