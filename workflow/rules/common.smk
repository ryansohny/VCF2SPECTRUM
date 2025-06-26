def get_mem_mb(wildcards, attempt):
    if attempt < 3:
        return attempt * 1024 * 8
    return attempt * 1024 * 16


def get_sample_vcf(wildcards):
    return config["samples"][wildcards.sample]["vcf"]


def get_all_samples():
    return list(config["samples"].keys())
