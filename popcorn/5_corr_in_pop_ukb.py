# -*- coding: utf-8 -*-
# caculate r2 for each trait in CAS cohort
import subprocess
import os
from multiprocessing import Pool, cpu_count

eur_munged_dir = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/merged_gwas/White_British/munged/"
# eas_munged_dir = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/CAS/gwas/munged/"
eur_ld_ref_path = "/data1/jiapl_group/lishuhua/software/general/ldsc/LD_SCORE/EUR_baselineLD/baselineLD."
eur_w_ld_path = "/data1/jiapl_group/lishuhua/software/general/ldsc/LD_SCORE/EUR_ldscores/LDscore."
# eas_ld_ref_path = "/data1/jiapl_group/lishuhua/software/general/ldsc/LD_SCORE/EAS_baselineLD/baselineLD."
# eas_w_ld_path = "/data1/jiapl_group/lishuhua/software/general/ldsc/LD_SCORE/EAS_ldscores/weights.EAS.hm3_noMHC."
eur_output_path = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/UKB/merged_gwas/White_British/corr/"
# eas_output_path = "/data1/jiapl_group/lishuhua/project/PRS_benchmark/real_data/CAS/gwas/corr/"

def run_ldsc(sumstats_file, ld_ref, weights_ref, out_prefix):
    cmd = [
        "python", "/data1/jiapl_group/lishuhua/software/general/ldsc/ldsc.py",
        "--rg", sumstats_file,
        "--ref-ld-chr", ld_ref,
        "--w-ld-chr", weights_ref,
        "--out", out_prefix,
    ]
    subprocess.call(cmd)

trait_dict = {
        'p48': 'waist',
        'p50': 'height',
        'p102': 'pulse',
        'p4079': 'dbp',
        'p4080': 'sbp',
        'p20116': 'smoke',
        'p20117': 'drink',
        'p21001': 'bmi',
        'p30000': 'wbc',
        'p30010': 'rbc',
        'p30020':'hb',
        'p30080': 'plt',
        'p30120': 'lymph',
        'p30130': 'mono',
        'p30140': 'neut',
        'p30150': 'eos',
        'p30620': 'alt',
        'p30650': 'ast',
        'p30670': 'bun',
        'p30690': 'cholesterol',
        'p30700': 'creatinine',
        'p30730': 'ggt',
        'p30740': 'glucose',
        'p30760': 'hdl',
        'p30780': 'ldl',
        'p30870': 'triglycerides',
        'p30880': 'ua'
    }
for trait, trait_name in trait_dict.iteritems():
    for other_trait, other_trait_name in trait_dict.iteritems():
        # print "Comparing {} with {}".format(trait_name, other_trait_name)
        # p30140_int.merged.sumstats.gz
        if trait == other_trait:
            continue
        # if {trait_name}_{other_trait_name}_corr.log or {other_trait_name}_{trait_name}_corr.log exists, skip
        if os.path.exists(os.path.join(eur_output_path, "{}_{}_corr.log".format(trait_name, other_trait_name))) or os.path.exists(os.path.join(eur_output_path, "{}_{}_corr.log".format(other_trait_name, trait_name))):
            print("Correlation file already exists for traits: {} and {}, skipping...".format(trait_name, other_trait_name))
            continue
        # alt_int.alt_int.sumstats.gz
        eur_sums_1 = os.path.join(eur_munged_dir, "{}_int.merged.sumstats.gz".format(trait))
        eur_sums_2 = os.path.join(eur_munged_dir, "{}_int.merged.sumstats.gz".format(other_trait))
        eur_sum_path = '{},{}'.format(eur_sums_1, eur_sums_2)
        eur_out_prefix = os.path.join(eur_output_path, "{}_{}_corr".format(trait_name, other_trait_name))
        run_ldsc(eur_sum_path, eur_ld_ref_path, eur_w_ld_path, eur_out_prefix)