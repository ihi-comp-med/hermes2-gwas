import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_qc_step1():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/qc_step1/data")
        expected_path = PurePosixPath(".tests/unit/qc_step1/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("results/qc/step1/StudyA/GWAS-QC1_Pheno1_EUR.tsv.gz results/qc/step1/StudyA/GWAS-QC1_Pheno1_EUR.tsv.gz.tbi results/qc/step1/StudyA/EXCLUDED-VARS-QC1_Pheno1_EUR.tsv.gz results/qc/step1/StudyA/VAR-STATS-QC1_Pheno1_EUR.tsv", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "results/qc/step1/StudyA/GWAS-QC1_Pheno1_EUR.tsv.gz results/qc/step1/StudyA/GWAS-QC1_Pheno1_EUR.tsv.gz.tbi results/qc/step1/StudyA/EXCLUDED-VARS-QC1_Pheno1_EUR.tsv.gz results/qc/step1/StudyA/VAR-STATS-QC1_Pheno1_EUR.tsv",
            "-f", 
            "-j1",
            "--target-files-omit-workdir-adjustment",
    
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
