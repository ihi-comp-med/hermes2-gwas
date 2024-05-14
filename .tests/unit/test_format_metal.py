import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_format_metal():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/format_metal/data")
        expected_path = PurePosixPath(".tests/unit/format_metal/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("results/meta/Pheno1_EUR/FORMAT-GWAS.tsv.gz results/meta/Pheno1_EUR/FORMAT-GWAS.tsv.gz.tbi results/meta/Pheno1_EUR/FORMAT-GWAS.tsv.gz.md5", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "results/meta/Pheno1_EUR/FORMAT-GWAS.tsv.gz results/meta/Pheno1_EUR/FORMAT-GWAS.tsv.gz.tbi results/meta/Pheno1_EUR/FORMAT-GWAS.tsv.gz.md5",
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
