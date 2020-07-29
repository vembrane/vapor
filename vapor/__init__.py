import os
import argparse
from multiprocessing import Pool
from pysam import VariantFile
import tempfile
import sys
import subprocess
import time

from functools import partial

def execute_single(filename, command, tmp_dir, region):
    #generate outfilename
    outfilename = os.path.join(tmp_dir, f"{'_'.join(map(str, region))}.vcf")
    #cmd = f"{command}| bcftools view -Ob > {tmp_dir}/{outfilename}.bcf"
    cmd = f"{command} > {outfilename}"
    #print(cmd)
    #| bcftools view -Ob > {tmp_dir}/{outfilename}.bcf
    sp = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE)
    pipe = sp.stdin

    with VariantFile(filename) as f:
        pipe.write(str(f.header).encode())
        for record in f.fetch(*region):
            pipe.write(str(record).encode())
    pipe.close()
    sp.wait()
    return outfilename


def execute(filename, command, stepsize, ncores, tmp_dir):
    # create regions
    regions = []
    with VariantFile(filename) as f:
        for name, contig in f.header.contigs.items():
            last = 0
            for x in list(range(0, contig.length, stepsize)):
                stop = min(last + stepsize - 1, contig.length)
                regions.append((name, last, stop))
                last = stop + 1

    # command = "SnpSift annotate {ss_args} {args.database_vcf} /dev/stdin"

    with tempfile.TemporaryDirectory(dir=tmp_dir) as tmp_dir:
        with Pool(ncores) as p:
            first = True
            execute_partial = partial(execute_single, filename, command, tmp_dir)
            for outfilename in p.imap(execute_partial, regions):
                #print("merge", outfilename, file=sys.stderr)
                if first:
                    cmd = f"cat {outfilename}"
                else:
                    cmd = f'cat {outfilename} | grep -v "#"'
                os.system(cmd)
                os.remove(outfilename)
                first = False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf")
    parser.add_argument("command")
    parser.add_argument("--stepsize", "-s", type=int, default=1000000)
    parser.add_argument("--threads", "-t", type=int, default=8)
    parser.add_argument(
        "--tmp", default=None, help="Define the base directory for the tmp files"
    )

    args = parser.parse_args()
    execute(args.vcf, args.command, args.stepsize, args.threads, args.tmp)
