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
    # generate outfilename
    outfilename = os.path.join(tmp_dir, f"{'_'.join(map(str, region))}.vcf")
    cmd = f"{command} > {outfilename}"
    sp = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE)
    pipe = sp.stdin

    with VariantFile(filename) as f:
        pipe.write(str(f.header).encode())
        for record in f.fetch(*region):
            pipe.write(str(record).encode())
    pipe.close()
    sp.wait()
    return outfilename


def create_regions(filename, stepsize):
    # smart region geration, large regions with no coverage are ignored
    with VariantFile(filename) as f:
        for name, contig in f.header.contigs.items():
            start = 0
            while start < contig.length:
                first_record = next(f.fetch(name, start), None)
                if not first_record:
                    break  # no more variants in contig
                start = first_record.pos
                stop = min(start + stepsize - 1, contig.length)
                yield (name, start, stop)
                start = stop + 1


def execute(filename, command, stepsize, ncores, tmp_dir):
    # create regions
    # regions = create_regions(filename, stepsize)
    regions = list(create_regions_advanced(filename, stepsize))
    print(regions)
    exit()

    # command = "SnpSift annotate {ss_args} {args.database_vcf} /dev/stdin"

    with tempfile.TemporaryDirectory(dir=tmp_dir) as tmp_dir:
        with Pool(ncores) as p:
            first = True
            execute_partial = partial(execute_single, filename, command, tmp_dir)
            for outfilename in p.imap(execute_partial, regions):
                # print("merge", outfilename, file=sys.stderr)
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
