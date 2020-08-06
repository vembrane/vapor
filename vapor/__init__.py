import os
import argparse
from multiprocessing import Pool, Queue
from pysam import VariantFile
import tempfile
import sys
import subprocess
import time
import re

from functools import partial
from itertools import chain, zip_longest , islice

import multiprocessing as mp
from threading import Thread

# ctx = mp.get_context()

# class VaporProcess(ctx.Process):
#     def __init__(self, command, header, queue, *args, **kwargs):
#         print("start VaporProcess", file=sys.stderr)
#         cmd = f"{command}"
#         self.sp = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE)
#         self.queue = queue
#         super().__init__(*args, **kwargs)

# ctx.Process = VaporProcess

def execute_single():
    command = execute_single.command
    inqueue = execute_single.inqueue
    outqueue = execute_single.outqueue

    sp = subprocess.Popen(f"{command}", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)    
    inpipe = sp.stdin
    outpipe = sp.stdout

    inpipe.write(execute_single.header.encode())

    def write():
        while True:
            records = inqueue.get()
            if records is None: # stop signal for write
                print("send close signal", file=sys.stderr)
                inpipe.close()
                break
            for record in records:
                inpipe.write(record.encode())

    def read():
        for line in outpipe:
            if line.startswith(b"#"):
                continue
            outqueue.put(line.decode())

    writer = Thread(target=write)
    writer.start()
    reader = Thread(target=read)
    reader.start()
    writer.join()
    print("writer closed", file=sys.stderr)
    reader.join()
    print("reader closed", file=sys.stderr)
    print("worker closed", file=sys.stderr)


def grouper(n, records, bnd_reg, mate_id_reg):
    bnds = dict()
    ret = []
    for i, record in enumerate(records):
        if bnd_reg.match(record) is not None:
            record_id = record.split("\t", 3)[2]
            mate_id_match = mate_id_reg.match(record)
            if mate_id_match:
                # if no mate_id is available, this part is skipped
                # and the record is just piped out normally
                mate_id = mate_id_match.group(1)
                if mate_id in bnds:
                    ret.append(bnds[mate_id])
                    del bnds[mate_id]
                else:
                    bnds[record_id] = record
                    continue
        ret.append(record)
        if (i-1) % n == 0:
            yield ret
            ret = []
    if ret:
        yield ret


def init_worker(function, command, header, inqueue, outqueue):
    function.header = header
    function.inqueue = inqueue
    function.outqueue = outqueue
    function.command = command
    function.header = header


def execute(filename, command, stepsize, ncores, tmp_dir):
    # generate header
    with open(filename) as f:
        header = []
        for line in f:
            if not line.startswith("#"):
                break
            header.append(line)
        header = "".join(header)

        bnd_reg = re.compile(".*SVTYPE=BND.*")
        mate_id_reg = re.compile(".*MATEID=(\w*).*")
        inqueue = mp.Queue()
        outqueue = mp.Queue()

        init_worker(execute_single, command, header, inqueue, outqueue)
        worker = []

        # add worker Process
        for i in range(ncores):
            print("add worker", file=sys.stderr)
            p = mp.Process(target=execute_single)
            p.start()
            worker.append(p)

        def read():
            print(header, end="")
            active_count = ncores
            while active_count > 0:
                data = outqueue.get()
                if not data:
                    active_count -= 1
                else:
                    print(data, end="")

        reader = Thread(target=read)
        reader.start()

        record_chunks = grouper(stepsize, chain([line], f), bnd_reg, mate_id_reg)
        for chunk in record_chunks:
            inqueue.put(chunk)

        for i in range(ncores):
            inqueue.put(None)

        reader.join()
        # with Pool(1, initializer=init_worker, initargs=(execute_single, command, header, inqueue, outqueue)) as p:
        #     list(p.imap(execute_single, records))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command")
    parser.add_argument("vcf")
    parser.add_argument("--stepsize", "-s", type=int, default=10000)
    parser.add_argument("--threads", "-t", type=int, default=4)
    parser.add_argument(
        "--tmp", default=None, help="Define the base directory for the tmp files"
    )

    args = parser.parse_args()
    execute(args.vcf, args.command, args.stepsize, args.threads, args.tmp)
