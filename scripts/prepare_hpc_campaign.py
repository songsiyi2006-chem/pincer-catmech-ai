"""Generate a portable ORCA/Slurm diagnostic bundle without submitting jobs."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pincer_catmech.hpc.bridge import prepare, resources


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--selection', type=Path, required=True)
    p.add_argument('--source-root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--nprocs', type=int, required=True)
    p.add_argument('--memory-mb', type=int, required=True)
    p.add_argument('--maxcore-mb', type=int, required=True)
    p.add_argument('--walltime', required=True)
    p.add_argument('--concurrency', type=int, required=True)
    a = p.parse_args()
    result = prepare(a.selection, a.source_root, a.output,
                     resources(a.nprocs, a.memory_mb, a.maxcore_mb, a.walltime, a.concurrency))
    print(json.dumps(dict(status=result['status'], tasks=len(result['tasks']),
                          output=str(a.output.resolve()), protocol_sha256=result['protocol_sha256']), indent=2))


if __name__ == '__main__':
    main()
