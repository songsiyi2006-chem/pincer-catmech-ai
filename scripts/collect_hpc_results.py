"""Collect an already returned directory; never executes jobs or extracts archives."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from pincer_catmech.hpc.bridge import main


if __name__ == '__main__':
    raise SystemExit(main(['collect', *sys.argv[1:]]))
