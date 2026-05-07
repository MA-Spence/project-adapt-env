#!/usr/bin/env python3
"""Demo training script that writes deterministic dummy metrics."""

import argparse
import json
from pathlib import Path

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics = {
        'accuracy': 0.61,
        'auroc': 0.67,
        'note': 'demo metrics only'
    }
    (output_dir / 'metrics.json').write_text(json.dumps(metrics, indent=2) + '\n', encoding='utf-8')
    print(f'Wrote metrics to {output_dir / "metrics.json"}')

if __name__ == '__main__':
    main()
