#!/usr/bin/env python3
"""Read demo metrics and print a small summary."""

import argparse
import json
from pathlib import Path

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--metrics', default='metrics.json')
    args = parser.parse_args()
    metrics = json.loads(Path(args.metrics).read_text(encoding='utf-8'))
    print('Demo validation summary')
    for key, value in metrics.items():
        print(f'{key}: {value}')

if __name__ == '__main__':
    main()
