#!/usr/bin/env python3
"""Validate only the requested number of refine-loop iterations."""

from __future__ import annotations

import argparse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, required=True)
    args = parser.parse_args()
    if args.iterations < 1:
        parser.error("--iterations must be a positive integer")
    print(args.iterations)


if __name__ == "__main__":
    main()
