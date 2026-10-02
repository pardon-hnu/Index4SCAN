"""Collect per-edge update times, excluding loading/construction timers.

Usage:
    python collect_update_time.py update.res
    python collect_update_time.py update.res --csv update_times.csv
    python collect_update_time.py "res/update/*.res"
"""

import argparse
import csv
import glob
import math
import re
import statistics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", help="Log files or wildcard patterns")
    parser.add_argument("--csv", help="Optional per-update CSV output")
    args = parser.parse_args()

    time_pattern = re.compile(
        r"\[UPDATE\]\s+update time\s*:\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    )
    edge_pattern = re.compile(
        r"now\s+(insert|remove|delete)\s+edge\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)"
    )
    files = []
    for pattern in args.files:
        matches = sorted(glob.glob(pattern))
        if not matches:
            parser.error("No matching file: " + pattern)
        for path in matches:
            if path not in files:
                files.append(path)

    if args.csv:
        import os
        if os.path.abspath(args.csv) in [os.path.abspath(path) for path in files]:
            parser.error("CSV output must not overwrite an input log")

    rows = []
    for path in files:
        times = []
        edge = ("", "", "")
        pending = False
        missing = 0
        with open(path, encoding="utf-8", errors="replace") as log:
            for line in log:
                match = edge_pattern.search(line)
                if match:
                    missing += int(pending)
                    edge = match.groups()
                    pending = True
                match = time_pattern.search(line)
                if match:
                    elapsed = float(match.group(1))
                    times.append(elapsed)
                    rows.append((path, len(times), *edge, match.group(1)))
                    edge = ("", "", "")
                    pending = False
        missing += int(pending)

        print("[FILE] " + path)
        print("[UPDATE COUNT] " + str(len(times)))
        if times:
            total = math.fsum(times)
            print("[TOTAL UPDATE TIME] {:.12g} s".format(total))
            print("[AVERAGE UPDATE TIME] {:.12g} s".format(total / len(times)))
            print("[MEDIAN UPDATE TIME] {:.12g} s".format(statistics.median(times)))
            print("[MIN UPDATE TIME] {:.12g} s".format(min(times)))
            print("[MAX UPDATE TIME] {:.12g} s".format(max(times)))
        if missing:
            print("[WARNING] {} edge update(s) have no recorded update time".format(missing))

    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8") as output:
            writer = csv.writer(output)
            writer.writerow(("file", "update_id", "operation", "u", "v", "time_seconds"))
            writer.writerows(rows)
        print("[CSV] " + args.csv)


if __name__ == "__main__":
    main()
