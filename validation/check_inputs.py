"""Verify supplied processed-report filenames, byte sizes and SHA-256 checksums."""
import argparse
import csv
import hashlib
import re
from pathlib import Path


def report_path(data_dir: Path, filename: str) -> Path | None:
    exact = data_dir / filename
    if exact.is_file():
        return exact
    # Some delivered filenames have a space immediately before ".csv".
    expected = Path(filename)
    matches = [path for path in data_dir.glob('*.csv')
               if path.stem.rstrip() == expected.stem and path.suffix.lower() == '.csv']
    return matches[0] if len(matches) == 1 else None


def check_g_headers(path: Path) -> bool:
    region_match = re.search(r'_(PMD|LPB|VHPC|PVH|BLA|DHPC|CEA|MSC|SSC|IL|PL)_DIA_', path.name)
    if region_match is None:
        return True  # The pooled-QC report has QC identifiers, not G1-G6.
    region = region_match.group(1)
    with path.open(encoding='utf-8-sig', newline='') as handle:
        columns = next(csv.reader(handle))
    groups = []
    for column in columns:
        if not column.endswith('.raw.PG.Quantity'):
            continue
        match = re.search(rf'_(G[1-6])_([1-4])_{region}(?:_\d{{14}})?\.raw\.PG\.Quantity$', column)
        if match is None:
            return False
        groups.append(match.groups())
    return len(groups) == 24 and set(groups) == {(f'G{g}', str(rep)) for g in range(1, 7) for rep in range(1, 5)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, required=True)
    args = parser.parse_args()
    manifest = Path(__file__).resolve().parents[1] / 'metadata/input_manifest.csv'
    failures = []
    with manifest.open(encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        path = report_path(args.data_dir, row['filename'])
        if path is None:
            failures.append(f"Missing report: {row['filename']}")
            continue
        with path.open('rb') as handle:
            actual = hashlib.file_digest(handle, 'sha256').hexdigest()
        if path.stat().st_size != int(row['bytes']) or actual != row['sha256']:
            failures.append(f"Changed report: {row['filename']}")
        elif not check_g_headers(path):
            failures.append(f"Expected 24 G1-G6 quantity columns: {path.name}")
    if failures:
        raise SystemExit('\n'.join(failures))
    print(f'{len(rows)} processed reports match the release input manifest.')

if __name__ == '__main__':
    main()
