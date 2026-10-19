"""Validated sample metadata reader for single and grouped Excel headers."""
from pathlib import Path
import re
import numpy as np
import pandas as pd

REQUIRED = {'raw_file', 'sample_type', 'region', 'group', 'injection_order'}
REGIONS = {'PL','IL','BLA','CEA','DHPC','VHPC','PVH','PMD','SSC','MSC','LPB'}


def read_sample_metadata(path: Path) -> pd.DataFrame:
    table = None
    for header in (0, 1):
        candidate = pd.read_csv(path, header=header) if Path(path).suffix.lower()=='.csv' else pd.read_excel(path, header=header)
        candidate.columns = [str(c).strip() for c in candidate.columns]
        if REQUIRED.issubset(candidate.columns):
            table = candidate
            break
    if table is None:
        raise ValueError(f'{path}: required metadata fields {sorted(REQUIRED)} were not found in header rows 1 or 2')
    raw = table['raw_file'].astype('string').str.strip()
    valid = raw.str.endswith('.raw', na=False)
    invalid = raw.notna() & raw.ne('') & ~valid
    if invalid.any():
        raise ValueError(f'{path}: invalid nonempty raw_file values: {raw[invalid].tolist()[:5]}')
    table = table.loc[valid].copy()
    if table.empty:
        raise ValueError(f'{path}: no raw sample records')
    table['raw_file'] = raw.loc[valid]
    table['sample_id'] = table['raw_file'].str.removesuffix('.raw')
    for c in ('sample_type', 'region', 'group'):
        table[c] = table[c].astype('string').str.strip()
    qc_blank_region = table['sample_type'].eq('QC') & (table['region'].isna() | table['region'].eq(''))
    table.loc[qc_blank_region, 'region'] = 'QC'
    order = pd.to_numeric(table['injection_order'], errors='coerce')
    if order.isna().any() or not np.isfinite(order).all() or (order <= 0).any() or (order != np.floor(order)).any():
        raise ValueError(f'{path}: injection_order must contain positive integers')
    table['injection_order'] = order.astype(int)
    if table['sample_id'].duplicated().any() or table['injection_order'].duplicated().any():
        raise ValueError(f'{path}: duplicate sample_id or injection_order')
    if not table['sample_type'].isin(['BIO', 'QC']).all():
        raise ValueError(f'{path}: sample_type must be BIO or QC')
    biological_fields = table.loc[table['sample_type'].eq('BIO'), ['region', 'group']]
    if (biological_fields.isna() | biological_fields.eq('')).any().any():
        raise ValueError(f'{path}: a biological sample lacks region or group')
    identities = []
    for row in table.itertuples():
        if row.sample_type == 'QC':
            if not re.fullmatch(r'Astral_MSP2600441_QC_\d+\.raw', row.raw_file) or row.group != 'QC' or row.region != 'QC':
                raise ValueError('QC filename/type/region/group mismatch')
            continue
        match = re.fullmatch(r'Astral_MSP2600441_\d+_(G[1-6])_([1-4])_([A-Z]+)(?:_\d{14})?\.raw', row.raw_file)
        if not match or row.region not in REGIONS or match.group(1) != row.group or match.group(3) != row.region:
            raise ValueError('Biological filename/type/region/group mismatch')
        identities.append((row.region, row.group, match.group(2)))
    if len(identities) != len(set(identities)):
        raise ValueError('Duplicate biological region/group/replicate')
    return table.reset_index(drop=True)


def match_metadata(table: pd.DataFrame, sample_ids) -> pd.DataFrame:
    expected = pd.Index(sample_ids)
    if expected.duplicated().any():
        raise ValueError('Duplicate sample IDs in the quantitative matrix')
    key = table.set_index('sample_id')
    missing = expected.difference(key.index)
    extra = key.index.difference(expected)
    if len(missing) or len(extra):
        raise ValueError(f'Metadata/matrix mismatch: missing={missing.tolist()[:5]}, extra={extra.tolist()[:5]}')
    return key.loc[expected]
