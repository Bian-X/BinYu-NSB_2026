"""Verify the source-only package without accessing study data or the network."""
from pathlib import Path
import hashlib,sys
root=Path(__file__).resolve().parents[1]
errors=[];count=0
for line in (root/'SOURCE_MANIFEST.sha256').read_text(encoding='utf-8').splitlines():
    if not line.strip():continue
    expected,name=line.split(None,1);path=(root/name.strip()).resolve()
    if not path.is_relative_to(root):raise ValueError('Manifest path escapes package root')
    count+=1
    if not path.is_file():errors.append('Missing: '+name);continue
    with path.open('rb') as f:actual=hashlib.file_digest(f,'sha256').hexdigest()
    if actual!=expected:errors.append('Changed: '+name)
if errors:print('\n'.join(errors));sys.exit(1)
print(f'{count} source files match SOURCE_MANIFEST.sha256')
