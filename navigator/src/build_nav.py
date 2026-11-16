"""Build the Relativistic Navigator mockup.

Usage: python3 -I build_nav.py <navigator.src.html> <stars.json> <out.html> [--standalone] [--stamp VERSION]

Injects the star list into the template. With --standalone, wraps the page in its own HTML skeleton so it runs as an
ordinary web page (needed for phone gyro and WebXR, which the claude.ai artifact viewer blocks).

With --stamp VERSION, writes a build stamp that the page's Debug panel shows: the artifact version it is published as,
the git commit and branch of the source, whether the source file is exactly that commit's (or has changes after it),
the build time (UTC) and the first 12 hex digits of the source's SHA-256. Commit the source first, then build with the
stamp, so that the stamp names the commit the page was built from.
"""
import datetime
import hashlib
import json
import os
import subprocess
import sys

src, data, out = sys.argv[1:4]
opts = sys.argv[4:]
standalone = '--standalone' in opts
version = opts[opts.index('--stamp') + 1] if '--stamp' in opts else None
raw = open(src, 'rb').read()
s = raw.decode('utf-8')
backend_path = os.path.join(os.path.dirname(os.path.abspath(src)), 'srkit-backend.js')
with open(backend_path, encoding='utf-8') as backend_file:
    backend = backend_file.read()
assert s.count('/*SRKIT_BACKEND*/') == 1
s = s.replace('/*SRKIT_BACKEND*/', backend.replace('</', '<\\/'))
j = open(data, encoding='utf-8').read().replace('</', '<\\/')
assert s.count('/*STARS_JSON*/') == 1 and s.count('/*BUILD_INFO*/') == 1
info = None
if version is not None:
    d = os.path.dirname(os.path.abspath(src))
    git = lambda *a: subprocess.run(['git', '-C', d, *a], capture_output=True, text=True).stdout.strip()
    commit, branch, top = git('rev-parse', '--short=12', 'HEAD'), git('rev-parse', '--abbrev-ref', 'HEAD'), git('rev-parse', '--show-toplevel')
    rel = os.path.relpath(os.path.abspath(src), top) if top else src
    at_head = subprocess.run(['git', '-C', d, 'show', f'HEAD:{rel}'], capture_output=True).stdout == raw if commit else False
    remote = git('config', '--get', 'remote.origin.url')
    info = {'version': version, 'commit': commit or None, 'branch': branch or None, 'clean': at_head,
            'built': datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC'),
            'src': hashlib.sha256(raw + backend.encode('utf-8') + j.encode('utf-8')).hexdigest()[:12], 'path': rel,
            'repo': remote.split('@')[-1].replace(':', '/').removesuffix('.git').split('github.com/')[-1] if 'github.com' in remote else None}
page = s.replace('/*STARS_JSON*/', j).replace('/*BUILD_INFO*/null', json.dumps(info).replace('</', '<\\/'))
if standalone:
    page = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">'
            '<style>:root{padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}'
            'body{margin:0}img{max-width:100%}[hidden]{display:none!important}</style></head><body>'
            + page + '</body></html>\n')
open(out, 'w', encoding='utf-8').write(page)
print(f'{out}: {len(page) / 1e6:.2f} MB{" (standalone)" if standalone else ""}' + (f'; stamp {json.dumps(info)}' if info else ''))
