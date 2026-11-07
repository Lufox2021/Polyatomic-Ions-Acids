"""Builds index.html from src/template.html and the src/*.js files."""
from pathlib import Path
root = Path(__file__).parent
t = (root / 'src/template.html').read_text(encoding='utf-8')
for tag, f in [('/*ENGINE*/', 'engine.js'), ('/*DATA*/', 'data.js'), ('/*UI*/', 'ui.js')]:
    t = t.replace(tag, (root / 'src' / f).read_text(encoding='utf-8'))
(root / 'index.html').write_text(t, encoding='utf-8')
print('index.html', len(t), 'bytes')
