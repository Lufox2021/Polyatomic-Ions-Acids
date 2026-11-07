# Polyatomic Ions & Acids

A coaching tool for general chemistry students: naming polyatomic anions, binary acids and oxoacids.

- **Structure → Name:** see an all-octet Lewis structure (with formal charges), type the name.
- **Name → Structure:** given a name, type the formula or draw the Lewis structure.
- **Review:** naming rules and reference tables.
- **Instructor:** upload a .csv or tab-delimited .txt to add or change species (saved in that browser only).

Live tool: https://lufox2021.github.io/Polyatomic-Ions-Acids/ (after GitHub Pages is turned on: Settings → Pages → deploy from `main`, `/ (root)`).

## Files
- `index.html` – the tool (what students load). Rebuild from `src/` with `python3 build.py`.
- `src/` – `engine.js`, `data.js`, `ui.js`, `template.html`
- `test/test_all.js` – checks every built-in species (`node test/test_all.js`).
