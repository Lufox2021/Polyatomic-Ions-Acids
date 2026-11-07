// ===== UI =====
const $ = id => document.getElementById(id);
const NS = 'http://www.w3.org/2000/svg';
const COLORS = {O:'#c0392b', N:'#2457c5', S:'#a87800', P:'#c35a00', Cl:'#1e7d46', Br:'#8a3b12', I:'#6b2fa0', F:'#1e7d46', H:'#55606f', C:'#1d2433', Cr:'#5a6b7f', Mn:'#7a4a9e', Se:'#a87800'};
const STORE_KEY = 'polyacid.uploads.v1';
function store(k, v){ try { if (v === undefined) return localStorage.getItem(k); if (v === null) localStorage.removeItem(k); else localStorage.setItem(k, v); } catch(e){ return null; } }

let LIST = [];
function loadSpecies(){
  const byId = new Map();
  RAW_SPECIES.forEach(r => byId.set(r.id, r));
  let up = [];
  try { up = JSON.parse(store(STORE_KEY) || '[]'); } catch(e){ up = []; }
  up.forEach(r => byId.set(r.id, r));
  const out = [];
  for (const r of byId.values()){ try { out.push(prepareSpecies(r)); } catch(e){ console.warn(r.id, e.message); } }
  LIST = linkSpecies(out);
  LIST.uploaded = up.length;
}

// ---------- SVG rendering ----------
function mk(tag, attrs, parent){ const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }
function angDist(a,b){ let d = Math.abs(a-b) % 360; return d > 180 ? 360 - d : d; }
function pickDirs(occupied, k, cands){
  const out = [], occ = occupied.slice();
  for (let n = 0; n < k; n++){
    let best = null, bs = -1;
    for (const c of cands){ if (out.includes(c)) continue; const s = occ.length ? Math.min(...occ.map(o => angDist(o,c))) : 999; if (s > bs + 0.01){ bs = s; best = c; } }
    out.push(best); occ.push(best);
  }
  return out;
}
// S: {atoms:[{el,x,y,lp,fc}], bonds:[{a,b,order}], charge}; coords in px (y down)
function renderStruct(svg, S, opt = {}){
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const g = mk('g', {}, svg);
  const bad = new Set(opt.bad || []);
  const R = 13;
  const bondAng = S.atoms.map(() => []);
  S.bonds.forEach(b => {
    const A = S.atoms[b.a], B = S.atoms[b.b];
    bondAng[b.a].push(Math.atan2(B.y - A.y, B.x - A.x)*180/Math.PI);
    bondAng[b.b].push(Math.atan2(A.y - B.y, A.x - B.x)*180/Math.PI);
  });
  for (const i of bad){ const a = S.atoms[i]; if (a) mk('circle', {cx:a.x, cy:a.y, r:22, fill:'#fdecea', stroke:'#b3261e', 'stroke-width':2}, g); }
  if (opt.sel != null && S.atoms[opt.sel]){ const a = S.atoms[opt.sel]; mk('circle', {cx:a.x, cy:a.y, r:20, fill:'#e8eefc', stroke:'#2457c5', 'stroke-width':1.5}, g); }
  S.bonds.forEach(b => {
    const A = S.atoms[b.a], B = S.atoms[b.b];
    const dx = B.x - A.x, dy = B.y - A.y, L = Math.hypot(dx,dy) || 1, ux = dx/L, uy = dy/L, px = -uy, py = ux;
    if (L < 2*R + 2) return;
    const x1 = A.x + ux*R, y1 = A.y + uy*R, x2 = B.x - ux*R, y2 = B.y - uy*R;
    const offs = b.order === 1 ? [0] : b.order === 2 ? [-3.6, 3.6] : [-5.5, 0, 5.5];
    for (const o of offs) mk('line', {x1:x1+px*o, y1:y1+py*o, x2:x2+px*o, y2:y2+py*o, stroke:'#1d2433', 'stroke-width':2.2, 'stroke-linecap':'round'}, g);
  });
  const lpC = [270, 0, 90, 180]; for (let d = 0; d < 360; d += 15) if (d % 90) lpC.push(d);
  S.atoms.forEach((a,i) => {
    const t = mk('text', {x:a.x, y:a.y, 'text-anchor':'middle', 'dominant-baseline':'central', 'font-size':21, 'font-weight':600, fill: COLORS[a.el] || '#1d2433', 'font-family':'system-ui,sans-serif'}, g);
    t.textContent = a.el;
    const occ = bondAng[i].map(x => (x + 360) % 360);
    const lpd = pickDirs(occ, a.lp || 0, lpC);
    const rr = a.el.length > 1 ? 19 : 16;
    lpd.forEach(d => {
      const t2 = d*Math.PI/180, cx = a.x + Math.cos(t2)*rr, cy = a.y + Math.sin(t2)*rr, qx = -Math.sin(t2)*4.2, qy = Math.cos(t2)*4.2;
      mk('circle', {cx:cx+qx, cy:cy+qy, r:2.3, fill:'#1d2433'}, g);
      mk('circle', {cx:cx-qx, cy:cy-qy, r:2.3, fill:'#1d2433'}, g);
    });
    if (a.fc){
      const d = pickDirs(occ.concat(lpd), 1, [315, 45, 225, 135, 270, 90, 0, 180])[0];
      const t3 = d*Math.PI/180, fx = a.x + Math.cos(t3)*26, fy = a.y + Math.sin(t3)*26;
      mk('circle', {cx:fx, cy:fy, r:9.5, fill:'#fff', stroke: a.fc > 0 ? '#2457c5' : '#b3261e', 'stroke-width':1.3}, g);
      const ft = mk('text', {x:fx, y:fy+0.5, 'text-anchor':'middle', 'dominant-baseline':'central', 'font-size':11.5, 'font-weight':700, fill: a.fc > 0 ? '#2457c5' : '#b3261e', 'font-family':'system-ui,sans-serif'}, g);
      ft.textContent = (Math.abs(a.fc) > 1 ? Math.abs(a.fc) : '') + (a.fc > 0 ? '+' : '−');
    }
  });
  if (!S.atoms.length && opt.empty){
    const t = mk('text', {x:350, y:195, 'text-anchor':'middle', 'font-size':16, fill:'#8a93a5'}, g); t.textContent = opt.empty;
  }
  let box = null;
  if (S.atoms.length){
    const xs = S.atoms.map(a => a.x), ys = S.atoms.map(a => a.y);
    box = {x0: Math.min(...xs) - 40, x1: Math.max(...xs) + 40, y0: Math.min(...ys) - 40, y1: Math.max(...ys) + 40};
    if (S.charge){
      const {x0,x1,y0,y1} = box, w = 8;
      mk('path', {d:`M${x0+w},${y0} H${x0} V${y1} H${x0+w}`, fill:'none', stroke:'#1d2433', 'stroke-width':2}, g);
      mk('path', {d:`M${x1-w},${y0} H${x1} V${y1} H${x1-w}`, fill:'none', stroke:'#1d2433', 'stroke-width':2}, g);
      const ct = mk('text', {x:x1+4, y:y0+4, 'font-size':18, 'font-weight':700, fill:'#1d2433', 'font-family':'system-ui,sans-serif'}, g);
      ct.textContent = chargeText(S.charge);
      box.x1 += 30;
    }
  }
  if (opt.temp) mk('line', Object.assign({stroke:'#2457c5', 'stroke-width':2, 'stroke-dasharray':'5 4'}, opt.temp), g);
  if (opt.fit && box){
    const pad = 10, w = box.x1 - box.x0 + 2*pad, h = box.y1 - box.y0 + 2*pad;
    const W = Math.max(w, opt.minW || 440), H = Math.max(h, opt.minH || 250);
    const cx = (box.x0 + box.x1)/2, cy = (box.y0 + box.y1)/2;
    svg.setAttribute('viewBox', `${cx - W/2} ${cy - H/2} ${W} ${H}`);
    svg.setAttribute('preserveAspectRatio', 'xMidYMid meet');
  }
}
function keyDraw(it, scale = 58){
  return {atoms: it.skel.atoms.map((a,i) => ({el:a.el, x:a.x*scale, y:-a.y*scale, lp:it.key.lp[i], fc:it.key.fc[i]})),
          bonds: it.skel.bonds.map(([a,b],i) => ({a, b, order: it.key.ord[i]})), charge: it.charge};
}

// ---------- shared panel state ----------
const CATS = [['all','All'],['anion','Polyatomic anions'],['binary','Binary acids'],['oxoacid','Oxoacids']];
const ST = {filter:{s2n: store('polyacid.f.s2n') || 'all', n2s: store('polyacid.f.n2s') || 'all'}, cur:{}, tries:{}, score:{s2n:{t:new Set(), s:new Set()}, n2s:{t:new Set(), s:new Set()}}, mode:'formula', promptName:''};
function pool(p){ const f = ST.filter[p]; return LIST.filter(x => f === 'all' || x.kind === f); }
function pick(p){ const P = pool(p); if (!P.length) return null; let c, n = 0; do { c = P[Math.floor(Math.random()*P.length)]; } while (P.length > 1 && c === ST.cur[p] && ++n < 20); return c; }
function drawFilters(p){
  const box = $('f-' + p); box.innerHTML = '';
  CATS.forEach(([k, lab]) => {
    const n = LIST.filter(x => k === 'all' || x.kind === k).length;
    const b = document.createElement('button'); b.innerHTML = `${lab}<span class="n">${n}</span>`;
    if (ST.filter[p] === k) b.className = 'on';
    b.onclick = () => { ST.filter[p] = k; store('polyacid.f.' + p, k); drawFilters(p); (p === 's2n' ? s2nLoad : n2sLoad)(pick(p)); };
    box.appendChild(b);
  });
}
function updScore(p){ const s = ST.score[p]; $(p + '-score').textContent = `This session: ${s.s.size} solved of ${s.t.size} tried`; }
function feedback(box, res, tries, rule, onShow, extraBtns){
  box.className = 'fb ' + (res.ok ? 'ok' : 'bad');
  let h = `<div>${res.msg}</div>`;
  if (!res.ok && tries >= 2 && rule) h += `<div class="rule"><b>Rule:</b> ${rule}</div>`;
  box.innerHTML = h;
  if (!res.ok && tries >= 3 && onShow){ const b = document.createElement('button'); b.className = 'plain'; b.textContent = 'Show me'; b.onclick = onShow; box.appendChild(b); }
  (extraBtns || []).forEach(([lab, fn, cls]) => { const b = document.createElement('button'); b.className = cls || 'plain'; b.textContent = lab; b.onclick = fn; box.appendChild(b); });
}
function clearFb(box){ box.className = 'fb'; box.innerHTML = ''; }

// ---------- Structure → Name ----------
function s2nLoad(it){
  if (!it) return;
  ST.cur.s2n = it; ST.tries.s2n = 0;
  renderStruct($('s2n-svg'), keyDraw(it), {fit:true});
  $('s2n-formula').textContent = $('s2n-showf').checked ? dispFormula(it) : '';
  $('s2n-in').value = ''; clearFb($('s2n-fb'));
}
function s2nCheck(){
  const it = ST.cur.s2n; if (!it) return;
  const r = diagnoseName(it, $('s2n-in').value, LIST);
  if (!$('s2n-in').value.trim()){ feedback($('s2n-fb'), r, 0); return; }
  ST.score.s2n.t.add(it.id);
  if (r.ok){ ST.score.s2n.s.add(it.id); feedback($('s2n-fb'), {ok:true, msg: r.msg + '<div class="small" style="margin-top:6px">' + explain(it, LIST) + '</div>'}, 0, null, null, [['Next structure →', () => s2nLoad(pick('s2n')), 'primary']]); }
  else { ST.tries.s2n++; feedback($('s2n-fb'), r, ST.tries.s2n, ruleFor(it), () => {
    const b = $('s2n-fb'); b.className = 'fb ok';
    b.innerHTML = `<div>${explain(it, LIST)}</div>` + (it.names.length > 1 ? `<div class="small">Also accepted: ${it.names.slice(1).join(', ')}</div>` : '');
    const nb = document.createElement('button'); nb.className = 'primary'; nb.textContent = 'Next structure →'; nb.onclick = () => s2nLoad(pick('s2n')); b.appendChild(nb);
  }); }
  updScore('s2n');
}

// ---------- Name → Structure ----------
function n2sLoad(it){
  if (!it) return;
  ST.cur.n2s = it; ST.tries.n2sf = 0; ST.tries.n2sd = 0;
  ST.promptName = it.names[Math.floor(Math.random()*it.names.length)];
  $('n2s-name').textContent = ST.promptName;
  $('n2s-fin').value = ''; $('n2s-reads').textContent = '';
  clearFb($('n2s-fb')); $('n2s-show').classList.add('hidden');
  edReset(); setMode(ST.mode);
}
function setMode(m){
  ST.mode = m;
  document.querySelectorAll('#n2s-mode button').forEach(b => b.classList.toggle('on', b.dataset.mode === m));
  $('n2s-fbox').classList.toggle('hidden', m !== 'formula');
  $('n2s-dbox').classList.toggle('hidden', m !== 'draw');
  clearFb($('n2s-fb'));
}
function n2sShow(){
  const it = ST.cur.n2s;
  $('n2s-show').classList.remove('hidden');
  renderStruct($('n2s-key'), keyDraw(it), {fit:true});
  $('n2s-keytext').innerHTML = `<div style="font-size:19px;margin-bottom:4px"><b>${dispFormula(it)}</b></div><div>${explain(it, LIST)}</div>` +
    (it.nBest > 1 ? `<div class="small">One of ${it.nBest} equivalent resonance structures is shown.</div>` : '');
}
function n2sFormulaCheck(){
  const it = ST.cur.n2s; if (!it) return;
  if (!$('n2s-fin').value.trim()) return;
  ST.score.n2s.t.add(it.id);
  const r = gradeFormula(it, $('n2s-fin').value, LIST);
  if (r.ok){ ST.score.n2s.s.add(it.id); feedback($('n2s-fb'), r, 0, null, null, [['Now draw its Lewis structure →', () => setMode('draw')], ['Next name →', () => n2sLoad(pick('n2s')), 'primary']]); }
  else { ST.tries.n2sf++; feedback($('n2s-fb'), r, ST.tries.n2sf, ruleFor(it), n2sShow); }
  updScore('n2s');
}
function n2sDrawCheck(){
  const it = ST.cur.n2s; if (!it) return;
  ST.score.n2s.t.add(it.id);
  const r = gradeLewis(it, {atoms: ED.atoms, bonds: ED.bonds, charge: ED.charge}, LIST);
  ED.bad = r.bad || []; edRender();
  if (r.ok){ ST.score.n2s.s.add(it.id); feedback($('n2s-fb'), r, 0, null, null, [['Next name →', () => n2sLoad(pick('n2s')), 'primary']]); }
  else { ST.tries.n2sd++; feedback($('n2s-fb'), r, ST.tries.n2sd, lewisRule(r.msg) || ruleFor(it), n2sShow); }
  updScore('n2s');
}
function lewisRule(msg){
  if (/valence electrons/.test(msg)) return 'Total valence electrons = sum of the valence electrons of all atoms, + 1 for each negative charge (− 1 for each positive charge).';
  if (/octet/.test(msg)) return 'Each atom other than H gets exactly 8 electrons (bonding pairs + lone pairs). H gets 2. Use single bonds to the central atom and formal charges rather than expanded octets.';
  if (/formal charge/i.test(msg)) return 'Formal charge = valence electrons − lone-pair electrons − number of bonds. The formal charges add up to the overall charge.';
  if (/skeleton|bonded to|bonded into/.test(msg)) return 'Skeleton: the central atom is bonded to every O. In acids and hydrogen anions, H atoms sit on O atoms.';
  return null;
}

// ---------- Lewis editor ----------
const ED = {atoms:[], bonds:[], charge:0, tool:'atom', el:'O', hist:[], bad:[], drag:null, move:null, temp:null};
const TOOLS = [
  ['atom', 'Atom', 'Click empty space to place the selected element. Click an atom to change it to the selected element.'],
  ['bond', 'Bond', 'Drag from one atom to another to bond them (or drag to empty space to add a bonded atom). Click a bond to cycle single → double → triple.'],
  ['lp+', 'Lone pair +', 'Click an atom to add a lone pair.'],
  ['lp-', 'Lone pair −', 'Click an atom to remove a lone pair.'],
  ['fc+', 'Charge +', 'Click an atom to raise its formal charge by 1.'],
  ['fc-', 'Charge −', 'Click an atom to lower its formal charge by 1.'],
  ['move', 'Move', 'Drag an atom to move it.'],
  ['erase', 'Erase', 'Click an atom or a bond to delete it.'],
];
const PAL = ['H','C','N','O','F','P','S','Cl','Br','I','Cr','Mn'];
function edSnap(){ ED.hist.push(JSON.stringify([ED.atoms, ED.bonds, ED.charge])); if (ED.hist.length > 150) ED.hist.shift(); }
function edReset(){ ED.atoms = []; ED.bonds = []; ED.charge = 0; ED.hist = []; ED.bad = []; $('ed-charge').value = '0'; edRender(); }
function edRender(){ renderStruct($('ed'), ED, {bad: ED.bad, temp: ED.temp, empty: 'Choose an element, then click here to place atoms.'}); }
function edPt(e){ const svg = $('ed'), p = svg.createSVGPoint(); p.x = e.clientX; p.y = e.clientY; return p.matrixTransform(svg.getScreenCTM().inverse()); }
function hitAtom(p){ let best = -1, bd = 20; ED.atoms.forEach((a,i) => { const d = Math.hypot(a.x - p.x, a.y - p.y); if (d < bd){ bd = d; best = i; } }); return best; }
function hitBond(p){
  let best = -1, bd = 8;
  ED.bonds.forEach((b,i) => {
    const A = ED.atoms[b.a], B = ED.atoms[b.b], dx = B.x - A.x, dy = B.y - A.y, L2 = dx*dx + dy*dy || 1;
    const t = Math.max(0, Math.min(1, ((p.x - A.x)*dx + (p.y - A.y)*dy)/L2));
    const d = Math.hypot(p.x - (A.x + t*dx), p.y - (A.y + t*dy));
    if (d < bd && t > 0.15 && t < 0.85){ bd = d; best = i; }
  });
  return best;
}
function findBond(a,b){ return ED.bonds.findIndex(x => (x.a === a && x.b === b) || (x.a === b && x.b === a)); }
function delAtom(i){ ED.bonds = ED.bonds.filter(b => b.a !== i && b.b !== i).map(b => ({a: b.a > i ? b.a-1 : b.a, b: b.b > i ? b.b-1 : b.b, order: b.order})); ED.atoms.splice(i,1); }
function addAtom(p){ ED.atoms.push({el: ED.el, x: Math.max(20, Math.min(680, p.x)), y: Math.max(20, Math.min(370, p.y)), lp:0, fc:0}); return ED.atoms.length - 1; }
function edDown(e){
  const svg = $('ed'); const p = edPt(e); const a = hitAtom(p); ED.bad = [];
  try { svg.setPointerCapture(e.pointerId); } catch(_){}
  const t = ED.tool;
  if (t === 'atom'){ edSnap(); if (a >= 0) ED.atoms[a].el = ED.el; else addAtom(p); }
  else if (t === 'bond'){
    if (a >= 0){ ED.drag = {from:a, x:p.x, y:p.y}; }
    else { const b = hitBond(p); if (b >= 0){ edSnap(); ED.bonds[b].order = ED.bonds[b].order % 3 + 1; } }
  }
  else if (t === 'lp+' || t === 'lp-' || t === 'fc+' || t === 'fc-'){
    if (a >= 0){ edSnap(); const A = ED.atoms[a];
      if (t === 'lp+') A.lp = Math.min(4, A.lp + 1); if (t === 'lp-') A.lp = Math.max(0, A.lp - 1);
      if (t === 'fc+') A.fc = Math.min(3, A.fc + 1); if (t === 'fc-') A.fc = Math.max(-3, A.fc - 1); }
  }
  else if (t === 'move'){ if (a >= 0){ edSnap(); ED.move = a; } }
  else if (t === 'erase'){
    if (a >= 0){ edSnap(); delAtom(a); }
    else { const b = hitBond(p); if (b >= 0){ edSnap(); ED.bonds.splice(b,1); } }
  }
  edRender();
}
function edMoveEv(e){
  if (!ED.drag && ED.move == null) return;
  const p = edPt(e);
  if (ED.drag){ const A = ED.atoms[ED.drag.from]; ED.temp = {x1:A.x, y1:A.y, x2:p.x, y2:p.y}; }
  if (ED.move != null){ ED.atoms[ED.move].x = Math.max(15, Math.min(685, p.x)); ED.atoms[ED.move].y = Math.max(15, Math.min(375, p.y)); }
  edRender();
}
function edUp(e){
  const p = edPt(e);
  if (ED.drag){
    const from = ED.drag.from, to = hitAtom(p);
    if (to >= 0 && to !== from){ edSnap(); const k = findBond(from, to); if (k >= 0) ED.bonds[k].order = ED.bonds[k].order % 3 + 1; else ED.bonds.push({a:from, b:to, order:1}); }
    else if (to < 0){
      const A = ED.atoms[from];
      if (Math.hypot(p.x - A.x, p.y - A.y) > 30){ edSnap(); const n = addAtom(p); ED.bonds.push({a:from, b:n, order:1}); }
    }
  }
  ED.drag = null; ED.temp = null; ED.move = null; edRender();
}
function edGiveAtoms(){
  const it = ST.cur.n2s; if (!it) return;
  edSnap(); ED.atoms = []; ED.bonds = []; ED.bad = [];
  const els = []; Object.keys(it.comp).forEach(e => { for (let k = 0; k < it.comp[e]; k++) els.push(e); });
  els.sort(() => Math.random() - 0.5);
  const n = els.length;
  els.forEach((el, i) => { const t = 2*Math.PI*i/n; ED.atoms.push({el, x: 350 + 150*Math.cos(t), y: 195 + 130*Math.sin(t), lp:0, fc:0}); });
  edRender();
}
function edInit(){
  const tb = $('tools');
  TOOLS.forEach(([k, lab, tip]) => {
    const b = document.createElement('button'); b.className = 'tool' + (k === ED.tool ? ' on' : ''); b.textContent = lab; b.dataset.k = k;
    b.onclick = () => { ED.tool = k; document.querySelectorAll('#tools .tool').forEach(x => x.classList.toggle('on', x.dataset.k === k)); $('tip').textContent = tip; };
    tb.appendChild(b);
  });
  $('tip').textContent = TOOLS[0][2];
  const pal = $('pal');
  const selEl = el => { ED.el = el; document.querySelectorAll('#pal button').forEach(x => x.classList.toggle('on', x.textContent === el)); if (ED.tool !== 'atom' && ED.tool !== 'bond') document.querySelector('#tools .tool[data-k="atom"]').click(); };
  PAL.forEach(el => { const b = document.createElement('button'); b.textContent = el; b.style.color = el === ED.el ? '' : (COLORS[el] || ''); b.onclick = () => selEl(el); if (el === ED.el) b.className = 'on'; pal.appendChild(b); });
  const other = document.createElement('input'); other.placeholder = 'other'; other.title = 'Type another element symbol';
  other.onchange = () => { const v = other.value.trim(); if (v in VAL){ ED.el = v; document.querySelectorAll('#pal button').forEach(x => x.classList.remove('on')); } else other.value = ''; };
  pal.appendChild(other);
  const svg = $('ed');
  svg.addEventListener('pointerdown', edDown); svg.addEventListener('pointermove', edMoveEv); svg.addEventListener('pointerup', edUp); svg.addEventListener('pointercancel', edUp);
  $('ed-charge').onchange = () => { edSnap(); ED.charge = +$('ed-charge').value; edRender(); };
  $('ed-undo').onclick = () => { const s = ED.hist.pop(); if (!s) return; const o = JSON.parse(s); ED.atoms = o[0]; ED.bonds = o[1]; ED.charge = o[2]; $('ed-charge').value = String(ED.charge); ED.bad = []; edRender(); };
  $('ed-clear').onclick = () => { edSnap(); ED.atoms = []; ED.bonds = []; ED.bad = []; edRender(); };
  $('ed-atoms').onclick = edGiveAtoms;
  $('ed-check').onclick = n2sDrawCheck;
  edRender();
}

// ---------- Review ----------
function reviewInit(){
  const anions = LIST.filter(x => x.kind === 'anion');
  const order = ['series','hanion','ide','thio','organic','metal'];
  anions.sort((a,b) => (order.indexOf(a.group) - order.indexOf(b.group)));
  let h = '<tr><th>Anion</th><th>Name</th><th>Acid</th><th>Acid name</th></tr>';
  anions.forEach(a => {
    const ac = a.acidItem;
    h += `<tr class="click" data-a="${a.id}" data-b="${ac ? ac.id : ''}"><td>${fmtU(a.formula)}</td><td>${a.names.join(' / ')}</td><td>${ac ? dispFormula(ac) : '—'}</td><td>${ac ? ac.names.join(' / ') : '—'}</td></tr>`;
  });
  $('rv-table').innerHTML = h;
  let hb = '<tr><th>Acid</th><th>Acid name</th><th>Anion</th><th>Anion name</th></tr>';
  LIST.filter(x => x.kind === 'binary').forEach(b => {
    const an = LIST.find(x => x.kind === 'anion' && x.formula === b.anionFormula);
    hb += `<tr class="click" data-a="${b.id}" data-b="${an ? an.id : ''}"><td>${dispFormula(b)}</td><td>${b.names.join(' / ')}</td><td>${b.anionFormula ? fmtU(b.anionFormula) : '—'}</td><td>${b.anionName || '—'}</td></tr>`;
  });
  $('rv-btable').innerHTML = hb;
  document.querySelectorAll('#tab-review tr.click').forEach(tr => tr.onclick = () => {
    document.querySelectorAll('#tab-review tr.sel').forEach(x => x.classList.remove('sel')); tr.classList.add('sel');
    const A = LIST.find(x => x.id === tr.dataset.a), B = LIST.find(x => x.id === tr.dataset.b);
    $('rv-hint').classList.add('hidden'); $('rv-pair').classList.remove('hidden');
    renderStruct($('rv-a'), keyDraw(A), {fit:true}); $('rv-acap').innerHTML = `${dispFormula(A)} · ${A.names[0]}`;
    if (B){ $('rv-b').parentNode.classList.remove('hidden'); renderStruct($('rv-b'), keyDraw(B), {fit:true}); $('rv-bcap').innerHTML = `${dispFormula(B)} · ${B.names[0]}`; }
    else $('rv-b').parentNode.classList.add('hidden');
  });
}

// ---------- Instructor ----------
const COLS = ['id','kind','formula','names','center','ligands','rule'];
function csvCell(v){ v = String(v == null ? '' : v); return /[",\n]/.test(v) ? '"' + v.replace(/"/g,'""') + '"' : v; }
function toCSV(rows){ return [COLS.join(',')].concat(rows.map(r => COLS.map(c => csvCell(r[c])).join(','))).join('\n') + '\n'; }
function rowOf(r){ return {id:r.id, kind:r.kind, formula:r.formula, names:r.names.join('; '), center: r.layout.special ? 'built-in' : r.layout.center, ligands: r.layout.special ? '' : r.layout.ligs.join(' '), rule: r.rule || ''}; }
function parseTable(text){
  text = text.replace(/^﻿/, '');
  const first = text.split(/\r?\n/)[0] || '', delim = first.includes('\t') ? '\t' : ',';
  const rows = []; let row = [], cell = '', q = false;
  for (let i = 0; i < text.length; i++){
    const c = text[i];
    if (q){ if (c === '"'){ if (text[i+1] === '"'){ cell += '"'; i++; } else q = false; } else cell += c; continue; }
    if (c === '"') q = true; else if (c === delim){ row.push(cell); cell = ''; }
    else if (c === '\n' || c === '\r'){ if (c === '\r' && text[i+1] === '\n') i++; row.push(cell); rows.push(row); row = []; cell = ''; }
    else cell += c;
  }
  if (cell || row.length){ row.push(cell); rows.push(row); }
  return rows.filter(r => r.some(x => x.trim()));
}
function download(name, text){
  try { const blob = new Blob([text], {type:'text/csv'}); const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click(); a.remove(); } catch(e){}
  $('in-copy').classList.remove('hidden'); $('in-copytext').value = text;
}
function handleUpload(text){
  const rows = parseTable(text), errs = $('in-errs'); errs.innerHTML = '';
  if (rows.length < 2){ $('in-status').textContent = 'The file has no data rows.'; return; }
  const head = rows[0].map(h => h.trim().toLowerCase());
  const miss = ['id','kind','formula','names'].filter(c => !head.includes(c));
  if (miss.length){ $('in-status').textContent = 'Missing column(s): ' + miss.join(', '); return; }
  const col = c => head.indexOf(c);
  const good = [], seen = new Set(); let nBad = 0;
  rows.slice(1).forEach((r, k) => {
    const get = c => col(c) >= 0 ? (r[col(c)] || '').trim() : '';
    const line = k + 2, li = document.createElement('li'), id = get('id');
    try {
      if (!id) throw new Error('missing id');
      if (seen.has(id)) throw new Error('duplicate id ' + id); seen.add(id);
      const kind = get('kind').toLowerCase();
      if (!['anion','binary','oxoacid'].includes(kind)) throw new Error('kind must be anion, binary or oxoacid');
      const names = get('names').split(';').map(s => s.trim().toLowerCase()).filter(Boolean);
      if (!names.length) throw new Error('missing names');
      const center = get('center'); let layout;
      if (!center || center.toLowerCase() === 'built-in'){
        const b = RAW_SPECIES.find(x => x.id === id);
        if (!b) throw new Error('"built-in" (or a blank center) only works for an id that is already built in');
        layout = b.layout;
      } else {
        const ligs = get('ligands').split(/\s+/).filter(Boolean);
        if (!(center in VAL)) throw new Error(`center "${center}" isn't an element the tool can draw`);
        if (ligs.length < 1 || ligs.length > 4) throw new Error('ligands: give 1 to 4 atoms around the center');
        for (const L0 of ligs){ const e = /H$/.test(L0) && L0.length > 1 ? L0.slice(0,-1) : L0; if (!(e in VAL)) throw new Error(`ligand "${L0}" isn't understood`); }
        layout = {center, ligs};
      }
      const raw = {id, kind, formula: get('formula'), names, layout};
      const base = RAW_SPECIES.find(x => x.id === id);
      if (base){ ['anionFormula','anionName','root','group','family','explain'].forEach(f => { if (base[f] != null) raw[f] = base[f]; }); if (base.kind !== kind) { delete raw.group; delete raw.explain; } }
      if (get('rule')) raw.rule = get('rule');
      prepareSpecies(raw);
      good.push(raw); li.className = 'ok'; li.textContent = `Row ${line} (${id}): OK — ${names[0]}`;
    } catch(e){ nBad++; li.className = 'bad'; li.textContent = `Row ${line}${id ? ' (' + id + ')' : ''}: ${e.message}`; }
    errs.appendChild(li);
  });
  if (good.length){
    let prev = []; try { prev = JSON.parse(store(STORE_KEY) || '[]'); } catch(e){}
    const m = new Map(prev.map(r => [r.id, r])); good.forEach(r => m.set(r.id, r));
    store(STORE_KEY, JSON.stringify([...m.values()]));
    refreshAll();
  }
  $('in-status').textContent = `${good.length} row(s) added or updated, ${nBad} row(s) skipped. Uploaded species are saved in this browser only.`;
}
function instrInit(){
  $('in-tpl').onclick = () => download('polyatomic-template.csv', toCSV([
    {id:'SeO4', kind:'anion', formula:'SeO4^2-', names:'selenate', center:'Se', ligands:'O O O O', rule:''},
    {id:'H2SeO4', kind:'oxoacid', formula:'H2SeO4', names:'selenic acid', center:'Se', ligands:'O OH O OH', rule:''},
    {id:'H2Se', kind:'binary', formula:'H2Se', names:'hydroselenic acid', center:'Se', ligands:'H H', rule:''},
  ]));
  $('in-cur').onclick = () => download('polyatomic-species.csv', toCSV(LIST.map(rowOf)));
  $('in-file').onchange = e => { const f = e.target.files[0]; if (!f) return; const rd = new FileReader(); rd.onload = () => handleUpload(String(rd.result)); rd.readAsText(f); e.target.value = ''; };
  $('in-remove').onclick = () => { store(STORE_KEY, null); refreshAll(); $('in-errs').innerHTML = ''; $('in-status').textContent = 'Uploaded species removed; the built-in list is back.'; };
}
function refreshAll(){
  loadSpecies(); drawFilters('s2n'); drawFilters('n2s'); reviewInit();
  s2nLoad(pick('s2n')); n2sLoad(pick('n2s'));
  $('in-status').textContent = LIST.uploaded ? `${LIST.uploaded} uploaded species in use in this browser.` : '';
}

// ---------- init ----------
function init(){
  document.querySelectorAll('#tabs button').forEach(b => b.onclick = () => {
    document.querySelectorAll('#tabs button').forEach(x => x.classList.toggle('on', x === b));
    document.querySelectorAll('.tab').forEach(t => t.classList.toggle('on', t.id === 'tab-' + b.dataset.tab));
  });
  $('s2n-next').onclick = () => s2nLoad(pick('s2n'));
  $('s2n-check').onclick = s2nCheck;
  $('s2n-in').addEventListener('keydown', e => { if (e.key === 'Enter') s2nCheck(); });
  $('s2n-showf').onchange = () => { const it = ST.cur.s2n; $('s2n-formula').textContent = $('s2n-showf').checked && it ? dispFormula(it) : ''; };
  $('n2s-next').onclick = () => n2sLoad(pick('n2s'));
  document.querySelectorAll('#n2s-mode button').forEach(b => b.onclick = () => setMode(b.dataset.mode));
  $('n2s-fcheck').onclick = n2sFormulaCheck;
  $('n2s-fin').addEventListener('keydown', e => { if (e.key === 'Enter') n2sFormulaCheck(); });
  $('n2s-fin').addEventListener('input', () => {
    const v = $('n2s-fin').value; if (!v.trim()){ $('n2s-reads').textContent = ''; return; }
    const p = parseFormula(v); $('n2s-reads').textContent = p.err ? p.err : 'Reads as: ' + readsAs(p);
  });
  edInit(); instrInit(); refreshAll(); updScore('s2n'); updScore('n2s');
}
init();
