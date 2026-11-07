// ===== Engine: structures, Lewis solver, graders =====
const VAL = {H:1,B:3,C:4,N:5,O:6,F:7,Si:4,P:5,S:6,Cl:7,As:5,Se:6,Br:7,Te:6,I:7,Cr:6,Mn:7,Mo:6,W:6};
const EN  = {H:2.20,B:2.04,C:2.55,N:3.04,O:3.44,F:3.98,Si:1.90,P:2.19,S:2.58,Cl:3.16,As:2.18,Se:2.55,Br:2.96,Te:2.10,I:2.66,Cr:1.66,Mn:1.55,Mo:2.16,W:2.36};
const OTHER_ELEMENTS = ['Li','Na','K','Rb','Cs','Be','Mg','Ca','Sr','Ba','Al','Fe','Co','Ni','Cu','Zn','Ag','Au','Hg','Pb','Sn','Ti','V','Ga','Ge','Sb','Bi','Ne','Ar','He','Kr','Xe'];
const SUBS = '₀₁₂₃₄₅₆₇₈₉', SUPS = '⁰¹²³⁴⁵⁶⁷⁸⁹';

// ---------- formatting ----------
function fmtU(f){
  const k = f.indexOf('^');
  let b = k < 0 ? f : f.slice(0,k), c = k < 0 ? null : f.slice(k+1);
  let s = b.replace(/\d/g, d => SUBS[+d]);
  if (c){
    const m = c.match(/^(\d*)([+-])$/);
    if (m) s += [...m[1]].map(x => SUPS[+x]).join('') + (m[2] === '-' ? '⁻' : '⁺');
  }
  return s;
}
function chargeText(q){ if (!q) return ''; const a = Math.abs(q); return (a > 1 ? a : '') + (q < 0 ? '−' : '+'); }
function fcText(q){ return q > 0 ? '+' + q : '−' + (-q); }
function dispFormula(it){ return fmtU(it.formula) + (it.kind !== 'anion' ? '(aq)' : ''); }

// ---------- skeleton builders (math coords, bond length 1, y up) ----------
function oxo(center, ligs, ang){
  const n = ligs.length;
  const A = ang || ({1:[180],2:[210,330],3:[90,210,330],4:[90,0,270,180]})[n];
  const atoms = [{el:center,x:0,y:0}], bonds = [];
  ligs.forEach((L,k) => {
    const t = A[k]*Math.PI/180, dx = Math.cos(t), dy = Math.sin(t);
    if (L.length > 1 && /H$/.test(L)){
      atoms.push({el:L.slice(0,-1), x:dx, y:dy}); const oi = atoms.length-1; bonds.push([0,oi]);
      atoms.push({el:'H', x:2*dx, y:2*dy}); bonds.push([oi, atoms.length-1]);
    } else { atoms.push({el:L, x:dx, y:dy}); bonds.push([0, atoms.length-1]); }
  });
  return {atoms, bonds};
}
function chain(els){ return {atoms: els.map((el,i)=>({el, x:i, y:0})), bonds: els.slice(1).map((_,i)=>[i,i+1])}; }
function manual(list, bonds){ return {atoms: list.map(([el,x,y])=>({el,x,y})), bonds}; }
const SPECIALS = {
  OH:   () => chain(['O','H']),
  CN:   () => chain(['C','N']),
  O2:   () => chain(['O','O']),
  SCN:  () => chain(['S','C','N']),
  HCN:  () => chain(['H','C','N']),
  acetate: () => manual([['C',-1,0],['C',0,0],['O',0.5,0.87],['O',0.5,-0.87],['H',-1,1],['H',-2,0],['H',-1,-1]],
                        [[0,1],[1,2],[1,3],[0,4],[0,5],[0,6]]),
  acetic:  () => manual([['C',-1,0],['C',0,0],['O',0.5,0.87],['O',0.5,-0.87],['H',-1,1],['H',-2,0],['H',-1,-1],['H',1.5,-0.87]],
                        [[0,1],[1,2],[1,3],[0,4],[0,5],[0,6],[3,7]]),
  oxalate: () => manual([['C',-0.5,0],['C',0.5,0],['O',-1,0.87],['O',-1,-0.87],['O',1,0.87],['O',1,-0.87]],
                        [[0,1],[0,2],[0,3],[1,4],[1,5]]),
  oxalic:  () => manual([['C',-0.5,0],['C',0.5,0],['O',-1,0.87],['O',-1,-0.87],['O',1,0.87],['O',1,-0.87],['H',-2,-0.87],['H',2,-0.87]],
                        [[0,1],[0,2],[0,3],[1,4],[1,5],[3,6],[5,7]]),
  Cr2O7:   () => manual([['Cr',-1,0],['O',0,0],['Cr',1,0],['O',-1,1],['O',-2,0],['O',-1,-1],['O',1,1],['O',2,0],['O',1,-1]],
                        [[0,1],[1,2],[0,3],[0,4],[0,5],[2,6],[2,7],[2,8]]),
};
function buildSkel(layout){
  if (layout.special){ if (!SPECIALS[layout.special]) throw new Error('unknown layout ' + layout.special); return SPECIALS[layout.special](); }
  return oxo(layout.center, layout.ligs, layout.ang);
}
function compOf(skel){ const c = {}; skel.atoms.forEach(a => c[a.el] = (c[a.el]||0) + 1); return c; }
function compEq(a,b){ const k = new Set([...Object.keys(a),...Object.keys(b)]); for (const e of k) if ((a[e]||0) !== (b[e]||0)) return false; return true; }
function compText(c){
  const order = Object.keys(c).sort((x,y) => (x==='C'?-2:x==='H'?-1:0) - (y==='C'?-2:y==='H'?-1:0) || x.localeCompare(y));
  return order.map(e => e + (c[e] > 1 ? c[e] : '')).join('');
}

// ---------- all-octet Lewis solver ----------
function solveOctet(skel, charge){
  const atoms = skel.atoms, B = skel.bonds;
  for (const a of atoms) if (!(a.el in VAL)) return null;
  const N = atoms.reduce((s,a) => s + VAL[a.el], 0) - charge;
  const free = []; B.forEach((b,i) => { if (atoms[b[0]].el !== 'H' && atoms[b[1]].el !== 'H') free.push(i); });
  const ord = B.map(() => 1); const all = [];
  function evalCur(){
    const bs = atoms.map(() => 0); B.forEach((b,i) => { bs[b[0]] += ord[i]; bs[b[1]] += ord[i]; });
    const lp = [], fc = [];
    for (let k = 0; k < atoms.length; k++){
      const el = atoms[k].el;
      if (el === 'H'){ if (bs[k] !== 1) return; lp.push(0); }
      else { if (bs[k] > 4) return; lp.push(4 - bs[k]); }
      fc.push(VAL[el] - 2*lp[k] - bs[k]);
    }
    const tot = 2*ord.reduce((s,o)=>s+o,0) + 2*lp.reduce((s,x)=>s+x,0);
    if (tot !== N) return;
    all.push({ord: ord.slice(), lp, fc,
      sumAbs: fc.reduce((s,f)=>s+Math.abs(f),0),
      en: fc.reduce((s,f,k)=>s + f*EN[atoms[k].el], 0)});
  }
  (function rec(i){ if (i === free.length){ evalCur(); return; } for (let o = 1; o <= 3; o++){ ord[free[i]] = o; rec(i+1); } ord[free[i]] = 1; })(0);
  if (!all.length) return null;
  let best = all[0];
  for (const r of all) if (r.sumAbs < best.sumAbs || (r.sumAbs === best.sumAbs && r.en < best.en - 1e-9)) best = r;
  const nBest = all.filter(r => r.sumAbs === best.sumAbs && Math.abs(r.en - best.en) < 1e-9).length;
  return {best, N, minAbs: best.sumAbs, bestEn: best.en, nBest};
}

// ---------- formula parsing ----------
function parseBody(s){
  if (!s) return {err:'Type a formula.'};
  let i = 0; const stack = [{}];
  const add = (o,e,n) => { o[e] = (o[e]||0) + n; };
  while (i < s.length){
    const c = s[i];
    if (c === ' '){ i++; continue; }
    if (c === '(' || c === '['){ stack.push({}); i++; continue; }
    if (c === ')' || c === ']'){
      if (stack.length < 2) return {err:'There is a closing parenthesis with no opening one.'};
      i++; let n = ''; while (/\d/.test(s[i]||'')) n += s[i++];
      const top = stack.pop(), k = n ? +n : 1;
      for (const e in top) add(stack[stack.length-1], e, top[e]*k);
      continue;
    }
    if (/[A-Z]/.test(c)){
      let e = c; i++;
      if (/[a-z]/.test(s[i]||'')){ e += s[i]; i++; }
      if (!(e in VAL) && !OTHER_ELEMENTS.includes(e)){
        if (e.length === 2 && (e[0] in VAL)) return {err:`"${e}" isn't an element symbol. (Element symbols are one capital letter, or a capital and a lowercase letter.)`};
        return {err:`"${e}" isn't an element symbol. Check capital letters: Cl, not CL or cl.`};
      }
      let n = ''; while (/\d/.test(s[i]||'')) n += s[i++];
      add(stack[stack.length-1], e, n ? +n : 1); continue;
    }
    if (/[a-z]/.test(c)) return {err:'Element symbols start with a capital letter (Cl, not cl).'};
    if (/\d/.test(c)) return {err:'Numbers go after an element symbol, as subscripts (H2SO4).'};
    return {err:`The character "${c}" can't be read in a formula.`};
  }
  if (stack.length > 1) return {err:'A parenthesis is not closed.'};
  if (!Object.keys(stack[0]).length) return {err:'Type a formula.'};
  return stack[0];
}
function parseFormula(raw){
  let s = String(raw||'').trim();
  s = s.replace(/[₀-₉]/g, c => SUBS.indexOf(c));
  s = s.replace(/([⁰¹²³⁴⁵⁶⁷⁸⁹]*)([⁺⁻])\s*$/, (m,d,sg) => '^' + [...d].map(c => SUPS.indexOf(c)).join('') + (sg === '⁺' ? '+' : '-'));
  s = s.replace(/[−–—]/g, '-').replace(/\((aq|g|l|s)\)/gi, '').replace(/\s+(ion|anion)$/i, '').trim();
  if (!s) return {err:'Type a formula.'};
  let charge = 0, body = s, alt = null, explicit = false, m;
  if ((m = s.match(/^(.*?)\s*\^\s*(\d*)\s*([+-])$/))){ body = m[1]; charge = (m[2] ? +m[2] : 1)*(m[3]==='-'?-1:1); explicit = true; }
  else if ((m = s.match(/^(.*\S)\s+(\d*)([+-])$/))){ body = m[1]; charge = (m[2] ? +m[2] : 1)*(m[3]==='-'?-1:1); explicit = true; }
  else if ((m = s.match(/^(.*?)([+-])$/))){
    body = m[1]; charge = m[2] === '-' ? -1 : 1;
    if (/[A-Za-z)\]]\d+$/.test(body)){ const d = +body.slice(-1); if (d > 1) alt = {body: body.slice(0,-1), charge: d*charge}; }
  }
  else if (/\^/.test(s)) return {err:'After ^, type the charge, e.g. ^2- or ^-.'};
  body = body.trim();
  const comp = parseBody(body); if (comp.err) return comp;
  const res = {comp, charge, body, explicit};
  if (alt){ const c2 = parseBody(alt.body); if (!c2.err) res.alt = {comp:c2, charge:alt.charge, body:alt.body}; }
  return res;
}
function readsAs(p){
  if (p.err) return '';
  return fmtU(compTextOrdered(p.body) + (p.charge ? '^' + (Math.abs(p.charge) > 1 ? Math.abs(p.charge) : '') + (p.charge < 0 ? '-' : '+') : ''));
}
function compTextOrdered(body){ return body.replace(/\s+/g,''); }

// ---------- species preparation ----------
function prepareSpecies(raw){
  const it = Object.assign({}, raw);
  const p = parseFormula(it.formula);
  if (p.err) throw new Error('formula: ' + p.err);
  it.charge = p.charge;
  it.skel = buildSkel(it.layout);
  it.comp = compOf(it.skel);
  if (!compEq(it.comp, p.comp)) throw new Error(`layout atoms (${compText(it.comp)}) don't match formula ${it.formula}`);
  const sol = solveOctet(it.skel, it.charge);
  if (!sol) throw new Error('no all-octet Lewis structure fits this skeleton and charge');
  Object.assign(it, {key: sol.best, N: sol.N, minAbs: sol.minAbs, bestEn: sol.bestEn, nBest: sol.nBest});
  if (it.kind === 'anion' && it.charge >= 0) throw new Error('an anion needs a negative charge');
  if (it.kind !== 'anion' && it.charge !== 0) throw new Error('acids must be neutral');
  const n0 = it.names[0];
  it.nO = it.comp.O || 0; it.nH = it.comp.H || 0;
  it.family = it.family || it.layout.center || it.layout.special;
  if (it.kind === 'anion') it.suf = /ate$/.test(n0) ? 'ate' : /ite$/.test(n0) ? 'ite' : /ide$/.test(n0) ? 'ide' : '';
  else it.suf = /ic acid$/.test(n0) ? 'ic' : /ous acid$/.test(n0) ? 'ous' : '';
  const core = n0.replace(/^(di|mono)?hydrogen\s*/,'');
  it.pre = /^per(?!ox)/.test(core) ? 'per' : /^hypo/.test(core) ? 'hypo' : '';
  if (!it.group) it.group = it.kind === 'binary' ? 'binary' : it.kind === 'oxoacid' ? 'oxoacid' : (it.nH ? 'hanion' : 'series');
  return it;
}
function linkSpecies(list){
  const nonH = c => { const o = Object.assign({}, c); delete o.H; return o; };
  for (const a of list){ a.parent = null; }
  for (const acid of list.filter(x => x.kind === 'oxoacid')){
    const cands = list.filter(x => x.kind === 'anion' && compEq(nonH(x.comp), nonH(acid.comp)) && x.nH - x.charge === acid.nH).sort((x,y) => x.nH - y.nH);
    if (cands.length){ acid.parent = cands[0]; if (!cands[0].acidItem) cands[0].acidItem = acid; }
  }
  for (const a of list.filter(x => x.kind === 'anion' && x.nH > 0)){
    a.hparent = list.find(x => x.kind === 'anion' && x !== a && x.nH === 0 && compEq(nonH(x.comp), nonH(a.comp))) || null;
  }
  return list;
}
function ateOf(it, list){
  return list.find(o => o.kind === 'anion' && o.family === it.family && o.suf === 'ate' && !o.pre && !o.nH) || null;
}

// ---------- name normalization & diagnosis ----------
function normName(s){
  return String(s||'').toLowerCase().replace(/\((aq)\)/g,'').replace(/\b(ion|anion)\b/g,'').replace(/[^a-z]/g,'');
}
const RULES = {
  series: 'Oxoanion series: memorize the -ate ion (NO₃⁻, SO₄²⁻, PO₄³⁻, CO₃²⁻, ClO₃⁻). One fewer O → -ite. One more O → per-…-ate. Two fewer O → hypo-…-ite. The charge is the same across a series, and Br and I follow Cl.',
  hanion: 'Adding H⁺ to an anion adds "hydrogen" to its name and raises the charge by 1 (two H⁺ → "dihydrogen"). The older prefix bi- means one H: bicarbonate = hydrogen carbonate.',
  ide: 'A few polyatomic anions end in -ide, like monatomic anions do: hydroxide OH⁻, cyanide CN⁻, peroxide O₂²⁻.',
  thio: 'Thio- means an S atom takes the place of an O: sulfate SO₄²⁻ → thiosulfate S₂O₃²⁻; cyanate OCN⁻ → thiocyanate SCN⁻.',
  organic: 'Two carbon-containing anions to memorize: acetate C₂H₃O₂⁻ (CH₃COO⁻) and oxalate C₂O₄²⁻.',
  metal: 'Transition-metal oxoanions to memorize: chromate CrO₄²⁻, dichromate Cr₂O₇²⁻, permanganate MnO₄⁻.',
  binary: 'Binary acids (H plus one nonmetal, dissolved in water): hydro- + root of the nonmetal + -ic + acid. HCl(aq) → hydrochloric acid; H₂S(aq) → hydrosulfuric acid.',
  hcn: 'HCN(aq) is named like a binary acid because its anion, cyanide, ends in -ide: hydro- + cyan + -ic acid.',
  oxoacid: 'Oxoacids: start from the anion and change the ending: -ate → -ic acid, -ite → -ous acid. Keep per- and hypo-; never use hydro-. Add one H for each negative charge. Sulfur keeps "ur" (sulfuric, sulfurous); phosphorus gives phosphoric.',
};
function ruleFor(it){ return it.rule || RULES[it.group] || RULES[it.kind === 'binary' ? 'binary' : it.kind === 'oxoacid' ? 'oxoacid' : 'series']; }
function cap(s){ return s.charAt(0).toUpperCase() + s.slice(1); }

function explain(it, list){
  const f = dispFormula(it), n0 = it.names[0];
  if (it.explain) return it.explain;
  if (it.group === 'series'){
    const ate = ateOf(it, list);
    if (!ate || ate === it) return `${f} is <b>${n0}</b>, one of the -ate ions to memorize.`;
    const d = it.nO - ate.nO, af = fmtU(ate.formula);
    if (d === 1)  return `${f} has one more O than ${ate.names[0]} (${af}), so it gets per-: <b>${n0}</b>.`;
    if (d === -1) return `${f} has one fewer O than ${ate.names[0]} (${af}), so the ending is -ite: <b>${n0}</b>.`;
    if (d === -2) return `${f} has two fewer O than ${ate.names[0]} (${af}), so it is hypo-…-ite: <b>${n0}</b>.`;
  }
  if (it.group === 'hanion' && it.hparent){
    const p = it.hparent;
    return `${fmtU(p.formula)} (${p.names[0]}) + ${it.nH === 1 ? '' : it.nH + ' '}H⁺ → ${f}: <b>${n0}</b>` + (it.names.length > 1 ? ` (also ${it.names.slice(1).join(', ')})` : '') + '.';
  }
  if (it.kind === 'binary' && it.anionName){
    return `${fmtU(it.anionFormula)} is ${it.anionName}. Dissolved in water, ${f} is an acid: hydro- + ${it.root} + -ic acid → <b>${n0}</b>.`;
  }
  if (it.kind === 'oxoacid' && it.parent){
    const p = it.parent, q = -p.charge;
    return `Parent anion: ${p.names[0]} (${fmtU(p.formula)}). ${p.suf === 'ite' ? '-ite → -ous' : '-ate → -ic'}, and ${q} H⁺ cancel${q === 1 ? 's' : ''} the ${chargeText(p.charge)} charge → ${f}: <b>${n0}</b>.`;
  }
  return `${f} is <b>${n0}</b>. ${ruleFor(it)}`;
}

function compareTo(it, o){
  if (o.family === it.family){
    if (o.nO !== it.nO) return `That one has ${o.nO} O atom${o.nO===1?'':'s'}; this one has ${it.nO}.`;
    if (o.nH !== it.nH) return `That one has ${o.nH} H atom${o.nH===1?'':'s'}; this one has ${it.nH}.`;
  }
  return 'Compare the formulas carefully.';
}

function diagnoseName(it, ans, list){
  const bad = msg => ({ok:false, msg});
  const a = normName(ans);
  if (!a) return bad('Type a name first.');
  const hit = it.names.find(n => normName(n) === a);
  if (hit){
    const others = it.names.filter(n => n !== hit);
    return {ok:true, msg:`Correct: <b>${hit}</b>.` + (others.length ? ` Also called ${others.join(' or ')}.` : '')};
  }
  const isAcid = it.kind !== 'anion', hasAcid = /acid$/.test(a), body = hasAcid ? a.slice(0,-4) : a;
  const other = list.find(o => o !== it && o.names.some(n => normName(n) === a));
  const binLike = it.kind === 'binary';
  if (isAcid && !hasAcid){
    if (binLike && /^hydrogen/.test(body)) return bad(`"${esc(ans.trim())}" is the name of the pure compound. In water, ${dispFormula(it)} is named as an acid: hydro- + root + -ic acid.`);
    if (binLike && it.anionName && normName(it.anionName) === body) return bad(`${cap(it.anionName)} is the anion ${fmtU(it.anionFormula)}. ${dispFormula(it)} is the acid, and acid names end in "acid".`);
    if (other && other.kind === 'anion') return bad(`${cap(other.names[0])} is the anion ${fmtU(other.formula)}. This species is neutral, with H attached; in water it's an acid, so its name ends in "acid".`);
    return bad('This is an acid: a neutral species with H in front, dissolved in water. Acid names end with the word "acid".');
  }
  if (!isAcid && hasAcid){
    if (other) return bad(`${cap(other.names[0])} is ${dispFormula(other)}, a neutral acid. This species has a ${chargeText(it.charge)} charge, so it's an anion, and anion names don't include "acid".`);
    return bad(`This species has a ${chargeText(it.charge)} charge, so it's an anion, not an acid. Anion names don't include "acid".`);
  }
  if (other) return bad(`${cap(other.names[0])} is ${dispFormula(other)}. ${compareTo(it, other)}`);
  if (isAcid){
    if (binLike && !/^hydro/.test(body)) return bad('Binary acids begin with hydro-: hydro- + root + -ic acid.');
    if (!binLike && it.group !== 'hcn' && /^hydro/.test(body)) return bad('Oxoacids don\'t use hydro-. That prefix is only for acids with no oxygen (binary acids and HCN).');
    if (it.group === 'hcn' && !/^hydro/.test(body)) return bad('HCN(aq) is named like a binary acid, so it begins with hydro-.');
    if (/(ate|ite|ide)$/.test(body)) return bad('In an acid name the anion ending changes: -ate → -ic, -ite → -ous, -ide → hydro-…-ic.');
    if (it.suf === 'ic' && /ous$/.test(body)) return bad(binLike ? 'Binary acids end in -ic.' : `The parent anion is ${it.parent ? it.parent.names[0] : 'an -ate ion'}, an -ate ion, so the acid ends in -ic, not -ous.`);
    if (it.suf === 'ous' && /ic$/.test(body)) return bad(`The parent anion is ${it.parent ? it.parent.names[0] : 'an -ite ion'}, an -ite ion, so the acid ends in -ous, not -ic.`);
  } else {
    if (it.suf === 'ide' && /(ate|ite)$/.test(body)) return bad('This one is an exception: only a few polyatomic anions end in -ide (hydroxide, cyanide, peroxide).');
    if (it.suf === 'ate' && /ite$/.test(body)) return bad(`-ite means one fewer O than the -ate ion. Count the O atoms: this one has ${it.nO}.`);
    if (it.suf === 'ite' && /ate$/.test(body)) return bad(`Count the O atoms: this one has ${it.nO}. Is that the -ate ion, or one fewer?`);
    if (it.suf !== 'ide' && /ide$/.test(body)) return bad('-ide is for monatomic anions (and hydroxide, cyanide, peroxide). An oxoanion ends in -ate or -ite.');
    if (it.nH > 0 && !/hydrogen|^bi/.test(body)) return bad(`This anion still has ${it.nH === 1 ? 'one H' : 'two H atoms'} attached, and the name has to say so (${it.nH === 1 ? 'hydrogen' : 'dihydrogen'} …).`);
    if (it.nH === 2 && !/dihydrogen/.test(body)) return bad('Two H atoms attached → dihydrogen.');
    if (it.nH === 1 && /dihydrogen/.test(body)) return bad('Only one H is attached → hydrogen, not dihydrogen.');
    if (it.nH === 0 && /hydrogen/.test(body)) return bad('This anion has no H atoms, so there is no "hydrogen" in its name.');
  }
  const core = body.replace(/^(di|mono)?hydrogen/,'').replace(/^bi/,'');
  const cPer = /^per(?!ox)/.test(core), cHypo = /^hypo/.test(core);
  const ate = ateOf(it.parent || it, list);
  const ateTxt = ate ? ` (${ate.names[0]}, ${fmtU(ate.formula)})` : '';
  if (it.pre === 'per' && !cPer) return bad(`This has one more O than the -ate ion${ateTxt}, so it takes the prefix per-.`);
  if (it.pre === 'hypo' && !cHypo) return bad(`This has two fewer O than the -ate ion${ateTxt}, so it takes the prefix hypo-.`);
  if (it.pre !== 'per' && cPer) return bad(`per- is only for one more O than the -ate ion${ateTxt}.`);
  if (it.pre !== 'hypo' && cHypo) return bad(`hypo- is only for two fewer O than the -ate ion${ateTxt}.`);
  if (/sulf(ic|ous)$/.test(body)) return bad('Sulfur acids keep the "ur": sulfuric, sulfurous, hydrosulfuric.');
  if (/phosph(ic|ous)$/.test(body)) return bad('The phosphorus acid is phosphoric acid.');
  return bad('Not quite. Check the root (it comes from the central element), any prefix, and the ending.');
}
function esc(s){ return String(s).replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'})[c]); }

// ---------- formula grading ----------
function gradeFormula(it, text, list){
  const bad = msg => ({ok:false, msg});
  const p = parseFormula(text);
  if (p.err) return bad(p.err);
  const same = (c,q) => compEq(c, it.comp) && q === it.charge;
  if (same(p.comp, p.charge)){
    return {ok:true, msg:`Correct: ${dispFormula(it)}.` + (it.kind !== 'anion' ? ' (Writing (aq) is optional here.)' : '')};
  }
  if (p.alt && same(p.alt.comp, p.alt.charge)) return bad(`Did you mean ${fmtU(it.formula)}? Put the charge after ^ (type ${it.formula}). As typed, the last digit reads as a subscript.`);
  const other = list.find(o => o !== it && compEq(o.comp, p.comp) && o.charge === p.charge);
  if (other) return bad(`${dispFormula(other)} is ${other.names[0]}. ${compareTo(it, other)}`);
  if (compEq(p.comp, it.comp)){
    if (it.charge === 0) return bad('Right atoms, but an acid is a neutral molecule: no charge.');
    if (p.charge === 0) return bad('Right atoms. An anion needs its charge: type it after ^, for example ^- or ^2-.');
    return bad('Right atoms, wrong charge.' + (it.kind === 'anion' && it.nH ? ' Each H⁺ added to an anion raises its charge by 1.' : ''));
  }
  const need = Object.keys(it.comp), have = Object.keys(p.comp);
  const missing = need.filter(e => !p.comp[e]), extra = have.filter(e => !it.comp[e]);
  if (missing.length || extra.length){
    let m = 'Check the elements.';
    if (missing.length) m += ` Missing: ${missing.join(', ')}.`;
    if (extra.length) m += ` Not in this species: ${extra.join(', ')}.`;
    return bad(m);
  }
  const diffs = need.filter(e => p.comp[e] !== it.comp[e]);
  if (diffs.length === 1 && diffs[0] === 'O'){
    return bad(`Right elements, but the number of O atoms is off (you have ${p.comp.O}).` + (it.suf === 'ite' || it.suf === 'ous' ? ' -ite / -ous means one fewer O than the -ate ion.' : it.pre ? ' per- and hypo- change the O count relative to the -ate ion.' : ''));
  }
  if (diffs.length === 1 && diffs[0] === 'H'){
    if (it.kind === 'oxoacid' || it.kind === 'binary') return bad(`Check the number of H atoms (you have ${p.comp.H}). An acid has one H for each negative charge on its anion.`);
    return bad(`Check the number of H atoms (you have ${p.comp.H}): hydrogen = 1 H, dihydrogen = 2 H.`);
  }
  if (it.charge !== p.charge) return bad('Check the subscripts and the charge.');
  return bad('Check the subscripts.');
}

// ---------- Lewis structure grading ----------
function iso(elsA, adjA, elsB, adjB){
  const n = elsA.length; if (n !== elsB.length) return null;
  const dA = adjA.map(s => s.size), dB = adjB.map(s => s.size);
  const order = [...Array(n).keys()].sort((x,y) => dA[y] - dA[x]);
  const map = Array(n).fill(-1), inv = Array(n).fill(-1);
  function bt(k){
    if (k === n) return true;
    const a = order[k];
    for (let b = 0; b < n; b++){
      if (inv[b] >= 0 || elsB[b] !== elsA[a] || dB[b] !== dA[a]) continue;
      let ok = true;
      for (const a2 of adjA[a]) if (map[a2] >= 0 && !adjB[b].has(map[a2])){ ok = false; break; }
      if (ok) for (const b2 of adjB[b]) if (inv[b2] >= 0 && !adjA[a].has(inv[b2])){ ok = false; break; }
      if (!ok) continue;
      map[a] = b; inv[b] = a;
      if (bt(k+1)) return true;
      map[a] = -1; inv[b] = -1;
    }
    return false;
  }
  return bt(0) ? map : null;
}
function atomDesc(d, i){
  const a = d.atoms[i];
  const nb = d.bonds.filter(b => b.a === i || b.b === i).map(b => d.atoms[b.a === i ? b.b : b.a].el);
  if (!nb.length) return `the unbonded ${a.el}`;
  const cnt = {}; nb.forEach(e => cnt[e] = (cnt[e]||0)+1);
  return `the ${a.el} bonded to ` + Object.keys(cnt).map(e => (cnt[e] > 1 ? cnt[e] + ' ' : '') + e).join(' and ');
}
function gradeLewis(it, d, list){
  const bad = (msg, atoms) => ({ok:false, msg, bad: atoms || []});
  const n = d.atoms.length;
  if (!n) return bad('Draw the structure first: place atoms, bond them, then add lone pairs and formal charges.');
  // 1. atoms
  const comp = {}; d.atoms.forEach(a => comp[a.el] = (comp[a.el]||0)+1);
  if (!compEq(comp, it.comp)){
    const other = list.find(o => o !== it && compEq(o.comp, comp));
    if (other) return bad(`Those atoms make ${dispFormula(other)} (${other.names[0]}). Check the formula of ${it.names[0]} first.`);
    const els = new Set([...Object.keys(comp), ...Object.keys(it.comp)]);
    const parts = [];
    for (const e of els){ const h = comp[e]||0, w = it.comp[e]||0; if (h !== w) parts.push(`${e}: you have ${h}${w ? '' : ', but there is no ' + e + ' in this species'}`); }
    return bad('Check your atoms. ' + parts.join('; ') + '. (Try typing the formula first if you\'re unsure.)');
  }
  // 2. connectivity
  const adj = d.atoms.map(() => new Set()), bs = d.atoms.map(() => 0);
  d.bonds.forEach(b => { adj[b.a].add(b.b); adj[b.b].add(b.a); bs[b.a] += b.order; bs[b.b] += b.order; });
  const seen = new Set([0]), q = [0];
  while (q.length){ const x = q.pop(); for (const y of adj[x]) if (!seen.has(y)){ seen.add(y); q.push(y); } }
  if (seen.size < n){ const un = [...Array(n).keys()].filter(i => !seen.has(i)); return bad('Every atom has to be bonded into one structure. Some atoms aren\'t connected.', un); }
  const hBad = [...Array(n).keys()].filter(i => d.atoms[i].el === 'H' && (adj[i].size !== 1 || bs[i] !== 1));
  if (hBad.length) return bad('H forms only one single bond.', hBad);
  // 3. skeleton
  const kAt = it.skel.atoms, kAdj = kAt.map(() => new Set());
  it.skel.bonds.forEach(([x,y]) => { kAdj[x].add(y); kAdj[y].add(x); });
  const map = iso(d.atoms.map(a => a.el), adj, kAt.map(a => a.el), kAdj);
  if (!map){
    const keyHonO = kAt.every((a,i) => a.el !== 'H' || [...kAdj[i]].every(j => kAt[j].el === 'O'));
    const hWrong = [...Array(n).keys()].filter(i => d.atoms[i].el === 'H' && [...adj[i]].some(j => d.atoms[j].el !== 'O'));
    if (keyHonO && hWrong.length) return bad('In oxoacids and hydrogen anions, each H is bonded to an O atom, not to the central atom.', hWrong);
    const keyOO = it.skel.bonds.some(([x,y]) => kAt[x].el === 'O' && kAt[y].el === 'O');
    const oo = d.bonds.filter(b => d.atoms[b.a].el === 'O' && d.atoms[b.b].el === 'O');
    if (!keyOO && oo.length) return bad('The O atoms in this species are not bonded to each other. Bond each O to the central atom.', [...new Set(oo.flatMap(b => [b.a,b.b]))]);
    let ci = 0; kAdj.forEach((s,i) => { if (s.size > kAdj[ci].size) ci = i; });
    const cEl = kAt[ci].el, cDeg = kAdj[ci].size;
    const dc = [...Array(n).keys()].filter(i => d.atoms[i].el === cEl);
    if (!dc.some(i => adj[i].size === cDeg)) return bad(`Check the skeleton. In this species a ${cEl} atom is bonded to ${cDeg} other atoms.`, dc);
    return bad('Check the skeleton: which atom is bonded to which? The atom that appears once (not O or H) is usually in the center.');
  }
  // 4. overall charge
  if (d.charge !== it.charge){
    if (it.charge === 0) return bad('Acids are neutral molecules. Set the overall charge to 0 (no brackets).');
    if (d.charge === 0) return bad('This is an ion. Set its overall charge (the box under the drawing) so it is drawn in brackets with its charge.');
    return bad(`The overall charge isn't right for ${it.names[0]}.`);
  }
  // 5. electrons
  const tot = 2*d.bonds.reduce((s,b) => s + b.order, 0) + 2*d.atoms.reduce((s,a) => s + a.lp, 0);
  if (tot !== it.N){
    const cnt = it.comp, parts = Object.keys(cnt).map(e => (cnt[e] > 1 ? cnt[e] + '×' : '') + e + '(' + VAL[e] + ')');
    const qtxt = it.charge ? ` + ${-it.charge} for the ${chargeText(it.charge)} charge` : '';
    return bad(`Your structure shows ${tot} valence electrons; it needs ${it.N}: ${parts.join(' + ')}${qtxt} = ${it.N}. (Each bond line = 2 e⁻, each lone pair = 2 e⁻.)`);
  }
  // 6. octets
  const eAround = i => 2*bs[i] + 2*d.atoms[i].lp;
  const octBad = [...Array(n).keys()].filter(i => eAround(i) !== (d.atoms[i].el === 'H' ? 2 : 8));
  if (octBad.length){
    const i0 = octBad[0], e0 = eAround(i0);
    let m = `Every atom except H needs an octet (8 e⁻). ${cap(atomDesc(d, i0))} has ${e0}.`;
    if (octBad.some(i => eAround(i) > 8)) m += ' This tool uses all-octet structures: no expanded octets. Use single bonds and formal charges instead of extra double bonds.';
    else m += ' Move lone pairs into bonds, or add lone pairs, until each atom has 8.';
    return bad(m, octBad);
  }
  // 7. formal charges
  const fc = d.atoms.map((a,i) => VAL[a.el] - 2*a.lp - bs[i]);
  const fcBad = [...Array(n).keys()].filter(i => fc[i] !== d.atoms[i].fc);
  if (fcBad.length){
    const i0 = fcBad[0];
    return bad(`Check the formal charges. ${cap(atomDesc(d, i0))} should be ${fc[i0] ? fcText(fc[i0]) : '0'} (you marked ${d.atoms[i0].fc ? fcText(d.atoms[i0].fc) : 'none'}). Formal charge = valence e⁻ − lone-pair e⁻ − bonds.` + (fcBad.length > 1 ? ` ${fcBad.length} atoms need fixing (circled).` : ''), fcBad);
  }
  // 8. best form
  const sumAbs = fc.reduce((s,f) => s + Math.abs(f), 0);
  if (sumAbs > it.minAbs) return bad(`This is a valid all-octet Lewis structure, but not the best one: another arrangement of bonds gives smaller formal charges (yours add up to ${sumAbs} in size; the best is ${it.minAbs}).`);
  const en = fc.reduce((s,f,i) => s + f*EN[d.atoms[i].el], 0);
  let msg = 'Correct Lewis structure!';
  if (en > it.bestEn + 1e-9) msg += ' (It\'s acceptable, but the preferred form puts the negative formal charge on the more electronegative atom.)';
  else if (it.nBest > 1) msg += ` ${cap(it.names[0])} has ${it.nBest} equivalent resonance structures; you drew one of them.`;
  return {ok:true, msg, bad:[]};
}

if (typeof module !== 'undefined') module.exports = {VAL, EN, fmtU, chargeText, dispFormula, parseFormula, readsAs, prepareSpecies, linkSpecies, solveOctet, diagnoseName, gradeFormula, gradeLewis, explain, ruleFor, normName, compEq, compText, buildSkel, ateOf, RULES};
