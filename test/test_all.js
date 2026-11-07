const fs = require('fs'), path = require('path');
const src = fs.readFileSync(path.join(__dirname,'../src/engine.js'),'utf8') + '\n' + fs.readFileSync(path.join(__dirname,'../src/data.js'),'utf8');
const ctx = {}; new Function('module', src + '\nmodule.exports.RAW = RAW_SPECIES; Object.assign(module.exports, {prepareSpecies, linkSpecies, diagnoseName, gradeFormula, gradeLewis, explain, fmtU, dispFormula, parseFormula});')(ctx.m = {exports:{}});
const E = ctx.m.exports;
let fails = 0; const fail = (...a) => { fails++; console.log('FAIL', ...a); };
const list = E.linkSpecies(E.RAW.map(r => { try { return E.prepareSpecies(r); } catch(e){ fail(r.id, e.message); } }).filter(Boolean));
console.log('species:', list.length, 'anions', list.filter(x=>x.kind==='anion').length, 'binary', list.filter(x=>x.kind==='binary').length, 'oxoacids', list.filter(x=>x.kind==='oxoacid').length);
const ids = new Set(); list.forEach(x => { if (ids.has(x.id)) fail('dup', x.id); ids.add(x.id); });
for (const it of list){
  // key structure as drawing
  const d = {atoms: it.skel.atoms.map((a,i) => ({el:a.el, x:a.x, y:a.y, lp: it.key.lp[i], fc: it.key.fc[i]})),
             bonds: it.skel.bonds.map(([a,b],i) => ({a, b, order: it.key.ord[i]})), charge: it.charge};
  const g = E.gradeLewis(it, d, list); if (!g.ok) fail(it.id, 'key lewis', g.msg);
  // shuffled atom order
  const perm = d.atoms.map((_,i)=>i).sort(()=>Math.random()-0.5), inv = []; perm.forEach((p,i)=>inv[p]=i);
  const d2 = {atoms: perm.map(p => d.atoms[p]), bonds: d.bonds.map(b => ({a:inv[b.a], b:inv[b.b], order:b.order})), charge: d.charge};
  if (!E.gradeLewis(it, d2, list).ok) fail(it.id, 'shuffled lewis');
  for (const n of it.names){ const r = E.diagnoseName(it, n, list); if (!r.ok) fail(it.id, 'name', n, r.msg); }
  const r2 = E.diagnoseName(it, it.names[0].toUpperCase() + ' ion', list); if (!r2.ok) fail(it.id, 'name caps');
  const f = E.gradeFormula(it, it.formula, list); if (!f.ok) fail(it.id, 'formula', f.msg);
  const fu = E.gradeFormula(it, E.dispFormula(it), list); if (!fu.ok) fail(it.id, 'unicode formula', fu.msg);
  // other species' names should all be rejected
  for (const o of list) if (o !== it) for (const n of o.names) if (E.diagnoseName(it, n, list).ok) fail(it.id, 'accepted other name', n);
  const ex = E.explain(it, list); if (!ex || /undefined/.test(ex)) fail(it.id, 'explain', ex);
  console.log(it.id.padEnd(8), E.dispFormula(it).padEnd(10), it.names[0].padEnd(24), 'N=' + it.N, 'FC', it.key.fc.join(','), 'res', it.nBest, it.parent ? '<' + it.parent.id : '');
}
const get = id => list.find(x => x.id === id);
const show = (id, ans) => console.log(`  [${id}] "${ans}" →`, E.diagnoseName(get(id), ans, list).msg);
show('ClO3','chlorite'); show('ClO3','chloric acid'); show('HClO2','chloric acid'); show('HClO2','chlorite acid');
show('HCl','hydrogen chloride'); show('HCl','chloric acid'); show('HCl','chloride'); show('H2SO4','hydrosulfuric acid');
show('H2SO4','sulfic acid'); show('HSO4','sulfate'); show('H2PO4','hydrogen phosphate'); show('IO4','iodate'); show('BrO','bromite');
show('OH','hydroxate'); show('NO3','nitric acid'); show('HNO3','nitrate');
const sf = (id, t) => console.log(`  [${id}] "${t}" →`, E.gradeFormula(get(id), t, list).msg);
sf('SO4','SO42-'); sf('SO4','SO4'); sf('SO4','SO3^2-'); sf('ClO2','ClO3-'); sf('H2SO4','H2SO4^2-'); sf('H3PO4','H2PO4'); sf('NO3','no3-'); sf('NO3','NO3 -'); sf('HC2H3O2','CH3COOH'); sf('C2H3O2','CH3COO^-');
// expanded-octet sulfate should be rejected with the octet message
{ const it = get('SO4'); const d = {atoms:[{el:'S',lp:0,fc:0}, {el:'O',lp:2,fc:0},{el:'O',lp:2,fc:0},{el:'O',lp:3,fc:-1},{el:'O',lp:3,fc:-1}], bonds:[{a:0,b:1,order:2},{a:0,b:2,order:2},{a:0,b:3,order:1},{a:0,b:4,order:1}], charge:-2};
  const g = E.gradeLewis(it, d, list); console.log('  expanded SO4:', g.msg); if (g.ok) fail('expanded accepted'); }
// SCN- S(-1) form
{ const it = get('SCN'); const d = {atoms:[{el:'S',lp:3,fc:-1},{el:'C',lp:0,fc:0},{el:'N',lp:1,fc:0}], bonds:[{a:0,b:1,order:1},{a:1,b:2,order:3}], charge:-1};
  console.log('  SCN S- form:', E.gradeLewis(it, d, list).msg); }
// H on central atom
{ const it = get('HNO2'); const d = {atoms:[{el:'H',lp:0,fc:0},{el:'N',lp:0,fc:0},{el:'O',lp:2,fc:0},{el:'O',lp:2,fc:0}], bonds:[{a:0,b:1,order:1},{a:1,b:2,order:1},{a:1,b:3,order:1}], charge:0};
  console.log('  HNO2 H on N:', E.gradeLewis(it, d, list).msg); }
console.log(fails ? fails + ' FAILURES' : 'ALL PASS');
