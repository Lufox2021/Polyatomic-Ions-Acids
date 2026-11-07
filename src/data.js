// ===== Built-in species =====
const RAW_SPECIES = [];
function A(id, formula, names, layout, extra){ RAW_SPECIES.push(Object.assign({id, kind:'anion', formula, names, layout}, extra||{})); }
function BA(id, formula, names, layout, anionFormula, anionName, root, extra){ RAW_SPECIES.push(Object.assign({id, kind:'binary', formula, names, layout, anionFormula, anionName, root}, extra||{})); }
function OA(id, formula, names, layout, extra){ RAW_SPECIES.push(Object.assign({id, kind:'oxoacid', formula, names, layout}, extra||{})); }
const L = (center, ligs, ang) => ({center, ligs, ang});

// --- polyatomic anions: oxoanion series ---
A('NO3', 'NO3^-',  ['nitrate'],  L('N',['O','O','O']));
A('NO2', 'NO2^-',  ['nitrite'],  L('N',['O','O']));
A('SO4', 'SO4^2-', ['sulfate'],  L('S',['O','O','O','O']));
A('SO3', 'SO3^2-', ['sulfite'],  L('S',['O','O','O']));
A('PO4', 'PO4^3-', ['phosphate'],L('P',['O','O','O','O']));
A('CO3', 'CO3^2-', ['carbonate'],L('C',['O','O','O']));
for (const [X, r] of [['Cl','chlor'],['Br','brom'],['I','iod']]){
  A(X+'O4', X+'O4^-', ['per'+r+'ate'], L(X,['O','O','O','O']));
  A(X+'O3', X+'O3^-', [r+'ate'],       L(X,['O','O','O']));
  A(X+'O2', X+'O2^-', [r+'ite'],       L(X,['O','O']));
  A(X+'O',  X+'O^-',  ['hypo'+r+'ite'],L(X,['O'],[0]));
}
// --- hydrogen anions ---
A('HCO3', 'HCO3^-',  ['hydrogen carbonate','bicarbonate'], L('C',['O','O','OH']));
A('HSO4', 'HSO4^-',  ['hydrogen sulfate','bisulfate'],     L('S',['O','O','O','OH']));
A('HSO3', 'HSO3^-',  ['hydrogen sulfite','bisulfite'],     L('S',['O','O','OH']));
A('HPO4', 'HPO4^2-', ['hydrogen phosphate','monohydrogen phosphate'], L('P',['O','O','O','OH']));
A('H2PO4','H2PO4^-', ['dihydrogen phosphate'],             L('P',['O','O','OH','OH']));
// --- other common anions ---
A('OH',   'OH^-',     ['hydroxide'],  {special:'OH'},  {group:'ide', family:'OH'});
A('CN',   'CN^-',     ['cyanide'],    {special:'CN'},  {group:'ide', family:'CN'});
A('O2',   'O2^2-',    ['peroxide'],   {special:'O2'},  {group:'ide', family:'O2'});
A('SCN',  'SCN^-',    ['thiocyanate'],{special:'SCN'}, {group:'thio', family:'SCN',
   explain:'SCN⁻ is cyanate (OCN⁻) with S in place of O: <b>thiocyanate</b>.'});
A('S2O3', 'S2O3^2-',  ['thiosulfate'],L('S',['S','O','O','O']), {group:'thio', family:'S2O3',
   explain:'S₂O₃²⁻ is sulfate (SO₄²⁻) with one O replaced by S: <b>thiosulfate</b>.'});
A('C2H3O2','C2H3O2^-',['acetate','ethanoate'], {special:'acetate'}, {group:'organic', family:'acetate',
   explain:'C₂H₃O₂⁻ (CH₃COO⁻) is <b>acetate</b>, one of the carbon-containing anions to memorize. (IUPAC name: ethanoate.)'});
A('C2O4', 'C2O4^2-',  ['oxalate'],    {special:'oxalate'}, {group:'organic', family:'oxalate',
   explain:'C₂O₄²⁻ is <b>oxalate</b>, one of the carbon-containing anions to memorize.'});
A('CrO4', 'CrO4^2-',  ['chromate'],   L('Cr',['O','O','O','O']), {group:'metal',
   explain:'CrO₄²⁻ is <b>chromate</b>, a transition-metal oxoanion to memorize.'});
A('Cr2O7','Cr2O7^2-', ['dichromate'], {special:'Cr2O7'}, {group:'metal', family:'Cr2O7',
   explain:'Cr₂O₇²⁻ is <b>dichromate</b>: two chromate units sharing one O.'});
A('MnO4', 'MnO4^-',   ['permanganate'],L('Mn',['O','O','O','O']), {group:'metal',
   explain:'MnO₄⁻ is <b>permanganate</b>, a transition-metal oxoanion to memorize.'});

// --- binary acids (named in water) ---
BA('HF',  'HF',  ['hydrofluoric acid'], L('F',['H']),  'F^-',  'fluoride', 'fluor');
BA('HCl', 'HCl', ['hydrochloric acid'], L('Cl',['H']), 'Cl^-', 'chloride', 'chlor');
BA('HBr', 'HBr', ['hydrobromic acid'],  L('Br',['H']), 'Br^-', 'bromide',  'brom');
BA('HI',  'HI',  ['hydroiodic acid','hydriodic acid'], L('I',['H']), 'I^-', 'iodide', 'iod');
BA('H2S', 'H2S', ['hydrosulfuric acid'],L('S',['H','H']), 'S^2-', 'sulfide', 'sulfur');
BA('HCN', 'HCN', ['hydrocyanic acid'],  {special:'HCN'}, 'CN^-', 'cyanide', 'cyan', {group:'hcn', family:'CN'});

// --- oxoacids ---
OA('HNO3',  'HNO3',  ['nitric acid'],     L('N',['O','O','OH']));
OA('HNO2',  'HNO2',  ['nitrous acid'],    L('N',['OH','O']));
OA('H2SO4', 'H2SO4', ['sulfuric acid'],   L('S',['O','OH','O','OH']));
OA('H2SO3', 'H2SO3', ['sulfurous acid'],  L('S',['O','OH','OH']));
OA('H3PO4', 'H3PO4', ['phosphoric acid'], L('P',['O','OH','OH','OH']));
OA('H2CO3', 'H2CO3', ['carbonic acid'],   L('C',['O','OH','OH']));
for (const [X, r] of [['Cl','chlor'],['Br','brom'],['I','iod']]){
  OA('H'+X+'O4', 'H'+X+'O4', ['per'+r+'ic acid'],  L(X,['O','O','O','OH']));
  OA('H'+X+'O3', 'H'+X+'O3', [r+'ic acid'],        L(X,['O','O','OH']));
  OA('H'+X+'O2', 'H'+X+'O2', [r+'ous acid'],       L(X,['OH','O']));
  OA('H'+X+'O',  'H'+X+'O',  ['hypo'+r+'ous acid'],L(X,['OH']));
}
OA('H2CrO4',  'H2CrO4',  ['chromic acid'],      L('Cr',['O','OH','O','OH']));
OA('HMnO4',   'HMnO4',   ['permanganic acid'],  L('Mn',['O','O','O','OH']));
OA('HC2H3O2', 'HC2H3O2', ['acetic acid','ethanoic acid'], {special:'acetic'}, {family:'acetate'});
OA('H2C2O4',  'H2C2O4',  ['oxalic acid'],       {special:'oxalic'}, {family:'oxalate'});
