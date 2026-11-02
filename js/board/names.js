/* names.js (run 32): chemical names and condensed formulas -> SMILES, for the reagent box and the Board's structure box.
   OCName.parse(text) -> canonical SMILES or null;  OCName.read(text) -> { smi, how } or null.
   Reads: SMILES; ~450 trivial names; functional-class names (isobutyl chloride, ethyl methyl ketone, ethylmagnesium
   bromide, sodium ethoxide, phenyl acetate, diethyl malonate); IUPAC substitutive names for chains, rings and benzene
   derivatives (2-bromo-2-methylpropane, 4-methylpent-3-en-2-one, 2,4,6-trinitrophenol, p-nitrotoluene, ethyl propanoate);
   condensed formulas ((CH3)2CHCH2Cl, CH3COOC2H5, C6H5CH2MgBr, p-CH3C6H4NO2, CH2(COOEt)2, (CH3CO)2O, C6H5N2+Cl-). */
(function (root) {
'use strict';
const RDK = () => (root.Chem && root.Chem.RD) || root.RDKitModule || null;
function canon(smi) { const RD = RDK(); if (!RD || !smi) return null; let m = null; try { m = RD.get_mol(smi); if (!m) return null; if (m.is_valid && !m.is_valid()) { m.delete(); return null; } const s = m.get_smiles(); m.delete(); return s || null; } catch (e) { try { if (m) m.delete(); } catch (e2) { } return null; } }

/* ---------------- tiny molecular graph + SMILES writer ---------------- */
const VAL = { C: 4, N: 3, O: 2, S: 2, P: 3, B: 3, F: 1, Cl: 1, Br: 1, I: 1, Mg: 2, Li: 1, Na: 1, K: 1, Ca: 2, Ba: 2, Ag: 1, Cu: 1, Zn: 2, Hg: 2, Cd: 2, Al: 3, Si: 4, Se: 2 };
const METALS = new Set(['Mg', 'Li', 'Na', 'K', 'Ca', 'Ba', 'Ag', 'Cu', 'Zn', 'Hg', 'Cd', 'Al']);
const IONIC = new Set(['Na', 'K', 'Li', 'Ca', 'Ba', 'Ag']); // bonds from O/N/S to these are drawn as ions
const ORG = new Set(['B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']);
const G = () => ({ a: [], b: [] });
function addA(g, el, o) { o = o || {}; g.a.push({ el, chg: o.chg || 0, h: o.h != null ? o.h : null }); return g.a.length - 1; }
function addB(g, i, j, o) { g.b.push([i, j, o || 1]); }
function bond(g, i, j) { return g.b.find(b => (b[0] === i && b[1] === j) || (b[0] === j && b[1] === i)); }
function bsum(g, i) { let s = 0; for (const b of g.b) if (b[0] === i || b[1] === i) s += b[2]; return s; }
function cval(el, chg) { if (!chg) return VAL[el] != null ? VAL[el] : 0; if (METALS.has(el)) return 0; if (el === 'N' || el === 'P') return chg > 0 ? 4 : 2; if (el === 'O' || el === 'S') return chg > 0 ? 3 : 1; if (el === 'C') return 3; if (el === 'B') return chg < 0 ? 4 : 2; return 0; }
function tok(g, i) { const a = g.a[i]; if (!a.chg && ORG.has(a.el)) return a.el; const h = METALS.has(a.el) ? 0 : Math.max(0, cval(a.el, a.chg) - bsum(g, i)); let s = '[' + a.el + (h ? 'H' + (h > 1 ? h : '') : ''); if (a.chg) s += (a.chg > 0 ? '+' : '-') + (Math.abs(a.chg) > 1 ? Math.abs(a.chg) : ''); return s + ']'; }
function toSmi(g) {
  const n = g.a.length, adj = Array.from({ length: n }, () => []); g.b.forEach(([i, j, o], k) => { adj[i].push([j, o, k]); adj[j].push([i, o, k]); });
  const seen = new Array(n).fill(false), used = new Set(), parts = []; const BS = o => o === 2 ? '=' : o === 3 ? '#' : '';
  for (let s0 = 0; s0 < n; s0++) { if (seen[s0]) continue; const ch = new Map(), cl = new Map();
    const dfs = x => { seen[x] = true; ch.set(x, []); if (!cl.has(x)) cl.set(x, []); for (const [y, o, k] of adj[x]) { if (used.has(k)) continue; used.add(k); if (!seen[y]) { ch.get(x).push([y, o]); dfs(y); } else { cl.get(x).push([k, o]); if (!cl.has(y)) cl.set(y, []); cl.get(y).push([k, o]); } } };
    dfs(s0); const num = new Map(); let next = 1; const free = [];
    const emit = x => { let s = tok(g, x); for (const [k, o] of cl.get(x)) { if (num.has(k)) { const d = num.get(k); s += BS(o) + (d < 10 ? d : '%' + d); free.push(d); num.delete(k); } else { const d = free.length ? free.sort((p, q) => p - q).shift() : next++; num.set(k, d); s += BS(o) + (d < 10 ? d : '%' + d); } } const c = ch.get(x); c.forEach(([y, o], i) => { const sub = BS(o) + emit(y); s += i < c.length - 1 ? '(' + sub + ')' : sub; }); return s; };
    parts.push(emit(s0)); }
  return parts.join('.');
}
/* parse a small Kekulé SMILES fragment into g; its first atom is bonded to `at` (order ord). Returns the new atom ids */
function frag(g, smi, at, ord) {
  let i = 0, prev = null, bo = 1; const st = [], rings = {}, ids = [];
  while (i < smi.length) { const c = smi[i];
    if (c === '(') { st.push(prev); i++; continue; } if (c === ')') { prev = st.pop(); i++; continue; }
    if (c === '=') { bo = 2; i++; continue; } if (c === '#') { bo = 3; i++; continue; } if (c === '-' || c === '/' || c === '\\') { i++; continue; }
    if (c === '.') { prev = null; at = null; i++; continue; }
    if (/\d/.test(c)) { if (rings[c]) { addB(g, prev, rings[c][0], Math.max(bo, rings[c][1])); delete rings[c]; } else rings[c] = [prev, bo]; bo = 1; i++; continue; }
    let el, chg = 0;
    if (c === '[') { const j = smi.indexOf(']', i); const t = smi.slice(i + 1, j); const mm = t.match(/^([A-Z][a-z]?)(?:H\d*)?([+-]+\d*)?$/); el = mm[1]; if (mm[2]) { const sg = mm[2][0] === '+' ? 1 : -1; const d = mm[2].replace(/[+-]/g, ''); chg = sg * (d ? +d : mm[2].length); } i = j + 1; }
    else if (smi.startsWith('Cl', i)) { el = 'Cl'; i += 2; } else if (smi.startsWith('Br', i)) { el = 'Br'; i += 2; } else { el = c; i++; }
    const id = addA(g, el, { chg }); ids.push(id);
    if (prev != null) addB(g, prev, id, bo); else if (at != null && ids.length === 1) addB(g, at, id, ord || 1);
    prev = id; bo = 1; }
  return ids;
}
/* attach a prefix template to atom `at`: templates starting with '=' make a double bond */
function attach(g, at, tpl) { if (tpl[0] === '=') return frag(g, tpl.slice(1), at, 2); return frag(g, tpl, at, 1); }
const outSmi = g => canon(toSmi(g));

/* ---------------- groups ---------------- */
const PH = 'C1=CC=CC=C1';
const ALK = { methyl: 'C', ethyl: 'CC', propyl: 'CCC', npropyl: 'CCC', isopropyl: 'C(C)C', ipropyl: 'C(C)C', butyl: 'CCCC', nbutyl: 'CCCC', isobutyl: 'CC(C)C', ibutyl: 'CC(C)C', secbutyl: 'C(C)CC', sbutyl: 'C(C)CC',
  tertbutyl: 'C(C)(C)C', tbutyl: 'C(C)(C)C', pentyl: 'CCCCC', npentyl: 'CCCCC', amyl: 'CCCCC', namyl: 'CCCCC', isopentyl: 'CCC(C)C', isoamyl: 'CCC(C)C', neopentyl: 'CC(C)(C)C', tertpentyl: 'C(C)(C)CC', tertamyl: 'C(C)(C)CC', tamyl: 'C(C)(C)CC',
  hexyl: 'CCCCCC', nhexyl: 'CCCCCC', heptyl: 'CCCCCCC', octyl: 'CCCCCCCC', nonyl: 'CCCCCCCCC', decyl: 'CCCCCCCCCC', cyclopropyl: 'C1CC1', cyclobutyl: 'C1CCC1', cyclopentyl: 'C1CCCC1', cyclohexyl: 'C1CCCCC1', cycloheptyl: 'C1CCCCCC1',
  vinyl: 'C=C', ethenyl: 'C=C', allyl: 'CC=C', propenyl: 'C=CC', isopropenyl: 'C(C)=C', propargyl: 'CC#C', ethynyl: 'C#C', crotyl: 'CC=CC', benzyl: 'C' + PH, phenyl: PH, tolyl: 'C1=CC=C(C)C=C1', ptolyl: 'C1=CC=C(C)C=C1', otolyl: 'C1=C(C)C=CC=C1', mtolyl: 'C1=CC(C)=CC=C1',
  phenethyl: 'CC' + PH, '2phenylethyl': 'CC' + PH, '1phenylethyl': 'C(C)' + PH, naphthyl: 'C1=CC=CC2=CC=CC=C12', '1naphthyl': 'C1=CC=CC2=CC=CC=C12', alphanaphthyl: 'C1=CC=CC2=CC=CC=C12', '2naphthyl': 'C1=CC2=CC=CC=C2C=C1', betanaphthyl: 'C1=CC2=CC=CC=C2C=C1',
  benzhydryl: 'C(' + PH + ')' + PH, diphenylmethyl: 'C(' + PH + ')' + PH, trityl: 'C(' + PH + ')(' + PH + ')' + PH, triphenylmethyl: 'C(' + PH + ')(' + PH + ')' + PH, pchlorophenyl: 'C1=CC=C(Cl)C=C1', pbromophenyl: 'C1=CC=C(Br)C=C1', pnitrophenyl: 'C1=CC=C([N+](=O)[O-])C=C1', pmethoxyphenyl: 'C1=CC=C(OC)C=C1', pmethylphenyl: 'C1=CC=C(C)C=C1' };
const ALKRE = Object.keys(ALK).sort((a, b) => b.length - a.length).join('|');
const ACYL = { formyl: 'C=O', acetyl: 'C(C)=O', propionyl: 'C(=O)CC', propanoyl: 'C(=O)CC', butyryl: 'C(=O)CCC', butanoyl: 'C(=O)CCC', isobutyryl: 'C(=O)C(C)C', valeryl: 'C(=O)CCCC', pentanoyl: 'C(=O)CCCC', pivaloyl: 'C(=O)C(C)(C)C', benzoyl: 'C(=O)' + PH,
  phenylacetyl: 'C(=O)C' + PH, chloroacetyl: 'C(=O)CCl', acryloyl: 'C(=O)C=C', crotonyl: 'C(=O)C=CC', cinnamoyl: 'C(=O)C=C' + PH, ptoluoyl: 'C(=O)C1=CC=C(C)C=C1', pnitrobenzoyl: 'C(=O)C1=CC=C([N+](=O)[O-])C=C1', ethanoyl: 'C(C)=O', methanoyl: 'C=O' };
const ACYLDI = { oxalyl: 'C(=O)C(=O)', malonyl: 'C(=O)CC(=O)', succinyl: 'C(=O)CCC(=O)', glutaryl: 'C(=O)CCCC(=O)', adipoyl: 'C(=O)CCCCC(=O)', phthaloyl: 'C(=O)C1=CC=CC=C1C(=O)', terephthaloyl: 'C(=O)C1=CC=C(C=C1)C(=O)' };
const ANION = { formate: 'C=O', methanoate: 'C=O', acetate: 'C(C)=O', ethanoate: 'C(C)=O', propionate: 'C(=O)CC', propanoate: 'C(=O)CC', butyrate: 'C(=O)CCC', butanoate: 'C(=O)CCC', isobutyrate: 'C(=O)C(C)C', valerate: 'C(=O)CCCC', pentanoate: 'C(=O)CCCC',
  benzoate: 'C(=O)' + PH, salicylate: 'C(=O)C1=CC=CC=C1O', cinnamate: 'C(=O)C=C' + PH, acrylate: 'C(=O)C=C', propenoate: 'C(=O)C=C', methacrylate: 'C(=O)C(C)=C', crotonate: 'C(=O)C=CC', chloroacetate: 'C(=O)CCl', bromoacetate: 'C(=O)CBr', iodoacetate: 'C(=O)CI',
  cyanoacetate: 'C(=O)CC#N', acetoacetate: 'C(=O)CC(C)=O', phenylacetate: 'C(=O)C' + PH, lactate: 'C(=O)C(C)O', pyruvate: 'C(=O)C(C)=O', trichloroacetate: 'C(=O)C(Cl)(Cl)Cl', trifluoroacetate: 'C(=O)C(F)(F)F', glycinate: 'C(=O)CN', anthranilate: 'C(=O)C1=CC=CC=C1N' };
const ANIONDI = { oxalate: 'C(=O)C(=O)', malonate: 'C(=O)CC(=O)', succinate: 'C(=O)CCC(=O)', glutarate: 'C(=O)CCCC(=O)', adipate: 'C(=O)CCCCC(=O)', phthalate: 'C(=O)C1=CC=CC=C1C(=O)', terephthalate: 'C(=O)C1=CC=C(C=C1)C(=O)', carbonate: 'C(=O)', maleate: 'C(=O)C=CC(=O)', fumarate: 'C(=O)C=CC(=O)', sulfate: 'S(=O)(=O)' };
const ALKOXIDE = { methoxide: 'C', ethoxide: 'CC', propoxide: 'CCC', npropoxide: 'CCC', isopropoxide: 'C(C)C', butoxide: 'CCCC', nbutoxide: 'CCCC', tertbutoxide: 'C(C)(C)C', tbutoxide: 'C(C)(C)C', isobutoxide: 'CC(C)C', secbutoxide: 'C(C)CC', phenoxide: PH, benzyloxide: 'C' + PH, tertpentoxide: 'C(C)(C)CC', tertamyloxide: 'C(C)(C)CC' };
const CATION = { sodium: ['[Na+]', 1], potassium: ['[K+]', 1], lithium: ['[Li+]', 1], silver: ['[Ag+]', 1], ammonium: ['[NH4+]', 1], calcium: ['[Ca+2]', 2], barium: ['[Ba+2]', 2], magnesium: ['[Mg+2]', 2] };

/* ---------------- trivial names ---------------- */
const NO2 = '[N+](=O)[O-]';
const TRIV_SRC = `methane=C;ethane=CC;propane=CCC;butane|nbutane=CCCC;isobutane=CC(C)C;pentane|npentane=CCCCC;isopentane=CCC(C)C;neopentane=CC(C)(C)C;hexane|nhexane=CCCCCC;heptane|nheptane=CCCCCCC;octane=CCCCCCCC;isooctane=CC(C)CC(C)(C)C;
ethylene|ethene=C=C;propylene|propene=C=CC;butene|1butene|butylene=C=CCC;2butene=CC=CC;isobutylene|isobutene=C=C(C)C;acetylene|ethyne=C#C;propyne|methylacetylene=CC#C;1butyne|ethylacetylene=C#CCC;2butyne|dimethylacetylene=CC#CC;allene|propadiene=C=C=C;butadiene|13butadiene=C=CC=C;isoprene=C=CC(C)=C;chloroprene=C=CC(Cl)=C;vinylacetylene=C=CC#C;
cyclopropane=C1CC1;cyclobutane=C1CCC1;cyclopentane=C1CCCC1;cyclohexane=C1CCCCC1;cycloheptane=C1CCCCCC1;cyclooctane=C1CCCCCCC1;cyclopentene=C1=CCCC1;cyclohexene=C1=CCCCC1;cyclopentadiene=C1=CCC=C1;13cyclohexadiene=C1=CC=CCC1;14cyclohexadiene=C1=CCC=CC1;cyclooctatetraene|cot=C1=CC=CC=CC=C1;methylenecyclohexane=C=C1CCCCC1;decalin=C1CCC2CCCCC2C1;tetralin=c1ccc2c(c1)CCCC2;
benzene|benzol=c1ccccc1;toluene|toluol=Cc1ccccc1;oxylene=Cc1ccccc1C;mxylene=Cc1cccc(C)c1;pxylene|xylene=Cc1ccc(C)cc1;mesitylene=Cc1cc(C)cc(C)c1;ethylbenzene=CCc1ccccc1;cumene=CC(C)c1ccccc1;styrene=C=Cc1ccccc1;naphthalene=c1ccc2ccccc2c1;anthracene=c1ccc2cc3ccccc3cc2c1;phenanthrene=c1ccc2c(c1)ccc1ccccc12;biphenyl|diphenyl=c1ccc(-c2ccccc2)cc1;diphenylmethane=c1ccc(Cc2ccccc2)cc1;triphenylmethane=c1ccc(C(c2ccccc2)c2ccccc2)cc1;stilbene|transstilbene=C(=C/c1ccccc1)\\c1ccccc1;cisstilbene=C(=C\\c1ccccc1)\\c1ccccc1;phenylacetylene=C#Cc1ccccc1;
chloroform=ClC(Cl)Cl;bromoform=BrC(Br)Br;iodoform=IC(I)I;carbontetrachloride|tetrachloromethane=ClC(Cl)(Cl)Cl;carbontetrabromide=BrC(Br)(Br)Br;methylenechloride|methylenedichloride=ClCCl;methylenebromide=BrCBr;methyleneiodide=ICI;vinylchloride=C=CCl;vinylbromide=C=CBr;allylchloride=C=CCCl;allylbromide=C=CCBr;allyliodide=C=CCI;benzylchloride=ClCc1ccccc1;benzylbromide=BrCc1ccccc1;chlorobenzene=Clc1ccccc1;bromobenzene=Brc1ccccc1;iodobenzene=Ic1ccccc1;fluorobenzene=Fc1ccccc1;
benzalchloride|benzylidenechloride|benzylidenedichloride=ClC(Cl)c1ccccc1;benzotrichloride=ClC(Cl)(Cl)c1ccccc1;ethylenedichloride|ethylenechloride=ClCCCl;ethylenedibromide|ethylenebromide=BrCCBr;ethylidenechloride|ethylidenedichloride=CC(Cl)Cl;ethylidenebromide|ethylidenedibromide=CC(Br)Br;chloral=O=CC(Cl)(Cl)Cl;chloralhydrate=OC(O)C(Cl)(Cl)Cl;chloropicrin=[O-][N+](=O)C(Cl)(Cl)Cl;phosgene|carbonylchloride=O=C(Cl)Cl;freon|freon12=FC(F)(Cl)Cl;
gammexane|bhc|lindane|benzenehexachloride|hexachlorocyclohexane=ClC1C(Cl)C(Cl)C(Cl)C(Cl)C1Cl;ddt=ClC(Cl)(Cl)C(c1ccc(Cl)cc1)c1ccc(Cl)cc1;westron=ClC(Cl)C(Cl)Cl;westrosol|trichloroethylene|trichloroethene=ClC=C(Cl)Cl;tetrachloroethylene=ClC(Cl)=C(Cl)Cl;
methanol|methylalcohol|woodspirit|carbinol=CO;ethanol|ethylalcohol|grainalcohol=CCO;propanol|npropanol|1propanol|npropylalcohol|propylalcohol=CCCO;isopropanol|isopropylalcohol|2propanol|ipa=CC(C)O;butanol|nbutanol|1butanol|nbutylalcohol|butylalcohol=CCCCO;2butanol|secbutylalcohol|sbutylalcohol|secbutanol=CCC(C)O;isobutanol|isobutylalcohol=CC(C)CO;tertbutanol|tbutanol|tertbutylalcohol|tbutylalcohol|2methyl2propanol=CC(C)(C)O;
pentanol|amylalcohol|namylalcohol=CCCCCO;isoamylalcohol|isopentylalcohol=CC(C)CCO;neopentylalcohol=CC(C)(C)CO;tertamylalcohol|tertpentylalcohol=CCC(C)(C)O;benzylalcohol=OCc1ccccc1;allylalcohol=C=CCO;propargylalcohol=C#CCO;ethyleneglycol|glycol=OCCO;propyleneglycol=CC(O)CO;glycerol|glycerine|glycerin=OCC(O)CO;pinacol=CC(C)(O)C(C)(C)O;cyclohexanol=OC1CCCCC1;cyclopentanol=OC1CCCC1;benzhydrol=OC(c1ccccc1)c1ccccc1;triphenylcarbinol=OC(c1ccccc1)(c1ccccc1)c1ccccc1;
phenol|carbolicacid=Oc1ccccc1;catechol=Oc1ccccc1O;resorcinol=Oc1cccc(O)c1;hydroquinone|quinol=Oc1ccc(O)cc1;ocresol=Cc1ccccc1O;mcresol=Cc1cccc(O)c1;pcresol|cresol=Cc1ccc(O)cc1;picricacid=Oc1c(${NO2})cc(${NO2})cc1${NO2};onitrophenol=Oc1ccccc1${NO2};mnitrophenol=Oc1cccc(${NO2})c1;pnitrophenol=Oc1ccc(${NO2})cc1;
1naphthol|alphanaphthol|anaphthol=Oc1cccc2ccccc12;2naphthol|betanaphthol|bnaphthol=Oc1ccc2ccccc2c1;salicylalcohol|saligenin=OCc1ccccc1O;pyrogallol=Oc1cccc(O)c1O;phloroglucinol=Oc1cc(O)cc(O)c1;guaiacol=COc1ccccc1O;ochlorophenol=Oc1ccccc1Cl;pchlorophenol=Oc1ccc(Cl)cc1;paminophenol=Nc1ccc(O)cc1;oaminophenol=Nc1ccccc1O;
anisole=COc1ccccc1;phenetole=CCOc1ccccc1;diphenylether=c1ccc(Oc2ccccc2)cc1;ethylmethylether|methylethylether=CCOC;mtbe=COC(C)(C)C;tetrahydrofuran=C1CCOC1;oxirane|epoxyethane=C1CO1;propyleneoxide|methyloxirane|12epoxypropane=CC1CO1;styreneoxide=C1OC1c1ccccc1;cyclohexeneoxide=C1CCC2OC2C1;epichlorohydrin=ClCC1CO1;
formaldehyde|methanal|formalin=C=O;acetaldehyde|ethanal=CC=O;propanal|propionaldehyde|propaldehyde=CCC=O;butanal|butyraldehyde|nbutyraldehyde=CCCC=O;isobutyraldehyde=CC(C)C=O;pentanal|valeraldehyde=CCCCC=O;pivaldehyde|trimethylacetaldehyde=CC(C)(C)C=O;benzaldehyde=O=Cc1ccccc1;salicylaldehyde=O=Cc1ccccc1O;phydroxybenzaldehyde=O=Cc1ccc(O)cc1;
anisaldehyde|panisaldehyde=COc1ccc(C=O)cc1;vanillin=COc1cc(C=O)ccc1O;cinnamaldehyde=O=C/C=C/c1ccccc1;acrolein|acrylaldehyde=C=CC=O;crotonaldehyde=CC=CC=O;glyoxal=O=CC=O;methylglyoxal|pyruvaldehyde=CC(=O)C=O;furfural=O=Cc1ccco1;phenylacetaldehyde=O=CCc1ccccc1;terephthalaldehyde=O=Cc1ccc(C=O)cc1;phthalaldehyde=O=Cc1ccccc1C=O;succinaldehyde=O=CCCC=O;glutaraldehyde=O=CCCCC=O;adipaldehyde=O=CCCCCC=O;pnitrobenzaldehyde=O=Cc1ccc(${NO2})cc1;mnitrobenzaldehyde=O=Cc1cccc(${NO2})c1;pchlorobenzaldehyde=O=Cc1ccc(Cl)cc1;ptolualdehyde=Cc1ccc(C=O)cc1;
acetone|dimethylketone=CC(C)=O;butanone|ethylmethylketone|methylethylketone|mek=CCC(C)=O;diethylketone=CCC(=O)CC;pinacolone=CC(=O)C(C)(C)C;acetophenone|methylphenylketone|phenylmethylketone=CC(=O)c1ccccc1;propiophenone|ethylphenylketone=CCC(=O)c1ccccc1;benzophenone|diphenylketone=O=C(c1ccccc1)c1ccccc1;cyclohexanone=O=C1CCCCC1;cyclopentanone=O=C1CCCC1;cyclobutanone=O=C1CCC1;cycloheptanone=O=C1CCCCCC1;
cyclohexenone=O=C1C=CCCC1;mesityloxide=CC(=O)C=C(C)C;diacetonealcohol=CC(=O)CC(C)(C)O;phorone=CC(C)=CC(=O)C=C(C)C;benzil=O=C(C(=O)c1ccccc1)c1ccccc1;benzoin=O=C(c1ccccc1)C(O)c1ccccc1;acetylacetone=CC(=O)CC(C)=O;biacetyl|diacetyl=CC(=O)C(C)=O;methylvinylketone|mvk=C=CC(C)=O;benzalacetone|benzylideneacetone=CC(=O)/C=C/c1ccccc1;dibenzalacetone=O=C(/C=C/c1ccccc1)/C=C/c1ccccc1;
chalcone|benzalacetophenone=O=C(/C=C/c1ccccc1)c1ccccc1;ketene=C=C=O;benzoquinone|pbenzoquinone|quinone=O=C1C=CC(=O)C=C1;acetonylacetone=CC(=O)CCC(C)=O;2methylcyclohexanone=CC1CCCCC1=O;camphor=CC1(C)C2CCC1(C)C(=O)C2;
formicacid=O=CO;aceticacid|vinegar=CC(=O)O;propionicacid=CCC(=O)O;butyricacid=CCCC(=O)O;isobutyricacid=CC(C)C(=O)O;valericacid=CCCCC(=O)O;pivalicacid=CC(C)(C)C(=O)O;caproicacid=CCCCCC(=O)O;palmiticacid=CCCCCCCCCCCCCCCC(=O)O;stearicacid=CCCCCCCCCCCCCCCCCC(=O)O;oleicacid=CCCCCCCC/C=C\\CCCCCCCC(=O)O;
benzoicacid=O=C(O)c1ccccc1;phenylaceticacid=O=C(O)Cc1ccccc1;salicylicacid|ohydroxybenzoicacid=O=C(O)c1ccccc1O;phydroxybenzoicacid=O=C(O)c1ccc(O)cc1;aspirin|acetylsalicylicacid=CC(=O)Oc1ccccc1C(=O)O;phthalicacid=O=C(O)c1ccccc1C(=O)O;isophthalicacid=O=C(O)c1cccc(C(=O)O)c1;terephthalicacid=O=C(O)c1ccc(C(=O)O)cc1;
oxalicacid=O=C(O)C(=O)O;malonicacid=O=C(O)CC(=O)O;succinicacid=O=C(O)CCC(=O)O;glutaricacid=O=C(O)CCCC(=O)O;adipicacid=O=C(O)CCCCC(=O)O;maleicacid=O=C(O)/C=C\\C(=O)O;fumaricacid=O=C(O)/C=C/C(=O)O;acrylicacid=C=CC(=O)O;methacrylicacid=C=C(C)C(=O)O;crotonicacid=C/C=C/C(=O)O;cinnamicacid=O=C(O)/C=C/c1ccccc1;
lacticacid=CC(O)C(=O)O;glycolicacid=O=C(O)CO;pyruvicacid=CC(=O)C(=O)O;tartaricacid=O=C(O)C(O)C(O)C(=O)O;citricacid=O=C(O)CC(O)(CC(=O)O)C(=O)O;malicacid=O=C(O)CC(O)C(=O)O;levulinicacid=CC(=O)CCC(=O)O;acetoaceticacid=CC(=O)CC(=O)O;chloroaceticacid|monochloroaceticacid=O=C(O)CCl;dichloroaceticacid=O=C(O)C(Cl)Cl;trichloroaceticacid=O=C(O)C(Cl)(Cl)Cl;trifluoroaceticacid|tfa=O=C(O)C(F)(F)F;bromoaceticacid=O=C(O)CBr;
panisicacid=COc1ccc(C(=O)O)cc1;pnitrobenzoicacid=O=C(O)c1ccc(${NO2})cc1;mnitrobenzoicacid=O=C(O)c1cccc(${NO2})c1;onitrobenzoicacid=O=C(O)c1ccccc1${NO2};ptoluicacid=Cc1ccc(C(=O)O)cc1;anthranilicacid=Nc1ccccc1C(=O)O;paba|paminobenzoicacid=Nc1ccc(C(=O)O)cc1;pchlorobenzoicacid=O=C(O)c1ccc(Cl)cc1;
propionylchloride=CCC(=O)Cl;butyrylchloride=CCCC(=O)Cl;oxalylchloride=O=C(Cl)C(=O)Cl;propionicanhydride=CCC(=O)OC(=O)CC;benzoicanhydride=O=C(OC(=O)c1ccccc1)c1ccccc1;phthalicanhydride=O=C1OC(=O)c2ccccc21;trifluoroaceticanhydride=O=C(OC(=O)C(F)(F)F)C(F)(F)F;
acetamide=CC(N)=O;formamide=NC=O;propionamide=CCC(N)=O;butyramide=CCCC(N)=O;benzamide=NC(=O)c1ccccc1;urea|carbamide=NC(N)=O;thiourea=NC(N)=S;phthalimide=O=C1NC(=O)c2ccccc21;succinimide=O=C1CCC(=O)N1;acetanilide=CC(=O)Nc1ccccc1;benzanilide=O=C(Nc1ccccc1)c1ccccc1;nmethylacetamide=CNC(C)=O;caprolactam=O=C1CCCCCN1;biuret=NC(=O)NC(N)=O;phenylurea=NC(=O)Nc1ccccc1;
ethylacetate=CCOC(C)=O;methylacetate=COC(C)=O;ethylformate=CCOC=O;methylformate=COC=O;methylsalicylate|oilofwintergreen=COC(=O)c1ccccc1O;ethylacetoacetate|acetoaceticester|eaa|acetoaceticethylester=CCOC(=O)CC(C)=O;methylacetoacetate=COC(=O)CC(C)=O;malonicester=CCOC(=O)CC(=O)OCC;methylmethacrylate|mma=C=C(C)C(=O)OC;vinylacetate=C=COC(C)=O;isopropenylacetate=C=C(C)OC(C)=O;butyrolactone|gammabutyrolactone=O=C1CCCO1;phthalide=O=C1OCc2ccccc21;coumarin=O=c1ccc2ccccc2o1;
acetonitrile|methylcyanide=CC#N;propionitrile|ethylcyanide=CCC#N;butyronitrile=CCCC#N;benzonitrile|phenylcyanide=N#Cc1ccccc1;benzylcyanide|phenylacetonitrile=N#CCc1ccccc1;malononitrile=N#CCC#N;succinonitrile=N#CCCC#N;adiponitrile=N#CCCCC#N;methylcarbylamine=[C-]#[N+]C;ethylcarbylamine=[C-]#[N+]CC;phenylcarbylamine=[C-]#[N+]c1ccccc1;mic=CN=C=O;
aminomethane=CN;aminoethane=CCN;dimethylamine=CNC;diethylamine=CCNCC;trimethylamine=CN(C)C;triethylamine=CCN(CC)CC;aniline|aminobenzene|benzenamine|phenylamine=Nc1ccccc1;nmethylaniline=CNc1ccccc1;nndimethylaniline|dimethylaniline=CN(C)c1ccccc1;nethylaniline=CCNc1ccccc1;nndiethylaniline|diethylaniline=CCN(CC)c1ccccc1;diphenylamine=c1ccc(Nc2ccccc2)cc1;triphenylamine=c1ccc(N(c2ccccc2)c2ccccc2)cc1;
otoluidine=Cc1ccccc1N;mtoluidine=Cc1cccc(N)c1;ptoluidine|toluidine=Cc1ccc(N)cc1;panisidine=COc1ccc(N)cc1;ophenylenediamine=Nc1ccccc1N;mphenylenediamine=Nc1cccc(N)c1;pphenylenediamine=Nc1ccc(N)cc1;ethylenediamine=NCCN;hexamethylenediamine=NCCCCCCN;ethanolamine=NCCO;pnitroaniline=Nc1ccc(${NO2})cc1;onitroaniline=Nc1ccccc1${NO2};mnitroaniline=Nc1cccc(${NO2})c1;pchloroaniline=Nc1ccc(Cl)cc1;pbromoaniline=Nc1ccc(Br)cc1;
sulfanilicacid=Nc1ccc(S(=O)(=O)O)cc1;sulfanilamide=Nc1ccc(S(N)(=O)=O)cc1;aniliniumchloride|anilinehydrochloride=[NH3+]c1ccccc1.[Cl-];hexamine|urotropine|hexamethylenetetramine=C1N2CN3CN1CN(C2)C3;pyrrole=c1cc[nH]c1;furan=c1ccoc1;thiophene=c1ccsc1;indole=c1ccc2[nH]ccc2c1;quinoline=c1ccc2ncccc2c1;imidazole=c1c[nH]cn1;piperidine=C1CCNCC1;morpholine=C1COCCN1;
nitrobenzene|nitrobenzol=O=[N+]([O-])c1ccccc1;nitromethane=C${NO2};nitroethane=CC${NO2};mdinitrobenzene=O=[N+]([O-])c1cccc(${NO2})c1;odinitrobenzene=O=[N+]([O-])c1ccccc1${NO2};pdinitrobenzene=O=[N+]([O-])c1ccc(${NO2})cc1;tnb=O=[N+]([O-])c1cc(${NO2})cc(${NO2})c1;tnt|trinitrotoluene=Cc1c(${NO2})cc(${NO2})cc1${NO2};
onitrotoluene=Cc1ccccc1${NO2};mnitrotoluene=Cc1cccc(${NO2})c1;pnitrotoluene=Cc1ccc(${NO2})cc1;pchloronitrobenzene=O=[N+]([O-])c1ccc(Cl)cc1;ochloronitrobenzene=O=[N+]([O-])c1ccccc1Cl;pfluoronitrobenzene=O=[N+]([O-])c1ccc(F)cc1;nitrosobenzene=O=Nc1ccccc1;azobenzene=c1ccc(/N=N/c2ccccc2)cc1;hydrazobenzene=c1ccc(NNc2ccccc2)cc1;phenylhydroxylamine|nphenylhydroxylamine=ONc1ccccc1;
benzenediazoniumchloride|phenyldiazoniumchloride|diazoniumchloride|diazoniumsalt=N#[N+]c1ccccc1.[Cl-];benzenediazonium=N#[N+]c1ccccc1;ptoluenediazoniumchloride=Cc1ccc([N+]#N)cc1.[Cl-];benzenesulfonicacid=O=S(=O)(O)c1ccccc1;ptoluenesulfonicacid=Cc1ccc(S(=O)(=O)O)cc1;benzenesulfonylchloride=O=S(=O)(Cl)c1ccccc1;
thiophenol=Sc1ccccc1;ethanethiol|ethylmercaptan=CCS;methanethiol|methylmercaptan=CS;diethylsulfide=CCSCC;dimethylsulfide=CSC;
glucose|dglucose|dextrose=O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO;fructose|dfructose|laevulose|levulose|fruitsugar=OCC(=O)[C@@H](O)[C@H](O)[C@H](O)CO;sucrose|canesugar=OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O;glyceraldehyde=O=CC(O)CO;dihydroxyacetone=O=C(CO)CO;ribose=O=CC(O)C(O)C(O)CO;sorbitol|glucitol=OCC(O)C(O)C(O)C(O)CO;gluconicacid=O=C(O)C(O)C(O)C(O)C(O)CO;saccharicacid|glucaricacid=O=C(O)C(O)C(O)C(O)C(O)C(=O)O;
glycine|aminoaceticacid=NCC(=O)O;alanine=CC(N)C(=O)O;phenylalanine=NC(Cc1ccccc1)C(=O)O;valine=CC(C)C(N)C(=O)O;leucine=CC(C)CC(N)C(=O)O;serine=NC(CO)C(=O)O;cysteine=NC(CS)C(=O)O;asparticacid=NC(CC(=O)O)C(=O)O;glutamicacid=NC(CCC(=O)O)C(=O)O;lysine=NCCCCC(N)C(=O)O;tyrosine=NC(Cc1ccc(O)cc1)C(=O)O;
calciumcarbide=[C-]#[C-].[Ca+2];sodiumacetylide=[C-]#C.[Na+];dimethylsulfoxide=CS(C)=O`;
const TRIV = new Map();
for (const ent of TRIV_SRC.replace(/\n/g, '').split(';')) { if (!ent.trim()) continue; const k = ent.lastIndexOf('='); const names = ent.slice(0, ent.indexOf('=') > 0 ? ent.indexOf('=') : k), smi = ent.slice(ent.indexOf('=') + 1); for (const n of names.split('|')) TRIV.set(n.trim(), smi.trim()); }
const MF = { ch4: 'C', c2h6: 'CC', c2h4: 'C=C', c2h2: 'C#C', c3h8: 'CCC', c3h6: 'C=CC', c3h4: 'CC#C', c4h10: 'CCCC', c6h6: 'c1ccccc1', c6h12: 'C1CCCCC1', c6h10: 'C1=CCCCC1', c7h8: 'Cc1ccccc1', c6h12o6: 'O=CC(O)C(O)C(O)C(O)CO', c12h22o11: 'OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O', cac2: '[C-]#[C-].[Ca+2]' };

const norm0 = s => String(s).toLowerCase().replace(/α/g, 'alpha').replace(/β/g, 'beta').replace(/γ/g, 'gamma').replace(/sulph/g, 'sulf').replace(/[\s\-–—,'’`()\[\].]+/g, '');

/* ---------------- functional-class names ---------------- */
const RX = (re) => new RegExp('^' + re + '$');
const A = '(' + ALKRE + ')';
const FC = [
  [RX(A + '(chloride|bromide|iodide|fluoride)'), m => ({ chloride: 'Cl', bromide: 'Br', iodide: 'I', fluoride: 'F' })[m[2]] + ALK[m[1]]],
  [RX(A + 'alcohol'), m => 'O' + ALK[m[1]]],
  [RX('tri' + A + 'amine'), m => 'N(' + ALK[m[1]] + ')(' + ALK[m[1]] + ')' + ALK[m[1]]],
  [RX('di' + A + 'amine'), m => 'N(' + ALK[m[1]] + ')' + ALK[m[1]]],
  [RX(A + A + 'amine'), m => 'N(' + ALK[m[1]] + ')' + ALK[m[2]]],
  [RX(A + 'amine'), m => 'N' + ALK[m[1]]],
  [RX(A + 'cyanide'), m => 'N#C' + ALK[m[1]]],
  [RX(A + '(?:isocyanide|carbylamine)'), m => '[C-]#[N+]' + ALK[m[1]]],
  [RX(A + 'isocyanate'), m => 'O=C=N' + ALK[m[1]]], [RX(A + 'isothiocyanate'), m => 'S=C=N' + ALK[m[1]]],
  [RX(A + 'nitrite'), m => 'O=NO' + ALK[m[1]]], [RX(A + 'nitrate'), m => '[O-][N+](=O)O' + ALK[m[1]]], [RX(A + 'azide'), m => '[N-]=[N+]=N' + ALK[m[1]]],
  [RX(A + '(?:thiol|mercaptan)'), m => 'S' + ALK[m[1]]], [RX(A + 'hydrogensulfate'), m => 'OS(=O)(=O)O' + ALK[m[1]]], [RX(A + 'tosylate'), m => 'CC1=CC=C(C=C1)S(=O)(=O)O' + ALK[m[1]]],
  [RX(A + 'magnesium(chloride|bromide|iodide)'), m => ({ chloride: 'Cl', bromide: 'Br', iodide: 'I' })[m[2]] + '[Mg]' + ALK[m[1]]],
  [RX(A + 'lithium'), m => '[Li]' + ALK[m[1]]],
  [RX('(?:lithium)?di' + A + '(?:cuprate|copperlithium|copper)'), m => ({ special: 'R2CuLi', smi: ALK[m[1]] })],
  [RX('di' + A + 'cadmium'), m => ({ special: 'R2Cd', smi: ALK[m[1]] })],
  [RX('di' + A + 'ether'), m => 'O(' + ALK[m[1]] + ')' + ALK[m[1]]], [RX(A + A + 'ether'), m => 'O(' + ALK[m[1]] + ')' + ALK[m[2]]],
  [RX('di' + A + 'ketone'), m => 'C(=O)(' + ALK[m[1]] + ')' + ALK[m[1]]], [RX(A + A + 'ketone'), m => 'C(=O)(' + ALK[m[1]] + ')' + ALK[m[2]]],
  [RX('di' + A + 'sulfide'), m => 'S(' + ALK[m[1]] + ')' + ALK[m[1]]], [RX(A + A + 'sulfide'), m => 'S(' + ALK[m[1]] + ')' + ALK[m[2]]],
  [RX('di' + A + '(' + Object.keys(ANIONDI).join('|') + ')'), m => 'O(' + ALK[m[1]] + ')' + ANIONDI[m[2]] + 'O' + ALK[m[1]]],
  [RX(A + '(' + Object.keys(ANION).sort((a, b) => b.length - a.length).join('|') + ')'), m => 'O(' + ALK[m[1]] + ')' + ANION[m[2]]],
  [RX('(' + Object.keys(CATION).join('|') + ')(' + Object.keys(ALKOXIDE).sort((a, b) => b.length - a.length).join('|') + ')'), m => '[O-]' + ALKOXIDE[m[2]] + '.' + CATION[m[1]][0]],
  [RX('(' + Object.keys(CATION).join('|') + ')(' + Object.keys(ANION).sort((a, b) => b.length - a.length).join('|') + ')'), m => { const [ion, z] = CATION[m[1]]; const one = '[O-]' + ANION[m[2]]; return z === 2 ? one + '.' + one + '.' + ion : one + '.' + ion; }],
  [RX('(' + Object.keys(CATION).join('|') + ')(' + Object.keys(ANIONDI).join('|') + ')'), m => { const [ion, z] = CATION[m[1]]; return '[O-]' + ANIONDI[m[2]] + '[O-].' + (z === 2 ? ion : ion + '.' + ion); }],
  [RX('(sodium|potassium|lithium)(acetylide|ethynide)'), m => '[C-]#C.' + CATION[m[1]][0]],
  [RX('(' + Object.keys(ACYL).sort((a, b) => b.length - a.length).join('|') + ')(chloride|bromide|iodide)'), m => ({ chloride: 'Cl', bromide: 'Br', iodide: 'I' })[m[2]] + ACYL[m[1]]],
  [RX('(' + Object.keys(ACYLDI).join('|') + ')(?:di)?(chloride|bromide)'), m => { const X = m[2] === 'bromide' ? 'Br' : 'Cl'; return X + ACYLDI[m[1]] + X; }],
  [RX('n' + A + '(acetamide|benzamide|formamide|ethanamide)'), m => 'N(' + ALK[m[1]] + ')' + ({ acetamide: 'C(C)=O', ethanamide: 'C(C)=O', benzamide: 'C(=O)' + PH, formamide: 'C=O' })[m[2]]],
  [RX('nn?di' + A + '(acetamide|formamide|benzamide)'), m => 'N(' + ALK[m[1]] + ')(' + ALK[m[1]] + ')' + ({ acetamide: 'C(C)=O', benzamide: 'C(=O)' + PH, formamide: 'C=O' })[m[2]]],
  [RX(A + 'idenechloride'), m => null],
];
function fclass(k) { for (const [re, fn] of FC) { const m = k.match(re); if (m) { const r = fn(m); if (r == null) continue; if (typeof r === 'object') { const c = canon(r.smi); return c ? { special: r.special, smi: c } : null; } const c = canon(r); if (c) return c; } } return null; }

/* ---------------- IUPAC substitutive names ---------------- */
const STEM = { meth: 1, eth: 2, prop: 3, but: 4, pent: 5, hex: 6, hept: 7, oct: 8, non: 9, dec: 10, undec: 11, dodec: 12 };
const PRE = Object.assign({}, ALK, { fluoro: 'F', chloro: 'Cl', bromo: 'Br', iodo: 'I', nitro: NO2, nitroso: 'N=O', amino: 'N', hydroxy: 'O', methoxy: 'OC', ethoxy: 'OCC', propoxy: 'OCCC', isopropoxy: 'OC(C)C', butoxy: 'OCCCC', tertbutoxy: 'OC(C)(C)C', phenoxy: 'O' + PH, benzyloxy: 'OC' + PH,
  acetoxy: 'OC(C)=O', oxo: '=O', imino: '=N', cyano: 'C#N', formyl: 'C=O', carboxy: 'C(=O)O', methoxycarbonyl: 'C(=O)OC', ethoxycarbonyl: 'C(=O)OCC', acetyl: 'C(C)=O', benzoyl: 'C(=O)' + PH, mercapto: 'S', sulfanyl: 'S', sulfo: 'S(=O)(=O)O', azido: 'N=[N+]=[N-]',
  methylidene: '=C', methylene: '=C', ethylidene: '=CC', isopropylidene: '=C(C)C', benzylidene: '=C' + PH, chloromethyl: 'CCl', bromomethyl: 'CBr', hydroxymethyl: 'CO', aminomethyl: 'CN', trifluoromethyl: 'C(F)(F)F', trichloromethyl: 'C(Cl)(Cl)Cl',
  dimethylamino: 'N(C)C', diethylamino: 'N(CC)CC', methylamino: 'NC', ethylamino: 'NCC', acetamido: 'NC(C)=O', acetylamino: 'NC(C)=O', hydroxyimino: '=NO', diazo: null });
const PREK = Object.keys(PRE).filter(k => PRE[k]).sort((a, b) => b.length - a.length);
const MULT = { di: 2, tri: 3, tetra: 4, penta: 5, hexa: 6, bis: 2, tris: 3 };
const ARP = { benzene: null, phenol: 'O', aniline: 'N', benzenamine: 'N', toluene: 'C', methylbenzene: 'C', benzoicacid: 'C(=O)O', benzenecarboxylicacid: 'C(=O)O', benzaldehyde: 'C=O', benzenecarbaldehyde: 'C=O', benzonitrile: 'C#N', benzenecarbonitrile: 'C#N',
  anisole: 'OC', methoxybenzene: 'OC', acetophenone: 'C(C)=O', styrene: 'C=C', benzamide: 'C(N)=O', benzoylchloride: 'C(Cl)=O', benzenesulfonicacid: 'S(=O)(=O)O', benzylalcohol: 'CO', phenylmethanol: 'CO', benzylchloride: 'CCl', benzylamine: 'CN', phenylacetic: null };
const ARPK = Object.keys(ARP).filter(k => k !== 'phenylacetic').sort((a, b) => b.length - a.length);
function k1norm(s) {
  return String(s).toLowerCase().replace(/α/g, 'alpha').replace(/β/g, 'beta').replace(/γ/g, 'gamma').replace(/sulph/g, 'sulf').replace(/[–—−]/g, '-').replace(/\s+/g, '').replace(/['’`]/g, '')
    .replace(/tert-?butyl/g, 'tertbutyl').replace(/sec-?butyl/g, 'secbutyl').replace(/(^|[^a-z])t-butyl/g, '$1tertbutyl').replace(/iso-?propyl/g, 'isopropyl').replace(/(^|[^a-z])i-propyl/g, '$1isopropyl').replace(/tert-?butoxy/g, 'tertbutoxy').replace(/tert-?(?:pentyl|amyl)/g, 'tertpentyl');
}
/* prefix list: "2-bromo-2-methyl" -> [{locs:[2],name:'bromo'},{locs:[2],name:'methyl'}] plus dangling locants */
function prefixes(s) {
  const out = []; let i = 0, dang = null;
  while (i < s.length) {
    let locs = null; const lm = s.slice(i).match(/^((?:\d+[a-z]?'?|n)(?:,(?:\d+[a-z]?|n))*)-/);
    if (lm) { locs = lm[1].split(','); i += lm[0].length; }
    if (i >= s.length) { dang = locs; break; }
    let mult = 1; const mm = s.slice(i).match(/^(di|tri|tetra|penta|hexa|bis|tris)(?=[a-z(])/);
    let name = null, tpl = null;
    const rest0 = s.slice(i + (mm ? mm[0].length : 0));
    if (rest0[0] === '(') { const j = rest0.indexOf(')'); if (j < 0) return null; const inner = rest0.slice(1, j); tpl = subst(inner); if (!tpl) return null; name = '(' + inner + ')'; if (mm) mult = MULT[mm[1]]; i += (mm ? mm[0].length : 0) + j + 1; }
    else { for (const k of PREK) if (rest0.startsWith(k)) { name = k; break; }
      if (!name && mm) { // "di" could be part of a name like "diazo"; try without the multiplier
        for (const k of PREK) if (s.slice(i).startsWith(k)) { name = k; break; } if (name) { i += name.length; tpl = PRE[name]; out.push({ locs, mult: 1, name, tpl }); if (s[i] === '-') i++; continue; } }
      if (!name) return null; if (mm) mult = MULT[mm[1]]; tpl = PRE[name]; i += (mm ? mm[0].length : 0) + name.length; }
    out.push({ locs, mult, name, tpl }); if (s[i] === '-') i++;
  }
  return { list: out, dang };
}
/* a substituent written out: "4-chlorophenyl", "1-methylethyl", "2-methylpropyl", "1,1-dimethylethyl" -> template, attach first */
function subst(inner) {
  const s = inner; let m = s.match(/^(.*?)(meth|eth|prop|but|pent|hex)yl$/);
  if (s.endsWith('phenyl')) { const pre = s.slice(0, -6); const g = G(); const r = ring6(g); const P = pre ? prefixes(pre) : { list: [] }; if (!P) return null; for (const p of P.list) { const locs = p.locs || []; for (let k = 0; k < p.mult; k++) { const L = +(locs[k] || 4); if (!(L >= 2 && L <= 6)) return null; attach(g, r[L - 1], p.tpl); } } return rootSmi(g, r[0]); }
  if (m) { const n = STEM[m[2]]; const g = G(); const c = []; for (let k = 0; k < n; k++) { c.push(addA(g, 'C')); if (k) addB(g, c[k - 1], c[k], 1); } const P = m[1] ? prefixes(m[1]) : { list: [] }; if (!P) return null; for (const p of P.list) { const locs = p.locs || []; for (let k = 0; k < p.mult; k++) { const L = +(locs[k] || 1); if (!(L >= 1 && L <= n)) return null; attach(g, c[L - 1], p.tpl); } } return rootSmi(g, c[0]); }
  if (PRE[s]) return PRE[s];
  return null;
}
function rootSmi(g, at) { // SMILES string whose first atom is `at`
  const n = g.a.length; const order = [at, ...[...Array(n).keys()].filter(x => x !== at)]; const g2 = G(); const map = new Map(); order.forEach(x => { const a = g.a[x]; map.set(x, addA(g2, a.el, { chg: a.chg })); }); for (const [i, j, o] of g.b) addB(g2, map.get(i), map.get(j), o); return toSmi(g2); }
function ring6(g) { const r = []; for (let k = 0; k < 6; k++) r.push(addA(g, 'C')); for (let k = 0; k < 6; k++) addB(g, r[k], r[(k + 1) % 6], k % 2 === 0 ? 2 : 1); return r; }
function naph(g) { // positions 1..8, 4a=8, 8a=9 (array index L-1; 4a -> 8, 8a -> 9)
  const p = []; for (let k = 0; k < 10; k++) p.push(addA(g, 'C')); const [c1, c2, c3, c4, c5, c6, c7, c8, c4a, c8a] = p;
  addB(g, c1, c2, 2); addB(g, c2, c3, 1); addB(g, c3, c4, 2); addB(g, c4, c4a, 1); addB(g, c4a, c5, 1); addB(g, c5, c6, 2); addB(g, c6, c7, 1); addB(g, c7, c8, 2); addB(g, c8, c8a, 1); addB(g, c8a, c1, 1); addB(g, c4a, c8a, 2); return p; }
const SUF = [ // [suffix, kind, count, keepE]
  ['dicarboxylicacid', 'xCOOH', 2, 1], ['carboxylicacid', 'xCOOH', 1, 1], ['dicarbaldehyde', 'xCHO', 2, 1], ['carbaldehyde', 'xCHO', 1, 1], ['carbonitrile', 'xCN', 1, 1], ['carbonylchloride', 'xCOCl', 1, 1], ['carboxamide', 'xCONH2', 1, 1], ['sulfonicacid', 'xSO3H', 1, 1],
  ['dioicacid', 'COOH', 2, 1], ['oicacid', 'COOH', 1, 0], ['dioylchloride', 'COCl', 2, 1], ['oylchloride', 'COCl', 1, 0], ['oylbromide', 'COBr', 1, 0], ['diamide', 'CONH2', 2, 1], ['amide', 'CONH2', 1, 0], ['dinitrile', 'CN', 2, 1], ['nitrile', 'CN', 1, 1],
  ['dial', 'CHO', 2, 1], ['al', 'CHO', 1, 0], ['trione', 'ONE', 3, 1], ['dione', 'ONE', 2, 1], ['one', 'ONE', 1, 0], ['triol', 'OH', 3, 1], ['diol', 'OH', 2, 1], ['ol', 'OH', 1, 0], ['triamine', 'NH2', 3, 1], ['diamine', 'NH2', 2, 1], ['amine', 'NH2', 1, 0], ['dithiol', 'SH', 2, 1], ['thiol', 'SH', 1, 1]];
function iupac(name, opt) {
  opt = opt || {}; let s = k1norm(name); let ez = null;
  s = s.replace(/^\((?:\d*[ezrs](?:,\d*[ezrs])*|[±+-])\)-?/, m0 => { ez = m0; return ''; }).replace(/^(?:cis|trans|meso|dl|d|l)-/, m0 => { ez = m0; return ''; });
  let omp = null; const om = s.match(/^(o|m|p|ortho|meta|para)-/); if (om) { omp = { o: 2, m: 3, p: 4 }[om[1][0]]; s = s.slice(om[0].length); }
  // aromatic parents
  for (const P of ARPK) {
    const re = new RegExp('^(.*?)' + P + '(?:-((?:\\d+)(?:,\\d+)*)-)?((?:di|tri)?(?:ol|amine|carboxylicacid|carbaldehyde|carbonitrile|sulfonicacid))?$'); const m = s.match(re);
    if (!m || (m[3] && P !== 'benzene')) continue;
    const res = (() => {
    const g = G(); const r = ring6(g); let used = new Set();
    const imp = ARP[P]; let nAt = null; if (imp) { const ids = attach(g, r[0], imp); used.add(1); if (/^N/.test(imp) || imp === 'C(N)=O') nAt = ids[imp === 'C(N)=O' ? 1 : 0]; }
    if (m[3]) { const sm = m[3].match(/^(di|tri)?(.*)$/); const cnt = sm[1] ? MULT[sm[1]] : 1; const locs = m[2] ? m[2].split(',').map(Number) : cnt === 1 ? [1] : cnt === 2 ? [1, omp ? omp - 1 + 1 : 2] : [1, 3, 5];
      const T = { ol: 'O', amine: 'N', carboxylicacid: 'C(=O)O', carbaldehyde: 'C=O', carbonitrile: 'C#N', sulfonicacid: 'S(=O)(=O)O' }[sm[2]]; if (locs.length !== cnt) return null;
      if (cnt === 2 && omp && !m[2]) locs[1] = omp; for (const L of locs) { if (!(L >= 1 && L <= 6)) return null; const ids = attach(g, r[L - 1], T); used.add(L); if (T === 'N') nAt = ids[0]; } }
    const pre = m[1] ? prefixes(m[1].replace(/-$/, '')) : { list: [] }; if (!pre) return null;
    let auto = imp ? (omp || 4) : null; let first = true;
    for (const p of pre.list) { let locs = (p.locs || []).filter(x => x !== 'n').map(Number);
      if ((p.locs || []).includes('n')) { if (nAt == null) { if (!(p.locs.length === 1 && ['propyl', 'butyl', 'pentyl', 'hexyl'].includes(p.name))) return null; locs = []; } else { for (let k = 0; k < p.mult; k++) attach(g, nAt, p.tpl); continue; } }
      if (!locs.length) { if (imp) { for (let k = 0; k < p.mult; k++) { let L = (pre.list.length === 1 && p.mult === 1) ? auto : null; if (L == null || used.has(L)) { L = [4, 2, 6, 3, 5].find(x => !used.has(x)); } locs.push(L); used.add(L); } }
        else { if (p.mult === 1) { const L = first ? 1 : (omp ? omp - 1 + 1 : [2, 3, 4, 5, 6].find(x => !used.has(x))); locs = [first ? 1 : (omp || L)]; } else if (p.mult === 2) locs = omp ? [1, omp] : [1, 4]; else if (p.mult === 3) locs = [1, 3, 5]; else locs = [1, 2, 3, 4, 5, 6].slice(0, p.mult); } }
      if (locs.length !== p.mult) return null; for (const L of locs) { if (!(L >= 1 && L <= 6)) return null; attach(g, r[L - 1], p.tpl); used.add(L); } first = false; }
    return outSmi(g); })();
    if (res) return res;
  }
  // naphthalene
  { const m = s.match(/^(.*?)naphthalene$/); if (m) { const g = G(); const p = naph(g); const pre = m[1] ? prefixes(m[1].replace(/-$/, '')) : { list: [] }; if (!pre) return null; for (const q of pre.list) { const locs = (q.locs || ['1']).map(Number); if (locs.length !== q.mult) { if (q.mult === 1 && !q.locs) { } else return null; } for (const L of locs) { if (!(L >= 1 && L <= 8)) return null; attach(g, p[L - 1], q.tpl); } } return outSmi(g); } }
  // chain or ring parent
  let suf = null, sufLocs = null; let rest = s;
  for (const [sx, kind, cnt, keepE] of SUF) { const re = new RegExp('^(.*?)(?:-((?:\\d+)(?:,\\d+)*)-)?' + sx + '$'); const m = rest.match(re); if (m && m[1]) { suf = { kind, cnt, sx }; sufLocs = m[2] ? m[2].split(',').map(Number) : null; rest = m[1]; if (keepE && rest.endsWith('e')) rest = rest.slice(0, -1); break; } }
  if (!suf) { if (!/(?:ane|ene|yne|diene|triene|diyne|enyne)$/.test(rest)) return null; rest = rest.replace(/e$/, ''); }
  // unsaturation
  const um = rest.match(/(a)?(?:-((?:\d+)(?:,\d+)*)-)?(an|en|yn|dien|trien|diyn)(?:-((?:\d+)(?:,\d+)*)-(yn|diyn))?$/); if (!um) return null;
  rest = rest.slice(0, rest.length - um[0].length); if (um[1] && !um[2]) { /* buta-diene style without locants */ }
  const sm = rest.match(/(cyclo)?(meth|eth|prop|but|pent|hex|hept|oct|non|undec|dodec|dec)$/); if (!sm) return null;
  const n = STEM[sm[2]], cyc = !!sm[1]; rest = rest.slice(0, rest.length - sm[0].length);
  const pre = rest ? prefixes(rest.replace(/-$/, '') + (rest.endsWith('-') ? '-' : '')) : { list: [], dang: null }; if (!pre) return null;
  let dang = pre.dang ? pre.dang.map(Number) : null;
  const g = G(); const c = []; for (let k = 0; k < n; k++) { c.push(addA(g, 'C')); if (k) addB(g, c[k - 1], c[k], 1); } if (cyc) { if (n < 3) return null; addB(g, c[n - 1], c[0], 1); }
  const at = L => (L >= 1 && L <= n) ? c[L - 1] : null; const nxt = L => cyc ? c[L % n] : at(L + 1);
  const setO = (L, o) => { const a = at(L), b = nxt(L); if (a == null || b == null) return false; const bd = bond(g, a, b); if (!bd) return false; bd[2] = o; return true; };
  // suffix locants: given, or dangling ("2-propanol"), or default
  const isTerm = suf && ['COOH', 'COCl', 'COBr', 'CONH2', 'CN', 'CHO'].includes(suf.kind);
  if (suf && !sufLocs && !isTerm && dang && dang.length === suf.cnt) { sufLocs = dang; dang = null; }
  // unsaturation locants
  const uk = um[3]; const ucnt = { an: 0, en: 1, yn: 1, dien: 2, trien: 3, diyn: 2 }[uk]; let ulocs = um[2] ? um[2].split(',').map(Number) : null;
  if (!ulocs && ucnt && dang && dang.length === ucnt) { ulocs = dang; dang = null; }
  if (dang) return null;
  if (ucnt && !ulocs) ulocs = ucnt === 1 ? [1] : ucnt === 2 ? (uk === 'dien' && !cyc ? [1, 3] : [1, 3]) : [1, 3, 5];
  if (ucnt && ulocs.length !== ucnt) return null;
  for (const L of ulocs || []) if (!setO(L, uk.includes('yn') ? 3 : 2)) return null;
  if (um[5]) { const l2 = um[4].split(',').map(Number); for (const L of l2) if (!setO(L, 3)) return null; }
  // suffix groups
  let nAt = null;
  if (suf) {
    const k = suf.kind;
    if (isTerm) { const ends = suf.cnt === 2 ? [1, n] : [1]; if (cyc) return null;
      for (const L of ends) { const a = at(L); if (k === 'CN') { const ids = frag(g, 'N', a, 3); } else if (k === 'CHO') frag(g, 'O', a, 2); else { frag(g, 'O', a, 2); if (k === 'COOH') { if (opt.esterR) frag(g, 'O' + opt.esterR, a, 1); else if (opt.salt) { const ids = frag(g, '[O-]', a, 1); } else frag(g, 'O', a, 1); } else if (k === 'COCl') frag(g, 'Cl', a, 1); else if (k === 'COBr') frag(g, 'Br', a, 1); else if (k === 'CONH2') nAt = frag(g, 'N', a, 1)[0]; } } }
    else { let locs = sufLocs || (k === 'ONE' ? (cyc ? [1, 2, 3].slice(0, suf.cnt) : suf.cnt === 1 ? [n >= 3 ? 2 : 1] : suf.cnt === 2 ? (n === 4 ? [2, 3] : [2, n - 1]) : [2, 3, 4]) : suf.cnt === 1 ? [1] : suf.cnt === 2 ? (cyc ? [1, 2] : [1, 2]) : [1, 2, 3]);
      if (locs.length !== suf.cnt) return null;
      for (const L of locs) { const a = at(L); if (a == null) return null;
        if (k === 'ONE') frag(g, 'O', a, 2); else if (k === 'OH') frag(g, 'O', a, 1); else if (k === 'SH') frag(g, 'S', a, 1); else if (k === 'NH2') nAt = frag(g, 'N', a, 1)[0];
        else if (k === 'xCOOH') { const ids = frag(g, 'C(=O)' + (opt.esterR ? 'O' + opt.esterR : opt.salt ? '[O-]' : 'O'), a, 1); } else if (k === 'xCHO') frag(g, 'C=O', a, 1); else if (k === 'xCN') frag(g, 'C#N', a, 1); else if (k === 'xCOCl') frag(g, 'C(=O)Cl', a, 1); else if (k === 'xCONH2') nAt = frag(g, 'C(=O)N', a, 1)[1]; else if (k === 'xSO3H') frag(g, 'S(=O)(=O)O', a, 1); } }
  }
  // prefixes
  let firstAlk = true;
  for (const p of pre.list) {
    let locs = p.locs ? p.locs.slice() : null;
    if (locs && locs.includes('n')) { if (nAt == null) { if (locs.length === 1 && ['propyl', 'butyl', 'pentyl', 'hexyl'].includes(p.name)) locs = null; else return null; } else { const nn = locs.filter(x => x === 'n').length; for (let k = 0; k < Math.max(nn, p.mult); k++) attach(g, nAt, p.tpl); continue; } }
    locs = locs ? locs.map(Number) : null;
    if (!locs) { const alk = /^[A-Z]/.test(p.tpl) && p.tpl[0] === 'C' && !['C#N', 'C=O', 'C(=O)O'].includes(p.tpl); const start = (alk && !cyc && n > 2) ? 2 : 1; locs = []; for (let k = 0; k < p.mult; k++) locs.push(p.mult > 1 && alk && !cyc ? start : start + (cyc ? k : 0)); if (!alk && !cyc) { locs = []; for (let k = 0; k < p.mult; k++) locs.push(1 + Math.floor(k / 3)); } }
    if (locs.length !== p.mult) return null;
    for (const L of locs) { const a = at(L); if (a == null) return null; attach(g, a, p.tpl); }
    firstAlk = false;
  }
  if (opt.salt) { frag(g, opt.salt, null); }
  // valence check
  for (let i = 0; i < g.a.length; i++) { const a = g.a[i]; if (a.el === 'C' && bsum(g, i) > 4) return null; if (a.el === 'O' && !a.chg && bsum(g, i) > 2) return null; if (a.el === 'N' && !a.chg && bsum(g, i) > 3) return null; }
  return outSmi(g);
}
/* esters / salts of IUPAC acids: "ethyl propanoate", "sodium 2-methylpropanoate", "methyl 3-oxobutanoate" */
function esterName(k1) {
  const s = k1norm(k1); const ALKS = Object.keys(ALK).sort((a, b) => b.length - a.length);
  if (/oate$/.test(s) || /carboxylate$/.test(s)) {
    for (const cat of Object.keys(CATION)) if (s.startsWith(cat)) { const acid = s.slice(cat.length).replace(/^-/, '').replace(/oate$/, 'oicacid').replace(/carboxylate$/, 'carboxylicacid'); const [ion] = CATION[cat]; const r = iupac(acid, { salt: ion }); if (r) return r; }
    for (const a of ALKS) if (s.startsWith(a)) { const acid = s.slice(a.length).replace(/^-/, '').replace(/oate$/, 'oicacid').replace(/carboxylate$/, 'carboxylicacid'); const r = iupac(acid, { esterR: ALK[a] }); if (r) return r; }
  }
  return null;
}

/* ---------------- condensed formulas ---------------- */
const CABBR = [['tert-c4h9', 'C(C)(C)C'], ['sec-c4h9', 'C(C)CC'], ['iso-c3h7', 'C(C)C'], ['n-c3h7', 'CCC'], ['i-c3h7', 'C(C)C'], ['n-c4h9', 'CCCC'], ['t-c4h9', 'C(C)(C)C'], ['i-c4h9', 'CC(C)C'], ['s-c4h9', 'C(C)CC'],
  ['c10h7', 'C1=CC=CC2=CC=CC=C12'], ['c6h11', 'C1CCCCC1'], ['c6h13', 'CCCCCC'], ['c5h11', 'CCCCC'], ['c6h5', PH], ['c4h9', 'CCCC'], ['c3h7', 'CCC'], ['c3h5', 'CC=C'], ['c2h5', 'CC'], ['c2h3', 'C=C'],
  ['sec-bu', 'C(C)CC'], ['tert-bu', 'C(C)(C)C'], ['t-bu', 'C(C)(C)C'], ['n-bu', 'CCCC'], ['i-bu', 'CC(C)C'], ['s-bu', 'C(C)CC'], ['i-pr', 'C(C)C'], ['n-pr', 'CCC'],
  ['co2me', 'C(=O)OC'], ['co2et', 'C(=O)OCC'], ['co2h', 'C(=O)O'], ['cooet', 'C(=O)OCC'], ['coome', 'C(=O)OC'], ['so3h', 'S(=O)(=O)O'], ['so2cl', 'S(=O)(=O)Cl'], ['so2nh2', 'S(N)(=O)=O'], ['nhac', 'NC(C)=O'], ['nme2', 'N(C)C'], ['net2', 'N(CC)CC'], ['nhme', 'NC'],
  ['oac', 'OC(C)=O'], ['ome', 'OC'], ['oet', 'OCC'], ['ots', 'OS(=O)(=O)C1=CC=C(C)C=C1'], ['no2', NO2], ['mgbr', '[Mg]Br'], ['mgcl', '[Mg]Cl'], ['mgi', '[Mg]I'], ['mgx', '[Mg]Br'],
  ['tbu', 'C(C)(C)C'], ['sbu', 'C(C)CC'], ['ibu', 'CC(C)C'], ['nbu', 'CCCC'], ['ipr', 'C(C)C'], ['npr', 'CCC'], ['bu', 'CCCC'], ['pr', 'CCC'], ['me', 'C'], ['et', 'CC'], ['ph', PH], ['bn', 'C' + PH], ['bz', 'C(=O)' + PH], ['ac', 'C(C)=O'], ['ts', 'S(=O)(=O)C1=CC=C(C)C=C1']];
const EL2 = ['cl', 'br', 'mg', 'li', 'na', 'ca', 'ag', 'cu', 'zn', 'hg', 'cd', 'ba'], EL1 = ['c', 'h', 'o', 'n', 's', 'p', 'f', 'i', 'k', 'b'];
const ELN = { c: 'C', h: 'H', o: 'O', n: 'N', s: 'S', p: 'P', f: 'F', i: 'I', k: 'K', b: 'B', cl: 'Cl', br: 'Br', mg: 'Mg', li: 'Li', na: 'Na', ca: 'Ca', ag: 'Ag', cu: 'Cu', zn: 'Zn', hg: 'Hg', cd: 'Cd', ba: 'Ba' };
function ctoks(t) { // -> tokens
  const T = []; let i = 0;
  while (i < t.length) { const c = t[i];
    if ('()'.includes(c)) { if (c === ')') { const m = t.slice(i + 1).match(/^\d+/); T.push({ k: ')', n: m ? +m[0] : 1 }); i += 1 + (m ? m[0].length : 0); } else { T.push({ k: '(' }); i++; } continue; }
    if (c === '=' || c === '#') { T.push({ k: 'bo', o: c === '=' ? 2 : 3 }); i++; continue; }
    if (c === '.') { T.push({ k: '.' }); i++; continue; }
    if (c === '+' || c === '-') { const nx = t[i + 1]; if (c === '-' && nx && /[a-z(]/.test(nx) && !(T.length && T[T.length - 1].k === 'chg')) { i++; continue; } T.push({ k: 'chg', v: c === '+' ? 1 : -1 }); i++; continue; }
    if (i === 0 || '(.'.includes(t[i - 1])) { let st = null; for (const [a, tpl] of [['hooc', 'C(=O)O'], ['ohc', 'C=O'], ['o2n', NO2], ['meooc', 'C(=O)OC'], ['etooc', 'C(=O)OCC'], ['h2nco', 'C(N)=O'], ['ncs', 'N=C=S']]) if (t.startsWith(a, i)) { st = [a, tpl]; break; }
      if (st) { T.push({ k: 'g', tpl: st[1], n: 1, name: st[0] }); i += st[0].length; continue; } }
    if (t.startsWith('c6h4', i)) { T.push({ k: 'ph2' }); i += 4; continue; }
    let hit = null; for (const [a, tpl] of CABBR) if (t.startsWith(a, i)) { const nx = t[i + a.length]; if (a.length <= 2 && nx && /[a-z]/.test(nx) && !['h', 'c', 'o', 'n', 's', 'b', 'i', 'f', 'p', 'k'].includes(nx)) continue; hit = [a, tpl]; break; }
    if (hit) { let j = i + hit[0].length; const m = t.slice(j).match(/^\d+/); T.push({ k: 'g', tpl: hit[1], n: m ? +m[0] : 1, name: hit[0] }); i = j + (m ? m[0].length : 0); continue; }
    let el = null; if (EL2.includes(t.slice(i, i + 2))) el = t.slice(i, i + 2); else if (EL1.includes(c)) el = c; if (!el) return null;
    i += el.length; let h = 0, n = 1;
    if (el === 'h') { const m = t.slice(i).match(/^\d+/); const hn = m ? +m[0] : 1; i += m ? m[0].length : 0; T.push({ k: 'H', n: hn }); continue; }
    if (['c', 'n', 'o', 's', 'p', 'b'].includes(el) && t[i] === 'h') { i++; const m = t.slice(i).match(/^\d+/); h = m ? +m[0] : 1; i += m ? m[0].length : 0; }
    const m2 = t.slice(i).match(/^\d+/); if (m2) { n = +m2[0]; i += m2[0].length; }
    T.push({ k: 'a', el: ELN[el], h, n }); }
  return T;
}
function condensed(raw) {
  let t = String(raw).replace(/\s+/g, '').replace(/[–—−]/g, '-').replace(/≡/g, '#').replace(/[₀-₉]/g, d => '₀₁₂₃₄₅₆₇₈₉'.indexOf(d)).replace(/⁺/g, '+').replace(/⁻/g, '-');
  let pos = 4; const om = t.match(/^(o|m|p|ortho|meta|para)-/i); if (om) { pos = { o: 2, m: 3, p: 4 }[om[1][0].toLowerCase()]; t = t.slice(om[0].length); }
  t = t.toLowerCase(); if (MF[t]) return canon(MF[t]);
  if (!/c(?![laurdsoe])|c(?=o|s)|ph|me|et|bu|pr|bn|bz/.test(t)) return null;
  const T = ctoks(t); if (!T || !T.length) return null;
  const g = G(); const free = [];
  const mk = (el, h, chg) => { const id = addA(g, el, { chg: chg || 0 }); free[id] = cval(el, chg || 0) - (h || 0); return id; };
  const link = (i, j, o) => { addB(g, i, j, o); free[i] -= o; free[j] -= o; };
  const placeTpl = (tpl) => { const ids = frag(g, tpl, null); for (const id of ids) free[id] = 0; const head = ids[0]; free[head] = 1; return { ids, head }; };
  function seq(T, a0, a1) { // parse T[a0..a1) into the graph; returns { first, last, att }
    let prev = null, pending = [], order = 1, preH = 0, salt = false, firstAtom = null, lastAtom = null, nextPort = null;
    const bondTo = (x, o) => { // bond new atom x to pending / prev
      if (pending.length) { for (const p of pending) { if (free[p.at] < 1 || free[x] < 1) return false; link(p.at, x, p.o || 1); } pending = []; return true; }
      if (prev != null && free[prev] < 1 && nextPort != null) { prev = nextPort; nextPort = null; }
      if (prev != null) { const o2 = o || 1; if (free[prev] < o2 || free[x] < o2) return false; link(prev, x, o2); if (free[prev] < 1 && nextPort != null) { prev = nextPort; nextPort = null; } return true; } return true; };
    for (let i = a0; i < a1; i++) { const tk = T[i];
      if (tk.k === 'bo') { order = tk.o; continue; }
      if (tk.k === '.') { prev = null; pending = []; order = 1; salt = false; continue; }
      if (tk.k === 'H') { preH += tk.n; continue; }
      if (tk.k === 'chg') { if (lastAtom == null) return null; const a = g.a[lastAtom]; a.chg += tk.v; free[lastAtom] = cval(a.el, a.chg) - (a._h || 0) - bsum(g, lastAtom); if (tk.v > 0) salt = true; continue; }
      if (tk.k === '(') { let d = 1, j = i + 1; for (; j < a1; j++) { if (T[j].k === '(') d++; else if (T[j].k === ')') { d--; if (!d) break; } } if (j >= a1) return null; const mult = T[j].n || 1; const leadBo = T[i + 1] && T[i + 1].k === 'bo' ? T[i + 1].o : 1;
        for (let r = 0; r < mult; r++) { const sub = seq(T, i + 1 + (leadBo > 1 ? 1 : 0), j); if (!sub) return null; const fr = sub.ids.filter(x => free[x] > 0);
          if (fr.length === 1 && free[fr[0]] >= 2 && prev != null && leadBo === 1 && !(sub.ids.length === 1 && g.a[sub.ids[0]].el === 'O' && g.a[prev].el === 'C' && free[prev] >= 2)) { // chain unit, e.g. (CH2)4
            if (!bondTo(sub.first, order)) return null; order = 1; prev = sub.last != null && free[sub.last] > 0 ? sub.last : fr[0]; continue; }
          if (fr.length === 1 && free[fr[0]] >= 2 && prev == null && pending.length === 0 && mult > 1) { pending.push({ at: fr[0], o: 1 }); continue; }
          const att = fr.length ? fr[0] : null; if (att == null) return null;
          if (prev == null) { if (leadBo > 1) return null; pending.push({ at: att, o: 1 }); continue; }
          let o = leadBo; if (o === 1 && sub.ids.length === 1 && ['O', 'S'].includes(g.a[att].el) && g.a[prev].el === 'C' && free[prev] >= 2 && free[att] >= 2) o = 2;
          if (free[prev] < o || free[att] < o) return null; link(prev, att, o); }
        i = j; continue; }
      if (tk.k === ')') return null;
      if (tk.k === 'ph2') { const ids = frag(g, PH, null); for (const id of ids) free[id] = 0; const r0 = ids[0], rp = ids[pos - 1]; free[r0] = 1; free[rp] = 1; if (firstAtom == null) firstAtom = r0; lastAtom = rp;
        if (prev == null && !pending.length) { prev = r0; nextPort = rp; continue; }
        if (!bondTo(r0, order)) return null; order = 1; prev = rp; continue; }
      if (tk.k === 'g') { for (let r = 0; r < tk.n; r++) { const { ids, head } = placeTpl(tk.tpl); lastAtom = ids[ids.length - 1]; if (firstAtom == null) firstAtom = head;
          if (tk.name === 'mgbr' || tk.name === 'mgcl' || tk.name === 'mgi' || tk.name === 'mgx') { if (prev == null) return null; if (!bondTo(head, 1)) return null; continue; }
          const fh = free[head]; if (prev == null && !pending.length) { if (r === 0 && tk.n === 1) { prev = head; } else pending.push({ at: head, o: 1 }); continue; }
          if (!bondTo(head, order)) return null; order = 1; if (free[head] > 0) prev = head; } continue; }
      if (tk.k === 'a') { let n = tk.n;
        // diazonium: Ar-N2(+)
        if (tk.el === 'N' && n === 2 && !tk.h && prev != null) { const n1 = mk('N', 0, 1), n2 = mk('N', 0); if (!bondTo(n1, 1)) return null; link(n1, n2, 3); prev = null; salt = true; lastAtom = n1; if (T[i + 1] && T[i + 1].k === 'chg' && T[i + 1].v > 0) i++; continue; }
        if (['F', 'Cl', 'Br', 'I'].includes(tk.el) && n > 1 && prev == null && !pending.length) { for (let r = 0; r < n; r++) { const x = mk(tk.el, 0, 0); if (firstAtom == null) firstAtom = x; lastAtom = x; pending.push({ at: x, o: 1 }); } continue; }
        if (n > 1 && !['F', 'Cl', 'Br', 'I', 'O', 'H'].includes(tk.el) && !(tk.el === 'C' && false)) { if (!(tk.el === 'N' && false)) return null; }
        for (let r = 0; r < n; r++) {
          const el = tk.el; const h = (tk.h || 0) + (r === 0 ? preH : 0); if (r === 0) preH = 0;
          let chg = 0; const nx = T[i + 1]; if (r === n - 1 && nx && nx.k === 'chg' && !(METALS.has(el))) { chg = 0; }
          // halide counter-ion after a cation
          if (['Cl', 'Br', 'I'].includes(el) && (salt || (prev != null && free[prev] < 1))) { const x = mk(el, 0, -1); lastAtom = x; if (T[i + 1] && T[i + 1].k === 'chg' && T[i + 1].v < 0) i++; continue; }
          const x = mk(el, h, 0); g.a[x]._h = h; lastAtom = x; if (firstAtom == null) firstAtom = x;
          // look-ahead charge on this atom ("CH3COO-", "CH3NH3+", "(CH3)4N+")
          if (r === n - 1 && T[i + 1] && T[i + 1].k === 'chg') { g.a[x].chg += T[i + 1].v; free[x] = cval(el, g.a[x].chg) - h; if (T[i + 1].v > 0) salt = true; i++; }
          if (prev == null && !pending.length) { prev = x; continue; }
          let o = order; order = 1;
          if (o === 1 && prev != null && !pending.length) { // carbonyl / nitrile / nitroso by valence
            if ((el === 'O' || el === 'S') && !h && !g.a[x].chg && free[prev] >= 2 && ['C', 'N', 'S'].includes(g.a[prev].el) && free[x] === 2 && !(g.a[prev].el === 'C' && g.b.some(b => b[2] === 2 && (b[0] === prev || b[1] === prev) && g.a[b[0] === prev ? b[1] : b[0]].el === 'O'))) { const nxt = T[i + 1]; const chainO = nxt && (nxt.k === 'a' || nxt.k === 'g' || nxt.k === '(') && free[prev] === 2 && g.a[prev].el === 'C' && !(nxt.k === 'a' && nxt.el === 'O') && false; if (!chainO) o = 2; }
            if (el === 'N' && !h && !g.a[x].chg && free[prev] >= 3 && g.a[prev].el === 'C' && !(T[i + 1] && (T[i + 1].k === 'a' || T[i + 1].k === 'g' || T[i + 1].k === '('))) o = 3; }
          if (!bondTo(x, o)) return null;
          if (free[x] > 0) prev = x;
        } continue; }
    }
    const ids = []; for (let k = 0; k < g.a.length; k++) ids.push(k);
    return { first: firstAtom, last: lastAtom, ids: ids.filter(k => k >= (firstAtom != null ? Math.min(firstAtom, ...[lastAtom].filter(x => x != null)) : 0)) };
  }
  const r = seq(T, 0, T.length); if (!r) return null;
  // implicit unsaturation between neighbours that both still have free valence
  for (let pass = 0; pass < 3; pass++) for (const b of g.b) { const [i, j] = b; const d = Math.min(free[i], free[j], 3 - b[2]); if (d > 0 && !METALS.has(g.a[i].el) && !METALS.has(g.a[j].el)) { b[2] += d; free[i] -= d; free[j] -= d; } }
  // isocyanide R-N=C with left-over valence -> R-[N+]#[C-]
  for (const b of g.b) { const [i, j] = b; for (const [n1, c1] of [[i, j], [j, i]]) if (g.a[n1].el === 'N' && g.a[c1].el === 'C' && free[c1] === 2 && free[n1] === 0 && b[2] === 2 && g.b.filter(x => x[0] === c1 || x[1] === c1).length === 1) { b[2] = 3; g.a[n1].chg = 1; g.a[c1].chg = -1; free[n1] = 0; free[c1] = 0; } }
  // ionic metal salts: O/N/S/C≡ - Na/K/Li/Ca/Ag
  for (let k = g.b.length - 1; k >= 0; k--) { const [i, j] = g.b[k]; for (const [m, x] of [[i, j], [j, i]]) if (IONIC.has(g.a[m].el) && ['O', 'N', 'S', 'C'].includes(g.a[x].el)) { g.b.splice(k, 1); g.a[x].chg -= 1; g.a[m].chg += 1; free[m] = 0; free[x] = 0; break; } }
  for (let k = 0; k < g.a.length; k++) { if (free[k] == null) continue; if (METALS.has(g.a[k].el)) continue; if (free[k] !== 0) { if (g.a[k].chg) continue; return null; } }
  return outSmi(g);
}

/* ---------------- main entry ---------------- */
const looksSmiles = s => /^[A-Za-z0-9@+\-\[\]\(\)=#$\/\\%.:*]+$/.test(s) && (/[a-z]\d|\d[a-z]|\[|@|[cnos]1/.test(s) || /[cnosp](?![a-z])/.test(s.replace(/cl|br/gi, '')) || !/H|[a-z]/.test(s.replace(/Cl|Br/g, '')));
function read(text, o) {
  o = o || {}; const raw = String(text || '').trim(); if (!raw || raw.length > 120) return null;
  const k0 = norm0(raw.replace(/^(?:\((?:\d*[ezrs](?:,\d*[ezrs])*|[±+-])\)-?|cis-|trans-|meso-|dl-|d-|l-)/i, '')).replace(/^(?:alpha|beta)(?=d|l)/, '');
  if (TRIV.has(k0)) { const s = canon(TRIV.get(k0)); if (s) return { smi: s, how: 'name' }; }
  if (MF[k0]) { const s = canon(MF[k0]); if (s) return { smi: s, how: 'formula' }; }
  const fc = fclass(k0); if (fc) return typeof fc === 'object' ? { smi: fc.smi, how: 'name', special: fc.special } : { smi: fc, how: 'name' };
  if (/[a-z]{3}/i.test(raw)) { let s = null; try { s = esterName(raw) || iupac(raw); } catch (e) { s = null; } if (s) return { smi: s, how: 'IUPAC name' }; }
  if (/^[a-z0-9()=#≡\-–+.,₀-₉⁺⁻]+$/i.test(raw.replace(/\s+/g, ''))) { let s = null; try { s = condensed(raw); } catch (e) { s = null; } if (s) return { smi: s, how: 'formula' }; }
  if (o.smiles !== false && looksSmiles(raw)) { const s = canon(raw); if (s) return { smi: s, how: 'SMILES' }; }
  return null;
}
const parse = t => { const r = read(t); return r ? r.smi : null; };
root.OCName = { parse, read, condensed, iupac, fclass, norm0, TRIV };
})(typeof window !== 'undefined' ? window : globalThis);
