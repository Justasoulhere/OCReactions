/* engine_r14.js: run 32. Chemistry gaps from the run-31 benchmark failures (Gambit + SKM) and the reagent matrix:
   organometallics + D2O, alkoxide + R-X, carbonate aldol, Gilman, Stork enamine, Darzens, barbituric acid, cumene
   process, Ca-salt distillation, haloform/decarboxylation on every site, full hydrolysis of peptides and glycosides,
   lactones, ring-closing alkoxides and malonates, alpha-halogenation, Br2 on activated rings, allylic and SO2Cl2
   chlorination, Rosenmund + Lindlar, alkyne + KMnO4, O3/H2O oxidative work-up, AgNO3 tests, SNAr with amines,
   aromatisation, iodolactonisation, carbylamine-negative amides. Recipes at the FRONT correct a wrong older answer. */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore; if (!C || !E) return;
const { BL, el, hTot } = E; const RP = () => root.RxParse;
const fresh = smi => { try { const m = C.fromSmilesMB(smi, BL); return m && m.atoms.size ? m : null; } catch (e) { return null; } };
const step = (L, R, text, cond) => ({ L: L.clone(), arrows: [], R: R.clone(), text, cond: cond || [] });
const out = o => Object.assign({ steps: [], minor: [], extra: [], why: [], db: [] }, o);
const M = (m, sm) => { try { return C.matchAll(m, sm, {}); } catch (e) { return []; } };
const RXC = new Map();
function rxn(sma, smis, max) {
  const RD = C.RD; let rx = RXC.get(sma); if (rx === undefined) { try { rx = RD.get_rxn(sma); } catch (e) { rx = null; } RXC.set(sma, rx); } if (!rx) return [];
  const L = new RD.MolList(); const ms = smis.map(s => RD.get_mol(s)); if (ms.some(x => !x)) { ms.forEach(x => x && x.delete()); return []; } ms.forEach(x => L.append(x));
  const res = []; let o = null; try { o = rx.run_reactants(L, max || 24); } catch (e) { o = null; }
  if (o) for (let i = 0; i < o.size(); i++) { const ml = o.get(i); const ps = []; let bad = false;
    for (let k = 0; k < ml.size(); k++) { const x = ml.at(k); let s = null; try { const q = RD.get_mol(x.get_smiles()); if (q) { s = q.get_smiles(); q.delete(); } } catch (e) { s = null; } x.delete && x.delete(); if (!s) bad = true; ps.push(s); }
    if (!bad && !res.some(r => r.join('.') === ps.join('.'))) res.push(ps); }
  ms.forEach(x => x.delete()); return res;
}
const canonS = s => { try { const q = C.RD.get_mol(s); if (!q) return null; const k = q.get_smiles(); q.delete(); return k; } catch (e) { return null; } };
const smiOf = m => { try { return C.canon(m); } catch (e) { return null; } };
const nC = s => (String(s).replace(/Cl|Ca|Cu|Cr|Cs|Cd/g, '').match(/[Cc]/g) || []).length;
const add = r => E.REC.push(r), front = r => E.REC.unshift(r);
const byC = arr => arr.slice().sort((a, b) => nC(b) - nC(a) || b.length - a.length);
function finish(m, prodSmi, o) { const R = fresh(prodSmi); if (!R) return null; return out(Object.assign({}, o, { product: R, steps: [step(m, R, o.text || o.title, o.cond)], extra: (o.extra || []).map(fresh).filter(Boolean), minor: (o.minor || []).map(fresh).filter(Boolean) })); }
const none = (title, why) => out({ title, product: null, none: why });
const positive = (title, why) => out({ title, product: null, none: why, positive: true, tags: ['test'] });
function applyAll(smi, T, lim) { let cur = smi; const extra = [], labels = [];
  for (const [sma, label] of T) for (let g = 0; g < (lim || 6); g++) { const p = rxn(sma, [cur]); if (!p.length) break; const frs = byC(p[0].join('.').split('.')); if (frs[0] === cur) break; cur = frs[0]; extra.push(...frs.slice(1)); if (!labels.includes(label)) labels.push(label); }
  return { cur, extra, labels }; }
const ctxWith = (c, add2, flags) => { const st = Object.assign({}, c.st, { sp: new Map(c.st.sp), flags: new Set([...c.st.flags, ...(flags || [])]) }); for (const [k, v] of Object.entries(add2 || {})) { if (v === null) st.sp.delete(k); else st.sp.set(k, Object.assign({ sp: k }, v)); } return RP().stageCtx(st); };
const recipe = id => E.REC.find(x => x.id === id);
function wrap(id, fn) { const r = recipe(id); if (!r) return; const run0 = r.run; r.run = function (m, c, all, k) { let o; try { o = fn.call(this, m, c, all, k, run0.bind(this)); } catch (e) { o = undefined; } return o === undefined ? run0.call(this, m, c, all, k) : o; }; }
const frags = s => String(s || '').split('.');
const mainFrag = s => byC(frags(s))[0];
const raw = c => String(c.raw || '').toLowerCase();
const isAlkoxide = c => c.any(['NaOEt', 'NaOMe', 'NaOR', 'KOtBu']);
const alkoxR = c => c.has('NaOEt') ? 'CC' : c.has('NaOMe') ? 'C' : c.has('KOtBu') ? 'C(C)(C)C' : ((c.get('NaOR') || {}).smi || 'CC');
const XOF = c => c.has('NaOCl') ? ((c.get('NaOCl') || {}).X || 'Cl') : c.has('Br2/NaOH') || c.has('Br2') ? 'Br' : c.has('I2') ? 'I' : c.has('Cl2') ? 'Cl' : 'I';
/* graph distance inside one molecule (bonds) */
function dist(m, a) { const d = new Map([[a, 0]]); const q = [a]; while (q.length) { const x = q.shift(); for (const y of m.nbrs(x)) if (!d.has(y)) { d.set(y, d.get(x) + 1); q.push(y); } } return d; }

/* ---------------- organometallics ---------------- */
front({ id: 'deuterate', name: 'Organometallic / carbanion + D₂O', fam: 'Organometallic', canon: 'D₂O', db: [],
  test: c => c.any(['D2O', 'ODm']),
  run(m) { const s = smiOf(m); if (!s) return null; const r = applyAll(s, [['[#6:1][Mg,Li,Na,K,Zn,Cd,Cu]>>[#6:1][2H]', 'C–M → C–D'], ['[#6-:1]>>[#6+0:1][2H]', 'C⁻ → C–D']]); if (!r.labels.length) return null;
    return finish(m, r.cur, { title: 'R–M + D₂O → R–D', text: 'The very basic carbon takes D⁺ from D₂O: the metal leaves as M–OD. This is how a deuterium label is put on a chosen carbon.', cond: ['D₂O'], why: ['Grignard reagents, organolithiums and acetylides are strong bases: any O–H/O–D protonates (deuterates) them'], tags: ['acidbase'] }); } });
/* Mg in a protic solvent: the Grignard reagent is destroyed as it forms */
front({ id: 'grig_protic', name: 'Mg + R–X in a protic solvent', fam: 'Organometallic', canon: 'Mg, H₂O / EtOH', db: [],
  test: c => c.has('Mg') && !c.has('RMgX') && (c.any(['H2O', 'H3O+']) || (c.has('ROH') && !c.flag('ether')) || (c.flag('aq'))),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[#6][Cl,Br,I]').length) return null; const r = applyAll(s, [['[#6:1][Cl,Br,I]>>[#6:1]', 'R–X → R–H']], 1); if (!r.labels.length) return null;
    return finish(m, r.cur, { title: 'R–MgX is destroyed by the protic solvent → R–H', text: 'Mg does insert into C–X, but water/alcohol at once protonates the carbanion-like carbon: R–MgX + H₂O → R–H + Mg(OH)X. Grignard reagents need a dry ether (Et₂O or THF).', why: ['Use dry ether or THF for Grignard reagents'] }); } });
/* alkoxide + R–X: Williamson for methyl and unhindered 1° halides (E2 is the answer for 2°/3° or t-BuO⁻) */
front({ id: 'rx_alkoxide', name: 'Alkoxide + alkyl halide (Williamson)', fam: 'Alkyl halide', canon: 'NaOEt, R–X', db: ['Williamson ether synthesis'],
  test: c => isAlkoxide(c) && !c.any(['RX', 'KOH', 'alcKOH', 'Br2', 'NaCN']) && !c.heat,
  run(m, c) { const s = smiOf(m); if (!s) return null; const hx = M(m, '[CX4:1][Cl,Br,I:2]'); if (hx.length !== 1) return null; const [cc] = hx[0];
    const heavyC = m.nbrs(cc).filter(n => el(m, n) === 'C'); const bulky = c.has('KOtBu'); const methyl = heavyC.length === 0; if (!methyl && (heavyC.length > 1 || bulky)) return null;
    if (!methyl && M(m, '[CX4]([CX4]([CH3])([CH3])[CH3])[Cl,Br,I]').length) return null; // neopentyl
    const R = alkoxR(c); const p = rxn('[CX4:1][Cl,Br,I]>>[C:1]O' + R, [s]); if (!p.length) return null; const prod = mainFrag(p[0].join('.'));
    const minor = []; if (!methyl) { const e = rxn('[CX4;!H0:2][CX4:1][Cl,Br,I]>>[C:2]=[C:1]', [s]); if (e.length) minor.push(mainFrag(e[0].join('.'))); }
    return finish(m, prod, { title: 'Williamson: alkoxide + R–X → ether (SN2)', text: 'RO⁻ attacks the carbon from the back and X⁻ leaves (SN2). Methyl and 1° halides give the ether; with 2°/3° halides the alkoxide acts as a base and elimination wins.', minor, why: ['Choose the halide to be methyl/1° and put the bulky group on the alkoxide'], tags: ['SN2'] }); } });

/* ---------------- carbonyl partners ---------------- */
front({ id: 'carbonate_aldol', name: 'Aldol with K₂CO₃ / Na₂CO₃', fam: 'Carbonyl', canon: 'aq. K₂CO₃', db: ['Aldol condensation'],
  test: c => c.has('NaHCO3') && /co3|carbonate/.test(raw(c)),
  run(m, c, all, k) { if (!M(m, '[CX3]=O').length || M(m, '[CX3](=O)[OX2H1,OX1-]').length) return null; const al = recipe('aldol'); if (!al) return null;
    const c2 = ctxWith(c, { aqKOH: { dil: 1 }, KOH: {}, NaHCO3: null }, ['dil']); let o = null; try { o = al.run(m, c2, all, k); } catch (e) { o = null; } if (!o || !o.product || o.none) return null;
    return Object.assign({}, o, { why: ['K₂CO₃ is a weak base: enough for an aldol (especially with HCHO), too weak to dehydrate the product'].concat(o.why || []) }); } });
front({ id: 'gilman', name: 'Gilman reagent (R₂CuLi)', fam: 'Organometallic', canon: '(CH₃)₂CuLi', db: ['Corey–House synthesis'],
  test: c => c.has('R2CuLi'),
  run(m, c) { const s = smiOf(m); if (!s) return null; const R = (c.get('R2CuLi') || {}).smi || 'C'; let p, what, text;
    if ((p = rxn('[CX3:1](=[O:2])[Cl,Br,I]>>' + R + '[C:1]=[O:2]', [s])).length) { what = 'acyl chloride → ketone'; text = 'The cuprate replaces Cl by R and stops at the ketone (it is too mild to add to the ketone).'; }
    else if ((p = rxn('[CX3:1]=[CX3:2][CX3:3]=[O:4]>>' + R + '[C:1][C:2][C:3]=[O:4]', [s])).length) { what = '1,4-(conjugate) addition'; text = 'Soft R₂CuLi adds R to the β-carbon of the α,β-unsaturated carbonyl; the enolate is protonated in the work-up (Grignard reagents add 1,2 instead).'; }
    else if ((p = rxn('[#6:1][Cl,Br,I]>>' + R + '[#6:1]', [s])).length) { what = 'coupling (Corey–House)'; text = 'R₂CuLi couples R with the halide carbon: R–X + R′₂CuLi → R–R′. Works for 1° alkyl, vinyl and aryl halides.'; }
    else return null;
    return finish(m, mainFrag(p[0].join('.')), { title: `Gilman reagent: ${what}`, text, cond: ['R₂CuLi', 'ether'], tags: ['C–C'] }); } });
front({ id: 'stork_enamine', name: 'Enamine from a 2° amine', fam: 'Carbonyl', canon: 'pyrrolidine, H⁺', db: ['Stork enamine'],
  test: c => c.any(['pyrrolidine', 'R2NH']) && !c.any(['NaBH4', 'LiAlH4', 'H2', 'RX', 'HCN']),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[CX4;!H0][CX3](=O)[#6]').length || M(m, '[CX3](=O)[O,N]').length) return null;
    const N = c.has('pyrrolidine') ? 'N1CCCC1' : 'N(' + ((c.get('R2NH') || {}).smi || 'C') + ')' + ((c.get('R2NH') || {}).smi || 'C');
    const p = rxn('[CX4;!H0:1][CX3:2](=O)[#6:3]>>[C:1]=[C:2]([#6:3])' + N, [s]); if (!p.length) return null;
    return finish(m, byC(p.map(x => mainFrag(x.join('.'))))[0], { title: 'Ketone + 2° amine → enamine', text: 'The amine adds to C=O, water is lost to give an iminium ion, and an α-H is removed: C=C–NR₂ (enamine). Its α-carbon is nucleophilic (Stork).', cond: ['H⁺ (cat.)', '−H₂O'], tags: ['C=N'] }); } });
front({ id: 'enamine_alk', name: 'Enamine + alkyl / acyl halide (Stork)', fam: 'Carbonyl', canon: 'CH₃I then H₃O⁺', db: ['Stork enamine'],
  test: c => c.any(['RX', 'AcCl', 'RCOCl']),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[CX3]=[CX3][NX3;H0;!$(N-C=O)]').length) return null; const rx = c.get('RX'); const R = rx ? (rx.smi || 'C') : c.has('AcCl') ? 'C(C)=O' : 'C(=O)' + ((c.get('RCOCl') || {}).smi || 'C');
    const p = rxn('[CX3:1]=[CX3:2][NX3;H0:3]>>' + R + '[C:1][C:2]=[N+:3]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Stork: enamine C-alkylation → iminium salt', text: 'The enamine α-carbon attacks R–X (SN2); the product is an iminium salt. Aqueous acid then hydrolyses it back to the ketone, now alkylated at the α-carbon.', tags: ['SN2', 'C–C'] }); } });
add({ id: 'iminium_hyd', name: 'Iminium hydrolysis', fam: 'Amine', canon: 'H₃O⁺', db: [],
  test: c => c.any(['H3O+', 'H2O']),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('[CX3:1]=[NX3+]>>[C:1]=O', [s]); if (!p.length) return null; return finish(m, mainFrag(p[0].join('.')), { title: 'Iminium ion + H₂O → C=O (+ R₂NH₂⁺)', text: 'Water adds to C=N⁺, the amine leaves and the carbonyl group is back.' }); } });
front({ id: 'darzens', name: 'Darzens glycidic ester', fam: 'Carbonyl', canon: 'ClCH₂COOEt, NaOEt', db: [], multi: true,
  test: c => c.any(['NaOEt', 'NaOMe', 'NaNH2', 'KOtBu', 'NaH']),
  run(m) { const s = smiOf(m); if (!s) return null; const fr = frags(s); if (fr.length !== 2) return null; const hx = fr.find(x => M(fresh(x), '[Cl,Br][CX4;!H0][CX3](=O)[OX2][#6]').length); const co = fr.find(x => x !== hx && M(fresh(x), '[CX3;!$(C(=O)[O,N])]=O').length); if (!hx || !co) return null;
    const p = rxn('[#6:1][CX3:2](=[O:3])[#6:4].[Cl,Br][CX4;!H0:5][C:6](=[O:7])[O:8][#6:9]>>[#6:1][C:2]1([#6:4])[O:3][C:5]1[C:6](=[O:7])[O:8][#6:9]', [co, hx]).concat(rxn('[#6:1][CX3H1:2]=[O:3].[Cl,Br][CX4;!H0:5][C:6](=[O:7])[O:8][#6:9]>>[#6:1][C:2]1[O:3][C:5]1[C:6](=[O:7])[O:8][#6:9]', [co, hx])); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Darzens: α-halo ester + C=O → α,β-epoxy (glycidic) ester', text: 'The base makes the enolate of the α-halo ester; it adds to C=O (aldol-like) and the alkoxide displaces the halide intramolecularly, closing the epoxide.', tags: ['C–C'] }); } });
front({ id: 'barbituric', name: 'Malonic ester + urea → barbituric acid', fam: 'Ester', canon: 'urea, NaOEt', db: [], multi: true,
  test: c => c.any(['NaOEt', 'NaOMe']),
  run(m) { const s = smiOf(m); if (!s) return null; const fr = frags(s); const u = fr.find(x => canonS(x) === canonS('NC(N)=O')); const es = fr.find(x => M(fresh(x), '[CX4]([CX3](=O)[OX2][#6])[CX3](=O)[OX2][#6]').length); if (!u || !es) return null;
    const p = rxn('[CX4:1]([CX3](=O)[OX2][#6])[CX3](=O)[OX2][#6].[NH2]C(=O)[NH2]>>[C:1]1C(=O)NC(=O)NC1=O', [es, u]); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Malonic ester + urea → barbituric acid', text: 'Each N–H of urea attacks an ester C=O (base catalysed), ethanol leaves twice and the six-membered ring closes.', extra: ['CCO'], tags: ['acyl substitution'] }); } });

/* ---------------- cumene process ---------------- */
add({ id: 'cumene_ox', name: 'Cumene + O₂ → hydroperoxide', fam: 'Aromatic', canon: 'O₂, Δ', db: ['Phenol from cumene'],
  test: c => c.has('O2') && !c.any(['FeTube', 'Cu', 'KMnO4', 'K2Cr2O7', 'Ag2O']),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('[c:1][CX4H1:2]([CH3:3])[CH3:4]>>[c:1][C:2]([C:3])([C:4])OO', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Cumene + O₂ → cumene hydroperoxide', text: 'Air oxidation (radical chain) takes the weak 3° benzylic C–H; the radical adds O₂ and abstracts H: Ar–C(CH₃)₂–O–O–H.', cond: ['O₂ (air)', 'radical initiator'], tags: ['radical'] }); } });
add({ id: 'hock', name: 'Hydroperoxide + acid → phenol + ketone', fam: 'Aromatic', canon: 'H₃O⁺', db: ['Phenol from cumene'],
  test: c => c.any(['H3O+', 'H2SO4', 'HCl']),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('[c:1][CX4:2]([#6:3])([#6:4])[OX2][OX2H1]>>[c:1][OH].[C:2](=O)([#6:3])[#6:4]', [s]); if (!p.length) return null; const fr = p[0];
    const ph = fr.find(x => /c/.test(x)) || fr[0]; const k2 = fr.find(x => x !== ph);
    return finish(m, ph, { title: 'Cumene hydroperoxide + H₃O⁺ → phenol + acetone', text: 'Acid protonates the O–OH; water leaves as the aryl group migrates from C to O (1,2-shift onto electron-poor O). Water adds to the oxocarbenium ion, which falls apart into phenol and the ketone.', extra: k2 ? [k2] : [], tags: ['rearrangement'] }); } });
/* calcium / barium carboxylates, dry distillation → ketones (formate gives an aldehyde) */
front({ id: 'ca_salt', name: 'Dry distillation of calcium carboxylates', fam: 'Acid derivative', canon: '(RCOO)₂Ca, Δ', db: ['Ketones from calcium salts'], multi: true,
  test: c => c.heat || c.flag('dry'),
  run(m) { const s = smiOf(m); if (!s || !/Ca|Ba/.test(s)) return null; const acs = [...new Set(frags(s).filter(x => /\[O-\]/.test(x) && M(fresh(x), '[CX3](=O)[OX1-]').length))]; if (!acs.length) return null;
    const fm = acs.find(x => canonS(x) === canonS('O=C[O-]')); const others = acs.filter(x => x !== fm); let prod = null, minor = [], text;
    if (others.length >= 2) { prod = rxn('[#6:1][CX3:2](=[O:3])[OX1-].[#6:4][CX3](=O)[OX1-]>>[#6:1][C:2](=[O:3])[#6:4]', [others[0], others[1]]); minor = [rxn('[#6:1][CX3:2](=[O:3])[OX1-].[#6:4][CX3](=O)[OX1-]>>[#6:1][C:2](=[O:3])[#6:4]', [others[0], others[0]]), rxn('[#6:1][CX3:2](=[O:3])[OX1-].[#6:4][CX3](=O)[OX1-]>>[#6:1][C:2](=[O:3])[#6:4]', [others[1], others[1]])].filter(x => x.length).map(x => x[0][0]); text = 'A mixture of two calcium salts gives the mixed ketone plus both symmetrical ketones.'; }
    else if (others.length === 1 && fm) { prod = rxn('[#6:1][CX3:2](=[O:3])[OX1-].[CX3H1](=O)[OX1-]>>[#6:1][CH1:2]=[O:3]', [others[0], fm]); text = 'With calcium formate one acyl group gets the formate H: the product is an aldehyde.'; }
    else if (others.length === 1) { prod = rxn('[#6:1][CX3:2](=[O:3])[OX1-].[#6:4][CX3](=O)[OX1-]>>[#6:1][C:2](=[O:3])[#6:4]', [others[0], others[0]]); text = '(RCOO)₂Ca → R–CO–R + CaCO₃ on strong heating.'; }
    else if (fm) { return finish(m, 'C=O', { title: 'Calcium formate, Δ → HCHO', text: '(HCOO)₂Ca → HCHO + CaCO₃.', extra: ['O=C([O-])[O-]'] }); }
    if (!prod || !prod.length) return null;
    return finish(m, prod[0][0], { title: 'Dry distillation of calcium salt → ketone (+ CaCO₃)', text, minor, extra: ['O=C([O-])[O-]'], tags: ['decarboxylation'] }); } });

/* ---------------- every site ---------------- */
front({ id: 'haloform_all', name: 'Haloform on every methyl ketone', fam: 'Carbonyl', canon: 'X₂, NaOH', db: ['Haloform reaction'],
  test: c => c.has('NaOCl') || c.has('Br2/NaOH') || (c.any(['I2', 'Br2', 'Cl2']) && c.any(['KOH', 'aqKOH'])),
  run(m, c) { if (M(m, '[CH3][CX3](=O)[#6]').length < 2) return null; const s = smiOf(m); const r = applyAll(s, [['[CH3][C:1](=[O:2])[#6:3]>>[O-][C:1](=[O:2])[#6:3]', 'haloform']], 6); if (!r.labels.length) return null;
    const X = XOF(c); const n = M(m, '[CH3][CX3](=O)[#6]').length;
    return finish(m, r.cur, { title: `Haloform on all ${n} CH₃CO groups → carboxylate + ${n} CHX₃`, text: 'Each CH₃–C=O is trihalogenated and cleaved by OH⁻: the ring keeps a –COO⁻ where each acetyl group was.', extra: Array(n).fill(X + 'C(' + X + ')' + X), tags: ['haloform'] }); } });
wrap('decarbox', function (m, c, all, k, run0) { if (!c.has('sodalime')) return undefined; let o = run0(m, c, all, k); if (!o || !o.product || o.none) return o; let p = o.product, n = 1;
  while (n < 6 && M(p, '[CX3](=O)[OX2H1,OX1-]').length) { let o2 = null; try { o2 = run0(p, c, all, k); } catch (e) { o2 = null; } if (!o2 || !o2.product) break; o.steps = (o.steps || []).concat(o2.steps || []); p = o2.product; n++; }
  if (n > 1) { o = Object.assign({}, o, { product: p, title: o.title + ` (all ${n} –COOH groups)` }); } return o; });
/* malonic / substituted esters: hydrolyse every ester, then decarboxylate. The arrowed r6 recipe goes first when it can. */
wrap('malonic_hyd', function (m, c, all, k, run0) { const kd = recipe('ketoester_decarb'); if (kd) { let ok = false; try { ok = kd.test(c, all, k); } catch (e) { ok = false; } if (ok) { let o0 = null; try { o0 = kd.run(m, c, all, k); } catch (e) { o0 = null; } if (o0 && o0.product && !o0.none && !M(o0.product, '[CX3](=O)[OX2][#6]').length) return o0; } }
  const o = run0(m, c, all, k); if (!o || !o.product || o.none) return o; const s = smiOf(o.product); if (!s || !M(o.product, '[CX3](=O)[OX2][#6]').length) return o;
  let r = applyAll(s, [['[C:1](=[O:2])[OX2][CX4]>>[C:1](=[O:2])O', 'ester hydrolysis']], 6); let cur = r.cur;
  const d = rxn('[CX3:1](=[O:5])([OX2H1])[CX4;!H0:2][CX3](=O)[OX2H1]>>[C:1](=[O:5])([OX2H1])[C:2]', [cur]); if (d.length) cur = mainFrag(d[0].join('.'));
  const R = fresh(cur); if (!R) return o; return Object.assign({}, o, { product: R, steps: (o.steps || []).concat([step(o.product, R, 'The remaining ester groups are hydrolysed too.', [])]) }); });
front({ id: 'malonic_hyd2', name: 'Substituted malonic ester: hydrolysis + decarboxylation', fam: 'Ester', canon: 'H₃O⁺, Δ', db: ['Malonic ester synthesis'],
  test: c => c.heat && c.any(['H3O+', 'H2SO4', 'HCl', 'KOH', 'aqKOH']),
  run(m) { const s = smiOf(m); if (!s || M(m, '[CX3](=O)[OX2][#6]').length < 1) return null; if (!M(m, '[CX4H0]([CX3](=O)[OX2,OX1-,OX2H1])[CX3](=O)[OX2,OX1-]').length && M(m, '[CX3](=O)[OX2][#6]').length < 3) return null;
    let cur = applyAll(s, [['[C:1](=[O:2])[OX2][CX4]>>[C:1](=[O:2])O', 'ester hydrolysis']], 8).cur;
    for (let g = 0; g < 3; g++) { const d = rxn('[CX3:1](=[O:5])([OX2H1])[CX4:2][CX3](=O)[OX2H1]>>[C:1](=[O:5])([OX2H1])[C:2]', [cur]).concat(rxn('[CX3:1](=[O:5])([#6])[CX4:2][CX3](=O)[OX2H1]>>[C:1](=[O:5])([#6])[C:2]', [cur])); if (!d.length) break; cur = mainFrag(d[0].join('.')); }
    if (cur === s) return null;
    return finish(m, cur, { title: 'Hydrolysis of every ester, then the gem-dicarboxylic (malonic) acid loses CO₂', text: 'Hot aqueous acid hydrolyses all the esters; a carbon carrying two –COOH groups loses one as CO₂ through a six-membered cyclic transition state.', extra: ['O=C=O', 'CCO'], tags: ['decarboxylation'] }); } });
/* peptides / esters + amides: hot aqueous acid cuts every C(=O)–N and C(=O)–O–C bond */
front({ id: 'full_hyd', name: 'Complete hydrolysis (peptides, esters)', fam: 'Biomolecule', canon: 'H₃O⁺, Δ', db: ['Hydrolysis of peptides'],
  test: c => c.heat && c.any(['H3O+', 'H2SO4', 'HCl']) && !c.any(['NaNO2', 'HNO2', 'KMnO4', 'K2Cr2O7', 'HNO3']),
  run(m) { const s = smiOf(m); if (!s) return null; const am = M(m, '[CX3](=O)!@[NX3;!$(N-N);!$(N=*)]').length, es = M(m, '[CX3](=O)!@[OX2][#6]').length; if (!am || am + es < 2) return null;
    let pool = frags(s), changed = true, g = 0; const T = ['[C:1](=[O:2])!@[NX3:3]>>[C:1](=[O:2])O.[N:3]', '[C:1](=[O:2])!@[OX2:3][#6:4]>>[C:1](=[O:2])O.[O:3][#6:4]'];
    while (changed && g++ < 12) { changed = false; const next = []; for (const f of pool) { let hit = null; for (const t of T) { const p = rxn(t, [f]); if (p.length) { hit = p[0]; break; } } if (hit) { next.push(...hit.join('.').split('.')); changed = true; } else next.push(f); } pool = next; }
    const uniq = [...new Set(pool.map(canonS).filter(Boolean))]; const main = byC(uniq)[0];
    return finish(m, main, { title: 'Complete acid hydrolysis: every amide (peptide) and ester bond is cut', text: 'Hot aqueous acid hydrolyses each C(=O)–NH (peptide) and C(=O)–O–C bond: the amino acids (and alcohols) are set free.', extra: uniq.filter(x => x !== main), tags: ['hydrolysis'] }); } });
/* glycosides / disaccharides: dilute acid cuts the glycosidic C–O, each ring opens to the open-chain sugar */
function openHemi(smi) { const t = fresh(smi); if (!t) return smi; for (let g = 0; g < 4; g++) { const h = M(t, '[CX4;R:1]([OX2H1:2])[OX2;R:3]')[0]; if (!h) break; const [c1, o5, o2] = h; const b = t.bondBetween(c1, o2); if (!b) break; t.removeBond(b.id); const b2 = t.bondBetween(c1, o5); if (b2) b2.o = 2; [c1, o5, o2].forEach(id => { const a = t.atoms.get(id); if (a) a.hFix = null; }); } return smiOf(t) || smi; }
front({ id: 'glycoside_hyd', name: 'Hydrolysis of glycosides / disaccharides', fam: 'Biomolecule', canon: 'dil. H₂SO₄ / H₃O⁺', db: ['Hydrolysis of sucrose'],
  test: c => c.any(['H3O+', 'H2SO4', 'HCl']) && (c.heat || c.has('H2O') || c.flag('dil')),
  run(m) { const s = smiOf(m); if (!s || M(m, '[CX4][OX2H1]').length < 3 || !M(m, '[CX4;R]([OX2;R])[OX2;!R][#6]').length) return null;
    let pool = frags(s), g = 0; while (g++ < 6) { let hit = false; const next = []; for (const f of pool) { const p = rxn('[C;R:1]([O;R:2])[O;!R:3][#6:4]>>[C:1]([O:2])O.[O:3][#6:4]', [f]); if (p.length) { next.push(...p[0].join('.').split('.')); hit = true; } else next.push(f); } pool = next; if (!hit) break; }
    const opened = [...new Set(pool.map(openHemi).map(canonS).filter(Boolean))]; const main = byC(opened)[0];
    return finish(m, main, { title: 'Glycosidic bond hydrolysed → monosaccharides (open-chain forms)', text: 'Aqueous acid protonates the glycosidic oxygen; the anomeric C–O bond breaks (oxocarbenium ion) and water adds. Each cyclic hemiacetal is in equilibrium with its open-chain aldose/ketose (sucrose → glucose + fructose, “invert sugar”).', extra: opened.filter(x => x !== main), tags: ['hydrolysis'] }); } });
/* γ/δ-hydroxy acids close to lactones in acid */
add({ id: 'lactonize', name: 'Hydroxy acid → lactone', fam: 'Acid derivative', canon: 'H₃O⁺', db: [],
  test: c => c.any(['H3O+', 'H2SO4']),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('[OX2H1:1][CX4:2][#6:3]~[#6:4][CX3:5](=[O:6])[OX2H1,OX1-]>>[O:1]1[C:2][#6:3]~[#6:4][C:5]1=[O:6]', [s]).concat(rxn('[OX2H1:1][CX4:2][#6:3]~[#6:7]~[#6:4][CX3:5](=[O:6])[OX2H1,OX1-]>>[O:1]1[C:2][#6:3]~[#6:7]~[#6:4][C:5]1=[O:6]', [s])); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'γ/δ-Hydroxy acid → lactone (intramolecular ester)', text: 'In acid the OH attacks its own –COOH: a five- (γ) or six-membered (δ) ring forms and water is lost. Five- and six-membered rings close easily.', tags: ['esterification'] }); } });
wrap('grignard_add', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || o.none) return o; const s = smiOf(o.product); if (!s || !M(o.product, '[CX4][Cl,Br,I]').length) return o;
  const T = ['[OX2H1,OX1-:1][C:2][C:3][C:4][C:5][Cl,Br,I]>>[O:1]1[C:2][C:3][C:4][C:5]1', '[OX2H1,OX1-:1][C:2][C:3][C:4][C:5][C:6][Cl,Br,I]>>[O:1]1[C:2][C:3][C:4][C:5][C:6]1', '[O:1]([Mg])[C:2][C:3][C:4][C:5][Cl,Br,I]>>[O:1]1[C:2][C:3][C:4][C:5]1', '[O:1]([Mg])[C:2][C:3][C:4][C:5][C:6][Cl,Br,I]>>[O:1]1[C:2][C:3][C:4][C:5][C:6]1'];
  for (const t of T) { const p = rxn(t, [s]); if (p.length) { const R = fresh(mainFrag(p[0].join('.'))); if (!R) break; return Object.assign({}, o, { product: R, title: o.title + ', then the alkoxide closes the ring (intramolecular SN2)', steps: (o.steps || []).concat([step(o.product, R, 'The new alkoxide O⁻ is 4–5 bonds from C–X: it displaces the halide inside the molecule and a cyclic ether (THF/THP ring) forms.', [])]) }); } } return o; });
/* active methylene with its own C–X: base closes the ring (malonic / acetoacetic) */
front({ id: 'am_ring', name: 'Intramolecular C-alkylation of an active methylene', fam: 'Ester', canon: 'NaOEt', db: ['Malonic ester synthesis'],
  test: c => c.any(['NaOEt', 'NaOMe', 'NaH', 'KOH', 'NaOR', 'KOtBu']),
  run(m) { const am = M(m, '[CX4;!H0:1]([$([CX3]=O),$(C#N)])[$([CX3]=O),$(C#N)]'); const hx = M(m, '[CX4:1][Cl,Br,I:2]'); if (!am.length || !hx.length) return null;
    for (const [a] of am) { const d = dist(m, a); for (const [b, x] of hx) { const L = d.get(b); if (L == null || L < 2 || L > 5) continue; const t = m.clone(); t.removeAtom(x); t.addBond(a, b, 1); [a, b].forEach(id => { const q = t.atoms.get(id); if (q) q.hFix = null; }); const smi = smiOf(t); if (!smi) continue;
      return finish(m, mainFrag(smi), { title: `Base closes a ${L + 1}-membered ring (intramolecular C-alkylation)`, text: 'The base removes the acidic H between the two C=O groups; the enolate carbon attacks the C–X of its own chain (SN2) and a carbocycle forms.', tags: ['SN2', 'C–C'] }); } }
    return null; } });

/* ---------------- halogenation ---------------- */
front({ id: 'alpha_x_ket', name: 'α-Halogenation of ketones (base or NBS)', fam: 'Carbonyl', canon: 'Br₂, NaOH / NBS', db: [],
  test: c => (c.any(['Br2', 'Cl2', 'I2', 'Br2/NaOH']) && c.any(['KOH', 'aqKOH', 'Br2/NaOH'])) || (c.has('NBS') && !c.flag('light') && !c.has('ROOR')),
  run(m, c) { const ks = M(m, '[CX4;!H0:1][CX3:2](=O)[#6]'); if (!ks.length || M(m, '[CX3](=O)[O,N]').length || M(m, '[CX3]=[CX3]').length) return null;
    if (!c.has('NBS') && M(m, '[CH3][CX3](=O)').length) return null; // methyl ketones: haloform
    const X = c.has('NBS') ? 'Br' : XOF(c); const alphas = [...new Set(ks.map(x => x[0]))]; const pick = alphas.sort((p, q) => c.has('NBS') ? hTot(m, p) - hTot(m, q) : hTot(m, q) - hTot(m, p))[0];
    const t = m.clone(); const xi = E.addOn ? E.addOn(t, pick, X) : null; if (xi == null) return null; const smi = smiOf(t); if (!smi) return null;
    return finish(m, smi, { title: c.has('NBS') ? 'NBS: α-bromination of the ketone' : 'Base-promoted α-halogenation (not a methyl ketone, so no haloform)', text: 'The enol/enolate attacks X–X (or NBS); X goes on the α-carbon. With OH⁻ each new X makes the remaining α-H more acidic, so excess reagent gives polyhalogenation at the same carbon.', tags: ['α-substitution'] }); } });
const DON = [['[NX3;H2,H1;!$(NC=O)]', 5], ['[OX2H1]', 4], ['[OX1-]', 6], ['[NX3;H0;!$(NC=O)]([#6])[#6]', 5], ['[NX3;H1]C(=O)', 3], ['[OX2][CX4]', 3], ['[OX2]C(=O)', 2]];
front({ id: 'br2_act', name: 'Br₂ (no catalyst) on an activated ring', fam: 'Aromatic', canon: 'Br₂ in CCl₄ / CS₂ / AcOH', db: ['Bromination of phenol'],
  test: c => c.has('Br2') && !c.any(['FeX3', 'HOBr', 'H2O', 'KOH', 'aqKOH', 'NaOCl', 'Br2/NaOH', 'redP', 'HVZ', 'NBS']) && !c.light && !c.flag('aq'),
  run(m) { const AR = E.aromSet(m); if (!AR.size || M(m, '[CX3]=[CX3]').length || M(m, '[CX2]#[CX2]').length) return null;
    let best = null; for (const [sm, sc] of DON) for (const mt of M(m, '[c:1]' + sm)) { const ca = mt[0]; if (!best || sc > best.sc) best = { sc, ca }; }
    if (!best) return null; const ring = C.rings(m, 6).find(r => r.length === 6 && r.includes(best.ca) && r.every(x => AR.has(x))); if (!ring) return null;
    const i0 = ring.indexOf(best.ca); const at = k => ring[(i0 + k + 6) % 6]; const cand = [at(3), at(1), at(-1)].filter(x => hTot(m, x) > 0 && el(m, x) === 'C'); if (!cand.length) return null;
    const pos = cand[0]; const t = m.clone(); E.addOn(t, pos, 'Br'); const smi = smiOf(t); if (!smi) return null;
    return finish(m, smi, { title: `Br₂ without a Lewis acid: the activated ring is brominated ${pos === at(3) ? 'para' : 'ortho'} to the donor`, text: 'An –OH/–OR/–NHCOR/–OCOR group pushes electron density into its ring, so Br₂ reacts without FeBr₃. In a non-polar solvent (CCl₄, CS₂) only one Br goes on, para (ortho when para is taken). In water, phenol and aniline give the 2,4,6-tribromo product.', extra: ['Br'], tags: ['SEAr'] }); } });
front({ id: 'allylic_x', name: 'Allylic halogenation (Cl₂ at high T / hν)', fam: 'Alkene', canon: 'Cl₂, 773 K', db: ['Allylic substitution'],
  test: c => c.any(['Cl2', 'Br2']) && (c.light || (c.T != null && c.T >= 600)) && !c.any(['FeX3', 'H2O']),
  run(m, c) { const s = smiOf(m); if (!s || E.aromSet(m).size || !M(m, '[CX3]=[CX3][CX4;!H0]').length) return null; const X = c.has('Cl2') ? 'Cl' : 'Br'; const p = rxn('[C:1]=[C:2][CX4;!H0:3]>>[C:1]=[C:2][C:3]' + X, [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: `Allylic substitution: ${X} replaces an allylic H`, text: 'At high temperature (or in light, low X₂ concentration) X• abstracts the weak allylic H; the resonance-stabilised allyl radical takes X from X₂. Addition to C=C is reversible at this temperature, so substitution wins.', extra: ['H' + X], tags: ['radical'] }); } });
front({ id: 'so2cl2', name: 'SO₂Cl₂ free-radical chlorination', fam: 'Radical', canon: 'SO₂Cl₂, hν / peroxide', db: [],
  test: c => c.has('SO2Cl2'),
  run(m, c, all, k) { const rx = recipe('radx'); if (!rx) return null; const c2 = ctxWith(c, { Cl2: {}, SO2Cl2: null }, ['light']); let o = null; try { o = rx.run(m, c2, all, k); } catch (e) { o = null; } if (!o || !o.product || o.none) return null;
    return Object.assign({}, o, { title: 'SO₂Cl₂: free-radical chlorination (side chain / C–H)', why: ['SO₂Cl₂ is a liquid Cl• source (SO₂ + Cl₂ in the flask); toluene gives benzyl chloride'].concat(o.why || []) }); } });
/* ---------------- reductions / oxidations ---------------- */
front({ id: 'rosen_lind', name: 'Rosenmund + Lindlar (both groups)', fam: 'Acid derivative', canon: 'H₂, Pd/BaSO₄', db: ['Rosenmund reduction'],
  test: c => c.has('H2') && (c.has('Lindlar') || /baso4|caco3|quinoline/.test(raw(c))),
  run(m) { if (!M(m, '[CX3](=O)Cl').length || !M(m, '[CX2]#[CX2]').length) return null; const r = applyAll(smiOf(m), [['[CX3:1](=[O:2])Cl>>[CH1:1]=[O:2]', 'COCl → CHO'], ['[C:1]#[C:2]>>[C:1]=[C:2]', 'C≡C → C=C']], 4);
    return finish(m, r.cur, { title: 'Poisoned Pd + H₂: –COCl → –CHO and C≡C → cis-C=C', text: 'The poisoned catalyst does both partial reductions: the acyl chloride stops at the aldehyde (Rosenmund) and the alkyne stops at the cis-alkene (Lindlar).', tags: ['red'] }); } });
front({ id: 'alkyne_kmno4', name: 'Alkyne + KMnO₄', fam: 'Alkyne', canon: 'KMnO₄', db: ["Baeyer's test"],
  test: c => c.has('KMnO4'),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[CX2]#[CX2]').length) return null; const mild = (c.flag('cold') || c.flag('dil') || c.flag('mild')) && !c.heat && !c.any(['H3O+', 'H2SO4', 'KOH']);
    let p; if (mild) p = rxn('[#6:1][C:2]#[C:3][#6:4]>>[#6:1][C:2](=O)[C:3](=O)[#6:4]', [s]); else { p = rxn('[#6:1][C:2]#[C:3][#6:4]>>[#6:1][C:2](=O)O.[#6:4][C:3](=O)O', [s]); if (!p.length) p = rxn('[#6:1][C:2]#[CH1]>>[#6:1][C:2](=O)O.O=C=O', [s]); if (!p.length && canonS(s) === 'C#C') p = [['O=C(O)C(=O)O']]; }
    if (!p || !p.length) return null; const fr = byC(p[0].join('.').split('.'));
    return finish(m, fr[0], { title: mild ? 'Cold dilute KMnO₄: alkyne → 1,2-dione (purple colour lost)' : 'KMnO₄ cleaves the C≡C → carboxylic acid(s) (Baeyer’s test positive)', text: 'Permanganate adds across the triple bond; under hot/acidic/alkaline conditions the C≡C is split and each carbon ends up as –COOH (a terminal ≡CH ends as CO₂). The purple colour disappears (brown MnO₂).', extra: fr.slice(1), tags: ['ox'] }); } });
wrap('ozone', function (m, c, all, k, run0) { const st = all || [c]; const red = st.some(x => x.any && x.any(['Zn', 'Me2S', 'PPh3', 'NaBH4', 'LiAlH4'])); const wet = st.some(x => x.any && x.any(['H2O', 'H3O+'])); if (red || !wet || st.some(x => x.has && x.has('H2O2'))) return undefined;
  const c2 = ctxWith(c, { H2O2: {} }); const o = run0(m, c2, all, k); if (o && o.product) o.why = ['No Zn/Me₂S: the H₂O₂ formed in the hydrolysis oxidises –CHO to –COOH (oxidative work-up)'].concat(o.why || []); return o; });
wrap('hydration', function (m, c, all, k, run0) { if (c.has('Hg2+')) return undefined; if (M(m, '[CX2]#[CX2]').length && !M(m, '[CX3]=[CX3]').length) return null; return undefined; });
front({ id: 'hemiacetal_test', name: 'Hemiacetal + Tollens / Fehling', fam: 'Carbonyl', canon: "Tollens'", db: [],
  test: c => c.any(['Tollens', 'Fehling']),
  run(m, c) { if (!M(m, '[CX4H1;!R]([OX2H1])[OX2][#6]').length || M(m, '[CX3H1]=O').length) return null;
    return positive(c.has('Tollens') ? 'Silver mirror: the hemiacetal opens to the aldehyde' : 'Red Cu₂O: the hemiacetal opens to the aldehyde', 'A hemiacetal is in equilibrium with its aldehyde + alcohol, so it reduces Tollens’/Fehling’s reagent like the aldehyde does.'); } });
front({ id: 'dehydro_arom', name: 'Dehydrogenation → aromatic ring', fam: 'Alkene', canon: 'Se or Pd/C, Δ', db: [],
  test: c => c.has('aromat') || (c.has('Pd') && c.heat && !c.has('H2')) || (c.has('S8') && c.heat),
  run(m) { const AR = E.aromSet(m); const ring = C.rings(m, 6).find(r => r.length === 6 && r.every(x => el(m, x) === 'C' && !AR.has(x)) && r.every(x => m.nbrs(x).filter(y => !r.includes(y) && el(m, y) !== 'H').length <= 1) && r.every(x => m.bondsOf(x).every(b => r.includes(m.other(b, x)) || b.o === 1))); if (!ring) return null;
    const t = m.clone(); for (let i = 0; i < 6; i++) { const b = t.bondBetween(ring[i], ring[(i + 1) % 6]); if (b) b.o = i % 2 === 0 ? 2 : 1; } ring.forEach(id => { t.atoms.get(id).hFix = null; }); const smi = smiOf(t); if (!smi || !C.isValid(t)) return null;
    return finish(m, smi, { title: 'Dehydrogenation (Se / Pd–C / S, Δ) → benzene ring', text: 'Heating with Se (or Pd/C, S) removes hydrogen from a six-membered carbocycle; the ring becomes aromatic (H₂Se / H₂ / H₂S given off).', tags: ['ox'] }); } });

/* ---------------- tests & small corrections ---------------- */
front({ id: 'agno3_rx', name: 'Alcoholic AgNO₃ + halide', fam: 'Alkyl halide', canon: 'alc. AgNO₃', db: [],
  test: c => c.has('AgNO3') && !c.any(['KOH', 'NH3', 'Tollens', 'aqKOH', 'NaNO2']),
  run(m) { const s = smiOf(m); if (!s) return null; if (/\[(Cl|Br|I)-\]/.test(s)) { const X = s.match(/\[(Cl|Br|I)-\]/)[1]; return positive(`White/yellow precipitate of Ag${X} at once`, `The halide here is an ion (X⁻), so Ag⁺ precipitates it immediately.`); }
    const sp3 = M(m, '[CX4:1][Cl,Br,I:2]'); if (!sp3.length) { if (M(m, '[#6][Cl,Br,I]').length) return none('No precipitate', 'Vinyl and aryl halides (C–X on sp² carbon, partial double bond) do not ionise, so no AgX precipitate.'); return null; }
    const cc = sp3[0][0]; if (m.nbrs(cc).filter(n => ['Cl', 'Br', 'I'].includes(el(m, n))).length >= 3) return none('No precipitate', 'CHCl₃/CCl₄: the C–Cl bonds are not ionised (no stable carbocation), so no AgCl forms.');
    const p = rxn('[CX4:1][Cl,Br,I]>>[C:1]O[N+](=O)[O-]', [s]); if (!p.length) return null; const deg = m.nbrs(cc).filter(n => el(m, n) === 'C').length; const fast = deg >= 3 || M(m, '[CX4]([Cl,Br,I])[$(C=C),$(c)]').length;
    return finish(m, mainFrag(p[0].join('.')), { title: `AgX precipitate${fast ? ' at once' : ' on warming'} + alkyl nitrate`, text: 'Ag⁺ pulls X⁻ off (SN1-like): AgX precipitates and NO₃⁻ bonds to the carbon (R–ONO₂). 3°, allylic and benzylic halides react fastest.', extra: ['[Ag+]'], tags: ['test'] }); } });
front({ id: 'ald_nh3', name: 'Aliphatic aldehyde + NH₃', fam: 'Carbonyl', canon: 'NH₃', db: [],
  test: c => c.has('NH3') && !c.any(['H2', 'Pd', 'NaBH4', 'HCN', 'LiAlH4', 'Na', 'RX', 'NaOCl', 'Br2/NaOH']) && !c.heat,
  run(m, c, all, k) { const s = canonS(smiOf(m)); if (s === 'C=O') return finish(m, 'C1N2CN3CN1CN(C2)C3', { title: 'HCHO + NH₃ → urotropine (hexamethylenetetramine)', text: '6 HCHO + 4 NH₃ → (CH₂)₆N₄ + 6 H₂O: the imine CH₂=NH trimerises and condenses into the cage.', tags: ['C=N'] });
    if (!M(m, '[CX3H1](=O)[CX4]').length || E.aromSet(m).size) return null; const r = recipe('ald_ammonia'); if (!r) return null; let o = null; try { o = r.run(m, c, all, k); } catch (e) { o = null; } return o && o.product ? o : null; } });
add({ id: 'acylx_hyd', name: 'Acyl halide + water', fam: 'Acid derivative', canon: 'H₂O', db: [],
  test: c => c.any(['H2O', 'H3O+']),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('[CX3:1](=[O:2])[Cl,Br,I]>>[C:1](=[O:2])O', [s]); if (!p.length) return null; const X = (M(m, '[CX3](=O)[Cl,Br,I]')[0] || []).map(i => el(m, i)).filter(e => ['Cl', 'Br', 'I'].includes(e))[0] || 'Cl'; return finish(m, mainFrag(p[0].join('.')), { title: 'Acyl halide + H₂O → carboxylic acid', text: 'Water attacks the very electrophilic C=O; X⁻ leaves (addition–elimination).', extra: [X], tags: ['acyl substitution'] }); } });
front({ id: 'na_roh_rx', name: 'Na + alcohol + alkyl halide → ether', fam: 'Alcohol', canon: 'Na then R–X', db: ['Williamson ether synthesis'], multi: true,
  test: c => c.any(['Na', 'K', 'NaH']),
  run(m) { const s = smiOf(m); if (!s) return null; const fr = frags(s); if (fr.length !== 2) return null; const oh = fr.find(x => M(fresh(x), '[CX4][OX2H1]').length && !M(fresh(x), '[Cl,Br,I]').length); const rx = fr.find(x => x !== oh && M(fresh(x), '[CX4;!$(C([#6])([#6])[#6])][Cl,Br,I]').length && !M(fresh(x), '[OX2H1]').length); if (!oh || !rx) return null;
    const p = rxn('[CX4:1][OX2H1:2].[CX4:3][Cl,Br,I]>>[C:1][O:2][C:3]', [oh, rx]); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Na makes the alkoxide from the alcohol, which then does SN2 on R–X (Williamson)', text: 'Sodium reacts with the O–H first (H₂↑); the alkoxide RO⁻Na⁺ attacks the alkyl halide. (Wurtz coupling needs the halide alone with Na.)', tags: ['SN2'] }); } });
front({ id: 'carbylamine_neg', name: 'Carbylamine test on amides / 2°–3° amines', fam: 'Amine', canon: 'CHCl₃, KOH', db: ['Carbylamine test'],
  test: c => c.has('CHCl3') && c.any(['KOH', 'alcKOH', 'aqKOH']),
  run(m) { if (M(m, '[NX3;H2;!$(N[CX3]=O);!$(N[SX4])][#6]').length || M(m, 'c[OX2H1,OX1-]').length || M(m, '[CX3]=[CX3]').length) return null; if (!M(m, '[CX3](=O)[NX3]').length && !M(m, '[NX3;H1,H0;!$(N[CX3]=O)]([#6])[#6]').length) return null;
    return none('Carbylamine test negative', 'Only 1° amines (–NH₂ on carbon) give an isocyanide with CHCl₃/KOH. Amides and 2°/3° amines do not.'); } });
front({ id: 'snar_nu', name: 'SNAr on an activated aryl halide', fam: 'Aromatic', canon: 'Nu⁻, Δ (o/p-NO₂)', db: ['Nucleophilic aromatic substitution'],
  test: c => c.any(['R2NH', 'RNH2', 'NH3', 'pyrrolidine', 'NaOMe', 'NaOEt', 'NaOR', 'NaSH', 'RS-']),
  run(m, c) { const s = smiOf(m); if (!s) return null; const act = '[c;$(c1ccc(cc1)[N+](=O)[O-]),$(c(:c[N+](=O)[O-])):1][F,Cl,Br,I]';
    const nu = c.has('R2NH') ? 'N(' + (c.get('R2NH').smi || 'C') + ')' + (c.get('R2NH').smi || 'C') : c.has('pyrrolidine') ? 'N1CCCC1' : c.has('RNH2') ? 'N' + (c.get('RNH2').smi || 'C') : c.has('NH3') ? 'N' : c.has('NaOMe') ? 'OC' : c.has('NaOEt') ? 'OCC' : c.has('NaOR') ? 'O' + (c.get('NaOR').smi || 'C') : c.has('NaSH') ? 'S' : 'S' + ((c.get('RS-') || {}).nu === 'SR' ? 'C' : 'C');
    const p = rxn(act + '>>[c:1]' + nu, [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'SNAr: the nucleophile replaces the halide ortho/para to NO₂', text: 'The nucleophile adds to the ring carbon carrying X (Meisenheimer complex, negative charge on the nitro group), then X⁻ leaves. F is the best leaving group here because the addition step is rate-limiting.', tags: ['SNAr'] }); } });
front({ id: 'acid_salt', name: 'Stronger acid + salt of a weaker acid', fam: 'Acid–base', canon: 'CH₃COONa', db: [],
  test: c => c.has('NaOAc'),
  run(m) { const so = M(m, '[SX4](=O)(=O)[OX2H1]'); if (!so.length) return null; const s = smiOf(m); const p = rxn('[SX4:1](=[O:2])(=[O:3])[OX2H1]>>[S:1](=[O:2])(=[O:3])[O-]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Sulphonic acid + CH₃COO⁻ → sulphonate + CH₃COOH', text: 'ArSO₃H (pKa ≈ −2) is far stronger than CH₃COOH (4.8), so it protonates acetate: the equilibrium lies on the side of the weaker acid.', extra: ['CC(=O)O'], tags: ['acidbase'] }); } });
front({ id: 'fc_ccl4', name: 'Benzene + CCl₄, AlCl₃', fam: 'Aromatic', canon: 'CCl₄, AlCl₃', db: [],
  test: c => c.has('FeX3') && /ccl4|tetrachloromethane|carbontetrachloride/.test(raw(c)) && !c.any(['RX', 'AcCl', 'RCOCl', 'BzCl']),
  run(m, c) { if (canonS(smiOf(m)) !== 'c1ccccc1') return null; const ex = c.eq === Infinity;
    return finish(m, ex ? 'ClC(Cl)(Cl)c1ccccc1' : 'ClC(c1ccccc1)(c1ccccc1)c1ccccc1', { title: ex ? 'Excess CCl₄: one FC alkylation → C₆H₅CCl₃' : 'FC alkylation three times → triphenylmethyl chloride', text: 'AlCl₃ makes ⁺CCl₃ from CCl₄; the ring attacks it. With excess benzene the product (ArCCl₂⁺ …) alkylates two more rings to give Ph₃C–Cl; with excess CCl₄ it stops at PhCCl₃.', tags: ['SEAr'] }); } });
front({ id: 'tetralin', name: 'Benzene + 1,4-dihalobutane → tetralin', fam: 'Aromatic', canon: 'Cl(CH₂)₄Cl, AlCl₃', db: [], multi: true,
  test: c => c.has('FeX3'),
  run(m) { const s = smiOf(m); if (!s) return null; const fr = frags(s); const ar = fr.find(x => canonS(x) === 'c1ccccc1'); const dh = fr.find(x => M(fresh(x), '[Cl,Br,I][CH2][CH2][CH2][CH2][Cl,Br,I]').length); if (!ar || !dh) return null;
    return finish(m, 'c1ccc2c(c1)CCCC2', { title: 'Two Friedel–Crafts alkylations → tetralin', text: 'The first C–Cl alkylates the ring; the second end then attacks the ortho position of the same ring (intramolecular FC), closing a six-membered ring.', tags: ['SEAr'] }); } });
front({ id: 'neopentyl_solv', name: 'Neopentyl halide solvolysis (rearranges)', fam: 'Alkyl halide', canon: 'H₂O, Δ', db: [],
  test: c => (c.has('H2O') || c.flag('aq') || c.has('AgNO3')) && !c.any(['KOH', 'aqKOH', 'alcKOH', 'NaOEt', 'NaCN', 'NaI']),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('[CH3:5][C:1]([CH3:6])([CH3:7])[CH2:2][Cl,Br,I]>>[CH3:5][C:1]([CH3:6])(O)[CH2:2][CH3:7]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Neopentyl halide: 1° cation → 1,2-CH₃ shift → 3° cation → 3° alcohol', text: 'SN2 is blocked by the t-butyl group, so on heating in water the slow SN1 path runs: the 1° carbocation rearranges at once by a methyl shift to the 3° tert-amyl cation, which water traps (some alkene too).', minor: ['CC=C(C)C'], tags: ['SN1', 'rearrangement'] }); } });
front({ id: 'sugar_hi_ring', name: 'Cyclic sugar + HI', fam: 'Biomolecule', canon: 'HI, Δ', db: [],
  test: c => c.has('HI'),
  run(m) { if (M(m, '[CX4;R]([OX2;R])[OX2H1]').length && M(m, '[CX4][OX2H1]').length >= 3 && !E.aromSet(m).size) { const n = [...m.atoms.values()].filter(a => a.el === 'C').length; if (n < 5 || n > 7) return null;
      return finish(m, 'C'.repeat(n), { title: 'Sugar + HI (prolonged heating) → n-alkane', text: 'The ring opens to the open-chain form; HI/heat replaces every O by H: all the carbons form one unbranched chain (glucose → n-hexane), showing the straight-chain skeleton.', tags: ['red'] }); } return null; } });
/* ester pyrolysis (syn elimination) needs strong heating (~500 °C), not a plain reflux */
wrap('ester_pyrolysis', function (m, c) { if ((c.T != null && c.T >= 573) || /pyroly|500|773|high\s*temp/.test(raw(c))) return undefined; return null; });
/* CCl₄ is not hydrolysed by aqueous base (no reaction in JEE terms) */
wrap('rx_base', function (m) { const s = canonS(smiOf(m)); if (s && /^(ClC\(Cl\)\(Cl\)Cl|BrC\(Br\)\(Br\)Br)$/.test(s)) return none('No reaction', 'CCl₄ has no C–H and no easy SN1/SN2 path: it is not hydrolysed by aqueous OH⁻ under these conditions.'); return undefined; });
})(typeof window !== 'undefined' ? window : globalThis);
