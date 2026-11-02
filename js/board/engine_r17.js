/* engine_r17.js: run 35. General rules from the blind sets B (SKM Blue Book JEE Advanced) and C (reagent x functional group):
   LiAlH4 on every group, tosylates as halides, oxidants block base side-reactions, Schmidt/Wolff, intramolecular FC acylation,
   naphthalene/heteroaromatic bromination, interhalogens, NCS, Ag coupling, thiols, Lucas 1° rule, HX shifts and more (qbench/run34_todo.md). */
(function (root) {
'use strict';
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
const isoOf = m => { try { return (root.Stereo && root.Stereo.iso(m)) || C.canon(m); } catch (e) { return smiOf(m); } };
const nC = s => (String(s).replace(/Cl|Ca|Cu|Cr|Cs|Cd/g, '').match(/[Cc]/g) || []).length;
const add = r => E.REC.push(r), front = r => E.REC.unshift(r);
const byC = arr => arr.slice().sort((a, b) => nC(b) - nC(a) || b.length - a.length);
function finish(m, prodSmi, o) { const R = fresh(prodSmi); if (!R) return null; return out(Object.assign({}, o, { product: R, steps: [step(m, R, o.text || o.title, o.cond)], extra: (o.extra || []).map(fresh).filter(Boolean), minor: (o.minor || []).map(fresh).filter(Boolean) })); }
const none = (title, why) => out({ title, product: null, none: why });
function applyAll(smi, T, lim) { let cur = smi; const extra = [], labels = [];
  for (const [sma, label] of T) for (let g = 0; g < (lim || 6); g++) { const p = rxn(sma, [cur]); if (!p.length) break; const frs = byC(p[0].join('.').split('.')); if (frs[0] === cur) break; cur = frs[0]; extra.push(...frs.slice(1)); if (!labels.includes(label)) labels.push(label); }
  return { cur, extra, labels }; }
const ctxWith = (c, add2, flags) => { const st = Object.assign({}, c.st, { sp: new Map(c.st.sp), flags: new Set([...c.st.flags, ...(flags || [])]) }); for (const [k, v] of Object.entries(add2 || {})) { if (v === null) st.sp.delete(k); else st.sp.set(k, Object.assign({ sp: k }, v)); } return RP().stageCtx(st); };
const recipe = id => E.REC.find(x => x.id === id);
function wrap(id, fn) { const r = recipe(id); if (!r) return; const run0 = r.run; r.run = function (m, c, all, k) { let o; try { o = fn.call(this, m, c, all, k, run0.bind(this)); } catch (e) { o = undefined; } return o === undefined ? run0.call(this, m, c, all, k) : o; }; }
function retest(id, fn) { const r = recipe(id); if (!r) return; const t0 = r.test; r.test = function (c, all, k) { return fn(c, all, k, (cc, aa, kk) => { try { return t0(cc, aa, kk); } catch (e) { return false; } }); }; }
const frags = s => String(s || '').split('.');
const mainFrag = s => byC(frags(s))[0];
const raw = c => String(c.raw || '').toLowerCase();
const isAlkoxide = c => c.any(['NaOEt', 'NaOMe', 'NaOR', 'KOtBu']);
const alkoxR = c => c.has('NaOEt') ? 'CC' : c.has('NaOMe') ? 'C' : c.has('KOtBu') ? 'C(C)(C)C' : ((c.get('NaOR') || {}).smi || 'CC');
const CARB = { HCHO: 'C=O', CH3CHO: 'CC=O', PhCHO: 'O=Cc1ccccc1', acetone: 'CC(C)=O', acrolein: 'C=CC=O' };
function partnerCarbonyl(c) { for (const [k, s] of Object.entries(CARB)) if (c.has(k)) return { k, smi: (c.get(k) || {}).smi || s }; const p = c.get('partner'); const ps = Array.isArray(p) ? p : p ? [p] : []; for (const x of ps) if (x.smi && M(fresh(x.smi), '[CX3H1,CX3H2]=O').length) return { k: 'partner', smi: x.smi }; return null; }



/* ===== r17 recipes (run 35): general rules from the blind sets B (SKM Blue Book JEE Advanced) and C (reagent x functional group) ===== */
const heat = c => !!(c.heat || c.flag('warm') || c.flag('reflux') || (c.T != null && c.T >= 323) || /\bwarm|reflux|boil|\bheat|Δ/i.test(c.raw || ''));
const OXID = ['KMnO4', 'K2Cr2O7', 'OsO4', 'O3', 'HIO4', 'mCPBA', 'Tollens', 'Fehling', 'PCC'];
const hasOx = c => c.any(OXID);
const sugarLike = m => M(m, '[CX4][OX2H1]').length >= 3 && (M(m, '[CX3]=O').length > 0 || M(m, '[CX4;R]([OX2;R])[OX2H1]').length > 0);
const nCn = (m, id) => m.nbrs(id).filter(n => el(m, n) === 'C').length;
const partners = c => [c.get('partner'), c.get('org'), c.get('RX')].flat().filter(x => x && x.smi);
const R17 = {};

/* LiAlH4 reduces every reducible group at once */
const LAH_RED = '[$([CX3](=O)[OX2][#6]),$([CX3](=O)[OX2H1]),$([CX3H1](=O)[#6]),$([CX3](=O)([#6])[#6]),$([CX3](=O)[NX3]),$([CX2]#[NX1]),$([CX3](=O)[Cl,Br])]';
wrap('hyd_red', function (m, c, all, k, run0) { let o = run0(m, c, all, k); if (!o || !o.product || o.none || !c.has('LiAlH4')) return o;
  for (let g = 0; g < 4; g++) { if (!M(o.product, LAH_RED).length) break; let o2 = null; try { o2 = run0(o.product, c, all, k); } catch (e) { o2 = null; }
    if (!o2 || !o2.product || o2.none || smiOf(o2.product) === smiOf(o.product)) break;
    o = Object.assign({}, o2, { steps: (o.steps || []).concat(o2.steps || []), extra: (o.extra || []).concat(o2.extra || []), title: 'LiAlH₄ reduces every C=O group (and C≡N) in the molecule', why: (o.why || []).concat(['LiAlH₄ is not selective: aldehyde, ketone, ester, acid and amide groups are all reduced in one go (NaBH₄ stops at aldehydes and ketones)']) }); }
  return o; });
/* nitrile + Grignard stops at the ketone, even with excess RMgX */
wrap('grignard_add', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || !/nitrile → ketone, then/.test(o.title || '') || !M(m, '[CX2]#[NX1]').length) return o;
  let o1 = null; try { o1 = run0(m, RP().stageCtx(Object.assign({}, c.st, { eq: 1 })), all, k); } catch (e) { o1 = null; }
  if (!o1 || !o1.product || /then Grignard/.test(o1.title || '')) return o;
  return Object.assign({}, o1, { title: 'Grignard + nitrile → ketone (the imine salt does not take a second RMgX, even in excess)', why: (o1.why || []).concat(['The R–C(=N⁻ MgX⁺)–R′ salt is negatively charged, so excess Grignard cannot add again: aqueous work-up gives the ketone']) }); });
/* alkyl tosylates / mesylates: the leaving group behaves like Br⁻ */
front({ id: 'sulfonate_lg', name: 'Alkyl tosylate / mesylate (leaves like Br⁻)', fam: 'Alkyl halide', canon: 'Nu⁻ or base', db: [],
  test: c => c.any(['NaI', 'NaCN', 'KOH', 'aqKOH', 'alcKOH', 'NaOEt', 'NaOMe', 'NaOR', 'KOtBu', 'NaN3', 'H2O', 'NH3', 'ROH', 'NaSH', 'AgCN', 'KNO2', 'LiAlH4']) || heat(c),
  run(m, c) { if (!M(m, '[CX4][OX2][SX4](=O)(=O)[#6]').length) return null; const s = smiOf(m); if (!s) return null;
    const p = rxn('[CX4:1][OX2][SX4](=O)(=O)[#6]>>[C:1]Br', [s]); if (!p.length) return null; const mb = fresh(mainFrag(p[0].join('.'))); if (!mb) return null;
    let o = null; try { o = root.EngineSolve.solve(mb, c.raw || ''); } catch (e) { o = null; }
    if (!o || o.status !== 'ok' || !o.product) return null; const ps = smiOf(o.product); if (!ps || /Br/.test(ps)) return null;
    return finish(m, ps, { title: 'Tosylate (–OTs) leaves like Br⁻: ' + String(o.title || '').replace(/^[^:]*:\s*/, ''), text: 'TsO⁻ is a very weak base (anion of a strong sulphonic acid), so –OTs is as good a leaving group as –Br: substitution and elimination go exactly as for the alkyl bromide.', extra: (o.extra || []).map(smiOf).filter(x => x && !/S\(=O\)/.test(x)), tags: o.tags || [] }); } });
/* aqueous OH⁻ substitutes (SN1 for 3°); elimination needs alcoholic KOH */
wrap('rx_base', function (m, c, all, k, run0) { const o = run0(m, c, all, k);
  if (!o || !o.product || !/E2/.test(o.title || '') || c.has('alcKOH') || !(c.has('aqKOH') || c.flag('aq') || c.has('H2O')) || M(m, '[CX4][Cl,Br,I]').length !== 1) return o;
  const s = smiOf(m); const p = s && rxn('[CX4:1][Cl,Br,I]>>[C:1]O', [s]); if (!p || !p.length) return o;
  return finish(m, mainFrag(p[0].join('.')), { title: 'Aqueous OH⁻ → substitution: R–X → R–OH', text: 'In water OH⁻ is heavily solvated and acts as a nucleophile: substitution (SN1 via the carbocation for a 3° halide, SN2 for 1°). Elimination needs alcoholic KOH.', tags: ['sub'] }); });
front({ id: 'acylx_oh', name: 'Acyl chloride + aqueous NaOH', fam: 'Acid derivative', canon: 'aq. NaOH', db: [],
  test: c => c.any(['KOH', 'aqKOH']) && !c.has('alcKOH') && !c.any(['NH3', 'ROH', 'org', 'partner', 'NaOCl', 'I2', 'Br2', 'Cl2', 'Br2/NaOH']),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX3](=O)[Cl,Br]').length || M(m, '[NX3;H1,H2;!$(NC=O)]').length || M(m, '[OX2H1]c').length) return null; const r = applyAll(s, [['[CX3:1](=[O:2])[Cl,Br]>>[C:1](=[O:2])[O-]', 'COCl → COO⁻']], 3); if (!r || r.cur === s) return null;
    return finish(m, r.cur, { title: 'Acyl chloride + aq. NaOH → carboxylate (fast hydrolysis)', text: 'OH⁻ adds to the C=O and Cl⁻ is expelled (addition–elimination); in base the acid stays as RCOO⁻ Na⁺ (acidify to get RCOOH).', tags: ['hyd'] }); } });
front({ id: 'epox_oh', name: 'Epoxide + OH⁻ → 1,2-diol', fam: 'Ether', canon: 'aq. NaOH', db: [],
  test: c => c.any(['KOH', 'aqKOH']) && !c.any(['RLi', 'R2CuLi', 'RMgX', 'partner', 'org', 'alcKOH', 'NaCN', 'NH3', 'NaOEt', 'NaOMe', 'RMgX', 'LiAlH4', 'org', 'partner']),
  run(m) { const s = smiOf(m); if (!s || m.frags().length > 1 || M(m, '[Li,Mg,Cu,Zn,Cd]').length || !M(m, '[CX4;r3]1[OX2;r3][CX4;r3]1').length) return null; const p = rxn('[C:1]1[O:2][C:3]1>>[C:1]([O:2])[C:3]O', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'OH⁻ opens the epoxide (SN2, anti) → 1,2-diol', text: 'The strained ring opens: OH⁻ attacks the less hindered carbon from the back, the ring O becomes –O⁻ and takes H⁺ from water.', tags: ['anti'] }); } });
front({ id: 'epox_h3o', name: 'Epoxide + H₃O⁺ → trans-1,2-diol', fam: 'Ether', canon: 'H₃O⁺', db: [],
  test: c => c.any(['H3O+', 'H2SO4']) && !c.any(['RLi', 'R2CuLi', 'RMgX', 'partner', 'org', 'ROH', 'HBr', 'HCl', 'HI', 'NaCN', 'NH3', 'KOH', 'aqKOH', 'LiAlH4', 'RMgX', 'HIO4', 'org', 'partner']),
  run(m) { const s = smiOf(m); if (!s || m.frags().length > 1 || M(m, '[Li,Mg,Cu,Zn,Cd]').length || !M(m, '[CX4;r3]1[OX2;r3][CX4;r3]1').length) return null; const p = rxn('[C:1]1[O:2][C:3]1>>[C:1]([O:2])[C:3]O', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'H₃O⁺ opens the epoxide → trans-1,2-diol (anti overall)', text: 'Protonation makes the ring O a good leaving group; water attacks the carbon from the back side, so the two OH groups end up trans.', tags: ['anti'] }); } });
front({ id: 'acid_alkoxide', name: 'Carboxylic acid + alkoxide (acid–base)', fam: 'Acid–base', canon: 'NaOEt', db: [],
  test: c => isAlkoxide(c) && !c.any(['org', 'partner', 'RX']),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX3](=O)[OX2H1]').length || M(m, '[CX3](=O)[OX2][#6]').length || M(m, '[CX4][Cl,Br,I]').length) return null; const p = rxn('[CX3:1](=[O:2])[OX2H1:3]>>[C:1](=[O:2])[O-:3]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Acid–base: RCOOH + RO⁻ → RCOO⁻ + ROH (no ester forms)', text: 'An alkoxide is a strong base: it just removes the acidic O–H proton. Esterification needs acid catalysis (Fischer) or an acyl chloride.', tags: ['acidbase'] }); } });
/* an oxidant (or borane) in the stage: no base / acid side-reactions on the substrate */
for (const id of ['aldol', 'acid_heat_dehyd', 'dehyd', 'acid_aldol', 'saponify', 'ester_sapon']) wrap(id, function (m, c) { if (hasOx(c) || c.has('BH3')) return null; return undefined; });
for (const id of ['pinacol2', 'pinacol', 'acid_heat_dehyd', 'dehyd', 'aldol']) wrap(id, function (m) { if (sugarLike(m)) return null; return undefined; });
front({ id: 'kmno4_cho', name: 'KMnO₄ oxidises aldehydes', fam: 'Carbonyl', canon: 'KMnO₄', db: [],
  test: c => c.has('KMnO4'),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX3H1](=O)[#6]').length || M(m, 'C=C').length || M(m, 'C#C').length || sugarLike(m) || M(m, 'c[CX4;!H0]').length) return null;
    const r = applyAll(s, [['[CX3H1:1](=[O:2])[#6:3]>>[C:1](=[O:2])(O)[*:3]', 'CHO → COOH']], 4); if (!r || r.cur === s) return null;
    return finish(m, r.cur, { title: 'KMnO₄: –CHO → –COOH (aldehydes reduce permanganate, hot or cold)', text: 'The aldehyde C–H is easily oxidised: even cold, dilute KMnO₄ is decolourised. Ketones resist this mild oxidation.', tags: ['ox'] }); } });
/* Fischer esterification wins over Friedel–Crafts alkylation for RCOOH + ROH, conc. H2SO4 */
wrap('fc_alk_drawn', function (m) { if (M(m, '[CX3](=O)[OX2H1]').length && M(m, '[CX4][OX2H1]').length) return null; return undefined; });
/* Lucas reagent: a 1° alcohol does not react at room temperature */
wrap('alc_hx', function (m, c) {
  if ((c.has('Lucas') || (c.has('ZnCl2') && c.has('HCl'))) && !heat(c)) { const oh = M(m, '[CX4:1][OX2H1]'); if (oh.length && oh.every(([cc]) => nCn(m, cc) <= 1) && !M(m, '[CX4]([OX2H1])[c,$(C=C)]').length)
    return none('Lucas test: no turbidity (1° alcohol)', 'At room temperature a 1° alcohol does not react with Lucas reagent (conc. HCl + anhyd. ZnCl₂): a 3° alcohol clouds at once, a 2° one in about 5 minutes, a 1° one only on heating.'); }
  return undefined; });
front({ id: 'alc_hx_shift', name: 'Alcohol + HX with a 1,2-shift', fam: 'Alcohol', canon: 'HX', db: [],
  test: c => c.any(['HCl', 'HBr', 'HI', 'Lucas']) && !c.has('ROOR'),
  run(m, c) { const s = smiOf(m); if (!s) return null; const X = c.has('HBr') ? 'Br' : c.has('HI') ? 'I' : 'Cl';
    const T = [['[CH3:5][CX4H0:1]([CH3:6])([CH3:7])[CH2:2][OX2H1]>>[CH3:6][C:1](' + X + ')([CH3:7])[CH2:2][CH3:5]', 'a methyl shift (neopentyl 1° → 3° cation)'],
      ['[CH3:5][CX4H0:1]([#6:6])([#6:7])[CX4H1:2]([#6:8])[OX2H1]>>[*:6][C:1](' + X + ')([*:7])[CH1:2]([*:8])[CH3:5]', 'a methyl shift (2° → 3° cation)'],
      ['[CX4H1:1]([CX4:6])([CX4:7])[CX4H1:2]([#6:8])[OX2H1]>>[*:6][C:1](' + X + ')([*:7])[CH2:2][*:8]', 'a hydride shift (2° → 3° cation)']];
    for (const [sma, lab] of T) { const p = rxn(sma, [s]); if (p.length) return finish(m, mainFrag(p[0].join('.')), { title: 'HX with ' + lab + ': the halide ends up on the 3° carbon', text: 'The protonated –OH leaves as water; the cation rearranges to the more stable 3° cation, which ' + X + '⁻ captures (SN1). The rearranged halide is the main product.', tags: ['rearr'] }); }
    return null; } });
/* alkynes */
front({ id: 'alkyne_hg_h2o', name: 'Alkyne + Hg²⁺, H₂O → ketone', fam: 'Alkyne', canon: 'Hg(OAc)₂, H₂O', db: ['Hydration of alkynes'],
  test: c => /hg\s*\(\s*oac\s*\)\s*2|hg\s*\(\s*ch3coo\s*\)\s*2|mercuric\s+acetate/i.test(c.raw || ''),
  run(m) { const s = smiOf(m); if (!s || !M(m, 'C#C').length || M(m, 'C=C').length) return null; let p = rxn('[CH1:1]#[C:2][#6:3]>>[CH3:1][C:2](=O)[*:3]', [s]); if (!p.length) p = rxn('[CH1:1]#[CH1:2]>>[CH3:1][CH1:2]=O', [s]); if (!p.length) p = rxn('[#6:3][C:1]#[C:2][#6:4]>>[*:3][CH2:1][C:2](=O)[*:4]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Hg²⁺-catalysed hydration of the alkyne → ketone (Markovnikov; the enol tautomerises)', text: 'Hg²⁺ activates the C≡C; water adds Markovnikov to give an enol, which tautomerises to the ketone (ethyne gives ethanal).', tags: ['add'] }); } });
front({ id: 'alkyne_hbr_rad', name: 'HBr / peroxide on a terminal alkyne', fam: 'Alkyne', canon: 'HBr, ROOR', db: [],
  test: c => c.has('HBr') && c.has('ROOR'),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CH1]#C').length || M(m, 'C=C').length) return null; const p = rxn('[CH1:1]#[C:2]>>[CH1:1](Br)=[CH1:2]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'HBr / peroxide: Br• adds to the terminal carbon (anti-Markovnikov) → 1-bromoalkene', text: 'With peroxide, Br• adds first, to the end carbon, giving the more stable vinyl radical; H–Br then gives H•. One equivalent gives the 1-bromoalkene (E/Z mixture).', tags: ['radical'] }); } });
front({ id: 'alkyne_brw', name: 'Bromine water + alkyne', fam: 'Alkyne', canon: 'Br₂ water', db: [],
  test: c => c.has('HOBr'),
  run(m) { const s = smiOf(m); if (!s || !M(m, 'C#C').length || M(m, 'C=C').length || M(m, 'c').length) return null; const p = rxn('[C:1]#[C:2]>>[C:1](Br)=[C:2]Br', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Bromine water is decolourised: Br₂ adds across C≡C → 1,2-dibromoalkene', text: 'Alkynes, like alkenes, add Br₂ (the unsaturation test): one Br₂ gives the (E)-1,2-dibromoalkene, more gives the 1,1,2,2-tetrabromide.', tags: ['add', 'test'] }); } });
/* catalytic H2: C=O, NO2, C=C, C#N go before an arene ring; Clemmensen leaves C=C alone */
wrap('arene_h2', function (m) { if (M(m, '[CX3]=O').length || M(m, '[N+](=O)[O-]').length || M(m, 'C=C').length || M(m, 'C#C').length || M(m, 'C#N').length) return null; return undefined; });
wrap('hx', function (m, c) { if (/zn\s*[-(]?\s*hg|amalgam|clemmensen/i.test(c.raw || '')) return null; return undefined; });
/* chloral hydrate; no SN2 hydrolysis of a CCl3 next to C=O */
front({ id: 'chloral_h2o', name: 'Chloral + water → chloral hydrate', fam: 'Carbonyl', canon: 'H₂O', db: [],
  test: c => c.has('H2O') && !c.any(['KOH', 'aqKOH', 'NaOCl', 'H2SO4', 'H3O+', 'HCl', 'Br2']),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX3H1](=O)C(Cl)(Cl)Cl').length) return null; const p = rxn('[CX3H1:1](=[O:2])[C:3]>>[C:1]([O:2])(O)[*:3]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Chloral + H₂O → chloral hydrate (a stable gem-diol)', text: 'Three Cl atoms pull electron density away from the C=O carbon, so water adds and the hydrate is stable (isolable crystals). The C–Cl bonds are not hydrolysed by water.', tags: ['add'] }); } });
wrap('rx_oh', function (m) { if (M(m, '[CX4](Cl)(Cl)Cl').length && M(m, '[CX3]=O').length) return null; return undefined; });
/* benzilic acid rearrangement needs a 1,2-diketone, not oxalic acid */
wrap('benzilic', function (m) { if (!M(m, '[CX3;!$(C[OX2H1]);!$(C[OX1-])](=O)[CX3;!$(C[OX2H1]);!$(C[OX1-])]=O').length) return null; return undefined; });
/* (run 36) excess HCHO + ketone: the r16 recipe stays (JEE official: the hydroxymethylated ketone, now with no α-H, is reduced by crossed Cannizzaro) */
/* nitration: at most three NO2, o/p to the activating group; benzene stops at m-dinitro */
wrap('eas', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || !/Nitration/.test(o.title || '')) return o; const nn = M(o.product, '[N+](=O)[O-]').length - M(m, '[N+](=O)[O-]').length; if (nn <= 3) return o;
  const s = smiOf(m); let ps = null;
  if (s === 'c1ccccc1') ps = 'O=[N+]([O-])c1cccc([N+](=O)[O-])c1';
  else { const p = rxn('[cH1:2]1[cH1:3][cH1:4][cH1:5][cH1:6][c:1]1[CX4,OX2,NX3:9]>>[c:2]1([N+](=O)[O-])[cH1:3][c:4]([N+](=O)[O-])[cH1:5][c:6]([N+](=O)[O-])[c:1]1[*:9]', [s]); if (p.length) ps = mainFrag(p[0].join('.')); }
  if (!ps) return o;
  return finish(m, ps, { title: s === 'c1ccccc1' ? 'Fuming HNO₃ / H₂SO₄, heat: benzene → m-dinitrobenzene' : 'Exhaustive nitration → 2,4,6-trinitro product (e.g. TNT)', text: 'Each NO₂ deactivates the ring, so substitution stops once the positions ortho and para to the activating group are filled: three NO₂ at most, never two side by side.', tags: ['eas'] }); });
/* aniline in the nitrating mixture: much meta (anilinium) as well as para */
wrap('eas', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || !/Nitration/.test(o.title || '') || !M(m, 'c[NX3H2]').length) return o; const s = smiOf(m);
  const p = rxn('[c:1]([NX3H2:9])[cH1:2][cH1:3]>>[c:1]([N:9])[c:2][c:3][N+](=O)[O-]', [s]); if (!p.length) return o; const R = fresh(mainFrag(p[0].join('.'))); if (!R) return o;
  o.minor = (o.minor || []).concat([R]); o.why = (o.why || []).concat(['About 47 % m-nitroaniline forms too: in the nitrating mixture much of the aniline is protonated to anilinium (–NH₃⁺ is meta-directing). p-Nitroaniline (51 %) is still the largest share; protect –NH₂ as acetanilide for clean para nitration.']); return o; });
/* radical chlorination does not go on from BHC; diaryl ethers resist HI */
wrap('radx', function (m) { if (M(m, 'ClC1C(Cl)C(Cl)C(Cl)C(Cl)C1Cl').length) return null; return undefined; });
wrap('ether_hi', function (m) { if (M(m, '[OX2](c)c').length && !M(m, '[OX2][CX4]').length) return none('No reaction: a diaryl ether is not cleaved by HI', 'Both C–O bonds are to sp² ring carbons (partial double-bond character): neither SN1 nor SN2 can happen at an aryl carbon, so Ph–O–Ph survives HI.'); return undefined; });
/* X2 + conjugated diene: 1,4 (thermodynamic) unless the temperature is low, the same rule as HBr */
wrap('diene_add', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || !/1,2-addition \(kinetic\)/.test(o.title || '') || !/Br₂|Cl₂/.test(o.title || '')) return o;
  if ((c.T != null && c.T < 273) || /low\s*temp|kinetic|cold|−\s*\d|-\s*\d+\s*°/i.test(c.raw || '')) return o;
  const alt = [...(o.extra || []), ...(o.minor || []).map(x => x && x.m || x)].find(x => x && x.atoms && M(x, '[Cl,Br][CX4][CX3]=[CX3][CX4][Cl,Br]').length); if (!alt) return o;
  return Object.assign({}, o, { product: alt, extra: (o.extra || []).filter(x => x !== alt), minor: [o.product], title: o.title.replace('1,2-addition (kinetic)', '1,4-addition (thermodynamic; at low temperature the 1,2-adduct is main)'), why: (o.why || []).concat(['No temperature given: the more stable internal alkene (the 1,4-adduct) is the main product, as with HBr']) }); });
/* intramolecular Friedel–Crafts acylation: indanone, tetralone, anthrone, anthraquinone */
front({ id: 'intra_fc_acyl', name: 'Intramolecular Friedel–Crafts acylation (ring closure)', fam: 'Aromatic', canon: 'H₃PO₄ / PPA, heat; or AlCl₃', db: [],
  test: c => (c.any(['H3PO4', 'H2SO4', 'HF']) && heat(c)) || c.has('FeX3') || /\bppa\b|polyphosphoric/i.test(c.raw || ''),
  run(m, c) { const s = smiOf(m); if (!s) return null; const acid = M(m, '[CX3](=O)[OX2H1]').length > 0;
    if (acid && !(c.any(['H3PO4', 'H2SO4', 'HF']) || /\bppa\b|polyphosphoric/i.test(c.raw || ''))) return null; if (!acid && !M(m, '[CX3](=O)[Cl,Br]').length) return null;
    const T = [['[cH1:1]:[c:2][CX4:3][CX4:4][CX3:5](=[O:6])[OX2H1,Cl,Br]>>[c:1]1:[c:2][C:3][C:4][C:5]1=[O:6]', 'a five-membered ring (1-indanone type)'],
      ['[cH1:1]:[c:2][CX4:3][CX4:4][CX4:7][CX3:5](=[O:6])[OX2H1,Cl,Br]>>[c:1]1:[c:2][C:3][C:4][C:7][C:5]1=[O:6]', 'a six-membered ring (1-tetralone type)'],
      ['[cH1:1]:[c:2][CX4:3][c:4]:[c:7][CX3:5](=[O:6])[OX2H1,Cl,Br]>>[c:1]1:[c:2][C:3][c:4]:[c:7][C:5]1=[O:6]', 'the central ring of anthrone'],
      ['[cH1:1]:[c:2][CX3:3](=[O:8])[c:4]:[c:7][CX3:5](=[O:6])[OX2H1,Cl,Br]>>[c:1]1:[c:2][C:3](=[O:8])[c:4]:[c:7][C:5]1=[O:6]', 'the central ring of anthraquinone']];
    for (const [sma, lab] of T) { const p = rxn(sma, [s]); if (p.length) return finish(m, mainFrag(p[0].join('.')), { title: 'Intramolecular Friedel–Crafts acylation closes ' + lab, text: 'The acylium ion (from –COOH + strong acid, or –COCl + AlCl₃) is attacked by the ring it is tethered to; five- and six-membered rings close easily.', tags: ['eas'] }); }
    return null; } });
/* strong KMnO4: an acyl side chain is cut back to –COOH */
wrap('kmno4hot', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (o && o.product && !o.none) return o; if (!c.has('KMnO4') || c.flag('mild')) return o; const s = smiOf(m); if (!s || !M(m, 'c[CX3](=O)[#6]').length) return o;
  const r = applyAll(s, [['[c:1][CX3](=O)[#6]>>[c:1]C(=O)O', 'ArCOR → ArCOOH']], 3); if (!r || r.cur === s) return o;
  return finish(m, r.cur, { title: 'Hot KMnO₄: the acyl side chain is cut back to –COOH (ArCOR → ArCOOH)', text: 'Strong oxidation attacks the benzylic position; whatever is attached to the ring through carbon ends up as –COOH.', tags: ['ox'] }); });
/* chromic acid / KMnO4 oxidise every 2° alcohol (and 1° → acid) */
wrap('kmno4hot', function (m, c, all, k, run0) { let o = run0(m, c, all, k); if (!o || !o.product || o.none || !c.any(['K2Cr2O7', 'KMnO4'])) return o;
  for (let g = 0; g < 3; g++) { if (!M(o.product, '[CX4;H1,H2]([#6])[OX2H1]').length) break; let o2 = null; try { o2 = run0(o.product, c, all, k); } catch (e) { o2 = null; } if (!o2 || !o2.product || o2.none || smiOf(o2.product) === smiOf(o.product)) break;
    o = Object.assign({}, o2, { steps: (o.steps || []).concat(o2.steps || []), title: 'Every 2° –OH → C=O (1° → –COOH): ' + String(o2.title || '') }); }
  return o; });
front({ id: 'ald_hypohalite', name: 'Aldehyde + NaOX (Br₂ / NaOH)', fam: 'Carbonyl', canon: 'Br₂, NaOH', db: [],
  test: c => c.any(['NaOCl', 'Br2/NaOH']) && !((c.get('NaOCl') || {}).X === 'I') && !c.has('I2') && !/\bi2\b|iod/i.test(c.raw || ''),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX3H1](=O)[#6]').length || M(m, '[CH3][CX3]=O').length || M(m, '[CH3][CX4H1][OX2H1]').length || M(m, '[NX3H2]C=O').length || M(m, '[CH2][CX3H1]=O').length) return null;
    const p = rxn('[CX3H1:1](=[O:2])[#6:3]>>[C:1](=[O:2])([O-])[*:3]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Hypohalite (Br₂ / NaOH) oxidises –CHO → –COO⁻', text: 'NaOBr formed in situ is an oxidant: an aldehyde with no CH₃–C=O unit is simply oxidised to the carboxylate (acidify to get RCOOH).', tags: ['ox'] }); } });
/* Schmidt and Wolff */
front({ id: 'schmidt', name: 'Schmidt reaction (HN₃, H₂SO₄)', fam: 'Amine', canon: 'HN₃, conc. H₂SO₄', db: [],
  test: c => c.has('HN3') || /\bhn3\b|hydrazoic/i.test(c.raw || ''),
  run(m) { const s = smiOf(m); if (!s) return null;
    if (M(m, '[CX3](=O)[OX2H1]').length) { const p = rxn('[#6:1][CX3](=O)[OX2H1]>>[*:1]N', [s]); if (p.length) return finish(m, mainFrag(p[0].join('.')), { title: 'Schmidt: R–COOH + HN₃ → R–NH₂ (+ CO₂ + N₂)', text: 'As in Hofmann and Curtius: the acyl azide loses N₂ while R migrates from C to N (isocyanate), which hydrolyses to the amine and CO₂. One carbon fewer.', extra: ['O=C=O'], tags: ['rearr'] }); }
    if (M(m, '[#6][CX3](=O)[#6]').length) { const cand = rxn('[#6:1][CX3:2](=[O:3])[#6:4]>>[*:1][NH1][C:2](=[O:3])[*:4]', [s]).map(x => mainFrag(x.join('.')));
      const sc = smi => { const P = fresh(smi); if (!P) return -9; const q = M(P, '[NX3H1]([#6:1])[CX3]=O'); if (!q.length) return -9; const cn = q[0][1]; return (E.aromSet(P).has(cn) ? 10 : 0) + P.nbrs(cn).filter(n => el(P, n) === 'C').length; };
      const best = cand.filter(Boolean).sort((a, b) => sc(b) - sc(a))[0]; if (best) return finish(m, best, { title: 'Schmidt: ketone + HN₃ → amide (the bulkier group moves to N; a cyclic ketone gives a lactam)', text: 'HN₃ adds to the protonated C=O; loss of N₂ goes with migration of the group anti to it (usually the larger one) from C to N, as in the Beckmann rearrangement.', tags: ['rearr'] }); }
    return null; } });
front({ id: 'wolff_rearr', name: 'Wolff rearrangement (α-diazoketone)', fam: 'Carbonyl', canon: 'Ag₂O, H₂O (or hν, heat)', db: ['Arndt–Eistert'],
  test: c => c.any(['Ag2O']) || c.light || heat(c),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[CX3](=O)[CX3]=[N+]=[N-]').length) return null;
    const nu = c.has('NH3') ? 'N' : c.has('ROH') && !c.has('H2O') ? 'OC' : (c.any(['H2O', 'H3O+']) || /h2o|water|moist/i.test(c.raw || '') ? 'O' : null);
    const p = rxn(nu ? '[#6:1][CX3:2](=[O:3])[CX3:4]=[N+]=[N-]>>[*:1][C:4][C:2](=[O:3])' + nu : '[#6:1][CX3:2](=[O:3])[CX3:4]=[N+]=[N-]>>[*:1][C:4]=[C:2]=[O:3]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Wolff rearrangement: R–CO–CHN₂ → ketene' + (nu === 'N' ? ' → amide' : nu === 'O' ? ' → acid (one CH₂ longer; a cyclic diazoketone contracts the ring)' : nu ? ' → ester' : ''), text: 'Ag₂O (or light / heat) removes N₂; the group on the C=O migrates to the carbene carbon, giving a ketene, which adds water (or ROH / NH₃).', extra: ['N#N'], tags: ['rearr'] }); } });
wrap('ag_salt', function (m, c) { if (c.any(['SOCl2', 'CH2N2']) || /ch2n2|diazomethane/i.test(c.raw || '')) return null; return undefined; });
/* calcium salt dry distillation (also crossed with calcium acetate / formate) */
front({ id: 'ca_dist', name: 'Dry distillation of a calcium carboxylate → ketone', fam: 'Acid derivative', canon: 'Ca(OH)₂ then dry distillation', db: [],
  test: c => /ca\s*\(\s*oh\s*\)\s*2|calcium|caco3/i.test(c.raw || '') && !/naoh|koh|soda\s*-?\s*lime/i.test(c.raw || '') && !c.has('sodalime') && (c.flag('dry') || heat(c)),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[CX3](=O)[OX2H1,OX1-]').length) return null; const r0 = raw(c);
    const R = /formate/.test(r0) ? 'H' : /acetate/.test(r0) && s !== 'CC(=O)O' ? 'C' : null; let p;
    if (R === 'H') p = rxn('[#6:1][CX3](=O)[OX2H1,OX1-]>>[*:1][CH1]=O', [s]); else if (R === 'C') p = rxn('[#6:1][CX3](=O)[OX2H1,OX1-]>>[*:1]C(=O)C', [s]);
    else p = rxn('[#6:1][CX3](=O)[OX2H1,OX1-].[#6:2][CX3](=O)[OX2H1,OX1-]>>[*:1]C(=O)[*:2]', [s, s]); if (!p || !p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: R === 'H' ? 'Dry distillation with calcium formate → aldehyde' : R === 'C' ? 'Dry distillation with calcium acetate → methyl ketone (crossed)' : 'Calcium salt, dry distillation → symmetrical ketone (+ CaCO₃)', text: '(RCOO)₂Ca on strong heating loses CaCO₃ and gives R–CO–R. With calcium formate one R is H (an aldehyde); with calcium acetate a methyl ketone forms (plus some acetone).', tags: ['decarb'] }); } });
/* amide hydrolysis needs heat (a cold base test leaves an amide alone); acetanilide + H3O+ at the end of a sequence comes off */
wrap('amide_rx', function (m, c) { if (!heat(c) && c.any(['KOH', 'aqKOH']) && !c.any(['H3O+', 'H2SO4', 'HCl', 'Br2/NaOH', 'NaOCl', 'Br2']) && !/hydroly|reflux|boil/i.test(c.raw || '')) return null; return undefined; });
front({ id: 'anilide_h3o', name: 'Acetanilide + H₃O⁺ → amine (protecting group off)', fam: 'Amine', canon: 'H₃O⁺, heat', db: [],
  test: c => c.any(['H3O+', 'H2SO4', 'HCl', 'KOH', 'aqKOH']) && (c.any(['H3O+']) || heat(c)) && !c.any(['HNO2', 'NO2+', 'HNO3', 'Br2', 'Cl2', 'FeX3', 'LiAlH4', 'oleum']),
  run(m) { const s = smiOf(m); if (!s || !M(m, 'c[NX3H1][CX3](=O)[CH3]').length) return null; const p = rxn('[c:1][NX3H1:2][CX3](=O)[CH3]>>[c:1][N:2]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Hydrolysis of the acetanilide: the –NHCOCH₃ protecting group comes off → –NH₂', text: 'Amide hydrolysis (acid or base, heat) gives back the free amine and acetic acid; this is the last step whenever –NH₂ was protected for a clean substitution.', extra: ['CC(=O)O'], tags: ['hyd'] }); } });
/* borane reduces aldehydes and ketones (no C=C present) */
front({ id: 'bh3_co', name: 'B₂H₆ reduces aldehydes and ketones', fam: 'Carbonyl', canon: 'B₂H₆ then H₂O', db: [],
  test: c => c.has('BH3'),
  run(m) { const s = smiOf(m); if (!s || M(m, 'C=C').length || M(m, 'C#C').length || M(m, '[CX3](=O)[OX2H1]').length || !M(m, '[CX3;$([CX3H1](=O)[#6]),$([CX3](=O)([#6])[#6])]=O').length) return null;
    const r = applyAll(s, [['[CX3;$([CX3H1](=O)[#6]),$([CX3](=O)([#6])[#6]):1]=[O:2]>>[C:1][O:2]', 'C=O → CH–OH']], 4); if (!r || r.cur === s) return null;
    return finish(m, r.cur, { title: 'B₂H₆ (BH₃) reduces the C=O → alcohol', text: 'Borane adds B–H across C=O (B to O, H to C); hydrolysis gives the alcohol. BH₃ reduces –COOH, aldehydes and ketones, but not esters or nitro groups.', tags: ['red'] }); } });
/* naphthalene and five-membered heteroaromatics: Br2 alone substitutes */
function naphScore(P, G) { const AS = new Set(M(P, 'c').map(x => x[0])); const xs = M(P, G || 'c[Br,Cl]'); if (!xs.length) return -99; const b = xs[0][0]; const fus = M(P, '[c;$(c(:c)(:c):c)]').map(x => x[0]);
  let sc = P.nbrs(b).some(n => fus.includes(n)) ? 3 : 0; const rings = C.rings(P, 6).filter(r => r.length === 6 && r.includes(b));
  for (const r of rings) { const i = r.indexOf(b); for (let j = 0; j < 6; j++) { const q = r[j]; if (q === b) continue; const d = Math.min(Math.abs(i - j), 6 - Math.abs(i - j)); const subs = P.nbrs(q).filter(n => !AS.has(n) && el(P, n) !== 'H');
      for (const n of subs) { const deact = M(P, '[c][$([CX3]=O),$([N+](=O)[O-]),$(S(=O)=O),$(C#N)]').some(x => x[0] === q); if (deact) sc -= 4; else if (['C', 'O', 'N'].includes(el(P, n))) sc += d === 3 ? 5 : d === 1 ? 2 : 0; } } }
  return sc; }
front({ id: 'polycyc_x2', name: 'Br₂ on naphthalene or thiophene / furan / pyrrole (no catalyst)', fam: 'Aromatic', canon: 'Br₂, CS₂ or CCl₄', db: [],
  test: c => (c.has('Br2') || c.has('Cl2')) && !c.has('FeX3') && !c.light && !c.has('HOBr') && !c.any(['KOH', 'aqKOH', 'redP']),
  run(m, c) { const s = smiOf(m); if (!s) return null; const X = c.has('Br2') ? 'Br' : 'Cl';
    if (M(m, '[cH1]1[s,o,nH1]ccc1').length && !M(m, 'C=C').length) { const p = rxn('[cH1:1]1:[s,o,nH1:2]:[c:3]:[c:4]:[c:5]1>>[c:1]1(' + X + '):[*:2]:[c:3]:[c:4]:[c:5]1', [s]); if (p.length) return finish(m, mainFrag(p[0].join('.')), { title: X + '₂ on a five-membered heteroaromatic ring → 2-substitution (next to the heteroatom)', text: 'Thiophene, furan and pyrrole are far more reactive than benzene (the heteroatom lone pair is part of the π system): the halogen enters next to the heteroatom, where the σ-complex has the most resonance forms. No Lewis acid is needed.', tags: ['eas'] }); }
    if (!M(m, 'c1ccc2ccccc2c1').length || M(m, 'C=C').length) return null;
    const cand = rxn('[cH1:1]>>[c:1]' + X, [s]).map(x => mainFrag(x.join('.'))).filter(Boolean); if (!cand.length) return null;
    const best = cand.map(x => [x, naphScore(fresh(x))]).sort((a, b) => b[1] - a[1])[0][0];
    return finish(m, best, { title: X + '₂ substitutes naphthalene without a catalyst, at an α-carbon (C-4 when C-1 carries an activating group)', text: 'Naphthalene is more reactive than benzene. Attack at C-1 (α) gives a σ-complex in which the other ring stays fully aromatic; an activating group at C-1 sends the halogen to C-4.', tags: ['eas'] }); } });
/* parser species added in rx_parse (run 35): NCS, BrCl/ICl/IBr, Ag powder */
front({ id: 'interhalogen', name: 'Interhalogen (BrCl, ICl, IBr) + alkene', fam: 'Alkene', canon: 'BrCl / ICl', db: [],
  test: c => c.any(['BrCl', 'ICl', 'IBr']),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, 'C=C').length) return null; const [X1, Y1] = c.has('BrCl') ? ['Br', 'Cl'] : c.has('ICl') ? ['I', 'Cl'] : ['I', 'Br'];
    const ps = rxn('[C:1]=[C:2]>>[C:1](' + X1 + ')[C:2]' + Y1, [s]).map(x => mainFrag(x.join('.'))).filter(Boolean); if (!ps.length) return null;
    const sc = smi => { const P = fresh(smi); if (!P) return -9; const q = M(P, '[CX4][' + Y1 + ']'); return q.length ? nCn(P, q[0][0]) : -9; };
    const best = ps.slice().sort((a, b) => sc(b) - sc(a))[0];
    return finish(m, best, { title: X1 + '–' + Y1 + ' adds Markovnikov-wise: ' + X1 + ' (the δ+ end) to the less substituted C, ' + Y1 + ' to the more substituted C', text: 'The less electronegative halogen is the electrophile (a halonium ion forms); ' + Y1 + '⁻ opens it at the carbon that carries more positive charge.', tags: ['add'] }); } });
front({ id: 'ncs', name: 'NCS: allylic / benzylic chlorination', fam: 'Radical', canon: 'NCS, hν', db: [],
  test: c => c.has('NCS'),
  run(m, c, all, k) { const r = E.REC.find(x => x.id !== 'nbs_alc' && x.id !== 'ncs' && /^nbs/i.test(x.id)); if (!r) return null; const c2 = ctxWith(c, { NBS: {}, NCS: null }); let o = null; try { o = r.run(m, c2, all, k); } catch (e) { o = null; } if (!o || !o.product) return null; const ps = smiOf(o.product); if (!ps) return null;
    return finish(m, ps.replace(/Br/g, 'Cl'), { title: 'NCS: allylic / benzylic chlorination (as NBS does with Br)', text: 'N-Chlorosuccinimide keeps the Cl₂ level tiny, so Cl• abstracts the weak allylic / benzylic H instead of adding to the C=C.', tags: ['radical'] }); } });
front({ id: 'ag_dehal', name: 'Ag powder, heat: CX₃ + CX₃ → C≡C', fam: 'Alkyl halide', canon: 'Ag, heat', db: [],
  test: c => c.has('Agm'),
  run(m) { const s = smiOf(m); if (!s) return null; const one = mainFrag(s); const mo = fresh(one); if (!mo || !M(mo, '[CX4]([Cl,Br,I])([Cl,Br,I])[Cl,Br,I]').length) return null;
    const p = rxn('[CX4:1]([Cl,Br,I])([Cl,Br,I])[Cl,Br,I].[CX4:2]([Cl,Br,I])([Cl,Br,I])[Cl,Br,I]>>[C:1]#[C:2]', [one, one]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Ag powder, heat: two CX₃ carbons lose all their halogen and join → C≡C', text: '6 Ag take the six halogens as AgX: 2 CHCl₃ (or 2 CHI₃) → HC≡CH; 2 CH₃CCl₃ → CH₃C≡CCH₃.', tags: ['coupling'] }); } });
/* thiols, cyclopentadiene: base gives the anion; thiolate + R–X → thioether */
front({ id: 'acidic_ch_sh', name: 'Base removes an acidic S–H / cyclopentadiene C–H', fam: 'Acid–base', canon: 'NaH / NaNH₂', db: [],
  test: c => c.any(['NaH', 'NaNH2', 'KOH', 'aqKOH', 'NaOEt', 'NaOMe']) && !partners(c).length,
  run(m) { const s = smiOf(m); if (!s) return null; let p = rxn('[SX2H1:1][#6:2]>>[S-:1][*:2]', [s]); let t = 'Thiol + base → thiolate RS⁻ (S–H is much more acidic than O–H)';
    if (!p.length && M(m, '[CH2]1C=CC=C1').length && M(m, 'C=C').length === 2) { p = rxn('[CH2:1]1[C:2]=[C:3][C:4]=[C:5]1>>[cH1-:1]1[cH1:2][cH1:3][cH1:4][cH1:5]1', [s]); t = 'Cyclopentadiene + strong base → cyclopentadienyl anion (aromatic, 6π electrons)'; }
    if (!p.length) return null; return finish(m, mainFrag(p[0].join('.')), { title: t, text: 'Acid–base: the conjugate base is stabilised (large, polarisable S; or aromaticity for C₅H₅⁻), so even moderate bases remove this proton.', tags: ['acidbase'] }); } });
front({ id: 'thiol_alk', name: 'Thiol / thiolate + alkyl halide → thioether', fam: 'Thiol', canon: 'NaH then R–X', db: [],
  test: c => partners(c).length > 0 || c.any(['org']),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[SX2H1,SX1-][#6]').length) return null; const fr = String(s).split('.'); const sh = fr.find(x => /S/.test(x)); let rx = fr.find(x => x !== sh && /Cl|Br|I/.test(x));
    if (!rx) { const q = partners(c).find(x => /Cl|Br|I/.test(x.smi)); rx = q && q.smi; } if (!sh || !rx) return null;
    const p = rxn('[SX2H1,SX1-:1].[CX4:2][Cl,Br,I]>>[S+0:1][C:2]', [sh, rx]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Thiolate + R–X (SN2) → thioether (sulphide)', text: 'RS⁻ is an excellent nucleophile and a weak base: clean SN2 on the alkyl halide (the sulphur analogue of the Williamson synthesis).', tags: ['sub'] }); } });
/* metals / small reagents */
front({ id: 'zn_nh4cl', name: 'Zn / NH₄Cl: nitroarene → N-arylhydroxylamine', fam: 'Nitro', canon: 'Zn, NH₄Cl', db: [],
  test: c => (c.has('Zn') && c.has('NH4Cl')) || c.has('Zn/NH4Cl'),
  run(m) { const s = smiOf(m); if (!s || !M(m, 'c[N+](=O)[O-]').length) return null; const p = rxn('[c:1][N+:2](=O)[O-]>>[c:1][N+0:2]O', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Zn / NH₄Cl (neutral): Ar–NO₂ → Ar–NHOH (N-phenylhydroxylamine)', text: 'In neutral solution the reduction stops at the hydroxylamine stage (Ar–NO₂ → Ar–N=O → Ar–NHOH); metal + acid (Sn/HCl) goes on to Ar–NH₂.', tags: ['red'] }); } });
front({ id: 'p_i2', name: 'Red P + I₂ (PI₃): R–OH → R–I', fam: 'Alcohol', canon: 'red P, I₂', db: [],
  test: c => c.has('I2') && c.has('redP') && !c.any(['HI', 'KOH', 'aqKOH', 'NaOCl']),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX4][OX2H1]').length || M(m, '[CX3](=O)[OX2H1]').length) return null; const r = applyAll(s, [['[CX4:1][OX2H1]>>[C:1]I', 'R–OH → R–I']], 4); if (!r || r.cur === s) return null;
    return finish(m, r.cur, { title: 'Red P + I₂ → PI₃ in situ: R–OH → R–I', text: '3 R–OH + PI₃ → 3 R–I + H₃PO₃ (iodides are made this way because PI₃ is not stored).', tags: ['sub'] }); } });
front({ id: 'hvz_cl', name: 'HVZ with Cl₂ / P', fam: 'Acid derivative', canon: 'Cl₂, red P', db: ['HVZ reaction'],
  test: c => c.has('Cl2') && c.has('redP'),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX4;!H0][CX3](=O)[OX2H1]').length) return null; const p = rxn('[CX4;!H0:1][CX3:2](=[O:3])[OX2H1:4]>>[C:1](Cl)[C:2](=[O:3])[O:4]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'HVZ (Cl₂, red P): α-chlorination of the acid', text: 'P turns some RCOOH into the acyl chloride, which enolises; the enol takes Cl at the α-carbon, and exchange gives back the α-chloro acid. One α-H per Cl₂.', tags: ['sub'] }); } });
front({ id: 'snar_h2o', name: 'Picryl chloride + warm water → picric acid', fam: 'Aromatic', canon: 'H₂O, warm', db: [],
  test: c => c.has('H2O') && !c.any(['KOH', 'aqKOH', 'NH3', 'NaOMe', 'NaOEt']) && heat(c),
  run(m) { const s = smiOf(m); if (!s || M(m, '[N+](=O)[O-]').length < 3) return null; const p = rxn('[c:1](Cl)(:[c:2][N+](=O)[O-]):[c:3][N+](=O)[O-]>>[c:1](O)(:[c:2][N+](=O)[O-]):[c:3][N+](=O)[O-]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Three o/p-NO₂ groups: warm water substitutes the Ar–Cl (SNAr) → picric acid', text: 'NCERT: 2,4,6-trinitrochlorobenzene + warm H₂O gives picric acid; 2,4-dinitrochlorobenzene needs aq. NaOH at 368 K, 4-nitrochlorobenzene 443 K.', tags: ['snar'] }); } });
front({ id: 'ether_pcl5', name: 'Ether + PCl₅ → alkyl chlorides', fam: 'Ether', canon: 'PCl₅, heat', db: [],
  test: c => c.has('PCl5'),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX4][OX2][CX4]').length || M(m, '[OX2;r3,r4,r5,r6]').length || M(m, '[OX2H1]').length || M(m, '[CX3]=O').length) return null; const p = rxn('[CX4:1][OX2][CX4:2]>>[C:1]Cl.[C:2]Cl', [s]); if (!p.length) return null; const fr = byC(p[0].join('.').split('.'));
    return finish(m, fr[0], { title: 'Ether + PCl₅ → two alkyl chlorides (+ POCl₃)', text: 'R–O–R′ + PCl₅ → R–Cl + R′–Cl + POCl₃ (no –OH is needed; the ether O attacks P first).', extra: fr.slice(1), tags: ['sub'] }); } });
front({ id: 'halo_acid_lactone', name: 'ω-Halo acid + base → lactone', fam: 'Acid derivative', canon: 'NaHCO₃ / base', db: [],
  test: c => c.any(['NaHCO3', 'KOH', 'aqKOH', 'NaOEt']),
  run(m) { const s = smiOf(m); if (!s) return null; for (const sma of ['[OX2H1,OX1-:1][CX3:2](=[O:3])[CX4:4][CX4:5][CX4:6][Cl,Br,I]>>[O+0:1]1[C:2](=[O:3])[C:4][C:5][C:6]1', '[OX2H1,OX1-:1][CX3:2](=[O:3])[CX4:4][CX4:5][CX4:6][CX4:7][Cl,Br,I]>>[O+0:1]1[C:2](=[O:3])[C:4][C:5][C:6][C:7]1']) { const p = rxn(sma, [s]); if (p.length) return finish(m, mainFrag(p[0].join('.')), { title: 'The carboxylate O⁻ displaces the halide inside the same molecule → lactone (5- or 6-ring)', text: 'Intramolecular SN2: a γ- (or δ-) halo carboxylate closes a five- (or six-) membered cyclic ester faster than any intermolecular reaction.', tags: ['sub'] }); } return null; } });
front({ id: 'dibal_cocl', name: 'DIBAL-H / LiAlH(OtBu)₃ on an acyl chloride → aldehyde', fam: 'Acid derivative', canon: 'DIBAL-H, −78 °C', db: [],
  test: c => c.any(['DIBAL']) || /lialh\s*\(\s*o-?t-?bu\s*\)\s*3/i.test(c.raw || ''),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX3](=O)Cl').length) return null; const p = rxn('[#6:1][CX3:2](=[O:3])Cl>>[*:1][CH1:2]=[O:3]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'A mild, bulky hydride reduces the acyl chloride only to the aldehyde', text: 'One hydride is delivered at low temperature; the aldehyde survives because the reagent is too weak / too hindered to reduce it further at −78 °C.', tags: ['red'] }); } });
front({ id: 'cyc_anhydride', name: '1,2-Dicarboxylic acid → cyclic anhydride', fam: 'Acid derivative', canon: '(CH₃CO)₂O or P₂O₅, heat', db: [],
  test: c => c.any(['Ac2O', 'P2O5']) || /\(ch3co\)2o|acetic\s+anhydride|p2o5|p4o10/i.test(c.raw || ''),
  run(m) { const s = smiOf(m); if (!s) return null; for (const sma of ['[CX3:1](=[O:2])([OX2H1:3])[c:4]:[c:5][CX3:6](=[O:7])[OX2H1]>>[C:1]1(=[O:2])[O:3][C:6](=[O:7])[c:5]:[c:4]1', '[CX3:1](=[O:2])([OX2H1:3])[CX4:4][CX4:5][CX3:6](=[O:7])[OX2H1]>>[C:1]1(=[O:2])[O:3][C:6](=[O:7])[C:5][C:4]1']) { const p = rxn(sma, [s]); if (p.length) return finish(m, mainFrag(p[0].join('.')), { title: 'The two neighbouring –COOH groups lose H₂O → cyclic anhydride (5-ring)', text: 'Phthalic or succinic acid closes to the five-membered anhydride on heating with a dehydrating agent (acetic anhydride, P₂O₅) or on its own above the melting point.', tags: ['dehyd'] }); } return null; } });
/* tests: Na / NaHCO3 gas evolution are tests */
wrap('alc_na', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (o && o.product) o.tags = (o.tags || []).concat(['test']); return o; });
wrap('acidbase', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (o && o.product && /effervescence|CO₂/.test(o.title || '')) o.tags = (o.tags || []).concat(['test']); return o; });
/* ---- run 35, pass 2: fixes found by the probe ---- */
/* NCS works like NBS (radical, allylic / benzylic), with Cl */
{ const r = recipe('ncs'); if (r) r.run = function (m, c, all, k) { const rx = recipe('radx'); if (!rx) return null; const s0 = smiOf(m); if (!s0 || /I/.test(s0.replace(/Ir|In/g, ''))) return null; const held = /Br/.test(s0); const m1 = held ? fresh(s0.replace(/Br/g, 'I')) : m; if (!m1) return null;
    const c2 = ctxWith(c, { NBS: {}, NCS: null }, ['light']); let o = null; try { o = rx.run(m1, c2, all, k); } catch (e) { o = null; } if (!o || !o.product) return null; const ps = smiOf(o.product); if (!ps || !/Br/.test(ps)) return null;
    return finish(m, held ? ps.replace(/Br/g, 'Cl').replace(/I/g, 'Br') : ps.replace(/Br/g, 'Cl'), { title: 'NCS, hν: allylic / benzylic chlorination (as NBS does with Br)', text: 'N-Chlorosuccinimide keeps the Cl₂ level tiny, so Cl• abstracts the weak allylic / benzylic H (resonance-stabilised radical) instead of adding to the C=C.', tags: ['radical'] }); }; }
/* nitrile + Grignard: the imine salt is the end of the stage, extra equivalents do nothing */
wrap('grignard_add', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (o && o.product && /nitrile → ketone/.test(o.title || '') && M(m, '[CX2]#[NX1]').length) o.stop = true; return o; });
/* nitration never goes past three NO2 on an activated ring, or two on benzene / a deactivated ring (also across repeats) */
const NO2 = '[c][N+](=O)[O-]';
const ringAct = P => M(P, 'c-[CX4,OX2H1,$([OX2][#6]),$([NX3;!$([N+]);!$(N=O)])]').length > 0;
wrap('eas', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || !/Nitration|nitration|TNT|dinitro/.test(o.title || '')) return o;
  const before = M(m, NO2).length, after = M(o.product, NO2).length, max = ringAct(m) ? 3 : 2; if (after <= max) return o; if (before >= max) return null; return o; });
/* alcohol + HX made in the flask (NaBr / H2SO4, HI, Lucas): substitution, not dehydration */
for (const id of ['dehyd', 'acid_heat_dehyd']) wrap(id, function (m, c) { if (c.any(['HBr', 'HI', 'Lucas']) || (c.has('HCl') && c.has('ZnCl2'))) return null; return undefined; });
/* Tollens / Fehling present: the aldehyde is oxidised, not turned into an aldehyde-ammonia */
wrap('ald_nh3', function (m, c) { if (c.any(['Tollens', 'AgNO3', 'Fehling', 'CuNH3'])) return null; return undefined; });
/* periodate: every C–C between two O-bearing carbons is cut; a carbon cut twice ends up as –COOH (glycerol → 2 HCHO + HCOOH) */
const NOTACID = '!$(C(=O)[O,N,S,Cl,Br,I])';
front({ id: 'hio4_all', name: 'HIO₄: cut every C–C between O-bearing carbons', fam: 'Oxidation', canon: 'HIO₄', db: ['Malaprade (periodate) cleavage'],
  test: c => c.has('HIO4'),
  run(m) { const s = smiOf(m); if (!s || M(m, '[CX4;R]([OX2H1])@[CX4;R][OX2H1]').length || M(m, '[#6;R]([OX2H1,OX1])@[#6;R]=O').length) return null; const T = [
      ['[CX4:1]([OX2H1:2])-[CX4:3][OX2H1:4]>>[C:1]=[O:2].[C:3]=[O:4]', '1,2-diol'],
      ['[CX3;' + NOTACID + ':1](=[O:2])-[CX4:3][OX2H1:4]>>[C:1](=[O:2])O.[C:3]=[O:4]', 'α-hydroxy carbonyl'],
      ['[CX3;' + NOTACID + ':1](=[O:2])-[CX3;' + NOTACID + ':3]=[O:4]>>[C:1](=[O:2])O.[C:3](=[O:4])O', '1,2-dicarbonyl']];
    let cur = s; const pieces = []; const labels = [];
    for (let g = 0; g < 16; g++) { let hit = false; for (const [sma, lab] of T) { const fr = cur.split('.'); let done = false;
        for (let i = 0; i < fr.length && !done; i++) { const p = rxn(sma, [fr[i]]); if (!p.length) continue; const rest = fr.filter((x, j) => j !== i); cur = rest.concat(p[0].join('.').split('.')).join('.'); if (!labels.includes(lab)) labels.push(lab); done = true; }
        if (done) { hit = true; break; } } if (!hit) break; }
    if (cur === s) return null; const fr = byC(cur.split('.').filter(Boolean)); const main = fr[0]; const ex = [...new Set(fr.slice(1))];
    const nC = ex.length + 1; if (labels.length === 1 && labels[0] === '1,2-diol' && nC <= 2 && !ex.length) { /* single ring diol: one dialdehyde */ }
    return finish(m, main, { title: 'HIO₄ cuts every C–C bond between two O-bearing carbons (' + labels.join(', ') + ')', text: 'Each cut gives each of the two carbons one more bond to O: CH₂OH → HCHO, a CHOH cut on both sides → HCOOH, a C=O next to the cut → –COOH. One HIO₄ is used per C–C bond broken.', extra: ex, tags: ['oxidation'] }); } });
/* conc. NaOH alone (Cannizzaro conditions) does nothing to a single ketone: no –CHO for Cannizzaro, and the ketone self-aldol lies back on the ketone side */
wrap('aldol', function (m, c) { if (!c.any(['KOH', 'aqKOH'])) return undefined; const conc = c.flag('conc') || /\bconc|concentrated/i.test(c.raw || ''); if (!conc || heat(c) || c.flag('dil')) return undefined;
  if (m.frags().length !== 1) return undefined; if (M(m, '[CX3H1](=O)').length || M(m, '[CH2]=O').length) return undefined; if (M(m, '[CX3](=O)([#6])[#6]').length !== 1) return undefined;
  return null; });
/* ---- run 36 ---- */
/* Wolff–Kishner (and other N-nucleophiles) on an enal: the Cannizzaro guard is only for Cannizzaro conditions */
wrap('cann_guard', function (m, c) { if (c.any(['NH2NH2', 'semicarbazide', 'NH2OH', 'PhNHNH2', 'DNP', 'NaCN', 'HCN'])) return null; return undefined; });
/* HI / red P: every C–O becomes C–H but an aromatic ring stays aromatic; ethers are reduced, not just cleaved */
wrap('hi_redp', function (m, c, all, k, run0) { if (!M(m, 'c').length && m.frags().length === 1 && !M(m, '[#6][OX2][#6]').length) return undefined; const s = smiOf(m); if (!s) return undefined;
  const r = applyAll(s, [['[CX3:1](=O)[OX2H1]>>[CH3:1]', 'COOH → CH₃'], ['[CX3:1]=[OX1]>>[C:1]', 'C=O → CH₂'], ['[#6:1][OX2H1]>>[*:1]', 'C–OH → C–H'], ['[#6:1][OX2][#6:2]>>[*:1].[*:2]', 'C–O–C → 2 C–H']], 8); if (!r || r.cur === s || /O/.test(r.cur)) return undefined;
  return finish(m, r.cur, { title: 'HI / red P, Δ: every C–O reduced to C–H (an aromatic ring is kept)', text: 'HI replaces each –OH (or ether O) by I, and HI/P then reduces each C–I to C–H; C=O goes to CH₂. Benzene rings are not hydrogenated.', extra: [...new Set(r.extra)], tags: ['reduction'] }); });
wrap('ether_hi', function (m, c) { if (c.has('redP')) return null; return undefined; });
/* acyl chloride + KCN → acyl cyanide (CH3COCl → CH3COCN; H3O+ then gives pyruvic acid) */
front({ id: 'acyl_cn', name: 'Acyl chloride + KCN → acyl cyanide', fam: 'Acid derivative', canon: 'KCN', db: [],
  test: c => c.has('NaCN') && !c.any(['HCN', 'H3O+', 'KOH', 'aqKOH', 'LiAlH4']),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX3](=O)Cl').length) return null; const p = rxn('[#6:3][CX3:1](=[O:2])Cl>>[*:3][C:1](=[O:2])C#N', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Acyl chloride + KCN → acyl cyanide (R–CO–CN)', text: 'CN⁻ adds to the acyl C=O and Cl⁻ leaves (addition–elimination). Hydrolysis of CH₃COCN gives pyruvic acid, CH₃COCOOH.', tags: ['sub'] }); } });
/* SeO2 on a substrate with no C=O: allylic C–H → allylic alcohol, benzylic CH3 → ArCHO */
front({ id: 'seo2_allyl', name: 'SeO₂: allylic / benzylic oxidation', fam: 'Oxidation', canon: 'SeO₂', db: [],
  test: c => c.has('SeO2'),
  run(m) { const s = smiOf(m); if (!s || M(m, '[CX3]=O').length) return null; let p = rxn('[C:1]=[C:2][CX4;H2,H3:3]>>[C:1]=[C:2][C:3]O', [s]), t = 'SeO₂: allylic C–H → C–OH (allylic oxidation; more SeO₂ gives the enal)';
    if (!p.length) { p = rxn('[c:1][CH3:2]>>[c:1][CH1:2]=O', [s]); t = 'SeO₂: benzylic CH₃ → CHO'; } if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: t, text: 'SeO₂ does an ene reaction with the C=C (or attacks the benzylic C–H), then a [2,3]-shift and hydrolysis put O on the carbon next to the double bond / ring.', tags: ['ox'] }); } });
/* ---- run 36: curly-arrow mechanisms for steps that used to morph (polar recipes) ---- */
const Mech = E.Mech || (root.EngineCore && root.EngineCore.Mech);
const neutK = mol => { try { return String(C.canon(mol)).replace(/[@\/\\]/g, '').replace(/\[O-\]/g, 'O').replace(/\[S-\]/g, 'S').replace(/\[NH3\+\]/g, 'N').replace(/\[NH2\+\]/g, 'N').replace(/\[NH\+\]/g, 'N').replace(/\[OH\+\]/g, 'O'); } catch (e) { return null; } };
const sameMain = (a, b) => { const ka = neutK(a), kb = neutK(b); if (!ka || !kb) return false; const fa = ka.split('.').map(x => { try { return C.canon(C.fromSmilesMB(x, BL)); } catch (e) { return x; } }); const fb = kb.split('.').map(x => { try { return C.canon(C.fromSmilesMB(x, BL)); } catch (e) { return x; } }); return fa.some(x => fb.includes(x)); };
const arrowless = o => o && o.steps && o.steps.length > 0 && o.steps.every(st => !(st.arrows || []).length);
const spOf = def => { try { return E.species(def); } catch (e) { return null; } };
const atomWhere = (sp, f) => { for (const a of sp.atoms.values()) if (f(a)) return a.id; return null; };
function baseDef(c) { return c.has('NaH') ? '[H-]' : c.has('NaNH2') ? '[NH2-]' : c.has('NaOEt') ? '[O-]CC' : c.has('NaOMe') ? '[O-]C' : c.has('KOtBu') ? '[O-]C(C)(C)C' : c.any(['KOH', 'aqKOH', 'alcKOH']) ? '[OH-]' : c.has('NaHCO3') ? '[O-]C(=O)O' : null; }
function addBase(mech, near, c, def0) { const sp = spOf(def0 || baseDef(c) || '[OH-]'); if (!sp) return null; const b = atomWhere(sp, a => (a.chg || 0) < 0); const ids = mech.add(sp, near, { att: b, dist: 1.55 * BL }); return ids.get(b); }
function protonOff(mech, at, c, text, def0) { const h = mech.H(at); if (h == null) return false; const b = addBase(mech, h, c, def0); if (b == null) return false; mech.step([[[b], [b, h]], [[h, at], [at]]], text); return true; }
function rxPartner(c) { for (const k of ['RX', 'partner', 'org']) { const v = c.get(k); if (!v) continue; let smi = null; if (k === 'RX' && v.smi && v.X) smi = v.X + v.smi; else if (v.smi && /Cl|Br|I/.test(v.smi)) smi = v.smi; if (!smi) continue;
    const sp = spOf(smi); if (!sp) continue; const x = atomWhere(sp, a => ['Cl', 'Br', 'I'].includes(a.el) && sp.nbrs(a.id).length === 1); if (x == null) continue; const cc = sp.nbrs(x)[0]; if (sp.atoms.get(cc).el !== 'C') continue; return { sp, cAt: cc, xAt: x }; } return null; }
function sn2Partner(mech, nu, rx, text) { const ids = mech.add(rx.sp, nu, { att: rx.cAt, dist: 1.75 * BL }); const cc = ids.get(rx.cAt), x = ids.get(rx.xAt); mech.step([[[nu], [nu, cc]], [[cc, x], [x]]], text); return true; }
function useMech(o, mech) { if (!mech || !mech.steps.length) return o; let P = null; try { P = mech.product(); } catch (e) { P = null; } if (!P || !sameMain(P, o.product)) return o; return Object.assign({}, o, { steps: mech.steps, product: P }); }
function arrowize(id, build) { wrap(id, function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || o.none || !arrowless(o)) return o; let mech = null; try { mech = build(m, c, o); } catch (e) { mech = null; } return mech ? useMech(o, mech) : o; }); }
/* alkoxide / phenoxide / thiolate + R–X (Williamson): deprotonate if needed, then SN2 */
const nuDefSN2 = c => { const v = c.get('NaOR'); if (v && v.smi) return '[O-]' + v.smi; return c.has('NaOEt') ? '[O-]CC' : c.has('NaOMe') ? '[O-]C' : c.has('KOtBu') ? '[O-]C(C)(C)C' : c.has('NaCN') ? '[C-]#N' : c.has('NaI') ? '[I-]' : c.has('NaSH') ? '[SH-]' : c.get('RS-') ? '[S-]C' : c.has('NaN3') ? '[N-]=[N+]=[N-]' : null; };
function sn2Reagent(m, c, txt) { const hx = M(m, '[CX4:1][Cl,Br,I:2]'); const nd = nuDefSN2(c); if (!hx.length || !nd) return null; const [cc, x] = hx[0]; const sp = spOf(nd); if (!sp) return null; const na = atomWhere(sp, a => (a.chg || 0) < 0); const mech = new Mech(m); const ids = mech.add(sp, cc, { att: na, dist: 1.75 * BL }); const nu = ids.get(na);
  mech.step([[[nu], [nu, cc]], [[cc, x], [x]]], txt || 'SN2: the nucleophile attacks the carbon from the back side as the halide leaves (one step, inversion).'); return mech; }
const nuSN2 = (m, c, txt) => { const rx = rxPartner(c); if (!rx && !M(m, '[OX1-,SX1-]').length) return sn2Reagent(m, c, txt); if (!rx) { const q0 = M(m, '[OX1-,SX1-]'); const hx = M(m, '[CX4:1][Cl,Br,I:2]'); if (!q0.length || !hx.length) return null; const mech0 = new Mech(m); mech0.step([[[q0[0][0]], [q0[0][0], hx[0][0]]], [[hx[0][0], hx[0][1]], [hx[0][1]]]], txt || 'SN2: the O⁻ attacks the carbon of R–X from the back side while the halide leaves (one step).'); return mech0; } const q = M(m, '[OX1-,SX1-]'); const q2 = q.length ? q : M(m, '[OX2H1,SX2H1;$(*[#6])]'); if (!q2.length) return null; const nu = q2[0][0]; const mech = new Mech(m);
  if (!q.length) { if (!protonOff(mech, nu, c, 'The base removes the O–H (or S–H) proton: the anion is the nucleophile.')) return null; }
  sn2Partner(mech, nu, rx, txt || 'SN2: the O⁻ attacks the carbon of R–X from the back side while the halide leaves (one step, inversion at that carbon).'); return mech; };
for (const id of ['williamson', 'rx_alkoxide', 'thiol_alk']) arrowize(id, (m, c) => nuSN2(m, c, id === 'thiol_alk' ? 'SN2: RS⁻ attacks the carbon of R–X from the back side; the halide leaves.' : null));
/* acyl chloride + OH⁻ / CN⁻: addition–elimination */
function acylSub(m, nuDef, txt1) { const q = M(m, '[CX3:1](=[O:2])[Cl,Br:3]'); if (!q.length) return null; const [cc, o2, x] = q[0]; const mech = new Mech(m); const sp = spOf(nuDef); if (!sp) return null; const na = atomWhere(sp, a => (a.chg || 0) < 0); const ids = mech.add(sp, cc, { att: na, dist: 1.7 * BL }); const nu = ids.get(na);
  mech.step([[[nu], [nu, cc]], [[cc, o2], [o2]]], txt1); mech.step([[[o2], [o2, cc]], [[cc, x], [x]]], 'The C=O re-forms and the halide leaves (addition–elimination).'); return { mech, nu }; }
arrowize('acylx_oh', (m, c) => { const r = acylSub(m, '[OH-]', 'OH⁻ adds to the C=O carbon; the π electrons move onto O (tetrahedral intermediate).'); if (!r) return null; protonOff(r.mech, r.nu, c, 'In base the acid at once loses H⁺ → carboxylate.', '[OH-]'); return r.mech; });
arrowize('acyl_cn', (m) => { const r = acylSub(m, '[C-]#N', 'CN⁻ adds to the C=O carbon; the π electrons move onto O.'); return r && r.mech; });
/* epoxides */
const ringC = (m, o) => m.nbrs(o).filter(n => el(m, n) === 'C');
arrowize('epox_oh', (m, c) => { const q = M(m, '[CX4;r3:1]1[OX2;r3:2][CX4;r3:3]1'); if (!q.length) return null; const [a0, o, b0] = q[0]; const [a, b] = nCn(m, a0) <= nCn(m, b0) ? [a0, b0] : [b0, a0]; const mech = new Mech(m);
  const sp = spOf('[OH-]'); const ids = mech.add(sp, a, { att: 1, dist: 1.7 * BL }); const nu = ids.get(1); mech.step([[[nu], [nu, a]], [[a, o], [o]]], 'OH⁻ attacks the less hindered ring carbon from the back side (SN2); the strained C–O bond breaks.');
  const w = spOf('[H]O[H]'); const wi = mech.add(w, o, { att: 1, dist: 1.5 * BL }); mech.step([[[o], [o, wi.get(1)]], [[wi.get(1), wi.get(2)], [wi.get(2)]]], 'The alkoxide takes H⁺ from water: the 1,2-diol.'); return mech; });
arrowize('epox_h3o', (m, c) => { const q = M(m, '[CX4;r3:1]1[OX2;r3:2][CX4;r3:3]1'); if (!q.length) return null; const [a0, o, b0] = q[0]; const [b, a] = nCn(m, a0) >= nCn(m, b0) ? [a0, b0] : [b0, a0]; const mech = new Mech(m);
  const hp = spOf('[H][OH2+]'); const hi = mech.add(hp, o, { att: 1, dist: 1.5 * BL }); mech.step([[[o], [o, hi.get(1)]], [[hi.get(1), hi.get(2)], [hi.get(2)]]], 'The epoxide O takes H⁺ from H₃O⁺ (now a good leaving group).');
  const w = spOf('[H]O[H]'); const wi = mech.add(w, b, { att: 2, dist: 1.7 * BL }); const ow = wi.get(2); mech.step([[[ow], [ow, b]], [[b, o], [o]]], 'Water attacks the more substituted carbon (it carries more + charge) from the back side; the C–O⁺ bond breaks.');
  E.deprotOn(mech, ow, 'Water removes H⁺: the trans-1,2-diol.'); return mech; });
/* acid–base proton transfers */
arrowize('acid_alkoxide', (m, c) => { const q = M(m, '[CX3](=O)[OX2H1:1]'); if (!q.length) return null; const mech = new Mech(m); return protonOff(mech, q[0][2] != null ? q[0][2] : q[0][q[0].length - 1], c, 'The alkoxide (a strong base) takes the acidic O–H proton: carboxylate + alcohol. No ester can form.') ? mech : null; });
arrowize('acidic_ch_sh', (m, c) => { let q = M(m, '[SX2H1][#6]'); let at = q.length ? q[0][0] : null; if (at == null) { q = M(m, '[CH2;r5]1C=CC=C1'); at = q.length ? q[0][0] : null; } if (at == null) return null; const mech = new Mech(m); return protonOff(mech, at, c, 'The base removes the acidic proton; the anion is stabilised (big polarisable S, or an aromatic 6π cyclopentadienyl ring).') ? mech : null; });
/* interhalogen: electrophilic end adds first (cation at the more substituted C), then the halide */
arrowize('interhalogen', (m, c) => { const q = M(m, 'C=C'); if (!q.length) return null; const def = c.has('BrCl') ? 'BrCl' : c.has('ICl') ? 'ICl' : 'IBr'; const mech = new Mech(m); const r = E.protonateAlkene(mech, q[0][0], q[0][1], def, 'The π bond attacks the less electronegative halogen (δ+ end); the other halogen leaves as X⁻. The cation forms on the more substituted carbon.');
  E.captureBy(mech, r.cp, r.X, 'The halide ion captures the carbocation (Markovnikov orientation).'); return mech; });
/* vic-dihalide + I⁻ → alkene */
arrowize('vic_nai', (m, c) => { const q = M(m, '[Cl,Br,I:3][CX4:1][CX4:2][Cl,Br,I:4]'); if (!q.length) return null; const [x1, c1, c2, x2] = q[0]; const mech = new Mech(m); const sp = spOf('[I-]'); const ids = mech.add(sp, x1, { att: 1, dist: 1.6 * BL }); const io = ids.get(1);
  mech.step([[[io], [io, x1]], [[x1, c1], [c1, c2]], [[c2, x2], [x2]]], 'I⁻ attacks one halogen; the C=C forms as the other halide leaves (anti elimination).'); return mech; });
/* conjugated diene + HX / X2: cation at C-2, allylic shift to C-4, capture */
arrowize('diene_add', (m, c, o) => { const q = M(m, '[C:1]=[C:2][C:3]=[C:4]'); if (!q.length) return null; const def = c.has('HBr') ? '[H]Br' : c.has('HCl') ? '[H]Cl' : c.has('HI') ? '[H]I' : c.has('Br2') ? 'BrBr' : c.has('Cl2') ? 'ClCl' : null; if (!def) return null;
  const tryV = (one4) => { const mech = new Mech(m); const [a, b, cc, d] = q[0]; const r = E.protonateAlkene(mech, a, b, def, def.startsWith('[H]') ? 'The terminal C=C takes H⁺: the cation forms at C-2, next to the other C=C (allylic, resonance-stabilised).' : 'The terminal C=C attacks Br₂: Br adds to C-1 and the allylic cation forms at C-2.');
    if (one4) { const cp = r.cp, far = cp === b ? d : a, mid = cp === b ? cc : b; mech.step([[[mid, far], [cp, mid]]], 'Resonance: the other C=C shifts, moving the + charge to C-4.'); E.captureBy(mech, far, r.X, 'The halide ion adds at C-4: the 1,4-adduct (the more stable, internal alkene).'); }
    else E.captureBy(mech, r.cp, r.X, 'The halide ion adds at C-2: the 1,2-adduct (formed faster, at low temperature).'); return mech; };
  for (const v of [true, false]) { const mech = tryV(v); let P = null; try { P = mech.product(); } catch (e) { P = null; } if (P && sameMain(P, o.product)) return mech; } return null; });
/* alcohol + HX with a 1,2-shift: protonate, lose water, shift, capture */
arrowize('alc_hx_shift', (m, c) => { const q = M(m, '[CX4:1][OX2H1:2]'); if (!q.length) return null; const [cc, o] = q[0]; const def = c.has('HBr') ? '[H]Br' : c.has('HI') ? '[H]I' : '[H]Cl'; const mech = new Mech(m); const sp = spOf(def); const ids = mech.add(sp, o, { att: 1, dist: 1.5 * BL }); const h = ids.get(1), x = ids.get(2);
  mech.step([[[o], [o, h]], [[h, x], [x]]], 'The –OH takes H⁺ from HX: –OH₂⁺ is a good leaving group.'); mech.step([[[cc, o], [o]]], 'Water leaves: a carbocation.');
  const sh = E.bestShift(mech.m, cc); let cp = cc; if (sh) cp = E.doShift(mech, cc, sh); E.captureBy(mech, cp, x, 'The halide ion captures the more stable (3°) carbocation.'); return mech; });
/* chloral hydrate */
arrowize('chloral_h2o', (m) => { const q = M(m, '[CX3H1:1](=[O:2])C(Cl)(Cl)Cl'); if (!q.length) return null; const [cc, o] = q[0]; const mech = new Mech(m); const w = spOf('[H]O[H]'); const wi = mech.add(w, cc, { att: 2, dist: 1.7 * BL }); const ow = wi.get(2), h = wi.get(1);
  mech.step([[[ow], [ow, cc]], [[cc, o], [o]]], 'Water adds to the very electron-poor C=O (the CCl₃ group pulls electrons away).'); mech.step([[[o], [o, h]], [[h, ow], [ow]]], 'Proton transfer from O⁺ to O⁻: the gem-diol, chloral hydrate.'); return mech; });
/* intramolecular SN2: ω-halo carboxylate → lactone */
arrowize('halo_acid_lactone', (m, c) => { const q = M(m, '[OX2H1,OX1-:1][CX3](=O)[CX4][CX4][CX4;$([CX4][Cl,Br,I]),$([CX4][CX4][Cl,Br,I])]'); if (!q.length) return null; const o = q[0][0]; const cx = M(m, '[CX4:1][Cl,Br,I:2]'); if (!cx.length) return null; const mech = new Mech(m);
  if (M(m, '[OX2H1][CX3]=O').length) { if (!protonOff(mech, o, c, 'The base removes the acidic O–H: carboxylate.')) return null; }
  mech.step([[[o], [o, cx[0][0]]], [[cx[0][0], cx[0][1]], [cx[0][1]]]], 'The carboxylate O⁻ attacks the C–X carbon of the same chain (intramolecular SN2): a five- or six-membered lactone closes.'); return mech; });
/* ---- run 36: ring substitutions get arrows too (EAS without a catalyst, SNAr) ---- */
function sitesFor(m, P, X, n) { const cand = M(m, '[cH1]').map(x => x[0]); const want = neutK(P); if (!want || n < 1 || n > 3) return null;
  const combos = []; const rec = (st, i) => { if (st.length === n) { combos.push(st.slice()); return; } for (let j = i; j < cand.length; j++) { st.push(cand[j]); rec(st, j + 1); st.pop(); } }; rec([], 0);
  for (const cb of combos) { const t = m.clone(); try { for (const id of cb) E.addOn(t, id, X); } catch (e) { continue; } const k = neutK(t); if (k && (k === want || want.split('.').includes(k))) return cb; } return null; }
function easMulti(m, o, X, eDef, eTxt, baseDef0) { const n = M(o.product, 'c' + (X === 'NO2' ? '[N+](=O)[O-]' : X)).length - M(m, 'c' + (X === 'NO2' ? '[N+](=O)[O-]' : X)).length; const sites = sitesFor(m, o.product, X === 'NO2' ? '[N+](=O)[O-]' : X, n); if (!sites) return null;
  const mech = new Mech(m); let i = 0;
  for (const p of sites) { const q = mech.m.nbrs(p).find(nb => mech.m.bondBetween(p, nb).o === 2 && el(mech.m, nb) === 'C'); if (q == null) return null; const sp = spOf(eDef); if (!sp) return null; const ea = atomWhere(sp, a => a.el === (X === 'NO2' ? 'N' : X)); const ids = mech.add(sp, p, { att: ea, dist: 1.7 * BL }); const e1 = ids.get(ea);
    const lv = X === 'NO2' ? null : atomWhere(sp, a => a.id !== ea && a.el === X); const arr = [[[p, q], [p, e1]]]; if (lv != null) arr.push([[e1, ids.get(lv)], [ids.get(lv)]]); if (X === 'NO2') { const od = atomWhere(sp, a => a.el === 'O' && sp.bondBetween && sp.bondBetween(ea, a.id) && sp.bondBetween(ea, a.id).o === 2); if (od != null) arr.push([[e1, ids.get(od)], [ids.get(od)]]); }
    mech.step(arr, i === 0 ? eTxt : 'Again at the next activated position: the ring π electrons attack the electrophile (σ-complex).'); if (lv != null) mech.dropFrag(ids.get(lv));
    const h = mech.H(p); if (h == null) return null; const b = addBase(mech, h, null, baseDef0 || '[Br-]'); if (b == null) return null; mech.step([[[b], [b, h]], [[p, h], [p, q]]], 'A base removes H⁺ from the sp³ ring carbon; the ring is aromatic again.'); mech.dropFrag(b); i++; }
  return mech; }
const brTxt = 'The very electron-rich ring polarises Br–Br and attacks it directly (no Lewis acid needed): a σ-complex; Br⁻ leaves.';
arrowize('br2_act', (m, c, o) => easMulti(m, o, c.has('Cl2') ? 'Cl' : 'Br', c.has('Cl2') ? 'ClCl' : 'BrBr', brTxt, c.has('HOBr') ? 'O' : null));
arrowize('polycyc_x2', (m, c, o) => easMulti(m, o, c.has('Br2') ? 'Br' : 'Cl', c.has('Br2') ? 'BrBr' : 'ClCl', 'The π system attacks Br₂ at the position whose σ-complex keeps the most aromatic stabilisation (C-2 of a five-membered heteroaromatic, C-1 of naphthalene).'));
arrowize('eas', (m, c, o) => { const t = o.title || ''; if (/Nitration|nitration|TNT|dinitro|picric/i.test(t)) return easMulti(m, o, 'NO2', '[N+](=O)=O', 'The ring π electrons attack the nitronium ion NO₂⁺ (from HNO₃ + H₂SO₄): a σ-complex.', '[O-]S(=O)(=O)O'); if (/Bromination|bromination/.test(t)) return easMulti(m, o, 'Br', 'BrBr', 'Br₂ is polarised by the Lewis acid (Br–Br–FeBr₃); the ring attacks the δ+ Br.'); if (/Chlorination|chlorination/.test(t)) return easMulti(m, o, 'Cl', 'ClCl', 'Cl₂ is polarised by the Lewis acid; the ring attacks the δ+ Cl.'); return null; });
/* SNAr: addition (Meisenheimer complex) then loss of X⁻ */
const nuDefOf = c => c.has('NaOMe') ? '[O-]C' : c.has('NaOEt') ? '[O-]CC' : (c.get('RS-') ? '[S-]' + (((c.get('RS-') || {}).txt || '').match(/c2h5|et/) ? 'CC' : 'C') : null) || (c.has('NaSH') ? '[SH-]' : null) || (c.any(['KOH', 'aqKOH']) ? '[OH-]' : null) || (c.has('NH3') ? 'N' : null) || (c.has('H2O') ? 'O' : null);
arrowize('snar_nu', (m, c, o) => { const nd = nuDefOf(c); if (!nd) return null; const xs = M(m, 'c[F,Cl,Br,I]'); for (const [ip, x] of xs) { const mech = new Mech(m); const q = mech.m.nbrs(ip).find(nb => mech.m.bondBetween(ip, nb).o === 2 && E.aromSet(mech.m).has(nb)); if (q == null) continue; const sp = spOf(nd); if (!sp) return null; const na = atomWhere(sp, a => (a.chg || 0) < 0) || atomWhere(sp, a => ['N', 'O', 'S'].includes(a.el)); const ids = mech.add(sp, ip, { att: na, dist: 1.7 * BL }); const nu = ids.get(na);
    mech.step([[[nu], [nu, ip]], [[ip, q], [q]]], 'The nucleophile adds to the carbon bearing the halogen; the negative charge spreads into the ring and onto the o/p-NO₂ groups (Meisenheimer complex).'); mech.step([[[q], [q, ip]], [[ip, x], [x]]], 'The ring re-aromatises and the halide leaves (addition–elimination, SNAr).');
    if ((sp.atoms.get(na).chg || 0) === 0) E.deprotOn(mech, nu, 'A base removes H⁺ from the new N–H / O–H.'); let P = null; try { P = mech.product(); } catch (e) { P = null; } if (P && sameMain(P, o.product)) return mech; } return null; });
/* ---- run 36, pass 2: general rules from the SKM Advanced misses ---- */
/* alkyl fluorides do E2 with strong base (poor leaving group → Hofmann, less substituted alkene) */
front({ id: 'rf_e2', name: 'Alkyl fluoride + strong base → Hofmann alkene', fam: 'Alkyl halide', canon: 'alc. KOH, heat', db: [],
  test: c => c.any(['alcKOH', 'NaOEt', 'NaOMe', 'KOtBu', 'NaNH2']) || (c.has('KOH') && heat(c) && !c.has('aqKOH')),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX4][F]').length || M(m, '[CX4][Cl,Br,I]').length) return null; const ps = rxn('[CX4;!H0:2][CX4:1]F>>[C:2]=[C:1]', [s]).map(x => mainFrag(x.join('.'))).filter(Boolean); if (!ps.length) return null;
    const sc = smi => { const P = fresh(smi); if (!P) return 99; const q = M(P, 'C=C'); return q.length ? nCn(P, q[0][0]) + nCn(P, q[0][1]) : 99; }; const best = ps.slice().sort((a, b) => sc(a) - sc(b))[0];
    return finish(m, best, { title: 'R–F + strong base: E2 gives the less substituted (Hofmann) alkene', text: 'F⁻ is a poor leaving group, so the E2 transition state is carbanion-like: the most acidic, least hindered β-H (on CH₃/CH₂) is removed.', tags: ['elim'] }); } });
/* Strecker: NH3 (or NH4Cl) + HCN/KCN with a carbonyl → α-amino nitrile (not the cyanohydrin) */
wrap('hcn', function (m, c) { if (c.any(['NH3', 'NH4Cl', 'RNH2'])) return null; return undefined; });
front({ id: 'strecker2', name: 'Strecker: C=O + NH₃ + HCN → α-amino nitrile', fam: 'Carbonyl', canon: 'NH₃, HCN', db: ['Strecker synthesis'],
  test: c => c.any(['NH3', 'NH4Cl']) && c.any(['HCN', 'NaCN']),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX3;$([CX3H1](=O)[#6]),$([CX3](=O)([#6])[#6])]=O').length) return null; const p = rxn('[CX3:1]=[OX1]>>[C:1](N)C#N', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Strecker: C=O + NH₃ + HCN → α-amino nitrile (acid hydrolysis then gives the α-amino acid)', text: 'NH₃ condenses with the C=O to an imine; CN⁻ adds to the C=N. Hydrolysis of –C≡N later gives –COOH.', tags: ['add'] }); } });
/* every OH reacts with Cu (573 K) / chromic acid / PCC, not just the first one */
front({ id: 'ox_all_oh', name: 'Oxidant / Cu on a molecule with several OH groups', fam: 'Alcohol', canon: 'Cu, 573 K / CrO₃', db: [],
  test: c => c.any(['K2Cr2O7', 'KMnO4', 'PCC']) || (c.has('Cu') && (heat(c) || (c.T && c.T > 500))),
  run(m, c) { const s = smiOf(m); const oh = M(m, '[CX4][OX2H1]'); if (!s || oh.length < 2 || M(m, 'C=C').length || M(m, 'C#C').length || M(m, '[CX3H1]=O').length || sugarLike(m)) return null; const prim = M(m, '[CX4;H1,H2]([#6])[OX2H1]').length; if (!prim) return null;
    const cu = c.has('Cu'), mild = c.has('PCC'); const T = [['[CX4H1:1]([#6:3])([#6:4])[OX2H1:2]>>[C:1]([*:3])([*:4])=[O:2]', '2° OH → C=O'], [mild || cu ? '[CX4H2:1]([#6:3])[OX2H1:2]>>[C:1]([*:3])=[O:2]' : '[CX4H2:1]([#6:3])[OX2H1:2]>>[C:1]([*:3])(=[O:2])O', mild || cu ? '1° OH → CHO' : '1° OH → COOH']];
    if (cu) T.push(['[CX4H0:1]([#6:5])([#6:6])([OX2H1])[CX4;!H0:2]>>[C:1]([*:5])([*:6])=[C:2]', '3° OH → alkene']);
    const r = applyAll(s, T, 6); if (!r || r.cur === s) return null;
    return finish(m, r.cur, { title: (cu ? 'Cu, 573 K: ' : mild ? 'PCC: ' : 'Chromic acid / KMnO₄: ') + r.labels.join(', '), text: cu ? 'Hot copper dehydrogenates 1° and 2° alcohols (H₂ lost) and dehydrates 3° alcohols, each –OH in the molecule independently.' : 'Each 1°/2° –OH is oxidised; a 3° –OH has no H on its carbon and survives.', tags: ['ox'] }); } });
/* acid aldol needs heat or a strong acid; Lucas reagent at room temperature does nothing to a ketone */
wrap('acid_aldol', function (m, c) { if ((c.has('Lucas') || c.has('ZnCl2')) && !heat(c)) return null; return undefined; });
/* chloroform: chloropicrin (conc. HNO3) and chloretone (acetone, KOH) */
front({ id: 'chcl3_x', name: 'CHCl₃ → chloropicrin / chloretone', fam: 'Alkyl halide', canon: 'conc. HNO₃ / acetone, KOH', db: [],
  test: c => c.has('HNO3') || (c.any(['KOH', 'aqKOH', 'alcKOH']) && (c.has('acetone') || /acetone|propanone/i.test(c.raw || ''))),
  run(m, c) { const s = smiOf(m); const fr = String(s || '').split('.'); if (!fr.includes('ClC(Cl)Cl')) return null;
    if (c.has('HNO3')) return finish(m, 'O=[N+]([O-])C(Cl)(Cl)Cl', { title: 'CHCl₃ + conc. HNO₃ → chloropicrin, CCl₃NO₂ (tear gas)', text: 'The C–H of chloroform is nitrated: Cl₃C–H + HO–NO₂ → Cl₃C–NO₂ + H₂O.', tags: ['sub'] });
    return finish(m, 'CC(C)(O)C(Cl)(Cl)Cl', { title: 'CHCl₃ + acetone, KOH → chloretone (1,1,1-trichloro-2-methylpropan-2-ol, a hypnotic)', text: 'OH⁻ removes the acidic C–H of CHCl₃; CCl₃⁻ adds to the acetone C=O and the alkoxide is protonated.', tags: ['add'] }); } });
/* acyl azide + heat → isocyanate → (water) amine: Curtius */
front({ id: 'curtius_drawn', name: 'Acyl azide, heat → isocyanate (Curtius)', fam: 'Acid derivative', canon: 'heat, H₂O', db: ['Curtius rearrangement'],
  test: c => heat(c) || c.any(['H2O', 'H3O+']),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[CX3](=O)N=[N+]=[N-]').length) return null; const wat = c.any(['H2O', 'H3O+', 'KOH', 'aqKOH']) || /h2o|water/i.test(c.raw || '');
    const p = rxn(wat ? '[#6:1][CX3](=O)N=[N+]=[N-]>>[*:1]N' : '[#6:1][CX3:2](=[O:3])N=[N+]=[N-]>>[*:1]N=[C:2]=[O:3]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: wat ? 'Curtius: R–CON₃, heat → R–N=C=O; water → R–NH₂ + CO₂' : 'Curtius: R–CON₃, heat → isocyanate R–N=C=O + N₂', text: 'N₂ leaves while R migrates from C to N (retention); the isocyanate adds water to a carbamic acid, which loses CO₂.', extra: wat ? ['O=C=O'] : ['N#N'], tags: ['rearr'] }); } });
/* K2Cr2O7 with distillation: the aldehyde boils off before it is oxidised further */
front({ id: 'ox_distil', name: '1° alcohol + K₂Cr₂O₇, distil → aldehyde', fam: 'Alcohol', canon: 'K₂Cr₂O₇, H₂SO₄, distil', db: [],
  test: c => c.has('K2Cr2O7') && /distil/i.test(c.raw || ''),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX4H2][OX2H1]').length) return null; const p = rxn('[CX4H2:1]([#6:3])[OX2H1:2]>>[C:1](=[O:2])[*:3]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'K₂Cr₂O₇ / H₂SO₄ with distillation: the aldehyde distils off as it forms (no further oxidation)', text: 'The aldehyde boils far lower than the alcohol (no H-bonding), so it leaves the hot flask before it can be oxidised to the acid.', tags: ['ox'] }); } });
/* cyclic anhydride + alcohol: mono-ester; excess alcohol with H+ (Fischer) → diester */
front({ id: 'anhyd_roh', name: 'Cyclic anhydride + alcohol', fam: 'Acid derivative', canon: 'ROH (H⁺)', db: [],
  test: c => c.any(['ROH']) || /ch3oh|methanol|etoh|ethanol|c2h5oh/i.test(c.raw || ''),
  run(m, c) { const q = M(m, '[CX3;R:1](=[O:2])[OX2;R:3][CX3;R:4]=[O:5]'); if (!q.length) return null; const [c1, , o3] = q[0]; const R = /eth|c2h5|etoh/i.test(c.raw || '') ? 2 : 1; const xs = (c.eq === Infinity || /excess/i.test(c.raw || '')) && c.any(['H3O+', 'H2SO4', 'HCl']);
    const t = m.clone(); const bd = t.bondBetween(c1, o3); if (!bd) return null; t.removeBond(bd.id); for (const id of [c1, o3]) { const a = t.atoms.get(id); if (a) a.hFix = null; }
    let on = E.addOn(t, c1, 'O'); let last = on; for (let i = 0; i < R; i++) last = E.addOn(t, last, 'C'); if (xs) { let l2 = E.addOn(t, o3, 'C'); if (R === 2) E.addOn(t, l2, 'C'); }
    const ps = smiOf(t); if (!ps) return null; const fo = finish(m, mainFrag(ps), { title: xs ? 'Cyclic anhydride + excess alcohol, H⁺ → diester' : 'Cyclic anhydride + alcohol → half ester (one –COOR, one –COOH)', text: 'The alcohol opens the anhydride (acyl substitution) to the mono-ester acid; with excess alcohol and acid the free –COOH is esterified too (Fischer).', tags: ['acyl'] }); if (fo) fo.stop = true; return fo; } });
/* a drawn alkoxide + R–X: Williamson beats Wurtz even when Na is written */
for (const id of ['wurtz', 'wurtz2']) wrap(id, function (m) { if (M(m, '[OX1-][CX4]').length && M(m, '[CX4][Cl,Br,I]').length) return null; return undefined; });
/* sulfonic acid + Na → salt + H2 */
front({ id: 'so3h_na', name: 'Sulphonic acid + Na', fam: 'Acid–base', canon: 'Na', db: [], test: c => c.has('Na') && !c.any(['NH3', 'Na/NH3', 'ROH']),
  run(m) { const s = smiOf(m); if (!s || !M(m, 'S(=O)(=O)[OX2H1]').length) return null; const p = rxn('[S:1](=[O:2])(=[O:3])[OX2H1:4]>>[S:1](=[O:2])(=[O:3])[O-:4]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'R–SO₃H + Na → R–SO₃⁻ Na⁺ + ½ H₂↑', text: 'Sulphonic acids are strong acids: sodium displaces H₂ at once.', tags: ['acidbase', 'test'] }); } });
/* aldehyde / ketone + 2 ArH, conc. H2SO4 → diarylmethane (DDT from chloral + chlorobenzene, bisphenol A) */
front({ id: 'diaryl_cond', name: 'C=O + 2 ArH (conc. H₂SO₄) → Ar₂C<', fam: 'Aromatic', canon: 'conc. H₂SO₄', db: [],
  test: c => c.any(['H2SO4', 'HCl']), multi: true,
  run(m) { const s = smiOf(m); const fr = String(s || '').split('.'); if (fr.length < 2) return null; const co = fr.find(x => M(fresh(x) || m, '[CX3;$([CX3H1](=O)C(Cl)(Cl)Cl),$([CX3](=O)([CH3])[CH3])]=O').length); const ar = fr.find(x => x !== co && /c1ccccc1|c1ccc/.test(x) && !/=O/.test(x)); if (!co || !ar) return null;
    const p = rxn('[CX3:1]=[OX1].[cH1:2]1[cH1:3][c:4]([Cl,OX2H1,CH3:9])[cH1:5][cH1:6][cH1:7]1.[cH1:12]1[cH1:13][c:14]([Cl,OX2H1,CH3:19])[cH1:15][cH1:16][cH1:17]1>>[C:1]([c:7]1[cH1:6][cH1:5][c:4]([*:9])[cH1:3][cH1:2]1)[c:17]1[cH1:16][cH1:15][c:14]([*:19])[cH1:13][cH1:12]1', [co, ar, ar]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Acid-catalysed condensation: the C=O joins two rings para to their Cl / OH (DDT from chloral, bisphenol A from acetone)', text: 'H⁺ activates the C=O; one ring attacks (Friedel–Crafts hydroxyalkylation), the OH is lost as water and the cation is attacked by a second ring.', tags: ['eas'] }); } });
/* Birch / enol ether hydrolysis with heat: the β,γ-enone moves into conjugation */
wrap('enol_ester_hyd', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || !heat(c)) return o; if (!M(o.product, '[CX3]=[CX3][CX4][CX3]=O').length || M(o.product, '[CX3]=[CX3][CX3]=O').length) return o; const s = smiOf(o.product); const p = s && rxn('[C:1]=[C:2][CX4;!H0:3][CX3:4]=[O:5]>>[C:1][C:2]=[C:3][C:4]=[O:5]', [s]); if (!p || !p.length) return o;
  const R = fresh(mainFrag(p[0].join('.'))); if (!R) return o; return Object.assign({}, o, { product: R, title: o.title + '; with heat the C=C moves into conjugation (α,β-enone)', steps: (o.steps || []).concat([step(o.product, R, 'Acid and heat move the C=C next to the C=O: the conjugated α,β-unsaturated ketone is more stable.', [])]) }); });
/* 2- and 4-halopyridines undergo SNAr (ring N stabilises the anion) */
front({ id: 'snar_pyr', name: '2-/4-Halopyridine + Nu⁻ (SNAr)', fam: 'Aromatic', canon: 'NaOCH₃, heat', db: [],
  test: c => c.any(['NaOMe', 'NaOEt', 'NH3', 'NaNH2', 'KOH', 'aqKOH', 'RS-']),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[$(c1([Cl,Br,F])nccc1),$(c1([Cl,Br,F])ccncc1)]').length) return null; const nu = c.has('NaOMe') ? 'OC' : c.has('NaOEt') ? 'OCC' : c.any(['NH3', 'NaNH2']) ? 'N' : c.get('RS-') ? 'SC' : 'O';
    let p = rxn('[c:1]([Cl,Br,F])1:[n:2]:[c:3]:[c:4]:[c:5]:[c:6]1>>[c:1](' + nu + ')1:[n:2]:[c:3]:[c:4]:[c:5]:[c:6]1', [s]); if (!p.length) p = rxn('[c:1]([Cl,Br,F])1:[c:2]:[c:3]:[n:4]:[c:5]:[c:6]1>>[c:1](' + nu + ')1:[c:2]:[c:3]:[n:4]:[c:5]:[c:6]1', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: '2- / 4-Halopyridine: SNAr (the ring N carries the negative charge of the Meisenheimer-type intermediate)', text: 'The nucleophile adds at the C–X carbon; the negative charge sits on the electronegative ring N; X⁻ then leaves. 3-Halopyridines do not react this way.', tags: ['snar'] }); } });
/* arene ozonolysis: both Kekulé forms contribute (o-xylene → glyoxal + methylglyoxal + biacetyl) */
wrap('arene_o3', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product) return o; const rings = C.aromaticRings(m); if (!rings.length) return o; const ring = rings[0]; const t = m.clone();
  for (let i = 0; i < 6; i++) { const b = t.bondBetween(ring[i], ring[(i + 1) % 6]); if (!b) return o; b.o = b.o === 2 ? 1 : 2; }
  let o2 = null; try { o2 = run0(t, c, all, k); } catch (e) { o2 = null; } if (!o2 || !o2.product) return o; const seen = new Set(); const ex = [];
  for (const x of [o.product, ...(o.extra || []), o2.product, ...(o2.extra || [])]) { let kk = null; try { kk = C.canon(x); } catch (e) { kk = null; } if (kk && !seen.has(kk)) { seen.add(kk); ex.push(x); } }
  return Object.assign({}, o, { extra: ex.filter(x => x !== o.product), why: (o.why || []).concat(['Both Kekulé structures react: the products of each set of three C=C appear (proof that benzene is a resonance hybrid)']) }); });
/* amide hydrolysis does not happen in the cold; trace acid alone does not hydrolyse an ester */
wrap('amide_rx', function (m, c) { if (c.cold || (c.T != null && c.T <= 283)) return null; return undefined; });
wrap('ester_hyd', function (m, c) { if ((c.flag('cat') || /trace/i.test(c.raw || '')) && !heat(c) && !c.has('H2O')) return null; if ((c.has('ROH') || /ch3oh|methanol|ethanol|etoh|c2h5oh/i.test(c.raw || '')) && !c.has('H2O') && !/h2o|water|aq/i.test(c.raw || '')) return null; return undefined; });
/* allene + aqueous acid → ketone (the central carbon ends up as C=O) */
front({ id: 'allene_h2o', name: 'Allene + H₃O⁺ → ketone', fam: 'Alkene', canon: 'H₂O, H⁺', db: [],
  test: c => c.any(['H3O+', 'H2SO4', 'Hg2+']) && (c.has('H2O') || c.has('H3O+') || /water|h2o/i.test(c.raw || '')),
  run(m) { const s = smiOf(m); if (!s || !M(m, 'C=C=C').length) return null; const p = rxn('[C:1]=[C:2]=[C:3]>>[C:1][C:2](=O)[C:3]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Allene + aqueous acid → ketone (H⁺ on a terminal C, water on the central C, enol → C=O)', text: 'Protonation at a terminal carbon gives a vinyl cation on the central carbon; water adds there and the enol tautomerises (propadiene → acetone).', tags: ['add'] }); } });
/* water + RMgX → R–H */
front({ id: 'rmgx_water', name: 'Grignard + water → alkane', fam: 'Organometallic', canon: 'H₂O', db: [], test: c => c.has('RMgX'),
  run(m, c) { const s = smiOf(m); if (s !== 'O') return null; const r = c.get('RMgX'); const R = r && r.smi ? r.smi : 'C'; return finish(m, R, { title: 'R–MgX + H₂O → R–H (the Grignard is destroyed by any acidic H)', text: 'The carbanion-like R takes H⁺ from water: alkane + Mg(OH)X.', tags: ['acidbase'] }); } });
/* Knoevenagel–Doebner: ArCHO + malonic acid, pyridine/piperidine, heat → cinnamic acid */
front({ id: 'knoev_doebner', name: 'Knoevenagel–Doebner condensation', fam: 'Carbonyl', canon: 'CH₂(COOH)₂, pyridine, heat', db: [], multi: true,
  test: c => /pyridine|piperidine/i.test(c.raw || '') || c.has('pyrrolidine') || c.flag('pyr'),
  run(m) { const s = smiOf(m); const fr = String(s || '').split('.'); const ald = fr.find(x => { const q = fresh(x); return q && M(q, '[CX3H1](=O)[#6]').length; }); const mal = fr.find(x => x !== ald && M(fresh(x) || m, 'OC(=O)[CH2]C(=O)O').length); if (!ald || !mal) return null;
    const p = rxn('[CX3H1:1](=O)[#6:2]>>[*:2][CH1:1]=[CH1]C(=O)O', [ald]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Knoevenagel–Doebner: ArCHO + malonic acid (base, heat) → α,β-unsaturated acid (CO₂ lost)', text: 'The base forms the malonate enolate, which adds to the C=O; dehydration and decarboxylation on heating give the (E)-acrylic acid (cinnamic acid from benzaldehyde).', extra: ['O=C=O'], tags: ['aldol'] }); } });
/* naphthalene nitration: same α / C-4 preference as bromination */
wrap('eas', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || !/Nitration/.test(o.title || '') || !M(m, 'c1ccc2ccccc2c1').length || M(m, '[N+](=O)[O-]').length) return o; const s = smiOf(m);
  const cand = rxn('[cH1:1]>>[c:1][N+](=O)[O-]', [s]).map(x => mainFrag(x.join('.'))).filter(Boolean); if (!cand.length) return o; const best = cand.map(x => [x, naphScore(fresh(x), 'c[N+](=O)[O-]')]).sort((a, b) => b[1] - a[1])[0][0];
  if (best === smiOf(o.product)) return o; return finish(m, best, { title: 'Nitration of the naphthalene: at C-1 (α), or C-4 when C-1 carries an activating group', text: 'Attack at an α-carbon gives the σ-complex that keeps one ring fully aromatic; an –OCH₃ / –CH₃ at C-1 directs to C-4 (para).', tags: ['eas'] }); });
wrap('aldol', function (m) { if (m.frags().length !== 1 || M(m, '[CX3H1]=O').length || M(m, '[CX3]=O').length !== 1) return undefined; if (M(m, '[CX3]=[CX3][CX3](=O)[#6]').length || M(m, '[CX3]=[CX3][CX4][CX3](=O)[#6]').length) return null; return undefined; });

/* ---- run 37: ketene / isocyanate + H–Nu (water, alcohol, ammonia, amine, acid), with curly arrows ---- */
const nuHOf = c => { const R = k => ((c.get(k) || {}).smi || 'C');
  if (c.has('NH3')) return ['[H]N([H])[H]', 'N', 'NH₃'];
  if (c.has('RNH2')) return ['[H]N([H])' + R('RNH2'), 'N' + R('RNH2'), 'R–NH₂'];
  if (c.has('R2NH')) return ['[H]N(' + R('R2NH') + ')' + R('R2NH'), 'N(' + R('R2NH') + ')' + R('R2NH'), 'R₂NH'];
  if (c.has('aniline')) return ['[H]N([H])c1ccccc1', 'Nc1ccccc1', 'aniline'];
  if (c.has('AcOH') && !c.any(['H2O', 'H3O+'])) return ['[H]OC(C)=O', 'OC(C)=O', 'CH₃COOH'];
  if (c.has('ROH') && !c.any(['H2O', 'H3O+'])) return ['[H]O' + R('ROH'), 'O' + R('ROH'), 'R–OH'];
  if (c.has('phenol') && !c.any(['H2O', 'H3O+'])) return ['[H]Oc1ccccc1', 'Oc1ccccc1', 'phenol'];
  if (c.any(['H2O', 'H3O+']) || /h2o|water|moist/i.test(c.raw || '')) return ['[H]O[H]', 'O', 'H₂O'];
  return null; };
front({ id: 'cumulene_nu', name: 'Ketene / isocyanate + H–Nu', fam: 'Acid derivative', canon: 'H₂O, ROH, NH₃ or RNH₂', db: ['Ketene', 'Isocyanate'],
  test: c => !!nuHOf(c),
  run(m, c) { const s = smiOf(m); if (!s) return null; const nu = nuHOf(c); if (!nu) return null;
    const kt = M(m, '[CX3]=[CX2]=[OX1]'), ic = M(m, '[#6][NX2]=[CX2]=[OX1]'); if (!kt.length && !ic.length) return null; const isK = kt.length > 0, wat = nu[1] === 'O';
    const p = isK ? rxn('[C:1]=[C:2]=[O:3]>>[C:1][C:2](=[O:3])' + nu[1], [s]) : wat ? rxn('[#6:1][N:2]=[C:3]=[O:4]>>[*:1][N:2].[C:3](=[O:4])=O', [s]) : rxn('[#6:1][N:2]=[C:3]=[O:4]>>[*:1][N:2][C:3](=[O:4])' + nu[1], [s]); if (!p.length) return null;
    const kind = isK ? (wat ? 'carboxylic acid' : /^N/.test(nu[1]) ? 'amide' : nu[1] === 'OC(C)=O' ? 'acid anhydride' : 'ester') : (wat ? 'carbamic acid → 1° amine + CO₂' : /^N/.test(nu[1]) ? 'substituted urea' : 'carbamate (urethane)');
    const o = finish(m, mainFrag(p[0].join('.')), { title: (isK ? 'Ketene' : 'Isocyanate') + ' + ' + nu[2] + ' → ' + kind, text: isK ? 'The sp carbon of a ketene is very electrophilic: H–Nu adds across the C=C, Nu on the carbonyl carbon and H on the other carbon.' : 'H–Nu adds across the N=C of the isocyanate (Nu on carbon, H on N).' + (wat ? ' The carbamic acid R–NH–COOH loses CO₂ at once.' : ''), extra: isK || !wat ? [] : ['O=C=O'], tags: ['addition'] });
    if (!o || !o.product) return o;
    try { const mech = new Mech(m); const [x, cc] = isK ? [kt[0][0], kt[0][1]] : [ic[0][1], ic[0][2]]; const sp = spOf(nu[0]); if (!sp) return o; const ids = mech.add(sp, cc, { att: 2, dist: 1.7 * BL }); const na = ids.get(2), h = ids.get(1);
      mech.step([[[na], [na, cc]], [[cc, x], [x]]], isK ? 'The lone pair of ' + nu[2] + ' attacks the sp carbon of the ketene; the C=C π pair moves onto the other carbon (an enolate-type carbanion).' : 'The lone pair of ' + nu[2] + ' attacks the isocyanate carbon; the C=N π pair moves onto N.');
      mech.step([[[x], [x, h]], [[h, na], [na]]], 'Proton transfer from the ' + (na && /^N/.test(nu[1]) ? 'N⁺–H' : 'O⁺–H') + ' to the ' + (isK ? 'carbon' : 'nitrogen') + ': the ' + (isK ? kind : wat ? 'carbamic acid' : kind) + '.');
      if (!isK && wat) { const h2 = ids.get(3); mech.step([[[x], [x, h2]], [[h2, na], [na]]], 'The N takes the O–H proton (zwitterion R–NH₂⁺–COO⁻).'); mech.step([[[na], [na, cc]], [[cc, x], [x]]], 'The carboxylate re-forms C=O and the C–N bond breaks: CO₂ leaves, the amine stays.'); }
      return useMech(o, mech); } catch (e) { return o; } } });
R17.nuSN2 = nuSN2; R17.rxPartner = rxPartner; R17.sameMain = sameMain;
R17.tk = { Mech, spOf, atomWhere, baseDef, addBase, protonOff, useMech, arrowize, arrowless, acylSub, sameMain, nuDefSN2 };
root.EngineR17 = R17;
})(typeof window !== 'undefined' ? window : globalThis);
