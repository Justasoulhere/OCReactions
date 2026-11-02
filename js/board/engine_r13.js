/* engine_r13.js: run 30. Gaps found by the real-question benchmark (Gambit bank + SKM Blue Book, build/t_qbench.js)
   and by a reagent-by-reagent sweep of the JEE Advanced syllabus (build/t_matrix.js).
   Mostly RDKit reaction templates (like r12). Recipes that correct a wrong older answer are put at the FRONT of the
   recipe list (registered last = tried first); gap-fillers are appended at the END. */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore; if (!C || !E) return;
const { BL, el, hTot } = E; const RP = () => root.RxParse;
const fresh = smi => { try { const m = C.fromSmilesMB(smi, BL); return m && m.atoms.size ? m : null; } catch (e) { return null; } };
const step = (L, R, text, cond) => ({ L: L.clone(), arrows: [], R: R.clone(), text, cond: cond || [] });
const out = o => Object.assign({ steps: [], minor: [], extra: [], why: [], db: [] }, o);
const M = (m, sm) => { try { return C.matchAll(m, sm, {}); } catch (e) { return []; } };
const AR = m => E.aromSet(m);
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
/* apply each template repeatedly on the main (most carbon) fragment; returns { cur, extra, labels } */
function applyAll(smi, T, lim) { let cur = smi; const extra = [], labels = [];
  for (const [sma, label] of T) for (let g = 0; g < (lim || 6); g++) { const p = rxn(sma, [cur]); if (!p.length) break; const frs = byC(p[0].join('.').split('.')); if (frs[0] === cur) break; cur = frs[0]; extra.push(...frs.slice(1)); if (!labels.includes(label)) labels.push(label); }
  return { cur, extra, labels }; }
const ctxWith = (c, add2, flags) => { const st = Object.assign({}, c.st, { sp: new Map(c.st.sp), flags: new Set([...c.st.flags, ...(flags || [])]) }); for (const [k, v] of Object.entries(add2 || {})) { if (v === null) st.sp.delete(k); else st.sp.set(k, Object.assign({ sp: k }, v)); } return RP().stageCtx(st); };
const recipe = id => E.REC.find(x => x.id === id);
/* wrap: fn runs first; when it returns null/undefined the old recipe runs */
function pre(id, fn) { const r = recipe(id); if (!r) return; const run0 = r.run; r.run = function (m, c, all, k) { let o = null; try { o = fn(m, c, all, k, run0.bind(this)); } catch (e) { o = null; } return o == null ? run0.call(this, m, c, all, k) : o; }; }
/* wrap: old recipe first; fn only when the old one has no product */
function post(id, fn) { const r = recipe(id); if (!r) return; const run0 = r.run; r.run = function (m, c, all, k) { const o = run0.call(this, m, c, all, k); if (o && o.product && !o.none) return o; let o2 = null; try { o2 = fn(m, c, all, k, o); } catch (e) { o2 = null; } return o2 || o; }; }
const isSugar = m => M(m, '[CX4][OX2H1]').length >= 3 && (M(m, '[CH1]=O').length + M(m, '[CX3](=O)([CX4][OH1])[CX4][OH1]').length) > 0 && !M(m, 'c').length;

/* ---------- carbocations: one 1,2-shift when it gives a more stable ion ---------- */
function shifted(t, c0) { // t is edited; returns the new centre (a 1,2-H shift or 1,2-methyl shift)
  const ar = AR(t); const base = E.cationScore(t, c0, ar); let best = null;
  for (const nb of t.nbrs(c0)) { if (el(t, nb) !== 'C' || ar.has(nb)) continue;
    if (hTot(t, nb) > 0) { const sc = E.cationScore(t, nb, ar); if (sc > base && (!best || sc > best.sc)) best = { sc, nb, me: null }; }
    else { const mes = t.nbrs(nb).filter(x => x !== c0 && el(t, x) === 'C' && !ar.has(x) && t.nbrs(x).filter(y => el(t, y) !== 'H').length === 1);
      if (mes.length) { const u = t.clone(); u.removeBond(u.bondBetween(nb, mes[0]).id); u.addBond(c0, mes[0], 1); const sc = E.cationScore(u, nb, AR(u)); if (sc > base && (!best || sc > best.sc)) best = { sc, nb, me: mes[0] }; } } }
  if (!best) return c0;
  if (best.me != null) { t.removeBond(t.bondBetween(best.nb, best.me).id); t.addBond(c0, best.me, 1); }
  [c0, best.nb].forEach(id => { const a = t.atoms.get(id); if (a) a.hFix = null; });
  return best.nb;
}
/* rooted SMILES of the alkyl group whose attachment point is atom `att` (first atom in the string) */
function rooted(t, att) { const ids = E.fragIds(t, att); const order = [att, ...ids.filter(id => id !== att)]; try { return C.toSmiles(t, { ids: order }).smi.replace(/\[(C|c|O|S|Cl|Br|I|F)H?\d*\]/g, '$1'); } catch (e) { return null; } }
/* alkyl group of R–X / R–OH after rearrangement of its carbocation; smi = rooted alkyl (attach first) */
function rearrAlkyl(rsmi) { const t = fresh(rsmi); if (!t) return { smi: rsmi, moved: false }; const c0 = [...t.atoms.keys()][0]; const c1 = shifted(t, c0); return { smi: c1 === c0 ? rsmi : (rooted(t, c1) || rsmi), moved: c1 !== c0 }; }

/* =====================================================================================================
   END-of-list gap fillers
   ===================================================================================================== */
/* Wolff–Kishner second stage: the hydrazone + KOH, heat (glycol) */
add({ id: 'wk_stage2', name: 'Wolff–Kishner (hydrazone + base, heat)', fam: 'Carbonyl', canon: 'KOH, ethylene glycol, Δ', db: ['Wolff–Kishner reduction'],
  test: c => c.any(['KOH', 'aqKOH', 'alcKOH', 'KOtBu', 'NaOEt']) && (c.heat || c.flag('glycol')),
  run(m) { const s = smiOf(m); if (!s) return null; const r = applyAll(s, [['[CX3:1]=[NX2][NX3H2]>>[C:1]', 'hydrazone → CH₂']]); if (!r.labels.length) return null;
    return finish(m, r.cur, { title: 'Wolff–Kishner: hydrazone + KOH, Δ → CH₂ (N₂ lost)', text: 'Base removes N–H, the carbanion takes H from the solvent, the second N–H goes, N₂ leaves and the carbanion is protonated: C=N–NH₂ becomes CH₂.', cond: ['KOH', 'glycol', 'Δ'], extra: ['N#N'], why: ['Basic conditions: use it when the molecule is acid-sensitive (Clemmensen is the acidic partner)'], tags: ['reduction'] }); } });
/* acid work-up of nitriles / isocyanides / imines / enamines (also used by solve() for "then H3O+") */
add({ id: 'isocyanide_hyd', name: 'Isocyanide hydrolysis', fam: 'Amine', canon: 'H₃O⁺', db: [],
  test: c => c.has('H3O+') || ((c.has('H2O') || c.flag('aq')) && c.heat),
  run(m) { const s = smiOf(m); if (!s || !/\[C-\]/.test(s)) return null; const p = rxn('[#6:1][N+]#[C-]>>[#6:1][NH2]', [s]); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Isocyanide + H₃O⁺ → 1° amine + HCOOH', text: 'Dilute acid adds water across the isocyanide carbon: R–N≡C → R–NH–CHO (formamide) → R–NH₂ + HCOOH. Nitriles need much harsher hydrolysis; isocyanides are NOT hydrolysed by base.', cond: ['H₃O⁺'], extra: ['O=CO'], why: ['Isocyanide → amine one carbon shorter than the isocyanide (the carbon leaves as formic acid)'], tags: ['hydrolysis'] }); } });
add({ id: 'imine_hyd', name: 'Imine hydrolysis', fam: 'Carbonyl', canon: 'H₃O⁺', db: [],
  test: c => c.has('H3O+') || ((c.has('H2O') || c.flag('aq')) && c.heat),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('[CX3:1]=[NX2;!$(N[OX2]);!$(N[NX3]);!$(N=N);!$(N#*):2]>>[C:1]=O.[N:2]', [s]); if (!p.length) return null; const frs = byC(p[0]);
    const main = M(m, '[CX3]=[NX2][#6]').length ? frs.find(x => /=O/.test(x)) || frs[0] : frs[0];
    return finish(m, main, { title: 'Imine + H₃O⁺ → carbonyl compound (+ amine)', text: 'Water adds to the protonated C=N (carbinolamine) and the amine leaves: the C=O is back.', cond: ['H₃O⁺'], extra: frs.filter(x => x !== main && nC(x) > 0), tags: ['hydrolysis'] }); } });
add({ id: 'enamine_hyd', name: 'Enamine hydrolysis', fam: 'Carbonyl', canon: 'H₃O⁺', db: [],
  test: c => c.has('H3O+'),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('([CX3:1]=[CX3:2]-[NX3H0;!$(N-C=O):3])>>([C:1]-[C:2]=O.[N:3])', [s]); if (!p.length) return null; const frs = byC(p[0].join('.').split('.')); const main = frs.find(x => /=O/.test(x)) || frs[0];
    return finish(m, main, { title: 'Enamine + H₃O⁺ → ketone (+ 2° amine)', text: 'The enamine C=C is protonated at carbon (iminium ion), water adds and the amine leaves: the ketone is regenerated (Stork enamine work-up).', cond: ['H₃O⁺'], extra: frs.filter(x => x !== main), tags: ['hydrolysis'] }); } });
/* C=N reductions (imine / oxime → amine) */
add({ id: 'cn_reduce', name: 'Imine / oxime reduction', fam: 'Amine', canon: 'H₂, Ni', db: [],
  test: c => c.any(['H2', 'NaBH4', 'LiAlH4']) || (c.has('Na') && (c.has('ROH') || c.flag('alc'))),
  run(m) { const s = smiOf(m); if (!s) return null; const r = applyAll(s, [['[CX3:1]=[NX2:2][OX2H1]>>[C:1][NH2:2]', 'oxime → 1° amine'], ['[CX3:1]=[NX2;!$(N[OX2]);!$(N[NX3]);!$(N=N):2]>>[C:1][N:2]', 'imine → amine']]); if (!r.labels.length) return null;
    return finish(m, r.cur, { title: `C=N reduced: ${r.labels.join(', ')}`, text: 'H adds across C=N exactly as across C=O: the carbon gets H and the nitrogen gets H.', extra: r.extra, why: ['Ketone + NH₃, then H₂/Ni = reductive amination'], tags: ['reduction'] }); } });
/* Strecker synthesis: aldehyde / ketone + NH4Cl + KCN → α-amino nitrile (→ α-amino acid with H3O+) */
add({ id: 'strecker', name: 'Strecker synthesis', fam: 'Amine', canon: 'NH₄Cl, KCN', db: ['Strecker synthesis'],
  test: c => c.has('NaCN') && c.any(['NH4Cl', 'NH3']),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('[CX3;$([CH1]([#6])=O),$(C([#6])([#6])=O),$([CH2]=O):1]=[OX1]>>[C:1]([NH2])C#N', [s]); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Strecker: C=O + NH₃ + HCN → α-amino nitrile', text: 'NH₃ (from NH₄Cl) forms the imine; CN⁻ adds to it. Acid hydrolysis of the –C≡N then gives the α-amino acid.', cond: ['NH₄Cl', 'KCN'], why: ['H₃O⁺, Δ next: α-amino acid'], tags: ['AdN'] }); } });
/* Gabriel: imide N–H + R–X (base) → N-alkyl imide */
add({ id: 'gabriel_alk', name: 'Gabriel: N-alkylation of phthalimide', fam: 'Amine', canon: 'KOH then R–X', db: ['Gabriel phthalimide synthesis'],
  test: c => c.has('RX') || c.has('partner'),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[NX3H1,NX3-,$([NX3][K,Na])]([CX3]=O)[CX3]=O').length) return null; const rx = c.get('RX'); const part = rx ? (rx.X || 'Cl') + rx.smi : (c.get('partner') || {}).smi; if (!part) return null;
    const p = rxn('[NX3:1]([CX3:3]=[O:4])[CX3:5]=[O:6].[CX4;!$(C(C)(C)C):2][Cl,Br,I]>>[N+0:1]([C:2])([C:3]=[O:4])[C:5]=[O:6]', [s.replace(/\.?\[K\+\]|\.?\[Na\+\]/g, '').replace('[N-]', '[NH]'), part]); if (!p.length) return null;
    return finish(m, byC(p[0])[0], { title: 'Gabriel: phthalimide anion + R–X → N-alkylphthalimide', text: 'The imide N–H is acidic (two C=O); its anion does SN2 on the alkyl halide. Hydrolysis (or NH₂NH₂) then frees the 1° amine.', why: ['Only 1° amines, and only from alkyl (not aryl) halides'], tags: ['SN2'] }); } });
/* Blanc chloromethylation */
add({ id: 'blanc', name: 'Chloromethylation', fam: 'Aromatic', canon: 'HCHO, HCl, ZnCl₂', db: [],
  test: c => c.has('HCHO') && c.any(['HCl', 'ZnCl2', 'Lucas']),
  run(m) { const s = smiOf(m); if (!s || !M(m, 'c').length) return null; const ps = rxn('[cH1:1]>>[c:1]CCl', [s]).map(x => x[0]); if (!ps.length) return null;
    const para = ps.find(x => M(fresh(x), '[!#1;!c]c1ccc(CCl)cc1').length) || ps[0];
    return finish(m, para, { title: 'Chloromethylation: Ar–H + HCHO + HCl → Ar–CH₂Cl', text: 'Protonated formaldehyde (with ZnCl₂) is the electrophile; the benzyl alcohol formed is converted to the chloride by HCl.', cond: ['HCHO', 'HCl', 'ZnCl₂'], tags: ['SEAr'] }); } });
/* aldehyde + NH3 (cold) → aldehyde ammonia */
add({ id: 'ald_ammonia', name: 'Aldehyde ammonia', fam: 'Carbonyl', canon: 'NH₃ (cold)', db: [],
  test: c => c.has('NH3') && (c.cold || c.flag('rt')),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('[CH1:1](=[O:2])[#6:3]>>[C:1]([OH1:2])([NH2])[#6:3]', [s]); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Aldehyde + NH₃ (cold) → aldehyde ammonia', text: 'NH₃ adds to C=O; in the cold the carbinolamine R–CH(OH)–NH₂ can be isolated (it loses water to the imine on warming).', why: ['HCHO is the exception: it gives urotropine'], tags: ['AdN'] }); } });
/* 1-halides rearrange with AlX3 alone (no arene) */
add({ id: 'rx_isom', name: 'Alkyl halide isomerisation (AlX₃)', fam: 'Halide', canon: 'AlBr₃', db: [],
  test: c => c.has('FeX3') && !c.any(['benzene', 'RX', 'AcCl', 'RCOCl', 'Br2', 'Cl2', 'I2', 'CO', 'partner']),
  run(m) { const s = smiOf(m); if (!s || M(m, 'c').length) return null; const q = M(m, '[CH2X4:1][Cl,Br,I:2]'); if (!q.length) return null;
    const t = m.clone(); const [cc, x] = q[0]; const X = el(t, x); t.removeAtom(x); const nc = shifted(t, cc); if (nc === cc) return null; const R = rooted(t, nc); if (!R) return null;
    return finish(m, X + R, { title: 'R–X + AlX₃ → rearranged (more substituted) halide', text: 'AlX₃ pulls the halide off; the 1° carbocation rearranges by a 1,2-shift and X⁻ comes back on the more stable carbon.', tags: ['rearrangement'] }); } });
/* SN1 / E1 solvolysis in ROH / H2O (no strong base): rearranged cation; heat → alkene */
add({ id: 'solvolysis', name: 'Solvolysis (SN1 / E1)', fam: 'Halide', canon: 'EtOH, Δ', db: [],
  test: c => (c.has('ROH') || c.has('H2O') || c.flag('alc') || c.flag('aq')) && !c.any(['KOH', 'aqKOH', 'alcKOH', 'NaOEt', 'NaOMe', 'NaOR', 'NaCN', 'NaI', 'AgNO3', 'Na', 'NaH', 'KOtBu', 'NH3', 'NaNH2', 'H3O+', 'H2SO4', 'HCl', 'HBr', 'HI']),
  run(m, c) { const q = M(m, '[CX4:1][Cl,Br,I:2]'); if (!q.length || M(m, '[CX4]([Cl,Br,I])[Cl,Br,I]').length) return null; const [cc, x] = q[0]; const ar = AR(m);
    const deg = E.cN(m, cc).length; const neo = deg === 1 && E.cN(m, E.cN(m, cc)[0]).length === 4; if (deg < 2 && !neo) return null;
    const t = m.clone(); t.removeAtom(x); const nc = shifted(t, cc); [cc, nc].forEach(id => { const a = t.atoms.get(id); if (a) a.hFix = null; });
    const s0 = smiOf(t); if (!s0) return null; let prod = null;
    if (c.heat) { // E1: Zaitsev alkene from the cation
      const T = t.clone(); const betas = T.nbrs(nc).filter(b => el(T, b) === 'C' && hTot(T, b) > 0 && !ar.has(b)); if (!betas.length) return null;
      const best = betas.map(b => { const u = T.clone(); u.bondBetween(nc, b).o = 2; [nc, b].forEach(id => { u.atoms.get(id).hFix = null; }); return { s: smiOf(u), sc: E.cN(u, nc).length + E.cN(u, b).length }; }).sort((p, q2) => q2.sc - p.sc)[0];
      prod = best && best.s; }
    else { const R = (c.get('ROH') || {}).smi || (c.has('H2O') || c.flag('aq') ? '' : 'CC'); const u = t.clone(); const p0 = u.atoms.get(nc); const o = u.addAtom('O', p0.x + BL, p0.y); u.addBond(nc, o, 1); let s1 = smiOf(u); if (R) { const v = rxn('[CX4:1][OH1:2]>>[C:1][O:2]' + R, [s1]); if (v.length) s1 = v.find(x => x[0]) [0]; } prod = s1; }
    if (!prod) return null;
    return finish(m, prod, { title: c.heat ? `Solvolysis, Δ: E1 → alkene${nc !== cc ? ' (after a 1,2-shift)' : ''}` : `Solvolysis: SN1 → ${c.has('H2O') || c.flag('aq') ? 'alcohol' : 'ether'}${nc !== cc ? ' (after a 1,2-shift)' : ''}`, text: `With no strong base or nucleophile, the C–X bond ionises (slowly) to a carbocation${nc !== cc ? ', which rearranges to the more stable one by a 1,2-shift' : ''}; the solvent then ${c.heat ? 'removes a β-H (E1, Saytzeff alkene)' : 'attacks it (SN1)'}.`, cond: c.heat ? ['Δ'] : [], why: ['Neopentyl halides cannot do SN2 (too hindered); their solvolysis always goes through the rearranged 3° cation'], tags: ['SN1', 'E1'] }); } });
/* alkane nitration (vapour phase) */
add({ id: 'alkane_nitration', name: 'Vapour-phase nitration of alkanes', fam: 'Alkane', canon: 'conc. HNO₃, 673 K', db: [],
  test: c => c.has('HNO3') && c.heat && !c.has('H2SO4'),
  run(m) { const s = smiOf(m); if (!s || /[^CH()0-9]/.test(s.replace(/C/g, '')) || !/^[C()0-9]+$/.test(s) || s.length < 2) return null;
    const all = new Set(); for (const sm of ['[CX4H1:1]>>[C:1][N+](=O)[O-]', '[CX4H2:1]>>[C:1][N+](=O)[O-]', '[CX4H3:1]>>[C:1][N+](=O)[O-]']) for (const p of rxn(sm, [s])) all.add(p[0]);
    const list = [...all]; if (!list.length) return null; const rank = x => (M(fresh(x), '[CX4H0]([N+])').length ? 3 : M(fresh(x), '[CX4H1]([N+])').length ? 2 : 1);
    list.sort((a, b) => rank(b) - rank(a)); return finish(m, list[0], { title: 'Alkane + HNO₃ (vapour, 673 K) → nitroalkanes (radical)', text: 'At high temperature HNO₃ gives •NO₂; H is abstracted (3° > 2° > 1°) and the alkyl radical takes NO₂. C–C bonds also break, so smaller nitroalkanes form too.', cond: ['HNO₃ (vapour)', '673 K'], minor: list.slice(1, 4), tags: ['radical'] }); } });
/* gem-dihalide / trihalide + water, heat */
add({ id: 'gem_hyd_acid', name: 'gem-Dihalide hydrolysis', fam: 'Halide', canon: 'H₂O, Δ', db: [],
  test: c => (c.has('H2O') || c.has('H3O+')) && c.heat,
  run(m) { const s = smiOf(m); if (!s) return null; const r = applyAll(s, [['[CX4:1]([Cl,Br,I])([Cl,Br,I])[Cl,Br,I]>>[C:1](=O)O', 'trihalide → acid'], ['[CX4:1]([Cl,Br,I])[Cl,Br,I]>>[C:1]=O', 'gem-dihalide → carbonyl']]); if (!r.labels.length) return null;
    return finish(m, r.cur, { title: `Hydrolysis: ${r.labels.join(', ')}`, text: 'Each C–X is replaced by OH; a carbon with two OH loses water (gem-diol → C=O).', tags: ['hydrolysis'] }); } });
/* desulphonation */
add({ id: 'desulph', name: 'Desulphonation', fam: 'Aromatic', canon: 'dil. H₂SO₄ / steam, Δ', db: ['Desulphonation'],
  test: c => (c.has('H2O') || c.has('H3O+') || (c.has('H2SO4') && c.flag('dil'))) && c.heat,
  run(m) { const s = smiOf(m); if (!s) return null; const r = applyAll(s, [['[c:1][SX4](=O)(=O)[OX2H1,OX1-]>>[cH1:1]', '–SO₃H removed']]); if (!r.labels.length) return null;
    return finish(m, r.cur, { title: 'Ar–SO₃H + H₂O (steam / dil. acid, Δ) → Ar–H', text: 'Sulphonation is reversible: H⁺ attacks the carbon carrying –SO₃H and SO₃ leaves.', extra: ['OS(=O)(=O)O'], why: ['Used to block a position temporarily'], tags: ['SEAr'] }); } });
/* reverse Friedel–Crafts: the drawn compound is the acyl / alkyl halide and benzene is the reagent */
add({ id: 'fc_reverse', name: 'Friedel–Crafts with benzene as the reagent', fam: 'Aromatic', canon: 'C₆H₆, AlCl₃', db: ['Friedel–Crafts acylation', 'Friedel–Crafts alkylation'],
  test: c => c.has('benzene') && c.any(['FeX3', 'BF3', 'H2SO4', 'H3PO4', 'HF']),
  run(m) { const s = smiOf(m); if (!s) return null;
    const a = rxn('[CX3:1](=[O:2])[Cl,Br]>>[C:1](=[O:2])c1ccccc1', [s]); if (a.length) return finish(m, a[0][0], { title: 'Friedel–Crafts acylation: R–COCl + C₆H₆ (AlCl₃) → aryl ketone', text: 'AlCl₃ makes the acylium ion R–C≡O⁺, which substitutes the benzene ring. No rearrangement, no polysubstitution.', extra: ['Cl'], tags: ['SEAr'] });
    const q = M(m, '[CX4:1][Cl,Br,I:2]'); if (q.length && !M(m, 'c').length) { const t = m.clone(); t.removeAtom(q[0][1]); const nc = shifted(t, q[0][0]); const R = rooted(t, nc); if (R) return finish(m, 'c1ccccc1' + R, { title: `Friedel–Crafts alkylation: R–X + C₆H₆ (AlCl₃)${nc !== q[0][0] ? ', the carbocation rearranges first' : ''}`, text: 'AlCl₃ ionises the C–X bond; the carbocation (after any 1,2-shift to a more stable one) substitutes the ring.', extra: ['Cl'], tags: ['SEAr'] }); }
    const o = M(m, '[CX4:1][OX2H1:2]'); if (o.length && !M(m, 'c').length) { const t = m.clone(); t.removeAtom(o[0][1]); const nc = shifted(t, o[0][0]); const R = rooted(t, nc); if (R) return finish(m, 'c1ccccc1' + R, { title: 'Alcohol + benzene (acid) → alkylbenzene', text: 'The protonated alcohol loses water; the carbocation (rearranged if possible) attacks the ring.', tags: ['SEAr'] }); }
    return null; } });
/* alkene + alcohol + acid → ether (two drawn molecules) */
add({ id: 'alkene_roh_ether', name: 'Alkene + alcohol (acid) → ether', fam: 'Alkene', canon: 'H₂SO₄', db: [], multi: true,
  test: c => c.any(['H2SO4', 'H3O+', 'H3PO4', 'BF3']),
  run(m) { const s = smiOf(m); if (!s) return null; const fr = s.split('.'); if (fr.length !== 2) return null; const al = fr.find(x => M(fresh(x), '[CX3]=[CX3]').length && !M(fresh(x), '[OH1]').length && !M(fresh(x), 'c').length), oh = fr.find(x => x !== al && M(fresh(x), '[CX4][OX2H1]').length && !M(fresh(x), 'c').length); if (!al || !oh) return null;
    const p = rxn('[CH2:1]=[C:2].[OX2H1:3][CX4:4]>>[CH3:1][C:2][O:3][C:4]', [al, oh]).concat(rxn('[C:1]=[C:2].[OX2H1:3][CX4:4]>>[C:1][C:2][O:3][C:4]', [al, oh])); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Alkene + alcohol, H⁺ → ether (Markovnikov)', text: 'H⁺ adds to the alkene to give the more stable carbocation; the alcohol oxygen attacks it and loses H⁺.', tags: ['AdE'] }); } });
/* aryl / alkyl FC alkylation from a drawn alkene or alcohol + arene */
add({ id: 'fc_alk_drawn', name: 'Friedel–Crafts alkylation (alkene / alcohol + arene, acid)', fam: 'Aromatic', canon: 'H₂SO₄ / HF / BF₃', db: ['Friedel–Crafts alkylation'], multi: true,
  test: c => c.any(['H2SO4', 'H3PO4', 'HF', 'BF3', 'H3O+', 'FeX3']),
  run(m, c, all, k) { const s = smiOf(m); if (!s) return null; const fr = s.split('.'); if (fr.length !== 2) return null;
    const ar = fr.find(x => /c/.test(x)); const other = fr.find(x => x !== ar); if (!ar || !other || /c/.test(other)) return null; const t = fresh(other); if (!t) return null; let cc = null;
    const al = M(t, '[CX3:1]=[CX3:2]'); const oh = M(t, '[CX4:1][OX2H1:2]'); const rx = M(t, '[CX4:1][Cl,Br,I:2]');
    if (al.length) { const [a, b] = al[0]; const ra = AR(t); cc = E.cationScore(t, a, ra) >= E.cationScore(t, b, ra) ? a : b; t.bondBetween(a, b).o = 1; [a, b].forEach(id => { t.atoms.get(id).hFix = null; }); }
    else if (oh.length) { cc = oh[0][0]; t.removeAtom(oh[0][1]); t.atoms.get(cc).hFix = null; }
    else if (rx.length) { cc = rx[0][0]; t.removeAtom(rx[0][1]); t.atoms.get(cc).hFix = null; }
    if (cc == null) return null; const nc = shifted(t, cc); const R = rooted(t, nc); if (!R) return null;
    const eas = recipe('eas'); const A = fresh(ar); if (!eas || !A) return null; const c2 = ctxWith(c, { FeX3: {}, RX: { smi: R, X: 'Cl', txt: 'carbocation' }, H2SO4: null, 'H3O+': null });
    let o = null; try { o = eas.run(A, c2, all, k); } catch (e) { o = null; } if (!o || !o.product) return null;
    return Object.assign({}, o, { title: `Friedel–Crafts alkylation by the ${al.length ? 'alkene' : oh.length ? 'alcohol' : 'halide'} (carbocation${nc !== cc ? ', rearranged' : ''})`, why: ['Alkenes and alcohols alkylate arenes with protic or Lewis acids: the electrophile is the carbocation (Markovnikov)'].concat(o.why || []) }); } });

/* =====================================================================================================
   FRONT recipes: they fix answers that older recipes got wrong, or reach reagents they missed.
   front() puts a recipe at the head of the list, so the LAST one registered here is tried first.
   ===================================================================================================== */
/* Na / ethanol: Bouveault–Blanc (esters), aldehydes, ketones, oximes, imines, alkyl halides (nitriles stay with Mendius) */
front({ id: 'na_etoh', name: 'Na / C₂H₅OH reductions', fam: 'Reduction', canon: 'Na, C₂H₅OH', db: ['Bouveault–Blanc reduction', 'Mendius reduction'],
  test: c => c.has('Na') && (c.has('ROH') || c.flag('alc')) && !c.any(['NH3', 'Na/NH3', 'RX', 'partner']) && !c.flag('liq'),
  run(m) { const s = smiOf(m); if (!s) return null; const has = sm => M(m, sm).length > 0;
    if (has('[CX2]#[NX1]') && !has('[CX3](=O)[OX2][#6]') && !has('[CX3;!$(C-[O,N,S])]=[OX1]') && !has('[CX3]=[NX2]')) return null;
    const r = applyAll(s, [['[CX3:1](=[O:2])[OX2:3][#6:4]>>[CH2:1][OH1:2].[OH1:3][#6:4]', 'ester → 1° alcohol (Bouveault–Blanc)'], ['[CX3;!$(C-[O,N,S,Cl,Br,I]):1]=[OX1:2]>>[C:1][OH1:2]', 'C=O → CH–OH'],
      ['[CX3:1]=[NX2:2][OX2H1]>>[C:1][NH2:2]', 'oxime → 1° amine'], ['[CX3:1]=[NX2;!$(N[OX2]);!$(N[NX3]);!$(N=N):2]>>[C:1][N:2]', 'imine → amine'], ['[CX2:1]#[NX1:2]>>[CH2:1][NH2:2]', 'nitrile → 1° amine (Mendius)'], ['[CX4:1][Cl,Br,I]>>[C:1]', 'R–X → R–H']]);
    if (!r.labels.length) return null;
    return finish(m, r.cur, { title: `Na / C₂H₅OH: ${r.labels.join(', ')}`, text: 'Sodium gives electrons one at a time and ethanol supplies the protons (dissolving-metal reduction). Esters are cut at the acyl–O bond: R–COOR′ → R–CH₂OH + R′OH (Bouveault–Blanc). C=O, C=N and C≡N are reduced; isolated C=C is not touched.', cond: ['Na', 'C₂H₅OH'], extra: r.extra.filter(x => nC(x) > 0), why: ['Na/EtOH does NOT reduce isolated C=C or C≡C (that needs Na/liq. NH₃ for alkynes)', 'Old alternative to LiAlH₄ for esters'], tags: ['reduction'] }); } });

/* haloform: ethanol is positive; CH3CO–OH / –Cl / –OR / –NH2 are negative; Br2 + NaOH works like NaOBr */
front({ id: 'haloform_fix', name: 'Haloform reaction', fam: 'Carbonyl', canon: 'I₂, NaOH', db: ['Haloform reaction'],
  test: c => c.has('NaOCl') || c.has('Br2/NaOH') || (c.any(['I2', 'Br2', 'Cl2']) && c.any(['KOH', 'aqKOH', 'NaHCO3'])),
  run(m, c, all, k) {
    const X = ((c.get('NaOCl') || {}).X) || (c.has('I2') ? 'I' : (c.has('Br2') || c.has('Br2/NaOH')) ? 'Br' : 'Cl'); const s = smiOf(m); if (!s) return null;
    const mk = M(m, '[CH3][CX3](=O)[#6,#1]').length + M(m, '[CH3][CH1]=O').length, mc = M(m, '[CH3][CX4H1]([OX2H1])[#6,#1]').length + M(m, '[CH3][CH2][OX2H1]').length;
    if (M(m, '[CX3](=O)[NH2]').length && (c.has('Br2/NaOH') || c.has('Br2'))) return null; // amides: Hofmann bromamide
    if (!mk && !mc) { if (M(m, '[CH3][CX3](=O)[OX2H1,OX1-,OX2,Cl,Br,NX3]').length) return none('Negative haloform test', 'CH₃CO– is present, but as an acid, ester, amide or acyl chloride: that C=O is not a ketone/aldehyde carbonyl, so the CH₃ is not trihalogenated. (Acid chlorides and esters are simply hydrolysed.)'); return null; }
    if (s === 'CCO' || s === 'OCC') return finish(m, 'O=C[O-]', { title: `Ethanol + ${X}₂ / NaOH → CH${X}₃ + HCOO⁻ (positive)`, text: `${X}₂/NaOH first oxidises CH₃CH₂OH to CH₃CHO; the CH₃ of CH₃CHO is trihalogenated and OH⁻ cleaves off CH${X}₃, leaving formate.`, cond: [`${X}₂`, 'NaOH'], extra: [X + 'C(' + X + ')' + X], why: ['Ethanol is the only 1° alcohol that gives the haloform test', 'Methanol does not'], tags: ['test'] });
    if (c.has('NaOCl') || c.has('I2')) return null; // older recipe (with arrows)
    let r = applyAll(s, [['[CH3:1][CX4H1:2]([OX2H1:3])[#6,#1:4]>>[CH3:1][C:2](=[O:3])[*:4]', 'oxidised to the methyl ketone']], 1); const t = r.cur;
    const p = rxn('[CH3][CX3:2](=[O:3])>>[OX1-][C:2]=[O:3]', [t]); if (!p.length) return null;
    return finish(m, byC(p[0])[0], { title: `Haloform (${X}₂ / NaOH): carboxylate (one C fewer) + CH${X}₃`, text: `OH⁻ makes the enolate of the CH₃–CO– group three times, each one taking ${X}; the C${X}₃ group then leaves when OH⁻ adds to C=O.`, cond: [`${X}₂`, 'NaOH'], extra: [X + 'C(' + X + ')' + X], tags: ['test'] }); } });

/* Hofmann bromamide must win over radical bromination when 'heat' is written */
front({ id: 'hofmann_first', name: 'Hofmann bromamide degradation', fam: 'Amide', canon: 'Br₂, NaOH', db: ['Hofmann bromamide degradation'],
  test: c => c.has('Br2/NaOH') || (c.has('Br2') && c.any(['KOH', 'aqKOH'])) || (c.has('NaOCl') && !c.has('I2')),
  run(m, c, all, k) { if (!M(m, '[CX3](=O)[NX3H2]').length) return null; const r = recipe('amide_rx') || recipe('nbromoamide'); if (!r) return null; const c2 = ctxWith(c, { 'Br2/NaOH': {} }, []); let o = null; try { o = r.run(m, c2, all, k); } catch (e) { o = null; } return o && o.product ? o : null; } });

/* Friedel–Crafts alkylation: 1° alkyl halides rearrange to the more stable carbocation first */
front({ id: 'fc_rearr', name: 'Friedel–Crafts alkylation (rearranging carbocation)', fam: 'Aromatic', canon: 'R–Cl, AlCl₃', db: ['Friedel–Crafts alkylation'],
  test: c => c.has('RX') && c.has('FeX3'),
  run(m, c, all, k) { if (!M(m, 'c').length) return null; const rx = c.get('RX'); if (!rx || !rx.smi) return null; const ra = rearrAlkyl(rx.smi); if (!ra.moved) return null;
    const eas = recipe('eas'); if (!eas) return null; const c2 = ctxWith(c, { RX: Object.assign({}, rx, { smi: ra.smi }) }); let o = null; try { o = eas.run(m, c2, all, k); } catch (e) { o = null; } if (!o || !o.product) return null;
    return Object.assign({}, o, { title: 'Friedel–Crafts alkylation: the 1° carbocation rearranges first', why: ['n-Propyl chloride gives isopropylbenzene; isobutyl chloride gives tert-butylbenzene (1,2-H shift); neopentyl gives tert-amyl (1,2-CH₃ shift)', 'This is why FC acylation + reduction is used to put on an n-alkyl chain'].concat(o.why || []) }); } });

/* nitrating mixture does not hydrolyse amides: nitrate the activated ring instead */
front({ id: 'amide_nitr_guard', name: 'Nitration of anilides', fam: 'Aromatic', canon: 'conc. HNO₃, conc. H₂SO₄', db: [],
  test: c => c.has('HNO3') && (c.has('H2SO4') || c.has('NO2+')) && !c.has('H2O'),
  run(m, c, all, k) { if (!M(m, '[CX3](=O)[NX3H1]c').length) return null; const eas = recipe('eas'); if (!eas) return null; let o = null; try { o = eas.run(m, c, all, k); } catch (e) { o = null; } return o && o.product ? o : null; } });

/* sugars + HI (prolonged heating) → n-alkane chain */
front({ id: 'sugar_hi', name: 'Glucose + HI', fam: 'Biomolecule', canon: 'HI, Δ', db: [],
  test: c => c.has('HI'),
  run(m) { if (!isSugar(m)) return null; const n = [...m.atoms.values()].filter(a => a.el === 'C').length; if (n < 4 || n > 7) return null;
    return finish(m, 'C'.repeat(n), { title: `${n === 6 ? 'Glucose / fructose' : 'Sugar'} + HI (prolonged heating) → n-${['', '', '', '', 'butane', 'pentane', 'hexane', 'heptane'][n]}`, text: 'HI removes every oxygen (reduces all –OH and the C=O): an unbranched alkane forms. This showed glucose has a straight chain of six carbons.', cond: ['HI', 'Δ'], tags: ['reduction'] }); } });

/* Tollens / Fehling: α-hydroxy ketones (fructose, benzoin), formic acid / formates, aryl hydroxylamines are positive */
front({ id: 'tollens_plus', name: "Tollens' / Fehling's (reducing compounds)", fam: 'Carbonyl', canon: "Tollens' reagent", db: ["Tollens' test", "Fehling's / Benedict's test"],
  test: c => c.any(['Tollens', 'Fehling']),
  run(m, c) { const s = smiOf(m); if (!s) return null; const fe = c.has('Fehling'); const name = fe ? "Fehling's" : "Tollens'"; const obs = fe ? 'red Cu₂O precipitate' : 'silver mirror';
    if (s === 'O=CO' || s === 'OC=O') return finish(m, 'O=C([O-])O', { title: `HCOOH + ${name}: positive (${obs})`, text: 'Formic acid has an aldehyde H–C=O: it is oxidised to carbonic acid / CO₂.', tags: ['test'] });
    if (M(m, '[CH1](=O)[OX2][#6]').length) return positive(`Formate ester + ${name}: positive (${obs})`, 'H–COOR still has the H–C=O unit, so it reduces the reagent.');
    if (M(m, 'c[NX3H1][OX2H1]').length && !fe) return positive(`Aryl hydroxylamine + ${name}: positive (${obs})`, 'Ar–NHOH is a reducing agent: it is oxidised (to Ar–N=O) and Ag⁺ is reduced to the mirror.');
    if (M(m, '[CH1]=O').length) return null; // aldehydes: the older recipe
    if (M(m, '[CX3](=O)[CX4;!H0][OX2H1]').length) {
      if (isSugar(m)) return finish(m, 'O=C(O)C(O)C(O)C(O)C(O)CO', { title: `Fructose + ${name}: positive (${obs})`, text: 'In the alkaline reagent the ketose isomerises through the enediol to aldoses (glucose / mannose), which are oxidised. So fructose is a reducing sugar even though it is a ketone.', tags: ['test'] });
      const p = rxn('[CX3:1](=[O:2])[CX4H1:3][OX2H1:4]>>[C:1](=[O:2])[C:3]=[O:4]', [s]); if (p.length) return finish(m, p[0][0], { title: `α-Hydroxy ketone + ${name}: positive (${obs}) → 1,2-diketone`, text: 'The –CH(OH)– next to C=O is oxidised (benzoin → benzil).', tags: ['test'] }); }
    return null; } });

/* Williamson: alcohol / phenol (+ base) + R–X, or an alkoxide / phenoxide + R–X */
front({ id: 'williamson', name: 'Williamson ether synthesis', fam: 'Ether', canon: 'NaH then CH₃I', db: ['Williamson ether synthesis'],
  test: c => c.has('RX'),
  run(m, c) { const s = smiOf(m); if (!s) return null; const rx = c.get('RX'); if (!rx || !rx.smi || rx.gen) return null; const base = c.any(['NaOEt', 'NaOMe', 'KOH', 'aqKOH', 'alcKOH', 'Na', 'NaH', 'NaHCO3', 'NaOR', 'KOtBu']);
    const anion = /\[O-\]|O\[Na\]|O\[K\]/.test(s); if (!base && !anion) return null; if (M(m, '[CX4][Cl,Br,I]').length && !anion) return null;
    const hasOH = M(m, '[OX2H1;$(O[#6]);!$(OC=O);!$(O[CX3]=[CX3])]').length || anion; if (!hasOH) return null;
    const part = (rx.X || 'Cl') + rx.smi; const deg = (() => { const q = fresh(part); if (!q) return 0; const x = [...q.atoms.values()].find(a => /Cl|Br|I/.test(a.el)); return x ? E.cN(q, q.nbrs(x.id)[0]).length : 0; })();
    if (deg >= 3) return null;
    const s1 = s.replace(/\.\[Na\+\]|\.\[K\+\]|\[Na\+\]\.|\[K\+\]\./g, '').replace(/O\[Na\]|O\[K\]/g, '[O-]');
    const p = rxn('[OX2H1,OX1-;$(O[#6]);!$(OC=O):1].[CX4:2][Cl,Br,I]>>[O+0:1][C:2]', [s1, part]); if (!p.length) return null; const pr = byC(p.map(x => x[0]))[0];
    return finish(m, pr, { title: `Williamson: ${M(m, 'c[OX2H1,OX1-]').length ? 'phenoxide' : 'alkoxide'} + R–X → ether`, text: 'The base turns –OH into –O⁻, which does SN2 on the alkyl halide.', why: ['Use a methyl / 1° halide: 3° halides give alkenes (E2) instead', 'Aryl ethers: make the phenoxide; aryl halides do not react'], tags: ['SN2'] }); } });

/* pinacol rearrangement: any acid; the more stable cation forms, then aryl > H > alkyl migrates */
front({ id: 'pinacol2', name: 'Pinacol–pinacolone rearrangement', fam: 'Alcohol', canon: 'H₂SO₄', db: ['Pinacol–pinacolone rearrangement'],
  test: c => c.any(['H2SO4', 'H3O+', 'H3PO4', 'pTsOH']) && !c.has('HNO3') && !c.has('HIO4'),
  run(m) { const d = M(m, '[CX4:1]([OX2H1:2])[CX4:3][OX2H1:4]'); if (!d.length) return null; const ar = AR(m); let best = null;
    for (const [a, oa, b, ob] of d) for (const [ca, oca, cb, ocb] of [[a, oa, b, ob], [b, ob, a, oa]]) { const sc = E.cationScore(m, ca, ar); if (!best || sc > best.sc) best = { ca, oca, cb, ocb, sc }; }
    const { ca, oca, cb, ocb } = best; const subs = m.nbrs(cb).filter(n => n !== ca && n !== ocb && el(m, n) !== 'H');
    const arom = subs.filter(n => ar.has(n)); let g = arom.length ? arom[0] : (hTot(m, cb) > 0 ? 'H' : (subs.sort((p, q) => E.cN(m, q).length - E.cN(m, p).length)[0]));
    if (g == null) return null; const t = m.clone(); t.removeAtom(oca);
    if (g !== 'H') { t.removeBond(t.bondBetween(cb, g).id); t.addBond(ca, g, 1); }
    t.bondBetween(cb, ocb).o = 2; [ca, cb, ocb].concat(g === 'H' ? [] : [g]).forEach(id => { const a = t.atoms.get(id); if (a) a.hFix = null; });
    const ps = smiOf(t); if (!ps) return null; const who = g === 'H' ? 'H' : ar.has(g) ? 'the aryl group' : 'an alkyl group';
    return finish(m, ps, { title: 'Pinacol–pinacolone rearrangement', text: `H⁺ takes the –OH whose loss gives the more stable carbocation; ${who} shifts from the neighbouring carbon (1,2-shift) and that carbon's –OH becomes C=O.`, cond: ['H⁺'], why: ['Order of migration: aryl > H > alkyl (bigger alkyl before CH₃)', 'Ring bonds can migrate too: ring expansion'], tags: ['rearrangement'] }); } });

/* alcoholic AgNO3: precipitate of AgX for halides that ionise (3°, benzylic, allylic > 2° > 1°); none for vinyl / aryl halides */
front({ id: 'alc_agno3', name: 'Alcoholic AgNO₃ test', fam: 'Halide', canon: 'alc. AgNO₃', db: [],
  test: c => c.has('AgNO3') && !c.has('NH3') && !c.has('Tollens') && !c.has('KOH'),
  run(m) { const s = smiOf(m); if (!s) return null;
    if (/\[Cl-\]|\[Br-\]|\[I-\]/.test(s)) return positive('Ionic halide + AgNO₃: precipitate', 'The halide is already an ion: AgX precipitates at once.');
    const q = M(m, '[CX4:1][Cl,Br,I:2]'); if (!q.length) { if (M(m, '[c,$(C=C)][Cl,Br,I]').length) return none('No precipitate', 'Aryl and vinyl halides do not ionise: the C–X bond has partial double-bond character and the cation would be unstable.'); return null; }
    if (M(m, '[CX4]([Cl,Br,I])([Cl,Br,I])[Cl,Br,I]').length || M(m, '[CX4]([Cl,Br,I])[Cl,Br,I]').length) return none('No precipitate', 'CCl₄ / CHCl₃ / CH₂Cl₂: the C–Cl bonds do not ionise, so no AgCl forms.');
    const [cc] = q[0]; const ar = AR(m); const deg = E.cN(m, cc).length; const benz = m.nbrs(cc).some(n => ar.has(n)); const allyl = M(m, '[CX4]([Cl,Br,I])[CX3]=[CX3]').length > 0;
    const how = benz || allyl || deg >= 3 ? 'at once (fast SN1: stable carbocation)' : deg === 2 ? 'slowly' : 'only on warming';
    return positive(`Alcoholic AgNO₃: precipitate of AgX ${how}`, `The C–X bond ionises (SN1) and Ag⁺ takes the halide. Order: 3°, benzylic, allylic > 2° > 1°; vinyl and aryl halides give no precipitate.`); } });

/* sugars with bromine water (aldoses only) and nitric acid */
front({ id: 'sugar_ox', name: 'Glucose: Br₂ water / HNO₃ oxidation', fam: 'Biomolecule', canon: 'Br₂ water', db: [],
  test: c => c.has('HOBr') || (c.has('Br2') && (c.has('H2O') || c.flag('aq'))) || (c.has('HNO3') && !c.has('H2SO4')),
  run(m, c) { if (!isSugar(m)) return null; const s = smiOf(m); const ald = M(m, '[CH1](=O)[CX4][OX2H1]').length > 0; const hno3 = c.has('HNO3') && !c.has('HOBr') && !c.has('Br2');
    if (!ald) return hno3 ? null : none('No reaction (ketose)', 'Bromine water oxidises only the –CHO of aldoses. Fructose (a ketose) does not decolourise it: this tells glucose from fructose.');
    if (hno3) { const r = applyAll(s, [['[CH1:1](=O)[CX4;$(C[OX2H1]):2]>>[C:1](=O)(O)[C:2]', 'CHO → COOH'], ['[CH2:1]([OX2H1:2])[CX4;$(C[OX2H1]):3]>>[C:1](=O)([OX2H1:2])[C:3]', 'CH₂OH → COOH']]); return finish(m, r.cur, { title: 'Glucose + HNO₃ → saccharic acid (glucaric acid)', text: 'Nitric acid oxidises both ends: the –CHO and the –CH₂OH become –COOH. So glucose has a primary –OH at the far end.', tags: ['oxidation'] }); }
    const p = rxn('[CH1:1](=O)[CX4;$(C[OX2H1]):2]>>[C:1](=O)(O)[C:2]', [s]); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Aldose + Br₂ water → aldonic acid (gluconic acid)', text: 'Bromine water is a mild oxidant: only –CHO becomes –COOH. The aldehyde group of glucose is shown this way.', tags: ['oxidation'] }); } });

/* sugars + ROH / dry HCl → glycoside (cyclic acetal), not an open-chain acetal */
front({ id: 'glycoside', name: 'Glycoside formation', fam: 'Biomolecule', canon: 'CH₃OH, dry HCl', db: [],
  test: c => c.any(['HCl', 'H3O+', 'H2SO4']) && (c.has('ROH') || c.flag('alc')),
  run(m, c) { if (!isSugar(m)) return null; const s = smiOf(m); const R = ((c.get('ROH') || {}).smi) || 'C';
    let p = rxn(`[CH1:1](=[O:2])[CX4:3][CX4:4][CX4:5][CX4:6]([OX2H1:7])>>[C:1]1(O${R})[C:3][C:4][C:5][C:6][O:7]1`, [s]);
    if (!p.length) p = rxn(`[CX3:1](=[O:2])([CX4:9])[CX4:3][CX4:4][CX4:5]([OX2H1:6])>>[C:1]1([C:9])(O${R})[C:3][C:4][C:5][O:6]1`, [s]);
    if (!p.length) return null;
    return finish(m, p[0][0], { title: `Sugar + ${R === 'C' ? 'CH₃OH' : 'ROH'}, dry HCl → ${R === 'C' ? 'methyl' : 'alkyl'} glycoside`, text: 'Only the anomeric (hemiacetal) –OH is replaced: the ring form reacts, giving the cyclic acetal (α and β anomers). One equivalent of alcohol, not two: more evidence for the ring structure.', tags: ['acetal'] }); } });

/* excess HIO4 (and every sugar): cut every C–C between carbons carrying OH / =O */
front({ id: 'hio4_all', name: 'Periodic acid, exhaustive cleavage', fam: 'Alcohol', canon: 'HIO₄ (excess)', db: [],
  test: c => c.has('HIO4'),
  run(m, c) { const s = smiOf(m); if (!s) return null; if (!(c.eq === Infinity || c.flag('excess') || isSugar(m) || (c.eq || 0) >= 2)) return null;
    const T = ['[CX4:1]([OX2H1:2])-[CX4:3]([OX2H1:4])>>[C:1]=[O:2].[C:3]=[O:4]', '[CX3:1](=[O:2])-[CX4:3]([OX2H1:4])>>[C:1](=[O:2])[OX2H1].[C:3]=[O:4]', '[CX3:1](=[O:2])-[CX3:3](=[O:4])>>[C:1](=[O:2])[OX2H1].[C:3](=[O:4])[OX2H1]'];
    let frs = [s], guard = 0, cut = 0;
    while (guard++ < 12) { let changed = false; const nx = []; for (const f of frs) { let done = false; for (const sm of T) { const p = rxn(sm, [f]); if (p.length) { nx.push(...p[0].join('.').split('.')); done = true; changed = true; cut++; break; } } if (!done) nx.push(f); } frs = nx; if (!changed) break; }
    if (!cut) return null; const order = byC(frs); const main = order[0];
    return finish(m, main, { title: `HIO₄ (excess): every 1,2-diol / α-hydroxy carbonyl bond is cut (${cut} cuts)`, text: 'Each C–C bond between two carbons that carry –OH or =O is cleaved through a cyclic periodate ester. Glucose takes 5 HIO₄: 5 HCOOH + 1 HCHO.', extra: order.slice(1), tags: ['oxidation'] }); } });

/* alkyl halides + HI / red P → alkane */
post('hi_redp', (m) => { const s = smiOf(m); if (!s) return null; const r = applyAll(s, [['[CX4:1][Cl,Br,I]>>[C:1]', 'R–X → R–H']]); if (!r.labels.length) return null; return finish(m, r.cur, { title: 'HI / red P: R–X → R–H', text: 'HI reduces the C–I bond (I₂ forms; red P turns the I₂ back into HI).', cond: ['HI', 'red P'], tags: ['reduction'] }); });

/* bromine water on phenols / anilines carrying –COOH or –SO3H ortho / para: ipso substitution → 2,4,6-tribromo */
front({ id: 'ipso_br2', name: 'Bromine water: ipso substitution', fam: 'Aromatic', canon: 'Br₂ water', db: [],
  test: c => c.has('HOBr') || (c.has('Br2') && (c.has('H2O') || c.flag('aq'))),
  run(m, c, all, k) { const s = smiOf(m); if (!s) return null; const T = ['[OX2H1,NX3H2:9][c:1]1[c:2][c:3][c:4]([CX3](=O)[OX2H1,OX1-])[c:5][c:6]1>>[*:9][c:1]1[c:2][c:3][cH1:4][c:5][c:6]1', '[OX2H1,NX3H2:9][c:1]1[c:2][c:3][c:4]([SX4](=O)(=O)[OX2H1,OX1-])[c:5][c:6]1>>[*:9][c:1]1[c:2][c:3][cH1:4][c:5][c:6]1', '[OX2H1,NX3H2:9][c:1][c:2][CX3](=O)[OX2H1,OX1-]>>[*:9][c:1][cH1:2]', '[OX2H1,NX3H2:9][c:1][c:2][SX4](=O)(=O)[OX2H1,OX1-]>>[*:9][c:1][cH1:2]'];
    const r = applyAll(s, T.map(t => [t, 'ipso'])); if (!r.labels.length) return null; const eas = recipe('eas'); const A = fresh(r.cur); if (!eas || !A) return null; let o = null; try { o = eas.run(A, c, all, k); } catch (e) { o = null; } if (!o || !o.product) return null;
    return Object.assign({}, o, { steps: [step(m, A, 'The –COOH / –SO₃H ortho or para to –OH is replaced by Br (ipso attack; CO₂ / H₂SO₄ leaves).', ['Br₂ water'])].concat(o.steps), title: 'Bromine water: ipso substitution, then 2,4,6-tribromination', why: ['Salicylic acid and p-hydroxybenzenesulphonic acid both give 2,4,6-tribromophenol'].concat(o.why || []) }); } });

/* Kolbe electrolysis of 1,2-dicarboxylates: succinate → ethene, maleate / fumarate → ethyne */
front({ id: 'kolbe_dicarb', name: 'Kolbe electrolysis of dicarboxylates', fam: 'Acid', canon: 'electrolysis', db: ['Kolbe electrolysis'],
  test: c => c.has('electrolysis'),
  run(m) { const s = smiOf(m); if (!s) return null; let p = rxn('[OX1-,OX2H1][CX3](=O)[CX3:1]=[CX3:2][CX3](=O)[OX1-,OX2H1]>>[C:1]#[C:2]', [s]); let what = 'ethyne (from maleate / fumarate)';
    if (!p.length) { p = rxn('[OX1-,OX2H1][CX3](=O)[CX4:1][CX4:2][CX3](=O)[OX1-,OX2H1]>>[C:1]=[C:2]', [s]); what = 'an alkene (from succinate)'; } if (!p.length) return null;
    return finish(m, byC(p[0])[0], { title: `Kolbe electrolysis of a 1,2-dicarboxylate → ${what}`, text: 'Both carboxylates lose an electron and CO₂ at the anode; the two radical centres on neighbouring carbons pair up as a new π bond.', extra: ['O=C=O'], tags: ['radical'] }); } });

/* malonic acid + P2O5 → carbon suboxide */
front({ id: 'suboxide', name: 'Malonic acid + P₂O₅', fam: 'Acid', canon: 'P₂O₅, Δ', db: [],
  test: c => c.has('P2O5'),
  run(m) { const s = smiOf(m); if (s !== 'O=C(O)CC(=O)O') return null; return finish(m, 'O=C=C=C=O', { title: 'Malonic acid + P₂O₅ → carbon suboxide (C₃O₂)', text: 'P₂O₅ removes two molecules of water from CH₂(COOH)₂: O=C=C=C=O.', extra: ['O'], tags: ['dehydration'] }); } });

/* lactone + NH3 / RNH2: the ring opens to a hydroxy amide */
front({ id: 'lactone_amm', name: 'Lactone ammonolysis', fam: 'Ester', canon: 'NH₃', db: [],
  test: c => c.any(['NH3', 'RNH2']),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[CX3;R](=O)[OX2;R][#6;R]').length) return null; const N = c.has('RNH2') ? 'N' + (c.get('RNH2').smi || 'C') : 'N';
    const p = rxn(`([CX3;R:1](=[O:2])-;@[OX2;R:3])>>([C:1](=[O:2])${N}.[OX2H1:3])`, [s]); if (!p.length) return null;
    return finish(m, byC(p[0].join('.').split('.'))[0], { title: 'Lactone + NH₃ → hydroxy amide (ring opens)', text: 'NH₃ attacks the ring C=O; the ring oxygen leaves but stays on the chain as –OH.', tags: ['acyl substitution'] }); } });

/* Tishchenko: aldehyde + Al(OEt)3 → ester */
front({ id: 'tishchenko', name: 'Tishchenko reaction', fam: 'Carbonyl', canon: 'Al(OC₂H₅)₃', db: ['Tishchenko reaction'],
  test: c => c.has('AlOEt') || (c.has('AlOiPr') && !c.flag('alc') && !c.has('ROH') && !c.has('acetone')),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CH1,CH2]=O').length || M(m, '[CX3](=O)[#6]').length - M(m, '[CH1](=O)[#6]').length > 0) return null;
    let p = rxn('([CH1:1](=[O:2])[CH2:3][CH2:4][CH1:5]=[O:6])>>([C:1]1(=[O:2])[CH2:3][CH2:4][CH2:5][O:6]1)', [s]); if (!p.length) p = rxn('([CH1:1](=[O:2])[CH2:3][CH2:4][CH2:7][CH1:5]=[O:6])>>([C:1]1(=[O:2])[CH2:3][CH2:4][CH2:7][CH2:5][O:6]1)', [s]);
    if (!p.length) p = rxn('[#6,#1:1][CH1:2]=[O:3].[#6,#1:4][CH1:5]=[O:6]>>[*:1][C:2](=[O:3])[O:6][CH2:5][*:4]', [s, s]); if (!p.length && s === 'C=O') p = [['COC=O']];
    if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Tishchenko: 2 RCHO + Al(OC₂H₅)₃ → RCOOCH₂R', text: 'Aluminium alkoxide joins two aldehyde molecules; a hydride moves from one to the other (a Cannizzaro-like disproportionation inside the complex) and the ester forms. CH₃CHO gives ethyl acetate.', why: ['With isopropanol present the same reagent does MPV reduction instead'], tags: ['hydride_transfer'] }); } });

/* Tiffeneau–Demjanov / neopentyl: 1° amine + HNO2 with a rearranging cation */
front({ id: 'demjanov', name: 'Demjanov rearrangement', fam: 'Amine', canon: 'NaNO₂, HCl', db: [],
  test: c => c.has('HNO2'),
  run(m) { const s = smiOf(m); if (!s || M(m, 'c[NX3]').length) return null;
    let p = rxn('([NX3H2:1][CH2:2][C;R:3]-;@[C;R:4])>>([OX2H1:1][C:3]-[C:2]-[C:4])', [s]); let what = 'ring expansion';
    if (!p.length) { p = rxn('[NX3H2:1][CH2:2][CX4H0:3][CH3:4]>>[OX2H1:1][C:3][CH2:2][CH3:4]', [s]); what = '1,2-methyl shift'; } if (!p.length) return null;
    return finish(m, byC(p[0].join('.').split('.'))[0], { title: `1° amine + HNO₂: N₂ leaves and the cation rearranges (${what})`, text: what === 'ring expansion' ? 'The CH₂–N₂⁺ loses N₂; a ring C–C bond moves to the CH₂ (the ring grows by one carbon, strain falls) and water traps the new cation (Tiffeneau–Demjanov).' : 'The neopentyl-type CH₂–N₂⁺ loses N₂; a methyl moves over to give a 3° carbocation, which water traps.', cond: ['NaNO₂', 'HCl'], extra: ['N#N'], tags: ['rearrangement'] }); } });

/* BH3 reduces –COOH (not esters, nitro, ketones with 1 eq) */
front({ id: 'bh3_acid', name: 'Diborane reduction of acids', fam: 'Acid', canon: 'B₂H₆ then H₂O', db: [],
  test: c => c.has('BH3') && !c.any(['H2O2', 'H2O2/OH', 'AcOH', 'NH2Cl']),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX3](=O)[OX2H1]').length || M(m, '[CX3]=[CX3]').length || M(m, 'C#C').length) return null; const p = rxn('[CX3:1](=[O:2])[OX2H1:3]>>[CH2:1][OX2H1:3]', [s]); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'B₂H₆: R–COOH → R–CH₂OH (selective)', text: 'Diborane, an electrophilic reducing agent, attacks the electron-rich carboxyl oxygen first and reduces –COOH quickly; esters, nitro groups and halides are left alone.', why: ['LiAlH₄ would reduce the ester / ketone too'], tags: ['reduction'] }); } });

/* 1° amide + POCl3 / SOCl2 / PCl5 → nitrile */
front({ id: 'amide_dehyd', name: 'Amide dehydration', fam: 'Amide', canon: 'POCl₃ / SOCl₂, Δ', db: [],
  test: c => c.any(['POCl3', 'SOCl2', 'PCl5']),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX3](=O)[NX3H2]').length || M(m, '[OX2H1]').length) return null; const p = rxn('[CX3:1](=O)[NX3H2:2]>>[C:1]#[N:2]', [s]); if (!p.length) return null;
    return finish(m, p[0][0], { title: '1° amide + POCl₃ / SOCl₂ / PCl₅ → nitrile', text: 'The amide oxygen attacks P (or S) and is pulled off; two N–H protons are lost: R–C≡N.', extra: ['O'], tags: ['dehydration'] }); } });

/* halohydrin / ω-halo alcohol + base → cyclic ether (epoxide, THF, THP) */
front({ id: 'halo_alc_ether', name: 'Intramolecular Williamson', fam: 'Ether', canon: 'NaOH', db: [],
  test: c => c.any(['KOH', 'aqKOH', 'NaH', 'Na', 'NaOEt', 'NaOMe', 'KOtBu']) && !c.has('RX'),
  run(m) { const s = smiOf(m); if (!s) return null; const T = [['([OX2H1:1][CX4:2][CX4:3][Cl,Br,I])>>([O:1]1[C:2][C:3]1)', 'epoxide'], ['([OX2H1:1][CX4:2][CX4:3][CX4:4][CX4:5][Cl,Br,I])>>([O:1]1[C:2][C:3][C:4][C:5]1)', 'tetrahydrofuran (5-ring)'], ['([OX2H1:1][CX4:2][CX4:3][CX4:4][CX4:5][CX4:6][Cl,Br,I])>>([O:1]1[C:2][C:3][C:4][C:5][C:6]1)', 'tetrahydropyran (6-ring)']];
    for (const [sm, what] of T) { const p = rxn(sm, [s]); if (p.length) return finish(m, byC(p[0].join('.').split('.'))[0], { title: `Halo-alcohol + base → cyclic ether (${what})`, text: 'The alkoxide made by the base attacks the C–X carbon of the same molecule (intramolecular SN2, backside): the ring closes.', why: ['Halohydrin + NaOH is how epoxides are made'], tags: ['SN2'] }); }
    return null; } });

/* Favorskii: cyclic α-halo ketone + OH⁻ / RO⁻ → ring-contracted acid / ester */
front({ id: 'favorskii', name: 'Favorskii rearrangement', fam: 'Carbonyl', canon: 'NaOH', db: [],
  test: c => c.any(['KOH', 'aqKOH', 'NaOMe', 'NaOEt', 'NaOR']) && !c.any(['Br2', 'I2', 'Cl2', 'NaOCl', 'Br2/NaOH', 'RX']),
  run(m, c) { const s = smiOf(m); if (!s) return null; const OR = c.has('NaOMe') ? 'OC' : c.has('NaOEt') ? 'OCC' : '[OH1]';
    const p = rxn(`([CX3;R:1](=[O:2])(-;@[CX4;R:3][Cl,Br,I])-;@[CX4;R;!H0:4])>>([C:1](=[O:2])(${OR})[C:4]-[C:3])`, [s]); if (!p.length) return null;
    return finish(m, byC(p[0].join('.').split('.'))[0], { title: 'Favorskii: cyclic α-halo ketone + base → ring-contracted acid / ester', text: 'Base removes the α′-H; the enolate closes a cyclopropanone by displacing the halide; OH⁻ (or RO⁻) adds to its C=O and the three-membered ring opens to the more stable carbanion. The ring is one carbon smaller.', why: ['2-Chlorocyclohexanone → cyclopentanecarboxylic acid'], tags: ['rearrangement'] }); } });

/* CHCl3, CCl4, CH2Cl2 and 1,2-dichloroethane do not react with plain water */
front({ id: 'polyhalo_water', name: 'Polyhalides and water', fam: 'Halide', canon: 'H₂O', db: [],
  test: c => c.has('H2O') && !c.any(['KOH', 'aqKOH', 'alcKOH', 'NaOEt', 'Ag2O', 'H3O+', 'AgNO3']) && !c.heat,
  run(m) { const s = smiOf(m); if (!s || M(m, '[!#6;!#1;!Cl;!Br;!I]').length || !M(m, '[CX4][Cl,Br,I]').length) return null; const nx = M(m, '[Cl,Br,I]').length; if (nx < 2) return null;
    return none('No reaction', 'Water alone does not hydrolyse these halides at room temperature (they are insoluble and the C–Cl bonds are strong). With aq. KOH they react: CHCl₃ → HCOO⁻, CH₂Cl₂ → HCHO.'); } });

/* "H⁺, Δ" (no water): dehydration of alcohols */
front({ id: 'acid_heat_dehyd', name: 'Dehydration with H⁺, Δ', fam: 'Alcohol', canon: 'H⁺, Δ', db: [],
  test: c => c.has('H3O+') && c.heat && !c.has('H2O') && !c.flag('aq') && !c.flag('dil'),
  run(m, c, all, k) { if (!M(m, '[CX4;!$(C[OX2H1][CX4][OX2H1])][OX2H1]').length || M(m, '[CX3](=O)[OX2,NX3]').length || M(m, 'C#N').length || M(m, '[CX4]([OX2])[OX2]').length) return null;
    const r = recipe('dehyd'); if (!r) return null; let o = null; try { o = r.run(m, ctxWith(c, { H2SO4: {} }, ['conc', 'heat']), all, k); } catch (e) { o = null; } return o && o.product ? o : null; } });

/* intramolecular Cannizzaro of non-enolisable dialdehydes (o-phthalaldehyde) */
front({ id: 'intra_cannizzaro', name: 'Intramolecular Cannizzaro', fam: 'Carbonyl', canon: 'conc. NaOH', db: [],
  test: c => c.any(['KOH', 'aqKOH']) && !c.flag('dil') && !c.any(['I2', 'Br2', 'NaOCl', 'CHCl3']),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('[CH1:1](=[O:2])[c:3]:[c:4][CH1:5]=[O:6]>>[C:1](=[O:2])([OX1-])[c:3]:[c:4][CH2:5][OX2H1:6]', [s]); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Intramolecular Cannizzaro: dialdehyde → hydroxy carboxylate', text: 'OH⁻ adds to one CHO; that carbon passes a hydride to the other CHO of the same molecule.', why: ['Acidification closes the γ-hydroxy acid to the lactone (phthalide)'], tags: ['hydride_transfer'] }); } });

/* LiAlH4 + nitroarene → azobenzene (aryl nitro groups are not reduced to amines by LiAlH4) */
front({ id: 'lah_arno2', name: 'LiAlH₄ + nitroarene', fam: 'Nitro', canon: 'LiAlH₄', db: [],
  test: c => c.has('LiAlH4'),
  run(m) { const s = smiOf(m); if (!s || !M(m, 'c[N+](=O)[O-]').length || M(m, '[CX3]=O').length || M(m, 'C#N').length) return null; const p = rxn('[c:1][N+](=O)[O-].[c:2][N+](=O)[O-]>>[c:1]N=N[c:2]', [s, s]); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Nitrobenzene + LiAlH₄ → azobenzene', text: 'LiAlH₄ reduces aryl nitro groups only part of the way: two Ar–N=O / Ar–NHOH units condense to Ar–N=N–Ar. (Aliphatic nitro compounds do give amines.)', tags: ['reduction'] }); } });

/* imines / oximes with H2: the older H2 recipe has no rule for C=N */
front({ id: 'cn_reduce_h2', name: 'Imine / oxime + H₂', fam: 'Amine', canon: 'H₂, Ni', db: [],
  test: c => c.has('H2') || c.has('NaBH4') || c.has('LiAlH4'),
  run(m, c, all, k) { if (!M(m, '[CX3]=[NX2;!$(N=N)]').length || M(m, '[CX3]=[CX3]').length || M(m, 'C#C').length || M(m, '[N+](=O)[O-]').length) return null; const r = recipe('cn_reduce'); return r ? r.run(m, c, all, k) : null; } });

/* Grignard formed in the presence of the drawn electrophile (Barbier) */
front({ id: 'barbier', name: 'Grignard reagent with a drawn partner', fam: 'Organometallic', canon: 'Mg, dry ether then H₃O⁺', db: [], multi: true,
  test: c => c.any(['Mg', 'Li']),
  run(m) { const s = smiOf(m); if (!s) return null; const fr = s.split('.'); if (fr.length !== 2) return null; const rx = fr.find(x => M(fresh(x), '[#6][Cl,Br,I]').length && !M(fresh(x), '[CX3]=O').length && !M(fresh(x), 'C#N').length); const el2 = fr.find(x => x !== rx); if (!rx || !el2) return null;
    const T = [['[#6:1][Cl,Br,I].[CX3;!$(C[OX2,NX3,Cl]):2]=[OX1:3]>>[#6:1][C:2][OX2H1:3]', 'adds to C=O'], ['[#6:1][Cl,Br,I].[CX2:2]#[NX1:3]>>[#6:1][C:2]=[NX2H1:3]', 'adds to C≡N (imine salt → ketone on work-up)'], ['[#6:1][Cl,Br,I].[CH2:2]1[OX2:3][C:4]1>>[#6:1][C:2][C:4][OX2H1:3]', 'opens the epoxide at the CH₂'], ['[#6:1][Cl,Br,I].[OX1]=[CX2]=[OX1]>>[#6:1]C(=O)O', 'adds to CO₂']];
    for (const [sm, what] of T) { const p = rxn(sm, [rx, el2]); if (p.length) { const pr = byC(p.map(x => x[0]))[0]; return finish(m, pr, { title: `R–MgX forms and ${what}`, text: 'Mg inserts into the C–X bond (R–MgX); its carbanion-like carbon attacks the electrophilic carbon of the other drawn molecule.', cond: ['Mg', 'dry ether'], tags: ['AdN'] }); } }
    return null; } });

/* terminal alkyne + R–X drawn together, NaNH2 → internal alkyne */
front({ id: 'acetylide_alk', name: 'Acetylide alkylation', fam: 'Alkyne', canon: 'NaNH₂ then R–X', db: [], multi: true,
  test: c => c.any(['NaNH2', 'Na/NH3', 'Na']),
  run(m) { const s = smiOf(m); if (!s) return null; const fr = s.split('.'); if (fr.length !== 2) return null; const yne = fr.find(x => M(fresh(x), '[CX2H1]#[CX2]').length), rx = fr.find(x => x !== yne && M(fresh(x), '[CX4;!$(C(C)(C)C)][Cl,Br,I]').length); if (!yne || !rx) return null;
    const p = rxn('[CX2:1]#[CX2H1:2].[CX4;!$(C(C)C):3][Cl,Br,I]>>[C:1]#[C:2][C:3]', [yne, rx]); if (!p.length) return null;
    return finish(m, byC(p.map(x => x[0]))[0], { title: 'Acetylide + R–X (SN2) → internal alkyne', text: 'NaNH₂ removes the acidic ≡C–H; the acetylide ion does SN2 on the 1° alkyl halide: a new C–C bond.', why: ['2°/3° halides give elimination with the strongly basic acetylide'], tags: ['SN2'] }); } });

/* active-methylene (malonic / acetoacetic ester) C-alkylation, and hydrolysis + decarboxylation */
front({ id: 'malonic_alk', name: 'Malonic / acetoacetic ester alkylation', fam: 'Ester', canon: 'NaOEt then R–X', db: ['Malonic ester synthesis'], multi: true,
  test: c => c.any(['NaOEt', 'NaOMe', 'NaH', 'KOH', 'NaOR']),
  run(m, c) { const s = smiOf(m); if (!s) return null; const fr = s.split('.'); let am = null, rx = null;
    if (fr.length === 2) { am = fr.find(x => M(fresh(x), '[CX3](=O)[CH2,CH1][CX3](=O)').length); rx = fr.find(x => x !== am && M(fresh(x), '[CX4][Cl,Br,I]').length); }
    if (!am || !rx) return null; const p = rxn('[CX3:1](=[O:5])[CH2,CH1:2][CX3:3](=[O:6]).[CX4:4][Cl,Br,I]>>[C:1](=[O:5])[C:2]([C:4])[C:3](=[O:6])', [am, rx]); if (!p.length) return null;
    return finish(m, byC(p.map(x => x[0]))[0], { title: 'Active methylene + base + R–X → C-alkylated ester', text: 'The CH₂ between two C=O groups is acidic; its enolate does SN2 on the alkyl halide.', why: ['Hydrolysis and heating then remove one –COOH (malonic) or give a ketone (acetoacetic)'], tags: ['SN2'] }); } });
add({ id: 'malonic_alk_rx', name: 'Malonic ester + R–X', fam: 'Ester', canon: 'NaOEt, R–X', db: ['Malonic ester synthesis'],
  test: c => c.has('RX') && c.any(['NaOEt', 'NaOMe', 'NaH', 'NaOR']),
  run(m, c) { const s = smiOf(m); const rx = c.get('RX'); if (!s || !rx || !rx.smi || !M(m, '[CX3](=O)[CH2,CH1][CX3](=O)').length) return null; const p = rxn('[CX3:1](=[O:5])[CH2,CH1:2][CX3:3](=[O:6]).[CX4:4][Cl,Br,I]>>[C:1](=[O:5])[C:2]([C:4])[C:3](=[O:6])', [s, (rx.X || 'Br') + rx.smi]); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Active methylene + base + R–X → C-alkylation', text: 'Enolate SN2 on the alkyl halide.', tags: ['SN2'] }); } });
front({ id: 'malonic_hyd', name: 'Malonic / acetoacetic ester: hydrolysis + decarboxylation', fam: 'Ester', canon: 'H₃O⁺, Δ', db: ['Malonic ester synthesis'],
  test: c => c.heat && c.any(['H3O+', 'H2SO4', 'HCl', 'KOH', 'aqKOH']),
  run(m) { const s = smiOf(m); if (!s) return null; let p = rxn('[CH3][CX3:1](=[O:5])[CX4;!H0:2][CX3](=O)[OX2][#6]>>[CH3][C:1](=[O:5])[C:2]', [s]); let what = 'acetoacetic ester → methyl ketone';
    if (!p.length) { p = rxn('[CX3:1](=[O:5])([OX2][#6])[CX4;!H0:2][CX3](=O)[OX2][#6]>>[C:1](=[O:5])([OX2H1])[C:2]', [s]); what = 'malonic ester → carboxylic acid'; }
    if (!p.length) { p = rxn('[CX3:1](=[O:5])([OX2H1])[CX4;!H0:2][CX3](=O)[OX2H1]>>[C:1](=[O:5])([OX2H1])[C:2]', [s]); what = 'substituted malonic acid → acid'; }
    if (!p.length) return null;
    return finish(m, byC(p[0])[0], { title: `Hydrolysis + decarboxylation: ${what}`, text: 'The esters are hydrolysed; the β-keto acid / malonic acid loses CO₂ through a six-membered cyclic transition state on heating.', extra: ['O=C=O'], tags: ['decarboxylation'] }); } });
})(typeof window !== 'undefined' ? window : globalThis);

;(function(){try{const kill=[],mv=[];for(const k of Object.getOwnPropertyNames(globalThis)){let v;try{v=globalThis[k]}catch(e){continue}if(!v||(typeof v!=='object'&&typeof v!=='function'))continue;for(const key of ['REC','recipes','RECIPES']){let R;try{R=v[key]}catch(e){continue}if(!Array.isArray(R))continue;for(let i=R.length-1;i>=0;i--)if(R[i]&&kill.includes(R[i].id))R.splice(i,1);const mo=[];for(let i=R.length-1;i>=0;i--)if(R[i]&&mv.includes(R[i].id))mo.unshift(R.splice(i,1)[0]);R.push(...mo)}}}catch(e){}})();
