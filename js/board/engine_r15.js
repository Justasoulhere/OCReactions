/* engine_r15.js: run 33. Fixes from the run-32 probes plus chemistry that unseen JEE questions need:
   catalytic hydrogenation of C=O, KMnO4 conventions, Cannizzaro past the aryl-halide guard, acid + heat
   decarboxylation, silver salts → Hunsdiecker, enol ester / enol ether hydrolysis, FeCl3 positive, Gabriel stops at the
   N-alkyl imide without a hydrolysis step, diene regiochemistry (1,4 default), Wittig Z, benzyne regiochemistry,
   α-haloketone + NaBH4 → epoxide, NCERT benzal-chloride route, gem-dihalide hydrolysis, O3 (1 eq), epoxide + PCl5,
   vic-dihalide + NaI, CCl4/KOH carboxylation, excess HCHO (crossed Cannizzaro), cyclic formals, nitrile Knoevenagel,
   D2O exchange, ¹³C-safe haloform, glycosides keep their stereo, E2 by-products. */
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

/* ---------------- catalytic hydrogenation: C=O too ---------------- */
wrap('h2', function (m, c, all, k, run0) {
  if (c.has('Na/NH3') || c.has('Lindlar')) return undefined;
  const co = M(m, '[CX3;$([CX3H1]),$([CX3]([#6])[#6]);!$(C(=O)[O,N,S,Cl,Br,I])]=[OX1]'); if (!co.length) return undefined;
  const r0 = raw(c); const pdOnly = /\bpd\b|pd\/c|pd-c|palladium/.test(r0) && !/\bni\b|raney|\bpt\b|pto2|platinum|\brh\b|nickel/.test(r0);
  const one = c.eq === 1 || /\b1\s*(?:eq|mol|equiv)/.test(r0);
  const hasCC = M(m, '[C;!a]=[C;!a]').length || M(m, '[CX2]#[CX2]').length;
  if (hasCC && (pdOnly || one)) return undefined;
  let base = m, pre = null; if (hasCC) { pre = run0(m, c, all, k); if (pre && pre.product) base = pre.product; }
  const s = smiOf(base); if (!s) return undefined;
  const r = applyAll(s, [['[CX3;$([CX3H1]),$([CX3]([#6])[#6]);!$(C(=O)[O,N,S,Cl,Br,I]):1]=[OX1:2]>>[C:1][O:2]', 'C=O → CH–OH']], 4); if (!r.labels.length) return undefined;
  const ald = M(m, '[CX3H1]=O').length, ket = co.length - ald;
  return finish(m, r.cur, { title: `Catalytic hydrogenation: ${hasCC ? 'C=C and ' : ''}C=O → ${ald && ket ? '1°/2°' : ald ? '1°' : '2°'} alcohol`, text: `H₂ on the metal surface adds across ${hasCC ? 'the C=C and ' : ''}the C=O: aldehydes give 1° alcohols, ketones 2° alcohols.${hasCC ? ' (With Pd/C or 1 eq H₂ only the C=C is reduced.)' : ''}`, cond: ['H₂', 'Ni/Pt/Pd'], why: ['NCERT: aldehydes/ketones are reduced by NaBH₄, LiAlH₄ or catalytic hydrogenation', 'Benzene rings and –COOH/esters stay'], tags: ['red'] });
});
/* ---------------- KMnO4: cold dilute (Baeyer) = diol; otherwise the strong oxidant ---------------- */
retest('diol', (c, all, k, t0) => c.has('OsO4') || (c.has('KMnO4') && (c.cold || c.flag('dil') || c.flag('mild') || /baeyer|0\s*°|273/.test(raw(c))) && !c.heat && !c.has('H2SO4')));
retest('kmno4hot', (c, all, k, t0) => t0(c, all, k) && !(c.has('KMnO4') && !c.has('K2Cr2O7') && !c.heat && (c.cold || c.flag('mild') || (c.flag('dil') && !c.any(['H3O+', 'H2SO4'])) || /baeyer/.test(raw(c)))));
/* aryl-halide "no reaction" must not block the rest of the molecule (Cannizzaro on m-ClC6H4CHO, esters, acids ...) */
wrap('aryl_snar', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (o && o.none && !o.product && (M(m, '[CX3]=O').length || M(m, '[OX2H1]').length || M(m, 'C#N').length)) return null; return o; });
/* β-keto / malonic acids lose CO2 on heating with acid as well */
retest('acid_heat', (c, all, k, t0) => t0(c, all, k) || (c.heat && c.st.sp.size > 0 && [...c.st.sp.keys()].every(x => ['H3O+', 'H2SO4', 'HCl', 'H2O'].includes(x))));
/* silver salts */
front({ id: 'ag_salt', name: 'Acid + AgOH / Ag₂O → silver carboxylate', fam: 'Acid', canon: 'AgOH', db: ['Hunsdiecker reaction'],
  test: c => c.any(['Ag2O', 'AgOH']) || (/agoh|ag2o/.test(raw(c))),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX3](=O)[OX2H1]').length) return null; const p = rxn('[#6:1][CX3:2](=[O:3])[OX2H1]>>[#6:1][C:2](=[O:3])O[Ag]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'RCOOH + AgOH → silver carboxylate (RCOOAg)', text: 'The acid neutralises the silver hydroxide/oxide: RCOO⁻Ag⁺ (dried before the Hunsdiecker step).', tags: ['acidbase'] }); } });
wrap('hunsdiecker', function (m, c, all, k, run0) { const s = smiOf(m); if (!s || !/\[Ag\+\]/.test(s)) return undefined; const cov = canonS(s.replace(/\.\[Ag\+\]/, '').replace('[O-]', 'O[Ag]')); if (!cov) return undefined; const m2 = fresh(cov); return m2 ? run0(m2, c, all, k) : undefined; });
/* enol esters / enol ethers in aqueous acid: hydrolysis, then the enol tautomerises */
front({ id: 'enol_ester_hyd', name: 'Enol ester / enol ether hydrolysis', fam: 'Ester', canon: 'H₃O⁺', db: [],
  test: c => c.any(['H3O+', 'H2SO4', 'HCl', 'HI', 'HBr', 'aqKOH']) || (c.has('H2O') && c.flag('aq')),
  run(m, c) { const s = smiOf(m); if (!s) return null; let p = rxn('[CX3:1]=[CX3:2][OX2][CX3](=O)[#6]>>[C:1][C:2]=O', [s]); let what = 'enol ester', ex = ['CC(=O)O'];
    if (!p.length) { if (c.has('aqKOH') && !c.any(['H3O+', 'HI'])) return null; p = rxn('[CX3:1]=[CX3:2][OX2][CX4:3]>>[C:1][C:2]=O.[C:3]O', [s]); what = 'enol ether'; ex = c.has('HI') ? ['CI'] : ['CO']; }
    if (!p.length) return null; const fr = byC(p[0].join('.').split('.'));
    return finish(m, fr[0], { title: `${what === 'enol ester' ? 'Enol ester' : 'Enol ether (vinyl ether)'} + aqueous acid → carbonyl compound`, text: `The ${what} is cleaved by aqueous acid: the vinyl alcohol formed is an enol and tautomerises at once to the aldehyde/ketone.`, extra: ex, tags: ['hydrolysis'] }); } });
/* FeCl3 colour test is a positive test */
wrap('phenol_rx', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (o && o.none && /Violet colour/.test(o.title || '')) { o.positive = true; o.tags = (o.tags || []).concat('test'); } return o; });
/* Gabriel: with base only (no hydrolysis step anywhere) the product is the N-alkylphthalimide */
wrap('gabriel', function (m, c, all, k, run0) { const later = (all || []).slice(k + 1); const hyd = later.some(s => s.any(['H3O+', 'NH2NH2', 'H2O', 'aqKOH', 'KOH', 'H2SO4', 'HCl'])) || (c.heat && c.any(['aqKOH', 'H2O'])); if (hyd) return undefined;
  const ga = recipe('gabriel_alk'); if (!ga) return undefined; let o = null; try { o = ga.run(m, c, all, k); } catch (e) { o = null; } return o && o.product ? o : undefined; });
/* conjugated dienes: protonate where the allylic cation is most substituted; no temperature given → 1,4 (thermodynamic) */
wrap('diene_add', function (m, c, all, k, run0) { const lowT = (c.T != null && c.T < 273) || /low\s*temp|kinetic|(?:−|-)\s*80|(?<![\d.])0\s*°|(?<!\d)273\s*k|\bcold\b/.test(raw(c)) || c.cold; const hotT = c.heat || (c.T != null && c.T >= 300);
  const s = smiOf(m); if (!s) return undefined; const hx = c.any(['HBr', 'HCl', 'HI']); if (!hx) return undefined; const X = c.has('HBr') ? 'Br' : c.has('HCl') ? 'Cl' : 'I';
    if (!M(m, '[CX3]=[CX3][CX3]=[CX3]').length) return undefined;
  /* choose the protonated end: the cation carbon (C2) with most C neighbours */
  const dm = M(m, '[CX3:1]=[CX3:2][CX3:3]=[CX3:4]'); let best = null; for (const [a, b, cc, d] of dm.flatMap(x => [x, x.slice().reverse()])) { const sc = m.nbrs(b).filter(x => el(m, x) === 'C').length * 10 + hTot(m, a); if (!best || sc > best.sc) best = { sc, a, b, cc, d }; }
  const t12 = m.clone(), t14 = m.clone(); const { a, b, cc, d } = best;
  for (const t of [t12, t14]) { t.bondBetween(a, b).o = 1; [a, b, cc, d].forEach(id => { t.atoms.get(id).hFix = null; }); }
  E.addOn(t12, b, X); t14.bondBetween(cc, d).o = 1; t14.bondBetween(b, cc).o = 2; E.addOn(t14, d, X);
  const s12 = smiOf(t12), s14 = smiOf(t14); if (!s12 || !s14) return undefined; const use14 = !lowT; const main = use14 ? s14 : s12, minor = use14 ? s12 : s14;
  return finish(m, main, { title: `H${X} + conjugated diene → ${use14 ? '1,4' : '1,2'}-addition (${use14 ? 'thermodynamic' : 'kinetic'})`, text: `H⁺ adds to the end that leaves the more substituted allylic cation (+ shared by C2 and C4). ${use14 ? `At room temperature or warmer the more stable product with the internal, more substituted C=C wins: ${X} on C4 (1,4).` : `At low temperature the faster product wins: ${X} on C2 (1,2).`}`, minor: [minor], cond: [use14 ? (hotT ? 'Δ / 40 °C' : 'rt') : 'low T'], why: ['Allylic cation: + charge on C2 and C4', 'Low T → 1,2 (kinetic); rt/40 °C → 1,4 (thermodynamic)'], tags: ['AdE'] });
});
/* Wittig: non-stabilised ylides give mainly the Z-alkene */
wrap('wittig', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || !c.has('Ph3P=CHCH3')) return o; const p = o.product; const cc = M(p, '[CX3H1:1]=[CX3H1:2]'); if (!cc.length || !E.setAlkeneGeom) return o;
  const [a, b] = cc[0]; try { const t = p.clone(); E.setAlkeneGeom(t, a, b, true); o.product = t; o.title = o.title + ' (non-stabilised ylide → mainly Z)'; o.why = ['Non-stabilised ylides (R = alkyl) give mostly the Z (cis) alkene; stabilised ylides (C=O, CN conjugated) give E'].concat(o.why || []); } catch (e) { } return o; });
/* benzyne: −I groups want the new carbanion next to them, alkyl (+I) groups want it far away */
wrap('benzyne', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || !(o.minor || []).length) return o;
  const AX = M(m, '[c:1][Cl,Br,I:2]'); if (!AX.length) return o; const [ca] = AX[0]; const subs = [];
  for (const sm of [['[c:1][CX4H3,CX4H2,CX4H1]', -1], ['[c:1][OX2]', 2], ['[c:1]F', 2], ['[c:1][Cl,Br]', 1], ['[c:1]C(F)(F)F', 2], ['[c:1][NX3]', 1]]) for (const mt of M(m, sm[0])) if (mt[0] !== ca) subs.push([mt[0], sm[1]]);
  if (!subs.length) return o; const ring = C.rings(m, 6).find(r => r.includes(ca)); if (!ring) return o; const pos = id => ring.indexOf(id); const rd = (x, y) => { const d = Math.abs(pos(x) - pos(y)); return Math.min(d, 6 - d); };
  const score = an => subs.reduce((sc, [s0, w]) => { const d = rd(an, s0); return sc + (w > 0 ? (d === 1 ? 2 * w : d === 2 ? w : 0) : (d === 1 ? 2 * w : d === 2 ? w : 0)); }, 0);
  const oh = m.nbrs(ca).filter(n => ring.includes(n) && hTot(m, n) > 0); let best = { sc: score(oh[0] != null ? oh[0] : ca), nh: ca };
  for (const h of oh) { const s1 = score(h); if (s1 > best.sc) best = { sc: s1, nh: ca }; const s2 = score(ca); if (s2 > best.sc) best = { sc: s2, nh: h }; }
  if (best.nh === ca) return o; const t = m.clone(); t.removeAtom(AX[0][1]); t.atoms.get(ca).hFix = null; E.addOn(t, best.nh, 'N'); const s = smiOf(t); if (!s) return o; const R = fresh(s); if (!R) return o;
  const old = o.product; o.product = R; o.minor = [old].concat((o.minor || []).filter(x => smiOf(x) !== s)); o.title = 'NaNH₂ / liq. NH₃: benzyne → the amine goes where the new carbanion is best placed';
  o.why = ['After NH₂⁻ adds to benzyne, the ring carbanion sits next to the other carbon: −I groups (OMe, F, Cl, CF₃) stabilise it when close; alkyl groups (+I) destabilise it, so it forms far from them'].concat(o.why || []); return o; });
/* NaBH4 on an α-halo ketone: the alkoxide closes an epoxide (JEE key) */
front({ id: 'haloketone_nabh4', name: 'α-Halo ketone + NaBH₄ → epoxide', fam: 'Carbonyl', canon: 'NaBH₄, MeOH', db: [],
  test: c => c.has('NaBH4'),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('[#6:1][CX3:2](=[O:3])[CX4:4][Cl,Br,I]>>[#6:1][C:2]1[O:3][C:4]1', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'NaBH₄ reduces the C=O; the alkoxide displaces the α-halide → epoxide', text: 'H⁻ adds to C=O giving an alkoxide right next to C–X; it attacks that carbon from the back (intramolecular SN2) and closes a three-membered ring.', tags: ['red', 'SN2'] }); } });
/* NCERT: ArCH3 + Cl2/hν (2 Cl2) then H2O, 373 K → ArCHO via the benzal chloride */
front({ id: 'benzal_route', name: 'Side-chain dichlorination → benzal chloride', fam: 'Radical', canon: 'Cl₂, hν then H₂O, 373 K', db: ['Benzaldehyde from toluene'],
  test: (c, all, k) => c.has('Cl2') && c.light && !!(all && all[k + 1]) && all[k + 1].any(['H2O']) && !all[k + 1].any(['KOH', 'aqKOH', 'alcKOH', 'NaCN', 'NH3']) && (all[k + 1].heat || /373/.test(String(all[k + 1].raw || ''))),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('[c:1][CH3:2]>>[c:1][CH1:2](Cl)Cl', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Side-chain chlorination (2 Cl₂, hν) → benzal chloride', text: 'NCERT route to benzaldehyde: two benzylic H are replaced by Cl (radical chain); hydrolysis of the gem-dichloride at 373 K gives –CHO.', extra: ['Cl'], tags: ['radical'] }); } });
front({ id: 'gem_hyd', name: 'gem-Dihalide hydrolysis → carbonyl', fam: 'Alkyl halide', canon: 'H₂O, Δ', db: [],
  test: c => c.any(['H2O', 'H3O+', 'aqKOH', 'KOH']) && !c.any(['alcKOH', 'NaNH2', 'Zn', 'KOtBu', 'NaOEt']),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('[CX4;!$(C([Cl,Br,I])([Cl,Br,I])[Cl,Br,I]):1]([Cl,Br,I])[Cl,Br,I]>>[C:1]=O', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'gem-Dihalide + H₂O → aldehyde / ketone', text: 'Both C–X are replaced by OH; a carbon carrying two –OH (gem-diol) loses water at once to give C=O.', extra: ['Cl'], tags: ['hydrolysis'] }); } });
/* ozonolysis with 1 eq O3 stops at the more reactive C=C (alkynes react much more slowly) */
wrap('ozone', function (m, c, all, k, run0) { const st = [c].concat((all || []).slice(k + 1, k + 2)); const one = st.some(x => x && (x.eq === 1 || /\b1\s*(?:eq|mol|equiv)/.test(String(x.raw || '').toLowerCase())));
  if (!one || !M(m, '[C;!a]=[C;!a]').length || !M(m, '[CX2]#[CX2]').length) return undefined; const s = smiOf(m); const p = rxn('[C;!a:1]=[C;!a:2]>>[C:1]=O.[C:2]=O', [s]); if (!p.length) return undefined; const fr = byC(p[0].join('.').split('.'));
  return finish(m, fr[0], { title: 'Ozonolysis (1 eq O₃): only the C=C is cleaved', text: 'O₃ attacks the electron-richer, more reactive C=C much faster than C≡C; with one equivalent the triple bond survives.', extra: fr.slice(1), tags: ['cleavage'] }); });
front({ id: 'epox_pcl5', name: 'Epoxide + PCl₅ → vicinal dichloride', fam: 'Ether', canon: 'PCl₅', db: [],
  test: c => c.has('PCl5'),
  run(m) { const s = smiOf(m); if (!s) return null; const p = rxn('[C:1]1[O][C:2]1>>[C:1](Cl)[C:2]Cl', [s]); if (!p.length) return null; return finish(m, mainFrag(p[0].join('.')), { title: 'Epoxide + PCl₅ → vicinal dichloride', text: 'PCl₅ replaces the ring C–O bonds by C–Cl (POCl₃ forms): a 1,2-dichloride.', extra: ['O=P(Cl)(Cl)Cl'], tags: ['substitution'] }); } });
front({ id: 'vic_nai', name: 'vic-Dihalide + NaI (acetone) → alkene', fam: 'Alkyl halide', canon: 'NaI, acetone', db: [],
  test: c => c.has('NaI') || c.has('Zn'),
  run(m, c) { const s = smiOf(m); if (!s || (c.has('Zn') && !c.has('NaI') && c.any(['HCl', 'H3O+', 'AcOH']))) return null; const p = rxn('[CX4:1]([Cl,Br,I])[CX4:2][Cl,Br,I]>>[C:1]=[C:2]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: `vic-Dihalide + ${c.has('NaI') ? 'NaI (acetone)' : 'Zn'} → alkene (dehalogenation)`, text: c.has('NaI') ? 'I⁻ attacks one halogen; the C–X electrons form the C=C as the other X⁻ leaves (anti elimination); I₂ is set free.' : 'Zn takes both halogens (ZnX₂) and the C=C forms again.', extra: [c.has('NaI') ? 'II' : '[Zn+2]'], tags: ['elimination'] }); } });
/* phenol / pyrrole + CCl4 + KOH → carboxylic acid (Reimer–Tiemann with CCl4) */
front({ id: 'ccl4_rt', name: 'Reimer–Tiemann with CCl₄ → carboxylic acid', fam: 'Phenol', canon: 'CCl₄, KOH', db: ['Reimer–Tiemann reaction'],
  test: c => c.any(['KOH', 'aqKOH', 'alcKOH']) && /ccl4|tetrachloromethane|carbon\s*tetrachloride/.test(raw(c)),
  run(m) { const s = smiOf(m); if (!s) return null; let p = rxn('[OX2H1:1][c:2][cH1:3]>>[OX2H1:1][c:2][c:3]C(=O)O', [s]); let what = 'salicylic acid type (o-COOH)';
    if (!p.length) { p = rxn('[nH1:1][cH1:2]>>[nH1:1][c:2]C(=O)O', [s]); what = 'pyrrole-2-carboxylic acid'; } if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: `CCl₄ + KOH on an activated ring → ${what}`, text: 'CCl₄/OH⁻ gives the electrophile ⁻CCl₃/:CCl₂-type species; it attacks the ring next to O⁻ (or N), and hydrolysis of the CCl₃ group gives –COOH.', tags: ['SEAr'] }); } });
/* excess HCHO / NaOH: every α-H → CH2OH, then crossed Cannizzaro reduces the C=O (pentaerythritol type) */
front({ id: 'hcho_excess', name: 'Excess HCHO, NaOH: hydroxymethylation + crossed Cannizzaro', fam: 'Carbonyl', canon: 'excess HCHO, NaOH, Δ', db: ['Crossed Cannizzaro'],
  test: c => c.has('HCHO') && c.any(['KOH', 'aqKOH', 'NaHCO3']) && (c.eq === Infinity || /excess|xs/.test(raw(c))),
  run(m) { const s = smiOf(m); if (!s || canonS(s) === 'C=O' || !M(m, '[CX4;!H0][CX3]=O').length) return null; const a = applyAll(s, [['[CX4;!H0:1][CX3:2]=[O:3]>>[C:1](CO)[C:2]=[O:3]', 'aldol']], 6); const r = rxn('[CX3;!$(C(=O)[O,N]):1]=[O:2]>>[C:1][O:2]', [a.cur]); if (!r.length) return null;
    return finish(m, mainFrag(r[0].join('.')), { title: 'Excess HCHO + NaOH: every α-H → CH₂OH (aldol), then crossed Cannizzaro reduces C=O', text: 'Each α-H is replaced by –CH₂OH in successive aldol additions with HCHO; with no α-H left, HCHO (the better hydride donor) transfers H⁻ to the C=O and leaves as HCOO⁻ (as in pentaerythritol from CH₃CHO).', extra: ['O=C[O-]'], tags: ['aldol', 'Cannizzaro'] }); } });
front({ id: 'diol_acetal', name: '1,2- / 1,3-Diol + carbonyl (H⁺) → cyclic acetal', fam: 'Alcohol', canon: 'HCHO, H⁺', db: ['Acetal formation'],
  test: c => c.any(['H3O+', 'H2SO4', 'HCl', 'pTsOH', 'TsOH']) && !!partnerCarbonyl(c),
  run(m, c) { const s = smiOf(m); if (!s || M(m, '[CX3]=O').length) return null; const pc = partnerCarbonyl(c); const ck = pc.k === 'HCHO' ? 'C' : pc.k === 'acetone' ? 'C(C)(C)' : pc.k === 'CH3CHO' ? 'C(C)' : null; if (!ck) return null;
    const p = rxn(`[OX2H1:1][CX4:2][CX4:3][OX2H1:4]>>[O:1]1[C:2][C:3][O:4]${ck}1`, [s]).concat(rxn(`[OX2H1:1][CX4:2][CX4:5][CX4:3][OX2H1:4]>>[O:1]1[C:2][C:5][C:3][O:4]${ck}1`, [s])); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Diol + carbonyl, H⁺ → cyclic acetal (1,3-dioxolane / 1,3-dioxane)', text: 'Both OH groups of the 1,2- or 1,3-diol add to the same C=O (acid catalysis, water lost): a five- or six-membered cyclic acetal.', extra: ['O'], tags: ['acetal'] }); } });
/* benzyl cyanide / active nitrile + aldehyde, base → conjugated nitrile */
front({ id: 'nitrile_knoev', name: 'Active nitrile + aldehyde (base) → α,β-unsaturated nitrile', fam: 'Nitrile', canon: 'NaOEt, ArCHO', db: [],
  test: c => c.any(['NaOEt', 'NaOMe', 'KOH', 'aqKOH', 'NaOR', 'pyrrolidine', 'R2NH']) && !!partnerCarbonyl(c),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[CX4;!H0]([$(c),$(C=O),$(C#N)])C#N').length) return null; const pc = partnerCarbonyl(c);
    const p = rxn('[CX4;!H0:1]([C:3]#[N:4])[c,$(C=O),$(C#N):5].[CX3H1,CX3H2:6]=O>>[C:1]([C:3]#[N:4])([*:5])=[C:6]', [s, pc.smi]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Knoevenagel-type condensation: the nitrile α-carbon adds to the aldehyde, then water is lost', text: 'The base removes the acidic H next to –CN (and the ring/C=O); the carbanion adds to the aldehyde and dehydration gives the conjugated C=C.', extra: ['O'], tags: ['aldol', 'C–C'] }); } });
/* enols / ketones with D2O (acid or base): every α-H exchanges */
front({ id: 'd2o_enol', name: 'Enol / carbonyl + D₂O → α-deuteration', fam: 'Carbonyl', canon: 'D₂O, H⁺/OH⁻', db: [],
  test: c => c.any(['D2O', 'ODm']),
  run(m) { const s = smiOf(m); if (!s || M(m, '[#6][Mg,Li,Na,K]').length) return null; let cur = s; const t = rxn('[CX3:1]=[CX3:2][OX2H1:3]>>[C:1][C:2]=[O:3]', [s]); const enol = t.length > 0; if (enol) cur = mainFrag(t[0].join('.'));
    if (!M(fresh(cur), '[CX4;!H0][CX3]=O').length) return null; const r = applyAll(cur, [['[CX4;!H0:1][CX3:2]=[O:3]>>[C:1]([2H])[C:2]=[O:3]', 'α-H → D']], 12); if (!r.labels.length) return null;
    return finish(m, r.cur, { title: `${enol ? 'The enol tautomerises; then e' : 'E'}very α-H exchanges with D (excess D₂O)`, text: 'Acid or base makes the enol/enolate again and again; each time the α-carbon picks up D⁺ from the large excess of D₂O, so all α-H end up as D.', tags: ['tautomerism'] }); } });
/* haloform: keep an isotope label on the CHX3 carbon */
for (const id of ['haloform', 'haloform_fix']) wrap(id, function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product) return o; const me = M(m, '[CH3][CX3]=O')[0] || M(m, '[CH3][CX4][OX2H1]')[0]; if (!me) return o; const iso = (m.atoms.get(me[0]) || {}).iso; if (!iso) return o;
  o.extra = (o.extra || []).map(x => { const s = smiOf(x); const mm = s && s.match(/^(Cl|Br|I)C\((Cl|Br|I)\)(Cl|Br|I)$/); return mm ? (fresh(`${mm[1]}[${iso}CH](${mm[1]})${mm[1]}`) || x) : x; }); return o; });
/* Fischer glycoside from the open-chain sugar keeps C2–C5; α-anomer main (anomeric effect), β minor */
wrap('glycoside', function (m, c, all, k, run0) { const s = isoOf(m); if (!s || !/@/.test(s)) return undefined; const R = c.has('NaOEt') ? 'CC' : 'C';
  const a = rxn(`[CH1:1](=[O])[C:2][C:3][C:4][C:5]([OX2H1:6])[CH2:7]>>${R}[O][C@@H1:1]1[C:2][C:3][C:4][C:5]([CH2:7])[O:6]1`, [s]), b = rxn(`[CH1:1](=[O])[C:2][C:3][C:4][C:5]([OX2H1:6])[CH2:7]>>${R}[O][C@H1:1]1[C:2][C:3][C:4][C:5]([CH2:7])[O:6]1`, [s]);
  if (!a.length) return undefined; return finish(m, mainFrag(a[0].join('.')), { title: 'Sugar + CH₃OH, dry HCl → methyl glycoside (α main, β minor)', text: 'The ring hemiacetal is protonated at the anomeric OH, water leaves (oxocarbenium ion) and CH₃OH adds. The configurations at C2–C5 are untouched; the α-anomer (axial OMe, anomeric effect) is the major product.', minor: b.length ? [mainFrag(b[0].join('.'))] : [], tags: ['acetal'] }); });
/* E2 with an alkoxide: the alcohol is a by-product; small alkoxides also give some ether with 3° halides */
wrap('rx_elim', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || !isAlkoxide(c)) return o; const R = alkoxR(c); const roh = fresh('O' + R); if (roh) o.extra = (o.extra || []).concat([roh]);
  if (!c.has('KOtBu') && M(m, '[CX4]([#6])([#6])([#6])[Cl,Br,I]').length) { const s = smiOf(m); const p = rxn('[CX4:1][Cl,Br,I]>>[C:1]O' + R, [s]); if (p.length) { const e = fresh(mainFrag(p[0].join('.'))); if (e) o.minor = (o.minor || []).concat([e]); } } return o; });
/* r14 corrections: iodine is the iodoform test (non-methyl ketones: negative), and acyl halides beat alkyl halides */
retest('alpha_x_ket', (c, all, k, t0) => t0(c, all, k) && !c.has('I2') && !(c.has('NaOCl') && (c.get('NaOCl') || {}).X === 'I'));
wrap('rx_alkoxide', function (m) { if (M(m, '[CX3](=O)[Cl,Br,I]').length || M(m, '[CX3](=O)O[CX3]=O').length) return null; return undefined; });
/* NH2–Z reagents and amines take the carbonyl first: no acid aldol */
wrap('acid_aldol', function (m, c) { if (c.any(['PhNHNH2', 'NH2OH', 'semicarbazide', 'DNP', 'NH2NH2', 'RNH2', 'NH3', 'aniline', 'R2NH', 'pyrrolidine'])) return null; return undefined; });
/* ethanol / acetaldehyde haloform: the iodoform is the product that matters (the other one is only HCOO⁻) */
for (const id of ['haloform', 'haloform_fix']) wrap(id, function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product) return o; const ps = smiOf(o.product); if (!/^(O=C\[O-\]|O=CO|\[O-\]C=O)$/.test(ps || '')) return o; const x = (o.extra || []).findIndex(e => /^(Cl|Br|I)/.test(smiOf(e) || '') && nC(smiOf(e)) === 1); if (x < 0) return o; const hx = o.extra[x]; o.extra = o.extra.filter((e, i) => i !== x).concat([o.product]); o.product = hx; return o; });
})(typeof window !== 'undefined' ? window : globalThis);
