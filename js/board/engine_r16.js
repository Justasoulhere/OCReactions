/* engine_r16.js: run 33b, generalisation. Written against questions the engine was never tuned on:
   reagent synonyms, multi-partner reactions (crossed aldol, mixed Wurtz, crossed Cannizzaro, crossed Claisen,
   Claisen–Schmidt, Kolbe electrolysis of mixtures, Wurtz–Fittig), and the last benchmark misses. */
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


/* ===== r16 recipes below ===== */

/* ---- residual benchmark fixes ---- */
/* drawn NH2–Z / amine partner: no acid self-aldol of the carbonyl, the partner condenses with it */
wrap('acid_aldol', function (m) { if (m.frags && m.frags().length >= 2 && M(m, '[NX3H2;!$(NC=O);$(N[#7,#8]),$(N[#6])]').length) return null; return undefined; });
/* Reimer–Tiemann: the p-isomer forms too (minor) */
wrap('phenol_rx', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || !c.has('CHCl3') || c.has('CCl4')) return o; const s = smiOf(m); if (!s) return o;
  const p = rxn('[OX2H1,OX1-:9][c:1]1[c:3][c:4][cH1:2][c:5][c:6]1>>[O:9][c:1]1[c:3][c:4][c:2](C=O)[c:5][c:6]1', [s]); if (!p.length) return o; const ps = canonS(mainFrag(p[0].join('.')).replace('[O-]', 'O')); const R = ps && fresh(ps); if (!R || ps === smiOf(o.product)) return o;
  o.minor = (o.minor || []).concat([R]); o.why = (o.why || []).concat(['Some p-hydroxybenzaldehyde forms as well (minor); the ortho product is main']); return o; });
/* EAS: o/p directors give the ortho isomer too; phenol + dilute HNO3 gives mainly o-nitrophenol */
wrap('eas', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || o.none) return o; const s = smiOf(m), ps = smiOf(o.product); if (!s || !ps) return o;
  const G = /^Dilute HNO₃|Nitration/.test(o.title || '') ? '[N+](=O)[O-]' : /Chlorination/.test(o.title || '') ? 'Cl' : /Bromination/.test(o.title || '') ? 'Br' : null; if (!G) return o;
  if (!(o.why || []).some(w => /ortho isomer forms too/.test(w)) && !/^Dilute HNO₃/.test(o.title || '')) return o;
  const orth = rxn(`[c:1](-[CX4,OX2,NX3,F,Cl,Br,I;!$(*=*):9]):[cH1:2]>>[c:1](-[*:9]):[c:2]${G}`, [s]).map(x => canonS(mainFrag(x.join('.')))).filter(x => x && x !== ps); if (!orth.length) return o;
  const oi = orth[0];
  if (/^Dilute HNO₃/.test(o.title || '') && M(m, 'c[OX2H1]').length) { const R = fresh(oi); if (!R) return o;
    return finish(m, oi, { title: 'Dilute HNO₃ (298 K) on phenol → o-nitrophenol (main) + p-nitrophenol', text: 'The –OH ring is so activated that dilute HNO₃ nitrates it without H₂SO₄. Ortho and para isomers both form; the ortho one is the larger share and is steam-volatile (intramolecular O–H···O=N bond), so the two are separated by steam distillation.', minor: [ps], why: ['Measured: ≈ 3–5 : 1 ortho : para with dilute HNO₃ at room temperature', 'p-Nitrophenol is not steam-volatile (intermolecular H-bonds)'], tags: ['SEAr'] }); }
  const R = fresh(oi); if (R) o.minor = (o.minor || []).concat([R]); return o; });
/* ethene with the strong (hot / plain alkaline) KMnO4: both ends become CO2 */
wrap('kmno4hot', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (o && o.product && !o.none) return o; if (!c.has('KMnO4') || !M(m, '[CH2]=[CH2]').length || smiOf(m) !== 'C=C') return o;
  return finish(m, 'O=C=O', { title: 'KMnO₄ (strong): ethene is cleaved completely → CO₂ + H₂O (the purple colour goes)', text: 'Each =CH₂ end is oxidised to HCOOH and then on to CO₂. (Cold, dilute alkaline KMnO₄, Baeyer’s reagent, stops at ethane-1,2-diol.)', extra: ['O'], tags: ['ox', 'cleavage'] }); });
/* SO2Cl2: side-chain chlorination is the JEE answer even when no initiator is written (two keys: plain SO2Cl2 → PhCH2Cl) */
/* CX3–C=O / chloral hydrate + OH⁻: haloform cleavage → CHX3 + carboxylate */
front({ id: 'cx3_cleave', name: 'Trihalomethyl carbonyl + OH⁻ → haloform', fam: 'Carbonyl', canon: 'NaOH (aq)', db: ['Haloform reaction'],
  test: c => c.any(['KOH', 'aqKOH', 'NaOCl', 'Br2/NaOH']) && !c.any(['alcKOH']),
  run(m) { const s = smiOf(m); if (!s) return null; let p = rxn('[CX4:1]([Cl,Br,I:3])([Cl,Br,I:4])([Cl,Br,I:5])[CX4H1:2]([OX2H1])[OX2H1]>>[CH1:1]([*:3])([*:4])[*:5].[C:2](=O)[O-]', [s]);
    if (!p.length) p = rxn('[CX4:1]([Cl,Br,I:3])([Cl,Br,I:4])([Cl,Br,I:5])[CX3:2](=[O:6])[#1,#6;!$([CX4](Cl)(Cl)Cl):7]>>[CH1:1]([*:3])([*:4])[*:5].[C:2](=[O:6])([O-])[*:7]', [s]);
    if (!p.length) p = rxn('[CX4:1]([Cl,Br,I:3])([Cl,Br,I:4])([Cl,Br,I:5])[CX3H1:2]=[O:6]>>[CH1:1]([*:3])([*:4])[*:5].[C:2](=[O:6])[O-]', [s]);
    if (!p.length) return null; const fr = p[0].join('.').split('.'); const hx = fr.find(x => /^(Cl|Br|I)C\((Cl|Br|I)\)(Cl|Br|I)$|^ClC\(Cl\)Cl$|^BrC\(Br\)Br$|^IC\(I\)I$/.test(canonS(x) || '')) || fr[0]; const rest = fr.filter(x => x !== hx);
    const big = rest.length && nC(rest[0]) > 1; const main = big ? rest[0] : hx; const ex = big ? [hx] : rest;
    return finish(m, main, { title: 'CX₃–C=O (or its hydrate) + OH⁻ → CHX₃ + carboxylate (haloform cleavage)', text: 'OH⁻ adds to the C=O; the tetrahedral intermediate throws out ⁻CX₃ (stabilised by three halogens), which takes H⁺: chloral (or its hydrate) gives CHCl₃ + HCOO⁻.', extra: ex, tags: ['haloform'] }); } });

})(typeof window !== 'undefined' ? window : globalThis);
