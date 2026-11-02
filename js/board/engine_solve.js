/* engine_solve.js: reagent text → stages → recipes; multi-step, work-ups, "no reaction" reasons, other ways.
   Also: several drawn molecules (multi-component recipes first, otherwise the biggest reacts and the rest
   join the reagents), equivalents ("2 eq", "excess", "2HBr") and a hint when more sites could react. */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore, RP = root.RxParse;
const WORK = new Set(['H3O+', 'H2O', 'NH4Cl']);
const HYD = ['[CX3]=[CX2]=[OX1]', '[#6][NX2]=[CX2]=[OX1]', '[CX2]#[NX1]', '[N+]#[C-]', '[CX3]=[NX2;!$(N-[OH1]);!$(N-N);!$(N=N)]', '[CX3]=[CX3]-[NX3;H0;!$(N-C=O)]', '[CX3]=[NX3+]', '[OX2H1][CX4][#6]~[#6][CX3](=O)[OX2H1,OX1-]', '[OX2H1][CX4][#6]~[#6]~[#6][CX3](=O)[OX2H1,OX1-]', '[CX3](=O)[Cl,Br,I]', '[c][CX4]([#6])([#6])[OX2][OX2H1]'];
const HYDW = ['[CX3]=[NX2;!$(N-[OH1]);!$(N-N);!$(N=N)]', '[CX3]=[NX3+]', '[CX3](=O)[Cl,Br,I]'];
const HYDR = new Set(['nitrile_rx', 'isocyanide_hyd', 'imine_hyd', 'enamine_hyd', 'iminium_hyd', 'lactonize', 'acylx_hyd', 'hock', 'epox_h3o', 'anilide_h3o', 'cumulene_nu']);
const isWork = st => st.sp.size > 0 && [...st.sp.keys()].every(k => WORK.has(k));
const METAL = new Set(['Mg', 'Li', 'Na', 'K', 'Al', 'B', 'Zn', 'Cu', 'Hg']);
const NOHOT = new Set(['grignard_add', 'grignard_sub', 'hyd_red', 'oxymerc', 'alkoxymerc', 'hydroboration']);
const STOICH = new Set(['hx', 'hbr_roor', 'x2', 'hydration', 'radx', 'rx_oh', 'rx_nu', 'rx_elim', 'ozone', 'h2', 'eas', 'alc_hx', 'alc_pcl', 'acylation', 'amine_rx', 'ether_hi', 'hyd_red', 'epox', 'diol', 'oxymerc', 'hydroboration', 'grignard_add', 'imine', 'ester_hyd', 'alc_ox_mild', 'kmno4hot']);
const CMET = new Set(['Mg', 'Li', 'Na', 'K', 'Zn', 'Cu', 'Cd']);
function neutral(m) {
  const t = m.clone(); let ch = false;
  for (const b of [...t.bonds.values()]) { const A = t.atoms.get(b.a), B = t.atoms.get(b.b); if ((['O', 'N'].includes(A.el) && METAL.has(B.el)) || (['O', 'N'].includes(B.el) && METAL.has(A.el))) { t.removeBond(b.id); ch = true; } }
  /* R–MgX / R–Li / R₂Cd left at the work-up: the C–metal bond takes H⁺ */
  for (const b of [...t.bonds.values()]) { const A = t.atoms.get(b.a), B = t.atoms.get(b.b); const c = A.el === 'C' && CMET.has(B.el) ? A : B.el === 'C' && CMET.has(A.el) ? B : null; if (c) { t.removeBond(b.id); c.hFix = null; c.chg = 0; c.rad = 0; ch = true; } }
  for (const a of t.atoms.values()) if (a.chg < 0 && ['O', 'N', 'S'].includes(a.el) && !t.nbrs(a.id).some(n => (t.atoms.get(n).chg || 0) > 0)) { a.chg = 0; a.hFix = null; ch = true; }
  if (!ch) return null;
  const fr = t.frags().map(ids => ({ ids, n: ids.filter(id => t.atoms.get(id).el === 'C').length })).sort((p, q) => q.n - p.n);
  return E.cleanProduct(t, fr[0].ids);
}
function recognised(st) { return [...st.sp.keys()].concat([...st.flags]); }
function stripStereo(s) { return s.replace(/[\/\\]/g, '').replace(/@+/g, ''); }
function key(m) { try { const c = C.canon(m); const mol = C.RD.get_mol(stripStereo(c)); const k = mol.get_smiles(); mol.delete(); return k; } catch (e) { return null; } }

/* ---- several molecules on the canvas ---- */
let REV = null;
function reverseAliases() { // canonical SMILES of a drawn small molecule -> parser alias entry
  if (REV) return REV; REV = new Map(); const pre = { RX: v => (v.X || 'Cl') + v.smi, RMgX: v => (v.X || 'Br') + '[Mg]' + v.smi, ROH: v => 'O' + v.smi, RNH2: v => 'N' + v.smi, RCOCl: v => 'ClC(=O)' + v.smi, R2NH: v => v.smi + 'N' + v.smi, RLi: v => '[Li]' + v.smi, NaOR: v => '[O-]' + v.smi }; const RSKIP = new Set(['R2CuLi', 'R2Cd', 'org']);
  const put = (smi, v) => { try { const mol = C.RD.get_mol(smi); if (!mol) return; const k = mol.get_smiles(); mol.delete(); if (!REV.has(k)) REV.set(k, v); } catch (e) { } };
  for (const s of RP.SP || []) if (s.smi && !pre[s.id]) put(s.smi, { sp: s.id, smi: s.smi, txt: s.al[0] });
  for (const [a, v] of RP.ALIAS) { if (!v.sp || RSKIP.has(v.sp)) continue; if (pre[v.sp] && v.smi && !v.gen) put(pre[v.sp](v), Object.assign({}, v, { txt: a })); else if (v.smi && !pre[v.sp]) put(v.smi, Object.assign({}, v, { txt: a })); }
  return REV;
}
const SPEC = { 'Br': 'Br2', 'ClCl': 'Cl2', 'II': 'I2', 'BrBr': 'Br2', 'Cl': 'HCl', 'Br_': 'HBr', 'I': 'HI', 'O': 'H2O', 'N': 'NH3', '[HH]': 'H2', 'O=C=O': 'CO2', '[C-]#[O+]': 'CO', 'ClC(Cl)Cl': 'CHCl3', 'NN': 'NH2NH2', 'NO': 'NH2OH', 'C1CO1': 'ethyleneoxide', 'OCCO': 'glycol' };
function heavyCount(m, ids) { return ids.filter(id => m.atoms.get(id).el !== 'H').length; }
function cCount(m, ids) { return ids.filter(id => m.atoms.get(id).el === 'C').length; }
function rootedR(m, ids, att, drop) { const keep = ids.filter(id => !drop.includes(id)); const order = [att, ...keep.filter(id => id !== att)]; try { return C.toSmiles(m, { ids: order }).smi.replace(/\[(C|O|S|Cl|Br|I|F)H?\d*\]/g, '$1'); } catch (e) { return null; } }
let SPEC2 = null; function spec2() { if (SPEC2) return SPEC2; SPEC2 = new Map(); const L = { 'NNc1ccccc1': 'PhNHNH2', 'NC(=O)NN': 'semicarbazide', 'NNc1ccc([N+](=O)[O-])cc1[N+](=O)[O-]': 'DNP', 'NN': 'NH2NH2', 'NO': 'NH2OH', 'N#C[Na]': 'NaCN', '[C-]#N.[Na+]': 'NaCN', '[C-]#N.[K+]': 'NaCN', 'C#N': 'HCN', 'C1CCNC1': 'pyrrolidine', 'C1CCNCC1': 'R2NH', 'Nc1ccccc1': 'aniline' };
  for (const [smi, sp] of Object.entries(L)) { try { const q = C.fromSmiles(smi, 40); const k = q && C.canon(q); if (k) SPEC2.set(k, sp); } catch (e) { } } return SPEC2; }
function classify(m, ids) { // -> { alias entry } to inject into the first stage, or null
  const sub = m.sub(ids); let k = null; try { k = C.canon(sub); } catch (e) { }
  if (k) { const s2 = spec2().get(k); if (s2 && s2 !== 'aniline') return { sp: s2, txt: 'drawn' }; const rv = reverseAliases().get(k); if (rv) return rv; const sp = SPEC[k]; if (sp) return { sp }; }
  const els = ids.map(id => m.atoms.get(id).el); const A_ = id => m.atoms.get(id);
  /* run 32: organometallics, alkoxides, diazonium salts, and halides with one sp3 C–X plus aryl halogens */
  const mg = ids.find(id => A_(id).el === 'Mg'); if (mg != null) { const cc = m.nbrs(mg).find(n => A_(n).el === 'C'); const xx = m.nbrs(mg).find(n => ['Cl', 'Br', 'I'].includes(A_(n).el)); if (cc != null) { const smi = rootedR(m, ids, cc, [mg, xx].filter(v => v != null)); if (smi) return { sp: 'RMgX', smi, X: xx != null ? A_(xx).el : 'Br', txt: 'drawn' }; } }
  const li = ids.find(id => A_(id).el === 'Li' && m.nbrs(id).some(n => A_(n).el === 'C')); if (li != null) { const cc = m.nbrs(li).find(n => A_(n).el === 'C'); const smi = rootedR(m, ids, cc, [li]); if (smi) return { sp: 'RLi', smi, txt: 'drawn' }; }
  const cu = ids.find(id => A_(id).el === 'Cu'); if (cu != null) { const cs = m.nbrs(cu).filter(n => A_(n).el === 'C'); if (cs.length) { const smi = rootedR(m, E.fragIds(m, cs[0]).filter(x => ids.includes(x)), cs[0], [cu, ...ids.filter(x => A_(x).el === 'Li')]); if (smi) return { sp: 'R2CuLi', smi: smi.split('.')[0], txt: 'drawn' }; } }
  { const om = ids.filter(id => A_(id).el === 'O' && A_(id).chg === -1 && m.nbrs(id).length === 1); const oo = ids.filter(id => A_(id).el === 'O');
    if (om.length === 1 && oo.length === 1 && els.every(e => ['C', 'H', 'O'].includes(e))) { const cc = m.nbrs(om[0])[0]; if (!E.isCarbonylC(m, cc)) { const smi = rootedR(m, ids, cc, [om[0]]); if (smi === 'C') return { sp: 'NaOMe', txt: 'drawn' }; if (smi === 'CC') return { sp: 'NaOEt', txt: 'drawn' }; if (smi === 'C(C)(C)C') return { sp: 'KOtBu', txt: 'drawn' }; const ar = E.aromSet(m).has(cc); if (smi) return { sp: 'NaOR', smi, aryl: ar, txt: 'drawn' }; } } }
  { const dz = ids.find(id => A_(id).el === 'N' && A_(id).chg === 1 && m.bondsOf(id).some(b => b.o === 3)); if (dz != null) { const cc = m.nbrs(dz).find(n => A_(n).el === 'C'); const n2 = m.nbrs(dz).find(n => A_(n).el === 'N'); if (cc != null && n2 != null) { const smi = rootedR(m, ids, cc, [dz, n2]); if (smi) return { sp: 'ArN2', smi, txt: 'drawn' }; } } }
  const hal = ids.filter(id => ['Cl', 'Br', 'I'].includes(m.atoms.get(id).el));
  if (hal.length > 1 && els.every(e => ['C', 'H', 'Cl', 'Br', 'I', 'F', 'O', 'N'].includes(e))) { const sp3 = hal.filter(x => { const cc = m.nbrs(x)[0]; return cc != null && A_(cc).el === 'C' && m.bondsOf(cc).every(b => b.o === 1) && m.nbrs(cc).filter(n => ['Cl', 'Br', 'I'].includes(A_(n).el)).length === 1; });
    if (sp3.length >= 1 && new Set(sp3.map(x => m.nbrs(x)[0])).size === sp3.length && !ids.some(id => A_(id).el === 'O' || A_(id).el === 'N')) { const x = sp3[0]; const cc = m.nbrs(x)[0]; const smi = rootedR(m, ids, cc, [x]); if (smi) return { sp: 'RX', smi, X: A_(x).el, txt: 'drawn' }; } }
  if (hal.length === 1 && els.every(e => ['C', 'H', 'Cl', 'Br', 'I'].includes(e))) { const x = hal[0]; const cc = m.nbrs(x)[0]; const sp3 = m.bondsOf(cc).every(b => b.o === 1); const smi = rootedR(m, ids, cc, [x]); if (smi && sp3) return { sp: 'RX', smi, X: m.atoms.get(x).el, txt: 'drawn' }; }
  const o = ids.filter(id => m.atoms.get(id).el === 'O');
  if (o.length === 1 && m.nbrs(o[0]).length === 1 && els.every(e => ['C', 'H', 'O'].includes(e)) && m.bondsOf(o[0])[0].o === 1) { const cc = m.nbrs(o[0])[0]; if (!E.isCarbonylC(m, cc)) { const smi = rootedR(m, ids, cc, [o[0]]); if (smi) return { sp: 'ROH', smi, txt: 'drawn' }; } }
  const n = ids.filter(id => m.atoms.get(id).el === 'N');
  if (n.length === 1 && els.every(e => ['C', 'H', 'N'].includes(e)) && m.nbrs(n[0]).length === 1) { const cc = m.nbrs(n[0])[0]; const smi = rootedR(m, ids, cc, [n[0]]); if (smi) return { sp: 'RNH2', smi, txt: 'drawn' }; }
  return null;
}
function splitInput(m) {
  const MSP = { Mg: 'Mg', Na: 'Na', K: 'Na', Li: 'Li', Zn: 'Zn', Sn: 'Sn', Fe: 'Fe', Cu: 'Cu' };
  const all = m.frags(); const metals = all.filter(ids => ids.length === 1 && MSP[m.atoms.get(ids[0]).el] && !m.atoms.get(ids[0]).chg);
  const fr = all.filter(ids => !(ids.length === 1 && (METAL.has(m.atoms.get(ids[0]).el) || MSP[m.atoms.get(ids[0]).el] || m.atoms.get(ids[0]).chg)));
  if (!fr.length || fr.length + metals.length < 2) return null;
  const RROLE = new Set(['PhNHNH2', 'NH2OH', 'semicarbazide', 'DNP', 'NH2NH2', 'ArN2', 'NaOEt', 'NaOMe', 'KOtBu', 'NaOR', 'RMgX', 'RLi', 'R2CuLi', 'HCN', 'NaCN']);
  const role = ids => { try { const v = classify(m, ids); return v && RROLE.has(v.sp) ? 1 : 0; } catch (e) { return 0; } };
  const rk = new Map(fr.map(ids => [ids, role(ids)]));
  fr.sort((p, q) => (rk.get(p) - rk.get(q)) || cCount(m, q) - cCount(m, p) || heavyCount(m, q) - heavyCount(m, p));
  const main = E.cleanProduct(m, fr[0]); const co = [], left = [];
  for (const ids of metals) { const e = m.atoms.get(ids[0]).el; co.push({ sp: MSP[e], txt: e }); }
  for (const ids of fr.slice(1)) { const v = classify(m, ids); if (v) co.push(v); else left.push(E.cleanProduct(m, ids)); }
  return { main, co, left };
}
/* several drawn molecules: when the biggest one does not react with the first reagent, try each of the others as the
   substrate (the rest become reagents) and keep the drawn partners for the later steps (aniline + β-naphthol, then NaOH) */
function altMain(m, ctx0, all, k) {
  const fr = m.frags(); if (fr.length < 2) return null;
  const order = fr.slice().sort((p, q) => cCount(m, q) - cCount(m, p) || heavyCount(m, q) - heavyCount(m, p));
  for (let i = 1; i < order.length; i++) {
    const ids = order[i]; if (cCount(m, ids) < 1) continue; const main = E.cleanProduct(m, ids);
    const st = Object.assign({}, ctx0.st, { sp: new Map(ctx0.st.sp), flags: new Set(ctx0.st.flags) }); const carry = [];
    for (let j = 0; j < order.length; j++) { if (j === i) continue; const v = classify(m, order[j]); if (v) { inject(st, v); carry.push(v); } else { const sub = E.cleanProduct(m, order[j]); let smi = null; try { smi = C.canon(sub); } catch (e) { } if (smi) { const v2 = { sp: 'partner', smi, txt: 'drawn' }; inject(st, v2); carry.push(v2); } } }
    const ctx = RP.stageCtx(st); const cands = E.REC.filter(r => { try { return r.test(ctx, all, k); } catch (e) { return false; } });
    const [o, rec] = runStage(main, ctx, all, k, cands);
    if (o && o.product && !o.none) return { o, rec, cands, main, ctx, carry, note: 'The smaller drawn molecule reacts with this reagent first; the other one is kept for the next step.' };
  }
  return null;
}
function inject(st, v) { const prev = st.sp.get(v.sp); if (!prev) st.sp.set(v.sp, v); else if (!Array.isArray(prev)) st.sp.set(v.sp, [prev, v]); else prev.push(v); }

function keeps(a, b) { const ha = [...a.atoms.values()].filter(x => x.el !== 'H'); if (!ha.length) return false; const kept = ha.filter(x => b.atoms.has(x.id)).length; return kept / ha.length >= 0.6; }
function runStage(cur, ctx, all, k, cands) {
  for (const r of cands) { let o = null; try { o = r.run(cur, ctx, all, k); } catch (e) { if (root.console) console.warn('recipe ' + r.id + ' failed:', e && e.message); o = null; } if (o) return [o, r]; }
  return [null, null];
}
/* run 32: organic co-reagents typed as names / formulas ('org' entries from the parser): reagent roles (R–X, R–MgX, RO⁻,
   ArN₂⁺ ...) join the stage; anything else in the first stage is added to the drawing as another molecule, in later
   stages it becomes a 'partner' */
function mergeInto(m, q) { const xs = [...m.atoms.values()].map(a => a.x), ys = [...m.atoms.values()].map(a => a.y); const qx = [...q.atoms.values()].map(a => a.x), qy = [...q.atoms.values()].map(a => a.y);
  const dx = m.atoms.size ? Math.max(...xs) - Math.min(...qx) + 3 * E.BL : 0, dy = m.atoms.size ? (Math.min(...ys) + Math.max(...ys)) / 2 - (Math.min(...qy) + Math.max(...qy)) / 2 : 0; m.merge(q, dx, dy); }
function resolveOrg(P, m) {
  let mm = null; const added = []; const drawn = new Set(); try { for (const ids of m.frags()) { const q = C.canon(E.cleanProduct(m, ids)); if (q) drawn.add(q); } } catch (e) { }
  if (!m.atoms.size && P.stages.length) { const st = P.stages[0]; const v = st.sp.get('org'); const list = v ? (Array.isArray(v) ? v : [v]) : [];
    if (list.length) { const first = list[0]; const q = (() => { try { return C.fromSmilesMB(first.smi, E.BL); } catch (e) { return null; } })(); if (q && q.atoms.size) { mm = q; added.push(first.txt); drawn.add(C.canon(q)); if (list.length === 1) st.sp.delete('org'); else st.sp.set('org', list.length === 2 ? list[1] : list.slice(1)); } }
    else if (st.sp.get('RMgX') && !mm) { const e = st.sp.get('RMgX'); const e1 = Array.isArray(e) ? e[0] : e; if (e1 && e1.smi) { const smi = '[Mg](' + (e1.X || 'Br') + ')' + e1.smi.replace(/^C1=CC=CC=C1$/, 'c1ccccc1'); const q = (() => { try { return C.fromSmilesMB(smi, E.BL); } catch (er) { return null; } })(); if (q && q.atoms.size) { mm = q; added.push(e1.txt || 'RMgX'); if (Array.isArray(e)) { e.shift(); if (!e.length) st.sp.delete('RMgX'); } else st.sp.delete('RMgX'); } } }
    if (!mm && !list.length) { for (const key of ['RX', 'ROH', 'RNH2']) { const e = st.sp.get(key); const e1 = Array.isArray(e) ? e[0] : e; if (e1 && e1.txt && e1.smi && root.OCName) { const r = root.OCName.read(e1.txt); if (r && r.smi) { const q = (() => { try { return C.fromSmilesMB(r.smi, E.BL); } catch (er) { return null; } })(); if (q && q.atoms.size) { mm = q; added.push(e1.txt); if (Array.isArray(e)) { e.shift(); if (!e.length) st.sp.delete(key); } else st.sp.delete(key); break; } } } } } }
  P.stages.forEach((st, k) => { const v = st.sp.get('org'); if (!v) return; st.sp.delete('org');
    for (const e of (Array.isArray(v) ? v : [v])) { let q = null; try { q = C.fromSmilesMB(e.smi, E.BL); } catch (er) { q = null; } if (!q || !q.atoms.size) continue; { let qs = null; try { qs = C.canon(q); } catch (er) { qs = null; } if (qs && drawn.has(qs)) continue; }
      const frs = q.frags().filter(ids => ids.some(id => q.atoms.get(id).el === 'C')).sort((a, b) => b.length - a.length); const ids = frs[0] || [...q.atoms.keys()];
      let cl = null; try { cl = classify(q, ids); } catch (er) { cl = null; }
      if (cl) { inject(st, Object.assign({}, cl, { txt: e.txt || cl.txt })); continue; }
      if (k === 0) { if (!mm) mm = m.clone(); let q2 = q; if (frs.length && frs[0].length < q.atoms.size) { try { q2 = E.cleanProduct(q, frs[0]); } catch (er) { q2 = q; } } mergeInto(mm, q2); added.push(e.txt); }
      else inject(st, { sp: 'partner', smi: e.smi, txt: e.txt }); } });
  return { m: mm || m, added, emptyStart: !m.atoms.size };
}
function mergeSt(a, b) { const sp = new Map(a.sp); for (const [k, v] of b.sp) { if (!sp.has(k)) sp.set(k, v); } return Object.assign({}, b, { raw: (a.raw || '').trim() + ', ' + (b.raw || '').trim(), sp, flags: new Set([...a.flags, ...b.flags]), T: b.T != null ? b.T : a.T, eq: b.eq != null ? b.eq : a.eq, unknown: [...(a.unknown || []), ...(b.unknown || [])] }); }
function solve(m, text, opt = {}) {
  const P = RP.parse(text || '');
  const RO = resolveOrg(P, m); m = RO.m;
  let stages = P.stages.filter((s, i) => s.sp.size || s.flags.size || s.T != null || (i === 0 && RO.added.length));
  if (!stages.length && m.frags().length >= 2 && P.stages.length) stages = [P.stages[0]];
  if (!stages.length) return { status: 'empty', unknown: P.unknown };
  /* "H₃O⁺ then heat" is one hot aqueous-acid step (ester hydrolysis + decarboxylation) */
  for (let i = 0; i + 1 < stages.length; i++) { const a = stages[i], b = stages[i + 1]; if (i > 0 && isWork(a) && !b.sp.size && b.flags.has('heat')) { stages.splice(i, 2, mergeSt(a, b)); i--; } }
  let all = stages.map(RP.stageCtx);
  let cur = m.clone(); const steps = [], titles = [], why = [], tags = [], used = []; let extra = [], minor = []; const extraAll = []; const notes = []; let more = null, eqInfo = null;
  if (RO.added.length) notes.push(m.atoms.size && RO.added.length && !RO.emptyStart ? `Read ${RO.added.map(x => '"' + x + '"').join(', ')} as ${RO.added.length > 1 ? 'molecules' : 'a molecule'} reacting with the drawn one` : `Nothing was drawn: "${RO.added[0]}" is used as the starting molecule`);
  let multi = cur.frags().length > 1; let carry = [];
  for (let k = 0; k < all.length; k++) {
    let ctx = all[k];
    /* "NaOH then I₂" on a methyl ketone / methyl carbinol is one haloform step */
    if (k + 1 < all.length && ctx.st.sp.size && [...ctx.st.sp.keys()].every(x => ['KOH', 'aqKOH'].includes(x)) && !ctx.heat && ['I2', 'Br2', 'Cl2'].some(x => all[k + 1].st.sp.has(x)) && (() => { try { return C.matchAll(cur, '[CH3][CX3]=O', {}).length + C.matchAll(cur, '[CH3][CX4H1][OX2H1]', {}).length > 0; } catch (e) { return false; } })()) { all = all.slice(0, k).concat([RP.stageCtx(mergeSt(ctx.st, all[k + 1].st))], all.slice(k + 2)); ctx = all[k]; }
    if (k > 0 && isWork(ctx.st) && !(ctx.heat || ctx.flag('warm') || ctx.flag('boil'))) {
      const n = neutral(cur); if (n) { steps.push({ L: cur.clone(), arrows: [], R: n.clone(), text: 'Work-up: the anion picks up H⁺ (and the metal salt washes away).', cond: ['H₃O⁺'] }); cur = n; }
      /* JEE convention: "then H₃O⁺" after a nitrile / isocyanide / imine means it is hydrolysed */
      const prevBase = k > 0 && all[k - 1] && all[k - 1].any(['KOH', 'aqKOH', 'alcKOH', 'NaOEt', 'NaOMe', 'NaOR', 'KOtBu', 'NaH', 'NaNH2', 'RMgX', 'RLi', 'LiAlH4', 'NaBH4', 'DIBAL', 'R2CuLi']);
      const protN = used.some(r => r && /^(acylation|acyl_amine|amine_acyl|ac2o_amine|acetyl)/.test(r.id)) ? ['[c][NX3H1][CX3](=O)[CH3]'] : [];
      if ((ctx.has('H3O+') ? HYD.concat(prevBase ? [] : ['[CX4;r3]1[OX2;r3][CX4;r3]1'], protN) : ctx.has('H2O') ? HYDW : []).some(sm => { try { return C.matchAll(cur, sm, {}).length > 0; } catch (e) { return false; } })) {
        const hst = Object.assign({}, ctx.st, { sp: new Map([...ctx.st.sp, ['H3O+', { sp: 'H3O+' }]]), flags: new Set([...ctx.st.flags, 'heat']) }); const hc = RP.stageCtx(hst);
        const hc2 = E.REC.filter(r => HYDR.has(r.id) && (() => { try { return r.test(hc, all, k); } catch (e) { return false; } })());
        const [ho, hr] = runStage(cur, hc, all, k, hc2);
        if (ho && ho.product && !ho.none) { steps.push(...ho.steps); titles.push(ho.title); why.push(...(ho.why || [])); used.push(hr); extra = ho.extra || []; extraAll.push(...extra); minor = ho.minor || []; cur = ho.product; }
      }
      continue; }
    if (!ctx.st.sp.size && !k && !ctx.st.flags.size && ctx.T == null && !multi) return { status: 'unknown', unknown: P.unknown, stage: ctx.raw };
    let o = null, rec = null, cands = [];
    const orig0 = k === 0 && multi ? cur.clone() : null; const ctx00 = k === 0 && multi ? RP.stageCtx(Object.assign({}, ctx.st, { sp: new Map(ctx.st.sp), flags: new Set(ctx.st.flags) })) : null;
    if (k === 0 && multi) {
      const mc = E.REC.filter(r => r.multi && (() => { try { return r.test(ctx, all, k); } catch (e) { return false; } })());
      [o, rec] = runStage(cur, ctx, all, k, mc);
      if (!o) { const sp = splitInput(cur); if (sp) { cur = sp.main; for (const v of sp.co) inject(ctx.st, v); const parts = sp.left.map(x => { try { return C.canon(x); } catch (e) { return null; } }).filter(Boolean).map(smi => ({ sp: 'partner', smi, txt: 'drawn' })); for (const v of parts) inject(ctx.st, v); carry = parts.slice(); all[0] = ctx = RP.stageCtx(ctx.st); if (sp.co.length) notes.push(`The other drawn molecule${sp.co.length > 1 ? 's were' : ' was'} read as reagent${sp.co.length > 1 ? 's' : ''} (${sp.co.map(v => v.txt === 'drawn' ? v.sp : v.txt).join(', ')}).`);  } }
      multi = false;
    }
    if (k > 0 && carry.length) { for (const v of carry) if (!ctx.st.sp.has(v.sp)) inject(ctx.st, v); all[k] = ctx = RP.stageCtx(ctx.st); }
    if (!o) { cands = E.REC.filter(r => { try { return r.test(ctx, all, k); } catch (e) { return false; } }); [o, rec] = runStage(cur, ctx, all, k, cands); }
    if (!o && orig0) { const alt = altMain(orig0, ctx00, all, k); if (alt) { o = alt.o; rec = alt.rec; cands = alt.cands; carry = alt.carry; cur = alt.main; all[0] = ctx = alt.ctx; if (alt.note) notes.push(alt.note); } }
    /* an anion from the previous step (carboxylate, phenoxide, alkoxide) is neutralised first if that is what the next reagent needs */
    if (!o && k > 0) { const n = neutral(cur); if (n) { const c2 = E.REC.filter(r => { try { return r.test(ctx, all, k); } catch (e) { return false; } }); const [o2, r2] = runStage(n, ctx, all, k, c2); if (o2 && o2.product && !o2.none) { steps.push({ L: cur.clone(), arrows: [], R: n.clone(), text: 'Acidify first: the anion picks up H⁺.', cond: ['H₃O⁺'] }); cur = n; o = o2; rec = r2; cands = c2; } } }
    /* run 32: a step that does nothing on its own is read together with the next one ("NaOEt then EtI", "NaOH then I₂", "Na then D₂O") */
    if (!o && k + 1 < all.length && !isWork(all[k + 1].st)) { const cm = RP.stageCtx(mergeSt(ctx.st, all[k + 1].st)); const c2 = E.REC.filter(r => { try { return r.test(cm, all, k); } catch (e) { return false; } });
      let [o2, r2] = runStage(cur, cm, all, k, c2); if ((!o2 || !o2.product || o2.none) && k > 0) { const n = neutral(cur); if (n) [o2, r2] = runStage(n, cm, all, k, c2); }
      if (o2 && o2.product && !o2.none) { notes.push(`Steps ${k + 1} and ${k + 2} (${ctx.raw.trim()} | ${all[k + 1].raw.trim()}) were read as one step`); all = all.slice(0, k).concat([cm], all.slice(k + 2)); ctx = cm; o = o2; rec = r2; cands = c2; } }
    if (!o && k > 0 && !P.unknown.length) { why.push(`Step ${k + 1} (${ctx.raw.trim()}) leaves this structure unchanged`); titles.push(`no change with ${ctx.raw.trim()}`); continue; }
    if (!o) {
      const sp = recognised(ctx.st);
      if (P.unknown.length) return { status: 'unknown', unknown: P.unknown, known: sp, stage: ctx.raw, reason: `I can't read "${P.unknown.join('", "')}".` };
      return { status: 'none', reason: cands.length ? `${ctx.raw.trim()} has no reaction with this structure in the library or the JEE syllabus.` : `No reaction of this structure with ${ctx.raw.trim()} that I know of.`, stage: ctx.raw, steps, titles, partial: k ? cur : null, known: sp };
    }
    if ((o.none || !o.product) && k > 0 && !(o.tags || []).includes('test') && !/\btest\b/i.test(o.title || '') && !o.positive && titles.length) { why.push(`Step ${k + 1} (${ctx.raw.trim()}): ${o.none || 'no reaction'}`); titles.push(`no change with ${ctx.raw.trim()}`); continue; }
    if (o.none || !o.product) return { status: 'none', reason: o.none || 'No reaction.', title: o.title, steps, titles, test: (o.tags || []).includes('test'), positive: !!o.positive };
    let title = o.title; const reps = [o];
    // equivalents: repeat the same stage on the product
    const eq = ctx.eq; if (!cands.length) cands = E.REC.filter(r => { try { return r.test(ctx, all, k); } catch (e) { return false; } });
    const lim = eq == null ? 1 : eq === Infinity ? 6 : Math.max(1, Math.min(6, Math.round(eq)));
    let p = o.product, pk = key(p), n = o.eqUsed || 1;
    while (n < lim && !o.stop) { const [o2, r2] = runStage(p, ctx, all, k, cands); if (!o2 || o2.none || !o2.product || !r2 || !(r2.stoich || STOICH.has(r2.id) || r2 === rec) || !keeps(p, o2.product)) break; const k2 = key(o2.product); if (!k2 || k2 === pk) break; reps.push(o2); p = o2.product; pk = k2; n++; if (o2.stop) break; }
    if (eq != null) eqInfo = { stage: k, asked: eq, used: n };
    else if (!opt.noMore && !o.stop && (rec.stoich || STOICH.has(rec.id))) { const [o2, r2] = runStage(p, ctx, all, k, cands); if (o2 && !o2.none && o2.product && key(o2.product) !== pk && r2 && (r2.stoich || STOICH.has(r2.id)) && keeps(p, o2.product)) more = { stage: k, title: o2.title }; }
    for (const x of reps) { steps.push(...x.steps); why.push(...(x.why || []).filter(Boolean)); tags.push(...(x.tags || [])); }
    if (reps.length > 1) { const ts = [...new Set(reps.map(x => x.title))]; title = ts.length === 1 ? `${ts[0]} (×${reps.length})` : ts.join(', then '); }
    if (n > 1) why.push(`${n} equivalents reacted`);
    if (eq != null && eq !== Infinity && n < lim) why.push(`Only ${n} equivalent${n > 1 ? 's' : ''} can react here; the rest stays unused`);
    titles.push(title); used.push(rec); const last = reps[reps.length - 1]; extra = last.extra || []; for (const x of reps) extraAll.push(...(x.extra || [])); minor = (reps.length > 1 ? last.minor : o.minor) || [];
    if (k === 0 && carry.length) { const hv = x => [...x.atoms.values()].filter(a => a.el !== 'H').length; if (hv(p) > hv(cur) + 1) carry = []; }
    cur = p;
    while (rec.eats && k + 1 < all.length && !(all[k + 1].heat && NOHOT.has(rec.id)) && (() => { try { return rec.eats(all[k + 1]); } catch (e) { return false; } })()) k++;
  }
  /* protecting group: a ketal made in an earlier stage comes off in a later aqueous-acid work-up */
  { const pi = used.findIndex(r => r && (r.id === 'acetal' || r.id === 'acetal_roh'));
    if (pi >= 0 && pi < used.length - 1 && all.slice(1).some(x => x.any(['H3O+', 'H2O']))) { const hy = E.REC.find(r => r.id === 'acetal_hyd');
      let o = null; try { o = hy && hy.run(cur, all[all.length - 1], all, all.length - 1); } catch (e) { o = null; }
      if (o && o.product && !o.none) { steps.push(...o.steps.map(s => Object.assign({}, s, { text: s.text }))); if (o.steps.length) steps[steps.length - 1] = Object.assign({}, steps[steps.length - 1], { text: (steps[steps.length - 1].text || '') + ' The aqueous acid work-up also removes the ketal protecting group: the C=O is back.' }); cur = o.product; titles.push('deprotection'); why.push('The ketal protected the C=O during the step that would otherwise have attacked it'); } } }
  let ok = true; try { ok = C.isValid(cur); } catch (e) { ok = false; }
  const exA = []; { const seen = new Set(); for (const x of extraAll) { let kk = null; try { kk = C.canon(x); } catch (e) { kk = null; } if (kk && !seen.has(kk)) { seen.add(kk); exA.push(x); } } }
  return { status: ok ? 'ok' : 'error', product: cur, extra, extraAll: exA, minor, steps, titles, title: titles.join(', then '), why: [...new Set(why.concat(notes))], tags, recipes: used, unknown: P.unknown, more, eq: eqInfo };
}
function ways(m, product, skip) {
  const want = key(product); if (!want) return [];
  const out = [], seen = new Set();
  for (const r of E.REC) {
    if (!r.canon || (skip && skip.includes(r))) continue;
    let res; try { res = solve(m, r.canon, { noMore: true }); } catch (e) { continue; } if (res.status !== 'ok') continue;
    if (res.recipes && res.recipes[0] !== r) continue;
    if (key(res.product) === want && !seen.has(r.canon)) { seen.add(r.canon); out.push({ name: r.name, reagent: r.canon, title: res.title, recipe: r }); }
  }
  return out;
}
root.EngineSolve = { solve, ways, key, neutral, splitInput };
})(typeof window !== 'undefined' ? window : globalThis);
