/* engine_solve.js: reagent text → stages → recipes; multi-step, work-ups, "no reaction" reasons, other ways.
   Also: several drawn molecules (multi-component recipes first, otherwise the biggest reacts and the rest
   join the reagents), equivalents ("2 eq", "excess", "2HBr") and a hint when more sites could react. */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore, RP = root.RxParse;
const WORK = new Set(['H3O+', 'H2O', 'NH4Cl']);
const isWork = st => st.sp.size > 0 && [...st.sp.keys()].every(k => WORK.has(k));
const METAL = new Set(['Mg', 'Li', 'Na', 'K', 'Al', 'B', 'Zn', 'Cu', 'Hg']);
const STOICH = new Set(['hx', 'hbr_roor', 'x2', 'hydration', 'radx', 'rx_oh', 'rx_nu', 'rx_elim', 'ozone', 'h2', 'eas', 'alc_hx', 'alc_pcl', 'acylation', 'amine_rx', 'ether_hi', 'hyd_red', 'epox', 'diol', 'oxymerc', 'hydroboration', 'grignard_add', 'imine', 'ester_hyd', 'alc_ox_mild', 'kmno4hot']);
function neutral(m) {
  const t = m.clone(); let ch = false;
  for (const b of [...t.bonds.values()]) { const A = t.atoms.get(b.a), B = t.atoms.get(b.b); if ((['O', 'N'].includes(A.el) && METAL.has(B.el)) || (['O', 'N'].includes(B.el) && METAL.has(A.el))) { t.removeBond(b.id); ch = true; } }
  for (const a of t.atoms.values()) if (a.chg < 0 && ['O', 'N', 'S'].includes(a.el)) { a.chg = 0; ch = true; }
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
  if (REV) return REV; REV = new Map(); const pre = { RX: v => (v.X || 'Cl') + v.smi, RMgX: v => (v.X || 'Br') + '[Mg]' + v.smi, ROH: v => 'O' + v.smi, RNH2: v => 'N' + v.smi, RCOCl: v => 'ClC(=O)' + v.smi, R2NH: v => v.smi + 'N' + v.smi };
  const put = (smi, v) => { try { const mol = C.RD.get_mol(smi); if (!mol) return; const k = mol.get_smiles(); mol.delete(); if (!REV.has(k)) REV.set(k, v); } catch (e) { } };
  for (const s of RP.SP || []) if (s.smi && !pre[s.id]) put(s.smi, { sp: s.id, smi: s.smi, txt: s.al[0] });
  for (const [a, v] of RP.ALIAS) { if (!v.sp) continue; if (pre[v.sp] && v.smi && !v.gen) put(pre[v.sp](v), Object.assign({}, v, { txt: a })); else if (v.smi && !pre[v.sp]) put(v.smi, Object.assign({}, v, { txt: a })); }
  return REV;
}
const SPEC = { 'Br': 'Br2', 'ClCl': 'Cl2', 'II': 'I2', 'BrBr': 'Br2', 'Cl': 'HCl', 'Br_': 'HBr', 'I': 'HI', 'O': 'H2O', 'N': 'NH3', '[HH]': 'H2', 'O=C=O': 'CO2', '[C-]#[O+]': 'CO', 'ClC(Cl)Cl': 'CHCl3', 'NN': 'NH2NH2', 'NO': 'NH2OH', 'C1CO1': 'ethyleneoxide', 'OCCO': 'glycol' };
function heavyCount(m, ids) { return ids.filter(id => m.atoms.get(id).el !== 'H').length; }
function cCount(m, ids) { return ids.filter(id => m.atoms.get(id).el === 'C').length; }
function rootedR(m, ids, att, drop) { const keep = ids.filter(id => !drop.includes(id)); const order = [att, ...keep.filter(id => id !== att)]; try { return C.toSmiles(m, { ids: order }).smi; } catch (e) { return null; } }
function classify(m, ids) { // -> { alias entry } to inject into the first stage, or null
  const sub = m.sub(ids); let k = null; try { k = C.canon(sub); } catch (e) { }
  if (k) { const rv = reverseAliases().get(k); if (rv) return rv; const sp = SPEC[k]; if (sp) return { sp }; }
  const els = ids.map(id => m.atoms.get(id).el);
  const hal = ids.filter(id => ['Cl', 'Br', 'I'].includes(m.atoms.get(id).el));
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
  fr.sort((p, q) => cCount(m, q) - cCount(m, p) || heavyCount(m, q) - heavyCount(m, p));
  const main = E.cleanProduct(m, fr[0]); const co = [], left = [];
  for (const ids of metals) { const e = m.atoms.get(ids[0]).el; co.push({ sp: MSP[e], txt: e }); }
  for (const ids of fr.slice(1)) { const v = classify(m, ids); if (v) co.push(v); else left.push(E.cleanProduct(m, ids)); }
  return { main, co, left };
}
function inject(st, v) { const prev = st.sp.get(v.sp); if (!prev) st.sp.set(v.sp, v); else if (!Array.isArray(prev)) st.sp.set(v.sp, [prev, v]); else prev.push(v); }

function keeps(a, b) { const ha = [...a.atoms.values()].filter(x => x.el !== 'H'); if (!ha.length) return false; const kept = ha.filter(x => b.atoms.has(x.id)).length; return kept / ha.length >= 0.6; }
function runStage(cur, ctx, all, k, cands) {
  for (const r of cands) { let o = null; try { o = r.run(cur, ctx, all, k); } catch (e) { if (root.console) console.warn('recipe ' + r.id + ' failed:', e && e.message); o = null; } if (o) return [o, r]; }
  return [null, null];
}
function solve(m, text, opt = {}) {
  const P = RP.parse(text || '');
  const stages = P.stages.filter(s => s.sp.size || s.flags.size || s.T != null);
  if (!stages.length) return { status: 'empty', unknown: P.unknown };
  let all = stages.map(RP.stageCtx);
  let cur = m.clone(); const steps = [], titles = [], why = [], tags = [], used = []; let extra = [], minor = []; const notes = []; let more = null, eqInfo = null;
  let multi = cur.frags().length > 1;
  for (let k = 0; k < all.length; k++) {
    let ctx = all[k];
    if (k > 0 && isWork(ctx.st)) { const n = neutral(cur); if (n) { steps.push({ L: cur.clone(), arrows: [], R: n.clone(), text: 'Work-up: the anion picks up H⁺ (and the metal salt washes away).', cond: ['H₃O⁺'] }); cur = n; } continue; }
    if (!ctx.st.sp.size && !k && !ctx.st.flags.size && ctx.T == null) return { status: 'unknown', unknown: P.unknown, stage: ctx.raw };
    let o = null, rec = null, cands = [];
    if (k === 0 && multi) {
      const mc = E.REC.filter(r => r.multi && (() => { try { return r.test(ctx, all, k); } catch (e) { return false; } })());
      [o, rec] = runStage(cur, ctx, all, k, mc);
      if (!o) { const sp = splitInput(cur); if (sp) { cur = sp.main; for (const v of sp.co) inject(ctx.st, v); all[0] = ctx = RP.stageCtx(ctx.st); if (sp.co.length) notes.push(`The other drawn molecule${sp.co.length > 1 ? 's were' : ' was'} read as reagent${sp.co.length > 1 ? 's' : ''} (${sp.co.map(v => v.txt === 'drawn' ? v.sp : v.txt).join(', ')}).`); if (sp.left.length) notes.push('One drawn molecule did not take part.'); } }
      multi = false;
    }
    if (!o) { cands = E.REC.filter(r => { try { return r.test(ctx, all, k); } catch (e) { return false; } }); [o, rec] = runStage(cur, ctx, all, k, cands); }
    if (!o) {
      const sp = recognised(ctx.st);
      if (P.unknown.length) return { status: 'unknown', unknown: P.unknown, known: sp, stage: ctx.raw, reason: `I can't read "${P.unknown.join('", "')}".` };
      return { status: 'none', reason: cands.length ? `${ctx.raw.trim()} has no reaction with this structure in the library or the JEE syllabus.` : `No reaction of this structure with ${ctx.raw.trim()} that I know of.`, stage: ctx.raw, steps, titles, partial: k ? cur : null, known: sp };
    }
    if (o.none || !o.product) return { status: 'none', reason: o.none || 'No reaction.', title: o.title, steps, titles, test: (o.tags || []).includes('test'), positive: !!o.positive };
    let title = o.title; const reps = [o];
    // equivalents: repeat the same stage on the product
    const eq = ctx.eq; if (!cands.length) cands = E.REC.filter(r => { try { return r.test(ctx, all, k); } catch (e) { return false; } });
    const lim = eq == null ? 1 : eq === Infinity ? 6 : Math.max(1, Math.min(6, Math.round(eq)));
    let p = o.product, pk = key(p), n = o.eqUsed || 1;
    while (n < lim) { const [o2, r2] = runStage(p, ctx, all, k, cands); if (!o2 || o2.none || !o2.product || !r2 || !(r2.stoich || STOICH.has(r2.id) || r2 === rec) || !keeps(p, o2.product)) break; const k2 = key(o2.product); if (!k2 || k2 === pk) break; reps.push(o2); p = o2.product; pk = k2; n++; }
    if (eq != null) eqInfo = { stage: k, asked: eq, used: n };
    else if (!opt.noMore && (rec.stoich || STOICH.has(rec.id))) { const [o2, r2] = runStage(p, ctx, all, k, cands); if (o2 && !o2.none && o2.product && key(o2.product) !== pk && r2 && (r2.stoich || STOICH.has(r2.id)) && keeps(p, o2.product)) more = { stage: k, title: o2.title }; }
    for (const x of reps) { steps.push(...x.steps); why.push(...(x.why || []).filter(Boolean)); tags.push(...(x.tags || [])); }
    if (reps.length > 1) { const ts = [...new Set(reps.map(x => x.title))]; title = ts.length === 1 ? `${ts[0]} (×${reps.length})` : ts.join(', then '); }
    if (n > 1) why.push(`${n} equivalents reacted`);
    if (eq != null && eq !== Infinity && n < lim) why.push(`Only ${n} equivalent${n > 1 ? 's' : ''} can react here; the rest stays unused`);
    titles.push(title); used.push(rec); const last = reps[reps.length - 1]; extra = last.extra || []; minor = (reps.length > 1 ? last.minor : o.minor) || [];
    cur = p;
    while (rec.eats && k + 1 < all.length && (() => { try { return rec.eats(all[k + 1]); } catch (e) { return false; } })()) k++;
  }
  /* protecting group: a ketal made in an earlier stage comes off in a later aqueous-acid work-up */
  { const pi = used.findIndex(r => r && (r.id === 'acetal' || r.id === 'acetal_roh'));
    if (pi >= 0 && pi < used.length - 1 && all.slice(1).some(x => x.any(['H3O+', 'H2O']))) { const hy = E.REC.find(r => r.id === 'acetal_hyd');
      let o = null; try { o = hy && hy.run(cur, all[all.length - 1], all, all.length - 1); } catch (e) { o = null; }
      if (o && o.product && !o.none) { steps.push(...o.steps.map(s => Object.assign({}, s, { text: s.text }))); if (o.steps.length) steps[steps.length - 1] = Object.assign({}, steps[steps.length - 1], { text: (steps[steps.length - 1].text || '') + ' The aqueous acid work-up also removes the ketal protecting group: the C=O is back.' }); cur = o.product; titles.push('deprotection'); why.push('The ketal protected the C=O during the step that would otherwise have attacked it'); } } }
  let ok = true; try { ok = C.isValid(cur); } catch (e) { ok = false; }
  return { status: ok ? 'ok' : 'error', product: cur, extra, minor, steps, title: titles.join(', then '), why: [...new Set(why.concat(notes))], tags, recipes: used, unknown: P.unknown, more, eq: eqInfo };
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
