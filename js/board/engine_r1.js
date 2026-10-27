/* engine_r1.js: recipes for alkenes, alkynes, radical halogenation, alkyl halides, alcohols and ethers */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore;
const { BL, el, P, heavyN, hExp, hTot, cN, bo, aromSet, Mech, cationScore, radicalScore, alkeneScore, freeDir, moveSet, rotSet, centroid, unit, rot } = E;
const REC = E.REC = E.REC || [];
const reg = r => REC.push(r);
const HX = { HCl: 'Cl', HBr: 'Br', HI: 'I' };
const SUBS = s => s.replace(/(\d)/g, d => '₀₁₂₃₄₅₆₇₈₉'[+d]);
/* ---------- small helpers ---------- */
function M(m, sm, o) { return C.matchAll(m, sm, o || {}); }
function uniq(list, keyf) { const s = new Set(), o = []; for (const x of list) { const k = keyf(x); if (!s.has(k)) { s.add(k); o.push(x); } } return o; }
function alkenes(m) { const AR = aromSet(m); return uniq(M(m, '[C:1]=[C:2]').filter(([a, b]) => !AR.has(a) && !AR.has(b)), p => [...p].sort().join()); }
function alkynes(m) { return uniq(M(m, '[C:1]#[C:2]'), p => [...p].sort().join()); }
const subtree = E.subtree = function (m, root, block) { const seen = new Set([root]), st = [root]; while (st.length) { const x = st.pop(); for (const n of m.nbrs(x)) if (n !== block && !seen.has(n)) { seen.add(n); st.push(n); } } return [...seen]; };
function moveGroupTo(m, g, from, to) { // rigidly move group rooted at g (detached from `from`) so g sits at to's free direction
  const ids = subtree(m, g, to); const d = freeDir(m, to, new Set(ids)); const tp = P(m, to);
  const len = el(m, g) === 'H' ? 0.82 * BL : BL; const target = [tp[0] + d[0] * len, tp[1] + d[1] * len];
  if (ids.length > 1) { const c = centroid(m, ids), gp = P(m, g); const v = [c[0] - gp[0], c[1] - gp[1]]; if (Math.hypot(v[0], v[1]) > 1e-3) rotSet(m, ids, gp[0], gp[1], Math.atan2(d[1], d[0]) - Math.atan2(v[1], v[0])); }
  const gp = P(m, g); moveSet(m, ids, target[0] - gp[0], target[1] - gp[1]);
}
E.moveGroupTo = moveGroupTo;
function setAlkeneGeom(m, a, b, cis) { // after alkyne -> alkene: zig-zag with cis or trans substituents
  const ra = heavyN(m, a).filter(x => x !== b)[0], rb = heavyN(m, b).filter(x => x !== a)[0];
  const pa = P(m, a), pb = P(m, b); const u = unit([pb[0] - pa[0], pb[1] - pa[1]]);
  const dA = rot(u, 2 * Math.PI / 3), dB = rot([-u[0], -u[1]], cis ? -2 * Math.PI / 3 : 2 * Math.PI / 3);
  const put = (r, at, d) => { if (r == null) return; const ids = subtree(m, r, at); const rp = P(m, r), ap = P(m, at); const cur = [rp[0] - ap[0], rp[1] - ap[1]]; rotSet(m, ids, ap[0], ap[1], Math.atan2(d[1], d[0]) - Math.atan2(cur[1], cur[0])); const rp2 = P(m, r); moveSet(m, ids, ap[0] + d[0] * BL - rp2[0], ap[1] + d[1] * BL - rp2[1]); };
  put(ra, a, dA); put(rb, b, dB);
}
E.setAlkeneGeom = setAlkeneGeom;
function dirTo(m, from, to) { const p = P(m, from), q = P(m, to); return unit([q[0] - p[0], q[1] - p[1]]); }
function opp(d) { return [-d[0], -d[1]]; }
function moveFragNear(mech, fragAtom, site, dist, dir) { const m = mech.m; const ids = E.fragIds(m, fragAtom); const d = dir || freeDir(m, site, new Set(ids)); const s = P(m, site), a = P(m, fragAtom); const tx = s[0] + d[0] * (dist || 1.7 * BL), ty = s[1] + d[1] * (dist || 1.7 * BL);
  if (ids.length > 1) { const c = centroid(m, ids); const v = [c[0] - a[0], c[1] - a[1]]; if (Math.hypot(v[0], v[1]) > 1e-3) rotSet(m, ids, a[0], a[1], Math.atan2(d[1], d[0]) - Math.atan2(v[1], v[0])); }
  const a2 = P(m, fragAtom); moveSet(m, ids, tx - a2[0], ty - a2[1]); }
E.moveFragNear = moveFragNear;
function deg(m, c) { const n = cN(m, c).length; return n >= 3 ? '3°' : n === 2 ? '2°' : n === 1 ? '1°' : 'methyl'; }
E.deg = deg;
/* cation rearrangement: returns {type, from, group} or null */
function bestShift(m, cp) {
  const AR = aromSet(m); const base = cationScore(m, cp, AR); let best = null;
  for (const n of cN(m, cp)) {
    if (AR.has(n) || m.bondsOf(n).some(b => b.o > 1)) continue;
    const cands = [];
    if (hTot(m, n) > 0) cands.push({ type: 'H', g: null });
    for (const g of heavyN(m, n)) if (g !== cp && (el(m, g) === 'C')) cands.push({ type: AR.has(g) ? 'aryl' : 'alkyl', g });
    for (const cd of cands) {
      const t = m.clone(); let gid = cd.g;
      if (cd.type === 'H') gid = E.expH(t, n);
      t.removeBond(t.bondBetween(n, gid).id); t.addBond(cp, gid, 1); t.atoms.get(cp).chg = 0; t.atoms.get(n).chg = 1;
      let sc = cationScore(t, n, aromSet(t));
      // ring expansion / strain relief
      const r0 = E.ringSizeOf(m, n, gid), r1 = E.ringSizeOf(t, cp, gid);
      if (cd.type !== 'H' && r0 && r1 && r1 > r0 && r0 <= 4) sc += r0 === 3 ? 3 : 2.2;
      const gain = sc - base + (cd.type === 'H' ? 0.15 : cd.type === 'aryl' ? 0.1 : 0);
      if (gain > 0.8 && (!best || gain > best.gain)) best = { type: cd.type, from: n, g: cd.g, gain };
    }
  }
  return best;
}
function doShift(mech, cp, sh) {
  const m = mech.m; let g = sh.g; if (sh.type === 'H') g = mech.H(sh.from);
  const lbl = sh.type === 'H' ? '1,2-hydride shift' : sh.type === 'aryl' ? '1,2-aryl shift' : (E.ringSizeOf(m, sh.from, g) ? 'Ring expansion' : '1,2-methyl shift');
  const fromDeg = deg(m, cp);
  mech.step([[[sh.from, g], [cp, g]]], `${lbl}: the positive charge moves to the more stable carbon.`, [], () => { if (!E.ringSizeOf(mech.m, cp, g)) moveGroupTo(mech.m, g, sh.from, cp); else E.tidyLongBonds(mech.m, 1.45); });
  return sh.from;
}
E.bestShift = bestShift; E.doShift = doShift;
/* carbocation from protonation of C=C: returns cation atom id; species = '[H]X' etc. */
function protonateAlkene(mech, a, b, spec, why) {
  const m = mech.m; const AR = aromSet(m);
  const sa = cationScore(m, a, AR), sb = cationScore(m, b, AR);
  const [hC, cp] = sa > sb ? [b, a] : [a, b];
  const mid = [(P(m, a)[0] + P(m, b)[0]) / 2, (P(m, a)[1] + P(m, b)[1]) / 2];
  const ids = mech.add(spec, hC, { dist: 1.65 * BL });
  const H = ids.get(1), X = ids.get(2);
  const t = Math.abs(sa - sb) < 0.5 ? 'The π bond grabs H⁺.' : `The π bond grabs H⁺. H adds to the carbon that leaves the ${deg(m, cp)} carbocation (more stable): Markovnikov.`;
  mech.step([[[a, b], [hC, H]], [[H, X], [X]]], why || t);
  return { cp, hC, X, tie: Math.abs(sa - sb) < 0.5 };
}
function captureBy(mech, cp, nu, text) { moveFragNear(mech, nu, cp, 1.65 * BL); mech.step([[[nu], [nu, cp]]], text); }
function deprotOn(mech, at, text, base, to) { // remove an H from atom `at` (O+ or N+) with water / base; `to` = atom that takes the bond pair (forms a π bond)
  const m = mech.m; const h = mech.H(at); if (h == null) return;
  const ids = mech.add(base || 'O', h, { dist: 1.5 * BL, dir: dirTo(m, at, h) }); const w = ids.get(1);
  mech.step([[[w], [w, h]], [[h, at], to != null ? [at, to] : [at]]], text || 'Water removes the proton.');
}
E.protonateAlkene = protonateAlkene; E.captureBy = captureBy; E.deprotOn = deprotOn;
function out(o) { return Object.assign({ steps: [], minor: [], extra: [], why: [], db: [] }, o); }
E.out = out;
function done(mech, o) { return out(Object.assign({ steps: mech.steps, product: mech.product() }, o)); }
E.done = done;
function prodWithout(m, fn) { const t = m.clone(); fn(t); return E.cleanProduct(t, t.frags().sort((p, q) => q.length - p.length)[0]); }
E.prodWithout = prodWithout;
function addOn(m, c, elx, o, chg) { const d = freeDir(m, c); const p = P(m, c); return m.addAtom(elx, p[0] + d[0] * BL, p[1] + d[1] * BL, { chg: chg || 0 }) && (() => { const id = m.nid - 1; m.addBond(c, id, o || 1); return id; })(); }
E.addOn = addOn;
function replaceAtom(m, id, elx) { m.atoms.get(id).el = elx; }

/* ===================== ALKENES ===================== */
reg({ id: 'hx', name: 'HX addition', fam: 'Alkene', canon: 'HBr', db: ['Markovnikov', 'HBr'],
  test: c => c.any(['HCl', 'HBr', 'HI']) && !(c.has('ROOR') && c.has('HBr')) && !c.light,
  run(m, c) {
    const X = c.has('HBr') ? 'Br' : c.has('HCl') ? 'Cl' : 'I';
    const al = alkenes(m), ay = alkynes(m); if (!al.length && !ay.length) return null;
    if (!al.length) { // alkyne: Markovnikov, 1 or 2 equivalents
      const [a, b] = ay[0]; const mech = new Mech(m);
      const r = protonateAlkene(mech, a, b, `[H]${X}`); captureBy(mech, r.cp, r.X, `${X}⁻ captures the vinyl cation.`);
      const two = c.flag('excess') || /2\s*(eq|mol|hx|hbr|hcl)|excess/i.test(c.raw);
      if (two) { const r2 = protonateAlkene(mech, a, b, `[H]${X}`, `Second HX adds the same way: the halogen already on the carbon stabilises the cation (lone-pair donation), so both X end up on one carbon.`); captureBy(mech, r2.cp, r2.X, `${X}⁻ captures it: a gem-dihalide.`); }
      return done(mech, { eqUsed: two ? 2 : 1, title: two ? 'HX × 2 on alkyne → gem-dihalide' : 'HX on alkyne → vinyl halide', why: ['Markovnikov: H on the carbon with more H'], tags: ['AdE'] });
    }
    // conjugated diene? 1,2 vs 1,4
    const [a, b] = al[0]; const mech = new Mech(m);
    const r = protonateAlkene(mech, a, b, `[H]${X}`);
    let cp = r.cp; const why = [];
    const sh = bestShift(mech.m, cp); if (sh) { cp = doShift(mech, cp, sh); why.push('Carbocation rearranges first (more stable cation)'); }
    captureBy(mech, cp, r.X, `${X}⁻ attacks the carbocation.`);
    const o = done(mech, { title: `Markovnikov addition of H${X}`, why: ['H⁺ adds first, to the carbon that gives the more stable carbocation', ...why], tags: ['AdE', 'Markovnikov'] });
    if (r.tie && !sh) o.why.push('Both carbons give similar cations: a mixture is possible');
    return o;
  } });
reg({ id: 'hbr_roor', name: 'HBr, peroxide (anti-Markovnikov)', fam: 'Alkene', canon: 'HBr, ROOR', db: ['Kharasch'],
  test: c => c.has('HBr') && (c.has('ROOR') || c.light),
  run(m, c) {
    const al = alkenes(m); if (!al.length) return null;
    const [a, b] = al[0]; const AR = aromSet(m); const mech = new Mech(m);
    const ra = radicalScore(m, a, AR), rb = radicalScore(m, b, AR); const [brC, rc] = ra > rb ? [b, a] : [a, b];
    // initiation shown compactly
    let ids = mech.add({ atoms: [{ el: 'O', rad: 1, x: 0 }, { el: 'C', x: 1 }], bonds: [[0, 1]] }, brC, { dist: 4.2 * BL });
    const ro = ids.get(1); ids = mech.add('[H]Br', ro, { dist: 1.2 * BL }); const h = ids.get(1), br = ids.get(2);
    mech.step([[[ro], [ro, h], 1], [[h, br], [ro, h], 1], [[h, br], [br], 1]], 'Initiation: an RO• radical (from the peroxide) pulls H• off HBr, leaving Br•.', ['ROOR']);
    mech.dropFrag(ro);
    moveFragNear(mech, br, brC, 1.6 * BL);
    mech.step([[[br], [br, brC], 1], [[a, b], [br, brC], 1], [[a, b], [rc], 1]], `Br• adds to the CH end so the radical sits on the more substituted carbon (${deg(mech.m, rc)} radical, more stable).`);
    ids = mech.add('[H]Br', rc, { dist: 1.6 * BL }); const h2 = ids.get(1), br2 = ids.get(2);
    mech.step([[[rc], [rc, h2], 1], [[h2, br2], [rc, h2], 1], [[h2, br2], [br2], 1]], 'The carbon radical takes H• from another HBr: product + new Br• (chain continues).');
    return done(mech, { title: 'Anti-Markovnikov addition of HBr (peroxide effect)', why: ['Radical chain: Br• adds first, to the less substituted carbon', 'Only HBr shows the peroxide effect (not HCl or HI)'], tags: ['AdR', 'anti-Markovnikov'] });
  } });
reg({ id: 'hydration', name: 'Acid-catalysed hydration', fam: 'Alkene', canon: 'dil. H₂SO₄', db: ['Acid-catalysed hydration of alkene', 'Kucherov'],
  test: (c, all) => (c.has('H3O+') || c.has('H2SO4') || c.has('H3PO4')) && !c.heat && !c.flag('conc') && !c.any(['HNO3', 'NaNO2', 'HNO2', 'KMnO4', 'K2Cr2O7']) || (c.has('Hg2+') && (c.has('H2SO4') || c.has('H3O+') || c.has('H2O')) && !(all || [c]).some(s => s.has('NaBH4'))),
  run(m, c) {
    const al = alkenes(m), ay = alkynes(m);
    if (!ay.length && c.has('Hg2+') && !c.any(['H2SO4', 'H3O+'])) return null;
    if (ay.length && (c.has('Hg2+') || !al.length)) {
      const [a, b] = ay[0]; const mech = new Mech(m); const AR = aromSet(m);
      const t = mech.m.clone(); // Markovnikov enol
      const sa = cationScore(t, a, AR), sb = cationScore(t, b, AR); const [hC, oC] = sa > sb ? [b, a] : [a, b];
      const enol = mech.m.clone(); enol.bondBetween(a, b).o = 2; const o = addOn(enol, oC, 'O');
      mech.note('Hg²⁺ / H⁺ add water across the triple bond (Markovnikov): the OH goes on the more substituted carbon, giving an enol.', mech.m.clone(), enol.clone());
      mech.m = enol; mech.sub.add(o);
      // keto-enol tautomerism (acid)
      let ids = mech.add('[H][OH2+]', hC, { dist: 1.65 * BL }); const h = ids.get(1), w = ids.get(2);
      mech.step([[[o], [o, oC]], [[oC, hC], [hC, h]], [[h, w], [w]]], 'Tautomerism 1: the O lone pair pushes, the enol C=C takes H⁺ on carbon (protonated ketone).');
      mech.dropFrag(w);
      deprotOn(mech, o, 'Tautomerism 2: water removes H⁺ from O, giving the C=O. The keto form is far more stable.');
      return done(mech, { title: 'Hydration of alkyne → ketone (Kucherov)', why: ['Markovnikov enol, then keto-enol tautomerism', 'Terminal alkynes give methyl ketones (ethyne gives ethanal)'], tags: ['AdE', 'tautomerism'] });
    }
    if (!al.length) return null;
    const [a, b] = al[0]; const mech = new Mech(m);
    const r = protonateAlkene(mech, a, b, '[H][OH2+]');
    let cp = r.cp; const why = [];
    const sh = bestShift(mech.m, cp); if (sh) { cp = doShift(mech, cp, sh); why.push('Carbocation rearranges (1,2-shift)'); }
    captureBy(mech, cp, r.X, 'Water attacks the carbocation.');
    deprotOn(mech, r.X, 'A second water removes H⁺ from the oxonium ion.');
    return done(mech, { title: 'Acid-catalysed hydration (Markovnikov)', why: ['Via carbocation: Markovnikov, rearrangement possible', ...why], tags: ['AdE'] });
  } });
reg({ id: 'oxymerc', name: 'Oxymercuration–demercuration', fam: 'Alkene', canon: 'Hg(OAc)₂, H₂O then NaBH₄', db: ['Oxymercuration–demercuration'], eats: s => s.any(['NaBH4', 'H2O', 'H3O+']),
  test: c => c.has('Hg2+') && (c.has('NaBH4') || !c.has('H2SO4')),
  run(m, c) {
    const al = alkenes(m); if (!al.length) return null; const [a, b] = al[0]; const AR = aromSet(m);
    const mech = new Mech(m); const sa = cationScore(m, a, AR), sb = cationScore(m, b, AR); const [more, less] = sa >= sb ? [a, b] : [b, a];
    const mid = [(P(m, a)[0] + P(m, b)[0]) / 2, (P(m, a)[1] + P(m, b)[1]) / 2];
    const ids = mech.add({ atoms: [{ el: 'Hg', chg: 1 }, { el: 'O', x: 1 }, { el: 'C', x: 1.6, y: -0.8 }, { el: 'O', x: 1.3, y: -1.7 }, { el: 'C', x: 2.6, y: -0.6 }], bonds: [[0, 1], [1, 2], [2, 3, 2], [2, 4]] }, less, { dist: 1.5 * BL });
    const hg = ids.get(1);
    mech.step([[[a, b], [more, hg]], [[hg], [hg, less]]], 'The π bond attacks Hg(OAc)⁺: a bridged mercurinium ion forms (no free carbocation, so no rearrangement).');
    const w = mech.add('O', more, { dist: 1.6 * BL }).get(1);
    mech.step([[[w], [w, more]], [[more, hg], [hg]]], 'Water opens the bridge at the more substituted carbon (it carries more δ+).');
    deprotOn(mech, w, 'Loss of H⁺ gives the organomercury alcohol.');
    const L = mech.m.clone(); const hgIds = E.fragIds(L, hg).filter(id => id !== less && !mech.sub.has(id));
    // demercuration: Hg(OAc) replaced by H
    const R = L.clone(); for (const id of E.subtree(R, hg, less)) R.removeAtom(id);
    mech.m = R; mech.note('NaBH₄ replaces HgOAc by H (demercuration).', L, R.clone());
    return done(mech, { title: 'Oxymercuration–demercuration (Markovnikov, no rearrangement)', why: ['Markovnikov alcohol', 'No carbocation, so no rearrangement'], tags: ['AdE'] });
  } });
function alkyneOnly() { return false; }
reg({ id: 'hydroboration', name: 'Hydroboration–oxidation', fam: 'Alkene', canon: 'BH₃·THF then H₂O₂/OH⁻', db: ['Hydroboration–oxidation'], eats: s => s.any(['H2O2', 'H2O2/OH', 'KOH', 'aqKOH', 'H2O']),
  test: c => c.has('BH3'),
  run(m, c, all) {
    const ox = all && all.some(s => s.has('H2O2') || s.has('H2O2/OH')) || c.has('H2O2') || c.has('H2O2/OH');
    const al = alkenes(m), ay = alkynes(m); if (!al.length && !ay.length) return null;
    const isYne = !al.length; const [a, b] = isYne ? ay[0] : al[0]; const AR = aromSet(m);
    // B goes to the less substituted (less hindered) carbon
    const ha = hTot(m, a), hb = hTot(m, b); const [bC, hC] = ha >= hb ? [a, b] : [b, a];
    const mech = new Mech(m);
    const ids = mech.add({ atoms: [{ el: 'B' }, { el: 'H', x: 0.15, y: 0.82 }], bonds: [[0, 1]] }, bC, { dist: 1.55 * BL }); const B = ids.get(1), H = ids.get(2);
    moveFragNear(mech, B, bC, 1.55 * BL, freeDir(mech.m, bC, new Set([B, H])));
    // put the B-H so that H is nearer the other carbon
    const pb = P(mech.m, B), ph = P(mech.m, hC); const hh = mech.m.atoms.get(H); hh.x = pb[0] + (ph[0] - pb[0]) * 0.45; hh.y = pb[1] + (ph[1] - pb[1]) * 0.45;
    mech.step([[[a, b], [bC, B]], [[B, H], [H, hC]]], 'One concerted step (4-centre TS): B bonds to the less hindered carbon while its H moves to the other carbon. Syn addition.');
    if (!ox) return done(mech, { title: 'Hydroboration (alkylborane)', why: ['Add H₂O₂ / OH⁻ next to get the alcohol'], tags: ['concerted'] });
    const oo = mech.add('[O-]O', B, { dist: 1.6 * BL }); const o1 = oo.get(1), o2 = oo.get(2);
    mech.step([[[o1], [o1, B]]], 'HOO⁻ adds to the empty p orbital of boron.');
    mech.step([[[bC, B], [bC, o1]], [[o1, o2], [o2]]], 'The alkyl group shifts from B to O with its electrons (retention), pushing out OH⁻.', [], () => {
      const m2 = mech.m; const pb2 = P(m2, bC); const d = unit([P(m2, o1)[0] - pb2[0], P(m2, o1)[1] - pb2[1]]); const o = m2.atoms.get(o1); o.x = pb2[0] + d[0] * BL; o.y = pb2[1] + d[1] * BL; const bb = m2.atoms.get(B); bb.x = o.x + d[0] * BL; bb.y = o.y + d[1] * BL; E.moveFragNear(mech, o2, B, 2.2 * BL); });
    mech.dropFrag(o2);
    const oh = mech.add('[OH-]', B, { dist: 1.5 * BL }).get(1);
    mech.step([[[oh], [oh, B]], [[B, o1], [o1]]], 'OH⁻ attacks boron and the B–O bond breaks: an alkoxide leaves.');
    mech.dropFrag(B);
    const w = mech.add('O', o1, { dist: 1.5 * BL }).get(1); const hw = E.expH(mech.m, w);
    mech.step([[[o1], [o1, hw]], [[hw, w], [w]]], 'Water protonates the alkoxide.');
    if (isYne) {
      // enol -> carbonyl
      deprotOn(mech, o1, 'The product is an enol. Tautomerism 1: OH⁻ removes the O–H proton (enolate).', '[OH-]');
      mech.dropFrag([...mech.m.atoms.keys()].find(id => !mech.sub.has(id) && el(mech.m, id) === 'O' && mech.m.heavyDegree(id) === 0 && id !== o1) || -1);
      const ids2 = mech.add('[H]O', hC, { dist: 1.6 * BL }); const h3 = ids2.get(1), w3 = ids2.get(2);
      mech.step([[[o1], [o1, bC]], [[bC, hC], [hC, h3]], [[h3, w3], [w3]]], 'Tautomerism 2: O⁻ pushes down to C=O while the C=C takes H⁺ from water on carbon. Terminal alkyne → aldehyde.'); mech.dropFrag(w3);
    }
    return done(mech, { title: isYne ? 'Hydroboration–oxidation of alkyne → carbonyl' : 'Hydroboration–oxidation (anti-Markovnikov, syn)', why: ['OH ends on the less substituted carbon', 'No carbocation: no rearrangement', 'Syn addition of H and OH'], tags: ['concerted', 'anti-Markovnikov'] });
  } });
reg({ id: 'x2', name: 'Halogen addition', fam: 'Alkene', canon: 'Br₂ in CCl₄', db: ['Bromination of alkene (halonium, anti)', 'Halohydrin formation'],
  test: c => c.any(['Br2', 'Cl2', 'I2', 'HOBr', 'HOCl']) && !c.light && !(c.T != null && c.T >= 450) && !c.has('FeX3') && !c.any(['KOH', 'aqKOH', 'alcKOH', 'Br2/NaOH', 'NaOCl', 'redP', 'HVZ']),
  run(m, c) {
    const al = alkenes(m), ay = alkynes(m); if (!al.length && !ay.length) return null;
    const X = c.any(['Br2', 'HOBr']) ? 'Br' : c.any(['Cl2', 'HOCl']) ? 'Cl' : 'I';
    const water = c.any(['HOBr', 'HOCl']) || c.has('H2O') || c.flag('aq');
    const [a, b] = al.length ? al[0] : ay[0]; const AR = aromSet(m);
    const mech = new Mech(m); const sa = cationScore(m, a, AR), sb = cationScore(m, b, AR); const [more, less] = sa >= sb ? [a, b] : [b, a];
    const ids = mech.add(`${X}${X}`, less, { dist: 1.55 * BL }); const x1 = ids.get(1), x2 = ids.get(2);
    mech.step([[[a, b], [less, x1]], [[x1], [x1, more]], [[x1, x2], [x2]]], `The π bond attacks ${X}₂; a bridged ${X === 'Br' ? 'bromonium' : X === 'Cl' ? 'chloronium' : 'iodonium'} ion forms and ${X}⁻ leaves.`);
    if (water) {
      mech.dropFrag(x2); const w = mech.add('O', more, { dist: 1.6 * BL, dir: unit([P(mech.m, more)[0] - P(mech.m, x1)[0], P(mech.m, more)[1] - P(mech.m, x1)[1]]) }).get(1);
      mech.step([[[w], [w, more]], [[more, x1], [x1]]], 'Water (solvent, in excess) opens the bridge from the back at the more substituted carbon.');
      deprotOn(mech, w, 'Loss of H⁺ gives the halohydrin.');
      return done(mech, { title: `Halohydrin formation (${X}₂/H₂O)`, why: ['OH on the more substituted carbon, X on the other', 'Anti addition'], tags: ['AdE', 'anti'] });
    }
    moveFragNear(mech, x2, more, 1.65 * BL, unit([P(mech.m, more)[0] - P(mech.m, x1)[0], P(mech.m, more)[1] - P(mech.m, x1)[1]]));
    mech.step([[[x2], [x2, more]], [[more, x1], [x1]]], `${X}⁻ opens the bridge from the opposite face: anti addition.`);
    if (!al.length && (c.flag('excess') || /2\s*eq|excess/i.test(c.raw))) {
      const t = mech.product(); const r2 = E.REC.find(r => r.id === 'x2').run(t, c); if (r2) { r2.steps = mech.steps.concat(r2.steps); r2.title = 'Excess X₂ on alkyne → tetrahalide'; return r2; }
    }
    return done(mech, { title: al.length ? `Anti addition of ${X}₂ (vicinal dihalide)` : `${X}₂ on alkyne (1 eq) → dihaloalkene`, why: ['Bridged halonium ion, then back-side attack: anti', 'Br₂ decolourises: test for unsaturation'], tags: ['AdE', 'anti'] });
  } });
reg({ id: 'h2', name: 'Catalytic hydrogenation', fam: 'Alkene', canon: 'H₂, Pd/C', db: ['Catalytic hydrogenation of alkenes', 'Lindlar reduction', 'Na/liq. NH3 reduction of internal alkyne'],
  test: c => (c.has('H2') && (c.has('Pd') || c.has('Lindlar'))) || c.has('Lindlar') || c.has('Na/NH3'),
  run(m, c) {
    const al = alkenes(m), ay = alkynes(m); const nitro = M(m, '[N+](=O)[O-]');
    if (c.has('Na/NH3') || c.has('Lindlar')) {
      if (!ay.length) return null; const [a, b] = ay[0]; const cis = c.has('Lindlar');
      const p = E.transform(m, t => { t.bondBetween(a, b).o = 2; setAlkeneGeom(t, a, b, cis); });
      const s = { L: m.clone(), arrows: [], R: p.clone(), text: cis ? 'H₂ adds from the catalyst surface on one face: syn addition, cis-alkene. The poisoned catalyst stops at the alkene.' : 'Na gives electrons one at a time (radical anion → vinyl radical → anion); the trans arrangement keeps the groups apart, so the trans-alkene forms.', cond: [] };
      return out({ title: cis ? 'Lindlar: alkyne → cis-alkene' : 'Na/liq. NH₃: alkyne → trans-alkene', product: E.cleanProduct(p, [...p.atoms.keys()]), steps: [s], why: [cis ? 'Syn addition on the catalyst surface' : 'Dissolving-metal reduction: anti, trans product'], tags: [cis ? 'syn' : 'anti'] });
    }
    if (!al.length && !ay.length) return null;
    const p = E.transform(m, t => { const A = aromSet(t); for (const bd of [...t.bonds.values()]) { if (bd.o >= 2 && el(t, bd.a) === 'C' && el(t, bd.b) === 'C' && !A.has(bd.a) && !A.has(bd.b)) bd.o = 1; } });
    return out({ title: 'Catalytic hydrogenation → alkane', product: p, steps: [{ L: m.clone(), arrows: [], R: p.clone(), text: 'H₂ adds across every C=C / C≡C on the metal surface (syn). Aromatic rings and C=O are left alone under these conditions.', cond: [] }], why: ['Syn addition', 'Benzene ring needs harsh conditions (high T, P)'], tags: ['surface'] });
  } });
reg({ id: 'ozone', name: 'Ozonolysis', fam: 'Alkene', canon: 'O₃ then Zn/H₂O', db: ['Ozonolysis'], eats: s => s.any(['Me2S', 'H2O2', 'Zn', 'H2O', 'H3O+']),
  test: c => c.has('O3'),
  run(m, c, all) {
    const al = alkenes(m), ay = alkynes(m); if (!al.length && !ay.length) return null;
    const oxw = (all || []).some(s => s.has('H2O2')) || c.has('H2O2');
    const p = E.transform(m, t => {
      for (const [a, b] of al.concat(ay)) {
        const pa = P(t, a), pb = P(t, b); t.removeBond(t.bondBetween(a, b).id);
        for (const [x, y] of [[a, pb], [b, pa]]) { const px = P(t, x); const d = unit([y[0] - px[0], y[1] - px[1]]); const o = t.addAtom('O', px[0] + d[0] * BL * 0.95, px[1] + d[1] * BL * 0.95); t.addBond(x, o, 2); if (oxw && t.hImp(x) > 0) { const d2 = freeDir(t, x, new Set([o])); const pp = P(t, x); const oh = t.addAtom('O', pp[0] + d2[0] * BL, pp[1] + d2[1] * BL); t.addBond(x, oh, 1); } }
        // pull the two halves apart
        const fa = E.fragIds(t, a), fb = E.fragIds(t, b); if (!fa.includes(b)) { const v = unit([pa[0] - pb[0], pa[1] - pb[1]]); moveSet(t, fa, v[0] * BL * 0.6, v[1] * BL * 0.6); moveSet(t, fb, -v[0] * BL * 0.6, -v[1] * BL * 0.6); }
      }
    });
    const fr = p.frags().map(ids => E.cleanProduct(p, ids)).sort((x, y) => y.atoms.size - x.atoms.size);
    return out({ title: oxw ? 'Ozonolysis, oxidative work-up (H₂O₂)' : 'Ozonolysis, reductive work-up (Zn/H₂O)', product: fr[0], extra: fr.slice(1), steps: [{ L: m.clone(), arrows: [], R: p.clone(), text: 'O₃ adds across C=C (molozonide → ozonide). Work-up splits it: every C=C becomes two C=O. ' + (oxw ? 'H₂O₂ then oxidises the aldehydes to acids.' : 'Zn/H₂O (or Me₂S) stops at aldehydes/ketones.'), cond: [] }], why: ['Cut each C=C, put =O on both carbons', 'Count the carbons to trace back the alkene'], tags: ['cleavage'] });
  } });
reg({ id: 'diol', name: 'Syn dihydroxylation', fam: 'Alkene', canon: 'cold dil. KMnO₄', db: ["Baeyer's reagent (syn dihydroxylation)", 'OsO4 dihydroxylation'],
  test: c => c.has('OsO4') || (c.has('KMnO4') && (c.cold || c.flag('dil') || c.flag('mild') || c.has('KOH') || c.has('aqKOH')) && !c.heat && !c.has('H2SO4')),
  run(m) {
    const al = alkenes(m); if (!al.length) return null; const [a, b] = al[0];
    const p = E.transform(m, t => { t.bondBetween(a, b).o = 1; addOn(t, a, 'O'); addOn(t, b, 'O'); });
    return out({ title: 'Syn dihydroxylation → vicinal diol', product: p, steps: [{ L: m.clone(), arrows: [], R: p.clone(), text: 'MnO₄⁻ (or OsO₄) adds as a cyclic ester on one face; hydrolysis frees both OH on the same side: syn diol. Purple colour vanishes (Baeyer test).', cond: ['cold', 'dil.'] }], why: ['Syn addition', 'Pink/purple KMnO₄ decolourises: test for C=C'], tags: ['syn'] });
  } });
reg({ id: 'kmno4hot', name: 'Hot KMnO₄ / K₂Cr₂O₇ oxidation', fam: 'Oxidation', canon: 'KMnO₄, H⁺, Δ', db: ['Hot KMnO4 cleavage of alkenes', 'Side-chain oxidation (KMnO4)', 'KMnO4 oxidation of alcohols', 'Acidic dichromate oxidation of alcohols'],
  test: c => (c.has('KMnO4') && !(c.cold || c.flag('mild'))) || c.has('K2Cr2O7'),
  run(m, c) {
    const AR = aromSet(m); const al = alkenes(m);
    const benzylic = M(m, '[c][CX4;!H0]').concat(M(m, '[c][CX4][#6]')); const side = M(m, '[c:1][CX4;!H0:2]');
    if (side.length || (AR.size && M(m, 'c[CX4]').length)) {
      if (!side.length) return out({ title: 'No side-chain oxidation', product: null, none: 'The side chain has no benzylic H (e.g. tert-butyl), so KMnO₄ cannot attack it.', why: ['Needs at least one benzylic H'] });
      const p = E.transform(m, t => { for (const [ar, cb] of side) { if (!t.atoms.has(cb)) continue; for (const id of E.subtree(t, cb, ar)) if (id !== cb) t.removeAtom(id); t.atoms.get(cb).el = 'C'; const o1 = addOn(t, cb, 'O', 2); addOn(t, cb, 'O', 1); } });
      return out({ title: 'Side-chain oxidation → benzoic acid', product: p, steps: [{ L: m.clone(), arrows: [], R: p.clone(), text: 'Every alkyl chain with a benzylic H is cut back to –COOH on the ring, whatever its length.', cond: ['Δ'] }], why: ['Needs a benzylic H', 'Chain length does not matter: always ArCOOH'], tags: ['ox'] });
    }
    if (al.length && c.has('KMnO4')) {
      const r = E.REC.find(x => x.id === 'ozone').run(m, c, [{ has: k => k === 'H2O2' }]);
      if (!r) return null;
      // terminal =CH2 -> CO2 ; HCOOH -> CO2
      const fix = mm => { if (!mm) return mm; const s = C.canon(mm); return s === 'O=CO' || s === 'OC=O' ? null : mm; };
      const prods = [r.product, ...r.extra].map(fix).filter(Boolean);
      r.title = 'Hot KMnO₄: oxidative cleavage of C=C'; r.product = prods[0]; r.extra = prods.slice(1); r.byLabel = ['CO₂'];
      r.steps[0].text = 'Hot acidic/alkaline KMnO₄ cleaves the C=C: =CHR → RCOOH, =CR₂ → ketone, =CH₂ → CO₂ + H₂O.'; r.why = ['=CH₂ ends become CO₂', '=CHR become acids, =CR₂ become ketones'];
      return r;
    }
    const alc = M(m, '[CX4:1][OX2H1:2]'); const ald = M(m, '[CX3H1:1](=O)');
    if (alc.length) {
      const [cc, o] = alc[0]; const n = hTot(m, cc);
      if (n === 0) return out({ title: 'No oxidation (3° alcohol)', product: null, none: '3° alcohols have no H on the carbinol carbon, so they are not oxidised (under mild conditions).', why: ['Needs H on the C–OH carbon'] });
      const p = E.transform(m, t => { t.bondBetween(cc, o).o = 2; if (n >= 2) addOn(t, cc, 'O', 1); });
      return out({ title: n >= 2 ? '1° alcohol → carboxylic acid' : '2° alcohol → ketone', product: p, steps: [{ L: m.clone(), arrows: [], R: p.clone(), text: n >= 2 ? 'Strong oxidants take a 1° alcohol all the way: aldehyde → acid (use PCC to stop at the aldehyde). Orange Cr₂O₇²⁻ turns green.' : 'A 2° alcohol stops at the ketone.', cond: [] }], why: ['1° → RCOOH, 2° → ketone, 3° → no reaction'], tags: ['ox'] });
    }
    if (ald.length) { const [cc] = ald[0]; const p = E.transform(m, t => { addOn(t, cc, 'O', 1); }); return out({ title: 'Aldehyde → carboxylic acid', product: p, steps: [{ L: m.clone(), arrows: [], R: p.clone(), text: 'The aldehyde C–H is oxidised to C–OH.', cond: [] }], why: ['Aldehydes oxidise easily; ketones resist'], tags: ['ox'] }); }
    return null;
  } });
reg({ id: 'epox', name: 'Epoxidation (peracid)', fam: 'Alkene', canon: 'mCPBA', db: ['Epoxidation', 'Prilezhaev'],
  test: c => c.has('mCPBA'),
  run(m) {
    const al = alkenes(m); if (!al.length) return null; const [a, b] = al[0]; const mech = new Mech(m);
    const pa = mech.add({ atoms: [{ el: 'O' }, { el: 'H', x: 0.55, y: -0.65 }, { el: 'O', x: 1.0, y: 0.2 }, { el: 'C', x: 1.9, y: -0.35 }, { el: 'O', x: 1.85, y: -1.35 }, { el: 'C', x: 2.85, y: 0.25 }], bonds: [[0, 1], [0, 2], [2, 3], [3, 4, 2], [3, 5]] }, a, { dist: 1.45 * BL, dir: (() => { const d = freeDir(m, a); return d; })() });
    const oa = pa.get(1), h = pa.get(2), ob = pa.get(3), cc = pa.get(4), oc = pa.get(5);
    mech.step([[[a, b], [a, oa]], [[oa, ob], [ob, cc]], [[cc, oc], [oc, h]], [[h, oa], [oa, b]]], 'One concerted "butterfly" step: the π bond takes the outer O, the weak O–O bond breaks, and the H moves to the carbonyl O.', ['RCO₃H']);
    return done(mech, { title: 'Epoxidation (syn, stereospecific)', why: ['Concerted, syn: cis-alkene → cis-epoxide'], tags: ['concerted'] });
  } });
/* radical substitution: Cl2/Br2 + hv, NBS */
reg({ id: 'radx', name: 'Free-radical halogenation', fam: 'Radical', canon: 'Cl₂, hν', db: ['Photochlorination of alkane', 'Photobromination of alkane (selective)', 'Wohl–Ziegler (NBS allylic/benzylic bromination)', 'Side-chain chlorination of toluene'],
  test: c => (c.any(['Cl2', 'Br2', 'SO2Cl2']) && (c.light || c.heat) && !c.has('FeX3')) || c.has('NBS'),
  run(m, c) {
    const AR = aromSet(m); const nbs = c.has('NBS'); const X = c.has('Cl2') || c.has('SO2Cl2') ? 'Cl' : 'Br';
    let sites = M(m, '[CX4;!H0:1]').map(x => x[0]);
    if (!sites.length) return null;
    const allyl = sites.filter(s => heavyN(m, s).some(n => AR.has(n) || m.bondsOf(n).some(b => b.o >= 2 && el(m, m.other(b, n)) === 'C')));
    if (nbs) { if (!allyl.length) return null; sites = allyl; }
    if (alkenes(m).length && !nbs && !c.heat) return null;
    const RATE = X === 'Br' ? [1, 1, 82, 1600] : [1, 1, 3.8, 5];
    const rate = s => hTot(m, s) * RATE[Math.min(3, cN(m, s).length)] * (allyl.includes(s) ? (X === 'Br' ? 5000 : 30) : 1);
    const grp = new Map();
    for (const s of sites) { const t = m.clone(); addOn(t, s, X); let k; try { k = C.canon(t); } catch (e) { k = 's' + s; } const g = grp.get(k) || { k, sites: [], w: 0, mol: t }; g.sites.push(s); g.w += rate(s); grp.set(k, g); }
    const G = [...grp.values()].sort((p, q) => q.w - p.w); const tot = G.reduce((a, g) => a + g.w, 0) || 1;
    const s = G[0].sites[0]; const pct = g => Math.round(100 * g.w / tot);
    const mixNote = G.length > 1 ? `Expected mixture: this isomer ≈ ${pct(G[0])}%, ${G.slice(1, 4).map(pct).join('%, ')}% for the other${G.length > 2 ? 's' : ''} (number of H × reactivity).` : '';
    const minorMols = G.slice(1, 4).map(g => E.cleanProduct(g.mol, E.fragIds(g.mol, g.sites[0])));
    const mech = new Mech(m);
    let ids = mech.add(`${X}${X}`, s, { dist: 3.4 * BL }); const x1 = ids.get(1), x2 = ids.get(2);
    mech.step([[[x1, x2], [x1], 1], [[x1, x2], [x2], 1]], `Initiation: light splits ${X}₂ into two ${X}• radicals (homolysis).`, [nbs ? 'NBS, hν' : 'hν']);
    mech.dropFrag(x2);
    const h = mech.H(s); moveFragNear(mech, x1, h, 1.1 * BL, unit([P(mech.m, h)[0] - P(mech.m, s)[0], P(mech.m, h)[1] - P(mech.m, s)[1]]));
    mech.step([[[x1], [x1, h], 1], [[s, h], [x1, h], 1], [[s, h], [s], 1]], `${X}• pulls off the H that leaves the most stable radical (${allyl.includes(s) ? (AR.size ? 'benzylic' : 'allylic') : deg(m, s)}).`);
    mech.dropFrag(x1);
    ids = mech.add(`${X}${X}`, s, { dist: 1.6 * BL }); const y1 = ids.get(1), y2 = ids.get(2);
    mech.step([[[s], [s, y1], 1], [[y1, y2], [s, y1], 1], [[y1, y2], [y2], 1]], `The carbon radical takes ${X} from ${nbs ? 'Br₂ (NBS keeps it at a trace level)' : X + '₂'}, regenerating ${X}•: chain propagation.`);
    return done(mech, { title: nbs ? 'NBS: allylic / benzylic bromination' : `Free-radical ${X === 'Cl' ? 'chlorination' : 'bromination'}`, minor: nbs ? [] : minorMols, why: [X === 'Br' ? 'Br• is very selective: 3° > 2° > 1° (≈1600 : 82 : 1)' : 'Cl• is reactive and less selective: mixtures, ratio ∝ (number of H) × (3° 5 : 2° 3.8 : 1° 1)', mixNote, 'Benzylic/allylic H react fastest'], tags: ['SR'] });
  } });
/* terminal alkynes: acetylide, alkylation, metal acetylides */
reg({ id: 'acetylide', name: 'Terminal alkyne + strong base', fam: 'Alkyne', canon: 'NaNH₂ then CH₃Br', db: ['Terminal alkyne + Na/liq. NH3', 'Wolf bottle (ammoniacal AgNO3)'], eats: s => s.has('RX'),
  test: c => c.any(['NaNH2', 'Na', 'Tollens', 'AgNO3', 'CuCl']) && !c.has('Na/NH3'),
  run(m, c, all) {
    const t = M(m, '[CH1:1]#[C:2]'); if (!t.length) return null; const [a, b] = t[0];
    const metal = c.any(['Tollens', 'AgNO3']) ? 'Ag' : c.has('CuCl') ? 'Cu' : null;
    if (metal) { const p = E.transform(m, tt => { addOn(tt, a, metal); }); return out({ title: `${metal === 'Ag' ? 'White' : 'Red'} precipitate: metal acetylide`, product: p, steps: [{ L: m.clone(), arrows: [], R: p.clone(), text: `Only terminal alkynes (acidic ≡C–H) react: ${metal === 'Ag' ? 'Tollens/ammoniacal AgNO₃ gives a white silver acetylide' : 'ammoniacal Cu₂Cl₂ gives a red copper acetylide'}. Test to tell terminal from internal alkynes.`, cond: [] }], why: ['sp C–H is acidic (pKa ≈ 25)'], tags: ['test'] }); }
    const mech = new Mech(m); const h = mech.H(a);
    const base = mech.add('[NH2-]', h, { dist: 1.4 * BL, dir: unit([P(mech.m, h)[0] - P(mech.m, a)[0], P(mech.m, h)[1] - P(mech.m, a)[1]]) }).get(1);
    mech.step([[[base], [base, h]], [[a, h], [a]]], 'NH₂⁻ removes the acidic terminal H (sp carbon holds the charge well).');
    mech.dropFrag(base);
    const rx = (all || []).map(s => s.get && s.get('RX')).find(Boolean) || c.get('RX');
    if (!rx) return done(mech, { title: 'Sodium acetylide', why: ['Add a 1° alkyl halide next to make a longer alkyne'], tags: ['acidbase'] });
    const sp = E.species((rx.smi || 'C') + rx.X); const ids = mech.add(sp, a, { dist: 1.6 * BL }); const cR = ids.get(1); const xR = ids.get(sp.atoms.size);
    mech.step([[[a], [a, cR]], [[cR, xR], [xR]]], 'The acetylide carbon does SN2 on the alkyl halide (1° works best; 2°/3° give elimination).');
    return done(mech, { title: 'Alkylation of acetylide → internal alkyne', why: ['SN2 needs a 1° (or methyl) halide'], tags: ['SN2'] });
  } });

/* ===================== ALKYL HALIDES ===================== */
function halides(m) { const AR = aromSet(m); return M(m, '[CX4:1][Cl,Br,I:2]').filter(([cc]) => !AR.has(cc)); }
function sn1ok(m, cc, lg) { let t = m; if (lg != null) { t = m.clone(); const b = t.bondBetween(cc, lg); if (b) t.removeBond(b.id); } return cationScore(t, cc, aromSet(t)) >= 2.6; } // lg: leaving group atom (not counted as a donor)
function betaHs(m, cc, x) { return cN(m, cc).filter(b => b !== x && hTot(m, b) > 0 && !aromSet(m).has(b)); }
function doSN2(m, cc, x, nuSpec, text, nuAtt) {
  const mech = new Mech(m); const back = unit([P(m, cc)[0] - P(m, x)[0], P(m, cc)[1] - P(m, x)[1]]);
  const ids = mech.add(nuSpec, cc, { dist: 1.65 * BL, dir: back, att: nuAtt }); const nu = ids.get(nuAtt || 1);
  mech.step([[[nu], [nu, cc]], [[cc, x], [x]]], text || 'Back-side attack: Nu in, leaving group out, one step (inversion).');
  return mech;
}
function doSN1(m, cc, x, nuSpec, nuName, deprot) {
  const mech = new Mech(m);
  mech.step([[[cc, x], [x]]], `The C–${el(m, x)} bond breaks by itself: a ${deg(m, cc)} carbocation (slow step).`);
  let cp = cc; const sh = bestShift(mech.m, cp); if (sh) cp = doShift(mech, cp, sh);
  const ids = mech.add(nuSpec, cp, { dist: 1.65 * BL }); const nu = ids.get(1);
  mech.step([[[nu], [nu, cp]]], `${nuName} attacks the flat carbocation from either face (racemisation).`);
  if (deprot) deprotOn(mech, nu, 'Loss of H⁺ gives the neutral product.');
  return mech;
}
/* stereo on rings from wedge/hash bonds: +1 = towards the viewer, -1 = away, 0 = not drawn */
function faceOf(m, c, s) {
  const bd = m.bondBetween(c, s); if (bd && bd.st && bd.a === c) return bd.st === 1 ? 1 : bd.st === 6 ? -1 : 0;
  const ringN = new Set((C.rings(m, 8).find(r => r.includes(c)) || []));
  const others = m.nbrs(c).filter(n => n !== s && !ringN.has(n));
  for (const n of others) { const b2 = m.bondBetween(c, n); if (b2 && b2.st && b2.a === c && others.length === 1) return b2.st === 1 ? -1 : b2.st === 6 ? 1 : 0; }
  return 0;
}
/* E2 needs H and leaving group anti-periplanar: on a cyclohexane ring that means trans (both axial) */
function antiOK(m, cc, x, b) {
  const ring = C.rings(m, 6).find(r => r.length === 6 && r.includes(cc) && r.includes(b)); if (!ring) return true;
  const fX = faceOf(m, cc, x); if (!fX) return true; if (hTot(m, b) >= 2) return true;
  const hs = m.nbrs(b).filter(n => el(m, n) === 'H'); for (const h of hs) { const f = faceOf(m, b, h); if (f) return f === -fX; }
  const subs = m.nbrs(b).filter(n => !ring.includes(n) && el(m, n) !== 'H'); if (subs.length !== 1) return true;
  const fR = faceOf(m, b, subs[0]); if (!fR) return true; return fR === fX; // the lone H is trans to X when R is cis to X
}
E.faceOf = faceOf; E.antiOK = antiOK;
function elim(m, cc, x, base, bulky, text, forceB) {
  const AR = aromSet(m); const bs0 = forceB != null ? [forceB].filter(b => hTot(m, b) > 0) : betaHs(m, cc, x); if (!bs0.length) return null;
  const antiL = forceB != null ? bs0 : bs0.filter(b => antiOK(m, cc, x, b)); const bs = antiL.length ? antiL : bs0;
  bs.sort((p, q) => bulky ? (hTot(m, q) - hTot(m, p)) || (alkeneScore(m, cc, p, AR) - alkeneScore(m, cc, q, AR)) : (alkeneScore(m, cc, q, AR) - alkeneScore(m, cc, p, AR)));
  const bC = bs[0]; const mech = new Mech(m); const h = mech.H(bC);
  const bid = mech.add(base, h, { dist: 1.45 * BL, dir: unit([P(mech.m, h)[0] - P(mech.m, bC)[0], P(mech.m, h)[1] - P(mech.m, bC)[1]]) }).get(1);
  mech.step([[[bid], [bid, h]], [[bC, h], [bC, cc]], [[cc, x], [x]]], text || 'E2, one step: the base takes a β-H, the C=C forms, and the halide leaves (anti-periplanar).');
  let antiForced = false; if (bs.length < bs0.length) { const best0 = bs0.slice().sort((p, q) => bulky ? (hTot(m, q) - hTot(m, p)) : (alkeneScore(m, cc, q, AR) - alkeneScore(m, cc, p, AR)))[0]; antiForced = best0 !== bC; }
  return { mech, bC, minorPossible: bs.length > 1, antiForced };
}
E.halides = halides; E.sn1ok = sn1ok; E.betaHs = betaHs; E.doSN2 = doSN2; E.doSN1 = doSN1; E.elimX = elim;
reg({ id: 'rx_oh', name: 'Hydrolysis of alkyl halide', fam: 'Alkyl halide', canon: 'aq. KOH', db: ['SN2 reaction', 'SN1 reaction', 'Moist / dry Ag2O with alkyl halide'],
  test: c => c.any(['aqKOH', 'Ag2O']) || (c.has('KOH') && !c.flag('alc') && !c.has('alcKOH')) || (c.has('H2O') && !c.any(['H2SO4', 'H3O+', 'Br2', 'Cl2', 'HOBr', 'HOCl'])),
  run(m, c) {
    const hs = halides(m); if (!hs.length) return null; const [cc, x] = hs[0];
    const weak = !c.any(['KOH', 'aqKOH']);
    if (sn1ok(m, cc, x) || (weak && cationScore(m, cc) >= 2)) {
      const mech = doSN1(m, cc, x, weak ? 'O' : '[OH-]', weak ? 'Water' : 'OH⁻', weak); mech.dropFrag(x);
      return done(mech, { title: `SN1 hydrolysis (${deg(m, cc)} halide)`, why: ['3°/benzylic/allylic halides ionise: SN1', 'Rearrangement possible, racemisation at a stereocentre'], tags: ['SN1'] });
    }
    const mech = doSN2(m, cc, x, '[OH-]', 'OH⁻ attacks from the back as the halide leaves: SN2, one step, inversion of configuration.');
    return done(mech, { title: `SN2 hydrolysis (${deg(m, cc)} halide) → alcohol`, why: ['1°/2° halide + strong nucleophile: SN2', 'Walden inversion'], tags: ['SN2'] });
  } });
reg({ id: 'rx_elim', name: 'Dehydrohalogenation (E2)', fam: 'Alkyl halide', canon: 'alc. KOH', db: ['E2 elimination', 'Saytzeff', 'Hofmann elimination'],
  test: c => c.any(['alcKOH', 'KOtBu', 'NaNH2']) || (c.has('KOH') && c.flag('alc')) || (c.any(['NaOEt', 'NaOMe', 'NaOR']) && (c.heat || true)),
  run(m, c) {
    const hs = halides(m); if (!hs.length) return null; const [cc, x] = hs[0];
    const bulky = c.has('KOtBu'); const alk = c.any(['NaOEt', 'NaOMe', 'NaOR']) && !c.has('alcKOH');
    if (alk && !bulky && deg(m, cc) !== '3°' && deg(m, cc) !== '2°' && !c.heat) return null; // 1° + alkoxide: Williamson
    const base = bulky ? '[O-]C(C)(C)C' : c.has('NaOEt') ? '[O-]CC' : c.has('NaNH2') ? '[NH2-]' : c.has('NaOMe') ? '[O-]C' : '[OH-]';
    // gem/vic dihalide + NaNH2 -> alkyne
    const r = elim(m, cc, x, base, bulky, bulky ? 'E2 with a bulky base: it reaches the least hindered β-H (Hofmann product).' : null);
    if (!r) return out({ title: 'No β-H: no elimination', product: null, none: 'This halide has no β-hydrogen, so E2 cannot happen (SN2 instead).', why: [] });
    if (r.antiForced) return done(r.mech, { title: 'E2: only the anti-periplanar β-H can leave → less substituted alkene', why: ['E2 needs H and X anti-periplanar: on a cyclohexane ring, trans and both axial', 'The H that would give the Saytzeff alkene is cis to the leaving group here, so that alkene cannot form', 'Classic case: menthyl chloride → 2-menthene only; neomenthyl chloride → mostly 3-menthene'], tags: ['E2', 'anti'] });
    const o = done(r.mech, { title: bulky ? 'E2 with t-BuOK → less substituted alkene (Hofmann)' : 'E2 → more substituted alkene (Saytzeff)', why: [bulky ? 'Bulky base: Hofmann (less substituted) alkene' : 'Small strong base, alcoholic solvent: elimination, Saytzeff (more substituted) alkene', 'trans alkene preferred over cis'], tags: ['E2'] });
    if (c.has('NaNH2')) {
      const v = M(o.product, '[CX3:1]([Cl,Br,I:2])=[#6]');
      if (v.length) { const [vc, vx] = v[0]; const r2 = elim(o.product, vc, vx, '[NH2-]', false, 'Second E2 on the vinyl halide: NH₂⁻ takes the β-H and the C=C becomes C≡C.'); if (r2) { const o2 = done(r2.mech, { title: 'Double E2 (NaNH₂) → alkyne', why: ['Vicinal or geminal dihalide + 2 NaNH₂: two E2 steps', 'The second E2 (from a vinyl halide) needs the very strong base NaNH₂'], tags: ['E2'] }); o2.steps = o.steps.concat(o2.steps); return o2; } }
      if (halides(o.product).length) { const r2 = E.REC.find(z => z.id === 'rx_elim').run(o.product, c); if (r2 && r2.product) { r2.steps = o.steps.concat(r2.steps); r2.title = 'Double E2 (NaNH₂) → alkyne'; return r2; } }
    }
    return o;
  } });
const NUC = [ // reagent id, nucleophile spec (atom 1 attacks), product name
  ['NaCN', '[C-]#N', 'Nitrile (–CN)', 'C attacks: alkyl cyanide (one C longer chain)'],
  ['AgCN', '[N]#[C-]', 'Isocyanide (–NC)', 'AgCN is covalent: the N end attacks, giving the isocyanide'],
  ['KNO2', '[O-]N=O', 'Alkyl nitrite (–O–N=O)', 'Ionic KNO₂: O attacks'], ['AgNO2', '[N](=O)[O-]', 'Nitroalkane (–NO₂)', 'Covalent AgNO₂: N attacks'],
  ['NaI', '[I-]', 'Alkyl iodide (Finkelstein)', 'NaI in dry acetone: NaCl/NaBr precipitate and drive it forward'],
  ['AgF', '[F-]', 'Alkyl fluoride (Swarts)', 'AgF / Hg₂F₂ swap Cl/Br for F'], ['NaSH', '[SH-]', 'Thiol', 'HS⁻ is a strong soft nucleophile'],
  ['NaN3', '[N-]=[N+]=[N-]', 'Alkyl azide', 'Azide does SN2'], ['NH3', 'N', 'Primary amine', 'NH₃ attacks; excess NH₃ limits over-alkylation'],
  ['NaOAc', '[O-]C(C)=O', 'Ester', 'Carboxylate O attacks']];
reg({ id: 'rx_nu', name: 'SN2 with nucleophiles', fam: 'Alkyl halide', canon: 'NaCN', db: ['Finkelstein reaction', 'Swarts reaction', 'Williamson ether synthesis'],
  test: c => NUC.some(n => c.has(n[0])) || c.any(['NaOEt', 'NaOMe', 'NaOR']) || c.has('RNH2'),
  run(m, c) {
    const hs = halides(m); if (!hs.length) return null; const [cc, x] = hs[0];
    let nu = NUC.find(n => c.has(n[0]));
    if (!nu && c.any(['NaOEt', 'NaOMe', 'NaOR'])) { if (deg(m, cc) === '3°') return null; const sm = c.has('NaOEt') ? 'CC' : c.has('NaOMe') ? 'C' : (c.get('NaOR').smi || 'CC'); nu = ['NaOR', '[O-]' + sm, 'Ether (Williamson)', 'Alkoxide does SN2; use a 1° halide, otherwise E2 wins']; }
    if (!nu) return null;
    if (nu[0] === 'NH3') { const mech = doSN2(m, cc, x, 'N', 'NH₃ attacks from the back: SN2.'); deprotOn(mech, mech.mainFrag().find(id => mech.m.atoms.get(id).el === 'N'), 'A second NH₃ removes H⁺ from the ammonium ion.', 'N'); return done(mech, { title: '1° amine (ammonolysis)', why: [nu[3]], tags: ['SN2'] }); }
    if (deg(m, cc) === '3°' && nu[0] !== 'NaI') return null;
    const mech = doSN2(m, cc, x, nu[1], `${nu[3]}. Back-side attack, inversion.`);
    return done(mech, { title: nu[2], why: [nu[3]], tags: ['SN2'] });
  } });
reg({ id: 'grignard_make', name: 'Grignard reagent formation', fam: 'Alkyl halide', canon: 'Mg, dry ether', db: ['Grignard reagent preparation'],
  test: c => c.has('Mg') && !c.has('RMgX'),
  run(m) { const hs = M(m, '[#6:1][Cl,Br,I:2]'); if (!hs.length) return null; const [cc, x] = hs[0];
    const p = E.transform(m, t => { const xp = P(t, x), cp = P(t, cc); t.removeBond(t.bondBetween(cc, x).id); const mg = t.addAtom('Mg', xp[0], xp[1]); t.addBond(cc, mg, 1); const d = unit([xp[0] - cp[0], xp[1] - cp[1]]); const xa = t.atoms.get(x); xa.x = xp[0] + d[0] * BL; xa.y = xp[1] + d[1] * BL; t.addBond(mg, x, 1); });
    return out({ title: 'Grignard reagent (R–MgX)', product: p, steps: [{ L: m.clone(), arrows: [], R: p.clone(), text: 'Mg inserts into the C–X bond (dry ether; water would destroy it). The carbon becomes strongly nucleophilic/basic.', cond: ['dry ether'] }], why: ['Must be dry: R–MgX + H₂O → R–H'], tags: ['organometallic'] }); } });
reg({ id: 'wurtz', name: 'Wurtz / Fittig coupling', fam: 'Alkyl halide', canon: 'Na, dry ether', db: ['Wurtz reaction', 'Fittig reaction'],
  test: c => c.has('Na') && (c.flag('dry') || c.flag('ether')),
  run(m) { const hs = M(m, '[#6:1][Cl,Br,I:2]'); if (!hs.length) return null; const [cc, x] = hs[0];
    const p = E.transform(m, t => { const cp = P(t, cc), xp = P(t, x); const d = unit([xp[0] - cp[0], xp[1] - cp[1]]); t.removeAtom(x); const ids = [...t.atoms.keys()]; const copy = t.sub(ids); const idm = t.merge(copy, 0, 0); const c2 = idm.get(cc); const all2 = [...idm.values()]; const pp = P(t, c2); rotSet(t, all2, pp[0], pp[1], Math.PI); const p2 = P(t, c2); moveSet(t, all2, cp[0] + d[0] * BL - p2[0], cp[1] + d[1] * BL - p2[1]); t.addBond(cc, c2, 1); });
    return out({ title: 'Wurtz coupling: R–X + X–R → R–R', product: p, steps: [{ L: m.clone(), arrows: [], R: p.clone(), text: 'Na gives electrons; two alkyl groups join (radical/carbanion route). Same halide → symmetric alkane with double the carbons. Aryl halides do this as the Fittig reaction.', cond: ['dry ether'] }], why: ['Doubles the carbon count', 'Mixed halides give mixtures'], tags: ['coupling'] }); } });

/* ===================== ALCOHOLS / ETHERS ===================== */
function alcohols(m) { const AR = aromSet(m); return M(m, '[CX4:1][OX2H1:2]').filter(([cc]) => !AR.has(cc)); }
reg({ id: 'dehyd', name: 'Acid dehydration', fam: 'Alcohol', canon: 'conc. H₂SO₄, Δ', db: ['Acid-catalysed dehydration of alcohol', 'E1 elimination'],
  test: c => (c.any(['H2SO4', 'H3PO4', 'Al2O3', 'KHSO4', 'pTsOH']) && (c.heat || c.flag('conc'))) || c.has('POCl3'),
  run(m, c) {
    const alc = alcohols(m); const diol = M(m, '[CX4:1]([OX2H1:2])[CX4:3][OX2H1:4]');
    if (diol.length && !c.has('POCl3')) return E.REC.find(r => r.id === 'pinacol').run(m, c);
    if (!alc.length) return null; const [cc, o] = alc[0];
    const mech = new Mech(m); const AR = aromSet(m);
    let ids = mech.add('[H][OH2+]', o, { dist: 1.5 * BL }); const h = ids.get(1), w = ids.get(2);
    mech.step([[[o], [o, h]], [[h, w], [w]]], 'The OH grabs H⁺: –OH₂⁺ is a good leaving group.'); mech.dropFrag(w);
    mech.step([[[cc, o], [o]]], `Water leaves: a ${deg(mech.m, cc)} carbocation (slow step).`); mech.dropFrag(o);
    let cp = cc; const sh = bestShift(mech.m, cp); const why = []; if (sh) { cp = doShift(mech, cp, sh); why.push('The carbocation rearranges before losing H⁺'); }
    const bs = cN(mech.m, cp).filter(b => hTot(mech.m, b) > 0 && !aromSet(mech.m).has(b)); if (!bs.length) return null;
    bs.sort((p, q) => alkeneScore(mech.m, cp, q, aromSet(mech.m)) - alkeneScore(mech.m, cp, p, aromSet(mech.m)));
    const bC = bs[0]; const hb = mech.H(bC); const wb = mech.add('O', hb, { dist: 1.4 * BL, dir: unit([P(mech.m, hb)[0] - P(mech.m, bC)[0], P(mech.m, hb)[1] - P(mech.m, bC)[1]]) }).get(1);
    mech.step([[[wb], [wb, hb]], [[bC, hb], [bC, cp]]], 'A base (H₂O / HSO₄⁻) removes a β-H: the more substituted alkene forms (Saytzeff).');
    return done(mech, { title: 'Dehydration (E1) → alkene', why: ['E1 via carbocation: 3° > 2° > 1°', 'Saytzeff alkene, rearrangement possible', ...why], tags: ['E1'] });
  } });
reg({ id: 'pinacol', name: 'Pinacol–pinacolone rearrangement', fam: 'Alcohol', canon: 'H⁺, Δ', db: ['Pinacol–pinacolone rearrangement'],
  test: c => c.any(['H2SO4', 'H3O+']) && (c.heat || c.flag('conc')),
  run(m) {
    const d = M(m, '[CX4:1]([OX2H1:2])[CX4:3][OX2H1:4]'); if (!d.length) return null;
    let [c1, o1, c2, o2] = d[0]; if (cationScore(m, c2) > cationScore(m, c1)) [c1, o1, c2, o2] = [c2, o2, c1, o1];
    const mech = new Mech(m);
    let ids = mech.add('[H][OH2+]', o1, { dist: 1.5 * BL }); const h = ids.get(1), w = ids.get(2);
    mech.step([[[o1], [o1, h]], [[h, w], [w]]], 'Protonate the OH whose loss gives the more stable carbocation.'); mech.dropFrag(w);
    mech.step([[[c1, o1], [o1]]], 'Water leaves: carbocation next to the second C–OH.'); mech.dropFrag(o1);
    const AR = aromSet(mech.m); const gs = heavyN(mech.m, c2).filter(g => g !== c1 && g !== o2 && el(mech.m, g) === 'C');
    gs.sort((p, q) => (AR.has(q) ? 2 : 0) + (hTot(mech.m, q) === 0 ? 0.5 : 0) - ((AR.has(p) ? 2 : 0) + (hTot(mech.m, p) === 0 ? 0.5 : 0)));
    let g = gs[0]; const typ = g == null ? 'H' : AR.has(g) ? 'aryl' : 'alkyl'; if (g == null) g = mech.H(c2);
    mech.step([[[c2, g], [c1, g]], [[o2], [o2, c2]]], `1,2-${typ === 'aryl' ? 'aryl' : typ === 'H' ? 'hydride' : 'alkyl'} shift, pushed by the O lone pair: the new cation is an oxocarbenium ion (very stable).`, [], () => moveGroupTo(mech.m, g, c2, c1));
    deprotOn(mech, o2, 'Loss of H⁺ from O gives the ketone (pinacolone).');
    return done(mech, { title: 'Pinacol–pinacolone rearrangement', why: ['Migratory aptitude: aryl > H > alkyl (roughly)', 'Driving force: C=O formation'], tags: ['1,2-shift'] });
  } });
reg({ id: 'alc_hx', name: 'Alcohol + HX', fam: 'Alcohol', canon: 'HBr', db: ['Lucas test'],
  test: c => c.any(['HBr', 'HI', 'Lucas']) || (c.has('HCl') && c.has('ZnCl2')) || (c.has('HCl') && c.flag('conc')),
  run(m, c) {
    const alc = alcohols(m); if (!alc.length) return null; const [cc, o] = alc[0];
    const X = c.has('HBr') ? 'Br' : c.has('HI') ? 'I' : 'Cl'; const lucas = c.has('ZnCl2') || c.has('Lucas');
    const mech = new Mech(m);
    let ids = mech.add(`[H]${X}`, o, { dist: 1.45 * BL }); const h = ids.get(1), xx = ids.get(2);
    mech.step([[[o], [o, h]], [[h, xx], [xx]]], lucas ? 'ZnCl₂/H⁺ activates the OH (Lucas reagent).' : `O grabs H⁺ from H${X}: –OH₂⁺ can leave as water.`);
    if (sn1ok(m, cc, o) || deg(m, cc) === '2°' && !c.has('HI')) {
      mech.step([[[cc, o], [o]]], `Water leaves: ${deg(m, cc)} carbocation.`); mech.dropFrag(o);
      let cp = cc; const sh = bestShift(mech.m, cp); if (sh) cp = doShift(mech, cp, sh);
      captureBy(mech, cp, xx, `${X}⁻ attacks the carbocation.`);
      return done(mech, { title: `SN1: alcohol → alkyl ${X === 'Cl' ? 'chloride' : X === 'Br' ? 'bromide' : 'iodide'}`, why: ['3°/2° alcohols react by SN1 (rearrangement possible)', lucas ? 'Lucas test: 3° turbid at once, 2° in ~5 min, 1° no turbidity at room T' : 'Reactivity HI > HBr > HCl'], tags: ['SN1'] });
    }
    moveFragNear(mech, xx, cc, 1.65 * BL, unit([P(mech.m, cc)[0] - P(mech.m, o)[0], P(mech.m, cc)[1] - P(mech.m, o)[1]]));
    mech.step([[[xx], [xx, cc]], [[cc, o], [o]]], `${X}⁻ pushes out water from the back: SN2 (1° alcohol).`);
    return done(mech, { title: 'SN2: 1° alcohol → alkyl halide', why: ['1° alcohols: SN2, no rearrangement', lucas ? 'Lucas test negative for 1° at room temperature' : ''], tags: ['SN2'] });
  } });
reg({ id: 'alc_pcl', name: 'SOCl₂ / PCl₅ / PBr₃', fam: 'Alcohol', canon: 'SOCl₂, pyridine', db: ['Darzens reaction (SOCl2, SNi)'],
  test: c => c.any(['SOCl2', 'PCl5', 'PCl3', 'PBr3', 'PI3']) && !c.any(['NH3']),
  run(m, c) {
    const X = c.has('PBr3') ? 'Br' : c.has('PI3') ? 'I' : 'Cl';
    const acid = M(m, '[CX3:1](=O)[OX2H1:2]');
    if (acid.length) { const [cc, , o] = acid[0]; const p = E.transform(m, t => { t.atoms.get(o).el = X; }); return out({ title: 'Carboxylic acid → acyl chloride', product: p, steps: [{ L: m.clone(), arrows: [], R: p.clone(), text: `The –OH of the acid is replaced by ${X}. With SOCl₂ the by-products SO₂ and HCl are gases (clean reaction).`, cond: [] }], why: ['Acyl chlorides are the most reactive acid derivatives'], tags: ['sub'] }); }
    const alc = M(m, '[#6:1][OX2H1:2]'); if (!alc.length) return null; const [cc, o] = alc[0];
    if (aromSet(m).has(cc)) return out({ title: 'Phenol: no easy replacement', product: null, none: 'The C–O bond of phenol has partial double-bond character, so PCl₅/SOCl₂ do not replace the OH easily.', why: [] });
    const p = E.transform(m, t => { t.atoms.get(o).el = X; });
    return out({ title: `Alcohol → alkyl ${X === 'Cl' ? 'chloride' : X === 'Br' ? 'bromide' : 'iodide'}`, product: p, steps: [{ L: m.clone(), arrows: [], R: p.clone(), text: c.has('SOCl2') ? 'R–OH + SOCl₂ → R–Cl + SO₂↑ + HCl↑ (Darzens). No carbocation, so no rearrangement; SNi with retention, inversion if pyridine is added.' : `${X} replaces OH without rearrangement.`, cond: [] }], why: ['No rearrangement (no free carbocation)'], tags: ['sub'] });
  } });
reg({ id: 'alc_ox_mild', name: 'Mild oxidation (PCC)', fam: 'Alcohol', canon: 'PCC, CH₂Cl₂', db: ['PCC / PDC / Collins oxidation', 'Cu / 300 °C dehydrogenation', 'MnO2 oxidation (allylic/benzylic)'],
  test: c => c.any(['PCC', 'CrO2Cl2']) || (c.has('Cu') && (c.heat || (c.T && c.T > 500))),
  run(m, c) {
    const alc = alcohols(m); if (!alc.length) return null; const [cc, o] = alc[0]; const n = hTot(m, cc);
    if (n === 0) { if (c.has('Cu')) { const r = E.REC.find(x => x.id === 'dehyd'); const o2 = r.run(m, { has: k => k === 'H2SO4', any: ks => ks.includes('H2SO4'), heat: true, flag: () => false, get: () => null, raw: '' }); if (o2) { o2.title = '3° alcohol + Cu/573 K → alkene (dehydration)'; return o2; } } return out({ title: 'No reaction (3° alcohol)', product: null, none: 'No H on the carbinol carbon, so PCC cannot oxidise a 3° alcohol.', why: ['1° → aldehyde, 2° → ketone, 3° → none'] }); }
    const mech = new Mech(m);
    const cr = mech.add({ atoms: [{ el: 'Cr' }, { el: 'O', x: 0.7, y: -0.7 }, { el: 'O', x: 0.7, y: 0.7 }, { el: 'O', x: 1, y: 0 }, { el: 'Cl', x: -0.1, y: 1 }], bonds: [[0, 1, 2], [0, 2, 2], [0, 3, 1], [0, 4, 1]] }, o, { dist: 1.6 * BL }).get(1);
    // chromate ester (shown as one step), then E2-like
    const L = mech.m.clone(); const mm = mech.m; const h0 = E.hExp(mm, o)[0] || null;
    mm.addBond(o, cr, 1); mm.removeBond(mm.bondsOf(cr).find(b => el(mm, mm.other(b, cr)) === 'Cl').id); for (const a of [...mm.atoms.values()]) if (a.el === 'Cl' && !mm.nbrs(a.id).length) mm.removeAtom(a.id);
    mech.note('The alcohol O attacks Cr: a chromate ester forms (Cl⁻ and H⁺ lost).', L, mm.clone());
    const h = mech.H(cc); const b = mech.add('n1ccccc1', h, { dist: 1.45 * BL, dir: unit([P(mech.m, h)[0] - P(mech.m, cc)[0], P(mech.m, h)[1] - P(mech.m, cc)[1]]) }).get(1);
    mech.step([[[b], [b, h]], [[cc, h], [cc, o]], [[o, cr], [cr]]], 'Base removes the H on the carbinol carbon; C=O forms and Cr(VI) leaves as Cr(IV).');
    return done(mech, { title: n >= 2 ? '1° alcohol → aldehyde (stops here)' : '2° alcohol → ketone', why: ['PCC is anhydrous: the aldehyde is not over-oxidised', '3° alcohols do not react'], tags: ['ox'] });
  } });
reg({ id: 'alc_na', name: 'Alcohol + Na', fam: 'Alcohol', canon: 'Na', db: [],
  test: c => (c.has('Na') || c.has('NaH')) && !c.flag('dry') && !c.has('NH3'),
  run(m) { const alc = M(m, '[#6:1][OX2H1:2]'); if (!alc.length) return null; const [, o] = alc[0];
    const p = E.transform(m, t => { t.atoms.get(o).chg = -1; }); return out({ title: 'Alkoxide + H₂↑', product: p, steps: [{ L: m.clone(), arrows: [], R: p.clone(), text: 'Na takes the acidic O–H: R–O⁻ Na⁺ + ½ H₂ (bubbles). Acidity: H₂O > 1° > 2° > 3° alcohol.', cond: [] }], why: ['Test for active H (H₂ gas)'], tags: ['acidbase'] }); } });
reg({ id: 'ether_hi', name: 'Ether cleavage (HI)', fam: 'Ether', canon: 'HI', db: ['Cleavage of ethers by HI'],
  test: c => c.any(['HI', 'HBr']),
  run(m, c) {
    const et = M(m, '[#6:1][OX2:2][#6:3]').filter(([a, , b]) => !E.isCarbonylC(m, a) && !E.isCarbonylC(m, b)); if (!et.length || alcohols(m).length) return null;
    let [a, o, b] = et[0]; const AR = aromSet(m); const X = c.has('HI') ? 'I' : 'Br';
    if (AR.has(a)) [a, b] = [b, a]; // b may be aryl
    const mech = new Mech(m);
    const ids = mech.add(`[H]${X}`, o, { dist: 1.45 * BL }); const h = ids.get(1), xx = ids.get(2);
    mech.step([[[o], [o, h]], [[h, xx], [xx]]], `The ether O takes H⁺ from H${X}.`);
    let target;
    if (AR.has(b)) target = a; else if (sn1ok(m, a, o) || sn1ok(m, b, o)) target = sn1ok(m, a, o) ? a : b; else target = cN(m, a).length <= cN(m, b).length ? a : b;
    const sn1 = sn1ok(m, target, o) && !AR.has(target);
    if (sn1) { mech.step([[[target, o], [o]]], `The O–C bond to the ${deg(m, target)} carbon breaks: stable carbocation.`); captureBy(mech, target, xx, `${X}⁻ captures it (SN1).`); }
    else { moveFragNear(mech, xx, target, 1.6 * BL, unit([P(mech.m, target)[0] - P(mech.m, o)[0], P(mech.m, target)[1] - P(mech.m, o)[1]])); mech.step([[[xx], [xx, target]], [[target, o], [o]]], `${X}⁻ attacks the less hindered carbon (SN2).`); }
    const f = E.fragOf(mech.m); const parts = mech.m.frags().filter(ids => ids.some(id => mech.sub.has(id))).map(ids => E.cleanProduct(mech.m, ids)).sort((p, q) => q.atoms.size - p.atoms.size);
    return out({ title: 'Ether cleavage', product: parts[0], extra: parts.slice(1), steps: mech.steps, why: [AR.has(b) ? 'Aryl–O bond never breaks: phenol + alkyl halide' : sn1 ? '3° alkyl: SN1, halide on the 3° carbon' : 'SN2 at the smaller alkyl: CH₃–X + the other alcohol'], tags: [sn1 ? 'SN1' : 'SN2'] });
  } });
})(typeof window !== 'undefined' ? window : globalThis);
