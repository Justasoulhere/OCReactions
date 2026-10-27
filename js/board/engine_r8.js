/* engine_r8.js: curly-arrow mechanisms, batch 2 (Oct 14). Same contract as engine_r7: each wrapper rebuilds the
   reaction as arrow steps on the same atoms; if its product differs from the plain recipe, or any step is not a
   valid structure, the plain version is kept. */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore, R7 = root.EngineR7;
if (!R7) return;
const { BL, el, P, unit, hTot, Mech } = E; const { wrap, makeDouble } = R7;
const M = (m, sm) => { try { return C.matchAll(m, sm, {}); } catch (e) { return []; } };
const relay = m => { const t = m.clone(); try { C.relayout(t, BL); } catch (e) { } return t; };
const clr = (m, ids) => ids.forEach(id => { const a = m.atoms.get(id); if (a) a.hFix = null; });
function distMap(m, a) { const d = new Map([[a, 0]]), q = [a]; while (q.length) { const x = q.shift(); for (const y of m.nbrs(x)) if (!d.has(y)) { d.set(y, d.get(x) + 1); q.push(y); } } return d; }
function ringOf6(m, a) { return C.rings(m, 6).find(r => r.length === 6 && r.includes(a)); }
function addAt(t, at, smi, o = 1) { const g = C.fromSmilesMB(smi, BL); for (const a of g.atoms.values()) a.hFix = null; const first = [...g.atoms.keys()][0]; const p = P(t, at); const idm = t.merge(g, p[0] + BL, p[1]); t.addBond(at, idm.get(first), o); clr(t, [at]); return idm.get(first); }
/* choose a Kekulé structure of the aromatic system in which every bond in `pairs` is double (fused rings allowed) */
function kek(m, pairs) {
  const arom = new Set(M(m, 'c').map(x => x[0])); if (!pairs.every(([a, b]) => arom.has(a) && arom.has(b) && m.bondBetween(a, b))) return false;
  const R6 = C.rings(m, 6).filter(r => r.length === 6 && r.every(id => arom.has(id)));
  const sys = new Set(); const stack = [pairs[0][0]]; while (stack.length) { const x = stack.pop(); if (sys.has(x)) continue; sys.add(x); for (const r of R6) if (r.includes(x)) for (const y of r) if (!sys.has(y)) stack.push(y); }
  const inRing = (a, b) => R6.some(r => { const i = r.indexOf(a), j = r.indexOf(b); return i >= 0 && j >= 0 && (Math.abs(i - j) === 1 || Math.abs(i - j) === 5); });
  const bonds = [...m.bonds.values()].filter(bd => sys.has(bd.a) && sys.has(bd.b) && inRing(bd.a, bd.b));
  const atoms = [...sys]; const adj = new Map(atoms.map(a => [a, []])); for (const bd of bonds) { adj.get(bd.a).push(bd); adj.get(bd.b).push(bd); }
  const need = new Set(atoms.filter(a => m.bondsOf(a).every(bd => bd.o !== 2 || bonds.includes(bd)))); // atoms whose double bond is inside the system
  const used = new Map(); const pick = new Set();
  const take = bd => { pick.add(bd); used.set(bd.a, bd); used.set(bd.b, bd); }; const drop = bd => { pick.delete(bd); used.delete(bd.a); used.delete(bd.b); };
  for (const [a, b] of pairs) { const bd = m.bondBetween(a, b); if (!bonds.includes(bd) || used.has(a) || used.has(b)) return false; take(bd); }
  const order = atoms.filter(a => need.has(a));
  const solve = i => { while (i < order.length && used.has(order[i])) i++; if (i >= order.length) return true; const a = order[i];
    for (const bd of adj.get(a)) { const o = bd.a === a ? bd.b : bd.a; if (used.has(o) || !need.has(o)) continue; take(bd); if (solve(i + 1)) return true; drop(bd); } return false; };
  if (!solve(0)) return false;
  for (const bd of bonds) bd.o = pick.has(bd) ? 2 : 1; return true;
}
const out = (mech, pre) => ({ steps: (pre || []).concat(mech.steps), product: mech.product() });
/* H⁺ onto atom `at` (from H₃O⁺); the water is removed afterwards */
function prot(mech, at, text, cond) { const ids = mech.add('[H][OH2+]', at, { dist: 1.5 * BL }); const h = ids.get(1), w = ids.get(2); mech.step([[[at], [at, h]], [[h, w], [w]]], text, cond); mech.dropFrag(w); return h; }
/* a base takes H⁺ from atom `at`; the pair stays on `at` (or forms a π bond to `to`) */
function deprot(mech, at, text, base, to, cond) { const m = mech.m; const h = mech.H(at); if (h == null) return false; const ph = P(m, h), pa = P(m, at);
  const w = mech.add(base || 'O', h, { dist: 1.45 * BL, dir: unit([ph[0] - pa[0], ph[1] - pa[1]]) }).get(1);
  mech.step([[[w], [w, h]], [[h, at], to != null ? [at, to] : [at]]], text, cond); mech.dropFrag(w); return true; }
/* a step drawn without arrows (tautomerism, work-up); fn edits a copy of the current structure */
function settle(mech, text, fn, cond) { const L = mech.m.clone(); const t = mech.m.clone(); fn(t); mech.steps.push({ L, arrows: [], R: t.clone(), text, cond: cond || [] }); mech.m = t; }
const dropH = (t, at) => { for (const h of t.nbrs(at).filter(n => el(t, n) === 'H')) t.removeAtom(h); };
const setO = (t, a, b, o) => { const x = t.bondBetween(a, b); if (x) x.o = o; clr(t, [a, b]); };
function groupPost(mech, g, from, to, ringCase, anchor) { return () => { if (ringCase) { const ids = E.fragIds(mech.m, anchor); const sub = mech.m.sub(ids); try { C.relayout(sub, BL); } catch (e) { } for (const a of mech.m.atoms.values()) { const q = sub.atoms.get(a.id); if (q) { a.x = q.x; a.y = q.y; } } } else { try { E.moveGroupTo(mech.m, g, from, to); } catch (e) { } E.tidyLongBonds(mech.m); } }; }

/* ---------------- alcohols and ethers ---------------- */
wrap('ether_413', m => {
  const al = M(m, '[CH2X4:1][OX2H1:2]').filter(([cc]) => !E.aromSet(m).has(cc)); if (!al.length) return null; const [cc, o] = al[0];
  const mech = new Mech(m);
  prot(mech, o, 'H⁺ from conc. H₂SO₄ protonates the –OH: –OH₂⁺ is a good leaving group.', ['conc. H₂SO₄']);
  const o2 = mech.add(m.clone(), cc, { att: o, dist: 1.75 * BL }).get(o);
  mech.step([[[o2], [o2, cc]], [[cc, o], [o]]], 'A second alcohol molecule attacks that carbon from the back side as water leaves (SN2).', ['413 K']); mech.dropFrag(o);
  deprot(mech, o2, 'HSO₄⁻ (or another alcohol) removes the H⁺: the ether.');
  return out(mech);
});
wrap('epox_roh_acid', (m, o, c) => {
  const ep = M(m, '[C:1]1[O:2][C:3]1')[0]; if (!ep) return null; let [a, ox, b] = ep; if (E.cN(m, b).length > E.cN(m, a).length) [a, b] = [b, a];
  const R = (c.get('ROH') || {}).smi || 'C'; const mech = new Mech(m);
  prot(mech, ox, 'H⁺ protonates the epoxide oxygen.', ['H⁺']);
  const n = mech.add('O' + R, a, { dist: 1.65 * BL }).get(1);
  mech.step([[[n], [n, a]], [[a, ox], [ox]]], 'The alcohol attacks the more substituted carbon (its C–O bond is weaker and it carries more δ+) from the back side as the ring opens: trans product.');
  deprot(mech, n, 'Loss of H⁺: the neutral alkoxy alcohol.');
  return out(mech);
});
wrap('diol_cyclic_ether', m => {
  const oh = M(m, '[CX4:1][OX2H1:2]');
  for (let i = 0; i < oh.length; i++) for (let j = i + 1; j < oh.length; j++) { const d = distMap(m, oh[i][1]).get(oh[j][0]);
    if (d === 4 || d === 5) { const [, oi] = oh[i], [cj, oj] = oh[j]; const mech = new Mech(m);
      prot(mech, oj, 'H⁺ protonates one –OH: it can now leave as water.', ['H⁺']);
      mech.step([[[oi], [oi, cj]], [[cj, oj], [oj]]], `The other –OH of the same molecule attacks that carbon from the back as water leaves (intramolecular SN2), closing a ${d + 1}-membered ring.`); mech.dropFrag(oj);
      deprot(mech, oi, 'Loss of H⁺: the cyclic ether.');
      return out(mech); } }
  return null;
});
wrap('intra_ether', (m, o) => {
  const AR = E.aromSet(m); let hit = null;
  for (const [a, b] of M(m, '[C:1]=[C:2]').filter(([a, b]) => !AR.has(a) && !AR.has(b))) for (const [cp, hC] of [[a, b], [b, a]]) { const ox = o.product.nbrs(cp).find(n => el(o.product, n) === 'O' && !m.bondBetween(cp, n)); if (ox != null) hit = { a, b, cp, hC, ox }; }
  if (!hit) return null; const { a, b, cp, hC, ox } = hit; const mech = new Mech(m);
  const ids = mech.add('[H][OH2+]', hC, { dist: 1.65 * BL }); const h = ids.get(1), w = ids.get(2);
  mech.step([[[a, b], [hC, h]], [[h, w], [w]]], 'The C=C π electrons take H⁺ so that the more stable carbocation forms (Markovnikov).', ['H⁺']); mech.dropFrag(w);
  mech.step([[[ox], [ox, cp]]], 'The molecule\'s own –OH is close enough to capture the cation: the ring closes.');
  deprot(mech, ox, 'Loss of H⁺: the cyclic ether.');
  return out(mech);
});
wrap('cumene', (m, o) => {
  if (!M(o.product, 'c[OX2H1]').length) return null;
  const cu = M(m, '[CH1X4:1]([CH3:2])([CH3:3])[c:4]')[0]; if (!cu) return null; const [ch, , , ar] = cu;
  const t0 = m.clone(); const oi = addAt(t0, ch, 'OO'); const ot = t0.nbrs(oi).find(n => n !== ch && el(t0, n) === 'O'); const L0 = relay(t0);
  const pre = [{ L: m.clone(), arrows: [], R: L0.clone(), text: 'Air (O₂) attacks the benzylic C–H by a radical chain: cumene hydroperoxide.', cond: ['O₂'] }];
  const mech = new Mech(L0);
  prot(mech, ot, 'Acid protonates the terminal –OH of the hydroperoxide.', ['H⁺']);
  mech.step([[[ch, ar], [ar, oi]], [[oi, ot], [ot]]], 'The phenyl group moves from carbon to oxygen as water leaves (a 1,2-shift onto O). Phenyl migrates in preference to methyl.', [], groupPost(mech, ar, ch, oi)); mech.dropFrag(ot);
  const w = mech.add('O', ch, { dist: 1.6 * BL }).get(1);
  mech.step([[[w], [w, ch]]], 'Water attacks the carbocation (stabilised by the O next to it).');
  deprot(mech, w, 'Loss of H⁺: a hemiketal.');
  const h = mech.H(w);
  mech.step([[[w], [w, ch]], [[ch, oi], [oi]]], 'The hemiketal falls apart: the –OH lone pair forms C=O as phenoxide leaves.');
  mech.step([[[oi], [oi, h]], [[h, w], [w]]], 'Proton transfer: phenol + acetone.');
  return out(mech, pre);
});

/* ---------------- carbonyl compounds ---------------- */
function carbonylOf(m) { return M(m, '[CX3:1]=[OX1:2]').find(([cc, ox]) => m.nbrs(cc).filter(n => n !== ox && el(m, n) !== 'H').every(n => el(m, n) === 'C')); }
wrap('acetal', (m, o, c) => {
  const k = carbonylOf(m); if (!k) return null; const [cc, ox] = k; const gly = c.has('glycol'); const R = gly ? null : ((c.get('ROH') || {}).smi || 'C');
  const mech = new Mech(m);
  prot(mech, ox, 'Acid protonates the C=O oxygen, which makes the carbonyl carbon much more electrophilic.', [gly ? 'H⁺' : 'dry HCl']);
  const ids = mech.add(gly ? 'OCCO' : 'O' + R, cc, { dist: 1.65 * BL }); const n1 = ids.get(1), n2 = gly ? ids.get(4) : null;
  mech.step([[[n1], [n1, cc]], [[cc, ox], [ox]]], 'An alcohol oxygen attacks the carbonyl carbon; the C=O π pair moves onto oxygen.');
  deprot(mech, n1, 'Loss of H⁺: the hemiacetal.');
  prot(mech, ox, 'The –OH is protonated: –OH₂⁺ can leave as water.');
  mech.step([[[n1], [n1, cc]], [[cc, ox], [ox]]], 'The OR oxygen pushes a lone pair in as water leaves: a resonance-stabilised oxocarbenium ion.'); mech.dropFrag(ox);
  const n3 = gly ? n2 : mech.add('O' + R, cc, { dist: 1.65 * BL }).get(1);
  mech.step([[[n3], [n3, cc]], [[cc, n1], [n1]]], gly ? 'The second –OH of the same glycol attacks: the five-membered ring (1,3-dioxolane) closes.' : 'A second alcohol molecule attacks the oxocarbenium carbon.');
  deprot(mech, n3, `Loss of H⁺: the ${gly ? 'cyclic ' : ''}${hTot(m, cc) ? 'acetal' : 'ketal'}. Every step is reversible: aqueous acid gives the carbonyl back.`);
  return out(mech);
});
wrap('baeyer_villiger', (m, o) => {
  const k = carbonylOf(m); if (!k) return null; const [cc, ox] = k; if (hTot(m, cc) > 0) return null;
  const mig = m.nbrs(cc).find(n => n !== ox && el(m, n) !== 'H' && !o.product.bondBetween(cc, n)); if (mig == null) return null;
  const ring = !!(E.ringSizeOf && E.ringSizeOf(m, cc, mig)); const mech = new Mech(m);
  const sp = E.species('OOC(C)=O'); const o1l = [...sp.atoms.values()].find(a => a.el === 'O' && sp.nbrs(a.id).length === 1 && el(sp, sp.nbrs(a.id)[0]) === 'O').id; const o2l = sp.nbrs(o1l)[0];
  const ids = mech.add(sp, cc, { att: o1l, dist: 1.7 * BL }); const o1 = ids.get(o1l), o2 = ids.get(o2l);
  mech.step([[[o1], [o1, cc]], [[cc, ox], [ox]]], 'The –OOH oxygen of the peracid adds to the carbonyl carbon; the C=O π pair moves onto oxygen.', ['CH₃CO₃H']);
  const h = mech.H(o1);
  mech.step([[[ox], [ox, h]], [[h, o1], [o1]]], 'Proton transfer gives the Criegee intermediate.');
  mech.step([[[ox], [ox, cc]], [[cc, mig], [mig, o1]], [[o1, o2], [o2]]], 'The O–H lone pair re-forms C=O while the group with the higher migratory aptitude moves from C to O and the weak O–O bond breaks: carboxylate leaves.', [], groupPost(mech, mig, cc, o1, ring, cc)); mech.dropFrag(o2);
  deprot(mech, ox, ring ? 'Loss of H⁺: the lactone (the ring has grown by one O).' : 'Loss of H⁺: the ester.');
  return out(mech);
});
wrap('conj_isom', m => {
  const h = M(m, '[CX3;!$(C=O)]=[CX3;!$(C=O)][CX4;!H0][CX3]=O')[0]; if (!h) return null; const [a, b, cc] = h; const mech = new Mech(m);
  deprot(mech, cc, 'OH⁻ removes the α-H, which sits between C=C and C=O and is especially acidic.', '[OH-]', null, ['OH⁻']);
  const ids = mech.add('[H]O', a, { dist: 1.6 * BL }); const h2 = ids.get(1), w2 = ids.get(2);
  mech.step([[[cc], [cc, b]], [[b, a], [a, h2]], [[h2, w2], [w2]]], 'The delocalised anion takes H⁺ back from water at the far (γ) carbon: the C=C moves into conjugation with C=O.'); mech.dropFrag(w2);
  return out(mech);
});
wrap('meyer_schuster', m => {
  const h = M(m, '[CX4;!H0]([OX2H1])C#C')[0]; if (!h) return null; const [c1, ox, c2, c3] = h; const mech = new Mech(m);
  prot(mech, ox, 'Acid protonates the –OH.', ['H⁺']);
  mech.step([[[c1, ox], [ox]]], 'Water leaves: a propargylic carbocation.'); mech.dropFrag(ox);
  mech.step([[[c2, c3], [c1, c2]]], 'Resonance: one π pair of the alkyne moves over, putting the + charge on the far carbon (an allenyl cation).');
  const w = mech.add('O', c3, { dist: 1.6 * BL }).get(1);
  mech.step([[[w], [w, c3]]], 'Water attacks that carbon.');
  deprot(mech, w, 'Loss of H⁺: an allenol (the enol of the product).');
  settle(mech, 'Keto–enol tautomerism: the allenol becomes the conjugated ketone.', t => { setO(t, c2, c3, 1); setO(t, c3, w, 2); dropH(t, w); });
  return out(mech);
});

/* ---------------- alkenes, alkynes ---------------- */
wrap('enyne_hx', (m, o, c) => {
  const q = M(m, '[CX3]=[CX3][CX2]#[CX2;H1]')[0]; if (!q) return null; const [, , c3, c4] = q; const X = c.has('HBr') ? 'Br' : 'Cl'; const mech = new Mech(m);
  const ids = mech.add('[H]' + X, c4, { dist: 1.65 * BL }); const h = ids.get(1), x = ids.get(2);
  mech.step([[[c3, c4], [c4, h]], [[h, x], [x]]], 'The triple bond takes H⁺ at its terminal carbon. The cation lands next to the C=C, where conjugation stabilises it.', [`H${X}`]);
  mech.step([[[x], [x, c3]]], `${X}⁻ captures the cation: a conjugated 2-halo-1,3-diene.`);
  return out(mech);
});
wrap('halolactone', (m, o, c) => {
  const acid = M(m, '[CX3](=O)[OX2H1,OX1-]')[0]; if (!acid) return null; const oa = acid[2]; const X = c.has('I2') ? 'I' : 'Br'; const AR = E.aromSet(m);
  let near = null, far = null; for (const [a, b] of M(m, '[C:1]=[C:2]').filter(([a, b]) => !AR.has(a) && !AR.has(b))) { if (o.product.bondBetween(oa, a)) { near = a; far = b; } else if (o.product.bondBetween(oa, b)) { near = b; far = a; } }
  if (near == null) return null; const mech = new Mech(m); const nm = X === 'I' ? 'iodonium' : 'bromonium';
  if (m.atoms.get(oa).chg !== -1) deprot(mech, oa, 'HCO₃⁻ removes the acid proton: the carboxylate.', '[O-]C(=O)O', null, ['NaHCO₃']);
  const ids = mech.add(X + X, far, { dist: 1.75 * BL }); const x1 = ids.get(1), x2 = ids.get(2);
  mech.step([[[near, far], [far, x1]], [[x1, x2], [x2]], [[x1], [x1, near]]], `The C=C π electrons attack ${X}₂; ${X}⁻ leaves and a lone pair on ${X} bridges back: a cyclic ${nm} ion.`, [`${X}₂`]); mech.dropFrag(x2);
  mech.step([[[oa], [oa, near]], [[near, x1], [x1]]], `The carboxylate of the same molecule opens the ${nm} ion from the back side, closing the lactone ring; ${X} stays on the other carbon (anti).`);
  return out(mech);
});

/* ---------------- acids and derivatives ---------------- */
function decarbSteps(mech, oa, cac, ca, ck, ok) { const h = mech.H(oa);
  mech.step([[[oa, h], [oa, cac]], [[cac, ca], [ca, ck]], [[ck, ok], [ok, h]]], 'Six-membered cyclic transition state: the C=O oxygen β to the acid takes the acid H, CO₂ leaves, and an enol forms. All in one step.', ['Δ']);
  mech.dropFrag(cac);
  settle(mech, 'Keto–enol tautomerism: the enol becomes the carbonyl compound.', t => { setO(t, ca, ck, 1); setO(t, ck, ok, 2); dropH(t, ok); }); }
wrap('acid_heat', (m, o) => {
  if (!/^β-Carbonyl acid/.test(o.title || '')) return null;
  for (const [cac, , oa] of M(m, '[CX3:1](=O)[OX2H1:2]')) { const ca = E.cN(m, cac).find(a => E.cN(m, a).some(b => b !== cac && E.isCarbonylC(m, b)) && hTot(m, a) > 0); if (ca == null) continue;
    const ck = E.cN(m, ca).find(b => b !== cac && E.isCarbonylC(m, b)); const ok = m.nbrs(ck).find(n => el(m, n) === 'O' && m.bondBetween(ck, n).o === 2);
    const mech = new Mech(m); decarbSteps(mech, oa, cac, ca, ck, ok); return out(mech); }
  return null;
});
wrap('ketoester_decarb', m => {
  const q = M(m, '[#6:1]C(=O)[CX4:2][CX3:3](=O)[OX2:4][#6:5]')[0]; if (!q) return null; const ck = q[1], ok = q[2], ca = q[3], cac = q[4], oe = q[6], r = q[7];
  if (M(m, '[CX3](=O)[OX2][#6]').length > 1) return null;
  const t0 = m.clone(); for (const id of E.subtree(m, r, oe)) t0.removeAtom(id); clr(t0, [oe]);
  const pre = [{ L: m.clone(), arrows: [], R: t0.clone(), text: 'Acid (or base, then acid) hydrolyses the ester: the β-keto acid.', cond: ['H₃O⁺'] }];
  const mech = new Mech(t0); decarbSteps(mech, oe, cac, ca, ck, ok); return out(mech, pre);
});
wrap('claisen_cond', (m, o) => {
  const enol = (mech, a) => { const h = mech.H(a); const ph = P(mech.m, h), pa = P(mech.m, a); const b = mech.add('[O-]CC', h, { dist: 1.45 * BL, dir: unit([ph[0] - pa[0], ph[1] - pa[1]]) }).get(1);
    mech.step([[[b], [b, h]], [[a, h], [a]]], 'Ethoxide removes an α-H: the ester enolate (drawn as its carbanion form).', ['NaOEt']); mech.dropFrag(b); };
  const q = M(m, '[CX4;!$(C(C=O)C=O):1][CX3:2](=[O:5])[OX2:3][CX4:4]').filter(([a]) => hTot(m, a) >= 2); if (!q.length) return null;
  if (/Dieckmann/.test(o.title || '')) {
    for (const [a, cc] of q) for (const [c2, oc2, o2] of M(m, '[CX3:1](=[O:3])[OX2:2][CX4]')) { if (c2 === cc) continue; const d = distMap(m, a).get(c2); if (d !== 4 && d !== 5) continue;
      const mech = new Mech(m); enol(mech, a);
      mech.step([[[a], [a, c2]], [[c2, oc2], [oc2]]], 'The enolate attacks the other ester C=O of the same molecule: the ring closes (tetrahedral alkoxide).');
      mech.step([[[oc2], [oc2, c2]], [[c2, o2], [o2]]], 'The C=O re-forms and the alkoxide leaves: the cyclic β-keto ester.'); mech.dropFrag(o2);
      return out(mech); }
    return null; }
  const [a, cc, oc, oe] = q[0]; const mech = new Mech(m); enol(mech, a);
  const ids = mech.add(m.clone(), a, { att: cc, dist: 1.75 * BL }); const c2 = ids.get(cc), o2 = ids.get(oc), oe2 = ids.get(oe);
  mech.step([[[a], [a, c2]], [[c2, o2], [o2]]], 'The enolate carbon attacks the C=O of a second ester molecule: a tetrahedral alkoxide.');
  mech.step([[[o2], [o2, c2]], [[c2, oe2], [oe2]]], 'The C=O re-forms and ethoxide leaves (acyl substitution): the β-keto ester. Its very acidic CH₂ is then deprotonated by EtO⁻, which drives the equilibrium; acid work-up returns it.'); mech.dropFrag(oe2);
  return out(mech);
});
wrap('nbromoamide', m => {
  const q = M(m, '[#6:1][CX3:2](=O)[NX3H1:4][Br,Cl:5]')[0]; if (!q) return null; const [r, cc, , n, x] = q; const mech = new Mech(m);
  deprot(mech, n, 'OH⁻ removes the N–H proton (made acidic by the C=O and the Br on N).', '[OH-]', null, ['OH⁻']);
  mech.step([[[cc, r], [r, n]], [[n, x], [x]], [[n], [n, cc]]], 'In one step the alkyl group moves from C to N as Br⁻ leaves, and the N lone pair forms C=N: an isocyanate, R–N=C=O. The migrating group keeps its configuration.', [], groupPost(mech, r, cc, n)); mech.dropFrag(x);
  return out(mech);
});

/* ---------------- phenols and arenes ---------------- */
/* phenoxide O⁻ pushes ring electrons onto `at` (ortho or para to the O) which bonds to electrophile atom `e`; `eArrow` is the electrophile's own arrow */
function phenoxideAttack(mech, oO, ca, at, e, eArrow, text, cond) {
  const m = mech.m; const ring = ringOf6(m, ca); if (!ring) return false; const i = ring.indexOf(ca);
  const nb = [ring[(i + 1) % 6], ring[(i + 5) % 6]];
  if (nb.includes(at)) { if (!kek(m, [[ca, at]]) && !makeDouble(m, ca, at)) return false; mech.step([[[oO], [oO, ca]], [[ca, at], [at, e]], eArrow], text, cond); return { kind: 'o', c2: at }; }
  const j = (i + 3) % 6; if (ring[j] !== at) return false; const c2 = ring[(i + 1) % 6], c3 = ring[(i + 2) % 6];
  if (!kek(m, [[ca, c2], [c3, at]]) && (!makeDouble(m, ca, c2) || m.bondBetween(c3, at).o !== 2)) return false;
  mech.step([[[oO], [oO, ca]], [[ca, c2], [c2, c3]], [[c3, at], [at, e]], eArrow], text, cond); return { kind: 'p', c2, c3 };
}
function rearomatise(mech, oO, ca, at, info, text, extra) {
  settle(mech, text, t => { setO(t, oO, ca, 1); t.atoms.get(oO).chg = 0; setO(t, ca, info.c2, 2); if (info.kind === 'p') { setO(t, info.c2, info.c3, 1); setO(t, info.c3, at, 2); } dropH(t, at); clr(t, [at]); if (extra) extra(t); }, null);
}
wrap('fries', (m, o) => {
  const q = M(m, '[c:1][O:2][C:3](=[O:5])[#6:4]')[0]; if (!q) return null; const [ca, ob, cc, oc] = q;
  const at = o.product.nbrs(cc).find(n => n !== oc && ringOf6(o.product, n) && C.aromaticRings(o.product).some(r => r.includes(n))); if (at == null) return null;
  const mech = new Mech(m);
  mech.step([[[oc], [oc, cc]], [[cc, ob], [ob]]], 'AlCl₃ (bound to the ester oxygens) pulls the acyl group off: an acylium ion R–C≡O⁺ and the aluminium phenoxide.', ['AlCl₃']);
  E.moveFragNear(mech, cc, at, 1.6 * BL);
  const info = phenoxideAttack(mech, ob, ca, at, cc, [[cc, oc], [oc]], 'The phenoxide ring attacks the acylium carbon (Friedel–Crafts acylation, inside the same molecule\'s ion pair), ortho or para to the O.');
  if (!info) return null;
  rearomatise(mech, ob, ca, at, info, 'Loss of H⁺ from the sp³ carbon and H⁺ on the oxygen (work-up) restore the aromatic ring: the hydroxyaryl ketone.');
  return out(mech);
});
wrap('kolbe2', (m, o, c) => {
  const q = M(m, '[c:1][OX1-,OX2H1:2]')[0]; if (!q) return null; const [ca, ox] = q; const acid = c.any(['H3O+', 'HCl', 'H2SO4']);
  const at = o.product.nbrs(ca).map(n => n).concat([]).find(n => o.product.nbrs(n).some(z => el(o.product, z) === 'C' && !m.atoms.has(z))); if (at == null) return null;
  const mech = new Mech(m);
  if (m.atoms.get(ox).chg !== -1) deprot(mech, ox, 'NaOH turns phenol into sodium phenoxide (much more reactive).', '[OH-]', null, ['NaOH']);
  const ids = mech.add('O=C=O', at, { att: 2, dist: 1.75 * BL }); const oA = ids.get(1), cC = ids.get(2);
  const info = phenoxideAttack(mech, ox, ca, at, cC, [[cC, oA], [oA]], 'The phenoxide ring is electron-rich at the ortho carbon; it attacks the carbon of CO₂ (a weak electrophile) at about 400 K and 4–7 atm.', ['CO₂', '400 K']);
  if (!info) return null;
  rearomatise(mech, ox, ca, at, info, acid ? 'Tautomerism re-aromatises the ring, then acid work-up gives salicylic acid.' : 'Tautomerism: the ortho H moves to the oxygen and the ring is aromatic again: sodium salicylate.', t => { if (acid) t.atoms.get(oA).chg = 0; });
  return out(mech);
});
wrap('azo_on_phenol', (m, o) => {
  const q = M(m, 'c[OX2H1,OX1-]')[0]; if (!q) return null; const [ca, ox] = q;
  const at = [...m.atoms.keys()].find(id => el(m, id) === 'C' && o.product.nbrs(id).some(n => el(o.product, n) === 'N' && !m.atoms.has(n))); if (at == null) return null;
  const mech = new Mech(m);
  if (m.atoms.get(ox).chg !== -1) deprot(mech, ox, 'In the mildly alkaline solution the phenol is present as phenoxide.', '[OH-]', null, ['OH⁻']);
  const sp = E.species('N#[N+]c1ccccc1'); const nTl = [...sp.atoms.values()].find(a => a.el === 'N' && !a.chg).id;
  const ids = mech.add(sp, at, { att: nTl, dist: 1.75 * BL }); const nT = ids.get(nTl), nI = [...ids.values()].find(id => el(mech.m, id) === 'N' && mech.m.atoms.get(id).chg === 1);
  const info = phenoxideAttack(mech, ox, ca, at, nT, [[nT, nI], [nI]], 'The diazonium ion is a weak electrophile, so it needs the very electron-rich phenoxide ring; its terminal N bonds to the ring carbon (para, or C1 of 2-naphthol).', ['0–5 °C']);
  if (!info) return null;
  rearomatise(mech, ox, ca, at, info, 'Loss of H⁺ from the ring carbon and H⁺ on the oxygen restore the aromatic ring: the azo dye.');
  return out(mech);
});
root.EngineR8 = { prot, deprot, settle, phenoxideAttack, rearomatise, kek, groupPost, setO, dropH, clr, addAt };
})(typeof window !== 'undefined' ? window : globalThis);
