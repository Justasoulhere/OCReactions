/* engine_r7.js: curly-arrow mechanisms for recipes that only had a before/after picture (Oct 13).
   Each wrapper rebuilds the reaction as arrow steps on the same atoms. If the arrow version disagrees with the
   plain recipe in any way (different product, invalid structure, exception) the plain version is kept. */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore;
const { BL, el, P, unit, hTot, Mech } = E;
const M = (m, sm) => { try { return C.matchAll(m, sm, {}); } catch (e) { return []; } };
const relay = m => { const t = m.clone(); try { C.relayout(t, BL); } catch (e) { } return t; };
const HALS = ['F', 'Cl', 'Br', 'I'];
function key(m) { try { const c = C.canon(m); const mol = C.RD.get_mol(c.replace(/[\/\\]/g, '').replace(/@+/g, '')); const k = mol.get_smiles(); mol.delete(); return k; } catch (e) { return null; } }
function validSteps(o) { return (o.steps || []).every(s => { try { return C.isValid(s.L) && (!s.R || C.isValid(s.R)); } catch (e) { return false; } }); }
function wrap(id, build) {
  const r = E.REC.find(x => x.id === id); if (!r) return; const run0 = r.run;
  r.run = function (m, c, all, k) {
    const o = run0.call(this, m, c, all, k); if (!o || !o.product || o.none) return o;
    const why = w => { if (root.__mechDebug) (root.__mechWhy = root.__mechWhy || {})[id] = w; };
    try { const o2 = build(m, o, c, all, k);
      if (!o2 || !o2.product) { if (!(o.steps || []).some(st => (st.arrows || []).length)) why('builder returned nothing'); }
      else if (key(o2.product) !== key(o.product)) why('product differs: ' + key(o2.product) + ' vs ' + key(o.product));
      else if (!validSteps(o2)) why('invalid step ' + o2.steps.findIndex(s => { try { return !(C.isValid(s.L) && (!s.R || C.isValid(s.R))); } catch (e) { return true; } }));
      else return Object.assign({}, o, { steps: o2.steps, product: o2.product, arrowed: true }); }
    catch (e) { why('threw ' + e.message + ' ' + String(e.stack || '').split('\n')[1]); }
    return o;
  };
}
function ringOf6(m, a, b) { return C.rings(m, 6).find(r => r.length === 6 && r.includes(a) && (b == null || r.includes(b))); }
function isolated(m, ring) { const R = C.rings(m, 8); return ring.every(id => R.filter(r => r.includes(id)).length === 1); }
function flipRing(m, ring) {
  const bs = []; for (let i = 0; i < 6; i++) { const b = m.bondBetween(ring[i], ring[(i + 1) % 6]); if (!b || (b.o !== 1 && b.o !== 2)) return false; bs.push(b); }
  for (let i = 0; i < 6; i++) if (bs[i].o === bs[(i + 1) % 6].o) return false;
  for (const b of bs) b.o = 3 - b.o; return true;
}
/* make the ring bond a–b a double bond by choosing the other Kekulé structure of an isolated benzene ring */
function makeDouble(m, a, b) { const bd = m.bondBetween(a, b); if (!bd) return false; if (bd.o === 2) return true; const ring = ringOf6(m, a, b); if (!ring || !isolated(m, ring) || !flipRing(m, ring)) return false; return m.bondBetween(a, b).o === 2; }
const copyXY = (to, from) => { for (const a of to.atoms.values()) { const q = from.atoms.get(a.id); if (q) { a.x = q.x; a.y = q.y; } } };
const clr = (m, ids) => ids.forEach(id => { const a = m.atoms.get(id); if (a) a.hFix = null; });

/* Diels–Alder: three arrows round the six-membered transition state */
wrap('diels_alder', (m, o, c) => {
  const D = root.EngineR6.dienes(m).filter(q => !m.bondBetween(q[0], q[3])); if (!D.length) return null; const [c1, c2, c3, c4] = D[0];
  const own = new Set(E.fragIds(m, c1)); const AR = E.aromSet(m); const t = m.clone();
  let d = M(t, '[C:1]=[C:2]').find(([a]) => !own.has(a) && !AR.has(a)) || M(t, '[C:1]#[C:2]').find(([a]) => !own.has(a));
  if (!d) { const sp = ['ethene', 'ethyne', 'acrolein', 'maleicanhydride', 'acrylate', 'acrylonitrile', 'quinone'].find(k => c.has(k)); if (!sp) return null;
    const g = C.fromSmilesMB(c.get(sp).smi, BL); const idm = t.merge(g, 0, 0); const dm = M(g, '[C:1]=[C:2]')[0] || M(g, '[C:1]#[C:2]')[0]; if (!dm) return null; d = [idm.get(dm[0]), idm.get(dm[1])]; }
  const [d1, d2] = d; const tri = t.bondBetween(d1, d2).o === 3; const six = [c1, c2, c3, c4, d1, d2];
  clr(t, six); const p = t.clone(); p.bondBetween(c1, c2).o = 1; p.bondBetween(c3, c4).o = 1; p.bondBetween(c2, c3).o = 2; p.bondBetween(d1, d2).o = tri ? 2 : 1; p.addBond(c1, d1, 1); p.addBond(c4, d2, 1);
  const R = relay(p); const L0 = t.clone(); copyXY(L0, R);
  const dien = E.fragIds(L0, d1); if (dien.includes(c1)) return null;
  const cen = ids => { let x = 0, y = 0; for (const id of ids) { const a = L0.atoms.get(id); x += a.x; y += a.y; } return [x / ids.length, y / ids.length]; };
  const A = cen([c1, c2, c3, c4]), B = cen([d1, d2]); const u = unit([B[0] - A[0], B[1] - A[1]]);
  for (const id of dien) { const a = L0.atoms.get(id); a.x += u[0] * 0.55 * BL; a.y += u[1] * 0.55 * BL; }
  const mech = new Mech(L0);
  mech.step([[[c1, c2], [c2, c3]], [[c3, c4], [c4, d2]], [[d1, d2], [d1, c1]]], 'One concerted step through a cyclic six-electron transition state: the diene (held s-cis) and the dienophile shift three π bonds at once. Two new σ bonds form at the ends and the diene\'s central bond becomes the new C=C.', ['heat'], () => copyXY(mech.m, R));
  return { steps: mech.steps, product: mech.product() };
});

/* Claisen rearrangement: [3,3] shift, then tautomerism back to the phenol */
wrap('claisen_rearr', (m) => {
  const q = M(m, '[c:5][O:1][CH2:2][CX3:3]=[CX3:4]')[0]; if (!q) return null; const [ca, o, c2, c3, c4] = q;
  const L0 = m.clone(); const ring = ringOf6(L0, ca); if (!ring) return null;
  const orth = L0.nbrs(ca).filter(n => ring.includes(n) && hTot(L0, n) > 0); if (!orth.length) return null;
  const ot = orth[0]; if (!makeDouble(L0, ca, ot)) return null; clr(L0, [ca, o, c2, c3, c4, ot]);
  const ts = L0.clone(); ts.addBond(c4, ot, 1); copyXY(L0, relay(ts));
  const mech = new Mech(L0);
  mech.step([[[c3, c4], [c4, ot]], [[o, c2], [c2, c3]], [[ot, ca], [ca, o]]], '[3,3]-Sigmatropic shift, one concerted step through a cyclic six-atom transition state: the O–CH₂ bond breaks while the far alkene carbon bonds to the ortho carbon. The π bonds move round the ring: a cyclohexadienone.', ['heat']);
  const L1 = mech.m.clone(); const t2 = mech.m.clone(); t2.bondBetween(ca, o).o = 1; t2.bondBetween(ca, ot).o = 2;
  for (const id of [ca, ot, o]) { const a = t2.atoms.get(id); a.hFix = null; a.chg = 0; } for (const h of t2.nbrs(ot).filter(n => el(t2, n) === 'H')) t2.removeAtom(h);
  mech.steps.push({ L: L1, arrows: [], R: t2.clone(), text: 'Tautomerism: the ortho H moves to the oxygen and the ring becomes aromatic again (o-allylphenol).', cond: [] });
  mech.m = t2; return { steps: mech.steps, product: mech.product() };
});

/* benzyne: elimination, addition, protonation */
wrap('benzyne', (m) => {
  const AX = M(m, '[c:1][Cl,Br,I:2]')[0]; if (!AX) return null; const [ca, x] = AX;
  const L0 = m.clone(); const ring = ringOf6(L0, ca); if (!ring) return null;
  const oh = L0.nbrs(ca).filter(n => ring.includes(n) && hTot(L0, n) > 0); if (!oh.length) return null;
  let o1 = oh.find(n => L0.bondBetween(ca, n).o === 2); if (o1 == null) { if (!makeDouble(L0, ca, oh[0])) return null; o1 = oh[0]; }
  const mech = new Mech(L0); const h = mech.H(o1); const ph = P(mech.m, h), po = P(mech.m, o1);
  const b = mech.add('[NH2-]', h, { dist: 1.45 * BL, dir: unit([ph[0] - po[0], ph[1] - po[1]]) }).get(1);
  mech.step([[[b], [b, h]], [[o1, h], [o1, ca]], [[ca, x], [x]]], 'NH₂⁻ removes the H next to the halogen; those C–H electrons form a second π bond as X⁻ leaves: benzyne (a strained, very reactive "triple bond" in the ring).', ['NaNH₂', 'liq. NH₃']);
  mech.dropFrag(b); mech.dropFrag(x);
  const n2 = mech.add('[NH2-]', ca, { dist: 1.6 * BL }).get(1);
  mech.step([[[n2], [n2, ca]], [[ca, o1], [o1]]], 'NH₂⁻ adds to one end of the benzyne; the other carbon keeps the electron pair: an aryl carbanion.');
  const ids = mech.add('[H]N', o1, { dist: 1.5 * BL }); const hh = ids.get(1), nn = ids.get(2);
  mech.step([[[o1], [o1, hh]], [[hh, nn], [nn]]], 'The aryl anion takes H⁺ from NH₃ (NH₂⁻ is regenerated).');
  mech.dropFrag(nn);
  return { steps: mech.steps, product: mech.product() };
});

/* SNAr: addition (Meisenheimer complex), then elimination */
wrap('aryl_snar', (m, o, c) => {
  const nuSpec = c.has('NaOMe') ? '[O-]C' : c.has('NaOEt') ? '[O-]CC' : c.has('NH3') ? null : '[OH-]'; if (!nuSpec) return null;
  const x = [...m.atoms.keys()].find(id => HALS.includes(el(m, id)) && !o.product.atoms.has(id)); if (x == null) return null; const ca = m.nbrs(x)[0];
  const L0 = m.clone(); const ring = ringOf6(L0, ca); if (!ring) return null; const nb = L0.nbrs(ca).filter(n => ring.includes(n));
  const nitroOn = id => L0.nbrs(id).some(s => el(L0, s) === 'N' && L0.atoms.get(s).chg > 0);
  let o1 = nb.find(n => L0.bondBetween(ca, n).o === 2 && nitroOn(n)) || nb.find(n => L0.bondBetween(ca, n).o === 2);
  if (o1 == null) { if (!makeDouble(L0, ca, nb[0])) return null; o1 = nb[0]; }
  const nNO2 = ring.filter(r => nitroOn(r)).length;
  const mech = new Mech(L0); const nu = mech.add(nuSpec, ca, { dist: 1.6 * BL }).get(1);
  mech.step([[[nu], [nu, ca]], [[ca, o1], [o1]]], `The nucleophile adds to the carbon carrying the halogen; a ring π pair is pushed onto the neighbouring carbon: the σ-complex (Meisenheimer complex).${nNO2 ? ' Its negative charge spreads onto the ortho/para –NO₂ oxygens, which is why those groups speed the reaction up.' : ' With no –NO₂ to hold the charge this needs about 623 K and 300 atm.'}`);
  mech.step([[[o1], [o1, ca]], [[ca, x], [x]]], 'The pair returns into the ring and the halide leaves: the ring is aromatic again (addition–elimination).');
  mech.dropFrag(x);
  return { steps: mech.steps, product: mech.product() };
});

/* Beckmann: protonation, concerted migration as water leaves, water addition, tautomerism */
wrap('beckmann', (m, o) => {
  const h = M(m, '[CX3]([#6])([#6])=N[OX2H1]')[0]; if (!h) return null; const [cc, r1, r2, n, ox] = h;
  const mig = [r1, r2].find(r => o.product.atoms.has(r) && o.product.bondBetween(r, n)); if (mig == null) return null;
  const ring = !!(E.ringSizeOf && E.ringSizeOf(m, cc, mig));
  const mech = new Mech(m);
  const ids = mech.add('[H][OH2+]', ox, { dist: 1.5 * BL }); const hp = ids.get(1), w0 = ids.get(2);
  mech.step([[[ox], [ox, hp]], [[hp, w0], [w0]]], 'Acid protonates the oxime –OH: –OH₂⁺ is a good leaving group.', ['H⁺']); mech.dropFrag(w0);
  mech.step([[[cc, mig], [mig, n]], [[n, ox], [ox]], [[n], [n, cc]]], 'The group anti to the leaving group moves from C to N exactly as water leaves (one step): a nitrilium ion, R–C≡N⁺–R′.', [],
    () => { if (ring) { const keep = new Set(E.fragIds(mech.m, cc)); const sub = mech.m.sub([...keep]); try { C.relayout(sub, BL); } catch (e) { } copyXY(mech.m, sub); } else { try { E.moveGroupTo(mech.m, mig, cc, n); } catch (e) { } E.tidyLongBonds(mech.m); } });
  mech.dropFrag(ox);
  const w = mech.add('O', cc, { dist: 1.6 * BL }).get(1);
  mech.step([[[w], [w, cc]], [[cc, n], [n]]], 'Water attacks the electrophilic carbon of the nitrilium ion; one π pair moves onto N.');
  const L1 = mech.m.clone(); const t2 = mech.m.clone(); t2.bondBetween(cc, n).o = 1; t2.bondBetween(cc, w).o = 2;
  for (const id of [cc, n, w]) { const a = t2.atoms.get(id); a.chg = 0; a.hFix = null; }
  mech.steps.push({ L: L1, arrows: [], R: t2.clone(), text: 'Loss of H⁺ gives the imidic acid, which tautomerises (H from O to N): the amide.', cond: [] });
  mech.m = t2; return { steps: mech.steps, product: mech.product() };
});
root.EngineR7 = { wrap, makeDouble, key };
})(typeof window !== 'undefined' ? window : globalThis);
