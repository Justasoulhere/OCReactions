/* engine_r11.js: curly-arrow mechanisms, batch 5 (Oct 16, run 25). Same contract as r7–r10: each wrapper rebuilds the
   plain recipe's result as arrow steps on the same atoms; if the product differs or a step is invalid, the plain one stays. */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore, R7 = root.EngineR7, R8 = root.EngineR8, R10 = root.EngineR10;
if (!R7 || !R8 || !R10) return;
const { BL, el, P, unit, hTot, Mech } = E; const { wrap, key } = R7; const { deprot } = R8; const { fin, put } = R10; const { prot } = R8;
const M = (m, sm) => { try { return C.matchAll(m, sm, {}); } catch (e) { return []; } };
const hasArrows = o => (o.steps || []).some(s => (s.arrows || []).length);
const W = (id, build) => wrap(id, (m, o, c, all, k) => { if (hasArrows(o)) return null; try { const r = build(m, o, c, all, k); if (root.__DBG) console.log('r11', id, r ? 'built ' + key(r.product) + ' want ' + key(o.product) : 'null'); return r; } catch (e) { if (root.__DBG) console.log('r11 ERR', id, e.message); return null; } });
const isEl = e => a => a.el === e;
const fish = (s, t) => [s, t, 1];
const jumpTo = (mech, o, text, cond) => { const L = mech.m.clone(); const R = o.product.clone(); mech.steps.push({ L, arrows: [], R: R.clone(), text, cond: cond || [] }); mech.m = R; return { steps: mech.steps, product: R }; };

/* ================= Wurtz / Fittig: radical coupling ================= */
W('wurtz2', (m, o, c) => {
  if (!c.has('Na') || c.has('RX') || (o.extra && o.extra.length > 1)) return null;
  const fr = m.frags(); if (fr.length !== 1) return null;
  const hx = M(m, '[#6:1][Cl,Br,I:2]'); if (hx.length !== 1) return null; const [cc, x] = hx[0];
  if (key(o.product) !== key(o.product)) return null;
  const mech = new Mech(m);
  /* a second molecule of the same halide, placed to the right */
  const smi = C.canon(m); const map2 = mech.add(smi, cc, { dist: 3.2 * BL }); if (!map2) return null;
  const sub2 = mech.m.sub([...map2.values()]); const h2 = M(sub2, '[#6:1][Cl,Br,I:2]'); if (h2.length !== 1) return null; const [cc2, x2] = h2[0];
  const na1 = put(mech, '[Na]', x, isEl('Na'), { dist: 1.6 * BL }).at, na2 = put(mech, '[Na]', x2, isEl('Na'), { dist: 1.6 * BL }).at;
  mech.step([fish([na1], [cc]), [[cc, x], [x]], fish([na2], [cc2]), [[cc2, x2], [x2]]],
    'Each Na atom gives one electron to a C–X bond (single-electron transfer in dry ether). The bond breaks: X⁻ goes off with Na⁺ as NaX, and an alkyl radical R• is left.', ['Na, dry ether']);
  for (const id of [x, x2]) mech.dropFrag(id); for (const id of [na1, na2]) if (mech.m.atoms.has(id)) mech.dropFrag(id);
  mech.step([fish([cc], [cc, cc2]), fish([cc2], [cc, cc2])], 'Two radicals pair their unpaired electrons: a new C–C bond, R–R. (In the alternative picture, R• takes a second electron from Na to give R⁻, which displaces X⁻ from another R–X by SN2: same product.)');
  return fin(mech, o);
});

/* ================= soda-lime decarboxylation ================= */
W('decarbox', (m, o, c) => {
  if (!c.has('sodalime') || c.has('electrolysis')) return null;
  let q = M(m, '[CX3:1](=[OX1:2])[OX1-:3]')[0]; let acid = false;
  if (!q) { q = M(m, '[CX3:1](=[OX1:2])[OX2H1:3]')[0]; acid = true; } if (!q) return null; const [c1, o2, o3] = q;
  const ca = m.nbrs(c1).find(n => el(m, n) === 'C' || el(m, n) === 'c'); if (ca == null) return null;
  const mech = new Mech(m);
  if (acid) { const h = mech.H(o3); const b = put(mech, '[OH-]', h, isEl('O'), { dist: 1.45 * BL }).at; mech.step([[[b], [b, h]], [[h, o3], [o3]]], 'NaOH first neutralises the acid: the sodium carboxylate RCOO⁻ Na⁺.', ['NaOH']); mech.dropFrag(b); }
  mech.step([[[o3], [o3, c1]], [[c1, ca], [ca]]], 'Strong heating with soda lime (NaOH + CaO): the O⁻ lone pair pushes in, the C–C bond to the carboxyl breaks and CO₂ leaves. A carbanion R⁻ is left behind.', ['NaOH, CaO, Δ']);
  mech.dropFrag(c1);
  const w = put(mech, 'O', ca, isEl('O'), { dist: 1.6 * BL }).at; const h = mech.H(w);
  mech.step([[[ca], [ca, h]], [[h, w], [w]]], 'The carbanion takes H⁺ from water (or NaOH): the alkane R–H, with one carbon fewer than the acid. CO₂ ends up as Na₂CO₃.'); mech.dropFrag(w);
  return fin(mech, o);
});

/* ================= Br₂ / HBr + conjugated diene: 1,2 vs 1,4 ================= */
W('diene_add', (m, o, c) => {
  const d = M(m, '[C:1]=[C:2]-[C:3]=[C:4]')[0]; if (!d) return null; const [a, b, cc, dd] = d;
  const X = c.has('Br2') ? 'Br' : c.has('Cl2') ? 'Cl' : c.has('HBr') ? 'HBr' : c.has('HCl') ? 'HCl' : c.has('HI') ? 'HI' : null; if (!X) return null;
  const P0 = o.product; const p14 = P0.bondBetween(b, cc) && P0.bondBetween(b, cc).o === 2;
  const mech = new Mech(m);
  if (X === 'Br' || X === 'Cl') {
    const s = put(mech, X + X, a, isEl(X), { dist: 1.7 * BL }); const e1 = s.at; const e2 = mech.m.nbrs(e1)[0];
    mech.step([[[a, b], [a, e1]], [[e1, e2], [e2]]], `The end C=C attacks ${X}₂ (Markovnikov end, C1): ${X}⁻ leaves and the cation forms on C2, where it is allylic: shared with C4 by resonance.`, [X + '₂']);
    if (p14) { mech.step([[[cc, dd], [b, cc]]], 'Resonance: the C3=C4 π electrons shift to C2–C3, so the positive charge now sits on C4.');
      mech.step([[[e2], [e2, dd]]], `${X}⁻ attacks C4: the 1,4-adduct (thermodynamic, more substituted internal C=C; favoured at higher temperature).`); }
    else mech.step([[[e2], [e2, b]]], `${X}⁻ attacks C2, right next to the first ${X}: the 1,2-adduct (kinetic, favoured at low temperature).`);
    return fin(mech, o);
  }
  const hal = X === 'HBr' ? 'Br' : X === 'HI' ? 'I' : 'Cl'; const s = put(mech, '[H]' + hal, a, a2 => a2.el === 'H', { dist: 1.6 * BL }); const h = s.at; const xx = mech.m.nbrs(h)[0];
  mech.step([[[a, b], [a, h]], [[h, xx], [xx]]], `H⁺ adds to the end carbon C1: the allylic cation, delocalised over C2 and C4.`, [X]);
  if (p14) { mech.step([[[cc, dd], [b, cc]]], 'Resonance moves the charge to C4.'); mech.step([[[xx], [xx, dd]]], `${hal}⁻ attacks C4: 1,4-addition (thermodynamic product).`); }
  else mech.step([[[xx], [xx, b]]], `${hal}⁻ attacks C2: 1,2-addition (kinetic product).`);
  return fin(mech, o);
});

/* ================= azo coupling: diazonium + phenol / aniline ================= */
W('diazo_rx', (m, o, c) => {
  if (!c.any(['phenol', 'aniline', 'naphthol'])) return null;
  const dz = M(m, '[c:1][N+:2]#[N:3]')[0]; if (!dz) return null; const [ca, n1, n2] = dz;
  const prod = o.product; const az = M(prod, '[c][N]=[N][c]')[0]; if (!az) return null;
  const nu = c.has('aniline') ? 'Nc1ccccc1' : c.has('naphthol') ? 'Oc1ccc2ccccc2c1' : 'Oc1ccccc1';
  const mech = new Mech(m);
  const g = C.fromSmilesMB(nu, BL); for (const a of g.atoms.values()) a.hFix = null;
  /* para carbon of the donor ring (for naphthol: C1, next to the OH) */
  const dn = [...g.atoms.values()].find(a => a.el === 'N' || a.el === 'O'); const ipso = g.nbrs(dn.id)[0];
  const dist = new Map([[ipso, 0]]); const q = [ipso]; while (q.length) { const x = q.shift(); for (const y of g.nbrs(x)) if (g.atoms.get(y).el === 'C' && !dist.has(y)) { dist.set(y, dist.get(x) + 1); q.push(y); } }
  const para = c.has('naphthol') ? [...dist].find(([id, d]) => d === 1 && g.nbrs(id).some(z => z !== ipso && g.nbrs(z).filter(w => g.atoms.get(w).el === 'C').length === 3))?.[0] : [...dist].find(([, d]) => d === 3)?.[0];
  if (para == null) return null;
  const map = mech.add(g, n2, { att: para, dist: 1.8 * BL }); const pc = map.get(para), dO = map.get(dn.id), ip = map.get(ipso);
  const ring = E.fragIds(mech.m, pc); const pb = mech.m.bondsOf(pc).find(b => b.o === 2 && ring.includes(b.a === pc ? b.b : b.a)); if (!pb) return null;
  const pnb = pb.a === pc ? pb.b : pb.a;
  mech.step([[[pc, pnb], [pc, n2]], [[n2, n1], [n1]]], `The ${c.has('aniline') ? '–NH₂' : '–OH'} group makes the ring electron-rich; its ${c.has('naphthol') ? 'C1' : 'para'} carbon attacks the terminal N of the diazonium ion (a weak electrophile, so only activated rings react). N≡N becomes N=N.`, [c.has('aniline') ? 'aniline, pH 4–5' : 'NaOH, pH 9–10, 0–5 °C']);
  const h = mech.H(pc);
  if (h == null) return null;
  const b = put(mech, '[OH-]', h, isEl('O'), { dist: 1.45 * BL }).at;
  mech.step([[[b], [b, h]], [[h, pc], [pc, pnb]]], 'The base removes H⁺ from the sp³ carbon: the ring is aromatic again. The azo dye Ar–N=N–Ar′ (orange / yellow, extended conjugation).'); mech.dropFrag(b);
  return fin(mech, o);
});

/* ================= terminal alkyne + ammoniacal AgNO₃ / Cu₂Cl₂ ================= */
W('acetylide_test', (m, o, c) => {
  const t = M(m, '[CH1;$(C#[#6]):1]'); if (!t.length || t.length > 2) return null;
  const metal = c.any(['CuNH3', 'CuCl']) ? 'Cu' : 'Ag';
  const mech = new Mech(m);
  for (const [c1] of t) {
    const h = mech.H(c1); const ph = P(mech.m, h), pc = P(mech.m, c1);
    const nN = mech.add('N', h, { dist: 1.45 * BL, dir: unit([ph[0] - pc[0], ph[1] - pc[1]]) }).get(1);
    mech.step([[[nN], [nN, h]], [[h, c1], [c1]]], 'NH₃ removes the terminal H⁺: the sp carbon (50% s-character) holds the negative charge well, so the ≡C–H is acidic enough (pKa ≈ 25).', [metal === 'Ag' ? "Tollens' / ammoniacal AgNO₃" : 'ammoniacal Cu₂Cl₂']);
    mech.dropFrag(nN);
    const mt = put(mech, '[' + metal + '+]', c1, isEl(metal), { dist: 1.6 * BL }).at;
    mech.step([[[c1], [c1, mt]]], metal === 'Ag' ? 'The acetylide binds Ag⁺: silver acetylide, a white precipitate. Internal alkynes (no ≡C–H) give no precipitate.' : 'The acetylide binds Cu⁺: red-brown copper(I) acetylide.');
  }
  return fin(mech, o);
});


/* ================= batch 6 (run 25) ================= */
const hid = (mech, h, base, at) => { const ph = P(mech.m, h), pa = P(mech.m, at); return put(mech, base, h, a => a.el === 'O' && (a.chg || 0) < 0, { dist: 1.45 * BL, dir: unit([ph[0] - pa[0], ph[1] - pa[1]]) }).at; };
/* 3° ether + dilute acid: SN1 */
W('ether3_hyd', (m, o) => {
  const t = M(m, '[CX4;H0:1]([#6])([#6])([#6])[OX2:2][#6:3]')[0]; if (!t) return null; const c3 = t[0], ox = t[4];
  const mech = new Mech(m);
  prot(mech, ox, 'H⁺ protonates the ether oxygen: a good leaving group (an alcohol) is made.', ['dil. H₂SO₄', 'H₂O']);
  mech.step([[[c3, ox], [ox]]], 'The C–O bond to the 3° carbon breaks on its own (SN1): a stable 3° carbocation and the other half leaves as an alcohol. The 1°/methyl side never ionises.');
  const w = put(mech, 'O', c3, isEl('O'), { dist: 1.65 * BL }).at;
  mech.step([[[w], [w, c3]]], 'Water attacks the flat carbocation.');
  deprot(mech, w, 'Loss of H⁺: the 3° alcohol.');
  return fin(mech, o);
});
/* enol ether + H3O+: C-protonation, then hemiacetal breakdown */
W('enol_ether', (m, o) => {
  const t = M(m, '[CX3:1]=[CX3:2][OX2:3][#6:4]')[0]; if (!t) return null; const [c1, c2, o3] = t;
  const mech = new Mech(m);
  const ids = mech.add('[H][OH2+]', c1, { dist: 1.5 * BL }); const h = ids.get(1), w0 = ids.get(2);
  mech.step([[[o3], [o3, c2]], [[c1, c2], [c1, h]], [[h, w0], [w0]]], 'The O lone pair pushes into the C=C, which takes H⁺ on carbon: a resonance-stabilised oxocarbenium ion (C=O⁺R).', ['dil. H₃O⁺']); mech.dropFrag(w0);
  const w = put(mech, 'O', c2, isEl('O'), { dist: 1.65 * BL }).at;
  mech.step([[[w], [w, c2]], [[c2, o3], [o3]]], 'Water attacks the carbon of C=O⁺; the π pair returns to oxygen.');
  deprot(mech, w, 'Loss of H⁺: a hemiacetal.');
  prot(mech, o3, 'H⁺ on the –OR oxygen…');
  mech.step([[[w], [w, c2]], [[c2, o3], [o3]]], '…the –OH lone pair pushes in and the alcohol leaves: a protonated carbonyl.');
  deprot(mech, w, 'Loss of H⁺: the aldehyde / ketone and the alcohol.');
  return fin(mech, o);
});
/* Weitz–Scheffer: HOO⁻ conjugate addition, then the enolate displaces OH⁻ */
W('weitz', (m, o) => {
  const t = M(m, '[CX3:1]=[CX3:2][CX3:3]=[OX1:4]')[0]; if (!t) return null; const [cb, ca, cc, ox] = t;
  const mech = new Mech(m);
  const s = put(mech, '[O-]O', cb, a => a.el === 'O' && a.chg < 0, { dist: 1.65 * BL }); const o1 = s.at, o2 = s.ids.find(x => x !== o1);
  mech.step([[[o1], [o1, cb]], [[cb, ca], [ca, cc]], [[cc, ox], [ox]]], 'OH⁻ + H₂O₂ gives HOO⁻, a soft α-effect nucleophile: it adds to the β carbon (1,4-addition) and the charge ends on the carbonyl O (an enolate).', ['H₂O₂', 'OH⁻']);
  mech.step([[[ox], [ox, cc]], [[cc, ca], [ca, o1]], [[o1, o2], [o2]]], 'The enolate swings back: the α carbon attacks the near peroxide O and the weak O–O bond breaks, OH⁻ leaving. A three-membered ring: the α,β-epoxy ketone.');
  mech.dropFrag(o2);
  return fin(mech, o);
});
/* chloral + OH⁻: CCl₃⁻ as leaving group */
W('chloral', (m, o) => {
  const t = M(m, '[CX3H1:1](=[OX1:2])[CX4:3](Cl)(Cl)Cl')[0]; if (!t) return null; const [c1, ox, c3] = t;
  const mech = new Mech(m);
  const w = put(mech, '[OH-]', c1, isEl('O'), { dist: 1.65 * BL }).at;
  mech.step([[[w], [w, c1]], [[c1, ox], [ox]]], 'OH⁻ adds to the very electrophilic C=O (three Cl pull electrons).', ['OH⁻']);
  mech.step([[[ox], [ox, c1]], [[c1, c3], [c3]]], 'The C=O re-forms and kicks out ⁻CCl₃: three Cl atoms stabilise that carbanion enough to be a leaving group (as in the haloform reaction).');
  const h = mech.H(w); if (h == null) return null;
  mech.step([[[c3], [c3, h]], [[h, w], [w]]], 'Proton transfer: ⁻CCl₃ takes the acidic H of HCOOH. Chloroform + formate.');
  return fin(mech, o);
});
/* R–MgX + R′–X  /  R–MgX + ClCN */
W('grignard_rx', (m, o) => {
  const g = M(m, '[#6:1][Mg:2]')[0]; if (!g) return null; const [rc, mg] = g; const mech = new Mech(m);
  const s = put(mech, 'CCBr', rc, (a, gg) => a.el === 'C' && gg.bondsOf(a.id).some(b => gg.atoms.get(gg.other(b, a.id)).el === 'Br')); const cx = s.at;
  const br = mech.m.nbrs(cx).find(n => el(mech.m, n) === 'Br');
  mech.step([[[rc, mg], [rc, cx]], [[cx, br], [br]]], 'The carbanion-like carbon of R–MgX attacks the C–X carbon from the back (SN2) and X⁻ leaves (as MgX₂): R–R′.', ['dry ether']);
  mech.dropFrag(mg);
  return fin(mech, o);
});
W('grignard_clcn', (m, o) => {
  const g = M(m, '[#6:1][Mg:2]')[0]; if (!g) return null; const [rc, mg] = g; const mech = new Mech(m);
  const s = put(mech, 'ClC#N', rc, isEl('C')); const cn = s.at; const cl = mech.m.nbrs(cn).find(n => el(mech.m, n) === 'Cl'), n = mech.m.nbrs(cn).find(x => el(mech.m, x) === 'N');
  mech.step([[[rc, mg], [rc, cn]], [[cn, n], [n]]], 'R⁻ (from R–MgX) adds to the nitrile carbon of Cl–C≡N; one π pair moves to N.', ['dry ether']);
  mech.step([[[n], [n, cn]], [[cn, cl], [cl]]], 'The pair on N comes back to re-form C≡N and Cl⁻ is expelled (Cl is the better leaving group): R–C≡N, one carbon longer.');
  mech.dropFrag(mg); mech.dropFrag(cl);
  return fin(mech, o);
});
/* IN₃: iodine electrophile, azide on the more substituted carbon */
W('in3_add', (m, o) => {
  const db = M(m, '[CX3:1]=[CX3:2]')[0]; if (!db) return null; let [a, b] = db;
  const deg = x => m.nbrs(x).filter(n => el(m, n) === 'C').length; if (deg(a) > deg(b)) [a, b] = [b, a];
  const mech = new Mech(m);
  const s = put(mech, 'IN=[N+]=[N-]', a, isEl('I'), { dist: 1.7 * BL }); const io = s.at; const n1 = mech.m.nbrs(io).find(n => el(mech.m, n) === 'N');
  mech.step([[[a, b], [a, io]], [[io, n1], [n1]]], 'I–N₃ is polarised Iδ⁺–N₃δ⁻: the π bond attacks iodine at the less substituted carbon (via the iodonium ion), and N₃⁻ leaves. Positive charge sits on the more substituted carbon.', ['IN₃']);
  mech.step([[[n1], [n1, b]]], 'Azide (N₃⁻) attacks the more substituted carbon, anti to iodine: Markovnikov-type, like halohydrin formation.');
  return fin(mech, o);
});
/* exhaustive methylation */
W('amine_rx', (m, o) => {
  if (!M(o.product, '[NX4+]([#6])([#6])([#6])[#6]').length) return null;
  const am = M(m, '[NX3;!$(NC=O):1]').filter(([n]) => hTot(m, n) > 0)[0]; if (!am) return null; const n = am[0];
  const mech = new Mech(m); const k = hTot(m, n);
  for (let i = 0; i <= k; i++) {
    const s = put(mech, 'CI', n, isEl('C'), { dist: 1.45 * BL }); const cm = s.at; const io = mech.m.nbrs(cm).find(x => el(mech.m, x) === 'I');
    mech.step([[[n], [n, cm]], [[cm, io], [io]]], i === 0 ? 'The N lone pair attacks CH₃ of CH₃I (SN2), iodide leaves.' : i < k ? 'The new amine is still nucleophilic (more so: alkyl groups push electrons): another SN2 on CH₃I.' : 'Last methylation: no N–H is left, so the quaternary ammonium iodide R–N⁺(CH₃)₃ I⁻ forms and stops here.', i === 0 ? ['CH₃I (excess)'] : []);
    mech.dropFrag(io);
    if (i < k) deprot(mech, n, 'A base (excess amine) takes H⁺ from N.', 'N');
  }
  return fin(mech, o);
});
/* Hofmann exhaustive methylation + elimination */
W('hofmann_elim', (m, o) => {
  const am = M(m, '[CX4:1][NX3H2:2]')[0]; if (!am) return null; const [ca, n] = am;
  const mech = new Mech(m);
  for (let i = 0; i < 3; i++) {
    const s = put(mech, 'CI', n, isEl('C'), { dist: 1.45 * BL }); const cm = s.at; const io = mech.m.nbrs(cm).find(x => el(mech.m, x) === 'I');
    mech.step([[[n], [n, cm]], [[cm, io], [io]]], i === 0 ? 'Exhaustive methylation: the N lone pair attacks CH₃I (SN2)…' : i === 1 ? '…again…' : '…and a third time: the quaternary ammonium ion R–N⁺(CH₃)₃.', i === 0 ? ['CH₃I (excess)'] : []);
    mech.dropFrag(io); if (i < 2) deprot(mech, n, 'Base takes H⁺ from N.', 'N');
  }
  const bs = mech.m.nbrs(ca).filter(x => el(mech.m, x) === 'C' && mech.H(x) != null);
  if (!bs.length) return null;
  const hc = x => mech.m.nbrs(x).filter(y => el(mech.m, y) === 'H').length + ((mech.m.atoms.get(x).hFix) || 0);
  const cb = bs.sort((x, y) => mech.m.nbrs(x).filter(z => el(mech.m, z) === 'C').length - mech.m.nbrs(y).filter(z => el(mech.m, z) === 'C').length)[0];
  const h = mech.H(cb); const b = hid(mech, h, '[OH-]', cb);
  mech.step([[[b], [b, h]], [[h, cb], [cb, ca]], [[ca, n], [n]]], 'Ag₂O/H₂O swaps I⁻ for OH⁻; on heating OH⁻ removes a β-H from the least hindered carbon (E2, bulky ⁺NMe₃ leaving group): Hofmann product, the less substituted alkene + N(CH₃)₃.', ['Ag₂O, H₂O', 'Δ']);
  mech.dropFrag(b); mech.dropFrag(n);
  return fin(mech, o);
});
/* nitrile α-carbanion + HCHO */
W('nitrile_hcho', (m, o) => {
  const t = M(m, '[CX4;!H0:1][CX2:2]#[NX1]')[0]; if (!t) return null; const ca = t[0];
  const mech = new Mech(m); const h = mech.H(ca); const b = hid(mech, h, '[OH-]', ca);
  mech.step([[[b], [b, h]], [[h, ca], [ca]]], 'OH⁻ removes an α-H: the carbanion is stabilised by the C≡N group.', ['NaOH']); mech.dropFrag(b);
  const s = put(mech, 'C=O', ca, isEl('C'), { dist: 1.75 * BL }); const cf = s.at, of = mech.m.nbrs(cf).find(x => el(mech.m, x) === 'O');
  mech.step([[[ca], [ca, cf]], [[cf, of], [of]]], 'The carbanion adds to HCHO; the C=O π pair goes to oxygen.', ['HCHO']);
  const ids = mech.add('[H]O', of, { dist: 1.5 * BL }); const hh = ids.get(1), w = ids.get(2);
  mech.step([[[of], [of, hh]], [[hh, w], [w]]], 'The alkoxide takes H⁺ from water (OH⁻ regenerated): the β-hydroxy nitrile.'); mech.dropFrag(w);
  return fin(mech, o);
});
/* diethyl oxalate + 1° amine (twice) */
W('oxalate_amine', (m, o) => {
  const es = M(m, '[CX3:1](=[OX1:2])[OX2:3][#6]'); if (es.length !== 2) return null;
  const mech = new Mech(m);
  for (const [c1, o2, o3] of es) {
    const s = put(mech, 'CCN', c1, isEl('N'), { dist: 1.7 * BL }); const n = s.at;
    mech.step([[[n], [n, c1]], [[c1, o2], [o2]]], 'The amine N attacks an ester carbonyl carbon (tetrahedral intermediate).', ['C₂H₅NH₂']);
    mech.step([[[o2], [o2, c1]], [[c1, o3], [o3]]], 'C=O re-forms and EtO⁻ leaves.');
    deprot(mech, n, 'EtO⁻ takes H⁺ from N⁺: the amide. (Then the same on the second ester.)', '[O-]CC');
    mech.dropFrag(o3);
  }
  return fin(mech, o);
});

/* cyclopropane + Br₂: strained σ bond acts like a π bond */
W('cyclopropane', (m, o) => {
  const t = M(m, '[CX4:1]1[CX4:2][CX4:3]1')[0]; if (!t) return null; const [c1, c2] = t;
  const mech = new Mech(m);
  const s = put(mech, 'BrBr', c1, isEl('Br'), { dist: 1.45 * BL }); const b1 = s.at, b2 = s.ids.find(x => x !== b1);
  mech.step([[[c1, c2], [c1, b1]], [[b1, b2], [b2]]], 'The bent, strained C–C σ bond of cyclopropane (about 60° angles, banana bonds with π-like character) attacks Br₂; Br⁻ leaves and the ring opens to a carbocation.', ['Br₂']);
  mech.step([[[b2], [b2, c2]]], 'Br⁻ captures the cation: the 1,3-dibromide. Larger rings (cyclobutane upwards) are not strained enough to do this.');
  return fin(mech, o);
});

/* ================= batch 7 (run 26) ================= */
const PPH3 = 'P(c1ccccc1)(c1ccccc1)c1ccccc1';
const CARB2 = { HCHO: 'C=O', CH3CHO: 'CC=O', acetone: 'CC(C)=O', PhCHO: 'O=Cc1ccccc1' };
const dO = (mm, x) => mm.nbrs(x).find(n => el(mm, n) === 'O' && mm.bondBetween(x, n).o === 2);
W('pph3_epox', (m, o) => {
  const ep = M(m, '[C:1]1[O:2][C:3]1')[0]; if (!ep) return null; let [c1, ox, c3] = ep;
  const deg = x => m.nbrs(x).filter(n => el(m, n) === 'C').length; if (deg(c1) > deg(c3)) [c1, c3] = [c3, c1];
  const mech = new Mech(m);
  const s = put(mech, PPH3, c1, isEl('P'), { dist: 1.8 * BL }); const pp = s.at;
  mech.step([[[pp], [pp, c1]], [[c1, ox], [ox]]], 'The P lone pair of PPh₃ attacks the less hindered ring carbon from the back (SN2) and the strained ring opens: a betaine (P⁺ … O⁻).', ['PPh₃']);
  mech.step([[[ox], [ox, pp]]], 'O⁻ bonds to P⁺ (after rotation about the C–C bond): the four-membered oxaphosphetane, as in the Wittig reaction.');
  mech.step([[[c1, pp], [c1, c3]], [[c3, ox], [ox, pp]]], 'The ring breaks up: C=C forms and Ph₃P=O leaves (the very strong P=O bond pays for it). The SN2 inversion means the alkene geometry flips: a cis-epoxide gives the trans-alkene.');
  return fin(mech, o);
});
W('mpv', (m, o, c) => {
  const mech = new Mech(m);
  const k = M(m, '[CX3:1]=[OX1:2]').filter(([x]) => m.nbrs(x).every(n => ['C', 'H', 'O'].includes(el(m, n)) || el(m, n) === 'c'))[0];
  if (k && key(o.product) !== key(m)) {
    const [cc, ox] = k;
    const s = put(mech, '[O-]C(C)C', cc, isEl('C'), { dist: 1.9 * BL }); const ci = s.at; const oi = mech.m.nbrs(ci).find(n => el(mech.m, n) === 'O'); const h = mech.H(ci); if (h == null) return null;
    mech.step([[[oi], [oi, ci]], [[ci, h], [h, cc]], [[cc, ox], [ox]]], 'Both oxygens sit on the same Al: a six-membered cyclic transition state. The isopropoxide C–H passes its hydride to the carbonyl carbon; the isopropoxide becomes acetone (distilled off, pulling the equilibrium).', ['Al(OiPr)₃', 'iPrOH']);
    mech.dropFrag(oi); prot(mech, ox, 'Work-up: the aluminium alkoxide is protonated: the alcohol. C=C, –NO₂ and ester groups are untouched.', ['H₃O⁺']);
    return fin(mech, o);
  }
  const al = M(m, '[CX4;H1:1]([#6])([#6])[OX2H1:4]')[0]; if (!al) return null; const cc = al[0], ox = al[3];
  deprot(mech, ox, 'The alcohol exchanges onto aluminium: an aluminium alkoxide (shown as the free alkoxide).', '[O-]C(C)(C)C', null, ['Al(OtBu)₃']);
  const s = put(mech, 'CC(C)=O', cc, (a, g) => a.el === 'C' && g.bondsOf(a.id).some(b => b.o === 2), { dist: 1.9 * BL }); const ca = s.at; const oa = dO(mech.m, ca); const h = mech.H(cc); if (h == null) return null;
  mech.step([[[ox], [ox, cc]], [[cc, h], [h, ca]], [[ca, oa], [oa]]], 'Six-membered cyclic transition state on Al: the carbinol C–H moves as a hydride to acetone (large excess); the alkoxide O forms the C=O. The ketone (Oppenauer oxidation).', ['acetone']);
  mech.dropFrag(oa);
  return fin(mech, o);
});
W('reformatsky', (m, o, c) => {
  const q = M(m, '[Br,Cl,I:1][CX4:2][CX3:3](=O)[OX2:4][#6]')[0]; if (!q) return null; const [x, ca] = q;
  const k = Object.keys(CARB2).find(z => c.has(z)); if (!k) return null;
  const mech = new Mech(m); let zn = null;
  R8.settle(mech, 'Zn inserts into the C–Br bond (dry ether/benzene): the organozinc Reformatsky reagent BrZn–CH₂COOR (it exists as a zinc enolate).', t => { t.removeBond(t.bondBetween(ca, x).id); const pa = P(t, ca), px = P(t, x); zn = t.addAtom('Zn', (pa[0] + px[0]) / 2, (pa[1] + px[1]) / 2, {}); t.addBond(ca, zn, 1); t.addBond(zn, x, 1); }, ['Zn']);
  const s = put(mech, CARB2[k], ca, (a, g) => a.el === 'C' && g.bondsOf(a.id).some(b => b.o === 2 && g.atoms.get(g.other(b, a.id)).el === 'O'), { dist: 1.8 * BL }); const ck = s.at, ok = dO(mech.m, ck);
  mech.step([[[ca, zn], [ca, ck]], [[ck, ok], [ok]]], 'The C–Zn carbon adds to the aldehyde/ketone C=O; the π pair moves to O. Organozinc is too weak to attack the ester, which survives.', [k === 'acetone' ? 'CH₃COCH₃' : k]);
  mech.dropFrag(zn); prot(mech, ok, 'H₃O⁺ work-up: the β-hydroxy ester.', ['H₃O⁺']);
  return fin(mech, o);
});
W('ddt', (m, o) => {
  const cl = M(m, '[Cl:1][c:2]')[0]; if (!cl) return null; const ip = cl[1];
  const dist = new Map([[ip, 0]]); const qq = [ip]; while (qq.length) { const x2 = qq.shift(); for (const y of m.nbrs(x2)) if ((el(m, y) === 'C' || el(m, y) === 'c') && !dist.has(y)) { dist.set(y, dist.get(x2) + 1); qq.push(y); } }
  const pc = [...dist].find(([, d]) => d === 3)?.[0]; if (pc == null) return null;
  const mech = new Mech(m);
  const s = put(mech, '[OH+]=CC(Cl)(Cl)Cl', pc, (a, g) => a.el === 'C' && g.bondsOf(a.id).some(b => b.o === 2), { dist: 1.8 * BL }); const ce = s.at; const oe = dO(mech.m, ce);
  const pb = mech.m.bondsOf(pc).find(b => b.o === 2); if (!pb) return null; const pn = pb.a === pc ? pb.b : pb.a;
  mech.step([[[pc, pn], [pc, ce]], [[ce, oe], [oe]]], 'Conc. H₂SO₄ protonates chloral; the para carbon of chlorobenzene (Cl is o/p-directing, para less hindered) attacks the carbonyl carbon.', ['chloral', 'conc. H₂SO₄']);
  deprot(mech, pc, 'Loss of H⁺ (to HSO₄⁻ / water) restores the ring: a benzylic alcohol, ArCH(OH)CCl₃.', 'O', pn);
  return jumpTo(mech, o, 'The –OH is protonated and leaves as water: a benzylic cation, which attacks a second chlorobenzene para to Cl. Loss of H⁺: DDT, p,p′-dichlorodiphenyltrichloroethane.');
});
W('na_nh3_terminal', (m, o) => {
  const t = M(m, '[CH1;$(C#[#6]):1]'); if (t.length !== 1) return null; const c1 = t[0][0];
  const mech = new Mech(m); const h = mech.H(c1); const ph = P(mech.m, h), pc = P(mech.m, c1);
  const nN = mech.add('[NH2-]', h, { dist: 1.45 * BL, dir: unit([ph[0] - pc[0], ph[1] - pc[1]]) }).get(1);
  mech.step([[[nN], [nN, h]], [[h, c1], [c1]]], 'Na in liquid NH₃ (with a trace of Fe) gives NaNH₂. Amide ion (conjugate acid NH₃, pKa ≈ 38) easily removes the terminal ≡C–H (pKa ≈ 25).', ['Na, liq. NH₃']);
  mech.dropFrag(nN);
  const na = put(mech, '[Na+]', c1, isEl('Na'), { dist: 1.6 * BL }).at;
  mech.step([[[c1], [c1, na]]], 'The acetylide pairs with Na⁺: sodium acetylide, a carbon nucleophile for making longer alkynes with R–X.');
  return fin(mech, o);
});
W('rx_base', (m, o) => {
  const t = M(m, '[CH1:1](Cl)(Cl)Cl')[0]; if (!t || m.atoms.size > 6 || !M(o.product, '[CH1](=O)[O-]').length) return null; const cm = t[0];
  const mech = new Mech(m); const h = mech.H(cm); const ph = P(mech.m, h), pc = P(mech.m, cm);
  const b = mech.add('[OH-]', h, { dist: 1.45 * BL, dir: unit([ph[0] - pc[0], ph[1] - pc[1]]) }).get(1);
  mech.step([[[b], [b, h]], [[h, cm], [cm]]], 'OH⁻ removes the H of chloroform (three Cl make it acidic): ⁻CCl₃.', ['aq. KOH']); mech.dropFrag(b);
  const cl = mech.m.nbrs(cm).find(n => el(mech.m, n) === 'Cl');
  mech.step([[[cm, cl], [cl]]], 'α-Elimination: Cl⁻ leaves: dichlorocarbene :CCl₂.', [], () => { const a = mech.m.atoms.get(cm); a.chg = 0; a.rad = 2; a.hFix = 0; }); mech.dropFrag(cl);
  return jumpTo(mech, o, 'OH⁻ adds to the carbene, and each remaining C–Cl is hydrolysed in turn (Cl⁻ out, OH in); the gem-diol loses water to HCOOH, which the base turns into formate HCOO⁻.');
});
W('wurtz2', (m, o, c) => {
  const hx = M(m, '[CX4:1][Cl,Br,I:2]'); if (hx.length !== 2 || m.frags().length !== 1) return null;
  const [[ca, xa], [cb, xb]] = hx; if (!o.product.bondBetween(ca, cb) || m.bondBetween(ca, cb)) return null;
  const metal = c.has('Na') ? 'Na' : c.has('Zn') ? 'Zn' : null; if (!metal) return null;
  const mech = new Mech(m);
  if (metal === 'Na') {
    const n1 = put(mech, '[Na]', xa, isEl('Na'), { dist: 1.6 * BL }).at, n2 = put(mech, '[Na]', xb, isEl('Na'), { dist: 1.6 * BL }).at;
    mech.step([fish([n1], [ca]), [[ca, xa], [xa]], fish([n2], [cb]), [[cb, xb], [xb]]], 'Each Na atom transfers one electron to a C–X bond: both halides leave as NaX, and a radical forms at each end of the chain.', ['Na, dry ether']);
    for (const id of [xa, xb, n1, n2]) mech.dropFrag(id);
  } else {
    const zn = put(mech, '[Zn]', xa, isEl('Zn'), { dist: 1.6 * BL }).at;
    mech.step([fish([zn], [ca]), [[ca, xa], [xa]], fish([zn], [cb]), [[cb, xb], [xb]]], 'Zn gives one electron to each C–X bond; ZnX₂ forms and a radical sits at each end.', ['Zn', 'Δ']);
    for (const id of [xa, xb, zn]) mech.dropFrag(id);
  }
  mech.step([fish([ca], [ca, cb]), fish([cb], [ca, cb])], 'The two ends of the same chain meet before they can find another molecule (intramolecular, effective concentration is high): the ring closes.', [], () => { const ids = E.fragIds(mech.m, ca); const sub = mech.m.sub(ids); try { C.relayout(sub, BL); } catch (e) { } for (const a of mech.m.atoms.values()) { const q = sub.atoms.get(a.id); if (q) { a.x = q.x; a.y = q.y; } } });
  return fin(mech, o);
});

/* ================= batch 8 (run 26) ================= */
const ringPara = (g, ip) => { const d = new Map([[ip, 0]]); const q = [ip]; while (q.length) { const x = q.shift(); for (const y of g.nbrs(x)) { const e = (g.atoms.get(y) || {}).el; if ((e === 'C' || e === 'c') && !d.has(y)) { d.set(y, d.get(x) + 1); q.push(y); } } } return [...d].find(([, k]) => k === 3)?.[0]; };
W('wurtz2', (m, o, c) => {
  if (!c.has('Na')) return null;
  const mech = new Mech(m); let A, B;
  const hx = M(m, '[#6:1][Cl,Br,I:2]');
  if (m.frags().length === 2 && hx.length === 2) [A, B] = hx;
  else if (m.frags().length === 1 && hx.length === 1 && c.has('RX')) { const rx = c.get('RX'); A = hx[0]; const g = C.fromSmilesMB((rx.X || 'Br') + (rx.smi || 'C'), BL); for (const a of g.atoms.values()) a.hFix = null;
    const xi = [...g.atoms.values()].find(a => a.el === (rx.X || 'Br')); const map = mech.add(g, A[0], { att: xi.id, dist: 3 * BL }); const xn = map.get(xi.id); B = [mech.m.nbrs(xn)[0], xn]; }
  else return null;
  const [ca, xa] = A, [cb, xb] = B;
  const n1 = put(mech, '[Na]', xa, isEl('Na'), { dist: 1.6 * BL }).at, n2 = put(mech, '[Na]', xb, isEl('Na'), { dist: 1.6 * BL }).at;
  mech.step([fish([n1], [ca]), [[ca, xa], [xa]], fish([n2], [cb]), [[cb, xb], [xb]]], 'Na gives one electron to each C–X bond: NaX forms and two different radicals (aryl and alkyl) are left.', ['Na, dry ether']);
  for (const id of [xa, xb, n1, n2]) mech.dropFrag(id);
  mech.step([fish([ca], [ca, cb]), fish([cb], [ca, cb])], 'Cross-coupling: the aryl and alkyl radicals pair up (Wurtz–Fittig). Some R–R and Ar–Ar also form, so the cross product is separated by distillation.');
  return fin(mech, o);
});
W('carbylamine_rev', (m, o, c) => {
  const t = M(m, '[CH1:1](Cl)(Cl)Cl')[0]; if (!t) return null; const cm = t[0];
  const mech = new Mech(m); const h = mech.H(cm); const ph = P(mech.m, h), pc = P(mech.m, cm);
  const b = mech.add('[OH-]', h, { dist: 1.45 * BL, dir: unit([ph[0] - pc[0], ph[1] - pc[1]]) }).get(1);
  mech.step([[[b], [b, h]], [[h, cm], [cm]]], 'OH⁻ removes the H of chloroform: ⁻CCl₃.', ['KOH']); mech.dropFrag(b);
  const cl = mech.m.nbrs(cm).find(n => el(mech.m, n) === 'Cl');
  mech.step([[[cm, cl], [cl]]], 'α-Elimination: dichlorocarbene :CCl₂.', [], () => { const a = mech.m.atoms.get(cm); a.chg = 0; a.rad = 2; a.hFix = 0; }); mech.dropFrag(cl);
  const am = (c.get('RNH2') || {}).smi ? 'N' + c.get('RNH2').smi : 'NCC';
  const s = put(mech, am, cm, isEl('N'), { dist: 1.8 * BL }); const n = s.at;
  mech.step([[[n], [n, cm]]], 'The 1° amine N lone pair attacks the empty orbital of the carbene carbon.', ['RNH₂'], () => { const a = mech.m.atoms.get(cm); a.rad = 0; a.chg = -1; a.hFix = 0; });
  return jumpTo(mech, o, 'Base removes the two N–H protons while both Cl⁻ leave (two α-eliminations): the isocyanide R–N⁺≡C⁻ (offensive smell). Only 1° amines do this.');
});
W('acid_heat', (m, o, c) => {
  const ac = M(m, '[CX3:1](=[OX1:2])[OX2H1:3]'); if (ac.length !== 1 || m.frags().length !== 1 || !M(o.product, '[CX3](=O)O[CX3]=O').length) return null; const [c1, o1, o3] = ac[0];
  const mech = new Mech(m);
  const g = C.fromSmilesMB(C.canon(m), BL); for (const a of g.atoms.values()) a.hFix = null;
  const q2 = M(g, '[CX3:1](=[OX1:2])[OX2H1:3]')[0]; const map = mech.add(g, c1, { att: q2[2], dist: 1.9 * BL }); const on = map.get(q2[2]);
  mech.step([[[on], [on, c1]], [[c1, o1], [o1]]], 'P₂O₅ is a powerful dehydrating agent: it activates the acid (it grabs water as it forms). The –OH oxygen of a second acid molecule attacks the C=O carbon.', [c.has('P2O5') ? 'P₂O₅' : 'Δ']);
  return jumpTo(mech, o, 'Proton transfers, then water leaves (taken up by P₂O₅ as H₃PO₄): the acid anhydride (RCO)₂O.');
});
W('amide_rx', (m, o, c) => {
  if (!c.has('P2O5') && !c.has('SOCl2')) return null;
  const q = M(m, '[CX3:1](=[OX1:2])[NX3H2:3]')[0]; if (!q || !M(o.product, 'C#N').length) return null; const [c1, ox, n] = q;
  const mech = new Mech(m);
  const s = put(mech, 'O=P(=O)[O-]', ox, isEl('P'), { dist: 1.7 * BL }); const pp = s.at; const op = mech.m.nbrs(pp).find(x => el(mech.m, x) === 'O' && mech.m.bondBetween(pp, x).o === 2);
  mech.step([[[ox], [ox, pp]], [[pp, op], [op]]], 'The amide oxygen (the more nucleophilic end, thanks to N lone-pair donation) attacks phosphorus: the O becomes part of a good leaving group.', ['P₂O₅', 'Δ']);
  const h1 = mech.H(n); const ph = P(mech.m, h1), pn = P(mech.m, n); const w = mech.add('O', h1, { dist: 1.45 * BL, dir: unit([ph[0] - pn[0], ph[1] - pn[1]]) }).get(1);
  mech.step([[[w], [w, h1]], [[h1, n], [n, c1]], [[c1, ox], [ox]]], 'Loss of an N–H proton: the N lone pair forms C=N as the C=O π pair moves to O: an imidate (O–P ester).'); mech.dropFrag(w);
  const h2 = mech.H(n); const ph2 = P(mech.m, h2), pn2 = P(mech.m, n); const w2 = mech.add('O', h2, { dist: 1.45 * BL, dir: unit([ph2[0] - pn2[0], ph2[1] - pn2[1]]) }).get(1);
  mech.step([[[w2], [w2, h2]], [[h2, n], [n, c1]], [[c1, ox], [ox]]], 'The second N–H goes the same way and the phosphate leaves from carbon: C≡N, the nitrile.'); mech.dropFrag(w2); mech.dropFrag(ox);
  return fin(mech, o);
});
W('benzilic', (m, o) => {
  const q = M(m, '[CX3H1:1](=[OX1:2])[CX3H1:3]=[OX1:4]')[0]; if (!q) return null; const [c1, o1, c2, o2] = q;
  const mech = new Mech(m);
  const w = put(mech, '[OH-]', c1, isEl('O'), { dist: 1.65 * BL }).at;
  mech.step([[[w], [w, c1]], [[c1, o1], [o1]]], 'OH⁻ adds to one aldehyde C=O.', ['NaOH']);
  const h = mech.H(c1); if (h == null) return null;
  mech.step([[[o1], [o1, c1]], [[c1, h], [h, c2]], [[c2, o2], [o2]]], 'The O⁻ pushes back down and the C–H on that carbon moves as a hydride to the neighbouring C=O (intramolecular Cannizzaro): one end is oxidised, the other reduced.');
  return jumpTo(mech, o, 'Proton transfer from –COOH to the alkoxide: the glycolate ion, HOCH₂COO⁻.');
});
W('phenolphthalein', (m, o, c) => {
  const q = M(m, '[CX3:1](=[OX1:2])[OX2][CX3]=O')[0]; if (!q) return null; const [c1, o1] = q;
  const mech = new Mech(m);
  prot(mech, o1, 'Conc. H₂SO₄ protonates one anhydride C=O: the carbon becomes a strong electrophile.', ['conc. H₂SO₄']);
  const g = C.fromSmilesMB('Oc1ccccc1', BL); for (const a of g.atoms.values()) a.hFix = null;
  const oa = [...g.atoms.values()].find(a => a.el === 'O'); const ip = g.nbrs(oa.id)[0]; const pa = ringPara(g, ip); if (pa == null) return null;
  const map = mech.add(g, c1, { att: pa, dist: 1.8 * BL }); const pc = map.get(pa);
  const pb = mech.m.bondsOf(pc).find(b => b.o === 2); if (!pb) return null; const pn = pb.a === pc ? pb.b : pb.a;
  mech.step([[[pc, pn], [pc, c1]], [[c1, o1], [o1]]], 'The para carbon of phenol (–OH is a strong o/p director) attacks the activated carbonyl carbon.', ['phenol']);
  deprot(mech, pc, 'Loss of H⁺: the ring is aromatic again.', 'O', pn);
  return jumpTo(mech, o, 'The new –OH is protonated and leaves as water; the cation formed is attacked by a second phenol (para). Loss of H⁺: phenolphthalein, the lactone form (colourless in acid, pink in base).');
});

/* ================= batch 9 (run 26) ================= */
W('acid_nh3', (m, o, c) => {
  const ac = M(m, '[CX3:1](=[OX1:2])[OX2H1:3]'); if (ac.length !== 2 || !M(o.product, '[CX3](=O)[NX3][CX3]=O').length) return null;
  const [c1, o1, o3] = ac[0]; const mech = new Mech(m);
  const h = mech.H(o3); const ph = P(mech.m, h), pa = P(mech.m, o3);
  const nN = mech.add('N', h, { dist: 1.45 * BL, dir: unit([ph[0] - pa[0], ph[1] - pa[1]]) }).get(1);
  mech.step([[[nN], [nN, h]], [[h, o3], [o3]]], 'Acid–base first: NH₃ takes H⁺ from a –COOH (with excess NH₃, from both): ammonium succinate.', ['NH₃']); mech.dropFrag(nN);
  return jumpTo(mech, o, 'Strong heating: NH₃ is released and attacks one acyl carbon (amide, loss of water); the amide N then attacks the second –COOH inside the molecule and a second water leaves. The five-membered cyclic imide (succinimide).', ['Δ']);
});
W('peracetyl', (m, o, c) => {
  const ohs = M(m, '[CX4][OX2H1:1]').map(x => x[1]); if (ohs.length < 3) return null; const ox = ohs[0];
  const mech = new Mech(m);
  const s = put(mech, 'CC(=O)OC(C)=O', ox, (a, g) => a.el === 'C' && g.bondsOf(a.id).some(b => b.o === 2), { dist: 1.8 * BL }); const ca = s.at;
  const oa = dO(mech.m, ca); const ol = mech.m.nbrs(ca).find(n => el(mech.m, n) === 'O' && mech.m.bondBetween(ca, n).o === 1);
  mech.step([[[ox], [ox, ca]], [[ca, oa], [oa]]], 'An –OH oxygen attacks a C=O of acetic anhydride (tetrahedral intermediate).', ['(CH₃CO)₂O']);
  mech.step([[[oa], [oa, ca]], [[ca, ol], [ol]]], 'The C=O re-forms and acetate leaves.'); mech.dropFrag(ol);
  deprot(mech, ox, 'Acetate takes the H⁺ from O: one ester (–OCOCH₃) made.', '[O-]C(C)=O');
  return jumpTo(mech, o, `The same happens at every –OH: all ${ohs.length} are acetylated (glucose gives a pentaacetate, proof of five –OH groups).`);
});
W('haworth', (m, o, c) => {
  const bz = C.canon(m) === 'c1ccccc1', anh = M(m, 'O=C1CCC(=O)O1').length > 0; if (!bz && !anh) return null;
  const mech = new Mech(m); let pc, cc;
  if (bz) { const g = C.fromSmilesMB('O=C1CCC(=O)O1', BL); for (const a of g.atoms.values()) a.hFix = null;
    const cA = [...g.atoms.values()].find(a => a.el === 'C' && g.bondsOf(a.id).some(b => b.o === 2));
    pc = [...m.atoms.keys()].find(id => el(m, id) === 'C' || el(m, id) === 'c'); const map = mech.add(g, pc, { att: cA.id, dist: 1.9 * BL }); cc = map.get(cA.id); }
  else { cc = M(m, '[CX3:1](=O)[OX2]')[0][0]; const g = C.fromSmilesMB('c1ccccc1', BL); for (const a of g.atoms.values()) a.hFix = null; const k0 = [...g.atoms.keys()][0];
    const map = mech.add(g, cc, { att: k0, dist: 1.9 * BL }); pc = map.get(k0); }
  const pb = mech.m.bondsOf(pc).find(b => b.o === 2); if (!pb) return null; const pn = pb.a === pc ? pb.b : pb.a;
  const oR = mech.m.nbrs(cc).find(n => el(mech.m, n) === 'O' && mech.m.bondBetween(cc, n).o === 1);
  const oD = dO(mech.m, cc);
  mech.step([[[pc, pn], [pc, cc]], [[cc, oR], [oR]]], 'AlCl₃ opens succinic anhydride to an acylium ion (O=C⁺–CH₂CH₂COO–AlCl₃⁻); benzene attacks it: Friedel–Crafts acylation (no rearrangement, no polysubstitution).', ['AlCl₃']);
  deprot(mech, pc, 'Loss of H⁺: β-benzoylpropionic acid.', 'O', pn);
  return jumpTo(mech, o, c.has('ZnHg') ? 'Clemmensen (Zn–Hg/HCl) reduces the ketone C=O to CH₂; conc. H₂SO₄ then closes the ring by an intramolecular Friedel–Crafts acylation: α-tetralone. (Clemmensen again, then Se/Δ, lead on to naphthalene.)' : 'Next steps of the Haworth sequence follow as shown.');
});

/* ================= batch 10 (run 27) ================= */
W('arndt', (m, o, c) => {
  const h = M(m, '[#6:1][CX3:2](=[OX1:3])[OX2H1,Cl:4]')[0]; if (!h) return null; const [r, cc, ox, lg] = h;
  const mech = new Mech(m);
  if (el(m, lg) === 'O') R8.settle(mech, 'SOCl₂ turns the acid into the acid chloride (–OH → –Cl; SO₂ and HCl leave).', t => { const a = t.atoms.get(lg); for (const x of t.nbrs(lg).filter(n => el(t, n) === 'H')) t.removeAtom(x); a.el = 'Cl'; a.hFix = null; a.chg = 0; }, ['SOCl₂']);
  const s = put(mech, '[CH2-][N+]#N', cc, a => a.el === 'C', { dist: 1.8 * BL }); const c2 = s.at; const n1 = mech.m.nbrs(c2).find(x => el(mech.m, x) === 'N');
  mech.step([[[c2], [c2, cc]], [[cc, ox], [ox]]], 'The carbon of diazomethane (CH₂⁻–N⁺≡N) attacks the acyl carbon.', ['CH₂N₂']);
  mech.step([[[ox], [ox, cc]], [[cc, lg], [lg]]], 'The C=O re-forms and Cl⁻ leaves (addition–elimination).'); mech.dropFrag(lg);
  const hh = mech.H(c2); if (hh == null) return null; const ph = P(mech.m, hh), pc = P(mech.m, c2);
  const b = mech.add('[CH2-][N+]#N', hh, { dist: 1.45 * BL, dir: unit([ph[0] - pc[0], ph[1] - pc[1]]) }).get(1);
  mech.step([[[b], [b, hh]], [[hh, c2], [c2]]], 'A second CH₂N₂ removes the acidic H next to C=O (it becomes CH₃Cl with the HCl): the diazoketone R–CO–CH⁻–N⁺≡N.'); mech.dropFrag(b);
  mech.step([[[c2, n1], [n1]], [[cc, r], [r, c2]], [[c2], [c2, cc]]], 'Wolff rearrangement (Ag₂O, or light/heat): N₂ leaves and, in the same moment, R moves from the carbonyl carbon to the neighbouring carbon. The result is a ketene, R–CH=C=O.', ['Ag₂O']);
  mech.dropFrag(n1);
  const w = put(mech, 'O', cc, isEl('O'), { dist: 1.65 * BL }).at;
  mech.step([[[w], [w, cc]], [[cc, c2], [c2]]], 'Water attacks the ketene carbon (sp, very electrophilic); the C=C pair moves to the CH carbon.', ['H₂O']);
  return jumpTo(mech, o, 'Proton transfers (O–H⁺ to the carbanion carbon): the carboxylic acid R–CH₂–COOH, one CH₂ longer than the start.');
});
W('sf4', (m, o) => {
  const q = M(m, '[CX3:1]=[OX1:2]')[0]; if (!q || !M(o.product, '[CX4](F)F').length) return null; const [cc, ox] = q;
  const mech = new Mech(m);
  const s = put(mech, 'FS(F)(F)F', ox, isEl('S'), { dist: 1.75 * BL }); const sS = s.at; const f1 = mech.m.nbrs(sS).find(x => el(mech.m, x) === 'F');
  mech.step([[[ox], [ox, sS]], [[sS, f1], [f1]]], 'The carbonyl O attacks sulfur of SF₄ and a fluoride leaves.', ['SF₄']);
  mech.step([[[f1], [f1, cc]], [[cc, ox], [ox]]], 'F⁻ adds to the carbon of C=O⁺–SF₃ (the π pair moves back to O).');
  const f2 = mech.m.nbrs(sS).find(x => el(mech.m, x) === 'F');
  mech.step([[[cc, ox], [ox]]], 'The C–O bond breaks: ⁻O–SF₃ leaves (it falls apart to O=SF₂ + F⁻) and a fluoro-carbocation forms.', [], () => {});
  mech.dropFrag(ox);
  const f3 = put(mech, '[F-]', cc, isEl('F'), { dist: 1.6 * BL }).at;
  mech.step([[[f3], [f3, cc]]], 'Fluoride captures the cation: the gem-difluoride (C=O → CF₂).');
  return fin(mech, o);
});

root.EngineR11 = {};
})(typeof window !== 'undefined' ? window : globalThis);
