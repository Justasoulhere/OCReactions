/* engine_r10.js: curly-arrow mechanisms, batch 3 (Oct 15). Same contract as engine_r7 / r8: each wrapper rebuilds the
   reaction as arrow steps on the same atoms; if its product differs from the plain recipe, or any step is not a valid
   structure, the plain version is kept. Recipes whose result already has arrows are left alone. */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore, R7 = root.EngineR7, R8 = root.EngineR8;
if (!R7 || !R8) return;
const { BL, el, P, unit, hTot, cN, Mech } = E; const { wrap, key, makeDouble } = R7; const { prot, deprot, settle, kek, groupPost, setO, dropH, clr } = R8;
const M = (m, sm) => { try { return C.matchAll(m, sm, {}); } catch (e) { return []; } };
const AR = m => E.aromSet(m);
const hasArrows = o => (o.steps || []).some(s => (s.arrows || []).length);
/* only build when the plain result has no arrows */
const W = (id, build) => wrap(id, (m, o, c, all, k) => hasArrows(o) ? null : build(m, o, c, all, k));
/* the fragment of the working structure that matches the plain product (falls back to the main fragment) */
function fin(mech, o) {
  const want = key(o.product);
  for (const ids of mech.m.frags()) { const p = E.cleanProduct(mech.m, ids); if (key(p) === want) return { steps: mech.steps, product: p }; }
  return { steps: mech.steps, product: mech.product() };
}
/* a reagent molecule with real 2D coordinates; hFix cleared so H counts follow valence */
function spec(smi) { const g = C.fromSmilesMB(smi, BL); if (!g) return null; for (const a of g.atoms.values()) a.hFix = null; return g; }
/* place a reagent next to `site`; `pick(atom, mol)` chooses the attacking atom; returns its id in the working structure */
function put(mech, smi, site, pick, opt) {
  const g = typeof smi === 'string' ? spec(smi) : smi; if (!g) return null;
  const at = [...g.atoms.values()].find(a => pick(a, g)); if (!at) return null;
  const map = mech.add(g, site, Object.assign({ att: at.id, dist: 1.75 * BL }, opt || {}));
  return { map, at: map.get(at.id), ids: [...map.values()] };
}
const isEl = e => a => a.el === e;
const firstKey = g => [...g.atoms.keys()][0];
const nb = (m, a, e) => m.nbrs(a).filter(n => el(m, n) === e);
const dblO = (m, c) => m.nbrs(c).find(n => el(m, n) === 'O' && m.bondBetween(c, n).o === 2);
function carbonyls(m) { return M(m, '[CX3:1]=[OX1:2]').filter(([c, o]) => m.nbrs(c).filter(n => n !== o && el(m, n) !== 'H').every(n => el(m, n) === 'C')); }
/* copy every substituent the plain product has gained on base atoms (except `skip`) onto t: used for "the same again" steps */
function graftFrom(t, prod, baseIds, skip) {
  const base = new Set(baseIds); const sk = new Set(skip || []);
  for (const a of baseIds) {
    if (sk.has(a) || !t.atoms.has(a) || !prod.atoms.has(a)) continue;
    for (const n of prod.nbrs(a)) {
      if (base.has(n)) continue;
      const map = new Map(); const q = [n]; const seen = new Set([n]);
      while (q.length) { const x = q.shift(); const pa = prod.atoms.get(x); map.set(x, t.addAtom(pa.el, pa.x, pa.y, { chg: pa.chg, rad: pa.rad }));
        for (const y of prod.nbrs(x)) if (!base.has(y) && !seen.has(y) && el(prod, y) !== 'H') { seen.add(y); q.push(y); } }
      for (const [x, nx] of map) for (const bd of prod.bondsOf(x)) { const y = prod.other(bd, x); if (map.has(y) && x < y) t.addBond(nx, map.get(y), bd.o); }
      t.addBond(a, map.get(n), prod.bondBetween(a, n).o); dropH(t, a); clr(t, [a]);
    }
  }
}
/* nucleophile lone pair onto the acyl carbon, then the C=O re-forms and `lg` leaves (addition–elimination) */
function addElim(mech, nu, c1, o2, lg, t1, t2, cond) {
  mech.step([[[nu], [nu, c1]], [[c1, o2], [o2]]], t1, cond);
  mech.step([[[o2], [o2, c1]], [[c1, lg], [lg]]], t2);
}
const chg = (m, a) => (m.atoms.get(a) || {}).chg || 0;

/* ================= Grignard reagents ================= */
const CARB = { HCHO: 'C=O', CH3CHO: 'CC=O', acetone: 'CC(C)=O', PhCHO: 'O=Cc1ccccc1' };
const carbC = (a, g) => a.el === 'C' && g.bondsOf(a.id).some(b => b.o === 2 && g.atoms.get(g.other(b, a.id)).el === 'O');
const POL = 'The C–Mg bond is strongly polarised (Cδ⁻–Mgδ⁺), so its carbon behaves as a carbanion';
W('grignard_sub', (m, o, c) => {
  const g = M(m, '[#6:1][Mg:2]')[0]; if (!g) return null; const [rc, mg] = g; const mech = new Mech(m);
  if (c.has('CO2')) {
    const s = put(mech, 'O=C=O', rc, isEl('C')); const cc = s.at, ox = nb(mech.m, cc, 'O')[0];
    mech.step([[[rc, mg], [rc, cc]], [[cc, ox], [ox]]], POL + ': it attacks the carbon of CO₂ and one C=O π pair moves onto oxygen. A magnesium carboxylate (RCOO⁻ ⁺MgX) forms.', ['dry ether']);
    mech.dropFrag(mg); prot(mech, ox, 'Work-up: H₃O⁺ protonates the carboxylate: the carboxylic acid, one carbon longer than R–X.', ['H₃O⁺']);
    return fin(mech, o);
  }
  const k = ['HCHO', 'CH3CHO', 'acetone', 'PhCHO'].find(x => c.has(x));
  if (k) {
    const s = put(mech, CARB[k], rc, carbC); const cc = s.at, ox = dblO(mech.m, cc);
    mech.step([[[rc, mg], [rc, cc]], [[cc, ox], [ox]]], POL + ': it adds to the carbonyl carbon and the C=O π electrons move onto oxygen (a magnesium alkoxide).', ['dry ether']);
    mech.dropFrag(mg); prot(mech, ox, `Work-up: H₃O⁺ protonates the alkoxide: a ${k === 'HCHO' ? '1°' : k === 'acetone' ? '3°' : '2°'} alcohol.`, ['H₃O⁺']);
    return fin(mech, o);
  }
  if (c.has('ethyleneoxide')) {
    const s = put(mech, 'C1CO1', rc, isEl('C')); const cc = s.at, ox = nb(mech.m, cc, 'O')[0];
    mech.step([[[rc, mg], [rc, cc]], [[cc, ox], [ox]]], POL + ': it attacks a ring carbon from the back (SN2) and the strained three-membered ring opens.', ['dry ether']);
    mech.dropFrag(mg); prot(mech, ox, 'Work-up: H₃O⁺ protonates the alkoxide: R–CH₂CH₂–OH, two carbons longer.', ['H₃O⁺']);
    return fin(mech, o);
  }
  const ids = mech.add('[H]O', rc, { dist: 1.6 * BL }); const h = ids.get(1), w = ids.get(2);
  mech.step([[[rc, mg], [rc, h]], [[h, w], [w]]], POL + ' and a very strong base: it takes H⁺ from water (or any O–H / N–H). R–H forms and Mg(OH)X is left.', []);
  mech.dropFrag(mg); mech.dropFrag(w);
  return fin(mech, o);
});
function grig(c) { const g = c.get('RMgX') || c.get('RLi') || c.get('R2CuLi'); if (!g) return null; const rs = g.smi || 'C'; const X = (c.get('RMgX') || {}).X || 'Br';
  const li = !c.get('RMgX') && c.get('RLi'); const sp = spec(rs + (li ? '[Li]' : '[Mg]' + X)); if (!sp) return null;
  const r = firstKey(sp); const mt = [...sp.atoms.values()].find(a => a.el === (li ? 'Li' : 'Mg')).id; return { sp, r, mt, li }; }
function grigAt(mech, c, site) { const G = grig(c); if (!G) return null; const map = mech.add(G.sp.clone ? G.sp.clone() : G.sp, site, { att: G.r, dist: 1.8 * BL }); return { r: map.get(G.r), mt: map.get(G.mt) }; }
W('grignard_add', (m, o, c) => {
  const mech = new Mech(m);
  /* acidic H first */
  const acidic = M(m, '[OX2H1:1]').map(x => x[0]).concat(M(m, '[NX3;!H0:1]').map(x => x[0]), M(m, '[CH1:1]#C').map(x => x[0]));
  if (acidic.length) {
    const x = acidic[0]; const h = mech.H(x); if (h == null) return null; const g = grigAt(mech, c, h); if (!g) return null;
    mech.step([[[g.r, g.mt], [g.r, h]], [[h, x], [x]]], 'R–MgX is a far stronger base than it is a nucleophile here: its carbanion carbon takes the acidic H. R–H forms and no addition happens.', ['dry ether']);
    mech.dropFrag(g.mt); return fin(mech, o);
  }
  const es = M(m, '[CX3:1](=[OX1:2])[OX2:3][#6:4]'), ax = M(m, '[CX3:1](=[OX1:2])[Cl,Br:3]');
  if (es.length || ax.length) {
    const [c1, o2, lg] = es.length ? es[0] : ax[0]; const one = c.cold || c.flag('one');
    let g = grigAt(mech, c, c1); if (!g) return null;
    mech.step([[[g.r, g.mt], [g.r, c1]], [[c1, o2], [o2]]], 'The Grignard carbon adds to the acyl carbon; the C=O π pair moves onto O: a tetrahedral alkoxide.', ['dry ether']); mech.dropFrag(g.mt);
    mech.step([[[o2], [o2, c1]], [[c1, lg], [lg]]], `The C=O re-forms and ${es.length ? 'the alkoxide OR′' : 'Cl⁻'} leaves: a ketone.${one ? '' : ' The ketone is more reactive than the starting ester, so it never survives.'}`); mech.dropFrag(lg);
    if (!one) {
      g = grigAt(mech, c, c1); if (!g) return null;
      mech.step([[[g.r, g.mt], [g.r, c1]], [[c1, o2], [o2]]], 'A second R–MgX adds to the ketone at once: the alkoxide of a 3° alcohol with two identical R groups.'); mech.dropFrag(g.mt);
      prot(mech, o2, 'Work-up: H₃O⁺ protonates the alkoxide: the 3° alcohol.', ['H₃O⁺']);
    }
    return fin(mech, o);
  }
  const ns = M(m, '[CX2:1]#[NX1:2]');
  if (ns.length) {
    const [c1, n2] = ns[0]; const g = grigAt(mech, c, c1); if (!g) return null;
    mech.step([[[g.r, g.mt], [g.r, c1]], [[c1, n2], [n2]]], 'The Grignard carbon adds to the nitrile carbon; one C≡N π pair moves onto N: the magnesium salt of an imine (it does not react again).', ['dry ether']); mech.dropFrag(g.mt);
    prot(mech, n2, 'Work-up: H₃O⁺ protonates the N: an imine, R₂C=NH.', ['H₃O⁺']);
    settle(mech, 'In aqueous acid the imine is hydrolysed (water adds to C=N, NH₃ is lost): the ketone.', t => { const a = t.atoms.get(n2); dropH(t, n2); a.el = 'O'; a.chg = 0; a.hFix = null; setO(t, c1, n2, 2); });
    return fin(mech, o);
  }
  const ep = M(m, '[C:1]1[O:2][C:3]1');
  if (ep.length) {
    const [c1, o2, c3] = ep[0]; const less = cN(m, c1).length <= cN(m, c3).length ? c1 : c3; const g = grigAt(mech, c, less); if (!g) return null;
    mech.step([[[g.r, g.mt], [g.r, less]], [[less, o2], [o2]]], 'The Grignard carbon attacks the less hindered ring carbon from the back (SN2); the strained epoxide ring opens.', ['dry ether']); mech.dropFrag(g.mt);
    prot(mech, o2, 'Work-up: H₃O⁺ protonates the alkoxide: an alcohol two carbons further along the chain.', ['H₃O⁺']);
    return fin(mech, o);
  }
  return null;
});

/* ================= C=O + H2N–Z ================= */
const NZ = [['NH2OH', 'NO'], ['PhNHNH2', 'NNc1ccccc1'], ['DNP', 'NNc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-]'], ['semicarbazide', 'NNC(N)=O'], ['NH2NH2', 'NN'], ['RNH2', null], ['aniline', 'Nc1ccccc1'], ['NH3', 'N']];
W('imine', (m, o, c) => {
  const cs = carbonyls(m); if (!cs.length) return null; const [cc, ox] = cs[0]; const z = NZ.find(q => c.has(q[0])); if (!z) return null;
  const smi = z[1] || ('N' + ((c.get('RNH2') || {}).smi || 'C'));
  const g = spec(smi); if (!g) return null;
  const pickN = (a, G) => a.el === 'N' && !a.chg && G.hImp(a.id) >= 2 && !G.nbrs(a.id).some(n => E.isCarbonylC(G, n));
  const mech = new Mech(m); const s = put(mech, g, cc, pickN, { dist: 1.7 * BL }); if (!s) return null; const n = s.at;
  mech.step([[[n], [n, cc]], [[cc, ox], [ox]]], 'The lone pair of the –NH₂ nitrogen attacks the carbonyl carbon (weak acid, pH ≈ 4.5, has made the C=O more electrophilic); the π pair moves onto O.', ['H⁺ (cat.)', 'pH ≈ 4.5']);
  prot(mech, ox, 'Proton transfer: the alkoxide O picks up H⁺…');
  deprot(mech, n, '…and the N loses its extra H⁺: the neutral tetrahedral carbinolamine.');
  prot(mech, ox, 'The –OH is protonated: –OH₂⁺ is a good leaving group.');
  mech.step([[[n], [n, cc]], [[cc, ox], [ox]]], 'The N lone pair pushes in to form C=N as water leaves: an iminium ion.'); mech.dropFrag(ox);
  deprot(mech, n, 'Loss of H⁺ from N: the C=N–Z product. At low pH the amine is protonated and cannot attack; at high pH the –OH cannot be protonated. Hence the optimum near pH 4.5.');
  return fin(mech, o);
});

/* ================= acid derivatives ================= */
function nucOf(c) { return c.has('NH3') ? 'N' : c.has('aniline') ? 'Nc1ccccc1' : c.has('RNH2') ? 'N' + (c.get('RNH2').smi || 'C') : c.has('R2NH') ? 'N(' + (c.get('R2NH').smi || 'C') + ')' + (c.get('R2NH').smi || 'C') : c.has('phenol') ? 'Oc1ccccc1' : c.has('ROH') ? 'O' + (c.get('ROH').smi || 'C') : c.has('NaOAc') ? '[O-]C(C)=O' : c.has('H2O') ? 'O' : null; }
const pickNu = (a, g) => (a.el === 'N' || a.el === 'O') && (a.chg < 0 || (g.hImp(a.id) >= 1 && !g.nbrs(a.id).some(n => E.isCarbonylC(g, n))));
W('acyl_sub', (m, o, c) => {
  const nuS = nucOf(c); if (!nuS) return null;
  const ax = M(m, '[CX3:1](=[OX1:2])[Cl,Br:3]'), an = M(m, '[CX3:1](=[OX1:2])[OX2:3][CX3:4]=[OX1:5]'), es = M(m, '[CX3:1](=[OX1:2])[OX2:3][#6:4]');
  const hit = ax[0] || an[0] || (/^N/.test(nuS) ? es[0] : null); if (!hit) return null; const [c1, o2, lg] = hit;
  const mech = new Mech(m); const s = put(mech, nuS, c1, pickNu, { dist: 1.7 * BL }); if (!s) return null; const n = s.at;
  const lgName = ax.length ? 'Cl⁻' : an.length ? 'the carboxylate ion' : 'the alkoxide';
  addElim(mech, n, c1, o2, lg, 'The nucleophile adds to the acyl carbon and the C=O π pair moves onto O: a tetrahedral intermediate.', `The C=O re-forms and ${lgName} leaves (addition–elimination).`, []);
  mech.dropFrag(lg);
  if (chg(mech.m, n) > 0) deprot(mech, n, `Loss of H⁺ (taken by ${/^N/.test(nuS) ? 'a second molecule of the amine' : 'Cl⁻ / pyridine'}): the neutral product.`);
  return fin(mech, o);
});
const ACYL = c => c.has('BzCl') ? 'ClC(=O)c1ccccc1' : c.has('RCOCl') ? 'ClC(=O)' + (c.get('RCOCl').smi || 'C') : c.has('Ac2O') ? 'CC(=O)OC(C)=O' : 'CC(=O)Cl';
W('acylation', (m, o, c) => {
  if (c.has('Hinsberg')) return null;
  const nh = M(m, '[NX3;!H0;!$(NC=O):1]'), oh = M(m, '[OX2H1;!$(OC=O):1]'); const site = nh.length ? nh[0][0] : oh.length ? oh[0][0] : null; if (site == null) return null;
  const g = spec(ACYL(c)); if (!g) return null;
  const acyC = [...g.atoms.values()].find(a => carbC(a, g)); if (!acyC) return null;
  const mech = new Mech(m); const map = mech.add(g, site, { att: acyC.id, dist: 1.75 * BL }); const c1 = map.get(acyC.id);
  const o2 = dblO(mech.m, c1); const lg = mech.m.nbrs(c1).find(x => ['Cl', 'Br'].includes(el(mech.m, x))) ?? mech.m.nbrs(c1).find(x => el(mech.m, x) === 'O' && x !== o2);
  if (lg == null) return null; const isAn = el(mech.m, lg) === 'O';
  addElim(mech, site, c1, o2, lg, `The ${nh.length ? 'N' : 'O'} lone pair attacks the acyl carbon of ${isAn ? 'the anhydride' : 'the acid chloride'}; the π pair moves onto O.`, `The C=O re-forms and ${isAn ? 'acetate' : 'Cl⁻'} leaves.`, [isAn ? '(CH₃CO)₂O' : 'RCOCl', 'pyridine']);
  mech.dropFrag(lg);
  deprot(mech, site, `The base (pyridine${isAn ? ' or acetate' : ''}) removes H⁺ from the ${nh.length ? 'N' : 'O'}: the ${nh.length ? 'amide' : 'ester'}. Pyridine also mops up the HCl, which pushes the equilibrium forward.`, 'c1ccncc1'.length && 'N');
  return fin(mech, o);
});
/* Fischer esterification: activate C=O, attack, proton shuffle, water leaves, deprotonate */
function fischerCore(mech, c1, o2, o3, nO, intro) {
  prot(mech, o2, intro || 'H⁺ (conc. H₂SO₄) protonates the C=O oxygen of the acid: the carbonyl carbon becomes much more electrophilic.', ['conc. H₂SO₄']);
  mech.step([[[nO], [nO, c1]], [[c1, o2], [o2]]], 'The alcohol oxygen attacks the carbonyl carbon: a tetrahedral intermediate.');
  deprot(mech, nO, 'Proton transfer: the new O loses H⁺…');
  prot(mech, o3, '…and an –OH is protonated, so it can leave as water.');
  mech.step([[[o2], [o2, c1]], [[c1, o3], [o3]]], 'An O lone pair pushes back in and water leaves: the protonated ester.'); mech.dropFrag(o3);
  deprot(mech, o2, 'Loss of H⁺: the ester. Every step is reversible (remove water or use excess alcohol); the ester O comes from the alcohol.');
}
W('fischer', (m, o, c) => {
  const ac = M(m, '[CX3:1](=[OX1:2])[OX2H1:3]'); const mech = new Mech(m);
  if (ac.length && c.has('ROH')) { const [c1, o2, o3] = ac[0]; const s = put(mech, 'O' + (c.get('ROH').smi || 'C'), c1, a => a.el === 'O', { dist: 1.75 * BL }); if (!s) return null; fischerCore(mech, c1, o2, o3, s.at); return fin(mech, o); }
  const al = M(m, '[CX4:1][OX2H1:2]'); if (!al.length || !c.has('AcOH')) return null; const nO = al[0][1];
  const g = spec('CC(=O)O'); const cA = [...g.atoms.values()].find(a => carbC(a, g)); const map = mech.add(g, nO, { att: cA.id, dist: 1.8 * BL }); const c1 = map.get(cA.id);
  const o2 = dblO(mech.m, c1), o3 = mech.m.nbrs(c1).find(x => el(mech.m, x) === 'O' && x !== o2);
  fischerCore(mech, c1, o2, o3, nO, 'H⁺ protonates the C=O oxygen of acetic acid.');
  return fin(mech, o);
});

/* ================= SOCl2 / PX3 / PCl5 ================= */
W('alc_pcl', (m, o, c) => {
  const X = c.has('PBr3') ? 'Br' : c.has('PI3') ? 'I' : 'Cl'; const so = c.has('SOCl2'), p5 = c.has('PCl5');
  const acid = M(m, '[CX3:1](=[OX1:2])[OX2H1:3]'); const alc = M(m, '[CX4:1][OX2H1:2]');
  const mech = new Mech(m);
  const rgS = so ? 'ClS(Cl)=O' : p5 ? 'ClP(Cl)(Cl)(Cl)Cl' : `${X}P(${X})${X}`; const cen = so ? 'S' : 'P';
  const [oS] = acid.length ? [acid[0][2]] : alc.length ? [alc[0][1]] : [null]; if (oS == null) return null;
  const s = put(mech, rgS, oS, isEl(cen), { dist: 1.75 * BL }); if (!s) return null; const ce = s.at;
  const xl = mech.m.nbrs(ce).find(n => el(mech.m, n) === X || el(mech.m, n) === 'Cl');
  mech.step([[[oS], [oS, ce]], [[ce, xl], [xl]]], `The –OH oxygen attacks ${so ? 'the S of SOCl₂' : p5 ? 'the P of PCl₅' : `the P of P${X}₃`} and ${X === 'Cl' || so || p5 ? 'Cl⁻' : X + '⁻'} is pushed out. The O is now part of a very good leaving group.`, [so ? 'SOCl₂' : p5 ? 'PCl₅' : `P${X}₃`]);
  const xe = el(mech.m, xl); mech.dropFrag(xl);
  deprot(mech, oS, `${xe}⁻ (or pyridine) removes the H⁺ from the O${so ? ': an alkyl chlorosulphite, R–O–S(=O)Cl' : ''}.`, so && c.has('pyridine') ? 'N' : `[${xe}-]`);
  const nuX = put(mech, `[${so || p5 ? 'Cl' : X}-]`, acid.length ? acid[0][0] : alc[0][0], () => true, { dist: 1.7 * BL }).at;
  if (acid.length) {
    const [c1, o2] = acid[0];
    addElim(mech, nuX, c1, o2, oS, 'Cl⁻ adds to the acyl carbon: a tetrahedral intermediate.', so ? 'The C=O re-forms and ⁻O–S(=O)Cl leaves; it falls apart into SO₂↑ and Cl⁻. Both by-products are gases, so the acyl chloride is easy to purify.' : 'The C=O re-forms and the O–P group leaves (as POCl₃ / H₃PO₃): the acyl chloride.');
    mech.dropFrag(oS); return fin(mech, o);
  }
  const cc = alc[0][0];
  mech.step([[[nuX], [nuX, cc]], [[cc, oS], [oS]]], so ? (c.has('pyridine') ? 'With pyridine, free Cl⁻ attacks the carbon from the back (SN2, inversion) as ⁻O–S(=O)Cl leaves; that piece falls apart into SO₂↑ + Cl⁻.' : 'Without pyridine the chlorosulphite ionises as a tight ion pair and its own Cl is delivered to the front face (SNi, retention); SO₂↑ and HCl↑ escape. The arrows show the C–Cl bond forming as the C–O bond breaks.') : `${X}⁻ attacks the carbon from the back (SN2) as the O–P group leaves. No free carbocation, so no rearrangement.`);
  mech.dropFrag(oS); return fin(mech, o);
});

/* ================= epoxides, acid–base ================= */
W('epox_open', (m, o, c) => {
  const ep = M(m, '[C:1]1[O:2][C:3]1')[0]; if (!ep) return null; const [c1, ox, c3] = ep; const acid = c.any(['H3O+', 'H2SO4', 'H2O']);
  const more = cN(m, c1).length >= cN(m, c3).length ? c1 : c3, less = more === c1 ? c3 : c1; const at = acid ? more : less;
  const mech = new Mech(m);
  if (acid) {
    prot(mech, ox, 'H⁺ protonates the epoxide oxygen.', ['H₃O⁺']);
    const w = put(mech, 'O', at, isEl('O'), { dist: 1.65 * BL }).at;
    mech.step([[[w], [w, at]], [[at, ox], [ox]]], 'Water attacks the more substituted carbon (that C–O bond is weaker and the carbon carries more δ+) from the back side as the ring opens: anti opening.');
    deprot(mech, w, 'Loss of H⁺: the trans-1,2-diol.');
    return fin(mech, o);
  }
  const nuS = c.has('NH3') ? 'N' : c.has('NaOMe') ? '[O-]C' : '[O-]CC';
  const nu = put(mech, nuS, at, a => a.el === (c.has('NH3') ? 'N' : 'O'), { dist: 1.65 * BL }).at;
  mech.step([[[nu], [nu, at]], [[at, ox], [ox]]], 'The strong nucleophile attacks the less hindered ring carbon from the back (SN2) and the ring opens.', [c.has('NH3') ? 'NH₃' : c.has('NaOMe') ? 'NaOMe' : 'NaOEt']);
  if (chg(mech.m, nu) > 0) deprot(mech, nu, 'Proton transfer from N…', 'N');
  const ids = mech.add(c.has('NH3') ? '[H][NH3+]' : c.has('NaOMe') ? '[H]OC' : '[H]OCC', ox, { dist: 1.5 * BL }); const h = ids.get(1), b = ids.get(2);
  mech.step([[[ox], [ox, h]], [[h, b], [b]]], c.has('NH3') ? '…to the alkoxide O: the amino alcohol.' : 'The alkoxide takes H⁺ from the alcohol solvent (the alkoxide is regenerated).'); mech.dropFrag(b);
  return fin(mech, o);
});
W('acidbase', (m, o, c) => {
  const ac = M(m, '[CX3](=[OX1])[OX2H1]').map(x => x[2]), so = M(m, '[SX4](=O)(=O)[OX2H1]').map(x => x[3]), ph = M(m, '[c][OX2H1]').map(x => x[1]);
  const bic = c.has('NaHCO3'); const sites = bic ? ac.concat(so) : ac.concat(ph, so); if (!sites.length || sites.length > 3) return null;
  const mech = new Mech(m);
  for (const x of sites) { const h = mech.H(x); if (h == null) return null; const ph0 = P(mech.m, h), px = P(mech.m, x);
    const s = put(mech, bic ? '[O-]C(=O)O' : '[OH-]', h, a => a.el === 'O' && a.chg < 0, { dist: 1.45 * BL, dir: unit([ph0[0] - px[0], ph0[1] - px[1]]) }); const b = s.at;
    mech.step([[[b], [b, h]], [[h, x], [x]]], bic ? 'HCO₃⁻ takes the acidic H⁺; the H₂CO₃ formed breaks down to CO₂↑ + H₂O (the brisk effervescence).' : 'OH⁻ takes the acidic H⁺; the O–H pair stays on oxygen. The anion is resonance-stabilised.', [bic ? 'NaHCO₃' : 'NaOH']);
    mech.dropFrag(b); }
  return fin(mech, o);
});

/* ================= haloform ================= */
W('haloform', (m, o, c) => {
  const X = ((c.get('NaOCl') || {}).X) || (c.has('I2') ? 'I' : /br/i.test(c.raw) ? 'Br' : 'Cl');
  const mech = new Mech(m); let mk = M(m, '[CH3:1][CX3:2]=[OX1:3]')[0];
  if (!mk) { const sa = M(m, '[CH3:1][CX4H1:2][OX2H1:3]')[0]; if (!sa) return null; mk = sa;
    settle(mech, `${X}₂ in alkali first oxidises the CH₃CH(OH)– group to CH₃CO– (a methyl ketone).`, t => { dropH(t, sa[2]); dropH(t, sa[1]); setO(t, sa[1], sa[2], 2); }, [`${X}₂, OH⁻`]); }
  const [me, cc, ox] = mk;
  const enol = () => { const h = mech.H(me); const ph = P(mech.m, h), pa = P(mech.m, me); const b = mech.add('[OH-]', h, { dist: 1.45 * BL, dir: unit([ph[0] - pa[0], ph[1] - pa[1]]) }).get(1);
    mech.step([[[b], [b, h]], [[h, me], [me, cc]], [[cc, ox], [ox]]], 'OH⁻ removes an α-H from the CH₃ next to C=O: the enolate ion.', ['OH⁻']); mech.dropFrag(b); };
  enol();
  const xs = mech.add(`${X}${X}`, me, { dist: 1.7 * BL }); const x1 = xs.get(1), x2 = xs.get(2);
  mech.step([[[ox], [ox, cc]], [[cc, me], [me, x1]], [[x1, x2], [x2]]], `The enolate carbon attacks ${X}₂ and ${X}⁻ leaves: one α-${X}. Each ${X} pulls electrons, so the remaining α-H get more acidic and the next two rounds are faster.`, [`${X}₂`]); mech.dropFrag(x2);
  settle(mech, `The same enolate + ${X}₂ steps happen twice more: –C${X}₃.`, t => { const p = P(t, me); for (let i = 0; i < 2; i++) { const d = E.freeDir(t, me); const xa = t.addAtom(X, p[0] + d[0] * BL, p[1] + d[1] * BL); t.addBond(me, xa, 1); } dropH(t, me); clr(t, [me]); });
  const w = put(mech, '[OH-]', cc, isEl('O'), { dist: 1.65 * BL }).at;
  mech.step([[[w], [w, cc]], [[cc, ox], [ox]]], 'OH⁻ adds to the C=O carbon: a tetrahedral intermediate.');
  mech.step([[[ox], [ox, cc]], [[cc, me], [me]]], `The C=O re-forms and ⁻C${X}₃ leaves: three ${X} atoms stabilise that carbanion enough to make it a leaving group.`);
  const h = mech.H(w); if (h == null) return null;
  mech.step([[[me], [me, h]], [[h, w], [w]]], `⁻C${X}₃ takes the acid's H⁺: the carboxylate (one C fewer) and CH${X}₃${X === 'I' ? ' (iodoform, yellow precipitate)' : ''}.`);
  return fin(mech, o);
});

/* ================= diazonium salts (aryl cation routes) ================= */
W('diazo_rx', (m, o, c) => {
  const dz = M(m, '[c:1][N+:2]#[N:3]')[0]; if (!dz) return null; const [ca, n1] = dz;
  const mode = c.has('NaI') ? 'I' : c.has('HBF4') ? 'F' : (!['CuCN', 'CuBr', 'CuCl', 'H3PO2', 'ROH', 'phenol', 'aniline', 'naphthol'].some(k => c.has(k)) && !c.has('Cu') && c.any(['H2O', 'H3O+'])) ? 'OH' : null;
  if (!mode) return null;
  const mech = new Mech(m);
  if (mode === 'F') {
    const s = put(mech, 'F[B-](F)(F)F', ca, isEl('F'), { dist: 2.1 * BL }); const f = s.at, b = mech.m.nbrs(f)[0];
    mech.step([[[ca, n1], [n1]]], 'The dry diazonium fluoroborate is heated: N₂ leaves and an aryl cation forms.', ['Δ']); mech.dropFrag(n1);
    mech.step([[[f, b], [f, ca]]], 'BF₄⁻ hands over a fluoride to the aryl cation; BF₃ is released (Balz–Schiemann).'); mech.dropFrag(b);
    return fin(mech, o);
  }
  mech.step([[[ca, n1], [n1]]], 'N₂ is an excellent leaving group: it departs and leaves an aryl cation (the empty orbital lies in the ring plane).', mode === 'OH' ? ['warm', 'H₂O'] : ['KI']); mech.dropFrag(n1);
  if (mode === 'I') { const i = put(mech, '[I-]', ca, isEl('I'), { dist: 1.7 * BL }).at; mech.step([[[i], [i, ca]]], 'I⁻ is a strong enough nucleophile to capture the cation directly, so no Cu(I) catalyst is needed: Ar–I.'); return fin(mech, o); }
  const w = put(mech, 'O', ca, isEl('O'), { dist: 1.65 * BL }).at;
  mech.step([[[w], [w, ca]]], 'Water attacks the aryl cation (SN1Ar).');
  deprot(mech, w, 'Loss of H⁺: phenol.');
  return fin(mech, o);
});

/* ================= Hofmann bromamide ================= */
W('amide_rx', (m, o, c) => {
  if (!c.has('Br2/NaOH')) return null;
  const am = M(m, '[#6:4][CX3:1](=[OX1:2])[NX3H2:3]')[0]; if (!am) return null; const [r, cc, ox, n] = am;
  const mech = new Mech(m);
  deprot(mech, n, 'OH⁻ removes one N–H proton of the amide (the C=O makes it acidic).', '[OH-]', null, ['Br₂', 'NaOH']);
  const xs = mech.add('BrBr', n, { dist: 1.6 * BL }); const b1 = xs.get(1), b2 = xs.get(2);
  mech.step([[[n], [n, b1]], [[b1, b2], [b2]]], 'The amide anion attacks Br₂: the N-bromoamide.'); mech.dropFrag(b2);
  deprot(mech, n, 'OH⁻ removes the second N–H (now even more acidic).', '[OH-]');
  mech.step([[[cc, r], [r, n]], [[n, b1], [b1]], [[n], [n, cc]]], 'In one step R moves from C to N as Br⁻ leaves, and the N lone pair forms C=N: the isocyanate R–N=C=O. The migrating group keeps its configuration.', [], groupPost(mech, r, cc, n)); mech.dropFrag(b1);
  const w = put(mech, '[OH-]', cc, isEl('O'), { dist: 1.65 * BL }).at;
  mech.step([[[w], [w, cc]], [[cc, n], [n]]], 'OH⁻ adds to the isocyanate carbon; the C=N π pair moves onto N.');
  const ids = mech.add('[H]O', n, { dist: 1.5 * BL }); const h = ids.get(1), hw = ids.get(2);
  mech.step([[[n], [n, h]], [[h, hw], [hw]]], 'N takes H⁺ from water: a carbamic acid (as its anion).'); mech.dropFrag(hw);
  settle(mech, 'The carbamate loses CO₂ (as CO₃²⁻ in the alkali) and N picks up H: R–NH₂, one carbon shorter than the amide.', t => { const keep = new Set(E.fragIds(t, n)); t.removeAtom(w); t.removeAtom(ox); t.removeAtom(cc); for (const id of [...t.atoms.keys()]) if (!keep.has(id)) continue; const a = t.atoms.get(n); a.chg = 0; a.hFix = null; dropH(t, n); });
  return fin(mech, o);
});

/* ================= Reimer–Tiemann ================= */
W('phenol_rx', (m, o, c) => {
  if (!(c.has('CHCl3') && !c.has('CCl4'))) return null;
  const ph = M(m, '[c:1][OX2H1,OX1-:2]')[0]; if (!ph) return null; const [ca, ox] = ph;
  const ring = C.aromaticRings(m).find(rr => rr.includes(ca)); if (!ring) return null;
  const at = m.nbrs(ca).find(n => ring.includes(n) && hTot(m, n) > 0 && o.product.nbrs(n).some(q => !m.atoms.has(q))); if (at == null) return null;
  const mech = new Mech(m);
  if (chg(m, ox) !== -1) deprot(mech, ox, 'NaOH turns phenol into the phenoxide ion.', '[OH-]', null, ['NaOH']);
  const ids = mech.add(spec('ClC(Cl)Cl'), at, { att: 2, dist: 2.6 * BL });
  const cm = [...ids.values()].find(id => el(mech.m, id) === 'C');
  const h = mech.H(cm); const b = mech.add('[OH-]', h, { dist: 1.45 * BL, dir: unit([P(mech.m, h)[0] - P(mech.m, cm)[0], P(mech.m, h)[1] - P(mech.m, cm)[1]]) }).get(1);
  mech.step([[[b], [b, h]], [[h, cm], [cm]]], 'OH⁻ removes the H of chloroform (made acidic by three Cl): ⁻CCl₃.', ['CHCl₃', 'NaOH']); mech.dropFrag(b);
  const cl = nb(mech.m, cm, 'Cl')[0];
  mech.step([[[cm, cl], [cl]]], 'α-Elimination: Cl⁻ leaves from the carbanion carbon, giving dichlorocarbene :CCl₂, an electron-deficient electrophile.', [], () => { const a = mech.m.atoms.get(cm); a.chg = 0; a.rad = 2; a.hFix = 0; });
  mech.dropFrag(cl);
  if (!kek(mech.m, [[ca, at]]) && !makeDouble(mech.m, ca, at)) return null;
  mech.step([[[ox], [ox, ca]], [[ca, at], [at, cm]]], 'The electron-rich phenoxide ring attacks the carbene carbon at the ortho position (the O⁻ pushes electrons into the ring).', [], () => { const a = mech.m.atoms.get(cm); a.rad = 0; a.chg = -1; a.hFix = 0; });
  settle(mech, 'Proton shift: the ortho H moves onto the CCl₂⁻ carbon and the ring becomes aromatic again (o-dichloromethylphenoxide).', t => { setO(t, ox, ca, 1); t.atoms.get(ox).chg = -1; setO(t, ca, at, 1); const ringNow = C.rings(t, 6).find(rr => rr.length === 6 && rr.includes(ca)); dropH(t, at); clr(t, [at]); const a = t.atoms.get(cm); a.chg = 0; a.hFix = null; a.rad = 0; for (const rr of [ringNow]) if (rr) { const i = rr.indexOf(ca); /* restore alternating double bonds */ const seq = []; for (let k = 0; k < 6; k++) seq.push(rr[(i + k) % 6]); for (let k = 0; k < 6; k++) { const bd = t.bondBetween(seq[k], seq[(k + 1) % 6]); if (bd) bd.o = k % 2 ? 1 : 2; } clr(t, seq); } });
  settle(mech, 'OH⁻ replaces both Cl atoms of the –CHCl₂ group; the gem-diol loses water: –CHO. Acid work-up turns the phenoxide back into the phenol: salicylaldehyde (ortho major).', t => { for (const x of nb(t, cm, 'Cl')) t.removeAtom(x); const p = P(t, cm), d = E.freeDir(t, cm); const oo = t.addAtom('O', p[0] + d[0] * BL, p[1] + d[1] * BL); t.addBond(cm, oo, 2); clr(t, [cm]); t.atoms.get(ox).chg = 0; clr(t, [ox]); }, ['H₃O⁺']);
  return fin(mech, o);
});

/* ================= LiAlH4 on esters and amides ================= */
function hydride(mech, site, lah) { const ids = mech.add({ atoms: [{ el: lah ? 'Al' : 'B', chg: -1, h: 3 }, { el: 'H', x: -0.82 }], bonds: [[0, 1]] }, site, { dist: 1.75 * BL, att: 2 }); return { b: ids.get(1), h: ids.get(2) }; }
W('hyd_red', (m, o, c) => {
  if (!c.has('LiAlH4')) return null;
  const es = M(m, '[CX3:1](=[OX1:2])[OX2:3][#6:4]'), am = M(m, '[CX3:1](=[OX1:2])[NX3:3]'), ac = M(m, '[CX3:1](=[OX1:2])[OX2H1:3]'), ni = M(m, '[CX2:1]#[NX1:2]');
  if (es.length + am.length + ac.length + ni.length !== 1) return null;
  const mech = new Mech(m);
  if (es.length || ac.length) {
    const [c1, o2, o3] = (es[0] || ac[0]);
    if (ac.length) { const hO = mech.H(o3); const H = hydride(mech, hO, true); mech.step([[[H.b, H.h], [H.h, hO]], [[hO, o3], [o3]]], 'The first H⁻ acts as a base on the acidic O–H: H₂↑ escapes and the carboxylate forms.', ['LiAlH₄']); mech.dropFrag(H.b); }
    let H = hydride(mech, c1, true);
    mech.step([[[H.b, H.h], [H.h, c1]], [[c1, o2], [o2]]], 'Hydride (H⁻) from AlH₄⁻ adds to the acyl carbon; the π pair moves onto O: a tetrahedral intermediate.', ac.length ? [] : ['LiAlH₄']); mech.dropFrag(H.b);
    mech.step([[[o2], [o2, c1]], [[c1, o3], [o3]]], ac.length ? 'The C=O re-forms and O²⁻ (held by aluminium) leaves: an aldehyde.' : 'The C=O re-forms and the alkoxide R′O⁻ leaves: an aldehyde, which is more reactive than the ester.'); mech.dropFrag(o3);
    H = hydride(mech, c1, true);
    mech.step([[[H.b, H.h], [H.h, c1]], [[c1, o2], [o2]]], 'A second H⁻ adds at once: the alkoxide of a 1° alcohol.'); mech.dropFrag(H.b);
    prot(mech, o2, 'Work-up: H₃O⁺ protonates the alkoxide: the 1° alcohol (RCH₂OH).', ['H₃O⁺']);
    return fin(mech, o);
  }
  if (am.length) {
    const [c1, o2, n3] = am[0]; let H = hydride(mech, c1, true);
    mech.step([[[H.b, H.h], [H.h, c1]], [[c1, o2], [o2]]], 'H⁻ adds to the amide carbon; the π pair moves onto O, which binds to aluminium.', ['LiAlH₄']);
    const al = H.b; mech.step([[[o2], [o2, al]]], 'The alkoxide O coordinates to AlH₃ (O–Al), turning it into a leaving group.');
    mech.step([[[n3], [n3, c1]], [[c1, o2], [o2]]], 'The N lone pair pushes in and the O–Al group leaves: an iminium ion (C=N⁺).'); mech.dropFrag(o2);
    H = hydride(mech, c1, true);
    mech.step([[[H.b, H.h], [H.h, c1]], [[c1, n3], [n3]]], 'A second H⁻ adds to the iminium carbon: the amine. The C=O has become CH₂.'); mech.dropFrag(H.b);
    return fin(mech, o);
  }
  const [c1, n2] = ni[0]; let H = hydride(mech, c1, true);
  mech.step([[[H.b, H.h], [H.h, c1]], [[c1, n2], [n2]]], 'H⁻ adds to the nitrile carbon; one π pair moves onto N (an imine anion bound to Al).', ['LiAlH₄']); mech.dropFrag(H.b);
  H = hydride(mech, c1, true);
  mech.step([[[H.b, H.h], [H.h, c1]], [[c1, n2], [n2]]], 'A second H⁻ adds to the C=N carbon: the N carries two negative charges' + ' (bound to aluminium).'); mech.dropFrag(H.b);
  prot(mech, n2, 'Work-up: water protonates the nitrogen…', ['H₂O']); prot(mech, n2, '…twice: the 1° amine R–CH₂NH₂.');
  return fin(mech, o);
});

/* ================= EAS: first substitution with arrows, then "the same again" ================= */
function easArrow(mech, p, eSmi, pickE, leave, text, cond) {
  const m = mech.m; const q = m.nbrs(p).find(n => m.bondBetween(p, n).o === 2 && AR(m).has(n)); if (q == null) return null;
  const s = put(mech, eSmi, p, pickE, { dist: 1.7 * BL }); if (!s) return null; const e = s.at; const arr = [[[p, q], [p, e]]];
  const lv = leave ? leave(mech.m, e) : null; if (lv) arr.push([[e, lv], [lv]]);
  mech.step(arr, text, cond); if (lv != null && E.fragIds(mech.m, lv).length < E.fragIds(mech.m, e).length) mech.dropFrag(lv);
  return { q, e };
}
function easBase(mech, p, q, base, text) { const h = mech.H(p); if (h == null) return false; const ph = P(mech.m, h), pp = P(mech.m, p);
  const b = mech.add(base, h, { dist: 1.45 * BL, dir: unit([ph[0] - pp[0], ph[1] - pp[1]]) }).get(1); mech.step([[[b], [b, h]], [[p, h], [p, q]]], text); mech.dropFrag(b); return true; }
W('eas', (m, o, c) => {
  const base = [...m.atoms.keys()]; const bs = new Set(base); const AS = AR(m);
  const sites = base.filter(a => AS.has(a) && o.product.atoms.has(a) && o.product.nbrs(a).some(n => !bs.has(n) && el(o.product, n) !== 'H'));
  if (!sites.length || base.some(a => !o.product.atoms.has(a) && el(m, a) !== 'H')) return null;
  const p = sites[0]; const gN = o.product.nbrs(p).find(n => !bs.has(n)); const ge = el(o.product, gN);
  const mech = new Mech(m); let r; const dbg = w => { if (root.__mechDebug) (root.__mechWhy = root.__mechWhy || {}).eas_d = w; };
  if (E.HAL.has(ge)) {
    const water = c.any(['HOBr', 'HOCl', 'H2O']) || c.flag('aq');
    r = easArrow(mech, p, ge + ge, isEl(ge), (t, e) => t.nbrs(e).find(n => el(t, n) === ge), `The very electron-rich ring (–OH / –NH₂ push electrons in) attacks ${ge}₂ directly, polarising it; ${ge}⁻ leaves. No Lewis acid is needed: an arenium ion (σ-complex).`, [water ? `${ge}₂ water` : `${ge}₂`]);
    if (!r) return null; if (!easBase(mech, p, r.q, 'O', 'Water removes H⁺ from the sp³ carbon: the ring is aromatic again.')) return null;
  } else if (ge === 'N') {
    r = easArrow(mech, p, '[N+](=O)=O', isEl('N'), (t, e) => t.nbrs(e).find(n => el(t, n) === 'O' && t.bondBetween(e, n).o === 2), 'The ring π electrons attack the nitronium ion NO₂⁺; one N=O π pair moves onto O: an arenium (σ-complex) ion.', ['conc. HNO₃']);
    if (!r) { dbg('easArrow N'); return null; } if (!easBase(mech, p, r.q, 'O', 'Water / HSO₄⁻ removes H⁺: the ring is aromatic again.')) return null;
  } else if (ge === 'S') {
    r = easArrow(mech, p, 'O=S(=O)=O', isEl('S'), (t, e) => t.nbrs(e).find(n => el(t, n) === 'O'), 'SO₃ (from oleum, or hot conc. H₂SO₄) is the electrophile: the ring attacks S and an S=O π pair moves onto O.', ['SO₃ / H₂SO₄']);
    if (!r) return null; if (!easBase(mech, p, r.q, '[O-]S(=O)(=O)O', 'HSO₄⁻ removes H⁺ from the sp³ carbon: the ring is aromatic again (the arenesulphonate).')) return null;
    const so = r.e; const oM = mech.m.nbrs(so).find(n => chg(mech.m, n) < 0); if (oM != null) prot(mech, oM, 'The –SO₃⁻ picks up H⁺: the sulphonic acid. Sulphonation is reversible (steam / dilute acid removes –SO₃H).');
  } else return null;
  if (sites.length > 1) settle(mech, `The same two steps repeat at every other free ortho/para position (${sites.length} substitutions in all).`, t => graftFrom(t, o.product, base, [p]));
  return fin(mech, o);
});


/* ================= batch 2 (Oct 15, run 19b) ================= */
/* :CCl2 from CHCl3 + OH⁻ near `site`; returns the carbene carbon id */
function relayFrag(mech, anchor) { return () => { const ids = E.fragIds(mech.m, anchor); const sub = mech.m.sub(ids); const old = new Map(ids.map(id => [id, P(mech.m, id)]));
  const badB = q => [...q.bonds.values()].some(b => { const x = q.atoms.get(b.a), y = q.atoms.get(b.b); const d = Math.hypot(x.x - y.x, x.y - y.y); return d > 1.6 * BL || d < 0.5 * BL; });
  { try { const { smi, order } = C.toSmiles(sub); const f = C.fromSmilesMB(smi, BL); const fk = [...f.atoms.keys()]; if (fk.length === order.length) order.forEach((id, i) => { const a = sub.atoms.get(id), q = f.atoms.get(fk[i]); if (a && q) { a.x = q.x; a.y = q.y; } }); } catch (e) { } }
  const hv = ids.filter(id => el(mech.m, id) !== 'H'); const c0 = E.centroid(mech.m, hv), c1 = E.centroid(sub, hv); let sxx = 0, sxy = 0;
  for (const id of hv) { const o = old.get(id), n = P(sub, id); const ox = o[0] - c0[0], oy = o[1] - c0[1], nx = n[0] - c1[0], ny = n[1] - c1[1]; sxx += nx * ox + ny * oy; sxy += nx * oy - ny * ox; }
  const t = Math.atan2(sxy, sxx); for (const id of ids) { const n = P(sub, id); const v = E.rot([n[0] - c1[0], n[1] - c1[1]], t); const a = mech.m.atoms.get(id); a.x = c0[0] + v[0]; a.y = c0[1] + v[1]; } }; }
function carbene(mech, site, dist) {
  const ids = mech.add(spec('ClC(Cl)Cl'), site, { att: 2, dist: dist || 2.6 * BL });
  const cm = [...ids.values()].find(id => el(mech.m, id) === 'C'); const h = mech.H(cm); if (h == null) return null;
  const ph = P(mech.m, h), pc = P(mech.m, cm); const b = mech.add('[OH-]', h, { dist: 1.45 * BL, dir: unit([ph[0] - pc[0], ph[1] - pc[1]]) }).get(1);
  mech.step([[[b], [b, h]], [[h, cm], [cm]]], 'OH⁻ removes the H of chloroform (made acidic by three Cl): ⁻CCl₃.', ['CHCl₃', 'KOH']); mech.dropFrag(b);
  const cl = nb(mech.m, cm, 'Cl')[0];
  mech.step([[[cm, cl], [cl]]], 'α-Elimination: Cl⁻ leaves the carbanion carbon: dichlorocarbene :CCl₂, an electron-deficient electrophile.', [], () => { const a = mech.m.atoms.get(cm); a.chg = 0; a.rad = 2; a.hFix = 0; });
  mech.dropFrag(cl); return cm;
}
W('wittig', (m, o, c) => {
  const cs = carbonyls(m); if (!cs.length) return null; const [cc, ox] = cs[0];
  const g = spec(c.has('Ph3P=CHCH3') ? '[CH-](C)[P+](c1ccccc1)(c1ccccc1)c1ccccc1' : '[CH2-][P+](c1ccccc1)(c1ccccc1)c1ccccc1'); if (!g) return null;
  const yc = [...g.atoms.values()].find(a => a.el === 'C' && a.chg === -1); const mech = new Mech(m);
  const map = mech.add(g, cc, { att: yc.id, dist: 1.8 * BL }); const cy = map.get(yc.id); const p = mech.m.nbrs(cy).find(n => el(mech.m, n) === 'P');
  mech.step([[[cy], [cy, cc]], [[cc, ox], [ox]]], 'The ylide carbon (a carbanion held by P⁺) attacks the carbonyl carbon; the π pair moves onto O: a betaine.', ['Ph₃P=CH₂']);
  mech.step([[[ox], [ox, p]]], 'O⁻ bonds to the positive phosphorus: a four-membered oxaphosphetane ring.', [], relayFrag(mech, cc));
  mech.step([[[cc, ox], [cc, cy]], [[cy, p], [p, ox]]], 'The ring falls apart in one step: the C=C forms exactly where C=O was, and Ph₃P=O leaves. The very strong P=O bond drives it.', [], () => { relayFrag(mech, cc)(); E.fixGeom(mech.m, mech.m.clone(), { formed: [], broken: [[cy, p]] }); });
  return fin(mech, o);
});
W('hvz', (m, o, c) => {
  const ac = M(m, '[CX3:1](=[OX1:2])[OX2H1:3]')[0]; if (!ac) return null; const [c1, o2, o3] = ac; const X = c.has('Cl2') ? 'Cl' : 'Br';
  const ca = m.nbrs(c1).find(n => el(m, n) === 'C' && hTot(m, n) > 0 && o.product.nbrs(n).some(q => el(o.product, q) === X && !m.atoms.has(q))); if (ca == null) return null;
  const mech = new Mech(m); let xa = null;
  settle(mech, `Red P + ${X}₂ give P${X}₃, which turns a little of the acid into the acyl ${X === 'Br' ? 'bromide' : 'chloride'}, R–CH₂–CO${X}.`, t => { const a = t.atoms.get(o3); dropH(t, o3); a.el = X; a.hFix = null; xa = o3; }, ['red P', `${X}₂`]);
  settle(mech, 'The acyl halide enolises far more easily than the acid itself (its α-H is more acidic): the enol.', t => { dropH(t, ca); clr(t, [ca]); setO(t, c1, ca, 2); setO(t, c1, o2, 1); });
  const xs = mech.add(X + X, ca, { dist: 1.7 * BL }); const x1 = xs.get(1), x2 = xs.get(2);
  mech.step([[[o2], [o2, c1]], [[c1, ca], [ca, x1]], [[x1, x2], [x2]]], `The enol C=C attacks ${X}₂ (the O lone pair pushes in) and ${X}⁻ leaves.`); mech.dropFrag(x2);
  deprot(mech, o2, `Loss of H⁺: the α-${X === 'Br' ? 'bromo' : 'chloro'} acyl halide.`, `[${X}-]`);
  settle(mech, `It swaps ${X} for OH with a fresh acid molecule (passing the acyl halide on, so P is only needed in catalytic amounts): the α-halo acid.`, t => { const a = t.atoms.get(xa); a.el = 'O'; a.hFix = null; clr(t, [xa]); }, ['H₂O']);
  return fin(mech, o);
});
W('alpha_x', (m, o, c) => {
  const X = c.has('Br2') ? 'Br' : c.has('Cl2') ? 'Cl' : 'I';
  const hs = M(m, '[CX3;$(C([#6])([#6])=O),$([CH1]([#6])=O)](=O)[CX4;!H0]'); if (!hs.length) return null;
  const [cc, ox, ca] = hs.sort((p, q) => hTot(m, p[2]) - hTot(m, q[2]))[0]; const mech = new Mech(m);
  prot(mech, ox, 'Acid protonates the C=O oxygen.', ['H⁺']);
  const h = mech.H(ca); const ph = P(mech.m, h), pa = P(mech.m, ca); const w = mech.add('O', h, { dist: 1.45 * BL, dir: unit([ph[0] - pa[0], ph[1] - pa[1]]) }).get(1);
  mech.step([[[w], [w, h]], [[h, ca], [ca, cc]], [[cc, ox], [ox]]], 'Water removes an α-H; the C=C forms and the C=O⁺ π pair moves onto O: the enol. This step is slow (rate-determining), so the rate does not depend on [X₂].'); mech.dropFrag(w);
  const xs = mech.add(X + X, ca, { dist: 1.7 * BL }); const x1 = xs.get(1), x2 = xs.get(2);
  mech.step([[[ox], [ox, cc]], [[cc, ca], [ca, x1]], [[x1, x2], [x2]]], `The enol C=C attacks ${X}₂ and ${X}⁻ leaves.`, [`${X}₂`]); mech.dropFrag(x2);
  deprot(mech, ox, `Loss of H⁺: the α-halo ketone. The ${X} pulls electrons from the C=O, so a second enolisation is slower: it stops at one ${X}.`);
  return fin(mech, o);
});
W('bisulfite', (m, o) => {
  const cs = carbonyls(m); if (!cs.length) return null; const [cc, ox] = cs[0]; const mech = new Mech(m);
  const g = spec('[O-]S(=O)O'); const sA = [...g.atoms.values()].find(a => a.el === 'S'); const map = mech.add(g, cc, { att: sA.id, dist: 1.8 * BL }); const S = map.get(sA.id);
  const om = mech.m.nbrs(S).find(n => chg(mech.m, n) < 0), oh = mech.m.nbrs(S).find(n => el(mech.m, n) === 'O' && mech.m.bondBetween(S, n).o === 1 && n !== om);
  mech.step([[[S], [S, cc]], [[om], [om, S]], [[cc, ox], [ox]]], 'The sulphur lone pair of HSO₃⁻ (S is the better nucleophile) attacks the carbonyl carbon; the π pair moves onto O.', ['NaHSO₃']);
  const h = mech.m.nbrs(oh).find(n => el(mech.m, n) === 'H') ?? mech.H(oh);
  mech.step([[[ox], [ox, h]], [[h, oh], [oh]]], 'Proton transfer from the –SO₃H to the alkoxide: the α-hydroxy sulphonate salt (crystalline). Dilute acid or base reverses it.');
  return fin(mech, o);
});
W('benzilic', (m, o) => {
  const h = M(m, '[CX3](=O)[CX3]=O')[0]; if (!h) return null; const [c1, o1, c2, o2] = h; if (hTot(m, c1) || hTot(m, c2)) return null;
  const R1 = m.nbrs(c1).find(n => n !== c2 && n !== o1); if (R1 == null) return null;
  for (const [a, oa, b, ob, R] of [[c1, o1, c2, o2, R1], [c2, o2, c1, o1, m.nbrs(c2).find(n => n !== c1 && n !== o2)]]) {
    const mech = new Mech(m); const w = put(mech, '[OH-]', a, isEl('O'), { dist: 1.65 * BL }).at;
    mech.step([[[w], [w, a]], [[a, oa], [oa]]], 'OH⁻ adds to one C=O carbon: a tetrahedral alkoxide.', ['OH⁻']);
    mech.step([[[oa], [oa, a]], [[a, R], [R, b]], [[b, ob], [ob]]], 'The C=O re-forms and pushes a group R across to the neighbouring carbonyl carbon (1,2-shift with its electron pair); that C=O π pair moves onto O.', [], relayFrag(mech, a));
    const hh = mech.H(w); if (hh == null) continue;
    mech.step([[[ob], [ob, hh]], [[hh, w], [w]]], 'Proton transfer from –COOH to the alkoxide: the α-hydroxy carboxylate (benzilic acid after acidification).');
    const r = fin(mech, o); if (key(r.product) === key(o.product)) return r;
  }
  return null;
});
W('diazotise', (m, o, c) => {
  const ar = M(m, '[c:1][NX3H2:2]')[0]; if (!ar || c.heat || (c.T && c.T > 300)) return null; const n = ar[1]; const mech = new Mech(m);
  const s = put(mech, 'N#[O+]', n, isEl('N'), { dist: 1.75 * BL }); const nn = s.at, ox = mech.m.nbrs(nn)[0];
  mech.step([[[n], [n, nn]], [[nn, ox], [ox]]], 'NaNO₂ + HCl give HNO₂, which loses water in acid to give the nitrosonium ion NO⁺. The amine N attacks it.', ['NaNO₂', 'HCl', '0–5 °C']);
  deprot(mech, n, 'Loss of H⁺: the N-nitrosamine, Ar–NH–N=O.');
  settle(mech, 'Tautomerism (H moves from N to O): the diazohydroxide, Ar–N=N–OH.', t => { dropH(t, n); clr(t, [n]); setO(t, n, nn, 2); setO(t, nn, ox, 1); });
  prot(mech, ox, 'Acid protonates the –OH: it can leave as water.');
  mech.step([[[n], [n, nn]], [[nn, ox], [ox]]], 'The N lone pair forms the third N–N bond as water leaves: the arenediazonium ion Ar–N≡N⁺. It is stable only at 0–5 °C (the ring delocalises the charge).'); mech.dropFrag(ox);
  return fin(mech, o);
});
W('nitrile_rx', (m, o, c) => {
  if (c.has('DIBAL') || c.has('SnCl2') || c.has('H2') || c.has('Na') || c.any(['KOH', 'aqKOH'])) return null;
  const ns = M(m, '[CX2:1]#[NX1:2]')[0]; if (!ns) return null; const [c1, n2] = ns; const mech = new Mech(m);
  prot(mech, n2, 'Acid protonates the nitrile N: the carbon becomes a strong electrophile.', ['H₃O⁺', 'Δ']);
  const w = put(mech, 'O', c1, isEl('O'), { dist: 1.65 * BL }).at;
  mech.step([[[w], [w, c1]], [[c1, n2], [n2]]], 'Water attacks the nitrile carbon; one π pair moves onto N.');
  deprot(mech, w, 'Loss of H⁺: an imidic acid, R–C(OH)=NH.');
  settle(mech, 'Tautomerism (H from O to N): the amide R–CONH₂.', t => { dropH(t, w); clr(t, [w, n2]); setO(t, c1, w, 2); setO(t, c1, n2, 1); });
  settle(mech, 'With heat the amide is hydrolysed in the same way (H⁺ on C=O, water adds, NH₃ leaves as NH₄⁺): the carboxylic acid.', t => { const a = t.atoms.get(n2); dropH(t, n2); a.el = 'O'; a.hFix = null; a.chg = 0; clr(t, [n2]); });
  return fin(mech, o);
});
W('amine_rx', (m, o, c) => {
  if (!c.has('CHCl3')) return null; const p = M(m, '[#6:1][NX3H2:2]').filter(x => !E.isCarbonylC(m, x[0]))[0]; if (!p) return null; const n = p[1];
  const mech = new Mech(m); const cm = carbene(mech, n, 2.4 * BL); if (cm == null) return null;
  mech.step([[[n], [n, cm]]], 'The amine N lone pair attacks the empty orbital of the carbene carbon.', [], () => { const a = mech.m.atoms.get(cm); a.rad = 0; a.chg = -1; a.hFix = 0; });
  settle(mech, 'Base removes H⁺ from N and Cl⁻ leaves from C, twice (two more α-eliminations): the isocyanide R–N⁺≡C⁻, foul-smelling.', t => { for (const x of nb(t, cm, 'Cl')) t.removeAtom(x); dropH(t, n); const a = t.atoms.get(n), b = t.atoms.get(cm); a.chg = 1; a.hFix = 0; b.chg = -1; b.hFix = 0; b.rad = 0; setO(t, n, cm, 3); t.atoms.get(n).hFix = 0; t.atoms.get(cm).hFix = 0; });
  return fin(mech, o);
});
W('acylation', (m, o, c) => {
  if (!c.has('Hinsberg')) return null;
  const am = M(m, '[NX3;!$(NC=O);!$(N=*);!$(N#*):1]').filter(([n]) => hTot(m, n) > 0 && !m.bondsOf(n).some(b => el(m, m.other(b, n)) === 'O')); if (!am.length) return null; const n = am[0][0];
  const g = spec('ClS(=O)(=O)c1ccccc1'); const sA = [...g.atoms.values()].find(a => a.el === 'S'); const mech = new Mech(m);
  const map = mech.add(g, n, { att: sA.id, dist: 1.8 * BL }); const S = map.get(sA.id); const cl = mech.m.nbrs(S).find(x => el(mech.m, x) === 'Cl');
  mech.step([[[n], [n, S]], [[S, cl], [cl]]], 'The amine N attacks the sulphur of benzenesulphonyl chloride and Cl⁻ leaves (substitution at S).', ["Hinsberg's reagent"]); mech.dropFrag(cl);
  deprot(mech, n, `Loss of H⁺: the sulphonamide. ${hTot(m, n) >= 2 ? 'From a 1° amine it still has an acidic N–H, so it dissolves in KOH.' : 'From a 2° amine no N–H is left: insoluble in KOH.'}`, '[OH-]');
  return fin(mech, o);
});

/* ================= batch 3 (Oct 15 evening, run 20) ================= */
/* last step drawn without arrows, landing exactly on the plain recipe's product */
function jump(mech, o, text, cond) { const L = mech.m.clone(); const R = o.product.clone(); mech.steps.push({ L, arrows: [], R: R.clone(), text, cond: cond || [] }); mech.m = R; return { steps: mech.steps, product: R }; }
const baseOn = (mech, h, from, base) => { const ph = P(mech.m, h), pa = P(mech.m, from); return mech.add(base || '[OH-]', h, { dist: 1.45 * BL, dir: unit([ph[0] - pa[0], ph[1] - pa[1]]) }).get(1); };
W('wolff', (m, o, c) => {
  if (!c.has('NH2NH2')) return null; const cs = carbonyls(m); if (cs.length !== 1) return null; const [cc, ox] = cs[0];
  const mech = new Mech(m); let nI, nT;
  settle(mech, 'Hydrazine first condenses with the C=O (addition, then loss of water): the hydrazone C=N–NH₂.', t => { const p0 = P(t, ox); t.removeAtom(ox); clr(t, [cc]); nI = t.addAtom('N', p0[0], p0[1]); t.addBond(cc, nI, 2); const d = E.freeDir(t, nI); nT = t.addAtom('N', p0[0] + d[0] * BL, p0[1] + d[1] * BL); t.addBond(nI, nT, 1); }, ['NH₂NH₂']);
  let h = mech.H(nT); let b = baseOn(mech, h, nT);
  mech.step([[[b], [b, h]], [[h, nT], [nT, nI]], [[nI, cc], [cc]]], 'Strong base (KOH in hot glycol) removes an N–H; the electrons shift to make N=N and push the C=N π pair onto carbon: a carbanion.', ['KOH', 'glycol, Δ']); mech.dropFrag(b);
  let ids = mech.add('[H]O', cc, { dist: 1.5 * BL }); let hh = ids.get(1), w = ids.get(2);
  mech.step([[[cc], [cc, hh]], [[hh, w], [w]]], 'The carbanion takes H⁺ from water: R₂CH–N=NH.'); mech.dropFrag(w);
  h = mech.H(nT); b = baseOn(mech, h, nT);
  mech.step([[[b], [b, h]], [[h, nT], [nT, nI]], [[nI, cc], [cc]]], 'OH⁻ removes the last N–H; the C–N bond breaks and N₂ gas leaves (the driving force): a carbanion again.'); mech.dropFrag(b); mech.dropFrag(nI);
  ids = mech.add('[H]O', cc, { dist: 1.5 * BL }); hh = ids.get(1); w = ids.get(2);
  mech.step([[[cc], [cc, hh]], [[hh, w], [w]]], 'Water protonates it: C=O has become CH₂.'); mech.dropFrag(w);
  return fin(mech, o);
});
W('rx_base', (m, o, c) => {
  const hs = E.halides(m); const by = new Map(); for (const [cc, x] of hs) { if (!by.has(cc)) by.set(cc, []); by.get(cc).push(x); }
  const g = [...by].find(([, xs]) => xs.length === 2); if (!g) return null; const [cc, [x1, x2]] = g;
  if (!o.product.nbrs(cc).some(n => el(o.product, n) === 'O' && o.product.bondBetween(cc, n).o === 2)) return null;
  const mech = new Mech(m); const w = put(mech, '[OH-]', cc, isEl('O'), { dist: 1.65 * BL }).at;
  mech.step([[[w], [w, cc]], [[cc, x1], [x1]]], 'OH⁻ replaces the first halogen (substitution).', ['aq. KOH']); mech.dropFrag(x1);
  mech.step([[[w], [w, cc]], [[cc, x2], [x2]]], 'The –OH lone pair pushes in and expels the second halide at once: a protonated carbonyl (the gem-diol route ends at the same place).'); mech.dropFrag(x2);
  deprot(mech, w, `Loss of H⁺: the ${hTot(m, cc) ? 'aldehyde' : 'ketone'}.`, '[OH-]');
  return fin(mech, o);
});
W('benzoin', (m, o) => {
  const q = M(m, '[c][CH1:1]=[O:2]')[0]; if (!q) return null; const [, cc, ox] = q; const mech = new Mech(m);
  const cn = put(mech, '[C-]#N', cc, isEl('C'), { dist: 1.7 * BL }).at;
  mech.step([[[cn], [cn, cc]], [[cc, ox], [ox]]], 'CN⁻ adds to the aldehyde carbon.', ['KCN', 'EtOH']);
  let ids = mech.add('[H]OCC', ox, { dist: 1.5 * BL }); const hO = ids.get(1), so = ids.get(2);
  mech.step([[[ox], [ox, hO]], [[hO, so], [so]]], 'The alkoxide is protonated by ethanol: a cyanohydrin.'); mech.dropFrag(so);
  const h = mech.H(cc); const b = baseOn(mech, h, cc, '[O-]CC');
  mech.step([[[b], [b, h]], [[h, cc], [cc]]], 'Now the old aldehyde H is acidic (the carbanion is stabilised by C≡N and the ring): EtO⁻ removes it. The carbonyl carbon has become a nucleophile (umpolung).'); mech.dropFrag(b);
  const g = spec('O=Cc1ccccc1'); const cA = [...g.atoms.values()].find(a => carbC(a, g)); const map = mech.add(g, cc, { att: cA.id, dist: 1.8 * BL }); const c2 = map.get(cA.id), o2 = dblO(mech.m, c2);
  mech.step([[[cc], [cc, c2]], [[c2, o2], [o2]]], 'The carbanion attacks a second benzaldehyde molecule.');
  mech.step([[[o2], [o2, hO]], [[hO, ox], [ox]]], 'Proton transfer from the first O–H to the new alkoxide.');
  mech.step([[[ox], [ox, cc]], [[cc, cn], [cn]]], 'The O⁻ pushes back to form C=O and CN⁻ is expelled (it is only a catalyst): benzoin, PhCO–CH(OH)Ph.'); mech.dropFrag(cn);
  return fin(mech, o);
});
W('aldol_dehyd', (m, o) => {
  const h0 = M(m, '[CX3:1](=[O:2])[CX4;!H0:3][CX4:4][OX2H1:5]')[0]; if (!h0) return null; const [cc, ox, ca, cb, oh] = h0; const mech = new Mech(m);
  const h = mech.H(ca); const b = baseOn(mech, h, ca);
  mech.step([[[b], [b, h]], [[h, ca], [ca, cc]], [[cc, ox], [ox]]], 'OH⁻ removes the acidic α-H (between the C=O and the β-carbon): the enolate.', ['OH⁻', 'Δ']); mech.dropFrag(b);
  mech.step([[[ox], [ox, cc]], [[cc, ca], [ca, cb]], [[cb, oh], [oh]]], 'The enolate pushes out OH⁻ from the β-carbon (E1cB: elimination from the conjugate base). The new C=C is conjugated with the C=O, which drives it.'); mech.dropFrag(oh);
  return fin(mech, o);
});
W('dibal', (m, o) => {
  const es = M(m, '[#6][CX3:1](=[O:2])[OX2:3][#6]')[0], ni = M(m, '[#6][CX2:1]#[NX1:2]')[0]; const hit = es ? [es[1], es[2], es[3]] : ni ? [ni[1], ni[2]] : null; if (!hit) return null;
  const [c1, x] = hit; const mech = new Mech(m);
  const g = spec('CC(C)C[AlH]CC(C)C'); const al = [...g.atoms.values()].find(a => a.el === 'Al'); for (const a of g.atoms.values()) a.hFix = a.el === 'Al' ? 1 : a.hFix;
  const map = mech.add(g, x, { att: al.id, dist: 1.8 * BL }); const A = map.get(al.id);
  mech.step([[[x], [x, A]]], `The Lewis-acidic aluminium binds the ${es ? 'C=O oxygen' : 'nitrile N'}, making the carbon more electrophilic.`, ['DIBAL-H', '−78 °C']);
  const hA = mech.H(A); mech.step([[[A, hA], [hA, c1]], [[c1, x], [x]]], 'Hydride moves from Al to the carbon (only one H⁻: DIBAL has just one).');
  return jump(mech, o, es ? 'At −78 °C this tetrahedral aluminium complex is stable, so a second hydride never adds. Aqueous work-up collapses it to the aldehyde and the alcohol.' : 'The imine–aluminium complex survives until work-up, where water hydrolyses C=N to the aldehyde.', ['H₂O']);
});
W('r2cd', (m, o, c) => {
  const q = M(m, '[CX3:1](=[O:2])[Cl:3]')[0]; if (!q) return null; const [c1, o2, x] = q; const R = (c.get('R2Cd') || {}).smi || 'C';
  const g = spec(R + '[Cd]' + R); if (!g) return null; const cd = [...g.atoms.values()].find(a => a.el === 'Cd'); const r0 = g.nbrs(cd.id)[0];
  const mech = new Mech(m); const map = mech.add(g, c1, { att: r0, dist: 1.8 * BL }); const rc = map.get(r0), Cd = map.get(cd.id);
  mech.step([[[rc, Cd], [rc, c1]], [[c1, o2], [o2]]], 'The C–Cd bond (much less polar than C–Mg) delivers R to the very reactive acyl carbon.', ['R₂Cd']); mech.dropFrag(Cd);
  mech.step([[[o2], [o2, c1]], [[c1, x], [x]]], 'The C=O re-forms and Cl⁻ leaves: the ketone. R₂Cd is too weak a nucleophile to attack the ketone, so it stops here (a Grignard would go on).'); mech.dropFrag(x);
  return fin(mech, o);
});
W('gattermann_koch', (m, o) => {
  const base = [...m.atoms.keys()]; const bs = new Set(base); const p = base.find(a => AR(m).has(a) && o.product.nbrs(a).some(n => !bs.has(n))); if (p == null) return null;
  const mech = new Mech(m); const r = easArrow(mech, p, '[CH+]=O', isEl('C'), null, 'CO + HCl with AlCl₃/CuCl give the formyl cation [HC=O]⁺; the ring attacks it.', ['CO', 'HCl', 'AlCl₃']); if (!r) return null;
  if (!easBase(mech, p, r.q, '[Cl-]', 'AlCl₄⁻ removes H⁺: the aldehyde.')) return null;
  return fin(mech, o);
});
W('stork', (m, o, c, all, k) => {
  const h = M(m, '[CX3;!H0:1]=[CX3:2][NX3;H0;!$(N[CX3]=O):3]')[0]; if (!h) return null; const [cb, ca, n] = h; const mech = new Mech(m);
  const s = put(mech, 'CI', cb, isEl('C'), { dist: 1.75 * BL }); const rc = s.at, x = mech.m.nbrs(rc).find(q => el(mech.m, q) === 'I');
  mech.step([[[n], [n, ca]], [[ca, cb], [cb, rc]], [[rc, x], [x]]], 'The N lone pair pushes through the C=C: the β-carbon attacks CH₃–I (SN2) and I⁻ leaves. An iminium ion forms.', ['CH₃I']); mech.dropFrag(x);
  return jump(mech, o, 'H₃O⁺ hydrolyses the iminium ion (water adds to C=N⁺, the amine leaves): the α-methyl ketone. Only one alkylation, and no strong base needed.', ['H₃O⁺']);
});
W('perkin', (m, o) => {
  const h0 = M(m, '[CX3H1:1](=[O:2])c')[0]; if (!h0) return null; const [cc, ox] = h0; const mech = new Mech(m);
  const g = spec('CC(=O)OC(C)=O'); const cm = [...g.atoms.values()].find(a => a.el === 'C' && g.hImp(a.id) === 3); const map = mech.add(g, cc, { att: cm.id, dist: 2.1 * BL }); const me = map.get(cm.id);
  const cA = mech.m.nbrs(me).find(q => el(mech.m, q) === 'C'); const oA = dblO(mech.m, cA);
  const hh = mech.H(me); const b = baseOn(mech, hh, me, '[O-]C(C)=O');
  mech.step([[[b], [b, hh]], [[hh, me], [me, cA]], [[cA, oA], [oA]]], 'Acetate (the base) removes an α-H of acetic anhydride: its enolate.', ['CH₃COONa', 'Δ']); mech.dropFrag(b);
  mech.step([[[oA], [oA, cA]], [[cA, me], [me, cc]], [[cc, ox], [ox]]], 'The enolate carbon attacks the aldehyde carbon (an aldol-type addition).');
  return jump(mech, o, 'The alkoxide is acetylated and acetic acid is eliminated (the new C=C is conjugated with the ring); hydrolysis of the anhydride end gives the β-aryl acrylic acid: cinnamic acid.', ['H₂O']);
});
W('pyrazole', (m, o) => {
  const h0 = M(m, '[CX3]([#6])(=O)[CX4;!H0][CX3]([#6])=O')[0]; if (!h0) return null; const [c1, , o1] = h0; const mech = new Mech(m);
  const s = put(mech, 'NN', c1, isEl('N'), { dist: 1.7 * BL }); const n = s.at;
  mech.step([[[n], [n, c1]], [[c1, o1], [o1]]], 'One –NH₂ of hydrazine attacks the first C=O.', ['NH₂NH₂']);
  return jump(mech, o, 'Proton transfers and loss of water give a hydrazone; the second –NH₂ closes onto the other C=O in the same way (5-membered ring), water is lost again, and tautomerism gives the aromatic pyrazole (6 π electrons).');
});

/* ================= batch 4 (Oct 16, run 21) ================= */
W('ozone', (m, o, c) => {
  const al = M(m, '[CX3:1]=[CX3:2]').filter(([a, b]) => !AR(m).has(a) && !AR(m).has(b)); if (al.length !== 1 || M(m, 'C#C').length) return null; const [a, b] = al[0];
  const mech = new Mech(m); const g = spec('[O-][O+]=O'); const ids = [...g.atoms.keys()]; const oT0 = ids.find(id => g.atoms.get(id).chg === -1), oC0 = ids.find(id => g.atoms.get(id).chg === 1), oE0 = ids.find(id => id !== oT0 && id !== oC0);
  const map = mech.add(g, a, { att: oT0, dist: 1.7 * BL }); const oT = map.get(oT0), oC = map.get(oC0), oE = map.get(oE0);
  mech.step([[[oT], [oT, a]], [[a, b], [b, oE]], [[oE, oC], [oC]]], 'Ozone adds across the C=C in one concerted step (1,3-dipolar cycloaddition): the molozonide, a 1,2,3-trioxolane.', ['O₃', '−78 °C'], relayFrag(mech, a));
  mech.step([[[oT, oC], [oC]], [[a, b], [a, oT]], [[oE], [oE, b]]], 'The unstable molozonide falls apart (retro-cycloaddition): a carbonyl compound and a carbonyl oxide (R₂C=O⁺–O⁻).', [], () => { E.fixGeom(mech.m, mech.m.clone(), { formed: [], broken: [[a, b]] }); });
  mech.step([[[oC], [oC, a]], [[a, oT], [oT, b]], [[b, oE], [oE]]], 'The two pieces recombine the other way round: the ozonide (1,2,4-trioxolane), with an O–O bridge and an ether O between the two carbons.', [], relayFrag(mech, a));
  return jump(mech, o, (c.has('H2O2') ? 'Oxidative work-up (H₂O₂) cleaves the ozonide; any aldehyde is oxidised further to the carboxylic acid.' : 'Reductive work-up (Zn/H₂O or Me₂S) cleaves the ozonide: two carbonyl compounds. The C=C is replaced by two C=O.'), [c.has('H2O2') ? 'H₂O₂' : 'Zn, H₂O']);
});
W('diol', (m, o, c) => {
  const al = M(m, '[CX3:1]=[CX3:2]').filter(([a, b]) => !AR(m).has(a) && !AR(m).has(b))[0]; if (!al) return null; const [a, b] = al; const os = c.has('OsO4');
  const g = spec(os ? 'O=[Os](=O)(=O)=O' : 'O=[Mn](=O)(=O)[O-]'); const mt = [...g.atoms.values()].find(x => x.el === (os ? 'Os' : 'Mn')).id; const dO = g.nbrs(mt).filter(n => g.bondBetween(mt, n).o === 2).slice(0, 2);
  const mech = new Mech(m); const map = mech.add(g, a, { att: mt, dist: 2.2 * BL }); const M0 = map.get(mt), o1 = map.get(dO[0]), o2 = map.get(dO[1]);
  mech.step([[[a, b], [b, o2]], [[M0, o1], [o1, a]], [[o2, M0], [M0]]], `${os ? 'OsO₄' : 'MnO₄⁻'} adds to one face of the C=C in a single concerted step: a five-membered cyclic ${os ? 'osmate' : 'manganate'} ester. Both C–O bonds form on the same side.`, [os ? 'OsO₄' : 'cold dil. KMnO₄'], relayFrag(mech, a));
  return jump(mech, o, `Hydrolysis breaks the two O–${os ? 'Os' : 'Mn'} bonds, so both –OH groups end up on the same side: the syn (cis) diol.${os ? '' : ' The purple colour disappears and brown MnO₂ forms (Baeyer test).'}`, [os ? 'NaHSO₃, H₂O' : 'H₂O']);
});
W('curtius_lossen', (m, o, c, all, k) => {
  const stages = (all || []).length ? all : [c]; if (!stages.some(s => s.has && s.has('NaN3'))) return null;
  const h0 = M(m, '[#6:1][CX3:2](=[O:3])[Cl,Br,I:4]')[0]; if (!h0) return null; const [r, c1, o2, x] = h0; const mech = new Mech(m);
  const g = spec('[N-]=[N+]=[N-]'); const n10 = [...g.atoms.values()].find(q => q.chg === -1).id; const map = mech.add(g, c1, { att: n10, dist: 1.7 * BL }); const n1 = map.get(n10); const n2 = mech.m.nbrs(n1)[0]; const n3 = mech.m.nbrs(n2).find(q => q !== n1);
  addElim(mech, n1, c1, o2, x, 'Azide (N₃⁻) adds to the acyl carbon.', 'The C=O re-forms and Cl⁻ leaves: the acyl azide, R–CO–N₃.', ['NaN₃']); mech.dropFrag(x);
  settle(mech, 'Drawn in its other resonance form, R–CO–N⁻–N⁺≡N, the N₂ leaving group is easy to see.', t => { setO(t, n1, n2, 1); setO(t, n2, n3, 3); const A = t.atoms.get(n1), B = t.atoms.get(n2), D = t.atoms.get(n3); A.chg = -1; B.chg = 1; D.chg = 0; A.hFix = 0; B.hFix = 0; D.hFix = 0; });
  mech.step([[[c1, r], [r, n1]], [[n1, n2], [n2]], [[n1], [n1, c1]]], 'On heating, R moves from C to N as N₂ leaves, all in one step: the isocyanate R–N=C=O. The migrating group keeps its configuration.', ['Δ'], groupPost(mech, r, c1, n1)); mech.dropFrag(n2);
  return jump(mech, o, 'Water adds to the isocyanate: a carbamic acid, R–NH–COOH, which loses CO₂ at once. R–NH₂ has one carbon fewer than the acid chloride.', ['H₂O']);
});
W('diazotise', (m, o, c) => {
  const al = M(m, '[CX4:1][NX3H2:2]')[0], sec = M(m, '[CX4][NX3H1:1][CX4]')[0]; if (M(m, '[c][NX3H2]').length) return null;
  const n = al ? al[1] : sec ? sec[1] : null; if (n == null) return null; const mech = new Mech(m);
  const s = put(mech, 'N#[O+]', n, isEl('N'), { dist: 1.75 * BL }); const nn = s.at, ox = mech.m.nbrs(nn)[0];
  mech.step([[[n], [n, nn]], [[nn, ox], [ox]]], 'HNO₂ (from NaNO₂ + HCl) gives the nitrosonium ion NO⁺; the amine N attacks it.', ['NaNO₂', 'HCl']);
  deprot(mech, n, al ? 'Loss of H⁺: the N-nitrosamine.' : 'Loss of H⁺: the N-nitrosamine, R₂N–N=O, a yellow oil (Liebermann nitroso test). 2° amines stop here: there is no second N–H.');
  if (!al) return fin(mech, o);
  const cc = al[0];
  settle(mech, 'Tautomerism (H from N to O): R–N=N–OH.', t => { dropH(t, n); clr(t, [n]); setO(t, n, nn, 2); setO(t, nn, ox, 1); });
  prot(mech, ox, 'Acid protonates the –OH.');
  mech.step([[[n], [n, nn]], [[nn, ox], [ox]]], 'Water leaves: the alkyl diazonium ion R–N≡N⁺.'); mech.dropFrag(ox);
  mech.step([[[cc, n], [n]]], 'Unlike the aryl case, nothing stabilises it: N₂ gas bubbles off at once (the test for 1° aliphatic amines) and a carbocation is left.'); mech.dropFrag(n);
  const w = put(mech, 'O', cc, isEl('O'), { dist: 1.65 * BL }).at;
  mech.step([[[w], [w, cc]]], 'Water captures the carbocation. (Rearranged alcohols, alkenes and alkyl chlorides also form, which is why this is not a useful synthesis.)');
  deprot(mech, w, 'Loss of H⁺: the alcohol.');
  return fin(mech, o);
});
W('vic_dehalo', (m, o, c) => {
  const h0 = M(m, '[CX4:1]([Cl,Br,I:2])[CX4:3][Cl,Br,I:4]')[0]; if (!h0) return null; const [c1, x1, c2, x2] = h0; const zn = c.has('Zn'); const mech = new Mech(m);
  const s = put(mech, zn ? '[Zn]' : '[I-]', x1, a => a.el === (zn ? 'Zn' : 'I'), { dist: 1.6 * BL }); const nu = s.at;
  mech.step([[[nu], [nu, x1]], [[x1, c1], [c1, c2]], [[c2, x2], [x2]]], zn ? 'Zn gives electrons to one halogen; the C–X electrons form the C=C as the second X⁻ leaves from the anti position (anti elimination).' : 'I⁻ attacks one halogen (it is the halogen, not the carbon, that is attacked); the C–X electrons form the C=C as the other X⁻ leaves anti-periplanar.', [zn ? 'Zn' : 'NaI, acetone']);
  return fin(mech, o);
});
W('rs_nu', (m, o, c) => {
  const R = /ch3|me/i.test((c.get('RS-') || {}).txt || '') ? 'C' : 'CC'; const hs = E.halides(m); const mech = new Mech(m);
  if (hs.length) { const [cc, x] = hs[0]; const sN = put(mech, '[S-]' + R, cc, isEl('S'), { dist: 1.7 * BL }).at; mech.step([[[sN], [sN, cc]], [[cc, x], [x]]], 'RS⁻ attacks the carbon from the back side as X⁻ leaves (SN2): a thioether. Soft nucleophile, soft sp³ carbon.', ['RS⁻']); mech.dropFrag(x); return fin(mech, o); }
  const q = M(m, '[CX3:1](=[O:2])[Cl,Br,I:3]')[0]; if (!q) return null; const [c1, o2, x] = q; const sN = put(mech, '[S-]' + R, c1, isEl('S'), { dist: 1.7 * BL }).at;
  addElim(mech, sN, c1, o2, x, 'RS⁻ adds to the acyl carbon.', 'The C=O re-forms and the halide leaves: the thioester.', ['RS⁻']); mech.dropFrag(x); return fin(mech, o);
});
W('acylx_alkoxide', (m, o, c) => {
  const q = M(m, '[CX3:1](=[O:2])[Cl,Br,I:3]')[0]; if (!q) return null; const [c1, o2, x] = q; const R = c.has('NaOMe') ? 'C' : c.has('NaOEt') ? 'CC' : ((c.get('NaOR') || {}).smi || 'C');
  const mech = new Mech(m); const oN = put(mech, '[O-]' + R, c1, isEl('O'), { dist: 1.7 * BL }).at;
  addElim(mech, oN, c1, o2, x, 'The alkoxide adds to the very reactive acyl C=O.', 'The C=O re-forms and X⁻ leaves: the ester. Acyl substitution is much faster than SN2 at any sp³ C–X in the molecule.', ['RO⁻']); mech.dropFrag(x); return fin(mech, o);
});
W('acid_nh3', (m, o, c) => {
  const q = M(m, '[CX3:1](=[O:2])[OX2H1:3]'); if (q.length !== 1) return null; const [c1, o2, o3] = q[0]; const mech = new Mech(m);
  const h = mech.H(o3); const nN = baseOn(mech, h, o3, 'N');
  mech.step([[[nN], [nN, h]], [[h, o3], [o3]]], 'Acid–base comes first: NH₃ takes the acidic H⁺. Ammonium carboxylate, RCOO⁻ NH₄⁺.', ['NH₃']);
  if (!c.heat) { mech.dropFrag(nN); return fin(mech, o); }
  mech.dropFrag(nN);
  return jump(mech, o, 'Strong heating: NH₃ is released again and attacks the acyl carbon; the tetrahedral intermediate loses water. RCONH₂, the amide.', ['Δ']);
});
W('acid_heat', (m, o, c) => {
  const ac = M(m, '[CX3:1](=[O:2])[OX2H1:3]'); if (ac.length !== 2) return null; const [ci, , oi] = ac[0], [cj, oj, oj2] = ac[1];
  if (!o.product.bondBetween(oi, cj) && !o.product.bondBetween(oj2, ci)) return null;
  const [nuO, tc, to, lo] = o.product.bondBetween(oi, cj) ? [oi, cj, oj, oj2] : [oj2, ci, ac[0][1], oi];
  const mech = new Mech(m);
  mech.step([[[nuO], [nuO, tc]], [[tc, to], [to]]], 'The –OH oxygen of one –COOH attacks the C=O carbon of the other (intramolecular, five- or six-membered ring).', [c.has('P2O5') ? 'P₂O₅' : 'Δ'], relayFrag(mech, tc));
  return jump(mech, o, 'Proton transfers, then loss of water from the tetrahedral carbon: the cyclic anhydride.');
});
W('intra_aldol', (m, o, c) => {
  const pairs = M(m, '[CX3:1]=[OX1:2]'); if (pairs.length < 2) return null; const pa = new Set([...o.product.atoms.keys()]);
  let hit = null; for (const [c1] of pairs) for (const ca of m.nbrs(c1).filter(x => el(m, x) === 'C' && hTot(m, x) > 0)) for (const [c2, o2] of pairs) if (c2 !== c1 && o.product.bondBetween(ca, c2) && !m.bondBetween(ca, c2)) hit = { c1, ca, c2, o2, o1: dblO(m, c1) };
  if (!hit) return null; const { c1, ca, c2, o2, o1 } = hit; const mech = new Mech(m);
  const h = mech.H(ca); const b = baseOn(mech, h, ca);
  mech.step([[[b], [b, h]], [[h, ca], [ca, c1]], [[c1, o1], [o1]]], 'OH⁻ removes an α-H: the enolate.', ['OH⁻']); mech.dropFrag(b);
  mech.step([[[o1], [o1, c1]], [[c1, ca], [ca, c2]], [[c2, o2], [o2]]], 'The enolate carbon attacks the other C=O of the same molecule: a five- or six-membered ring closes (the strain-free ring wins).', [], relayFrag(mech, ca));
  return jump(mech, o, o.product.bondBetween(ca, c2).o === 2 ? 'Protonation, then E1cB dehydration on heating: the cyclic α,β-unsaturated ketone.' : 'Water protonates the alkoxide: the cyclic β-hydroxy ketone.', c.heat ? ['Δ'] : []);
});

/* ================= Cope elimination (new recipes, FRONT) ================= */
const front = r => E.REC.unshift(r);
const amineOxide = m => M(m, '[NX4+:1][OX1-:2]').filter(([n]) => m.nbrs(n).filter(x => el(m, x) === 'C').length === 3);
front({ id: 'cope_ox', name: '3° amine → amine oxide', fam: 'Amine', canon: 'H₂O₂', db: ['Cope elimination'],
  test: c => c.any(['H2O2', 'mCPBA']) && !c.any(['KOH', 'aqKOH', 'NaOH']),
  run(m) {
    const am = M(m, '[NX3;H0;!$(N[#6]=[O,S,N]);!$(N-a);!$(N=*);!$(N#*):1]').filter(([n]) => m.nbrs(n).filter(x => el(m, x) === 'C').length === 3); if (!am.length) return null; const n = am[0][0];
    const mech = new Mech(m); const ids = mech.add('OO', n, { dist: 1.7 * BL }); const o1 = ids.get(1), o2 = ids.get(2);
    mech.step([[[n], [n, o1]], [[o1, o2], [o2]]], 'The N lone pair attacks one O of H₂O₂ (or of a peracid); the weak O–O bond breaks and OH⁻ (RCOO⁻) leaves.', ['H₂O₂']); mech.dropFrag(o2);
    deprot(mech, o1, 'Loss of H⁺ from the O: the amine oxide R₃N⁺–O⁻.', '[OH-]');
    return E.done(mech, { title: '3° amine + H₂O₂ → amine oxide (N-oxide)', why: ['Only 3° amines give stable N-oxides', 'Heating the N-oxide gives the Cope elimination'], tags: ['ox'] });
  } });
front({ id: 'cope_elim', name: 'Cope elimination', fam: 'Amine', canon: 'Δ (amine oxide)', db: ['Cope elimination'],
  test: c => c.heat || (c.T != null && c.T >= 370) || c.flag('heat'),
  run(m) {
    const ao = amineOxide(m); if (!ao.length) return null; const [n, ox] = ao[0];
    let best = null;
    for (const ca of m.nbrs(n).filter(x => el(m, x) === 'C')) for (const cb of m.nbrs(ca).filter(x => el(m, x) === 'C' && x !== n)) {
      if (hTot(m, cb) < 1 || m.bondsOf(cb).some(b => b.o > 1) || AR(m).has(cb)) continue;
      const sc = hTot(m, cb) * 2 - E.alkeneScore(m, ca, cb, AR(m)) * 0.5; if (!best || sc > best.sc) best = { ca, cb, sc }; }
    if (!best) return E.out({ title: 'No β-H: no Cope elimination', product: null, none: 'The amine oxide needs a β-hydrogen (on a carbon next to the C–N carbon) for the cyclic syn elimination.', why: [] });
    const { ca, cb } = best; const mech = new Mech(m); const h = mech.H(cb, new Set([n])); if (h == null) return null;
    mech.step([[[ox], [ox, h]], [[cb, h], [cb, ca]], [[ca, n], [n]]], 'One concerted step through a planar five-membered cyclic transition state: the N–O⁻ oxygen pulls off a β-H on the same side (syn), the C–H electrons form the C=C and the C–N bond breaks. No external base is needed.', ['Δ (≈ 100–150 °C)']);
    const fr = mech.m.frags(); const hn = fr.find(f => f.includes(n)), al = fr.find(f => f.includes(ca));
    const res = E.done(mech, { title: 'Cope elimination: amine oxide → alkene + R₂N–OH', extra: hn ? [E.cleanProduct(mech.m, hn)] : [], why: ['Syn (cis) elimination, intramolecular (Ei)', 'Gives mainly the less substituted (Hofmann-type) alkene', 'By-product: an N,N-dialkylhydroxylamine'], tags: ['Ei', 'syn'] });
    if (res && al) res.product = E.cleanProduct(mech.m, al); return res;
  } });

const BASES10 = ['KOH', 'aqKOH', 'alcKOH', 'NaOEt', 'NaOMe', 'KOtBu', 'NaNH2', 'LDA', 'NaH'];
front({ id: 'enolate_alkyl', name: 'Enolate alkylation', fam: 'Carbonyl', canon: 'base, CH₃I', db: ['Enolate alkylation'],
  test: c => c.has('RX') && c.any(BASES10),
  run(m, c) {
    if (M(m, '[CX3](=O)[CX4][CX3]=O').length || E.halides(m).length) return null;
    const cs = carbonyls(m); if (!cs.length) return null; const [cc, ox] = cs[0];
    const al = m.nbrs(cc).filter(n => el(m, n) === 'C' && !AR(m).has(n) && hTot(m, n) > 0 && m.bondsOf(n).every(b => b.o === 1)).sort((p, q) => hTot(m, q) - hTot(m, p)); if (!al.length) return null; const ca = al[0];
    const R = (c.get('RX') || {}).smi || 'C'; const X = (c.get('RX') || {}).X || 'I';
    const mech = new Mech(m); const h = mech.H(ca); const b = baseOn(mech, h, ca);
    mech.step([[[b], [b, h]], [[h, ca], [ca, cc]], [[cc, ox], [ox]]], 'The base removes an α-H: the enolate.', ['OH⁻']); mech.dropFrag(b);
    const s = put(mech, R + X, ca, a => a.el === 'C' && s0(a), { dist: 1.75 * BL }); function s0() { return true; } if (!s) return null; const rc = s.at; const x = mech.m.nbrs(rc).find(q => el(mech.m, q) === X); if (x == null) return null;
    mech.step([[[ox], [ox, cc]], [[cc, ca], [ca, rc]], [[rc, x], [x]]], `The enolate carbon attacks the alkyl halide (SN2) and ${X}⁻ leaves: an α-alkyl carbonyl. (With a simple aldehyde, aldol competes; in practice LDA at −78 °C is used.)`, ['R–X']); mech.dropFrag(x);
    return E.done(mech, { title: 'Enolate alkylation → α-alkyl carbonyl', why: ['The α-carbon is the nucleophile', 'SN2: works best with CH₃–X and 1° R–X'], tags: ['enolate', 'SN2'] });
  } });
front({ id: 'claisen_formyl', name: 'Crossed Claisen with ethyl formate', fam: 'Carbonyl', canon: 'HCOOC₂H₅, NaOEt', db: ['Crossed Claisen with ethyl formate'],
  test: c => c.has('HCOOEt') && c.any(['NaOEt', 'NaOMe', 'NaH', 'KOtBu', 'NaNH2']),
  run(m) {
    const cs = carbonyls(m); if (!cs.length) return null; const [cc, ox] = cs[0];
    const al = m.nbrs(cc).filter(n => el(m, n) === 'C' && !AR(m).has(n) && hTot(m, n) > 0 && m.bondsOf(n).every(b => b.o === 1)).sort((p, q) => hTot(m, p) - hTot(m, q)); if (!al.length) return null; const ca = al[0];
    const mech = new Mech(m); const h = mech.H(ca); const b = baseOn(mech, h, ca, '[O-]CC');
    mech.step([[[b], [b, h]], [[h, ca], [ca, cc]], [[cc, ox], [ox]]], 'EtO⁻ removes an α-H of the ketone: the enolate.', ['NaOEt']); mech.dropFrag(b);
    const g = spec('O=COCC'); const cF = [...g.atoms.values()].find(a => carbC(a, g)); const map = mech.add(g, ca, { att: cF.id, dist: 1.8 * BL }); const cf = map.get(cF.id), of = dblO(mech.m, cf), oe = mech.m.nbrs(cf).find(q => el(mech.m, q) === 'O' && q !== of);
    mech.step([[[ox], [ox, cc]], [[cc, ca], [ca, cf]], [[cf, of], [of]]], 'The enolate attacks the C=O of ethyl formate (it has no α-H of its own, so it can only be the electrophile).');
    mech.step([[[of], [of, cf]], [[cf, oe], [oe]]], 'The C=O re-forms and EtO⁻ leaves: the α-formyl ketone (a 1,3-dicarbonyl).'); mech.dropFrag(oe);
    return E.done(mech, { title: 'Crossed Claisen: ketone + HCOOEt → α-formyl ketone', why: ['Ethyl formate has no α-H: a clean crossed Claisen', 'The product (a β-keto aldehyde) is deprotonated by EtO⁻, which drives the equilibrium'], tags: ['enolate', 'cond'] });
  } });
front({ id: 'chugaev', name: 'Chugaev elimination', fam: 'Alcohol', canon: 'Na, CS₂, CH₃I, Δ', db: ['Chugaev elimination'],
  test: c => /cs2|cs₂/i.test(c.raw || '') && c.has('RX'),
  run(m) {
    const al = M(m, '[CX4:1][OX2H1:2]').filter(([cc]) => !AR(m).has(cc)); if (!al.length) return null; const [cc, o] = al[0];
    const bs = E.betaHs(m, cc, o).sort((p, q) => hTot(m, q) - hTot(m, p)); if (!bs.length) return E.out({ title: 'No β-H: no Chugaev elimination', product: null, none: 'The xanthate needs a β-H for the cyclic syn elimination.', why: [] });
    const bC = bs[0]; const mech = new Mech(m); let cs, sT, sD;
    settle(mech, 'Na gives the alkoxide; it adds to CS₂ and the S⁻ is methylated by CH₃I: the xanthate ester R–O–C(=S)–SCH₃.', t => { dropH(t, o); const p0 = P(t, o), d = E.freeDir(t, o); cs = t.addAtom('C', p0[0] + d[0] * BL, p0[1] + d[1] * BL); t.addBond(o, cs, 1); clr(t, [o]); sD = t.addAtom('S', p0[0] + d[0] * 2 * BL, p0[1] + d[1] * 2 * BL - BL * .4); t.addBond(cs, sD, 2); sT = t.addAtom('S', p0[0] + d[0] * 2 * BL, p0[1] + d[1] * 2 * BL + BL * .6); t.addBond(cs, sT, 1); const me = t.addAtom('C', p0[0] + d[0] * 3 * BL, p0[1] + d[1] * 3 * BL + BL * .6); t.addBond(sT, me, 1); E.tidyLongBonds(t); }, ['Na', 'CS₂', 'CH₃I']);
    const h = mech.H(bC, new Set([o]));
    mech.step([[[sD], [sD, h]], [[bC, h], [bC, cc]], [[cc, o], [o]]], 'On heating, the C=S sulphur pulls off a β-H on the same side through a six-membered cyclic transition state (syn, Ei); the C=C forms and the C–O bond breaks. No rearrangement.', ['Δ']);
    const fr = mech.m.frags(); const alk = fr.find(f => f.includes(cc)); const R = E.cleanProduct(mech.m, alk);
    return { title: 'Chugaev elimination: alcohol → alkene (syn, via xanthate)', product: R, steps: mech.steps, extra: [C.fromSmilesMB('O=C=S', BL), C.fromSmilesMB('CS', BL)], why: ['Syn (cis) elimination, no carbocation, so no rearrangement', 'By-products COS and CH₃SH'], tags: ['Ei'], arrowed: true };
  } });
front({ id: 'grignard_s', name: 'Grignard + S / O₂', fam: 'Organometallic', canon: 'S then H₃O⁺', db: ['Grignard + S / O2'],
  test: c => c.has('S8'),
  run(m, c) {
    const g = M(m, '[#6:1][Mg:2]')[0]; if (!g) return null; const [rc, mg] = g; const X = 'S'; const mech = new Mech(m);
    const s = put(mech, X === 'S' ? '[S]' : 'O=O', rc, isEl(X), { dist: 1.7 * BL }); const a = s.at;
    if (X === 'S') { mech.m.atoms.get(a).hFix = 0;
      mech.step([[[rc, mg], [rc, a]]], 'The Grignard carbon attacks sulphur: R–S–MgX (a thiolate salt).', ['S']); }
    else { const b2 = mech.m.nbrs(a)[0]; mech.step([[[rc, mg], [rc, a]], [[a, b2], [b2]]], 'The Grignard carbon attacks O₂; after a second R–MgX and work-up, R–O–MgX gives the alcohol.', ['O₂']); mech.dropFrag(b2); }
    mech.dropFrag(mg);
    return jump(mech, { product: (() => { const t = mech.m.clone(); const at = t.atoms.get(a); at.chg = 0; at.hFix = null; at.rad = 0; return E.cleanProduct(t, E.fragIds(t, rc)); })() }, `Work-up with H₃O⁺ protonates it: ${X === 'S' ? 'the thiol R–SH' : 'the alcohol R–OH'}.`, ['H₃O⁺']) && E.done(mech, { title: X === 'S' ? 'Grignard + sulphur → thiol' : 'Grignard + O₂ → alcohol', why: [X === 'S' ? 'Sulphur inserts into the C–Mg bond' : 'Usually avoided: it is why Grignards are kept under N₂'], tags: ['organometallic'] });
  } });
front({ id: 'liebermann', name: "Liebermann's nitroso test", fam: 'Phenol', canon: 'NaNO₂, conc. H₂SO₄, then NaOH', db: ["Liebermann nitroso test"],
  test: c => c.has('HNO2') && c.has('H2SO4'),
  run(m) {
    if (key(m) !== key(C.fromSmilesMB('Oc1ccccc1', BL))) return null;
    const R = C.fromSmilesMB('O=C1C=CC(=Nc2ccc(O)cc2)C=C1', BL); for (const a of R.atoms.values()) a.hFix = null;
    return E.out({ title: "Liebermann's nitroso test: red, then deep blue/green in NaOH", product: R, steps: [{ L: m.clone(), arrows: [], R: R.clone(), text: 'HNO₂ nitrosates phenol at para (p-nitrosophenol, which tautomerises to the quinone monoxime). In conc. H₂SO₄ this condenses with a second phenol to the red indophenol; NaOH turns it deep blue (its anion).', cond: ['NaNO₂', 'conc. H₂SO₄'] }], why: ['Test for phenols (and for 2° amines via their N-nitrosamines)', 'Red → blue/green on adding alkali'], tags: ['test'] });
  } });

/* Tollens' / Fehling's take priority over imine (NH₃) and aldol (OH⁻) when they are named: the reagent mixture contains both */
front({ id: 'tollens_first', name: "Tollens' / Fehling's test", fam: 'Carbonyl', canon: "Tollens' reagent", db: ["Tollens' test", "Fehling's / Benedict's test"],
  test: c => c.any(['Tollens', 'Fehling']),
  run(m, c, all, k) { if (M(m, '[CH1]#[#6]').length) return null; const r = E.REC.find(x => x.id === 'tollens'); return r ? r.run(m, c, all, k) : null; } });

/* ---------- carbocation rearrangement (1,2-H, alkyl, aryl shifts; ring expansion / contraction) ---------- */
front({ id: 'cation_shift', name: 'Carbocation rearrangement', fam: 'Intermediate', canon: 'rearrange', db: ['1,2-hydride shift', '1,2-methyl shift', '1,2-aryl shift', 'Ring expansion', 'Cyclopropylcarbinyl–cyclobutyl rearrangement', 'Ring contraction'],
  test: c => c.has('rearr'),
  run(m) {
    const cat = [...m.atoms.values()].find(a => a.el === 'C' && a.chg === 1); if (!cat) return null; const c0 = cat.id; const A0 = AR(m);
    const base = E.cationScore(m, c0, A0); let best = null;
    for (const b of m.nbrs(c0).filter(n => el(m, n) === 'C' && !A0.has(n))) {
      const hasOH = m.nbrs(b).some(n => el(m, n) === 'O' && hTot(m, n) > 0);
      const groups = m.nbrs(b).filter(g => g !== c0 && el(m, g) !== 'O' && el(m, g) !== 'N');
      const cand = groups.map(g => ({ g, kind: el(m, g) === 'H' ? 'H' : A0.has(g) ? 'aryl' : 'C' }));
      if (hTot(m, b) > 0 && !cand.some(x => x.kind === 'H')) cand.push({ g: null, kind: 'H' });
      for (const x of cand) {
        const t = m.clone(); let g = x.g; if (x.kind === 'H') { g = E.expH(t, b); if (g == null) continue; }
        const ring = x.kind !== 'H' && E.ringSizeOf(t, b, g) && !E.ringSizeOf(t, c0, b) ? E.ringSizeOf(t, b, g) : 0; // migrating bond is in a ring that does not contain c0 → expansion
        const ringC = x.kind !== 'H' && E.ringSizeOf(t, c0, b) && E.ringSizeOf(t, b, g) === E.ringSizeOf(t, c0, b) ? E.ringSizeOf(t, c0, b) : 0; // c0 and b in the same ring and g in it → contraction
        t.removeBond(t.bondBetween(b, g).id); t.addBond(c0, g, 1); t.atoms.get(c0).chg = 0; t.atoms.get(b).chg = 1; clr(t, [c0, b]);
        let sc = E.cationScore(t, b, AR(t)) + (ring === 3 || ring === 4 ? 3 : 0) + (x.kind === 'H' ? 0.3 : x.kind === 'aryl' ? 0.2 : 0) - (ringC && !hasOH ? 3 : 0);
        if (hasOH) sc = 6 + (x.kind === 'aryl' ? 2 : ringC ? 1.8 : x.kind === 'H' ? 1.5 : 1); // the new cation is oxygen-stabilised either way: migratory aptitude decides (aryl > ring bond > H > alkyl)
        if (sc > base + 0.5 && (!best || sc > best.sc)) best = { sc, b, x, ring, ringC, hasOH };
      }
    }
    if (!best) return E.out({ title: 'No better carbocation by a 1,2-shift', product: null, none: 'No neighbouring H, alkyl or aryl shift gives a more stable carbocation here, so it does not rearrange.', why: [] });
    const { b, x, ring, ringC, hasOH } = best; const mech = new Mech(m); let g = x.g; if (x.kind === 'H') g = mech.H(b); if (g == null) return null;
    const what = x.kind === 'H' ? '1,2-hydride shift' : x.kind === 'aryl' ? '1,2-aryl shift (via a bridged phenonium ion)' : ring ? `ring expansion (${ring}- to ${ring + 1}-membered ring)` : ringC ? `ring contraction (${ringC}- to ${ringC - 1}-membered ring)` : '1,2-alkyl (methyl) shift';
    mech.step([[[b, g], [g, c0]]], `${what[0].toUpperCase() + what.slice(1)}: the ${x.kind === 'H' ? 'H' : 'group'} moves with its bonding pair to the neighbouring cation carbon, so the + charge moves to the carbon it left. ${ring ? 'Relieving the ring strain is the driving force.' : 'The new carbocation is more stable.'}`, [], x.kind === 'H' ? null : relayFrag(mech, c0));
    if (ringC && hasOH) { const o = mech.m.nbrs(b).find(n => el(mech.m, n) === 'O'); mech.step([[[o], [o, b]]], 'The –OH lone pair stabilises the new cation: a protonated carbonyl (semipinacol). Loss of H⁺ gives the ring-contracted aldehyde or ketone.'); }
    return E.done(mech, { title: `Carbocation rearrangement: ${what}`, why: ['Carbocations rearrange whenever a 1,2-shift gives a more stable ion (3° > 2° > 1°, resonance, less strain)', 'Migratory aptitude roughly: H ≈ aryl > alkyl'], tags: ['1,2-shift', 'rearr'] });
  } });
/* carbocation + nucleophile */
front({ id: 'cation_capture', name: 'Carbocation + nucleophile', fam: 'Intermediate', canon: 'Nu⁻', db: [],
  test: c => c.any(['NaI', 'AgF', 'NaSH', 'NaN3']),
  run(m, c) {
    const cat = [...m.atoms.values()].find(a => a.el === 'C' && a.chg === 1); if (!cat) return null;
    const sp = ['NaI', 'AgF', 'NaSH', 'NaN3'].find(k => c.has(k)); const nu = { NaI: '[I-]', AgF: '[F-]', NaSH: '[SH-]', NaN3: '[N-]=[N+]=[N-]' }[sp];
    const mech = new Mech(m); const s = put(mech, nu, cat.id, a => a.chg === -1, { dist: 1.7 * BL }); if (!s) return null;
    mech.step([[[s.at], [s.at, cat.id]]], 'The nucleophile gives a lone pair to the empty p orbital of the carbocation.');
    return E.done(mech, { title: 'Carbocation captured by the nucleophile', why: ['The last step of SN1'], tags: ['sub'] });
  } });
/* ---------- tautomerism: 1,3-H shift (keto→enol, imine→enamine, nitroso→oxime, nitro→aci-nitro); dienone→phenol ---------- */
front({ id: 'tautomer', name: 'Tautomerism', fam: 'Isomerism', canon: 'tautomerise', db: ['Acid-catalysed keto–enol tautomerism', 'Enol-favoured tautomerism', 'Nitro–aci-nitro tautomerism', 'Nitroso–oxime tautomerism', 'Imine–enamine tautomerism'],
  test: c => c.has('taut'),
  run(m) {
    const A0 = AR(m);
    /* cyclohexadienone → phenol */
    const dn = M(m, '[CX3:1]1(=[OX1:2])[#6]=[#6][CX4;!H0:3][#6]=[#6]1')[0] || M(m, '[CX3:1]1(=[OX1:2])[CX4;!H0:3][#6]=[#6][#6]=[#6]1')[0];
    if (dn) { const [cc, ox, ch] = dn; const ring = C.rings(m, 6).find(r => r.length === 6 && r.includes(cc) && r.includes(ch)); if (ring) { const mech = new Mech(m);
      prot(mech, ox, 'Acid (or base) catalysis: the C=O oxygen is protonated.', ['H⁺ (cat.)']);
      settle(mech, 'Loss of the H⁺ from the sp³ ring carbon gives a fully conjugated ring of six π electrons: the phenol. Aromaticity makes the enol form overwhelmingly favoured here.', t => { dropH(t, ch); const i = ring.indexOf(cc); for (let q = 0; q < 6; q++) { const bd = t.bondBetween(ring[(i + q) % 6], ring[(i + q + 1) % 6]); if (bd) bd.o = q % 2 ? 2 : 1; } setO(t, cc, ox, 1); t.atoms.get(ox).chg = 0; clr(t, [...ring, ox]); });
      return E.done(mech, { title: 'Tautomerism: cyclohexadienone → phenol (aromatic enol wins)', why: ['Normally the keto form is far more stable; aromaticity reverses that'], tags: ['tautomerism'] }); } }
    /* 1,3-shift: H–Cα–Y=A → Cα=Y–A–H */
    const pats = [['[CX4;!H0:1][CX3:2]=[OX1:3]', 'keto → enol', 'The keto form is normally far more stable; the enol is the minor tautomer (but it is the reactive one).'],
      ['[CX4;!H0:1][CX3:2]=[NX2:3]', 'imine → enamine', 'Imine and enamine interconvert; with an N–H the imine is usually favoured.'],
      ['[CX4;!H0:1][NX2:2]=[OX1:3]', 'nitroso → oxime', 'The oxime is much more stable: a C-nitroso compound with an α-H tautomerises to it.'],
      ['[CX4;!H0:1][N+:2](=[OX1:3])[O-]', 'nitro → aci-nitro (nitronic acid)', 'The nitro form is favoured; the aci form explains why 1°/2° nitroalkanes dissolve in alkali.']];
    for (const [sm, what, note] of pats) { const q = M(m, sm).find(([a]) => !A0.has(a)); if (!q) continue; const [ca, y, a] = q; const mech = new Mech(m);
      prot(mech, a, 'Acid catalysis: H⁺ adds to the ' + (el(m, a) === 'O' ? 'oxygen' : 'nitrogen') + '.', ['H⁺ (cat.)']);
      const h = mech.H(ca); const b = baseOn(mech, h, ca, 'O');
      mech.step([[[b], [b, h]], [[h, ca], [ca, y]], [[y, a], [a]]], `Water removes the H⁺ from the α-carbon; its electrons form the C=${el(m, y)} π bond and the old π pair moves onto the heteroatom: the ${what.split('→')[1].trim()}. ${note}`); mech.dropFrag(b);
      return E.done(mech, { title: `Tautomerism: ${what}`, why: ['Tautomers differ only in the position of an H and a π bond (a 1,3-shift)', 'Acid- or base-catalysed; they are in equilibrium'], tags: ['tautomerism'] }); }
    return null;
  } });
/* azo coupling with N,N-dimethylaniline (para) */
front({ id: 'azo_dma', name: 'Azo coupling with N,N-dimethylaniline', fam: 'Amine', canon: 'C₆H₅N(CH₃)₂', db: ['Azo coupling'],
  test: c => c.has('DMA'),
  run(m) {
    const dz = M(m, '[c:1][N+:2]#[N:3]')[0]; if (!dz) return null; const [ca, n1, n2] = dz; const mech = new Mech(m);
    const g = spec('CN(C)c1ccccc1'); const R6 = C.rings(g, 6)[0]; const nA = [...g.atoms.values()].find(a => a.el === 'N').id; const ipso = g.nbrs(nA).find(x => R6.includes(x)); const i = R6.indexOf(ipso); const para = R6[(i + 3) % 6];
    const map = mech.add(g, n2, { att: para, dist: 1.75 * BL }); const pp = map.get(para), nn = map.get(nA), ip = map.get(ipso);
    const r2 = C.rings(mech.m, 6).find(r => r.includes(pp)); const j = r2.indexOf(ip); const c2 = r2[(j + 1) % 6], c3 = r2[(j + 2) % 6];
    if (!kek(mech.m, [[ip, c2], [c3, pp]])) return null;
    mech.step([[[nn], [nn, ip]], [[ip, c2], [c2, c3]], [[c3, pp], [pp, n2]], [[n2, n1], [n1]]], 'The –NMe₂ lone pair pushes electrons through the ring; the para carbon attacks the terminal N of the diazonium ion (a weak electrophile, so it needs this very activated ring). Mildly acidic solution.', ['0–5 °C', 'pH 4–5']);
    settle(mech, 'Loss of the para H⁺ restores the aromatic ring and the N⁺=C becomes N–C again: the azo dye (here methyl orange, as its sulphonate).', t => { dropH(t, pp); setO(t, nn, ip, 1); t.atoms.get(nn).chg = 0; setO(t, ip, c2, 2); setO(t, c2, c3, 1); setO(t, c3, pp, 2); clr(t, [pp, nn, ip, c2, c3]); });
    return E.done(mech, { title: 'Azo coupling with N,N-dimethylaniline → para azo dye', why: ['Para to –NMe₂', 'Coloured: extended conjugation through –N=N–'], tags: ['SEAr'] });
  } });

/* acetal / ketal hydrolysis (deprotection) */
front({ id: 'acetal_hyd', name: 'Acetal / ketal hydrolysis', fam: 'Carbonyl', canon: 'H₃O⁺', db: ['Acetal / ketal formation', 'Acetal protection (Grignard)'],
  test: c => c.any(['H3O+', 'HCl', 'H2SO4']) && (c.has('H2O') || c.has('H3O+') || c.flag('dil') || c.flag('aq')) && !c.any(['ROH', 'glycol']) && !c.flag('glycol') && !c.flag('dry'),
  run(m) {
    const q = M(m, '[CX4:1]([OX2:2][#6:3])[OX2:4][#6:5]').find(([cc, o1, r1, o2]) => m.nbrs(cc).filter(n => el(m, n) === 'O').length === 2 && !E.isCarbonylC(m, r1)); if (!q) return null;
    const [cc, o1, , o2] = q; const mech = new Mech(m);
    prot(mech, o1, 'Aqueous acid protonates one of the ether oxygens.', ['H₃O⁺']);
    mech.step([[[o2], [o2, cc]], [[cc, o1], [o1]]], 'The other O pushes its lone pair in and the protonated O leaves as an alcohol: an oxocarbenium ion.'); 
    const w = put(mech, 'O', cc, isEl('O'), { dist: 1.65 * BL }).at;
    mech.step([[[w], [w, cc]], [[cc, o2], [o2]]], 'Water attacks the oxocarbenium carbon.');
    deprot(mech, w, 'Loss of H⁺: the hemiacetal.');
    prot(mech, o2, 'The second OR oxygen is protonated…');
    mech.step([[[w], [w, cc]], [[cc, o2], [o2]]], '…and the –OH lone pair pushes it out: a protonated carbonyl.');
    deprot(mech, w, 'Loss of H⁺: the aldehyde or ketone is back. Every step is the reverse of acetal formation: excess water drives it this way. (Acetals are stable to base and Grignards, which is why they protect C=O.)');
    const fr = mech.m.frags(); const main = fr.find(f => f.includes(cc)); const R = E.cleanProduct(mech.m, main); const ex = fr.filter(f => f !== main).map(f => E.cleanProduct(mech.m, f)).filter(x => [...x.atoms.values()].some(a => a.el === 'C'));
    return { title: 'Acetal / ketal hydrolysis (H₃O⁺) → carbonyl + alcohol (deprotection)', product: R, extra: ex, steps: mech.steps, why: ['Acetals survive base, Grignards and hydrides; aqueous acid removes them', 'Reverse of acetal formation'], tags: ['hydrolysis'], arrowed: true };
  } });

/* ---------- radicals ---------- */
front({ id: 'rad_couple', name: 'Radical coupling (termination)', fam: 'Radical', canon: 'couple', db: [],
  test: c => c.has('couple'),
  run(m) {
    const r = [...m.atoms.values()].find(a => a.el === 'C' && a.rad === 1); if (!r) return null; const mech = new Mech(m);
    const g = m.clone(); const map = mech.add(g, r.id, { att: r.id, dist: 1.9 * BL }); const r2 = map.get(r.id);
    mech.step([[[r.id], [r.id, r2], 1], [[r2], [r.id, r2], 1]], 'Two radicals meet and each gives its unpaired electron to the new σ bond (two fishhook arrows): termination. Two ethyl radicals give butane.');
    return E.done(mech, { title: 'Radical coupling (chain termination)', why: ['Termination removes radicals: R• + R• → R–R', 'Explains the small amount of ethane in methane chlorination'], tags: ['radical'] });
  } });
front({ id: 'allyl_h_abs', name: 'Allylic H-abstraction by Br•', fam: 'Radical', canon: 'Br•', db: [],
  test: c => c.has('Brrad'),
  run(m) {
    const AL = M(m, '[CX4;!H0:1][CX3:2]=[CX3:3]')[0]; if (!AL) return null; const [s0, a, b] = AL; const mech = new Mech(m);
    const h = mech.H(s0); const x = mech.add({ atoms: [{ el: 'Br', rad: 1, h: 0 }] }, h, { dist: 1.1 * BL, dir: unit([P(mech.m, h)[0] - P(mech.m, s0)[0], P(mech.m, h)[1] - P(mech.m, s0)[1]]) }).get(1);
    mech.step([[[x], [x, h], 1], [[s0, h], [x, h], 1], [[s0, h], [s0], 1]], 'Br• pulls off an allylic H (fishhook arrows: one electron from the radical and one from the C–H bond make H–Br). The allylic C–H is weak because the radical left behind is resonance-stabilised.', ['Br•']); mech.dropFrag(x);
    return E.done(mech, { title: 'Allylic H-abstraction → allyl radical (resonance-stabilised)', why: ['Allylic and benzylic C–H bonds are the weakest: those radicals are delocalised', 'This is the propagation step of NBS allylic bromination'], tags: ['radical'] });
  } });
front({ id: 'f2_explosive', name: 'Alkane + F₂', fam: 'Radical', canon: 'F₂', db: [],
  test: c => c.has('F2'),
  run(m) {
    if (![...m.atoms.values()].every(a => a.el === 'C' || a.el === 'H') || M(m, 'C=C').length || M(m, 'C#C').length || AR(m).size) return null;
    const R = C.fromSmilesMB('[C]', BL); R.atoms.forEach(a => { a.hFix = 0; });
    return E.out({ title: 'Alkane + F₂: explosive, the molecule is destroyed (C + HF)', product: R, extra: [C.fromSmilesMB('F', BL)], steps: [{ L: m.clone(), arrows: [], R: R.clone(), text: 'Fluorination is so exothermic (strong H–F and C–F bonds, weak F–F) that the reaction runs away: C–C and C–H bonds are broken and carbon plus HF result. Controlled fluorination needs dilution with N₂ or reagents such as CoF₃ or SbF₃ (Swarts).', cond: ['F₂'] }], why: ['Reactivity: F₂ ≫ Cl₂ > Br₂ > I₂ (I₂ needs an oxidant such as HIO₃)'], tags: ['radical'] });
  } });
/* ---------- Sandmeyer: SET from Cu(I), aryl radical, ligand transfer from Cu(II) (fishhooks) ---------- */
W('diazo_rx', (m, o, c) => {
  const z = ['CuCl', 'CuBr', 'CuCN'].find(k => c.has(k)); if (!z) return null;
  const dz = M(m, '[c:1][N+:2]#[N:3]')[0]; if (!dz) return null; const [ca, n1, n2] = dz; const Y = { CuCl: 'Cl', CuBr: 'Br', CuCN: 'C#N' }[z];
  const mech = new Mech(m);
  settle(mech, `Cu(I) gives one electron to the diazonium ion (single-electron transfer): the neutral diazenyl radical Ar–N=N•, and Cu(II)${z === 'CuCN' ? '(CN)₂' : Y + '₂'}.`, t => { setO(t, n1, n2, 2); const A = t.atoms.get(n1), B = t.atoms.get(n2); A.chg = 0; B.chg = 0; B.rad = 1; A.hFix = 0; B.hFix = 0; }, [z]);
  mech.step([[[n2], [n1, n2], 1], [[ca, n1], [n1, n2], 1], [[ca, n1], [ca], 1]], 'The C–N bond breaks homolytically: N₂ gas leaves and an aryl radical Ar• is formed.');
  mech.dropFrag(n1);
  const g = spec(Y === 'C#N' ? 'N#C[Cu]C#N' : `${Y}[Cu]${Y}`); const cu = [...g.atoms.values()].find(a => a.el === 'Cu').id; const y0 = g.nbrs(cu)[0];
  const map = mech.add(g, ca, { att: y0, dist: 1.5 * BL }); const Cu = map.get(cu), yA = map.get(y0);
  mech.step([[[ca], [ca, yA], 1], [[Cu, yA], [ca, yA], 1], [[Cu, yA], [Cu], 1]], `The aryl radical takes ${Y === 'C#N' ? 'a CN' : 'a ' + Y} from Cu(II), which is reduced back to Cu(I): Ar–${Y === 'C#N' ? 'CN' : Y}. Cu(I) is regenerated, so it acts as a catalyst.`);
  mech.dropFrag(Cu);
  return fin(mech, o);
});
front({ id: 'diazo_nitro', name: 'Diazonium + NaNO₂ / Cu → nitroarene', fam: 'Amine', canon: 'NaNO₂, Cu, Δ', db: [],
  test: c => c.has('HNO2') && c.has('Cu'),
  run(m) {
    const dz = M(m, '[c:1][N+:2]#[N:3]')[0]; if (!dz) return null; const [ca, n1] = dz; const t = m.clone(); t.removeAtom(dz[2]); t.removeAtom(n1); clr(t, [ca]);
    const p0 = P(t, ca), d = E.freeDir(t, ca); const n = t.addAtom('N', p0[0] + d[0] * BL, p0[1] + d[1] * BL, { chg: 1 }); t.addBond(ca, n, 1); const oa = t.addAtom('O', p0[0] + d[0] * 1.6 * BL, p0[1] + d[1] * 1.6 * BL - BL * .6); t.addBond(n, oa, 2); const ob = t.addAtom('O', p0[0] + d[0] * 1.6 * BL, p0[1] + d[1] * 1.6 * BL + BL * .6, { chg: -1 }); t.addBond(n, ob, 1);
    return E.out({ title: 'Diazonium + NaNO₂ / Cu → nitroarene', product: t, steps: [{ L: m.clone(), arrows: [], R: t.clone(), text: 'The diazonium fluoroborate is warmed with NaNO₂ and copper powder: by the same Cu-mediated radical route as Sandmeyer, –N₂⁺ is replaced by –NO₂. This is a way to put –NO₂ where direct nitration would not.', cond: ['NaNO₂', 'Cu', 'Δ'] }], why: ['Sandmeyer-type replacement'], tags: ['radical'] });
  } });
/* ---------- thiocarboxylate + acyl halide: O-acylation (HSAB) ---------- */
front({ id: 'thioacid_acyl', name: 'Thiocarboxylate + acyl halide', fam: 'Acid derivative', canon: 'CH₃COCl', db: [],
  test: c => c.any(['AcCl', 'RCOCl', 'BzCl']),
  run(m, c) {
    const q = M(m, '[CX3:1](=[O:2])[SX1-:3]')[0]; if (!q) return null; const [c1, o1, s1] = q; const mech = new Mech(m);
    settle(mech, 'Thiocarboxylate is ambident: drawn in its other resonance form, the negative charge sits on oxygen (R–C(=S)–O⁻).', t => { setO(t, c1, o1, 1); setO(t, c1, s1, 2); t.atoms.get(o1).chg = -1; t.atoms.get(s1).chg = 0; clr(t, [o1, s1]); });
    const g = spec(c.has('BzCl') ? 'O=C(Cl)c1ccccc1' : 'CC(=O)Cl'); const aC = [...g.atoms.values()].find(a => carbC(a, g)); const map = mech.add(g, o1, { att: aC.id, dist: 1.8 * BL }); const c2 = map.get(aC.id); const o2 = dblO(mech.m, c2); const x = mech.m.nbrs(c2).find(n => el(mech.m, n) === 'Cl');
    addElim(mech, o1, c2, o2, x, 'The hard acyl carbon is attacked by the hard end of the nucleophile, the oxygen (HSAB: hard–hard).', 'The C=O re-forms and the halide leaves: the O-acyl thioanhydride R–C(=S)–O–COR.', [c.has('BzCl') ? 'C₆H₅COCl' : 'CH₃COX']); mech.dropFrag(x);
    return E.done(mech, { title: 'Thiocarboxylate + acyl halide → O-acylation (hard–hard)', why: ['With CH₃I (a soft sp³ carbon) the same ion reacts at S instead (soft–soft)', 'Ambident nucleophiles: the site depends on the hardness of the electrophile'], tags: ['hsab', 'AdN-E'] });
  } });
/* ---------- homophthalaldehyde + NH₃ → isoquinoline ---------- */
front({ id: 'isoquinoline', name: '1,5-Dicarbonyl + NH₃ → (iso)quinoline / pyridine ring', fam: 'Carbonyl', canon: 'NH₃, H⁺', db: [],
  test: c => c.has('NH3'),
  run(m) {
    const q = M(m, '[CX3H1:1](=[O:2])[CH2:3][c:4]:[c:5][CX3H1:6]=[O:7]')[0]; if (!q) return null; const [c1, o1, ch, , , c2, o2] = q; const mech = new Mech(m);
    const n = put(mech, 'N', c2, isEl('N'), { dist: 1.7 * BL }).at;
    mech.step([[[n], [n, c2]], [[c2, o2], [o2]]], 'NH₃ attacks the aryl aldehyde.', ['NH₃', 'H⁺']);
    settle(mech, 'Proton transfers and loss of water: the imine Ar–CH=NH.', t => { t.removeAtom(o2); dropH(t, n); setO(t, c2, n, 2); clr(t, [c2, n]); t.atoms.get(n).chg = 0; });
    mech.step([[[n], [n, c1]], [[c1, o1], [o1]]], 'The imine N is close to the other C=O and attacks it: the six-membered ring closes.', [], relayFrag(mech, n));
    return jump(mech, { product: (() => { const t = mech.m.clone(); t.removeAtom(o1); dropH(t, n); dropH(t, ch); clr(t, [n, c1, ch]); t.atoms.get(n).chg = 0; setO(t, n, c1, 1); setO(t, c1, ch, 2); return E.cleanProduct(t, E.fragIds(t, c1)); })() }, 'Loss of water, then the CH₂ loses an H: the new ring has six π electrons, so it is aromatic. Isoquinoline.') && E.done(mech, { title: 'Homophthalaldehyde + NH₃ → isoquinoline', why: ['Two C=O groups 1,5 to each other + NH₃ close a six-membered N ring', 'Aromaticity drives the final dehydration'], tags: ['cyclisation'] });
  } });

front({ id: 'schiff', name: "Schiff's test", fam: 'Carbonyl', canon: "Schiff's reagent", db: ["Schiff's test"],
  test: c => c.has('Schiff'),
  run(m) {
    const ald = M(m, '[CX3H1:1](=O)[#6,#1]').concat(M(m, '[CH2:1]=O'));
    const sugar = M(m, '[CX4][OX2H1]').length >= 3;
    if (ald.length && !sugar) return E.out({ title: "Schiff's test positive: pink-magenta colour", product: null, none: 'The aldehyde adds to the decolourised (SO₂-bleached) rosaniline and restores the conjugated dye: the colourless reagent turns pink-magenta. Ketones do not do this.', positive: true, why: ['Aldehydes positive, ketones negative', 'Glucose is negative: in solution it is almost entirely the cyclic hemiacetal, so there is too little free –CHO'], tags: ['test'] });
    if (ald.length && sugar) return E.out({ title: "Schiff's test negative (glucose)", product: null, none: 'Glucose exists mostly as the cyclic hemiacetal in solution, so it does not restore the colour of Schiff\'s reagent, even though its open form is an aldehyde. One of the observations that pointed to the ring structure.', why: [], tags: ['test'] });
    if (M(m, '[#6][CX3](=O)[#6]').length) return E.out({ title: "Schiff's test negative", product: null, none: 'Ketones do not restore the colour of Schiff\'s reagent: only aldehydes do.', why: [], tags: ['test'] });
    return null;
  } });

root.EngineR10 = { fin, spec, put, graftFrom };
})(typeof window !== 'undefined' ? window : globalThis);
