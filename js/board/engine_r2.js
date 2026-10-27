/* engine_r2.js: carbonyls, acids and derivatives, aromatic substitution, amines, diazonium salts, phenols */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore;
const { BL, el, P, heavyN, hTot, cN, aromSet, Mech, unit, freeDir, out, done } = E;
const REC = E.REC, reg = r => REC.push(r);
const M = (m, s) => C.matchAll(m, s);
const dir = (m, a, b) => { const p = P(m, a), q = P(m, b); return unit([q[0] - p[0], q[1] - p[1]]); };
const ws = s => s.any(['H3O+', 'H2O', 'NH4Cl']) && s.st.sp.size <= 2;
function ord(smi) { const g = C.fromSmilesMB(smi, BL); if (!g) return null; for (const a of g.atoms.values()) { a.hFix = null; if (a.el === 'H') a.hFix = null; } return { m: g, first: [...g.atoms.keys()][0] }; }
function graft(t, at, smi, o = 1, d) { const g = ord(smi); const map = E.placeSpecies(t, g.m, at, { att: g.first, dist: BL, dir: d || freeDir(t, at) }); const f = map.get(g.first); t.addBond(at, f, o); return { map, first: f }; }
function cut(t, keep, rootId) { for (const id of E.subtree(t, rootId, keep)) t.removeAtom(id); }
function swap(t, at, old, smi, o = 1) { const d = dir(t, at, old); cut(t, at, old); return graft(t, at, smi, o, d); }
function frags(t) { return t.frags().map(ids => ({ ids, n: ids.filter(id => el(t, id) === 'C').length, h: ids.length })).sort((p, q) => q.n - p.n || q.h - p.h).map(f => E.cleanProduct(t, f.ids)); }
function tx(L, product, o) { return out(Object.assign({ product, steps: [{ L: L.clone(), arrows: [], R: product.clone(), text: o.text, cond: o.cond || [] }] }, o)); }
function txF(L, t, o) { const fr = frags(t); const main = o.pick ? o.pick(fr) : fr[0]; return tx(L, main, Object.assign({ extra: fr.filter(f => f !== main) }, o)); }
const none = (title, why) => out({ title, product: null, none: why, why: [] });
/* perception */
function carbonyls(m) { return M(m, '[CX3:1]=[OX1:2]').filter(([c, o]) => heavyN(m, c).filter(n => n !== o).every(n => el(m, n) === 'C')); }
const isAld = (m, c) => hTot(m, c) >= 1;
function alphaH(m, c) { const AR = aromSet(m); return cN(m, c).filter(a => !AR.has(a) && !E.isCarbonylC(m, a) && hTot(m, a) > 0 && m.bondsOf(a).every(b => b.o === 1)).sort((p, q) => hTot(m, q) - hTot(m, p)); }
const esters = m => M(m, '[CX3:1](=[OX1:2])[OX2:3][#6:4]');
const acids = m => M(m, '[CX3:1](=[OX1:2])[OX2H1:3]');
const carboxylates = m => M(m, '[CX3:1](=[OX1:2])[OX1-:3]');
const acylX = m => M(m, '[CX3:1](=[OX1:2])[Cl,Br:3]');
const anhyd = m => M(m, '[CX3:1](=[OX1:2])[OX2:3][CX3:4]=[OX1:5]');
const amides = m => M(m, '[CX3:1](=[OX1:2])[NX3:3]');
const nitriles = m => M(m, '[CX2:1]#[NX1:2]');
const epox = m => M(m, '[C:1]1[O:2][C:3]1');
const arOH = m => M(m, '[c:1][OX2H1,OX1-:2]');
const arNH2 = m => M(m, '[c:1][NX3H2:2]');
const nitroAr = m => M(m, '[c:1][N+:2](=[O:3])[O-:4]');
const nitroAny = m => M(m, '[#6:1][N+:2](=[O:3])[O-:4]');
const diazo = m => M(m, '[c:1][N+:2]#[N:3]');
const RGS = c => { const g = c.get('RMgX') || c.get('RLi') || c.get('R2CuLi'); return g ? (g.smi || 'C') : null; };
const CARB = { HCHO: 'C=O', CH3CHO: 'CC=O', acetone: 'CC(C)=O', PhCHO: 'O=CC1=CC=CC=C1' };
const partnerKey = c => ['HCHO', 'CH3CHO', 'acetone', 'PhCHO'].find(k => c.has(k));
function protO(mech, o, text) { const ids = mech.add('[H][OH2+]', o, { dist: 1.5 * BL }); const h = ids.get(1), w = ids.get(2); mech.step([[[o], [o, h]], [[h, w], [w]]], text || 'Work-up: the alkoxide takes H⁺ from H₃O⁺.', ['H₃O⁺']); mech.dropFrag(w); }
function waterH(mech, o, text) { const ids = mech.add('[H]O', o, { dist: 1.5 * BL }); const h = ids.get(1), w = ids.get(2); mech.step([[[o], [o, h]], [[h, w], [w]]], text || 'Water protonates the alkoxide (OH⁻ is regenerated).'); mech.dropFrag(w); }
function deprotC(mech, a, base, text) { const h = mech.H(a); const b = mech.add(base || '[OH-]', h, { dist: 1.45 * BL, dir: dir(mech.m, a, h) }).get(1); mech.step([[[b], [b, h]], [[a, h], [a]]], text); mech.dropFrag(b); }
/* ================= carbonyl compounds ================= */
reg({ id: 'hyd_red', name: 'Hydride reduction (NaBH₄ / LiAlH₄)', fam: 'Carbonyl', canon: 'NaBH₄, EtOH', db: ['Meerwein–Ponndorf–Verley reduction'], eats: ws,
  test: c => c.any(['NaBH4', 'LiAlH4']) && !c.has('Hg2+'),
  run(m, c) {
    const lah = c.has('LiAlH4'), cs = carbonyls(m), ac = acylX(m);
    if (cs.length) {
      const [cc, o] = cs[0]; const mech = new Mech(m);
      const ids = mech.add({ atoms: [{ el: lah ? 'Al' : 'B', chg: -1, h: 3 }, { el: 'H', x: -0.82 }], bonds: [[0, 1]] }, cc, { dist: 1.75 * BL, att: 2 });
      const b = ids.get(1), h = ids.get(2);
      mech.step([[[b, h], [h, cc]], [[cc, o], [o]]], `Hydride (H⁻) from ${lah ? 'AlH₄⁻' : 'BH₄⁻'} adds to the carbonyl carbon; the C=O π electrons move onto O.`, [lah ? 'LiAlH₄' : 'NaBH₄']);
      mech.dropFrag(b); protO(mech, o, lah ? 'Work-up: H₃O⁺ protonates the alkoxide.' : 'The solvent (or H₃O⁺) protonates the alkoxide.');
      const kind = isAld(m, cc) ? (hTot(m, cc) >= 2 ? 'HCHO → methanol' : 'aldehyde → 1° alcohol') : 'ketone → 2° alcohol';
      return done(mech, { title: `${lah ? 'LiAlH₄' : 'NaBH₄'}: ${kind}`, why: ['H⁻ adds to the carbonyl carbon, then O picks up H⁺', lah ? 'LiAlH₄ also reduces esters, acids, amides and nitriles' : 'NaBH₄ leaves esters, acids, amides and C=C alone'], tags: ['AdN'] });
    }
    if (!lah) {
      if (ac.length) { const t = m.clone(); t.removeAtom(ac[0][2]); t.bondBetween(ac[0][0], ac[0][1]).o = 1; return tx(m, t, { title: 'NaBH₄: acyl chloride → 1° alcohol', text: 'Acyl chlorides are reactive enough for NaBH₄: reduced to the 1° alcohol.', why: [] }); }
      if (esters(m).length || acids(m).length || amides(m).length) return none('No reaction with NaBH₄', 'NaBH₄ is too mild for esters, acids and amides. LiAlH₄ reduces them.');
      return null;
    }
    const t = m.clone(); const did = [];
    for (const [c1, o2, o3] of esters(t)) { if (!t.bondBetween(c1, o3)) continue; t.removeBond(t.bondBetween(c1, o3).id); t.bondBetween(c1, o2).o = 1; did.push('ester → two alcohols'); }
    for (const [c1, o2, o3] of acids(t)) { if (!t.atoms.has(o3)) continue; t.removeAtom(o3); t.bondBetween(c1, o2).o = 1; did.push('acid → 1° alcohol'); }
    for (const [c1, o2, x] of acylX(t)) { if (!t.atoms.has(x)) continue; t.removeAtom(x); t.bondBetween(c1, o2).o = 1; did.push('acyl chloride → 1° alcohol'); }
    for (const [, o2] of amides(t)) { if (!t.atoms.has(o2)) continue; t.removeAtom(o2); did.push('amide → amine (C=O becomes CH₂)'); }
    for (const [c1, n2] of nitriles(t)) { const bd = t.bondBetween(c1, n2); if (bd) { bd.o = 1; did.push('nitrile → 1° amine'); } }
    for (const [c1, o2, c3] of epox(t)) { if (!t.bondBetween(c1, o2)) continue; const less = cN(t, c1).length <= cN(t, c3).length ? c1 : c3; t.removeBond(t.bondBetween(less, o2).id); did.push('epoxide → alcohol (H⁻ hits the less hindered carbon)'); }
    for (const [c1, n2] of nitroAny(t).filter(x => !aromSet(t).has(x[0]))) { for (const o of heavyN(t, n2).filter(x => el(t, x) === 'O')) t.removeAtom(o); t.atoms.get(n2).chg = 0; did.push('nitroalkane → amine'); }
    if (!did.length) return null;
    return txF(m, t, { title: 'LiAlH₄ reduction', text: 'LiAlH₄ delivers H⁻ to every polar multiple bond it can reach: ' + [...new Set(did)].join('; ') + '. Then H₃O⁺ work-up.', why: [...new Set(did)], tags: ['red'] });
  } });
function rmgx(rs, X) { const g = ord(rs); const t = g.m, a = g.first; const d = freeDir(t, a), p = P(t, a); const mg = t.addAtom('Mg', p[0] + d[0] * BL, p[1] + d[1] * BL); t.addBond(a, mg, 1); const x = t.addAtom(X || 'Br', p[0] + d[0] * 2 * BL, p[1] + d[1] * 2 * BL); t.addBond(mg, x, 1); return { m: t, c: a, mg }; }
reg({ id: 'grignard_add', name: 'Grignard reagent addition', fam: 'Carbonyl', canon: 'CH₃MgBr then H₃O⁺', db: ['Grignard addition to carbonyl', 'Grignard + ester / acyl halide (excess)', 'Grignard + nitrile → ketone', 'Grignard + epoxide', 'Grignard as base (active H)'], eats: ws,
  test: c => c.any(['RMgX', 'RLi']),
  run(m, c) {
    const rs = RGS(c), X = (c.get('RMgX') || {}).X || 'Br';
    if (M(m, '[OX2H1]').length || M(m, '[NX3;!H0]').length || M(m, '[CH1]#C').length) { const p = ord(rs).m; return tx(m, p, { title: 'Grignard acts as a base: R–H forms', text: 'This molecule has an acidic H (O–H, N–H, COOH or ≡C–H). R–MgX is a very strong base, so it just takes that proton: the alkane R–H is the product, not an addition product.', why: ['Active H reacts first (the Zerewitinoff method counts them)', 'Protect –OH groups before using a Grignard'], tags: ['acidbase'] }); }
    const cs = carbonyls(m);
    if (cs.length) {
      const [cc, o] = cs[0]; const g = rmgx(rs, X); const mech = new Mech(m);
      const map = E.placeSpecies(mech.m, g.m, cc, { att: g.c, dist: 1.8 * BL }); const rc = map.get(g.c), mg = map.get(g.mg);
      mech.step([[[rc, mg], [rc, cc]], [[cc, o], [o]]], 'The C–Mg bond is polarised (Cδ⁻–Mgδ⁺): that carbon attacks the carbonyl carbon and the π electrons move onto O.', ['dry ether']);
      mech.dropFrag(mg); protO(mech, o);
      const what = isAld(m, cc) ? (hTot(m, cc) >= 2 ? 'HCHO → 1° alcohol' : 'aldehyde → 2° alcohol') : 'ketone → 3° alcohol';
      return done(mech, { title: 'Grignard + ' + what, why: ['HCHO gives 1°, other aldehydes 2°, ketones 3° alcohols', 'Keep everything dry until the H₃O⁺ work-up'], tags: ['AdN'] });
    }
    const es = esters(m), ax = acylX(m);
    if (es.length || ax.length) {
      const one = c.cold || c.flag('one'); const t = m.clone();
      if (es.length) { const [c1, o2, o3] = es[0]; t.removeBond(t.bondBetween(c1, o3).id); if (one) graft(t, c1, rs); else { t.bondBetween(c1, o2).o = 1; graft(t, c1, rs); graft(t, c1, rs); } }
      else { const [c1, o2, x] = ax[0]; t.removeAtom(x); if (one) graft(t, c1, rs); else { t.bondBetween(c1, o2).o = 1; graft(t, c1, rs); graft(t, c1, rs); } }
      return txF(m, t, { title: one ? 'Grignard (1 eq) + acyl compound → ketone' : 'Grignard (excess) + ester → 3° alcohol', text: one ? 'One R adds and the leaving group (OR′ or Cl) is expelled: a ketone.' : 'The first R adds and pushes out OR′ to give a ketone, which is more reactive than the ester, so a second R adds at once. After H₃O⁺: a 3° alcohol with two identical R groups.', why: ['Two equivalents of R–MgX are used up', 'Esters of HCOOH give 2° alcohols'], tags: ['AdN-E'] });
    }
    const ns = nitriles(m);
    if (ns.length) { const t = m.clone(); const [c1, n2] = ns[0]; t.atoms.get(n2).el = 'O'; t.bondBetween(c1, n2).o = 2; graft(t, c1, rs); return tx(m, t, { title: 'Grignard + nitrile → ketone', text: 'R adds to the C≡N carbon to give an imine salt; H₃O⁺ hydrolyses the C=N to C=O.', why: ['Stops at the ketone: the imine salt does not react again'], tags: ['AdN'] }); }
    const ep = epox(m);
    if (ep.length) { const t = m.clone(); const [c1, o2, c3] = ep[0]; const less = cN(t, c1).length <= cN(t, c3).length ? c1 : c3; t.removeBond(t.bondBetween(less, o2).id); graft(t, less, rs); return tx(m, t, { title: 'Grignard + epoxide → alcohol (2 C longer)', text: 'R attacks the less hindered epoxide carbon (SN2) and opens the ring; H₃O⁺ gives the alcohol.', why: ['Attack at the less substituted carbon'], tags: ['SN2'] }); }
    if (M(m, 'O=C=O').length) { const t = m.clone(); const [cc, o] = M(t, '[C:1]=[O:2]')[0]; t.bondBetween(cc, o).o = 1; graft(t, cc, rs); return tx(m, t, { title: 'Grignard + CO₂ → carboxylic acid', text: 'R adds to CO₂ (carboxylation); H₃O⁺ gives R–COOH, one carbon longer.', why: [] }); }
    return null;
  } });
reg({ id: 'grignard_sub', name: 'A Grignard reagent reacting', fam: 'Carbonyl', canon: 'CO₂ then H₃O⁺', db: ['Grignard + CO2'], eats: ws,
  test: c => c.any(['CO2', 'HCHO', 'CH3CHO', 'acetone', 'PhCHO', 'ethyleneoxide', 'H2O', 'H3O+', 'ROH']),
  run(m, c) {
    const g = M(m, '[#6:1][Mg:2]'); if (!g.length) return null; const [rc, mg] = g[0]; const t = m.clone(); cut(t, rc, mg);
    if (c.has('CO2')) { graft(t, rc, 'C(=O)O'); return tx(m, t, { title: 'Grignard + CO₂ → carboxylic acid', text: 'The nucleophilic carbon adds to CO₂; H₃O⁺ work-up gives the acid with one more carbon.', why: [] }); }
    const k = partnerKey(c);
    if (k) { graft(t, rc, { HCHO: 'CO', CH3CHO: 'C(C)O', acetone: 'C(C)(C)O', PhCHO: 'C(O)C1=CC=CC=C1' }[k]); return tx(m, t, { title: `Grignard + ${k === 'HCHO' ? 'HCHO → 1°' : k === 'acetone' ? 'ketone → 3°' : 'aldehyde → 2°'} alcohol`, text: 'R adds to the carbonyl carbon; H₃O⁺ protonates the alkoxide.', why: [] }); }
    if (c.has('ethyleneoxide')) { graft(t, rc, 'CCO'); return tx(m, t, { title: 'Grignard + ethylene oxide → 1° alcohol (2 C longer)', text: 'R opens the epoxide ring; H₃O⁺ gives R–CH₂CH₂–OH.', why: [] }); }
    return tx(m, t, { title: 'Grignard + water / alcohol → alkane', text: 'Any proton source destroys R–MgX: R–H forms.', why: ['That is why Grignard work must be dry'] });
  } });
reg({ id: 'hcn', name: 'Cyanohydrin formation', fam: 'Carbonyl', canon: 'HCN', db: ['Cyanohydrin formation (HCN)'],
  test: c => c.has('HCN') || (c.has('NaCN') && c.any(['H3O+', 'HCl', 'H2SO4'])),
  run(m) {
    const cs = carbonyls(m); if (!cs.length) return null; const [cc, o] = cs[0]; const mech = new Mech(m);
    const cn = mech.add('[C-]#N', cc, { dist: 1.7 * BL }).get(1);
    mech.step([[[cn], [cn, cc]], [[cc, o], [o]]], 'CN⁻ (from HCN with a trace of base) adds to the carbonyl carbon.');
    const ids = mech.add('[H]C#N', o, { dist: 1.5 * BL }); const h = ids.get(1), c2 = ids.get(2);
    mech.step([[[o], [o, h]], [[h, c2], [c2]]], 'The alkoxide takes H⁺ from HCN, regenerating CN⁻.'); mech.dropFrag(c2);
    return done(mech, { title: 'Cyanohydrin', why: ['Nucleophilic addition; aldehydes react faster than ketones', 'Hydrolysis of the CN gives an α-hydroxy acid'], tags: ['AdN'] });
  } });
const NZ = [['NH2OH', 'NO', 'oxime'], ['PhNHNH2', 'NNC1=CC=CC=C1', 'phenylhydrazone'], ['DNP', 'NNC1=CC=C(C=C1[N+](=O)[O-])[N+](=O)[O-]', '2,4-DNP derivative (orange-red precipitate)'], ['semicarbazide', 'NNC(N)=O', 'semicarbazone'], ['NH2NH2', 'NN', 'hydrazone'], ['RNH2', null, 'imine (Schiff base)'], ['aniline', 'NC1=CC=CC=C1', 'imine (Schiff base)'], ['NH3', 'N', 'imine']];
reg({ id: 'imine', name: 'Addition of NH₂–Z (C=N)', fam: 'Carbonyl', canon: 'NH₂OH, H⁺', db: ['Imine formation (carbonyl + H2N–Z)'],
  test: c => NZ.some(z => c.has(z[0])) && !(c.has('NH2NH2') && (c.any(['KOH', 'aqKOH', 'NaOEt']) || c.heat || c.flag('glycol'))),
  run(m, c) {
    const cs = carbonyls(m); if (!cs.length) return null; const [cc, o] = cs[0]; const z = NZ.find(z => c.has(z[0]));
    const smi = z[1] || ('N' + ((c.get('RNH2') || {}).smi || 'C'));
    const t = m.clone(); const d = dir(t, cc, o); t.removeAtom(o); graft(t, cc, smi, 2, d);
    return tx(m, t, { title: `C=O → C=N–Z: ${z[2]}`, text: 'NH₂–Z adds to C=O (tetrahedral carbinolamine), then water is lost: C=N–Z. Best near pH 4.5, where C=O is activated but the nitrogen is not fully protonated.', why: ['Addition, then elimination of H₂O', 'Used to identify carbonyl compounds'], tags: ['AdN-E'] });
  } });
reg({ id: 'wolff', name: 'Wolff–Kishner / Clemmensen', fam: 'Carbonyl', canon: 'Zn–Hg, conc. HCl', db: ['Wolff–Kishner reduction', 'Clemmensen reduction'],
  test: c => (c.has('NH2NH2') && (c.any(['KOH', 'aqKOH', 'NaOEt']) || c.heat || c.flag('glycol'))) || (c.has('ZnHg')) || (c.has('redP') && c.has('HI')),
  run(m, c) {
    const cs = carbonyls(m); if (!cs.length) return null; const t = m.clone(); for (const [, o] of cs) t.removeAtom(o);
    const wk = c.has('NH2NH2');
    return tx(m, t, { title: wk ? 'Wolff–Kishner: C=O → CH₂' : c.has('ZnHg') ? 'Clemmensen: C=O → CH₂' : 'HI / red P: C=O → CH₂', text: wk ? 'Hydrazone formation, then base and heat drive off N₂: the carbonyl carbon becomes CH₂. Basic conditions, so use it for acid-sensitive molecules.' : 'Zn–Hg with conc. HCl removes the carbonyl oxygen completely: C=O → CH₂. Acidic conditions, so use it for base-sensitive molecules.', why: ['Carbon skeleton unchanged'], tags: ['red'] });
  } });
function aldolMech(m, c, sub, partner, part) {
  const mech = new Mech(m); let nuA, ce, oe;
  const AH = alphaH(m, sub[0]); const pc = partner ? carbonyls(partner)[0] : null; const pAH = pc ? alphaH(partner, pc[0]) : [];
  if (partner && !AH.length && pAH.length) {
    const map = E.placeSpecies(mech.m, partner, sub[0], { att: pAH[0], dist: 2.2 * BL }); nuA = map.get(pAH[0]); ce = sub[0]; oe = sub[1];
  } else if (partner && AH.length) {
    const map = E.placeSpecies(mech.m, partner, AH[0], { att: pc[0], dist: 2.2 * BL }); nuA = AH[0]; ce = map.get(pc[0]); oe = map.get(pc[1]);
  } else if (!partner && AH.length) {
    const map = E.placeSpecies(mech.m, m.clone(), AH[0], { att: sub[0], dist: 2.3 * BL }); nuA = AH[0]; ce = map.get(sub[0]); oe = map.get(sub[1]);
  } else return null;
  deprotC(mech, nuA, '[OH-]', 'OH⁻ removes an α-H: the enolate (shown as its carbanion form) forms.');
  mech.step([[[nuA], [nuA, ce]], [[ce, oe], [oe]]], 'The enolate carbon attacks the carbonyl carbon of the other molecule.');
  waterH(mech, oe, 'Water protonates the alkoxide: a β-hydroxy carbonyl compound (the aldol).');
  const aryl = heavyN(mech.m, ce).some(n => aromSet(mech.m).has(n));
  let cond = false;
  if ((c.heat || aryl) && hTot(mech.m, nuA) > 0) {
    const h = mech.H(nuA); const b = mech.add('[OH-]', h, { dist: 1.45 * BL, dir: dir(mech.m, nuA, h) }).get(1);
    mech.step([[[b], [b, h]], [[nuA, h], [nuA, ce]], [[ce, oe], [oe]]], 'On heating (or when the new C=C would be conjugated with a ring), OH⁻ removes the second α-H and OH⁻ is expelled (E1cB): an α,β-unsaturated carbonyl.', ['Δ']);
    mech.dropFrag(b); mech.dropFrag(oe); cond = true;
  }
  return { mech, cond };
}
reg({ id: 'aldol', name: 'Aldol reaction / condensation', fam: 'Carbonyl', canon: 'dil. NaOH', db: ['Aldol reaction / condensation', 'Claisen–Schmidt condensation', 'Crossed aldol with HCHO'],
  test: c => c.any(['KOH', 'aqKOH', 'NaOEt', 'NaOMe']) && !c.any(['alcKOH', 'CHCl3', 'Br2/NaOH', 'NaOCl', 'CO2', 'I2', 'Br2', 'Cl2']),
  run(m, c) {
    const cs = carbonyls(m); if (!cs.length) return null; const k = partnerKey(c);
    const partner = k ? C.fromSmiles(CARB[k], BL) : null;
    const r = aldolMech(m, c, cs[0], partner);
    if (!r) { if (!partner && c.flag('dil') && cs.every(([cc]) => !alphaH(m, cc).length)) return none('No aldol', 'No α-hydrogen, so no enolate can form. (With conc. NaOH an aldehyde like this does Cannizzaro instead.)'); return null; }
    const crossed = !!partner;
    return done(r.mech, { title: (crossed ? (k === 'PhCHO' || aromSet(m).size ? 'Claisen–Schmidt / crossed ' : 'Crossed ') : '') + (r.cond ? 'aldol condensation → α,β-unsaturated carbonyl' : 'aldol → β-hydroxy carbonyl'), why: ['Needs an α-H on the molecule that forms the enolate', crossed ? 'Crossed aldol is clean when one partner has no α-H (HCHO, PhCHO)' : 'Heat dehydrates the aldol (E1cB)'], tags: ['enolate'] });
  } });
reg({ id: 'cannizzaro', name: 'Cannizzaro reaction', fam: 'Carbonyl', canon: 'conc. NaOH', db: ['Cannizzaro reaction', 'Crossed Cannizzaro'],
  test: c => c.any(['KOH', 'aqKOH']) && !c.flag('dil') && !c.any(['alcKOH', 'CHCl3', 'Br2/NaOH', 'NaOCl', 'CO2', 'I2', 'Br2']),
  run(m, c) {
    const cs = carbonyls(m).filter(([cc]) => isAld(m, cc) && !alphaH(m, cc).length); if (!cs.length) return null;
    const [cc, o] = cs[0]; const crossed = c.has('HCHO') && hTot(m, cc) < 2; const mech = new Mech(m);
    let other;
    if (crossed) other = E.placeSpecies(mech.m, C.fromSmiles('C=O', BL), cc, { att: 1, dist: 2.4 * BL });
    else other = E.placeSpecies(mech.m, m.clone(), cc, { att: cc, dist: 2.6 * BL });
    // attacked molecule A (oxidised) and hydride acceptor B (reduced)
    const [aC, aO, bC, bO] = crossed ? (() => { const c1 = [...other.values()].find(id => el(mech.m, id) === 'C'); const o1 = [...other.values()].find(id => el(mech.m, id) === 'O'); return [c1, o1, cc, o]; })() : [cc, o, other.get(cc), other.get(o)];
    const oh = mech.add('[OH-]', aC, { dist: 1.6 * BL }).get(1);
    mech.step([[[oh], [oh, aC]], [[aC, aO], [aO]]], `OH⁻ adds to the carbonyl carbon of ${crossed ? 'HCHO (the more reactive aldehyde)' : 'one aldehyde molecule'}.`);
    const h = mech.H(aC);
    mech.step([[[aO], [aO, aC]], [[aC, h], [h, bC]], [[bC, bO], [bO]]], 'The O⁻ pushes back down; the C–H leaves as a hydride and adds to the second aldehyde (the slow step).');
    const ho = E.hExp(mech.m, oh)[0] || mech.H(oh);
    mech.step([[[bO], [bO, ho]], [[ho, oh], [oh]]], 'Proton transfer: the acid gives H⁺ to the alkoxide. Result: a carboxylate ion and an alcohol.');
    const alc = E.cleanProduct(mech.m, E.fragIds(mech.m, bC)), carbx = E.cleanProduct(mech.m, E.fragIds(mech.m, aC));
    return out({ title: crossed ? 'Crossed Cannizzaro: HCHO is oxidised, the other aldehyde reduced' : 'Cannizzaro: alcohol + carboxylate', product: alc, extra: [carbx], steps: mech.steps, why: ['Only aldehydes with no α-H', 'One molecule is oxidised, the other reduced (disproportionation)'], tags: ['hydride_transfer'] });
  } });
reg({ id: 'tollens', name: "Tollens' / Fehling's test", fam: 'Carbonyl', canon: "Tollens' reagent", db: ["Tollens' test", "Fehling's / Benedict's test"],
  test: c => c.any(['Tollens', 'Fehling']),
  run(m, c) {
    const cs = carbonyls(m); if (!cs.length) return null; const ald = cs.filter(([cc]) => isAld(m, cc)); const fe = c.has('Fehling');
    if (!ald.length) return none(fe ? "Negative Fehling's test" : "Negative Tollens' test", 'Ketones are not oxidised by these mild reagents (α-hydroxy ketones are the exception).');
    const [cc] = ald[0];
    if (fe && heavyN(m, cc).some(n => aromSet(m).has(n))) return none("Negative Fehling's test", 'Aromatic aldehydes (like benzaldehyde) do not reduce Fehling\'s solution, though they do give the silver mirror with Tollens\'.');
    const t = m.clone(); t.addBond(cc, t.addAtom('O', P(t, cc)[0] + freeDir(t, cc)[0] * BL, P(t, cc)[1] + freeDir(t, cc)[1] * BL, { chg: -1 }), 1);
    return tx(m, t, { title: fe ? "Fehling's: red Cu₂O precipitate, aldehyde → carboxylate" : "Tollens': silver mirror, aldehyde → carboxylate", text: fe ? 'Cu²⁺ (blue) is reduced to Cu₂O (red-brown precipitate); the aldehyde is oxidised to its carboxylate.' : 'Ag(NH₃)₂⁺ is reduced to metallic silver (the mirror); the aldehyde is oxidised to its carboxylate.', why: ['Aldehydes are oxidised, ketones are not'], tags: ['test'] });
  } });
reg({ id: 'haloform', name: 'Haloform reaction', fam: 'Carbonyl', canon: 'I₂, NaOH', db: ['Haloform reaction'],
  test: c => c.has('NaOCl') || (c.has('I2') && c.any(['KOH', 'aqKOH', 'NaHCO3'])),
  run(m, c) {
    const X = ((c.get('NaOCl') || {}).X) || (c.has('I2') ? 'I' : /br/i.test(c.raw) ? 'Br' : 'Cl');
    let mk = M(m, '[CH3:1][CX3:2]=[OX1:3]'); const t = m.clone();
    if (!mk.length) { const sa = M(t, '[CH3:1][CX4H1:2][OX2H1:3]'); if (!sa.length) return none('Negative haloform test', 'Needs CH₃CO– or CH₃CH(OH)– in the molecule.'); t.bondBetween(sa[0][1], sa[0][2]).o = 2; mk = [sa[0]]; }
    const [me, cc] = mk[0]; const d = dir(t, cc, me); t.removeAtom(me); const op = P(t, cc); t.addBond(cc, t.addAtom('O', op[0] + d[0] * BL, op[1] + d[1] * BL, { chg: -1 }), 1);
    const hf = C.fromSmiles(`${X}C(${X})${X}`, BL);
    return out({ title: `Haloform: carboxylate (one C fewer) + CH${X === 'I' ? 'I₃ (yellow precipitate)' : X + '₃'}`, product: t, extra: [hf], steps: [{ L: m.clone(), arrows: [], R: t.clone(), text: `The CH₃ next to C=O is halogenated three times (each H gets more acidic), then OH⁻ attacks the C=O and ⁻C${X}₃ leaves. Products: R–COO⁻ and CH${X}₃.`, cond: [] }], why: ['Test for CH₃CO– and CH₃CH(OH)– groups', 'Ethanol and ethanal also give iodoform'], tags: ['test'] });
  } });
reg({ id: 'wittig', name: 'Wittig reaction', fam: 'Carbonyl', canon: 'Ph₃P=CH₂', db: ['Wittig reaction'],
  test: c => c.any(['Ph3P=CH2', 'Ph3P=CHCH3']),
  run(m, c) { const cs = carbonyls(m); if (!cs.length) return null; const [cc, o] = cs[0]; const t = m.clone(); const d = dir(t, cc, o); t.removeAtom(o); graft(t, cc, c.has('Ph3P=CHCH3') ? 'CC' : 'C', 2, d); return tx(m, t, { title: 'Wittig: C=O → C=C', text: 'The ylide carbon adds to C=O; a four-membered oxaphosphetane forms and splits into the alkene and Ph₃P=O (the strong P=O bond drives it).', why: ['The C=C appears exactly where the C=O was'], tags: ['cyclicTS'] }); } });
/* ================= acids and derivatives ================= */
reg({ id: 'ester_hyd', name: 'Ester hydrolysis', fam: 'Acid derivative', canon: 'aq. NaOH, Δ', db: ['Ester hydrolysis (BAC2)', 'Saponification'],
  test: c => c.any(['KOH', 'aqKOH', 'H3O+', 'H2SO4', 'HCl', 'H2O']) && !c.any(['alcKOH', 'LiAlH4', 'NaBH4', 'RMgX']),
  run(m, c) {
    const es = esters(m); if (!es.length) return null; const [c1, o2, o3, c4] = es[0];
    if (c.any(['KOH', 'aqKOH'])) {
      const mech = new Mech(m); const oh = mech.add('[OH-]', c1, { dist: 1.65 * BL }).get(1);
      mech.step([[[oh], [oh, c1]], [[c1, o2], [o2]]], 'OH⁻ adds to the ester C=O: a tetrahedral intermediate.');
      mech.step([[[o2], [o2, c1]], [[c1, o3], [o3]]], 'The C=O re-forms and the alkoxide (R′O⁻) is expelled.');
      const h = E.hExp(mech.m, oh)[0] || mech.H(oh);
      mech.step([[[o3], [o3, h]], [[h, oh], [oh]]], 'R′O⁻ takes the acidic H of the acid: carboxylate + alcohol. This last step is irreversible, so saponification goes to completion.');
      const acyl = E.cleanProduct(mech.m, E.fragIds(mech.m, c1)), alc = E.cleanProduct(mech.m, E.fragIds(mech.m, o3));
      return out({ title: 'Saponification (BAC2): carboxylate + alcohol', product: acyl, extra: E.fragIds(mech.m, o3).includes(c1) ? [] : [alc], steps: mech.steps, why: ['Acyl–oxygen bond breaks (B: base, AC: acyl cleavage, 2: bimolecular)', 'Irreversible under base'], tags: ['BAC2'] });
    }
    const t = m.clone(); t.removeBond(t.bondBetween(c1, o3).id); graftOH(t, c1);
    return txF(m, t, { pick: fr => fr.find(f => M(f, 'C(=O)[OH]').length) || fr[0], title: 'Acid hydrolysis of ester: acid + alcohol', text: 'H⁺ activates the C=O, water adds, and the alcohol leaves (AAC2). Reversible: excess water pushes it forward.', why: ['Reverse of Fischer esterification'], tags: ['AAC2'] });
  } });
function graftOH(t, c) { const p = P(t, c), d = freeDir(t, c); const o = t.addAtom('O', p[0] + d[0] * BL, p[1] + d[1] * BL); t.addBond(c, o, 1); return o; }
const NUCS = c => c.has('NH3') ? ['N', 'amide'] : c.has('aniline') ? ['NC1=CC=CC=C1', 'anilide'] : c.has('RNH2') ? ['N' + (c.get('RNH2').smi || 'C'), 'N-alkyl amide'] : c.has('R2NH') ? ['N(' + (c.get('R2NH').smi || 'C') + ')' + (c.get('R2NH').smi || 'C'), 'N,N-dialkyl amide'] : c.has('phenol') ? ['OC1=CC=CC=C1', 'phenyl ester'] : c.has('ROH') ? ['O' + (c.get('ROH').smi || 'C'), 'ester'] : c.has('NaOAc') ? ['OC(C)=O', 'anhydride'] : c.has('H2O') ? ['O', 'carboxylic acid'] : null;
reg({ id: 'acyl_sub', name: 'Acyl substitution of acid chlorides / anhydrides', fam: 'Acid derivative', canon: 'NH₃', db: ['Nucleophilic acyl substitution'],
  test: c => !!NUCS(c),
  run(m, c) {
    const ax = acylX(m), an = anhyd(m); const nu = NUCS(c);
    if (!ax.length && !an.length) {
      const es = esters(m); if (!es.length || !/^N/.test(nu[0])) return null; const [c1, , o3] = es[0]; const t2 = m.clone(); swap(t2, c1, o3, nu[0]);
      return txF(m, t2, { title: `Ammonolysis of ester → ${nu[1]}`, text: 'NH₃ (or the amine) adds to the ester C=O and the alkoxide leaves: amide + alcohol. Slower than with acyl chlorides.', why: ['Esters react with NH₃ / amines to give amides'], tags: ['AdN-E'] });
    }
    const t = m.clone();
    if (ax.length) { const [c1, , x] = ax[0]; swap(t, c1, x, nu[0]); } else { const [c1, , o3] = an[0]; const d = dir(t, c1, o3); t.removeBond(t.bondBetween(c1, o3).id); graft(t, c1, nu[0], 1, d); }
    return txF(m, t, { title: `Acyl substitution → ${nu[1]}`, text: 'The nucleophile adds to C=O (tetrahedral intermediate) and the leaving group (Cl⁻ or RCOO⁻) is expelled. Acyl chloride > anhydride > ester > amide in reactivity.', why: ['Addition–elimination at the acyl carbon'], tags: ['AdN-E'] });
  } });
reg({ id: 'fischer', name: 'Esterification', fam: 'Acid derivative', canon: 'CH₃OH, conc. H₂SO₄', db: ['Fischer esterification'],
  test: c => c.has('ROH') || (c.has('AcOH') && c.any(['H2SO4', 'H3O+'])),
  run(m, c) {
    const ac = acids(m);
    if (ac.length && c.has('ROH')) { const t = m.clone(); const [, , o3] = ac[0]; graft(t, o3, c.get('ROH').smi || 'C'); return tx(m, t, { title: 'Fischer esterification', text: 'Acid-catalysed: the alcohol adds to the protonated C=O and water leaves. The O of the ester comes from the alcohol.', why: ['Reversible: remove water or use excess alcohol'], tags: ['AAC2'] }); }
    const al = M(m, '[CX4:1][OX2H1:2]'); if (al.length && c.has('AcOH')) { const t = m.clone(); graft(t, al[0][1], 'C(C)=O'); return tx(m, t, { title: 'Esterification of the alcohol', text: 'The alcohol O attacks the protonated acid; water is lost.', why: [] }); }
    return null;
  } });
reg({ id: 'acylation', name: 'Acylation (–OH / –NH)', fam: 'Acid derivative', canon: '(CH₃CO)₂O', db: ['Acetylation of amines', 'Acetylation'],
  test: c => c.any(['AcCl', 'Ac2O', 'BzCl', 'RCOCl', 'Hinsberg']),
  run(m, c) {
    const nh = M(m, '[NX3;!H0;!$(NC=O):1]'), oh = M(m, '[OX2H1;!$(OC=O):1]');
    const hins = c.has('Hinsberg');
    if (c.has('FeX3') && !hins && M(m, '[c][NX3;!$(NC=O)]').length) return none('No Friedel–Crafts on aniline', 'The –NH₂ lone pair binds AlCl₃ (a Lewis acid–base salt). N becomes positive and strongly deactivates the ring, so Friedel–Crafts fails. (Without AlCl₃, CH₃COCl would simply acetylate the N.)');
    if (hins) {
      const am = M(m, '[NX3;!$(NC=O);!$(N=*);!$(N#*):1]').filter(([n]) => !m.bondsOf(n).some(b => el(m, m.other(b, n)) === 'O'));
      if (!am.length) return null; const [n] = am[0]; const nh2 = hTot(m, n);
      if (!nh2) return none("No reaction with Hinsberg's reagent", '3° amines have no N–H, so they do not form a sulphonamide (they stay insoluble in KOH, and dissolve in acid).');
      const t = m.clone(); graft(t, n, 'S(=O)(=O)C1=CC=CC=C1');
      return tx(m, t, { title: nh2 >= 2 ? 'Hinsberg: 1° amine → N-alkyl sulphonamide (soluble in KOH)' : 'Hinsberg: 2° amine → N,N-dialkyl sulphonamide (insoluble in KOH)', text: 'The amine N attacks the S of benzenesulphonyl chloride; HCl is lost.', why: ['1°: the sulphonamide still has an acidic N–H and dissolves in KOH', '2°: no N–H left, insoluble', '3°: no reaction'], tags: ['test'] });
    }
    const sm = c.has('BzCl') ? 'C(=O)C1=CC=CC=C1' : c.has('RCOCl') ? 'C(=O)' + (c.get('RCOCl').smi || 'C') : 'C(C)=O';
    const site = nh.length ? nh[0][0] : oh.length ? oh[0][0] : null; if (site == null) return null;
    const t = m.clone(); graft(t, site, sm);
    return tx(m, t, { title: nh.length ? 'N-Acylation → amide' : 'O-Acylation → ester', text: 'The N or O lone pair attacks the acyl carbon; Cl⁻ (or acetate) leaves. Pyridine mops up the HCl.', why: nh.length ? ['Acetylation tames –NH₂ (used to control EAS of aniline)'] : ['Phenols and alcohols both acylate'], tags: ['AdN-E'] });
  } });
reg({ id: 'amide_rx', name: 'Amide reactions', fam: 'Acid derivative', canon: 'Br₂, NaOH', db: ['Hofmann bromamide degradation'],
  test: c => c.any(['Br2/NaOH', 'P2O5', 'SOCl2', 'H3O+', 'KOH', 'aqKOH', 'H2SO4']),
  run(m, c) {
    const am = M(m, '[#6:4][CX3:1](=[OX1:2])[NX3H2:3]'); const anyAm = amides(m); if (!anyAm.length) return null;
    if (c.has('Br2/NaOH') && am.length) { const [r4, c1, o2, n3] = am[0]; const t = m.clone(); const pc = P(t, c1); t.removeAtom(o2); t.removeAtom(c1); const na = t.atoms.get(n3); na.x = pc[0]; na.y = pc[1]; t.addBond(r4, n3, 1); return tx(m, t, { title: 'Hofmann bromamide: amide → 1° amine (one C fewer)', text: 'N-bromination, loss of HBr gives a nitrene-like species, R migrates from C to N (isocyanate), and hydrolysis releases CO₂: R–NH₂.', why: ['Only 1° amides (–CONH₂)', 'The amine has one carbon fewer'], tags: ['1,2-shift'] }); }
    if (c.any(['P2O5', 'SOCl2']) && am.length) { const [, c1, o2, n3] = am[0]; const t = m.clone(); t.removeAtom(o2); t.bondBetween(c1, n3).o = 3; return tx(m, t, { title: 'Dehydration: amide → nitrile', text: 'P₂O₅ (or SOCl₂) removes H₂O from –CONH₂: –C≡N.', why: [] }); }
    if (c.any(['H3O+', 'H2SO4', 'KOH', 'aqKOH'])) { const [c1, o2, n3] = anyAm[0]; const t = m.clone(); const d = dir(t, c1, n3); t.removeBond(t.bondBetween(c1, n3).id); const p = P(t, c1); const o = t.addAtom('O', p[0] + d[0] * BL, p[1] + d[1] * BL, { chg: c.any(['KOH', 'aqKOH']) ? -1 : 0 }); t.addBond(c1, o, 1); return txF(m, t, { title: 'Amide hydrolysis', text: 'Amides need heat with acid or base; the C–N bond breaks to give the acid (or carboxylate) and NH₃ / amine.', why: ['Slowest of the acid derivatives'], tags: ['hydrolysis'] }); }
    return null;
  } });
reg({ id: 'nitrile_rx', name: 'Nitrile reactions', fam: 'Acid derivative', canon: 'H₃O⁺, Δ', db: ['Stephen reduction'],
  test: c => c.any(['H3O+', 'H2SO4', 'HCl', 'DIBAL', 'SnCl2', 'H2', 'Na', 'KOH', 'aqKOH']),
  run(m, c) {
    const ns = nitriles(m); if (!ns.length) return null; const [c1, n2] = ns[0]; const t = m.clone();
    if (c.has('DIBAL') || c.has('SnCl2')) { t.atoms.get(n2).el = 'O'; t.bondBetween(c1, n2).o = 2; return tx(m, t, { title: c.has('SnCl2') ? 'Stephen reduction: nitrile → aldehyde' : 'DIBAL-H: nitrile → aldehyde', text: 'Partial reduction to the imine, then hydrolysis gives the aldehyde.', why: [] }); }
    if (c.has('H2') || c.has('Na')) { t.bondBetween(c1, n2).o = 1; return tx(m, t, { title: 'Nitrile → 1° amine', text: 'H₂/Ni (or Na/ethanol, the Mendius reduction) adds 2 H₂ across C≡N: R–CH₂NH₂.', why: [] }); }
    const base = c.any(['KOH', 'aqKOH']); t.atoms.get(n2).el = 'O'; t.bondBetween(c1, n2).o = 2; const p = P(t, c1), d = freeDir(t, c1); t.addBond(c1, t.addAtom('O', p[0] + d[0] * BL, p[1] + d[1] * BL, { chg: base ? -1 : 0 }), 1);
    return tx(m, t, { title: base ? 'Nitrile hydrolysis (base) → carboxylate' : 'Nitrile hydrolysis → carboxylic acid', text: 'Water adds across C≡N (via the amide); full hydrolysis gives the acid and NH₄⁺ (or NH₃ in base).', why: ['Alkyl halide → nitrile → acid adds one carbon'], tags: ['hydrolysis'] });
  } });
reg({ id: 'decarbox', name: 'Decarboxylation / Kolbe electrolysis', fam: 'Acid derivative', canon: 'soda lime, Δ', db: ['Soda-lime decarboxylation', 'Kolbe electrolysis'],
  test: c => c.any(['sodalime', 'electrolysis']),
  run(m, c) {
    const ac = acids(m).concat(carboxylates(m)); if (!ac.length) return null; const [c1] = ac[0]; const rC = heavyN(m, c1).find(n => el(m, n) === 'C');
    if (rC == null) return null; const t = m.clone(); cut(t, rC, c1);
    if (c.has('electrolysis')) { const half = t.clone(); const idm = t.merge(half, 0, 0); const r2 = idm.get(rC); const all = [...idm.values()]; const p = P(t, rC), d = freeDir(t, rC); E.rotSet(t, all, P(t, r2)[0], P(t, r2)[1], Math.PI); const q = P(t, r2); E.moveSet(t, all, p[0] + d[0] * BL - q[0], p[1] + d[1] * BL - q[1]); t.addBond(rC, r2, 1); E.tidyLongBonds(t); return tx(m, t, { title: 'Kolbe electrolysis: 2 RCOO⁻ → R–R', text: 'At the anode RCOO⁻ loses an electron and CO₂; two R• radicals combine.', why: ['Even-carbon symmetric alkanes'], tags: ['radical'] }); }
    return tx(m, t, { title: 'Soda-lime decarboxylation → alkane (one C fewer)', text: 'Heating the sodium salt with NaOH/CaO removes CO₂ as Na₂CO₃: R–COONa → R–H.', why: ['One carbon fewer'], tags: ['cleavage'] });
  } });
reg({ id: 'hvz', name: 'Hell–Volhard–Zelinsky', fam: 'Acid derivative', canon: 'Br₂, red P', db: ['HVZ reaction'],
  test: c => c.has('HVZ') || (c.any(['Br2', 'Cl2']) && c.has('redP')),
  run(m, c) { const ac = acids(m); if (!ac.length) return null; const [c1] = ac[0]; const a = alphaH(m, c1)[0]; if (a == null) return none('No HVZ reaction', 'There is no α-H next to the –COOH.'); const t = m.clone(); const X = c.has('Cl2') ? 'Cl' : 'Br'; const p = P(t, a), d = freeDir(t, a); t.addBond(a, t.addAtom(X, p[0] + d[0] * BL, p[1] + d[1] * BL), 1); return tx(m, t, { title: `HVZ: α-${X === 'Br' ? 'bromo' : 'chloro'} acid`, text: 'Red P makes a little acyl halide, which enolises; the enol takes X₂ at the α-carbon.', why: ['Needs an α-H'], tags: ['enol'] }); } });
/* ================= aromatic ================= */
function subsOn(m, ring) { const out = []; for (const a of ring) for (const s of heavyN(m, a)) if (!ring.includes(s)) out.push({ a, s }); return out; }
function gtype(m, a, s) {
  const e = el(m, s), at = m.atoms.get(s);
  if (e === 'O') { if (at.chg < 0) return { op: true, act: 5.5, n: 'O⁻' }; return heavyN(m, s).some(n => n !== a && E.isCarbonylC(m, n)) ? { op: true, act: 1.5, n: '–OCOR' } : { op: true, act: 4, n: hTot(m, s) ? '–OH' : '–OR' }; }
  if (e === 'N') { if (at.chg > 0) return { op: false, act: -4, n: m.bondsOf(s).some(b => el(m, m.other(b, s)) === 'O') ? '–NO₂' : (m.bondsOf(s).some(b => b.o === 3) ? '–N₂⁺' : '–NR₃⁺') }; if (heavyN(m, s).some(n => n !== a && E.isCarbonylC(m, n))) return { op: true, act: 2, n: '–NHCOR' }; return { op: true, act: 4.5, n: '–NH₂ / –NR₂' }; }
  if (E.HAL.has(e)) return { op: true, act: -1, n: '–' + e };
  if (e === 'C') { if (E.isCarbonylC(m, s)) return { op: false, act: -2.6, n: '–CHO / –COR / –COOH' }; if (m.bondsOf(s).some(b => b.o === 3)) return { op: false, act: -3, n: '–CN' }; if (heavyN(m, s).filter(n => E.HAL.has(el(m, n))).length >= 3) return { op: false, act: -3, n: '–CX₃' }; return { op: true, act: 1, n: 'alkyl' }; }
  if (e === 'S') return m.bondsOf(s).some(b => b.o === 2) ? { op: false, act: -3, n: '–SO₃H' } : { op: true, act: 2, n: '–SR' };
  return { op: true, act: 0, n: e };
}
function sites(m, ring, all) {
  const subs = subsOn(m, ring).map(x => Object.assign(x, gtype(m, x.a, x.s)));
  const free = ring.filter(a => hTot(m, a) > 0); if (!free.length) return { subs, list: [] };
  const rel = (p, q) => { const d = Math.abs(ring.indexOf(p) - ring.indexOf(q)); return Math.min(d, 6 - d); };
  if (!subs.length) return { subs, list: [free[0]] };
  const sc = p => { let s = 0; for (const x of subs) { const r = rel(p, x.a); const good = x.op ? (r === 1 || r === 3) : r === 2; s += (good ? 1 : -1) * (Math.abs(x.act) + .5) + (r === 1 ? -0.35 * heavyN(m, x.s).length : 0) + (r === 3 ? .45 : 0); } return s; };
  const ranked = free.map(p => [p, sc(p)]).sort((a, b) => b[1] - a[1]);
  if (all) { const best = Math.max(...subs.map(x => x.act)); const lead = subs.find(x => x.act === best); return { subs, list: free.filter(p => [1, 3].includes(rel(p, lead.a))) }; }
  return { subs, list: [ranked[0][0]], ranked, rel };
}
function easStep(mech, p, eSpec, eAtt, text, cond, extraArrows) {
  const m = mech.m; const q = m.nbrs(p).find(n => m.bondBetween(p, n).o === 2 && aromSet(m).has(n));
  const ids = mech.add(eSpec, p, { dist: 1.7 * BL, att: eAtt || 1 }); const E1 = ids.get(eAtt || 1);
  mech.step([[[p, q], [p, E1]], ...(extraArrows ? extraArrows(ids) : [])], text, cond);
  return { q, E1 };
}
function easFinish(mech, p, q, base, text) { const h = mech.H(p); const b = mech.add(base, h, { dist: 1.45 * BL, dir: dir(mech.m, p, h) }).get(1); mech.step([[[b], [b, h]], [[p, h], [p, q]]], text || 'A base removes H⁺ from the sp³ carbon; the ring becomes aromatic again.'); mech.dropFrag(b); }
reg({ id: 'eas', name: 'Electrophilic aromatic substitution', fam: 'Aromatic', canon: 'HNO₃, H₂SO₄', db: ['Nitration of benzene', 'Halogenation of benzene', 'Sulphonation of benzene', 'Friedel–Crafts alkylation', 'Friedel–Crafts acylation'],
  test: c => c.any(['NO2+', 'HNO3', 'oleum', 'FeX3', 'HOBr', 'HOCl', 'Br2', 'Cl2']) || (c.has('H2SO4') && (c.flag('conc') || c.heat)),
  run(m, c) {
    const rings = C.aromaticRings(m); if (!rings.length) return null; const ring = rings[0];
    const S0 = sites(m, ring); const act = S0.subs.length ? Math.max(...S0.subs.map(x => x.act)) : 0;
    const anilineN = S0.subs.some(x => x.n === '–NH₂ / –NR₂'); const deact = S0.subs.length && S0.subs.every(x => x.act <= -2);
    const why = S0.subs.length ? S0.subs.map(x => `${x.n}: ${x.op ? 'ortho/para' : 'meta'} director, ${x.act > 0 ? 'activating' : 'deactivating'}`) : ['Benzene: any position'];
    const X = c.any(['Br2', 'HOBr']) ? 'Br' : 'Cl';
    const halo = c.any(['Br2', 'Cl2', 'HOBr', 'HOCl']);
    if (halo && !c.has('FeX3')) {
      const modAct = S0.subs.some(x => x.n === '–NHCOR') && (c.any(['AcOH']) || /acoh|ch3cooh|acetic/i.test(c.raw || ''));
      if (act < 3.5 && !modAct) return none('No reaction', `Plain ${X}₂ does not substitute this ring: a Lewis acid (FeBr₃/AlCl₃) is needed${alkeneHint(m)}.`);
      const water = c.any(['HOBr', 'HOCl', 'H2O']) || c.flag('aq');
      const t = m.clone(); const S = sites(t, ring, water);
      const list = water ? S.list : (S.list.length ? [S.list[0]] : []); if (!list.length) return null;
      for (const p of list) { const pp = P(t, p), d = freeDir(t, p); t.addBond(p, t.addAtom(X, pp[0] + d[0] * BL, pp[1] + d[1] * BL), 1); }
      return tx(m, t, { title: water ? `${X}₂ water: every free ortho/para position (${list.length}× ${X})` : `${X}₂ in CS₂ / low T: mono-${X === 'Br' ? 'bromination' : 'chlorination'} (para major)`, text: water ? `–OH / –NH₂ make the ring so electron-rich that ${X}₂ water substitutes all free ortho and para positions at once (e.g. 2,4,6-tribromophenol, white precipitate).` : `In a non-polar solvent at low temperature the reaction can be stopped at one ${X}; para is major, ortho minor.`, why, tags: ['SEAr'] });
    }
    const nitr = c.any(['NO2+', 'HNO3']); const fc = c.has('FeX3') && c.any(['RX', 'AcCl', 'BzCl', 'Ac2O', 'RCOCl']);
    if (fc && (deact || anilineN)) return none('No reaction', anilineN ? 'The –NH₂ lone pair binds the Lewis acid (AlCl₃): the ring becomes strongly deactivated, so Friedel–Crafts / Lewis-acid halogenation fails.' : 'Friedel–Crafts fails on strongly deactivated rings such as nitrobenzene.');
    if (!S0.list.length) return null; const p = S0.list[0];
    const mech = new Mech(m); let title = '', base = 'O';
    if (nitr) {
      if (act >= 3.5 && !c.has('H2SO4') && !c.flag('conc')) { title = 'Dilute HNO₃: o- and p-nitro product'; }
      if (act >= 3.5 && (c.flag('conc') || c.has('H2SO4')) && S0.subs.some(x => x.n === '–OH')) { const t = m.clone(); for (const pp of sites(t, ring, true).list) graft(t, pp, '[N+](=O)[O-]'); return tx(m, t, { title: 'Conc. HNO₃ on phenol → picric acid (2,4,6-trinitrophenol)', text: 'The strongly activated ring is nitrated at every ortho and para position.', why, tags: ['SEAr'] }); }
      const r = easStep(mech, p, '[N+](=O)=O', 1, 'The ring π electrons attack the nitronium ion NO₂⁺ (made from HNO₃ + H₂SO₄): an arenium (σ-complex) ion.', ['HNO₃/H₂SO₄'], ids => [[[ids.get(1), ids.get(2)], [ids.get(2)]]]);
      easFinish(mech, p, r.q, 'O', 'HSO₄⁻ / H₂O removes H⁺: the ring is aromatic again.');
      title = title || 'Nitration';
      if (anilineN && c.has('H2SO4')) why.push('In conc. acid aniline is protonated (–NH₃⁺, a meta director): a lot of m-nitroaniline also forms, and the ring gets oxidised');
    } else if (halo) {
      const r = easStep(mech, p, `[${X}+]`, 1, `${X}₂ and Fe${X}₃ give ${X}⁺ (as ${X}–${X}–Fe${X}₃); the ring attacks it.`, [`${X}₂, Fe${X}₃`]);
      easFinish(mech, p, r.q, `[${X}-]`, `Fe${X}₄⁻ removes H⁺ and Fe${X}₃ is regenerated.`); title = `${X === 'Br' ? 'Bromination' : 'Chlorination'} (Lewis acid)`;
    } else if (fc) {
      if (c.any(['AcCl', 'Ac2O', 'BzCl', 'RCOCl'])) { const acy = c.has('BzCl') ? '[C+](=O)C1=CC=CC=C1' : c.has('RCOCl') ? '[C+](=O)' + (c.get('RCOCl').smi || 'C') : '[C+](=O)C'; const r = easStep(mech, p, acy, 1, 'AlCl₃ pulls Cl⁻ off the acyl chloride: the acylium ion R–C≡O⁺ is attacked by the ring.', ['AlCl₃']); easFinish(mech, p, r.q, '[Cl-]', 'AlCl₄⁻ removes H⁺. The ketone product is deactivated, so no second acylation.'); title = 'Friedel–Crafts acylation → aryl ketone'; why.push('No rearrangement, no polysubstitution'); }
      else { let rs = (c.get('RX') || {}).smi || 'C'; let note = ''; if (/^CCC+$/.test(rs)) { rs = 'C(C)' + rs.slice(2); note = 'The 1° carbocation rearranges (1,2-H shift) to the 2° one first, so the ring gets the branched group.'; } const cat = rs === 'C' ? '[CH3+]' : '[CH+]' + rs.slice(1).replace(/^\((\w+)\)/, '($1)'); const spec = rs === 'C' ? '[CH3+]' : rs.startsWith('C(C)(C)') ? '[C+](C)(C)' + rs.slice(7) : rs.startsWith('C(C)') ? '[CH+](C)' + rs.slice(4) : '[CH2+]' + rs.slice(1); const r = easStep(mech, p, spec, 1, 'AlCl₃ takes Cl⁻ from the alkyl halide: the carbocation is attacked by the ring.' + (note ? ' ' + note : ''), ['AlCl₃']); easFinish(mech, p, r.q, '[Cl-]', 'AlCl₄⁻ removes H⁺.'); title = 'Friedel–Crafts alkylation'; why.push('Alkyl groups activate the ring: polyalkylation is a problem'); if (note) why.push('Carbocations rearrange before attacking'); }
    } else if (c.any(['oleum', 'H2SO4'])) { const t = m.clone(); graft(t, p, 'S(=O)(=O)O'); return tx(m, t, { title: 'Sulphonation → arenesulphonic acid', text: 'SO₃ (from oleum or hot conc. H₂SO₄) is the electrophile. Sulphonation is reversible: steam removes –SO₃H again.', why, tags: ['SEAr'] }); }
    else return null;
    const o = done(mech, { title: title + (S0.subs.length ? ' (major isomer)' : ''), why, tags: ['SEAr'] });
    if (S0.subs.length && S0.ranked && S0.ranked.length > 1) { const alt = S0.ranked.find(([q]) => q !== p && S0.rel(q, S0.subs[0].a) !== S0.rel(p, S0.subs[0].a)); if (alt && S0.subs[0].op) o.why.push('The ortho isomer forms too (minor, more hindered)'); }
    return o;
  } });
function alkeneHint(m) { return M(m, 'C=C').filter(([a]) => !aromSet(m).has(a)).length ? '' : ''; }
reg({ id: 'nitro_red', name: 'Reduction of nitro compounds', fam: 'Aromatic', canon: 'Sn, HCl', db: ['Zinin reduction (selective)'],
  test: c => c.any(['Sn/HCl', 'Zn/HCl', 'SnCl2', 'Zn/NH4Cl']) || (c.has('H2') && c.has('Pd')) || (c.any(['Sn', 'Fe', 'Zn']) && c.has('HCl')),
  run(m, c) { const nx = nitroAny(m); if (!nx.length) return null; const t = m.clone(); const hydroxyl = c.has('Zn/NH4Cl');
    for (const [, n] of nx) { const os = heavyN(t, n).filter(x => el(t, x) === 'O'); t.removeAtom(os[0]); if (!hydroxyl) t.removeAtom(os[1]); else t.bondBetween(n, os[1]).o = 1; const na = t.atoms.get(n); na.chg = 0; if (hydroxyl) t.atoms.get(os[1]).chg = 0; }
    return tx(m, t, { title: hydroxyl ? 'Zn/NH₄Cl: nitro → hydroxylamine' : 'Nitro → amine', text: hydroxyl ? 'In neutral medium the reduction stops at Ar–NHOH.' : 'Sn/HCl (or Fe/HCl, H₂/Pd) reduces –NO₂ all the way to –NH₂ (aniline after a base work-up).', why: ['Fe/HCl is preferred industrially (FeCl₂ formed regenerates HCl)'], tags: ['red'] }); } });
reg({ id: 'diazotise', name: 'Amines + nitrous acid', fam: 'Amine', canon: 'NaNO₂, HCl, 0–5 °C', db: ['Diazotisation of aniline', 'Diazotisation of aliphatic 1° amine', 'Liebermann nitroso test'],
  test: c => c.has('HNO2'),
  run(m, c) {
    const ar = arNH2(m), al = M(m, '[CX4:1][NX3H2:2]'), sec = M(m, '[#6:1][NX3H1:2][#6:3]').filter(x => !E.isCarbonylC(m, x[0]) && !E.isCarbonylC(m, x[2]));
    if (ar.length) {
      const [, n] = ar[0]; const t = m.clone();
      if (c.heat || (c.T && c.T > 300)) { t.atoms.get(n).el = 'O'; return tx(m, t, { title: 'Warm: diazonium salt → phenol', text: 'Above 5 °C the diazonium salt decomposes in water: N₂ is lost and the ring takes –OH.', why: ['Diazotisation must be done at 0–5 °C'], tags: ['SN1Ar'] }); }
      t.atoms.get(n).chg = 1; graft(t, n, 'N', 3);
      return tx(m, t, { title: 'Diazotisation → arenediazonium salt (0–5 °C)', text: 'NaNO₂ + HCl make HNO₂ → NO⁺; the amine N attacks NO⁺, then tautomerism and loss of water give Ar–N≡N⁺ Cl⁻. Stable only when cold (resonance with the ring).', why: ['Only aromatic 1° amines give stable diazonium salts', 'Use it at once: Sandmeyer, coupling, etc.'], tags: ['diazotisation'] });
    }
    if (al.length) { const [cc, n] = al[0]; const t = m.clone(); t.atoms.get(n).el = 'O'; return tx(m, t, { title: 'Aliphatic 1° amine + HNO₂ → alcohol + N₂↑', text: 'The aliphatic diazonium ion is unstable: N₂ bubbles off at once and the carbocation is captured by water (rearranged products and alkenes also form).', why: ['N₂ evolution: test for aliphatic 1° amines'], tags: ['SN1'] }); }
    if (sec.length) { const [, n] = sec[0]; const t = m.clone(); graft(t, n, 'N=O'); return tx(m, t, { title: '2° amine → N-nitrosamine (yellow oil)', text: 'The N–H is replaced by N–N=O (Liebermann nitroso test).', why: [], tags: ['test'] }); }
    const dm = M(m, '[c:1][NX3H0:2]([CH3])[CH3]'); if (dm.length) { const rings = C.aromaticRings(m); const S = sites(m, rings[0]); const t = m.clone(); graft(t, S.list[0], 'N=O'); return tx(m, t, { title: '3° aromatic amine → p-nitroso compound', text: 'NO⁺ is a weak electrophile, but the –NMe₂ ring is activated enough: substitution at para.', why: [], tags: ['SEAr'] }); }
    return null;
  } });
const DZ = [['CuCN', 'C#N', 'Sandmeyer: Ar–CN'], ['CuBr', 'Br', 'Sandmeyer: Ar–Br'], ['CuCl', 'Cl', 'Sandmeyer: Ar–Cl'], ['NaI', 'I', 'KI: Ar–I'], ['HBF4', 'F', 'Balz–Schiemann: Ar–F'], ['H3PO2', null, 'Deamination: Ar–H'], ['ROH', null, 'Deamination: Ar–H']];
reg({ id: 'diazo_rx', name: 'Reactions of diazonium salts', fam: 'Amine', canon: 'CuCN', db: ['Sandmeyer reaction', 'Gattermann reaction', 'Balz–Schiemann reaction', 'Deamination of diazonium', 'Azo coupling'],
  test: c => DZ.some(z => c.has(z[0])) || c.any(['phenol', 'aniline', 'naphthol', 'H2O', 'H3O+', 'Cu']),
  run(m, c) {
    const dz = diazo(m); if (!dz.length) return null; const [ca, n1] = dz[0]; const t = m.clone(); const d = dir(t, ca, n1); cut(t, ca, n1);
    const z = DZ.find(z => c.has(z[0]));
    if (z) { if (z[1]) graft(t, ca, z[1], 1, d); return tx(m, t, { title: z[2], text: z[0].startsWith('Cu') ? 'Cu(I) salts transfer the nucleophile through an aryl radical; N₂ is lost.' : z[0] === 'HBF4' ? 'The diazonium fluoroborate is isolated and heated: Ar–F + N₂ + BF₃.' : z[0] === 'NaI' ? 'I⁻ is a strong enough nucleophile on its own, no Cu needed.' : 'H₃PO₂ (or ethanol) replaces –N₂⁺ by H.', why: ['N₂ is an excellent leaving group'], tags: ['radical'] }); }
    if (c.has('Cu') && c.has('HCl')) { graft(t, ca, 'Cl', 1, d); return tx(m, t, { title: 'Gattermann: Ar–Cl', text: 'Copper powder with HCl does the Sandmeyer job.', why: [] }); }
    const cp = ['phenol', 'aniline', 'naphthol'].find(k => c.has(k));
    if (cp) { const smi = { phenol: 'OC1=CC=CC=C1', aniline: 'NC1=CC=CC=C1', naphthol: 'OC1=CC=C2C=CC=CC2=C1' }[cp]; const pm = C.fromSmiles(smi, BL); const prs = C.aromaticRings(pm); const S = sites(pm, prs[0]); let target = S.list[0];
      if (cp === 'naphthol') { /* 2-naphthol couples at C1: next to C–OH and to the ring fusion */ const oA = [...pm.atoms.values()].find(a => a.el === 'O'); const cO = pm.nbrs(oA.id)[0]; const rs = C.rings(pm, 6); const fus = new Set([...pm.atoms.keys()].filter(id => rs.filter(r => r.includes(id)).length > 1));
        const c1 = pm.nbrs(cO).find(x => x !== oA.id && pm.nbrs(x).some(y => fus.has(y)) && !fus.has(x)); if (c1 != null) target = c1; }
      const g = graft(t, ca, 'N=N', 1, d); const n2 = [...g.map.values()].find(id => id !== g.first);
      const map = E.placeSpecies(t, pm, n2, { att: target, dist: BL }); t.addBond(n2, map.get(target), 1); E.tidyLongBonds(t);
      return tx(m, t, { title: `Azo coupling with ${cp} → ${cp === 'aniline' ? 'yellow' : 'orange'} azo dye`, text: `The diazonium ion is a weak electrophile; it attacks the activated ring at para (${cp === 'phenol' ? 'mildly alkaline' : cp === 'aniline' ? 'mildly acidic' : 'alkaline'} solution).`, why: ['–N=N– joins two rings: coloured compounds'], tags: ['SEAr'] }); }
    if (c.any(['H2O', 'H3O+'])) { graft(t, ca, 'O', 1, d); return tx(m, t, { title: 'Warm water: Ar–OH', text: 'Warming the diazonium solution: N₂ leaves, water attacks (SN1Ar): phenol.', why: [] }); }
    return null;
  } });
reg({ id: 'phenol_rx', name: 'Phenol reactions', fam: 'Phenol', canon: 'CHCl₃, NaOH', db: ['Reimer–Tiemann reaction', 'Kolbe–Schmitt reaction', 'Neutral FeCl3 test'],
  test: c => c.any(['CHCl3', 'CO2', 'Zn', 'FeX3', 'CCl4']),
  run(m, c) {
    const ph = arOH(m); if (!ph.length) return null; const [ca, o] = ph[0]; const ring = C.aromaticRings(m).find(r => r.includes(ca)) || [...aromSet(m)];
    const ortho = m.nbrs(ca).filter(n => ring.includes(n) && hTot(m, n) > 0);
    if (c.has('CHCl3') || c.has('CCl4')) { if (!ortho.length) return null; const t = m.clone(); t.atoms.get(o).chg = 0; graft(t, ortho[0], c.has('CCl4') ? 'C(=O)O' : 'C=O'); return tx(m, t, { title: c.has('CCl4') ? 'Reimer–Tiemann with CCl₄ → salicylic acid' : 'Reimer–Tiemann → salicylaldehyde (ortho –CHO)', text: 'CHCl₃ + OH⁻ give dichlorocarbene (:CCl₂), the electrophile. The phenoxide ring attacks it at ortho; hydrolysis of the –CHCl₂ gives –CHO.', why: ['Ortho product major'], tags: ['carbene'] }); }
    if (c.has('CO2')) { if (!ortho.length) return null; const t = m.clone(); t.atoms.get(o).chg = 0; graft(t, ortho[0], 'C(=O)O'); return tx(m, t, { title: 'Kolbe–Schmitt → salicylic acid', text: 'Sodium phenoxide is carboxylated by CO₂ (≈400 K, 4–7 atm) at the ortho position; acidification gives salicylic acid.', why: ['Na⁺ salt: ortho; K⁺ salt gives more para'], tags: ['SEAr'] }); }
    if (c.has('Zn')) { const t = m.clone(); t.removeAtom(o); return tx(m, t, { title: 'Zinc dust distillation → arene', text: 'Heating with Zn dust removes the –OH oxygen: phenol → benzene.', why: [] }); }
    if (c.has('FeX3')) return none('Violet colour (positive FeCl₃ test)', 'Phenols (and enols) give a violet/blue complex with neutral FeCl₃: that is the test for a phenolic –OH.');
    return null;
  } });
reg({ id: 'amine_rx', name: 'Amine reactions', fam: 'Amine', canon: 'CHCl₃, alc. KOH', db: ['Carbylamine reaction (isocyanide test)', 'Hofmann ammonolysis'],
  test: c => (c.has('CHCl3') && c.any(['KOH', 'alcKOH', 'aqKOH'])) || c.has('RX'),
  run(m, c) {
    const p = M(m, '[#6:1][NX3H2:2]').filter(x => !E.isCarbonylC(m, x[0]));
    if (c.has('CHCl3')) { if (!p.length) { if (M(m, '[NX3;!$(NC=O)]').length) return none('Negative carbylamine test', 'Only 1° amines (–NH₂) give the foul-smelling isocyanide.'); return null; } const [, n] = p[0]; const t = m.clone(); t.atoms.get(n).chg = 1; const pp = P(t, n), d = freeDir(t, n); const cc = t.addAtom('C', pp[0] + d[0] * BL, pp[1] + d[1] * BL, { chg: -1 }); t.addBond(n, cc, 3); return tx(m, t, { title: 'Carbylamine test → isocyanide (foul smell)', text: 'CHCl₃ + KOH give :CCl₂; the 1° amine N adds to it and two HCl are lost: R–N≡C.', why: ['Test for 1° amines (aliphatic and aromatic)'], tags: ['test'] }); }
    const am = M(m, '[NX3;!$(NC=O);!$(N=*);!$(N#*):1]').filter(([n]) => hTot(m, n) > 0 || true); if (!am.length || M(m, '[CX4][Cl,Br,I]').length) return null;
    const [n] = am[0]; const rs = (c.get('RX') || {}).smi || 'C'; const t = m.clone(); const excess = c.flag('excess');
    if (excess) { const k = hTot(t, n); for (let i = 0; i < k + 1; i++) graft(t, n, rs); t.atoms.get(n).chg = 1; return tx(m, t, { title: 'Exhaustive alkylation → quaternary ammonium salt', text: 'With excess R–X every N–H is replaced and a fourth R adds: R₄N⁺ X⁻.', why: [], tags: ['SN2'] }); }
    graft(t, n, rs); return tx(m, t, { title: 'N-Alkylation (SN2)', text: 'The amine lone pair attacks the alkyl halide; a mixture of 2°, 3° amines and the quaternary salt usually forms.', why: ['Hard to stop at one alkylation'], tags: ['SN2'] });
  } });
reg({ id: 'epox_open', name: 'Epoxide ring opening', fam: 'Ether', canon: 'H₃O⁺', db: ['Acid hydrolysis of ethers'],
  test: c => c.any(['H3O+', 'H2SO4', 'NaOMe', 'NaOEt', 'NH3', 'H2O']),
  run(m, c) { const ep = epox(m); if (!ep.length) return null; const [c1, o2, c3] = ep[0]; const t = m.clone(); const acid = c.any(['H3O+', 'H2SO4', 'H2O']);
    const more = cN(t, c1).length >= cN(t, c3).length ? c1 : c3, less = more === c1 ? c3 : c1; const at = acid ? more : less;
    t.removeBond(t.bondBetween(at, o2).id); graft(t, at, acid ? 'O' : c.has('NH3') ? 'N' : c.has('NaOMe') ? 'OC' : 'OCC');
    return tx(m, t, { title: acid ? 'Acid opening → trans-1,2-diol' : 'Base / nucleophile opening at the less hindered carbon', text: acid ? 'H⁺ on the ring O; water attacks the more substituted carbon (more δ+) from the back: anti opening.' : 'A strong nucleophile does SN2 at the less hindered ring carbon.', why: [acid ? 'Acid: more substituted carbon' : 'Base: less substituted carbon'], tags: [acid ? 'SN1' : 'SN2'] }); } });
reg({ id: 'acidbase', name: 'Acid–base', fam: 'Acid–base', canon: 'NaHCO₃', db: ['Acid–base extraction'],
  test: c => c.any(['NaHCO3', 'KOH', 'aqKOH']),
  run(m, c) {
    const ac = acids(m), so = M(m, '[SX4](=O)(=O)[OX2H1:1]'), ph = arOH(m).filter(([, o]) => m.atoms.get(o).chg === 0);
    const t = m.clone();
    if (c.has('NaHCO3')) { if (!ac.length && !so.length) return none('No reaction with NaHCO₃', ph.length ? 'Phenol is too weak an acid to release CO₂ from NaHCO₃ (that is how you tell it from a carboxylic acid).' : 'Nothing here is acidic enough.'); for (const x of ac) t.atoms.get(x[2]).chg = -1; for (const x of so) t.atoms.get(x[0]).chg = -1; return tx(m, t, { title: 'Brisk effervescence: carboxylate + CO₂↑', text: 'Carboxylic (and sulphonic) acids are stronger than H₂CO₃, so they release CO₂ from NaHCO₃.', why: ['Test for –COOH'], tags: ['acidbase'] }); }
    if (!ac.length && !ph.length && !so.length) { if (M(m, '[CX4][OX2H1]').length) return none('No reaction with aq. NaOH', 'Alcohols are weaker acids than water, so NaOH cannot deprotonate them.'); return null; }
    for (const x of ac) t.atoms.get(x[2]).chg = -1; for (const x of ph) t.atoms.get(x[1]).chg = -1; for (const x of so) t.atoms.get(x[0]).chg = -1;
    return tx(m, t, { title: ac.length ? 'Neutralisation → carboxylate salt' : 'Phenol + NaOH → sodium phenoxide', text: ac.length ? 'R–COOH + NaOH → R–COO⁻Na⁺ + H₂O.' : 'Phenol is acidic enough (pKa ≈ 10) to react with NaOH; the phenoxide is resonance-stabilised.', why: ['Acidity: RCOOH > phenol > water > alcohol'], tags: ['acidbase'] });
  } });
reg({ id: 'rosenmund', name: 'Rosenmund reduction', fam: 'Acid derivative', canon: 'H₂, Pd/BaSO₄', db: ['Rosenmund reduction'],
  test: c => c.has('H2') && (c.has('Lindlar') || /baso4/i.test(c.raw || '')),
  run(m, c) {
    const ax = acylX(m); if (!ax.length) return null; const [c1, , x] = ax[0]; const t = m.clone(); t.removeAtom(x);
    return tx(m, t, { title: 'Rosenmund: acyl chloride → aldehyde', text: 'H₂ on Pd poisoned with BaSO₄ (and quinoline) replaces Cl by H and stops at the aldehyde, because the poisoned catalyst cannot reduce the C=O further.', why: ['Poisoned Pd: stops at the aldehyde', 'Formaldehyde cannot be made this way (HCOCl is unstable)'], tags: ['red'] });
  } });
})(typeof window !== 'undefined' ? window : globalThis);
