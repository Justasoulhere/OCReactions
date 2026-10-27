/* engine_r4.js: second batch from the coverage audit (Oct 12, late). FRONT-registered, so these specific
   reactions win over the broad recipes (aldol, Cannizzaro, SN2, hydration...) that used to take them. */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore;
const { BL, el, P, unit, out } = E;
const REC = E.REC; const front = r => REC.unshift(r);
const M = (m, sm) => C.matchAll(m, sm, {});
const tidy = m => { const t = m.clone(); try { C.relayout(t, BL); } catch (e) { } return t; };
const fresh = smi => C.fromSmilesMB(smi, BL);
const step = (L, R, text, cond) => ({ L: L.clone(), arrows: [], R: R.clone(), text, cond: cond || [] });
const later = (c, all, k) => [c, ...((all || []).slice((k || 0) + 1))];
const anyLater = (c, all, k, ks) => later(c, all, k).some(s => s.any(ks));
const OHm = ['KOH', 'aqKOH'];
const hT = (m, id) => E.hTot(m, id);
const linear = (t, from, to, elx) => { const p0 = P(t, from), p1 = P(t, to); const d = unit([p0[0] - p1[0], p0[1] - p1[1]]); const id = t.addAtom(elx, p0[0] + d[0] * BL, p0[1] + d[1] * BL); t.addBond(from, id, 1); return id; };
const mainOf = (t, keep) => { const f = t.frags().find(x => x.includes(keep)) || t.frags()[0]; return E.cleanProduct(t, f); };
const othersOf = (t, keep) => t.frags().filter(x => !x.includes(keep)).map(f => E.cleanProduct(t, f)).filter(x => x && x.atoms.size > 1 || (x && [...x.atoms.values()].some(a => a.el === 'C')));
function attachR(t, at, Rs) { const sp = E.species(Rs); const pa = P(t, at); const d = E.freeDir(t, at); const idm = t.merge(sp, pa[0] + d[0] * BL, pa[1] + d[1] * BL); t.addBond(at, idm.get(1), 1); return idm.get(1); }

/* terminal alkyne + Na / liq. NH3: deprotonation wins over reduction */
front({ id: 'na_nh3_terminal', name: 'Terminal alkyne + Na/NH₃', fam: 'Alkyne', canon: 'Na, liq. NH₃', db: [],
  test: c => c.has('Na/NH3') || (c.has('Na') && c.has('NH3')),
  run(m) { const h = M(m, '[CH1]#[#6]')[0]; if (!h) return null; const t = m.clone(); linear(t, h[0], h[1], 'Na');
    return out({ title: 'Terminal alkyne + Na in liquid NH₃ → sodium acetylide (+ ½ H₂)', product: t, steps: [step(m, t, 'The ≡C–H is acidic (sp carbon), so Na removes it before any reduction: R–C≡C⁻ Na⁺.', ['Na, liq. NH₃'])], why: ['Internal alkynes are reduced to trans-alkenes instead', 'Terminal alkynes form acetylides first'], tags: ['acidbase'] }); } });

/* nucleophilic (Weitz-Scheffer) epoxidation of enones */
front({ id: 'weitz', name: 'Nucleophilic epoxidation of an enone', fam: 'Carbonyl', canon: 'H₂O₂, NaOH', db: [],
  test: c => c.has('H2O2/OH') || (c.has('H2O2') && c.any(OHm)),
  run(m) { const h = M(m, '[CX3]=[CX3][CX3]=O')[0]; if (!h) return null; const [cb, ca] = h;
    const t = m.clone(); t.bondBetween(cb, ca).o = 1; const pb = P(t, cb), pa = P(t, ca); const d = E.freeDir(t, cb); const o = t.addAtom('O', (pb[0] + pa[0]) / 2 + d[0] * BL * .7, (pb[1] + pa[1]) / 2 + d[1] * BL * .7); t.addBond(cb, o, 1); t.addBond(ca, o, 1); const R = tidy(t);
    return out({ title: 'H₂O₂/OH⁻ on an enone → α,β-epoxy ketone', product: R, steps: [step(m, R, 'HOO⁻ adds conjugately (1,4) to the β-carbon; the enolate then closes the ring, pushing out OH⁻.', ['H₂O₂, OH⁻'])], why: ['Electron-poor C=C: needs a nucleophilic oxidant, not a peracid', 'mCPBA would attack an electron-rich C=C instead'], tags: ['AdN'] }); } });

/* chloral + OH-: haloform-type cleavage */
front({ id: 'chloral', name: 'Chloral + base', fam: 'Carbonyl', canon: 'NaOH', db: ['Haloform reaction'],
  test: c => c.any(OHm),
  run(m) { if (!M(m, '[CX3H1](=O)C(Cl)(Cl)Cl').length) return null; const A1 = fresh('ClC(Cl)Cl'), B1 = fresh('O=C[O-]');
    return out({ title: 'Chloral + OH⁻ → chloroform + formate', product: A1, extra: [B1], steps: [step(m, A1, 'OH⁻ adds to C=O; the tetrahedral intermediate pushes out CCl₃⁻ (stabilised by three Cl), which takes H⁺: CHCl₃.', ['OH⁻'])], why: ['Haloform-type C–C cleavage', 'Used to make pure chloroform'], tags: ['AdN'] }); } });

/* benzilic acid rearrangement / intramolecular Cannizzaro of 1,2-dicarbonyls */
front({ id: 'benzilic', name: 'Benzilic acid rearrangement', fam: 'Carbonyl', canon: 'conc. KOH', db: ['Benzilic acid rearrangement', 'Cannizzaro reaction'],
  test: c => c.any(OHm) && !c.flag('dil'),
  run(m) {
    const h = M(m, '[CX3](=O)[CX3]=O')[0]; if (!h) return null; let [c1, o1, c2, o2] = h;
    const ald = id => hT(m, id) > 0; if (ald(c2) && !ald(c1)) { [c1, o1, c2, o2] = [c2, o2, c1, o1]; }
    const t = m.clone(); const R1 = m.nbrs(c1).find(n => n !== c2 && n !== o1);
    if (ald(c1)) { t.bondBetween(c2, o2).o = 1; } else { if (R1 == null) return null; t.removeBond(t.bondBetween(R1, c1).id); t.addBond(R1, c2, 1); t.bondBetween(c2, o2).o = 1; }
    const p1 = P(t, c1); const d = E.freeDir(t, c1); const om = t.addAtom('O', p1[0] + d[0] * BL, p1[1] + d[1] * BL, { chg: -1 }); t.addBond(c1, om, 1); const R = tidy(t);
    const intra = ald(c1);
    return out({ title: intra ? 'Intramolecular Cannizzaro: 1,2-dialdehyde → α-hydroxy carboxylate' : 'Benzilic acid rearrangement: 1,2-diketone → α-hydroxy carboxylate', product: R, steps: [step(m, R, intra ? 'OH⁻ adds to one CHO; its H shifts as hydride to the other CHO (internal Cannizzaro).' : 'OH⁻ adds to one C=O; a group migrates with its electrons to the neighbouring C=O (1,2-shift), giving the α-hydroxy acid salt.', ['OH⁻'])], why: ['Needs a 1,2-dicarbonyl', 'Benzil → benzilic acid is the classic example'], tags: ['rearrangement'] });
  } });

/* base dehydration of an aldol (E1cB) */
front({ id: 'aldol_dehyd', name: 'Aldol dehydration', fam: 'Carbonyl', canon: 'NaOH, Δ', db: ['Aldol reaction / condensation'],
  test: c => c.any([...OHm, 'NaOEt']),
  run(m) { const h = M(m, '[CX3](=O)[CX4;!H0][CX4][OX2H1]')[0]; if (!h) return null; const [, , ca, cb, oh] = h;
    const t = m.clone(); t.removeAtom(oh); t.bondBetween(ca, cb).o = 2;
    return out({ title: 'Aldol dehydration (E1cB) → α,β-unsaturated carbonyl', product: t, steps: [step(m, t, 'OH⁻ removes the acidic α-H; the enolate expels OH⁻ from the β-carbon (E1cB): the C=C forms in conjugation with C=O.', ['OH⁻'])], why: ['Conjugation drives the dehydration', 'Heat usually needed'], tags: ['E1cB'] }); } });

/* enol ether hydrolysis */
front({ id: 'enol_ether', name: 'Enol ether hydrolysis', fam: 'Ether', canon: 'H₃O⁺', db: [],
  test: c => c.any(['H3O+', 'H2SO4']) && !c.any(['HI', 'HBr', 'HCl']),
  run(m) { const h = M(m, '[CX3]=[CX3][OX2][CX4]')[0]; if (!h) return null; const [cb, ca, o, r] = h;
    const t = m.clone(); t.removeBond(t.bondBetween(o, r).id); t.bondBetween(cb, ca).o = 1; t.bondBetween(ca, o).o = 2; E.addOn(t, r, 'O');
    const A1 = tidy(mainOf(t, ca)); const B1 = othersOf(t, ca);
    return out({ title: 'Enol ether + H₃O⁺ → carbonyl + alcohol', product: A1, extra: B1, steps: [step(m, A1, 'H⁺ adds to the C=C carbon (giving an oxocarbenium ion); water adds; the hemiacetal falls apart to C=O + ROH.', ['H₃O⁺'])], why: ['Vinyl ethers are acid-sensitive', 'Markovnikov protonation on carbon'], tags: ['hydrolysis'] }); } });

/* vicinal dihalide + NaI or Zn: dehalogenation */
front({ id: 'vic_dehalo', name: 'Dehalogenation of vicinal dihalides', fam: 'Alkyl halide', canon: 'NaI, acetone', db: [],
  test: c => c.any(['NaI', 'Zn']) && !c.any(['HCl', 'H3O+']),
  run(m, c) { const h = M(m, '[CX4]([Cl,Br,I])[CX4][Cl,Br,I]')[0]; if (!h) return null; const [c1, x1, c2, x2] = h;
    const t = m.clone(); t.removeAtom(x1); t.removeAtom(x2); t.bondBetween(c1, c2).o = 2;
    return out({ title: `Vicinal dihalide + ${c.has('Zn') ? 'Zn' : 'NaI (acetone)'} → alkene`, product: t, steps: [step(m, t, c.has('Zn') ? 'Zn removes both halogens (anti elimination): ZnX₂ + alkene.' : 'I⁻ attacks one halogen; the other leaves anti-periplanar: alkene + IX + X⁻.', [c.has('Zn') ? 'Zn' : 'NaI, acetone'])], why: ['Anti elimination: stereospecific', 'Used to protect/regenerate a C=C'], tags: ['E2'] }); } });

/* Hofmann elimination (exhaustive methylation, Ag2O, heat) */
front({ id: 'hofmann_elim', name: 'Hofmann elimination', fam: 'Amine', canon: 'excess CH₃I then Ag₂O, H₂O, Δ', db: ['Hofmann elimination'],
  eats: s => s.any(['Ag2O', 'H2O']),
  test: (c, all, k) => anyLater(c, all, k, ['Ag2O']),
  run(m) {
    const n = [...m.atoms.values()].find(a => a.el === 'N' && !M(m, '[NX3][CX3]=O').some(x => x[0] === a.id)); if (!n) return null;
    const cs = m.nbrs(n.id).filter(x => el(m, x) === 'C'); const ca = cs.sort((p, q) => m.heavyDegree(q) - m.heavyDegree(p))[0]; if (ca == null || m.heavyDegree(ca) < 2) return null;
    const bs = m.nbrs(ca).filter(x => x !== n.id && el(m, x) === 'C' && hT(m, x) > 0 && m.bondBetween(ca, x).o === 1); if (!bs.length) return null;
    const cb = bs.sort((p, q) => hT(m, q) - hT(m, p))[0];
    const t = m.clone(); for (const id of E.subtree(t, n.id, ca)) t.removeAtom(id); t.bondBetween(ca, cb).o = 2;
    return out({ title: 'Hofmann elimination → less substituted alkene + N(CH₃)₃', product: t, extra: [fresh('CN(C)C')], steps: [step(m, t, 'Excess CH₃I makes the quaternary ammonium salt; Ag₂O/H₂O swaps in OH⁻; on heating OH⁻ removes the most accessible β-H (E2) and N(CH₃)₃ leaves.', ['CH₃I', 'Ag₂O, Δ'])], why: ['Hofmann rule: least substituted alkene (bulky leaving group, acidic primary H)', 'Opposite of Saytzeff'], tags: ['E2'] });
  } });

/* acyl halide + alkoxide: ester (beats SN2 on an alkyl halide elsewhere) */
front({ id: 'acylx_alkoxide', name: 'Acyl halide + alkoxide', fam: 'Acid derivative', canon: 'CH₃ONa', db: [],
  test: c => c.any(['NaOMe', 'NaOEt', 'NaOR']),
  run(m, c) { const h = M(m, '[CX3](=O)[Cl,Br,I]')[0]; if (!h) return null; const x = h[2];
    const R = c.has('NaOMe') ? 'C' : c.has('NaOEt') ? 'CC' : ((c.get('NaOR') || {}).smi || 'C');
    const t = m.clone(); const a = t.atoms.get(x); a.el = 'O'; a.hFix = null; attachR(t, x, R); const Pm = tidy(t);
    return out({ title: 'Acyl halide + alkoxide → ester', product: Pm, steps: [step(m, Pm, 'RO⁻ adds to the very reactive acyl C=O and X⁻ leaves (addition–elimination). Acyl substitution is far faster than SN2 at any alkyl halide in the molecule.', [R === 'C' ? 'CH₃O⁻' : 'RO⁻'])], why: ['Acid chlorides are the most reactive acid derivatives'], tags: ['AdN-E'] }); } });

/* intramolecular aldol of 1,4- / 1,5-dicarbonyls */
front({ id: 'intra_aldol', name: 'Intramolecular aldol', fam: 'Carbonyl', canon: 'NaOH, Δ', db: ['Aldol reaction / condensation'],
  test: c => c.any([...OHm, 'NaOEt']) && !c.any(['RX', 'I2', 'Br2']),
  run(m, c) {
    const cos = M(m, '[CX3;$([CX3]([#6])[#6]),$([CX3H1][#6])]=O').map(x => x[0]); if (cos.length < 2) return null;
    const dist = (a) => { const d = new Map([[a, 0]]), q = [a]; while (q.length) { const x = q.shift(); for (const y of m.nbrs(x)) if (!d.has(y)) { d.set(y, d.get(x) + 1); q.push(y); } } return d; };
    let best = null;
    for (const c1 of cos) for (const a of m.nbrs(c1)) { if (el(m, a) !== 'C' || hT(m, a) < 1 || m.bondsOf(a).some(b => b.o > 1)) continue; const d = dist(a);
      for (const c2 of cos) { if (c2 === c1) continue; const ring = (d.get(c2) || 0) + 1; if (ring !== 5 && ring !== 6) continue; const sc = (ring === 6 ? 2 : 1) + (hT(m, a) >= 2 ? .5 : 0); if (!best || sc > best[0]) best = [sc, a, c2, ring]; } }
    if (!best) return null; const [, a, c2] = best; const o2 = m.nbrs(c2).find(x => el(m, x) === 'O' && m.bondBetween(c2, x).o === 2);
    const t = m.clone(); t.addBond(a, c2, 1); t.bondBetween(c2, o2).o = 1;
    const dehyd = c.heat && hT(m, a) >= 2; if (dehyd) { t.removeAtom(o2); t.bondBetween(a, c2).o = 2; }
    const R = tidy(t);
    return out({ title: dehyd ? `Intramolecular aldol condensation → ${best[3]}-ring enone` : `Intramolecular aldol → ${best[3]}-membered ring`, product: R, steps: [step(m, R, `The enolate of one C=O attacks the other C=O in the same molecule, closing a ${best[3]}-membered ring${dehyd ? '; heat then removes water to give the conjugated enone' : ''}.`, ['OH⁻'])], why: ['5- and 6-membered rings form fastest', 'Hexane-2,5-dione → 3-methylcyclopent-2-enone'], tags: ['enolate'] });
  } });

/* Curtius (acyl azide) and Lossen (hydroxamic acid): R-CO-X → R-NH2 */
front({ id: 'curtius_lossen', name: 'Curtius / Lossen rearrangement', fam: 'Acid derivative', canon: 'NaN₃ then Δ, H₂O', db: ['Curtius rearrangement', 'Lossen rearrangement'],
  eats: s => s.any(['H2O', 'H3O+', 'KOH', 'aqKOH']),
  test: (c, all, k) => anyLater(c, all, k, ['NaN3']) || (c.has('NH2OH') && anyLater(c, all, k, OHm)),
  run(m, c, all, k) { const h = M(m, '[#6][CX3](=O)[Cl,Br,I,OX2H1]')[0]; if (!h) return null; const [r, cc] = h; const cur = anyLater(c, all, k, ['NaN3']);
    const t = m.clone(); for (const id of E.subtree(t, cc, r)) t.removeAtom(id); E.addOn(t, r, 'N'); const R = tidy(t);
    return out({ title: cur ? 'Curtius rearrangement: R–COCl → R–NH₂ (one C fewer)' : 'Lossen rearrangement: hydroxamic acid → R–NH₂ (one C fewer)', product: R, extra: [fresh('O=C=O')], steps: [step(m, R, cur ? 'Azide gives the acyl azide; on heating N₂ leaves and R migrates to N: isocyanate R–N=C=O; water hydrolyses it to the amine + CO₂.' : 'NH₂OH gives the hydroxamic acid (acylated on O); base triggers R migration to N with loss of the carboxylate: isocyanate, then amine + CO₂.', [cur ? 'NaN₃, Δ' : 'NH₂OH, OH⁻'])], why: ['Same isocyanate route as Hofmann bromamide', 'Product has one carbon fewer'], tags: ['rearrangement'] }); } });

/* HI / red P: complete reduction (glucose → n-hexane) */
front({ id: 'hi_redp', name: 'HI / red P reduction', fam: 'Reduction', canon: 'HI, red P, Δ', db: [],
  test: c => c.has('HI') && c.has('redP'),
  run(m) { const at = [...m.atoms.values()]; if (!at.some(a => a.el === 'O') || at.some(a => !['C', 'O', 'H'].includes(a.el))) return null;
    const t = m.clone(); for (const a of [...t.atoms.values()]) if (a.el === 'O') t.removeAtom(a.id); for (const b of t.bonds.values()) b.o = 1; for (const a of t.atoms.values()) a.hFix = null;
    if (t.frags().length > 1) return null; const R = tidy(t);
    return out({ title: 'HI / red P, Δ: every C–O reduced → alkane', product: R, steps: [step(m, R, 'Prolonged heating with HI and red P replaces every –OH and C=O oxygen by H.', ['HI, red P, Δ'])], why: ['Glucose → n-hexane: proves the straight 6-carbon chain'], tags: ['reduction'] }); } });

/* exhaustive acetylation of polyols / sugars */
front({ id: 'peracetyl', name: 'Acetylation of every –OH', fam: 'Acid derivative', canon: '(CH₃CO)₂O', db: ['Acetylation of –OH / –NH₂'],
  test: c => c.any(['Ac2O', 'AcCl', 'RCOCl']),
  run(m) { const ohs = M(m, '[CX4][OX2H1]').map(x => x[1]); if (ohs.length < 3) return null;
    const t = m.clone(); for (const o of ohs) { const c1 = E.addOn(t, o, 'C'); E.addOn(t, c1, 'O', 2); E.addOn(t, c1, 'C'); } const R = tidy(t);
    return out({ title: `Acetylation of all ${ohs.length} –OH groups (${ohs.length === 5 ? 'pentaacetate' : 'polyacetate'})`, product: R, steps: [step(m, R, 'Each –OH attacks acetic anhydride; acetate leaves. With excess Ac₂O every –OH is acetylated.', ['(CH₃CO)₂O'])], why: ['Glucose pentaacetate: proves five –OH groups'], tags: ['acylation'] }); } });

/* formaldehyde + ammonia: urotropine */
front({ id: 'urotropine', name: 'Formaldehyde + NH₃', fam: 'Carbonyl', canon: 'NH₃', db: [],
  test: c => c.has('NH3') && !c.any(['Na', 'Li', 'Na/NH3']),
  run(m) { if (m.atoms.size !== 2 || !M(m, '[CH2]=O').length) return null; const R = fresh('C1N2CN3CN1CN(C2)C3');
    return out({ title: 'HCHO + NH₃ → urotropine (hexamethylenetetramine)', product: R, extra: [fresh('O')], steps: [step(m, R, '6 HCHO + 4 NH₃ condense (losing 6 H₂O) into the cage compound (CH₂)₆N₄.', ['NH₃'])], why: ['Formaldehyde does not give a simple imine with NH₃', 'Urotropine is used as a urinary antiseptic'], tags: ['condensation'] }); } });

/* Perkin reaction */
front({ id: 'perkin', name: 'Perkin reaction', fam: 'Carbonyl', canon: '(CH₃CO)₂O, CH₃COONa, Δ', db: ['Perkin reaction'],
  test: c => c.has('Ac2O') && c.any(['NaOAc', 'KOAc']),
  run(m) { const h = M(m, '[CX3H1](=O)c')[0]; if (!h) return null; const [c0, o] = h;
    const t = m.clone(); t.removeAtom(o); const c2 = E.addOn(t, c0, 'C', 2); const c3 = E.addOn(t, c2, 'C'); E.addOn(t, c3, 'O', 2); E.addOn(t, c3, 'O'); const R = tidy(t);
    return out({ title: 'Perkin reaction: ArCHO → cinnamic acid', product: R, steps: [step(m, R, 'Acetate removes an α-H of the anhydride; the enolate adds to ArCHO; dehydration and hydrolysis give Ar–CH=CH–COOH.', ['(CH₃CO)₂O, CH₃COONa, Δ'])], why: ['Aromatic aldehydes only (no α-H of their own)', 'Product: β-aryl acrylic acid'], tags: ['condensation'] }); } });

/* Wurtz-Fittig cross coupling */
front({ id: 'wurtz_fittig', name: 'Wurtz–Fittig reaction', fam: 'Aromatic', canon: 'CH₃Br, Na, dry ether', db: ['Fittig reaction'],
  test: c => c.has('Na') && c.has('RX') && (c.flag('dry') || c.flag('ether')),
  run(m, c) { const h = M(m, 'c[Cl,Br,I]')[0]; if (!h) return null; const [ar, x] = h; const Rs = (c.get('RX') || {}).smi || 'C';
    const t = m.clone(); t.removeAtom(x); attachR(t, ar, Rs); const R = tidy(t);
    return out({ title: 'Wurtz–Fittig: Ar–X + R–X + Na → Ar–R', product: R, steps: [step(m, R, 'Na couples the aryl and alkyl groups (dry ether). Side products: R–R and Ar–Ar.', ['Na, dry ether'])], why: ['Mixed coupling: Wurtz (R–R) + Fittig (Ar–Ar) side products'], tags: ['coupling'] }); } });

/* DIBAL-H at low temperature */
front({ id: 'dibal', name: 'DIBAL-H reduction', fam: 'Reduction', canon: 'DIBAL-H, −78 °C then H₂O', db: [],
  eats: s => s.any(['H2O', 'H3O+']),
  test: c => c.has('DIBAL'),
  run(m) { let h = M(m, '[#6][CX3](=O)[OX2][#6]')[0];
    if (h) { const [, cc, , o2] = h; const t = m.clone(); const leave = E.subtree(t, o2, cc); const alc = E.cleanProduct(m, leave); for (const id of leave) t.removeAtom(id);
      return out({ title: 'DIBAL-H (−78 °C): ester → aldehyde', product: t, extra: alc ? [alc] : [], steps: [step(m, t, 'One hydride adds; at −78 °C the tetrahedral aluminium intermediate is stable and only collapses to the aldehyde on work-up.', ['DIBAL-H, −78 °C'])], why: ['Stops at the aldehyde (LiAlH₄ would give the alcohol)'], tags: ['reduction'] }); }
    h = M(m, '[#6]C#N')[0]; if (!h) return null; const t = m.clone(); const n = t.atoms.get(h[2]); n.el = 'O'; t.bondBetween(h[1], h[2]).o = 2;
    return out({ title: 'DIBAL-H: nitrile → aldehyde (after hydrolysis)', product: t, steps: [step(m, t, 'One hydride gives the imine–aluminium complex; water hydrolyses it to the aldehyde.', ['DIBAL-H'])], why: ['Like Stephen reduction, nitrile → aldehyde'], tags: ['reduction'] }); } });

/* Stork enamine alkylation */
front({ id: 'stork', name: 'Stork enamine alkylation', fam: 'Carbonyl', canon: 'CH₃I then H₃O⁺', db: [],
  eats: s => s.any(['H3O+', 'H2O']),
  test: (c, all, k) => later(c, all, k).some(s => s.has('RX')),
  run(m, c, all, k) { const h = M(m, '[CX3;!H0]=[CX3][NX3;H0;!$(N[CX3]=O)]')[0]; if (!h) return null; const [cb, ca, n] = h; const Rs = (later(c, all, k).map(s => s.get('RX')).find(Boolean) || {}).smi || 'C';
    const t = m.clone(); for (const id of E.subtree(t, n, ca)) t.removeAtom(id); t.bondBetween(cb, ca).o = 1; E.addOn(t, ca, 'O', 2); attachR(t, cb, Rs); const R = tidy(t);
    return out({ title: 'Stork enamine alkylation → α-alkyl ketone', product: R, steps: [step(m, R, 'The enamine β-carbon (nucleophilic, like an enolate) attacks R–X; the iminium ion is hydrolysed by H₃O⁺ back to the ketone.', ['R–X', 'H₃O⁺'])], why: ['Mono-alkylation without strong base'], tags: ['enamine'] }); } });
})(typeof window !== 'undefined' ? window : globalThis);
