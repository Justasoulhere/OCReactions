/* engine_r3.js: recipes found missing by the library coverage audit (Oct 12). Registered at the FRONT so they
   win over the broader recipes they specialise (hydroboration follow-ups, X2 in alcohol, alkoxymercuration,
   active-methylene alkylation, Gabriel, Arndt-Eistert, acyloin). */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore;
const { BL, el, P, aromSet, Mech, cationScore, freeDir, unit, out, done, deprotOn } = E;
const REC = E.REC; const front = r => REC.unshift(r);
const M = (m, sm) => C.matchAll(m, sm, {});
function alkenes(m) { const AR = aromSet(m); const s = new Set(), o = []; for (const p of M(m, '[C:1]=[C:2]')) { if (AR.has(p[0]) || AR.has(p[1])) continue; const k = [...p].sort().join(); if (!s.has(k)) { s.add(k); o.push(p); } } return o; }
const tidy = m => { const t = m.clone(); try { C.relayout(t, BL); } catch (e) { } return t; };
const fresh = smi => C.fromSmilesMB(smi, BL);
const noteStep = (L, R, text, cond) => ({ L: L.clone(), arrows: [], R: R.clone(), text, cond: cond || [] });
const later = (c, all, k) => [c, ...((all || []).slice((k || 0) + 1))];
function alcoholOf(c) {
  const r = c.get('ROH'); if (r && r.smi) return r.smi;
  const t = (c.raw || '').toLowerCase().replace(/\s+/g, '');
  if (/meoh|ch3oh|methanol/.test(t)) return 'C'; if (/etoh|c2h5oh|ethanol/.test(t)) return 'CC'; return null;
}
const ALC_NAME = { C: 'CH₃OH', CC: 'C₂H₅OH' };
const BASES = ['KOH', 'aqKOH', 'alcKOH', 'NaOEt', 'NaOR', 'NaH', 'NaNH2', 'tBuOK', 'K2CO3'];

/* 1 · hydroboration, then something other than H2O2/OH- */
front({ id: 'hb_follow', name: 'Hydroboration then protonolysis / amination / halogenation / coupling', fam: 'Alkene', canon: 'BH₃ then CH₃COOH', db: ['Hydroboration–protonolysis'],
  eats: s => s.any(['AcOH', 'NH2Cl', 'Cl2', 'Br2', 'I2', 'AgNO3', 'KOH', 'aqKOH', 'H2O']),
  test: (c, all, k) => { const cs = later(c, all, k); return c.has('BH3') && cs.some(s => s.any(['AcOH', 'NH2Cl', 'Cl2', 'Br2', 'I2', 'AgNO3'])) && !cs.some(s => s.any(['H2O2', 'H2O2/OH'])); },
  run(m, c, all, k) {
    const base = REC.find(r => r.id === 'hydroboration'); const hb = base.run(m, c, null); if (!hb || !hb.product) return null;
    const Pm = hb.product; const B = [...Pm.atoms.values()].find(a => a.el === 'B'); if (!B) return null; const cb = Pm.nbrs(B.id).find(n => Pm.atoms.get(n).el === 'C');
    const cs = later(c, all, k); const has = x => cs.some(s => s.has(x));
    let R = Pm.clone(), text, title, why;
    if (has('AcOH')) { R.removeAtom(B.id); text = 'CH₃COOH cleaves the C–B bond with retention: boron is replaced by H.'; title = 'Hydroboration–protonolysis: alkene → alkane'; why = ['Net syn addition of H₂ without a metal catalyst', 'Useful when other groups must survive hydrogenation']; }
    else if (has('NH2Cl')) { const a = R.atoms.get(B.id); a.el = 'N'; a.hFix = null; text = 'NH₂Cl replaces boron by NH₂ with retention.'; title = 'Hydroboration–amination: anti-Markovnikov 1° amine'; why = ['NH₂ ends on the less substituted carbon']; }
    else if (has('AgNO3')) {
      R.removeAtom(B.id); const R0 = R.clone(); const idm = R.merge(R0, 0, 3 * BL); R.addBond(cb, idm.get(cb), 1); R = tidy(R);
      text = 'AgNO₃ / OH⁻ couples two alkyl groups taken from boron: R–R.'; title = 'Hydroboration then AgNO₃/OH⁻: coupling to R–R'; why = ['Two anti-Markovnikov alkyl groups join'];
    } else { const X = has('Br2') ? 'Br' : has('Cl2') ? 'Cl' : 'I'; const a = R.atoms.get(B.id); a.el = X; a.hFix = null; text = `${X}₂ (with base) replaces boron by ${X}.`; title = `Hydroboration–halogenation: anti-Markovnikov alkyl ${X === 'Br' ? 'bromide' : X === 'Cl' ? 'chloride' : 'iodide'}`; why = ['X ends on the less substituted carbon']; }
    return out({ title, product: R, steps: [...hb.steps, noteStep(Pm, R, text)], why, tags: ['concerted', 'anti-Markovnikov'] });
  } });

/* 2 · X2 in an alcohol solvent: haloether */
front({ id: 'haloether', name: 'Haloether formation (X₂ in alcohol)', fam: 'Alkene', canon: 'Br₂, CH₃OH', db: ['Halohydrin formation'],
  test: c => c.any(['Br2', 'Cl2', 'I2']) && !!alcoholOf(c) && !c.has('H2O') && !c.flag('aq') && !c.any([...BASES, 'FeX3']) && !c.light,
  run(m, c) {
    const al = alkenes(m); if (!al.length) return null;
    const X = c.has('Br2') ? 'Br' : c.has('Cl2') ? 'Cl' : 'I'; const R = alcoholOf(c);
    const [a, b] = al[0]; const AR = aromSet(m); const mech = new Mech(m);
    const sa = cationScore(m, a, AR), sb = cationScore(m, b, AR); const [more, less] = sa >= sb ? [a, b] : [b, a];
    const ids = mech.add(X + X, less, { dist: 1.55 * BL }); const x1 = ids.get(1), x2 = ids.get(2);
    mech.step([[[a, b], [less, x1]], [[x1], [x1, more]], [[x1, x2], [x2]]], `The π bond attacks ${X}₂: a bridged halonium ion forms and ${X}⁻ leaves.`);
    mech.dropFrag(x2);
    const dir = unit([P(mech.m, more)[0] - P(mech.m, x1)[0], P(mech.m, more)[1] - P(mech.m, x1)[1]]);
    const w = mech.add('O' + R, more, { dist: 1.6 * BL, dir }).get(1);
    mech.step([[[w], [w, more]], [[more, x1], [x1]]], 'The alcohol is the solvent, so it, not X⁻, opens the bridge from the back at the more substituted carbon.');
    deprotOn(mech, w, 'Loss of H⁺ gives the haloether.', 'O' + R);
    return done(mech, { title: `Haloether formation (${X}₂ in ${ALC_NAME[R] || 'ROH'})`, why: ['Solvent alcohol captures the halonium ion', 'OR on the more substituted carbon, X on the other', 'Anti addition'], tags: ['AdE', 'anti'] });
  } });

/* 3 · alkoxymercuration-demercuration: ether */
front({ id: 'alkoxymerc', name: 'Alkoxymercuration–demercuration', fam: 'Alkene', canon: 'Hg(OAc)₂, C₂H₅OH then NaBH₄', db: ['Oxymercuration–demercuration'], eats: s => s.any(['NaBH4']),
  test: c => c.has('Hg2+') && !!alcoholOf(c) && !c.has('H2O') && !c.any(['H2SO4', 'H3O+']),
  run(m, c) {
    const al = alkenes(m); if (!al.length) return null; const [a, b] = al[0]; const AR = aromSet(m); const R = alcoholOf(c);
    const mech = new Mech(m); const sa = cationScore(m, a, AR), sb = cationScore(m, b, AR); const [more, less] = sa >= sb ? [a, b] : [b, a];
    const ids = mech.add({ atoms: [{ el: 'Hg', chg: 1 }, { el: 'O', x: 1 }, { el: 'C', x: 1.6, y: -0.8 }, { el: 'O', x: 1.3, y: -1.7 }, { el: 'C', x: 2.6, y: -0.6 }], bonds: [[0, 1], [1, 2], [2, 3, 2], [2, 4]] }, less, { dist: 1.5 * BL });
    const hg = ids.get(1);
    mech.step([[[a, b], [more, hg]], [[hg], [hg, less]]], 'The π bond attacks Hg(OAc)⁺: a bridged mercurinium ion (no free carbocation, so no rearrangement).');
    const w = mech.add('O' + R, more, { dist: 1.6 * BL }).get(1);
    mech.step([[[w], [w, more]], [[more, hg], [hg]]], 'The alcohol opens the bridge at the more substituted carbon.');
    deprotOn(mech, w, 'Loss of H⁺ gives the organomercury ether.', 'O' + R);
    const L = mech.m.clone(); const Rm = L.clone(); for (const id of E.subtree(Rm, hg, less)) Rm.removeAtom(id);
    mech.m = Rm; mech.note('NaBH₄ replaces HgOAc by H (demercuration).', L, Rm.clone());
    return done(mech, { title: 'Alkoxymercuration–demercuration → ether (Markovnikov)', why: ['OR on the more substituted carbon', 'No carbocation, so no rearrangement'], tags: ['AdE'] });
  } });

/* 4 · active methylene (1,3-dicarbonyl) C-alkylation */
const rxOf = (c, all, k) => later(c, all, k).map(s => s.get('RX')).find(Boolean);
front({ id: 'am_alkyl', name: 'Alkylation of an active methylene compound', fam: 'Carbonyl', canon: 'NaOEt, C₂H₅Br', db: ['Acetoacetic ester synthesis', 'Malonic ester synthesis'], eats: s => s.has('RX'),
  test: (c, all, k) => c.any(BASES) && !!rxOf(c, all, k),
  run(m, c, all, k) {
    const h = M(m, '[CX3](=O)[CX4;H1,H2][CX3](=O)')[0]; if (!h) return null; const cc = h[2];
    const rx = rxOf(c, all, k) || {}; const Rs = rx.smi || 'C', X = rx.X || 'I';
    const mech = new Mech(m);
    deprotOn(mech, cc, 'The base removes the very acidic α-H between the two C=O groups: a resonance-stabilised enolate.', '[OH-]');
    for (const g of mech.others()) mech.drop(g);
    const ids = mech.add(X + Rs, cc, { dist: 1.7 * BL, att: 2 }); const rc = ids.get(2), xa = ids.get(1);
    mech.step([[[cc], [cc, rc]], [[rc, xa], [xa]]], `The enolate carbon attacks the alkyl halide (SN2) and ${X}⁻ leaves.`);
    return done(mech, { title: 'C-alkylation of an active methylene (enolate + R–X)', why: ['α-H between two C=O is very acidic', 'The enolate alkylates on carbon; best with 1° R–X'], tags: ['enolate', 'SN2'] });
  } });

/* 5 · Gabriel phthalimide synthesis */
front({ id: 'gabriel', name: 'Gabriel phthalimide synthesis', fam: 'Amine', canon: 'KOH, C₂H₅Br then H₃O⁺', db: ['Gabriel phthalimide synthesis'],
  eats: s => s.any(['H3O+', 'H2O', 'KOH', 'aqKOH', 'NH2NH2', 'RX']),
  test: (c, all, k) => c.any(BASES) && !!rxOf(c, all, k),
  run(m, c, all, k) {
    const h = M(m, 'O=C1[NH]C(=O)c2ccccc12')[0]; if (!h) return null; const n = h[2];
    const rx = rxOf(c, all, k) || {}; const Rs = rx.smi || 'C', X = rx.X || 'Br';
    const mech = new Mech(m);
    deprotOn(mech, n, 'KOH removes the N–H proton (flanked by two C=O, so acidic): the phthalimide anion.', '[OH-]');
    for (const g of mech.others()) mech.drop(g);
    const ids = mech.add(X + Rs, n, { dist: 1.7 * BL, att: 2 }); const rc = ids.get(2), xa = ids.get(1);
    mech.step([[[n], [n, rc]], [[rc, xa], [xa]]], 'The imide anion attacks the alkyl halide (SN2): N-alkylphthalimide. It cannot alkylate twice.');
    for (const g of mech.others()) mech.drop(g);
    const cs = later(c, all, k); const hz = cs.some(s => s.has('NH2NH2')); const hyd = hz || cs.some(s => s.any(['H2O', 'H3O+', 'aqKOH', 'KOH', 'H2SO4']));
    if (!hyd) return done(mech, { title: 'N-Alkylphthalimide (Gabriel, before hydrolysis)', why: ['Hydrolyse (H₂O/OH⁻ or NH₂NH₂) to free the 1° amine'], tags: ['SN2'] });
    const Lp = mech.m.clone(); const amine = fresh('N' + Rs); const side = fresh(hz ? 'O=C1NNC(=O)c2ccccc12' : 'O=C(O)c1ccccc1C(=O)O');
    mech.steps.push(noteStep(Lp, amine, hz ? 'Hydrazine cleaves both C–N bonds: free 1° amine + phthalhydrazide.' : 'Hydrolysis breaks both amide bonds: free 1° amine + phthalic acid.'));
    return out({ title: 'Gabriel phthalimide synthesis → 1° amine', product: amine, extra: [side], steps: mech.steps, why: ['Gives pure 1° amines only', 'Fails for aryl halides: no SN2 on Ar–X'], tags: ['SN2'] });
  } });

/* 6 · Arndt-Eistert homologation */
front({ id: 'arndt', name: 'Arndt–Eistert homologation', fam: 'Acid derivative', canon: 'SOCl₂ then CH₂N₂ then Ag₂O, H₂O', db: ['Arndt–Eistert synthesis'],
  eats: s => s.any(['CH2N2', 'Ag2O', 'H2O', 'H3O+', 'SOCl2', 'AgNO3']),
  test: (c, all, k) => { const cs = later(c, all, k); return cs.some(s => s.has('CH2N2')) && cs.some(s => s.any(['Ag2O', 'Ag+', 'AgNO3'])); },
  run(m) {
    const h = M(m, '[#6][CX3](=O)[OX2H1,Cl]')[0]; if (!h) return null; const [r, cc, , lg] = h;
    const steps = []; let cur = m.clone();
    if (el(cur, lg) === 'O') { const t = cur.clone(); const a = t.atoms.get(lg); a.el = 'Cl'; a.hFix = null; steps.push(noteStep(cur, t, 'SOCl₂ turns the acid into the acid chloride.', ['SOCl₂'])); cur = t; }
    const dz = cur.clone(); const pl = P(dz, lg); dz.removeAtom(lg); const c2 = dz.addAtom('C', pl[0], pl[1]); dz.addBond(cc, c2, 1);
    const pc = P(dz, cc); const d = unit([pl[0] - pc[0], pl[1] - pc[1]]);
    const n1 = dz.addAtom('N', pl[0] + d[0] * BL, pl[1] + d[1] * BL, { chg: 1 }); dz.addBond(c2, n1, 2);
    const n2 = dz.addAtom('N', pl[0] + 2 * d[0] * BL, pl[1] + 2 * d[1] * BL, { chg: -1 }); dz.addBond(n1, n2, 2);
    steps.push(noteStep(cur, dz, 'CH₂N₂ attacks the acid chloride: a diazoketone R–CO–CH=N₂ (a second CH₂N₂ mops up the HCl).', ['CH₂N₂']));
    const kt = dz.clone(); kt.removeAtom(n1); kt.removeAtom(n2); kt.removeBond(kt.bondBetween(r, cc).id); kt.addBond(r, c2, 1); kt.bondBetween(cc, c2).o = 2;
    const ktT = tidy(kt); steps.push(noteStep(dz, ktT, 'Ag₂O catalyses the Wolff rearrangement: N₂ leaves and R shifts to the neighbouring carbon, giving a ketene.', ['Ag₂O']));
    const ac = ktT.clone(); ac.bondBetween(cc, c2).o = 1; const po = P(ac, cc); const dd = freeDir(ac, cc); const oh = ac.addAtom('O', po[0] + dd[0] * BL, po[1] + dd[1] * BL); ac.addBond(cc, oh, 1);
    const acT = tidy(ac); steps.push(noteStep(ktT, acT, 'Water adds across the ketene C=C: the acid, now one CH₂ longer.', ['H₂O']));
    return out({ title: 'Arndt–Eistert: acid → next higher acid (+CH₂)', product: acT, steps, why: ['Adds one CH₂ between R and COOH', 'Key step: Wolff rearrangement (diazoketone → ketene)'], tags: ['rearrangement'] });
  } });

/* 7 · acyloin condensation */
front({ id: 'acyloin', name: 'Acyloin condensation', fam: 'Acid derivative', canon: 'Na, xylene then H₂O', db: ['Acyloin condensation'],
  eats: s => s.any(['H2O', 'H3O+']),
  test: c => c.has('Na') && !c.any(['NH3', 'Na/NH3', 'ROH']),
  run(m) {
    const h = M(m, '[#6][CX3](=O)[OX2][#6]')[0]; if (!h) return null; const [, cc, o1, o2] = h;
    const leave = E.subtree(m, o2, cc); const acyl = m.clone(); for (const id of leave) acyl.removeAtom(id);
    const dk = acyl.clone(); const idm = dk.merge(acyl, 0, 2.2 * BL); dk.addBond(cc, idm.get(cc), 1); const dkT = tidy(dk);
    const ac = dkT.clone(); ac.bondBetween(idm.get(cc), idm.get(o1)).o = 1; const acT = tidy(ac);
    const alc = E.cleanProduct(m, leave);
    return out({ title: 'Acyloin condensation → α-hydroxy ketone', product: acT, extra: alc && alc.atoms.size ? [alc] : [],
      steps: [noteStep(m, dkT, 'Na gives one electron to each of two ester C=O groups; the radical anions couple and lose 2 RO⁻: a 1,2-diketone.', ['Na']), noteStep(dkT, acT, 'More Na reduces it to an enediolate; water work-up and tautomerism give the acyloin.', ['H₂O'])],
      why: ['Esters + Na in dry ether/xylene (no alcohol, unlike Bouveault–Blanc)', 'Product: R–CO–CH(OH)–R'], tags: ['SET'] });
  } });

/* 8 · 3° alcohol + HCl works without ZnCl2 (SN1, fast) */
front({ id: 'alc3_hcl', name: '3° alcohol + HCl', fam: 'Alcohol', canon: 'HCl', db: ['Lucas test'],
  test: c => c.has('HCl') && !c.has('ZnCl2') && !c.flag('conc'),
  run(m, c) { if (!M(m, '[CX4;H0]([#6])([#6])([#6])[OX2H1]').length) return null; const r = REC.find(x => x.id === 'alc_hx').run(m, c); if (r) { r.why = ['3° alcohols react with HCl alone at room temperature (SN1)', ...(r.why || [])]; } return r; } });

/* 9 · anionic nucleophile (alkoxide, thiolate, carboxylate, phenoxide) + R–X */
front({ id: 'nu_rx', name: 'Nucleophile + alkyl halide (SN2)', fam: 'Alkyl halide', canon: 'CH₃I', db: ['Williamson ether synthesis', 'SN2 reaction'],
  test: c => c.has('RX'),
  run(m, c) {
    const h = M(m, '[O-,S-;X1;!$([O-][#7+]);!$([O-][#16])]')[0]; if (!h) return null; const nu = h[0];
    const rx = c.get('RX') || {}; const Rs = rx.smi || 'C', X = rx.X || 'I';
    const mech = new Mech(m);
    const ids = mech.add(X + Rs, nu, { dist: 1.7 * BL, att: 2 }); const rc = ids.get(2), xa = ids.get(1);
    const carb = M(m, '[O-][CX3]=O').length > 0, S_ = el(m, nu) === 'S';
    mech.step([[[nu], [nu, rc]], [[rc, xa], [xa]]], `The ${S_ ? 'thiolate' : carb ? 'carboxylate' : 'alkoxide'} attacks the carbon of R–X from the back (SN2) and ${X}⁻ leaves.`);
    return done(mech, { title: S_ ? 'SN2: thiolate + R–X → thioether' : carb ? 'SN2: carboxylate + R–X → ester' : 'Williamson synthesis: alkoxide + R–X → ether', why: ['SN2: best with CH₃ and 1° halides', '3° halides eliminate instead'], tags: ['SN2'] });
  } });

/* 10 · cyclopropane ring opening (X2, HX, H2/Ni) */
front({ id: 'cyclopropane', name: 'Cyclopropane ring opening', fam: 'Alkane', canon: 'Br₂', db: [],
  test: c => c.any(['Br2', 'Cl2', 'HBr', 'HI', 'HCl', 'H2']) && !c.light,
  run(m, c) {
    const h = M(m, '[CX4;r3]1[CX4;r3][CX4;r3]1')[0]; if (!h) return null;
    const hc = id => E.hTot(m, id); const pairs = [[h[0], h[1]], [h[1], h[2]], [h[0], h[2]]].sort((p, q) => (hc(q[0]) + hc(q[1])) - (hc(p[0]) + hc(p[1])));
    const isHX = !c.has('H2') && !c.any(['Br2', 'Cl2']);
    let [a, b] = pairs[0];
    if (isHX) { const sub = [...h].sort((p, q) => hc(p) - hc(q))[0]; const part = h.filter(x => x !== sub).sort((p, q) => hc(q) - hc(p))[0]; a = part; b = sub; }
    const t = m.clone(); t.removeBond(t.bondBetween(a, b).id);
    let title, text;
    if (c.has('H2')) { title = 'Hydrogenolysis of cyclopropane → alkane'; text = 'Ring strain (60° angles) lets H₂/Ni open the ring.'; }
    else if (!isHX) { const X = c.has('Br2') ? 'Br' : 'Cl'; E.addOn(t, a, X); E.addOn(t, b, X); title = `Cyclopropane + ${X}₂ → 1,3-di${X === 'Br' ? 'bromo' : 'chloro'}alkane`; text = `The strained ring behaves like a weak C=C: ${X}₂ opens it and adds to both ends (1,3-addition).`; }
    else { const X = c.has('HBr') ? 'Br' : c.has('HI') ? 'I' : 'Cl'; E.addOn(t, b, X); title = `Cyclopropane + H${X}: ring opening (Markovnikov-like)`; text = `H goes to the carbon with more H; ${X} to the more substituted end.`; }
    const R = tidy(t);
    return out({ title, product: R, steps: [noteStep(m, R, text)], why: ['Angle strain (Baeyer): cyclopropane reacts like an alkene', 'Cyclobutane and larger rings do not open this easily'], tags: ['ring-opening'] });
  } });

/* 11 · Frankland (Zn) and Ullmann (Cu) couplings via the Wurtz recipe */
front({ id: 'zn_cu_coupling', name: 'Frankland / Ullmann coupling', fam: 'Alkyl halide', canon: 'Cu, Δ', db: ['Ullmann reaction', 'Frankland reaction'],
  test: c => c.any(['Zn', 'Cu']) && !c.any(['HCl', 'Sn/HCl', 'H2SO4', 'H3O+', 'CuCN', 'CuCl', 'CuBr']),
  run(m, c) {
    const aryl = M(m, 'c[Cl,Br,I]').length > 0; const alkyl = M(m, '[CX4][Cl,Br,I]').length > 0;
    if (c.has('Cu') && !aryl) return null; if (c.has('Zn') && !alkyl) return null;
    const r = REC.find(x => x.id === 'wurtz').run(m); if (!r) return null;
    r.title = c.has('Cu') ? 'Ullmann coupling: 2 Ar–I + Cu → Ar–Ar' : 'Frankland coupling: 2 R–X + Zn → R–R';
    if (r.steps[0]) r.steps[0].text = c.has('Cu') ? 'Copper (heat) joins two aryl groups: biaryl + CuI.' : 'Zinc couples two alkyl groups, like Na in the Wurtz reaction.';
    return r;
  } });

/* 12 · acid-catalysed α-halogenation of ketones (one X, via the enol) */
front({ id: 'alpha_x', name: 'α-Halogenation (acid)', fam: 'Carbonyl', canon: 'Br₂, CH₃COOH', db: [],
  test: c => c.any(['Br2', 'Cl2', 'I2']) && c.any(['AcOH', 'H3O+']) && !c.any(BASES),
  run(m, c) {
    const hs = M(m, '[CX3;$(C([#6])([#6])=O),$([CH1]([#6])=O)](=O)[CX4;!H0]'); if (!hs.length) return null;
    const best = hs.sort((p, q) => E.hTot(m, p[2]) - E.hTot(m, q[2]))[0]; const [cc, o, ca] = best;
    const X = c.has('Br2') ? 'Br' : c.has('Cl2') ? 'Cl' : 'I';
    const en = m.clone(); en.bondBetween(cc, o).o = 1; en.bondBetween(cc, ca).o = 2;
    const pr = m.clone(); E.addOn(pr, ca, X);
    return out({ title: `Acid-catalysed α-${X === 'Br' ? 'bromination' : X === 'Cl' ? 'chlorination' : 'iodination'} (one ${X} only)`, product: pr,
      steps: [noteStep(m, en, 'H⁺ catalyses enolisation; the more substituted enol forms (rate-determining).', ['H⁺']), noteStep(en, pr, `The enol C=C attacks ${X}₂; loss of H⁺ gives the α-halo ketone. The new ${X} slows a second enolisation, so it stops at one.`, [`${X}₂`])],
      why: ['Rate is independent of [X₂]: enol formation is slow', 'In base, every α-H is replaced (haloform for methyl ketones)'], tags: ['enol'] });
  } });

/* 13 · 3° alkyl ether + dilute acid: SN1 cleavage */
front({ id: 'ether3_hyd', name: 'Acid hydrolysis of a 3° alkyl ether', fam: 'Ether', canon: 'dil. H₂SO₄, H₂O', db: ['Acid hydrolysis of ethers'],
  test: c => c.any(['H3O+', 'H2SO4']) && (c.has('H2O') || c.flag('dil') || c.flag('aq')) && !c.any(['HI', 'HBr', 'HCl']),
  run(m) {
    const h = M(m, '[CX4;H0]([#6])([#6])([#6])[OX2][CX4]')[0]; if (!h) return null; const cq = h[0], o = h[4];
    const t = m.clone(); t.removeBond(t.bondBetween(cq, o).id); E.addOn(t, cq, 'O');
    const fr = t.frags(); const fa = fr.find(f => f.includes(cq)), fb = fr.find(f => f.includes(o));
    const A1 = tidy(E.cleanProduct(t, fa)), B1 = fb ? E.cleanProduct(t, fb) : null;
    return out({ title: '3° ether + dilute acid → 3° alcohol + alcohol (SN1)', product: A1, extra: B1 ? [B1] : [], steps: [noteStep(m, A1, 'H⁺ protonates the ether O; the C–O bond to the 3° carbon breaks (SN1) because that gives the stable 3° carbocation; water captures it.', ['H₃O⁺'])], why: ['3° alkyl ethers cleave easily (stable carbocation)', '1°/2° ethers need HI or HBr'], tags: ['SN1'] });
  } });

/* 14 · isomerisation of n-alkanes (anhydrous AlCl3 / HCl) */
front({ id: 'alkane_isom', name: 'Isomerisation of alkanes', fam: 'Alkane', canon: 'anhyd. AlCl₃, HCl, Δ', db: [],
  test: c => c.has('FeX3') && !c.any(['Br2', 'Cl2', 'RX', 'RCOCl', 'CH3COCl']) && /alcl3/i.test((c.raw || '').replace(/[₀-₉]/g, d => '0123456789'['₀₁₂₃₄₅₆₇₈₉'.indexOf(d)])),
  run(m) {
    const at = [...m.atoms.values()]; if (at.some(a => a.el !== 'C') || at.length < 4 || [...m.bonds.values()].some(b => b.o !== 1) || E.aromSet(m).size || C.rings(m).length) return null;
    const ends = at.filter(a => m.degree(a.id) === 1); if (ends.length !== 2) return null;
    const path = [ends[0].id]; while (path.length < at.length) { const nx = m.nbrs(path[path.length - 1]).find(n => !path.includes(n)); if (nx == null) break; path.push(nx); }
    const t = m.clone(); t.removeBond(t.bondBetween(path[0], path[1]).id); t.addBond(path[0], path[2], 1); const R = tidy(t);
    return out({ title: 'Isomerisation: n-alkane → branched isomer', product: R, steps: [noteStep(m, R, 'Anhydrous AlCl₃ / HCl (heat) generates carbocations; methyl shifts give the more branched, more stable isomer.', ['AlCl₃, HCl'])], why: ['Branched alkanes are more stable', 'Used to raise the octane number'], tags: ['rearrangement'] });
  } });

/* 15 · terminal alkyne + ammoniacal AgNO3 / Cu2Cl2: metal acetylide (test for ≡C–H) */
front({ id: 'acetylide_test', name: 'Terminal alkyne test (Ag⁺ / Cu⁺ acetylide)', fam: 'Alkyne', canon: 'ammoniacal AgNO₃', db: ["Tollens' test"],
  test: c => c.any(['Tollens', 'CuNH3']) || (c.any(['AgNO3', 'CuCl', 'Ag+']) && c.has('NH3')),
  run(m, c) {
    const h = M(m, '[CH1]#[#6]')[0]; if (!h) return null; const Cu = c.any(['CuNH3', 'CuCl']);
    const t = m.clone(); for (const [c0] of M(m, '[CH1;$(C#[#6])]')) { const c1 = m.nbrs(c0).find(n => m.bondBetween(c0, n).o === 3); const p0 = P(t, c0), p1 = P(t, c1); const d = unit([p0[0] - p1[0], p0[1] - p1[1]]); const mid = t.addAtom(Cu ? 'Cu' : 'Ag', p0[0] + d[0] * BL, p0[1] + d[1] * BL); t.addBond(c0, mid, 1); t.atoms.get(c0).hFix = null; }
    return out({ title: Cu ? 'Terminal alkyne + ammoniacal Cu₂Cl₂ → copper(I) acetylide (red ppt)' : 'Terminal alkyne + ammoniacal AgNO₃ → silver acetylide (white ppt)', product: t, steps: [noteStep(m, t, `The acidic ≡C–H (sp carbon, 50% s character) is replaced by ${Cu ? 'Cu⁺' : 'Ag⁺'}.`)], why: ['Positive only for terminal alkynes (≡C–H)', 'Distinguishes but-1-yne from but-2-yne'], tags: ['test'] });
  } });

/* 16 · no Cannizzaro for α,β-unsaturated aldehydes that have γ-H (they enolise) */
front({ id: 'cann_guard', name: 'Cannizzaro guard', fam: 'Carbonyl', canon: '', db: [],
  test: c => c.any(['KOH', 'aqKOH']) && !c.flag('dil'),
  run(m) { if (!M(m, '[CH1](=O)[#6]=[#6][CX4;!H0]').length) return null; return { none: 'This aldehyde has γ-hydrogens through its C=C (vinylogous enolisation), so it condenses with itself in base instead of doing Cannizzaro.', title: 'No Cannizzaro' }; } });
})(typeof window !== 'undefined' ? window : globalThis);
