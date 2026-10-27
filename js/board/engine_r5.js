/* engine_r5.js: third audit batch (Oct 13). FRONT-registered. Intramolecular captures, rearrangements and
   classic named preparations that the broad recipes got wrong. */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore;
const { BL, el, P, unit, out } = E;
const REC = E.REC; const front = r => REC.unshift(r);
const M = (m, sm) => C.matchAll(m, sm, {});
const tidy = m => { const t = m.clone(); try { C.relayout(t, BL); } catch (e) { } return t; };
const fresh = smi => C.fromSmilesMB(smi, BL);
const step = (L, R, text, cond) => ({ L: L.clone(), arrows: [], R: R.clone(), text, cond: cond || [] });
const hT = (m, id) => E.hTot(m, id);
const ACID = ['H3O+', 'H2SO4', 'H3PO4', 'pTsOH', 'KHSO4'];
function distMap(m, a) { const d = new Map([[a, 0]]), q = [a]; while (q.length) { const x = q.shift(); for (const y of m.nbrs(x)) if (!d.has(y)) { d.set(y, d.get(x) + 1); q.push(y); } } return d; }
function alkenes(m) { const AR = E.aromSet(m); return M(m, '[C:1]=[C:2]').filter(([a, b]) => !AR.has(a) && !AR.has(b)); }

/* acid + alkene with a pendant OH that can close a 5/6 ring: cyclic ether */
front({ id: 'intra_ether', name: 'Intramolecular capture by –OH', fam: 'Alkene', canon: 'H⁺', db: [],
  test: c => c.any(ACID) && !c.any(['HBr', 'HCl', 'HI', 'HNO3']),
  run(m) {
    const ohs = M(m, '[CX4][OX2H1]').map(x => x[1]); if (!ohs.length) return null; const AR = E.aromSet(m);
    for (const [a, b] of alkenes(m)) {
      const sa = E.cationScore(m, a, AR), sb = E.cationScore(m, b, AR); const cp = sa >= sb ? a : b;
      for (const o of ohs) { const d = distMap(m, o).get(cp); if (d === 4 || d === 5) {
        const t = m.clone(); t.bondBetween(a, b).o = 1; t.addBond(o, cp, 1); const R = tidy(t); const n = d + 1;
        return out({ title: `Intramolecular: alkene + –OH (H⁺) → ${n}-membered cyclic ether`, product: R, steps: [step(m, R, `H⁺ adds to the C=C (Markovnikov) to give the more stable carbocation; the molecule's own –OH is close enough to capture it, closing a ${n}-membered ring; loss of H⁺.`, ['H⁺'])], why: ['Intramolecular capture beats water when a 5- or 6-ring can form', 'Markovnikov: O ends on the more substituted carbon'], tags: ['AdE'] }); } }
    }
    return null;
  } });

/* conjugated enyne + HX: X adds to the triple bond (chloroprene) */
front({ id: 'enyne_hx', name: 'Vinylacetylene + HX', fam: 'Alkyne', canon: 'HCl', db: [],
  test: c => c.any(['HCl', 'HBr']) && !c.any(['ROOR']),
  run(m, c) { const h = M(m, '[CX3]=[CX3][CX2]#[CX2;H1]')[0]; if (!h) return null; const [, , c3, c4] = h; const X = c.has('HBr') ? 'Br' : 'Cl';
    const t = m.clone(); t.bondBetween(c3, c4).o = 2; E.addOn(t, c3, X); const R = tidy(t);
    return out({ title: `Vinylacetylene + H${X} → 2-${X === 'Cl' ? 'chloro' : 'bromo'}buta-1,3-diene${X === 'Cl' ? ' (chloroprene)' : ''}`, product: R, steps: [step(m, R, `H${X} adds to the triple bond (Markovnikov), keeping the C=C conjugated: the product is a 1,3-diene.`, [`H${X}`])], why: ['Chloroprene → neoprene rubber', 'Addition keeps the diene conjugated'], tags: ['AdE'] }); } });

/* halolactonisation */
front({ id: 'halolactone', name: 'Halolactonisation', fam: 'Alkene', canon: 'I₂, NaHCO₃', db: [],
  test: c => c.any(['I2', 'Br2']) && c.any(['NaHCO3', 'NaOCl', 'KOH']),
  run(m, c) { const acid = M(m, '[CX3](=O)[OX2H1,OX1-]')[0]; if (!acid) return null; const oa = acid[2]; const X = c.has('I2') ? 'I' : 'Br'; const AR = E.aromSet(m);
    const opts = []; const da = distMap(m, oa); for (const [a, b] of alkenes(m)) for (const [near, far] of [[a, b], [b, a]]) { const d = da.get(near); if (d === 4 || d === 5) opts.push([d, a, b, near, far]); }
    opts.sort((p, q) => p[0] - q[0]);
    for (const [d, a, b, near, far] of opts) { { {
      const t = m.clone(); const oo = t.atoms.get(oa); oo.chg = 0; oo.hFix = null; t.bondBetween(a, b).o = 1; t.addBond(oa, near, 1); E.addOn(t, far, X); const R = tidy(t); const n = d + 1;
      return out({ title: `${X === 'I' ? 'Iodo' : 'Bromo'}lactonisation → ${n === 5 ? 'γ' : 'δ'}-lactone`, product: R, steps: [step(m, R, `${X}₂ forms a halonium ion on the C=C; the carboxylate (from NaHCO₃) opens it from the back, closing a ${n}-membered lactone ring with ${X} on the other carbon.`, [`${X}₂, NaHCO₃`])], why: ['Intramolecular opening of the halonium ion', '5-membered (γ) lactones form fastest'], tags: ['AdE', 'anti'] }); } } }
    return null;
  } });

/* Beckmann rearrangement */
front({ id: 'beckmann', name: 'Beckmann rearrangement', fam: 'Carbonyl', canon: 'conc. H₂SO₄', db: ['Beckmann rearrangement'],
  test: c => c.any(['H2SO4', 'PCl5', 'H3PO4', 'SOCl2', 'pTsOH', 'H3O+']),
  run(m) { const h = M(m, '[CX3]([#6])([#6])=N[OX2H1]')[0]; if (!h) return null; const [cc, r1, r2, n, o] = h;
    const AR = E.aromSet(m); const score = r => (AR.has(r) ? 10 : 0) + E.subtree(m, r, cc).length; const [mig, stay] = score(r1) >= score(r2) ? [r1, r2] : [r2, r1];
    const t = m.clone(); t.removeAtom(o); t.removeBond(t.bondBetween(cc, mig).id); t.addBond(n, mig, 1); t.bondBetween(cc, n).o = 1; E.addOn(t, cc, 'O', 2); const R = tidy(t);
    const ring = E.ringSizeOf && E.ringSizeOf(m, cc, mig);
    return out({ title: ring ? 'Beckmann rearrangement: cyclic ketoxime → lactam (ring grows by N)' : 'Beckmann rearrangement: ketoxime → N-substituted amide', product: R, steps: [step(m, R, 'Acid turns the oxime –OH into a leaving group; the group anti to it migrates from C to N as water leaves (nitrilium ion); water adds and tautomerism gives the amide.', ['H⁺'])], why: ['The group anti to –OH migrates (usually the bulkier / aryl group)', 'Cyclohexanone oxime → caprolactam → nylon-6'], tags: ['rearrangement'] }); } });

/* Knorr pyrazole: 1,3-dicarbonyl + hydrazine */
front({ id: 'pyrazole', name: '1,3-Dicarbonyl + hydrazine', fam: 'Carbonyl', canon: 'NH₂NH₂', db: [],
  test: c => c.has('NH2NH2') && !c.any(['KOH', 'aqKOH', 'glycol']),
  run(m) { const h = M(m, '[CX3]([#6])(=O)[CX4;!H0][CX3]([#6])=O')[0]; if (!h) return null; const [c1, , o1, c2, c3, , o3] = h;
    const t = m.clone(); t.removeAtom(o1); t.removeAtom(o3); const p1 = P(t, c1), p3 = P(t, c3); const d1 = E.freeDir(t, c1), d3 = E.freeDir(t, c3);
    const n1 = t.addAtom('N', p1[0] + d1[0] * BL, p1[1] + d1[1] * BL); const n2 = t.addAtom('N', p3[0] + d3[0] * BL, p3[1] + d3[1] * BL);
    t.addBond(c1, n1, 2); t.addBond(n1, n2, 1); t.addBond(n2, c3, 1); t.bondBetween(c2, c3).o = 2; const R = tidy(t);
    return out({ title: '1,3-Dicarbonyl + hydrazine → pyrazole', product: R, steps: [step(m, R, 'Hydrazine condenses with both C=O groups (two C=N formations, –2 H₂O); tautomerism makes the aromatic pyrazole ring.', ['NH₂NH₂'])], why: ['Both carbonyls are 1,3 to each other: a 5-membered ring with N–N closes', 'The ring is aromatic (6 π electrons)'], tags: ['condensation'] }); } });

/* β,γ- to α,β-unsaturated carbonyl in base */
front({ id: 'conj_isom', name: 'Isomerisation into conjugation', fam: 'Carbonyl', canon: 'OH⁻', db: [],
  test: c => c.any(['KOH', 'aqKOH', 'NaOEt', 'H3O+']),
  run(m) { const h = M(m, '[CX3;!$(C=O)]=[CX3;!$(C=O)][CX4;!H0][CX3]=O')[0]; if (!h) return null; const [a, b, cc] = h;
    const t = m.clone(); t.bondBetween(a, b).o = 1; t.bondBetween(b, cc).o = 2;
    return out({ title: 'β,γ-Unsaturated → α,β-unsaturated carbonyl (conjugation)', product: t, steps: [step(m, t, 'Base removes the α-H between C=C and C=O; the delocalised enolate takes H⁺ back at the γ-carbon, moving the C=C into conjugation with C=O.', ['OH⁻'])], why: ['Conjugated C=C–C=O is more stable'], tags: ['tautomerism'] }); } });

/* Meyer-Schuster: propargylic alcohol + acid → enone */
front({ id: 'meyer_schuster', name: 'Meyer–Schuster rearrangement', fam: 'Alkyne', canon: 'H₂SO₄', db: [],
  test: c => c.any(ACID) && !c.has('Hg2+'),
  run(m) { const h = M(m, '[CX4;!H0]([OX2H1])C#C')[0]; if (!h) return null; const [c1, o, c2, c3] = h;
    const t = m.clone(); t.removeAtom(o); t.bondBetween(c1, c2).o = 2; t.bondBetween(c2, c3).o = 1; E.addOn(t, c3, 'O', 2); const R = tidy(t);
    return out({ title: 'Meyer–Schuster: propargylic alcohol → α,β-unsaturated ketone', product: R, steps: [step(m, R, 'H⁺ helps –OH leave; water returns at the far alkyne carbon (1,3-shift) giving an allenol, which tautomerises to the conjugated C=O.', ['H⁺'])], why: ['Formal 1,3-shift of OH across the C≡C'], tags: ['rearrangement'] }); } });

/* DDT: chlorobenzene + chloral, conc. H2SO4 */
front({ id: 'ddt', name: 'DDT preparation', fam: 'Aromatic', canon: 'chloral, conc. H₂SO₄', db: [],
  test: c => c.has('chloral') && c.has('H2SO4'),
  run(m) { if (!M(m, 'Clc1ccccc1').length || m.atoms.size !== 7) return null; const R = fresh('Clc1ccc(C(c2ccc(Cl)cc2)C(Cl)(Cl)Cl)cc1');
    return out({ title: 'Chlorobenzene + chloral (conc. H₂SO₄) → DDT', product: R, steps: [step(m, R, 'Protonated chloral is the electrophile; it substitutes para on two chlorobenzene rings (loss of water).', ['chloral, conc. H₂SO₄'])], why: ['p,p′-Dichlorodiphenyltrichloroethane', 'Banned: non-biodegradable, accumulates in fat'], tags: ['EAS'] }); } });

/* phenolphthalein: phthalic anhydride + phenol */
front({ id: 'phenolphthalein', name: 'Phenolphthalein preparation', fam: 'Phenol', canon: 'phenol, conc. H₂SO₄', db: [],
  test: c => c.has('phenol') && c.has('H2SO4'),
  run(m) { if (!M(m, 'O=C1OC(=O)c2ccccc12').length) return null; const R = fresh('O=C1OC(c2ccc(O)cc2)(c2ccc(O)cc2)c2ccccc21');
    return out({ title: 'Phthalic anhydride + 2 phenol (conc. H₂SO₄) → phenolphthalein', product: R, steps: [step(m, R, 'One C=O of the anhydride, protonated, attacks two phenol rings at the para position (Friedel–Crafts type), losing water.', ['phenol, conc. H₂SO₄'])], why: ['Indicator: colourless in acid, pink in base (pH 8.3–10)'], tags: ['EAS'] }); } });

/* azo coupling with the phenol / naphthol as the drawn structure */
front({ id: 'azo_on_phenol', name: 'Azo coupling (phenol/naphthol + diazonium)', fam: 'Phenol', canon: 'C₆H₅N₂⁺Cl⁻, NaOH', db: ['Coupling reaction (azo dyes)'],
  test: c => c.has('ArN2'),
  run(m) { const oh = M(m, 'c[OX2H1,OX1-]')[0] || M(m, 'c[NX3;H2,H1]')[0]; if (!oh) return null; const [co] = oh; const AR = E.aromSet(m);
    const d = distMap(m, co); const fusion = new Set(M(m, '[c;R2]').map(x => x[0])); const arom = new Set(M(m, 'c').map(x => x[0]));
    const ringA = C.rings(m, 6).find(r => r.length === 6 && r.includes(co)) || []; const cands = M(m, '[cH]').map(x => x[0]).filter(x => ringA.includes(x) && (d.get(x) === 1 || d.get(x) === 3));
    const sc = x => (m.nbrs(x).some(y => fusion.has(y)) ? 3 : 0) + (d.get(x) === 3 ? 1 : 0); const at = cands.sort((p, q) => sc(q) - sc(p))[0]; if (at == null) return null;
    const t = m.clone(); const n1 = E.addOn(t, at, 'N'); const n2 = E.addOn(t, n1, 'N', 2); const sp = E.species('c1ccccc1'); const pn = P(t, n2), dd = E.freeDir(t, n2); const idm = t.merge(sp, pn[0] + dd[0] * BL, pn[1] + dd[1] * BL); t.addBond(n2, idm.get(1), 1);
    for (const a of t.atoms.values()) if (a.el === 'O' && a.chg === -1) a.chg = 0; const R = tidy(t);
    return out({ title: fusion.size ? 'Naphthol + diazonium (OH⁻, 0–5 °C) → orange-red azo dye' : 'Phenol + diazonium (OH⁻, 0–5 °C) → p-hydroxyazobenzene (orange)', product: R, steps: [step(m, R, 'The phenoxide ring is very electron-rich; the weak electrophile ArN₂⁺ attacks the para position (or position 1 of 2-naphthol).', ['ArN₂⁺, OH⁻, 0–5 °C'])], why: ['Mildly basic for phenols (phenoxide), mildly acidic for amines', 'Coupling goes para; if para is blocked, ortho'], tags: ['EAS'] }); } });
})(typeof window !== 'undefined' ? window : globalThis);
