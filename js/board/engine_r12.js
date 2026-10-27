/* engine_r12.js: general recipes for reactions the engine could only answer from the library (run 28).
   Written as RDKit reaction templates; appended to the END of the recipe list so they never pre-empt older recipes. */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore; if (!C || !E) return;
const { BL } = E;
const fresh = smi => C.fromSmilesMB(smi, BL);
const step = (L, R, text, cond) => ({ L: L.clone(), arrows: [], R: R.clone(), text, cond: cond || [] });
const out = o => Object.assign({ steps: [], minor: [], extra: [], why: [], db: [] }, o);
const M = (m, sm) => { try { return C.matchAll(m, sm, {}); } catch (e) { return []; } };
const RXC = new Map();
function run(sma, smis) {
  const RD = C.RD; let rx = RXC.get(sma); if (!rx) { try { rx = RD.get_rxn(sma); } catch (e) { rx = null; } RXC.set(sma, rx); } if (!rx) return [];
  const L = new RD.MolList(); const ms = smis.map(s => RD.get_mol(s)); if (ms.some(x => !x)) return []; ms.forEach(x => L.append(x));
  const res = []; let o = null; try { o = rx.run_reactants(L, 24); } catch (e) { o = null; }
  if (o) for (let i = 0; i < o.size(); i++) { const ml = o.get(i); const ps = []; let bad = false;
    for (let k = 0; k < ml.size(); k++) { const x = ml.at(k); let s = null; try { const q = RD.get_mol(x.get_smiles()); if (q) { s = q.get_smiles(); q.delete(); } } catch (e) { s = null; } x.delete && x.delete(); if (!s) bad = true; ps.push(s); }
    if (!bad && !res.some(r => r.join('.') === ps.join('.'))) res.push(ps); }
  ms.forEach(x => x.delete()); return res;
}
const smiOf = m => { try { return C.canon(m); } catch (e) { return null; } };
const nC = s => (String(s).replace(/Cl|Ca|Cu|Cr/g, '').match(/[Cc]/g) || []).length;
const add = r => E.REC.push(r);
const finish = (m, prodSmi, o) => { const R = fresh(prodSmi); if (!R || !R.atoms.size) return null; return out(Object.assign({ product: R, steps: [step(m, R, o.text, o.cond)] }, o, { extra: (o.extra || []).map(fresh).filter(Boolean), minor: (o.minor || []).map(fresh).filter(Boolean) })); };

add({ id: 'acetal_roh', name: 'Acetal / ketal formation', fam: 'Carbonyl', canon: 'ROH (2 eq), dry HCl', db: ['Acetal formation'],
  test: c => c.any(['H3O+', 'HCl', 'H2SO4', 'TsOH']) && c.flag('alc') && c.flag('dry'),
  run(m, c) { const s = smiOf(m); if (!s || M(m, '[CX3](=O)[O,N,Cl,Br]').length) return null;
    const R = /meoh|methanol|ch3oh/i.test(c.raw) ? 'C' : 'CC';
    const p = run(`[CX3;!$(C-[O,N,S,Cl,Br]):1]=[OX1:2]>>[C:1](O${R})O${R}`, [s]); if (!p.length) return null;
    const ald = M(m, '[CX3H1]=O').length > 0;
    return finish(m, p[0][0], { title: `${ald ? 'Aldehyde' : 'Ketone'} + ${R === 'C' ? 'MeOH' : 'EtOH'} (dry H⁺) → ${ald ? 'acetal' : 'ketal'}`, text: 'Dry acid protonates C=O; one alcohol adds (hemiacetal), the –OH is protonated and leaves as water, and the oxocarbenium ion takes a second alcohol. Water must be kept out (it reverses it).', cond: ['ROH', 'dry H⁺'], extra: ['O'], why: ['Acetals survive base, Grignard and hydride reagents: used to protect C=O', 'Aqueous acid gives the carbonyl back'], tags: ['AdN'] }); } });

add({ id: 'hydroxy_acid_heat', name: 'Hydroxy acids on heating', fam: 'Acid', canon: 'Δ', db: ['Effect of heat on hydroxy acids'],
  test: c => c.heat,
  run(m) { const s = smiOf(m); if (!s) return null;
    const a = run('[OX2H1:1][CX4:2][CX3:3](=[O:4])[OX2H1:5].[OX2H1:6][CX4:7][CX3:8](=[O:9])[OX2H1:10]>>[C:2]1[C:3](=[O:4])[O:6][C:7][C:8](=[O:9])[O:1]1', [s, s]);
    if (a.length) return finish(m, a[0][0], { title: 'α-Hydroxy acid + Δ → lactide (cyclic diester)', text: 'Two molecules esterify each other (each –OH with the other –COOH): a six-membered cyclic diester, the lactide. Two water molecules leave.', cond: ['Δ'], extra: ['O'], why: ['α: lactide · β: α,β-unsaturated acid · γ, δ: lactone'], tags: ['ester'] });
    const b = run('[OX2H1:1][CX4:2][CX4;!H0:3][CX3:4](=[O:5])[OX2H1:6]>>[C:2]=[C:3][C:4](=[O:5])[O:6]', [s]);
    if (b.length) return finish(m, b[0][0], { title: 'β-Hydroxy acid + Δ → α,β-unsaturated acid', text: 'The α-H is acidic (next to C=O) and the β-OH leaves: dehydration gives the conjugated α,β-unsaturated acid.', cond: ['Δ'], extra: ['O'], why: ['Conjugation with C=O drives the loss of water', 'α: lactide · γ, δ: lactone'], tags: ['elim'] });
    const g = run('[OX2H1:1][CX4:2][CX4:3][CX4:4][CX3:5](=[O:6])[OX2H1:7]>>[O:1]1[C:2][C:3][C:4][C:5]1=[O:6]', [s]);
    if (g.length) return finish(m, g[0][0], { title: 'γ-Hydroxy acid + Δ → γ-lactone (5-ring)', text: 'The –OH and –COOH of the same molecule esterify: a strain-free five-membered cyclic ester (γ-butyrolactone type). It forms so easily that the free γ-hydroxy acid is hard to keep.', cond: ['Δ'], extra: ['O'], why: ['γ → 5-ring lactone, δ → 6-ring lactone'], tags: ['ester'] });
    const d = run('[OX2H1:1][CX4:2][CX4:3][CX4:4][CX4:8][CX3:5](=[O:6])[OX2H1:7]>>[O:1]1[C:2][C:3][C:4][C:8][C:5]1=[O:6]', [s]);
    if (d.length) return finish(m, d[0][0], { title: 'δ-Hydroxy acid + Δ → δ-lactone (6-ring)', text: 'Intramolecular esterification: a six-membered lactone.', cond: ['Δ'], extra: ['O'], why: ['γ → 5-ring, δ → 6-ring lactone'], tags: ['ester'] });
    return null; } });

add({ id: 'ca_salt', name: 'Dry distillation of calcium carboxylates', fam: 'Acid', canon: 'dry distillation (Δ)', db: ['Dry distillation of calcium salts'],
  test: c => c.heat,
  run(m) { const s = smiOf(m); if (!s || !/Ca/.test(s)) return null;
    const p = run('[#6:6][C:1](=[O:2])[O][Ca][O][C:3](=[O:4])[#6:5]>>[#6:6][C:1](=[O:2])[#6:5]', [s]);
    const f = run('[C:1](=[O:2])([#6:6])[O][Ca][O][CH1:3]=[O:4]>>[#6:6][CH1:1]=[O:2]', [s]);
    const all = [...new Set(p.concat(f).map(x => x[0]))]; if (!all.length) return null;
    return finish(m, all[0], { title: 'Calcium carboxylate, dry distillation → ketone (+ CaCO₃)', text: 'Heating the dry calcium salt: one carboxylate passes its alkyl group to the other carbonyl carbon and CaCO₃ is left behind. (RCOO)₂Ca → R–CO–R; with calcium formate the product is an aldehyde.', cond: ['dry distillation', 'Δ'], extra: ['O=C([O-])[O-]'], minor: all.slice(1), why: ['Mixed salts give a mixture of three carbonyls', 'Calcium formate + another calcium salt → aldehyde'], tags: ['decarbox'] }); } });

E.REC.unshift({ id: 'hunsdiecker', name: 'Hunsdiecker reaction', fam: 'Acid', canon: 'Br₂, CCl₄, Δ', db: ['Hunsdiecker reaction'],
  test: c => c.any(['Br2', 'Cl2']),
  run(m, c) { const s = smiOf(m); if (!s || !/Ag/.test(s)) return null; const X = c.has('Br2') ? 'Br' : 'Cl';
    const p = run(`[#6:1][CX3](=O)[O][Ag]>>[#6:1]${X}`, [s]); if (!p.length) return null;
    return finish(m, p[0][0], { title: `Hunsdiecker: RCOOAg + ${X}₂ → R–${X} (one C fewer)`, text: `The silver salt reacts with ${X}₂ (dry CCl₄, heat): an acyl hypohalite RCOO–${X} forms, its O–${X} bond breaks homolytically, the RCOO• radical loses CO₂, and R• takes ${X} from the next molecule (radical chain).`, cond: [X + '₂', 'CCl₄', 'Δ'], extra: ['O=C=O', '[Ag+]'], why: ['Product has one carbon fewer than the acid', 'Order of yield: 1° > 2° > 3°', 'With I₂ (2 : 1) the ester forms instead (Simonini)'], tags: ['radical'] }); } });

add({ id: 'pinacol_mg', name: 'Pinacol formation (Mg–Hg)', fam: 'Carbonyl', canon: 'Mg–Hg, then H₂O', db: ['Pinacol formation'],
  test: c => c.has('MgHg'),
  run(m, c) { const s = smiOf(m); if (!s) return null;
    const T = '[#6:1][CX3:2](=[O:3])[#6,#1:4].[#6:5][CX3:6](=[O:7])[#6,#1:8]>>[#6:1][C:2]([OH:3])([*:4])[C:6]([OH:7])([*:5])[*:8]';
    const self = run(T, [s, s]); if (!self.length) return null; let main = self[0][0], minor = [];
    if (c.has('acetone') && s !== 'CC(C)=O') { const x = run(T, [s, 'CC(C)=O']); if (x.length) { minor = [main]; main = x[0][0]; } }
    return finish(m, main, { title: 'Ketone + Mg–Hg, then H₂O → pinacol (1,2-diol)', text: 'Mg gives one electron to each C=O (single-electron transfer): two ketyl radical anions pair up at carbon, giving a magnesium 1,2-dialkoxide. Water then protonates both O: the pinacol.', cond: ['Mg–Hg', 'H₂O'], minor, why: ['Reductive coupling: two carbonyl carbons become one C–C bond', 'Pinacol + H⁺ → pinacolone (pinacol–pinacolone rearrangement)'], tags: ['radical'] }); } });

add({ id: 'ester_pyrolysis', name: 'Ester pyrolysis (syn elimination)', fam: 'Ester', canon: 'Δ (~500 °C)', db: ['Pyrolysis of esters'],
  test: c => c.heat,
  run(m, c) { const s = smiOf(m); if (!s) return null;
    const p = run('[#6,#1:1][CX3:2](=[O:3])[OX2:4][CX4:5][CX4;!H0:6]>>[*:1][C:2](=[O:3])[OH:4].[C:5]=[C:6]', [s]); if (!p.length) return null;
    const alk = p.map(x => x[1]); const sub = q => { const mm = C.fromSmilesMB(q, BL); return M(mm, '[CX3]=[CX3]').map(([a, b]) => mm.nbrs(a).concat(mm.nbrs(b)).filter(n => mm.atoms.get(n).el === 'C').length)[0] || 0; };
    const ord = [...new Set(alk)].sort((x, y) => sub(x) - sub(y));
    return finish(m, ord[0], { title: 'Ester pyrolysis → alkene + carboxylic acid', text: 'Strong heating: a six-membered cyclic transition state in which the C=O oxygen pulls off a β-H on the same side while the C–O bond breaks (syn, Ei elimination). No base, no carbocation, no rearrangement.', cond: ['Δ'], extra: [p[0][0]], minor: ord.slice(1), why: ['Syn elimination; the less substituted alkene usually predominates (more β-H on CH₃ groups)'], tags: ['elim'] }); } });

add({ id: 'acyl_azide_heat', name: 'Acyl azide → isocyanate (Curtius)', fam: 'Acid derivative', canon: 'Δ', db: ['Curtius rearrangement'],
  test: c => c.heat && !c.any(['H2O', 'H3O+', 'KOH', 'NaOH', 'aqKOH']),
  run(m) { const s = smiOf(m); if (!s) return null;
    const p = run('[#6:1][CX3:2](=[O:3])N=[N+]=[N-]>>[#6:1]N=[C:2]=[O:3]', [s]); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Acyl azide + Δ → isocyanate (Curtius rearrangement)', text: 'N₂ leaves and, at the same moment, R moves from the carbonyl carbon to nitrogen (with retention of configuration): R–N=C=O. Water would then give R–NH₂ + CO₂.', cond: ['Δ'], extra: ['N#N'], why: ['Isocyanate + H₂O → carbamic acid → amine + CO₂', 'Isocyanate + ROH → carbamate (urethane)'], tags: ['rearr'] }); } });

add({ id: 'tosylation', name: 'Tosylation (TsCl, pyridine)', fam: 'Alcohol', canon: 'TsCl, pyridine', db: ['Tosylates'],
  test: c => c.has('TsCl'),
  run(m) { const s = smiOf(m); if (!s) return null; const Ts = 'S(=O)(=O)c1ccc(C)cc1';
    let p = run(`[CX4:1][OX2H1:2]>>[C:1][O:2]${Ts}`, [s]); let kind = 'alcohol';
    if (!p.length) { p = run(`[NX3;H2,H1;!$(NC=O);!$(NS=O):1]>>[N:1]${Ts}`, [s]); kind = 'amine'; }
    if (!p.length) return null;
    return finish(m, p[0][0], kind === 'alcohol'
      ? { title: 'R–OH + TsCl (pyridine) → R–OTs (tosylate)', text: 'The alcohol O attacks sulfur of TsCl and Cl⁻ leaves; pyridine takes the H⁺. The C–O bond is not touched, so the configuration at carbon is retained. –OTs is an excellent leaving group (as good as a halide).', cond: ['TsCl', 'pyridine'], extra: ['Cl'], why: ['Converts a poor leaving group (–OH) into a very good one (–OTs) without a carbocation', 'Then SN2 with Nu⁻ gives inversion', 'A 3° tosylate ionises easily: TsO⁻ leaves and the 3° carbocation forms'], tags: ['sub'] }
      : { title: 'Amine + TsCl → sulfonamide', text: 'The amine N attacks sulfur and Cl⁻ leaves: an N-tosyl sulfonamide (Hinsberg-type).', cond: ['TsCl'], extra: ['Cl'], why: ['1° sulfonamide has an acidic N–H (dissolves in alkali); 2° does not; 3° amines do not react'], tags: ['sub'] }); } });


/* ================= batch 2 (run 29): gaps found by t_general ================= */
E.REC.unshift({ id: 'victor_meyer', name: 'Victor Meyer test (1°/2°/3° alcohol)', fam: 'Alcohol', canon: 'P/I₂, AgNO₂, HNO₂, NaOH', db: ['Victor Meyer test'],
  test: c => (c.has('redP') || c.has('PI3') || c.has('I2')) && c.has('AgNO2') && c.has('HNO2'),
  run(m) { const s = smiOf(m); if (!s) return null;
    const p1 = run('[CX4H2:1]([#6:2])[OX2H1]>>[*:2]/[C:1](=N/[O-])[N+](=O)[O-]', [s]);
    if (p1.length) return finish(m, p1[0][0], { title: 'Victor Meyer: 1° alcohol → nitrolic acid → blood-red salt', text: 'P/I₂ gives R–CH₂–I, AgNO₂ gives the nitroalkane R–CH₂–NO₂, HNO₂ nitrosates the CH₂ (two α-H): the nitrolic acid R–C(=NOH)–NO₂. NaOH turns it into its blood-red sodium salt.', cond: ['P/I₂', 'AgNO₂', 'HNO₂', 'NaOH'], why: ['1°: red · 2°: blue (pseudonitrole) · 3°: colourless', 'The colour depends on how many α-H the nitroalkane has'], tags: ['test'] });
    const p2 = run('[CX4H1:1]([#6:2])([#6:3])[OX2H1]>>[*:2][C:1]([*:3])(N=O)[N+](=O)[O-]', [s]);
    if (p2.length) return finish(m, p2[0][0], { title: 'Victor Meyer: 2° alcohol → pseudonitrole (blue)', text: 'R₂CH–I → R₂CH–NO₂; HNO₂ nitrosates the single α-H: R₂C(NO)(NO₂), a pseudonitrole. It gives a blue colour (it has no α-H left, so NaOH makes no salt).', cond: ['P/I₂', 'AgNO₂', 'HNO₂', 'NaOH'], why: ['1°: red · 2°: blue · 3°: colourless'], tags: ['test'] });
    if (M(m, '[CX4H0]([#6])([#6])([#6])[OX2H1]').length) return out({ title: 'Victor Meyer: 3° alcohol → colourless', product: null, none: 'The 3° nitroalkane R₃C–NO₂ has no α-H, so HNO₂ cannot nitrosate it: no colour (colourless).', why: ['1°: red · 2°: blue · 3°: colourless'], tags: ['test'] });
    return null; } });

add({ id: 'acid_aldol', name: 'Acid-catalysed aldol condensation', fam: 'Carbonyl', canon: 'H⁺ (dry HCl / conc. H₂SO₄)', db: ['Aldol condensation'],
  test: c => c.any(['H3O+', 'HCl', 'H2SO4']) && !c.flag('alc') && !c.flag('cat') && (c.flag('conc') || c.flag('dry') || c.heat || c.any(['acetone', 'CH3CHO', 'HCHO'])) && !c.any(['NaOH', 'KOH', 'aqKOH', 'NaBH4', 'LiAlH4', 'Zn', 'ZnHg', 'HCN', 'NaCN']),
  run(m, c) { const s = smiOf(m); if (!s || m.frags().length !== 1 || M(m, '[CX3](=O)[O,N,Cl,Br]').length) return null;
    if (!M(m, '[CX3](=O)[CX4;!H0]').length) return null;
    if (c.has('H2SO4') && c.flag('conc')) { const t = run('[CH3:1][CX3:2](=O)[#6:3].[CH3:4][CX3:5](=O)[#6:6].[CH3:7][CX3:8](=O)[#6:9]>>[c:2]1([*:3])[cH:4][c:5]([*:6])[cH:7][c:8]([*:9])[cH:1]1', [s, s, s]);
      if (t.length) return finish(m, t[0][0], { title: 'Methyl ketone + conc. H₂SO₄ → 1,3,5-trisubstituted benzene', text: 'Three molecules undergo successive acid-catalysed aldol condensations (enol attacks protonated C=O, water leaves each time); the last one closes a six-membered ring, which loses water to become aromatic. Acetone gives mesitylene.', cond: ['conc. H₂SO₄', 'Δ'], extra: ['O'], why: ['3 CH₃COCH₃ → 1,3,5-trimethylbenzene + 3 H₂O'], tags: ['aldol'] }); }
    const one = run('[CX3:1](=[O:2])[CX4;!H0:3].[CX3;!$(C-[O,N,Cl]):4]=[O:5]>>[C:1](=[O:2])[C:3]=[C:4]', [s, s]); if (!one.length) return null;
    const pick = one.map(x => x[0]); const main1 = pick.sort((a, b) => a.length - b.length)[0];
    if (c.has('HCl') && c.flag('dry')) { const two = run('[CX3:1](=[O:2])[CX4;!H0:3].[CX3;!$(C-[O,N,Cl]):4]=[O:5]>>[C:1](=[O:2])[C:3]=[C:4]', [main1, s]);
      const ok2 = two.map(x => x[0]).filter(x => x !== main1 && nC(x) === nC(main1) + nC(s) && !/C=CC=C/.test(x.replace(/[()]/g, '')) );
      if (ok2.length) return finish(m, ok2[0], { title: 'Ketone + dry HCl → double aldol condensation (phorone type)', text: 'The enol of one molecule attacks a protonated C=O, water leaves (α,β-unsaturated ketone: mesityl oxide from acetone); a third molecule condenses on the other α-carbon: phorone.', cond: ['dry HCl'], extra: ['O'], minor: [main1], why: ['Acetone: mesityl oxide (1 condensation) and phorone (2)'], tags: ['aldol'] }); }
    return finish(m, main1, { title: 'Acid-catalysed aldol condensation → α,β-unsaturated carbonyl', text: 'Acid makes the enol; the enol carbon attacks a protonated C=O of a second molecule; under acid the β-hydroxy carbonyl loses water at once (E1, conjugated product).', cond: ['H⁺'], extra: ['O'], why: ['Under acid the dehydrated product forms directly', 'CH₃CHO → crotonaldehyde · acetone → mesityl oxide'], tags: ['aldol'] }); } });

add({ id: 'alkane_sulfonation', name: 'Sulfonation of alkanes (3° C–H)', fam: 'Alkane', canon: 'fuming H₂SO₄, Δ', db: ['Sulphonation of alkanes'],
  test: c => c.has('H2SO4') && c.heat,
  run(m) { const s = smiOf(m); if (!s || !/^[C()]+$/.test(s)) return null;
    const p = run('[CX4H1:1]([CH3:2])([#6:3])[#6:4]>>[C:1]([*:2])([*:3])([*:4])S(=O)(=O)O', [s]);
    if (!p.length) return out({ title: 'No sulfonation', product: null, none: 'Only alkanes with a 3° C–H (and long-chain alkanes, slowly) are sulfonated; this one has no 3° hydrogen.', why: [], tags: [] });
    return finish(m, p[0][0], { title: 'Alkane + fuming H₂SO₄ (Δ) → alkanesulfonic acid', text: 'Fuming sulfuric acid (SO₃) replaces the 3° hydrogen by –SO₃H (free-radical-like substitution at the weakest C–H).', cond: ['H₂SO₄ (fuming)', 'Δ'], extra: ['O'], why: ['Ease: 3° > 2° > 1°; lower alkanes with no 3° H do not react'], tags: ['sub'] }); } });

add({ id: 'orthoester_hyd', name: 'Orthoester hydrolysis', fam: 'Ether', canon: 'H₂O, H⁺', db: [],
  test: c => c.any(['H3O+', 'H2O']),
  run(m) { const s = smiOf(m); if (!s) return null; const p = run('[CX4:1]([OX2:2][#6:5])([OX2][#6])([OX2][#6])>>[C:1](=O)[O:2][*:5]', [s]); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Orthoester + H₂O (H⁺) → ester + 2 alcohol', text: 'Like an acetal: H⁺ on one –OR, the alcohol leaves, the dialkoxy-carbocation takes water; repeat once more and the ester remains (more water/heat would hydrolyse it to the acid).', cond: ['H₂O', 'H⁺'], extra: ['CO'], why: ['Orthoesters are stable to base, hydrolysed by aqueous acid'], tags: ['hyd'] }); } });

add({ id: 'simonini', name: 'Simonini (Birnbaum–Simonini) reaction', fam: 'Acid', canon: 'I₂ (RCOOAg : I₂ = 2 : 1)', db: [],
  test: c => c.has('I2') && !c.any(['Br2', 'Cl2', 'NaOH', 'KOH']),
  run(m) { const s = smiOf(m); if (!s || !/Ag/.test(s)) return null; const p = run('[#6:1][CX3:2](=[O:3])[O][Ag].[#6:4][CX3](=O)[O][Ag]>>[*:1][C:2](=[O:3])O[*:4]', [s, s]); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Simonini: 2 RCOOAg + I₂ → RCOOR (ester, one alkyl lost as CO₂)', text: 'With I₂ (half the silver salt), the R–I made by the Hunsdiecker route meets a second carboxylate: it is esterified, RCOO–R.', cond: ['I₂'], extra: ['O=C=O', '[Ag+]'], why: ['RCOOAg : I₂ = 1 : 1 → R–I (Hunsdiecker); 2 : 1 → ester (Simonini)'], tags: ['radical'] }); } });

add({ id: 'd_exchange', name: 'α-H exchange (OD⁻/D₂O)', fam: 'Carbonyl', canon: 'NaOD, D₂O', db: ['Keto–enol tautomerism'],
  test: c => c.has('D2O'),
  run(m) { const a = M(m, '[CX4;!H0:1][CX3]=O').map(x => x[0]); if (!a.length) return out({ title: 'No exchange', product: null, none: 'No α-hydrogen next to a carbonyl group, so nothing exchanges with D₂O.', why: [], tags: [] });
    const t = m.clone(); let n = 0; for (const id of new Set(a)) { let g = 0; while (t.hImp(id) > 0 && g++ < 4) { const d = E.freeDir(t, id); const at = t.atoms.get(id); const h = t.addAtom('H', at.x + d[0] * BL * .82, at.y + d[1] * BL * .82); t.addBond(id, h, 1); if (at.hFix != null) at.hFix = Math.max(0, at.hFix - 1); } for (const h of E.hExp(t, id)) { t.atoms.get(h).iso = 2; n++; } }
    let R = null; try { R = fresh(C.canon(t)); } catch (e) { R = null; } if (!R) return null;
    return out({ title: `Every α-H exchanged for D (${n} D)`, product: R, steps: [step(m, R, 'OD⁻ removes an α-H (the enolate is stabilised by C=O); the enolate takes D⁺ from D₂O. Repeated until every α-H is D. The aldehyde H and other non-α H do not exchange.', ['NaOD', 'D₂O'])], why: ['Counting D atoms counts the α-hydrogens', 'Evidence for enol / enolate formation'], tags: ['tautomer'] }); } });

add({ id: 'cation_e1', name: 'Carbocation loses H⁺ (E1)', fam: 'Carbocation', canon: 'Δ', db: [],
  test: c => c.heat,
  run(m) { const s = smiOf(m); if (!s || !/\+\]/.test(s)) return null; const p = run('[C+:1][CX4;!H0:2]>>[C+0:1]=[C:2]', [s]); if (!p.length) return null;
    const sub = q => { const mm = fresh(q); return M(mm, '[CX3]=[CX3]').map(([a, b]) => mm.nbrs(a).concat(mm.nbrs(b)).filter(x => mm.atoms.get(x).el === 'C').length)[0] || 0; };
    const ord = [...new Set(p.map(x => x[0]))].sort((x, y) => sub(y) - sub(x));
    return finish(m, ord[0], { title: 'Carbocation → alkene (loss of β-H⁺, Zaitsev)', text: 'A base (water, the counter-ion) removes a β-H; its pair forms the C=C. The more substituted alkene is preferred.', cond: ['Δ'], minor: ord.slice(1), why: ['E1 final step'], tags: ['elim'] }); } });

add({ id: 'alkene_dimer', name: 'Acid-catalysed alkene dimerisation', fam: 'Alkene', canon: 'H₂SO₄ (60%), Δ', db: [],
  test: c => c.has('H2SO4'),
  run(m) { const s = smiOf(m); if (!s) return null; const p = run('[CH2:1]=[C:2]([#6:3])[#6:4].[CH2:5]=[C:6]([CH3:7])[#6:8]>>[CH3:1][C:2]([*:3])([*:4])[CH2:5][C:6](=[CH2:7])[*:8]', [s, s]); if (!p.length) return null;
    return finish(m, p[0][0], { title: 'Isobutylene-type alkene + H₂SO₄ → dimer (diisobutylene)', text: 'H⁺ adds to the CH₂ end: a 3° carbocation. It is attacked by the C=C of a second molecule (new C–C bond, a new 3° cation), which then loses H⁺ from a CH₃: the branched C₈ alkene. Hydrogenation gives "iso-octane".', cond: ['H₂SO₄', 'Δ'], why: ['Isobutylene → 2,4,4-trimethylpent-1-ene (major) and -2-ene'], tags: ['cation'] }); } });

E.REC.unshift({ id: 'tosyl_inversion', name: 'Inversion of an alcohol via its tosylate', fam: 'Alcohol', canon: 'TsCl, then AcO⁻, then OH⁻', db: [],
  test: c => c.has('TsCl') && c.any(['NaOAc', 'AcO-']) && c.any(['KOH', 'NaOH']),
  run(m) { const s = smiOf(m); if (!s) return null; const q = M(m, '[CX4:1][OX2H1]'); if (q.length !== 1) return null;
    const iso = (() => { try { return (root.Stereo && root.Stereo.iso(m)) || s; } catch (e) { return s; } })();
    const n = (iso.match(/@/g) || []).length; const inv = n && (iso.match(/\[C@@?H?\]/g) || []).length === 1 ? iso.replace(/@@|@/, x => x === '@@' ? '@' : '@@') : iso;
    return finish(m, inv, { title: 'TsCl, then AcO⁻ (SN2), then OH⁻ → alcohol with inverted configuration', text: 'TsCl makes the tosylate with retention (the C–O bond is untouched); acetate displaces –OTs by SN2 from the back (inversion); hydrolysis of the acetate cuts the acyl–O bond, so the inverted carbon keeps its O. Net: the alcohol, configuration inverted.', cond: ['TsCl', 'AcO⁻', 'OH⁻'], why: ['Retention, inversion, retention → overall inversion', 'Compare SOCl₂ alone (retention) and with pyridine (inversion)'], tags: ['SN2', 'stereo'] }); } });

/* with an organometallic in the same stage, H₃O⁺ is only the work-up: ketal hydrolysis must not run first */
{ const ah = E.REC.find(r => r.id === 'acetal_hyd'); if (ah) { const t0 = ah.test; ah.test = (c, all, k) => !c.any(['RMgX', 'RLi', 'LiAlH4', 'R2CuLi', 'DIBAL']) && t0(c, all, k); } }

/* ---------- curly arrows for the new recipes ---------- */
const R7 = root.EngineR7, R10 = root.EngineR10, R8 = root.EngineR8;
if (R7 && R10 && R8) {
  const { Mech, P, unit, el } = E; const { fin, put } = R10; const { deprot } = R8;
  const hasA = o => (o.steps || []).some(s => (s.arrows || []).length);
  const W = (id, b) => R7.wrap(id, (m, o, c, all, k) => { if (hasA(o)) return null; try { return b(m, o, c); } catch (e) { if (root.__DBG) console.log('r12 ERR', id, e.message); return null; } });
  const isEl = e => a => a.el === e; const fish = (a, b) => [a, b, 1];
  W('tosylation', (m, o) => {
    const q = M(m, '[CX4][OX2H1:1]')[0] || M(m, '[NX3;H2,H1;!$(NC=O):1]')[0]; if (!q) return null; const nu = q[q.length - 1];
    const mech = new Mech(m); const sp = put(mech, 'Cc1ccc(S(=O)(=O)Cl)cc1', nu, isEl('S'), { dist: 1.8 * BL }); const sS = sp.at; const cl = mech.m.nbrs(sS).find(x => el(mech.m, x) === 'Cl');
    mech.step([[[nu], [nu, sS]], [[sS, cl], [cl]]], 'The lone pair attacks the sulfur of TsCl and Cl⁻ leaves (substitution at S; the C–O bond is untouched).', ['TsCl']); mech.dropFrag(cl);
    const h = mech.H(nu); if (h == null) return null; const ph = P(mech.m, h), pn = P(mech.m, nu);
    const py = put(mech, 'c1ccncc1', h, isEl('N'), { dist: 1.45 * BL, dir: unit([ph[0] - pn[0], ph[1] - pn[1]]) }).at;
    mech.step([[[py], [py, h]], [[h, nu], [nu]]], 'Pyridine takes the H⁺ (pyridinium chloride forms): the tosylate (or sulfonamide).'); mech.dropFrag(py); return fin(mech, o);
  });
  W('ester_pyrolysis', (m, o) => {
    const q = M(m, '[CX3:1](=[OX1:2])[OX2:3][CX4:4][CX4;!H0:5]'); const want = R7.key(o.product);
    for (const [c2, oc, oe, ca, cb] of q) {
      const mech = new Mech(m); const h = mech.H(cb); if (h == null) continue;
      mech.step([[[oc], [oc, h]], [[h, cb], [cb, ca]], [[ca, oe], [oe, c2]], [[c2, oc], [oc]]], 'One concerted, six-membered cyclic step (Ei, syn): the C=O oxygen takes the β-H on the same side, the C–H pair becomes the C=C, the C–O bond breaks and the carboxyl group re-forms as the acid.', ['Δ']);
      const r = fin(mech, o); if (r && R7.key(r.product) === want) return r; }
    return null;
  });

  const { prot } = R8;
  const jumpTo = (mech, o, text, cond) => { const L = mech.m.clone(); const R = o.product.clone(); mech.steps.push({ L, arrows: [], R: R.clone(), text, cond: cond || [] }); mech.m = R; return { steps: mech.steps, product: R }; };
  W('acetal_roh', (m, o, c) => {
    const q = M(m, '[CX3:1]=[OX1:2]')[0]; if (!q) return null; const [cc, ox] = q; const R = /meoh|methanol|ch3oh/i.test(c.raw) ? 'C' : 'CC';
    const mech = new Mech(m);
    prot(mech, ox, 'Dry acid protonates the carbonyl O: the carbon becomes much more electrophilic.', ['dry H⁺']);
    const oa = put(mech, 'O' + R, cc, isEl('O'), { dist: 1.65 * BL }).at;
    mech.step([[[oa], [oa, cc]], [[cc, ox], [ox]]], 'An alcohol O attacks the carbonyl carbon; the π pair goes to O.', [R === 'C' ? 'CH₃OH' : 'C₂H₅OH']);
    deprot(mech, oa, 'Loss of H⁺: the hemiacetal.', 'O' + R);
    prot(mech, ox, 'The –OH is protonated: –OH₂⁺ is a good leaving group.');
    mech.step([[[oa], [oa, cc]], [[cc, ox], [ox]]], 'The OR lone pair pushes in and water leaves: a resonance-stabilised oxocarbenium ion.'); mech.dropFrag(ox);
    const ob = put(mech, 'O' + R, cc, isEl('O'), { dist: 1.65 * BL }).at;
    mech.step([[[ob], [ob, cc]], [[cc, oa], [oa]]], 'A second alcohol molecule attacks the oxocarbenium carbon.');
    deprot(mech, ob, 'Loss of H⁺: the acetal (ketal). Removing the water drives the equilibrium.', 'O' + R);
    return fin(mech, o);
  });
  W('hydroxy_acid_heat', (m, o) => {
    if (!M(o.product, '[CX3;R](=O)[OX2;R]').length || M(o.product, '[CX3;R](=O)[OX2;R][CX4;R][CX3;R](=O)').length) return null;
    const q = M(m, '[CX3:1](=[OX1:2])[OX2H1]')[0]; const oh = M(m, '[CX4][OX2H1:1]')[0]; if (!q || !oh) return null; const [c1, o1] = q; const on = oh[1];
    const mech = new Mech(m);
    prot(mech, o1, 'A trace of acid (or the acid itself) protonates the C=O of –COOH.', ['Δ']);
    mech.step([[[on], [on, c1]], [[c1, o1], [o1]]], 'The –OH at the far end of the same chain attacks the carbonyl carbon: a five- or six-membered ring closes (intramolecular, strain-free).', []);
    return jumpTo(mech, o, 'Proton transfers, then water leaves and the C=O re-forms: the lactone (cyclic ester).');
  });
  W('hunsdiecker', (m, o, c) => {
    const q = M(m, '[#6:1][CX3:2](=[OX1:3])[OX2:4][Ag:5]')[0]; if (!q) return null; const [r, cc, od, os, ag] = q; const X = c.has('Br2') ? 'Br' : 'Cl';
    const mech = new Mech(m);
    const s2 = put(mech, X + X, os, isEl(X), { dist: 1.7 * BL }); const x1 = s2.at, x2 = s2.ids.find(z => z !== x1);
    mech.step([[[os, ag], [os, x1]], [[x1, x2], [x2]]], `The carboxylate O attacks ${X}₂; AgX precipitates (it drives the step): an acyl hypohalite RCOO–${X}.`, [X + '₂', 'CCl₄']);
    mech.dropFrag(ag); mech.dropFrag(x2);
    mech.step([fish([os, x1], [os]), fish([os, x1], [x1])], `Heat breaks the weak O–${X} bond homolytically: a carboxyl radical RCOO• and ${X}•.`, ['Δ']);
    mech.step([fish([os], [cc, os]), fish([cc, r], [cc, os]), fish([cc, r], [r])], 'The carboxyl radical loses CO₂ at once: an alkyl radical R•, one carbon shorter.');
    mech.dropFrag(cc);
    mech.step([fish([r], [r, x1]), fish([x1], [r, x1])], `R• takes ${X} (from ${X}₂ in the chain; shown pairing with ${X}•): R–${X}.`);
    return fin(mech, o);
  });
  W('acyl_azide_heat', (m, o) => {
    const q = M(m, '[#6:1][CX3:2](=[OX1:3])[N:4]=[N+:5]=[N-:6]')[0]; if (!q) return null; const [r, c2, o3, n1, n2, n3] = q;
    const mech = new Mech(m);
    mech.step([[[c2, r], [r, n1]], [[n1, n2], [n2]], [[n1, n2], [n1, c2]], [[n3], [n2, n3]]], 'Concerted on heating: R moves from the carbonyl carbon to the nitrogen (retaining its configuration) as N₂ leaves; the N lone pair becomes the C=N. The isocyanate R–N=C=O.', ['Δ']);
    mech.dropFrag(n2);
    return fin(mech, o);
  });
  W('pinacol_mg', (m, o, c) => {
    const q = M(m, '[CX3:1]=[OX1:2]')[0]; if (!q || m.frags().length !== 1) return null; const [c1, o1] = q;
    const mech = new Mech(m); const other = (c.has('acetone') && C.canon(m) !== 'CC(C)=O') ? 'CC(C)=O' : C.canon(m);
    const map2 = mech.add(other, c1, { dist: 3 * BL }); const sub = mech.m.sub([...map2.values()]); const q2 = M(sub, '[CX3:1]=[OX1:2]')[0]; if (!q2) return null; const [c2, o2] = q2;
    const mg = put(mech, '[Mg]', o1, isEl('Mg'), { dist: 1.9 * BL }).at;
    mech.step([fish([mg], [o1]), fish([c1, o1], [c1]), fish([c1, o1], [o1]), fish([mg], [o2]), fish([c2, o2], [c2]), fish([c2, o2], [o2])], 'Mg (activated by Hg) gives one electron to each C=O: two ketyl radical anions, C•–O⁻, held by Mg²⁺.', ['Mg–Hg'], () => { for (const x of [o1, o2]) { const a = mech.m.atoms.get(x); if (a) a.rad = 0; } E.tidyLongBonds(mech.m); });
    mech.dropFrag(mg);
    mech.step([fish([c1], [c1, c2]), fish([c2], [c1, c2])], 'The two carbon radicals pair: a new C–C bond and a magnesium 1,2-dialkoxide.');
    for (const ox of [o1, o2]) { const ids = mech.add('[H]O', ox, { dist: 1.5 * BL }); const hh = ids.get(1), w = ids.get(2); mech.step([[[ox], [ox, hh]], [[hh, w], [w]]], 'Water protonates the alkoxide.', ox === o1 ? ['H₂O'] : []); mech.dropFrag(w); }
    return fin(mech, o);
  });
}
root.EngineR12 = { run };
})(typeof window !== 'undefined' ? window : globalThis);
