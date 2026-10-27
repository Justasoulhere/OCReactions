/* engine_r9.js: stereo-aware eliminations and small additions (Oct 14). FRONT-registered. */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore;
const { BL, el, P, unit, hTot, out, done, Mech } = E;
const REC = E.REC; const front = r => REC.unshift(r);
const M = (m, sm) => { try { return C.matchAll(m, sm, {}); } catch (e) { return []; } };
const AR = m => E.aromSet(m);
function pyridine() { const sp = E.species('c1ccncc1'); const n = [...sp.atoms.values()].find(a => a.el === 'N').id; return { sp, n }; }
/* POCl₃ / pyridine: –OH → –OPOCl₂, then E2 (no carbocation, so no rearrangement; anti-periplanar H needed) */
front({ id: 'pocl3_e2', name: 'POCl₃ / pyridine dehydration (E2)', fam: 'Alcohol', canon: 'POCl₃, pyridine', db: [],
  test: c => c.has('POCl3'),
  run(m) {
    const al = M(m, '[CX4:1][OX2H1:2]').filter(([cc]) => !AR(m).has(cc)); if (!al.length) return null; const [cc, o] = al[0];
    const bs0 = E.betaHs(m, cc, o); if (!bs0.length) return out({ title: 'No β-H: no elimination', product: null, none: 'This alcohol has no β-hydrogen, so POCl₃/pyridine cannot dehydrate it.' });
    const anti = bs0.filter(b => E.antiOK(m, cc, o, b)); const bs = (anti.length ? anti : bs0).slice().sort((p, q) => E.alkeneScore(m, cc, q, AR(m)) - E.alkeneScore(m, cc, p, AR(m)));
    const best0 = bs0.slice().sort((p, q) => E.alkeneScore(m, cc, q, AR(m)) - E.alkeneScore(m, cc, p, AR(m)))[0]; const bC = bs[0]; const forced = bC !== best0;
    const mech = new Mech(m);
    const pc = E.species('O=P(Cl)(Cl)Cl'); const pl = [...pc.atoms.values()].find(a => a.el === 'P').id; const cl = pc.nbrs(pl).find(n => el(pc, n) === 'Cl');
    const ids = mech.add(pc, o, { att: pl, dist: 1.7 * BL }); const p = ids.get(pl), x = ids.get(cl);
    mech.step([[[o], [o, p]], [[p, x], [x]]], 'The –OH oxygen attacks phosphorus and Cl⁻ leaves.', ['POCl₃']); mech.dropFrag(x);
    let py = pyridine(); let n = mech.add(py.sp, mech.H(o), { att: py.n, dist: 1.45 * BL }).get(py.n); const h0 = mech.m.nbrs(o).find(z => el(mech.m, z) === 'H');
    mech.step([[[n], [n, h0]], [[h0, o], [o]]], 'Pyridine removes the O–H proton: –O–POCl₂, an excellent leaving group (no carbocation is formed).', ['pyridine']); mech.dropFrag(n);
    const h = mech.H(bC); py = pyridine(); const ph = P(mech.m, h), pb = P(mech.m, bC);
    n = mech.add(py.sp, h, { att: py.n, dist: 1.45 * BL, dir: unit([ph[0] - pb[0], ph[1] - pb[1]]) }).get(py.n);
    mech.step([[[n], [n, h]], [[bC, h], [bC, cc]], [[cc, o], [o]]], 'E2 in one step: pyridine takes the β-H that is anti-periplanar to the leaving group, the C=C forms and ⁻OPOCl₂ leaves.');
    mech.dropFrag(n); mech.dropFrag(o);
    return done(mech, { title: forced ? 'POCl₃ / pyridine: E2 with the anti-periplanar H → less substituted alkene' : 'POCl₃ / pyridine: E2 dehydration (no rearrangement)', why: ['–OH becomes –OPOCl₂ without a carbocation: no rearrangement (unlike conc. H₂SO₄)', 'E2: the β-H must be anti-periplanar (trans-diaxial on a cyclohexane ring)', forced ? 'Here the H that would give the more substituted alkene is cis to the O, so the other alkene forms' : 'Saytzeff product when the anti H allows it', 'Pyridine is the base and also takes up HCl'], tags: ['E2', 'anti'] });
  } });

const tidy = m => { const t = m.clone(); try { C.relayout(t, BL); } catch (e) { } return t; };
const fresh = smi => C.fromSmilesMB(smi, BL);
const clr = (t, ids) => ids.forEach(id => { const a = t.atoms.get(id); if (a) a.hFix = null; });
function addAt(t, at, smi, o = 1) { const g = C.fromSmilesMB(smi, BL); for (const a of g.atoms.values()) a.hFix = null; const first = [...g.atoms.keys()][0]; const p = P(t, at); const idm = t.merge(g, p[0] + BL, p[1]); t.addBond(at, idm.get(first), o); clr(t, [at]); return idm.get(first); }
const step = (L, R, text, cond) => ({ L: L.clone(), arrows: [], R: R.clone(), text, cond: cond || [] });
function fragsOf(t) { return t.frags().map(ids => ({ ids, n: ids.filter(id => el(t, id) === 'C').length, h: ids.length })).sort((p, q) => q.n - p.n || q.h - p.h).map(f => E.cleanProduct(t, f.ids)); }
function simple(L, t, o) { const R = tidy(t); let main = R, extra = o.extra || []; if (R.frags().length > 1) { const all = fragsOf(R); main = o.pick ? (all.find(o.pick) || all[0]) : all[0]; extra = extra.concat(all.filter(f => f !== main)); }
  return out(Object.assign({}, o, { product: main, extra, steps: [step(L, R, o.text, o.cond)] })); }
const carbonyl = m => M(m, '[CX3:1]=[OX1:2]').find(([cc, ox]) => m.nbrs(cc).filter(n => n !== ox && el(m, n) !== 'H').every(n => el(m, n) === 'C'));
function markov(m) { const ar = AR(m); const al = M(m, '[C:1]=[C:2]').filter(([a, b]) => !ar.has(a) && !ar.has(b)); if (!al.length) return null; const [a, b] = al[0]; const sa = E.cationScore(m, a, ar), sb = E.cationScore(m, b, ar); return sa >= sb ? { cp: a, hC: b } : { cp: b, hC: a }; }

/* Zinin: (NH₄)₂S / Na₂S reduces one –NO₂ only */
front({ id: 'zinin', name: 'Zinin reduction', fam: 'Reduction', canon: '(NH₄)₂S', db: ['Zinin reduction (selective)'],
  test: c => c.has('Zinin'),
  run(m) { const nn = M(m, '[c:1][N+:2](=[O:3])[O-:4]'); if (!nn.length) return null; const [, n, o1, o2] = nn[0]; const t = m.clone(); t.removeAtom(o1); t.removeAtom(o2); const a = t.atoms.get(n); a.chg = 0; a.hFix = null;
    return simple(m, t, { title: nn.length > 1 ? 'Zinin reduction: only one –NO₂ → –NH₂' : 'Zinin reduction: Ar–NO₂ → Ar–NH₂', text: 'Sulphide / hydrosulphide is a mild reducing agent: it reduces one nitro group to –NH₂ and stops, because the new –NH₂ makes the ring electron-rich and the second –NO₂ harder to reduce.', why: ['m-Dinitrobenzene → m-nitroaniline (Sn/HCl would reduce both)'], tags: ['reduction'] }); } });
/* Meerwein–Ponndorf–Verley and Oppenauer */
front({ id: 'mpv', name: 'Meerwein–Ponndorf–Verley / Oppenauer', fam: 'Reduction', canon: 'Al(OiPr)₃, isopropanol', db: ['Meerwein–Ponndorf–Verley reduction'],
  test: c => c.has('AlOiPr'),
  run(m, c) { const k = carbonyl(m);
    if (k) { const [cc, ox] = k; const t = m.clone(); t.bondBetween(cc, ox).o = 1; clr(t, [cc, ox]); const o = simple(m, t, { title: 'MPV reduction: C=O → CH–OH (Al(OiPr)₃)', text: 'Aluminium isopropoxide passes a hydride from its isopropoxide C–H to the carbonyl carbon through a six-membered cyclic transition state; isopropanol becomes acetone, which is distilled off.', why: ['Very selective: C=C, –NO₂ and esters are untouched', 'Reverse reaction (with acetone, Al(OtBu)₃) is the Oppenauer oxidation'], tags: ['reduction'] }); o.extra = (o.extra || []).concat([fresh('CC(C)=O')]); return o; }
    const al = M(m, '[CX4;H1:1]([#6])([#6])[OX2H1:2]')[0]; if (!al) return null; const cc = al[0], ox = al[3]; const t = m.clone(); t.bondBetween(cc, ox).o = 2; clr(t, [cc, ox]);
    const o = simple(m, t, { title: 'Oppenauer oxidation: 2° alcohol → ketone', text: 'The same cyclic hydride transfer run backwards: the alcohol gives its C–H hydride to acetone (in excess), becoming the ketone.', why: ['Reverse of the MPV reduction'], tags: ['oxidation'] }); o.extra = (o.extra || []).concat([fresh('CC(C)O')]); return o; } });
/* osazone */
front({ id: 'osazone', name: 'Osazone formation', fam: 'Carbonyl', canon: 'PhNHNH₂ (excess)', db: ['Osazone formation'],
  test: c => c.has('PhNHNH2') && !(c.eq != null && c.eq < 3),
  run(m) { const q = M(m, '[CX3:1](=[O:2])[CX4;!H0:3][OX2H1:4]').sort((p, q2) => hTot(m, q2[2]) - hTot(m, p[2]))[0]; if (!q) return null; const [c1, o1, c2, o2] = q; // ketoses: the CH₂OH end (C1) reacts
    const t = m.clone(); t.removeAtom(o1); t.removeAtom(o2); clr(t, [c1, c2]); addAt(t, c1, 'NNc1ccccc1', 2); addAt(t, c2, 'NNc1ccccc1', 2);
    const o = simple(m, t, { title: 'Excess PhNHNH₂ → osazone (C1 and C2 both C=N–NHPh)', text: 'The first PhNHNH₂ forms the phenylhydrazone; a second one oxidises the neighbouring CH–OH to C=O (and is itself cut to aniline + NH₃); a third condenses with that new C=O. Three equivalents in all.', why: ['Glucose, fructose and mannose give the same osazone (they differ only at C1/C2)', 'Only C1 and C2 react: the reaction stops there'], tags: ['condensation'] });
    o.extra = (o.extra || []).concat([fresh('Nc1ccccc1'), fresh('N')]); return o; } });
/* periodic acid / lead tetra-acetate: 1,2-diols cleave */
front({ id: 'hio4', name: 'HIO₄ cleavage of 1,2-diols', fam: 'Oxidation', canon: 'HIO₄', db: ['Malaprade (periodate) cleavage'],
  test: c => c.has('HIO4'),
  run(m) { const q = M(m, '[CX4:1]([OX2H1:2])[CX4:3][OX2H1:4]'); if (!q.length) return null; const t = m.clone(); let n = 0;
    for (const [a, oa, b, ob] of q) { const bd = t.bondBetween(a, b); if (!bd) continue; t.removeBond(bd.id); for (const [x, ox] of [[a, oa], [b, ob]]) { const bx = t.bondBetween(x, ox); if (bx) { bx.o = 2; clr(t, [x, ox]); } } n++; }
    if (!n) return null;
    return simple(m, t, { title: 'HIO₄: 1,2-diol → two carbonyl compounds', text: 'Periodic acid forms a cyclic periodate ester with the two neighbouring –OH groups; it breaks apart, cutting the C–C bond: each carbon becomes C=O.', why: ['Ethane-1,2-diol → 2 HCHO', 'Used to count adjacent –OH groups in sugars (one HIO₄ per C–C cut)', 'Pb(OAc)₄ does the same'], tags: ['oxidation'] }); } });
/* NOCl and IN₃: electrophile to the less substituted carbon, anion to the more substituted */
front({ id: 'nocl_add', name: 'NOCl addition (Tilden reagent)', fam: 'Alkene', canon: 'NOCl', db: [],
  test: c => c.has('NOCl'),
  run(m) { const mk = markov(m); if (!mk) return null; const { cp, hC } = mk; const mech = new Mech(m);
    const sp = E.species('O=NCl'); const nl = [...sp.atoms.values()].find(a => a.el === 'N').id, cl = [...sp.atoms.values()].find(a => a.el === 'Cl').id;
    const ids = mech.add(sp, hC, { att: nl, dist: 1.65 * BL }); const n = ids.get(nl), x = ids.get(cl);
    mech.step([[[cp, hC], [hC, n]], [[n, x], [x]]], 'The C=C π electrons attack the N of NOCl (N=O⁺ is the electrophile); Cl⁻ leaves. N goes to the carbon that leaves the more stable carbocation.', ['NOCl']);
    E.captureBy(mech, cp, x, 'Cl⁻ captures the carbocation (Markovnikov).');
    return done(mech, { title: 'NOCl adds Markovnikov: –N=O on the CH₂ end, Cl on the more substituted C', why: ['NO⁺ is the electrophile (like H⁺ in HCl)', 'The nitroso chloride often dimerises (blue → colourless)'], tags: ['AdE'] }); } });
front({ id: 'in3_add', name: 'IN₃ addition', fam: 'Alkene', canon: 'IN₃', db: [],
  test: c => c.has('IN3'),
  run(m) { const mk = markov(m); if (!mk) return null; const { cp, hC } = mk; const t = m.clone(); const bd = t.bondBetween(cp, hC); bd.o = 1; clr(t, [cp, hC]); addAt(t, hC, 'I'); addAt(t, cp, 'N=[N+]=[N-]');
    return simple(m, t, { title: 'IN₃: I on the CH₂ end, N₃ on the more substituted carbon (anti)', text: 'I⁺ forms a cyclic iodonium ion; azide (N₃⁻) opens it at the carbon that carries more + charge, from the back side.', why: ['Same pattern as halohydrin formation: electrophile (I) + nucleophile (N₃⁻), Markovnikov, anti'], tags: ['AdE', 'anti'] }); } });
/* SF₄ / DAST: C=O → CF₂, C–OH → C–F */
front({ id: 'sf4', name: 'SF₄ / DAST fluorination', fam: 'Carbonyl', canon: 'SF₄', db: [],
  test: c => c.has('SF4'),
  run(m) { const k = carbonyl(m); const t = m.clone();
    if (k) { const [cc, ox] = k; t.removeAtom(ox); clr(t, [cc]); addAt(t, cc, 'F'); addAt(t, cc, 'F'); const o = simple(m, t, { title: 'SF₄: C=O → CF₂ (gem-difluoride)', text: 'SF₄ replaces the carbonyl oxygen by two fluorines; SOF₂ is the by-product.', why: ['DAST does the same under milder conditions'], tags: ['substitution'] }); o.extra = (o.extra || []).concat([fresh('O=S(F)F')]); return o; }
    const al = M(m, '[CX4:1][OX2H1:2]')[0]; if (!al) return null; t.removeAtom(al[1]); clr(t, [al[0]]); addAt(t, al[0], 'F');
    return simple(m, t, { title: 'SF₄ / DAST: R–OH → R–F', text: 'The –OH is replaced by F.', why: [], tags: ['substitution'] }); } });
/* PPh₃ deoxygenates epoxides */
front({ id: 'pph3_epox', name: 'Epoxide + PPh₃ → alkene', fam: 'Ether', canon: 'PPh₃', db: [],
  test: c => c.has('PPh3'),
  run(m) { const ep = M(m, '[C:1]1[O:2][C:3]1')[0]; if (!ep) return null; const [a, ox, b] = ep; const t = m.clone(); t.removeAtom(ox); t.bondBetween(a, b).o = 2; clr(t, [a, b]);
    const o = simple(m, t, { title: 'PPh₃: epoxide → alkene (deoxygenation)', text: 'PPh₃ opens the epoxide (SN2 at carbon), the alkoxide attacks P, and a four-membered ring breaks into the alkene and Ph₃P=O (the strong P=O bond drives it).', why: ['The double bond geometry is inverted (cis-epoxide → trans-alkene): two inversions… overall anti elimination'], tags: ['elimination'] }); o.extra = (o.extra || []).concat([fresh('O=P(c1ccccc1)(c1ccccc1)c1ccccc1')]); return o; } });
/* NBS without light oxidises alcohols */
front({ id: 'nbs_alc', name: 'NBS oxidation of alcohols', fam: 'Alcohol', canon: 'NBS', db: [],
  test: c => c.has('NBS') && !c.light && !c.has('ROOR'),
  run(m) { if (M(m, '[C;!a]=[C;!a]').length) return null; const al = M(m, '[CX4;!H0:1][OX2H1:2]').filter(([cc]) => !AR(m).has(cc))[0]; if (!al) return null; const [cc, ox] = al;
    const t = m.clone(); t.bondBetween(cc, ox).o = 2; clr(t, [cc, ox]); const o = simple(m, t, { title: 'NBS (no light): alcohol → carbonyl compound', text: 'In the dark NBS acts as a source of Br⁺ (an oxidant), not of Br• radicals: the CH–OH is oxidised to C=O; succinimide and HBr are by-products.', why: ['With light / peroxide and an allylic or benzylic H, NBS brominates instead'], tags: ['oxidation'] }); o.extra = (o.extra || []).concat([fresh('O=C1CCC(=O)N1')]); return o; } });
/* Grignard as the drawn structure */
front({ id: 'grignard_clcn', name: 'Grignard + ClCN', fam: 'Organometallic', canon: 'ClCN', db: ['Grignard + ClCN / (CN)2'],
  test: c => c.has('ClCN'),
  run(m) { const g = M(m, '[#6:1][Mg:2]')[0]; if (!g) return null; const [r, mg] = g; const t = m.clone(); for (const id of E.fragIds(m, mg).filter(id => id !== r && !E.subtree(m, r, mg).includes(id))) t.removeAtom(id); clr(t, [r]); addAt(t, r, 'C#N');
    return simple(m, t, { title: 'R–MgX + ClCN → R–C≡N (nitrile)', text: 'The Grignard carbon attacks the cyanide carbon of ClCN (or cyanogen) and Cl⁻ (or CN⁻) leaves: a nitrile with one more carbon.', why: ['A way to add one carbon as –CN'], tags: ['organometallic'] }); } });
front({ id: 'grignard_rx', name: 'Grignard + alkyl halide', fam: 'Organometallic', canon: 'R–X', db: ['Grignard + alkyl halide'],
  test: c => c.has('RX') && !c.any(['Na', 'Mg', 'NaOEt', 'KOH']),
  run(m, c) { const g = M(m, '[#6:1][Mg:2]')[0]; if (!g) return null; const [r, mg] = g; const rx = c.get('RX') || {}; const t = m.clone(); for (const id of E.fragIds(m, mg).filter(id => id !== r && !E.subtree(m, r, mg).includes(id))) t.removeAtom(id); clr(t, [r]); addAt(t, r, rx.smi || 'C');
    return simple(m, t, { title: 'R–MgX + R′–X → R–R′', text: 'The carbanion-like Grignard carbon displaces the halide from the alkyl halide (works well for reactive allyl, benzyl or methyl halides).', why: ['Side reaction: Wurtz-type coupling and halogen exchange'], tags: ['organometallic'] }); } });
/* thiolates: SN2 at alkyl halides first, otherwise acyl substitution */
front({ id: 'rs_nu', name: 'Thiolate nucleophile', fam: 'Alkyl halide', canon: 'C₂H₅S⁻', db: [],
  test: c => c.has('RS-'),
  run(m, c) { const R = /ch3|me/i.test((c.get('RS-') || {}).txt || '') ? 'C' : 'CC'; const hs = E.halides(m); const t = m.clone();
    if (hs.length) { const [cc, x] = hs[0]; t.removeAtom(x); clr(t, [cc]); addAt(t, cc, 'S' + R); return simple(m, t, { title: 'RS⁻ + R′–X → thioether (SN2)', text: 'Thiolate is a very strong, soft nucleophile: it does SN2 at the sp³ C–X (even when an acyl halide is present elsewhere).', why: ['Soft nucleophile, soft (sp³ carbon) electrophile'], tags: ['SN2'] }); }
    const q = M(m, '[CX3:1](=O)[Cl,Br,I:2]')[0]; if (!q) return null; const [cc, , x] = q; t.removeAtom(x); clr(t, [cc]); addAt(t, cc, 'S' + R);
    return simple(m, t, { title: 'RS⁻ + acyl halide → thioester', text: 'Thiolate adds to the C=O (tetrahedral intermediate) and the halide leaves.', why: [], tags: ['acyl substitution'] }); } });

/* ---------------- boranes ---------------- */
/* BH₃ alone: unhindered alkenes go on to R₃B, trisubstituted stop at R₂BH (disiamylborane), tetrasubstituted at RBH₂ (thexylborane); 1,5-cyclooctadiene closes to 9-BBN */
(function () { const r = REC.find(x => x.id === 'hydroboration'); if (!r) return; const run0 = r.run;
  r.run = function (m, c, all, k) {
    const o = run0.call(this, m, c, all, k); if (!o || !o.product || !all || all.length > 1 || c.st.sp.size > 1 || M(m, 'C#C').length) return o;
    const ar = AR(m); const al = M(m, '[C:1]=[C:2]').filter(([a, b]) => !ar.has(a) && !ar.has(b)); if (!al.length) return o;
    const cod = M(m, '[C:1]1=[C:2][C:3][C:4][C:5]=[C:6][C:7][C:8]1')[0];
    let t = null, title = '', why = [];
    if (cod) { const [c1, c2, , , c5, c6] = cod; t = m.clone(); t.bondBetween(c1, c2).o = 1; t.bondBetween(c5, c6).o = 1; clr(t, [c1, c2, c5, c6]); const b = addAt(t, c1, 'B'); t.addBond(b, c5, 1);
      title = 'BH₃ + cycloocta-1,5-diene → 9-BBN'; why = ['The second B–H adds across the other C=C of the same ring: a bicyclic borane (9-borabicyclo[3.3.1]nonane)', '9-BBN is a bulky, very regioselective hydroborating agent']; }
    else { const [a, b] = al[0]; const subs = [a, b].reduce((n, x) => n + E.cN(m, x).length, 0) - 2; const n = subs >= 4 ? 1 : subs === 3 ? 2 : 3; if (n === 1) return o;
      const bc = hTot(m, a) >= hTot(m, b) ? a : b, oc = bc === a ? b : a; t = m.clone(); t.bondBetween(a, b).o = 1; clr(t, [a, b]); const B = addAt(t, bc, 'B');
      for (let i = 1; i < n; i++) { const u = m.clone(); u.bondBetween(a, b).o = 1; clr(u, [a, b]); const idm = t.merge(u, 320 * i, 0); t.addBond(B, idm.get(bc), 1); }
      title = n === 3 ? 'BH₃ + 3 alkene → trialkylborane R₃B' : 'BH₃ + 2 hindered alkene → R₂BH (disiamylborane type)';
      why = n === 3 ? ['Each B–H adds to another alkene the same way (anti-Markovnikov, syn): three alkyl groups on B', 'Then H₂O₂ / OH⁻ turns every C–B into C–OH'] : ['The alkene is too crowded for a third one to add: the dialkylborane stops at R₂BH', 'Disiamylborane hydroborates terminal alkynes only once (→ aldehyde after oxidation)']; }
    const R = tidy(t); const steps = (o.steps || []).concat([step(o.product, R, cod ? 'The remaining B–H adds intramolecularly to the second C=C, closing the bicyclic ring.' : 'The remaining B–H bonds add to further alkene molecules the same way.', [])]);
    return Object.assign({}, o, { product: R, steps, title, why: why.concat(o.why || []) });
  }; })();
/* a drawn organoborane + the second reagent (two arrows on the Board) */
front({ id: 'borane_follow', name: 'Organoborane reactions', fam: 'Alkene', canon: 'H₂O₂, NaOH', db: ['Hydroboration–oxidation'],
  test: c => c.any(['H2O2/OH', 'H2O2', 'AcOH', 'NH2Cl']) || (c.any(['Br2', 'I2']) && c.any(['KOH', 'NaOMe'])),
  run(m, c) { const bs = [...m.atoms.values()].filter(a => a.el === 'B').map(a => a.id); if (!bs.length) return null;
    const cs = bs.flatMap(b => m.nbrs(b).filter(n => el(m, n) === 'C')); if (!cs.length) return null; const t = m.clone(); for (const b of bs) t.removeAtom(b);
    const kind = c.any(['H2O2/OH', 'H2O2']) ? ['O', 'alcohol', 'H₂O₂ / OH⁻: each C–B becomes C–OH with retention (HOO⁻ adds to B, then the alkyl group moves from B to O).'] : c.has('AcOH') ? [null, 'alkane', 'CH₃COOH replaces each B by H (protonolysis, retention).'] : c.has('NH2Cl') ? ['N', 'amine', 'NH₂Cl turns C–B into C–NH₂.'] : [c.has('I2') ? 'I' : 'Br', 'alkyl halide', 'X₂ / base turns C–B into C–X.'];
    for (const x of cs) { clr(t, [x]); if (kind[0]) addAt(t, x, kind[0]); }
    const o = simple(m, t, { title: `Organoborane → ${kind[1]}`, text: kind[2], why: ['The B ended up on the less substituted carbon, so the group does too (anti-Markovnikov overall)', 'Syn addition then retention: the H and the new group are cis'], tags: ['hydroboration'] });
    const seen = new Set([C.canon(o.product)]); o.extra = (o.extra || []).filter(x => { const k = C.canon(x); if (seen.has(k)) return false; seen.add(k); return true; }); return o; } });

/* ---------------- alcohols, sugars ---------------- */
/* 1,2-diols + HI: diiodide loses I₂ → alkene; with excess HI the alkene adds HI again (glycerol → allyl iodide / 2-iodopropane) */
front({ id: 'diol_hi', name: '1,2-Diol + HI', fam: 'Alcohol', canon: 'HI, Δ', db: [],
  test: c => c.has('HI') && !c.has('redP'),
  run(m, c) { if (!M(m, '[CX4]([OX2H1])[CX4][OX2H1]').length) return null; const ex = c.eq === Infinity || c.flag('excess') || /excess|xs/i.test(c.raw || '');
    const t = m.clone(); for (const [cc, o] of M(t, '[CX4:1][OX2H1:2]')) { t.removeAtom(o); clr(t, [cc]); addAt(t, cc, 'I'); }
    let guard = 0, changed = true; const notes = [];
    while (changed && guard++ < 6) { changed = false;
      const v = M(t, '[CX4:1]([I:2])[CX4:3][I:4]')[0]; if (v) { const [a, ia, b, ib] = v; t.removeAtom(ia); t.removeAtom(ib); t.bondBetween(a, b).o = 2; clr(t, [a, b]); notes.push('the vicinal diiodide loses I₂ → C=C'); changed = true; continue; }
      if (ex) { const al = M(t, '[C:1]=[C:2]')[0]; if (al) { const [a, b] = al; const ar = AR(t); const sa = E.cationScore(t, a, ar), sb = E.cationScore(t, b, ar); const cp = sa >= sb ? a : b; t.bondBetween(a, b).o = 1; clr(t, [a, b]); addAt(t, cp, 'I'); notes.push('HI adds (Markovnikov)'); changed = true; } } }
    const glycerol = M(m, '[CX4]([OX2H1])[CX4]([OX2H1])[CX4][OX2H1]').length > 0;
    return simple(m, t, { title: glycerol ? (ex ? 'Glycerol + excess HI → 2-iodopropane' : 'Glycerol + HI → allyl iodide') : (ex ? '1,2-Diol + excess HI → alkyl iodide' : '1,2-Diol + HI → alkene (via the diiodide)'), text: `Each –OH becomes –I; ${[...new Set(notes)].join(', then ')}.`, why: ['A vicinal diiodide is unstable: it loses I₂ to give the alkene', 'Glycerol: small amount of HI → allyl iodide; excess → isopropyl iodide', 'Ethylene glycol + excess HI → ethyl iodide'], tags: ['substitution'] }); } });
/* glucose: no bisulphite adduct (cyclic hemiacetal) */
front({ id: 'sugar_nahso3', name: 'Glucose + NaHSO₃', fam: 'Biomolecule', canon: 'NaHSO₃', db: [],
  test: c => c.has('NaHSO3'),
  run(m) { if (M(m, '[CX4][OX2H1]').length < 4 || !carbonyl(m)) return null;
    return out({ title: 'No adduct with NaHSO₃', product: null, none: 'Glucose does not form a bisulphite adduct (and gives no Schiff test): in solution almost all of it is the cyclic hemiacetal, with very little free –CHO. This was evidence for the ring structure.' }); } });

/* ---------------- nitrogen ---------------- */
front({ id: 'carbylamine_rev', name: 'Carbylamine reaction (chloroform drawn)', fam: 'Amine', canon: 'C₂H₅NH₂, alc. KOH', db: ['Carbylamine reaction (isocyanide test)'],
  test: c => (c.has('RNH2') || c.has('aniline')) && c.any(['KOH', 'alcKOH', 'aqKOH']),
  run(m, c) { if (C.canon(m) !== C.canon(fresh('ClC(Cl)Cl'))) return null; const R = c.has('aniline') ? 'c1ccccc1' : (c.get('RNH2').smi || 'C');
    const P0 = fresh(R + '[N+]#[C-]'); return out({ title: 'Carbylamine reaction: R–NH₂ + CHCl₃ + 3 KOH → R–N≡C (foul smell)', product: P0, steps: [step(m, P0, 'Base turns CHCl₃ into dichlorocarbene (:CCl₂); the amine N attacks it and loses two HCl: the isocyanide.', ['alc. KOH'])], why: ['Test for 1° amines only (aliphatic and aromatic): 2° and 3° amines do not give it'], tags: ['carbene'] }); } });
front({ id: 'oxalate_amine', name: "Hofmann's method (diethyl oxalate)", fam: 'Amine', canon: 'C₂H₅NH₂', db: [],
  test: c => c.any(['RNH2', 'R2NH', 'NH3', 'aniline']),
  run(m, c) { const q = M(m, '[O]=[C:1]([OX2:3][#6])[C:2](=[O])[OX2:4][#6]')[0]; if (!q) return null; const c1 = q[1], o1 = q[2], c2 = q[4], o2 = q[6];
    const t = m.clone(); const R = c.has('RNH2') ? 'N' + (c.get('RNH2').smi || 'C') : c.has('aniline') ? 'Nc1ccccc1' : c.has('NH3') ? 'N' : 'N(' + (c.get('R2NH').smi || 'C') + ')' + (c.get('R2NH').smi || 'C'); const two = !c.has('R2NH');
    const cut = (cc, oo) => { for (const id of E.subtree(m, oo, cc)) if (t.atoms.has(id)) t.removeAtom(id); clr(t, [cc]); addAt(t, cc, R); };
    cut(c1, o1); if (two) cut(c2, o2);
    return simple(m, t, { title: two ? "Diethyl oxalate + 1° amine → N,N′-dialkyloxamide (solid)" : "Diethyl oxalate + 2° amine → dialkyl oxamic ester (liquid)", text: two ? 'The 1° amine (two N–H) replaces both OEt groups: a crystalline oxamide.' : 'A 2° amine (one N–H) replaces only one OEt group; a 3° amine does not react.', why: ["Hofmann's method of separating 1°, 2° and 3° amines", 'Each product is hydrolysed back to its amine afterwards'], tags: ['acyl substitution'] }); } });
front({ id: 'nitrile_hcho', name: 'Nitrile carbanion + HCHO', fam: 'Acid derivative', canon: 'HCHO, NaOH', db: [],
  test: c => c.has('HCHO') && c.any(['KOH', 'aqKOH', 'NaOEt']),
  run(m) { const q = M(m, '[CX4;!H0:1][C:2]#[N:3]')[0]; if (!q) return null; const t = m.clone(); clr(t, [q[0]]); addAt(t, q[0], 'CO');
    return simple(m, t, { title: 'α-Carbanion of the nitrile + HCHO → β-hydroxy nitrile', text: 'Base removes an α-H (stabilised by –C≡N); the carbanion adds to formaldehyde, and the alkoxide picks up H⁺.', why: ['Aldol-type addition: –CN makes the α-H acidic, like C=O does'], tags: ['AdN'] }); } });

/* ---------------- Haworth synthesis ---------------- */
front({ id: 'haworth', name: 'Haworth synthesis', fam: 'Aromatic', canon: 'benzene, AlCl₃', db: ['Haworth synthesis'],
  test: c => c.has('FeX3') || c.has('ZnHg') || c.has('H2SO4'),
  run(m, c) { const anh = M(m, 'O=C1CCC(=O)O1').length > 0, bz = C.canon(m) === C.canon(fresh('c1ccccc1'));
    if (!(anh && c.has('benzene')) && !(bz && c.has('succanh'))) return null; if (!c.has('FeX3')) return null;
    const red = c.has('ZnHg') || c.has('NH2NH2'), cyc = c.has('H2SO4') || c.has('H3PO4');
    const smi = red ? (cyc ? 'O=C1CCCc2ccccc21' : 'O=C(O)CCCc1ccccc1') : 'O=C(O)CCC(=O)c1ccccc1'; const P0 = fresh(smi);
    const parts = ['Friedel–Crafts acylation with succinic anhydride (AlCl₃) → 4-oxo-4-phenylbutanoic acid'].concat(red ? ['Clemmensen (Zn–Hg/HCl) → 4-phenylbutanoic acid'] : []).concat(red && cyc ? ['conc. H₂SO₄: intramolecular acylation → α-tetralone'] : []);
    return out({ title: red && cyc ? 'Haworth synthesis → α-tetralone' : red ? 'Haworth: → 4-phenylbutanoic acid' : 'Benzene + succinic anhydride (AlCl₃) → β-benzoylpropionic acid', product: P0, steps: [step(m, P0, parts.join('; then ') + '.', ['AlCl₃'])], why: ['Then Clemmensen again, Se/Δ (dehydrogenation) → naphthalene', 'Each step is a reaction you already know: FC acylation, Clemmensen, FC acylation again'], tags: ['SEAr'] }); } });
})(typeof window !== 'undefined' ? window : globalThis);
