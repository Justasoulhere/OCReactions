(function (root) {
'use strict';
'use strict';
const C = root.Chem, E = root.EngineCore; if (!C || !E) return;
const { BL, el, hTot } = E; const RP = () => root.RxParse;
const fresh = smi => { try { const m = C.fromSmilesMB(smi, BL); return m && m.atoms.size ? m : null; } catch (e) { return null; } };
const step = (L, R, text, cond) => ({ L: L.clone(), arrows: [], R: R.clone(), text, cond: cond || [] });
const out = o => Object.assign({ steps: [], minor: [], extra: [], why: [], db: [] }, o);
const M = (m, sm) => { try { return C.matchAll(m, sm, {}); } catch (e) { return []; } };
const RXC = new Map();
function rxn(sma, smis, max) {
  const RD = C.RD; let rx = RXC.get(sma); if (rx === undefined) { try { rx = RD.get_rxn(sma); } catch (e) { rx = null; } RXC.set(sma, rx); } if (!rx) return [];
  const L = new RD.MolList(); const ms = smis.map(s => RD.get_mol(s)); if (ms.some(x => !x)) { ms.forEach(x => x && x.delete()); return []; } ms.forEach(x => L.append(x));
  const res = []; let o = null; try { o = rx.run_reactants(L, max || 24); } catch (e) { o = null; }
  if (o) for (let i = 0; i < o.size(); i++) { const ml = o.get(i); const ps = []; let bad = false;
    for (let k = 0; k < ml.size(); k++) { const x = ml.at(k); let s = null; try { const q = RD.get_mol(x.get_smiles()); if (q) { s = q.get_smiles(); q.delete(); } } catch (e) { s = null; } x.delete && x.delete(); if (!s) bad = true; ps.push(s); }
    if (!bad && !res.some(r => r.join('.') === ps.join('.'))) res.push(ps); }
  ms.forEach(x => x.delete()); return res;
}
const canonS = s => { try { const q = C.RD.get_mol(s); if (!q) return null; const k = q.get_smiles(); q.delete(); return k; } catch (e) { return null; } };
const smiOf = m => { try { return C.canon(m); } catch (e) { return null; } };
const isoOf = m => { try { return (root.Stereo && root.Stereo.iso(m)) || C.canon(m); } catch (e) { return smiOf(m); } };
const nC = s => (String(s).replace(/Cl|Ca|Cu|Cr|Cs|Cd/g, '').match(/[Cc]/g) || []).length;
const add = r => E.REC.push(r), front = r => E.REC.unshift(r);
const byC = arr => arr.slice().sort((a, b) => nC(b) - nC(a) || b.length - a.length);
function finish(m, prodSmi, o) { const R = fresh(prodSmi); if (!R) return null; return out(Object.assign({}, o, { product: R, steps: [step(m, R, o.text || o.title, o.cond)], extra: (o.extra || []).map(fresh).filter(Boolean), minor: (o.minor || []).map(fresh).filter(Boolean) })); }
const none = (title, why) => out({ title, product: null, none: why });
function applyAll(smi, T, lim) { let cur = smi; const extra = [], labels = [];
  for (const [sma, label] of T) for (let g = 0; g < (lim || 6); g++) { const p = rxn(sma, [cur]); if (!p.length) break; const frs = byC(p[0].join('.').split('.')); if (frs[0] === cur) break; cur = frs[0]; extra.push(...frs.slice(1)); if (!labels.includes(label)) labels.push(label); }
  return { cur, extra, labels }; }
const ctxWith = (c, add2, flags) => { const st = Object.assign({}, c.st, { sp: new Map(c.st.sp), flags: new Set([...c.st.flags, ...(flags || [])]) }); for (const [k, v] of Object.entries(add2 || {})) { if (v === null) st.sp.delete(k); else st.sp.set(k, Object.assign({ sp: k }, v)); } return RP().stageCtx(st); };
const recipe = id => E.REC.find(x => x.id === id);
function wrap(id, fn) { const r = recipe(id); if (!r) return; const run0 = r.run; r.run = function (m, c, all, k) { let o; try { o = fn.call(this, m, c, all, k, run0.bind(this)); } catch (e) { o = undefined; } return o === undefined ? run0.call(this, m, c, all, k) : o; }; }
function retest(id, fn) { const r = recipe(id); if (!r) return; const t0 = r.test; r.test = function (c, all, k) { return fn(c, all, k, (cc, aa, kk) => { try { return t0(cc, aa, kk); } catch (e) { return false; } }); }; }
const frags = s => String(s || '').split('.');
const mainFrag = s => byC(frags(s))[0];
const raw = c => String(c.raw || '').toLowerCase();
const isAlkoxide = c => c.any(['NaOEt', 'NaOMe', 'NaOR', 'KOtBu']);
const alkoxR = c => c.has('NaOEt') ? 'CC' : c.has('NaOMe') ? 'C' : c.has('KOtBu') ? 'C(C)(C)C' : ((c.get('NaOR') || {}).smi || 'CC');
const CARB = { HCHO: 'C=O', CH3CHO: 'CC=O', PhCHO: 'O=Cc1ccccc1', acetone: 'CC(C)=O', acrolein: 'C=CC=O' };
function partnerCarbonyl(c) { for (const [k, s] of Object.entries(CARB)) if (c.has(k)) return { k, smi: (c.get(k) || {}).smi || s }; const p = c.get('partner'); const ps = Array.isArray(p) ? p : p ? [p] : []; for (const x of ps) if (x.smi && M(fresh(x.smi), '[CX3H1,CX3H2]=O').length) return { k: 'partner', smi: x.smi }; return null; }



/* ===== r17 recipes (run 35): general rules from the blind sets B (SKM Blue Book JEE Advanced) and C (reagent x functional group) ===== */
const heat = c => !!(c.heat || c.flag('warm') || c.flag('reflux') || (c.T != null && c.T >= 323) || /\bwarm|reflux|boil|\bheat|Δ/i.test(c.raw || ''));
const OXID = ['KMnO4', 'K2Cr2O7', 'OsO4', 'O3', 'HIO4', 'mCPBA', 'Tollens', 'Fehling', 'PCC'];
const hasOx = c => c.any(OXID);
const sugarLike = m => M(m, '[CX4][OX2H1]').length >= 3 && (M(m, '[CX3]=O').length > 0 || M(m, '[CX4;R]([OX2;R])[OX2H1]').length > 0);
const nCn = (m, id) => m.nbrs(id).filter(n => el(m, n) === 'C').length;
const partners = c => [c.get('partner'), c.get('org'), c.get('RX')].flat().filter(x => x && x.smi);

/* engine_r18.js: run 37. Fixes found by sweeping every reagent over every functional group (prof/), so the Reagents view can show a full JEE profile. */
const R18 = {};
const ringSubsHetero = m => M(m, 'c[!#6;!#1]').length > 0 || M(m, 'c[CX3]=O').length > 0 || M(m, 'c[CX2]#N').length > 0;
const isTert = m => M(m, '[CX4H0]([#6])([#6])([#6])[Cl,Br,I]').length > 0;
const anhydrousROH = c => c.has('ROH') && !c.has('H2O') && !/h2o|water|\baq/i.test(c.raw || '');

/* isocyanides are reduced to N-methyl 2° amines (LiAlH4, H2/Ni, Na/EtOH) */
front({ id: 'isocy_red', name: 'Isocyanide reduction → 2° amine', fam: 'Amine', canon: 'LiAlH₄ or H₂/Ni', db: [],
  test: c => c.has('LiAlH4') || (c.has('H2') && c.has('Pd') && !c.has('Lindlar')) || (c.has('Na') && c.has('ROH')),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[#6][N+]#[C-]').length) return null; const p = rxn('[#6:1][N+:2]#[C-:3]>>[*:1][N+0:2][C+0:3]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Isocyanide + 4[H] → N-methyl 2° amine (R–NH–CH₃)', text: 'The terminal isocyanide carbon becomes a CH₃ on nitrogen: R–N≡C → R–NH–CH₃. (A nitrile R–C≡N gives the 1° amine R–CH₂NH₂ instead.)', tags: ['reduction'] }); } });
/* NaBH4 does not reduce oximes / hydrazones; a poisoned catalyst (Lindlar) leaves C≡N and C=N alone */
wrap('cn_reduce_h2', function (m, c) { if (c.has('Lindlar')) return null; if (c.has('NaBH4') && M(m, '[CX3]=N[OX2,NX3]').length) return null; return undefined; });
wrap('nitrile_rx', function (m, c) { if (c.has('Lindlar') || c.has('Na/NH3') || c.has('HNO2')) return null; if (c.has('Na') && !c.has('ROH') && !c.any(['H2O', 'H3O+', 'KOH', 'aqKOH', 'LiAlH4', 'H2'])) return null; if (anhydrousROH(c) && c.has('H3O+') && !c.any(['LiAlH4', 'H2', 'Na'])) return null; return undefined; });
/* Na / liq. NH3 is a dissolving-metal reduction: NH3 is only the solvent */
for (const id of ['rx_nu', 'epox_open', 'imine', 'acyl_sub', 'ald_nh3']) wrap(id, function (m, c) { if (c.has('Na/NH3')) return null; return undefined; });
wrap('birch', function (m, c) { if (M(m, 'c[N+](=O)[O-]').length || M(m, 'c[N+]#N').length || M(m, 'c[OX2H1]').length || M(m, 'c[CX3]=O').length || M(m, 'c[Cl,Br,I]').length) return null; return undefined; });
front({ id: 'arx_na_nh3', name: 'Aryl halide + Na / liq. NH₃ → arene', fam: 'Aromatic', canon: 'Na, liq. NH₃', db: [], test: c => c.has('Na/NH3'),
  run(m) { const s = smiOf(m); if (!s || !M(m, 'c[Cl,Br,I]').length) return null; const r = applyAll(s, [['[c:1][Cl,Br,I]>>[c:1]', 'C–X → C–H']], 4);
    return finish(m, r.cur, { title: 'Na / liq. NH₃: Ar–X → Ar–H (the ring itself is not reduced first)', text: 'Electron transfer from Na breaks the C–X bond (aryl radical, then anion); NH₃ protonates it.', tags: ['reduction'] }); } });
wrap('acetylide', function (m, c) { if (c.has('ROH') && !c.has('NaNH2')) return null; return undefined; });
/* Zn in aqueous NH4Cl / HCl is a reducing agent, not a coupling metal; Zn-dust distillation needs heat; Ullmann needs Ar–I */
wrap('zn_cu_coupling', function (m, c) { if (c.any(['NH4Cl', 'H2O', 'HCl', 'H3O+'])) return null; if (c.has('Cu') && !M(m, 'c[I]').length && !M(m, '[CX4][Cl,Br,I]').length) return null; return undefined; });
wrap('phenol_rx', function (m, c) { if (c.has('Zn') && !heat(c) && !c.any(['KOH', 'aqKOH', 'CHCl3', 'CO2']) && (c.any(['HCl', 'NH4Cl', 'H2O', 'H3O+', 'AcOH', 'NaOH']) || !/dust/i.test(c.raw || ''))) return null; return undefined; });
/* NaHSO3 after OsO4 / KMnO4 is the reductive work-up, not a carbonyl reagent */
wrap('bisulfite', function (m, c, all, k) { if (k > 0 && all.slice(0, k).some(s => s && s.any(['OsO4', 'KMnO4', 'O3']))) return null; return undefined; });
/* ozonolysis of an alkyne: 1,2-dicarbonyl (reductive work-up) or acids (water / H2O2) */
front({ id: 'alkyne_o3', name: 'Ozonolysis of an alkyne', fam: 'Alkyne', canon: 'O₃ then Zn/H₂O', db: [], eats: s => s.any(['Me2S', 'H2O', 'H2O2', 'Zn']),
  test: c => c.has('O3'),
  run(m, c, all, k) { const s = smiOf(m); if (!s || !M(m, 'C#C').length || M(m, 'C=C').length) return null; const nx = all && all[k + 1]; const red = nx ? nx.any(['Me2S', 'Zn']) : true; const ox = nx && nx.any(['H2O2']);
    const p = ox || !red ? rxn('[C:1]#[C:2]>>[C:1](=O)O.[C:2](=O)O', [s]) : rxn('[C:1]#[C:2]>>[C:1](=O)[C:2]=O', [s]); if (!p.length) return null; const fr = p[0].join('.').split('.');
    return finish(m, byC(fr)[0], { title: ox || !red ? 'Alkyne + O₃, then H₂O / H₂O₂ → carboxylic acids (HC≡ end → CO₂)' : 'Alkyne + O₃, then Zn / H₂O → 1,2-dicarbonyl compound', text: 'The ozonide of a C≡C keeps both carbons joined as a 1,2-dicarbonyl under reductive work-up; oxidative work-up cleaves it to two –COOH groups.', extra: byC(fr).slice(1), tags: ['oxidation'] }); } });
wrap('ozone', function (m) { if (M(m, 'C#C').length && !M(m, 'C=C').length) return null; return undefined; });
wrap('arene_o3', function (m) { if (ringSubsHetero(m)) return null; return undefined; });
/* diazonium salts are not EAS substrates (they couple / substitute instead) */
for (const id of ['eas', 'br2_act']) wrap(id, function (m) { if (M(m, '[N+]#N').length) return null; return undefined; });
/* aniline (any Br2) and phenol in water / base are so activated that every free o/p position is brominated */
front({ id: 'tribromo', name: 'Aniline / phenoxide + Br₂ → 2,4,6-tribromo', fam: 'Aromatic', canon: 'Br₂', db: [],
  test: c => c.any(['Br2', 'HOBr', 'Br2/NaOH']) && !c.any(['AcOH', 'Ac2O']),
  run(m, c, all, k) { const an = M(m, 'c[NX3H2]').length > 0; const ph = M(m, 'c[OX2H1,OX1-]').length > 0 && (c.any(['HOBr', 'KOH', 'aqKOH', 'Br2/NaOH', 'H2O']) || /water|aq/i.test(c.raw || '')); if (!an && !ph) return null; if (c.has('HOBr') && !c.any(['KOH', 'aqKOH', 'Br2/NaOH'])) return null;
    const r = recipe('eas'); if (!r) return null; let o = null; try { o = r.run(m, ctxWith(c, { HOBr: {}, Br2: null, FeX3: null, CCl4: null, 'Br2/NaOH': null, KOH: null, aqKOH: null }), all, k); } catch (e) { o = null; } if (!o || !o.product) return null;
    return Object.assign({}, o, { title: an ? 'Aniline + Br₂ (any solvent): 2,4,6-tribromoaniline at once (white ppt). Protect as acetanilide for mono-bromination.' : 'Phenoxide + Br₂ (aqueous base): every free o/p position → 2,4,6-tribromophenol' }); } });
/* radical halogenation: not on N–H amines, isocyanides or oximes; Cl2 adds to C≡C and dienes */
wrap('radx', function (m, c) { if (M(m, '[NX3;!H0]').length || M(m, '[N+]#[C-]').length || M(m, '[CX3]=N[OX2H1]').length) return null; if (c.has('Cl2') && (M(m, 'C#C').length || M(m, 'C=CC=C').length)) return null; return undefined; });
/* nitrating mixture: alcohols give nitrate esters (glyceryl trinitrate); aliphatic aldehydes are oxidised */
front({ id: 'alc_nitrate', name: 'Alcohol + HNO₃ / H₂SO₄ → alkyl nitrate', fam: 'Alcohol', canon: 'conc. HNO₃, conc. H₂SO₄', db: [],
  test: c => c.has('HNO3') && c.has('H2SO4') && !(c.T != null && c.T > 400),
  run(m) { const s = smiOf(m); if (!s || M(m, 'c').length || !M(m, '[CX4][OX2H1]').length) return null; const r = applyAll(s, [['[CX4:1][OX2H1:2]>>[C:1][O:2][N+](=O)[O-]', 'R–OH → R–O–NO₂']], 4);
    return finish(m, r.cur, { title: 'R–OH + HNO₃ (H₂SO₄) → alkyl nitrate R–O–NO₂ (glycerol → nitroglycerine)', text: 'The alcohol O attacks NO₂⁺ (an O-nitration, an ester of nitric acid); no C–N bond forms.', tags: ['ester'] }); } });
front({ id: 'hno3_cho', name: 'HNO₃ oxidises an aliphatic aldehyde', fam: 'Aldehyde', canon: 'conc. HNO₃', db: [], test: c => c.has('HNO3'),
  run(m) { const s = smiOf(m); if (!s || M(m, 'c').length || !M(m, '[CX3H1](=O)[#6]').length) return null; const p = rxn('[CX3H1:1](=[O:2])[#6:3]>>[C:1](=[O:2])(O)[*:3]', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Conc. HNO₃ oxidises –CHO → –COOH', text: 'Nitric acid is a strong oxidant: the aldehyde C–H becomes C–OH.', tags: ['oxidation'] }); } });
wrap('acid_aldol', function (m, c) { if (c.has('HNO3')) return null; if (!(c.flag('conc') || c.flag('dry') || /conc|dry\s*hcl|hcl\s*gas/i.test(c.raw || ''))) return null; return undefined; });
/* Beckmann needs PCl5 / SOCl2 / P2O5 / conc. H2SO4; dilute acid simply hydrolyses the oxime back */
wrap('beckmann', function (m, c) { if (c.any(['PCl5', 'SOCl2', 'P2O5', 'PPA', 'POCl3'])) return undefined; if (c.has('H2SO4') && !c.flag('dil') && !c.flag('aq') && !/dil|aq/i.test(c.raw || '')) return undefined; return null; });
/* plain "H2SO4" in JEE means conc.; only dil./aq. acid hydrolyses the oxime back */
front({ id: 'oxime_hyd', name: 'Oxime / hydrazone + dilute acid → carbonyl', fam: 'Imine', canon: 'H₃O⁺', db: [],
  test: c => (c.any(['H3O+']) || (c.any(['H2SO4', 'HCl']) && (c.flag('dil') || c.flag('aq')))) && !c.flag('conc') && !c.any(['PCl5', 'SOCl2', 'P2O5', 'LiAlH4', 'H2', 'NaBH4', 'Na']),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[CX3]=N[OX2H1]').length || anhydrousROH(c)) return null; const p = rxn('[C:1]=[N][OX2H1]>>[C:1]=O', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Oxime + dilute acid, warm → carbonyl compound + NH₂OH (hydrolysis, the reverse of oxime formation)', text: 'Water adds to the protonated C=N; NH₂OH leaves and the C=O returns. (Beckmann rearrangement needs PCl₅ / conc. H₂SO₄.)', extra: ['NO'], tags: ['hydrolysis'] }); } });
wrap('diazo_rx', function (m, c) { if (c.has('KOtBu')) return null; return undefined; });
wrap('acyloin', function (m) { if (M(m, '[CX3](=O)O[CX3]=O').length) return null; return undefined; });
/* SN2 limits at a 3° carbon: no Finkelstein; NH3 / amines eliminate instead */
wrap('rx_nu', function (m, c) { if (!isTert(m)) return undefined; if (c.has('NaI')) return none('No Finkelstein at a 3° carbon', 'SN2 needs back-side attack, and a 3° carbon is too crowded. NaI / acetone works for 1° > 2° halides.');
  if (c.any(['NH3', 'RNH2', 'R2NH']) && !c.has('Na/NH3')) { const s = smiOf(m); const p = s && rxn('[CX4;!H0:1][CX4H0:2][Cl,Br,I]>>[C:1]=[C:2]', [s]); if (p && p.length) return finish(m, mainFrag(p[0].join('.')), { title: '3° halide + NH₃: elimination (E2) → alkene, not the amine', text: 'NH₃ cannot reach the crowded 3° carbon (no SN2); acting as a base it removes a β-H.', tags: ['elim'] }); }
  return undefined; });
wrap('gilman', function (m) { if (isTert(m)) return none('Corey–House fails for 3° halides', 'R₂CuLi couples with methyl, 1° and 2° halides (SN2-like). A 3° halide eliminates instead.'); return undefined; });
/* PCl5 converts a C=O into a gem-dichloride */
front({ id: 'co_pcl5', name: 'Aldehyde / ketone + PCl₅ → gem-dichloride', fam: 'Carbonyl', canon: 'PCl₅', db: [], test: c => c.has('PCl5'),
  run(m) { const s = smiOf(m); if (!s || M(m, '[OX2H1]').length || !M(m, '[CX3;$([CX3H1](=O)[#6]),$([CX3](=O)([#6])[#6])]=O').length) return null; const p = rxn('[CX3;$([CX3H1](=O)[#6]),$([CX3](=O)([#6])[#6]):1]=[OX1]>>[C:1](Cl)Cl', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Aldehyde / ketone + PCl₅ → gem-dichloride (+ POCl₃)', text: 'The C=O oxygen is replaced by two Cl (CH₃CHO → CH₃CHCl₂).', extra: ['O=P(Cl)(Cl)Cl'], tags: ['sub'] }); } });
/* excess HI: both halves of an ether end up as iodides */
wrap('ether_hi', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || !c.has('HI') || !(c.flag('excess') || /excess|xs/i.test(c.raw || '')) || c.has('redP')) return o; if (M(m, 'c[OX2][CX4]').length) return o;
  const conv = x => { const s0 = smiOf(x); if (!s0) return x; const r = applyAll(s0, [['[CX4:1][OX2H1]>>[C:1]I', 'R–OH → R–I']], 4); return fresh(r.cur) || x; };
  return Object.assign({}, o, { product: conv(o.product), extra: (o.extra || []).map(conv), title: 'Excess HI, Δ: ether → two alkyl iodides (the alcohol formed first is converted too)' }); });
/* alkene + alcohol / H+ → ether (Markovnikov); transesterification */
front({ id: 'alkene_roh', name: 'Alkene + ROH, H⁺ → ether', fam: 'Alkene', canon: 'ROH, H⁺', db: [], test: c => anhydrousROH(c) && c.any(['H3O+', 'H2SO4']),
  run(m, c) { const q = M(m, '[C;!a:1]=[C;!a:2]'); if (!q.length || M(m, '[CX3]=O').length || M(m, '[OX2H1]').length) return null; const [a, b] = q[0]; const R = ((c.get('ROH') || {}).smi) || 'C'; const cp = nCn(m, a) >= nCn(m, b) ? a : b; const t = m.clone(); t.bondBetween(a, b).o = 1; for (const id of [a, b]) t.atoms.get(id).hFix = null;
    let last = E.addOn(t, cp, 'O'); const rs = fresh(R); const nR = rs ? [...rs.atoms.values()].filter(x => x.el === 'C').length : 1; for (let i = 0; i < Math.min(nR, 4); i++) last = E.addOn(t, last, 'C'); const ps = smiOf(t); if (!ps) return null;
    return finish(m, ps, { title: 'Alkene + alcohol, H⁺ → ether (Markovnikov: OR goes to the carbon that holds the carbocation)', text: 'H⁺ adds to give the more stable carbocation; the alcohol O attacks it and loses H⁺.', tags: ['add'] }); } });
front({ id: 'transester', name: 'Transesterification', fam: 'Ester', canon: 'R′OH, H⁺ or R′O⁻', db: [], test: c => c.has('ROH') && (c.any(['H3O+', 'H2SO4']) || isAlkoxide(c)) && !c.has('H2O'),
  run(m, c) { const s = smiOf(m); const R = ((c.get('ROH') || {}).smi) || 'C'; if (!s || !M(m, '[CX3](=O)[OX2][CX4]').length || M(m, '[CX3](=O)[OX2][CX3]=O').length) return null; const p = rxn('[#6:5][CX3:1](=[O:2])[OX2:3][CX4:4]>>[*:5][C:1](=[O:2])O' + R + '.[O:3][*:4]', [s]); if (!p.length) return null; const fr = p[0].join('.').split('.'); const main = fr[0]; if (canonS(main) === canonS(s)) return null;
    return finish(m, main, { title: 'Transesterification: the ester swaps its alkoxy group for the solvent alcohol', text: 'Acid (or alkoxide) catalysed acyl substitution: R′OH attacks the C=O, the old OR leaves as alcohol. Excess alcohol drives it.', extra: fr.slice(1), tags: ['acyl'] }); } });
for (const id of ['amide_rx', 'isocyanide_hyd']) wrap(id, function (m, c) { if (anhydrousROH(c) && !c.any(['LiAlH4', 'H2', 'Br2/NaOH', 'KOH', 'aqKOH', 'P2O5'])) return null; return undefined; });
wrap('hydration', function (m, c) { if (anhydrousROH(c)) return null; return undefined; });
/* NaNO2 / HCl: the dilute cold HCl is only there to make HNO2 */
for (const id of ['hx', 'diene_add', 'alc3_hcl', 'ester_hyd', 'alc_hx']) wrap(id, function (m, c) { if (c.has('HNO2')) return null; return undefined; });
/* HI / red P on a non-aromatic O-compound: remove every O (rings stay closed) */
wrap('hi_redp', function (m, c) { if (M(m, 'c').length || !M(m, '[#8]').length || [...m.atoms.values()].some(a => !['C', 'O', 'H'].includes(a.el))) return undefined; const t = m.clone(); for (const a of [...t.atoms.values()]) if (a.el === 'O') t.removeAtom(a.id); for (const b of t.bonds.values()) b.o = 1; for (const a of t.atoms.values()) { a.hFix = null; a.chg = 0; }
  const fr = String(smiOf(t) || '').split('.').filter(Boolean); if (!fr.length) return undefined;
  return finish(m, byC(fr)[0], { title: 'HI / red P, Δ: every C–O becomes C–H (C=O → CH₂, ethers and esters are cut at the O)', text: 'Prolonged heating with HI and red phosphorus removes every oxygen; HI reduces each C–I to C–H and P recycles the I₂.', extra: byC(fr).slice(1), tags: ['reduction'] }); });
/* DIBAL-H also reduces aldehydes and ketones (to alcohols); BH3 reduces amides and nitriles to amines */
front({ id: 'dibal_co', name: 'DIBAL-H + aldehyde / ketone', fam: 'Reduction', canon: 'DIBAL-H', db: [], test: c => c.has('DIBAL'),
  run(m) { const s = smiOf(m); if (!s || M(m, '[CX3](=O)[OX2,Cl,NX3]').length || M(m, 'C#N').length || !M(m, '[CX3;$([CX3H1](=O)[#6]),$([CX3](=O)([#6])[#6])]=O').length) return null; const r = applyAll(s, [['[CX3;$([CX3H1](=O)[#6]),$([CX3](=O)([#6])[#6]):1]=[OX1:2]>>[C:1][O:2]', 'C=O → CH–OH']], 4);
    return finish(m, r.cur, { title: 'DIBAL-H: aldehyde / ketone → alcohol (it stops at the aldehyde only for esters and nitriles at −78 °C)', text: 'Al–H delivers hydride to the C=O carbon; work-up protonates the alkoxide.', tags: ['reduction'] }); } });
front({ id: 'bh3_n', name: 'BH₃: amide / nitrile → amine', fam: 'Reduction', canon: 'B₂H₆ then H₃O⁺', db: [], test: c => c.has('BH3') && !c.has('H2O2'),
  run(m) { const s = smiOf(m); if (!s || M(m, 'C=C').length || M(m, 'C#C').length) return null; const r = applyAll(s, [['[CX3:1](=[OX1])[NX3:2]>>[C:1][N:2]', 'amide C=O → CH₂'], ['[C:1]#[N:2]>>[C:1][N:2]', 'C≡N → CH₂NH₂']], 4); if (r.cur === s) return null;
    return finish(m, r.cur, { title: 'BH₃ (diborane): ' + r.labels.join(', ') + ' → amine', text: 'Borane is an electrophilic hydride: it coordinates the electron-rich amide O / nitrile N and delivers H to carbon.', tags: ['reduction'] }); } });

wrap('cn_reduce', function (m, c) { if (c.has('Lindlar')) return null; if (c.has('NaBH4') && M(m, '[CX3]=N[OX2,NX3]').length) return null; return undefined; });
wrap('arene_h2', function (m) { if (M(m, 'c[Cl,Br,I]').length || M(m, '[N+]#N').length) return null; return undefined; });
wrap('bisulfite', function (m, c) { if (c.any(['OsO4', 'KMnO4', 'O3'])) return null; return undefined; });
/* Cl2 / Br2 add to C≡C and dienes even in light (radical substitution needs an sp3 C–H and no π bond to attack) */
front({ id: 'x2_light_add', name: 'X₂ (light) + alkyne / diene → addition', fam: 'Alkene', canon: 'Cl₂', db: [], test: c => c.any(['Cl2', 'Br2']) && c.light && !c.has('FeX3'),
  run(m, c, all, k) { if (!M(m, 'C#C').length && !M(m, 'C=CC=C').length) return null; for (const id of ['diene_add', 'x2']) { const r = recipe(id); if (!r) continue; let o = null; try { o = r.run(m, ctxWith(c, {}, []), all, k); } catch (e) { o = null; } if (!o) { const c2 = ctxWith(c, {}, []); c2.light = false; try { o = r.run(m, c2, all, k); } catch (e) { o = null; } } if (o && o.product) return o; } return null; } });
/* diazomethane: RCOOH → methyl ester (+ N2); with an acid chloride → diazoketone (Arndt–Eistert) is handled elsewhere */
front({ id: 'ch2n2_acid', name: 'RCOOH + CH₂N₂ → methyl ester', fam: 'Acid', canon: 'CH₂N₂', db: [], test: c => (c.has('CH2N2') || /ch2n2|diazomethane/i.test(c.raw || '')) && !c.any(['SOCl2', 'PCl5', 'PCl3', 'Ag2O']),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX3](=O)[OX2H1]').length) return null; const r = applyAll(s, [['[CX3:1](=[O:2])[OX2H1:3]>>[C:1](=[O:2])[O:3]C', 'COOH → COOCH₃']], 3);
    return finish(m, r.cur, { title: 'R–COOH + CH₂N₂ → R–COOCH₃ + N₂↑ (clean methyl ester)', text: 'The acid protonates diazomethane (CH₃–N₂⁺); the carboxylate attacks the CH₃ in SN2 and N₂ leaves.', extra: ['N#N'], tags: ['ester'] }); } });

/* MnO2 oxidises only allylic / benzylic alcohols */
for (const id of ['alc_ox_mild', 'ox_all_oh']) wrap(id, function (m, c) { if (!c.has('MnO2') && !/mno2/i.test(c.raw || '')) return undefined; if (M(m, '[OX2H1][CX4;$(C[CX3]=[CX3]),$(Cc)]').length) return undefined; return none('MnO₂ needs an allylic or benzylic –OH', 'Active MnO₂ oxidises only allylic / benzylic alcohols (to the α,β-unsaturated or aryl carbonyl); a saturated alcohol is untouched. Use PCC or K₂Cr₂O₇ for those.'); });
wrap('solvolysis', function (m, c) { if (c.has('Al(OiPr)3') || /al\s*\(\s*o-?i-?pr\)\s*3|isopropoxide/i.test(c.raw || '')) return null; return undefined; });
/* Hofmann elimination needs a quaternary ammonium salt (or CH3I first) */
wrap('hofmann_elim', function (m, c) { if (M(m, '[NX4+]').length || c.any(['RX']) || /ch3i|methyl iodide/i.test(c.raw || '')) return undefined; return null; });
/* no Grignard from a molecule that already carries C=O / acidic H / acyl halide */
wrap('grignard_make', function (m) { if (M(m, '[CX3]=O').length || M(m, '[OX2H1,NX3;!H0]').length) return null; return undefined; });
/* an acid does not react with NaN3 alone (Schmidt needs HN3 + conc. H2SO4) */
wrap('curtius_lossen', function (m, c, all, k) { if (M(m, '[CX3](=O)[OX2H1]').length && !M(m, '[CX3](=O)[Cl,Br]').length && !(all || []).slice(0, k + 1).some(s => s && s.any(['SOCl2', 'PCl5', 'PCl3', 'NH2OH']))) return null; return undefined; });
/* acyl halides precipitate AgCl at once with AgNO3 */
wrap('agno3_rx', function (m) { if (!M(m, '[CX3](=O)[Cl,Br,I]').length) return undefined; return Object.assign(none('Positive: AgX at once', 'An acyl halide is attacked instantly by the solvent (water / alcohol); the X⁻ released precipitates with Ag⁺ at once.'), { positive: true }); });

/* ---- run 38 ---- */
/* an amine / phenol / alcohol partner is the nucleophile in Schotten–Baumann: the NaOH only mops up HCl */
wrap('acylx_oh', function (m, c) { if (c.any(['phenol', 'aniline', 'RNH2', 'R2NH', 'ROH'])) return null; return undefined; });
/* nitrating mixture / anhydrous acid has no water to hydrolyse an ester; Baeyer's reagent (cold, dilute) leaves esters and epoxides alone */
wrap('ester_hyd', function (m, c) { if (c.has('HNO3')) return null; if (c.has('KMnO4') && (c.flag('cold') || c.flag('mild'))) return null; return undefined; });
wrap('epox_oh', function (m, c) { if (c.has('KMnO4')) return null; return undefined; });
/* conc. HNO3 oxidises BOTH ends of an aldose (saccharic / glucaric acid) */
wrap('hno3_cho', function (m) { if (M(m, '[CX4][OX2H1]').length >= 3) return null; return undefined; });
/* E1cB dehydration must not leave an enol (that is the α-OH of a sugar, not an aldol) */
wrap('aldol_dehyd', function (m, c, all, k, run0) { let o; try { o = run0(m, c, all, k); } catch (e) { o = null; } if (!o || !o.product) return o === undefined ? null : o; if (M(o.product, '[OX2H1][CX3]=[CX3]').length) return null; return o; });
/* cold dilute KMnO4 on an internal alkyne → 1,2-diketone (no C–C cleavage) */
front({ id: 'alkyne_mild', name: 'Internal alkyne + cold dil. KMnO₄ → 1,2-diketone', fam: 'Alkyne', canon: 'cold dil. KMnO₄', db: [],
  test: c => c.has('KMnO4') && (c.flag('cold') || c.flag('mild')) && !heat(c),
  run(m) { const s = smiOf(m); if (!s || M(m, 'C=C').length || !M(m, '[#6][CX2]#[CX2][#6]').length) return null; const r = applyAll(s, [['[#6:3][C:1]#[C:2][#6:4]>>[*:3][C:1](=O)[C:2](=O)[*:4]', 'C≡C → CO–CO']], 2); if (!r || r.cur === s) return null;
    return finish(m, r.cur, { title: 'Internal alkyne + cold, dilute KMnO₄ → 1,2-diketone (the C–C bond survives)', text: 'Under cold, near-neutral conditions MnO₄⁻ adds twice across the C≡C (syn, cyclic manganate esters); the tetra-ol loses water to the 1,2-dione. Hot / acidic KMnO₄ cleaves it to acids.', tags: ['oxidation'] }); } });
/* phenol / aniline + Br2 water: an SO3H or COOH at an o/p position is replaced too (ipso attack, SO3 / CO2 lost) */
wrap('tribromo', function (m, c, all, k, run0) { const s = smiOf(m); if (!s) return undefined; const G = '[$([SX4](=O)(=O)[OX2H1,OX1-]),$([CX3](=O)[OX2H1,OX1-])]';
  if (!M(m, '[OX2H1,OX1-,NX3H2]c:c-' + G).length && !M(m, '[OX2H1,OX1-,NX3H2]c:c:c:c-' + G).length) return undefined;
  const T = [['[OX2H1,OX1-,NX3H2:9][c:1]:[c:2]-[SX4](=O)(=O)[OX2H1,OX1-]>>[*:9][c:1]:[c:2]', 'ipso -SO3H'], ['[OX2H1,OX1-,NX3H2:9][c:1]:[c:2]:[c:3]:[c:4]-[SX4](=O)(=O)[OX2H1,OX1-]>>[*:9][c:1]:[c:2]:[c:3]:[c:4]', 'ipso -SO3H'],
    ['[OX2H1,OX1-,NX3H2:9][c:1]:[c:2]-[CX3](=O)[OX2H1,OX1-]>>[*:9][c:1]:[c:2]', 'ipso -COOH'], ['[OX2H1,OX1-,NX3H2:9][c:1]:[c:2]:[c:3]:[c:4]-[CX3](=O)[OX2H1,OX1-]>>[*:9][c:1]:[c:2]:[c:3]:[c:4]', 'ipso -COOH']];
  const r = applyAll(s, T, 3); if (!r || r.cur === s) return undefined; const m2 = fresh(r.cur); if (!m2) return undefined; let o; try { o = run0(m2, c, all, k); } catch (e) { o = null; } if (!o || !o.product) return undefined;
  const lost = r.labels.some(x => /SO3H/.test(x)) ? 'SO₃H' : 'COOH'; return finish(m, smiOf(o.product), { title: o.title.replace(/:.*$/, '') + ': the –' + lost + ' at an o/p position is replaced by Br too (ipso substitution)', text: 'The ring is so activated that Br⁺ also attacks the carbon carrying –' + lost + '; that group leaves as ' + (lost === 'COOH' ? 'CO₂' : 'SO₃') + ' instead of H⁺, so the product is the 2,4,6-tribromo compound.', extra: [lost === 'COOH' ? 'O=C=O' : 'O=S(=O)=O'], tags: ['eas'] }); });

/* ---- run 38: curly arrows for polar steps that still morphed ---- */
const TK = (root.EngineR17 && root.EngineR17.tk) || null;
if (TK) {
const { Mech, spOf, atomWhere, sameMain, arrowless } = TK;
const fragsOf = mech => [mech.mainFrag()].concat(mech.others());
function useMech2(o, mech) { if (!mech || !mech.steps.length) return o; for (const ids of fragsOf(mech)) { let P = null; try { P = mech.product(ids); } catch (e) { P = null; } if (P && sameMain(P, o.product)) { const bare = mech.steps.filter(st => !st.arrows.length).length; if (!arrowless(o) && bare >= o.steps.filter(st => !(st.arrows || []).length).length) return o; return Object.assign({}, o, { steps: mech.steps, product: P }); } } return o; }
const someBare = o => o && o.steps && o.steps.some(st => !(st.arrows || []).length);
function arrow2(id, build) { wrap(id, function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product || o.none || !someBare(o)) return o; let mech = null; try { mech = build(m, c, o); } catch (e) { mech = null; } return mech ? useMech2(o, mech) : o; }); }
const addSp = (mech, def, near, pick, dist) => { const sp = spOf(def); if (!sp) return null; const a = atomWhere(sp, pick); const ids = mech.add(sp, near, { att: a, dist: (dist || 1.7) * BL }); return { ids, at: ids.get(a), sp }; };
const isO = a => a.el === 'O', isN = a => a.el === 'N', neg = a => (a.chg || 0) < 0;
/* H3O+ gives H+ to a lone pair on `at` */
function protOn(mech, at, text, cond) { const h = addSp(mech, '[H][OH2+]', at, a => a.el === 'H', 1.5); if (!h) return false; const H = h.at, O = h.ids.get(atomWhere(h.sp, isO)); mech.step([[[at], [at, H]], [[H, O], [O]]], text, cond || []); mech.dropFrag(O); return true; }
/* a base (water / OH-) takes H+ from `at`; to = atom that takes the pair as a π bond; also = extra arrow */
function deprot(mech, at, text, def, to, also) { const h = mech.H(at); if (h == null) return false; const b = addSp(mech, def || 'O', h, def ? neg : isO, 1.5); if (!b) return false; const A = [[[b.at], [b.at, h]], [[h, at], to != null ? [at, to] : [at]]]; if (also) A.push(also); mech.step(A, text); mech.dropFrag(b.at); return true; }
const BASE = c => c.any(['KOH', 'aqKOH', 'alcKOH', 'NaOEt', 'NaOMe', 'Ba(OH)2', 'Ca(OH)2']) || /naoh|koh|oh-|alkali|base/i.test(c.raw || '');
/* ---------- nucleophilic acyl substitution by water (ester, anhydride, amide) ---------- */
function acylHyd(m, c, o, pat, lgName, lgN) { const q = M(m, pat); if (q.length !== 1) return null; const [c1, o2, x] = q[0]; const mech = new Mech(m); const base = BASE(c) && !c.any(['H3O+', 'H2SO4', 'HCl']);
  if (base) { const nu = addSp(mech, '[OH-]', c1, isO); if (!nu) return null;
    mech.step([[[nu.at], [nu.at, c1]], [[c1, o2], [o2]]], 'OH⁻ adds to the carbonyl carbon; the C=O π pair moves onto O: a tetrahedral intermediate.', ['OH⁻']);
    mech.step([[[o2], [o2, c1]], [[c1, x], [x]]], 'The O⁻ pushes back down to re-form C=O and ' + lgName + ' leaves (addition–elimination).');
    const h = mech.H(nu.at); if (h == null) return mech; mech.step([[[x], [x, h]], [[h, nu.at], [nu.at]]], 'Proton transfer: the strongly basic leaving group takes H⁺ from the new –COOH: carboxylate. This step is irreversible, so base hydrolysis goes to completion.'); return mech; }
  if (!protOn(mech, o2, 'The C=O oxygen takes H⁺: the carbonyl carbon becomes much more electrophilic.', ['H₃O⁺'])) return null;
  const w = addSp(mech, '[H]O[H]', c1, isO); if (!w) return null;
  mech.step([[[w.at], [w.at, c1]], [[c1, o2], [o2]]], 'Water (weak nucleophile) attacks the activated carbonyl carbon; the π pair moves onto O.');
  deprot(mech, w.at, 'Water removes H⁺ from the attacking O: a neutral tetrahedral intermediate.');
  if (!protOn(mech, x, 'The ' + lgName + ' takes H⁺, turning it into a good (neutral) leaving group.')) return mech;
  mech.step([[[o2], [o2, c1]], [[c1, x], [x]]], 'An O–H lone pair re-forms C=O and ' + (lgN || 'the neutral leaving group') + ' leaves.');
  deprot(mech, o2, 'Water removes H⁺ from C=OH⁺: the carboxylic acid (H₃O⁺ is regenerated: acid is a catalyst).'); return mech; }
arrow2('ester_hyd', (m, c, o) => { if (/transest/i.test(o.title || '')) return null; return M(m, '[CX3](=[OX1])[OX2][CX3]=O').length ? acylHyd(m, c, o, '[CX3:1](=[OX1:2])[OX2:3][CX3]=O', 'carboxylate (RCOO⁻)', 'a carboxylic acid molecule') : acylHyd(m, c, o, '[CX3:1](=[OX1:2])[OX2:3][#6]', 'alkoxide (RO⁻)', 'the alcohol R′–OH'); });
arrow2('amide_rx', (m, c, o) => { if (!/hydrol/i.test(o.title || '')) return null; return acylHyd(m, c, o, '[CX3:1](=[OX1:2])[NX3:3]', 'amide ion (⁻NH₂)', 'NH₃ (as NH₄⁺ in acid)'); });
arrow2('acylx_hyd', (m, c, o) => acylHyd(m, c, o, '[CX3:1](=[OX1:2])[Cl,Br:3]', 'halide (Cl⁻)', 'HCl'));
/* ---------- nitrile hydrolysis: C≡N → amide (→ acid) ---------- */
arrow2('nitrile_rx', (m, c, o) => { const q = M(m, '[CX2:1]#[NX1:2]'); if (q.length !== 1) return null; const [c1, n] = q[0]; const toAmide = M(o.product, 'C(=O)[NX3H2]').length > 0, toAcid = M(o.product, '[CX3](=O)[OX2H1,OX1-]').length > 0; if (!toAmide && !toAcid) return null; const mech = new Mech(m); const base = BASE(c) && !c.any(['H3O+', 'H2SO4', 'HCl']);
  if (base) { const nu = addSp(mech, '[OH-]', c1, isO); if (!nu) return null; mech.step([[[nu.at], [nu.at, c1]], [[c1, n], [n]]], 'OH⁻ adds to the nitrile carbon; one π pair of C≡N moves onto N.', ['OH⁻']);
    const w = addSp(mech, '[H]O[H]', n, a => a.el === 'H', 1.5); const wh = w.ids.get(1), wo = w.ids.get(2); mech.step([[[n], [n, wh]], [[wh, wo], [wo]]], 'The N⁻ takes H⁺ from water: an imidic acid R–C(OH)=NH.'); mech.dropFrag(wo);
    const h = mech.H(nu.at); const b = addSp(mech, '[OH-]', h, isO, 1.5); mech.step([[[b.at], [b.at, h]], [[h, nu.at], [nu.at, c1]], [[c1, n], [n]]], 'Tautomerism: OH⁻ removes the O–H proton, the C=O forms and the C=N pair moves onto N.'); mech.dropFrag(b.at);
    const w2 = addSp(mech, '[H]O[H]', n, a => a.el === 'H', 1.5); mech.step([[[n], [n, w2.ids.get(1)]], [[w2.ids.get(1), w2.ids.get(2)], [w2.ids.get(2)]]], 'N⁻ takes H⁺ from water: the amide R–CONH₂.'); mech.dropFrag(w2.ids.get(2));
    if (toAmide && !toAcid) return mech;
    const nu2 = addSp(mech, '[OH-]', c1, isO); mech.step([[[nu2.at], [nu2.at, c1]], [[c1, nu.at], [nu.at]]], 'Second stage, amide hydrolysis: OH⁻ adds to the amide C=O (tetrahedral intermediate).');
    mech.step([[[nu.at], [nu.at, c1]], [[c1, n], [n]]], 'C=O re-forms and ⁻NH₂ leaves.'); const h2 = mech.H(nu2.at); mech.step([[[n], [n, h2]], [[h2, nu2.at], [nu2.at]]], '⁻NH₂ takes the acid proton: NH₃↑ and the carboxylate (acidify to get RCOOH).'); return mech; }
  if (!protOn(mech, n, 'The nitrile N takes H⁺: the C≡N⁺–H carbon is now electrophilic.', ['H₃O⁺'])) return null;
  const w = addSp(mech, '[H]O[H]', c1, isO); mech.step([[[w.at], [w.at, c1]], [[c1, n], [n]]], 'Water attacks the nitrile carbon; a π pair moves onto N.');
  deprot(mech, w.at, 'Water removes H⁺ from O: the imidic acid R–C(OH)=NH.');
  if (!protOn(mech, n, 'Tautomerism 1: the imine N takes H⁺.')) return mech;
  deprot(mech, w.at, 'Tautomerism 2: water removes the O–H proton; the C=O forms and the C=N⁺ pair moves onto N: the amide.', null, c1, [[c1, n], [n]]);
  if (toAmide && !toAcid) return mech;
  const o2 = w.at; if (!protOn(mech, o2, 'Amide hydrolysis: the amide C=O takes H⁺.')) return mech;
  const w2 = addSp(mech, '[H]O[H]', c1, isO); mech.step([[[w2.at], [w2.at, c1]], [[c1, o2], [o2]]], 'Water attacks the carbonyl carbon (tetrahedral intermediate).');
  deprot(mech, w2.at, 'Water removes H⁺ from the new O.'); if (!protOn(mech, n, 'The –NH₂ takes H⁺ (NH₃ is a good leaving group, NH₂⁻ is not).')) return mech;
  mech.step([[[o2], [o2, c1]], [[c1, n], [n]]], 'C=O re-forms and NH₃ leaves (in acid it ends up as NH₄⁺).'); deprot(mech, o2, 'Water removes H⁺: the carboxylic acid.'); return mech; });
/* ---------- metal hydrides: H⁻ to C=O, epoxide, acyl halide, nitrile ---------- */
const hydSp = (kind) => kind === 'B' ? { atoms: [{ el: 'B', chg: 0, h: 2 }, { el: 'H', x: -0.82 }], bonds: [[0, 1]] } : kind === 'Bm' ? { atoms: [{ el: 'B', chg: -1, h: 3 }, { el: 'H', x: -0.82 }], bonds: [[0, 1]] } : kind === 'D' ? { atoms: [{ el: 'Al', chg: 0, h: 2 }, { el: 'H', x: -0.82 }], bonds: [[0, 1]] } : { atoms: [{ el: 'Al', chg: -1, h: 3 }, { el: 'H', x: -0.82 }], bonds: [[0, 1]] };
function hydAdd(mech, site, kind, pi, text, cond) { const ids = mech.add(hydSp(kind), site, { dist: 1.75 * BL, att: 2 }); const b = ids.get(1), h = ids.get(2); const A = [[[b, h], [h, site]]]; if (pi) A.push(pi); mech.step(A, text, cond || []); mech.dropFrag(b); return h; }
const hName = k => k === 'B' ? 'BH₃' : k === 'Bm' ? 'BH₄⁻' : k === 'D' ? '(i-Bu)₂Al–H' : 'AlH₄⁻';
const hKind = c => c.has('LiAlH4') ? 'L' : c.has('NaBH4') ? 'Bm' : (c.has('DIBAL') || /dibal/i.test(c.raw || '')) ? 'D' : (c.any(['BH3', 'B2H6']) || /bh3|b2h6|borane|diborane/i.test(c.raw || '')) ? 'B' : null;
function hydCO(m, c) { const k = hKind(c); if (!k) return null; const q = M(m, '[CX3:1]=[OX1:2]'); if (q.length !== 1 || M(m, '[CX3](=O)[O,N,Cl,Br]').length) return null; const [c1, o2] = q[0]; const mech = new Mech(m);
  hydAdd(mech, c1, k, [[c1, o2], [o2]], k === 'B' ? 'BH₃ (a Lewis acid) coordinates the C=O oxygen and delivers H⁻ to the carbonyl carbon; the π pair moves onto O (O–B bond).' : 'H⁻ from ' + hName(k) + ' attacks the carbonyl carbon; the π pair moves onto O: an alkoxide.', [hName(k)]);
  protOn(mech, o2, 'Work-up: the alkoxide takes H⁺ from water / H₃O⁺: the alcohol.', ['H₂O']); return mech; }
for (const id of ['bh3_co', 'dibal_co']) arrow2(id, hydCO);
arrow2('hyd_red', (m, c, o) => { const k = hKind(c) || 'L'; const ep = M(m, '[CX4;r3:1]1[OX2;r3:2][CX4;r3:3]1'), ax = M(m, '[CX3:1](=[OX1:2])[Cl,Br:3]'); const mech = new Mech(m);
  if (ep.length === 1 && !ax.length) { const [a0, oo, b0] = ep[0]; const [a] = nCn(m, a0) <= nCn(m, b0) ? [a0, b0] : [b0, a0]; hydAdd(mech, a, k, [[a, oo], [oo]], 'H⁻ attacks the less hindered epoxide carbon from the back side (SN2); the strained C–O bond breaks.', [hName(k)]); protOn(mech, oo, 'Work-up: the alkoxide takes H⁺: the alcohol (OH stays on the more substituted carbon).', ['H₃O⁺']); return mech; }
  if (ax.length === 1 && !ep.length) { const [c1, o2, x] = ax[0]; hydAdd(mech, c1, k, [[c1, o2], [o2]], 'H⁻ adds to the acyl carbon; the π pair moves onto O.', [hName(k)]); mech.step([[[o2], [o2, c1]], [[c1, x], [x]]], 'C=O re-forms and Cl⁻ leaves: an aldehyde (more reactive than the acid chloride).');
    hydAdd(mech, c1, k, [[c1, o2], [o2]], 'A second H⁻ adds at once: the alkoxide of the 1° alcohol.'); protOn(mech, o2, 'Work-up: H₃O⁺ protonates the alkoxide: R–CH₂OH.', ['H₃O⁺']); return mech; }
  return null; });
/* ---------- diazonium: N₂ leaves (aryl cation), Nu captures; azo coupling ---------- */
arrow2('diazo_rx', (m, c, o) => { const q = M(m, '[c:1][N+:2]#[N:3]'); if (q.length !== 1) return null; const [ca, n1, n2] = q[0]; const mech = new Mech(m);
  const P = o.product; const ph = M(P, 'c[OX2H1]').length && !M(m, 'c[OX2H1]').length, ai = M(P, 'cI').length > M(m, 'cI').length, af = M(P, 'cF').length > M(m, 'cF').length;
  const azo = M(P, 'cN=Nc').length > 0;
  if (azo) return null;
  if (!ph && !ai && !af) return null;
  mech.step([[[ca, n1], [n1]]], 'On warming the C–N bond breaks heterolytically: N₂ (a superb leaving group) escapes and an aryl cation is left (the + sits in an sp² orbital in the ring plane).', [ph ? 'H₂O, warm' : ai ? 'KI' : 'HBF₄, Δ']); mech.dropFrag(n1);
  if (ph) { const w = addSp(mech, '[H]O[H]', ca, isO); mech.step([[[w.at], [w.at, ca]]], 'Water captures the aryl cation.'); deprot(mech, w.at, 'Water removes H⁺: phenol.'); return mech; }
  const nu = addSp(mech, ai ? '[I-]' : '[F-]', ca, neg); mech.step([[[nu.at], [nu.at, ca]]], ai ? 'I⁻ (soft, strong nucleophile) captures the aryl cation: iodoarene (no Cu needed).' : 'F⁻ from BF₄⁻ captures the aryl cation: fluoroarene (Balz–Schiemann; BF₃ escapes).'); return mech; });
/* ---------- CH₂N₂ + RCOOH ---------- */
arrow2('ch2n2_acid', (m) => { const q = M(m, '[CX3](=O)[OX2H1:1]'); if (q.length !== 1) return null; const oh = q[0][2]; const mech = new Mech(m); const h = mech.H(oh); if (h == null) return null;
  const d = addSp(mech, '[CH2-][N+]#N', h, a => a.el === 'C', 1.5); if (!d) return null; const cC = d.at, nP = d.ids.get(atomWhere(d.sp, a => a.el === 'N' && (a.chg || 0) > 0));
  mech.step([[[cC], [cC, h]], [[h, oh], [oh]]], 'The carbanion carbon of diazomethane (CH₂⁻–N⁺≡N form) takes the acidic proton: RCOO⁻ + CH₃–N₂⁺.', ['CH₂N₂']);
  mech.step([[[oh], [oh, cC]], [[cC, nP], [nP]]], 'SN2: the carboxylate O attacks the CH₃; N₂ (best leaving group known) escapes: the methyl ester.'); return mech; });
/* ---------- acid-catalysed addition of ROH / transesterification ---------- */
arrow2('alkene_roh', (m, c, o) => { const q = M(m, '[CX3:1]=[CX3:2]'); if (!q.length) return null; const v = c.get('ROH'); const R = (v && v.smi) || 'C'; const mech = new Mech(m); const r = E.protonateAlkene(mech, q[0][0], q[0][1], '[H][OH2+]'); if (!r) return null; mech.dropFrag(r.X);
  const al = addSp(mech, 'O' + R, r.cp, isO); if (!al) return null; mech.step([[[al.at], [al.at, r.cp]]], 'The alcohol O (lone pair) captures the carbocation: an oxonium ion.'); deprot(mech, al.at, 'Another ROH / water removes H⁺: the ether.'); return mech; });
arrow2('transester', (m, c, o) => { const q = M(m, '[CX3:1](=[OX1:2])[OX2:3][#6]'); if (q.length !== 1) return null; const [c1, o2, x] = q[0]; const v = c.get('ROH'); const R = (v && v.smi) || 'C'; const mech = new Mech(m);
  if (!protOn(mech, o2, 'The ester C=O takes H⁺ (acid catalyst).', ['H⁺'])) return null; const al = addSp(mech, 'O' + R, c1, isO); mech.step([[[al.at], [al.at, c1]], [[c1, o2], [o2]]], 'The new alcohol attacks the carbonyl carbon.');
  deprot(mech, al.at, 'Loss of H⁺ from the incoming O.'); protOn(mech, x, 'The old –OR′ takes H⁺.'); mech.step([[[o2], [o2, c1]], [[c1, x], [x]]], 'C=O re-forms; R′OH leaves (excess new alcohol drives this equilibrium).'); deprot(mech, o2, 'Loss of H⁺: the new ester.'); return mech; });
/* ---------- redox steps: highlight view, labelled ---------- */
for (const id of ['h2', 'hi_redp', 'kmno4hot', 'rx_reduce', 'nitro_red', 'na_etoh', 'birch', 'alc_na', 'arene_h2', 'alkyne_kmno4', 'kmno4_cho', 'isocy_red', 'tollens_first', 'tollens_plus', 'ox_all_oh', 'zn_cu_coupling', 'grignard_make', 'phenol_quinone', 'arx_na_nh3', 'zn_nh4cl', 'wurtz2', 'acyloin', 'rosenmund', 'lah_arno2', 'phenol_rx', 'alkyne_mild', 'sugar_ox', 'ald_hypohalite', 'hno3_cho', 'seo2', 'seo2_allyl', 'hio4_all'])
  wrap(id, function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (o && o.steps && arrowless(o)) o.steps.forEach(s => { s.redox = true; }); return o; });
R18.useMech2 = useMech2;
}

/* ---- run 38: two different partners that can both play both roles give every combination (statistical mixture) ---- */
const flat = v => Array.isArray(v) ? v : [v];
function partnersOf(c) { const out = []; for (const [k, v0] of c.st.sp) for (const v of flat(v0)) { if (!v || !v.smi) continue; let s = null;
    if (k === 'partner') s = v.smi; else if (k === 'ROH') s = 'O' + v.smi; else if (k === 'RX') s = (v.X || 'Br') + v.smi; else if (k === 'RNH2') s = 'N' + v.smi; else if (/[Cc]/.test(v.smi) && !['RMgX', 'RLi', 'R2CuLi', 'NaOR', 'ArN2'].includes(k)) s = v.smi;
    const q = s && canonS(s); if (q) out.push({ key: k, smi: q }); } return out; }
const MIXT = {
  aldol: { ok: s => { const q = fresh(s); return q && M(q, '[CX4;!H0][CX3;!$(C[O,N,Cl,Br])]=O').length && M(q, '[CX3]=O').length === 1 && !M(q, 'C=CC=O').length && !M(q, 'C#C').length; },
    t: [['[CX4;!H0:1][CX3;!$(C[O,N,Cl,Br]):2]=[O:3].[CX3;!$(C[O,N,Cl,Br]):4]=[O:5]>>[O:5][C:4][C:1][C:2]=[O:3]'], ['[CX4;!H0;!H1:1][CX3;!$(C[O,N,Cl,Br]):2]=[O:3].[CX3;!$(C[O,N,Cl,Br]):4]=[O]>>[C:4]=[C:1][C:2]=[O:3]']],
    pick: c => heat(c) ? 1 : 0, why: 'Both carbonyls have α-H, so each forms an enolate and each is attacked: four crossed / self aldols (a mixture, so this is a poor synthesis).' },
  decarbox: { gate: c => c.has('electrolysis') || /kolbe|electroly/i.test(c.raw || ''), ok: s => { const q = fresh(s); return q && M(q, '[CX3](=O)[OX2H1,OX1-]').length === 1; },
    t: [['[#6:1][CX3](=O)[OX2H1,OX1-].[#6:2][CX3](=O)[OX2H1,OX1-]>>[*:1][*:2]']], pick: () => 0, why: 'Kolbe radicals R• and R′• combine at random: R–R, R–R′ and R′–R′ all form.' },
  ether_413: { ok: s => { const q = fresh(s); return q && M(q, '[CX4][OX2H1]').length === 1 && M(q, '[CX4;H2,H3][OX2H1]').length === 1; },
    t: [['[CX4:1][OX2H1:2].[CX4:3][OX2H1]>>[C:1][O:2][C:3]']], pick: () => 0, why: 'Either alcohol can be protonated and either can attack: a mixture of three ethers (only symmetrical ethers are made this way).' },
  wurtz2: { ok: s => { const q = fresh(s); return q && M(q, '[CX4][Cl,Br,I]').length === 1 && !M(q, 'c').length; },
    t: [['[CX4:1][Cl,Br,I].[CX4:2][Cl,Br,I]>>[C:1][C:2]']], pick: () => 0, why: 'Two different R–X give R–R, R–R′ and R′–R′ (a mixture: Wurtz is used only for symmetrical alkanes).' } };
for (const [id, D] of Object.entries(MIXT)) wrap(id, function (m, c, all, k, run0) { if (D.gate && !D.gate(c)) return undefined; const a = smiOf(m); if (!a || !D.ok(a)) return undefined;
  const pt = partnersOf(c).find(p => p.smi !== a && D.ok(p.smi)); if (!pt) return undefined; const b = pt.smi; const T = D.t[D.pick(c)] || D.t[0]; const prods = [];
  for (const [x, y] of [[a, a], [b, b], [a, b], [b, a]]) for (const sma of T) for (const pr of rxn(sma, [x, y])) { const main = byC(pr.join('.').split('.'))[0]; const q = canonS(main); if (q && !prods.includes(q)) prods.push(q); }
  if (prods.length < 2) return undefined; let o; try { o = run0(m, c, all, k); } catch (e) { o = null; } const mainS = o && o.product ? smiOf(o.product) : null;
  const list = mainS && prods.includes(mainS) ? [mainS].concat(prods.filter(x => x !== mainS)) : prods; const base = o && o.product ? o : finish(m, list[0], { title: 'Mixture', text: D.why, tags: [] });
  if (!base) return undefined; const P0 = mainS && prods.includes(mainS) ? base.product : fresh(list[0]);
  return Object.assign({}, base, { product: P0, title: (base.title || '') + ' (all ' + list.length + ' combinations form)', minor: list.slice(1).map(fresh).filter(Boolean), why: (base.why || []).concat([D.why]) }); });
/* crossed Claisen: the ester with no α-H can only be the acceptor */
wrap('claisen_cond', function (m, c) { const a = smiOf(m); const pt = partnersOf(c).map(p => p.smi); if (!a || !pt.length) return undefined; const fr = [a].concat(pt);
  const hasEs = s => { const q = fresh(s); return q && M(q, '[CX3](=O)[OX2][#6]').length > 0; }, hasA = s => { const q = fresh(s); return q && M(q, '[CX4;!H0;!H1][CX3](=O)[OX2][#6]').length > 0; };
  const acc = fr.find(s => hasEs(s) && !hasA(s)), don = fr.find(s => s !== acc && hasA(s)); if (!acc || !don) return undefined;
  const p = rxn('[CX4;!H0;!H1:1][CX3:2](=[O:5])[OX2:3][#6:6].[CX3:7](=[O:8])[OX2][#6]>>[C:1]([C:2](=[O:5])[O:3][*:6])[C:7]=[O:8]', [don, acc]); if (!p.length) return undefined; const fr2 = byC(p[0].join('.').split('.'));
  return finish(m, fr2[0], { title: 'Crossed Claisen condensation → β-keto ester (the ester without α-H is only the acceptor)', text: 'Ethoxide forms the enolate of the ester that has α-H; it attacks the C=O of the ester with no α-H (that one cannot form an enolate), and EtO⁻ leaves. Self-Claisen of the α-H ester is the minor side reaction.', tags: ['condensation'] }); });

/* ---- run 38b: SKM Advanced gaps, as general rules ---- */
/* Robinson annulation: enolisable carbonyl + methyl vinyl ketone-type enone, base → Michael, then intramolecular aldol condensation */
front({ id: 'robinson', name: 'Robinson annulation', fam: 'Carbonyl', canon: 'base, Δ (with CH₂=CH–CO–CH₃)', db: ['Robinson annulation', 'Michael addition'],
  test: c => c.any(['KOH', 'aqKOH', 'NaOEt', 'NaOMe', 'KOtBu']),
  run(m, c) { const a = smiOf(m); if (!a) return null; const fr = a.split('.').concat(partnersOf(c).map(p => p.smi)); const en = fr.find(s => { const q = fresh(s); return q && M(q, '[CH2]=[CH1][CX3](=O)[CH3]').length; }); if (!en) return null;
    const dn = fr.find(s => s !== en && (() => { const q = fresh(s); return q && M(q, '[CX4;!H0][CX3;!$(C[O,N,Cl,Br])]=O').length && !M(q, 'C=CC=O').length; })()); if (!dn) return null;
    const p = rxn('[CX4;!H0:1][CX3;!$(C[O,N,Cl,Br]):2]=[O].[CH2:4]=[CH1:5][CX3:6](=[O:7])[CH3:8]>>[C:1]1[C:4][C:5][C:6](=[O:7])[C:8]=[C:2]1', [dn, en]); if (!p.length) return null; const prods = byC(p.map(x => x.join('.')).map(x => byC(x.split('.'))[0]));
    const mich = rxn('[CX4;!H0:1][CX3;!$(C[O,N,Cl,Br]):2]=[O:3].[CH2:4]=[CH1:5][CX3:6](=[O:7])[CH3:8]>>[C:1]([C:4][C:5][C:6](=[O:7])[C:8])[C:2]=[O:3]', [dn, en]);
    const o = finish(m, prods[0], { title: 'Robinson annulation → cyclohex-2-enone (Michael addition, then intramolecular aldol condensation)', text: 'The enolate of the carbonyl compound adds 1,4 to the enone (Michael). The CH₃ of the new methyl ketone then forms an enolate that attacks the other C=O intramolecularly (6-ring), and heat dehydrates the aldol to the conjugated cyclohexenone.', tags: ['enolate'] });
    if (o && mich.length) { const mi = fresh(byC(mich[0].join('.').split('.'))[0]); if (mi) o.steps = [step(m, mi, 'Michael addition: the enolate (from the α-C–H) adds to the β-carbon of the enone; protonation gives a 1,5-dicarbonyl compound.'), step(mi, o.product, 'Intramolecular aldol: the enolate of the CH₃–C=O attacks the other carbonyl (a 6-membered ring closes); E1cB loss of water gives the α,β-unsaturated ketone.')]; }
    return o; } });
/* benzidine rearrangement: hydrazobenzene + acid → 4,4′-diaminobiphenyl ([5,5]-sigmatropic) */
front({ id: 'benzidine', name: 'Benzidine rearrangement', fam: 'Amine', canon: 'H⁺ (HCl), warm', db: ['Benzidine rearrangement'],
  test: c => c.any(['H3O+', 'HCl', 'H2SO4']) || /h\+|acid/i.test(c.raw || ''),
  run(m) { const s = smiOf(m); if (!s || !M(m, 'c[NX3H1][NX3H1]c').length) return null; const p = rxn('[cH1:5]1[cH1:6][c:1]([NH1:2][NH1:3][c:4]2[cH1:7][cH1:8][cH1:9][cH1:10][cH1:11]2)[cH1:12][cH1:13][cH1:14]1>>N[c:1]1[cH1:6][cH1:5][c:14](-[c:9]2[cH1:8][cH1:7][c:4](N)[cH1:11][cH1:10]2)[cH1:13][cH1:12]1', [s]);
    if (!p.length) return null; return finish(m, mainFrag(p[0].join('.')), { title: 'Benzidine rearrangement: hydrazobenzene + H⁺ → 4,4′-diaminobiphenyl (benzidine)', text: 'Both N are protonated; the N–N bond breaks as a new C–C bond forms between the two para carbons in one concerted [5,5]-sigmatropic shift; loss of the two ring H⁺ restores both rings.', tags: ['rearrangement'] }); } });
/* Prévost: I2 + silver benzoate (dry) → trans-1,2-dibenzoate; hydrolysis → trans (anti) diol */
front({ id: 'prevost', name: 'Prévost reaction', fam: 'Alkene', canon: 'I₂, PhCOOAg (dry), then H₃O⁺', db: ['Prévost reaction'],
  test: c => c.has('I2') && /silver\s*benzoate|phcooag|c6h5cooag|agobz|silver\s*acetate|ch3cooag/i.test(c.raw || ''),
  run(m, c, all, k) { const s = smiOf(m); if (!s || !M(m, '[CX3]=[CX3]').length) return null; const hyd = all && all.slice(k + 1).some(x => x && x.any(['H3O+', 'KOH', 'aqKOH', 'H2O']));
    const tmpl = hyd ? '[C:1]=[C:2]>>[C:1](O)[C:2]O' : '[C:1]=[C:2]>>[C:1](OC(=O)c1ccccc1)[C:2]OC(=O)c1ccccc1'; const p = rxn(tmpl, [s]); if (!p.length) return null; let prod = mainFrag(p[0].join('.'));
    const iso = rxn('[C;R:1]=[C;R:2]>>[C@@H:1](O)[C@H:2]O', [s]); if (hyd && iso.length) prod = mainFrag(iso[0].join('.'));
    return finish(m, prod, { title: hyd ? 'Prévost: I₂ / silver benzoate, then hydrolysis → trans-1,2-diol (anti)' : 'Prévost: I₂ / silver benzoate → trans-1,2-dibenzoate', text: 'I⁺ forms a cyclic iodonium ion; benzoate opens it anti, then the neighbouring ester carbonyl displaces I⁻ (a cyclic benzoxonium ion); a second benzoate opens that anti again. Both O end up trans; hydrolysis gives the trans diol (compare cold KMnO₄ / OsO₄: cis).', tags: ['addition'] }); } });
/* Houben–Hoesch / Gattermann-type: arene + R–C≡N, HCl, Lewis acid → ketimine salt → H2O → aryl ketone */
front({ id: 'hoesch', name: 'Houben–Hoesch acylation', fam: 'Aromatic', canon: 'R–C≡N, HCl, ZnCl₂ / AlCl₃; then H₂O', db: ['Houben–Hoesch reaction'],
  test: c => c.has('HCl') && c.any(['AlCl3', 'ZnCl2', 'FeX3']) && /ch3cn|c2h5cn|nitrile|r-?cn|acetonitrile|phcn/i.test(c.raw || ''),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, 'c[H,#1]').length && !M(m, '[cH]').length) return null; const raw = (c.raw || '').toLowerCase(); const R = /c2h5cn|ch3ch2cn|propanenitrile/.test(raw) ? 'CC' : /phcn|c6h5cn|benzonitrile/.test(raw) ? 'c1ccccc1' : 'C';
    let o = null; const r = recipe('eas'); if (r) { try { o = r.run(m, ctxWith(c, { AcCl: { acyl: R === 'C' ? 'CC(=O)' : R + 'C(=O)', smi: (R === 'C' ? 'CC(=O)' : R + 'C(=O)') + 'Cl' }, HCl: null }), [], 0); } catch (e) { o = null; } }
    if (!o || !o.product) { const p = rxn('[cH:1]>>[c:1]C(=O)' + R, [s]); if (!p.length) return null; o = { product: fresh(mainFrag(p[0].join('.'))) }; }
    return finish(m, smiOf(o.product), { title: 'Houben–Hoesch: Ar–H + R–C≡N / HCl / Lewis acid → ketimine, then H₂O → aryl ketone Ar–CO–R', text: 'HCl and the Lewis acid turn the nitrile into an electrophilic nitrilium / iminium species R–C≡N⁺H; the ring attacks it (EAS) to give a ketimine salt Ar–C(R)=NH₂⁺, which water hydrolyses to the ketone. (Works best on electron-rich rings such as phenols.)', tags: ['eas'] }); } });
/* nitration of a phenol / aniline sulphonic acid replaces the SO3H groups too (route to picric acid) */
wrap('eas', function (m, c, all, k, run0) { if (!c.has('HNO3') || !M(m, '[OX2H1]c').length || !M(m, 'c[SX4](=O)(=O)[OX2H1]').length) return undefined; const s = smiOf(m); const r = applyAll(s, [['[c:1][SX4](=O)(=O)[OX2H1]>>[c:1]', 'ipso']], 3); const m2 = fresh(r.cur); if (!m2) return undefined; let o; try { o = run0(m2, c, all, k); } catch (e) { o = null; } if (!o || !o.product) return undefined;
  return finish(m, smiOf(o.product), { title: 'Nitration of a phenolsulphonic acid: NO₂⁺ replaces the –SO₃H groups as well (ipso) → ' + (M(o.product, 'c[N+](=O)[O-]').length >= 3 ? 'picric acid (2,4,6-trinitrophenol)' : 'nitrophenol'), text: 'Sulphonating first keeps the ring from being oxidised by conc. HNO₃; the –SO₃H groups are then displaced by NO₂⁺ (ipso attack, SO₃ lost). This is the industrial route to picric acid.', extra: ['O=S(=O)=O'], tags: ['eas'] }); });
/* SNAr: an amine present is the nucleophile, the base only removes HX */
wrap('aryl_snar', function (m, c) { if (c.any(['RNH2', 'R2NH', 'NH3', 'aniline'])) return null; return undefined; });
/* cyclic imide + Br2 / NaOH → OH⁻ opens the imide, then Hofmann on the amide half (phthalimide → anthranilic acid) */
wrap('amide_rx', function (m, c, all, k, run0) { if (!((c.has('Br2') && c.any(['KOH', 'aqKOH'])) || c.has('Br2/NaOH') || c.has('NaOCl')) || !M(m, '[CX3;R](=O)[NX3H1;R][CX3;R]=O').length) return undefined; const o = run0(m, c, all, k); if (!o || !o.product) return o;
  const r = applyAll(smiOf(o.product), [['[#6:1][CX3](=O)[NX3H2]>>[*:1]N', 'Hofmann']], 1); if (!r || r.cur === smiOf(o.product)) return o;
  return finish(m, r.cur, { title: 'Cyclic imide + Br₂ / NaOH → amino acid (Hofmann): phthalimide → anthranilic acid', text: 'OH⁻ opens the imide to the amide-carboxylate; the amide half then does the Hofmann bromamide reaction (N-bromoamide → isocyanate → amine + CO₂). Acidify for the free amino acid.', extra: ['O=C=O'], tags: ['hofmann'] }); });
/* malonic / acetoacetic ester + 1,n-dihalide (n = 2–5) → the second alkylation closes a 3- to 6-ring */
wrap('malonic_alk', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product) return o; const s = smiOf(o.product); if (!s || !M(o.product, '[CX4;H1]([CX3]=O)([CX3]=O)[CH2][Cl,Br,I]').length && !M(o.product, '[CX4;H1]([CX3]=O)([CX3]=O)[CH2][CH2][Cl,Br,I]').length && !M(o.product, '[CX4;H1]([CX3]=O)([CX3]=O)[CH2][CH2][CH2][Cl,Br,I]').length && !M(o.product, '[CX4;H1]([CX3]=O)([CX3]=O)[CH2][CH2][CH2][CH2][Cl,Br,I]').length) return o;
  for (const n of [2, 3, 4, 5]) { const ch = [...Array(n)].map((_, i) => '[CH2:' + (10 + i) + ']').join(''); const prodCh = [...Array(n)].map((_, i) => '[C:' + (10 + i) + ']').join(''); const p = rxn('[CX4;H1:1]([CX3:2]=[O:3])([CX3:4]=[O:5])' + ch + '[Cl,Br,I]>>[C:1]1([C:2]=[O:3])([C:4]=[O:5])' + prodCh + '1', [s]); if (p.length) { const pr = byC(p[0].join('.').split('.'))[0];
      return finish(m, pr, { title: 'Malonic ester + 1,' + n + '-dihalide, 2 NaOEt → cycloalkane-1,1-diester (' + (n + 1) + '-membered ring)', text: 'The first enolate alkylation puts the haloalkyl chain on the α-carbon; that carbon still has one acidic H, so a second enolate forms and displaces the other halide intramolecularly, closing the ring. Hydrolysis + decarboxylation then gives the cycloalkanecarboxylic acid.', extra: ['[Na+].[Br-]'], tags: ['enolate'] }); } }
  return o; });
/* alkyne + BH3, then CH3COOH: syn addition, protonolysis with retention → cis-alkene */
wrap('hb_follow', function (m, c, all, k, run0) { if (!M(m, '[#6][CX2]#[CX2][#6]').length || M(m, 'C=C').length) return undefined; if (!(all || [c]).slice(k).some(x => x && (x.any(['AcOH']) || /ch3cooh|acetic|propanoic|rcooh/i.test(x.raw || '')))) return undefined; const s = smiOf(m);
  const p = rxn('[#6:3][C:1]#[C:2][#6:4]>>[*:3]/[C:1]=[C:2]\\[*:4]', [s]); if (!p.length) return undefined; return finish(m, mainFrag(p[0].join('.')), { title: 'Alkyne + BH₃ (syn), then CH₃COOH → cis-alkene', text: 'B–H adds syn across one π bond (vinylborane); the carboxylic acid replaces B by H with retention, so both H are on the same side: the cis (Z) alkene. (Na / liq. NH₃ gives trans.)', tags: ['reduction'] }); });

/* ---- run 39: SKM Advanced gaps, general rules ---- */
/* break ring bonds safely (reaction SMARTS with '.' cannot open rings reliably): cut a–b, cap each end with an atom or H */
function cutCap(mol, cuts) { const t = mol.clone(); for (const [a, b, capA, capB] of cuts) { const bd = t.bondBetween(a, b); if (!bd) return null; t.removeBond(bd.id);
    for (const [at, cap] of [[a, capA], [b, capB]]) { const A = t.atoms.get(at); A.hFix = null; if (cap) { const nid = t.addAtom(cap, A.x + 20, A.y + 20, {}); t.addBond(at, typeof nid === 'number' ? nid : nid.id, 1); } } }
  try { return C.canon(t); } catch (e) { return null; } }

const partnerSmis = c => partnersOf(c).map(p => p.smi);
/* SOBr2 / "SOX2": same as SOCl2 but the new halogen is Br */
for (const id of ['alc_pcl', 'socl2_alc', 'alc_socl2']) wrap(id, function (m, c, all, k, run0) { const v = c.get('SOCl2'); if (!v || v.X !== 'Br') return undefined; const o = run0(m, c, all, k); if (!o || !o.product) return o; if (M(m, '[CX4]Cl').length) return o; const r = applyAll(smiOf(o.product), [['[CX4:1]Cl>>[C:1]Br', 'Br']], 4); return finish(m, r.cur, { title: (o.title || '').replace(/Cl/g, 'Br').replace(/chloride/g, 'bromide'), text: 'SOBr₂ works like SOCl₂: the alcohol becomes a chlorosulphite-type ester and Br⁻ displaces it (SO₂ escapes).', extra: ['O=S=O'], tags: o.tags || [] }); });
/* tert-butyl hypochlorite: radical chlorination of C–H like SO2Cl2 */
front({ id: 'tbuocl', name: 't-BuOCl radical chlorination', fam: 'Alkane', canon: 't-BuOCl, hν or Δ', db: [],
  test: c => c.has('tBuOCl'),
  run(m, c, all, k) { const r = recipe('radx'); if (!r) return null; let o = null; try { o = r.run(m, ctxWith(c, { Cl2: {}, tBuOCl: null }, ['light']), all, k); } catch (e) { o = null; } if (!o || !o.product) return null;
    return Object.assign({}, o, { title: 't-BuOCl, Δ / hν: free-radical chlorination (like Cl₂ / hν or SO₂Cl₂)', extra: (o.extra || []).concat([fresh('CC(C)(C)O')].filter(Boolean)) }); } });
/* CrO3 alone (anhydrous) stops at the aldehyde (NCERT); with H2SO4 / water (Jones) the 1° alcohol goes to the acid */
front({ id: 'cro3_anhyd', name: 'CrO₃ (anhydrous) oxidation', fam: 'Alcohol', canon: 'CrO₃ (anhydrous)', db: [],
  test: c => /cro3/i.test(c.raw || '') && !c.any(['H2SO4', 'H3O+', 'H2O', 'AcOH', 'KOH']) && !/jones|aq|dil/i.test(c.raw || ''),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX4;H1,H2][OX2H1]').length || M(m, '[CX3H1]=O').length) return null; const r = applyAll(s, [['[CX4;H2:1][OX2H1:2]>>[C:1]=[O:2]', '1° → CHO'], ['[CX4;H1;$(C([#6])[#6]):1][OX2H1:2]>>[C:1]=[O:2]', '2° → C=O']], 4); if (r.cur === s) return null;
    return finish(m, r.cur, { title: 'CrO₃ (anhydrous) oxidises 1° alcohol → aldehyde, 2° → ketone', text: 'Without water the aldehyde cannot form its hydrate, so it is not oxidised further (NCERT: "CrO₃ in anhydrous medium"). With acid / water (Jones) the 1° alcohol goes on to RCOOH.', tags: ['oxidation'] }); } });
/* MnO2 also oxidises propargylic and cyclopropylmethyl alcohols (π-like neighbours) */
for (const id of ['alc_ox_mild', 'ox_all_oh']) wrap(id, function (m, c, all, k, run0) { if (!c.has('MnO2') && !/mno2/i.test(c.raw || '')) return undefined; const q = M(m, '[OX2H1][CX4;H1,H2;$(CC#C),$(C[CX4;r3])]'); if (!q.length) return undefined; const s = smiOf(m);
  const r = applyAll(s, [['[OX2H1:2][CX4;H2;$(CC#C),$(C[CX4;r3]):1]>>[O:2]=[C:1]', 'CHO'], ['[OX2H1:2][CX4;H1;$(CC#C),$(C[CX4;r3]):1]>>[O:2]=[C:1]', 'C=O']], 2); if (r.cur === s) return undefined;
  return finish(m, r.cur, { title: 'MnO₂ oxidises an activated alcohol (propargylic / cyclopropylmethyl) → carbonyl', text: 'Active MnO₂ needs a π-type neighbour that stabilises the radical at the carbinol carbon: allylic, benzylic, propargylic, or the bent bonds of a cyclopropane ring.', tags: ['oxidation'] }); });
/* Oppenauer: Al(OtBu)3 / Al(OiPr)3 with excess acetone oxidises 1°/2° alcohols (C=C untouched) */
front({ id: 'oppenauer', name: 'Oppenauer oxidation', fam: 'Alcohol', canon: 'Al(OtBu)₃, acetone', db: ['Oppenauer oxidation'],
  test: c => c.has('AlOiPr') && (c.has('acetone') || /otbu|o-t-bu|tert-?butoxide|t-?butoxide|quinone/i.test(c.raw || '')) && !/isopropanol|propan-2-ol|ipr?oh/i.test(c.raw || ''),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX4;H1,H2][OX2H1]').length) return null; const r = applyAll(s, [['[CX4;H2:1][OX2H1:2]>>[C:1]=[O:2]', 'CHO'], ['[CX4;H1;$(C([#6])[#6]):1][OX2H1:2]>>[C:1]=[O:2]', 'C=O']], 2); if (r.cur === s) return null;
    return finish(m, r.cur, { title: 'Oppenauer oxidation: alcohol → aldehyde / ketone (C=C left alone)', text: 'The aluminium alkoxide passes a hydride from the carbinol C–H to acetone through a six-membered ring (the reverse of MPV reduction); excess acetone drives it. C=C bonds are not touched.', extra: ['CC(C)O'], tags: ['oxidation'] }); } });
/* Ag2O (moist / ammoniacal) oxidises an aldehyde to the acid, like Tollens */
front({ id: 'ag2o_cho', name: 'Ag₂O oxidises –CHO → –COOH', fam: 'Aldehyde', canon: 'moist Ag₂O', db: [],
  test: c => c.has('Ag2O') && !c.any(['RX', 'Cl-']),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX3H1](=O)[#6,#1]').length) return null; const r = applyAll(s, [['[CX3H1:1](=[O:2])>>[C:1](=[O:2])O', 'CHO → COOH']], 3); if (r.cur === s) return null;
    return finish(m, r.cur, { title: 'Ag₂O oxidises the aldehyde → carboxylic acid (silver mirror)', text: 'Ag(I) is a mild oxidant that attacks only the aldehyde C–H (the same chemistry as Tollens’ reagent); C=C and alcohols are untouched.', extra: ['[Ag]'], tags: ['oxidation'] }); } });
/* quaternary ammonium halide + moist Ag2O → R4N+ OH- + AgX (white ppt with Cl-) */
front({ id: 'quat_ag2o', name: 'R₄N⁺X⁻ + moist Ag₂O → R₄N⁺OH⁻', fam: 'Amine', canon: 'moist Ag₂O', db: [],
  test: c => c.has('Ag2O'),
  run(m, c) { if (!M(m, '[NX4+]').length) return null; const s = smiOf(m); const p = s.split('.').filter(x => !/^\[(Cl|Br|I)-\]$/.test(x)).join('.') + '.[OH-]'; const X = /Br-|\[Br-\]/.test(s + (c.raw || '')) ? 'Br' : /I-|\[I-\]/.test(s + (c.raw || '')) ? 'I' : 'Cl';
    const o = finish(m, p, { title: 'Quaternary ammonium halide + moist Ag₂O (AgOH) → R₄N⁺ OH⁻ + Ag' + X + '↓ (' + (X === 'Cl' ? 'white' : 'pale yellow') + ' ppt)', text: 'Ag⁺ removes the halide as insoluble Ag' + X + '; the strongly basic quaternary hydroxide stays in solution (heat it for Hofmann elimination).', extra: ['[Ag+].[' + X + '-]'], tags: ['test'] }); if (o) o.positive = true; return o; } });
/* CHCl3 + O2 (air) / O3 in light → phosgene + HCl (why chloroform is stored in dark bottles) */
front({ id: 'chcl3_air', name: 'CHCl₃ + air, light → COCl₂', fam: 'Halide', canon: 'O₂ (air), light', db: [],
  test: c => c.any(['O2', 'O3']) && (c.light || /light|sun|hν|hv|hnu/i.test(c.raw || '')),
  run(m) { const s = smiOf(m); if (s !== 'ClC(Cl)Cl') return null; return finish(m, 'O=C(Cl)Cl', { title: 'CHCl₃ + O₂ (air), light → COCl₂ (phosgene, very poisonous) + HCl', text: 'Light-driven radical oxidation of the C–H. That is why chloroform is kept in dark bottles filled to the top, with 1% ethanol added to destroy any phosgene.', extra: ['Cl'], tags: ['oxidation'] }); } });
/* V2O5, air, ~773 K: naphthalene → phthalic anhydride, benzene → maleic anhydride (NCERT) */
front({ id: 'v2o5_ox', name: 'V₂O₅ catalytic air oxidation of arenes', fam: 'Aromatic', canon: 'O₂, V₂O₅, 773 K', db: [],
  test: c => /v2o5/i.test(c.raw || ''),
  run(m) { const s = smiOf(m); if (s === 'c1ccc2ccccc2c1') return finish(m, 'O=C1OC(=O)c2ccccc21', { title: 'Naphthalene + O₂, V₂O₅, ~773 K → phthalic anhydride', text: 'One ring is burnt away (to CO₂ + H₂O); the two carbons left on the other ring become the cyclic anhydride.', extra: ['O=C=O', 'O'], tags: ['oxidation'] });
    if (s === 'c1ccccc1') return finish(m, 'O=C1C=CC(=O)O1', { title: 'Benzene + O₂, V₂O₅, ~773 K → maleic anhydride', text: 'Catalytic air oxidation opens the ring and loses two carbons as CO₂.', extra: ['O=C=O', 'O'], tags: ['oxidation'] }); return null; } });
/* Friedel–Crafts with a cyclic anhydride (succinic, phthalic): ring opens to a keto acid */
front({ id: 'fc_anhydride', name: 'Friedel–Crafts acylation with a cyclic anhydride', fam: 'Aromatic', canon: 'cyclic anhydride, AlCl₃', db: [],
  test: c => c.has('FeX3') || /alcl3/i.test(c.raw || ''),
  run(m, c) { const a = smiOf(m); if (!a) return null; const fr = a.split('.').concat(partnerSmis(c)); if (c.has('benzene') && !fr.includes('c1ccccc1')) fr.push('c1ccccc1'); const an = fr.find(x => { const q = fresh(x); return q && M(q, '[CX3;R](=O)[OX2;R][CX3;R]=O').length; }); const ar = fr.find(x => x !== an && /^c1ccccc1$|^Cc1ccccc1$/.test(x)); if (!an || !ar) return null;
    const ring = fresh(an); const q = M(ring, '[CX3;R:1](=[O:2])[OX2;R:3][CX3;R:4]=[O:5]'); if (!q.length) return null; const [c1, , o3] = q[0]; const acl = cutCap(ring, [[c1, o3, 'Cl', null]]); if (!acl) return null;
    const p = rxn('[CX3:1](=[O:2])Cl.[cH:6]>>[c:6][C:1]=[O:2]', [acl, ar]); if (!p.length) return null; const prod = byC(p[0].join('.').split('.'))[0]; const pm = fresh(prod); if (!pm || !M(pm, 'c[CX3](=O)').length || !M(pm, '[CX3](=O)[OX2H1]').length) return null;
    return finish(m, prod, { title: 'Friedel–Crafts acylation with a cyclic anhydride → keto acid (the ring opens)', text: 'AlCl₃ opens the anhydride to an acylium ion still joined to –COO–AlCl₃⁻; the arene attacks the acylium carbon (no rearrangement). Work-up gives Ar–CO–…–COOH.', tags: ['eas'] }); } });
/* acids that are β-keto / malonic after one loss of CO2 keep decarboxylating on heating */
wrap('acid_heat', function (m, c, all, k, run0) { let o = run0(m, c, all, k); if (!o || !o.product) return o; for (let g = 0; g < 3; g++) { const P = o.product; if (!M(P, '[CX3](=O)[CX4][CX3](=O)[OX2H1]').length && !M(P, '[OX2H1][CX3](=O)[CX4][CX3](=O)[OX2H1]').length) break; let o2 = null; try { o2 = run0(P, c, all, k); } catch (e) { o2 = null; } if (!o2 || !o2.product || smiOf(o2.product) === smiOf(P)) break;
    o = Object.assign({}, o2, { title: 'Repeated decarboxylation: each β-keto / malonic –COOH loses CO₂ on heating', steps: (o.steps || []).concat(o2.steps || []), extra: (o.extra || []).concat(o2.extra || []) }); } return o; });
/* carbonates take three RMgX: R3C–OH + the alcohols / diol */
front({ id: 'carbonate_rmgx', name: 'Carbonate ester + excess RMgX → R₃C–OH', fam: 'Ester', canon: 'excess RMgX then H₃O⁺', db: [],
  test: c => c.has('RMgX'),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[#8X2]~[#6X3](=[#8X1])~[#8X2]').length) return null; const v = c.get('RMgX'); const R = flat(v)[0].smi || 'C';
    const q = M(m, '[#8X2:1]~[#6X3:2](=[#8X1:4])~[#8X2:3]')[0]; const t = m.clone(); const [o1, c2, o4, o3] = [q[0], q[1], q[2], q[3]]; for (const [a, b] of [[c2, o1], [c2, o3]]) { const bd = t.bondBetween(a, b); if (bd) t.removeBond(bd.id); } t.removeAtom(o4); t.removeAtom(c2); for (const at of [o1, o3]) { const A = t.atoms.get(at); if (A) A.hFix = null; }
    let left = null; try { left = C.canon(t); } catch (e) { left = null; } if (!left) return null; const fr = left.split('.'); const main = canonS('OC(' + R + ')(' + R + ')' + R);
    return finish(m, canonS('OC(' + R + ')(' + R + ')' + R) || main, { title: 'Carbonate + excess R–MgX → R₃C–OH (three R add to the carbonate carbon)', text: 'The first R⁻ gives an ester (an alkoxide leaves), the second a ketone (another alkoxide leaves), the third the 3° alkoxide. Work-up: R₃C–OH plus the two alcohols (or the diol of a cyclic carbonate).', extra: fr.filter(x => canonS(x) !== canonS('OC(' + R + ')(' + R + ')' + R)), tags: ['grignard'] }); } });
/* R–X + tertiary amine → quaternary ammonium salt (SN2; works on the drawn partner too) */
front({ id: 'quat_make', name: 'R–X + R₃N → R₄N⁺X⁻', fam: 'Amine', canon: 'R₃N', db: [],
  test: c => partnerSmis(c).some(x => { const q = fresh(x); return q && M(q, '[NX3H0;!$(NC=O);!$(N-a);!$(N=*);!$(N#*)]([CX4])([CX4])[CX4]').length; }),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[CX4;H2,H3][Cl,Br,I]').length) return null; const am = partnerSmis(c).find(x => { const q = fresh(x); return q && M(q, '[NX3H0;!$(NC=O);!$(N-a)]([CX4])([CX4])[CX4]').length; });
    const p = rxn('[CX4;H2,H3:1][Cl,Br,I:2].[NX3H0;!$(NC=O):3]>>[C:1][N+:3].[*-:2]', [s, am]); if (!p.length) return null; const fr = byC(p[0].join('.').split('.'));
    return finish(m, fr[0], { title: 'R–X + tertiary amine → quaternary ammonium salt (SN2)', text: 'The amine lone pair attacks the carbon of R–X from the back side; X⁻ becomes the counter-ion. A 3° amine stops here (no N–H left to lose).', extra: fr.slice(1), tags: ['sn2'] }); } });
/* Gabriel: K-phthalimide picks the better leaving group (I > Br > Cl) when the substrate has two */
front({ id: 'gabriel_k', name: 'Gabriel: K-phthalimide + R–X', fam: 'Amine', canon: 'K-phthalimide then NH₂NH₂', db: ['Gabriel phthalimide synthesis'],
  test: c => partnerSmis(c).some(x => /\[N-\]|N\[K\]|\[K\]/.test(x) && /C\(=O\)/.test(x)) || /phthalimide/i.test(c.raw || ''),
  run(m) { const s = smiOf(m); if (!s) return null; const hs = M(m, '[CX4;H2,H3][I]').length ? 'I' : M(m, '[CX4;H2,H3][Br]').length ? 'Br' : M(m, '[CX4;H2,H3][Cl]').length ? 'Cl' : null; if (!hs) return null;
    const p = rxn('[CX4;H2,H3:1][' + hs + ']>>[C:1]N1C(=O)c2ccccc2C1=O', [s]); if (!p.length) return null; const prod = byC(p[0].join('.').split('.'))[0];
    return finish(m, prod, { title: 'Gabriel step 1: phthalimide anion + R–' + hs + ' → N-alkylphthalimide (SN2' + (M(m, '[CX4][Cl,Br,I]').length > 1 ? ', at the C–' + hs + ': the better leaving group' : '') + ')', text: 'The imide N⁻ (stabilised by both C=O) is a good nucleophile but cannot alkylate twice. Hydrazine (or KOH, H₂O) then frees the 1° amine.', extra: ['[K+].[' + hs + '-]'], tags: ['sn2'] }); } });

/* NBS without acid is a radical (allylic / benzylic) brominating agent; α-bromination of a ketone needs acid catalysis */
wrap('alpha_x_ket', function (m, c) { if (!c.has('NBS') || c.any(['H3O+', 'HCl', 'H2SO4', 'pTsOH', 'TsOH', 'AcOH'])) return undefined; if (M(m, '[CX4;!H0;$(Cc),$(CC=C)]').length) return null; return undefined; });
/* Cu2Cl2 / NH4Cl (no ammonia): acetylene dimerises to vinylacetylene (Nieuwland); ammoniacal Cu2Cl2 gives the red acetylide */
wrap('acetylide', function (m, c) { if (c.has('CuCl') && c.has('NH4Cl') && !c.has('NH3')) return null; return undefined; });
front({ id: 'nieuwland', name: 'Acetylene + Cu₂Cl₂ / NH₄Cl → vinylacetylene', fam: 'Alkyne', canon: 'Cu₂Cl₂, NH₄Cl', db: [],
  test: c => c.has('CuCl') && c.has('NH4Cl') && !c.has('NH3'),
  run(m) { const s = smiOf(m); if (s !== 'C#C') return null; return finish(m, 'C=CC#C', { title: 'HC≡CH + Cu₂Cl₂ / NH₄Cl → CH₂=CH–C≡CH (vinylacetylene, dimerisation)', text: 'In the slightly acidic NH₄Cl medium the copper(I) acetylide is not precipitated; it adds to a second acetylene instead (Nieuwland). With ammoniacal Cu₂Cl₂ the red copper acetylide precipitates. HCl then adds to vinylacetylene to give chloroprene.', tags: ['addition'] }); } });
/* Friedel–Crafts with an alkene under the industrial conditions (high pressure, heat): the cumene process */
front({ id: 'fc_alk_press', multi: true, name: 'Arene + alkene, high pressure, heat (cumene process)', fam: 'Aromatic', canon: 'H₃PO₄, ~520 K, high pressure', db: [],
  test: c => c.flag('pressure') && c.heat && !c.any(['H2SO4', 'H3PO4', 'HF', 'BF3', 'H3O+', 'FeX3']),
  run(m, c, all, k) { const r = recipe('fc_alk_drawn'); if (!r) return null; let o = null; try { o = r.run(m, ctxWith(c, { H3PO4: {} }), all, k); } catch (e) { o = null; } if (!o || !o.product) return null; return Object.assign({}, o, { title: 'Arene + alkene at high pressure and heat (acid catalyst, H₃PO₄) → alkylbenzene: the cumene process' }); } });
/* Michael addition of a malonate to an enone, then intramolecular Claisen onto the methyl ketone → cyclohexane-1,3-dione ester (dimedone route) */
front({ id: 'michael_claisen', name: 'Malonate + enone (NaOEt) → Michael, then intramolecular Claisen', fam: 'Carbonyl', canon: 'CH₂(COOEt)₂, NaOEt', db: ['Michael addition'],
  test: c => c.any(['NaOEt', 'NaOMe']) && partnerSmis(c).some(x => /^CCOC\(=O\)CC\(=O\)OCC$/.test(x) || /^COC\(=O\)CC\(=O\)OC$/.test(x)),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[CH3][CX3](=O)[CH1]=[CX3]').length) return null; const mal = partnerSmis(c).find(x => /C\(=O\)CC\(=O\)O/.test(x));
    const p = rxn('[CH3:1][CX3:2](=[O:3])[CH1:4]=[CX3:5].[CH2:6]([C:7](=[O:8])[O][CH2][CH3])[C:10](=[O:11])[O:12][#6:13]>>[C:6]1([C:10](=[O:11])[O:12][*:13])[C:5][C:4][C:2](=[O:3])[C:1][C:7]1=[O:8]', [s, mal]); if (!p.length) return null; const prod = byC(p[0].join('.').split('.'))[0];
    const mi = rxn('[CH3:1][CX3:2](=[O:3])[CH1:4]=[CX3:5].[CH2:6]([C:7](=[O:8])[O:9][#6:14])[C:10](=[O:11])[O:12][#6:13]>>[CH3:1][C:2](=[O:3])[C:4][C:5][C:6]([C:7](=[O:8])[O:9][*:14])[C:10](=[O:11])[O:12][*:13]', [s, mal]);
    const o = finish(m, prod, { title: 'Malonate + α,β-unsaturated ketone, NaOEt → Michael adduct, then intramolecular Claisen → cyclohexane-1,3-dione ester', text: 'The malonate enolate adds 1,4 to the enone (Michael). The CH₃ of the methyl ketone then forms an enolate that attacks one ester C=O (6-ring, EtO⁻ leaves). Acid + heat hydrolyses and decarboxylates the remaining ester: mesityl oxide gives dimedone.', extra: ['CCO'], tags: ['enolate'] });
    if (o && mi.length) { const mm = fresh(byC(mi[0].join('.').split('.'))[0]); if (mm) o.steps = [step(m, mm, 'Michael addition: the malonate enolate attacks the β-carbon of the enone; protonation gives the 1,5-dicarbonyl adduct.'), step(mm, o.product, 'Intramolecular Claisen: the enolate of the CH₃–C=O attacks an ester carbonyl, EtO⁻ leaves and the six-membered 1,3-dione ring closes.')]; } return o; } });
/* DIBAL-H on a lactone → lactol ⇌ hydroxy aldehyde (Tollens positive) */
wrap('dibal', function (m, c, all, k, run0) { const q = M(m, '[CX3;R:1](=[OX1:2])[OX2;R:3][CX4;R:4]'); if (!q.length || q.length > 1) return undefined; const [c1, , o3] = q[0]; const smi = cutCap(m, [[c1, o3, null, null]]); if (!smi) return undefined;
  return finish(m, smi, { title: 'DIBAL-H (−78 °C) on a lactone → hydroxy aldehyde (lactol in equilibrium); Tollens positive', text: 'One H⁻ adds to the lactone C=O; at −78 °C the tetrahedral aluminium intermediate does not collapse. Work-up releases the lactol, which opens to the ω-hydroxy aldehyde. LiAlH₄ would go on to the diol.', tags: ['reduction'] }); });
/* Perkin-type condensation of ArCHO with a cyclic anhydride that has α-CH2 (NaOAc, heat) → arylidene anhydride */
front({ id: 'perkin_cyclic', name: 'ArCHO + cyclic anhydride (NaOAc, Δ) → arylidene anhydride', fam: 'Carbonyl', canon: 'NaOAc, Δ', db: ['Perkin reaction'],
  test: c => c.any(['NaOAc', 'KOAc']) || /naoac|ch3coona|sodium\s*acetate/i.test(c.raw || ''),
  run(m, c) { const a = smiOf(m); if (!a) return null; const fr = a.split('.').concat(partnerSmis(c)); const ald = fr.find(x => { const q = fresh(x); return q && M(q, '[CX3H1](=O)c').length; }); const an = fr.find(x => { const q = fresh(x); return q && M(q, '[CX3;R](=O)[OX2;R][CX3;R](=O)').length && M(q, '[CH2;R][CX3;R](=O)[OX2;R]').length; }); if (!ald || !an) return null;
    const p = rxn('[CX3H1:1](=O)[c:2].[CH2;R:3][CX3;R:4](=[O:5])[OX2;R:6]>>[c:2][CH1:1]=[C:3][C:4](=[O:5])[O:6]', [ald, an]); if (!p.length) return null; const prod = byC(p[0].join('.').split('.'))[0];
    return finish(m, prod, { title: 'Perkin-type condensation: ArCHO + cyclic anhydride (NaOAc, Δ) → arylidene anhydride', text: 'Acetate removes an α-H of the anhydride; the enolate adds to ArCHO and loses water (E1cB), giving Ar–CH= on the α-carbon. With acetic anhydride the same reaction gives cinnamic acid after hydrolysis.', extra: ['O'], tags: ['condensation'] }); } });
/* carbonyl + a molecule with two OH 1,2 or 1,3 apart (glycol, catechol, o-hydroxybenzyl alcohol), H+ → cyclic acetal */
add({ id: 'cyc_acetal', name: 'Carbonyl + diol (H⁺, −H₂O) → cyclic acetal', fam: 'Carbonyl', canon: 'diol, H⁺', db: ['Acetal formation'],
  test: c => c.any(['H3O+', 'H2SO4', 'HCl', 'pTsOH', 'TsOH']) || /h\+|acid/i.test(c.raw || ''),
  run(m, c) { const a = smiOf(m); if (!a) return null; const fr = a.split('.').concat(partnerSmis(c)); const co = fr.find(x => { const q = fresh(x); return q && M(q, '[CX3;$([CH1](=O)[#6]),$([CH2]=O),$(C(=O)([#6])[#6])]=O').length && !M(q, '[OX2H1]').length; }); if (!co) return null;
    const T = ['[OX2H1:2][c:3]:[c:4][OX2H1:5]', '[OX2H1:2][c:3]:[c:4][CX4:6][OX2H1:5]', '[OX2H1:2][CX4:3][CX4:4][OX2H1:5]', '[OX2H1:2][CX4:3][CX4:4][CX4:6][OX2H1:5]'];
    const P = ['[*:3]:[*:4]', '[*:3]:[*:4][*:6]', '[*:3][*:4]', '[*:3][*:4][*:6]'];
    for (const dio of fr) { if (dio === co) continue; const q = fresh(dio); if (!q) continue; for (let i = 0; i < T.length; i++) { if (!M(q, T[i]).length) continue; const p = rxn('[CX3:1]=[O].' + T[i] + '>>[C:1]1[O:2]' + P[i] + '[O:5]1', [co, dio]); if (p.length) { const prod = byC(p[0].join('.').split('.'))[0];
      return finish(m, prod, { title: 'Carbonyl + ' + ['catechol', 'o-hydroxybenzyl alcohol', '1,2-diol', '1,3-diol'][i] + ' (H⁺, water removed) → cyclic acetal (' + (i === 0 || i === 2 ? '5' : '6') + '-membered ring)', text: 'H⁺ activates the C=O; one OH adds (hemiacetal), the OH of the hemiacetal is protonated and leaves as water, and the second OH of the same molecule closes the ring. Removing water drives it; aqueous acid reverses it.', extra: ['O'], tags: ['acetal'] }); } } }
    return null; } });
/* tertiary propargylic alcohol + SOCl2 (no base): SNi′ → chloroallene */
front({ id: 'propargyl_sni', name: '3° propargylic alcohol + SOCl₂ → chloroallene (SNi′)', fam: 'Alcohol', canon: 'SOCl₂', db: [],
  test: c => c.has('SOCl2') && !c.flag('pyr') && !(c.get('SOCl2') || {}).X,
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX4;H0]([OX2H1])([#6])([#6])C#[CH1]').length) return null; const p = rxn('[CX4;H0:1]([OX2H1])C#[CH1:2]>>[C:1]=C=[CH1:2]Cl', [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: '3° propargylic alcohol + SOCl₂ → chloroallene (SNi′: Cl is delivered to the far alkyne carbon)', text: 'The chlorosulphite ester forms; a 3° carbon resists SN2, so the Cl of the –OS(O)Cl group is delivered intramolecularly to the end of the C≡C through a six-membered ring while SO₂ leaves: C=C=C–Cl.', extra: ['O=S=O', 'Cl'], tags: ['substitution'] }); } });
/* high-temperature halogenation of a simple alkene (no solvent): allylic substitution instead of addition */
wrap('allylic_x', function (m, c, all, k, run0) { return undefined; });
front({ id: 'allylic_x_heat', name: 'Alkene + X₂, heat (no solvent) → allylic substitution', fam: 'Alkene', canon: 'Cl₂ / Br₂, high T', db: ['Allylic substitution'],
  test: c => c.any(['Cl2', 'Br2']) && c.heat && c.st.sp.size === 1 && !c.flag('ccl4') && !c.flag('aq') && !/ccl4|ch2cl2|water|h2o|acoh|cs2/i.test(c.raw || ''),
  run(m, c, all, k) { const s = smiOf(m); if (!s || E.aromSet(m).size || M(m, 'C=CC=C').length || M(m, '[#8,#7]').length || !M(m, '[CX3]=[CX3][CX4;!H0]').length) return null; const r = recipe('allylic_x'); if (!r) return null; let o = null; try { o = r.run(m, Object.assign(ctxWith(c, {}), { T: 773 }), all, k); } catch (e) { o = null; } if (!o || !o.product) { const X = c.has('Cl2') ? 'Cl' : 'Br'; const p = rxn('[C:1]=[C:2][CX4;!H0:3]>>[C:1]=[C:2][C:3]' + X, [s]); if (!p.length) return null; return finish(m, mainFrag(p[0].join('.')), { title: 'Alkene + ' + X + '₂ at high temperature → allylic substitution', text: 'At high T the addition is reversible and X• abstracts the weak allylic C–H; the allylic radical takes X from X₂.', extra: ['H' + X], tags: ['radical'] }); }
    return Object.assign({}, o, { title: 'Alkene + X₂ at high temperature (no solvent) → allylic substitution (radical), not addition' }); } });
/* SN2 at an activated (benzylic / allylic) C–F with an alkoxide */
front({ id: 'rf_act_sn2', name: 'Benzylic / allylic fluoride + RO⁻ → ether', fam: 'Halide', canon: 'NaOR', db: [],
  test: c => c.any(['NaOEt', 'NaOMe', 'NaOR']),
  run(m, c) { const s = smiOf(m); if (!s || !M(m, '[CX4;H2,H1;$(Cc),$(CC=C)]F').length || M(m, '[CX4][Cl,Br,I]').length) return null; const v = c.get('NaOR'); const R = c.has('NaOEt') ? 'CC' : c.has('NaOMe') ? 'C' : (v && v.smi) || 'CC'; const p = rxn('[CX4;$(Cc),$(CC=C):1]F>>[C:1]O' + R, [s]); if (!p.length) return null;
    return finish(m, mainFrag(p[0].join('.')), { title: 'Benzylic fluoride + RO⁻ → benzyl ether (SN2; F⁻ is a poor leaving group, but the benzylic transition state is stabilised)', text: 'Reactivity R–I > R–Br > R–Cl ≫ R–F; only an activated (benzylic / allylic) C–F is displaced by a strong nucleophile.', extra: ['[Na+].[F-]'], tags: ['sn2'] }); } });
/* 1,3-dicarbonyl + 2 eq. strong base (LiNH2 / NaNH2 / LDA) → dianion: the terminal CH3 reacts with an ester (γ-acylation) */
front({ id: 'dianion_acyl', name: 'β-Diketone dianion + ester → γ-acylation', fam: 'Carbonyl', canon: '2 LiNH₂, then RCOOR′', db: [],
  test: c => c.any(['LiNH2', 'NaNH2', 'LDA']) || /linh2|nanh2|lda/i.test(c.raw || ''),
  run(m, c) { const a = smiOf(m); if (!a) return null; const fr = a.split('.').concat(partnerSmis(c)); const dk = fr.find(x => { const q = fresh(x); return q && M(q, '[CH3][CX3](=O)[CH2][CX3](=O)[#6]').length; }); const es = fr.find(x => x !== dk && (() => { const q = fresh(x); return q && M(q, '[CX3](=O)[OX2][CH3,CH2]').length; })()); if (!dk || !es) return null;
    const p = rxn('[CH3:1][CX3:2](=[O:3])[CH2:4][CX3:5](=[O:6]).[CX3:7](=[O:8])[OX2][#6]>>[C:7](=[O:8])[CH2:1][C:2](=[O:3])[C:4][C:5]=[O:6]', [dk, es]); if (!p.length) return null; const prod = byC(p[0].join('.').split('.'))[0];
    return finish(m, prod, { title: '1,3-Diketone + 2 equiv. LiNH₂ → dianion; the terminal CH₂⁻ is acylated by the ester (γ-acylation) → 1,3,5-triketone', text: 'The first equivalent removes the most acidic CH₂ proton (between the C=O groups); the second removes a terminal CH₃ proton. The second carbanion is the less stabilised, more reactive one, so the electrophile attacks there.', extra: ['CO'], tags: ['enolate'] }); } });
/* anthrone / anthraquinone + HI / red P → anthracene (the 9-anthrol enol is aromatic) */
wrap('hi_redp', function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product) return o; if (smiOf(o.product) !== 'c1ccc2c(c1)Cc1ccccc1C2') return o; return finish(m, 'c1ccc2cc3ccccc3cc2c1', { title: 'HI / red P on anthrone → anthracene', text: 'Anthrone exists partly as its aromatic enol, 9-anthrol; HI / P removes the C–OH as for a phenol-type C–O and the middle ring ends fully aromatic: anthracene (not 9,10-dihydroanthracene).', tags: ['reduction'] }); });

/* a Grignard / hydride step may swallow the next stage only when it is a plain work-up (water, acid, NH4Cl), not a new reagent such as Hg(OAc)2 / H2O */
{ const WORK = new Set(['H2O', 'H3O+', 'NH4Cl', 'HCl', 'H2SO4', 'AcOH', 'ether']);
  for (const r of E.REC) { if (!r || typeof r.eats !== 'function' || !/grignard|rmgx|rli|gilman|hyd_red|lah|nabh4|dibal/i.test(r.id)) continue; const e0 = r.eats; r.eats = function (st) { let ok = false; try { ok = e0.call(this, st); } catch (e) { ok = false; } if (!ok) return false; try { const keys = st.st ? [...st.st.sp.keys()] : []; if (keys.some(x => !WORK.has(x))) return false; } catch (e) { } return true; }; } }

/* aldoxime + Ac2O, heat → nitrile (the O-acetyl oxime eliminates AcOH); any OH groups are acetylated too (Wohl degradation, step 2) */
wrap('peracetyl', function (m, c, all, k, run0) { if (!M(m, '[CX3H1]=N[OX2H1]').length || !heat(c)) return undefined; const o = run0(m, c, all, k); if (!o || !o.product) return o;
  const r = applyAll(smiOf(o.product), [['[CX3H1:1]=[N:2][OX2][CX3](=O)C>>[C:1]#[N:2]', 'nitrile'], ['[CX3H1:1]=[N:2][OX2H1]>>[C:1]#[N:2]', 'nitrile']], 1);
  return finish(m, r.cur, { title: 'Aldoxime + Ac₂O, Δ → nitrile (dehydration); every –OH is acetylated as well', text: 'The oxime O is acetylated; AcO⁻ then leaves as the C–H is removed (an E2-type elimination), giving C≡N. With a sugar oxime this is step 2 of the Wohl degradation: an acetylated cyanohydrin.', extra: ['CC(=O)O'], tags: ['elimination'] }); });
front({ id: 'oxime_ac2o', name: 'Aldoxime + Ac₂O, Δ → nitrile', fam: 'Imine', canon: 'Ac₂O, Δ', db: [],
  test: c => c.has('Ac2O') && heat(c),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[CX3H1]=N[OX2H1]').length || M(m, '[OX2H1][CX4]').length) return null; const r = applyAll(s, [['[CX3H1:1]=[N:2][OX2H1]>>[C:1]#[N:2]', 'nitrile']], 1); if (r.cur === s) return null;
    return finish(m, r.cur, { title: 'Aldoxime + Ac₂O, Δ → nitrile (dehydration)', text: 'O-Acetylation turns the oxime OH into a leaving group; loss of AcOH across C–H / N–OAc gives the nitrile.', extra: ['CC(=O)O'], tags: ['elimination'] }); } });
/* acetylated (or free) cyanohydrin + base / moist Ag2O → the esters are cleaved and the cyanohydrin falls apart: carbonyl + HCN (Wohl, step 3) */
front({ id: 'wohl_end', name: 'Cyanohydrin collapse (Wohl degradation)', fam: 'Biomolecule', canon: 'NaOMe or moist Ag₂O', db: ['Wohl degradation'],
  test: c => c.any(['Ag2O', 'NaOMe', 'NH3', 'KOH', 'aqKOH']),
  run(m, c) { const s = smiOf(m); if (!s) return null; const ac = M(m, '[CX4]([OX2][CX3](=O)[CH3])C#N').length > 0; const free = M(m, '[CX4]([OX2H1])C#N').length > 0; if (!ac && !(free && c.has('Ag2O'))) return null;
    const r = applyAll(s, [['[OX2:1][CX3](=O)[CH3]>>[O:1]', 'OAc → OH']], 8); const r2 = applyAll(r.cur, [['[CX4:1]([OX2H1:2])C#N>>[C:1]=[O:2]', 'cyanohydrin → C=O']], 1); if (r2.cur === r.cur) return null;
    return finish(m, r2.cur, { title: 'Cyanohydrin (acetylated) + base / moist Ag₂O → carbonyl + HCN: the aldose is one carbon shorter (Wohl degradation)', text: 'The acetates are removed; the free cyanohydrin O⁻ pushes out CN⁻ (the reverse of HCN addition), so C-2 becomes the new –CHO. Net: aldohexose → aldopentose.', extra: ['C#N'], tags: ['degradation'] }); } });
/* nitrile + H2 over poisoned Pd (Pd/BaSO4): partial reduction to the imine; water / acid then gives the aldehyde (modern Kiliani–Fischer) */
front({ id: 'nitrile_pdbaso4', name: 'Nitrile + H₂, Pd/BaSO₄ → imine → aldehyde', fam: 'Nitrile', canon: 'H₂, Pd/BaSO₄; H₃O⁺', db: ['Kiliani–Fischer synthesis'],
  eats: s => s.any(['H3O+', 'HCl', 'H2O']) && [...s.st.sp.keys()].every(x => ['H3O+', 'HCl', 'H2O'].includes(x)),
  test: c => c.has('H2') && /baso4/i.test(c.raw || ''),
  run(m, c, all, k) { const s = smiOf(m); if (!s || !M(m, '[CX2]#N').length) return null; const hyd = all && all.slice(k + 1).some(x => x && x.any(['H3O+', 'HCl', 'H2O']));
    const r = applyAll(s, [hyd ? ['[CX2:1]#[N]>>[CH1:1]=O', 'CHO'] : ['[CX2:1]#[N:2]>>[CH1:1]=[NH1:2]', 'imine']], 1); if (r.cur === s) return null;
    return finish(m, r.cur, { title: hyd ? 'Nitrile + H₂, Pd/BaSO₄ (poisoned) → imine, then H₃O⁺ → aldehyde (Kiliani–Fischer: the aldose gains one carbon)' : 'Nitrile + H₂, Pd/BaSO₄ (poisoned) → imine', text: 'The poisoned catalyst adds only one H₂ to C≡N, stopping at the imine R–CH=NH; water hydrolyses it to R–CHO + NH₄⁺. Applied to the cyanohydrin of an aldose this lengthens the chain by one carbon (two epimers form).', extra: hyd ? ['[NH4+]'] : [], tags: ['reduction'] }); } });

/* Oppenauer conditions (Al alkoxide + excess acetone / quinone, no isopropanol) cannot reduce a ketone: nothing donates H⁻ net */
wrap('mpv', function (m, c) { if ((c.has('acetone') || /excess\s+acetone|quinone/i.test(c.raw || '')) && !/isopropanol|propan-2-ol|ipr?oh/i.test(c.raw || '') && !M(m, '[CX4;!H0][OX2H1]').length) return none('No net reaction', 'Al(OR)₃ with excess acetone are Oppenauer (oxidation) conditions: a ketone has no carbinol C–H to give, and the acetone ⇌ isopropoxide exchange is degenerate.'); return undefined; });
/* haloform-type cleavages: the carboxylate is the product, CHX3 the by-product */
for (const id of ['haloform', 'cx3_cleave', 'haloform_fix']) wrap(id, function (m, c, all, k, run0) { const o = run0(m, c, all, k); if (!o || !o.product) return o; const ps = smiOf(o.product); if (!/^(ClC\(Cl\)Cl|BrC\(Br\)Br|IC\(I\)I)$/.test(ps)) return o; const ex = (o.extra || []).map(x => smiOf(x)); const cb = ex.findIndex(x => x && x !== 'O=C[O-]' && x !== 'O=CO' && /C\(=O\)\[O-\]|C\(=O\)O$|\[O-\]C=O/.test(x)); if (cb < 0) return o; /* formate (from ethanol / acetaldehyde): the iodoform itself is the product that matters */
  const extra = (o.extra || []).slice(); const P = extra.splice(cb, 1)[0]; extra.unshift(o.product); return Object.assign({}, o, { product: P, extra }); });
/* heteroaryl carboxylic acids (pyridine, pyrrole, furan, thiophene) lose CO2 on strong heating */
front({ id: 'hetaryl_decarb', name: 'Heteroaryl –COOH, strong heat → heteroarene + CO₂', fam: 'Acid', canon: 'strong heat (soda lime / Cu)', db: [],
  test: c => c.heat && (c.st.sp.size === 0 || c.any(['sodalime', 'Cu', 'CaO'])),
  run(m) { const s = smiOf(m); if (!s || !M(m, '[$(c1ccncc1),$(c1ccccn1),$(c1cccnc1),$(c1cc[nH]c1),$(c1ccoc1),$(c1ccsc1),$(c1c[nH]cc1),$(c1cocc1)][CX3](=O)[OX2H1]').length) return null; const r = applyAll(s, [['[c:1][CX3](=O)[OX2H1]>>[cH1:1]', '-CO2']], 3); if (r.cur === s) return null;
    return finish(m, r.cur, { title: 'Heteroaromatic carboxylic acid, strong heat → loses CO₂ (pyridine-, pyrrole-, furan-carboxylic acids)', text: 'The ring N (or O) stabilises the carbanion / ylide left when CO₂ departs (2- fastest, via the ylide; 3- and 4- need stronger heating, or soda lime / Cu).', extra: ['O=C=O'], tags: ['decarboxylation'] }); } });

root.EngineR18 = R18;
})(typeof window !== 'undefined' ? window : globalThis);
