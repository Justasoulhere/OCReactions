/* engine_r6.js: JEE-syllabus audit (Oct 13). FRONT-registered, so these win over the broad recipes.
   Covers: KOH substitution vs elimination by conditions, gem/vic/tri-halides, dihalide → alkyne without allenes,
   Wurtz family (mixtures, 3° failure, intramolecular), aryl halides (Dow, SNAr, benzyne), conjugated dienes,
   Diels–Alder, thermal rearrangements, aromatics and a long tail of named carbonyl / acid reactions. */
(function (root) {
'use strict';
const C = root.Chem, E = root.EngineCore;
const { BL, el, P, out, done, hTot, cN, heavyN, aromSet, Mech } = E;
const REC = E.REC; const front = r => REC.unshift(r);
const M = (m, sm) => { try { return C.matchAll(m, sm, {}); } catch (e) { return []; } };
const tidy = m => { const t = m.clone(); try { C.relayout(t, BL); } catch (e) { } return t; };
const fresh = smi => C.fromSmilesMB(smi, BL);
const step = (L, R, text, cond) => ({ L: L.clone(), arrows: [], R: R.clone(), text, cond: cond || [] });
const ACID = ['H3O+', 'H2SO4', 'H3PO4', 'pTsOH', 'KHSO4', 'HCl'];
const isCO = (m, c) => E.isCarbonylC(m, c);
const clr = (t, ...ids) => ids.forEach(id => { const a = t.atoms.get(id); if (a) a.hFix = null; });
function addAt(t, at, smi, o = 1) { const g = C.fromSmilesMB(smi, BL); for (const a of g.atoms.values()) a.hFix = null; const first = [...g.atoms.keys()][0]; const p = P(t, at); const idm = t.merge(g, p[0] + BL, p[1]); t.addBond(at, idm.get(first), o); clr(t, at); return idm.get(first); }
function setO(t, a, b, o) { const x = t.bondBetween(a, b); if (x) x.o = o; else if (o) t.addBond(a, b, o); clr(t, a, b); }
function cutB(t, a, b) { const x = t.bondBetween(a, b); if (x) t.removeBond(x.id); clr(t, a, b); }
function fragsOf(t) { return t.frags().map(ids => ({ ids, n: ids.filter(id => el(t, id) === 'C').length, h: ids.length })).sort((p, q) => q.n - p.n || q.h - p.h).map(f => E.cleanProduct(t, f.ids)); }
function simple(L, t, o) { const R = tidy(t); const fr = R.frags(); let main = R, extra = o.extra || [];
  if (fr.length > 1) { const all = fragsOf(R); main = o.pick ? (all.find(o.pick) || all[0]) : all[0]; extra = extra.concat(all.filter(f => f !== main)); }
  return out(Object.assign({}, o, { product: main, extra, steps: [step(L, R, o.text, o.cond)].concat(o.more || []) })); }
function distMap(m, a) { const d = new Map([[a, 0]]), q = [a]; while (q.length) { const x = q.shift(); for (const y of m.nbrs(x)) if (!d.has(y)) { d.set(y, d.get(x) + 1); q.push(y); } } return d; }
function key(m) { try { return C.canon(m); } catch (e) { return null; } }
function ringOf(m, a) { return C.aromaticRings(m).find(r => r.includes(a)); }
function ringPos(ring, a, b) { const i = ring.indexOf(a), j = ring.indexOf(b); if (i < 0 || j < 0) return -1; const d = Math.abs(i - j); return Math.min(d, ring.length - d); }
function subsOnRing(m, ring) { const S = new Set(ring); const o = []; for (const r of ring) for (const n of m.nbrs(r)) if (!S.has(n) && el(m, n) !== 'H') o.push([r, n]); return o; }
const NO2 = (m, a) => el(m, a) === 'N' && m.atoms.get(a).chg > 0 && heavyN(m, a).filter(x => el(m, x) === 'O').length === 2;
const OHB = ['KOH', 'aqKOH', 'alcKOH'];

/* ===================== aryl halides ===================== */
front({ id: 'aryl_snar', name: 'Aryl halide + nucleophile (Dow / SNAr)', fam: 'Aryl halide', canon: 'NaOH, 623 K, 300 atm', db: ['Dow process', 'Nucleophilic aromatic substitution'],
  test: c => c.any(['KOH', 'aqKOH', 'NH3', 'NaOMe', 'NaOEt', 'alcKOH']) && !c.any(['NaNH2', 'CHCl3', 'Na', 'Br2/NaOH', 'NaOCl', 'I2', 'CO2', 'ArN2', 'Na/NH3', 'Br2', 'Cl2']),
  run(m, c) {
    if (E.halides(m).length) return null;
    const AX = M(m, '[c:1][F,Cl,Br,I:2]'); if (!AX.length) return null;
    let best = null;
    for (const [ca, x] of AX) { const ring = ringOf(m, ca); if (!ring) continue; let n = 0;
      for (const [r, s] of subsOnRing(m, ring)) if (NO2(m, s)) { const d = ringPos(ring, ca, r); if (d === 1 || d === 3) n++; }
      if (!best || n > best.n) best = { ca, x, n }; }
    if (!best) return null; const { ca, x, n } = best;
    const strong = (c.T != null && c.T >= 573) || c.flag('pressure');
    if (!n && !strong) return out({ title: 'No reaction under these conditions', product: null, none: 'Aryl halides resist nucleophilic substitution: C–X has partial double-bond character (resonance with the ring), the carbon is sp², and the π cloud repels the nucleophile. Chlorobenzene needs NaOH at 623 K and 300 atm (Dow process), unless –NO₂ groups sit ortho/para to the halogen.' });
    const nu = c.has('NH3') ? ['N', 'amine', 'NH₃'] : c.has('NaOMe') ? ['OC', 'methyl aryl ether', 'CH₃O⁻'] : c.has('NaOEt') ? ['OCC', 'ethyl aryl ether', 'C₂H₅O⁻'] : ['O', 'phenol (as phenoxide in the base)', 'OH⁻'];
    const t = m.clone(); t.removeAtom(x); addAt(t, ca, nu[0]);
    const T = ['about 623 K and 300 atm', 'about 443 K', 'about 368 K', 'warm water is enough'][Math.min(n, 3)];
    return simple(m, t, { title: n ? `SNAr (${n} –NO₂ ortho/para) → ${nu[1]}` : `Dow process: Ar–Cl + NaOH (623 K, 300 atm) → ${nu[1]}`,
      text: n ? `${nu[2]} adds to the carbon carrying the halogen; the negative charge of this σ-complex (Meisenheimer complex) spreads onto the oxygens of the ortho/para –NO₂ group(s). The halide then leaves and the ring is aromatic again: addition–elimination.` : `Under forcing conditions ${nu[2]} replaces Cl on the ring. Acidifying the sodium phenoxide gives phenol.`,
      why: [`Conditions: ${T}`, 'Each –NO₂ ortho/para to X stabilises the anionic intermediate: 1 → 443 K, 2 → 368 K, 3 → warm water', '–NO₂ meta to X does not help: the charge cannot reach it', nu[0] === '[O-]' ? 'H₃O⁺ work-up turns the phenoxide into the phenol' : ''], tags: ['SNAr'] });
  } });
front({ id: 'benzyne', name: 'Benzyne (elimination–addition)', fam: 'Aryl halide', canon: 'NaNH₂, liq. NH₃', db: [],
  test: c => c.has('NaNH2') || (c.has('NH3') && /k?nh2-|amide/i.test(c.raw || '')),
  run(m) {
    if (E.halides(m).length || M(m, '[CX3]([Cl,Br,I])=[CX3]').length) return null;
    const AX = M(m, '[c:1][Cl,Br,I:2]'); if (!AX.length) return null; const [ca, x] = AX[0]; const ring = ringOf(m, ca);
    const oh = m.nbrs(ca).filter(n => ring && ring.includes(n) && hTot(m, n) > 0);
    if (!oh.length) return out({ title: 'No ortho-H: no benzyne', product: null, none: 'Benzyne forms by removing the H ortho to the halogen; with both ortho positions substituted nothing happens.' });
    const mk = at => { const t = m.clone(); t.removeAtom(x); clr(t, ca); addAt(t, at, 'N'); return tidy(t); };
    const ipso = mk(ca); const cines = oh.map(mk).filter(p => key(p) !== key(ipso)); const uniq = [...new Map(cines.map(p => [key(p), p])).values()];
    const R = ipso;
    return out({ title: 'NaNH₂ / liq. NH₃: benzyne route → aniline', product: R, minor: uniq,
      steps: [step(m, R, 'NH₂⁻ pulls off the H ortho to the halogen; X⁻ leaves, giving benzyne (a strained, very reactive "triple bond" in the ring). NH₂⁻ adds to either end of it and the aryl anion takes H⁺ from NH₃.', ['NaNH₂', 'liq. NH₃'])],
      why: ['Elimination–addition, not SNAr: works with no –NO₂ present', uniq.length ? 'Substituted rings give two products: –NH₂ where X was (ipso) and on the neighbouring carbon (cine)' : 'Unsubstituted: both ends of benzyne give the same aniline', '¹⁴C labelling proved the symmetric intermediate'], tags: ['benzyne'] });
  } });

/* ===================== alkyl halides ===================== */
front({ id: 'rx_reduce', name: 'Alkyl halide → alkane (reduction)', fam: 'Alkyl halide', canon: 'Zn, HCl', db: [],
  test: c => c.any(['LiAlH4', 'Sn/HCl', 'Zn/HCl']) || (c.has('H2') && c.has('Pd')) || (c.has('Zn') && c.any(['HCl', 'AcOH', 'H3O+'])),
  run(m, c) {
    const hs = E.halides(m); if (!hs.length) return null;
    if (M(m, '[CX3]=[OX1]').length || M(m, '[N+](=O)[O-]').length || M(m, 'C#N').length || M(m, 'N=[N+]=[N-]').length) return null;
    if (c.has('H2') && (M(m, '[C;!a]=[C;!a]').length || M(m, 'C#C').length)) return null;
    const t = m.clone(); hs.forEach(([cc, x]) => { if (t.atoms.has(x)) t.removeAtom(x); clr(t, cc); });
    const how = c.has('LiAlH4') ? ['LiAlH₄', 'H⁻ from LiAlH₄ displaces the halide (SN2 by hydride).'] : c.has('H2') ? ['H₂ / Pd', 'Hydrogenolysis: H₂ on Pd breaks C–X and adds H.'] : ['Zn / HCl', 'Zn and acid supply electrons and H⁺ ("nascent hydrogen"): C–X → C–H.'];
    return simple(m, t, { title: `${how[0]}: R–X → R–H (alkane)`, text: how[1], why: ['Reduction of alkyl halides is a way to make alkanes with the same carbon skeleton'], tags: ['reduction'] });
  } });

/* KOH + alkyl halide: the conditions decide (aqueous → substitution, alcoholic / concentrated / hot → elimination) */
front({ id: 'rx_base', name: 'Alkyl halide + KOH (substitution vs elimination)', fam: 'Alkyl halide', canon: 'conc. KOH, Δ', stoich: true, db: ['SN2 reaction', 'E2 elimination', 'Saytzeff rule'],
  test: c => c.any(OHB) && !c.any(['CHCl3', 'Br2/NaOH', 'NaOCl', 'H2O2/OH', 'Br2', 'Cl2', 'I2', 'NaOEt', 'NaOMe', 'KOtBu', 'NaNH2', 'phenol', 'aniline', 'naphthol', 'ArN2', 'CO2', 'HCHO', 'CH3CHO', 'acetone', 'PhCHO', 'RX', 'Ac2O', 'AcCl']),
  run(m, c) {
    const hs = E.halides(m); if (!hs.length) return null;
    const byC = new Map(); for (const [cc, x] of hs) { if (!byC.has(cc)) byC.set(cc, []); byC.get(cc).push(x); }
    const wet = c.has('aqKOH') || c.flag('aq') || c.flag('dil') || c.has('H2O'); const alc = (c.has('alcKOH') || c.flag('alc')) && !wet; const aq = !alc && wet; // aqueous-alcoholic KOH: substitution
    const conc = c.flag('conc'), hot = c.heat;
    for (const [cc, xs] of byC) if (xs.length >= 3 && cN(m, cc).some(n => isCO(m, n))) return null; // CX₃–C=O: haloform-type cleavage (chloral, trihalo ketones)
    for (const [cc, xs] of byC) if (xs.length >= 3) {
      const t = m.clone(); xs.forEach(x => t.removeAtom(x)); clr(t, cc); addAt(t, cc, 'O', 2); const o2 = addAt(t, cc, 'O', 1); t.atoms.get(o2).chg = -1;
      return simple(m, t, { title: hTot(m, cc) ? 'CHX₃ + KOH → formate ion (HCOO⁻)' : '1,1,1-Trihalide + KOH → carboxylate', text: 'Each C–X is replaced by OH; a carbon carrying three –OH groups loses water to the acid, and the base turns it into its salt.', why: ['Three halogens on one carbon: the carbon ends up at the acid oxidation level', 'Chloroform + aq. KOH → potassium formate'], tags: ['SN'] });
    }
    for (const [cc, xs] of byC) if (xs.length === 2) {
      if (!alc) { const t = m.clone(); xs.forEach(x => t.removeAtom(x)); clr(t, cc); addAt(t, cc, 'O', 2); const ald = hTot(m, cc) > 0;
        return simple(m, t, { title: `gem-Dihalide + aq. KOH → ${ald ? 'aldehyde' : 'ketone'}`, text: 'Both halogens are replaced by –OH; the gem-diol is unstable and loses water to give the C=O compound.', why: ['Terminal gem-dihalide (–CHX₂) → aldehyde; internal → ketone', 'Two –OH on one carbon dehydrate spontaneously'], tags: ['SN'] }); }
      const r = E.elimX(m, cc, xs[0], '[OH-]', false, 'E2: alcoholic OH⁻ takes a β-H, the C=C forms and one halide leaves: a vinyl halide.');
      if (!r) return null;
      return done(r.mech, { title: 'gem-Dihalide + alc. KOH → vinyl halide', why: ['The second E2 (vinyl halide → alkyne) needs the stronger base NaNH₂ (or 2 eq / hot conc. alc. KOH)'], tags: ['E2'] });
    }
    const vic = hs.find(([a]) => hs.some(([b]) => b !== a && m.bondBetween(a, b)));
    if (vic) { const [a, xa] = vic; const [b, xb] = hs.find(([b2]) => b2 !== a && m.bondBetween(a, b2));
      if (!alc) { const t = m.clone(); [[a, xa], [b, xb]].forEach(([cc, x]) => { t.removeAtom(x); clr(t, cc); addAt(t, cc, 'O'); });
        return simple(m, t, { title: 'vic-Dihalide + aq. KOH → 1,2-diol (glycol)', text: 'Each C–X undergoes SN2 hydrolysis: OH⁻ replaces both halogens.', why: ['1,2-Dibromoethane + aq. KOH → ethylene glycol', 'With alc. KOH the same compound eliminates instead (vinyl halide)'], tags: ['SN2'] }); }
      const r = E.elimX(m, a, xa, '[OH-]', false, 'E2: alcoholic OH⁻ takes the H on the carbon that also carries a halogen; the C=C forms and one halide leaves.', hTot(m, b) ? b : undefined) || E.elimX(m, a, xa, '[OH-]', false);
      if (!r) return null;
      return done(r.mech, { title: 'vic-Dihalide + alc. KOH → vinyl halide', why: ['Then NaNH₂ (second E2) → alkyne', 'One equivalent of alc. KOH removes only one HX'], tags: ['E2'] });
    }
    const [cc, x] = hs[0]; const d = E.deg(m, cc); const hasB = E.betaHs(m, cc, x).length > 0;
    let mode = 'S';
    if (hasB) { if (alc) mode = 'E'; else if (aq) mode = d === '3°' && hot ? 'E' : 'S'; else if (conc || hot || d === '3°') mode = 'E'; }
    const subst = () => { const o = E.REC.find(r => r.id === 'rx_oh').run(m, c); return o && o.product ? o : null; };
    if (mode === 'S') { const o = subst(); if (!o) return null;
      if (hasB && d !== '1°' && d !== 'methyl') { const e = E.elimX(m, cc, x, '[OH-]', false); if (e) o.minor = [e.mech.product()]; }
      o.why = (o.why || []).concat([aq ? 'Aqueous / dilute KOH: OH⁻ acts as a nucleophile → substitution (alcohol)' : 'Plain KOH is read as aqueous: substitution (alcohol)', 'Alcoholic or hot concentrated KOH acts as a base → elimination (alkene)']); return o; }
    const r = E.elimX(m, cc, x, '[OH-]', false, `E2: ${alc ? 'alcoholic' : 'concentrated'} OH⁻ acts as a base. It takes a β-H while the halide leaves, forming the C=C (anti-periplanar, one step).`);
    if (!r) return subst();
    if (r.antiForced) return done(r.mech, { title: 'E2: only the anti-periplanar β-H can leave → less substituted alkene', why: ['E2 needs H and X anti-periplanar: on a cyclohexane ring, trans and both axial', 'The H that would give the Saytzeff alkene is cis to the leaving group here', 'Classic case: menthyl chloride → 2-menthene only; neomenthyl chloride → mostly 3-menthene (Saytzeff)'], tags: ['E2', 'anti'] });
    const lab = alc ? 'alc. KOH' : conc ? 'conc. KOH' : hot ? 'hot KOH' : '3° halide + OH⁻'; const o = done(r.mech, { title: `${lab}: E2 → more substituted alkene (Saytzeff)`, why: [!alc && !conc && !hot ? '3° halide + strong base: E2 wins (SN2 is blocked at a 3° carbon; SN1 needs a weak nucleophile such as water)' : alc ? 'Alcoholic KOH (ethoxide/hydroxide in ethanol): strong base, poor solvation → elimination' : 'High base concentration and heat favour E2 over SN2', 'Saytzeff: the more substituted (more stable) alkene is major', 'Aqueous / dilute KOH would substitute instead (alcohol)'], tags: ['E2'] });
    if (d !== '3°') { const s = subst(); if (s) o.minor = [s.product]; }
    return o;
  } });
/* NaNH₂ (or 2 eq alc. KOH): gem/vic dihalide or vinyl halide → alkyne. The second H always comes from the other alkene carbon (no allenes). */
front({ id: 'alkyne_from_dihalide', name: 'Dihalide / vinyl halide → alkyne', fam: 'Alkyl halide', canon: 'NaNH₂', stoich: true, db: ['Dehydrohalogenation'],
  test: c => c.has('NaNH2') || (c.has('alcKOH') && ((c.eq != null && c.eq >= 2) || (c.heat && c.flag('conc')))),
  run(m, c) {
    const base = c.has('NaNH2') ? '[NH2-]' : '[OH-]'; const bn = base === '[NH2-]' ? 'NH₂⁻' : 'OH⁻';
    const second = (mm, steps0) => { const v = M(mm, '[CX3:1]([Cl,Br,I:2])=[CX3:3]').filter(([, , b]) => hTot(mm, b) > 0)[0]; if (!v) return null; const [a, x, b] = v;
      const r = E.elimX(mm, a, x, base, false, `Second E2 on the vinyl halide: ${bn} takes the H on the other alkene carbon (not an allylic H), so C=C becomes C≡C.`, b);
      if (!r) return null; const o = done(r.mech, { title: 'Double dehydrohalogenation → alkyne', why: ['The second H comes from the other alkene carbon; taking an allylic H would give an allene instead', 'The vinyl-halide E2 is slow: NaNH₂ (very strong base) is the usual reagent', 'Terminal alkynes are then deprotonated by excess NaNH₂; water work-up returns the alkyne'], tags: ['E2'] });
      o.steps = (steps0 || []).concat(o.steps); return o; };
    if (M(m, '[CX3:1]([Cl,Br,I:2])=[CX3:3]').some(([, , b]) => hTot(m, b) > 0)) return second(m, []);
    const hs = E.halides(m); if (!hs.length) return null; let first = null;
    for (const [cc, x] of hs) { const gem = hs.filter(h => h[0] === cc).length >= 2; const partner = hs.find(([b]) => b !== cc && m.bondBetween(cc, b));
      if (gem) { first = E.elimX(m, cc, x, base, false, `First E2: ${bn} takes a β-H and one halide leaves → vinyl halide.`); if (first) break; }
      if (partner && hTot(m, partner[0]) > 0) { first = E.elimX(m, cc, x, base, false, `First E2: ${bn} takes the H on the carbon that carries the other halogen → vinyl halide.`, partner[0]); if (first) break; } }
    if (!first) return null;
    return second(first.mech.product(), first.mech.steps) || done(first.mech, { title: 'E2 → vinyl halide', why: [], tags: ['E2'] });
  } });

/* Wurtz family: one halide → R–R; two different halides → mixture; aryl + alkyl → Wurtz–Fittig; 3° → fails; dihalides → rings */
function rxFrags(m) { const AR = aromSet(m); const o = [];
  for (const ids of m.frags()) { const S = new Set(ids); const h = M(m, '[#6:1][Cl,Br,I:2]').filter(([cc]) => S.has(cc) && (AR.has(cc) || m.bondsOf(cc).every(b => b.o === 1))); if (h.length) o.push({ m: m.sub(ids), hs: h, ar: AR.has(h[0][0]) }); }
  return o; }
function couple(A, a1, x1, B, a2, x2) { const t = A.clone(); t.removeAtom(x1); const u = B.clone(); u.removeAtom(x2); const idm = t.merge(u, 300, 0); t.addBond(a1, idm.get(a2), 1); clr(t, a1, idm.get(a2)); return tidy(t); }
front({ id: 'wurtz2', multi: true, name: 'Wurtz / Fittig / Wurtz–Fittig', fam: 'Alkyl halide', canon: 'Na, dry ether', db: ['Wurtz reaction', 'Fittig reaction', 'Wurtz–Fittig reaction'],
  test: c => (c.has('Mg') && !c.has('RMgX')) || (c.has('Na') && !c.any(['NH3', 'Na/NH3', 'ROH', 'H2O', 'H3O+', 'liq'])) || (c.has('Zn') && !c.any(['HCl', 'Sn/HCl', 'H2SO4', 'H3O+', 'AcOH', 'CuCN', 'CuCl', 'CuBr', 'HBr', 'H2O', 'NaI'])),
  run(m, c) {
    const Na = c.has('Na'); const Mg = c.has('Mg') && !Na && !c.has('Zn'); const F = rxFrags(m); if (!F.length) return null; const metal = Na ? 'Na' : Mg ? 'Mg' : 'Zn';
    if (F.length === 1 && !c.has('RX')) {
      const f = F[0]; const sp3 = f.hs.filter(([cc]) => !aromSet(f.m).has(cc)); const sites = [...new Map(sp3.map(h => [h[0], h])).values()];
      if (sites.length >= 2) { let best = null;
        for (let i = 0; i < sites.length; i++) for (let j = i + 1; j < sites.length; j++) { const d = distMap(f.m, sites[i][0]).get(sites[j][0]); if (d >= 2 && d <= 5 && (!best || d < best.d)) best = { d, A: sites[i], B: sites[j] }; }
        if (best && Mg && best.d !== 2) best = null;
        if (best) { const t = f.m.clone(); t.removeAtom(best.A[1]); t.removeAtom(best.B[1]); t.addBond(best.A[0], best.B[0], 1); clr(t, best.A[0], best.B[0]); const n = best.d + 1;
          const nm = ['', '', '', 'cyclopropane', 'cyclobutane', 'cyclopentane', 'cyclohexane'][n];
          return simple(m, t, { title: `Intramolecular coupling (${metal}) → ${nm}`, text: `The two C–X ends of the same chain couple with each other: ${metal} removes both halogens and joins the carbons, closing a ${n}-membered ring.`, why: [n === 3 ? '1,3-Dihalide + Zn (or Na) → cyclopropane (Freund / Gustavson)' : `1,${n}-dihalide → ${n}-membered ring`, 'Intramolecular reaction beats coupling with a second molecule at these ring sizes'], tags: ['coupling'] }); }
        if (sites.length >= 2 && sites.some(([a]) => sites.some(([b]) => b !== a && f.m.bondBetween(a, b)))) return null; // vicinal: dehalogenation recipe
      }
      if (!Na) return null; // Zn / Mg alone: only the ring closures
      const [cc, x] = f.hs[0];
      if (!f.ar && E.deg(f.m, cc) === '3°') {
        const bs = E.betaHs(f.m, cc, x); const t1 = f.m.clone(); t1.removeAtom(x); if (bs.length) setO(t1, cc, bs[0], 2); const t2 = f.m.clone(); t2.removeAtom(x); clr(t2, cc);
        const A = tidy(t1), B = tidy(t2);
        return out({ title: 'Wurtz fails for 3° halides → alkene + alkane', product: A, extra: [B], steps: [step(m, A, 'The bulky 3° radicals/carbanions do not couple; instead one takes an H from the other (disproportionation): an alkene and an alkane form.', ['Na, dry ether'])], why: ['Wurtz works for 1° (and some 2°) halides only', '3° halides: elimination / disproportionation'], tags: ['coupling'] });
      }
      const R = couple(f.m, cc, x, f.m, cc, x);
      return out({ title: f.ar ? 'Fittig reaction: 2 Ar–X + 2 Na → Ar–Ar' : 'Wurtz reaction: 2 R–X + 2 Na → R–R', product: R, steps: [step(m, R, `Na transfers electrons to the C–X bond (dry ether); two ${f.ar ? 'aryl' : 'alkyl'} groups join, NaX precipitates. One halide gives a symmetric product with twice the carbons.`, ['Na, dry ether'])], why: ['Doubles the carbon count: good only for symmetric (even-carbon) products', 'Water must be absent (it would destroy the organosodium)'], tags: ['coupling'] });
    }
    if (!Na) return null;
    const pieces = F.map(f => ({ m: f.m, cc: f.hs[0][0], x: f.hs[0][1], ar: f.ar }));
    if (pieces.length < 2 && c.has('RX')) { const rx = c.get('RX'); const g = fresh((rx.X || 'Br') + (rx.smi || 'C')); const ids = [...g.atoms.keys()]; pieces.push({ m: g, cc: ids[1], x: ids[0], ar: false, text: 1 }); }
    if (pieces.length < 2) return null;
    const [A, B] = pieces; const AB = couple(A.m, A.cc, A.x, B.m, B.cc, B.x), AA = couple(A.m, A.cc, A.x, A.m, A.cc, A.x), BB = couple(B.m, B.cc, B.x, B.m, B.cc, B.x);
    const wf = A.ar !== B.ar; const side = [AA, BB].filter(p => key(p) !== key(AB));
    return out({ title: wf ? 'Wurtz–Fittig: Ar–X + R–X + 2 Na → Ar–R' : A.ar ? 'Fittig with two aryl halides → mixture' : 'Wurtz with two different halides → mixture of three alkanes', product: AB, minor: side,
      steps: [step(m, AB, wf ? 'Na couples the aryl and the alkyl group (dry ether). The cross product Ar–R is the useful one; some R–R (Wurtz) and Ar–Ar (Fittig) also form.' : 'All three couplings compete: R–R′ (shown) plus R–R and R′–R′. Their boiling points are close, so separating them is hard.', ['Na, dry ether'])],
      why: wf ? ['Cross product is major: the aryl and alkyl organosodiums react differently, so the mixed coupling wins', 'Side products: R–R and Ar–Ar'] : ['Two different alkyl halides give three alkanes: not a clean way to make odd-carbon alkanes'], tags: ['coupling'] });
  } });

/* ===================== alkenes, dienes, alkynes, alkanes ===================== */
function dienes(m) { const AR = aromSet(m); return M(m, '[C:1]=[C:2]-[C:3]=[C:4]').filter(q => q.every(a => !AR.has(a))); }
front({ id: 'diene_add', name: 'Conjugated diene: 1,2- vs 1,4-addition', fam: 'Alkene', canon: 'HBr, 40 °C', stoich: true, db: ['1,2- and 1,4-addition'],
  test: c => c.any(['HBr', 'HCl', 'HI', 'Br2', 'Cl2']) && !c.light && !c.has('ROOR') && !c.has('FeX3') && !c.any(['H2O', 'HOBr', 'HOCl', 'ROH', 'AcOH']),
  run(m, c) {
    const D = dienes(m); if (!D.length) return null;
    let [c1, c2, c3, c4] = D[0]; if (hTot(m, c4) > hTot(m, c1)) [c1, c2, c3, c4] = [c4, c3, c2, c1];
    const hot = c.heat || (c.T != null && c.T >= 300); const hx = c.any(['HBr', 'HCl', 'HI']);
    const X = c.has('HBr') || c.has('Br2') ? 'Br' : c.has('HCl') || c.has('Cl2') ? 'Cl' : 'I';
    const make = mode => { const t = m.clone(); clr(t, c1, c2, c3, c4); setO(t, c1, c2, 1); if (!hx) addAt(t, c1, X);
      if (mode === 12) addAt(t, c2, X); else { setO(t, c3, c4, 1); setO(t, c2, c3, 2); addAt(t, c4, X); } return tidy(t); };
    const p12 = make(12), p14 = make(14); const main = hot ? p14 : p12, minor = hot ? p12 : p14;
    const why = [hot ? '≈40 °C: equilibrium → 1,4-product (more substituted internal C=C, thermodynamic)' : 'Low temperature: 1,2-product forms faster (kinetic control)', 'The allylic cation shares its + charge between C2 and C4', 'Both products form; the temperature decides which is major'];
    if (hx) { try { const mech = new Mech(m); const r = E.protonateAlkene(mech, c1, c2, `[H]${X}`);
        if (r && r.cp === c2) { if (hot) { mech.step([[[c3, c4], [c2, c3]]], 'Resonance: the allylic cation spreads its + charge onto C4.'); E.captureBy(mech, c4, r.X, `${X}⁻ attacks C4: 1,4-addition.`); }
          else E.captureBy(mech, c2, r.X, `${X}⁻ attacks C2, right next to where H⁺ added: 1,2-addition.`);
          const o = done(mech, { title: `H${X} + conjugated diene → ${hot ? '1,4' : '1,2'}-addition (${hot ? 'thermodynamic' : 'kinetic'})`, why, tags: ['AdE'] }); o.minor = [minor]; return o; } } catch (e) { }
    }
    return out({ title: `${hx ? 'H' + X : X + '₂'} + conjugated diene → ${hot ? '1,4' : '1,2'}-addition (${hot ? 'thermodynamic' : 'kinetic'})`, product: main, minor: [minor],
      steps: [step(m, main, `${hx ? 'H⁺' : X + '⁺'} adds to C1, giving an allylic cation whose + charge is shared by C2 and C4. ${X}⁻ attacks ${hot ? 'C4 (1,4-addition; the C=C moves to C2–C3)' : 'C2 (1,2-addition)'}.`, [hot ? '40 °C' : 'low T'])], why, tags: ['AdE'] });
  } });
const DPH = ['ethene', 'ethyne', 'acrolein', 'maleicanhydride', 'acrylate', 'acrylonitrile', 'quinone'];
front({ id: 'diels_alder', multi: true, name: 'Diels–Alder [4+2] cycloaddition', fam: 'Alkene', canon: 'heat', db: ['Diels–Alder reaction'],
  test: c => c.heat || c.st.sp.size === 0 || DPH.some(k => c.has(k)),
  run(m, c) {
    const D = dienes(m).filter(q => !m.bondBetween(q[0], q[3])); if (!D.length) return null;
    const [c1, c2, c3, c4] = D[0]; const own = new Set(E.fragIds(m, c1)); const AR = aromSet(m);
    let t = m.clone(); let d = M(t, '[C:1]=[C:2]').find(([a, b]) => !own.has(a) && !AR.has(a)) || M(t, '[C:1]#[C:2]').find(([a]) => !own.has(a));
    const sp = DPH.find(k => c.has(k));
    if (!d && sp) { const g = fresh(c.get(sp).smi); const idm = t.merge(g, 0, 0); const dm = M(g, '[C:1]=[C:2]')[0] || M(g, '[C:1]#[C:2]')[0]; if (dm) d = [idm.get(dm[0]), idm.get(dm[1])]; }
    if (!d) return null; const [d1, d2] = d; const tri = t.bondBetween(d1, d2).o === 3;
    setO(t, c1, c2, 1); setO(t, c3, c4, 1); setO(t, c2, c3, 2); setO(t, d1, d2, tri ? 2 : 1); t.addBond(c1, d1, 1); t.addBond(c4, d2, 1); clr(t, d1, d2);
    return simple(m, t, { title: `Diels–Alder [4+2] → ${tri ? 'cyclohexa-1,4-diene' : 'cyclohexene'} ring`, text: 'One concerted step through a cyclic transition state: the s-cis diene and the dienophile form two new σ bonds at once, and the diene\'s central bond becomes the new C=C.',
      why: ['Electron-poor dienophiles (–CHO, –COOR, –CN, anhydride) react fastest', 'Suprafacial, stereospecific: cis groups on the dienophile stay cis', 'The diene must be able to adopt the s-cis conformation'], tags: ['pericyclic'] });
  } });
front({ id: 'aromatise', name: 'Aromatisation (catalytic reforming)', fam: 'Alkane', canon: 'Cr₂O₃ / V₂O₅, 773 K, 10–20 atm', db: [],
  test: c => c.has('aromat'),
  run(m) {
    if ([...m.atoms.values()].some(a => a.el !== 'C' && a.el !== 'H') || [...m.bonds.values()].some(b => b.o !== 1) || C.rings && C.rings(m).length) return null;
    const ends = [...m.atoms.keys()].filter(a => el(m, a) === 'C'); if (ends.length < 6) return out({ title: 'Too short', product: null, none: 'Aromatisation needs a chain of six or more carbons (n-hexane → benzene).' });
    const far = a => { const d = distMap(m, a); let b = a; for (const [k, v] of d) if (v > d.get(b)) b = k; return [b, d]; };
    const [u] = far(ends[0]); const [v, du] = far(u); if (du.get(v) < 5) return out({ title: 'Chain too short', product: null, none: 'Needs six carbons in a row to close a benzene ring.' });
    const path = [v]; while (path[path.length - 1] !== u) { const cur = path[path.length - 1]; path.push(m.nbrs(cur).find(n => du.get(n) === du.get(cur) - 1)); }
    const r = path.slice(0, 6); const t = m.clone(); t.addBond(r[0], r[5], 1); setO(t, r[0], r[1], 2); setO(t, r[2], r[3], 2); setO(t, r[4], r[5], 2); r.forEach(a => clr(t, a));
    if (!C.isValid(t)) return null;
    return simple(m, t, { title: 'Aromatisation: n-alkane → arene', text: 'Over Cr₂O₃ / V₂O₅ / Mo₂O₃ at about 773 K and 10–20 atm, six carbons of the chain cyclise and lose hydrogen (dehydrocyclisation): n-hexane → benzene, n-heptane → toluene.', why: ['Also called reforming', 'Loses 4 H₂'], tags: ['alkane'] });
  } });
front({ id: 'alkyne_trimer', name: 'Cyclic polymerisation of alkynes', fam: 'Alkyne', canon: 'red-hot iron tube, 873 K', db: [],
  test: c => c.has('FeTube'),
  run(m) {
    const ay = M(m, '[C:1]#[C:2]'); if (!ay.length) return null; let [a, b] = ay[0]; if (hTot(m, b) > hTot(m, a)) [a, b] = [b, a];
    const t = m.clone(); const c2 = t.merge(m, 200, 0), c3 = t.merge(m, 400, 0);
    const A = [a, c2.get(a), c3.get(a)], B = [b, c2.get(b), c3.get(b)];
    for (let i = 0; i < 3; i++) { setO(t, A[i], B[i], 2); t.addBond(B[i], A[(i + 1) % 3], 1); clr(t, A[i], B[i]); }
    return simple(m, t, { title: 'Alkyne → benzene ring (red-hot iron tube, 873 K)', text: 'Three alkyne molecules join head-to-tail into a six-membered aromatic ring: ethyne → benzene, propyne → mesitylene (1,3,5-trimethylbenzene).', why: ['Cyclic trimerisation'], tags: ['alkyne'] });
  } });

/* ===================== alcohols, phenols, ethers ===================== */
front({ id: 'ether_413', name: 'Intermolecular dehydration → ether', fam: 'Alcohol', canon: 'conc. H₂SO₄, 413 K', db: ['Preparation of ethers by dehydration'],
  test: c => c.any(['H2SO4', 'H3PO4']) && c.T != null && c.T > 380 && c.T < 430,
  run(m) {
    const al = M(m, '[CH2X4:1][OX2H1:2]').filter(([cc]) => !aromSet(m).has(cc)); if (!al.length) return null; const [cc, o] = al[0];
    const t = m.clone(); const u = m.clone(); u.removeAtom(o); const idm = t.merge(u, 200, 0); t.addBond(o, idm.get(cc), 1); clr(t, o, idm.get(cc));
    return simple(m, t, { title: 'conc. H₂SO₄ at 413 K: 2 R–OH → R–O–R (ether)', text: 'One alcohol is protonated; a second alcohol molecule attacks its carbon (SN2) and water leaves. Loss of H⁺ gives the ether.', why: ['413 K: ether; 443 K: alkene (elimination wins at higher T)', 'Works for 1° alcohols; 2°/3° give alkenes'], tags: ['SN2'] });
  } });
front({ id: 'phenol_quinone', name: 'Phenol → p-benzoquinone', fam: 'Phenol', canon: 'Na₂Cr₂O₇, H₂SO₄', db: [],
  test: c => c.has('K2Cr2O7'),
  run(m) {
    const ph = M(m, '[c:1][OX2H1:2]'); if (!ph.length) return null; const [ca, o] = ph[0]; const ring = ringOf(m, ca); if (!ring || ring.length !== 6) return null;
    const i = ring.indexOf(ca); const R = k => ring[(i + k) % 6]; if (hTot(m, R(3)) < 1) return null;
    const t = m.clone(); setO(t, R(0), o, 2); setO(t, R(0), R(1), 1); setO(t, R(1), R(2), 2); setO(t, R(2), R(3), 1); setO(t, R(3), R(4), 1); setO(t, R(4), R(5), 2); setO(t, R(5), R(0), 1); addAt(t, R(3), 'O', 2);
    return simple(m, t, { title: 'Phenol + chromic acid → p-benzoquinone', text: 'The very electron-rich phenol ring is oxidised at the para position: two C=O groups in a cross-conjugated ring.', why: ['Phenols are easily oxidised (they darken in air for the same reason)'], tags: ['oxidation'] });
  } });
front({ id: 'cumene', name: 'Cumene process', fam: 'Phenol', canon: 'O₂, then H₃O⁺', db: ['Cumene process'],
  eats: s => s.any(['H3O+', 'H2SO4', 'H2O', 'HCl']),
  test: c => c.has('O2'),
  run(m, c, all, k) {
    const cu = M(m, '[CH1X4:1]([CX4:2])([CX4:3])[c:4]'); if (!cu.length) return null; const [ch, , , ar] = cu[0];
    const acidNext = (all && all[k + 1] && all[k + 1].any(['H3O+', 'H2SO4', 'H2O', 'HCl'])) || c.any(['H2SO4', 'H3O+', 'HCl']);
    const t = m.clone(); const o1 = addAt(t, ch, 'OO'); const hp = tidy(t);
    if (!acidNext) return out({ title: 'Air oxidation → cumene hydroperoxide', product: hp, steps: [step(m, hp, 'O₂ (air) attacks the benzylic C–H by a radical chain: cumene hydroperoxide.', ['O₂'])], why: ['Next: dilute acid rearranges it to phenol + acetone'], tags: ['radical'] });
    const u = m.clone(); cutB(u, ch, ar); const oa = addAt(u, ar, 'O'); addAt(u, ch, 'O', 2);
    return simple(m, u, { pick: f => [...f.atoms.values()].some(a => a.el === 'O') && C.aromaticRings(f).length, title: 'Cumene process → phenol + acetone', text: 'Air turns cumene into its hydroperoxide; acid then protonates the –O–OH, the phenyl group migrates from C to O as water leaves, and the cation adds water and falls apart into phenol and acetone.', more: [step(hp, tidy(u), 'H⁺: phenyl migrates to O (water leaves), then hydrolysis → phenol + acetone.', ['H₃O⁺'])], why: ['Industrial phenol; acetone is the co-product', 'Phenyl migrates (higher migratory aptitude than methyl)'], tags: ['rearrangement'] });
  } });
front({ id: 'claisen_rearr', name: 'Claisen rearrangement', fam: 'Ether', canon: 'heat (≈200 °C)', db: [],
  test: c => c.heat && c.st.sp.size === 0,
  run(m) {
    const q = M(m, '[c:5][O:1][CH2:2][CX3:3]=[CX3:4]'); if (!q.length) return null; const [ca, o, c2, c3, c4] = q[0]; const ring = ringOf(m, ca);
    const orth = m.nbrs(ca).filter(n => ring && ring.includes(n) && hTot(m, n) > 0); if (!orth.length) return null;
    const t = m.clone(); cutB(t, o, c2); setO(t, c2, c3, 2); setO(t, c3, c4, 1); t.addBond(c4, orth[0], 1); clr(t, orth[0], ca, o);
    return simple(m, t, { title: 'Claisen rearrangement → o-allylphenol', text: '[3,3]-Sigmatropic shift through a cyclic six-atom transition state: the O–CH₂ bond breaks while the far end of the allyl group bonds to the ortho carbon; the dienone then tautomerises back to the aromatic phenol.', why: ['Allyl group arrives through its terminal carbon (inversion of the allyl unit)', 'Both ortho positions blocked → it migrates on to para'], tags: ['pericyclic'] });
  } });
front({ id: 'fries', name: 'Fries rearrangement', fam: 'Phenol', canon: 'AlCl₃, heat', db: [],
  test: c => c.has('FeX3') && !c.any(['RX', 'AcCl', 'BzCl', 'Ac2O', 'RCOCl', 'Br2', 'Cl2', 'CO', 'HCl']),
  run(m, c) {
    const q = M(m, '[c:1][O:2][C:3](=O)[#6:4]'); if (!q.length) return null; const [ca, o, cc] = q[0]; const ring = ringOf(m, ca); if (!ring) return null;
    const i = ring.indexOf(ca); const para = ring[(i + 3) % 6]; const orth = [ring[(i + 1) % 6], ring[(i + 5) % 6]].filter(n => hTot(m, n) > 0);
    const mk = at => { const t = m.clone(); cutB(t, o, cc); t.addBond(cc, at, 1); clr(t, at, o); return tidy(t); };
    const hot = (c.T != null && c.T >= 413) || /high\s*t|hight|160|165|170/i.test(c.raw || ''); const pP = hTot(m, para) ? mk(para) : null, pO = orth.length ? mk(orth[0]) : null;
    const main = hot ? (pO || pP) : (pP || pO), minor = [hot ? pP : pO].filter(x => x && x !== main);
    if (!main) return null;
    return out({ title: `Fries rearrangement → ${main === pO ? 'o' : 'p'}-hydroxyaryl ketone`, product: main, minor, steps: [step(m, main, 'AlCl₃ binds the ester oxygen; the acylium ion (R–C≡O⁺) is released and re-attacks the activated ring (Friedel–Crafts acylation) ortho or para to the O.', ['AlCl₃'])], why: ['Below about 413 K (or just "Δ"): para major; high temperature (≈433 K): ortho major, held by chelation with AlCl₃', 'Both isomers form; steam distillation separates the o-isomer (intramolecular H-bond, volatile)'], tags: ['rearrangement'] });
  } });
front({ id: 'epox_roh_acid', name: 'Epoxide + alcohol (acid)', fam: 'Ether', canon: 'CH₃OH, H⁺', db: [],
  test: c => c.has('ROH') && c.any(['H3O+', 'H2SO4', 'HCl', 'pTsOH']),
  run(m, c) {
    const ep = M(m, '[C:1]1[O:2][C:3]1'); if (!ep.length) return null; let [a, o, b] = ep[0]; if (cN(m, b).length > cN(m, a).length) [a, b] = [b, a];
    const R = (c.get('ROH').smi || 'C'); const t = m.clone(); cutB(t, a, o); addAt(t, a, 'O' + R); clr(t, o);
    return simple(m, t, { title: 'Acid-catalysed epoxide opening by ROH → alkoxy alcohol', text: 'H⁺ protonates the epoxide oxygen. The C–O bond to the more substituted carbon is longer and carries more + charge, so the alcohol attacks there (SN1-like, from the back side: trans).', why: ['Acid: nucleophile at the more substituted carbon', 'Base (RO⁻): nucleophile at the less hindered carbon (SN2)'], tags: ['SN'] });
  } });
front({ id: 'diol_cyclic_ether', name: '1,4-/1,5-diol → cyclic ether', fam: 'Alcohol', canon: 'conc. H₂SO₄, Δ', db: [],
  test: c => c.any(['H2SO4', 'H3PO4']) && (c.heat || c.flag('conc')),
  run(m) {
    const oh = M(m, '[CX4:1][OX2H1:2]'); for (let i = 0; i < oh.length; i++) for (let j = i + 1; j < oh.length; j++) { const d = distMap(m, oh[i][1]).get(oh[j][0]);
      if (d === 4 || d === 5) { const t = m.clone(); t.removeAtom(oh[j][1]); t.addBond(oh[i][1], oh[j][0], 1); clr(t, oh[i][1], oh[j][0]);
        return simple(m, t, { title: `Diol + acid → ${d + 1}-membered cyclic ether`, text: 'One –OH is protonated and leaves as water; the other –OH of the same molecule attacks that carbon and closes the ring (intramolecular SN2).', why: ['Butane-1,4-diol → tetrahydrofuran (THF)', 'Five- and six-membered rings close easily'], tags: ['SN2'] }); } }
    return null;
  } });

/* ===================== aldehydes and ketones ===================== */
function aldKet(m) { return M(m, '[CX3:1](=[OX1:2])').filter(([cc, o]) => heavyN(m, cc).filter(n => n !== o).every(n => el(m, n) === 'C')); }
front({ id: 'bisulfite', name: 'Sodium bisulphite addition', fam: 'Carbonyl', canon: 'NaHSO₃', db: [],
  test: c => c.has('NaHSO3'),
  run(m) { const k = aldKet(m); if (!k.length) return null; const [cc, o] = k[0]; const ald = hTot(m, cc) > 0; const nC = cN(m, cc).length;
    if (!ald && !cN(m, cc).some(n => hTot(m, n) === 3 && cN(m, n).length === 1) && !(E.ringSizeOf && E.ringSizeOf(m, cc))) return out({ title: 'No adduct', product: null, none: 'Only aldehydes, methyl ketones and small cyclic ketones form bisulphite adducts: bulkier ketones are too hindered.' });
    const t = m.clone(); setO(t, cc, o, 1); const s = addAt(t, cc, 'S(=O)(=O)[O-]');
    return simple(m, t, { title: 'NaHSO₃ adduct (crystalline, water-soluble)', text: 'The S of HSO₃⁻ attacks the carbonyl carbon; proton transfer gives the α-hydroxy sulphonate salt.', why: ['Used to separate/purify aldehydes and methyl ketones: dilute acid or base gives the carbonyl back'], tags: ['AdN'] }); } });
front({ id: 'acetal', name: 'Acetal / ketal formation', fam: 'Carbonyl', canon: 'excess ROH, dry HCl', db: ['Acetal formation'],
  test: c => (c.has('ROH') && c.any(['HCl', 'H3O+', 'H2SO4', 'pTsOH'])) || (c.has('glycol') && c.any(['HCl', 'H3O+', 'H2SO4', 'pTsOH'])),
  run(m, c) { if (M(m, 'C(=O)[OX2H1]').length || M(m, 'C1OC1').length) return null; const k = aldKet(m); if (!k.length) return null; const [cc, o] = k[0]; const ald = hTot(m, cc) > 0;
    const t = m.clone(); t.removeAtom(o); clr(t, cc);
    if (c.has('glycol')) { const o1 = addAt(t, cc, 'O'); const c1 = addAt(t, o1, 'C'); const c2 = addAt(t, c1, 'C'); const o2 = addAt(t, c2, 'O'); t.addBond(o2, cc, 1); }
    else { const R = c.get('ROH').smi || 'C'; addAt(t, cc, 'O' + R); addAt(t, cc, 'O' + R); }
    return simple(m, t, { title: c.has('glycol') ? `Cyclic ${ald ? 'acetal' : 'ketal'} (1,3-dioxolane) with ethylene glycol` : `${ald ? 'Acetal' : 'Ketal'} with excess alcohol (dry HCl)`, text: 'Acid protonates C=O; one alcohol adds (hemiacetal); the new –OH is protonated and leaves as water; a second alcohol O attacks the oxocarbenium ion and loses H⁺.', why: ['Reversible: aqueous acid gives the carbonyl back, so acetals protect C=O groups', 'Stable to bases, Grignard reagents and hydrides'], tags: ['AdN'] }); } });
front({ id: 'benzoin', name: 'Benzoin condensation', fam: 'Carbonyl', canon: 'KCN, ethanol', db: ['Benzoin condensation'],
  test: c => c.has('NaCN') && !c.any(['H3O+', 'HCl', 'H2SO4', 'HCN']),
  run(m) { if (E.halides(m).length) return null; const q = M(m, '[c:3][CH1:1]=[O:2]'); if (!q.length) return null; const [, cc, o] = q[0];
    const t = m.clone(); const idm = t.merge(m, 250, 0); setO(t, cc, o, 1); t.addBond(cc, idm.get(cc), 1); clr(t, idm.get(cc));
    return simple(m, t, { title: 'Benzoin condensation: 2 ArCHO → ArCH(OH)COAr', text: 'CN⁻ adds to one aldehyde; the cyanohydrin anion (C–H now acidic) attacks a second aldehyde molecule; CN⁻ is expelled again (it is a catalyst).', why: ['Aromatic aldehydes only (no α-H)', 'Cyanide umpolung: the carbonyl carbon becomes a nucleophile'], tags: ['condensation'] }); } });
front({ id: 'baeyer_villiger', name: 'Baeyer–Villiger oxidation', fam: 'Carbonyl', canon: 'peracid (CH₃CO₃H / mCPBA)', db: ['Baeyer–Villiger oxidation'],
  test: c => c.has('mCPBA'),
  run(m) { if (M(m, '[C;!a]=[C;!a]').length) return null; const k = aldKet(m); if (!k.length) return null; const [cc, o] = k[0];
    const AR = aromSet(m); const g = heavyN(m, cc).filter(n => n !== o);
    const apt = n => { const nC = cN(m, n).length; if (AR.has(n)) return 2.5; return nC >= 3 ? 4 : nC === 2 ? 3 : nC === 1 ? 2 : 1; };
    const t = m.clone(); let mig;
    if (hTot(m, cc) > 0) { const oo = addAt(t, cc, 'O'); return simple(m, t, { title: 'Peracid oxidises the aldehyde → carboxylic acid', text: 'H migrates (highest migratory aptitude) in the Criegee intermediate: RCHO → RCOOH.', why: [], tags: ['oxidation'] }); }
    mig = g.sort((p, q) => apt(q) - apt(p))[0]; cutB(t, cc, mig); const ox = t.addAtom('O', (P(t, cc)[0] + P(t, mig)[0]) / 2, (P(t, cc)[1] + P(t, mig)[1]) / 2); t.addBond(cc, ox, 1); t.addBond(ox, mig, 1);
    return simple(m, t, { title: E.ringSizeOf && E.ringSizeOf(m, cc) ? 'Baeyer–Villiger: cyclic ketone → lactone (ring grows by O)' : 'Baeyer–Villiger: ketone → ester', text: 'The peracid adds to C=O (Criegee intermediate); the group with the higher migratory aptitude moves from C to O as the carboxylate leaves.', why: ['Migratory aptitude: H > 3° > 2° ≈ cyclohexyl ≈ aryl > 1° > CH₃', 'Configuration of the migrating carbon is retained'], tags: ['oxidation'] }); } });
front({ id: 'seo2', name: 'SeO₂ oxidation (α-CH₂ → C=O)', fam: 'Carbonyl', canon: 'SeO₂', db: [],
  test: c => c.has('SeO2'),
  run(m) { const k = aldKet(m); for (const [cc, o] of k) { const a = cN(m, cc).find(n => hTot(m, n) >= 2 && !isCO(m, n)); if (a == null) continue;
      const t = m.clone(); addAt(t, a, 'O', 2);
      return simple(m, t, { title: 'SeO₂: α-CH₂ / α-CH₃ → C=O (1,2-dicarbonyl)', text: 'Selenium dioxide oxidises the carbon next to C=O (via the enol): CH₃CHO → OHC–CHO (glyoxal), CH₃COCH₃ → CH₃COCHO.', why: ['Riley oxidation'], tags: ['oxidation'] }); }
    return null; } });
front({ id: 'etard', name: 'Étard / CrO₃–Ac₂O → ArCHO', fam: 'Aromatic', canon: 'CrO₂Cl₂, then H₃O⁺', db: ['Etard reaction'],
  eats: s => s.any(['H3O+', 'H2O', 'H2SO4']),
  test: c => c.has('CrO2Cl2') || (c.has('K2Cr2O7') && c.has('Ac2O')),
  run(m, c, all, k) { const q = M(m, '[c:2][CH3:1]'); if (!q.length) return null; const [, me] = q[0]; const t = m.clone(); addAt(t, me, 'O', 2); const ac = c.has('Ac2O');
    const hyd = all && all[k + 1] && all[k + 1].any(['H3O+', 'H2O', 'H2SO4']);
    if (ac && !hyd) { const u = m.clone(); addAt(u, me, 'OC(C)=O'); addAt(u, me, 'OC(C)=O'); return simple(m, u, { title: 'CrO₃ / (CH₃CO)₂O → benzylidene diacetate', text: 'The CH₃ is oxidised but trapped as the diacetate, so it cannot go on to the acid.', why: ['Next: H₃O⁺ hydrolyses it to the aldehyde'], tags: ['oxidation'] }); }
    return simple(m, t, { title: ac ? 'CrO₃ / Ac₂O then H₃O⁺ → aromatic aldehyde' : 'Étard reaction: Ar–CH₃ → Ar–CHO', text: ac ? 'Benzylidene diacetate forms first (protects the aldehyde level), then acid hydrolysis releases ArCHO.' : 'Chromyl chloride forms a brown chromium complex (Étard complex) with the methyl group; hydrolysis gives the aldehyde without over-oxidation.', why: ['KMnO₄ / K₂Cr₂O₇ alone would go all the way to the acid'], tags: ['oxidation'] }); } });
front({ id: 'gattermann_koch', name: 'Gattermann–Koch formylation', fam: 'Aromatic', canon: 'CO, HCl, anhyd. AlCl₃', db: ['Gattermann–Koch reaction'],
  test: c => c.has('CO') && (c.has('FeX3') || c.has('HCl') || /alcl3|cucl/i.test(c.raw || '')),
  run(m) { const rings = C.aromaticRings(m); if (!rings.length) return null; const ring = rings[0]; const sub = subsOnRing(m, ring);
    if (sub.some(([, s]) => NO2(m, s) || isCO(m, s))) return out({ title: 'No reaction', product: null, none: 'Gattermann–Koch needs an unreactive-to-activated ring: deactivated rings (–NO₂, –CHO) do not react.' });
    let at = ring.find(a => hTot(m, a) > 0); if (sub.length) { const i = ring.indexOf(sub[0][0]); const p = ring[(i + 3) % 6]; if (hTot(m, p)) at = p; }
    const t = m.clone(); addAt(t, at, 'C=O');
    return simple(m, t, { title: 'Gattermann–Koch: arene + CO/HCl → ArCHO', text: 'CO + HCl with AlCl₃ (CuCl) give the formyl cation [HC≡O]⁺, which substitutes the ring (para to an activating group).', why: ['Benzene → benzaldehyde, toluene → p-tolualdehyde'], tags: ['SEAr'] }); } });

/* ===================== acids and derivatives ===================== */
front({ id: 'r2cd', name: 'Acid chloride + R₂Cd → ketone', fam: 'Acid derivative', canon: '(CH₃)₂Cd', db: [],
  test: c => c.has('R2Cd'),
  run(m, c) { const q = M(m, '[CX3:1](=O)[Cl:2]'); if (!q.length) return null; const [cc, , x] = q[0]; const t = m.clone(); t.removeAtom(x); addAt(t, cc, c.get('R2Cd').smi || 'C');
    return simple(m, t, { title: 'R₂Cd + acid chloride → ketone', text: 'Dialkylcadmium is a mild organometallic: it replaces –Cl by R but does not attack the ketone formed (a Grignard would go on to the 3° alcohol).', why: ['2 RCOCl + R′₂Cd → 2 RCOR′ + CdCl₂'], tags: ['acyl substitution'] }); } });
front({ id: 'acid_nh3', name: 'Carboxylic acid + NH₃ → amide', fam: 'Acid', canon: 'NH₃, heat', db: [],
  test: c => c.has('NH3') && !c.any(['Na', 'Li', 'Na/NH3', 'NaNH2']),
  run(m, c) { const q = M(m, '[CX3:1](=O)[OX2H1:2]'); if (!q.length || E.halides(m).length) return null; const [cc, , o] = q[0];
    if (!c.heat) { const t = m.clone(); t.atoms.get(o).chg = -1; clr(t, o); const R = tidy(t); return out({ title: 'Acid + NH₃ → ammonium carboxylate (salt)', product: R, extra: [fresh('[NH4+]')], steps: [step(m, R, 'Acid–base first: NH₃ takes the proton.', ['NH₃'])], why: ['Heating the ammonium salt drives off water → amide'], tags: ['acid–base'] }); }
    const acs = M(m, '[CX3:1](=O)[OX2H1:2]'); for (let i = 0; i < acs.length; i++) for (let j = i + 1; j < acs.length; j++) { const d = distMap(m, acs[i][0]).get(acs[j][0]);
      if (d === 3 || d === 4) { const t = m.clone(); t.removeAtom(acs[i][2]); t.removeAtom(acs[j][2]); const n = addAt(t, acs[i][0], 'N'); t.addBond(n, acs[j][0], 1); clr(t, acs[j][0]);
        return simple(m, t, { title: `Dicarboxylic acid + NH₃, strong heat → cyclic imide (${d + 2}-ring)`, text: 'The diammonium salt loses water on heating: first the amide, then the second –COOH closes onto the N–H (intramolecular), giving the imide.', why: ['Succinic acid → succinimide; phthalic acid → phthalimide (Gabriel reagent)'], tags: ['condensation'] }); } }
    const t = m.clone(); t.removeAtom(o); addAt(t, cc, 'N');
    return simple(m, t, { title: 'RCOOH + NH₃, heat → RCONH₂ (amide)', text: 'NH₃ first forms the ammonium salt RCOO⁻NH₄⁺; strong heating removes water to give the amide.', why: ['Further heating with P₂O₅ dehydrates the amide to the nitrile'], tags: ['acyl substitution'] }); } });
front({ id: 'acid_heat', name: 'Heating acids: anhydrides and decarboxylation', fam: 'Acid', canon: 'heat', db: [],
  test: c => c.has('P2O5') || (c.heat && c.st.sp.size === 0),
  run(m, c) {
    const ac = M(m, '[CX3:1](=O)[OX2H1:2]').map(x => [x[0], x[2]]); if (!ac.length) return null;
    for (const [cc, o] of ac) { const bk = cN(m, cc).find(a => cN(m, a).some(b => b !== cc && isCO(m, b)) && hTot(m, a) > 0 && m.bondsOf(a).every(x => x.o === 1));
      if (bk != null && !c.has('P2O5')) { const t = m.clone(); for (const n of heavyN(m, cc)) if (n !== bk) t.removeAtom(n); t.removeAtom(cc); clr(t, bk);
        const R = tidy(t); return out({ title: 'β-Carbonyl acid, heat → loses CO₂', product: R, extra: [fresh('O=C=O')], steps: [step(m, R, 'A cyclic six-membered transition state: the C=O oxygen takes the acid H while CO₂ leaves, giving an enol that tautomerises.', ['heat'])], why: ['Malonic acid → acetic acid; acetoacetic acid → acetone', 'Needs a C=O β to the –COOH'], tags: ['decarboxylation'] }); } }
    for (let i = 0; i < ac.length; i++) for (let j = i + 1; j < ac.length; j++) { const d = distMap(m, ac[i][0]).get(ac[j][0]);
      if (d === 3 || d === 4) { const t = m.clone(); t.removeAtom(ac[j][1]); t.addBond(ac[i][1], ac[j][0], 1); clr(t, ac[i][1], ac[j][0]);
        return simple(m, t, { title: `Dicarboxylic acid, heat → cyclic anhydride (${d + 2}-ring)`, text: 'The two –COOH groups of one molecule lose water between them, closing a five- (or six-) membered anhydride ring.', why: ['Succinic, maleic, phthalic acids → their anhydrides', 'Fumaric (trans) cannot do this'], tags: ['dehydration'] }); } }
    if (c.has('P2O5')) { const [cc, o] = ac[0]; const t = m.clone(); const u = m.clone(); const uo = u.nbrs(cc).find(n => el(u, n) === 'O' && u.bondBetween(cc, n).o === 1); u.removeAtom(uo); const idm = t.merge(u, 250, 0); t.addBond(o, idm.get(cc), 1); clr(t, o, idm.get(cc));
      return simple(m, t, { title: '2 RCOOH + P₂O₅ → acid anhydride', text: 'P₂O₅ removes one water from two acid molecules: (RCO)₂O.', why: [], tags: ['dehydration'] }); }
    return null;
  } });
front({ id: 'ketoester_decarb', name: 'β-Keto ester / malonic ester: hydrolysis + decarboxylation', fam: 'Acid derivative', canon: 'H₃O⁺, heat', db: [],
  test: c => c.any(['H3O+', 'H2SO4', 'KOH', 'aqKOH', 'HCl']) && c.heat,
  run(m) { const q = M(m, '[#6:1]C(=O)[CX4:2][CX3:3](=O)[OX2:4][#6:5]'); if (!q.length) return null; const a = q[0][3], ec = q[0][4];
    const t = m.clone(); for (const n of heavyN(m, ec)) if (n !== a) for (const id of E.subtree(m, n, ec)) if (t.atoms.has(id)) t.removeAtom(id); t.removeAtom(ec); clr(t, a);
    for (const [cc, , o2, r] of M(t, '[CX3:1](=O)[OX2:2][CX4:3]')) { for (const id of E.subtree(t, r, o2)) if (t.atoms.has(id)) t.removeAtom(id); clr(t, o2); }
    const R = tidy(t); return out({ title: 'Ester hydrolysis then decarboxylation (−CO₂)', product: R, extra: [fresh('O=C=O')], steps: [step(m, R, 'Acid (or base then acid) hydrolyses the ester to the β-keto / malonic acid, which loses CO₂ on heating through a six-membered cyclic transition state.', ['H₃O⁺', 'Δ'])], why: ['Acetoacetic ester → methyl ketone; malonic ester → substituted acetic acid'], tags: ['decarboxylation'] }); } });
front({ id: 'claisen_cond', name: 'Claisen condensation', fam: 'Acid derivative', canon: 'NaOEt, then H₃O⁺', db: ['Claisen condensation'],
  test: c => c.any(['NaOEt', 'NaOMe', 'NaH']) && !c.any(['RX', 'PhCHO', 'HCHO', 'CH3CHO', 'acetone']),
  run(m) { if (E.halides(m).length) return null; const q = M(m, '[CX4;!$(C(C=O)C=O):1][CX3:2](=O)[OX2:3][CX4:4]').filter(([a]) => hTot(m, a) >= 2); if (!q.length) return null;
    const es = M(m, '[CX3:2](=O)[OX2:3][CX4:4]');
    if (es.length >= 2) { for (const [a, cc] of q) for (const [c2, , o2] of es) { if (c2 === cc) continue; const d = distMap(m, a).get(c2); if (d === 4 || d === 5) { const t = m.clone(); for (const id of E.subtree(m, o2, c2)) if (t.atoms.has(id)) t.removeAtom(id); t.addBond(a, c2, 1); clr(t, a, c2);
        return simple(m, t, { title: `Dieckmann condensation → ${d + 1}-membered cyclic β-keto ester`, text: 'Intramolecular Claisen: the enolate of one ester attacks the other ester carbonyl of the same molecule; alkoxide leaves.', why: ['1,6-diesters → 5-rings, 1,7-diesters → 6-rings'], tags: ['condensation'] }); } } }
    const [a, cc, , o, r] = q[0]; const t = m.clone(); const u = m.clone(); for (const id of E.subtree(u, o, cc)) u.removeAtom(id); const idm = t.merge(u, 250, 0); t.addBond(a, idm.get(cc), 1); clr(t, a, idm.get(cc));
    return simple(m, t, { title: 'Claisen condensation → β-keto ester', text: 'Ethoxide makes the ester enolate; it attacks the C=O of a second ester (acyl substitution) and EtO⁻ leaves. The very acidic β-keto ester is deprotonated (this drives the equilibrium); acid work-up gives it back.', why: ['Ethyl acetate → ethyl acetoacetate', 'Needs two α-H (the product must be deprotonated)'], tags: ['condensation'] }); } });

/* ===================== nitrogen ===================== */
front({ id: 'diazo_gattermann', name: 'Gattermann reaction (Cu / HX)', fam: 'Diazonium', canon: 'Cu, HBr', db: ['Gattermann reaction'],
  test: c => c.has('Cu') && c.any(['HBr', 'HCl']),
  run(m, c) { const q = M(m, '[c:1][N+:2]#[N:3]'); if (!q.length) return null; const [ca, n1, n2] = q[0]; const X = c.has('HBr') ? 'Br' : 'Cl'; const t = m.clone(); t.removeAtom(n2); t.removeAtom(n1); addAt(t, ca, X);
    return simple(m, t, { title: `Gattermann: ArN₂⁺ + Cu/H${X} → Ar–${X}`, text: `Copper powder with H${X} replaces –N₂⁺ by ${X} (N₂ is lost).`, why: ['Sandmeyer uses Cu₂X₂ instead and gives better yields'], tags: ['diazonium'] }); } });
front({ id: 'phthalimide_cleave', name: 'N-Alkylphthalimide → 1° amine', fam: 'Amine', canon: 'NaOH (aq), heat', db: ['Gabriel phthalimide synthesis'],
  test: c => c.any(['KOH', 'aqKOH', 'H3O+', 'NH2NH2', 'H2SO4']),
  run(m, c) { const q = M(m, 'O=[C:1][N:2]([CX4:4])[C:3]=O'); if (!q.length) return null; const [, c1, n, , c3] = q[0];
    const t = m.clone(); cutB(t, c1, n); cutB(t, c3, n); const hz = c.has('NH2NH2');
    if (hz) { const a = addAt(t, c1, 'N'); const b = addAt(t, c3, 'N'); t.addBond(a, b, 1); } else { const a = addAt(t, c1, c.has('H3O+') || c.has('H2SO4') ? 'O' : '[O-]'); const b = addAt(t, c3, c.has('H3O+') || c.has('H2SO4') ? 'O' : '[O-]'); }
    return simple(m, t, { pick: f => [...f.atoms.values()].some(a => a.el === 'N') && !C.aromaticRings(f).length, title: 'Hydrolysis of N-alkylphthalimide → 1° amine (Gabriel, last step)', text: hz ? 'Hydrazine attacks both C=O groups (Ing–Manske); phthalhydrazide is released with the free amine.' : 'Both amide bonds are hydrolysed: the 1° amine is released with the phthalate.', why: ['Gabriel gives pure 1° amines (no over-alkylation)', 'Aryl amines cannot be made this way (ArX does not do SN2)'], tags: ['hydrolysis'] }); } });
front({ id: 'azide_red', name: 'Alkyl azide → 1° amine', fam: 'Amine', canon: 'LiAlH₄', db: [],
  test: c => c.any(['LiAlH4', 'Sn/HCl']) || (c.has('H2') && c.has('Pd')),
  run(m) { const q = M(m, '[#6:1][N:2]=[N+:3]=[N-:4]'); if (!q.length) return null; const [cc, n1, n2, n3] = q[0]; const t = m.clone(); t.removeAtom(n3); t.removeAtom(n2); clr(t, n1);
    return simple(m, t, { title: 'R–N₃ reduction → R–NH₂', text: 'The azide is reduced; N₂ is lost and the 1° amine remains.', why: ['RX + NaN₃ then LiAlH₄: a clean route to 1° amines'], tags: ['reduction'] }); } });

/* ===================== arenes ===================== */
front({ id: 'bhc', name: 'Benzene + Cl₂, UV → BHC (addition)', fam: 'Aromatic', canon: 'Cl₂, UV light', db: [],
  test: c => c.has('Cl2') && c.light && !c.has('FeX3'),
  run(m) { const rings = C.aromaticRings(m); if (!rings.length || [...m.atoms.values()].filter(a => a.el !== 'H').length !== 6) return null; const ring = rings[0];
    const t = m.clone(); for (let i = 0; i < 6; i++) setO(t, ring[i], ring[(i + 1) % 6], 1); ring.forEach(a => addAt(t, a, 'Cl'));
    return simple(m, t, { title: 'Benzene + 3 Cl₂ (UV) → BHC (C₆H₆Cl₆)', text: 'Under UV light and without a Lewis acid, chlorine radicals add across all three π bonds: benzene hexachloride (gammexane/lindane is one stereoisomer).', why: ['With FeCl₃ instead: substitution (chlorobenzene)', 'Side-chain C–H (toluene) would be substituted instead'], tags: ['radical'] }); } });
front({ id: 'arene_h2', name: 'Arene hydrogenation', fam: 'Aromatic', canon: 'H₂, Ni, 473 K', db: [],
  test: c => c.has('H2') && c.has('Pd') && !c.has('Lindlar') && (c.heat || c.flag('pressure') || (c.T != null && c.T >= 373) || (c.eq != null && c.eq >= 3)),
  run(m) { const rings = C.aromaticRings(m); if (!rings.length) return null; const t = m.clone();
    for (const b of t.bonds.values()) { const A = t.atoms.get(b.a), B = t.atoms.get(b.b); if (A.el === 'C' && B.el === 'C' && b.o > 1) { b.o = 1; clr(t, b.a, b.b); } }
    return simple(m, t, { title: 'H₂ / Ni, heat → cyclohexane ring', text: 'Benzene adds three H₂ only under forcing conditions (Ni, ≈473 K, pressure): the aromatic stabilisation has to be paid for.', why: ['Alkenes hydrogenate at room temperature; arenes need heat and pressure'], tags: ['reduction'] }); } });
front({ id: 'arene_o3', name: 'Ozonolysis of benzene', fam: 'Aromatic', canon: 'O₃, then Zn/H₂O', db: [],
  eats: s => s.any(['Me2S', 'H2O']),
  test: c => c.has('O3'),
  run(m) { if (M(m, '[C;!a]=[C;!a]').length || M(m, 'C#C').length) return null; const rings = C.aromaticRings(m); if (!rings.length) return null; const ring = rings[0];
    const t = m.clone(); const dbl = []; for (let i = 0; i < 6; i++) { const b = t.bondBetween(ring[i], ring[(i + 1) % 6]); if (b && b.o === 2) dbl.push([ring[i], ring[(i + 1) % 6]]); }
    if (dbl.length !== 3) return null; for (const [a, b] of dbl) { cutB(t, a, b); addAt(t, a, 'O', 2); addAt(t, b, 'O', 2); }
    return simple(m, t, { title: 'Benzene + 3 O₃ (triozonide), then Zn/H₂O → 3 glyoxal', text: 'Each Kekulé C=C is cut like an alkene double bond: benzene → three OHC–CHO.', why: ['Proof of three "double bonds" in the Kekulé structure', 'o-Xylene gives glyoxal, methylglyoxal and dimethylglyoxal (both Kekulé forms)'], tags: ['oxidation'] }); } });
front({ id: 'birch', name: 'Birch reduction', fam: 'Aromatic', canon: 'Na, liq. NH₃, EtOH', db: ['Birch reduction'],
  test: c => c.has('Na/NH3') || (c.any(['Na', 'Li']) && c.has('NH3')),
  run(m) { if (M(m, 'C#C').length) return null; const rings = C.aromaticRings(m); if (!rings.length) return null; const ring = rings[0]; const sub = subsOnRing(m, ring);
    let i0 = 0; if (sub.length) { const [r, s] = sub[0]; const ewg = isCO(m, s) || NO2(m, s) || M(m, 'C#N').some(x => x.includes(s)); i0 = ring.indexOf(r) + (ewg ? 0 : 1); }
    const R = k => ring[(i0 + k) % 6]; const t = m.clone();
    for (let k = 0; k < 6; k++) setO(t, R(k), R(k + 1), 1); setO(t, R(1), R(2), 2); setO(t, R(4), R(5), 2);
    return simple(m, t, { title: 'Birch reduction → cyclohexa-1,4-diene', text: 'Na in liquid NH₃ gives solvated electrons; the ring takes one electron, the radical anion is protonated by ethanol, a second electron and proton follow. The two new sp³ carbons sit para to each other (non-conjugated diene).', why: ['Electron-donating groups (–OCH₃, –CH₃) end up on a C=C', 'Electron-withdrawing groups (–COOH) end up on an sp³ carbon'], tags: ['reduction'] }); } });

/* Reformatsky: α-bromo ester + Zn + aldehyde/ketone → β-hydroxy ester */
const CARB = { acetone: 'C(C)(C)O', CH3CHO: 'C(C)O', HCHO: 'CO', PhCHO: 'C(O)c1ccccc1' };
front({ id: 'reformatsky', name: 'Reformatsky reaction', fam: 'Acid derivative', canon: 'Zn, CH₃COCH₃ then H₃O⁺', db: ['Reformatsky reaction'],
  eats: s => s.any(['H3O+', 'H2O', 'NH4Cl']),
  test: c => c.has('Zn') && Object.keys(CARB).some(k => c.has(k)),
  run(m, c) { const q = M(m, '[Br,Cl,I:1][CX4:2][CX3:3](=O)[OX2:4][#6]'); if (!q.length) return null; const [x, ca] = q[0]; const k = Object.keys(CARB).find(z => c.has(z));
    const t = m.clone(); t.removeAtom(x); clr(t, ca); addAt(t, ca, CARB[k]);
    return simple(m, t, { title: 'Reformatsky: α-bromo ester + Zn + C=O → β-hydroxy ester', text: 'Zn inserts into the C–Br bond, giving a zinc enolate (the Reformatsky reagent). It is too weak to attack the ester, but it adds to the aldehyde/ketone C=O; H₃O⁺ work-up gives the β-hydroxy ester.', why: ['Organozinc: less reactive than a Grignard, so the ester group survives', 'Heating with acid then dehydrates it to the α,β-unsaturated ester'], tags: ['AdN'] }); } });
/* Hofmann bromamide, isolated step: N-bromoamide + OH⁻ → isocyanate */
front({ id: 'nbromoamide', name: 'N-Bromoamide + base → isocyanate', fam: 'Acid derivative', canon: 'OH⁻', db: ['Hofmann bromamide degradation'],
  test: c => c.any(['KOH', 'aqKOH', 'alcKOH']),
  run(m) { const q = M(m, '[#6:1][CX3:2](=O)[NX3H1:4][Br,Cl:5]'); if (!q.length) return null; const [r, cc, o, n, x] = q[0];
    const t = m.clone(); t.removeAtom(x); cutB(t, cc, r); t.addBond(r, n, 1); t.bondBetween(cc, n).o = 2; clr(t, cc, n, r);
    return simple(m, t, { title: 'N-Bromoamide + OH⁻ → isocyanate (Hofmann rearrangement step)', text: 'OH⁻ removes the N–H; the anion loses Br⁻ while the alkyl group moves from C to N (nitrene-like, concerted): R–N=C=O. Water then hydrolyses it to R–NH₂ + CO₂.', why: ['The migrating group keeps its configuration', 'Overall: RCONH₂ + Br₂ + 4 OH⁻ → RNH₂ + CO₃²⁻'], tags: ['rearrangement'] }); } });
/* Kolbe–Schmitt: the salt first (salicylate), acid only on work-up */
front({ id: 'kolbe2', name: 'Kolbe–Schmitt reaction', fam: 'Phenol', canon: 'NaOH, CO₂, 400 K, then H₃O⁺', db: ['Kolbe–Schmitt reaction'],
  test: c => c.has('CO2') && !c.any(['RMgX', 'RLi']),
  run(m, c) { const q = M(m, '[c:1][OX1-,OX2H1:2]'); if (!q.length) return null; const [ca, o] = q[0]; const ring = ringOf(m, ca); if (!ring) return null;
    const orth = m.nbrs(ca).filter(n => ring.includes(n) && hTot(m, n) > 0); if (!orth.length) return null;
    if (m.atoms.get(o).chg !== -1 && !c.any(['KOH', 'aqKOH', 'Na'])) return null;
    const acid = c.any(['H3O+', 'HCl', 'H2SO4']); const t = m.clone(); const oo = t.atoms.get(o); oo.chg = 0; oo.hFix = null; addAt(t, orth[0], acid ? 'C(=O)O' : 'C(=O)[O-]');
    return simple(m, t, { title: acid ? 'Kolbe–Schmitt → salicylic acid' : 'Kolbe–Schmitt → sodium salicylate', text: 'The phenoxide ring attacks CO₂ (weak electrophile) at the ortho position (≈400 K, 4–7 atm); the H moves to the phenoxide O. The product is the salicylate salt; acidifying gives salicylic acid.', why: ['Na⁺ phenoxide: ortho (chelation); K⁺ phenoxide gives more para', 'Then (CH₃CO)₂O → aspirin'], tags: ['SEAr'] }); } });

root.EngineR6 = { dienes, rxFrags };
})(typeof window !== 'undefined' ? window : globalThis);
