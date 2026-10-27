/* engine_core.js: mechanism builder for the sketch board.
   A working copy of the drawn molecule is edited by curly arrows (electron bookkeeping decides
   charges and bond orders), reagents are placed next to the reacting atoms, and every step is
   recorded as {L, arrows, R, text}. Recipes (engine_recipes.js) use these helpers. */
(function (root) {
'use strict';
const C = root.Chem;
const BL = 40;
const HAL = new Set(['F', 'Cl', 'Br', 'I']);
const el = (m, id) => m.atoms.get(id).el;
const P = (m, id) => { const a = m.atoms.get(id); return [a.x, a.y]; };
const nbrs = (m, id) => m.nbrs(id);
const heavyN = (m, id) => m.nbrs(id).filter(n => m.atoms.get(n).el !== 'H');
const hExp = (m, id) => m.nbrs(id).filter(n => m.atoms.get(n).el === 'H');
const hTot = (m, id) => m.hImp(id) + hExp(m, id).length;
const cN = (m, id) => m.nbrs(id).filter(n => m.atoms.get(n).el === 'C');
const bo = (m, a, b) => { const x = m.bondBetween(a, b); return x ? x.o : 0; };
const unit = v => { const n = Math.hypot(v[0], v[1]) || 1; return [v[0] / n, v[1] / n]; };
const rot = (v, t) => [v[0] * Math.cos(t) - v[1] * Math.sin(t), v[0] * Math.sin(t) + v[1] * Math.cos(t)];

function aromSet(m) { const s = new Set(); for (const r of C.aromaticRings(m)) for (const id of r) s.add(id); return s; }
function fragOf(m) { const f = new Map(); m.frags().forEach((ids, k) => ids.forEach(id => f.set(id, k))); return f; }
function fragIds(m, id) { const f = fragOf(m), k = f.get(id); return [...f.keys()].filter(x => f.get(x) === k); }

/* free direction at an atom (unit vector), avoiding listed atoms */
function freeDir(m, id, exclude) {
  const a = m.atoms.get(id); const ex = exclude || new Set();
  const dirs = m.nbrs(id).filter(n => !ex.has(n)).map(n => { const b = m.atoms.get(n); return Math.atan2(b.y - a.y, b.x - a.x); });
  if (!dirs.length) {
    // isolated atom: point away from the nearest other atom
    let best = null, bd = 1e9; for (const o of m.atoms.values()) { if (o.id === id || ex.has(o.id)) continue; const d = Math.hypot(o.x - a.x, o.y - a.y); if (d < bd) { bd = d; best = o; } }
    return best ? unit([a.x - best.x, a.y - best.y]) : [1, 0];
  }
  const ang = C.freeAngle(dirs, a, m); return [Math.cos(ang), Math.sin(ang)];
}
/* translate + rotate a set of atoms rigidly */
function moveSet(m, ids, dx, dy) { for (const id of ids) { const a = m.atoms.get(id); a.x += dx; a.y += dy; } }
function rotSet(m, ids, cx, cy, t) { for (const id of ids) { const a = m.atoms.get(id); const v = rot([a.x - cx, a.y - cy], t); a.x = cx + v[0]; a.y = cy + v[1]; } }
function centroid(m, ids) { let x = 0, y = 0; for (const id of ids) { const a = m.atoms.get(id); x += a.x; y += a.y; } return [x / ids.length, y / ids.length]; }

/* make one implicit H explicit (for arrows that use it) */
function expH(m, id, avoid) {
  const ex = hExp(m, id); if (ex.length) return ex[0];
  if (m.hImp(id) < 1) return null;
  const d = freeDir(m, id, avoid); const a = m.atoms.get(id); const hf = a.hFix;
  const h = m.addAtom('H', a.x + d[0] * BL * 0.82, a.y + d[1] * BL * 0.82);
  m.addBond(id, h, 1); if (hf != null) a.hFix = Math.max(0, hf - 1);
  return h;
}
/* species: tiny molecules in local coordinates (atom 0 attaches, rest extends along +x) */
function species(def) {
  // def: 'H-Br' style chain or SMILES; returns Mol with coordinates
  const sm = new C.Mol();
  if (typeof def === 'object') { def.atoms.forEach((a, i) => sm.addAtom(a.el, (a.x || 0) * BL, (a.y || 0) * BL, { id: i + 1, chg: a.chg || 0, rad: a.rad || 0, hFix: a.h != null ? a.h : null })); (def.bonds || []).forEach(([i, j, o]) => sm.addBond(i + 1, j + 1, o || 1)); return sm; }
  const p = C.parseSmiles(def);
  if (p.atoms.some(a => a.arom) || p.atoms.length > 4) {
    const f = C.fromSmiles(def, BL); if (!f) return null;
    // fromSmiles drops explicit H; keep as is
    const ids = [...f.atoms.keys()]; const c0 = P(f, ids[0]); moveSet(f, ids, -c0[0], -c0[1]); return f;
  }
  // chain layout along +x with zig-zag
  const pos = [[0, 0]]; const placed = [0]; const adjL = p.atoms.map(() => []);
  p.bonds.forEach(b => { adjL[b.a].push(b.b); adjL[b.b].push(b.a); });
  const q = [0]; const seen = new Set([0]); const ang = [0];
  while (q.length) {
    const x = q.shift(); let k = 0;
    for (const y of adjL[x]) { if (seen.has(y)) continue; seen.add(y); const base = x === 0 ? 0 : ang[x]; const t = base + (k === 0 ? (x === 0 ? 0 : (p.atoms.length > 2 ? (x % 2 ? -0.5 : 0.5) : 0)) : (k === 1 ? 1.2 : -1.2)); ang[y] = t; pos[y] = [pos[x][0] + Math.cos(t) * BL * (p.atoms[y].el === 'H' || p.atoms[x].el === 'H' ? 0.82 : 1), pos[x][1] + Math.sin(t) * BL * (p.atoms[y].el === 'H' || p.atoms[x].el === 'H' ? 0.82 : 1)]; q.push(y); k++; }
  }
  p.atoms.forEach((a, i) => { const pp = pos[i] || [i * BL, 0]; sm.addAtom(a.el, pp[0], pp[1], { id: i + 1, chg: a.chg, rad: a.rad, iso: a.iso }); });
  p.bonds.forEach(b => sm.addBond(b.a + 1, b.b + 1, Math.round(b.o)));
  return sm;
}
/* place species so its atom `att` (local id, default first) sits at distance d from site atom along its free direction */
function placeSpecies(m, sp, site, opt = {}) {
  const att = opt.att || [...sp.atoms.keys()][0];
  const d = opt.dist != null ? opt.dist : 1.75 * BL;
  const dir = opt.dir || freeDir(m, site, opt.avoid);
  const s = P(m, site); const target = [s[0] + dir[0] * d, s[1] + dir[1] * d];
  const idm = m.merge(sp, 0, 0, 'sp');
  const ids = [...idm.values()]; const a0 = idm.get(att);
  // orient: species extends away from the site
  if (ids.length > 1) {
    const c = centroid(m, ids); const ap = P(m, a0); const v = [c[0] - ap[0], c[1] - ap[1]];
    if (Math.hypot(v[0], v[1]) > 1e-3) { const t = Math.atan2(dir[1], dir[0]) - Math.atan2(v[1], v[0]); rotSet(m, ids, ap[0], ap[1], t); }
  }
  const ap = P(m, a0); moveSet(m, ids, target[0] - ap[0], target[1] - ap[1]);
  // nudge away if it collides with the substrate
  for (let k = 0; k < 6; k++) {
    let dmin = 1e9; for (const id of ids) { const a = m.atoms.get(id); for (const o of m.atoms.values()) { if (ids.includes(o.id)) continue; dmin = Math.min(dmin, Math.hypot(o.x - a.x, o.y - a.y)); } }
    if (dmin >= 0.95 * BL) break; moveSet(m, ids, dir[0] * 0.3 * BL, dir[1] * 0.3 * BL);
  }
  const map = new Map(); for (const [k, v] of idm) map.set(k, v); return map;
}

/* ---------------- electron bookkeeping ---------------- */
function applyArrows(m, arrows) {
  const dnb = new Map(), dbo = new Map(), drad = new Map();
  const inc = (M, k, v) => M.set(k, (M.get(k) || 0) + v);
  for (const ar of arrows) {
    const e = ar.f ? 1 : 2;
    if (ar.s.length === 1) { inc(dnb, ar.s[0], -e); if (ar.f) inc(drad, ar.s[0], -1); } else inc(dbo, C.key(ar.s[0], ar.s[1]), -e / 2);
    if (ar.t.length === 1) { inc(dnb, ar.t[0], e); if (ar.f) inc(drad, ar.t[0], 1); } else inc(dbo, C.key(ar.t[0], ar.t[1]), e / 2);
  }
  const dsum = new Map(); const formed = [], broken = [];
  for (const [k, v] of dbo) {
    if (Math.abs(v) < 1e-9) continue;
    const [a, b] = k.split('|').map(Number); const bd = m.bondBetween(a, b); const o = (bd ? bd.o : 0) + v;
    if (o < 0.25) { if (bd) { m.removeBond(bd.id); broken.push([a, b]); } }
    else if (bd) bd.o = Math.round(o * 2) / 2; else { m.addBond(a, b, Math.max(1, Math.round(o))); formed.push([a, b]); }
    inc(dsum, a, v); inc(dsum, b, v);
  }
  for (const id of new Set([...dnb.keys(), ...dsum.keys(), ...drad.keys()])) {
    const at = m.atoms.get(id); if (!at) continue;
    at.chg = Math.round((at.chg || 0) - (dnb.get(id) || 0) - (dsum.get(id) || 0));
    if (drad.has(id)) at.rad = Math.max(0, (at.rad || 0) + drad.get(id));
  }
  return { formed, broken };
}
/* keep drawings tidy after a step: attach incoming pieces, push leaving groups away */
function fixGeom(m, before, ch) {
  const fb = fragOf(before);
  for (const [a, b] of ch.formed) {
    if (fb.get(a) === fb.get(b)) continue;
    const fa = fragOf(m); const F = [...fa.keys()].filter(x => fa.get(x) === fa.get(a));
    const pa = F.filter(x => fb.get(x) === fb.get(a)), pb = F.filter(x => fb.get(x) === fb.get(b));
    const ha = pa.filter(x => el(m, x) !== 'H').length, hb = pb.filter(x => el(m, x) !== 'H').length;
    const [mov, p, q] = (hb < ha || (hb === ha && pb.length <= pa.length)) ? [pb, a, b] : [pa, b, a];
    const exclude = new Set(mov); const pp = P(m, p), qp0 = P(m, q);
    const len = (el(m, q) === 'H' || el(m, p) === 'H') ? 0.82 * BL : BL;
    let d = unit([qp0[0] - pp[0], qp0[1] - pp[1]]);
    const nd = m.nbrs(p).filter(n => !exclude.has(n)).map(n => unit([m.atoms.get(n).x - pp[0], m.atoms.get(n).y - pp[1]]));
    if (nd.some(v => v[0] * d[0] + v[1] * d[1] > 0.7) || Math.hypot(qp0[0] - pp[0], qp0[1] - pp[1]) < 1) d = freeDir(m, p, exclude);
    const target = [pp[0] + d[0] * len, pp[1] + d[1] * len];
    if (mov.length > 1) {
      const c = centroid(m, mov); const qp = P(m, q); const v = [c[0] - qp[0], c[1] - qp[1]];
      if (Math.hypot(v[0], v[1]) > 1e-3) rotSet(m, mov, qp[0], qp[1], Math.atan2(d[1], d[0]) - Math.atan2(v[1], v[0]));
    }
    const qp = P(m, q); moveSet(m, mov, target[0] - qp[0], target[1] - qp[1]);
  }
  const fa = fragOf(m);
  for (const [a, b] of ch.broken) {
    if (!m.atoms.has(a) || !m.atoms.has(b) || fa.get(a) === fa.get(b)) continue;
    const A = [...fa.keys()].filter(x => fa.get(x) === fa.get(a)), Bs = [...fa.keys()].filter(x => fa.get(x) === fa.get(b));
    const [mov, from, to] = Bs.length <= A.length ? [Bs, a, b] : [A, b, a];
    const v = unit([m.atoms.get(to).x - m.atoms.get(from).x, m.atoms.get(to).y - m.atoms.get(from).y]);
    moveSet(m, mov, v[0] * 0.75 * BL, v[1] * 0.75 * BL);
    for (let k = 0; k < 5; k++) { let dmin = 1e9; for (const id of mov) { const x = m.atoms.get(id); for (const o of m.atoms.values()) { if (mov.includes(o.id)) continue; dmin = Math.min(dmin, Math.hypot(o.x - x.x, o.y - x.y)); } } if (dmin > 0.95 * BL) break; moveSet(m, mov, v[0] * 0.35 * BL, v[1] * 0.35 * BL); }
  }
}
/* ring closures that came out stretched get a local relayout (RDKit), aligned to the old drawing */
function tidyLongBonds(m, lim) {
  lim = lim || 1.7;
  let bad = false; for (const b of m.bonds.values()) { const d = Math.hypot(m.atoms.get(b.a).x - m.atoms.get(b.b).x, m.atoms.get(b.a).y - m.atoms.get(b.b).y); if (d > lim * BL || d < 0.45 * BL) { bad = true; break; } }
  if (!bad) return false;
  for (const ids of m.frags()) {
    if (ids.length < 3) continue;
    let fragBad = false; for (const b of m.bonds.values()) { if (!ids.includes(b.a)) continue; const d = Math.hypot(m.atoms.get(b.a).x - m.atoms.get(b.b).x, m.atoms.get(b.a).y - m.atoms.get(b.b).y); if (d > lim * BL || d < 0.45 * BL) fragBad = true; }
    if (!fragBad) continue;
    const sub = m.sub(ids); const old = new Map(ids.map(id => [id, P(m, id)]));
    try { C.relayout(sub, BL); } catch (e) { continue; }
    // align new coordinates to the old ones (rigid fit on heavy atoms)
    const hv = ids.filter(id => el(m, id) !== 'H');
    const c0 = centroid(m, hv), c1 = centroid(sub, hv);
    let sxx = 0, sxy = 0;
    for (const id of hv) { const o = old.get(id), n = P(sub, id); const ox = o[0] - c0[0], oy = o[1] - c0[1], nx = n[0] - c1[0], ny = n[1] - c1[1]; sxx += nx * ox + ny * oy; sxy += nx * oy - ny * ox; }
    const t = Math.atan2(sxy, sxx);
    for (const id of ids) { const n = P(sub, id); const v = rot([n[0] - c1[0], n[1] - c1[1]], t); const a = m.atoms.get(id); a.x = c0[0] + v[0]; a.y = c0[1] + v[1]; }
  }
  return true;
}

/* ---------------- mechanism builder ---------------- */
class Mech {
  constructor(m) { this.m = m.clone(); this.steps = []; this.sub = new Set([...m.atoms.keys()].filter(id => m.atoms.get(id).el !== 'H')); this.notes = []; }
  H(id, avoid) { return expH(this.m, id, avoid); }
  add(def, site, opt) { const sp = (typeof def === 'string' || (def.atoms && !(def.atoms instanceof Map))) ? species(def) : def; return placeSpecies(this.m, sp, site, opt || {}); }
  drop(ids) { for (const id of ids) this.m.removeAtom(id); }
  dropFrag(id) { if (this.m.atoms.has(id)) this.drop(fragIds(this.m, id)); }
  step(arrows, text, cond, post) {
    const L = this.m.clone();
    const ar = arrows.map(a => ({ s: a[0], t: a[1], f: !!a[2] }));
    const ch = applyArrows(this.m, ar); fixGeom(this.m, L, ch); if (post) post(); else tidyLongBonds(this.m);
    this.steps.push({ L, arrows: ar, R: this.m.clone(), text: text || '', cond: cond || [] });
    return ch;
  }
  note(text, L, R) { this.steps.push({ L: L || this.m.clone(), arrows: [], R: R || null, text, cond: [] }); }
  mainFrag() {
    const f = fragOf(this.m); const cnt = new Map();
    for (const id of this.sub) if (this.m.atoms.has(id)) cnt.set(f.get(id), (cnt.get(f.get(id)) || 0) + 1);
    let best = null, bc = -1; for (const [k, v] of cnt) if (v > bc) { bc = v; best = k; }
    return [...f.keys()].filter(id => f.get(id) === best);
  }
  product(ids) { return cleanProduct(this.m, ids || this.mainFrag()); }
  others() { const keep = new Set(this.mainFrag()); const f = fragOf(this.m); const groups = new Map(); for (const [id, k] of f) if (!keep.has(id)) { if (!groups.has(k)) groups.set(k, []); groups.get(k).push(id); } return [...groups.values()]; }
}
function cleanProduct(m, ids) {
  const sub = m.sub(ids);
  for (const a of [...sub.atoms.values()]) {
    if (a.el !== 'H' || a.chg) continue; const nb = sub.nbrs(a.id); if (nb.length !== 1) continue;
    const p = sub.atoms.get(nb[0]); if (p.el === 'H') continue; sub.removeAtom(a.id);
  }
  return sub;
}
/* a direct (no-arrow) transformation on a copy: fn edits the copy */
function transform(m, fn) { const c = m.clone(); fn(c); return c; }

/* ---------------- structure perception ---------------- */
function isCarbonylC(m, id) { return el(m, id) === 'C' && m.bondsOf(id).some(b => b.o === 2 && ['O', 'S'].includes(el(m, m.other(b, id)))); }
function hasPi(m, id, exceptId) { return m.bondsOf(id).some(b => b.o >= 2 && m.other(b, id) !== exceptId); }
function substitution(m, id) { return cN(m, id).length; }
function cationScore(m, c, AR) {
  AR = AR || aromSet(m);
  if (AR.has(c)) return -6;
  let s = 0;
  for (const n of heavyN(m, c)) {
    const e = el(m, n);
    if (e === 'C') {
      s += 1;
      if (AR.has(n)) s += 2.2;
      else if (isCarbonylC(m, n)) s -= 1.8;
      else if (m.bondsOf(n).some(b => b.o >= 2 && m.other(b, n) !== c && el(m, m.other(b, n)) === 'C')) s += 1.6;
      // cyclopropylmethyl
      for (const r of C.rings(m, 3)) if (r.length === 3 && r.includes(n) && !r.includes(c)) s += 1.5;
    } else if ((e === 'O' || e === 'N') && !m.bondsOf(n).some(b => b.o === 2)) s += 3.2;
    else if (HAL.has(e)) s += 0.2;
  }
  if (m.bondsOf(c).some(b => b.o >= 2)) s -= 4; // vinyl
  return s;
}
function radicalScore(m, c, AR) {
  AR = AR || aromSet(m); if (AR.has(c)) return -6;
  let s = 0; for (const n of heavyN(m, c)) { if (el(m, n) === 'C') { s += 1; if (AR.has(n)) s += 2.4; else if (m.bondsOf(n).some(b => b.o >= 2 && m.other(b, n) !== c)) s += 1.8; } }
  return s;
}
function alkeneScore(m, a, b, AR) {
  AR = AR || aromSet(m);
  let s = cN(m, a).filter(x => x !== b).length + cN(m, b).filter(x => x !== a).length;
  for (const x of [a, b]) for (const n of heavyN(m, x)) { if (n === a || n === b) continue; if (AR.has(n) || m.bondsOf(n).some(bd => bd.o >= 2 && m.other(bd, n) !== a && m.other(bd, n) !== b)) s += 1.5; }
  return s;
}
const ringsCache = new WeakMap();
function ringSizeOf(m, a, b) { let best = 0; for (const r of C.rings(m, 8)) { if (b == null ? r.includes(a) : (r.includes(a) && r.includes(b))) { if (!best || r.length < best) best = r.length; } } return best; }

root.EngineCore = { BL, HAL, el, P, nbrs, heavyN, hExp, hTot, cN, bo, unit, rot, aromSet, fragOf, fragIds, freeDir, moveSet, rotSet, centroid, expH, species, placeSpecies, applyArrows, fixGeom, tidyLongBonds, Mech, cleanProduct, transform, isCarbonylC, hasPi, substitution, cationScore, radicalScore, alkeneScore, ringSizeOf };
})(typeof window !== 'undefined' ? window : globalThis);
