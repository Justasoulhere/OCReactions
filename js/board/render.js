/* render.js: ChemDraw-style skeletal structures (own layout, Nimbus Sans labels) + curved-arrow engine.
   Units are board pixels; bond length ST.BL. Output is SVG markup plus the geometry the editor and
   the arrow engine need (atom label boxes, glyph boxes, bond line segments). */
(function (root) {
'use strict';
const FM = root.FMET;
const ST = { BL: 40, fs: 19.5, lw: 1.6, dsp: 0.17, marg: 2.6 };
const K = () => ST.fs / 1000;
const adv = ch => (FM.g[ch] || FM.g['n'])[0];
function textW(s, fs) { let w = 0; for (const ch of s) w += adv(ch); return w * fs / 1000; }
function gbox(ch, x, base, fs) { const g = FM.g[ch]; if (!g || g[3] <= g[1]) return null; const k = fs / 1000; return [x + g[1] * k, base - g[4] * k, x + g[3] * k, base - g[2] * k]; }
const MINUS = '\u2212';
function chgText(c) { if (!c) return ''; const n = Math.abs(c); return (n > 1 ? n : '') + (c > 0 ? '+' : MINUS); }
const ELCLS = { O: 'eo', N: 'en', S: 'es', F: 'ex', Cl: 'ex', Br: 'ebr', I: 'ei', P: 'ep', B: 'eb', H: 'eh', C: 'ec' };
const METALS = new Set(['Li', 'Na', 'K', 'Mg', 'Al', 'Ca', 'Zn', 'Cu', 'Ag', 'Hg', 'Fe', 'Sn', 'Pd', 'Pt', 'Ni', 'Cr', 'Mn', 'Ba', 'Cs', 'Rb', 'Pb', 'Co', 'Ti', 'Si', 'Se']);
const elcls = el => ELCLS[el] || (METALS.has(el) ? 'em' : 'ec');
const esc = t => String(t).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const f1 = v => (Math.round(v * 10) / 10).toString();

/* ---------- which atoms are drawn, H folding ---------- */
function prepare(m, keepH) {
  // explicit H atoms with exactly one heavy neighbour fold into that atom unless an arrow uses them
  const fold = new Map(); const hidden = new Set();
  for (const a of m.atoms.values()) {
    if (a.el !== 'H' || (keepH && keepH.has(a.id)) || a.chg) continue;
    const nb = m.nbrs(a.id); if (nb.length !== 1) continue;
    const p = m.atoms.get(nb[0]); if (p.el === 'H') continue;
    if (m.bondBetween(a.id, nb[0]).o !== 1) continue;
    hidden.add(a.id); fold.set(p.id, (fold.get(p.id) || 0) + 1);
  }
  return { hidden, fold };
}
function hCount(m, a, fold) { return m.hImp(a) + (fold.get(a.id) || 0); }
function needsLabel(m, a, vis) {
  if (a.el !== 'C') return true;
  const deg = m.nbrs(a.id).filter(n => vis(n)).length;
  return deg === 0;
}
/* ---------- label layout ---------- */
function labelLayout(m, a, h, vis, fs) {
  const sym = a.el === '*' ? 'R' : a.el; const k = fs / 1000; const base = a.y + FM.cap * k / 2;
  // H side: right unless bonds point mostly right
  let sx = 0, n = 0;
  for (const nb of m.nbrs(a.id)) { if (!vis(nb)) continue; const o = m.atoms.get(nb); const d = Math.hypot(o.x - a.x, o.y - a.y) || 1; sx += (o.x - a.x) / d; n++; }
  const hLeft = h > 0 && n > 0 && sx > 0.35;
  const parts = [];
  const x0 = a.x - adv(sym[0]) * k / 2;
  const ss = fs * 0.7, subDy = fs * 0.3, supDy = -fs * 0.42;
  const symW = textW(sym, fs);
  if (!hLeft) {
    let x = x0; parts.push({ t: sym, x, y: base, fs, core: true }); x += symW;
    if (h > 0) { parts.push({ t: 'H', x, y: base, fs }); x += textW('H', fs); if (h > 1) { parts.push({ t: String(h), x, y: base + subDy, fs: ss }); x += textW(String(h), ss); } }
    const ct = chgText(a.chg) + (a.rad ? '\u2022'.repeat(a.rad) : '');
    if (ct) { parts.push({ t: ct, x: x + 0.5, y: base + supDy, fs: ss, chg: true }); }
  } else {
    let x = x0; parts.push({ t: sym, x, y: base, fs, core: true });
    const ct = chgText(a.chg) + (a.rad ? '\u2022'.repeat(a.rad) : '');
    if (ct) parts.push({ t: ct, x: x + symW + 0.5, y: base + supDy, fs: ss, chg: true });
    let xl = x0;
    if (h > 1) { const w = textW(String(h), ss); xl -= w; parts.push({ t: String(h), x: xl, y: base + subDy, fs: ss }); }
    if (h > 0) { const w = textW('H', fs); xl -= w; parts.push({ t: 'H', x: xl, y: base, fs }); }
  }
  // glyph boxes
  const glyphs = [], core = [];
  for (const p of parts) {
    let x = p.x;
    for (const ch of p.t) { const b = gbox(ch, x, p.y, p.fs); if (b) { glyphs.push(b); if (p.core) core.push(b); } x += adv(ch) * p.fs / 1000; }
  }
  let bx = [1e9, 1e9, -1e9, -1e9]; for (const b of glyphs) { bx = [Math.min(bx[0], b[0]), Math.min(bx[1], b[1]), Math.max(bx[2], b[2]), Math.max(bx[3], b[3])]; }
  return { parts, glyphs, core, box: bx, hLeft };
}
/* ---------- geometry helpers ---------- */
function rayHit(c, u, b, pad) {
  const x0 = b[0] - pad, y0 = b[1] - pad, x1 = b[2] + pad, y1 = b[3] + pad; let tmin = -1e9, tmax = 1e9;
  const lo = [x0, y0], hi = [x1, y1], cc = [c[0], c[1]], uu = [u[0], u[1]];
  for (let k = 0; k < 2; k++) {
    if (Math.abs(uu[k]) < 1e-9) { if (cc[k] < lo[k] || cc[k] > hi[k]) return null; }
    else { let t1 = (lo[k] - cc[k]) / uu[k], t2 = (hi[k] - cc[k]) / uu[k]; if (t1 > t2) [t1, t2] = [t2, t1]; tmin = Math.max(tmin, t1); tmax = Math.min(tmax, t2); }
  }
  if (tmax < Math.max(tmin, 0)) return null; return tmax;
}
function rayUnion(c, u, boxes, pad) { let r = 0; for (const b of boxes) { const h = rayHit(c, u, b, pad); if (h != null) r = Math.max(r, h); } return r; }
function segClip(p, q, labP, labQ) {
  // shorten segment p->q where it exits label glyph boxes at either end
  const d = [q[0] - p[0], q[1] - p[1]]; const L = Math.hypot(d[0], d[1]) || 1; const u = [d[0] / L, d[1] / L];
  let t0 = 0, t1 = L;
  if (labP) t0 = rayUnion(p, u, labP.glyphs, ST.marg);
  if (labQ) t1 = L - rayUnion(q, [-u[0], -u[1]], labQ.glyphs, ST.marg);
  if (t1 - t0 < 2) return null;
  return [p[0] + u[0] * t0, p[1] + u[1] * t0, p[0] + u[0] * t1, p[1] + u[1] * t1];
}
/* ---------- full layout ---------- */
function layout(m, opt = {}) {
  const fs = opt.fs || ST.fs; const BL = opt.BL || ST.BL;
  const { hidden, fold } = prepare(m, opt.keepH);
  const vis = id => !hidden.has(id);
  const atoms = new Map(); const labels = new Map();
  for (const a of m.atoms.values()) {
    if (!vis(a.id)) continue;
    const h = a.el === 'H' ? 0 : hCount(m, a, fold);
    const lab = needsLabel(m, a, vis) ? labelLayout(m, a, a.el === 'C' && !m.nbrs(a.id).some(vis) ? h : h, vis, fs) : null;
    atoms.set(a.id, { id: a.id, x: a.x, y: a.y, el: a.el, chg: a.chg, rad: a.rad, h, lab });
    if (lab) labels.set(a.id, lab);
  }
  // rings for double-bond placement
  const ringsL = root.Chem.rings(m, 8).filter(r => r.every(vis));
  const ringOf = new Map();
  for (const r of ringsL) for (let i = 0; i < r.length; i++) { const k = root.Chem.key(r[i], r[(i + 1) % r.length]); const prev = ringOf.get(k); if (!prev || r.length < prev.length || (r.length === prev.length && aromaticish(m, r) && !aromaticish(m, prev))) ringOf.set(k, r); }
  const bonds = [];
  for (const b of m.bonds.values()) {
    if (!vis(b.a) || !vis(b.b)) continue;
    const A = atoms.get(b.a), B = atoms.get(b.b); const p = [A.x, A.y], q = [B.x, B.y];
    const d = [q[0] - p[0], q[1] - p[1]]; const L = Math.hypot(d[0], d[1]) || 1; const u = [d[0] / L, d[1] / L]; const nrm = [-u[1], u[0]];
    const lines = [];
    const main = segClip(p, q, A.lab, B.lab);
    const off = (s, shrink) => {
      const pp = [p[0] + nrm[0] * s, p[1] + nrm[1] * s], qq = [q[0] + nrm[0] * s, q[1] + nrm[1] * s];
      let seg = segClip(pp, qq, A.lab, B.lab); if (!seg) return null;
      if (shrink) {
        const sa = A.lab ? 0 : shrink, sb = B.lab ? 0 : shrink;
        seg = [seg[0] + u[0] * sa, seg[1] + u[1] * sa, seg[2] - u[0] * sb, seg[3] - u[1] * sb];
      }
      return seg;
    };
    const o = Math.round(b.o);
    let wedge = null;
    if (o === 1 && (b.st === 1 || b.st === 6) && main) {
      // stereo bond drawn from the stereocentre (b.a): solid wedge or hashes widening toward b.b
      wedge = { t: b.st, p: [main[0], main[1]], q: [main[2], main[3]], w: 0.15 * BL };
    } else if (o === 1 || b.o === 1.5) { if (main) lines.push(main); }
    else if (o === 2) {
      const D = ST.dsp * BL;
      const ring = ringOf.get(root.Chem.key(b.a, b.b));
      let mode = 0;
      if (ring) {
        let cx = 0, cy = 0; for (const id of ring) { const t = m.atoms.get(id); cx += t.x; cy += t.y; } cx /= ring.length; cy /= ring.length;
        mode = ((cx - p[0]) * nrm[0] + (cy - p[1]) * nrm[1]) > 0 ? 1 : -1;
      } else {
        const na = m.nbrs(b.a).filter(x => x !== b.b && vis(x)), nb = m.nbrs(b.b).filter(x => x !== b.a && vis(x));
        let s = 0;
        for (const x of na) { const t = m.atoms.get(x); s += Math.sign((t.x - p[0]) * nrm[0] + (t.y - p[1]) * nrm[1]); }
        for (const x of nb) { const t = m.atoms.get(x); s += Math.sign((t.x - q[0]) * nrm[0] + (t.y - q[1]) * nrm[1]); }
        if ((na.length === 0 || nb.length === 0) && (na.length >= 2 || nb.length >= 2 || (na.length === 0 && nb.length === 0))) mode = 0;
        else if (A.lab && B.lab) mode = 0;
        else mode = s > 0 ? 1 : s < 0 ? -1 : 0;
      }
      if (mode === 0) { const l1 = off(D / 2), l2 = off(-D / 2); if (l1) lines.push(l1); if (l2) lines.push(l2); }
      else { if (main) lines.push(main); const l2 = off(mode * D, ring ? 0.16 * BL : 0.11 * BL); if (l2) lines.push(l2); }
    } else if (o >= 3) {
      const D = ST.dsp * BL * 0.95; if (main) lines.push(main); const l1 = off(D), l2 = off(-D); if (l1) lines.push(l1); if (l2) lines.push(l2);
    }
    bonds.push({ id: b.id, a: b.a, b: b.b, o: b.o, lines, wedge });
  }
  // charges / radicals on unlabelled carbons
  const marks = [];
  for (const A of atoms.values()) {
    if (A.lab || (!A.chg && !A.rad)) continue;
    const dirs = m.nbrs(A.id).filter(vis).map(n => { const t = m.atoms.get(n); return Math.atan2(t.y - A.y, t.x - A.x); });
    const ang = bestGap(dirs, -Math.PI / 4);
    const r = 0.36 * BL; const x = A.x + Math.cos(ang) * r, y = A.y + Math.sin(ang) * r;
    const t = chgText(A.chg) + (A.rad ? '\u2022'.repeat(A.rad) : ''); const ss = fs * 0.78; const w = textW(t, ss);
    const base = y + FM.cap * ss / 1000 / 2 - (A.chg ? 0 : 0);
    const gl = []; let xx = x - w / 2; for (const ch of t) { const bb = gbox(ch, xx, base, ss); if (bb) gl.push(bb); xx += adv(ch) * ss / 1000; }
    marks.push({ id: A.id, t, x: x - w / 2, y: base, fs: ss, glyphs: gl });
  }
  // Avoid hetero-label collisions (especially fused OH labels in cyclic products).
  // Keep the structure fixed; nudge only the later label vertically.
  const shiftLabel = (L, dy) => { for (const q of L.parts) q.y += dy; for (const q of L.glyphs) { q[1] += dy; q[3] += dy; } for (const q of L.core) { q[1] += dy; q[3] += dy; } L.box[1] += dy; L.box[3] += dy; };
  const labs = [...labels.entries()];
  for (let i = 0; i < labs.length; i++) for (let j = i + 1; j < labs.length; j++) {
    const A = labs[i][1], B = labs[j][1]; const ox = Math.min(A.box[2], B.box[2]) - Math.max(A.box[0], B.box[0]); const oy = Math.min(A.box[3], B.box[3]) - Math.max(A.box[1], B.box[1]);
    if (ox > 2 && oy > 2) { const dy = (B.box[1] >= A.box[1] ? 24 : -24); shiftLabel(B, dy); }
  }
  // bounds
  let x0 = 1e9, y0 = 1e9, x1 = -1e9, y1 = -1e9;
  const grow = (x, y) => { x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y); };
  for (const A of atoms.values()) { grow(A.x, A.y); if (A.lab) { grow(A.lab.box[0], A.lab.box[1]); grow(A.lab.box[2], A.lab.box[3]); } }
  for (const mk of marks) for (const g of mk.glyphs) { grow(g[0], g[1]); grow(g[2], g[3]); }
  if (!atoms.size) { x0 = y0 = 0; x1 = y1 = 0; }
  return { m, atoms, labels, bonds, marks, hidden, fold, box: [x0, y0, x1, y1], fs, BL };
}
function aromaticish(m, r) { let d = 0; for (let i = 0; i < r.length; i++) { const b = m.bondBetween(r[i], r[(i + 1) % r.length]); if (b && b.o === 2) d++; } return r.length === 6 && d === 3; }
function bestGap(dirs, pref) {
  if (!dirs.length) return pref;
  if (dirs.length === 1) { const d = dirs[0]; const c1 = d + 2 * Math.PI / 3, c2 = d - 2 * Math.PI / 3; return Math.cos(c1 - pref) > Math.cos(c2 - pref) ? c1 : c2; }
  const s = dirs.map(x => (x + 4 * Math.PI) % (2 * Math.PI)).sort((p, q) => p - q);
  let best = 0, bg = -1;
  for (let i = 0; i < s.length; i++) { const a0 = s[i], a1 = i + 1 < s.length ? s[i + 1] : s[0] + 2 * Math.PI; const g = a1 - a0; if (g > bg + 1e-6) { bg = g; best = a0 + g / 2; } }
  return best;
}
/* ---------- SVG ---------- */
function svgMol(L, opt = {}) {
  const out = [];
  const cls = opt.cls || 'mol';
  out.push(`<g class="${cls}">`);
  let k = 0;
  for (const b of L.bonds) {
    for (const s of b.lines) out.push(`<path class="bd" data-b="${b.id}" style="--k:${k++}" pathLength="1" d="M${f1(s[0])},${f1(s[1])}L${f1(s[2])},${f1(s[3])}"/>`);
    const w = b.wedge; if (!w) continue;
    const d = [w.q[0] - w.p[0], w.q[1] - w.p[1]]; const len = Math.hypot(d[0], d[1]) || 1; const n = [-d[1] / len, d[0] / len];
    if (w.t === 1) out.push(`<polygon class="wd" data-b="${b.id}" points="${f1(w.p[0])},${f1(w.p[1])} ${f1(w.q[0] + n[0] * w.w)},${f1(w.q[1] + n[1] * w.w)} ${f1(w.q[0] - n[0] * w.w)},${f1(w.q[1] - n[1] * w.w)}"/>`);
    else { const N = Math.max(4, Math.round(len / 5)); let dd = ''; for (let i = 1; i <= N; i++) { const t = i / N; const cx = w.p[0] + d[0] * t, cy = w.p[1] + d[1] * t, hw = w.w * t + 0.4; dd += `M${f1(cx + n[0] * hw)},${f1(cy + n[1] * hw)}L${f1(cx - n[0] * hw)},${f1(cy - n[1] * hw)}`; } out.push(`<path class="bd hs" data-b="${b.id}" d="${dd}"/>`); }
  }
  for (const A of L.atoms.values()) {
    if (!A.lab) continue;
    const c = elcls(A.el);
    out.push(`<g class="lab ${c}" data-a="${A.id}">` + A.lab.parts.map(p => `<text x="${f1(p.x)}" y="${f1(p.y)}" style="font-size:${f1(p.fs)}px">${esc(p.t)}</text>`).join('') + '</g>');
  }
  for (const mk of L.marks) out.push(`<text class="mk" data-a="${mk.id}" x="${f1(mk.x)}" y="${f1(mk.y)}" style="font-size:${f1(mk.fs)}px">${esc(mk.t)}</text>`);
  out.push('</g>');
  return out.join('');
}

/* ================= curved arrows (port of mechlib2) ================= */
const TH_PTS = [[0.3, 205], [0.6, 182], [1.0, 142], [2.0, 98], [3.0, 78], [4.5, 64], [7.0, 56]];
function theta(Ld) {
  const x = Ld / ST.BL; let d;
  if (x <= TH_PTS[0][0]) d = TH_PTS[0][1]; else if (x >= TH_PTS[TH_PTS.length - 1][0]) d = TH_PTS[TH_PTS.length - 1][1];
  else for (let i = 0; i < TH_PTS.length - 1; i++) { const [x0, y0] = TH_PTS[i], [x1, y1] = TH_PTS[i + 1]; if (x0 <= x && x <= x1) { d = y0 + (y1 - y0) * (x - x0) / (x1 - x0); break; } }
  return d * Math.PI / 180;
}
const sub = (a, b) => [a[0] - b[0], a[1] - b[1]], add = (a, b) => [a[0] + b[0], a[1] + b[1]], mul = (a, s) => [a[0] * s, a[1] * s];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1], nrm = a => Math.hypot(a[0], a[1]);
const unit = a => { const n = nrm(a) || 1; return [a[0] / n, a[1] / n]; };
class Arc {
  constructor(S, T, th, side) {
    const v = sub(T, S); const L = Math.max(1e-6, nrm(v)); const u = mul(v, 1 / L); const n = [-u[1], u[0]];
    th = Math.min(Math.max(th, 30 * Math.PI / 180), 240 * Math.PI / 180);
    const R = L / (2 * Math.sin(th / 2)); const dd = R * Math.cos(th / 2); const M = mul(add(S, T), 0.5);
    const C = sub(M, mul(n, side * dd)); const B = add(M, mul(n, side * (R - dd)));
    const a0 = Math.atan2(S[1] - C[1], S[0] - C[0]), a1 = Math.atan2(T[1] - C[1], T[0] - C[0]);
    const dp = ((a1 - a0) % (2 * Math.PI) + 2 * Math.PI) % (2 * Math.PI); let best = null;
    for (const da of [dp, dp - 2 * Math.PI]) { const am = a0 + da / 2; const P = add(C, mul([Math.cos(am), Math.sin(am)], R)); const e = nrm(sub(P, B)); if (!best || e < best[0]) best = [e, da]; }
    Object.assign(this, { S, T, C, R, a0, da: best[1], L });
  }
  pt(t) { const a = this.a0 + this.da * t; return add(this.C, mul([Math.cos(a), Math.sin(a)], this.R)); }
  tan(t) { const a = this.a0 + this.da * t; return mul([-Math.sin(a), Math.cos(a)], this.da > 0 ? 1 : -1); }
  length() { return Math.abs(this.da) * this.R; }
  tAt(s) { return Math.max(0, Math.min(1, s / Math.max(1e-6, this.length()))); }
  d(t0, t1) { const P0 = this.pt(t0), P1 = this.pt(t1); const sw = Math.abs(this.da) * (t1 - t0); return `M${f1(P0[0])},${f1(P0[1])}A${f1(this.R)},${f1(this.R)} 0 ${sw > Math.PI ? 1 : 0} ${this.da > 0 ? 1 : 0} ${f1(P1[0])},${f1(P1[1])}`; }
}
function segDist(p, a, b) { const ab = sub(b, a); const l2 = dot(ab, ab); if (l2 < 1e-9) return nrm(sub(p, a)); const t = Math.max(0, Math.min(1, dot(sub(p, a), ab) / l2)); return nrm(sub(add(a, mul(ab, t)), p)); }
class Scene {
  constructor(L) {
    this.L = L; this.m = L.m; this.xy = new Map(); for (const A of L.atoms.values()) this.xy.set(A.id, [A.x, A.y]);
    this.gly = new Map(); for (const [id, lab] of L.labels) this.gly.set(id, lab.glyphs);
    for (const mk of L.marks) { if (!this.gly.has(mk.id)) this.gly.set(mk.id, []); }
    this.blines = new Map(); for (const b of L.bonds) this.blines.set(root.Chem.key(b.a, b.b), b);
    this.bonds = L.bonds.map(b => [b.a, b.b]);
    this.pts = [];
  }
  labeled(i) { return this.L.labels.has(i); }
  rad(i, u, pad) { const g = this.L.labels.get(i); return g ? rayUnion(this.xy.get(i), u, g.glyphs, pad) : 0; }
  coreSplit(i) { const g = this.L.labels.get(i); return g ? [g.core, g.glyphs.filter(b => !g.core.includes(b))] : [[], []]; }
  bond(i, j) { return this.m.bondBetween(i, j); }
  nbrDirs(i) { const c = this.xy.get(i); return this.m.nbrs(i).filter(n => this.xy.has(n)).map(n => unit(sub(this.xy.get(n), c))); }
}
function endpoints(sc, src, tgt, fish, arc) {
  const u0 = arc.tan(0), v0 = mul(arc.tan(1), -1); let dots = null, S, T;
  if (src.length === 1) {
    const i = src[0], c = sc.xy.get(i), a = sc.m.atoms.get(i);
    if (sc.labeled(i)) {
      const r = sc.rad(i, u0, 1.0);
      const wants = !fish && (a.el !== 'C' || a.chg < 0) && a.el !== 'H';
      if (wants) {
        const [core, rest] = sc.coreSplit(i);
        const rc = rayUnion(c, u0, core, 1.0);
        const blocked = rest.some(b => { const h = rayHit(c, u0, b, 2.2); return h != null && h < rc + 13; });
        const rr = blocked ? r : rc;
        dots = [add(c, mul(u0, rr + 3.6)), [-u0[1], u0[0]], blocked];
        S = add(c, mul(u0, rr + 8.2));
      } else S = add(c, mul(u0, r + 3.5));
    } else {
      const wants = !fish && a.chg < 0;
      if (wants) { dots = [add(c, mul(u0, 6)), [-u0[1], u0[0]], false]; S = add(c, mul(u0, 10.5)); }
      else S = add(c, mul(u0, 5.0));
    }
  } else {
    const [i, j] = src; const a = sc.xy.get(i), b = sc.xy.get(j); const M = mul(add(a, b), 0.5); const e = unit(sub(b, a)); let nb = [-e[1], e[0]];
    if (dot(nb, u0) < 0) nb = mul(nb, -1);
    const bd = sc.blines.get(root.Chem.key(i, j)); let M2 = M;
    if (bd && bd.o >= 2 && bd.lines.length) {
      const mids = bd.lines.map(s => [(s[0] + s[2]) / 2, (s[1] + s[3]) / 2]); let far = mids[0], fv = -1e9;
      for (const x of mids) { const v = dot(sub(x, M), nb); if (v > fv) { fv = v; far = x; } }
      if (fv > 0.5) M2 = far;
    }
    S = add(M2, mul(nb, 2.6));
  }
  if (tgt.length === 1) {
    const i = tgt[0], c = sc.xy.get(i);
    T = add(c, mul(v0, sc.labeled(i) ? sc.rad(i, v0, 1.0) + 3.2 : 4.5));
  } else {
    const [i, j] = tgt; const bd = sc.bond(i, j);
    if (bd) {
      const a = sc.xy.get(i), b = sc.xy.get(j); const M = mul(add(a, b), 0.5); const e = unit(sub(b, a)); let nb = [-e[1], e[0]];
      if (dot(nb, v0) < 0) nb = mul(nb, -1); T = add(M, mul(nb, 3.6));
    } else {
      const sset = new Set(src); let y = null;
      if (fish) { T = mul(add(sc.xy.get(i), sc.xy.get(j)), 0.5); }
      else {
        if (sset.has(i) && !sset.has(j)) y = j; else if (sset.has(j) && !sset.has(i)) y = i;
        if (y == null) T = mul(add(sc.xy.get(i), sc.xy.get(j)), 0.5);
        else { const c = sc.xy.get(y); T = add(c, mul(v0, sc.labeled(y) ? sc.rad(y, v0, 1.0) + 3.6 : 5.0)); }
      }
    }
  }
  return [S, T, dots];
}
function anchors(sc, src, tgt, fish) {
  const S0 = src.length === 1 ? sc.xy.get(src[0]) : mul(add(sc.xy.get(src[0]), sc.xy.get(src[1])), 0.5);
  let T0;
  if (tgt.length === 1) T0 = sc.xy.get(tgt[0]);
  else {
    const [i, j] = tgt;
    if (sc.bond(i, j) || fish) T0 = mul(add(sc.xy.get(i), sc.xy.get(j)), 0.5);
    else { const s = new Set(src); if (s.has(i) && !s.has(j)) T0 = sc.xy.get(j); else if (s.has(j) && !s.has(i)) T0 = sc.xy.get(i); else T0 = mul(add(sc.xy.get(i), sc.xy.get(j)), 0.5); }
  }
  return [S0, T0];
}
function geom(sc, src, tgt, fish, side, ml) {
  const [S0, T0] = anchors(sc, src, tgt, fish);
  let arc = new Arc(S0, T0, theta(nrm(sub(T0, S0))) * ml, side), dots = null;
  for (let k = 0; k < 2; k++) { const [S, T, d] = endpoints(sc, src, tgt, fish, arc); dots = d; arc = new Arc(S, T, theta(nrm(sub(T, S))) * ml, side); }
  return [arc, dots];
}
function score(sc, arc, src, tgt, dots, ml) {
  const own = new Set([...src, ...tgt]);
  let pen = ml === 1 ? 0 : 1.4; const N = 36;
  for (let k = 0; k <= N; k++) {
    const t = k / N, p = arc.pt(t), mid = t > 0.14 && t < 0.86;
    for (const [i, c] of sc.xy) {
      if (own.has(i) && !mid) continue;
      const g = sc.gly.get(i);
      if (g && g.length) {
        let d = 99; for (const b of g) { const dx = Math.max(b[0] - 3 - p[0], 0, p[0] - b[2] - 3), dy = Math.max(b[1] - 3 - p[1], 0, p[1] - b[3] - 3); d = Math.min(d, Math.hypot(dx, dy)); }
        if (d < 3.0) pen += (3.0 - d) * 2.2 + 1.2;
      } else { const d = nrm(sub(p, c)); if (d < 6.5) pen += (6.5 - d) * 0.9; }
    }
    for (const [i, j] of sc.bonds) {
      if (!mid && (own.has(i) || own.has(j))) continue;
      const d = segDist(p, sc.xy.get(i), sc.xy.get(j)); if (d < 3.2) pen += (3.2 - d) * 1.1;
    }
    for (const q of sc.pts) { const d = nrm(sub(p, q)); if (d < 9) pen += (9 - d) * 0.9; }
  }
  if (dots) {
    if (dots[2]) pen += 7.0;
    const u = arc.tan(0);
    for (const v of sc.nbrDirs(src[0])) { const ang = Math.acos(Math.max(-1, Math.min(1, dot(u, v)))) * 180 / Math.PI; if (ang < 38) pen += (38 - ang) * 0.35; }
  }
  const B = arc.pt(0.5); let crowd = 0; for (const c of sc.xy.values()) if (nrm(sub(B, c)) < 0.9 * ST.BL) crowd++;
  return pen + crowd * 0.8;
}
function head(T, tdir, hl, hw, notch, fish, outward) {
  let p = [-tdir[1], tdir[0]];
  let pts;
  if (fish) { if (dot(p, outward) < 0) p = mul(p, -1); pts = [T, add(sub(T, mul(tdir, hl)), mul(p, hw * 1.05)), sub(T, mul(tdir, hl * notch))]; }
  else pts = [T, add(sub(T, mul(tdir, hl)), mul(p, hw)), sub(T, mul(tdir, hl * notch)), sub(sub(T, mul(tdir, hl)), mul(p, hw))];
  return pts.map(q => f1(q[0]) + ',' + f1(q[1])).join(' ');
}
const AVARS = ['a1', 'a2', 'a3', 'a4'];
/* arrow: {s:[id] | [id,id], t:[id] | [id,id], f:bool}. Returns svg string + sample points */
function drawArrow(sc, ar, i) {
  for (const k of [...ar.s, ...ar.t]) if (!sc.xy.has(k)) return ['', []];
  let best = null;
  for (const ml of [1.0, 0.78, 1.28]) for (const side of [1, -1]) {
    const [arc, dots] = geom(sc, ar.s, ar.t, ar.f, side, ml); const s = score(sc, arc, ar.s, ar.t, dots, ml);
    if (!best || s < best[0]) best = [s, arc, dots];
  }
  const [, arc, dots] = best; const Lc = arc.length();
  const hl = Math.min(10.5, 2.5 + 0.55 * arc.R, 0.45 * Lc + 2), hw = hl * 0.37, notch = 0.62;
  const tb = arc.tAt(Lc - hl); const T = arc.T; const Pb = arc.pt(tb); const tdir = unit(sub(T, Pb));
  const outward = unit(sub(T, arc.C)); const tEnd = arc.tAt(Lc - hl * notch);
  const v = AVARS[i % AVARS.length]; const o = [];
  const th0 = arc.tAt(7.5), th1 = arc.tAt(Lc - hl - 1.5);
  if (th1 > th0 + 0.04) o.push(`<path class="arx" d="${arc.d(th0, th1)}"/>`);
  o.push(`<path class="arr" style="--i:${i};stroke:var(--${v})" pathLength="1" d="${arc.d(0, tEnd)}"/>`);
  o.push(`<polygon class="arh" style="--i:${i};fill:var(--${v})" points="${head(T, tdir, hl, hw, notch, ar.f, outward)}"/>`);
  if (dots) { const [D, pn] = dots; for (const sg of [1, -1]) { const q = add(D, mul(pn, 2.7 * sg)); o.push(`<circle class="lp" style="--i:${i}" cx="${f1(q[0])}" cy="${f1(q[1])}" r="1.6"/>`); } }
  const pts = []; for (let k = 0; k <= 24; k++) pts.push(arc.pt(k / 24));
  return [o.join(''), pts];
}
function drawArrows(L, arrows) {
  const sc = new Scene(L); const out = []; let bx = L.box.slice();
  arrows.forEach((a, i) => { const [g, pts] = drawArrow(sc, a, i); out.push(g); sc.pts.push(...pts); for (const p of pts) { bx[0] = Math.min(bx[0], p[0] - 6); bx[1] = Math.min(bx[1], p[1] - 6); bx[2] = Math.max(bx[2], p[0] + 6); bx[3] = Math.max(bx[3], p[1] + 6); } });
  return { svg: `<g class="arrows">${out.join('')}</g>`, box: bx };
}
/* plain reaction arrow (straight, barbed) */
function rxnArrow(x1, y, x2, cls = 'rxa') {
  const hl = 13, hw = 4.8;
  return `<path class="${cls}" pathLength="1" d="M${f1(x1)},${f1(y)}L${f1(x2 - hl * 0.62)},${f1(y)}"/><polygon class="${cls}h" points="${f1(x2)},${f1(y)} ${f1(x2 - hl)},${f1(y - hw)} ${f1(x2 - hl * 0.62)},${f1(y)} ${f1(x2 - hl)},${f1(y + hw)}"/>`;
}
/* condensed formula text for small byproducts, e.g. H2O -> H₂O using tspans */
function formulaSpans(m, ids) {
  const cnt = {}; let chg = 0;
  for (const id of ids) { const a = m.atoms.get(id); cnt[a.el] = (cnt[a.el] || 0) + 1; chg += a.chg; if (a.el !== 'H') cnt.H = (cnt.H || 0) + m.hImp(a); }
  const order = []; if (cnt.C) { order.push('C'); if (cnt.H) order.push('H'); }
  const rest = Object.keys(cnt).filter(e => !order.includes(e) && cnt[e] > 0).sort((p, q) => (p === 'H' ? -1 : q === 'H' ? 1 : p < q ? -1 : 1));
  // inorganic convention: H first (H2O, HBr) except for metals
  const els = [...order, ...rest];
  const parts = []; for (const e of els) { if (!cnt[e]) continue; parts.push([e, 0]); if (cnt[e] > 1) parts.push([String(cnt[e]), 1]); }
  if (chg) parts.push([chgText(chg), 2]);
  return parts;
}
root.Render = { ST, layout, svgMol, drawArrows, rxnArrow, textW, esc, f1, chgText, formulaSpans, bestGap, elcls };
})(typeof window !== 'undefined' ? window : globalThis);
