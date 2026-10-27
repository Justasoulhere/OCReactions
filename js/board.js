/* Board v2 · chemistry whiteboard.
   Draw a structure, drag a reaction arrow, type reagents: products and an animated mechanism appear live.
   Model: one editable Mol (the starting structure) + a chain of reaction arrows, each solved on the
   previous species. The engine keeps atom ids, so products and mechanism steps move atom-by-atom
   from where they were instead of popping in. */
(() => {
'use strict';
const A = window.App; if (!A) return;
const { $, $$, esc, toast } = A;
const Chem = window.Chem, Render = window.Render, Solve = window.EngineSolve, Core = window.EngineCore;
const BL = 40, FS = 19, TAU = Math.PI * 2, D30 = Math.PI / 6;
const clamp = (v, a, b) => Math.max(a, Math.min(b, v));
const lerp = (a, b, t) => a + (b - a) * t;
const easeOut = t => 1 - Math.pow(1 - t, 3);
const easeIO = t => t < .5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
const f1 = v => Math.round(v * 10) / 10;
const motion = () => { try { return localStorage.getItem('ocr3-motion') !== 'reduced' && localStorage.getItem('ocr3-bxanim') !== 'off'; } catch (e) { return true; } };
const raf = fn => requestAnimationFrame(fn);

/* ---------------- chemistry text (H2SO4 → H₂SO₄ as runs) ---------------- */
const SUBD = '₀₁₂₃₄₅₆₇₈₉', SUPD = '⁰¹²³⁴⁵⁶⁷⁸⁹';
function chemRuns(s) {
  const t = String(s == null ? '' : s); const runs = []; let cur = '', mode = 0;
  const push = (ch, m) => { if (m !== mode) { if (cur) runs.push([cur, mode]); cur = ''; mode = m; } cur += ch; };
  for (let i = 0; i < t.length; i++) {
    const ch = t[i], prev = t[i - 1] || ' ', next = t[i + 1] || ' ';
    let k;
    if ((k = SUBD.indexOf(ch)) >= 0) { push(String(k), 1); continue; }
    if ((k = SUPD.indexOf(ch)) >= 0) { push(String(k), 2); continue; }
    if (ch === '⁺') { push('+', 2); continue; } if (ch === '⁻') { push('\u2212', 2); continue; }
    if (/[0-9]/.test(ch) && /[A-Za-z)\]]/.test(prev) && !/[a-z]{2}$/.test(t.slice(Math.max(0, i - 3), i)) && mode !== 2) { push(ch, 1); continue; }
    if (/[0-9]/.test(ch) && mode === 1 && /[0-9]/.test(prev)) { push(ch, 1); continue; }
    if ((ch === '+' || ch === '-') && /[A-Za-z0-9)\]]/.test(prev) && /[\s,;)/]/.test(next) && /[A-Z0-9)\]]/.test(prev)) { push(ch === '-' ? '\u2212' : '+', 2); continue; }
    push(ch, 0);
  }
  if (cur) runs.push([cur, mode]);
  return runs;
}
function chemHTML(s) { return chemRuns(s).map(([t, m]) => m === 1 ? `<sub>${esc(t)}</sub>` : m === 2 ? `<sup>${esc(t)}</sup>` : esc(t)).join(''); }
function chemSVG(s) {
  let out = '', base = 0;
  for (const [t, m] of chemRuns(s)) { const want = m === 1 ? .32 : m === 2 ? -.42 : 0; const dy = want - base; base = want; out += `<tspan${dy ? ` dy="${f1(dy * 100) / 100}em"` : ''}${m ? ' class="sm"' : ''}>${esc(t)}</tspan>`; }
  if (base) out += `<tspan dy="${f1(-base * 100) / 100}em"></tspan>`;
  return out;
}
const textW = (s, fs) => { try { return Render.textW(String(s), fs); } catch (e) { return String(s).length * fs * .55; } };
const normRg = s => A.norm(s || '').replace(/\s+/g, '');

/* ---------------- state ---------------- */
const S = {
  el: null, can: null, svg: null, cam: null, Lg: {}, // DOM + layer <g>s
  mol: null,
  tool: 'select', bondType: 'single', ringType: 'benz', atomLabel: 'O', chargeMode: '+',
  sel: { a: new Set(), b: new Set() }, hover: null, ghost: null,
  view: { x: 0, y: 0, k: 1.5 }, camT: 0,
  undo: [], redo: [],
  chain: [], active: -1, nextId: 1,
  g: null, space: false, rd: false,
  cache: new Map(), waysCache: new Map(),
  flights: new Map(), morph: null, animT: 0,
  editing: null, prodHover: -1, scheme: null, L0: null,
};

/* ---------------- molecule helpers ---------------- */
/* "Wedges & hashes" off in the View menu: draw stereo bonds as plain single bonds (the molecule itself keeps them) */
const flatW = new WeakMap();
function plainSt(m) { if (!S.stOpt || S.stOpt.wedge !== false || !m) return m; let f = flatW.get(m); if (!f || m === S.mol) { f = m.clone(); for (const b of f.bonds.values()) if (b.st === 1 || b.st === 6) b.st = 0; flatW.set(m, f); } return f; }
function lay(m) { return Render.layout(plainSt(m), { BL, fs: FS }); }
const LC = new WeakMap();
function layC(m) { let L = LC.get(m); if (!L) { L = lay(m); LC.set(m, L); } return L; }
function atomDirs(m, id, skip) { const a = m.atoms.get(id); return m.nbrs(id).filter(n => n !== skip).map(n => { const b = m.atoms.get(n); return Math.atan2(b.y - a.y, b.x - a.x); }); }
function bondProps(t) { return t === 'double' ? { o: 2, st: 0 } : t === 'triple' ? { o: 3, st: 0 } : t === 'wedge' ? { o: 1, st: 1 } : t === 'hash' ? { o: 1, st: 6 } : { o: 1, st: 0 }; }
function sproutAngle(m, id, type) {
  const a = m.atoms.get(id), nb = m.nbrs(id);
  if (!nb.length) return -D30;
  const dirs = atomDirs(m, id);
  if (nb.length === 1) {
    const b = m.bondBetween(id, nb[0]); const o = bondProps(type).o;
    if (o === 3 || b.o === 3 || (o === 2 && b.o === 2)) return dirs[0] + Math.PI;
  }
  return Chem.freeAngle(dirs, a, m);
}
function makeBond(a, b, type) { const p = bondProps(type); const id = S.mol.addBond(a, b, p.o, { st: p.st }); const bd = S.mol.bonds.get(id); if (bd) bd.st = p.st; return id; }
function modifyBond(bd, type) {
  if (type === 'single') { if (bd.st) { bd.st = 0; bd.o = 1; } else bd.o = bd.o >= 3 ? 1 : Math.round(bd.o) + 1; }
  else if (type === 'double') { bd.st = 0; bd.o = bd.o === 2 ? 1 : 2; }
  else if (type === 'triple') { bd.st = 0; bd.o = bd.o === 3 ? 1 : 3; }
  else { const st = type === 'wedge' ? 1 : 6; if (bd.o === 1 && bd.st === st) { const t = bd.a; bd.a = bd.b; bd.b = t; } else { bd.o = 1; bd.st = st; } }
}
function mergeAtoms(keep, drop) {
  const m = S.mol; if (keep === drop || !m.atoms.has(keep) || !m.atoms.has(drop)) return;
  for (const b of m.bondsOf(drop)) { const o = m.other(b, drop); if (o === keep) continue; if (!m.bondBetween(keep, o)) { const id = m.addBond(keep, o, b.o); const nb = m.bonds.get(id); nb.st = b.st; if (b.st && b.a === drop) { nb.a = keep; nb.b = o; } } }
  m.removeAtom(drop);
}
function nearAtom(x, y, r, skip) { let best = null, d = r; for (const a of S.mol.atoms.values()) { if (skip && skip.has(a.id)) continue; const q = Math.hypot(a.x - x, a.y - y); if (q < d) { d = q; best = a; } } return best; }

/* ring geometry: free placement, spiro on an atom, or fused on a bond */
function ringGeom(n, at) {
  const m = S.mol, pts = [];
  if (at.k === 'free') {
    const R = BL / (2 * Math.sin(Math.PI / n)); const off = -Math.PI / 2 + (n === 4 || n === 8 ? Math.PI / n : 0);
    for (let k = 0; k < n; k++) pts.push({ x: at.x + R * Math.cos(off + k * TAU / n), y: at.y + R * Math.sin(off + k * TAU / n) });
  } else if (at.k === 'atom') {
    const a = m.atoms.get(at.id); const phi = m.degree(at.id) ? Chem.freeAngle(atomDirs(m, at.id), a, m) : -Math.PI / 2;
    const R = BL / (2 * Math.sin(Math.PI / n)); const cx = a.x + R * Math.cos(phi), cy = a.y + R * Math.sin(phi); const st = phi + Math.PI;
    for (let k = 0; k < n; k++) pts.push({ x: cx + R * Math.cos(st + k * TAU / n), y: cy + R * Math.sin(st + k * TAU / n), id: k ? undefined : at.id });
  } else {
    const bd = m.bonds.get(at.id); const a = m.atoms.get(bd.a), b = m.atoms.get(bd.b);
    const l = Math.hypot(b.x - a.x, b.y - a.y) || BL; const R = l / (2 * Math.sin(Math.PI / n)), ap = l / (2 * Math.tan(Math.PI / n));
    const mx = (a.x + b.x) / 2, my = (a.y + b.y) / 2, ux = (b.x - a.x) / l, uy = (b.y - a.y) / l, nx = -uy, ny = ux;
    let s = 0; for (const o of m.atoms.values()) { if (o.id === a.id || o.id === b.id) continue; const dx = o.x - mx, dy = o.y - my; const dd = Math.hypot(dx, dy); if (dd > 3 * BL) continue; s += Math.sign(dx * nx + dy * ny) / Math.max(.5, dd / BL); }
    let side = s > 0 ? -1 : 1; if (at.side) side = at.side;
    const cx = mx + nx * side * ap, cy = my + ny * side * ap; const aa = Math.atan2(a.y - cy, a.x - cx);
    const tryDir = d => { const t = aa + d * (n - 1) * TAU / n; return Math.hypot(cx + R * Math.cos(t) - b.x, cy + R * Math.sin(t) - b.y); };
    const dir = tryDir(1) < tryDir(-1) ? 1 : -1;
    for (let k = 0; k < n; k++) pts.push({ x: cx + R * Math.cos(aa + dir * k * TAU / n), y: cy + R * Math.sin(aa + dir * k * TAU / n), id: k === 0 ? bd.a : k === n - 1 ? bd.b : undefined });
  }
  return pts;
}
function placeRing(n, arom, at) {
  const pts = ringGeom(n, at); const m = S.mol; const created = [];
  const ids = pts.map(p => { if (p.id != null) return p.id; const ex = nearAtom(p.x, p.y, BL * .3); if (ex) return ex.id; const id = m.addAtom('C', p.x, p.y); created.push(id); return id; });
  for (let k = 0; k < n; k++) { const a = ids[k], b = ids[(k + 1) % n]; if (a !== b && !m.bondBetween(a, b)) m.addBond(a, b, 1); }
  if (arom && n === 6) {
    const dbl = id => m.bondsOf(id).some(b => b.o === 2);
    for (let k = 0; k < n; k++) { const bd = m.bondBetween(ids[k], ids[(k + 1) % n]); if (bd && bd.o === 1 && !bd.st && !dbl(ids[k]) && !dbl(ids[(k + 1) % n])) bd.o = 2; }
  }
  return created;
}
/* zig-zag chain from a start point toward the cursor (30° snapping) */
function chainPts(sx, sy, cx, cy, fromId) {
  const d = Math.hypot(cx - sx, cy - sy); let th = Math.round(Math.atan2(cy - sy, cx - sx) / D30) * D30;
  const step = BL * Math.cos(D30); const n = clamp(Math.round(d / step), 1, 40);
  const ux = Math.cos(th), uy = Math.sin(th), nx = -uy, ny = ux;
  let side = -1;
  if (fromId != null && S.mol.degree(fromId)) { const nb = S.mol.atoms.get(S.mol.nbrs(fromId)[0]); const p1 = [sx + ux * step - nx * BL / 2, sy + uy * step - ny * BL / 2], p2 = [sx + ux * step + nx * BL / 2, sy + uy * step + ny * BL / 2]; side = Math.hypot(p1[0] - nb.x, p1[1] - nb.y) > Math.hypot(p2[0] - nb.x, p2[1] - nb.y) ? -1 : 1; }
  const pts = []; for (let k = 1; k <= n; k++) { const o = k % 2 ? side * BL / 2 : 0; pts.push({ x: sx + ux * step * k + nx * o, y: sy + uy * step * k + ny * o }); }
  return pts;
}

/* ---------------- labels: elements + group abbreviations ---------------- */
const ABBR = {
  oh: 'O', sh: 'S', nh2: 'N', cn: 'C#N', nc: '[N+]#[C-]', no2: '[N+](=O)[O-]', no: 'N=O', cho: 'C=O', cooh: 'C(=O)O', co2h: 'C(=O)O', 'coo-': 'C(=O)[O-]', coona: 'C(=O)O[Na]',
  cocl: 'C(=O)Cl', conh2: 'C(=O)N', coome: 'C(=O)OC', co2me: 'C(=O)OC', cooch3: 'C(=O)OC', cooet: 'C(=O)OCC', co2et: 'C(=O)OCC', cooc2h5: 'C(=O)OCC', coch3: 'C(=O)C', come: 'C(=O)C', ac: 'C(=O)C',
  oac: 'OC(C)=O', ococh3: 'OC(C)=O', ome: 'OC', och3: 'OC', oet: 'OCC', oc2h5: 'OCC', oph: 'Oc1ccccc1', me: 'C', et: 'CC', c2h5: 'CC', pr: 'CCC', npr: 'CCC', c3h7: 'CCC', ipr: 'C(C)C',
  bu: 'CCCC', nbu: 'CCCC', tbu: 'C(C)(C)C', sbu: 'C(C)CC', ibu: 'CC(C)C', ph: 'c1ccccc1', c6h5: 'c1ccccc1', bn: 'Cc1ccccc1', bz: 'C(=O)c1ccccc1', vinyl: 'C=C', allyl: 'CC=C',
  cf3: 'C(F)(F)F', ccl3: 'C(Cl)(Cl)Cl', so3h: 'S(=O)(=O)O', so2cl: 'S(=O)(=O)Cl', mgbr: '[Mg]Br', mgcl: '[Mg]Cl', mgi: '[Mg]I', 'n2+': '[N+]#N', n2: '[N+]#N', nhme: 'NC', nhch3: 'NC', nme2: 'N(C)C',
  nhac: 'NC(C)=O', ts: 'S(=O)(=O)c1ccc(C)cc1', ots: 'OS(=O)(=O)c1ccc(C)cc1', ch2oh: 'CO', ch2cl: 'CCl', ch2br: 'CBr', ona: 'O[Na]', ok: 'O[K]', 'o-': '[O-]', ccoh: 'C(C)=O', cono: 'C(=O)N'
};
const ABBR_SHOW = ['OH', 'NH2', 'NO2', 'CN', 'CHO', 'COOH', 'COCl', 'CONH2', 'COOMe', 'OMe', 'OEt', 'OAc', 'Me', 'Et', 'iPr', 'tBu', 'Ph', 'Bn', 'Ac', 'CF3', 'CCl3', 'SO3H', 'MgBr', 'N2+', 'NMe2', 'NHAc', 'OTs', 'ONa', 'NC'];
const ELEMS = ['C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'Si', 'Li', 'Na', 'K', 'Mg', 'Al', 'Zn', 'Cu', 'Hg', 'Sn', 'Se'];
function parseEl(t) {
  const mm = t.match(/^([A-Za-z][a-z]?)(?:H(\d*))?(\d?[+-]{1,2}|[+-]\d)?$/); if (!mm) return null;
  const el = mm[1][0].toUpperCase() + (mm[1][1] || '').toLowerCase(); if (!Chem.Z[el]) return null;
  let chg = 0; const c = mm[3]; if (c) { const sg = c.includes('-') ? -1 : 1; const dg = c.replace(/[+-]/g, ''); chg = sg * (dg ? +dg : (c.match(/[+-]/g) || []).length); }
  return { el, chg };
}
function groupMol(text) {
  const k = text.trim().toLowerCase().replace(/[₀-₉]/g, ch => SUBD.indexOf(ch)).replace('⁺', '+').replace('⁻', '-');
  if (ABBR[k]) return Chem.fromSmilesMB(ABBR[k], BL);
  return null;
}
function applyLabel(id, text, at) {
  text = String(text || '').trim(); if (!text) return false;
  const m = S.mol; let a = id != null ? m.atoms.get(id) : null;
  const G = groupMol(text); const el = G ? null : parseEl(text);
  let fragM = null;
  if (!G && !el) { try { fragM = Chem.fromSmilesMB(text, BL); } catch (e) { fragM = null; } if (!fragM || !fragM.atoms.size) { toast(`“${text}” isn't an element or a group I know`); return false; } }
  pushUndo();
  if (!a) { const nid = m.addAtom('C', at.x, at.y); a = m.atoms.get(nid); id = nid; }
  if (el) { a.el = el.el; a.chg = el.chg; a.hFix = null; a.rad = 0; return true; }
  expandGroup(id, G || fragM); return true;
}
function expandGroup(id, G) {
  const m = S.mol, a = m.atoms.get(id); const ga = [...G.atoms.values()]; const g0 = ga[0];
  let rot = 0;
  if (ga.length > 1) {
    let cx = 0, cy = 0; for (const t of ga.slice(1)) { cx += t.x; cy += t.y; } cx /= ga.length - 1; cy /= ga.length - 1;
    const gAng = Math.atan2(cy - g0.y, cx - g0.x); const tAng = m.degree(id) ? Chem.freeAngle(atomDirs(m, id), a, m) : 0; rot = tAng - gAng;
  }
  const c = Math.cos(rot), s = Math.sin(rot);
  a.el = g0.el; a.chg = g0.chg || 0; a.rad = g0.rad || 0; a.hFix = Chem.valList(g0.el, g0.chg) ? null : g0.hFix;
  const idm = new Map([[g0.id, id]]);
  for (const t of ga.slice(1)) { const dx = t.x - g0.x, dy = t.y - g0.y; idm.set(t.id, m.addAtom(t.el, a.x + dx * c - dy * s, a.y + dx * s + dy * c, { chg: t.chg, rad: t.rad, hFix: t.hFix })); }
  for (const b of G.bonds.values()) { const bid = m.addBond(idm.get(b.a), idm.get(b.b), Math.round(b.o) || 1); m.bonds.get(bid).st = 0; }
}

/* ---------------- analysis: formula, mass, DBE, valence ---------------- */
const MASS = { H: 1.008, B: 10.81, C: 12.011, N: 14.007, O: 15.999, F: 18.998, Na: 22.99, Mg: 24.305, Al: 26.982, Si: 28.085, P: 30.974, S: 32.06, Cl: 35.45, K: 39.098, Ca: 40.078, Br: 79.904, I: 126.904, Li: 6.94, Zn: 65.38, Cu: 63.546, Hg: 200.59, Ag: 107.868, Sn: 118.71, Fe: 55.845, Cr: 51.996, Mn: 54.938, Se: 78.971 };
function analysis(m) {
  if (!m || !m.atoms.size) return null;
  const cnt = {}; let chg = 0;
  for (const a of m.atoms.values()) { cnt[a.el] = (cnt[a.el] || 0) + 1; chg += a.chg || 0; if (a.el !== 'H') cnt.H = (cnt.H || 0) + m.hImp(a); }
  const keys = Object.keys(cnt).filter(k => cnt[k] > 0); const order = [];
  if (cnt.C) { order.push('C'); if (cnt.H) order.push('H'); } order.push(...keys.filter(k => !order.includes(k)).sort());
  const formula = order.map(k => k + (cnt[k] > 1 ? cnt[k] : '')).join('') + (chg ? (Math.abs(chg) > 1 ? Math.abs(chg) : '') + (chg > 0 ? '+' : '−') : '');
  let mw = 0, known = true; for (const k of keys) { if (MASS[k] == null) known = false; else mw += MASS[k] * cnt[k]; }
  let dbe = null; const X = (cnt.F || 0) + (cnt.Cl || 0) + (cnt.Br || 0) + (cnt.I || 0);
  if (keys.every(k => 'C H N O S F Cl Br I'.split(' ').includes(k))) dbe = ((2 * (cnt.C || 0) + 2 + (cnt.N || 0) - (cnt.H || 0) - X) / 2);
  const frags = m.frags().length;
  const bad = [...m.atoms.keys()].filter(id => { try { return m.valenceError(id); } catch (e) { return false; } });
  return { formula, mw: known ? mw : null, dbe, frags, bad, chg };
}
function badText(m, bad) {
  if (!bad.length) return '';
  const a = m.atoms.get(bad[0]); const n = m.bondSum(a.id);
  const nm = { C: 'Carbon', N: 'Nitrogen', O: 'Oxygen', S: 'Sulfur', F: 'Fluorine', Cl: 'Chlorine', Br: 'Bromine', I: 'Iodine', H: 'Hydrogen', P: 'Phosphorus', B: 'Boron' }[a.el] || a.el;
  return `${nm}${a.chg ? (a.chg > 0 ? ' cation' : ' anion') : ''} with ${n} bond${n === 1 ? '' : 's'}${bad.length > 1 ? ` (+${bad.length - 1} more)` : ''}`;
}

/* ---------------- history ---------------- */
/* an arrow reacts one source: { f: [atom ids] } = the drawn compound(s) holding those atoms, { a: arrowId } = that arrow's product */
function arrowObj(reagent, id, src) { return { id: id || ('r' + (S.nextId++)), reagent: reagent || '', src: src ? JSON.parse(JSON.stringify(src)) : null, res: null, key: '', parts: [], x1: 0, x2: 0, y: 0, srcMol: null, conn: null }; }
function snapshot() { return { mol: S.mol.clone(), chain: S.chain.map(c => ({ id: c.id, reagent: c.reagent, src: c.src ? JSON.parse(JSON.stringify(c.src)) : null, forceLib: !!c.forceLib })), active: S.active }; }
function pushUndo() { S.undo.push(snapshot()); if (S.undo.length > 120) S.undo.shift(); S.redo.length = 0; }
function restoreSnap(s) {
  S.mol = s.mol; S.chain = s.chain.map(c => { const a = arrowObj(c.reagent, c.id, c.src); a.forceLib = c.forceLib; return a; }); S.active = clamp(s.active, -1, S.chain.length - 1);
  S.sel.a.clear(); S.sel.b.clear(); closeEditors(true); changed({ noFlight: true });
}
function undo() { const s = S.undo.pop(); if (!s) { toast('Nothing to undo'); return; } S.redo.push(snapshot()); restoreSnap(s); }
function redo() { const s = S.redo.pop(); if (!s) { toast('Nothing to redo'); return; } S.undo.push(snapshot()); restoreSnap(s); }

/* every edit funnels through here: redraw now, re-solve shortly */
function changed(o = {}) { S._noFlight = !!o.noFlight; drawAll(); scheduleSolve(o.fast ? 60 : 170); scheduleSmiles(); }

/* ---------------- solving the chain ---------------- */
let solveTimer = 0;
function scheduleSolve(ms = 170) { clearTimeout(solveTimer); solveTimer = setTimeout(solveAll, ms); }
function molKey(m) { const { smi, order } = Chem.toSmiles(m); return smi + '#' + order.join(','); }
function solveOne(src, rg) {
  let sk = ''; if (window.Stereo && [...src.bonds.values()].some(b => b.st === 1 || b.st === 6)) { try { sk = window.Stereo.iso(src) || ''; } catch (e) { sk = ''; } }
  const key = molKey(src) + '|' + sk + '|' + normRg(rg);
  let r = S.cache.get(key);
  if (!r) {
    try { r = Solve.solve(src, rg); } catch (e) { r = { status: 'error', reason: 'The engine tripped on this one: ' + (e && e.message) }; }
    S.cache.set(key, r); if (S.cache.size > 400) S.cache.delete(S.cache.keys().next().value);
  }
  return [r, key];
}
/* lanes: which compound / product each arrow acts on */
const arrowById = id => S.chain.find(a => a.id === id);
function fragHit(src, fr) { const want = new Set(src.f); return fr.filter(c => c.some(id => want.has(id))); }
function fragAt(id) { if (id == null || !S.mol.atoms.has(id)) return null; return S.mol.frags().find(c => c.includes(id)) || null; }
/* the arrow's source if none was set: the compound the user is working on */
function defaultSrc(i) {
  const prev = i > 0 ? S.chain[i - 1] : null; if (prev) return { a: prev.id };
  return { f: [...S.mol.atoms.keys()] };
}
function descendants(id) { const out = new Set([id]); let grew = true; while (grew) { grew = false; for (const a of S.chain) if (a.src && a.src.a != null && out.has(a.src.a) && !out.has(a.id)) { out.add(a.id); grew = true; } } return out; }
function pruneOrphans() {
  /* arrows whose compound was erased disappear together with everything downstream of them */
  if (!S.mol.atoms.size) return false; const fr = S.mol.frags(); let gone = new Set();
  for (const a of S.chain) if (a.src && a.src.f && !fragHit(a.src, fr).length) for (const d of descendants(a.id)) gone.add(d);
  for (const a of S.chain) if (a.src && a.src.a != null && !arrowById(a.src.a)) for (const d of descendants(a.id)) gone.add(d);
  if (!gone.size) return false; S.chain = S.chain.filter(a => !gone.has(a.id)); S.active = Math.min(S.active, S.chain.length - 1); return true;
}
function solveAll() {
  if (!S.rd || !S.el) return;
  pruneOrphans();
  layoutScheme();
  const before = S.chain.map(ar => ar.res);
  const srcPos = S.chain.map(ar => srcWorld(ar));
  const an = analysis(S.mol); const fr = S.mol.atoms.size ? S.mol.frags() : [];
  S.chain.forEach((ar, i) => { if (!ar.src) ar.src = defaultSrc(i); });
  for (let i = 0; i < S.chain.length; i++) {
    const ar = S.chain[i]; const rg = (ar.reagent || '').trim(); let src = null; ar.srcMol = null;
    if (ar.src.a != null) {
      const p = arrowById(ar.src.a); const pi = S.chain.indexOf(p);
      if (!p || pi > i || !p.res || p.res.status !== 'ok') { ar.res = { status: 'nosrc', parent: pi }; ar.key = ''; continue; }
      src = (p.st && p.st.disp) || p.res.product;
    } else {
      const hit = fragHit(ar.src, fr);
      if (!hit.length) { ar.res = { status: 'nomol' }; ar.key = ''; continue; }
      const ids = hit.flat(); ar.src.f = ids;
      const bad = an ? an.bad.filter(id => ids.includes(id)) : [];
      if (bad.length) { ar.res = { status: 'invalid', reason: badText(S.mol, bad) }; ar.key = ''; continue; }
      src = ids.length === S.mol.atoms.size ? S.mol : S.mol.sub(ids);
    }
    ar.srcMol = src;
    if (!rg) { ar.res = { status: 'empty' }; ar.key = ''; continue; }
    let [r, key] = solveOne(src, rg);
    if (ar.forceLib || (r.status !== 'ok' && r.status !== 'empty')) { const fk = 'nb|' + key; let fb = S.cache.get(fk); if (fb === undefined) { try { fb = libFallback(src, rg); } catch (e) { fb = null; } S.cache.set(fk, fb); } if (fb && fb.noRx) { if (r.status === 'none') r = Object.assign({}, r, { agree: fb.id }); } else if (fb) r = fb;
      if (r.status !== 'ok' && !(fb && fb.noRx) && !(r.status === 'none' && r.positive)) { const ak = 'an|' + key; let an = S.cache.get(ak); if (an === undefined) { try { an = analogyFallback(src, rg); } catch (e) { an = null; } S.cache.set(ak, an); } if (an) r = an; } }
    if (r.status === 'ok' && r.product && !r._al) { try { alignTo(r.product, src, r.extra); } catch (e) { } try { const St = window.Stereo; if (St) for (const m of [r.product, ...(r.extra || [])].filter(x => x && x.bonds)) { const c = St.cip(m); for (const b of m.bonds.values()) if ((b.st === 1 || b.st === 6) && !c.all.has(b.a)) b.st = 0; } } catch (e) { } r._al = true; }
    ar.res = r; ar.key = key; ar.st = r.status === 'ok' ? stereoFor(ar) : null;
  }
  layoutScheme();
  if (!S._noFlight && motion()) S.chain.forEach((ar, i) => {
    const r = ar.res; if (!r || r.status !== 'ok' || r === before[i]) return;
    const prevSmi = before[i] && before[i].product ? Chem.toSmiles(before[i].product).smi : null;
    if (prevSmi && prevSmi === Chem.toSmiles(r.product).smi) return;
    if (r.fromLib) return; S.flights.set(ar.id, { from: srcPos[i], rg: [(ar.x1 + ar.x2) / 2, ar.y - 12], t0: performance.now(), dur: 760 });
  });
  S._noFlight = false;
  if (S.flights.size) loop();
  drawScheme(); drawHL(); placeOverlays(); renderPanel(); updateStatus();
  ensureSchemeVisible();
}
/* which stereoisomers this arrow gives (cached per reaction + drawn configuration) */
function stereoFor(ar) {
  const St = window.Stereo; const r = ar.res; if (!St || !r || r.status !== 'ok' || !ar.srcMol || r.fromLib) return null;
  let sk = ''; try { sk = St.iso(ar.srcMol) || ''; } catch (e) { sk = ''; }
  const par = ar.src && ar.src.a != null ? arrowById(ar.src.a) : null; const racSrc = !!(par && par.st && (par.st.kind === 'racemic' || par.st.racSrc));
  const k = ar.key + '|' + sk + '|' + racSrc; if (!S.stCache) S.stCache = new Map(); if (S.stCache.has(k)) return S.stCache.get(k);
  let o = null; try { o = St.outcome(ar.srcMol, r, ar.reagent); } catch (e) { o = null; }
  if (o) {
    if (racSrc && o.kind === 'single') { const m2 = St.mirror(o.list[0].mol); o = Object.assign({}, o, { kind: 'racemic', racSrc: true, list: [o.list[0], { mol: m2, smi: St.iso(m2), labels: St.cip(m2) }], note: 'The starting material is racemic, so this product forms as a racemic pair too: each enantiomer of the start gives one enantiomer of the product.' }); }
    const one = o.list && o.list[0] && o.list[0].mol; const hasW = one && [...one.bonds.values()].some(b => b.st === 1 || b.st === 6);
    if (hasW && (o.kind === 'meso' || o.kind === 'single' || o.kind === 'racemic')) o.disp = one;
    if (o.kind === 'achiral' && !o.note) o = null;
  }
  S.stCache.set(k, o); if (S.stCache.size > 300) S.stCache.delete(S.stCache.keys().next().value); return o;
}
const ST_BADGE = { meso: 'meso', racemic: 'racemic (±)', single: 'one stereoisomer', diastereomers: n => n + ' stereoisomers', 'pair-meso': 'meso + (±) pair' };
function stBadge(st) { if (!st || !st.kind) return ''; const b = ST_BADGE[st.kind]; if (!b) return ''; return typeof b === 'function' ? b(st.list.length) : b; }
/* rigid 2D fit of a product (and its by-products) onto the source drawing: same atom ids → least motion.
   A mirror fit is allowed when it is clearly better; wedges then swap with hashes so the molecule itself is unchanged. */
function fit2D(P, Q) {
  const n = P.length; if (n < 2) return null; let px = 0, py = 0, qx = 0, qy = 0;
  for (let k = 0; k < n; k++) { px += P[k][0]; py += P[k][1]; qx += Q[k][0]; qy += Q[k][1]; } px /= n; py /= n; qx /= n; qy /= n;
  const best = [];
  for (const refl of [1, -1]) { let sxx = 0, sxy = 0; for (let k = 0; k < n; k++) { const ax = (P[k][0] - px) * refl, ay = P[k][1] - py, bx = Q[k][0] - qx, by = Q[k][1] - qy; sxx += ax * bx + ay * by; sxy += ax * by - ay * bx; }
    const th = Math.atan2(sxy, sxx), c = Math.cos(th), s = Math.sin(th); let err = 0;
    for (let k = 0; k < n; k++) { const ax = (P[k][0] - px) * refl, ay = P[k][1] - py; const x = ax * c - ay * s + qx, y = ax * s + ay * c + qy; err += (x - Q[k][0]) ** 2 + (y - Q[k][1]) ** 2; }
    best.push({ refl, c, s, err, px, py, qx, qy }); }
  const [r0, r1] = best; return r1.err < r0.err * .55 - 1 ? r1 : r0;
}
function alignTo(prod, src, extra) {
  const ids = [...prod.atoms.keys()].filter(id => src.atoms.has(id) && prod.atoms.get(id).el !== 'H'); if (ids.length < 2) return;
  const P = ids.map(id => { const a = prod.atoms.get(id); return [a.x, a.y]; }), Q = ids.map(id => { const a = src.atoms.get(id); return [a.x, a.y]; });
  const T = fit2D(P, Q); if (!T) return; const { refl, c, s, px, py, qx, qy } = T;
  const id0 = Math.abs(c - 1) < 1e-3 && refl === 1; let moved = 0; for (let k = 0; k < ids.length; k++) moved += Math.hypot(P[k][0] - Q[k][0], P[k][1] - Q[k][1]);
  if (id0 && moved / ids.length < 2) return; // already drawn on the source geometry
  for (const m of [prod, ...(extra || []).filter(e => e && e.atoms)]) {
    for (const a of m.atoms.values()) { const ax = (a.x - px) * refl, ay = a.y - py; a.x = ax * c - ay * s + qx; a.y = ax * s + ay * c + qy; }
    if (refl < 0) for (const b of m.bonds.values()) { if (b.st === 1) b.st = 6; else if (b.st === 6) b.st = 1; }
  }
}
/* world positions of an arrow's source atoms (drawn compound or the parent's main product) */
function srcWorld(ar) {
  const mp = new Map(); if (!ar) return mp;
  if (ar.src && ar.src.a != null) { const p = arrowById(ar.src.a); const q = p && p.parts && p.parts.find(x => x.main); if (q) for (const a of q.m.atoms.values()) mp.set(a.id, [a.x + q.dx, a.y + q.dy]); return mp; }
  const want = ar.src && ar.src.f ? new Set(ar.src.f) : null;
  for (const a of S.mol.atoms.values()) if (!want || want.has(a.id)) mp.set(a.id, [a.x, a.y]);
  return mp;
}
function speciesWorld(i) { if (i < 0) return srcWorld({ src: null }); const ar = S.chain[i]; const mp = new Map(); const p = ar && ar.parts && ar.parts.find(q => q.main); if (p) for (const a of p.m.atoms.values()) mp.set(a.id, [a.x + p.dx, a.y + p.dy]); return mp; }

/* ---------------- scheme layout (world coordinates) ---------------- */
function shortTitle(t) { t = String(t || '').replace(/\s*\(major isomer\)/, ''); return t.length > 44 ? t.slice(0, 42).replace(/[\s,;]+\S*$/, '') + '…' : t; }
function layoutScheme() {
  const has = S.mol && S.mol.atoms.size;
  S.L0 = has ? lay(S.mol) : null;
  const fr = has ? S.mol.frags() : [];
  /* one box per drawn compound (atom positions + room for labels) */
  S.fragBox = fr.map(ids => { let x0 = 1e9, y0 = 1e9, x1 = -1e9, y1 = -1e9; for (const id of ids) { const a = S.mol.atoms.get(id); x0 = Math.min(x0, a.x); y0 = Math.min(y0, a.y); x1 = Math.max(x1, a.x); y1 = Math.max(y1, a.y); } return { ids, box: [x0 - 18, y0 - 16, x1 + 18, y1 + 16] }; });
  const b0 = S.L0 ? S.L0.box : [-30, -20, 30, 20];
  const box = b0.slice(); const occ = S.fragBox.map(f => f.box);
  const grow = b => { box[0] = Math.min(box[0], b[0]); box[1] = Math.min(box[1], b[1]); box[2] = Math.max(box[2], b[2]); box[3] = Math.max(box[3], b[3]); };
  const hitAny = (r, skip) => occ.some(o => o !== skip && r[0] < o[2] && r[2] > o[0] && r[1] < o[3] && r[3] > o[1]);
  const kids = new Map();
  S.chain.forEach((ar, i) => {
    if (!ar.src) ar.src = defaultSrc(i);
    let sb = null;
    if (ar.src.a != null) { const p = arrowById(ar.src.a); const q = p && p.parts && p.parts.find(x => x.main); sb = q ? q.box : p ? (p.slot || p.box) : null; }
    else { const want = new Set(ar.src.f); const bs = S.fragBox.filter(f => f.ids.some(id => want.has(id))).map(f => f.box); if (bs.length) sb = [Math.min(...bs.map(b => b[0])), Math.min(...bs.map(b => b[1])), Math.max(...bs.map(b => b[2])), Math.max(...bs.map(b => b[3]))]; }
    if (!sb) sb = b0;
    const kk = ar.src.a != null ? 'a' + ar.src.a : 'f' + Math.min(...ar.src.f); const k = kids.get(kk) || 0; kids.set(kk, k + 1);
    const r = ar.res; const lw = textW(ar.reagent || 'reagent', 15.5); const tw = r && r.status === 'ok' ? textW(shortTitle(r.title), 11.5) : 0;
    const len = clamp(Math.max(lw, tw) + 46, 136, 360);
    /* measure the product side first, then find a free row */
    const sp = r && r.status === 'ok' && r.product ? [(ar.st && ar.st.disp) || r.product, ...(r.extra || []).filter(e => e && e.atoms && e.atoms.size)] : null;
    const Ls = sp ? sp.map(m => layC(m)) : null;
    let wide = 0, up = 36, dn = 36;
    if (Ls) { Ls.forEach((L, j) => { const b = L.box; wide += (j ? 34 : 0) + (b[2] - b[0]) + 18; up = Math.max(up, (b[3] - b[1]) / 2 + 4); dn = Math.max(dn, (b[3] - b[1]) / 2 + 4); }); }
    else wide = (r && (r.status === 'none' || r.status === 'unknown' || r.status === 'invalid' || r.status === 'error') ? 178 : 112) + 30;
    const cy = (sb[1] + sb[3]) / 2; const x1 = sb[2] + 30; let y = k ? sb[3] + up + 26 : cy;
    const region = yy => [x1 - 6, yy - Math.max(up, 36), x1 + len + 30 + wide, yy + Math.max(dn, 24)];
    let guard = 0; while (hitAny(region(y), sb) && guard++ < 80) y += 18;
    ar.x1 = x1; ar.x2 = x1 + len; ar.y = y; ar.parts = []; ar.slot = null;
    ar.conn = Math.abs(y - cy) > 1 ? (k ? { x: (sb[0] + sb[2]) / 2, y0: sb[3] - 2, y1: y, x1 } : { x: sb[2] + 10, y0: cy, y1: y, x1 }) : null;
    let px = ar.x2 + 30;
    if (Ls) {
      Ls.forEach((L, j) => {
        if (j) { ar.parts.push({ plus: true, x: px + 4, y }); px += 34; }
        const b = L.box; const dx = px - b[0], dy = y - (b[1] + b[3]) / 2;
        const pb = [b[0] + dx, b[1] + dy, b[2] + dx, b[3] + dy];
        ar.parts.push({ m: sp[j], L, dx, dy, box: pb, main: j === 0 }); px = pb[2] + 18; grow(pb);
      });
      const mb = ar.parts.find(p => p.main).box; ar.badge = stBadge(ar.st); ar.box = [ar.x1, Math.min(mb[1], y - 30), px, Math.max(mb[3] + (ar.badge ? 20 : 0), y + 22)];
    } else {
      const w = wide - 30; ar.slot = [px, y - 36, px + w, y + 36]; grow(ar.slot); ar.box = [ar.x1, y - 36, ar.slot[2], y + 36];
    }
    if (ar.conn) grow([ar.conn.x - 4, Math.min(ar.conn.y0, y), ar.x1, Math.max(ar.conn.y0, y)]);
    occ.push(ar.box); grow([ar.x1, y - 34, ar.x2, y + 22]);
  });
  S.scheme = { box, cy: (b0[1] + b0[3]) / 2 };
}

/* ---------------- camera ---------------- */
function cw() { return S.can.clientWidth; } function ch() { return S.can.clientHeight; }
function toWorld(cx, cy) { const r = S.can.getBoundingClientRect(); return [(cx - r.left - cw() / 2 - S.view.x) / S.view.k, (cy - r.top - ch() / 2 - S.view.y) / S.view.k]; }
function toScreen(wx, wy) { return [wx * S.view.k + cw() / 2 + S.view.x, wy * S.view.k + ch() / 2 + S.view.y]; }
function applyCam() {
  if (!S.cam) return; const v = S.view;
  S.cam.setAttribute('transform', `translate(${f1(cw() / 2 + v.x)},${f1(ch() / 2 + v.y)}) scale(${v.k.toFixed(4)})`);
  const g = 24 * v.k; let gs = g; while (gs < 14) gs *= 2;
  S.can.style.backgroundSize = `${gs}px ${gs}px`; S.can.style.backgroundPosition = `${f1(cw() / 2 + v.x)}px ${f1(ch() / 2 + v.y)}px`;
  const z = S.el.querySelector('[data-z]'); if (z) z.textContent = Math.round(v.k / 1.5 * 100) + '%';
  placeOverlays();
}
function animateCam(to, ms = 520) {
  cancelAnimationFrame(S.camT);
  if (!motion()) { Object.assign(S.view, to); applyCam(); return; }
  const from = { ...S.view }, t0 = performance.now();
  const tick = now => { const t = clamp((now - t0) / ms, 0, 1), e = easeIO(t); S.view.x = lerp(from.x, to.x, e); S.view.y = lerp(from.y, to.y, e); S.view.k = lerp(from.k, to.k, e); applyCam(); if (t < 1) S.camT = raf(tick); };
  S.camT = raf(tick);
}
function contentBox() {
  layoutScheme(); const b = S.scheme.box.slice();
  if (!S.mol.atoms.size && !S.chain.length) return null; return b;
}
function fitView(anim = true, maxK = 1.75) {
  const b = contentBox(); if (!b) { const to = { x: 0, y: 0, k: 1.5 }; anim ? animateCam(to) : (Object.assign(S.view, to), applyCam()); return; }
  const pad = 90, W = cw(), H = ch(); const w = b[2] - b[0], h = b[3] - b[1];
  const k = clamp(Math.min((W - pad * 2) / Math.max(w, 1), (H - pad * 2) / Math.max(h, 1)), .35, maxK);
  const to = { k, x: -((b[0] + b[2]) / 2) * k, y: -((b[1] + b[3]) / 2) * k };
  anim ? animateCam(to) : (Object.assign(S.view, to), applyCam());
}
function ensureSchemeVisible() {
  if (!S.scheme || S.g) return; const b = S.scheme.box; const [x0, y0] = toScreen(b[0], b[1]), [x1, y1] = toScreen(b[2], b[3]);
  const m = 36; if (x0 < m || y0 < m || x1 > cw() - m || y1 > ch() - 52) fitView(true, S.view.k);
}
function zoomAt(f, sx, sy) {
  const k0 = S.view.k, k1 = clamp(k0 * f, .3, 4); if (k1 === k0) return;
  const px = sx - cw() / 2 - S.view.x, py = sy - ch() / 2 - S.view.y;
  S.view.x -= px * (k1 / k0 - 1); S.view.y -= py * (k1 / k0 - 1); S.view.k = k1; applyCam();
}

/* ---------------- drawing ---------------- */
function drawAll() { layoutScheme(); drawMol(); drawScheme(); drawHL(); placeOverlays(); updateStatus(); emptyState(); }
function drawMol() {
  const g = S.Lg.mol; if (!g) return;
  if (S.morph) return; // tidy animation owns this layer
  g.innerHTML = S.L0 ? Render.svgMol(S.L0, { cls: 'mol' }) : '';
  const an = analysis(S.mol); let w = '';
  if (an) for (const id of an.bad) { const a = S.mol.atoms.get(id); w += `<circle class="bx-bad" cx="${f1(a.x)}" cy="${f1(a.y)}" r="15"/>`; }
  S.Lg.warn.innerHTML = w;
  if (S.Lg.st) S.Lg.st.innerHTML = S.mol && S.mol.atoms.size ? stMarks(S.mol, 0, 0, true) : '';
}
function arrowSVG(x1, y, x2, cls) { return `<g class="${cls}">${Render.rxnArrow(x1, y, x2, 'rxa')}</g>`; }
function drawScheme() {
  const g = S.Lg.sch; if (!g) return; let o = '';
  const now = performance.now();
  S.chain.forEach((ar, i) => {
    const act = i === S.active; const editing = S.editing && S.editing.kind === 'reagent' && S.editing.id === ar.id;
    let x1 = ar.x1, x2 = ar.x2, y = ar.y;
    if (ar.intro) { const t = clamp((now - ar.intro.t0) / 420, 0, 1), e = easeOut(t); x1 = lerp(ar.intro.x1, x1, e); x2 = lerp(ar.intro.x2, x2, e); y = lerp(ar.intro.y, y, e); if (t >= 1) ar.intro = null; }
    const r = ar.res; const st = r ? r.status : 'empty';
    if (ar.conn && !ar.intro) { const c = ar.conn; o += `<path class="bx-conn${act ? ' on' : ''}" d="M${f1(c.x)},${f1(c.y0)} L${f1(c.x)},${f1(c.y1 - Math.sign(c.y1 - c.y0) * 8)} Q${f1(c.x)},${f1(c.y1)} ${f1(c.x + 8)},${f1(c.y1)} L${f1(c.x1 - 3)},${f1(c.y1)}"/>`; }
    o += `<g class="bx-ar${act ? ' on' : ''} st-${st}" data-ar="${i}"><rect class="hit" x="${f1(x1 - 8)}" y="${f1(y - 34)}" width="${f1(x2 - x1 + 16)}" height="52" rx="8"/>`;
    o += arrowSVG(x1, y, x2, 'line');
    if (!editing) o += ar.reagent ? `<text class="rg" x="${f1((x1 + x2) / 2)}" y="${f1(y - 11)}">${chemSVG(ar.reagent)}</text>` : `<text class="rg ph" x="${f1((x1 + x2) / 2)}" y="${f1(y - 11)}">type a reagent</text>`;
    if (st === 'ok' && r.title) o += `<text class="tt" x="${f1((x1 + x2) / 2)}" y="${f1(y + 21)}">${esc(shortTitle(r.title))}</text>`;
    o += '</g>';
    // product side
    const fl = S.flights.get(ar.id);
    if (ar.parts.length) {
      ar.parts.forEach((p, j) => {
        if (p.plus) { o += `<text class="bx-plus" x="${f1(p.x + 9)}" y="${f1(p.y + 7)}">+</text>`; return; }
        if (fl) { o += flightSVG(p, fl, now, ar.id, p.main); return; }
        o += `<g class="bx-sp${S.prodHover === i && p.main ? ' hv' : ''}${p.main ? ' main' : ' extra'}" data-sp="${i}" transform="translate(${f1(p.dx)},${f1(p.dy)})">${Render.svgMol(p.L, { cls: 'mol' })}</g>`;
        o += stMarks(p.m, p.dx, p.dy);
        if (p.main && ar.badge && !fl) o += `<g class="bx-stb k-${ar.st.kind}"><text x="${f1((p.box[0] + p.box[2]) / 2)}" y="${f1(p.box[3] + 15)}">${esc(ar.badge)}</text></g>`;
      });
    } else if (ar.slot) {
      const [sx0, sy0, sx1, sy1] = ar.slot; const cx = (sx0 + sx1) / 2;
      const msg = { empty: ['?', ''], nomol: ['', 'draw a structure'], nosrc: ['', 'no product before this'], unknown: ['?', 'unknown reagent'], none: ['×', 'no reaction'], invalid: ['!', 'fix the structure'], error: ['!', 'engine error'] }[st] || ['…', ''];
      o += `<g class="bx-slot st-${st}"><rect x="${f1(sx0)}" y="${f1(sy0)}" width="${f1(sx1 - sx0)}" height="${f1(sy1 - sy0)}" rx="12"/>${msg[0] ? `<text class="big" x="${f1(cx)}" y="${f1((sy0 + sy1) / 2 + (msg[1] ? 1 : 8))}">${msg[0]}</text>` : ''}${msg[1] ? `<text class="sm" x="${f1(cx)}" y="${f1((sy0 + sy1) / 2 + (msg[0] ? 22 : 5))}">${msg[1]}</text>` : ''}</g>`;
    }
  });
  if (S.g && S.g.k === 'arrow' && S.g.moved) { const q = S.g; o += `<g class="bx-ar ghost">${Render.rxnArrow(q.x, q.y, Math.max(q.x + 20, q.cx), 'rxa')}</g>`; }
  g.innerHTML = o;
}
/* product atoms glide from where they were in the source; atoms that come from the reagent drift in from its label */
function flightSVG(p, fl, now, id, main) {
  const t = clamp((now - fl.t0) / fl.dur, 0, 1), e = easeIO(t);
  if (t >= 1 && main) S.flights.delete(id);
  const m = p.m.clone(); const fresh = new Set(); const [rx, ry] = fl.rg || [p.box[0], p.box[1]];
  /* fresh atoms start bunched at the reagent label (kept in shape, shrunk), so a group arrives as a group */
  let cx = 0, cy = 0, n = 0; for (const a of m.atoms.values()) if (!fl.from.get(a.id)) { cx += a.x; cy += a.y; n++; } if (n) { cx /= n; cy /= n; }
  const ef = easeOut(clamp((t - .12) / .88, 0, 1));
  for (const a of m.atoms.values()) {
    const tx = a.x + p.dx, ty = a.y + p.dy; const f = fl.from.get(a.id);
    if (f) { a.x = lerp(f[0], tx, e); a.y = lerp(f[1], ty, e); }
    else { const sx = rx + (a.x - cx) * .35, sy = ry + (a.y - cy) * .35; a.x = lerp(sx, tx, ef); a.y = lerp(sy, ty, ef); fresh.add(a.id); }
  }
  const L = lay(m); let svg = Render.svgMol(L, { cls: 'mol' });
  if (fresh.size) svg = svg.replace(/data-a="(\d+)"/g, (s, q) => fresh.has(+q) ? s + ' data-new="1"' : s).replace(/data-b="(\d+)"/g, (s, q) => { const b = m.bonds.get(+q); return b && (fresh.has(b.a) || fresh.has(b.b)) ? s + ' data-new="1"' : s; });
  return `<g class="bx-sp flying${main ? '' : ' extra'}" style="--e:${(main ? e : ef).toFixed(3)}">${svg}</g>`;
}
function loop() { cancelAnimationFrame(S.animT); const tick = () => { drawScheme(); if (S.flights.size || S.chain.some(a => a.intro)) S.animT = raf(tick); else { drawScheme(); placeOverlays(); } }; S.animT = raf(tick); }

/* highlights: hover, selection, ghosts, marquee (under/over the molecule) */
function drawHL() {
  const g = S.Lg.hl, top = S.Lg.ghost; if (!g) return; const m = S.mol; let o = '', t = '';
  const k = S.view.k;
  const bondCap = (b, cls) => { const a = m.atoms.get(b.a), c = m.atoms.get(b.b); return `<line class="${cls}" x1="${f1(a.x)}" y1="${f1(a.y)}" x2="${f1(c.x)}" y2="${f1(c.y)}"/>`; };
  for (const id of S.sel.b) { const b = m.bonds.get(id); if (b) o += bondCap(b, 'bx-selb'); }
  for (const id of S.sel.a) { const a = m.atoms.get(id); if (a) o += `<circle class="bx-sela" cx="${f1(a.x)}" cy="${f1(a.y)}" r="10.5"/>`; }
  const h = S.hover;
  if (h && !S.g) {
    if (h.k === 'atom') { const a = m.atoms.get(h.id); if (a) o += `<circle class="bx-hva" cx="${f1(a.x)}" cy="${f1(a.y)}" r="12"/>`; }
    else if (h.k === 'bond') { const b = m.bonds.get(h.id); if (b) o += bondCap(b, 'bx-hvb'); }
  }
  const G = S.g;
  if (G && G.k === 'move' && G.merge) { const a = m.atoms.get(G.merge[1]); if (a) o += `<circle class="bx-merge" cx="${f1(a.x)}" cy="${f1(a.y)}" r="14"/>`; }
  if (G && G.k === 'marq' && G.moved) { const x = Math.min(G.x0, G.x1), y = Math.min(G.y0, G.y1); t += `<rect class="bx-marq" x="${f1(x)}" y="${f1(y)}" width="${f1(Math.abs(G.x1 - G.x0))}" height="${f1(Math.abs(G.y1 - G.y0))}" style="stroke-width:${(1.2 / k).toFixed(3)}"/>`; }
  const gh = S.ghost;
  if (gh) {
    if (gh.k === 'bond') t += `<line class="bx-gb" x1="${f1(gh.x1)}" y1="${f1(gh.y1)}" x2="${f1(gh.x2)}" y2="${f1(gh.y2)}"/>${gh.dot ? `<circle class="bx-gd" cx="${f1(gh.x2)}" cy="${f1(gh.y2)}" r="3.2"/>` : ''}`;
    if (gh.k === 'ring') t += `<polygon class="bx-gr" points="${gh.pts.map(p => f1(p.x) + ',' + f1(p.y)).join(' ')}"/>`;
    if (gh.k === 'chain') { t += `<polyline class="bx-gc" points="${[[gh.x, gh.y], ...gh.pts.map(p => [p.x, p.y])].map(p => f1(p[0]) + ',' + f1(p[1])).join(' ')}"/>`; const l = gh.pts[gh.pts.length - 1]; t += `<g class="bx-cnt" transform="translate(${f1(l.x + 12 / k)},${f1(l.y - 12 / k)}) scale(${(1 / k).toFixed(3)})"><rect x="-2" y="-14" width="${String(gh.pts.length).length * 8 + 14}" height="20" rx="10"/><text x="5" y="0">${gh.pts.length}</text></g>`; }
  }
  g.innerHTML = o; top.innerHTML = t;
}

/* ---------------- picking ---------------- */
function pickAtom(wx, wy, px = 13, skip) { return nearAtom(wx, wy, Math.max(9, px / S.view.k), skip); }
function pickBond(wx, wy) {
  const m = S.mol; let best = null, bd = Math.max(5.5, 8 / S.view.k);
  for (const b of m.bonds.values()) {
    const a = m.atoms.get(b.a), c = m.atoms.get(b.b); const vx = c.x - a.x, vy = c.y - a.y, l = vx * vx + vy * vy || 1;
    const t = clamp(((wx - a.x) * vx + (wy - a.y) * vy) / l, .1, .9); const q = Math.hypot(wx - (a.x + t * vx), wy - (a.y + t * vy));
    if (q < bd) { bd = q; best = b; }
  }
  return best;
}
function pick(wx, wy) {
  const a = pickAtom(wx, wy); if (a) return { k: 'atom', id: a.id };
  const b = pickBond(wx, wy); if (b) return { k: 'bond', id: b.id };
  for (let i = 0; i < S.chain.length; i++) {
    const ar = S.chain[i];
    if (wx >= ar.x1 - 8 && wx <= ar.x2 + 8 && wy >= ar.y - 34 && wy <= ar.y + 18) return { k: 'arrow', i };
    const p = ar.parts && ar.parts.find(q => q.main); if (p && wx >= p.box[0] - 10 && wx <= p.box[2] + 10 && wy >= p.box[1] - 10 && wy <= p.box[3] + 10) return { k: 'prod', i };
    if (ar.slot && wx >= ar.slot[0] && wx <= ar.slot[2] && wy >= ar.slot[1] && wy <= ar.slot[3]) return { k: 'arrow', i };
  }
  return null;
}
function moveSet() { const s = new Set(S.sel.a); for (const id of S.sel.b) { const b = S.mol.bonds.get(id); if (b) { s.add(b.a); s.add(b.b); } } return [...s].filter(id => S.mol.atoms.has(id)); }
function pruneSel() { for (const id of [...S.sel.a]) if (!S.mol.atoms.has(id)) S.sel.a.delete(id); for (const id of [...S.sel.b]) if (!S.mol.bonds.has(id)) S.sel.b.delete(id); }
function dropLoneCarbons(ids) { for (const id of ids) { const a = S.mol.atoms.get(id); if (a && a.el === 'C' && !a.chg && !S.mol.degree(id)) S.mol.removeAtom(id); } }
function deleteSelection() {
  const atoms = [...S.sel.a], bonds = [...S.sel.b]; if (!atoms.length && !bonds.length) return false;
  pushUndo(); const touched = new Set();
  for (const id of bonds) { const b = S.mol.bonds.get(id); if (b) { touched.add(b.a); touched.add(b.b); S.mol.removeBond(id); } }
  for (const id of atoms) { for (const n of S.mol.nbrs(id)) touched.add(n); S.mol.removeAtom(id); }
  dropLoneCarbons(touched); S.sel.a.clear(); S.sel.b.clear(); S.hover = null; changed(); return true;
}
function deleteHovered() {
  const h = S.hover; if (!h) return false;
  if (h.k === 'atom' || h.k === 'bond') { pushUndo(); if (h.k === 'atom') { const nb = S.mol.nbrs(h.id); S.mol.removeAtom(h.id); dropLoneCarbons(nb); } else { const b = S.mol.bonds.get(h.id); S.mol.removeBond(h.id); if (b) dropLoneCarbons([b.a, b.b]); } S.hover = null; pruneSel(); changed(); return true; }
  if (h.k === 'arrow') { removeArrow(h.i); return true; }
  return false;
}
function eraseAt(wx, wy) {
  const G = S.g; const h = pick(wx, wy); if (!h) return;
  if (h.k === 'arrow' || h.k === 'prod') { if (!G.did && !G.moved) removeArrow(h.i); return; }
  if (!G.did) { pushUndo(); G.did = true; }
  if (h.k === 'atom') { const nb = S.mol.nbrs(h.id); S.mol.removeAtom(h.id); dropLoneCarbons(nb); }
  else { const b = S.mol.bonds.get(h.id); S.mol.removeBond(h.id); if (b) dropLoneCarbons([b.a, b.b]); }
  pruneSel(); S.hover = null; layoutScheme(); drawMol(); drawScheme(); drawHL(); updateStatus();
}
function selectAll() { S.sel.a = new Set(S.mol.atoms.keys()); S.sel.b = new Set(S.mol.bonds.keys()); drawHL(); updateStatus(); }

/* ---------------- arrows ---------------- */
/* which compound an arrow starts from: a product box it starts on, else the nearest drawn compound;
   a selection that spans several compounds makes them react together */
function selFrags() { if (!S.sel.a.size) return null; const fr = S.mol.frags().filter(c => c.some(id => S.sel.a.has(id))); return fr.length ? fr : null; }
function srcNear(wx, wy, reach = 120) {
  let best = null, bd = reach;
  for (const ar of S.chain) { const q = ar.parts && ar.parts.find(x => x.main); if (!q) continue; const b = q.box; const d = Math.hypot(Math.max(b[0] - wx, 0, wx - b[2]), Math.max(b[1] - wy, 0, wy - b[3])); if (d < bd * .7) { bd = d / .7; best = { a: ar.id }; } }
  for (const f of S.fragBox || []) { const b = f.box; const d = Math.hypot(Math.max(b[0] - wx, 0, wx - b[2]), Math.max(b[1] - wy, 0, wy - b[3])); if (d < bd) { bd = d; best = { f: [...f.ids] }; } }
  const sf = selFrags(); if (sf && (!best || (best.f && sf.some(c => c.includes(best.f[0]))))) return { f: sf.flat() };
  return best;
}
function autoSrc() {
  const sf = selFrags(); if (sf) return { f: sf.flat() };
  const act = S.chain[S.active]; if (act && act.res && act.res.status === 'ok') return { a: act.id };
  const fr = S.mol.frags(); const used = new Set(S.chain.filter(a => a.src && a.src.f).flatMap(a => a.src.f));
  const free = fr.find(c => !c.some(id => used.has(id))); if (free) return { f: free };
  const last = S.chain[S.chain.length - 1]; if (last && last.res && last.res.status === 'ok') return { a: last.id };
  return fr.length ? { f: fr[0] } : null;
}
function addArrow(reagent, from, src) {
  if (!S.mol.atoms.size && !(src && src.a != null)) { toast('Draw a structure first, then add the arrow'); return; }
  src = src || autoSrc(); if (!src) { toast('Draw a structure first, then add the arrow'); return; }
  if (src.a != null) { const p = arrowById(src.a); const pi = S.chain.indexOf(p); if (!p || !(p.res && p.res.status === 'ok')) { if (pi >= 0) { setActive(pi); openReagentEditor(pi); } toast('Finish this arrow first: it has no product yet'); return; } }
  pushUndo(); const ar = arrowObj(reagent || '', null, src); S.chain.push(ar); S.active = S.chain.length - 1;
  layoutScheme();
  if (from && motion()) { ar.intro = { x1: from.x1, x2: from.x2, y: from.y, t0: performance.now() }; loop(); }
  drawAll(); renderPanel(); ensureSchemeVisible();
  if (!reagent) setTimeout(() => openReagentEditor(S.chain.indexOf(ar)), from && motion() ? 260 : 0); else scheduleSolve(30);
}
function removeArrow(i) {
  const ar = S.chain[i]; if (!ar) return; pushUndo(); const gone = descendants(ar.id); const n = gone.size;
  S.chain = S.chain.filter(a => !gone.has(a.id)); S.active = Math.min(S.active, S.chain.length - 1); closeEditors(true);
  toast(n > 1 ? `Removed arrow ${i + 1} and the ${n - 1} that follow from it` : `Removed arrow ${i + 1}`); changed({ noFlight: true });
}
function setActive(i) { if (S.active === i) return; S.active = i; drawScheme(); renderPanel(); }

/* ---------------- pointer ---------------- */
function onDown(e) {
  if (e.target.closest('.bx-edit,.bx-pill,.bx-empty button,.bx-help')) return;
  if (S.editing) closeEditors();
  S.can.focus({ preventScroll: true });
  const [wx, wy] = toWorld(e.clientX, e.clientY);
  const base = { sx: e.clientX, sy: e.clientY, moved: false };
  if (e.button === 1 || e.button === 2 || S.space) { S.g = { ...base, k: 'pan', vx: S.view.x, vy: S.view.y }; S.can.setPointerCapture(e.pointerId); S.can.classList.add('panning'); e.preventDefault(); return; }
  if (e.button !== 0) return;
  cancelAnimationFrame(S.camT);
  let h = pick(wx, wy); const T = S.tool;
  S.can.setPointerCapture(e.pointerId);
  /* our own double-click: the first click may already have redrawn the SVG, so the browser's
     dblclick never fires. Undo whatever the first click drew, then open the label editor. */
  const now = performance.now(), pc = S.lastDown;
  S.lastDown = { t: now, x: e.clientX, y: e.clientY, undoLen: S.undo.length };
  if (pc && now - pc.t < 430 && Math.hypot(e.clientX - pc.x, e.clientY - pc.y) < 8 && T !== 'erase' && T !== 'arrow' && !(h && (h.k === 'arrow' || h.k === 'prod'))) {
    S.lastDown = null;
    if (S.undo.length > pc.undoLen) { const snap = S.undo[pc.undoLen]; S.undo.length = pc.undoLen; S.redo.length = 0; restoreSnap(snap); h = pick(wx, wy); }
    const target = h && h.k === 'atom' ? h.id : null;
    if (!h || h.k === 'atom') S.g = { ...base, k: 'click', fn: () => openLabelEditor(target, target == null ? { x: wx, y: wy } : null) };
    return;
  }
  if (h && h.k === 'arrow') { if (T === 'erase') { removeArrow(h.i); return; } S.g = { ...base, k: 'click', fn: () => { setActive(h.i); openReagentEditor(h.i); } }; return; }
  if (h && h.k === 'prod' && T !== 'arrow') { S.g = { ...base, k: 'click', fn: () => setActive(h.i) }; return; }
  if (T === 'select') {
    if (h) {
      const set = h.k === 'atom' ? S.sel.a : S.sel.b; const was = set.has(h.id);
      if (!was && !e.shiftKey) { S.sel.a.clear(); S.sel.b.clear(); }
      if (!was) set.add(h.id);
      const ids = moveSet(); const orig = new Map(ids.map(id => { const a = S.mol.atoms.get(id); return [id, [a.x, a.y]]; }));
      S.g = { ...base, k: 'move', ids, orig, wx0: wx, wy0: wy, hit: h, was, shift: e.shiftKey }; drawHL(); updateStatus(); return;
    }
    S.g = { ...base, k: 'marq', x0: wx, y0: wy, x1: wx, y1: wy, shift: e.shiftKey }; return;
  }
  if (T === 'bond') {
    if (h && h.k === 'bond') { S.g = { ...base, k: 'click', fn: () => { pushUndo(); modifyBond(S.mol.bonds.get(h.id), S.bondType); changed(); } }; return; }
    const a = h ? S.mol.atoms.get(h.id) : null; S.g = { ...base, k: 'bond', from: a ? a.id : null, x: a ? a.x : wx, y: a ? a.y : wy }; return;
  }
  if (T === 'chain') { const a = h && h.k === 'atom' ? S.mol.atoms.get(h.id) : null; S.g = { ...base, k: 'chain', from: a ? a.id : null, x: a ? a.x : wx, y: a ? a.y : wy }; return; }
  if (T === 'ring') {
    S.g = { ...base, k: 'click', fn: () => { const n = S.ringType === 'benz' ? 6 : +S.ringType; pushUndo(); placeRing(n, S.ringType === 'benz', ringAnchor(h, wx, wy)); S.ghost = null; changed(); } }; return;
  }
  if (T === 'atom') {
    S.g = { ...base, k: 'click', fn: () => { if (h && h.k === 'atom') { if (applyLabel(h.id, S.atomLabel)) changed(); } else if (!h) { if (applyLabel(null, S.atomLabel, { x: wx, y: wy })) changed(); } } }; return;
  }
  if (T === 'charge') {
    const neg = e.shiftKey;
    S.g = { ...base, k: 'click', fn: () => { if (!h || h.k !== 'atom') return; pushUndo(); const a = S.mol.atoms.get(h.id); if (S.chargeMode === 'rad') a.rad = a.rad ? 0 : 1; else { const sg = (S.chargeMode === '-') !== neg ? -1 : 1; a.chg = clamp((a.chg || 0) + sg, -2, 2); } a.hFix = null; changed(); } }; return;
  }
  if (T === 'erase') { S.g = { ...base, k: 'erase', did: false }; eraseAt(wx, wy); return; }
  if (T === 'arrow') { const b = S.scheme ? S.scheme.box : [0, 0, 0, 0]; S.g = { ...base, k: 'arrow', x: wx, y: wy, cx: wx, cy: wy }; return; }
}
function ringAnchor(h, wx, wy) { return h && h.k === 'atom' ? { k: 'atom', id: h.id } : h && h.k === 'bond' ? { k: 'bond', id: h.id } : { k: 'free', x: wx, y: wy }; }
function bondEnd(G, wx, wy, e) {
  const t = pickAtom(wx, wy, 15, G.from != null ? new Set([G.from]) : null);
  if (t && (G.from == null || t.id !== G.from)) return { x: t.x, y: t.y, id: t.id };
  let th = Math.atan2(wy - G.y, wx - G.x); const step = e.shiftKey ? D30 / 2 : D30; if (!e.altKey) th = Math.round(th / step) * step;
  return { x: G.x + BL * Math.cos(th), y: G.y + BL * Math.sin(th), id: null };
}
let hovRaf = 0, lastPt = null;
function onMove(e) {
  const G = S.g; lastPt = [e.clientX, e.clientY];
  if (!G) { if (!hovRaf) hovRaf = raf(() => { hovRaf = 0; if (lastPt && !S.g) { const [x, y] = toWorld(lastPt[0], lastPt[1]); hoverAt(x, y); } }); return; }
  const [wx, wy] = toWorld(e.clientX, e.clientY);
  if (!G.moved && Math.hypot(e.clientX - G.sx, e.clientY - G.sy) > 4) G.moved = true;
  if (!G.moved) return;
  if (G.k === 'pan') { S.view.x = G.vx + e.clientX - G.sx; S.view.y = G.vy + e.clientY - G.sy; applyCam(); return; }
  if (G.k === 'move') {
    if (!G.undone) { pushUndo(); G.undone = true; }
    const dx = wx - G.wx0, dy = wy - G.wy0;
    for (const [id, [x, y]] of G.orig) { const a = S.mol.atoms.get(id); if (a) { a.x = x + dx; a.y = y + dy; } }
    G.merge = null; const still = new Set(G.ids);
    if (G.ids.length < S.mol.atoms.size) for (const id of G.ids) { const a = S.mol.atoms.get(id); const t = nearAtom(a.x, a.y, BL * .3, still); if (t) { G.merge = [id, t.id]; break; } }
    layoutScheme(); drawMol(); drawScheme(); drawHL(); placeOverlays(); return;
  }
  if (G.k === 'marq') { G.x1 = wx; G.y1 = wy; drawHL(); return; }
  if (G.k === 'bond') { const end = bondEnd(G, wx, wy, e); S.ghost = { k: 'bond', x1: G.x, y1: G.y, x2: end.x, y2: end.y, dot: end.id == null }; drawHL(); return; }
  if (G.k === 'chain') { S.ghost = { k: 'chain', x: G.x, y: G.y, pts: chainPts(G.x, G.y, wx, wy, G.from) }; drawHL(); return; }
  if (G.k === 'erase') { eraseAt(wx, wy); return; }
  if (G.k === 'arrow') { G.cx = wx; G.cy = wy; drawScheme(); return; }
}
function onUp(e) {
  const G = S.g; if (!G) return; S.g = null; S.can.classList.remove('panning');
  try { S.can.releasePointerCapture(e.pointerId); } catch (_) { }
  const [wx, wy] = toWorld(e.clientX, e.clientY); const m = S.mol;
  if (G.k === 'click') { if (!G.moved) G.fn(); }
  else if (G.k === 'move') {
    if (G.moved) { if (G.merge) { mergeAtoms(G.merge[1], G.merge[0]); S.sel.a.delete(G.merge[0]); pruneSel(); changed(); } else { drawAll(); } }
    else { const set = G.hit.k === 'atom' ? S.sel.a : S.sel.b; if (G.shift) { if (G.was) set.delete(G.hit.id); } else { S.sel.a.clear(); S.sel.b.clear(); set.add(G.hit.id); } drawHL(); updateStatus(); }
  } else if (G.k === 'marq') {
    if (!G.shift) { S.sel.a.clear(); S.sel.b.clear(); }
    if (G.moved) {
      const x0 = Math.min(G.x0, G.x1), x1 = Math.max(G.x0, G.x1), y0 = Math.min(G.y0, G.y1), y1 = Math.max(G.y0, G.y1);
      for (const a of m.atoms.values()) if (a.x >= x0 && a.x <= x1 && a.y >= y0 && a.y <= y1) S.sel.a.add(a.id);
      for (const b of m.bonds.values()) if (S.sel.a.has(b.a) && S.sel.a.has(b.b)) S.sel.b.add(b.id);
    }
    drawHL(); updateStatus();
  } else if (G.k === 'bond') {
    S.ghost = null; pushUndo();
    if (!G.moved) {
      if (G.from != null) { const a = m.atoms.get(G.from); const th = sproutAngle(m, G.from, S.bondType); const px = a.x + BL * Math.cos(th), py = a.y + BL * Math.sin(th); const ex = nearAtom(px, py, BL * .3, new Set([G.from])); makeBond(G.from, ex ? ex.id : m.addAtom('C', px, py), S.bondType); }
      else { const a = m.addAtom('C', G.x, G.y), b = m.addAtom('C', G.x + BL * Math.cos(-D30), G.y + BL * Math.sin(-D30)); makeBond(a, b, S.bondType); }
    } else {
      const end = bondEnd(G, wx, wy, e); const from = G.from != null ? G.from : m.addAtom('C', G.x, G.y);
      const to = end.id != null ? end.id : ((nearAtom(end.x, end.y, BL * .3, new Set([from])) || {}).id ?? m.addAtom('C', end.x, end.y));
      if (from !== to) { const ex = m.bondBetween(from, to); if (ex) modifyBond(ex, S.bondType); else makeBond(from, to, S.bondType); }
    }
    changed();
  } else if (G.k === 'chain') {
    S.ghost = null; pushUndo();
    let prev = G.from != null ? G.from : m.addAtom('C', G.x, G.y);
    let pts;
    if (G.moved) pts = chainPts(G.x, G.y, wx, wy, G.from);
    else { const a = m.atoms.get(prev); const th = m.degree(prev) ? sproutAngle(m, prev, 'single') : -D30; pts = [{ x: a.x + BL * Math.cos(th), y: a.y + BL * Math.sin(th) }]; }
    for (const p of pts) { const ex = nearAtom(p.x, p.y, BL * .3, new Set([prev])); const id = ex ? ex.id : m.addAtom('C', p.x, p.y); if (id !== prev && !m.bondBetween(prev, id)) m.addBond(prev, id, 1); prev = id; }
    changed();
  } else if (G.k === 'erase') { if (G.did) changed(); }
  else if (G.k === 'arrow') { drawScheme(); const src = srcNear(G.x, G.y) || srcNear(wx, wy); addArrow('', G.moved ? { x1: G.x, x2: Math.max(G.x + 30, G.cx), y: G.y } : null, src || undefined); }
  const [hx, hy] = toWorld(e.clientX, e.clientY); hoverAt(hx, hy);
}
function hoverAt(wx, wy) {
  if (!S.can) return;
  const h = pick(wx, wy); const prev = S.hover; S.hover = h;
  let ph = h && h.k === 'prod' ? h.i : -1;
  if (ph < 0 && S.prodHover >= 0 && !h) { const ar = S.chain[S.prodHover]; const q = ar && ar.parts && ar.parts.find(x => x.main); if (q && wx >= q.box[0] - 30 && wx <= q.box[2] + 30 && wy >= q.box[1] - 70 / Math.max(.5, S.view.k) && wy <= q.box[3] + 10) ph = S.prodHover; }
  if (ph !== S.prodHover) { if (ph < 0 && S.prodHover >= 0) { S.pillHold = S.prodHover; clearTimeout(S.pillT); S.pillT = setTimeout(() => { if (!S.overPill) { S.pillHold = null; placeOverlays(); } }, 650); } S.prodHover = ph; if (ph >= 0) S.pillHold = null; drawScheme(); placeOverlays(); }
  S.ghost = null; const T = S.tool;
  if (T === 'ring' && (!h || h.k === 'atom' || h.k === 'bond')) { const n = S.ringType === 'benz' ? 6 : +S.ringType; try { S.ghost = { k: 'ring', pts: ringGeom(n, ringAnchor(h, wx, wy)) }; } catch (e) { S.ghost = null; } }
  if ((T === 'bond' || T === 'chain') && h && h.k === 'atom') { const a = S.mol.atoms.get(h.id); const th = sproutAngle(S.mol, h.id, T === 'bond' ? S.bondType : 'single'); S.ghost = { k: 'bond', x1: a.x, y1: a.y, x2: a.x + BL * Math.cos(th), y2: a.y + BL * Math.sin(th), dot: true }; }
  const cur = !h ? (T === 'select' ? 'default' : T === 'erase' ? 'cell' : 'crosshair') : (h.k === 'arrow' || h.k === 'prod') ? 'pointer' : T === 'select' ? 'move' : T === 'erase' ? 'cell' : 'crosshair';
  S.can.style.cursor = S.space ? 'grab' : cur;
  drawHL(); if (!prev !== !h || (prev && h && (prev.k !== h.k || prev.id !== h.id))) updateStatus();
}
function onWheel(e) {
  e.preventDefault(); const r = S.can.getBoundingClientRect();
  if (!e.ctrlKey && Math.abs(e.deltaX) > Math.abs(e.deltaY) * 1.2) { S.view.x -= e.deltaX; applyCam(); return; }
  if (e.shiftKey) { S.view.x -= e.deltaY; applyCam(); return; }
  const d = e.deltaMode === 1 ? e.deltaY * 16 : e.deltaY;
  zoomAt(Math.exp(-d * (e.ctrlKey ? .01 : .0016)), e.clientX - r.left, e.clientY - r.top);
}
function onDbl(e) {
  const [wx, wy] = toWorld(e.clientX, e.clientY); const h = pick(wx, wy);
  if (h && h.k === 'atom') openLabelEditor(h.id); else if (!h && S.tool !== 'erase') openLabelEditor(null, { x: wx, y: wy });
}

/* ---------------- keyboard ---------------- */
const boardOn = () => S.el && S.el.classList.contains('on') && $('#viewer').hidden && $('#cmdp').hidden;
const TOOLKEYS = { v: 'select', b: 'bond', z: 'chain', r: 'ring', t: 'atom', q: 'charge', e: 'erase', a: 'arrow' };
const ELKEYS = { c: 'C', n: 'N', o: 'O', s: 'S', f: 'F', p: 'P', h: 'H', l: 'Cl', b: 'Br', i: 'I' };
function onKey(e) {
  if (!boardOn()) return;
  const t = e.target; if (t && t.matches && t.matches('input,textarea,select,[contenteditable]')) return;
  const k = e.key, lk = k.length === 1 ? k.toLowerCase() : k, mod = e.ctrlKey || e.metaKey;
  if (mod) {
    if (lk === 'z') { e.preventDefault(); e.stopImmediatePropagation(); e.shiftKey ? redo() : undo(); return; }
    if (lk === 'y') { e.preventDefault(); redo(); return; }
    if (lk === 'a') { e.preventDefault(); selectAll(); return; }
    if (lk === 'c') { e.preventDefault(); copySmiles(); return; }
    if (lk === 'd') { e.preventDefault(); duplicateSel(); return; }
    if (lk === 'k' && e.shiftKey) { e.preventDefault(); e.stopImmediatePropagation(); tidy(); return; }
    if (k === '0') { e.preventDefault(); fitView(); return; }
    if (k === '=' || k === '+') { e.preventDefault(); zoomAt(1.2, cw() / 2, ch() / 2); return; }
    if (k === '-') { e.preventDefault(); zoomAt(1 / 1.2, cw() / 2, ch() / 2); return; }
    return;
  }
  if (e.altKey) return;
  if (k === ' ') { if (!S.space) { S.space = true; S.can.style.cursor = 'grab'; } e.preventDefault(); return; }
  if (k === 'Escape') { if (S.g) { S.g = null; S.ghost = null; } if (helpOpen()) { toggleHelp(false); return; } if (S.sel.a.size || S.sel.b.size) { S.sel.a.clear(); S.sel.b.clear(); } else if (S.tool !== 'select') setTool('select'); drawHL(); updateStatus(); return; }
  if (k === 'Delete' || k === 'Backspace') { e.preventDefault(); if (!deleteSelection()) deleteHovered(); return; }
  if (k === '?') { toggleHelp(); return; }
  if (lk === 'm' && !(S.hover && S.hover.k === 'atom')) { e.preventDefault(); togglePanel(); return; }
  if (lk === 'n' && !(S.hover && S.hover.k === 'atom')) { const ar = S.chain[S.active]; if (ar && ar.res && ar.res.status === 'ok') { e.preventDefault(); addArrow('', null, { a: ar.id }); } return; }
  if (k === '*') { e.preventDefault(); setStOpt('star'); toast(S.stOpt.star ? 'Stereocentres starred' : 'Stereocentre stars off'); return; }
  if (k.startsWith('Arrow')) {
    const ids = moveSet(); const d = e.shiftKey ? 10 : 2; const v = { ArrowLeft: [-d, 0], ArrowRight: [d, 0], ArrowUp: [0, -d], ArrowDown: [0, d] }[k];
    e.preventDefault();
    if (ids.length) { if (!e.repeat) pushUndo(); for (const id of ids) { const a = S.mol.atoms.get(id); a.x += v[0]; a.y += v[1]; } drawAll(); return; }
    if (k === 'ArrowLeft' || k === 'ArrowRight') mvStep(k === 'ArrowRight' ? 1 : -1);
    return;
  }
  const h = S.hover;
  if (h && h.k === 'atom') {
    const a = S.mol.atoms.get(h.id);
    if (ELKEYS[lk]) { e.preventDefault(); pushUndo(); a.el = ELKEYS[lk]; a.hFix = null; changed(); return; }
    if (k === '+' || k === '=' || k === '-') { e.preventDefault(); pushUndo(); a.chg = clamp((a.chg || 0) + (k === '-' ? -1 : 1), -2, 2); a.hFix = null; changed(); return; }
    if (k === 'Enter') { e.preventDefault(); openLabelEditor(h.id); return; }
    if (/^[1-3]$/.test(k)) { e.preventDefault(); pushUndo(); const ty = ['single', 'double', 'triple'][+k - 1]; const th = sproutAngle(S.mol, h.id, ty); const nid = S.mol.addAtom('C', a.x + BL * Math.cos(th), a.y + BL * Math.sin(th)); makeBond(h.id, nid, ty); changed(); return; }
  }
  if (h && h.k === 'bond') {
    const b = S.mol.bonds.get(h.id);
    if (/^[1-3]$/.test(k)) { e.preventDefault(); pushUndo(); b.o = +k; b.st = 0; changed(); return; }
    if (lk === 'w' || lk === 'h') { e.preventDefault(); pushUndo(); modifyBond(b, lk === 'w' ? 'wedge' : 'hash'); changed(); return; }
  }
  if (S.tool === 'ring' && /^[1-8]$/.test(k) && k !== '2') { S.ringType = k === '1' ? 'benz' : k; renderCtx(); hoverAgain(); return; }
  if (S.tool === 'bond' && /^[1-3]$/.test(k)) { S.bondType = ['single', 'double', 'triple'][+k - 1]; renderCtx(); return; }
  if (S.tool === 'bond' && (lk === 'w' || lk === 'h')) { S.bondType = lk === 'w' ? 'wedge' : 'hash'; renderCtx(); return; }
  if (lk === 'f') { fitView(); return; }
  if (TOOLKEYS[lk]) { setTool(TOOLKEYS[lk]); return; }
}
function onKeyUp(e) { if (e.key === ' ') { S.space = false; if (S.can) S.can.style.cursor = ''; } }
function hoverAgain() { if (lastPt) { const [x, y] = toWorld(lastPt[0], lastPt[1]); hoverAt(x, y); } }

/* ---------------- inline editors (HTML over the canvas) ---------------- */
const POPULAR = ['HBr', 'Br2 in CCl4', 'H2, Pd', 'alc. KOH', 'aq. KOH', 'conc. H2SO4, heat', 'NaBH4', 'LiAlH4', 'PCC', 'KMnO4, H+, heat', 'HNO3, H2SO4', 'Br2, FeBr3', 'CH3MgBr then H3O+', 'O3 then Zn/H2O', 'BH3.THF then H2O2/NaOH', 'HBr, ROOR', 'NaNH2', 'mCPBA', 'dil. NaOH', 'I2, NaOH'];
const EXTRA_RG = ['H3O+', 'Zn-Hg, HCl', 'NH2NH2, KOH, glycol, heat', 'NaNO2, HCl, 0-5 C', 'SOCl2', 'PCl5', 'conc. NaOH', 'HgSO4, dil. H2SO4', 'Br2, hv', 'Cl2, hv', 'CHCl3, KOH', 'CHCl3, NaOH', 'CuCN', 'Sn/HCl', 't-BuOK', 'NaCN', 'AgCN', 'HCN', 'NH2OH', 'H2, Lindlar', 'Br2/H2O', 'cold dil. KMnO4', 'K2Cr2O7, H2SO4', 'NaHCO3', 'Na', "Tollens' reagent", "Fehling's solution", 'soda lime, heat', 'Lucas reagent', 'HI', 'CH3COCl, AlCl3', 'NH3'];
let RGC = null;
function rgCands() {
  if (RGC) return RGC; const seen = new Set(); RGC = [];
  const add = (v, src, n) => { const k = normRg(v); if (!k || seen.has(k)) return; seen.add(k); RGC.push({ v, k, src, n: n || 0 }); };
  POPULAR.forEach(v => add(v, 'pop')); (Core && Core.REC || window.EngineCore.REC || []).forEach(r => r.canon && add(r.canon, 'eng')); EXTRA_RG.forEach(v => add(v, 'pop'));
  (A.RGL || []).forEach(e => add(e.label, 'lib', e.ids.size));
  return RGC;
}
const SEP = /(\s*(?:,|;|\bthen\b|→|->)\s*)/i;
function lastSeg(s) { const parts = s.split(SEP); return { head: parts.slice(0, -1).join(''), seg: parts[parts.length - 1] || '' }; }
function rgMatches(text) {
  const all = rgCands(); const q = normRg(text); const { head, seg } = lastSeg(text); const qs = normRg(seg);
  if (!q) return POPULAR.slice(0, 8).map(v => ({ v, label: v, src: 'pop' }));
  const out = [], used = new Set();
  const push = (c, full) => { const v = full ? c.v : head + c.v; const k = normRg(v); if (used.has(k) || k === q) return; used.add(k); out.push({ v, label: c.v, src: c.src, n: c.n, full }); };
  for (const c of all) if (c.k.startsWith(q)) push(c, true);
  if (head && qs) for (const c of all) if (c.k.startsWith(qs)) push(c, false);
  if (out.length < 7 && q.length > 1) for (const c of all) if (c.k.includes(q)) push(c, true);
  return out.slice(0, 7);
}
function closeEditors(silent) {
  const ed = S.editing; if (!ed) return; S.editing = null;
  if (ed.kind === 'reagent') { const ar = S.chain.find(a => a.id === ed.id); if (ar && ar.reagent !== ed.orig) { S.undo.push(ed.snap); S.redo.length = 0; } }
  ed.el.remove(); if (!silent) { drawScheme(); placeOverlays(); }
}
function openReagentEditor(i) {
  const ar = S.chain[i]; if (!ar) return; closeEditors(true); if (S.active !== i) { S.active = i; renderPanel(); }
  const box = document.createElement('div'); box.className = 'bx-edit bx-rge';
  box.innerHTML = `<div class="f"><input spellcheck="false" autocomplete="off" placeholder="reagent, conditions…" aria-label="Reagent for arrow ${i + 1}"><kbd>↵</kbd></div><div class="bx-sug" role="listbox"></div>`;
  S.can.append(box); const ed = S.editing = { kind: 'reagent', id: ar.id, el: box, sel: -1, orig: ar.reagent, snap: snapshot() };
  const inp = box.querySelector('input'), sug = box.querySelector('.bx-sug'); inp.value = ar.reagent;
  const list = () => {
    const ms = rgMatches(inp.value); ed.items = ms; ed.sel = Math.min(ed.sel, ms.length - 1);
    sug.innerHTML = ms.length ? ms.map((x, j) => `<button type="button" class="it${j === ed.sel ? ' on' : ''}" data-j="${j}" role="option"><span>${chemHTML(x.label)}</span>${x.src === 'lib' ? `<em>${x.n} in library</em>` : x.src === 'eng' ? '<em>engine</em>' : ''}</button>`).join('') : '';
    sug.hidden = !ms.length;
  };
  const upd = () => { ar.reagent = inp.value; ar.forceLib = false; layoutScheme(); drawScheme(); placeOverlays(); scheduleSolve(140); list(); updateRgState(); };
  const accept = j => { const x = ed.items && ed.items[j]; if (!x) return false; inp.value = x.v; ed.sel = -1; upd(); return true; };
  inp.addEventListener('input', () => { ed.sel = -1; upd(); });
  inp.addEventListener('keydown', e => {
    if (e.key === 'ArrowDown' || e.key === 'ArrowUp') { e.preventDefault(); const n = (ed.items || []).length; if (!n) return; let j = ed.sel + (e.key === 'ArrowDown' ? 1 : -1); if (j >= n) j = -1; if (j < -1) j = n - 1; ed.sel = j; list(); return; }
    if (e.key === 'Tab') { e.preventDefault(); accept(ed.sel >= 0 ? ed.sel : 0); return; }
    if (e.key === 'Enter') { e.preventDefault(); if (ed.sel >= 0 && accept(ed.sel)) return; scheduleSolve(0); closeEditors(); S.can.focus({ preventScroll: true }); return; }
    if (e.key === 'Escape') { e.preventDefault(); closeEditors(); S.can.focus({ preventScroll: true }); return; }
    e.stopPropagation();
  });
  sug.addEventListener('pointerdown', e => { const b = e.target.closest('[data-j]'); if (!b) return; e.preventDefault(); accept(+b.dataset.j); inp.focus(); });
  inp.addEventListener('blur', () => setTimeout(() => { if (S.editing === ed && !box.contains(document.activeElement)) closeEditors(); }, 140));
  list(); drawScheme(); placeOverlays(); updateRgState(); inp.focus(); inp.setSelectionRange(inp.value.length, inp.value.length);
}
function updateRgState() { const ed = S.editing; if (!ed || ed.kind !== 'reagent') return; const ar = S.chain.find(a => a.id === ed.id); const st = ar && ar.res ? ar.res.status : 'empty'; ed.el.dataset.st = st; }
function openLabelEditor(id, at) {
  closeEditors(true); const a = id != null ? S.mol.atoms.get(id) : null; if (id != null && !a) return;
  const box = document.createElement('div'); box.className = 'bx-edit bx-lbe';
  box.innerHTML = `<input spellcheck="false" autocomplete="off" placeholder="OH, NO₂…" aria-label="Atom label"><div class="bx-chips"></div>`;
  S.can.append(box); const ed = S.editing = { kind: 'label', id, at: at || null, el: box };
  const inp = box.querySelector('input'), chips = box.querySelector('.bx-chips');
  inp.value = a && a.el !== 'C' ? a.el : '';
  const list = () => { const q = inp.value.trim().toLowerCase(); const all = [...ELEMS, ...ABBR_SHOW]; const ms = (q ? all.filter(x => x.toLowerCase().startsWith(q)) : ABBR_SHOW.slice(0, 10)).slice(0, 10); ed.items = ms; chips.innerHTML = ms.map((x, j) => `<button type="button" data-j="${j}"><span>${chemHTML(x)}</span></button>`).join(''); chips.hidden = !ms.length; };
  const apply = v => { const ok = applyLabel(id, v, at); closeEditors(true); if (ok) changed(); else drawAll(); S.can.focus({ preventScroll: true }); };
  inp.addEventListener('input', list);
  inp.addEventListener('keydown', e => { e.stopPropagation(); if (e.key === 'Enter') { e.preventDefault(); apply(inp.value); } else if (e.key === 'Escape') { e.preventDefault(); closeEditors(); S.can.focus({ preventScroll: true }); } else if (e.key === 'Tab') { e.preventDefault(); if (ed.items && ed.items[0]) { inp.value = ed.items[0]; list(); } } });
  chips.addEventListener('pointerdown', e => { const b = e.target.closest('[data-j]'); if (!b) return; e.preventDefault(); apply(ed.items[+b.dataset.j]); });
  inp.addEventListener('blur', () => setTimeout(() => { if (S.editing === ed && !box.contains(document.activeElement)) closeEditors(); }, 140));
  list(); placeOverlays(); inp.focus(); inp.select();
}
function placeOverlays() {
  if (!S.can) return; const ed = S.editing;
  if (ed && ed.kind === 'reagent') { const ar = S.chain.find(a => a.id === ed.id); if (ar) { const [sx, sy] = toScreen((ar.x1 + ar.x2) / 2, ar.y); const w = ed.el.offsetWidth || 240; ed.el.style.transform = `translate(${f1(clamp(sx - w / 2, 8, cw() - w - 8))}px,${f1(sy - 50)}px)`; } }
  if (ed && ed.kind === 'label') { const p = ed.id != null ? S.mol.atoms.get(ed.id) : ed.at; if (p) { const [sx, sy] = toScreen(p.x, p.y); ed.el.style.transform = `translate(${f1(sx - 60)}px,${f1(sy - 18)}px)`; } }
  const pill = S.el.querySelector('.bx-pill'); const i = S.prodHover >= 0 ? S.prodHover : S.pillHold != null ? S.pillHold : -1;
  if (pill) {
    const ar = S.chain[i]; const p = ar && ar.parts && ar.parts.find(q => q.main);
    if (p && !S.g) { const [sx, sy] = toScreen((p.box[0] + p.box[2]) / 2, p.box[1]); pill.hidden = false; pill.dataset.i = i; const nk = S.chain.filter(a => a.src && a.src.a === ar.id).length; pill.querySelector('[data-pa=next]').innerHTML = nk ? '+ Branch' : '+ Next arrow'; pill.style.transform = `translate(${f1(sx - pill.offsetWidth / 2)}px,${f1(Math.max(6, sy - 44))}px)`; }
    else pill.hidden = true;
  }
}

/* ---------------- structure-wide actions ---------------- */
function loadSmiles(smi, o = {}) {
  if (!S.rd) { S.pending = () => loadSmiles(smi, o); return false; }
  let m = null; try { m = Chem.fromSmilesMB(smi, BL); } catch (e) { m = null; }
  if (!m || !m.atoms.size) { toast('Could not read that structure'); return false; }
  const xs = [...m.atoms.values()]; const cx = (Math.min(...xs.map(a => a.x)) + Math.max(...xs.map(a => a.x))) / 2, cy = (Math.min(...xs.map(a => a.y)) + Math.max(...xs.map(a => a.y))) / 2;
  for (const a of xs) { a.x -= cx; a.y -= cy; }
  if (o.undo !== false && (S.mol.atoms.size || S.chain.length)) pushUndo();
  if (o.add && S.mol.atoms.size) { addCompound(m); return true; }
  S.mol = m; S.sel.a.clear(); S.sel.b.clear(); closeEditors(true);
  S.chain = []; (o.reagents || []).forEach((r, j) => S.chain.push(arrowObj(r, null, j ? { a: S.chain[j - 1].id } : { f: [...m.atoms.keys()] }))); S.active = S.chain.length ? 0 : -1;
  changed({ noFlight: !o.flight, fast: true }); fitView(o.anim !== false, 1.75);
  return true;
}
function tidy() {
  if (!S.mol.atoms.size) return; pushUndo();
  const from = new Map([...S.mol.atoms.values()].map(a => [a.id, [a.x, a.y]]));
  /* each compound (or just the selected ones) is redrawn in place, turned to match how it was drawn; wedges are re-derived so R/S survives */
  const back = S.mol.clone(); let ok = true; const sf = selFrags(); const frs = sf || S.mol.frags(); const St = window.Stereo;
  for (const ids of frs) {
    const sub = S.mol.sub(ids); const hadSt = [...sub.bonds.values()].some(b => b.st === 1 || b.st === 6); let smi = null; if (hadSt && St) { try { smi = St.iso(sub); } catch (e) { smi = null; } }
    try { Chem.relayout(sub, BL); } catch (e) { ok = false; break; }
    if ([...sub.atoms.values()].some(a => !isFinite(a.x) || !isFinite(a.y))) { ok = false; break; }
    const hv = ids.filter(id => sub.atoms.get(id).el !== 'H'); const T = fit2D(hv.map(id => [sub.atoms.get(id).x, sub.atoms.get(id).y]), hv.map(id => [back.atoms.get(id).x, back.atoms.get(id).y]));
    if (T) for (const a of sub.atoms.values()) { const ax = (a.x - T.px) * T.refl, ay = a.y - T.py; a.x = ax * T.c - ay * T.s + T.qx; a.y = ax * T.s + ay * T.c + T.qy; }
    let fin = sub; if (smi) { try { fin = St.wedgeTo(sub, smi) || sub; } catch (e) { fin = sub; } }
    for (const a of fin.atoms.values()) { const t = S.mol.atoms.get(a.id); t.x = a.x; t.y = a.y; }
    for (const b of fin.bonds.values()) { const t = S.mol.bonds.get(b.id); if (t) { t.st = b.st; t.a = b.a; t.b = b.b; } }
  }
  if (!ok) { S.mol = back; S.undo.pop(); toast('Fix the red-ringed atoms first, then tidy'); drawAll(); return; }
  if (!motion()) { changed({ noFlight: true }); return; }
  const to = new Map([...S.mol.atoms.values()].map(a => [a.id, [a.x, a.y]])); const t0 = performance.now();
  const tick = now => { const t = clamp((now - t0) / 480, 0, 1), e = easeIO(t); for (const a of S.mol.atoms.values()) { const p = from.get(a.id), q = to.get(a.id); a.x = lerp(p[0], q[0], e); a.y = lerp(p[1], q[1], e); } layoutScheme(); drawMol(); drawScheme(); drawHL(); if (t < 1) raf(tick); else changed({ noFlight: true }); };
  raf(tick);
}
/* a new compound below everything on the board (ids renumbered), e.g. a product copied out to start a new lane */
function addCompound(m0) {
  const b = S.scheme ? S.scheme.box : S.mol.bbox(); const mb = m0.bbox(); const dx = (b[0] + 40) - mb.x0, dy = (b[3] + 70) - mb.y0;
  const idm = S.mol.merge(m0, dx, dy); S.sel.a = new Set(idm.values()); S.sel.b.clear(); for (const bd of S.mol.bonds.values()) if (S.sel.a.has(bd.a) && S.sel.a.has(bd.b)) S.sel.b.add(bd.id);
  closeEditors(true); changed({ noFlight: true }); setTimeout(() => fitView(true), 30); return idm;
}
/* Flip turns the drawing over in the plane. Done naively that silently draws the enantiomer, so wedges and
   hashes swap: the molecule stays the same. Mirror (below) is the one that makes the enantiomer. */
function flip(axis) {
  const ids = moveSet().length ? moveSet() : [...S.mol.atoms.keys()]; if (!ids.length) return; pushUndo();
  const bb = S.mol.bbox(ids); const cx = (bb.x0 + bb.x1) / 2, cy = (bb.y0 + bb.y1) / 2; const set = new Set(ids);
  for (const id of ids) { const a = S.mol.atoms.get(id); if (axis === 'h') a.x = 2 * cx - a.x; else a.y = 2 * cy - a.y; }
  let sw = 0; for (const b of S.mol.bonds.values()) if (set.has(b.a) && set.has(b.b) && (b.st === 1 || b.st === 6)) { b.st = b.st === 1 ? 6 : 1; sw++; }
  if (sw) toast('Flipped · wedges and hashes swapped so it is still the same molecule'); changed();
}
function mirrorSel() {
  const ids = moveSet().length ? moveSet() : [...S.mol.atoms.keys()]; const set = new Set(ids); let n = 0;
  const bs = [...S.mol.bonds.values()].filter(b => set.has(b.a) && set.has(b.b) && (b.st === 1 || b.st === 6));
  if (!bs.length) { toast('No wedges or hashes here, so the mirror image is the same drawing'); return; }
  pushUndo(); for (const b of bs) { b.st = b.st === 1 ? 6 : 1; n++; } toast('Mirror image: every wedge ↔ hash (the enantiomer)'); changed();
}
/* SMILES keeps R/S and E/Z when the drawing has them */
function isoSmi(m) { const St = window.Stereo; if (!St) return null; try { const c = St.cip(m); return c.atoms.size || c.bonds.size ? St.iso(m) : null; } catch (e) { return null; } }
function copySmiles() {
  if (!S.mol.atoms.size) { toast('Nothing to copy yet'); return; }
  let s = ''; try { s = isoSmi(S.mol) || Chem.canon(S.mol) || Chem.toSmiles(S.mol).smi; } catch (e) { s = Chem.toSmiles(S.mol).smi; }
  (navigator.clipboard ? navigator.clipboard.writeText(s) : Promise.reject()).then(() => toast('Copied SMILES ' + s)).catch(() => toast('SMILES: ' + s));
}
function duplicateSel() {
  const fr = selFrags() || (S.mol.frags().length === 1 ? S.mol.frags() : null); if (!fr) { toast('Select a compound to duplicate'); return; }
  pushUndo(); const sub = S.mol.sub(fr.flat()); addCompound(sub); toast('Duplicated · draw an arrow from the copy to try another reagent');
}
function clearArrows() { if (!S.chain.length) { toast('No arrows to clear'); return; } pushUndo(); S.chain = []; S.active = -1; closeEditors(true); changed({ noFlight: true }); toast('Arrows cleared, structures kept · Ctrl+Z brings them back'); }
function clearBoard() { if (!S.mol.atoms.size && !S.chain.length) return; pushUndo(); S.mol = new Chem.Mol(); S.chain = []; S.active = -1; S.sel.a.clear(); S.sel.b.clear(); closeEditors(true); changed({ noFlight: true }); toast('Board cleared · Ctrl+Z brings it back'); }
let smiT = 0;
function scheduleSmiles() { clearTimeout(smiT); smiT = setTimeout(updateSmilesField, 220); }
function updateSmilesField() {
  const f = S.el && S.el.querySelector('.bx-smi input'); if (!f || document.activeElement === f) return;
  let s = ''; if (S.rd && S.mol.atoms.size) { try { s = isoSmi(S.mol) || Chem.canon(S.mol) || ''; } catch (e) { s = ''; } }
  f.value = s; f.parentElement.classList.toggle('bad', !!S.mol.atoms.size && !s);
}

/* ---------------- mechanism viewer: atom-level continuity between steps ---------------- */
const MV = { res: null, T: [], i: 0, shown: null, box: null, playing: false, timer: 0, anim: 0, big: false };
/* every frame is turned to match the one before it (same atom ids), so a step never spins or mirrors the drawing */
function alignClone(m, ref) { if (!m) return m; const c = m.clone(); if (ref) { try { alignTo(c, ref); } catch (e) { } } return c; }
function mvTargets(r, src) {
  const T = []; let ref = src || null;
  (r.steps || []).forEach((s, j) => { const mol = alignClone(s.L, ref); T.push({ mol, arrows: s.arrows || [], text: s.text || '', cond: s.cond || [], j }); ref = mol; });
  T.push({ mol: alignClone(r.product, ref), arrows: [], text: '', final: true, j: T.length });
  T.forEach((t, j) => { const st = (r.steps || [])[j]; if (st && st.R) t.R = alignClone(st.R, t.mol); });
  return T;
}
function mvSvg() { return S.el && S.el.querySelector('.bx-mvs'); }
function mvSet(r, src) { if (MV.res === r) return; mvStop(); MV.res = r; MV.T = mvTargets(r, src); MV.i = 0; MV.shown = null; MV.box = null; mvPaint(0, false); }
function mvStep(d) { if (MV.res && MV.res.fromLib) { NM.i += d; nmPaint(); return; } if (!MV.res) return; mvStop(); mvGo(MV.i + d); }
function mvGo(i, anim = true) { if (!MV.res) return; i = clamp(i, 0, MV.T.length - 1); if (i === MV.i && MV.shown) return; const from = MV.i; MV.i = i; mvPaint(i, anim, from); }
function mvStop() { MV.playing = false; clearTimeout(MV.timer); const b = S.el && S.el.querySelector('[data-mv=play]'); if (b) b.classList.remove('on'); }
function mvPlay() {
  if (!MV.res) return; if (MV.playing) { mvStop(); return; }
  MV.playing = true; S.el.querySelector('[data-mv=play]').classList.add('on');
  if (MV.i >= MV.T.length - 1) mvGo(0, false);
  const next = () => { if (!MV.playing) return; if (MV.i >= MV.T.length - 1) { mvStop(); return; } mvGo(MV.i + 1); MV.timer = setTimeout(next, 2300); };
  MV.timer = setTimeout(next, 1400);
}
const keepOf = T => { const k = new Set(); for (const a of (T && T.arrows) || []) { for (const id of a.s) k.add(id); for (const id of a.t) k.add(id); } return k; };
function layK(m, keep) { return Render.layout(m, { BL, fs: FS, keepH: keep }); }
function stBox(T) { const L = layK(T.mol, keepOf(T)); const ar = T.arrows.length ? Render.drawArrows(L, T.arrows) : null; return { L, ar, box: ar ? ar.box : L.box }; }
function setVB(b) { const s = mvSvg(); if (!s) return; const p = 26; s.setAttribute('viewBox', `${f1(b[0] - p)} ${f1(b[1] - p)} ${f1(b[2] - b[0] + p * 2)} ${f1(b[3] - b[1] + p * 2)}`); }
function mvStatic(T, drawIn) {
  const s = mvSvg(); if (!s) return; const { L, ar, box } = stBox(T);
  s.classList.remove('draw'); s.style.removeProperty('--e');
  s.querySelector('.mvg').innerHTML = `<g class="mv-st">${Render.svgMol(L, { cls: 'mol' })}${ar ? ar.svg : ''}</g>`;
  if (drawIn && motion()) { void s.getBoundingClientRect(); s.classList.add('draw'); }
  setVB(box); MV.shown = T.mol; MV.box = box;
}
function sigs(m) {
  const b = new Set(), a = new Map();
  for (const x of m.bonds.values()) b.add(`${Math.min(x.a, x.b)}-${Math.max(x.a, x.b)}-${x.o}-${x.st || 0}`);
  for (const x of m.atoms.values()) a.set(x.id, `${x.el}|${x.chg || 0}|${x.rad || 0}|${x.el === 'H' ? '' : m.hTot(x.id)}|${m.degree(x.id)}`);
  return { b, a };
}
function tagged(svg, m, other) {
  return svg.replace(/data-b="(\d+)"/g, (s, q) => { const x = m.bonds.get(+q); const k = x && `${Math.min(x.a, x.b)}-${Math.max(x.a, x.b)}-${x.o}-${x.st || 0}`; return x && other.b.has(k) ? s : s + ' data-x="1"'; })
    .replace(/data-a="(\d+)"/g, (s, q) => { const x = m.atoms.get(+q); const k = x && `${x.el}|${x.chg || 0}|${x.rad || 0}|${x.el === 'H' ? '' : m.hTot(x.id)}|${m.degree(x.id)}`; return x && other.a.get(+q) === k ? s : s + ' data-x="1"'; });
}
function morph(from, to, dur, fadeSvg, box0, box1, done, kf, kt) {
  const s = mvSvg(); if (!s) { done(); return; }
  const sf = sigs(from), stg = sigs(to); const t0 = performance.now(); s.classList.remove('draw');
  /* atoms that arrive grow out of the atom they bond to (or drift in from outside); atoms that leave drift away */
  const cen = m => { let x = 0, y = 0, n = 0; for (const a of m.atoms.values()) { x += a.x; y += a.y; n++; } return n ? [x / n, y / n] : [0, 0]; };
  const [cxF, cyF] = cen(from), [cxT, cyT] = cen(to); const enter = new Map(), leave = new Map();
  for (const a of to.atoms.values()) { if (from.atoms.has(a.id)) continue; const nb = to.nbrs(a.id).find(n => from.atoms.has(n));
    if (nb != null) { const p = from.atoms.get(nb); enter.set(a.id, [p.x + (a.x - to.atoms.get(nb).x) * .3, p.y + (a.y - to.atoms.get(nb).y) * .3]); }
    else { const dx = a.x - cxT, dy = a.y - cyT, l = Math.hypot(dx, dy) || 1; enter.set(a.id, [a.x + dx / l * 30, a.y + dy / l * 30]); } }
  for (const a of from.atoms.values()) { if (to.atoms.has(a.id)) continue; const dx = a.x - cxF, dy = a.y - cyF, l = Math.hypot(dx, dy) || 1; leave.set(a.id, [dx / l * 26, dy / l * 26]); }
  const frame = now => {
    if (!mvSvg()) return;
    const t = clamp((now - t0) / dur, 0, 1), e = easeIO(t);
    const mf = from.clone(), mt = to.clone();
    for (const a of mt.atoms.values()) { const p = from.atoms.get(a.id); if (p) { a.x = lerp(p.x, a.x, e); a.y = lerp(p.y, a.y, e); } else { const s0 = enter.get(a.id); if (s0) { a.x = lerp(s0[0], a.x, e); a.y = lerp(s0[1], a.y, e); } } }
    for (const a of mf.atoms.values()) { const q = to.atoms.get(a.id); if (q) { a.x = lerp(a.x, q.x, e); a.y = lerp(a.y, q.y, e); } else { const d = leave.get(a.id); if (d) { a.x += d[0] * e; a.y += d[1] * e; } } }
    const so = tagged(Render.svgMol(layK(mf, kf), { cls: 'mol' }), from, stg), sn = tagged(Render.svgMol(layK(mt, kt), { cls: 'mol' }), to, sf);
    s.style.setProperty('--e', e.toFixed(3));
    s.querySelector('.mvg').innerHTML = `<g class="mv-old">${so}</g><g class="mv-new">${sn}</g>${fadeSvg ? `<g class="mv-fade" style="opacity:${Math.max(0, 1 - t * 2.6).toFixed(3)}">${fadeSvg}</g>` : ''}`;
    setVB(box0.map((v, k) => lerp(v, box1[k], e)));
    if (t < 1) MV.anim = raf(frame); else done();
  };
  cancelAnimationFrame(MV.anim); MV.anim = raf(frame);
}
function mvPaint(i, anim, from) {
  const T = MV.T[i]; if (!T) return; mvText(i);
  if (!anim || !motion() || !MV.shown || !mvSvg()) { cancelAnimationFrame(MV.anim); mvStatic(T, !!MV.shown); return; }
  const fin = stBox(T).box; const seq = [];
  const prevT = MV.T[from]; const R = from >= 0 && from < MV.res.steps.length && i === from + 1 ? (prevT.R || MV.res.steps[from].R) : null;
  const fade = prevT && prevT.arrows.length ? stBox(prevT).ar.svg : '';
  const kPrev = keepOf(prevT), kT = keepOf(T);
  if (R) { const rb = layK(R, kT).box; seq.push([MV.shown, R, 680, fade, MV.box, rb, kPrev, kT]); seq.push([R, T.mol, 380, '', rb, fin, kT, kT]); }
  else seq.push([MV.shown, T.mol, 520, fade, MV.box, fin, kPrev, kT]);
  const run = k => { if (k >= seq.length) { mvStatic(T, true); return; } const q = seq[k]; morph(q[0], q[1], q[2], q[3], q[4] || fin, q[5], () => run(k + 1), q[6], q[7]); };
  run(0);
}
function mvText(i) {
  const T = MV.T[i]; const box = S.el.querySelector('.bx-mvt'); if (!box || !T) return; const n = MV.T.length - 1;
  const dots = S.el.querySelector('.bx-dots'); if (dots) dots.innerHTML = MV.T.map((t, j) => `<button type="button" class="${j === i ? 'on' : ''}${j < i ? ' past' : ''}${t.final ? ' fin' : ''}" data-mvi="${j}" aria-label="${t.final ? 'Product' : 'Step ' + (j + 1)}"></button>`).join('');
  S.el.querySelector('[data-mv=prev]').disabled = i <= 0; S.el.querySelector('[data-mv=next]').disabled = i >= n;
  if (T.final) { const an = analysis(T.mol); box.innerHTML = `<div class="k">Product</div><p>${an ? `${chemHTML(an.formula)}${an.mw ? ` · ${an.mw.toFixed(2)} g/mol` : ''}${an.dbe != null ? ` · DBE ${an.dbe}` : ''}` : ''}</p>`; return; }
  box.innerHTML = `<div class="k">Step ${i + 1} of ${n}${T.cond && T.cond.length ? ` <span class="cond">${T.cond.map(chemHTML).join(' · ')}</span>` : ''}</div><p>${esc(T.text || (T.arrows.length ? 'Electron movement shown by the curved arrows.' : ''))}</p>`;
}

/* ---------------- the right panel ---------------- */
const TRY = [['Propene', 'CC=C', 'HBr'], ['Toluene', 'Cc1ccccc1', 'Br2, FeBr3'], ['Acetone', 'CC(C)=O', 'CH3MgBr then H3O+'], ['Ethanal', 'CC=O', 'dil. NaOH'], ['Cyclohexanol', 'OC1CCCCC1', 'conc. H2SO4, heat'], ['2-Bromobutane', 'CCC(C)Br', 'alc. KOH'], ['Benzene', 'c1ccccc1', 'CH3COCl, AlCl3'], ['Propyne', 'CC#C', 'HgSO4, dil. H2SO4']];
function panelIndex() { if (!S.chain.length) return -1; return S.active >= 0 && S.active < S.chain.length ? S.active : S.chain.length - 1; }
function lev(a, b) { const m = a.length, n = b.length; if (!m) return n; if (!n) return m; let p = Array.from({ length: n + 1 }, (_, j) => j); for (let i = 1; i <= m; i++) { const c = [i]; for (let j = 1; j <= n; j++) c[j] = Math.min(p[j] + 1, c[j - 1] + 1, p[j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1)); p = c; } return p[n]; }
function closest(text) { const q = normRg(text); return rgCands().map(c => [lev(q, c.k), c]).filter(([d, c]) => d <= Math.max(2, q.length / 3)).sort((x, y) => x[0] - y[0]).slice(0, 5).map(x => x[1].v); }
function renderPanel() {
  const P = S.el && S.el.querySelector('.bx-panel'); if (!P) return;
  const i = panelIndex(); const tabs = P.querySelector('.bx-tabs'); const body = P.querySelector('.bx-pbody');
  tabs.innerHTML = S.chain.map((ar, j) => `<button type="button" class="bx-tab${j === i ? ' on' : ''} st-${ar.res ? ar.res.status : 'empty'}" data-tab="${j}"><b>${j + 1}</b><span>${ar.reagent ? chemHTML(ar.reagent) : '<i>no reagent</i>'}</span></button>`).join('');
  tabs.hidden = S.chain.length < 2;
  if (i < 0) { MV.res = null; if (body.dataset.k !== 'idle') { body.dataset.k = 'idle'; body.innerHTML = idleHTML(); } return; }
  const ar = S.chain[i], r = ar.res, st = r ? r.status : 'empty';
  const coK = coKey(ar);
  if (st !== 'ok') { MV.res = null; const k = st + '|' + ar.reagent + '|' + i + '|' + coK; if (body.dataset.k === k) return; body.dataset.k = k; body.innerHTML = stateHTML(ar, i, st, r) + coHTML(ar, i); return; }
  const bk = ar.key + '|' + coK + '|' + i + '|' + (ar.st ? ar.st.kind + ar.st.list.length : '');
  if (body.dataset.k === bk && MV.res === r && (mvSvg() || r.fromLib)) return;
  body.dataset.k = bk; mvStop(); MV.res = null;
  if (r.fromLib) { body.innerHTML = libBody(ar, i, r); MV.res = r; deferExtras(i); return; }
  const minor = (r.minor || []).filter(m => m && m.atoms);
  const ign = (r.unknown || []).filter(Boolean);
  body.innerHTML = `<div class="bx-res"><div class="eyebrow">Arrow ${i + 1} · ${chemHTML(ar.reagent)}</div><h3>${esc(r.title || 'Product')}</h3>${ign.length ? `<p class="bx-ign">Ignored “${ign.map(chemHTML).join('”, “')}”: not a reagent I can read.</p>` : ''}${eqHTML(i, r)}</div>
    <div class="bx-mv"><svg class="bx-mvs" preserveAspectRatio="xMidYMid meet" aria-label="Mechanism"><g class="mvg board-mol"></g></svg>
      <div class="bx-mvc"><button type="button" class="ib" data-mv="prev" aria-label="Previous step">${ICON.prev}</button><div class="bx-dots"></div><button type="button" class="ib" data-mv="next" aria-label="Next step">${ICON.next}</button><span class="sp"></span><button type="button" class="ib" data-mv="play" aria-label="Play mechanism">${ICON.play}</button><button type="button" class="ib" data-mv="big" aria-label="Enlarge">${ICON.big}</button></div></div>
    <div class="bx-mvt"></div>
    ${r.why && r.why.length ? `<ul class="bx-why">${r.why.map(x => `<li>${esc(x)}</li>`).join('')}</ul>` : ''}
    ${stereoHTML(ar, i)}${coHTML(ar, i)}
    ${minor.length ? `<section class="bx-sec"><h4>Also formed</h4><div class="bx-minor">${minor.map(m => `<div class="bx-mini">${molCard(m)}</div>`).join('')}</div></section>` : ''}
    <section class="bx-sec bx-ways"><h4>Other ways to the same product</h4><div class="bx-chipset"><span class="faint">checking…</span></div></section>
    <section class="bx-sec bx-lib"><h4>In the library</h4><div class="bx-nl"><span class="faint">looking…</span></div></section>`;
  mvSet(r, ar.srcMol); deferExtras(i);
}
/* ---------------- stereochemistry of this arrow ---------------- */
const cipText = lab => { const v = [...lab.atoms.values()]; return v.length ? '(' + v.join(',') + ')' : ''; };
const ST_KIND = { meso: 'meso', racemic: 'racemic mixture (±)', single: 'single stereoisomer', diastereomers: 'diastereomers', 'pair-meso': 'meso + racemic pair', undrawn: 'configuration not drawn', achiral: 'no stereocentres' };
const DESC = new WeakMap();
function descOf(m) { let d = DESC.get(m); if (!d) { try { d = window.Stereo.describe(m); } catch (e) { d = null; } DESC.set(m, d); } return d; }
function startLine(m) {
  const d = descOf(m); if (!d) return ''; const ez = [...d.labels.bonds.values()];
  if (!d.n && !ez.length) return 'Start: achiral, no stereocentres.';
  const parts = [];
  if (d.n) parts.push(`${d.n} stereocentre${d.n > 1 ? 's' : ''}${d.assigned ? ' ' + cipText(d.labels) : ''}`);
  if (ez.length) parts.push(ez.map(v => v === 'Z' ? 'Z (cis)' : 'E (trans)').join(', ') + ' C=C');
  const what = d.meso ? 'meso, optically inactive' : d.chiral ? 'chiral, optically active' : d.unassigned ? 'configuration not drawn (use wedges to set it)' : d.n ? 'achiral' : '';
  return `Start: ${parts.join(' · ')}${what ? ' · ' + what : ''}.`;
}
function stereoHTML(ar, i) {
  const St = window.Stereo; if (!St || !ar.srcMol) return ''; const o = ar.st; const sl = startLine(ar.srcMol); const d = descOf(ar.srcMol);
  if (!o && (!d || (!d.n && !d.labels.bonds.size))) return '';
  const list = o && o.list ? o.list.slice(0, 4) : [];
  const cards = list.length > 1 || (o && (o.kind === 'meso' || o.kind === 'single')) ? `<div class="bx-isos">${list.map((x, k) => `<div class="bx-iso"><div class="m">${molCard(x.mol)}</div><div class="l"><b>${esc(cipText(x.labels) || (o.kind === 'meso' ? 'meso' : 'isomer ' + (k + 1)))}</b>${list.length > 1 && o.kind === 'racemic' ? `<span>${k ? 'enantiomer' : '50%'}</span>` : ''}<button type="button" class="bx-fbtn" data-fis="${i}:${k}">Fischer</button></div></div>`).join('')}</div>` : '';
  const sf = d && d.assigned ? `<button type="button" class="bx-fbtn" data-fis="${i}:-1">Fischer of the start</button>` : '';
  return `<section class="bx-sec bx-stereo"><h4>Stereochemistry${o && ST_KIND[o.kind] ? ` <span class="k k-${o.kind}">${ST_KIND[o.kind]}</span>` : ''}</h4><p class="ex">${esc(sl)}${sf ? ' ' + sf : ''}</p>${o && o.note ? `<p class="ex">${esc(o.note)}</p>` : ''}${cards}<div class="bx-fis" hidden></div></section>`;
}
function showFischer(i, k) {
  const ar = S.chain[i]; if (!ar) return; const box = S.el.querySelector('.bx-stereo .bx-fis'); if (!box) return;
  const m = k < 0 ? ar.srcMol : ar.st && ar.st.list[k] && ar.st.list[k].mol; if (!m) return;
  let F = null; try { F = window.Stereo.fischer(m); } catch (e) { F = null; }
  box.hidden = false;
  box.innerHTML = F && F.rows ? `<div class="hd"><b>Fischer projection${k < 0 ? ' · start' : ''}</b><span>horizontal bonds point toward you, vertical ones away</span></div>${window.Stereo.fischerSVG(F, { labels: true })}<button type="button" class="btn sm" data-fisdraw="${i}:${k}">Draw it on the board</button>` : `<p class="ex">${esc((F && F.reason) || 'No Fischer projection for this one.')}</p>`;
  S.fisT = { i, k, F };
}
/* which drawn compounds take part (only for arrows that start on the drawing) */
function coKey(ar) { return ar.src && ar.src.f ? (S.fragBox || []).map(f => ar.src.f.includes(f.ids[0]) ? 1 : 0).join('') : ''; }
function coHTML(ar, i) {
  if (!ar.src || !ar.src.f || !S.fragBox || S.fragBox.length < 2) return '';
  const chips = S.fragBox.map((f, k) => { const on = f.ids.some(id => ar.src.f.includes(id)); const an = analysis(S.mol.sub(f.ids)); return `<button type="button" class="chip${on ? ' on' : ''}" data-co="${i}:${k}" aria-pressed="${on}"><span>${an ? chemHTML(an.formula) : 'compound ' + (k + 1)}</span></button>`; }).join('');
  return `<section class="bx-sec bx-co"><h4>Compounds in this reaction</h4><p class="ex">Each arrow reacts the compound it starts from. Add another drawn compound to make them react together.</p><div class="bx-chipset">${chips}</div></section>`;
}
function toggleCo(i, k) {
  const ar = S.chain[i]; const f = S.fragBox && S.fragBox[k]; if (!ar || !f || !ar.src || !ar.src.f) return;
  const on = f.ids.some(id => ar.src.f.includes(id)); const next = on ? ar.src.f.filter(id => !f.ids.includes(id)) : ar.src.f.concat(f.ids);
  if (!next.length) { toast('An arrow needs at least one compound'); return; }
  pushUndo(); ar.src = { f: next }; S.cache.delete(ar.key); changed({ noFlight: true });
}
/* ---------------- amount of reagent (equivalents) ---------------- */
function eqHTML(i, r) {
  if (!r || (!r.more && !r.eq)) return '';
  const asked = r.eq ? r.eq.asked : 1; const cur = asked === Infinity ? 'x' : String(asked);
  const opts = [['1', '1 eq'], ['2', '2 eq'], ['x', 'Excess']]; if (!opts.some(o => o[0] === cur)) opts.splice(2, 0, [cur, cur + ' eq']);
  const used = r.eq ? r.eq.used : 1; let note = '';
  if (!r.eq && r.more) note = `More of it reacts again: ${esc(shortTitle(r.more.title))}`;
  else if (r.eq && asked !== Infinity && used < asked) note = `Only ${used} equivalent${used > 1 ? 's' : ''} can react here; the rest stays unused`;
  else if (r.eq && used > 1 && asked === Infinity) note = `${used} equivalents reacted`;
  return `<div class="bx-eq"><span class="lbl">Amount</span><div class="seg bx-eqseg" role="radiogroup" aria-label="Amount of reagent">${opts.map(([v, t]) => `<button type="button" role="radio" aria-checked="${v === cur}" class="${v === cur ? 'on' : ''}" data-eq="${i}:${v}">${t}</button>`).join('')}</div>${note ? `<p class="bx-eqn">${note}</p>` : ''}</div>`;
}
const EQ_SPECIES = 'H2|HBr|HCl|HI|Br2|Cl2|I2|HNO3|NaOH|KOH|NaNH2|CH3I|MeI|CH3MgBr|CH3MgI|O3|NaBH4|LiAlH4|Na|NaOMe|NaOEt|CH3Cl|CH3Br';
function stripEq(t) {
  return String(t || '').replace(/\(?\s*\d+(?:\.\d+)?\s*(?:eq|equiv|equivalents?|moles?|mol)\b\.?\s*\)?/gi, ' ')
    .replace(/\(?\s*\b(?:in\s+)?(?:large\s+)?(?:excess|xs)\b\s*\)?/gi, ' ')
    .replace(new RegExp('(^|[\\s,(+/])[2-6]\\s*(?=(?:' + EQ_SPECIES + ')\\b)', 'g'), '$1')
    .replace(/\(\s*\)/g, ' ').replace(/\s{2,}/g, ' ').replace(/\s+,/g, ',').replace(/,\s*,/g, ',').replace(/^[,\s]+|[,\s]+$/g, '').trim();
}
function withEq(t, v) { const s = stripEq(t); if (v === 1) return s; const tag = v === Infinity ? '(excess)' : `(${v} eq)`; const m = s.match(/^((?:(?:alc|aq|conc|dil|liq|dry|cold|hot|anhyd)\.?\s+)*[^\s,;\/]+)/i); return m ? `${m[1]} ${tag}${s.slice(m[1].length)}` : `${s} ${tag}`; }
function setAmount(i, v) {
  const ar = S.chain[i]; if (!ar) return; const r = ar.res || {}; const RP = window.RxParse; const orig = ar.reagent || '';
  const k = r.eq ? r.eq.stage : r.more ? r.more.stage : 0;
  const norm = RP.normalize(orig); const base = norm.length === orig.length ? orig : norm;
  const P = RP.parse(orig); const live = P.stages.filter(st => st.sp.size || st.flags.size || st.T != null); const target = live[k] || P.stages[0];
  let next = withEq(base, v);
  if (P.stages.length > 1 && target) { let pos = 0, span = null; for (const st of P.stages) { const j = norm.indexOf(st.raw, pos); if (j < 0) break; if (st === target) { span = [j, j + st.raw.length]; break; } pos = j + st.raw.length; }
    if (span) next = base.slice(0, span[0]) + withEq(base.slice(span[0], span[1]), v) + base.slice(span[1]); }
  if (next === orig) return;
  pushUndo(); ar.reagent = next; ar.forceLib = false; closeEditors(true); layoutScheme(); drawAll(); scheduleSolve(20);
  toast(v === 1 ? 'One equivalent' : v === Infinity ? 'Excess reagent' : `${v} equivalents`);
}
function molCard(m) { const L = layC(m); const b = L.box; const p = 18; return `<svg class="board-mol" viewBox="${f1(b[0] - p)} ${f1(b[1] - p)} ${f1(b[2] - b[0] + 2 * p)} ${f1(b[3] - b[1] + 2 * p)}">${Render.svgMol(L, { cls: 'mol' })}</svg>`; }
function rgLinks(text) {
  const seen = new Set(); const out = [];
  for (const seg of String(text || '').split(SEP).filter((x, j) => j % 2 === 0)) { const k = normRg(seg); const e = k && A.RGI.get(k); if (!e || seen.has(k)) continue; seen.add(k); out.push(`<button type="button" class="bx-link" data-rgpage="${esc(k)}"><span>${chemHTML(e.label)}</span><em>${e.ids.size} reaction${e.ids.size === 1 ? '' : 's'} in the library</em><i>›</i></button>`); }
  return out.join('');
}
function idleHTML() {
  return `<div class="bx-idle"><div class="eyebrow">How the board works</div><ol class="bx-how"><li><b>Draw</b> a structure: <kbd>B</kbd> bonds, <kbd>R</kbd> rings, <kbd>Z</kbd> chains, hover an atom and type <kbd>O</kbd> <kbd>N</kbd> <kbd>L</kbd>…</li><li><b>Drag an arrow</b> from it with <kbd>A</kbd>, or press the arrow button.</li><li><b>Type the reagent.</b> The product, every mechanism step and the other routes appear as you type.</li></ol>
  <div class="eyebrow">Try one</div><div class="bx-try">${TRY.map(([n, s, rg], j) => `<button type="button" class="row" data-try="${j}"><span class="t">${esc(n)}</span><span class="r">${chemHTML(rg)}</span></button>`).join('')}</div>
  <p class="faint sm">Press <kbd>?</kbd> for every shortcut.</p></div>`;
}
function agreeHTML(id) { const rx = A.BY.get(id); if (!rx) return ''; return `<button type="button" class="bx-agree" data-nmopen="${id}"><span class="k">The library agrees</span><span class="l">${esc(A.chLine(rx))} · ${A.fmt(A.title(rx))}</span></button>`; }
function stateHTML(ar, i, st, r) {
  const chips = list => `<div class="bx-chipset">${list.map(v => `<button type="button" class="chip" data-rg="${esc(v)}"><span>${chemHTML(v)}</span></button>`).join('')}</div>`;
  const head = `<div class="bx-res"><div class="eyebrow">Arrow ${i + 1}</div>`;
  if (st === 'empty') return `${head}<h3>Type a reagent</h3></div><p class="ex">Click the arrow and type, or pick one. The product appears the moment the engine recognises it.</p>${chips(POPULAR.slice(0, 14))}`;
  if (st === 'nomol') return `${head}<h3>Nothing to react yet</h3></div><p class="ex">Draw a starting structure on the left of the arrow.</p>`;
  if (st === 'nosrc') return `${head}<h3>Waiting on arrow ${(r.parent >= 0 ? r.parent : i - 1) + 1}</h3></div><p class="ex">That arrow has no product yet, so this one has nothing to act on.</p>`;
  if (st === 'invalid') return `${head}<h3>Fix the structure first</h3></div><p class="ex">${esc(r.reason)}. Atoms with impossible valence are ringed in red on the board.</p>`;
  if (st === 'unknown') { const near = closest(ar.reagent); const known = (r.known || []).filter(k => !/^(metal|heat|dil|conc)$/.test(k)); return `${head}<h3>I can't read “${chemHTML(ar.reagent)}”</h3></div><p class="ex">${r.unknown && r.unknown.length ? `Not recognised: <b>${r.unknown.map(chemHTML).join(', ')}</b>.` : ''}${known.length ? ` Read so far: ${known.map(chemHTML).join(', ')}.` : ''} ${near.length ? 'Did you mean one of these?' : 'Try one of these:'}</p>${chips(near.length ? near : POPULAR.slice(0, 8))}`; }
  if (st === 'none' && r.positive) return `${head}<h3>${chemHTML(r.title || 'Positive test')}</h3></div><p class="ex">${chemHTML(r.reason || '')}</p><p class="ex">This is a colour test: the structure is not changed in a way worth drawing, so there is no product here.</p>${r.agree ? agreeHTML(r.agree) : ''}`;
  if (st === 'none') { const lk = rgLinks(ar.reagent); return `${head}<h3>No reaction</h3></div><p class="ex">${chemHTML(r.reason || 'Nothing happens with this reagent.')}</p>${r.test ? '<p class="ex">That is a test: the useful fact is that it gives no positive result here.</p>' : ''}${r.agree ? agreeHTML(r.agree) : ''}${lk ? `<section class="bx-sec"><h4>What it does react with</h4>${lk}</section>` : ''}`; }
  return `${head}<h3>Engine error</h3></div><p class="ex">${esc(r.reason || 'Something went wrong.')}</p>`;
}
/* other ways + matching reactions from the library, computed after the result settles */
let exT = 0, NIDX = null;
const COREACT = { Br: 'HBr', Cl: 'HCl', I: 'HI', BrBr: 'Br2', ClCl: 'Cl2', II: 'I2', O: 'H2O', OO: 'H2O2', 'C#N': 'HCN', '[C-]#N': 'HCN', N: 'NH3', '[OH-]': 'KOH', 'O=[N+]([O-])O': 'HNO3', 'O=S(=O)(O)O': 'H2SO4', 'C[Mg]Br': 'CH3MgBr', 'C[Mg]I': 'CH3MgI', NN: 'NH2NH2', NO: 'NH2OH', 'O=O': 'O2', '[O-][O+]=O': 'O3', 'ClC(Cl)Cl': 'CHCl3', 'O=CO': 'HCOOH', 'C=O': 'HCHO', 'CC(=O)Cl': 'CH3COCl', 'O=C=O': 'CO2' };
const WEAK = new Set(['H3O+', 'H2O', 'CCl4', 'heat', 'ether', 'EtOH']);
function canonS(smi) { try { const m = Chem.RD.get_mol(String(smi).replace(/[\/\\]/g, '').replace(/@+/g, '')); if (!m) return null; const c = m.get_smiles(); m.delete(); return c; } catch (e) { return null; } }
function fpOf(smi) { try { const m = Chem.RD.get_mol(smi); if (!m) return null; const f = m.get_morgan_fp(JSON.stringify({ radius: 2, nBits: 512 })); m.delete(); return f; } catch (e) { return null; } }
function tani(a, b) { if (!a || !b) return 0; let i = 0, u = 0; for (let k = 0; k < a.length; k++) { const x = a.charCodeAt(k) === 49, y = b.charCodeAt(k) === 49; if (x && y) i++; if (x || y) u++; } return u ? i / u : 0; }
function speciesOf(text) { const out = new Set(); try { for (const st of window.RxParse.parse(text || '').stages) { for (const k of st.sp.keys()) out.add(k); for (const f of st.flags) out.add(f); } } catch (e) { } return out; }
const heavy = s => { const t = String(s).replace(/\[[^\]]*\]/g, 'X'); return (t.match(/Cl|Br|[BCNOPSFIX]|[bcnops]/g) || []).length; };
/* the library’s main product: the heaviest carbon product that is not a reagent by-product (succinimide, Ph₃P=O, acids from peracids...) */
const hasCarbon = x => /C|c/.test(String(x).replace(/Cl|Ca|Cd|Cu|Cr|Co|Cs/g, ''));
let BYP = null;
function pickProd(list) { if (!BYP) BYP = new Set(['O=C1CCC(=O)N1', 'O=P(c1ccccc1)(c1ccccc1)c1ccccc1', 'CC(=O)O', 'CC(=O)[O-]', 'O=C(O)c1cccc(Cl)c1', 'O=C=O', 'O=C=S', 'CCS', 'C=O' + '#never'].map(canonS).filter(Boolean));
  const s = list.slice().sort((a, b) => heavy(b) - heavy(a)); const i = s.findIndex(x => hasCarbon(x) && !BYP.has(canonS(x))); if (i <= 0) return s; return [s[i], ...s.filter((_, j) => j !== i)]; }
function libIndex() {
  if (NIDX) return NIDX; NIDX = [];
  for (const r of A.R) {
    const parts = (r.rxn || '').split('>'); if (parts.length < 3) continue;
    const pm = pickMain(parts[0].split('.').filter(Boolean)); const left = [pm.main, ...pm.rest].filter(Boolean); const right = pickProd(parts[2].split('.').filter(Boolean));
    const sp = speciesOf([...(r.rl || []), ...(r.cl || [])].join(', '));
    const mid = !(r.rl || []).filter(Boolean).length ? parts[1].split('.').filter(Boolean) : [];
    for (const x of left.slice(1).concat(mid)) { const c = canonS(x); const nm = c && ((typeof CO_NAMES !== 'undefined' && CO_NAMES[c]) || COREACT[c]); if (nm) speciesOf(nm).forEach(k => sp.add(k)); }
    NIDX.push({ r, sp, main: left[0] ? canonS(left[0]) : null, prod: right[0] ? canonS(right[0]) : null, prods: right.filter(x => heavy(x) >= 1 && hasCarbon(x)).map(canonS).filter(Boolean), prodSmi: right[0] || null, extraSmi: right.slice(1).filter(x => heavy(x) >= 2 && /[Cc]/.test(x)), fp: null });
  }
  return NIDX;
}
/* When the engine has no general rule, but this exact structure + reagent is a reaction in the library,
   use the library’s product (and show the library’s own mechanism). */
const FLAGW = new Set(['heat', 'cold', 'dil', 'conc', 'aq', 'alc', 'dry', 'ether', 'excess', 'one', 'pyr', 'acetone', 'ccl4', 'polar', 'pressure', 'glycol', 'dark', 'cat', 'basecat', 'protic', 'aprotic', 'moist', 'sealed', 'rt', 'lewis', 'metal']);
const weak = x => WEAK.has(x) || FLAGW.has(x);
function libFallback(src, rg) {
  const k = canonS(Chem.canon(src) || ''); if (!k) return null; const U = speciesOf(rg); if (!U.size) return null;
  let best = null;
  for (const n of libIndex()) {
    if (n.main !== k || !n.prodSmi) continue;
    let sc = 0, hit = 0; for (const x of U) if (n.sp.has(x)) { hit++; sc += weak(x) ? .25 : 1; }
    let miss = 0; for (const x of n.sp) if (!U.has(x) && !weak(x)) miss++;
    if (!(sc >= 1 || (hit && hit === U.size && !miss))) continue;
    const s2 = sc - miss * .3; if (!best || s2 > best[0]) best = [s2, n, miss];
  }
  if (!best) return null; const n = best[1];
  if (n.prod && n.prod === n.main) { if (best[2] === 0 && !/@/.test(n.prodSmi || '')) return { noRx: true, id: n.r.id }; return null; }
  let Pm = null; try { Pm = Chem.fromSmilesMB(n.prodSmi, BL); } catch (e) { } if (!Pm || !Pm.atoms.size) return null;
  const extra = n.extraSmi.map(s => { try { return Chem.fromSmilesMB(s, BL); } catch (e) { return null; } }).filter(x => x && x.atoms.size);
  return { status: 'ok', product: Pm, extra, minor: [], steps: [], title: A.title(n.r), why: [], tags: [], fromLib: n.r.id, unknown: [] };
}
/* last resort: a rule learned from a library reaction with the same reagents and the same kind of reactive site */
function analogyFallback(src, rg) {
  const L = window.Learn; if (!L || !window.LEARNED) return null; const U = speciesOf(rg); if (!U.size) return null;
  let smi = null; try { smi = Chem.canon(src); } catch (e) { smi = null; } if (!smi) return null;
  let p = null; try { p = L.predict(smi, U); } catch (e) { p = null; } if (!p) return null;
  const rx = A.BY.get(p.rule.id); if (!rx) return null;
  let Pm = null; try { Pm = Chem.fromSmilesMB(p.prod, BL); } catch (e) { Pm = null; } if (!Pm || !Pm.atoms.size) return null;
  const minor = p.alts.map(s => { try { return Chem.fromSmilesMB(s, BL); } catch (e) { return null; } }).filter(x => x && x.atoms.size);
  return { status: 'ok', product: Pm, extra: [], minor, steps: [], title: 'By analogy: ' + A.title(rx).replace(/<[^>]+>/g, ''), why: [], tags: [], fromLib: rx.id, analog: { n: p.rules.length, ids: p.rules.map(t => t.id) }, unknown: [] };
}
const NM = { id: null, i: 0 };
/* the library’s pre-drawn steps are a whole row (start → product); the Board shows each step's starting structure with its
   arrows, cropped and large like its own mechanisms, then the final product as the last view */
function cropStep(svg, part) {
  const m = svg.match(new RegExp(`<svg class='${part}' x='(-?[\\d.]+)' y='(-?[\\d.]+)' width='([\\d.]+)' height='([\\d.]+)'`)); if (!m) return svg.replace("class='step'", "class='step anim'");
  let [x, y, w, h] = m.slice(1).map(Number); let x2 = x + w, y2 = y + h;
  if (part === 'ls') { const g = svg.match(/<g transform='translate\((-?[\d.]+),(-?[\d.]+)\)'>([\s\S]*?)<\/g>/);
    if (g) { const tx = +g[1], ty = +g[2]; for (const z of g[3].matchAll(/(-?\d+(?:\.\d+)?),(-?\d+(?:\.\d+)?)/g)) { const X = +z[1] + tx, Y = +z[2] + ty; x = Math.min(x, X); x2 = Math.max(x2, X); y = Math.min(y, Y); y2 = Math.max(y2, Y); } } }
  const pad = 10; const vb = `${(x - pad).toFixed(1)} ${(y - pad).toFixed(1)} ${(x2 - x + 2 * pad).toFixed(1)} ${(y2 - y + 2 * pad).toFixed(1)}`;
  return svg.replace(/<svg viewBox='[^']*' class='step'>/, `<svg viewBox='${vb}' class='step anim nm-${part}'>`);
}
function nmViews(id) { const ms = A.ART.ms[id] || []; const v = ms.map(s => ({ svg: cropStep(s, 'ls'), fin: false })); if (ms.length) v.push({ svg: cropStep(ms[ms.length - 1], 'ps'), fin: true }); return v; }
function nmPaint() {
  const V = nmViews(NM.id); const box = S.el.querySelector('.bx-nm'); if (!box || !V.length) return; NM.i = clamp(NM.i, 0, V.length - 1);
  box.innerHTML = V[NM.i].svg; S.el.querySelectorAll('[data-nmi]').forEach((b, j) => { b.classList.toggle('on', j === NM.i); b.classList.toggle('past', j < NM.i); });
  const rx = A.BY.get(NM.id); const cl = rx && rx.mcl && rx.mcl[NM.i] ? rx.mcl[NM.i].filter(Boolean).join(', ') : '';
  const t = S.el.querySelector('.bx-nmt'); if (t) t.textContent = V[NM.i].fin ? 'Product' : `Step ${NM.i + 1} of ${V.length - 1}${cl ? ' · ' + cl : ''}`;
}
function libBody(ar, i, r) {
  const rx = A.BY.get(r.fromLib); const ms = A.ART.ms[rx.id] || []; NM.id = rx.id; NM.i = 0; const V = nmViews(rx.id);
  if (r.analog) return `<div class="bx-res"><div class="eyebrow">Arrow ${i + 1} · ${chemHTML(ar.reagent)} · by analogy</div><h3>Predicted by analogy</h3><p class="bx-src">The engine has no general rule here and the library don't have this exact structure, so the board applied the pattern of <b>${A.fmt(A.title(rx))}</b> (${esc(A.chLine(rx))}): same reagents, same reactive site${r.analog.n > 1 ? `, and ${r.analog.n - 1} more reaction${r.analog.n > 2 ? 's' : ''} in the library agree` : ''}. Worth a check. Below is how that reaction goes.</p></div>
    <div class="bx-mv"><div class="bx-nm">${V.length ? V[0].svg : (A.ART.rx[rx.id] || '')}</div>
    <div class="bx-mvc">${V.length ? `<button type="button" class="ib" data-nm="-1" aria-label="Previous step">${ICON.prev}</button><div class="bx-dots">${V.map((v, j) => `<button type="button" data-nmi="${j}" class="${j ? '' : 'on'}${v.fin ? ' fin' : ''}" aria-label="${v.fin ? 'Product' : 'Step ' + (j + 1)}"></button>`).join('')}</div><button type="button" class="ib" data-nm="1" aria-label="Next step">${ICON.next}</button>` : ''}<span class="sp"></span><button type="button" class="btn sm" data-nmopen="${rx.id}">Open that reaction</button></div></div>
    ${(r.minor || []).length ? `<section class="bx-sec"><h4>Other rules suggest</h4><div class="bx-minor">${r.minor.map(m => `<div class="bx-mini">${molCard(m)}</div>`).join('')}</div></section>` : ''}`;
  return `<div class="bx-res"><div class="eyebrow">Arrow ${i + 1} · ${chemHTML(ar.reagent)} · from the library</div><h3>${A.fmt(A.title(rx))}</h3><p class="bx-src">${esc(A.chLine(rx))}. ${ar.forceLib ? 'You chose the library’s recorded product over the engine’s answer (Ctrl+Z switches back).' : 'The engine has no general rule for this one, so the product comes straight from the library.'}</p></div>
    <div class="bx-mv"><div class="bx-nm">${V.length ? V[0].svg : (A.ART.rx[rx.id] || '')}</div>
    <div class="bx-mvc">${V.length ? `<button type="button" class="ib" data-nm="-1" aria-label="Previous step">${ICON.prev}</button><div class="bx-dots">${V.map((v, j) => `<button type="button" data-nmi="${j}" class="${j ? '' : 'on'}${v.fin ? ' fin' : ''}" aria-label="${v.fin ? 'Product' : 'Step ' + (j + 1)}"></button>`).join('')}</div><button type="button" class="ib" data-nm="1" aria-label="Next step">${ICON.next}</button>` : `<span class="faint sm" style="padding:0 6px">no mechanism recorded</span>`}<span class="sp"></span><span class="faint sm bx-nmt">${V.length ? `Step 1 of ${V.length - 1}` : ''}</span><button type="button" class="ib" data-mv="big" aria-label="Enlarge">${ICON.big}</button><button type="button" class="ib" data-nmopen="${rx.id}" aria-label="Open in the viewer" title="Open in the viewer">↗</button></div></div>
    <section class="bx-sec bx-ways"><h4>Other ways to the same product</h4><div class="bx-chipset"><span class="faint">checking…</span></div></section>
    <section class="bx-sec bx-lib"><h4>In the library</h4><div class="bx-nl"><span class="faint">looking…</span></div></section>`;
}
function deferExtras(i) {
  clearTimeout(exT);
  exT = setTimeout(() => {
    const ar = S.chain[i]; if (!ar || !ar.res || ar.res.status !== 'ok' || panelIndex() !== i) return;
    const r = ar.res; const src = ar.srcMol || S.mol;
    if (!r.fromLib) { let fb = null; try { fb = libFallback(src, ar.reagent); } catch (e) { fb = null; }
      if (fb && fb.noRx) { const rx = A.BY.get(fb.id); const box = S.el.querySelector('.bx-res'); if (rx && box && !box.querySelector('.bx-diff')) box.insertAdjacentHTML('beforeend', `<div class="bx-diff solo"><div class="t"><b>The library records no reaction here</b><span>${esc(A.chLine(rx))} · ${A.fmt(A.title(rx))}</span></div><button type="button" class="btn sm" data-nmopen="${fb.id}">Open in the library</button></div>`); }
      else if (fb) { const kk = x => canonS(Chem.canon(x) || ''); const mine = new Set([r.product, ...(r.extra || [])].filter(Boolean).map(kk)); const theirs = kk(fb.product);
        const ent = libIndex().find(e => e.r.id === fb.fromLib); const anyMine = ent && ent.prods && ent.prods.some(p => mine.has(p));
        if (theirs && !mine.has(theirs) && !anyMine) { const rx = A.BY.get(fb.fromLib); const box = S.el.querySelector('.bx-res'); if (box && !box.querySelector('.bx-diff')) box.insertAdjacentHTML('beforeend', `<div class="bx-diff"><div class="t"><b>The library records a different product here</b><span>${esc(A.chLine(rx))} · ${A.fmt(A.title(rx))}</span></div><div class="m">${molCard(fb.product)}</div>${rx && rx.dn ? `<p class="why">${A.fmt(rx.dn)}</p>` : ''}<button type="button" class="btn sm" data-uselib="${i}">Use the library’s product</button></div>`); } } }
    let ways = S.waysCache.get(ar.key);
    if (!ways) { try { ways = Solve.ways(src, r.product).filter(w => normRg(w.reagent) !== normRg(ar.reagent) && w.title !== r.title).slice(0, 8); } catch (e) { ways = []; } S.waysCache.set(ar.key, ways); }
    const wb = S.el.querySelector('.bx-ways .bx-chipset');
    if (wb) wb.innerHTML = ways.length ? ways.map(w => `<button type="button" class="chip" data-rg="${esc(w.reagent)}" title="${esc(w.title || w.name || '')}"><span>${chemHTML(w.reagent)}</span></button>`).join('') : '<span class="faint">No other reagent in the engine gives exactly this product.</span>';
    let top = S.libCache && S.libCache.get(ar.key);
    if (!top) {
      const pk = canonS(Chem.canon(r.product) || ''); const sk = canonS(Chem.canon(src) || ''); const sfp = sk ? fpOf(sk) : null;
      const U = speciesOf(ar.reagent); const scored = [];
      for (const n of libIndex()) {
        let ov = 0; for (const k of U) if (n.sp.has(k)) ov += WEAK.has(k) ? .25 : 1;
        const sameP = pk && n.prod === pk, sameS = sk && n.main === sk;
        if (!ov && !sameP) continue;
        if (n.fp == null && n.main) n.fp = fpOf(n.main) || '';
        const sim = tani(sfp, n.fp);
        const score = (sameP ? 3 : 0) + ov + sim * 1.6 + (sameS ? .6 : 0);
        if (!sameP && (ov < 1 || sim < .08)) continue;
        scored.push([score, n.r, sameP && ov >= 1 ? 'same product and reagent' : sameP ? 'same product' : sameS ? 'same reagent, same substrate' : `same reagent · ${Math.round(sim * 100)}% similar substrate`]);
      }
      scored.sort((a, b) => b[0] - a[0]); top = scored.slice(0, 4);
      (S.libCache || (S.libCache = new Map())).set(ar.key, top);
    }
    S.noteIds = top.map(x => x[1].id);
    const nb = S.el.querySelector('.bx-lib .bx-nl');
    if (nb) { nb.innerHTML = top.length ? top.map(([, rx, why]) => A.miniHTML(rx, why)).join('') : '<span class="faint">Nothing in the library uses this reagent on a structure like this.</span>'; A.hydrate(nb); }
  }, 380);
}
function updateStatus() {
  if (!S.el) return; const st = S.el.querySelector('.bx-status'); if (!st) return;
  const an = analysis(S.mol); const sel = S.sel.a.size;
  st.querySelector('.an').innerHTML = an ? `<b>${chemHTML(an.formula)}</b>${an.mw ? `<span>${an.mw.toFixed(2)} g/mol</span>` : ''}${an.dbe != null && an.dbe >= 0 && Number.isInteger(an.dbe) ? `<span>DBE ${an.dbe}</span>` : ''}${an.frags > 1 ? `<span>${an.frags} fragments</span>` : ''}${an.bad.length ? `<span class="warn">${esc(badText(S.mol, an.bad))}</span>` : ''}${sel ? `<span class="acc">${sel} selected</span>` : ''}` : '<span>Empty board</span>';
  st.querySelector('.hint').innerHTML = hintText();
}
function hintText() {
  const h = S.hover;
  if (h && h.k === 'atom') return '<kbd>C</kbd><kbd>N</kbd><kbd>O</kbd><kbd>L</kbd>… replace · <kbd>+</kbd><kbd>−</kbd> charge · <kbd>1</kbd>–<kbd>3</kbd> grow · <kbd>↵</kbd> label · <kbd>Del</kbd>';
  if (h && h.k === 'bond') return '<kbd>1</kbd><kbd>2</kbd><kbd>3</kbd> order · <kbd>W</kbd> wedge · <kbd>H</kbd> hash · <kbd>Del</kbd> remove';
  if (h && h.k === 'arrow') return 'Click to edit the reagent · <kbd>Del</kbd> removes this arrow and the ones that follow from it';
  if (h && h.k === 'prod') return 'Product of this arrow · add the next arrow, branch off, or copy it out as a new compound';
  return ({ select: 'Drag atoms to move · drag empty space to box-select · <kbd>Shift</kbd> adds · <kbd>Space</kbd>-drag pans', bond: 'Click an atom to grow · drag for direction (30°, <kbd>Shift</kbd> 15°) · click a bond to change it', chain: 'Drag to draw a zig-zag chain; the bubble counts new carbons', ring: 'Click empty space, an atom (spiro) or a bond (fused) · <kbd>1</kbd> benzene, <kbd>3</kbd>–<kbd>8</kbd> sizes', atom: `Click an atom to make it ${chemHTML(S.atomLabel)} · double-click any atom to type a label`, charge: 'Click an atom to add charge · <kbd>Shift</kbd>-click the opposite sign', erase: 'Click or scrub across atoms and bonds', arrow: 'Drag from the structure to add a reaction arrow, then type the reagent' })[S.tool] || '';
}

/* ---------------- icons ---------------- */
const sv = (d, extra = '') => `<svg viewBox="0 0 20 20" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"${extra}>${d}</svg>`;
const ICON = {
  select: sv('<path d="M5.2 3.2 15 9.4l-4.4 1.1L8.2 15z" fill="currentColor" fill-opacity=".14"/>'),
  bond: sv('<path d="M4 15.5 16 4.5"/><circle cx="4" cy="15.5" r="1.4" fill="currentColor"/><circle cx="16" cy="4.5" r="1.4" fill="currentColor"/>'),
  chain: sv('<path d="M2.5 13.5 6.3 7l3.8 6.5L13.9 7l3.8 6.5"/>'),
  ring: sv('<path d="M10 2.8 16.2 6.4v7.2L10 17.2 3.8 13.6V6.4z"/><path d="M10 5.6 13.8 7.8M13.8 12.2 10 14.4M6.2 12.2V7.8" opacity=".55"/>'),
  atom: sv('<circle cx="10" cy="10" r="7.2"/><path d="M7.4 13.4V6.6l5.2 6.8V6.6" stroke-width="1.5"/>'),
  charge: sv('<circle cx="10" cy="10" r="7.2"/><path d="M10 6.6v6.8M6.6 10h6.8"/>'),
  erase: sv('<path d="m8.4 16.2-4.6-4.6a1.6 1.6 0 0 1 0-2.3l5.7-5.7a1.6 1.6 0 0 1 2.3 0l4.6 4.6a1.6 1.6 0 0 1 0 2.3l-5.7 5.7"/><path d="M8.2 16.4h8.3M6.4 6.9l6.9 6.9"/>'),
  arrow: sv('<path d="M2.8 10h12.6"/><path d="M12.4 6.2 16.6 10l-4.2 3.8" fill="currentColor" fill-opacity=".9"/>'),
  undo: sv('<path d="M7.4 5.2 3.8 8.8l3.6 3.6"/><path d="M4.2 8.8h7.4a4.6 4.6 0 0 1 0 9.2H9"/>'),
  redo: sv('<path d="m12.6 5.2 3.6 3.6-3.6 3.6"/><path d="M15.8 8.8H8.4a4.6 4.6 0 0 0 0 9.2H11"/>'),
  clear: sv('<path d="M4.5 6h11M8 6V4.2h4V6M6 6l.8 10.2h6.4L14 6"/>'),
  tidy: sv('<path d="M10 3v2.4M10 14.6V17M3 10h2.4M14.6 10H17M5.1 5.1l1.7 1.7M13.2 13.2l1.7 1.7M5.1 14.9l1.7-1.7M13.2 6.8l1.7-1.7"/>'),
  fliph: sv('<path d="M10 3v14" stroke-dasharray="1.6 2.2"/><path d="M7.6 6 3.2 14h4.4zM12.4 6l4.4 8h-4.4z"/>'),
  flipv: sv('<path d="M3 10h14" stroke-dasharray="1.6 2.2"/><path d="M6 7.6 14 3.2v4.4zM6 12.4l8 4.4v-4.4z"/>'),
  copy: sv('<rect x="6.6" y="6.6" width="9.6" height="9.6" rx="2"/><path d="M13.4 4.6V4.2a1.4 1.4 0 0 0-1.4-1.4H5a1.4 1.4 0 0 0-1.4 1.4V12A1.4 1.4 0 0 0 5 13.4h.4"/>'),
  fit: sv('<path d="M3.5 7.5v-4h4M16.5 7.5v-4h-4M3.5 12.5v4h4M16.5 12.5v4h-4"/>'),
  help: sv('<circle cx="10" cy="10" r="7.2"/><path d="M7.9 7.8a2.2 2.2 0 1 1 3 2.1c-.6.3-.9.8-.9 1.4v.5"/><circle cx="10" cy="14.2" r=".5" fill="currentColor"/>'),
  prev: sv('<path d="M12.4 4.8 7.2 10l5.2 5.2"/>'), next: sv('<path d="M7.6 4.8 12.8 10l-5.2 5.2"/>'),
  play: sv('<path d="M6.6 4.6v10.8L15.4 10z" fill="currentColor" fill-opacity=".9"/>'), big: sv('<path d="M11.8 3.6h4.6v4.6M8.2 16.4H3.6v-4.6M16.4 3.6l-5.2 5.2M3.6 16.4l5.2-5.2"/>'),
  mirror: sv('<path d="M10 2.5v15" stroke-dasharray="1.6 2.2"/><path d="M7.6 6.2 3.4 10l4.2 3.8z" fill="currentColor" fill-opacity=".9"/><path d="M12.4 6.2 16.6 10l-4.2 3.8" />'),
  fischer: sv('<path d="M10 2.8v14.4M4 7h12M4 13h12"/>'), side: sv('<rect x="3" y="4" width="14" height="12" rx="2"/><path d="M12.4 4v12"/>'), stereo: sv('<path d="M10 3.2l1.9 4.1 4.4.5-3.3 3 .9 4.4L10 13l-3.9 2.2.9-4.4-3.3-3 4.4-.5z"/>'),
  plus: sv('<path d="M10 4.5v11M4.5 10h11"/>'), minus: sv('<path d="M4.5 10h11"/>'), tpl: sv('<path d="M3.6 5.2h12.8M3.6 10h12.8M3.6 14.8h8"/>'),
};
const ringIcon = n => { if (n === 'benz') return sv('<path d="M10 2.8 16.2 6.4v7.2L10 17.2 3.8 13.6V6.4z"/><circle cx="10" cy="10" r="3.3"/>'); const k = +n; const R = 7.2; let p = ''; for (let i = 0; i < k; i++) { const t = -Math.PI / 2 + (k === 4 || k === 8 ? Math.PI / k : 0) + i * TAU / k; p += (i ? 'L' : 'M') + f1(10 + R * Math.cos(t)) + ',' + f1(10.4 + R * Math.sin(t)); } return sv(`<path d="${p}z"/>`); };
const bondIcon = t => sv({ single: '<path d="M4 15.5 16 4.5"/>', double: '<path d="M3 13.8 14.4 3.4M5.6 16.6 17 6.2"/>', triple: '<path d="M2.6 13 13.6 3M4.3 14.9 15.6 4.5M6 16.8 17 6.8"/>', wedge: '<path d="M3.6 16.4 15 3.6l2 2z" fill="currentColor"/>', hash: '<path d="M4.6 15.3l.6.7M7 12.6l1.3 1.3M9.4 10l1.9 1.9M11.8 7.4l2.6 2.5M14.2 4.8l3 3"/>' }[t]);

const TOOLS = [['select', 'Select & move', 'V'], ['bond', 'Bond', 'B'], ['chain', 'Chain', 'Z'], ['ring', 'Ring', 'R'], ['atom', 'Atom label', 'T'], ['charge', 'Charge', 'Q'], ['erase', 'Eraser', 'E'], ['arrow', 'Reaction arrow', 'A']];
const TEMPLATES = [
  ['Alkenes & alkynes', [['Propene', 'CC=C'], ['2-Methylpropene', 'CC(C)=C'], ['But-2-ene', 'CC=CC'], ['Cyclohexene', 'C1CCC=CC1'], ['Buta-1,3-diene', 'C=CC=C'], ['Propyne', 'CC#C'], ['But-1-yne', 'CCC#C']]],
  ['Halides, alcohols, ethers', [['Bromoethane', 'CCBr'], ['2-Bromobutane', 'CCC(C)Br'], ['tert-Butyl bromide', 'CC(C)(C)Br'], ['Ethanol', 'CCO'], ['Propan-2-ol', 'CC(C)O'], ['tert-Butanol', 'CC(C)(C)O'], ['Cyclohexanol', 'OC1CCCCC1'], ['tert-Butyl methyl ether', 'COC(C)(C)C']]],
  ['Carbonyls & acid derivatives', [['Methanal', 'C=O'], ['Ethanal', 'CC=O'], ['Propanone', 'CC(C)=O'], ['Benzaldehyde', 'O=Cc1ccccc1'], ['Acetophenone', 'CC(=O)c1ccccc1'], ['Ethanoic acid', 'CC(=O)O'], ['Ethanoyl chloride', 'CC(=O)Cl'], ['Ethyl ethanoate', 'CC(=O)OCC'], ['Ethanamide', 'CC(N)=O'], ['Ethanenitrile', 'CC#N']]],
  ['Aromatics & amines', [['Benzene', 'c1ccccc1'], ['Toluene', 'Cc1ccccc1'], ['Phenol', 'Oc1ccccc1'], ['Aniline', 'Nc1ccccc1'], ['Nitrobenzene', 'O=[N+]([O-])c1ccccc1'], ['Benzoic acid', 'O=C(O)c1ccccc1'], ['Chlorobenzene', 'Clc1ccccc1'], ['Benzenediazonium', 'N#[N+]c1ccccc1'], ['Ethanamine', 'CCN'], ['Isobutane', 'CC(C)C']]],
];
const QUICK_ATOMS = ['C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'OH', 'NH2', 'NO2', 'CN', 'CHO', 'COOH', 'OMe', 'Ph'];

function setTool(t) {
  if (S.tool === t) return; S.tool = t; S.ghost = null; closeEditors(true);
  $$('.bx-rail [data-tool]', S.el).forEach(b => b.classList.toggle('on', b.dataset.tool === t));
  renderCtx(); hoverAgain(); updateStatus(); drawHL();
}
function seg(items, cur, attr) { return `<div class="bx-seg">${items.map(([v, html, tip]) => `<button type="button" class="${String(v) === String(cur) ? 'on' : ''}" data-${attr}="${v}" title="${esc(tip || '')}">${html}</button>`).join('')}</div>`; }
function renderCtx() {
  const o = S.el.querySelector('.bx-opts'); if (!o) return; const T = S.tool; let h = `<span class="bx-tn">${esc(TOOLS.find(x => x[0] === T)[1])}</span>`;
  if (T === 'select') h += `<div class="bx-acts"><button type="button" class="bx-b" data-act="tidy" title="Tidy structure · Ctrl+Shift+K">${ICON.tidy}<span>Tidy</span></button><button type="button" class="bx-b" data-act="fliph" title="Flip horizontally">${ICON.fliph}<span>Flip</span></button><button type="button" class="bx-b" data-act="flipv" title="Flip vertically">${ICON.flipv}<span>Flip</span></button><button type="button" class="bx-b" data-act="mirror" title="Mirror image: swap every wedge and hash (makes the enantiomer)">${ICON.mirror}<span>Mirror</span></button><button type="button" class="bx-b" data-act="fischer" title="Redraw the compound as a Fischer projection (wedges horizontal, hashes vertical)">${ICON.fischer}<span>Fischer</span></button><button type="button" class="bx-b" data-act="dup" title="Duplicate the selected compound · Ctrl+D">${ICON.copy}<span>Duplicate</span></button><button type="button" class="bx-b" data-act="clrarrows" title="Remove every arrow, keep the structures"><span>Clear arrows</span></button><button type="button" class="bx-b" data-act="all" title="Select all · Ctrl+A"><span>Select all</span></button></div>`;
  if (T === 'bond') h += seg([['single', bondIcon('single'), 'Single · 1'], ['double', bondIcon('double'), 'Double · 2'], ['triple', bondIcon('triple'), 'Triple · 3'], ['wedge', bondIcon('wedge'), 'Wedge · W'], ['hash', bondIcon('hash'), 'Hash · H']], S.bondType, 'bt');
  if (T === 'ring') h += seg([['benz', ringIcon('benz'), 'Benzene · 1'], ['3', ringIcon(3), '3 · 3'], ['4', ringIcon(4), '4 · 4'], ['5', ringIcon(5), '5 · 5'], ['6', ringIcon(6), '6 · 6'], ['7', ringIcon(7), '7 · 7'], ['8', ringIcon(8), '8 · 8']], S.ringType, 'rt');
  if (T === 'atom') h += seg(QUICK_ATOMS.map(x => [x, `<span>${chemHTML(x)}</span>`, x]), S.atomLabel, 'al');
  if (T === 'charge') h += seg([['+', ICON.plus, 'Positive (Shift-click: negative)'], ['-', ICON.minus, 'Negative'], ['rad', '<span class="dot">•</span>', 'Radical']], S.chargeMode, 'cm');
  if (T === 'chain') h += '<span class="faint">drag from an atom or empty space</span>';
  if (T === 'erase') h += '<span class="faint">click, or hold and scrub</span>';
  if (T === 'arrow') h += `<div class="bx-quick">${POPULAR.slice(0, 9).map(v => `<button type="button" class="chip" data-rg="${esc(v)}"><span>${chemHTML(v)}</span></button>`).join('')}</div>`;
  o.innerHTML = h;
}
function helpOpen() { const h = S.el && S.el.querySelector('.bx-help'); return h && !h.hidden; }
function toggleHelp(v) { const h = S.el.querySelector('.bx-help'); h.hidden = v == null ? !h.hidden : !v; }
function emptyState() { const e = S.el && S.el.querySelector('.bx-empty'); if (e) e.hidden = !!(S.mol && S.mol.atoms.size) || !S.rd; const l = S.el && S.el.querySelector('.bx-load'); if (l) l.hidden = S.rd; }
function useReagentChip(v) {
  const i = panelIndex(); const ar = S.chain[i];
  const sf = selFrags(); if (sf && !(ar && ar.src && ar.src.f && sf.every(c => c.some(id => ar.src.f.includes(id))))) { addArrow(v, null, { f: sf.flat() }); return; }
  if (ar && (!ar.res || ar.res.status !== 'ok' || S.tool !== 'arrow')) { const snap = snapshot(); ar.reagent = v; S.undo.push(snap); S.redo.length = 0; closeEditors(true); layoutScheme(); drawAll(); scheduleSolve(20); return; }
  addArrow(v);
}

/* ---------------- stereo marks: ★ on stereocentres, R/S, E/Z (saved) ---------------- */
function loadStOpt() { try { return Object.assign({ star: true, rs: false, ez: false, wedge: true }, JSON.parse(localStorage.getItem('ocr3-stereo') || '{}')); } catch (e) { return { star: true, rs: false, ez: false, wedge: true }; } }
const getLS = (k, d) => { try { const v = localStorage.getItem(k); return v == null ? d : v; } catch (e) { return d; } };
const setLS = (k, v) => { try { localStorage.setItem(k, v); } catch (e) { } };
function togglePanel(v) { const on = v == null ? getLS('ocr3-bxpanel', 'on') === 'off' : v; setLS('ocr3-bxpanel', on ? 'on' : 'off'); S.el.querySelector('.bx').classList.toggle('nopanel', !on); applyCam(); const p = S.el.querySelector('.bx-stpop'); if (p) p.innerHTML = stPopHTML(); const b = S.el.querySelector('[data-act=panel]'); if (b) b.classList.toggle('on', on); if (on) renderPanel(); }
S.stOpt = loadStOpt();
function stPopHTML() { const O = S.stOpt; const tg = (attr, k, on, t, sub) => `<button type="button" class="tg${on ? ' on' : ''}" data-${attr}="${k}" role="switch" aria-checked="${!!on}"><i></i><span>${t}<em>${sub}</em></span></button>`;
  const row = (k, t, sub) => tg('sto', k, O[k], t, sub);
  return `<div class="gh">Board</div>${tg('vw', 'anim', getLS('ocr3-bxanim', 'on') !== 'off', 'Animations', 'products glide in, mechanism steps morph')}${tg('vw', 'panel', getLS('ocr3-bxpanel', 'on') !== 'off', 'Mechanism sidebar', 'the panel on the right · M')}
  <div class="gh">Stereochemistry</div>${row('wedge', 'Wedges & hashes', 'off draws stereo bonds as plain lines')}${row('star', 'Stereocentres *', 'star every chiral carbon · *')}${row('rs', 'R / S labels', 'CIP configuration of each centre')}${row('ez', 'E / Z labels', 'on stereogenic double bonds')}`; }
function setView(k) { if (k === 'anim') { const off = getLS('ocr3-bxanim', 'on') !== 'off'; setLS('ocr3-bxanim', off ? 'off' : 'on'); toast(off ? 'Board animations off' : 'Board animations on'); } if (k === 'panel') togglePanel(); const p = S.el.querySelector('.bx-stpop'); if (p) p.innerHTML = stPopHTML(); }
function setStOpt(k) { S.stOpt[k] = S.stOpt[k] === false ? true : !S.stOpt[k]; if (k === 'wedge') { layoutScheme(); drawAll(); } try { localStorage.setItem('ocr3-stereo', JSON.stringify(S.stOpt)); } catch (e) { } const p = S.el.querySelector('.bx-stpop'); if (p) p.innerHTML = stPopHTML(); drawMol(); drawScheme(); }
function fischerSel() {
  const fr = selFrags() || (S.mol.frags().length === 1 ? S.mol.frags() : null); if (!fr) { toast('Select the compound to redraw first'); return; }
  const ids = fr[0]; const sub = S.mol.sub(ids); let F = null; try { F = window.Stereo.fischer(sub); } catch (e) { F = null; }
  if (!F || !F.mol || !F.ok) { toast((F && F.reason) || 'No Fischer projection for this one'); return; }
  pushUndo(); placeInto(ids, F.mol); changed({ noFlight: true }); toast('Fischer projection: wedges across, hashes up and down');
}
/* swap a compound's drawing for another drawing of it (same atom ids, maybe extra explicit H), keeping its place */
function placeInto(ids, g) {
  const b0 = S.mol.bbox(ids); const gb = g.bbox(); const dx = (b0.x0 + b0.x1) / 2 - (gb.x0 + gb.x1) / 2, dy = (b0.y0 + b0.y1) / 2 - (gb.y0 + gb.y1) / 2;
  for (const a of g.atoms.values()) { const t = S.mol.atoms.get(a.id); if (t) { t.x = a.x + dx; t.y = a.y + dy; t.hFix = a.hFix; } else S.mol.addAtom(a.el, a.x + dx, a.y + dy, { id: a.id, hFix: a.hFix, chg: a.chg }); }
  for (const b of g.bonds.values()) { let t = S.mol.bondBetween(b.a, b.b); if (!t) { const id = S.mol.addBond(b.a, b.b, b.o); t = S.mol.bonds.get(id); } t.st = b.st; t.a = b.a; t.b = b.b; t.o = b.o; }
}
function fischerDraw(i, k) {
  const ar = S.chain[i]; if (!ar) return; const m = k < 0 ? ar.srcMol : ar.st && ar.st.list[k] && ar.st.list[k].mol; if (!m) return;
  let F = null; try { F = window.Stereo.fischer(m); } catch (e) { F = null; } if (!F || !F.mol) { toast('No Fischer projection for this one'); return; }
  pushUndo();
  if (k < 0 && ar.src && ar.src.f) { placeInto([...m.atoms.keys()].filter(id => S.mol.atoms.has(id)), F.mol); changed({ noFlight: true }); toast('Redrawn as a Fischer projection'); }
  else { addCompound(F.mol.clone()); toast('Fischer projection added as a new compound'); }
}
const CIPC = new WeakMap();
function stMarks(m, dx, dy, live) {
  const O = S.stOpt; if (!O || (!O.star && !O.rs && !O.ez)) return ''; const St = window.Stereo; if (!St || !S.rd || !m || !m.atoms.size) return '';
  let c = live ? null : CIPC.get(m); if (!c) { try { c = St.cip(m); } catch (e) { c = null; } if (c && !live) CIPC.set(m, c); } if (!c) return '';
  let o = '';
  for (const id of c.all) { const a = m.atoms.get(id); if (!a) continue; const ang = Chem.freeAngle(atomDirs(m, id), a, m); const lab = c.atoms.get(id); const ux = Math.cos(ang), uy = Math.sin(ang);
    if (O.star) o += `<text class="bx-star" x="${f1(a.x + dx + ux * 12)}" y="${f1(a.y + dy + uy * 12 + 5)}">*</text>`;
    if (O.rs && lab) o += `<text class="bx-cip" x="${f1(a.x + dx + ux * (O.star ? 25 : 18))}" y="${f1(a.y + dy + uy * (O.star ? 25 : 18) + 4)}">${esc(lab)}</text>`; }
  if (O.ez) for (const [bid, v] of c.bonds) { const b = m.bonds.get(bid); if (!b) continue; const p = m.atoms.get(b.a), q = m.atoms.get(b.b); const mx = (p.x + q.x) / 2, my = (p.y + q.y) / 2; let nx = -(q.y - p.y), ny = q.x - p.x; const l = Math.hypot(nx, ny) || 1; nx /= l; ny /= l;
    let sc = 0; for (const t of m.atoms.values()) sc += Math.sign((t.x - mx) * nx + (t.y - my) * ny); if (sc > 0) { nx = -nx; ny = -ny; }
    o += `<text class="bx-ez" x="${f1(mx + dx + nx * 15)}" y="${f1(my + dy + ny * 15 + 4)}">${esc(v)}</text>`; }
  return o ? `<g class="bx-stm">${o}</g>` : '';
}

function mount(el) {
  S.el = el; S.mol = new Chem.Mol();
  el.classList.add('bx-view');
  el.innerHTML = `<div class="bx">
  <nav class="bx-rail" aria-label="Drawing tools">${TOOLS.map(([k, n, key]) => `<button type="button" class="bx-t${k === S.tool ? ' on' : ''}" data-tool="${k}" title="${n} · ${key}" aria-label="${n}">${ICON[k]}<kbd>${key}</kbd></button>`).join('')}<span class="gap"></span><button type="button" class="bx-t" data-act="undo" title="Undo · Ctrl+Z" aria-label="Undo">${ICON.undo}</button><button type="button" class="bx-t" data-act="redo" title="Redo · Ctrl+Y" aria-label="Redo">${ICON.redo}</button><button type="button" class="bx-t" data-act="clear" title="Clear board" aria-label="Clear board">${ICON.clear}</button></nav>
  <section class="bx-main">
    <header class="bx-ctx"><div class="bx-opts"></div><div class="bx-right">
      <div class="bx-menu"><button type="button" class="bx-b" data-act="stmenu" aria-haspopup="true" title="View: animations, sidebar, stereo marks">${ICON.stereo}<span>View</span></button><div class="bx-pop bx-stpop" hidden>${stPopHTML()}</div></div>
      <div class="bx-menu"><button type="button" class="bx-b" data-act="tpl" aria-haspopup="true">${ICON.tpl}<span>Templates</span></button><div class="bx-pop bx-tplpop" hidden>${TEMPLATES.map(([g, items]) => `<div class="grp"><div class="gh">${esc(g)}</div>${items.map(([n, s]) => `<button type="button" data-smi="${esc(s)}">${esc(n)}</button>`).join('')}</div>`).join('')}</div></div>
      <label class="bx-smi" title="SMILES of the drawn structure. Paste one and press Enter."><span>SMILES</span><input spellcheck="false" autocomplete="off" placeholder="paste SMILES, Enter"></label>
      <button type="button" class="bx-b ic" data-act="copy" title="Copy SMILES · Ctrl+C" aria-label="Copy SMILES">${ICON.copy}</button>
      <button type="button" class="bx-b ic" data-act="help" title="Shortcuts · ?" aria-label="Shortcuts">${ICON.help}</button><button type="button" class="bx-b ic on" data-act="panel" title="Show / hide the mechanism sidebar · M" aria-label="Toggle mechanism sidebar">${ICON.side}</button></div></header>
    <div class="bx-can" tabindex="0" aria-label="Chemistry board">
      <svg class="bx-sheet" aria-hidden="true"><g class="cam"><g class="l-hl"></g><g class="l-mol board-mol"></g><g class="l-st"></g><g class="l-warn"></g><g class="l-sch board-mol"></g><g class="l-ghost"></g></g></svg>
      <div class="bx-pill" hidden><button type="button" data-pa="next">+ Next arrow</button><button type="button" data-pa="start" title="Copy this product onto the board as a new compound with its own arrows">Copy out</button></div>
      <div class="bx-empty" hidden><div class="in"><h3>Blank board</h3><p>Press <kbd>B</kbd> and click to draw bonds, <kbd>R</kbd> for rings, or start from a structure:</p><div class="bx-chipset">${['Propene', 'Toluene', 'Ethanal', 'Propanone', 'Cyclohexanol', 'Benzoic acid', 'Aniline', 'Propyne'].map(n => { const t = TEMPLATES.flatMap(g => g[1]).find(x => x[0] === n); return `<button type="button" class="chip" data-smi="${esc(t[1])}">${esc(n)}</button>`; }).join('')}</div></div></div>
      <div class="bx-load"><span></span>Loading the chemistry engine…</div>
      <div class="bx-help" hidden><div class="hd"><h3>Board shortcuts</h3><button type="button" class="bx-b ic" data-act="help" aria-label="Close">×</button></div><div class="cols">
        <div><h4>Tools</h4><p><kbd>V</kbd> select · <kbd>B</kbd> bond · <kbd>Z</kbd> chain · <kbd>R</kbd> ring · <kbd>T</kbd> label · <kbd>Q</kbd> charge · <kbd>E</kbd> erase · <kbd>A</kbd> arrow</p><p>Bond tool: <kbd>1</kbd><kbd>2</kbd><kbd>3</kbd> <kbd>W</kbd> <kbd>H</kbd> pick the type. Ring tool: <kbd>1</kbd> benzene, <kbd>3</kbd>–<kbd>8</kbd>.</p></div>
        <div><h4>Hover an atom</h4><p><kbd>C</kbd> <kbd>N</kbd> <kbd>O</kbd> <kbd>S</kbd> <kbd>F</kbd> <kbd>P</kbd> <kbd>H</kbd>, <kbd>L</kbd> Cl, <kbd>B</kbd> Br, <kbd>I</kbd> iodine · <kbd>+</kbd> <kbd>−</kbd> charge · <kbd>1</kbd>–<kbd>3</kbd> grow a bond · <kbd>↵</kbd> or double-click to type OH, NO2, COOH…</p><h4>Hover a bond</h4><p><kbd>1</kbd> <kbd>2</kbd> <kbd>3</kbd> order · <kbd>W</kbd> wedge · <kbd>H</kbd> hash</p></div>
        <div><h4>Board</h4><p><kbd>Ctrl</kbd><kbd>Z</kbd>/<kbd>Y</kbd> undo/redo · <kbd>Ctrl</kbd><kbd>A</kbd> select all · <kbd>Del</kbd> delete · arrows nudge · <kbd>Ctrl</kbd><kbd>Shift</kbd><kbd>K</kbd> tidy · <kbd>Ctrl</kbd><kbd>C</kbd> copy SMILES · <kbd>F</kbd> fit · wheel zooms · <kbd>Space</kbd>/right-drag pans</p><p>Mechanism: <kbd>←</kbd> <kbd>→</kbd> step through when nothing is selected.</p><p><kbd>M</kbd> sidebar · <kbd>N</kbd> next arrow from the active product · <kbd>Ctrl</kbd><kbd>D</kbd> duplicate · <kbd>*</kbd> stereocentres</p></div></div></div>
    </div>
    <footer class="bx-status"><div class="an"></div><div class="hint"></div><div class="zm"><button type="button" class="bx-b ic" data-act="zout" aria-label="Zoom out">${ICON.minus}</button><span data-z>100%</span><button type="button" class="bx-b ic" data-act="zin" aria-label="Zoom in">${ICON.plus}</button><button type="button" class="bx-b ic" data-act="fit" title="Fit · F" aria-label="Fit to content">${ICON.fit}</button></div></footer>
  </section>
  <aside class="bx-panel"><div class="bx-ph"><h2>Mechanism</h2><div class="bx-tabs" hidden></div></div><div class="bx-pbody"></div></aside></div>`;
  S.can = el.querySelector('.bx-can'); S.svg = el.querySelector('.bx-sheet'); S.cam = el.querySelector('.cam');
  S.Lg = { hl: el.querySelector('.l-hl'), mol: el.querySelector('.l-mol'), st: el.querySelector('.l-st'), warn: el.querySelector('.l-warn'), sch: el.querySelector('.l-sch'), ghost: el.querySelector('.l-ghost') };
  renderCtx(); if (getLS('ocr3-bxpanel', 'on') === 'off') { el.querySelector('.bx').classList.add('nopanel'); const pb = el.querySelector('[data-act=panel]'); if (pb) pb.classList.remove('on'); } applyCam(); renderPanel(); updateStatus();
  const pillEl = el.querySelector('.bx-pill'); pillEl.addEventListener('pointerenter', () => { S.overPill = true; }); pillEl.addEventListener('pointerleave', () => { S.overPill = false; S.pillHold = null; placeOverlays(); });
  const can = S.can;
  can.addEventListener('pointerdown', onDown); can.addEventListener('pointermove', onMove); can.addEventListener('pointerup', onUp);
  can.addEventListener('pointercancel', () => { S.g = null; S.ghost = null; drawHL(); });
  can.addEventListener('pointerleave', () => { if (!S.g) { S.hover = null; S.ghost = null; if (S.prodHover >= 0 && !document.querySelector('.bx-pill:hover')) { S.prodHover = -1; drawScheme(); placeOverlays(); } drawHL(); updateStatus(); } });
  can.addEventListener('wheel', onWheel, { passive: false }); can.addEventListener('contextmenu', e => e.preventDefault());
  can.addEventListener('paste', e => { const t = (e.clipboardData && e.clipboardData.getData('text') || '').trim(); if (t && !/\s/.test(t)) { e.preventDefault(); if (loadSmiles(t)) toast('Pasted structure · Ctrl+Z to undo'); } });
  document.addEventListener('keydown', e => { if (e.ctrlKey && e.shiftKey && e.key.toLowerCase() === 'k' && boardOn() && !(e.target.matches && e.target.matches('input,textarea'))) { e.preventDefault(); e.stopImmediatePropagation(); tidy(); } }, true);
  document.addEventListener('keydown', onKey); document.addEventListener('keyup', onKeyUp);
  window.addEventListener('blur', () => { S.space = false; });
  new ResizeObserver(() => { applyCam(); }).observe(can);
  el.addEventListener('click', e => {
    const t = e.target.closest('[data-tool],[data-act],[data-bt],[data-rt],[data-al],[data-cm],[data-rg],[data-rgpage],[data-nm],[data-nmi],[data-nmopen],[data-uselib],[data-smi],[data-try],[data-tab],[data-mv],[data-mvi],[data-pa],[data-eq],[data-fis],[data-fisdraw],[data-co],[data-sto],[data-vw]'); if (!t || !el.contains(t)) { if (!e.target.closest('.bx-menu')) el.querySelectorAll('.bx-pop').forEach(p => p.hidden = true); return; }
    const d = t.dataset;
    if (d.tool) { setTool(d.tool); return; }
    if (d.bt) { S.bondType = d.bt; renderCtx(); hoverAgain(); return; }
    if (d.rt) { S.ringType = d.rt; renderCtx(); hoverAgain(); return; }
    if (d.al) { S.atomLabel = d.al; renderCtx(); updateStatus(); return; }
    if (d.cm) { S.chargeMode = d.cm; renderCtx(); return; }
    if (d.smi != null) { el.querySelectorAll('.bx-pop').forEach(p => p.hidden = true); loadSmiles(d.smi, { add: e.shiftKey }); return; }
    if (d.sto) { setStOpt(d.sto); return; }
    if (d.vw) { setView(d.vw); return; }
    if (d.fis) { const [ii, k] = d.fis.split(':').map(Number); showFischer(ii, k); return; }
    if (d.fisdraw) { const [ii, k] = d.fisdraw.split(':').map(Number); fischerDraw(ii, k); return; }
    if (d.co) { const [ii, k] = d.co.split(':').map(Number); toggleCo(ii, k); return; }
    if (d.rg != null) { useReagentChip(d.rg); return; }
    if (d.rgpage) { A.go('reagents'); A.VIEWS.reagents.select(d.rgpage); return; }
    if (d.uselib != null) { const ar = S.chain[+d.uselib]; if (ar) { pushUndo(); ar.forceLib = true; scheduleSolve(0); toast('Using the library’s product · Ctrl+Z to go back'); } return; }
    if (d.nm) { NM.i += +d.nm; nmPaint(); return; } if (d.nmi != null) { NM.i = +d.nmi; nmPaint(); return; } if (d.nmopen) { A.openViewer(d.nmopen, [d.nmopen], { tab: (A.ART.ms[d.nmopen] || []).length ? 'mech' : 'reaction' }); return; }
    if (d.try != null) { const x = TRY[+d.try]; loadSmiles(x[1], { reagents: [x[2]], flight: true }); return; }
    if (d.eq) { const [ii, v] = d.eq.split(':'); setAmount(+ii, v === 'x' ? Infinity : +v); return; }
    if (d.tab != null) { S.active = +d.tab; drawScheme(); renderPanel(); return; }
    if (d.mvi != null) { mvStop(); mvGo(+d.mvi); return; }
    if (d.mv) { if (d.mv === 'prev') mvStep(-1); else if (d.mv === 'next') mvStep(1); else if (d.mv === 'play') mvPlay(); else if (d.mv === 'big') { MV.big = !MV.big; el.querySelector('.bx').classList.toggle('mv-big', MV.big); } return; }
    if (d.pa) { const i = +el.querySelector('.bx-pill').dataset.i; const ar = S.chain[i]; if (!ar) return; if (d.pa === 'next') addArrow('', null, { a: ar.id }); else { const p = ar.parts.find(q => q.main); pushUndo(); S.prodHover = -1; addCompound(p.m.clone()); toast('Copied out as a new compound · draw an arrow from it'); } return; }
    const a = d.act;
    if (a === 'undo') undo(); else if (a === 'redo') redo(); else if (a === 'clear') clearBoard(); else if (a === 'tidy') tidy(); else if (a === 'fliph') flip('h'); else if (a === 'flipv') flip('v'); else if (a === 'all') selectAll(); else if (a === 'mirror') mirrorSel(); else if (a === 'panel') togglePanel(); else if (a === 'clrarrows') clearArrows(); else if (a === 'dup') duplicateSel(); else if (a === 'fischer') fischerSel();
    else if (a === 'stmenu') { const p = el.querySelector('.bx-stpop'); el.querySelector('.bx-tplpop').hidden = true; p.hidden = !p.hidden; }
    else if (a === 'copy') copySmiles(); else if (a === 'help') toggleHelp(); else if (a === 'fit') fitView(); else if (a === 'zin') zoomAt(1.25, cw() / 2, ch() / 2); else if (a === 'zout') zoomAt(.8, cw() / 2, ch() / 2);
    else if (a === 'tpl') { const p = el.querySelector('.bx-tplpop'); el.querySelector('.bx-stpop').hidden = true; p.hidden = !p.hidden; }
  });
  const smi = el.querySelector('.bx-smi input');
  smi.addEventListener('keydown', e => { e.stopPropagation(); if (e.key === 'Enter') { e.preventDefault(); const v = smi.value.trim(); if (v && loadSmiles(v)) { smi.blur(); } } if (e.key === 'Escape') { smi.blur(); updateSmilesField(); } });
  smi.addEventListener('blur', () => setTimeout(updateSmilesField, 50));
  A.bindCards(el.querySelector('.bx-panel'), () => S.noteIds || []);
  A.rdkit().then(RD => { Chem.setRDKit(RD); S.rd = true; emptyState(); const p = S.pending; S.pending = null; if (p) p(); else loadSmiles('CC=C', { undo: false, anim: false }); }).catch(() => { toast('Chemistry engine could not load'); });
  drawAll();
}
function canOpen(r) { return !!(r && r.rxn && r.rxn.includes('>')); }
/* co-reactants written on the left of the arrow in the library (HBr in "CC=C.Br") become reagent text */
const CO_NAMES = { '[Br]': 'Br•', 'CC(=O)I': 'CH3COI', '[F-]': 'F-', 'CN(C)c1ccccc1': 'N,N-dimethylaniline', 'OP=O': 'H3PO2', '[Ag+]': 'Tollens', '[Cu+2]': 'Fehling', '[S]': 'sulphur', 'CCOC=O': 'HCOOC2H5', 'c1ccccc1': 'benzene', 'O=C1CCC(=O)O1': 'succinic anhydride', 'OS(=O)[O-]': 'NaHSO3', 'O=S([O-])O': 'NaHSO3', '[N+](=O)=O': 'NO2+', 'O=[N+]=O': 'NO2+', 'ClC#N': 'ClCN', 'N#CC#N': '(CN)2', 'IN=[N+]=[N-]': 'IN3', '[N-]=[N+]=NI': 'IN3', 'FS(F)(F)F': 'SF4', 'O=C1CCC(=O)N1Br': 'NBS', 'c1ccc(P(c2ccccc2)c2ccccc2)cc1': 'PPh3', 'O[I](=O)(=O)=O': 'HIO4', 'O=I(=O)(=O)O': 'HIO4', 'CC[S-]': 'C2H5S-', 'C[S-]': 'CH3S-', 'CCBr': 'C2H5Br', 'CCI': 'C2H5I', 'CCCl': 'C2H5Cl', 'CC(=O)[O-]': 'CH3COO-', '[S-2]': '(NH4)2S', 'CC(C)O[Al](OC(C)C)OC(C)C': 'Al(OiPr)3', 'Cl[O-]': 'bleaching powder', 'Cl[Sb](Cl)(Cl)(Cl)Cl': 'SbCl5', 'CC(C)=O': 'acetone', 'CC=O': 'CH3CHO', 'C=O': 'HCHO', 'O=Cc1ccccc1': 'C6H5CHO', 'CCN': 'C2H5NH2', 'CN': 'CH3NH2', 'Nc1ccccc1': 'aniline', 'B': 'BH3', 'NNc1ccccc1': 'PhNHNH2', 'O=NCl': 'NOCl', 'FF': 'F2', 'CC(=O)OO': 'CH3CO3H', '[Mg]': 'Mg', 'NCl': 'NH2Cl', 'C1CO1': 'ethylene oxide', '[CH3][Cu-][CH3]': '(CH3)2CuLi', 'C[Cu-]C': '(CH3)2CuLi', '[Cl][Fe]([Cl])[Cl]': 'FeCl3', 'Cl[Fe](Cl)Cl': 'FeCl3', 'C[CH2][Cd][CH2]C': '(C2H5)2Cd', 'CC[Cd]CC': '(C2H5)2Cd', 'O=CC(Cl)(Cl)Cl': 'chloral', 'Oc1ccccc1': 'phenol', 'N#[N+]c1ccccc1': 'C6H5N2+', '[H+]': 'H+', Br: 'HBr', Cl: 'HCl', I: 'HI', BrBr: 'Br2', ClCl: 'Cl2', II: 'I2', O: 'H2O', '[OH-]': 'OH-', N: 'NH3', 'C#N': 'HCN', '[C-]#N': 'KCN', 'O=S(=O)(O)O': 'H2SO4', 'O=[N+]([O-])O': 'HNO3', OO: 'H2O2', CI: 'CH3I', CBr: 'CH3Br', CCl: 'CH3Cl', 'C[Mg]Br': 'CH3MgBr', 'C[Mg]I': 'CH3MgI', 'CC[Mg]Br': 'C2H5MgBr', 'O=C=O': 'CO2', 'C=O': 'HCHO', 'CC(=O)Cl': 'CH3COCl', 'ClC(Cl)Cl': 'CHCl3', 'O=O': 'O2', 'NN': 'NH2NH2', NO: 'NH2OH', '[Na]': 'Na', '[BH4-]': 'NaBH4', 'CC[O-]': 'NaOEt', 'C[O-]': 'NaOMe', 'CC(C)(C)[O-]': 't-BuOK', 'CCO': 'EtOH', 'CO': 'MeOH', 'CC(=O)OC(C)=O': '(CH3CO)2O', 'O=C(Cl)c1ccccc1': 'PhCOCl', 'Cc1ccc(S(=O)(=O)Cl)cc1': 'TsCl', '[Cl-]': 'Cl-', '[Br-]': 'Br-', '[I-]': 'I-', 'ClP(Cl)(Cl)(Cl)Cl': 'PCl5', 'O=S(Cl)Cl': 'SOCl2', 'ClP(Cl)Cl': 'PCl3', 'O=Cc1ccccc1': 'PhCHO', '[NH2-]': 'NaNH2', 'O=N[O-]': 'NaNO2', 'BrN1C(=O)CCC1=O': 'NBS', 'O=C1CCC(=O)N1Br': 'NBS' };
/* the drawn structure: the heaviest reactant that is not itself a known reagent (PhNHNH₂, NBS, peracid...) */
const isCoName = x => !!(CO_NAMES[canonS(x)] || CO_NAMES[x]);
function pickMain(list) { const s = list.slice().sort((a, b) => heavy(b) - heavy(a)); const hasC = x => /C|c/.test(x.replace(/Cl|Ca|Cd|Cu|Cr|Co|Cs/g, '')); const i = s.findIndex(x => hasC(x) && !isCoName(x)); const k = i < 0 ? 0 : i; return { main: s[k] || '', rest: s.filter((_, j) => j !== k) }; }
function boardSetup(r) {
  const sides = (r.rxn || '').split('>'); const { main, rest } = pickMain(sides[0].split('.').filter(Boolean));
  const mid = !(r.rl || []).filter(Boolean).length && sides.length > 2 ? sides[1].split('.').filter(Boolean) : [];
  const co = rest.concat(mid).map(x => CO_NAMES[canonS(x)] || CO_NAMES[x]).filter(Boolean);
  const lab = [...(r.rl || []), ...(r.cl || [])].filter(Boolean); const have = speciesOf(lab.join(', '));
  const co2 = co.filter(nm => { const sp = [...speciesOf(nm)]; return !sp.length || !sp.every(k => have.has(k)); }); // "H⁺" already in "H⁺ (cat.)"
  let rg = [...new Set([...co2, ...lab])].join(', ');
  if (!rg && new Set(sides[0].split('.').filter(Boolean).map(x => canonS(x) || x)).size === 1) rg = /\[(C|CH|CH2|c)\+\]/.test(main) ? 'rearrange' : /\[(CH2|CH|CH3)\]/.test(main) ? 'couple' : 'tautomerise';
  return { main, rg };
}
function openReaction(r) {
  const go0 = () => { const { main, rg } = boardSetup(r); loadSmiles(main, { reagents: [rg || ''], flight: true }); if (!rg) setTimeout(() => openReagentEditor(0), 300); };
  if (S.rd) go0(); else S.pending = go0;
}
function openReactionOld(r) {
  const smi = (r.rxn || '').split('>')[0]; const rg = (r.rl || []).filter(Boolean).join(', ') || (r.cl || []).filter(Boolean).join(', ');
  const go = () => { loadSmiles(smi, { reagents: [rg || ''], flight: true }); if (!rg) setTimeout(() => openReagentEditor(0), 300); };
  if (S.rd) go(); else S.pending = go;
}
A.VIEWS.board = { mount, show: () => { if (S.el) { applyCam(); drawAll(); } }, hide: () => { closeEditors(true); mvStop(); } };
window.Board = { canOpen, openReaction, _S: S, _scr: (x, y) => { const r = S.can.getBoundingClientRect(); const [a, b] = toScreen(x, y); return [a + r.left, b + r.top]; } };
})();
