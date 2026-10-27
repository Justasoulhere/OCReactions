/* chem.js: small chemistry core for the sketch board.
   Graph model with implicit hydrogens by valence, SMILES writer (exact bracket atoms, Kekule),
   RDKit.js bridge (matching, canonical SMILES, coordinates) and a graph-rewrite reaction engine. */
(function (root) {
'use strict';
const SYM = [null, 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba'];
const EXTRA = { 78: 'Pt', 79: 'Au', 80: 'Hg', 82: 'Pb', 58: 'Ce' };
const Z = {}; SYM.forEach((s, i) => { if (s) Z[s] = i }); for (const k in EXTRA) Z[EXTRA[k]] = +k;
const symOf = z => EXTRA[z] || SYM[z] || '*';
const BASEVAL = { H: [1], B: [3], C: [4], N: [3], O: [2], F: [1], Si: [4], P: [3, 5], S: [2, 4, 6], Cl: [1], Br: [1], I: [1, 3, 5, 7], Se: [2, 4, 6] };
const ISO = { C: { 1: [3], '-1': [3] }, N: { 1: [4], '-1': [2] }, O: { 1: [3], '-1': [1] }, S: { 1: [3, 5], '-1': [1] }, P: { 1: [4], '-1': [2] }, B: { '-1': [4], 1: [2] },
  F: { 1: [0, 1, 2], '-1': [0] }, Cl: { 1: [0, 1, 2], '-1': [0] }, Br: { 1: [0, 1, 2], '-1': [0] }, I: { 1: [0, 1, 2], '-1': [0] }, H: { 1: [0], '-1': [0] }, Se: { 1: [3], '-1': [1] } };
function valList(el, chg) {
  if (!chg) return BASEVAL[el] || null;
  const t = ISO[el]; return t && t[chg] ? t[chg] : null;
}
const key = (a, b) => a < b ? a + '|' + b : b + '|' + a;

class Mol {
  constructor() { this.atoms = new Map(); this.bonds = new Map(); this.adj = new Map(); this.nid = 1; this.nbid = 1; }
  addAtom(el, x = 0, y = 0, p = {}) {
    const id = p.id != null ? p.id : this.nid++; if (id >= this.nid) this.nid = id + 1;
    const a = { id, el, x, y, chg: p.chg || 0, rad: p.rad || 0, iso: p.iso || 0, hFix: p.hFix != null ? p.hFix : null, tag: p.tag || null, isNew: !!p.isNew };
    this.atoms.set(id, a); this.adj.set(id, new Set()); return id;
  }
  addBond(a, b, o = 1, p = {}) {
    if (a === b || !this.atoms.has(a) || !this.atoms.has(b)) return null;
    const ex = this.bondBetween(a, b); if (ex) { ex.o = o; return ex.id; }
    const id = p.id != null ? p.id : this.nbid++; if (id >= this.nbid) this.nbid = id + 1;
    this.bonds.set(id, { id, a, b, o, st: p.st || 0 }); this.adj.get(a).add(id); this.adj.get(b).add(id); return id;
  }
  removeBond(id) { const b = this.bonds.get(id); if (!b) return; this.adj.get(b.a).delete(id); this.adj.get(b.b).delete(id); this.bonds.delete(id); }
  removeAtom(id) { if (!this.atoms.has(id)) return; for (const bid of [...this.adj.get(id)]) this.removeBond(bid); this.adj.delete(id); this.atoms.delete(id); }
  bondBetween(a, b) { const s = this.adj.get(a); if (!s) return null; for (const bid of s) { const bd = this.bonds.get(bid); if (bd.a === b || bd.b === b) return bd; } return null; }
  bondsOf(a) { return [...(this.adj.get(a) || [])].map(id => this.bonds.get(id)); }
  nbrs(a) { return this.bondsOf(a).map(b => b.a === a ? b.b : b.a); }
  other(b, a) { return b.a === a ? b.b : b.a; }
  degree(a) { return (this.adj.get(a) || new Set()).size; }
  heavyDegree(a) { return this.nbrs(a).filter(n => this.atoms.get(n).el !== 'H').length; }
  bondSum(a) { let v = 0; for (const b of this.bondsOf(a)) v += b.o; return v; }
  hImp(a) {
    const at = typeof a === 'object' ? a : this.atoms.get(a); if (at.hFix != null) return at.hFix;
    const vl = valList(at.el, at.chg); if (!vl) return 0;
    const v = this.bondSum(at.id) + (at.rad || 0);
    for (const t of vl) if (t >= v) return t - v; return 0;
  }
  hTot(a) { return this.hImp(a) + this.nbrs(a).filter(n => this.atoms.get(n).el === 'H').length; }
  valenceError(a) {
    const at = this.atoms.get(a); const vl = valList(at.el, at.chg); if (!vl || at.hFix != null) return false;
    return this.bondSum(a) + (at.rad || 0) > Math.max(...vl);
  }
  clone() {
    const m = new Mol(); m.nid = this.nid; m.nbid = this.nbid;
    for (const a of this.atoms.values()) { m.atoms.set(a.id, { ...a }); m.adj.set(a.id, new Set()); }
    for (const b of this.bonds.values()) { m.bonds.set(b.id, { ...b }); m.adj.get(b.a).add(b.id); m.adj.get(b.b).add(b.id); }
    return m;
  }
  frags() {
    const seen = new Set(), out = [];
    for (const id of this.atoms.keys()) {
      if (seen.has(id)) continue; const comp = [], st = [id]; seen.add(id);
      while (st.length) { const x = st.pop(); comp.push(x); for (const n of this.nbrs(x)) if (!seen.has(n)) { seen.add(n); st.push(n); } }
      out.push(comp);
    }
    return out;
  }
  sub(ids) {
    const S = new Set(ids), m = new Mol(); m.nid = this.nid; m.nbid = this.nbid;
    for (const id of S) { m.atoms.set(id, { ...this.atoms.get(id) }); m.adj.set(id, new Set()); }
    for (const b of this.bonds.values()) if (S.has(b.a) && S.has(b.b)) { m.bonds.set(b.id, { ...b }); m.adj.get(b.a).add(b.id); m.adj.get(b.b).add(b.id); }
    return m;
  }
  merge(o, dx = 0, dy = 0, tag = null) {
    const idm = new Map();
    for (const a of o.atoms.values()) { const id = this.addAtom(a.el, a.x + dx, a.y + dy, { ...a, id: undefined, tag: tag || a.tag }); idm.set(a.id, id); }
    for (const b of o.bonds.values()) this.addBond(idm.get(b.a), idm.get(b.b), b.o, { st: b.st });
    return idm;
  }
  bbox(ids) {
    let x0 = 1e9, y0 = 1e9, x1 = -1e9, y1 = -1e9;
    for (const id of ids || this.atoms.keys()) { const a = this.atoms.get(id); x0 = Math.min(x0, a.x); y0 = Math.min(y0, a.y); x1 = Math.max(x1, a.x); y1 = Math.max(y1, a.y); }
    return { x0, y0, x1, y1 };
  }
  toJSON() { return { atoms: [...this.atoms.values()].map(a => ({ id: a.id, el: a.el, x: +a.x.toFixed(2), y: +a.y.toFixed(2), chg: a.chg || undefined, rad: a.rad || undefined, hFix: a.hFix == null ? undefined : a.hFix })), bonds: [...this.bonds.values()].map(b => [b.a, b.b, b.o]) }; }
  static fromJSON(j) { const m = new Mol(); for (const a of j.atoms) m.addAtom(a.el, a.x, a.y, a); for (const [a, b, o] of j.bonds) m.addBond(a, b, o); return m; }
}

/* ---------- rings & aromaticity (Kekule 6-rings with alternating bonds) ---------- */
function rings(m, maxSize = 8) {
  // smallest cycles through each bond (good enough for drawing and aromaticity)
  const out = [], seen = new Set();
  for (const b of m.bonds.values()) {
    // BFS from b.a to b.b without using b
    const prev = new Map([[b.a, null]]), q = [b.a]; let found = false;
    while (q.length && !found) {
      const x = q.shift();
      for (const bd of m.bondsOf(x)) {
        if (bd.id === b.id) continue; const y = m.other(bd, x);
        if (prev.has(y)) continue; prev.set(y, x);
        if (y === b.b) { found = true; break; }
        q.push(y);
      }
    }
    if (!found) continue;
    const path = []; let x = b.b; while (x != null) { path.push(x); x = prev.get(x); }
    if (path.length > maxSize) continue;
    const k = [...path].sort((p, q) => p - q).join(',');
    if (!seen.has(k)) { seen.add(k); out.push(path); }
  }
  return out;
}
function aromaticRings(m) {
  const res = [];
  for (const r of rings(m, 6)) {
    if (r.length !== 6) continue;
    let ok = true, dbl = 0;
    for (let i = 0; i < 6; i++) {
      const a = m.atoms.get(r[i]); if (!['C', 'N'].includes(a.el)) { ok = false; break; }
      const b = m.bondBetween(r[i], r[(i + 1) % 6]); if (!b) { ok = false; break; }
      if (b.o === 2) dbl++;
    }
    if (!ok) continue;
    // each ring atom must carry exactly one double bond (in or out of ring makes it non-aromatic if exocyclic)
    let alt = true;
    for (const id of r) { const d = m.bondsOf(id).filter(b => b.o === 2); if (d.length !== 1) { alt = false; break; } const o = m.other(d[0], id); if (!r.includes(o)) { alt = false; break; } }
    if (alt) res.push(r);
  }
  return res;
}

/* ---------- SMILES writer: every atom bracketed with exact H count ---------- */
function atomToken(m, a, opt) {
  const h = (opt.implicitH === false) ? 0 : m.hImp(a);
  let s = '[' + (a.iso ? a.iso : '') + a.el;
  if (h) s += 'H' + (h > 1 ? h : '');
  if (a.chg > 0) s += '+' + (a.chg > 1 ? a.chg : ''); else if (a.chg < 0) s += '-' + (a.chg < -1 ? -a.chg : '');
  if (opt.maps && opt.maps.get(a.id)) s += ':' + opt.maps.get(a.id);
  return s + ']';
}
const BSYM = o => o === 2 ? '=' : o === 3 ? '#' : '';
function toSmiles(m, opt = {}) {
  const ids = opt.ids ? opt.ids.filter(id => m.atoms.has(id)) : [...m.atoms.keys()];
  const S = new Set(ids), visited = new Set(), order = [], parts = [];
  for (const start of ids) {
    if (visited.has(start)) continue;
    // pass 1: DFS tree + ring closures
    const children = new Map(), closures = new Map(), usedB = new Set(), seq = [];
    const dfs1 = (x) => {
      visited.add(x); seq.push(x); children.set(x, []); closures.set(x, closures.get(x) || []);
      for (const b of m.bondsOf(x)) {
        const y = m.other(b, x); if (!S.has(y) || usedB.has(b.id)) continue;
        usedB.add(b.id);
        if (!visited.has(y)) { children.get(x).push([y, b]); dfs1(y); }
        else { closures.get(x).push(b); if (!closures.has(y)) closures.set(y, []); closures.get(y).push(b); }
      }
    };
    dfs1(start);
    // pass 2: emit
    const ringNum = new Map(); const free = [];
    let next = 1;
    const takeNum = () => free.length ? free.sort((p, q) => p - q).shift() : next++;
    const emit = (x) => {
      order.push(x); let s = atomToken(m, m.atoms.get(x), opt);
      for (const b of closures.get(x)) {
        if (ringNum.has(b.id)) { const n = ringNum.get(b.id); s += BSYM(b.o) + (n < 10 ? n : '%' + n); ringNum.delete(b.id); free.push(n); }
        else { const n = takeNum(); ringNum.set(b.id, n); s += BSYM(b.o) + (n < 10 ? n : '%' + n); }
      }
      const ch = children.get(x);
      ch.forEach(([y, b], k) => { const sub = BSYM(b.o) + emit(y); s += (k < ch.length - 1) ? '(' + sub + ')' : sub; });
      return s;
    };
    parts.push(emit(start));
  }
  return { smi: parts.join('.'), order };
}

/* ---------- SMILES reader (for templates and co-reactants) ---------- */
const ORGANIC2 = ['Cl', 'Br'];
function parseSmiles(s) {
  const atoms = [], bonds = []; let prev = null, bo = null; const stack = [], ringsO = {};
  let i = 0;
  while (i < s.length) {
    const c = s[i];
    if (c === '(') { stack.push(prev); i++; continue; }
    if (c === ')') { prev = stack.pop(); i++; continue; }
    if (c === '.') { prev = null; i++; continue; }
    if (c === '-') { bo = 1; i++; continue; } if (c === '=') { bo = 2; i++; continue; } if (c === '#') { bo = 3; i++; continue; }
    if (c === ':') { bo = 1.5; i++; continue; } if (c === '/' || c === '\\') { i++; continue; }
    if (/[0-9]/.test(c) || c === '%') {
      let num; if (c === '%') { num = s.slice(i + 1, i + 3); i += 3; } else { num = c; i++; }
      if (ringsO[num]) { const r = ringsO[num]; bonds.push({ a: r.atom, b: prev, o: bo || r.o || 1 }); delete ringsO[num]; }
      else ringsO[num] = { atom: prev, o: bo };
      bo = null; continue;
    }
    const a = { el: null, chg: 0, map: 0, rad: 0, arom: false, h: null, iso: 0 };
    if (c === '[') {
      const j = s.indexOf(']', i); let t = s.slice(i + 1, j); i = j + 1;
      let m1 = t.match(/^(\d+)/); if (m1) { a.iso = +m1[1]; t = t.slice(m1[1].length); }
      m1 = t.match(/^([A-Z][a-z]?|[a-z][a-z]?|\*)/); let el = m1 ? m1[1] : '*';
      if (el.length === 2 && !Z[el[0].toUpperCase() + el[1]] && Z[el[0].toUpperCase()]) el = el[0];
      t = t.slice(el.length);
      if (el !== '*' && el[0] === el[0].toLowerCase()) { a.arom = true; el = el[0].toUpperCase() + el.slice(1); }
      a.el = el;
      m1 = t.match(/^H(\d*)/); if (m1) { a.h = m1[1] ? +m1[1] : 1; t = t.slice(m1[0].length); }
      m1 = t.match(/^(\++|-+)(\d*)/); if (m1) { const sg = m1[1][0] === '+' ? 1 : -1; a.chg = sg * (m1[2] ? +m1[2] : m1[1].length); t = t.slice(m1[0].length); }
      m1 = t.match(/^(\^+)/); if (m1) { a.rad = m1[1].length; t = t.slice(m1[0].length); }
      m1 = t.match(/:(\d+)$/); if (m1) a.map = +m1[1];
    } else {
      if (ORGANIC2.includes(s.slice(i, i + 2))) { a.el = s.slice(i, i + 2); i += 2; }
      else if (/[BCNOPSFI]/.test(c)) { a.el = c; i++; }
      else if (/[cnospb]/.test(c)) { a.el = c.toUpperCase(); a.arom = true; i++; }
      else if (c === '*') { a.el = '*'; i++; }
      else { i++; continue; }
    }
    atoms.push(a); const idx = atoms.length - 1;
    if (prev !== null) bonds.push({ a: prev, b: idx, o: bo || (a.arom && atoms[prev].arom ? 1.5 : 1) });
    prev = idx; bo = null;
  }
  return { atoms, bonds };
}
function smartsMaps(sm) {
  const maps = []; let i = 0;
  while (i < sm.length) {
    const c = sm[i];
    if (c === '[') { let d = 1, j = i + 1; while (j < sm.length && d) { if (sm[j] === '[') d++; else if (sm[j] === ']') d--; j++; } const t = sm.slice(i, j); const mm = t.match(/:(\d+)\]$/); maps.push(mm ? +mm[1] : 0); i = j; continue; }
    if (c === '%') { i += 3; continue; }
    if (sm.startsWith('Cl', i) || sm.startsWith('Br', i)) { maps.push(0); i += 2; continue; }
    if (/[BCNOPSFIcnospbAa*]/.test(c)) { maps.push(0); i++; continue; }
    i++;
  }
  return maps;
}

/* ---------- RDKit bridge ---------- */
let RD = null; const QC = new Map();
function setRDKit(r) { RD = r; }
function qmol(sm) { let q = QC.get(sm); if (!q) { q = RD.get_qmol(sm); QC.set(sm, q); } return q; }
function rdMol(m, opt = {}) {
  const { smi, order } = toSmiles(m, opt);
  const mol = RD.get_mol(smi, JSON.stringify({ removeHs: false, setAromaticity: opt.kekule ? false : true }));
  return { mol, order, smi };
}
function matchAll(m, smarts, opt = {}) {
  const { mol, order } = rdMol(m, opt); if (!mol) return [];
  let res; try { res = JSON.parse(mol.get_substruct_matches(qmol(smarts), JSON.stringify({ uniquify: opt.uniquify !== false ? true : false, maxMatches: 400 }))); } catch (e) { res = []; }
  mol.delete(); if (!Array.isArray(res)) res = [];
  return res.map(r => r.atoms.map(k => order[k]));
}
function canon(m, ids) {
  const { smi } = toSmiles(m, { ids });
  const mol = RD.get_mol(smi); if (!mol) return null;
  const c = mol.get_smiles(); mol.delete(); return c;
}
function isValid(m) {
  const { smi } = toSmiles(m); const mol = RD.get_mol(smi); if (!mol) return false; mol.delete(); return true;
}
function fromSmiles(smi, scale = 1) {
  const mol = RD.get_mol(smi); if (!mol) return null;
  const k = RD.get_mol(mol.get_smiles(JSON.stringify({ kekuleSmiles: true })), JSON.stringify({ setAromaticity: false }));
  mol.delete(); if (!k) return null;
  k.set_new_coords(); const j = JSON.parse(k.get_json()).molecules[0]; const xy = k.get_coords(); k.delete();
  return fromJSONMol(j, xy, scale);
}
/* coordinates + wedges from RDKit's molblock, exact H counts from its JSON */
function fromSmilesMB(smi, scale = 40) {
  const mol = RD.get_mol(smi); if (!mol) return null;
  try { mol.set_new_coords(); } catch (e) { }
  const mb = mol.get_molblock(); let js = null; try { js = JSON.parse(mol.get_json()).molecules[0]; } catch (e) { }
  mol.delete();
  const L = mb.split('\n'); const na = +L[3].slice(0, 3), nb = +L[3].slice(3, 6);
  const m = new Mol(); const ids = []; const k = scale / 1.5;
  for (let i = 0; i < na; i++) { const t = L[4 + i]; const x = +t.slice(0, 10), y = +t.slice(10, 20); const e = t.slice(31, 34).trim(); ids.push(m.addAtom(e, x * k, -y * k)); }
  for (let i = 0; i < nb; i++) { const t = L[4 + na + i]; const a = +t.slice(0, 3), b = +t.slice(3, 6), o = +t.slice(6, 9), st = +t.slice(9, 12); m.addBond(ids[a - 1], ids[b - 1], o === 4 ? 1.5 : o, { st: o === 1 ? st : 0 }); }
  for (let i = 4 + na + nb; i < L.length; i++) {
    const t = L[i]; const tag = t.slice(0, 6);
    if (tag === 'M  CHG' || tag === 'M  RAD' || tag === 'M  ISO') { const n = +t.slice(6, 9); for (let j = 0; j < n; j++) { const idx = +t.slice(10 + 8 * j, 13 + 8 * j), v = +t.slice(14 + 8 * j, 17 + 8 * j); const a = m.atoms.get(ids[idx - 1]); if (tag === 'M  CHG') a.chg = v; else if (tag === 'M  RAD') a.rad = v === 2 ? 1 : 2; else a.iso = v; } }
  }
  if (js) js.atoms.forEach((a, i) => { const at = m.atoms.get(ids[i]); if (!at) return; const want = a.impHs || 0; if (m.hImp(at) !== want) at.hFix = want; });
  return m;
}
function fromJSONMol(j, xy, scale) {
  const m = new Mol(); const ids = [];
  j.atoms.forEach((a, i) => { const el = symOf(a.z == null ? 6 : a.z); ids.push(m.addAtom(el, xy[i][0] * scale, -xy[i][1] * scale, { chg: a.chg || 0, rad: a.nRad || 0, iso: a.isotope || 0 })); });
  j.bonds.forEach(b => m.addBond(ids[b.atoms[0]], ids[b.atoms[1]], b.bo == null ? 1 : b.bo));
  return m;
}
/* fresh 2D coordinates for some atoms, optionally aligned to the unchanged part */
function relayout(m, scale, keepIds) {
  const { smi, order } = toSmiles(m);
  const mol = RD.get_mol(smi, JSON.stringify({ removeHs: false, setAromaticity: false })); if (!mol) return false;
  let aligned = false;
  if (keepIds && keepIds.length >= 2) {
    const sub = m.sub(keepIds.filter(id => m.atoms.has(id)));
    const { smi: s2, order: o2 } = toSmiles(sub);
    const t = RD.get_mol(s2, JSON.stringify({ removeHs: false, setAromaticity: false }));
    if (t) {
      const mb = molblock(sub, o2, scale); const t2 = RD.get_mol(mb, JSON.stringify({ removeHs: false, setAromaticity: false }));
      if (t2) { try { const r = JSON.parse(mol.generate_aligned_coords(t2, JSON.stringify({ useCoordGen: false, acceptFailure: true, allowRGroups: false }))); aligned = r && r.atoms && r.atoms.length > 0; } catch (e) { } t2.delete(); }
      t.delete();
    }
  }
  if (!aligned) mol.set_new_coords();
  const xy = mol.get_coords(); mol.delete();
  order.forEach((id, k) => { const a = m.atoms.get(id); a.x = xy[k][0] * scale; a.y = -xy[k][1] * scale; });
  return aligned;
}
function molblock(m, order, scale) {
  const idx = new Map(order.map((id, k) => [id, k + 1]));
  const L = ['', '  board', '']; const bonds = [...m.bonds.values()].filter(b => idx.has(b.a) && idx.has(b.b));
  L.push(String(order.length).padStart(3) + String(bonds.length).padStart(3) + '  0  0  0  0  0  0  0  0999 V2000');
  for (const id of order) { const a = m.atoms.get(id); L.push((a.x / scale).toFixed(4).padStart(10) + (-a.y / scale).toFixed(4).padStart(10) + '    0.0000 ' + a.el.padEnd(3) + ' 0  0  0  0  0  0  0  0  0  0  0  0'); }
  for (const b of bonds) L.push(String(idx.get(b.a)).padStart(3) + String(idx.get(b.b)).padStart(3) + String(b.o === 1.5 ? 4 : b.o).padStart(3) + '  0');
  const ch = order.map(id => m.atoms.get(id)).map((a, k) => [k + 1, a.chg]).filter(x => x[1]);
  if (ch.length) L.push('M  CHG' + String(ch.length).padStart(3) + ch.map(([k, c]) => String(k).padStart(4) + String(c).padStart(4)).join(''));
  L.push('M  END'); return L.join('\n');
}

/* ---------- geometry helpers for placing new atoms ---------- */
function placeNew(m, BL) {
  const pending = [...m.atoms.values()].filter(a => a.isNew).map(a => a.id);
  let guard = 0;
  while (pending.length && guard++ < 200) {
    let progressed = false;
    for (let k = 0; k < pending.length; k++) {
      const id = pending[k]; const a = m.atoms.get(id);
      const anchors = m.nbrs(id).filter(n => !m.atoms.get(n).isNew);
      if (!anchors.length) continue;
      const an = m.atoms.get(anchors[0]);
      const dirs = m.nbrs(an.id).filter(n => n !== id && !m.atoms.get(n).isNew).map(n => Math.atan2(m.atoms.get(n).y - an.y, m.atoms.get(n).x - an.x));
      const ang = freeAngle(dirs, an, m);
      a.x = an.x + Math.cos(ang) * BL; a.y = an.y + Math.sin(ang) * BL; a.isNew = false;
      pending.splice(k, 1); k--; progressed = true;
    }
    if (!progressed) { const a = m.atoms.get(pending.shift()); a.isNew = false; }
  }
}
function freeAngle(dirs, an, m) {
  if (!dirs.length) return -Math.PI / 6;
  if (dirs.length === 1) {
    // zig-zag: 120 degrees from the existing bond, choose side away from its own neighbour
    const d = dirs[0]; const c1 = d + 2 * Math.PI / 3, c2 = d - 2 * Math.PI / 3;
    return score(c1) >= score(c2) ? c1 : c2;
  }
  const s = [...dirs].map(x => (x + 2 * Math.PI) % (2 * Math.PI)).sort((p, q) => p - q);
  let best = 0, bg = -1;
  for (let i = 0; i < s.length; i++) { const a0 = s[i], a1 = i + 1 < s.length ? s[i + 1] : s[0] + 2 * Math.PI; if (a1 - a0 > bg) { bg = a1 - a0; best = a0 + (a1 - a0) / 2; } }
  return best;
  function score(ang) { // distance to nearby atoms
    const px = an.x + Math.cos(ang) * 30, py = an.y + Math.sin(ang) * 30; let dmin = 1e9;
    for (const o of m.atoms.values()) { if (o.id === an.id) continue; dmin = Math.min(dmin, Math.hypot(o.x - px, o.y - py)); }
    return dmin;
  }
}

root.Chem = { Mol, Z, symOf, valList, key, rings, aromaticRings, toSmiles, parseSmiles, smartsMaps, setRDKit, rdMol, matchAll, canon, isValid, fromSmiles, fromSmilesMB, fromJSONMol, relayout, molblock, placeNew, freeAngle, get RD() { return RD; } };
})(typeof window !== 'undefined' ? window : globalThis);
