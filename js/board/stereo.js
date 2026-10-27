/* stereo.js: stereochemistry for the Board (Oct 16, run 23).
   Everything goes through RDKit with a stereo-aware molblock (wedges + 2D coordinates), so R/S, E/Z and
   stereocentre detection follow the drawing exactly. Atom ids of our Mol are kept everywhere. */
(function (root) {
'use strict';
const C = () => root.Chem; const RD = () => root.Chem.RD;
const BL = 40;

/* molblock that keeps wedges (st 1 = wedge, 6 = hash, drawn from bond.a), charges, radicals and isotopes */
function mb(m, opt = {}) {
  const order = opt.order || [...m.atoms.keys()]; const idx = new Map(order.map((id, k) => [id, k + 1]));
  const bonds = [...m.bonds.values()].filter(b => idx.has(b.a) && idx.has(b.b));
  const L = ['', '  stereo', '', String(order.length).padStart(3) + String(bonds.length).padStart(3) + '  0  0  0  0  0  0  0  0999 V2000'];
  for (const id of order) { const a = m.atoms.get(id); L.push((a.x / BL).toFixed(4).padStart(10) + (-a.y / BL).toFixed(4).padStart(10) + '    0.0000 ' + a.el.padEnd(3) + ' 0  0  0  0  0  0  0  0  0  0  0  0'); }
  for (const b of bonds) { const st = opt.noStereo ? 0 : (b.o === 1 && (b.st === 1 || b.st === 6) ? b.st : 0); L.push(String(idx.get(b.a)).padStart(3) + String(idx.get(b.b)).padStart(3) + String(b.o === 1.5 ? 4 : Math.round(b.o)).padStart(3) + String(st).padStart(3)); }
  const ch = order.map((id, k) => [k + 1, m.atoms.get(id)]).filter(x => x[1].chg); if (ch.length) for (let i = 0; i < ch.length; i += 8) { const p = ch.slice(i, i + 8); L.push('M  CHG' + String(p.length).padStart(3) + p.map(([k, a]) => String(k).padStart(4) + String(a.chg).padStart(4)).join('')); }
  const rd = order.map((id, k) => [k + 1, m.atoms.get(id)]).filter(x => x[1].rad); if (rd.length) L.push('M  RAD' + String(rd.length).padStart(3) + rd.map(([k, a]) => String(k).padStart(4) + String(a.rad === 2 ? 1 : 2).padStart(4)).join(''));
  const iso = order.map((id, k) => [k + 1, m.atoms.get(id)]).filter(x => x[1].iso); if (iso.length) L.push('M  ISO' + String(iso.length).padStart(3) + iso.map(([k, a]) => String(k).padStart(4) + String(a.iso).padStart(4)).join(''));
  L.push('M  END'); return { text: L.join('\n'), order };
}
function rdOf(m, opt) { const { text, order } = mb(m, opt); let mol = null; try { mol = RD().get_mol(text, JSON.stringify({ removeHs: false })); } catch (e) { mol = null; } if (mol && !mol.is_valid()) { mol.delete(); mol = null; } return { mol, order }; }
const frag1 = m => m.atoms.size > 0;

/* R/S and E/Z labels from the drawing: { atoms: Map(id → 'R'|'S'|'r'|'s'), bonds: Map(bondId → 'E'|'Z'), all: Set(potential centres), allB: Set(stereogenic C=C) } */
function cip(m) {
  const out = { atoms: new Map(), bonds: new Map(), all: new Set(), allB: new Set() }; if (!frag1(m)) return out;
  const { mol, order } = rdOf(m); if (!mol) return out;
  try { const t = JSON.parse(mol.get_stereo_tags());
    for (const [i, c] of t.CIP_atoms || []) { const v = String(c).replace(/[()]/g, ''); out.all.add(order[i]); if (v !== '?') out.atoms.set(order[i], v); }
    for (const [i, j, c] of t.CIP_bonds || []) { const b = m.bondBetween(order[i], order[j]); if (!b) continue; const v = String(c).replace(/[()]/g, ''); out.allB.add(b.id); if (v !== '?') out.bonds.set(b.id, v); }
  } catch (e) { }
  mol.delete(); return out;
}
/* every potential stereocentre (assigned or not) and every stereogenic C=C */
function centres(m) { const c = cip(m); return { atoms: c.all, bonds: c.allB }; }
/* canonical isomeric SMILES of the drawing */
function noH(s) { if (!s || s.indexOf('[H]') < 0) return s; const q = RD().get_mol(s); if (!q) return s; const t = q.get_smiles(); q.delete(); return t; }
function iso(m, ids) { const sub = ids ? m.sub(ids) : m; const { mol } = rdOf(sub); if (!mol) return null; const s = mol.get_smiles(); mol.delete(); return noH(s); }
function flat(m, ids) { const sub = ids ? m.sub(ids) : m; const { mol } = rdOf(sub, { noStereo: true }); if (!mol) return null; const s = mol.get_smiles(); mol.delete(); return s ? s.replace(/[\/\\]/g, '').replace(/@+/g, '') : s; }
/* same atoms + coordinates as `tpl`, wedges chosen by RDKit so that the drawing shows stereoisomer `smi` */
function wedgeTo(tpl, smi) {
  const { text, order } = mb(tpl, { noStereo: true }); const t = RD().get_mol(text, JSON.stringify({ removeHs: false })); const q = RD().get_mol(noH(smi));
  if (!t || !q) { t && t.delete(); q && q.delete(); return null; }
  let map = null; try { map = JSON.parse(q.generate_aligned_coords(t, JSON.stringify({ acceptFailure: false, allowRGroups: false }))); } catch (e) { map = null; }
  if (!map || !map.atoms || map.atoms.length !== order.length) { t.delete(); q.delete(); return null; }
  const L = q.get_molblock().split('\n'); const na = +L[3].slice(0, 3), nb = +L[3].slice(3, 6);
  const qi2id = new Map(); map.atoms.forEach((qi, k) => qi2id.set(qi, order[k])); // template atom k ↔ query atom qi
  const out = tpl.clone(); for (const b of out.bonds.values()) if (b.st === 1 || b.st === 6) b.st = 0;
  for (let i = 0; i < nb; i++) { const s = L[4 + na + i]; const a = +s.slice(0, 3) - 1, b = +s.slice(3, 6) - 1, st = +s.slice(9, 12);
    if (st !== 1 && st !== 6) continue; const ia = qi2id.get(a), ib = qi2id.get(b); const bd = out.bondBetween(ia, ib); if (!bd) continue; bd.st = st; if (bd.a !== ia) { bd.b = bd.a; bd.a = ia; } }
  t.delete(); q.delete(); return out;
}
/* mirror image: every wedge ↔ hash (coordinates unchanged), and E/Z is unaffected */
function mirror(m) { const o = m.clone(); for (const b of o.bonds.values()) { if (b.st === 1) b.st = 6; else if (b.st === 6) b.st = 1; } return o; }
/* chirality of a drawn compound */
function describe(m) {
  const c = centres(m); const lab = cip(m); const n = c.atoms.size;
  const assigned = [...c.atoms].filter(id => lab.atoms.has(id)); const s = iso(m); const flipped = iso(mirror(m));
  const chiral = assigned.length > 0 && s && flipped && s !== flipped;
  const meso = n >= 2 && assigned.length === n && s === flipped;
  const total = Math.pow(2, n + c.bonds.size);
  return { centres: c.atoms, dbl: c.bonds, labels: lab, n, assigned: assigned.length, chiral, meso, total, unassigned: n - assigned.length };
}

/* ---------------- carrying a configuration through a reaction ----------------
   The engine redraws products, so a wedge copied onto new coordinates can mean the other enantiomer.
   Handedness is kept as the sign of a signed volume over the centre's neighbours (wedge z = +1, hash −1),
   then the wedges are chosen again on the new drawing. One swapped neighbour: inversion for SN2, else retention. */
function vec3(m, c, n) { const a = m.atoms.get(c), b = m.atoms.get(n); const bd = m.bondBetween(c, n); let z = 0; if (bd && bd.a === c) { if (bd.st === 1) z = 1; else if (bd.st === 6) z = -1; } const dx = b.x - a.x, dy = b.y - a.y, l = Math.hypot(dx, dy) || 1; return [dx / l, dy / l, z]; }
const det3 = (a, b, c) => a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
const sub3 = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
function handed(m, c, order) { const v = order.map(n => vec3(m, c, n)); if (!v.some(x => x[2])) return 0; const d = v.length === 3 ? det3(v[0], v[1], v[2]) : det3(sub3(v[1], v[0]), sub3(v[2], v[0]), sub3(v[3], v[0])); return Math.abs(d) < 1e-6 ? 0 : Math.sign(d); }
function subMode(res, rg) {
  const text = res ? [res.title || '', ...(res.why || []), ...(res.steps || []).map(s => s.text || '')].join(' ') : '';
  const pyr = /pyridine|\bpy\b|Et3N|base/i.test(rg || '');
  if (/\bSNi\b|retention/i.test(text)) return pyr ? 'inv' : 'ret';
  if (/\bSN2\b|inversion|back ?side|Walden/i.test(text)) return 'inv';
  if (/carbocation|\bSN1\b|planar/i.test(text)) return 'rac';
  return 'ret';
}
/* RDKit chirality tag of each atom with the neighbour ids in the order the tag refers to (bond order of the JSON) */
function tags(m, noStereo) {
  const { mol, order } = rdOf(m, { noStereo }); if (!mol) return null; const j = JSON.parse(mol.get_json()); mol.delete();
  const M = j.molecules[0]; const nb = M.atoms.map(() => []); (M.bonds || []).forEach(b => { nb[b.atoms[0]].push(order[b.atoms[1]]); nb[b.atoms[1]].push(order[b.atoms[0]]); });
  const t = new Map(); M.atoms.forEach((a, i) => t.set(order[i], { st: a.stereo || null, nb: nb[i], i })); return { j, t, order };
}
function parity(a, b) { /* parity of the permutation taking list a to list b (same members) */ const p = a.map(x => b.indexOf(x)); if (p.some(x => x < 0)) return null; let n = 0; const seen = p.map(() => false); for (let i = 0; i < p.length; i++) { if (seen[i]) continue; let k = i, len = 0; while (!seen[k]) { seen[k] = true; k = p[k]; len++; } n += len - 1; } return n % 2; }
function carry(src, prod, res, rg) {
  if (!src || !prod) return 0; const sc = cip(src); if (!sc.atoms.size) return 0;
  const inv = subMode(res, rg) === 'inv'; const ts = tags(src); const tp = tags(prod, true); if (!ts || !tp) return 0; let n = 0;
  for (const c of sc.atoms.keys()) {
    const A = ts.t.get(c), B = tp.t.get(c); if (!A || !B || (A.st !== 'cw' && A.st !== 'ccw')) continue;
    if (A.nb.length !== B.nb.length || A.nb.length < 3 || src.hTot(c) !== prod.hTot(c)) continue;
    const gone = A.nb.filter(x => !B.nb.includes(x)), added = B.nb.filter(x => !A.nb.includes(x)); if (gone.length > 1 || gone.length !== added.length) continue;
    const mapped = A.nb.map(x => x === gone[0] ? added[0] : x); const par = parity(mapped, B.nb); if (par == null) continue;
    let st = A.st; if (par) st = st === 'cw' ? 'ccw' : 'cw'; if (gone.length && inv) st = st === 'cw' ? 'ccw' : 'cw';
    tp.j.molecules[0].atoms[B.i].stereo = st; n++;
  }
  if (!n) return 0; delete tp.j.molecules[0].conformers;
  let smi = null; try { const q = RD().get_mol(JSON.stringify(tp.j)); if (q) { smi = q.get_smiles(); q.delete(); } } catch (e) { smi = null; }
  if (!smi) return 0; const w = wedgeTo(prod, smi); if (!w) return 0;
  for (const b of prod.bonds.values()) { const x = w.bonds.get(b.id); if (x) { b.st = x.st; b.a = x.a; b.b = x.b; } }
  return n;
}
/* ---------------- stereospecific E2: anti-periplanar H and leaving group ----------------
   Newman view along Cα→Cβ. The leaving group X on Cα sits at 0°, the H removed from Cβ at 180°.
   Each carbon's chirality fixes where its other two groups sit (±120° / 180°∓120°); groups at 120° and 60°
   end up cis on the new C=C, and so do 240° and 300°. The product's double bond is then drawn that way. */
function withH(m, ids) { const d = m.clone(); for (const c of ids) { let h = d.hImp(c); while (h-- > 0) { const [ux, uy] = freeDirAt(d, c); const a = d.atoms.get(c); const id = d.addAtom('H', a.x + ux * BL * .7, a.y + uy * BL * .7); d.addBond(c, id, 1); if (a.hFix != null) a.hFix = Math.max(0, a.hFix - 1); } } return d; }
function tagFrom(T, c, first) { const A = T.t.get(c); if (!A || (A.st !== 'cw' && A.st !== 'ccw') || A.nb.length !== 4) return null; const want = [first, ...A.nb.filter(x => x !== first)]; const par = parity(want, A.nb); let st = A.st; if (par) st = st === 'cw' ? 'ccw' : 'cw'; return { st, rest: want.slice(1) }; }
function e2Geometry(src, prod, res) {
  if (!src || !prod) return false; const text = [res.title || '', ...(res.why || []), ...(res.steps || []).map(s => s.text || '')].join(' ');
  if (!/\bE2\b|anti-?periplanar|elimination|alc\.? ?KOH|dehydrohalogenation|vinyl halide/i.test(text)) return false;
  const sc = cip(src); if (sc.atoms.size === 0) return false; let did = false;
  for (const pb of prod.bonds.values()) {
    if (pb.o !== 2) continue; const sb = src.bondBetween(pb.a, pb.b); if (!sb || sb.o !== 1) continue;
    if (C().rings(prod, 8).some(r => r.includes(pb.a) && r.includes(pb.b))) continue;
    /* which end lost the heavy leaving group, which lost an H */
    let al = null, be = null, X = null;
    for (const [p, q] of [[pb.a, pb.b], [pb.b, pb.a]]) { const lost = src.nbrs(p).filter(n => !prod.atoms.has(n) || !prod.bondBetween(p, n)).filter(n => src.atoms.get(n).el !== 'H'); if (lost.length === 1 && src.hTot(q) - prod.hTot(q) === 1) { al = p; be = q; X = lost[0]; break; } }
    if (al == null) continue;
    const sH = withH(src, [al, be]); const T = tags(sH); if (!T) continue;
    const Hs = sH.nbrs(be).filter(n => sH.atoms.get(n).el === 'H' && !src.atoms.has(n)).concat(src.nbrs(be).filter(n => src.atoms.get(n).el === 'H'));
    if (Hs.length !== 1) continue; /* CH₂: either H can go, no stereospecific outcome */
    const ta = tagFrom(T, al, be), tb = tagFrom(T, be, al); if (!ta || !tb) continue;
    /* α, viewed from the front: from β the order X→a1→a2 is ta.st; the viewer sees it reversed */
    const ra = [X, ...ta.rest.filter(x => x !== X)]; const pa = parity(ra, ta.rest); let sa = ta.st; if (pa) sa = sa === 'cw' ? 'ccw' : 'cw';
    const ang = new Map(); ang.set(ra[1], sa === 'cw' ? 120 : 240); ang.set(ra[2], sa === 'cw' ? 240 : 120);
    const H = Hs[0]; const rb = [H, ...tb.rest.filter(x => x !== H)]; const pbp = parity(rb, tb.rest); let sB = tb.st; if (pbp) sB = sB === 'cw' ? 'ccw' : 'cw';
    ang.set(rb[1], sB === 'cw' ? 60 : 300); ang.set(rb[2], sB === 'cw' ? 300 : 60);
    /* pick one heavy group on each side that survives in the product */
    const ga = ra.slice(1).find(x => prod.atoms.has(x) && prod.bondBetween(al, x)), gb = rb.slice(1).find(x => prod.atoms.has(x) && prod.bondBetween(be, x));
    if (ga == null || gb == null) continue;
    const upA = ang.get(ga) === 120, upB = ang.get(gb) === 60; const wantCis = upA === upB;
    const A0 = prod.atoms.get(al), B0 = prod.atoms.get(be); const side = id => { const t = prod.atoms.get(id); return Math.sign((B0.x - A0.x) * (t.y - A0.y) - (B0.y - A0.y) * (t.x - A0.x)); };
    const isCis = side(ga) === side(gb);
    if (isCis !== wantCis) { /* reflect β's side across the double-bond axis */
      const seen = new Set([be]); const q = [be]; while (q.length) { const x = q.shift(); for (const y of prod.nbrs(x)) if (y !== al && !seen.has(y)) { seen.add(y); q.push(y); } } seen.delete(be);
      const dx = B0.x - A0.x, dy = B0.y - A0.y, l2 = dx * dx + dy * dy || 1;
      for (const id of seen) { const t = prod.atoms.get(id); const k = ((t.x - A0.x) * dx + (t.y - A0.y) * dy) / l2; const px = A0.x + k * dx, py = A0.y + k * dy; t.x = 2 * px - t.x; t.y = 2 * py - t.y; }
      for (const b of prod.bonds.values()) if (seen.has(b.a) && seen.has(b.b) && (b.st === 1 || b.st === 6)) b.st = b.st === 1 ? 6 : 1;
    }
    did = true;
  }
  return did;
}
/* ---------------- reaction outcome: which stereoisomers form ---------------- */
const ANTI = /\banti(?![- ]?mark)\b|halonium|bromonium|trans-dihal|anti addition|trans-1,2-diol|back side/i, SYN = /\bsyn\b|cis-|concerted|cyclic ester|osmate|manganate|same face|same side|Lindlar/i;
function addedTo(src, prod) {
  /* alkene/alkyne carbons that became sp3 (or lost a π bond) and what they gained */
  const out = [];
  for (const b of src.bonds.values()) { if (b.o < 2) continue; const pb = prod.bondBetween(b.a, b.b); if (!pb || pb.o >= b.o) continue;
    const ga = prod.nbrs(b.a).filter(n => !src.bondBetween(b.a, n)), gb = prod.nbrs(b.b).filter(n => !src.bondBetween(b.b, n));
    const ha = prod.hTot(b.a) - src.hTot(b.a), hb = prod.hTot(b.b) - src.hTot(b.b); out.push({ a: b.a, b: b.b, ga, gb, ha, hb, ring: pb.o }); }
  return out;
}
function freeDirAt(m, id, skip) { const a = m.atoms.get(id); const dirs = m.nbrs(id).filter(n => !(skip && skip.has(n))).map(n => { const q = m.atoms.get(n); return Math.atan2(q.y - a.y, q.x - a.x); }); const ang = dirs.length ? C().freeAngle(dirs, a, m) : -Math.PI / 2; return [Math.cos(ang), Math.sin(ang)]; }
/* product drawn on the source's 2D geometry; the groups added to the alkene carbons are wedged (syn) or wedge/hash (anti) */
function addDrawing(src, prod, add, mode) {
  const d = prod.clone(); for (const b of d.bonds.values()) if (b.st) b.st = 0;
  for (const a of d.atoms.values()) { const s = src.atoms.get(a.id); if (s) { a.x = s.x; a.y = s.y; } }
  const placed = new Set([...d.atoms.keys()].filter(id => src.atoms.has(id)));
  const put = (c, g, st) => { /* group root g on carbon c */
    const [ux, uy] = freeDirAt(d, c, new Set([g, ...d.nbrs(c).filter(n => !placed.has(n))])); const ca = d.atoms.get(c); const ga = d.atoms.get(g);
    const dx = ca.x + ux * BL - ga.x, dy = ca.y + uy * BL - ga.y; const sub = C().subtree ? null : null;
    // move the whole new group rigidly (atoms not in the source reachable from g without crossing c)
    const seen = new Set([g]); const q = [g]; while (q.length) { const x = q.shift(); for (const y of d.nbrs(x)) if (y !== c && !seen.has(y) && !placed.has(y)) { seen.add(y); q.push(y); } }
    for (const id of seen) { const t = d.atoms.get(id); t.x += dx; t.y += dy; placed.add(id); }
    const bd = d.bondBetween(c, g); bd.st = st; if (bd.a !== c) { bd.b = bd.a; bd.a = c; } };
  const hAt = (c, st) => { const [ux, uy] = freeDirAt(d, c); const ca = d.atoms.get(c); if (ca.hFix != null) ca.hFix = Math.max(0, ca.hFix - 1); const h = d.addAtom('H', ca.x + ux * BL * .8, ca.y + uy * BL * .8); d.addBond(c, h, 1, { st }); const bd = d.bondBetween(c, h); bd.st = st; bd.a = c; bd.b = h; placed.add(h); };
  for (const x of add) {
    const sa = 1, sb = mode === 'anti' ? 6 : 1; let did = false;
    if (x.ga.length) { put(x.a, x.ga[0], sa); did = true; } else if (x.ha > 0) { hAt(x.a, sa); did = true; }
    if (x.gb.length) { put(x.b, x.gb[0], did ? sb : sa); } else if (x.hb > 0) { hAt(x.b, did ? sb : sa); }
  }
  return d;
}
function uniqIso(list) { const seen = new Map(); for (const m of list) { const k = iso(m); if (k && !seen.has(k)) seen.set(k, m); } return [...seen.entries()]; }
/* all stereoisomers by setting every unassigned stereocentre both ways (≤ 4 centres) */
function enumerate(prod, only) {
  const c = centres(prod); const lab = cip(prod); const free = [...c.atoms].filter(id => !lab.atoms.has(id) && (!only || only.includes(id))).slice(0, 4);
  if (!free.length) return [[iso(prod), prod]];
  const { mol, order } = rdOf(prod); if (!mol) return [[iso(prod), prod]]; const j = JSON.parse(mol.get_json()); mol.delete(); delete j.molecules[0].conformers; j.molecules[0].atoms.forEach(a => { delete a.stereo; });
  const idx = new Map(order.map((id, k) => [id, k])); const out = new Map();
  for (let mask = 0; mask < (1 << free.length); mask++) {
    const jj = JSON.parse(JSON.stringify(j)); free.forEach((id, k) => { jj.molecules[0].atoms[idx.get(id)].stereo = mask & (1 << k) ? 'cw' : 'ccw'; });
    let q = null; try { q = RD().get_mol(JSON.stringify(jj)); } catch (e) { q = null; } if (!q) continue; const s = q.get_smiles(); q.delete(); if (out.has(s)) continue;
    const w = wedgeTo(prod, s); if (w) out.set(s, w); else root.__stDbg = 'wedgeTo failed ' + s;
  }
  return [...out.entries()];
}
function stereoNote(list, kind) {
  const n = list.length;
  if (kind === 'meso') return 'meso compound: it has stereocentres but is superimposable on its mirror image, so it is optically inactive.';
  if (kind === 'racemic') return 'Both enantiomers form in equal amounts: a racemic mixture (±), optically inactive.';
  if (kind === 'diastereomers') return `${n} stereoisomers form (diastereomers and their enantiomers).`;
  if (kind === 'pair-meso') return 'A meso compound and a racemic pair form.';
  return '';
}
/* outcome(src, res): { list: [{ mol, smi, labels }], kind, note } or null when nothing stereo happens */
function outcome(src, res, rg) {
  if (!res || !res.product) return null; const prod = res.product; if (src && !res._carried) { try { carry(src, prod, res, rg); } catch (e) { } try { if (e2Geometry(src, prod, res)) res._e2 = true; } catch (e) { } res._carried = true; }
  const pc = centres(prod); if (!pc.atoms.size && !pc.bonds.size) return null;
  let plab = cip(prod);
  const text0 = [res.title || '', ...(res.why || []), ...(res.steps || []).map(s => s.text || '')].join(' ');
  if (src && subMode(res, rg) === 'rac') {
    /* a stereocentre that went through a flat cation loses its configuration */
    const flatIds = [...plab.atoms.keys()].filter(id => src.atoms.has(id) && prod.nbrs(id).some(n => !src.bondBetween(id, n)));
    if (flatIds.length) { const p2 = prod.clone(); for (const b of p2.bonds.values()) if ((b.st === 1 || b.st === 6) && (flatIds.includes(b.a))) b.st = 0; return outcome(src, Object.assign({}, res, { product: p2, _flat: flatIds }), rg); }
  }
  const unassigned = [...pc.atoms].filter(id => !plab.atoms.has(id));
  if (!unassigned.length) { const d = describe(prod); const ez = [...plab.bonds.values()]; const ezN = ez.length ? 'Gives the ' + ez.map(v => v === 'Z' ? 'Z (cis)' : 'E (trans)').join(', ') + ' alkene' + (ez.length > 1 ? 's' : '') + ' only.' : '';
    return { list: [{ mol: prod, smi: iso(prod), labels: plab }], kind: d.meso ? 'meso' : d.chiral ? 'single' : 'achiral', note: [d.meso ? stereoNote([], 'meso') : d.chiral ? 'A single stereoisomer: the configuration is set by the starting material.' : '', ezN].filter(Boolean).join(' ') }; }
  /* centres that were already there, undrawn, in the starting material stay undrawn */
  const sc = src ? cip(src) : null; const fresh = unassigned.filter(id => (res._flat && res._flat.includes(id)) || !(sc && sc.all.has(id) && !sc.atoms.has(id)));
  if (!fresh.length) return { list: [{ mol: prod, smi: iso(prod), labels: plab }], kind: 'undrawn', note: 'The starting material has no configuration drawn at its stereocentre' + (unassigned.length > 1 ? 's' : '') + ', so the product has none either. Draw a wedge or hash on the start to follow it through.' };
  const text = [res.title || '', ...(res.why || []), ...(res.steps || []).map(s => s.text || '')].join(' ');
  const tags = (res.tags || []).join(' ');
  const mode = /\banti(?![- ]?mark)\b/i.test(tags) || ANTI.test(text) ? 'anti' : /\bsyn\b/.test(tags) || SYN.test(text) ? 'syn' : null;
  let list = null;
  const add = src ? addedTo(src, prod).filter(x => unassigned.includes(x.a) || unassigned.includes(x.b)) : [];
  if (mode && add.length) {
    const d1 = addDrawing(src, prod, add, mode); const L = uniqIso([d1, mirror(d1)]);
    if (L.length) list = L.map(([s, m]) => ({ mol: wedgeTo(prod, s) || m, smi: s }));
  }
  if (!list) list = enumerate(prod, fresh).map(([s, m]) => ({ mol: m, smi: s }));
  list.forEach(x => { x.labels = cip(x.mol); });
  let kind = 'single';
  if (list.length === 1) { const d = describe(list[0].mol); kind = d.meso ? 'meso' : 'single'; }
  else if (list.length === 2 && iso(mirror(list[0].mol)) === list[1].smi) kind = 'racemic';
  else { const mesoN = list.filter(x => iso(mirror(x.mol)) === x.smi).length; kind = mesoN && list.length === 3 ? 'pair-meso' : 'diastereomers'; }
  const how = mode && add.length ? (mode === 'anti' ? 'Anti addition: the two new groups end up on opposite faces.' : 'Syn addition: both new groups add to the same face.') : unassigned.length ? 'The new stereocentre forms from a flat (planar) intermediate or face, so attack from either side is equally likely.' : '';
  return { list, kind, mode, note: [how, stereoNote(list, kind)].filter(Boolean).join(' ') };
}

/* ---------------- Fischer projection ---------------- */
const GROUP = [['CHO', '[CX3H1](=O)'], ['COOH', 'C(=O)[OH]'], ['COO⁻', 'C(=O)[O-]'], ['CH₂OH', '[CH2][OH]'], ['CN', 'C#N'], ['CONH₂', 'C(=O)N'], ['COOCH₃', 'C(=O)OC'], ['CH₃', '[CH3]'], ['C₂H₅', '[CH2][CH3]']];
function label(m, root, from) {
  /* condensed label for the substituent at `root` (coming from `from`) */
  const el = m.atoms.get(root).el; const ids = []; const st = [root]; const seen = new Set([from, root]);
  while (st.length) { const x = st.pop(); ids.push(x); for (const y of m.nbrs(x)) if (!seen.has(y)) { seen.add(y); st.push(y); } }
  const h = m.hTot(root);
  if (ids.length === 1) { if (el === 'C') return h === 3 ? 'CH₃' : h ? 'CH' + (h > 1 ? '₂' : '') : 'C'; return h ? (el === 'O' || el === 'S' ? el + 'H' : el === 'N' ? (h === 2 ? 'NH₂' : 'NH') : el + 'H' + (h > 1 ? '₂' : '')) : el + (m.atoms.get(root).chg < 0 ? '⁻' : ''); }
  const sub = m.sub(ids); let smi = null; try { smi = C().canon(sub); } catch (e) { }
  const T = { 'C=O': 'CHO', 'O=CO': 'COOH', 'O=C[O-]': 'COO⁻', 'OC': el === 'C' ? 'CH₂OH' : 'OCH₃', 'C#N': 'CN', 'NC=O': 'CONH₂', 'CC': 'C₂H₅', 'COC=O': 'COOCH₃', 'CCC': 'C₃H₇', 'CC(C)C': 'CH(CH₃)₂', 'c1ccccc1': 'C₆H₅', 'CO': el === 'C' ? 'CH₂OH' : 'OCH₃', 'CN': el === 'C' ? 'CH₂NH₂' : 'NHCH₃', 'CBr': 'CH₂Br', 'CCl': 'CH₂Cl', 'CC=O': 'CH₂CHO', 'CC(=O)O': 'CH₂COOH', 'OCC': 'OC₂H₅' };
  if (smi && T[smi]) return T[smi];
  const cnt = {}; for (const id of ids) { const a = m.atoms.get(id); cnt[a.el] = (cnt[a.el] || 0) + 1; cnt.H = (cnt.H || 0) + m.hImp(id); }
  const sb = n => n > 1 ? String(n).split('').map(d => '₀₁₂₃₄₅₆₇₈₉'[+d]).join('') : '';
  return Object.keys(cnt).sort((a, b) => (a === 'C' ? -2 : a === 'H' ? -1 : 0) - (b === 'C' ? -2 : b === 'H' ? -1 : 0) || a.localeCompare(b)).map(k => k + sb(cnt[k])).join('');
}
/* longest carbon chain through the stereocentres, top = the more oxidised end */
function heavyOut(m, id, nb) { const seen = new Set([nb, id]); const st = [id]; let n = 0; while (st.length) { const x = st.pop(); n++; for (const y of m.nbrs(x)) if (!seen.has(y) && m.atoms.get(y).el !== 'H') { seen.add(y); st.push(y); } } return n; }
function fischerChain(m, sc) {
  const C_ = [...m.atoms.keys()].filter(id => m.atoms.get(id).el === 'C' && !C().rings(m, 8).some(r => r.includes(id)));
  if (!sc.length || sc.some(id => !C_.includes(id))) return null;
  const set = new Set(C_); let best = null;
  const ends = C_.filter(id => m.nbrs(id).filter(n => set.has(n)).length <= 1);
  for (const s of ends) { const prev = new Map([[s, null]]); const q = [s]; while (q.length) { const x = q.shift(); for (const y of m.nbrs(x)) if (set.has(y) && !prev.has(y)) { prev.set(y, x); q.push(y); } }
    for (const e of ends) { if (e === s || !prev.has(e)) continue; const path = []; let x = e; while (x != null) { path.unshift(x); x = prev.get(x); }
      if (!sc.every(id => path.includes(id))) continue; if (!best || path.length > best.length) best = path; } }
  if (!best) { if (sc.length === 1) { const c = sc[0]; const nb = m.nbrs(c).filter(n => set.has(n)); if (nb.length >= 1) best = nb.length >= 2 ? [nb[0], c, nb[1]] : [nb[0], c]; } }
  if (!best) return null;
  const ox = id => m.nbrs(id).filter(n => m.atoms.get(n).el !== 'C' && m.atoms.get(n).el !== 'H').reduce((s, n) => s + m.bondBetween(id, n).o, 0);
  /* only one carbon beyond the outermost stereocentres: what lies further out is written as a group (C₂H₅, CH₂OH…) */
  const i0 = best.findIndex(id => sc.includes(id)), i1 = best.length - 1 - [...best].reverse().findIndex(id => sc.includes(id));
  best = best.slice(Math.max(0, i0 - 1), i1 + 2);
  if (ox(best[best.length - 1]) > ox(best[0]) || (ox(best[best.length - 1]) === ox(best[0]) && heavyOut(m, best[best.length - 1], best[best.length - 2]) > heavyOut(m, best[0], best[1]))) best.reverse();
  return best;
}
function fischer(m) {
  const d = describe(m); const sc = [...d.centres].filter(id => d.labels.atoms.has(id));
  if (!sc.length) return { ok: false, reason: d.centres.size ? 'Give the stereocentres a configuration first (wedge or hash bonds), then the Fischer projection can be drawn.' : 'No stereocentres, so there is nothing to project.' };
  const chain = fischerChain(m, sc); if (!chain) return { ok: false, reason: 'Fischer projections need the stereocentres on an open (acyclic) carbon chain.' };
  /* build a Fischer drawing: vertical chain hashed, horizontal groups wedged; then fix left/right per centre by comparing R/S */
  const f = m.clone(); for (const b of f.bonds.values()) b.st = 0; const pos = new Map(); const step = BL * 1.25;
  chain.forEach((id, k) => pos.set(id, [0, k * step]));
  const horiz = new Map();
  for (const c of chain) { if (!sc.includes(c)) continue; const k = chain.indexOf(c); const side = f.nbrs(c).filter(n => !chain.includes(n)); let hs = f.hTot(c) - f.nbrs(c).filter(n => f.atoms.get(n).el === 'H').length;
    const L = [...side]; while (L.length < 2 && hs > 0) { const h = f.addAtom('H', 0, 0); f.addBond(c, h, 1); f.atoms.get(c).hFix = f.atoms.get(c).hFix != null ? Math.max(0, f.atoms.get(c).hFix - 1) : null; L.push(h); hs--; }
    horiz.set(c, L.slice(0, 2)); }
  const place = flipSet => {
    const g = f.clone(); for (const b of g.bonds.values()) b.st = 0; const fixed = new Set();
    for (const [id, [x, y]] of pos) { const a = g.atoms.get(id); a.x = x; a.y = y; fixed.add(id); }
    for (const [c, [l, r]] of horiz) { const [cx, cy] = pos.get(c); const flip = flipSet.has(c); const [L, R] = flip ? [r, l] : [l, r];
      const mv = (root, x) => { const a = g.atoms.get(root); const dx = x - a.x, dy = cy - a.y; const seen = new Set([root]); const q = [root]; while (q.length) { const z = q.shift(); for (const w of g.nbrs(z)) if (w !== c && !seen.has(w) && !pos.has(w)) { seen.add(w); q.push(w); } } for (const z of seen) { const t = g.atoms.get(z); t.x += dx; t.y += dy; } };
      for (const [root, x] of [[L, cx - BL], [R, cx + BL]]) { const a = g.atoms.get(root); a.x = x; a.y = cy; fixed.add(root); }
      for (const n of [L, R]) { const bd = g.bondBetween(c, n); bd.st = 1; if (bd.a !== c) { bd.b = bd.a; bd.a = c; } }
      for (const n of g.nbrs(c).filter(n => chain.includes(n))) { const bd = g.bondBetween(c, n); bd.st = 6; if (bd.a !== c) { bd.b = bd.a; bd.a = c; } } }
    /* everything else (end groups, the rest of each horizontal group) is grown out from what is placed */
    for (const a of g.atoms.values()) a.isNew = !fixed.has(a.id);
    /* end groups point away from the chain: seed the first atom beyond each end straight up / down */
    const ends = [[chain[0], chain[1], -1], [chain[chain.length - 1], chain[chain.length - 2], 1]];
    for (const [e, nb, dir] of ends) { if (e == null) continue; const ea = g.atoms.get(e); for (const n of g.nbrs(e)) { if (fixed.has(n) || n === nb) continue; const a = g.atoms.get(n); if (!a.isNew) continue; const k = g.nbrs(e).filter(x => !fixed.has(x) && x !== nb).indexOf(n); const ang = dir < 0 ? -Math.PI / 2 + (k ? (k % 2 ? -1 : 1) * Math.PI / 3 : 0) : Math.PI / 2 + (k ? (k % 2 ? 1 : -1) * Math.PI / 3 : 0); a.x = ea.x + Math.cos(ang) * BL; a.y = ea.y + Math.sin(ang) * BL; a.isNew = false; } }
    try { C().placeNew(g, BL); } catch (e) { }
    for (const a of g.atoms.values()) a.isNew = false;
    return g; };
  let flips = new Set(); let g = place(flips); let lab = cip(g).atoms;
  for (const c of sc) if (lab.get(c) !== d.labels.atoms.get(c)) flips.add(c);
  g = place(flips); lab = cip(g).atoms;
  const okAll = sc.every(c => lab.get(c) === d.labels.atoms.get(c));
  /* rows for rendering */
  const rows = chain.map((id, k) => { const isC = horiz.has(id); const [l, r] = isC ? (flips.has(id) ? [horiz.get(id)[1], horiz.get(id)[0]] : horiz.get(id)) : [null, null];
    return { id, center: isC, left: l != null ? label(f, l, id) : null, right: r != null ? label(f, r, id) : null, lab: d.labels.atoms.get(id) || '' }; });
  const top = chain[0], bot = chain[chain.length - 1];
  const endLabel = (id, nb) => { if (horiz.has(id)) return null; const keep = new Set([id]); const st = [id]; const seen = new Set([nb, id]); while (st.length) { const x = st.pop(); for (const y of f.nbrs(x)) if (!seen.has(y) && !chain.includes(y)) { seen.add(y); keep.add(y); st.push(y); } } return label(f, id, nb); };
  return { ok: okAll, rows, top: endLabel(top, chain[1]), bottom: endLabel(bot, chain[chain.length - 2]), chain, flips, describe: d, mol: g };
}
function fischerSVG(F, opt = {}) {
  if (!F || !F.rows) return ''; const step = 46, w = 64; const rows = F.rows; let y = 30, o = ''; const cx = 120;
  const T = (x, y, s, cls) => `<text class="${cls || 'fl'}" x="${x}" y="${y}">${String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;')}</text>`;
  const pos = [];
  rows.forEach((r, k) => {
    const yy = y + k * step; pos.push(yy);
    if (r.center) { o += `<line class="fv" x1="${cx}" y1="${yy - step / 2 + 6}" x2="${cx}" y2="${yy + step / 2 - 6}"/><line class="fh" x1="${cx - w}" y1="${yy}" x2="${cx + w}" y2="${yy}"/>`;
      o += T(cx - w - 8, yy + 5, r.left, 'fl l') + T(cx + w + 8, yy + 5, r.right, 'fl r'); if (opt.labels && r.lab) o += T(cx + 14, yy - 8, r.lab, 'fcip'); }
    else if (k === 0) o += T(cx, yy + 4, F.top || 'C', 'fl c'); else if (k === rows.length - 1) o += T(cx, yy + 6, F.bottom || 'C', 'fl c');
    else o += `<line class="fv" x1="${cx}" y1="${yy - step / 2 + 6}" x2="${cx}" y2="${yy + step / 2 - 6}"/>${T(cx + 10, yy + 4, 'CH₂', 'fl')}`;
  });
  for (let k = 0; k < rows.length - 1; k++) { const a = rows[k], b = rows[k + 1]; if (a.center || b.center) continue; o += `<line class="fv" x1="${cx}" y1="${pos[k] + 10}" x2="${cx}" y2="${pos[k + 1] - 12}"/>`; }
  if (!rows[0].center) o += `<line class="fv" x1="${cx}" y1="${pos[0] + 10}" x2="${cx}" y2="${pos[1] - step / 2 + 6}"/>`;
  const n = rows.length; if (!rows[n - 1].center && n > 1) o += `<line class="fv" x1="${cx}" y1="${pos[n - 2] + step / 2 - 6}" x2="${cx}" y2="${pos[n - 1] - 12}"/>`;
  const H = y + (rows.length - 1) * step + 30;
  return `<svg class="bx-fischer" viewBox="0 0 240 ${H}" role="img" aria-label="Fischer projection">${o}</svg>`;
}

root.Stereo = { carry, handed, mb, cip, centres, iso, flat, wedgeTo, mirror, describe, outcome, enumerate, fischer, fischerSVG, label };
})(typeof window !== 'undefined' ? window : globalThis);
