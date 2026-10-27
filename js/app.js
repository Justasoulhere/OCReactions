/* app.js · OC Reactions v3 core: data helpers, Reactions list, viewer, Reagents index, command palette, theme */
(() => {
'use strict';
const D = window.DB, ART = window.ART;
const R = D.R, BY = new Map(R.map(r => [r.id, r]));
const $ = (s, el = document) => el.querySelector(s), $$ = (s, el = document) => [...el.querySelectorAll(s)];
const esc = s => String(s == null ? '' : s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const SUBD = '₀₁₂₃₄₅₆₇₈₉', SUPD = '⁰¹²³⁴⁵⁶⁷⁸⁹';
function fmt(s) {
  let out = '', mode = 0;
  const open = m => m === 1 ? '<sub>' : m === 2 ? '<sup>' : '', close = m => m === 1 ? '</sub>' : m === 2 ? '</sup>' : '';
  for (const ch of String(s == null ? '' : s)) {
    let m = 0, c = ch, i;
    if ((i = SUBD.indexOf(ch)) >= 0) { m = 1; c = String(i); } else if ((i = SUPD.indexOf(ch)) >= 0) { m = 2; c = String(i); } else if (ch === '⁺') { m = 2; c = '+'; } else if (ch === '⁻') { m = 2; c = '\u2212'; }
    if (m !== mode) { out += close(mode) + open(m); mode = m; }
    out += esc(c);
  }
  return out + close(mode);
}
function norm(s) { let o = ''; for (const ch of String(s == null ? '' : s)) { let i = SUBD.indexOf(ch); if (i >= 0) { o += i; continue; } i = SUPD.indexOf(ch); if (i >= 0) { o += i; continue; } o += ch === '⁺' ? '+' : (ch === '⁻' || ch === '\u2212' || ch === '–' || ch === '—') ? '-' : ch; } return o.toLowerCase(); }
const debounce = (f, ms) => { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => f(...a), ms); }; };
const toast = (msg, ms = 2200) => { const t = $('#toast'); t.innerHTML = msg; t.classList.add('on'); clearTimeout(toast._t); toast._t = setTimeout(() => t.classList.remove('on'), ms); };

/* ---------- vocabulary ---------- */
const CLS = { sub: ['Substitution', 250], elim: ['Elimination', 52], add: ['Addition', 145], ox: ['Oxidation', 27], red: ['Reduction', 300], rearr: ['Rearrangement', 345], cond: ['Condensation', 192], hydrolysis: ['Hydrolysis', 222], acidbase: ['Acid/base', 96], radical: ['Radical', 8], pericyclic: ['Pericyclic', 275], test: ['Test', 122], cyclisation: ['Cyclisation', 205], coupling: ['Coupling', 322], exchange: ['Exchange', 240, .35], isomerisation: ['Isomerisation', 168, .45], cleavage: ['Cleavage', 16], prep: ['Preparation', 232, .3], polymerisation: ['Polymerisation', 70, .4] };
const kcol = k => { const c = CLS[k] || ['', 232, .3]; return `oklch(var(--kl) calc(var(--kc) * ${c[2] || 1}) ${c[1]})`; };
const clsOf = r => (r.tags.cls && r.tags.cls.length ? r.tags.cls : ['prep']);
function barCSS(r) { const cs = clsOf(r); if (cs.length === 1) return kcol(cs[0]); const n = cs.length; return 'linear-gradient(' + cs.map((c, i) => `${kcol(c)} ${(i * 100 / n).toFixed(2)}% ${((i + 1) * 100 / n).toFixed(2)}%`).join(',') + ')'; }
const PRETTY = { markovnikov: 'Markovnikov', antimarkovnikov: 'anti-Markovnikov', saytzeff: 'Saytzeff', hofmann: 'Hofmann', ortho_para: 'ortho/para', '1_2': '1,2-addition', '1_4': '1,4-addition', less_hindered: 'less hindered site', more_substituted: 'more substituted', ortho_major: 'ortho major', para_if_ortho_blocked: 'para if ortho is blocked', Na_ortho: 'Na⁺ salt: ortho', K_para: 'K⁺ salt: para', cyclicTS: 'cyclic TS', sigma_complex: 'σ-complex', hydride_transfer: 'hydride transfer', SN2p: "SN2'", AdN_E: 'AdN-E', tetrahedral: 'tetrahedral intermediate', oxocarbenium: 'oxocarbenium ion', acylium: 'acylium ion', halonium: 'halonium ion', mercurinium: 'mercurinium ion', carbocation: 'carbocation', carbanion: 'carbanion', radical: 'radical', stereospecific: 'stereospecific', stereoselective: 'stereoselective' };
const pv = v => PRETTY[v] || String(v).replace(/_/g, ' ');
const FGL = Object.fromEntries(D.FG.map(n => [n.k, n.label]));
const ORDER = ['Reaction Mechanism', 'Electrophilic Addition', 'Free Radicals', 'Substitution & Elimination', 'Name Reactions (Aldol & more)', 'Biomolecules', 'Carbonyl Compounds', 'Grignard, Oxidation & Reduction', 'Conformational Isomerism', 'General Organic Chemistry', 'Geometrical Isomerism', 'Ester Hydrolysis', 'Structural Isomerism', 'Nomenclature', 'Stereoisomerism', 'Aromatic Compounds', 'Carboxylic Acids & Derivatives', 'Alkyl Halides', 'Alkyl Halides (Free Radical & Electrophilic Addition)'];
const chIdx = c => { const i = ORDER.indexOf(c); return i < 0 ? 98 : i; };
const prim = r => r.chs[0] || { chapter: 'Other' };
const chName = l => l.chapter;
const chText = l => chName(l);
const chLine = r => r.chs.length ? chText(r.chs[0]) + (r.chs.length > 1 ? `  +${r.chs.length - 1}` : '') : '';
R.forEach((r, i) => { const l = prim(r); r._k = chIdx(l.chapter) * 1e6 + i; });
function title(r) {
  if (r.name) return r.name;
  const [a, b] = r.fg || []; const mt = (r.tags.mt || [])[0];
  if (a && b && a !== b) return `${FGL[a]} → ${FGL[b]}` + (mt ? ` (${pv(mt)})` : '');
  return (FGL[a] || 'Reaction') + (mt ? `: ${pv(mt)}` : (clsOf(r)[0] ? `: ${CLS[clsOf(r)[0]][0].toLowerCase()}` : ''));
}
const natW = svg => { const m = /width='(\d+)'/.exec(svg || ''); return m ? +m[1] : 400; };

/* ---------- reagent index ---------- */
const RGI = new Map();
function addRg(label, id, inEx) {
  for (let part of String(label).split(/,\s+|\s+then\s+|;\s*/)) {
    part = part.replace(/^\(?(?:i{1,3}|iv|v|\d)[.)]\s*/i, '').trim(); if (!part || part.length > 42) continue;
    const k = norm(part).replace(/\s+/g, ''); let e = RGI.get(k);
    if (!e) RGI.set(k, e = { k, label: part, ids: new Set(), main: new Set() });
    e.ids.add(id); if (!inEx) e.main.add(id);
  }
}
for (const r of R) { for (const l of r.rl) addRg(l, r.id, false); for (const e of r.exl) for (const l of e) addRg(l, r.id, true); }
const RGL = [...RGI.values()].sort((a, b) => b.ids.size - a.ids.size || a.label.localeCompare(b.label));
const rgKeysOf = new Map(); for (const e of RGL) for (const id of e.ids) { if (!rgKeysOf.has(id)) rgKeysOf.set(id, []); rgKeysOf.get(id).push(e); }

/* ---------- state ---------- */
const S = { view: 'reactions', q: '', ch: new Set(), cls: new Set(), mt: '', chf: '', named: false, mech: false, sub: '', subHits: null, group: 'none', fa: '', fb: '', rg: RGL.length ? RGL[0].k : '' };
const hay = r => r._h || (r._h = norm([r.name, title(r), ...r.rl, ...r.cl, ...r.rel, ...Object.values(r.tags).flat().map(pv), ...r.chs.map(chText), FGL[r.fg[0]], FGL[r.fg[1]], ...clsOf(r).map(c => CLS[c] ? CLS[c][0] : c)].join(' ')));
function filtered() {
  const toks = norm(S.q).split(/\s+/).filter(Boolean);
  return R.filter(r => {
    if (S.ch.size && !r.chs.some(l => S.ch.has(chName(l)))) return false;
    if (S.chf && !r.chs.some(l => chText(l) === S.chf)) return false;
    if (S.cls.size && !clsOf(r).some(c => S.cls.has(c))) return false;
    if (S.mt && !(r.tags.mt || []).includes(S.mt)) return false;
    if (S.fa && !(r.fgs ? r.fgs[0] : [r.fg[0]]).includes(S.fa)) return false;
    if (S.fb && !(r.fgs ? r.fgs[1] : [r.fg[1]]).includes(S.fb)) return false;
    if (S.named && !r.name) return false;
    if (S.mech && !r.ns) return false;
    if (S.sub && S.subHits && !S.subHits.has(r.id)) return false;
    if (toks.length) { const h = hay(r); if (!toks.every(t => h.includes(t))) return false; }
    return true;
  }).sort((a, b) => a._k - b._k);
}

/* ---------- RDKit (lazy) ---------- */
let RDP = null;
function rdkit() {
  if (RDP) return RDP;
  RDP = new Promise((res, rej) => {
    const add = (src, cb) => { const s = document.createElement('script'); s.src = src; s.onload = cb; s.onerror = () => rej(new Error('load ' + src)); document.head.append(s); };
    add('js/vendor/rdkit_wasm.js', () => add('js/vendor/RDKit_minimal.js', () => {
      const bin = atob(window.RDKIT_WASM_B64); const u = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i);
      window.initRDKitModule({ instantiateWasm: (info, receive) => { WebAssembly.instantiate(u, info).then(r => receive(r.instance, r.module)); return {}; }, locateFile: f => f })
        .then(m => { window.RDK = m; try { m.disable_logging && m.disable_logging(); } catch (e) { } res(m); }).catch(rej);
    }));
  });
  return RDP;
}

/* ---------- micro-interactions ---------- */
/* press line: stays while the pointer is held down on a control and still inside it; clears on release,
   on cancel, when the pointer is dragged off the element, or if the window loses focus */
const PRESSABLE = '.btn,.chip,.seg button,.kbar,.row,.card,.mini,.pi';
let pressEl = null, pressSafety = 0;
function endPress() { if (pressEl) pressEl.classList.remove('pressing'); pressEl = null; clearTimeout(pressSafety); }
document.addEventListener('pointerdown', e => {
  if (e.button != null && e.button > 0) return;
  const el = e.target.closest && e.target.closest(PRESSABLE);
  endPress(); if (!el) return; const r = el.getBoundingClientRect(); el.style.setProperty('--px', (e.clientX - r.left) + 'px');
  el.classList.remove('pressing'); void el.offsetWidth; el.classList.add('pressing'); pressEl = el;
  pressSafety = setTimeout(endPress, 90000);
}, true);
document.addEventListener('pointermove', e => {
  if (!pressEl) return;
  if (!pressEl.isConnected) { endPress(); return; }
  const r = pressEl.getBoundingClientRect(), pad = 2;
  if (e.clientX < r.left - pad || e.clientX > r.right + pad || e.clientY < r.top - pad || e.clientY > r.bottom + pad) endPress();
}, { passive: true, capture: true });
document.addEventListener('pointerup', endPress, true); document.addEventListener('pointercancel', endPress, true);
window.addEventListener('blur', endPress);
/* a held mouse that moves would otherwise start a native image/SVG drag (which cancels the pointer) */
document.addEventListener('dragstart', e => { if (e.target.closest && e.target.closest(PRESSABLE) && !e.target.closest('[draggable="true"]')) e.preventDefault(); }, true);
let latt = null;
document.addEventListener('pointermove', e => {
  if (latt) return; latt = requestAnimationFrame(() => { latt = null; });
  const c = e.target.closest && e.target.closest('.card'); if (!c) return;
  const r = c.getBoundingClientRect(); c.style.setProperty('--mx', (e.clientX - r.left) + 'px'); c.style.setProperty('--my', (e.clientY - r.top) + 'px');
}, { passive: true });
const FR = $('#fring'); let frEl = null;
function placeRing(el, anim) {
  if (!el || !el.isConnected) { FR.classList.remove('on'); frEl = null; return; }
  const r = el.getBoundingClientRect(); const pad = 3; const br = parseFloat(getComputedStyle(el).borderTopLeftRadius) || 7;
  FR.style.width = (r.width + pad * 2) + 'px'; FR.style.height = (r.height + pad * 2) + 'px';
  FR.style.transform = `translate(${r.left - pad}px,${r.top - pad}px)`; FR.style.borderRadius = Math.min(br + pad, 22) + 'px';
  if (anim) { FR.classList.remove('on'); void FR.offsetWidth; FR.classList.add('on'); }
}
document.addEventListener('focusin', e => {
  const el = e.target; if (!(el instanceof Element) || el === document.body || !el.matches(':focus-visible') || el.closest('.gcan,.bcan,.bx-can,.bx-smi') || el.matches('.card,.mini')) { FR.classList.remove('on'); frEl = null; return; }
  FR.classList.toggle('still', !frEl); frEl = el; placeRing(el, true);
});
document.addEventListener('focusout', () => setTimeout(() => { const a = document.activeElement; if (!a || a === document.body) { FR.classList.remove('on'); frEl = null; } }, 0));
window.addEventListener('scroll', () => frEl && placeRing(frEl, false), true);
window.addEventListener('resize', () => { frEl && placeRing(frEl, false); inkTabs(); });

/* ---------- views & routing ---------- */
const VIEWS = {}; const MOUNTED = new Set();
const viewEl = name => { let el = $('#v-' + name); if (!el) { el = document.createElement('div'); el.id = 'v-' + name; el.className = 'v'; $('#view').append(el); } return el; };
function inkTabs() { const a = $(`#tabs a[data-v="${S.view}"]`), ink = $('#tabs .tab-ink'); if (!a || !ink) return; ink.style.left = (a.offsetLeft + 9) + 'px'; ink.style.width = (a.offsetWidth - 18) + 'px'; }
function go(name) {
  if (!VIEWS[name]) name = 'reactions';
  const prev = S.view; S.view = name;
  for (const a of $$('#tabs a')) a.classList.toggle('on', a.dataset.v === name);
  inkTabs();
  for (const k in VIEWS) { const el = $('#v-' + k); if (el) el.classList.toggle('on', k === name); }
  const el = viewEl(name); el.classList.add('on');
  if (!MOUNTED.has(name)) { MOUNTED.add(name); VIEWS[name].mount(el); }
  if (VIEWS[name].show) VIEWS[name].show(el);
  if (prev !== name && VIEWS[prev] && VIEWS[prev].hide) VIEWS[prev].hide();
  if (location.hash !== '#' + name) history.replaceState(null, '', '#' + name);
}
window.addEventListener('hashchange', () => go(location.hash.slice(1)));

/* ---------- cards ---------- */
const IO = new IntersectionObserver(es => { for (const e of es) if (e.isIntersecting) { const el = e.target; const a = el.querySelector('.art'); if (a && !a.firstChild) a.innerHTML = ART.rx[el.dataset.id] || ''; IO.unobserve(el); } }, { rootMargin: '700px 0px' });
const hydrate = root => { for (const c of root.querySelectorAll('.card,.mini')) IO.observe(c); };
function cardHTML(r) {
  const t = r.tags;
  const badges = [...(t.mt || []).slice(0, 2), ...(t.int || []).filter(x => x !== 'none').slice(0, 1)].map(x => `<span class="tagx">${esc(pv(x))}</span>`).join('');
  return `<article class="card" tabindex="0" data-id="${r.id}" style="--bar:${barCSS(r)}"><i class="bar"></i><div class="art"></div><h3>${fmt(title(r))}</h3><div class="chs">${esc(chLine(r))}</div><div class="tg">${badges}${r.ns ? `<span class="tagx">${r.ns} step${r.ns > 1 ? 's' : ''}</span>` : ''}</div></article>`;
}
function miniHTML(r, note) { return `<div class="mini" tabindex="0" data-id="${r.id}"><div class="art"></div><h4>${fmt(title(r))}</h4><p>${esc(note || chLine(r))}</p></div>`; }
function bindCards(root, idsFn) {
  const open = el => { const ids = idsFn(); openViewer(el.dataset.id, ids, { from: el.getBoundingClientRect() }); };
  root.addEventListener('click', e => { const c = e.target.closest('.card,.mini'); if (c && root.contains(c)) open(c); });
  root.addEventListener('keydown', e => { if (e.key !== 'Enter') return; const c = e.target.closest('.card,.mini'); if (c) { e.preventDefault(); open(c); } });
}

/* ---------- Reactions view ---------- */
const SEARCH_ICON = '<svg viewBox="0 0 20 20" aria-hidden="true"><circle cx="8.5" cy="8.5" r="5.5" fill="none" stroke="currentColor" stroke-width="1.7"/><path d="m13 13 4.2 4.2" stroke="currentColor" stroke-width="1.7" stroke-linecap="round"/></svg>';
const bondSw = (id, on, label) => `<label class="sw"><button class="bond${on ? ' on' : ''}" id="${id}" role="switch" aria-checked="${on}"><i class="a1"></i><i class="l1"></i><i class="l2"></i><i class="a2"></i></button><span>${label}</span></label>`;
const CHAPS = (() => { const m = new Map(); for (const r of R) for (const l of r.chs) { const k = chName(l); m.set(k, (m.get(k) || 0) + 1); } return [...m].sort((a, b) => chIdx(a[0]) - chIdx(b[0])); })();
const CHLIST = [];
/* From / To: type a group, pick from live suggestions (with how many reactions each would give) */
const FGLIST = () => D.FG.filter(f => R.some(r => r.fg.includes(f.k) || (r.fgs && (r.fgs[0].includes(f.k) || r.fgs[1].includes(f.k))))).sort((a, b) => a.label.localeCompare(b.label));
function fgCombo(id, lab) { const v = S[id]; return `<div class="fgc" data-fgc="${id}"><input class="in" id="${id}" autocomplete="off" spellcheck="false" placeholder="${lab}: any group" value="${v ? esc(FGL[v] || v) : ''}" aria-label="${lab} group" role="combobox" aria-expanded="false"><button type="button" class="fgx" data-fgx="${id}" aria-label="Clear" ${v ? '' : 'hidden'}>×</button><div class="fgl" role="listbox" hidden></div></div>`; }
function fgCount(side, k) { const other = side === 'fa' ? S.fb : S.fa; return R.filter(r => { const a = r.fgs ? r.fgs[0] : [r.fg[0]], b = r.fgs ? r.fgs[1] : [r.fg[1]]; const mine = side === 'fa' ? a : b, oth = side === 'fa' ? b : a; return mine.includes(k) && (!other || oth.includes(other)); }).length; }
function bindFgCombo(root, id, onPick) {
  const box = root.querySelector(`[data-fgc="${id}"]`); if (!box) return; const inp = box.querySelector('input'), list = box.querySelector('.fgl'), x = box.querySelector('.fgx'); let items = [], sel = -1;
  const show = () => { const q = norm(inp.value); const all = FGLIST().map(f => ({ f, n: fgCount(id, f.k) })).filter(o => o.n > 0);
    const sc = o => { const l = norm(o.f.label), k = o.f.k; if (!q) return 1; if (l.startsWith(q) || k.startsWith(q)) return 3; if (l.split(/[\s/(,-]+/).some(w => w.startsWith(q))) return 2; return l.includes(q) || k.includes(q) ? 1 : 0; };
    items = all.map(o => [sc(o), o]).filter(z => z[0] > 0).sort((a, b) => b[0] - a[0] || b[1].n - a[1].n).slice(0, 9).map(z => z[1]); sel = Math.min(sel, items.length - 1);
    list.innerHTML = items.length ? items.map((o, j) => `<button type="button" class="${j === sel ? 'on' : ''}" data-j="${j}" role="option"><span>${esc(o.f.label)}</span><em>${o.n}</em></button>`).join('') : '<div class="none">No group matches</div>';
    list.hidden = false; inp.setAttribute('aria-expanded', 'true'); };
  const pick = o => { S[id] = o ? o.f.k : ''; inp.value = o ? o.f.label : ''; x.hidden = !o; list.hidden = true; inp.setAttribute('aria-expanded', 'false'); onPick(); };
  inp.addEventListener('focus', () => { inp.select(); show(); });
  inp.addEventListener('input', () => { sel = 0; show(); });
  inp.addEventListener('keydown', e => { if (e.key === 'ArrowDown' || e.key === 'ArrowUp') { e.preventDefault(); if (list.hidden) show(); sel = (sel + (e.key === 'ArrowDown' ? 1 : -1) + items.length) % Math.max(1, items.length); show(); }
    else if (e.key === 'Enter') { e.preventDefault(); if (!inp.value.trim()) pick(null); else if (items[sel >= 0 ? sel : 0]) pick(items[sel >= 0 ? sel : 0]); inp.blur(); }
    else if (e.key === 'Escape') { list.hidden = true; inp.value = S[id] ? FGL[S[id]] : ''; inp.blur(); } });
  inp.addEventListener('blur', () => setTimeout(() => { list.hidden = true; if (!inp.value.trim() && S[id]) pick(null); else inp.value = S[id] ? FGL[S[id]] : ''; }, 150));
  list.addEventListener('pointerdown', e => { const b = e.target.closest('[data-j]'); if (!b) return; e.preventDefault(); pick(items[+b.dataset.j]); inp.blur(); });
  x.addEventListener('click', () => pick(null));
}
const FGOPT = v => D.FG.filter(f => R.some(r => r.fg.includes(f.k) || (r.fgs && (r.fgs[0].includes(f.k) || r.fgs[1].includes(f.k))))).sort((a, b) => a.label.localeCompare(b.label)).map(f => `<option value="${f.k}"${v === f.k ? ' selected' : ''}>${esc(f.label)}</option>`).join('');
const MTS = (() => { const m = new Map(); for (const r of R) for (const x of r.tags.mt || []) m.set(x, (m.get(x) || 0) + 1); return [...m].sort((a, b) => b[1] - a[1]); })();
const CLSN = (() => { const m = new Map(); for (const r of R) for (const c of clsOf(r)) m.set(c, (m.get(c) || 0) + 1); return m; })();
let LIST = [];
VIEWS.reactions = {
  mount(el) {
    el.classList.add('rv');
    el.innerHTML = `<aside class="rail" id="rail"></aside><section class="res" id="res"><div class="bar-top" id="bartop"><div class="search">${SEARCH_ICON}<input class="in" id="q" placeholder="Search names, reagents, mechanisms, groups" autocomplete="off" spellcheck="false"></div>
      <select class="in" id="grp" style="width:auto" aria-label="Group by"><option value="none">No grouping</option><option value="chapter">Group by chapter</option><option value="class">Group by class</option><option value="mt">Group by mechanism</option></select><span class="count" id="cnt"></span></div>
      <div class="pills" id="pills"></div><div class="grid" id="grid"></div></section>`;
    this.rail(); const q = $('#q', el); q.value = S.q;
    q.addEventListener('input', debounce(() => { S.q = q.value; this.list(); }, 140));
    $('#grp', el).value = S.group; $('#grp', el).addEventListener('change', e => { S.group = e.target.value; this.list(); });
    const res = $('#res', el); res.addEventListener('scroll', () => $('#bartop').classList.toggle('stuck', res.scrollTop > 4), { passive: true });
    bindCards($('#grid', el), () => LIST.map(r => r.id));
    $('#pills', el).addEventListener('click', e => { const b = e.target.closest('button[data-k]'); if (!b) return; const [k, v] = [b.dataset.k, b.dataset.v]; if (k === 'ch') S.ch.delete(v); else if (k === 'cls') S.cls.delete(v); else if (k === 'sub') { S.sub = ''; S.subHits = null; } else if (k === 'conv') { S.fa = ''; S.fb = ''; } else S[k] = k === 'named' || k === 'mech' ? false : ''; this.rail(); this.list(); });
    this.list();
  },
  rail() {
    const el = $('#rail');
    el.innerHTML = `<div><h4>Chapter</h4><div class="rows">${CHAPS.map(([c, n]) => `<button class="row${S.ch.has(c) ? ' on' : ''}" data-ch="${esc(c)}"><span class="t">${esc(c)}</span><span class="n">${n}</span></button>`).join('')}</div></div>
      <div><h4>Class</h4><div class="chips">${Object.keys(CLS).filter(k => CLSN.has(k)).map(k => `<button class="chip${S.cls.has(k) ? ' on' : ''}" data-cls="${k}" style="--k:${kcol(k)}"><i class="sw8"></i>${CLS[k][0]}<span class="n">${CLSN.get(k)}</span></button>`).join('')}</div></div>
      <div><h4>Mechanism type</h4><select class="in" id="mt"><option value="">Any</option>${MTS.map(([m, n]) => `<option value="${esc(m)}"${S.mt === m ? ' selected' : ''}>${esc(pv(m))} (${n})</option>`).join('')}</select></div>
      <div><h4>Conversion</h4><div style="display:flex;flex-direction:column;gap:6px">${fgCombo('fa', 'From')}${fgCombo('fb', 'To')}</div></div>
      <div><h4>Show only</h4><div style="display:flex;flex-direction:column;gap:10px">${bondSw('f-named', S.named, 'Named reactions')}${bondSw('f-mech', S.mech, 'With a mechanism')}</div></div>
      <div><h4>Reactant contains</h4><input class="in mono" id="sub" placeholder="SMILES or SMARTS, e.g. C=O" value="${esc(S.sub)}" spellcheck="false" autocomplete="off"><div class="faint" style="font-size:12px;margin-top:6px" id="substate">Press Enter to filter by a substructure.</div></div>
      <div><button class="btn sm" id="clr">Clear all filters</button></div>`;
    el.onclick = e => {
      const b = e.target.closest('[data-ch],[data-cls],#clr,.bond'); if (!b) return;
      if (b.dataset.ch) { S.ch.has(b.dataset.ch) ? S.ch.delete(b.dataset.ch) : S.ch.add(b.dataset.ch); b.classList.toggle('on'); }
      else if (b.dataset.cls) { S.cls.has(b.dataset.cls) ? S.cls.delete(b.dataset.cls) : S.cls.add(b.dataset.cls); b.classList.toggle('on'); }
      else if (b.id === 'clr') { Object.assign(S, { q: '', mt: '', fa: '', fb: '', chf: '', named: false, mech: false, sub: '', subHits: null }); S.ch.clear(); S.cls.clear(); $('#q').value = ''; this.rail(); }
      else if (b.classList.contains('bond')) { const on = !b.classList.contains('on'); b.classList.toggle('on', on); b.setAttribute('aria-checked', on); if (b.id === 'f-named') S.named = on; else S.mech = on; }
      this.list();
    };
    $('#mt', el).onchange = e => { S.mt = e.target.value; this.list(); };
    bindFgCombo(el, 'fa', () => this.list()); bindFgCombo(el, 'fb', () => this.list());
    const sub = $('#sub', el);
    sub.onkeydown = e => { if (e.key === 'Enter') { e.preventDefault(); this.setSub(sub.value.trim()); } };
    sub.onchange = () => this.setSub(sub.value.trim());
  },
  setSub(v) {
    if (v === S.sub) return; S.sub = v; S.subHits = null; const st = $('#substate');
    if (!v) { st.textContent = 'Press Enter to filter by a substructure.'; this.list(); return; }
    st.textContent = 'Loading the chemistry engine…'; this.list();
    rdkit().then(RD => {
      let q = null; try { q = RD.get_qmol(v); } catch (e) { }
      if (!q || !q.is_valid()) { st.textContent = 'Not a valid SMILES/SMARTS.'; S.subHits = new Set(); this.list(); return; }
      const hits = new Set();
      for (const r of R) { const L = r.rxn.split('>')[0]; let m = null; try { m = RD.get_mol(L); if (m && m.get_substruct_match(q) !== '{}') hits.add(r.id); } catch (e) { } finally { m && m.delete(); } }
      q.delete(); S.subHits = hits; st.textContent = `${hits.size} reactions have it in a reactant.`; this.list();
    }).catch(() => { st.textContent = 'Could not load the chemistry engine.'; });
  },
  list() {
    LIST = filtered(); const grid = $('#grid'); if (!grid) return;
    const pend = S.sub && !S.subHits;
    $('#cnt').textContent = pend ? 'matching structures…' : `${LIST.length} of ${R.length}`;
    const pills = [];
    for (const c of S.ch) pills.push(['ch', c, c]); for (const c of S.cls) pills.push(['cls', c, CLS[c][0]]);
    if (S.mt) pills.push(['mt', S.mt, pv(S.mt)]); if (S.fa || S.fb) pills.push(['conv', '', (FGL[S.fa] || 'Any group') + ' → ' + (FGL[S.fb] || 'any group')]); if (S.chf) pills.push(['chf', S.chf, S.chf]); if (S.named) pills.push(['named', '', 'Named']); if (S.mech) pills.push(['mech', '', 'Has mechanism']); if (S.sub) pills.push(['sub', S.sub, 'Contains ' + S.sub]);
    $('#pills').innerHTML = pills.map(([k, v, l]) => `<span class="pill">${esc(l)}<button data-k="${k}" data-v="${esc(v)}" aria-label="Remove filter">✕</button></span>`).join('');
    $('#pills').style.paddingBottom = pills.length ? '6px' : '0';
    if (!LIST.length) { grid.innerHTML = `<div class="empty" style="grid-column:1/-1"><b>Nothing matches</b>Try fewer words, or clear a filter.</div>`; return; }
    let html = '';
    if (S.group === 'none') html = LIST.map(cardHTML).join('');
    else {
      const key = r => S.group === 'chapter' ? chName(prim(r)) : S.group === 'class' ? CLS[clsOf(r)[0]][0] : pv((r.tags.mt || ['other'])[0]);
      const G = new Map(); for (const r of LIST) { const k = key(r); if (!G.has(k)) G.set(k, []); G.get(k).push(r); }
      let groups = [...G]; if (S.group === 'class' || S.group === 'mt') groups.sort((a, b) => b[1].length - a[1].length);
      LIST = groups.flatMap(g => g[1]);
      html = groups.map(([k, rs]) => `<div class="ghead"><h2>${esc(k)}</h2><span>${rs.length}</span></div>` + rs.map(cardHTML).join('')).join('');
    }
    grid.innerHTML = html; hydrate(grid);
  },
  show() { inkTabs(); }
};

/* ---------- viewer ---------- */
const VW = { ids: [], i: 0, tab: 'reaction', step: 0, el: $('#viewer'), lastFocus: null };
const FACT = [['mt', 'Mechanism'], ['int', 'Intermediate'], ['st', 'Stereochemistry'], ['rg', 'Regioselectivity'], ['cond', 'Conditions'], ['sel', 'Selectivity'], ['rate', 'Rate'], ['order', 'Order'], ['rds', 'Slow step'], ['kie', 'Isotope effect'], ['rev', 'Reversibility'], ['drive', 'Driving force'], ['rearr', 'Rearrangement'], ['ctrl', 'Control'], ['taut', 'Tautomerism'], ['hsab', 'HSAB'], ['effect', 'Effect'], ['obs', 'Observation'], ['use', 'Used for'], ['need', 'Needs'], ['only', 'Only'], ['count', 'Count'], ['name_alt', 'Also called'], ['rxn_none', 'No reaction when'], ['note', 'Remarks']];
function rgChipsHTML(r) {
  const seen = new Set(); const out = [];
  for (const l of r.rl || []) for (let part of String(l).split(/,\s+|\s+then\s+|;\s*/)) {
    part = part.replace(/^\(?(?:i{1,3}|iv|v|\d)[.)]\s*/i, '').trim(); const k = norm(part).replace(/\s+/g, ''); const e = RGI.get(k);
    if (!e || seen.has(k)) continue; seen.add(k); out.push(`<button type="button" class="chip" data-act="rg" data-k="${esc(k)}" title="${e.ids.size} reactions use ${esc(e.label)}"><span>${fmt(e.label)}</span><span class="n">${e.ids.size}</span></button>`);
  }
  return out.join('');
}
function factsHTML(r) {
  const t = r.tags; const out = [];
  for (const [k, lab] of FACT) {
    let v = t[k]; if (!v || !v.length) continue; v = v.filter(x => !(k === 'int' && x === 'none') && !(k === 'rearr' && x === 'none'));
    if (!v.length) continue;
    const txt = v.map(x => k === 'cond' ? (D.CL[x] || pv(x)) : pv(x)).join(k === 'note' ? '. ' : ', ');
    out.push(`<div class="fact"><h5>${lab}</h5><p>${fmt(txt)}</p></div>`);
  }
  const [a, b] = r.fg; out.unshift(`<div class="fact"><h5>Conversion</h5><p>${esc(FGL[a])}${a !== b ? ' → ' + esc(FGL[b]) : ''}</p></div>`);
  if (r.chs.length > 1) out.push(`<div class="fact"><h5>Also in</h5><p>${r.chs.slice(1).map(l => esc(chText(l))).join('<br>')}</p></div>`);
  return `<div class="facts">${out.join('')}</div>`;
}
function related(r) {
  const same = R.filter(x => x !== r && x.fg[0] === r.fg[0] && x.fg[1] === r.fg[1]).sort((a, b) => a._k - b._k);
  const mine = (rgKeysOf.get(r.id) || []).filter(e => e.ids.size > 1 && e.ids.size <= 14);
  const sh = new Map(); for (const e of mine) for (const id of e.ids) if (id !== r.id && !sh.has(id)) sh.set(id, e.label);
  const next = R.filter(x => x !== r && x.fg[0] === r.fg[1] && x.fg[0] !== x.fg[1]).sort((a, b) => (b.name ? 1 : 0) - (a.name ? 1 : 0) || a._k - b._k).slice(0, 18);
  const prev = R.filter(x => x !== r && x.fg[1] === r.fg[0] && x.fg[0] !== x.fg[1]).sort((a, b) => (b.name ? 1 : 0) - (a.name ? 1 : 0) || a._k - b._k).slice(0, 12);
  return { same, shared: [...sh].map(([id, lab]) => [BY.get(id), lab]), next, prev };
}
function vtabCounts(r) { return { mech: r.ns, ex: r.ex.length }; }
function openViewer(id, ids, opt = {}) {
  VW.ids = ids && ids.length ? ids.slice() : [id]; VW.i = Math.max(0, VW.ids.indexOf(id)); if (VW.ids[VW.i] !== id) { VW.ids = [id]; VW.i = 0; }
  VW.tab = opt.tab || 'reaction'; VW.step = 0; VW.prevStep = null; VW.lastFocus = document.activeElement;
  const was = !VW.el.hidden; VW.el.hidden = false; renderViewer();
  const panel = $('.vw', VW.el);
  if (opt.from && !was && panel.animate) {
    const f = panel.getBoundingClientRect(), o = opt.from;
    panel.animate([{ transform: `translate(${o.left - f.left}px,${o.top - f.top}px) scale(${Math.max(.05, o.width / f.width)},${Math.max(.05, o.height / f.height)})`, opacity: .35 }, { transform: 'none', opacity: 1 }], { duration: 360, easing: 'cubic-bezier(.2,.75,.2,1)' });
  }
}
function closeViewer() { VW.el.hidden = true; VW.el.innerHTML = ''; if (VW.lastFocus && VW.lastFocus.isConnected) VW.lastFocus.focus({ preventScroll: true }); }
function renderViewer() {
  const r = BY.get(VW.ids[VW.i]); const c = vtabCounts(r);
  const cls = clsOf(r).map(k => `<span class="ks" style="--k:${kcol(k)}"><i></i>${CLS[k] ? CLS[k][0] : k}</span>`).join('');
  const tabs = [['reaction', 'Reaction', ''], ['mech', 'Mechanism', c.mech ? c.mech : ''], ['ex', 'Examples', c.ex ? c.ex : ''], ['rel', 'Related', '']];
  const canBoard = !!(window.Board && window.Board.canOpen && window.Board.canOpen(r));
  VW.el.innerHTML = `<div class="vw" role="dialog" aria-modal="true" aria-label="${esc(title(r))}">
    <div class="vw-h"><div class="tt"><h1>${fmt(title(r))}</h1><div class="sub"><span>${esc(chLine(r))}</span>${cls}</div></div>
      ${canBoard ? '<button class="btn sm" data-act="board">Open on Board</button>' : ''}<button class="btn ghost sm" data-act="close" aria-label="Close">Close <kbd>Esc</kbd></button></div>
    <div class="vtabs" role="tablist">${tabs.map(([k, l, n], i) => `<button role="tab" data-tab="${k}" class="${VW.tab === k ? 'on' : ''}" aria-selected="${VW.tab === k}">${l}${n !== '' ? `<span class="n">${n}</span>` : ''}</button>`).join('')}<span class="ink"></span></div>
    <div class="vw-b" id="vwb"></div>
    <div class="vw-f"><button class="btn sm" data-act="prev"${VW.i ? '' : ' disabled'}>← Previous</button><span class="sp"></span><span>${VW.i + 1} of ${VW.ids.length}</span><span class="sp"></span><span class="faint">${VW.tab === 'mech' ? '← → steps · Shift ← → reactions' : '← → reactions · 1 to 4 tabs'}</span><button class="btn sm" data-act="next"${VW.i < VW.ids.length - 1 ? '' : ' disabled'}>Next →</button></div></div>`;
  renderTab(); requestAnimationFrame(inkV);
}
function inkV() { const b = $('.vtabs button.on', VW.el), ink = $('.vtabs .ink', VW.el); if (!b || !ink) return; ink.style.left = (b.offsetLeft + 6) + 'px'; ink.style.width = (b.offsetWidth - 12) + 'px'; }
function setTab(t) { VW.tab = t; if(t !== 'mech') VW.prevStep = null; for (const b of $$('.vtabs button', VW.el)) { b.classList.toggle('on', b.dataset.tab === t); b.setAttribute('aria-selected', b.dataset.tab === t); } inkV(); renderTab(); const f = $('.vw-f .faint', VW.el); if (f) f.textContent = t === 'mech' ? '← → steps · Shift ← → reactions' : '← → reactions · 1 to 4 tabs'; }
function renderTab() {
  const r = BY.get(VW.ids[VW.i]); const b = $('#vwb'); b.scrollTop = 0;
  b.classList.toggle('mode-mech', VW.tab === 'mech'); b.classList.toggle('mode-rx', VW.tab === 'reaction');
  if (VW.tab === 'reaction') {
    const svg = ART.rx[r.id] || ''; const ex = ART.ex[r.id] || []; const canBoard = !!(window.Board && window.Board.canOpen && window.Board.canOpen(r) && window.BOARDCOV && window.BOARDCOV.has(r.id)); const bmode = window.BOARDCOV && window.BOARDCOV.get(r.id);
    const rg = rgChipsHTML(r);
    b.innerHTML = `<div class="rxv"><div class="rxv-main"><div class="rxv-stage" style="--nw:${Math.round(Math.min(natW(svg) * 1.7, 1000))}px">${svg}</div>
      ${ex.length ? `<section class="rxv-ex"><h3>From class<span>${ex.length}</span><button class="lnk" data-tab="ex">${ex.length > 3 ? 'All ' + ex.length : 'Open'} ›</button></h3><div class="rxv-exs">${ex.slice(0, 3).map(s => `<button type="button" class="rxv-exi" data-tab="ex">${s}</button>`).join('')}</div></section>` : ''}</div>
      <aside class="rxv-side">${factsHTML(r)}${rg ? `<div class="fact"><h5>Reagents</h5><div class="rxv-rg">${rg}</div></div>` : ''}
      <div class="rxv-acts">${r.ns ? `<button type="button" class="btn sm" data-tab="mech">Watch the mechanism<span class="faint">${r.ns} step${r.ns > 1 ? 's' : ''}</span></button>` : ''}${canBoard ? `<button type="button" class="btn sm" data-act="board">Try it on the Board<span class="faint">${bmode === 'e' ? 'engine reproduces it' : bmode === 'd' ? 'the library differs: see why' : 'product from the library'}</span></button>` : ''}</div></aside></div>`;
  }
  else if (VW.tab === 'mech') renderStep();
  else if (VW.tab === 'ex') {
    const ex = ART.ex[r.id] || [];
    b.innerHTML = ex.length ? `<div style="display:flex;align-items:center;gap:10px;margin-bottom:12px"><span class="mute">Class questions recorded for this reaction.</span><span class="sp" style="flex:1"></span>${window.Study ? '<button class="btn sm" data-act="quizex">Quiz me on these</button>' : ''}</div><div class="exg">${ex.map(s => `<div class="exi">${s}</div>`).join('')}</div>` : '<div class="empty"><b>No examples yet</b>None were recorded for this reaction.</div>';
  } else {
    const rel = related(r);
    const sec = (h, list, note) => list.length ? `<section><h3>${h}<span>${list.length}</span></h3><div class="minis">${list.map(x => Array.isArray(x) ? miniHTML(x[0], 'shares ' + x[1]) : miniHTML(x, note && note(x))).join('')}</div></section>` : '';
    b.innerHTML = `<div class="rel">${sec('Same conversion', rel.same)}${sec('What can come next', rel.next, x => `${FGL[x.fg[0]]} → ${FGL[x.fg[1]]}`)}${sec('What leads here', rel.prev, x => `${FGL[x.fg[0]]} → ${FGL[x.fg[1]]}`)}${sec('Shares a reagent', rel.shared.slice(0, 18))}</div>` || '';
    if (!b.textContent.trim()) b.innerHTML = '<div class="empty"><b>Nothing related</b></div>';
    const ids = $$('.mini', b).map(m => m.dataset.id); hydrate(b);
    b.onclick = e => { const m = e.target.closest('.mini'); if (m) { VW.ids = [...new Set([r.id, ...ids])]; VW.i = VW.ids.indexOf(m.dataset.id); VW.tab = 'reaction'; renderViewer(); } };
    return;
  }
  b.onclick = null;
}
function renderStep() {
  const r = BY.get(VW.ids[VW.i]); const b = $('#vwb'); const n = r.ns; const cam = localStorage.getItem('ocr3-camera') || 'slide';
  if (!n) { b.innerHTML = '<div class="empty"><b>No mechanism for this one</b>It is a fact to remember rather than a mechanism to push.</div>'; return; }
  const k = Math.min(VW.step, n - 1); const steps = ART.ms[r.id] || []; const svg = steps[k] || '';
  const oldK = VW.prevStep != null ? Math.min(VW.prevStep, n - 1) : null;
  const nav = `<div class="mnav"><span class="stepno">Step ${k + 1} of ${n}</span><div class="dots">${Array.from({ length: n }, (_, i) => `<button data-step="${i}" class="${i === k ? 'on' : ''}${i < k ? ' past' : ''}" aria-label="Step ${i + 1}"></button>`).join('')}</div>
    <button class="btn sm" data-act="replay">Replay</button><button class="btn sm" data-act="sprev"${k ? '' : ' disabled'}>←</button><button class="btn sm" data-act="snext"${k < n - 1 ? '' : ' disabled'}>→</button></div>`;
  let st = $('.mstage', b);
  const moving = cam === 'slide' && oldK != null && oldK !== k && st && SET.motion() && st.querySelector('.mv-layer.cur');
  if (!moving) {
    b.innerHTML = `${r.std ? '<div class="banner">Standard mechanism, not traced from the library.</div>' : ''}<div class="mstage"><div class="mv-layer cur">${svg.replace("class='step'", "class='step anim'")}</div></div>${nav}`;
    return;
  }
  b.querySelector('.mnav').outerHTML = nav;
  for (const x of [...st.querySelectorAll('.mv-layer:not(.cur)')]) x.remove();
  const old = st.querySelector('.mv-layer.cur'); old.classList.remove('cur'); old.classList.add('old');
  old.getAnimations({ subtree: true }).forEach(a => a.finish());
  const nw = document.createElement('div'); nw.className = 'mv-layer cur entering'; nw.innerHTML = svg; st.append(nw);
  const fwd = k > oldK;
  const oSel = fwd ? '.ps' : '.ls', nSel = fwd ? '.ls' : '.ps';
  const oEl = old.querySelector(oSel), nEl = nw.querySelector(nSel);
  const sr = st.getBoundingClientRect(); const oB = oEl && oEl.getBBox ? box(oEl) : null, nB = nEl ? box(nEl) : null;
  function box(el) { const q = el.getBoundingClientRect(); return { x: q.left - sr.left + q.width / 2, y: q.top - sr.top + q.height / 2, w: q.width, h: q.height }; }
  const ease = 'cubic-bezier(.45,.05,.2,1)', D = 620;
  if (oB && nB && nB.w > 4 && oB.w > 4) {
    const sc = Math.max(.55, Math.min(1.8, (oB.h || 1) / (nB.h || 1)));
    const dx = oB.x - nB.x, dy = oB.y - nB.y;
    nw.style.transformOrigin = `${nB.x}px ${nB.y}px`; old.style.transformOrigin = `${oB.x}px ${oB.y}px`;
    nw.animate([{ transform: `translate(${dx}px,${dy}px) scale(${sc})` }, { transform: 'none' }], { duration: D, easing: ease });
    old.animate([{ transform: 'none', opacity: 1 }, { transform: `translate(${-dx}px,${-dy}px) scale(${1 / sc})`, opacity: 0 }], { duration: D, easing: ease, fill: 'forwards' });
    // the shared species stays solid; everything else in the old frame fades out, the new frame's other parts fade in
    old.querySelectorAll(`.step>:not(${oSel})`).forEach(el => el.animate([{ opacity: 1 }, { opacity: 0 }], { duration: 240, fill: 'forwards' }));
    oEl.animate([{ opacity: 1 }, { opacity: 1, offset: .55 }, { opacity: 0 }], { duration: D, fill: 'forwards' });
    nEl.animate([{ opacity: 0 }, { opacity: 0, offset: .4 }, { opacity: 1 }], { duration: D, easing: ease });
    nw.querySelectorAll(`.step>:not(${nSel})`).forEach(el => el.animate([{ opacity: 0 }, { opacity: 0, offset: .6 }, { opacity: 1 }], { duration: D + 120, easing: 'ease-out' }));
  } else {
    nw.animate([{ opacity: 0, transform: `translateX(${fwd ? 40 : -40}px)` }, { opacity: 1, transform: 'none' }], { duration: 380, easing: ease });
    old.animate([{ opacity: 1 }, { opacity: 0 }], { duration: 240, fill: 'forwards' });
  }
  clearTimeout(VW.mvT);
  VW.mvT = setTimeout(() => { old.remove(); nw.classList.remove('entering'); const sv = nw.querySelector('svg.step'); if (sv) sv.classList.add('anim'); }, D - 40);
}
function stepGo(d) { const r = BY.get(VW.ids[VW.i]); if (!r.ns) return; const k = Math.max(0, Math.min(r.ns - 1, VW.step + d)); if (k !== VW.step) { VW.prevStep = VW.step; VW.step = k; renderStep(); VW.prevStep = null; } }
function vwGo(d) { const i = VW.i + d; if (i < 0 || i >= VW.ids.length) return; VW.i = i; VW.step = 0; renderViewer(); }
VW.el.addEventListener('click', e => {
  if (e.target === VW.el) { closeViewer(); return; }
  const t = e.target.closest('[data-tab]'); if (t) { setTab(t.dataset.tab); return; }
  const s = e.target.closest('[data-step]'); if (s) { VW.prevStep = VW.step; VW.step = +s.dataset.step; renderStep(); VW.prevStep = null; return; }
  const a = e.target.closest('[data-act]'); if (!a) return; const act = a.dataset.act; const r = BY.get(VW.ids[VW.i]);
  if (act === 'close') closeViewer(); else if (act === 'prev') vwGo(-1); else if (act === 'next') vwGo(1);
  else if (act === 'sprev') stepGo(-1); else if (act === 'snext') stepGo(1); else if (act === 'replay') { VW.prevStep = null; const st = $('.mstage', VW.el); if (st) st.remove(); renderStep(); }
  else if (act === 'board') { closeViewer(); go('board'); window.Board.openReaction(r); }
  else if (act === 'quizex') { closeViewer(); go('study'); window.Study.quizReaction(r.id); }
  else if (act === 'rg') { const k = a.dataset.k; closeViewer(); go('reagents'); VIEWS.reagents.select(k); }
});

/* ---------- Reagents view ---------- */
VIEWS.reagents = {
  mount(el) {
    el.classList.add('gv');
    el.innerHTML = `<aside class="glist"><div class="hd"><div class="search">${SEARCH_ICON}<input class="in" id="rgq" placeholder="Find a reagent" autocomplete="off" spellcheck="false"></div></div><div class="rows" id="rgl"></div></aside><section class="gmain" id="rgm"></section>`;
    const q = $('#rgq', el); q.addEventListener('input', debounce(() => this.rows(q.value), 100));
    $('#rgl', el).addEventListener('click', e => { const b = e.target.closest('[data-rg]'); if (!b) return; S.rg = b.dataset.rg; for (const x of $$('#rgl .row')) x.classList.toggle('on', x === b); this.main(); });
    bindCards($('#rgm', el), () => $$('#rgm .card').map(c => c.dataset.id));
    this.rows(''); this.main();
  },
  rows(q) { const nq = norm(q).replace(/\s+/g, ''); const list = RGL.filter(e => !nq || e.k.includes(nq)); $('#rgl').innerHTML = list.map(e => `<button class="row${e.k === S.rg ? ' on' : ''}" data-rg="${esc(e.k)}"><span class="t">${fmt(e.label)}</span><span class="n">${e.ids.size}</span></button>`).join('') || '<div class="empty">No reagent by that name.</div>'; },
  main() {
    const e = RGI.get(S.rg); const m = $('#rgm'); if (!e) { m.innerHTML = ''; return; }
    const main = [...e.main].map(id => BY.get(id)).sort((a, b) => a._k - b._k), exOnly = [...e.ids].filter(id => !e.main.has(id)).map(id => BY.get(id)).sort((a, b) => a._k - b._k);
    m.innerHTML = `<div class="ghero"><h1>${fmt(e.label)}</h1><span>${e.ids.size} reaction${e.ids.size > 1 ? 's' : ''}</span></div><div class="grid">${main.map(cardHTML).join('')}${exOnly.length ? `<div class="ghead"><h2>Also used in class questions</h2><span>${exOnly.length}</span></div>` + exOnly.map(cardHTML).join('') : ''}</div>`;
    m.scrollTop = 0; hydrate(m);
  },
  select(k) { S.rg = k; if (MOUNTED.has('reagents')) { this.rows($('#rgq').value = ''); this.main(); const b = $(`#rgl [data-rg="${CSS.escape(k)}"]`); b && b.scrollIntoView({ block: 'center' }); } }
};

/* ---------- command palette ---------- */
const CP = { el: $('#cmdp'), items: [], i: 0 };
const FGSYN = { alkane: 'alkane|alkanes|paraffin', alkene: 'alkene|alkenes|olefin', alkyne: 'alkyne|alkynes|acetylene', alcohol: 'alcohol|alcohols', ether: 'ether|ethers', phenol: 'phenol|phenols', aldehyde: 'aldehyde|aldehydes', ketone: 'ketone|ketones', acid: 'acid|acids|carboxylic acid|carboxylic', ester: 'ester|esters|lactone', amide: 'amide|amides|lactam', amine: 'amine|amines|aniline', nitrile: 'nitrile|nitriles|cyanide|alkyl cyanide', alkylhalide: 'alkyl halide|alkyl halides|haloalkane|halide|halides|alkyl chloride|alkyl bromide', arylhalide: 'aryl halide|haloarene|halobenzene|chlorobenzene|bromobenzene', arene: 'arene|benzene|aromatic|arenes', diazonium: 'diazonium|diazonium salt', nitroarene: 'nitrobenzene|nitroarene', nitro: 'nitro|nitroalkane', organometallic: 'grignard|organometallic|grignard reagent', epoxide: 'epoxide|oxirane', imine: 'imine|oxime|hydrazone|schiff base', acylhalide: 'acyl chloride|acid chloride|acyl halide', anhydride: 'anhydride|anhydrides', azo: 'azo|azo dye', isocyanide: 'isocyanide|carbylamine', enamine: 'enamine', thiol: 'thiol|thiols|sulphide|sulfide', sulfonic: 'sulphonic acid|sulfonic acid', borane: 'borane|organoborane', isocyanate: 'isocyanate', carbocation: 'carbocation', sulfonate: 'tosylate|sulfonate' };
function fgOf(s) { s = s.trim().toLowerCase().replace(/s$/, ''); let best = null; for (const k in FGSYN) for (const w of FGSYN[k].split('|')) { const ww = w.replace(/s$/, ''); if (ww === s) return k; if (!best && s.length >= 3 && ww.startsWith(s)) best = k; } return best; }
function routes(A, B) {
  const direct = R.filter(r => r.fg[0] === A && r.fg[1] === B);
  const two = [];
  for (const r1 of R) { if (r1.fg[0] !== A || r1.fg[1] === A || r1.fg[1] === B) continue; for (const r2 of R) if (r2.fg[0] === r1.fg[1] && r2.fg[1] === B) two.push([r1, r2]); }
  two.sort((p, q) => ((q[0].name ? 1 : 0) + (q[1].name ? 1 : 0)) - ((p[0].name ? 1 : 0) + (p[1].name ? 1 : 0)));
  const seen = new Set(); const uniq = []; for (const p of two) { const k = title(p[0]) + '|' + title(p[1]); if (!seen.has(k)) { seen.add(k); uniq.push(p); } }
  return { direct, two: uniq.slice(0, 8) };
}
function cpItems(q) {
  const out = []; const nq = norm(q).trim();
  const add = (sec, ic, lb, ht, run) => out.push({ sec, ic, lb, ht, run });
  const m = /^(.+?)\s*(?:\bto\b|→|->|=>|>)\s*(.+)$/i.exec(q.trim());
  if (m) {
    const A = fgOf(m[1]), B = fgOf(m[2]);
    if (A && B) {
      const { direct, two } = routes(A, B);
      for (const r of direct.slice(0, 8)) add(`${FGL[A]} to ${FGL[B]}`, '→', fmt(title(r)), esc(chLine(r)), () => openViewer(r.id, direct.map(x => x.id)));
      for (const [r1, r2] of two) add('Two-step routes', '⇢', `<span class="route"><b>${esc(FGL[A])}</b><i>→</i>${fmt(title(r1))}<i>→</i><b>${esc(FGL[r1.fg[1]])}</b><i>→</i>${fmt(title(r2))}<i>→</i><b>${esc(FGL[B])}</b></span>`, '', () => openViewer(r1.id, [r1.id, r2.id]));
      if (!direct.length && !two.length) add('Routes', '·', `No route from ${esc(FGL[A])} to ${esc(FGL[B])} in the library`, '', () => { });
    }
  }
  const VIEWLIST = [['Reactions', 'reactions'], ['Reagents', 'reagents'], ['Graph', 'graph'], ['Study', 'study'], ['Board', 'board']];
  const modes = window.Study ? window.Study.MODES : [];
  if (!nq) {
    for (const [l, v] of VIEWLIST) add('Go to', '↗', l, '', () => go(v));
    for (const md of modes) add('Study', '◆', md.label, md.hint, () => { go('study'); window.Study.setMode(md.k); });
    add('Try', '?', 'Type “alkene to alcohol” for routes', '', () => { CP.inp.value = 'alkene to alcohol'; cpRender(); });
    add('Actions', '◐', 'Switch theme', '', toggleTheme);
    return out;
  }
  const toks = nq.split(/\s+/);
  const score = (s) => { const n = norm(s); if (n.startsWith(nq)) return 3; if (toks.every(t => n.includes(t))) return 2; return 0; };
  const rs = R.map(r => [r, Math.max(score(title(r)) + (r.name ? .5 : 0), toks.every(t => hay(r).includes(t)) ? 1 : 0)]).filter(x => x[1] > 0).sort((a, b) => b[1] - a[1] || a[0]._k - b[0]._k);
  const rids = rs.map(x => x[0].id);
  for (const [r] of rs.slice(0, 8)) add('Reactions', '⌬', fmt(title(r)), esc(chLine(r)), () => openViewer(r.id, rids));
  if (rs.length > 8) add('Reactions', '⋯', `Show all ${rs.length} matches in the list`, '', () => { S.q = q; go('reactions'); const i = $('#q'); if (i) { i.value = q; VIEWS.reactions.list(); } });
  for (const e of RGL.filter(e => e.k.includes(nq.replace(/\s+/g, ''))).slice(0, 5)) add('Reagents', '⚗', fmt(e.label), `${e.ids.size} reactions`, () => { go('reagents'); VIEWS.reagents.select(e.k); });
  for (const [c] of CHAPS.filter(([c]) => toks.every(t => norm(c).includes(t))).slice(0, 3)) add('Filter', '▦', 'Chapter: ' + esc(c), '', () => { S.ch.add(c); go('reactions'); VIEWS.reactions.rail(); VIEWS.reactions.list(); });
  for (const k of Object.keys(CLS).filter(k => CLSN.has(k) && norm(CLS[k][0]).startsWith(nq)).slice(0, 2)) add('Filter', '▦', 'Class: ' + CLS[k][0], '', () => { S.cls.add(k); go('reactions'); VIEWS.reactions.rail(); VIEWS.reactions.list(); });
  for (const [mt] of MTS.filter(([mt]) => norm(pv(mt)).startsWith(nq)).slice(0, 3)) add('Filter', '▦', 'Mechanism: ' + esc(pv(mt)), '', () => { S.mt = mt; go('reactions'); VIEWS.reactions.rail(); VIEWS.reactions.list(); });
  for (const [l, v] of VIEWLIST.filter(([l]) => norm(l).startsWith(nq))) add('Go to', '↗', l, '', () => go(v));
  for (const md of modes.filter(md => toks.every(t => norm(md.label + ' ' + md.hint).includes(t)))) add('Study', '◆', md.label, md.hint, () => { go('study'); window.Study.setMode(md.k); });
  if ('theme dark light'.includes(nq)) add('Actions', '◐', 'Switch theme', '', toggleTheme);
  return out;
}
function cpRender() {
  CP.items = cpItems(CP.inp.value); CP.i = 0; let sec = '', html = '';
  CP.items.forEach((it, i) => { if (it.sec !== sec) { sec = it.sec; html += `<div class="sec">${esc(sec)}</div>`; } html += `<div class="pi${i === 0 ? ' on' : ''}" data-i="${i}" role="option"><span class="ic">${it.ic}</span><span class="lb">${it.lb}</span>${it.ht ? `<span class="ht">${it.ht}</span>` : ''}</div>`; });
  CP.lst.innerHTML = html || '<div class="empty">No matches.</div>';
}
function cpMove(d) { if (!CP.items.length) return; CP.i = (CP.i + d + CP.items.length) % CP.items.length; $$('.pi', CP.lst).forEach((p, i) => p.classList.toggle('on', i === CP.i)); const on = $('.pi.on', CP.lst); on && on.scrollIntoView({ block: 'nearest' }); }
function cpRun(i) { const it = CP.items[i]; if (!it) return; closeCP(); it.run(); }
function openCP(q = '') {
  CP.el.hidden = false; CP.el.innerHTML = `<div class="cp" role="dialog" aria-label="Command palette"><input id="cpi" placeholder="Search reactions, reagents, chapters… or “alkene to alcohol”" autocomplete="off" spellcheck="false"><div class="lst" id="cpl" role="listbox"></div><div class="ft"><span><kbd>↑</kbd> <kbd>↓</kbd> move</span><span><kbd>Enter</kbd> open</span><span><kbd>Esc</kbd> close</span></div></div>`;
  CP.inp = $('#cpi'); CP.lst = $('#cpl'); CP.inp.value = q; cpRender(); CP.inp.focus();
  CP.inp.addEventListener('input', cpRender);
  CP.inp.addEventListener('keydown', e => { if (e.key === 'ArrowDown') { e.preventDefault(); cpMove(1); } else if (e.key === 'ArrowUp') { e.preventDefault(); cpMove(-1); } else if (e.key === 'Enter') { e.preventDefault(); cpRun(CP.i); } else if (e.key === 'Escape') { e.preventDefault(); closeCP(); } });
  CP.lst.addEventListener('click', e => { const p = e.target.closest('.pi'); if (p) cpRun(+p.dataset.i); });
  CP.lst.addEventListener('pointermove', e => { const p = e.target.closest('.pi'); if (p && +p.dataset.i !== CP.i) { CP.i = +p.dataset.i; $$('.pi', CP.lst).forEach((x, i) => x.classList.toggle('on', i === CP.i)); } });
}
function closeCP() { CP.el.hidden = true; CP.el.innerHTML = ''; }
CP.el.addEventListener('click', e => { if (e.target === CP.el) closeCP(); });
$('#kbar').addEventListener('click', () => openCP());

/* ---------- Settings page: custom controls, not a modal ---------- */
const THEMES = [['dark','Slate','tinted graphite + chalk','oklch(.235 .013 175)'],['graphite','Graphite','neutral grey, no tint','oklch(.235 .004 260)'],['oled','OLED Black','true black for dark rooms','oklch(.09 0 0)'],['midnight','Midnight','blue-black + electric pink','oklch(.17 .035 260)'],['dune','Dune','warm umber + sand','oklch(.235 .02 62)'],['plum','Plum','deep violet dusk','oklch(.235 .026 318)'],['forest','Forest','pine green + moss','oklch(.235 .027 150)'],['light','Paper','warm paper + ink','oklch(.95 .01 95)'],['snow','Snow','cool white + slate ink','oklch(.965 .005 250)']];
const LIGHT_THEMES = new Set(['light','snow']);
const ACCENTS = [['auto','Theme',85],['amber','Amber',85],['lime','Lime',128],['mint','Mint',165],['cyan','Cyan',205],['blue','Blue',255],['violet','Violet',295],['pink','Pink',345],['coral','Coral',28]];
function applyAccent(k) { const r = document.documentElement.style; const a = ACCENTS.find(x => x[0] === k); const vars = ['--acc','--acct','--accbg','--onacc','--latt'];
  if (!a || k === 'auto') { vars.forEach(v => r.removeProperty(v)); return; } const h = a[2]; const lt = LIGHT_THEMES.has(document.documentElement.dataset.theme);
  if (lt) { r.setProperty('--acc', `oklch(.68 .16 ${h})`); r.setProperty('--acct', `oklch(.47 .13 ${h})`); r.setProperty('--accbg', `oklch(.74 .15 ${h} / .22)`); r.setProperty('--onacc', `oklch(.99 .01 ${h})`); r.setProperty('--latt', `oklch(.6 .13 ${h} / .4)`); }
  else { r.setProperty('--acc', `oklch(.82 .15 ${h})`); r.setProperty('--acct', `oklch(.86 .13 ${h})`); r.setProperty('--accbg', `oklch(.82 .15 ${h} / .14)`); r.setProperty('--onacc', `oklch(.22 .05 ${h})`); r.setProperty('--latt', `oklch(.82 .15 ${h} / .32)`); } }
const SET = { accent: () => localStorage.getItem('ocr3-accent') || 'auto', theme: () => localStorage.getItem('ocr3-theme') || 'dark', camera: () => localStorage.getItem('ocr3-camera') || 'slide', repel: () => localStorage.getItem('ocr3-repel') || 'subtle', motion: () => localStorage.getItem('ocr3-motion') !== 'reduced' };
function customDrop(id, label, value, options) { return `<div class="setting-control"><span class="setting-label">${label}</span><div class="drop" data-drop="${id}"><button class="drop-trigger" type="button"><span>${options.find(x=>x[0]===value)?.[1]||value}</span><b>⌄</b></button><div class="drop-menu">${options.map(x=>`<button type="button" data-drop-value="${x[0]}" class="${x[0]===value?'on':''}"><span>${x[1]}</span><small>${x[2]||''}</small></button>`).join('')}</div></div></div>`; }
function saveSetting(k,v){ localStorage.setItem('ocr3-'+k,v); if(k==='motion') document.documentElement.dataset.motion=v; if(k==='theme') setTheme(v); if(k==='accent') applyAccent(v); if(k==='camera') toast(v==='slide'?'Camera movement enabled':'Instant step changes enabled'); if(k==='motion') window.dispatchEvent(new Event('motionchange')); }
VIEWS.settings = { mount(el){ el.classList.add('settings-page'); el.innerHTML=`<div class="settings-shell"><div class="settings-intro"><div class="eyebrow">OC Reactions / Preferences</div><h1>Make it yours.</h1><p>Quiet controls for the things you stare at repeatedly: colour, motion, and how the chemistry moves.</p></div><section class="settings-section"><div class="section-head"><h2>Appearance</h2><span>the surface of the app</span></div><div class="theme-grid">${THEMES.map(x=>`<button class="theme-card ${SET.theme()===x[0]?'on':''}" data-theme-choice="${x[0]}"><i style="background:${x[3]}"></i><b>${x[1]}</b><small>${x[2]}</small></button>`).join('')}</div><div class="setting-control accent-wrap"><span class="setting-label">Accent colour</span><div class="accent-row">${ACCENTS.map(([k,n,h])=>`<button type="button" class="accent-dot${k==='auto'?' auto':''}${SET.accent()===k?' on':''}" data-accent="${k}" style="--sw:oklch(.78 .15 ${h})" title="${n}" aria-label="${n} accent"><span>${n}</span></button>`).join('')}</div></div></section><section class="settings-section"><div class="section-head"><h2>Movement</h2><span>continuity over novelty</span></div>${customDrop('camera','Mechanism next step',SET.camera(),[['slide','Move the camera','bring the next step in from the right'],['cut','Instant change','no transition, maximum speed']])}${customDrop('repel','Graph hover behaviour',SET.repel(),[['subtle','Subtle drift','nearby nodes yield a little'],['push','Pointer repulsion','slow pointer movement pushes nodes away'],['still','Stationary','hover labels only']])}${customDrop('motion','Animation amount',SET.motion()?'full':'reduced',[['full','Full motion','all signature transitions'],['reduced','Reduced motion','respectful, minimal movement']])}</section><section class="settings-section"><div class="section-head"><h2>Board</h2><span>the workbench</span></div><div class="setting-note"><b>Live solving is always on.</b><span>Type into an arrow and the product updates automatically. Use the Board tab for the full editor.</span></div></section><div class="settings-foot">Settings are saved in this browser only.<button class="btn sm ghost" data-reset-settings>Reset preferences</button></div></div>`; el.onclick=e=>{const th=e.target.closest('[data-theme-choice]');if(th){saveSetting('theme',th.dataset.themeChoice);$$('[data-theme-choice]',el).forEach(x=>x.classList.toggle('on',x===th));}const ac=e.target.closest('[data-accent]');if(ac){saveSetting('accent',ac.dataset.accent);$$('[data-accent]',el).forEach(x=>x.classList.toggle('on',x===ac));toast('Accent: '+ac.title);}const d=e.target.closest('[data-drop]');if(d&&e.target.closest('.drop-trigger')){e.stopPropagation();$$('.drop.open',el).forEach(x=>x!==d&&x.classList.remove('open'));d.classList.toggle('open');}const opt=e.target.closest('[data-drop-value]');if(opt){const drop=opt.closest('.drop'),k=drop.dataset.drop;saveSetting(k,opt.dataset.dropValue);drop.classList.remove('open');drop.querySelector('.drop-trigger span').textContent=opt.querySelector('span').textContent;$$('[data-drop-value]',drop).forEach(x=>x.classList.toggle('on',x===opt));}if(e.target.closest('[data-reset-settings]')){['theme','camera','repel','motion','accent'].forEach(k=>localStorage.removeItem('ocr3-'+k));setTheme('dark');applyAccent('auto');this.mount(el);}}; }, show(el){ $$('.theme-card',el).forEach(x=>x.classList.toggle('on',x.dataset.themeChoice===SET.theme())); }};
$('#settings').addEventListener('click',()=>go('settings'));
/* ---------- theme ---------- */
function setTheme(t) { document.documentElement.dataset.theme = t; try { if (LIGHT_THEMES.has(t)) localStorage.setItem('ocr3-lastlight', t); else localStorage.setItem('ocr3-lastdark', t); } catch (e) { } try { localStorage.setItem('ocr3-theme', t); } catch (e) { } const b = $('#theme'); b.classList.toggle('on', !LIGHT_THEMES.has(t)); b.setAttribute('aria-checked', !LIGHT_THEMES.has(t)); $('.theme-sw span').textContent = (THEMES.find(x => x[0] === t) || THEMES[0])[1]; applyAccent(SET.accent()); window.dispatchEvent(new Event('themechange')); }
function toggleTheme() { const t = document.documentElement.dataset.theme; if (LIGHT_THEMES.has(t)) setTheme(localStorage.getItem('ocr3-lastdark') || 'dark'); else { try { localStorage.setItem('ocr3-lastdark', t); } catch (e) { } setTheme(localStorage.getItem('ocr3-lastlight') || 'light'); } }
$('#theme').addEventListener('click', toggleTheme);

/* ---------- keyboard ---------- */
document.addEventListener('keydown', e => {
  const inField = e.target.matches && e.target.matches('input,textarea,select,[contenteditable]');
  if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === 'k') { e.preventDefault(); CP.el.hidden ? openCP() : closeCP(); return; }
  if (!CP.el.hidden) return;
  if (!VW.el.hidden) {
    if (e.key === 'Escape') { e.preventDefault(); closeViewer(); return; }
    if (inField) return;
    if (e.key === 'ArrowRight' || e.key === 'ArrowLeft') { e.preventDefault(); const d = e.key === 'ArrowRight' ? 1 : -1; if (VW.tab === 'mech' && !e.shiftKey) stepGo(d); else vwGo(d); return; }
    if (/^[1-4]$/.test(e.key)) { setTab(['reaction', 'mech', 'ex', 'rel'][+e.key - 1]); return; }
    return;
  }
  if (inField) return;
  if (e.key === '/' && S.view === 'reactions') { e.preventDefault(); $('#q') && $('#q').focus(); }
});

window.App = { D, R, BY, ART, S, $, $$, esc, fmt, norm, title, chLine, chText, chName, prim, clsOf, kcol, CLS, FGL, pv, RGI, RGL, rgKeysOf, VIEWS, MOUNTED, go, openViewer, toast, rdkit, filtered, cardHTML, miniHTML, hydrate, bindCards, natW, debounce, bondSw, inkTabs };
setTheme(document.documentElement.dataset.theme || 'dark'); document.documentElement.dataset.motion = localStorage.getItem('ocr3-motion') === 'reduced' ? 'reduced' : 'full';
window.addEventListener('DOMContentLoaded', () => { go((location.hash || '#reactions').slice(1)); if (document.fonts) document.fonts.ready.then(inkTabs); });
})();
