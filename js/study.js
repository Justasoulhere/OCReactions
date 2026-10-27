/* study.js · spaced repetition over class questions, products, reagents, names and mechanism steps */
(() => {
'use strict';
const A = window.App; if (!A) return;
const { R, BY, ART, $, $$, esc, fmt, title, chLine, natW } = A;
const MODES = [
  { k: 'class', label: 'Class questions', hint: 'predict the product, from your class examples' },
  { k: 'product', label: 'Predict the product', hint: 'the main reaction of each card' },
  { k: 'reagent', label: 'Recall the reagents', hint: 'what goes over the arrow' },
  { k: 'name', label: 'Name the reaction', hint: 'named reactions only' },
  { k: 'step', label: 'Next mechanism step', hint: 'what the arrows make' }];
const KEY = 'ocr3-srs', IV = [0, 1, 3, 7, 16, 35, 80];
let SRS = {}; try { SRS = JSON.parse(localStorage.getItem(KEY) || '{}'); } catch (e) { SRS = {}; }
const save = () => { try { localStorage.setItem(KEY, JSON.stringify(SRS)); } catch (e) { } };
const T = { mode: 'class', scope: 'all', cur: null, shown: false, extra: 0, ses: [], t0: 0 };
const ivl = ms => ms < 36e5 ? Math.max(1, Math.round(ms / 6e4)) + 'm' : ms < 864e5 ? Math.round(ms / 36e5) + 'h' : Math.round(ms / 864e5) + 'd';
function preview(it, g) { const s = SRS[it.key] || { b: 0 }; if (g === 1) return 60e3; if (g === 2) return 10 * 60e3; const b = Math.min(6, s.b + (g === 3 ? 1 : 2)); return IV[b] * 864e5; }
function pool() {
  let rs = R;
  if (T.scope === 'filters') rs = A.filtered();
  else if (T.scope.startsWith('r:')) rs = [BY.get(T.scope.slice(2))].filter(Boolean);
  else if (T.scope.startsWith('ch:')) rs = R.filter(r => r.chs.some(l => A.chName(l) === T.scope.slice(3)));
  const out = [];
  for (const r of rs) {
    if (T.mode === 'class') (ART.ex[r.id] || []).forEach((s, k) => s && out.push({ key: `ex|${r.id}|${k}`, r, k }));
    else if (T.mode === 'product') { if (r.rxn.split('>')[2]) out.push({ key: `rx|${r.id}`, r }); }
    else if (T.mode === 'reagent') { if (r.rl.length || r.cl.length) out.push({ key: `rg|${r.id}`, r }); }
    else if (T.mode === 'name') { if (r.name) out.push({ key: `nm|${r.id}`, r }); }
    else for (let k = 0; k < r.ns; k++) out.push({ key: `st|${r.id}|${k}`, r, k });
  }
  return out;
}
function hash(s) { let h = 2166136261; for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 16777619); } return h >>> 0; }
function pick() {
  const now = Date.now(), items = pool(); const day = new Date().toDateString();
  const due = items.filter(it => SRS[it.key] && SRS[it.key].due <= now).sort((a, b) => SRS[a.key].due - SRS[b.key].due);
  const fresh = items.filter(it => !SRS[it.key]).sort((a, b) => hash(day + a.key) - hash(day + b.key));
  const learned = items.filter(it => SRS[it.key] && SRS[it.key].b >= 3).length;
  T.stats = { due: due.length, fresh: fresh.length, learned, total: items.length };
  if (due.length) return due[0]; if (fresh.length) return fresh[0];
  if (T.extra && items.length) return items[hash(day + T.extra) % items.length];
  return null;
}
function grade(g) {
  const it = T.cur; if (!it) return; const s = SRS[it.key] || { b: 0, n: 0 }; const now = Date.now();
  if (g === 1) { s.b = 0; s.due = now + 60e3; } else if (g === 2) { s.due = now + 10 * 60e3; } else if (g === 3) { s.b = Math.min(6, s.b + 1); s.due = now + IV[s.b] * 864e5; } else { s.b = Math.min(6, s.b + 2); s.due = now + IV[s.b] * 864e5; }
  s.n = (s.n || 0) + 1; s.last = now; SRS[it.key] = s; save(); if (T.extra) T.extra++;
  T.ses.push({ g, ms: now - (T.t0 || now) }); next();
}
let EL = null;
function chapters() { const m = new Map(); for (const r of R) for (const l of r.chs) m.set(A.chName(l), 1); return [...m.keys()]; }
function shell() {
  EL.innerHTML = `<div class="wrap"><div class="bar"><div class="seg" id="smode">${MODES.map(m => `<button data-m="${m.k}" class="${m.k === T.mode ? 'on' : ''}" title="${esc(m.hint)}">${m.label}</button>`).join('')}</div>
    <select class="in" id="sscope" style="width:auto;margin-left:auto"><option value="all">All reactions</option><option value="filters">Reactions tab filters</option>${chapters().map(c => `<option value="ch:${esc(c)}">${esc(c)}</option>`).join('')}${T.scope.startsWith('r:') ? `<option value="${T.scope}">This reaction only</option>` : ''}</select></div>
    <div class="qcard" id="qcard"></div><div class="sfoot"><div class="stats" id="sstats"></div><div class="ses" id="sses"></div></div></div>`;
  $('#sscope', EL).value = T.scope;
  $('#smode', EL).onclick = e => { const b = e.target.closest('[data-m]'); if (b) setMode(b.dataset.m); };
  $('#sscope', EL).onchange = e => { T.scope = e.target.value; T.extra = 0; T.ses = []; next(); };
}
function stageHTML(it) {
  const r = it.r;
  if (T.mode === 'class') { const s = ART.ex[r.id][it.k]; return { q: 'What forms?', cls: 'hidep', svg: s, a: `<span class="mute" style="font:14px var(--font)">from</span> ${fmt(title(r))}` }; }
  if (T.mode === 'product') return { q: 'What forms?', cls: 'hidep', svg: ART.rx[r.id], a: fmt(title(r)) };
  if (T.mode === 'reagent') return { q: 'Which reagents and conditions?', cls: 'hider', svg: ART.rx[r.id], a: fmt([...r.rl, ...r.cl].join(' · ')) || '' };
  if (T.mode === 'name') return { q: 'Name this reaction', cls: '', svg: ART.rx[r.id], a: fmt(r.name) };
  return { q: `What do these arrows give? Step ${it.k + 1} of ${r.ns}`, cls: 'hidep', svg: (ART.ms[r.id] || [])[it.k], a: fmt(title(r)), step: true };
}
function next() {
  T.cur = pick(); T.shown = false; const q = $('#qcard', EL); const st = T.stats;
  $('#sstats', EL).innerHTML = `<span><b>${st.due}</b>due</span><span><b>${st.fresh}</b>new</span><span><b>${st.learned}</b>learned</span><span><b>${st.total}</b>in this set</span>`;
  if (!T.cur) { q.innerHTML = `<div class="empty" style="margin:auto"><b>${st.total ? 'All caught up' : 'Nothing in this set'}</b>${st.total ? 'Nothing is due right now.' : 'Pick another mode or scope.'}${st.total ? '<div style="margin-top:14px"><button class="btn" id="smore">Keep practising anyway</button></div>' : ''}</div>`; const m = $('#smore', q); if (m) m.onclick = () => { T.extra = 1; next(); }; return; }
  const h = stageHTML(T.cur); const nw = Math.round(natW(h.svg) * 1.3);
  q.innerHTML = `<div class="q">${esc(h.q)}</div><div class="stage ${h.cls}" id="sstage" style="--nw:${nw}px">${h.step ? (h.svg || '') : (h.svg || '')}</div><div class="ans" id="sans"></div>
    <div class="qfoot"><button class="btn" id="sopen">Open reaction</button><span class="sp"></span><button class="btn acc" id="sshow">Show answer <kbd>Space</kbd></button>
    <div class="grade" id="sgrade" hidden>${[['Again', 1], ['Hard', 2], ['Good', 3], ['Easy', 4]].map(([l, g]) => `<button class="btn g${g}${g === 3 ? ' acc' : ''}" data-g="${g}"><span>${l}</span><em>${ivl(preview(T.cur, g))}</em><kbd>${g}</kbd></button>`).join('')}</div></div>`;
  T.t0 = Date.now(); sesHTML();
  $('#sshow', q).onclick = reveal; $('#sopen', q).onclick = () => A.openViewer(T.cur.r.id, [T.cur.r.id], { tab: T.mode === 'step' ? 'mech' : T.mode === 'class' ? 'ex' : 'reaction' });
  $('#sgrade', q).onclick = e => { const b = e.target.closest('[data-g]'); if (b) grade(+b.dataset.g); };
}
function reveal() {
  if (!T.cur || T.shown) return; T.shown = true; const h = stageHTML(T.cur); const s = $('#sstage', EL);
  s.classList.remove('hidep', 'hider'); s.classList.add('revealed'); $('#sans', EL).innerHTML = h.a; $('#sans', EL).classList.add('shown'); $('#sshow', EL).hidden = true; $('#sgrade', EL).hidden = false;
}
function sesHTML() {
  const el = $('#sses', EL); if (!el) return; const n = T.ses.length;
  if (!n) { el.innerHTML = '<span class="faint">Space shows the answer · 1 to 4 grades it</span>'; return; }
  const ok = T.ses.filter(x => x.g >= 3).length; const avg = T.ses.reduce((a, x) => a + x.ms, 0) / n / 1000;
  el.innerHTML = `<span><b>${n}</b>reviewed</span><span><b>${Math.round(ok / n * 100)}%</b>recalled</span><span><b>${avg < 60 ? avg.toFixed(1) + 's' : (avg / 60).toFixed(1) + 'm'}</b>per card</span><div class="ticks">${T.ses.slice(-40).map(x => `<i class="g${x.g}"></i>`).join('')}</div>`;
}
function setMode(k) { T.mode = k; T.extra = 0; T.ses = []; if (EL) { for (const b of $$('#smode button', EL)) b.classList.toggle('on', b.dataset.m === k); next(); } }
A.VIEWS.study = {
  mount(el) { EL = el; el.classList.add('sv'); shell(); next(); },
  show() { if (EL) next(); }
};
document.addEventListener('keydown', e => {
  if (A.S.view !== 'study' || !$('#viewer').hidden || !$('#cmdp').hidden || (e.target.matches && e.target.matches('input,select,textarea'))) return;
  if ((e.key === ' ' || e.key === 'Enter') && !T.shown) { e.preventDefault(); reveal(); }
  else if (T.shown && /^[1-4]$/.test(e.key)) grade(+e.key);
});
window.Study = { MODES, setMode, quizReaction(id) { T.scope = 'r:' + id; T.mode = 'class'; T.extra = 1; if (EL) shell(); if (EL) next(); } };
})();
