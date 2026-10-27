/* graph.js · reactions-only graph: pointer pushes dots aside, a fast pointer can catch and drag one */
(() => {
'use strict';
const A = window.App; if (!A) return;
const { R, BY, ART, D, $, $$, esc, fmt, title, chLine, clsOf, CLS, FGL, RGL } = A;
function oklchRGB(L, C, h, a = 1) {
  const hr = h * Math.PI / 180, a_ = C * Math.cos(hr), b_ = C * Math.sin(hr);
  const l = (L + .3963377774 * a_ + .2158037573 * b_) ** 3, m = (L - .1055613458 * a_ - .0638541728 * b_) ** 3, s = (L - .0894841775 * a_ - 1.291485548 * b_) ** 3;
  const f = x => { x = Math.max(0, Math.min(1, x)); return Math.round((x <= .0031308 ? 12.92 * x : 1.055 * Math.pow(x, 1 / 2.4) - .055) * 255); };
  return `rgba(${f(4.0767416621 * l - 3.3077115913 * m + .2309699292 * s)},${f(-1.2684380046 * l + 2.6097574011 * m - .3413193965 * s)},${f(-.0041960863 * l - .7034186147 * m + 1.707614701 * s)},${a})`;
}
const cs = () => getComputedStyle(document.documentElement);
function cssColor(name, alpha) { const v = cs().getPropertyValue(name).trim(); const m = /oklch\(\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*(?:\/\s*([\d.]+))?\s*\)/.exec(v); return m ? oklchRGB(+m[1], +m[2], +m[3], alpha != null ? alpha : (m[4] ? +m[4] : 1)) : v; }
let G = null;
A.VIEWS.graph = {
  mount(el) {
    el.classList.add('grv');
    el.innerHTML = `<div class="gcan" id="gcan" tabindex="-1"><canvas id="gcv"></canvas><div class="legend"><span><i></i>same conversion</span><span><i class="d"></i>shared reagent</span><span><i class="a"></i>next step</span></div><div class="hint" id="ghint">Hover to reveal labels, click to select, drag to reposition. Slow movement gently repels nearby dots. Double-click for a neighbourhood. Scroll to zoom.</div></div><aside class="gside" id="gside"></aside>`;
    build(el);
  },
  show() { if (G) { G.resize(); G.sim.alpha(Math.max(G.sim.alpha(), .2)).restart(); } },
  hide() { if (G) G.sim.stop(); }
};
function build(el) {
  const cv = $('#gcv', el), box = $('#gcan', el), side = $('#gside', el), ctx = cv.getContext('2d');
  const POS = Object.fromEntries(D.FG.map(n => [n.k, n]));
  const nodes = R.map((r, i) => { const p = POS[r.fg[0]] || { x: 6, y: 5 }, q = POS[r.fg[1]] || p; return { i, r, x: ((p.x * .6 + q.x * .4) - 6) * 64 + (Math.random() - .5) * 50, y: ((p.y * .6 + q.y * .4) - 5) * 64 + (Math.random() - .5) * 50, cls: clsOf(r), rad: r.name ? 5.8 : 4.8 }; });
  const byId = new Map(nodes.map(n => [n.r.id, n]));
  const links = [], seen = new Set();
  const link = (a, b, t) => { const k = a < b ? a + '|' + b : b + '|' + a; if (a === b || seen.has(k)) return; seen.add(k); links.push({ source: a, target: b, t }); };
  const conv = new Map(); for (const n of nodes) { const k = n.r.fg.join('>'); if (!conv.has(k)) conv.set(k, []); conv.get(k).push(n); }
  for (const g of conv.values()) { g.sort((a, b) => a.r._k - b.r._k); for (let i = 1; i < g.length; i++) link(g[i - 1].i, g[i].i, 'c'); }
  for (const e of RGL) { if (e.main.size < 2 || e.main.size > 6) continue; const ids = [...e.main].map(id => byId.get(id).i); for (let i = 0; i < ids.length; i++) for (let j = i + 1; j < ids.length; j++) link(ids[i], ids[j], 'r'); }
  const adj = nodes.map(() => new Set()); for (const l of links) { adj[l.source].add(l.target); adj[l.target].add(l.source); }
  const nextOf = n => R.filter(x => x !== n.r && x.fg[0] === n.r.fg[1] && x.fg[0] !== x.fg[1]).map(x => byId.get(x.id)).slice(0, 40);
  let W = 0, H = 0, dpr = 1; const tf = { x: 0, y: 0, k: 1 };
  const P = { x: 0, y: 0, lx: 0, ly: 0, wx: 0, wy: 0, inside: false, speed: 0, dt: 16, lt: 0, down: null };
  let hover = null, sel = null, drag = null, vis = null, depth = 0, showNext = true, nextList = [], idleT = 0;
  let C = {};
  function colors() {
    const kl = parseFloat(cs().getPropertyValue('--kl')) || .79, kc = parseFloat(cs().getPropertyValue('--kc')) || .1;
    C = { dot: cssColor('--dot'), link: cssColor('--line2'), acc: cssColor('--acc'), ink: cssColor('--ink'), ink2: cssColor('--ink2'), bg: cssColor('--bg'), mute: cssColor('--mute'), k: {} };
    for (const k in CLS) C.k[k] = oklchRGB(kl, kc * (CLS[k][2] || 1), CLS[k][1]);
  }
  colors(); window.addEventListener('themechange', () => { colors(); draw(); if (sel) sideSel(); else sideIdle(); });
  // A small, local pointer force: nearby dots make room, but the graph never runs away.
  function fPointer() {
    if (!P.inside || drag || P.down) return;
    const speed = Math.hypot(P.x - P.lx, P.y - P.ly) / Math.max(1, P.dt || 16); P.lx = P.x; P.ly = P.y;
    const slow = Math.max(0, 1 - speed / 1.8); if (slow < .12) return;
    const rw = 34 / tf.k, r2 = rw * rw;
    for (const n of nodes) { if (n === hover) continue; const dx=n.x-P.wx,dy=n.y-P.wy,d2=dx*dx+dy*dy; if(d2>r2||d2<.01)continue; const d=Math.sqrt(d2),f=(1-d/rw)*.045*slow; n.vx+=dx/d*f;n.vy+=dy/d*f; n.vx=Math.max(-1.2,Math.min(1.2,n.vx));n.vy=Math.max(-1.2,Math.min(1.2,n.vy)); }
  }
  /* arrangement: free blob, or clusters by chapter / class (cluster centres packed with a tiny simulation) */
  let arr = localStorage.getItem('ocr3-garr') || 'free'; let groups = [];
  const groupKey = (n, mode) => mode === 'chapter' ? ((n.r.chs[0] || {}).chapter || 'Other') : mode === 'class' ? (n.cls[0] || 'other') : '';
  function layoutGroups(mode) {
    groups = []; if (mode === 'free') { for (const n of nodes) { n.gx = 0; n.gy = 0; n.g = null; } return; }
    const m = new Map(); for (const n of nodes) { const k = groupKey(n, mode); if (!m.has(k)) m.set(k, []); m.get(k).push(n); }
    groups = [...m.entries()].map(([k, list]) => ({ k, list, r: 26 + Math.sqrt(list.length) * 13.5, label: mode === 'class' ? (CLS[k] ? CLS[k][0] : k) : k }));
    groups.sort((a, b) => b.list.length - a.list.length);
    groups.forEach((g, i) => { const t = i * 2.39996, rr = 46 * Math.sqrt(i); g.x = Math.cos(t) * rr; g.y = Math.sin(t) * rr * .78; });
    const gs = d3.forceSimulation(groups).force('c', d3.forceCollide(g => g.r + 34).iterations(3)).force('x', d3.forceX(0).strength(.06)).force('y', d3.forceY(0).strength(.085)).stop();
    for (let i = 0; i < 320; i++) gs.tick();
    for (const g of groups) for (const n of g.list) { n.gx = g.x; n.gy = g.y; n.g = g; }
  }
  layoutGroups(arr);
  const sim = d3.forceSimulation(nodes)
    .force('charge', d3.forceManyBody().strength(() => arr === 'free' ? -30 : -14).distanceMax(300))
    .force('link', d3.forceLink(links).distance(l => l.t === 'c' ? 26 : 52).strength(l => arr === 'free' ? (l.t === 'c' ? .5 : .1) : (l.source.g === l.target.g ? (l.t === 'c' ? .35 : .05) : 0)))
    .force('x', d3.forceX(n => n.gx || 0).strength(() => arr === 'free' ? .04 : .16)).force('y', d3.forceY(n => n.gy || 0).strength(() => arr === 'free' ? .055 : .16))
    .force('collide', d3.forceCollide(n => n.rad + 3.2))
    .force('pointer', fPointer).alphaDecay(.028).on('tick', draw);
  G = { sim, resize };
  function setArr(mode) {
    arr = mode; try { localStorage.setItem('ocr3-garr', mode); } catch (e) { }
    layoutGroups(mode);
    sim.force('x').x(n => n.gx || 0).strength(() => arr === 'free' ? .04 : .16); sim.force('y').y(n => n.gy || 0).strength(() => arr === 'free' ? .055 : .16);
    sim.force('link').strength(l => arr === 'free' ? (l.t === 'c' ? .5 : .1) : (l.source.g === l.target.g ? (l.t === 'c' ? .35 : .05) : 0));
    sim.force('charge').strength(() => arr === 'free' ? -30 : -14);
    sim.alpha(.9).restart(); fitAfter(1300);
  }
  let fitT = 0;
  function fitAfter(ms) { clearTimeout(fitT); fitT = setTimeout(() => fitAll(true), ms); }
  function fitAll(anim) {
    if (!W) return; let x0 = 1e9, y0 = 1e9, x1 = -1e9, y1 = -1e9;
    for (const n of nodes) { if (vis && !vis.has(n.i)) continue; x0 = Math.min(x0, n.x); y0 = Math.min(y0, n.y); x1 = Math.max(x1, n.x); y1 = Math.max(y1, n.y); }
    if (x0 > x1) return; const pad = 70; const k = Math.max(.25, Math.min(3, Math.min((W - pad * 2) / Math.max(1, x1 - x0), (H - pad * 2 - 40) / Math.max(1, y1 - y0))));
    const tx = W / 2 - (x0 + x1) / 2 * k, ty = H / 2 + 10 - (y0 + y1) / 2 * k;
    if (!anim) { tf.x = tx; tf.y = ty; tf.k = k; draw(); return; }
    const a = { ...tf }, t0 = performance.now();
    const st = t => { const p = Math.min(1, (t - t0) / 560), e = p < .5 ? 4 * p * p * p : 1 - Math.pow(-2 * p + 2, 3) / 2; tf.x = a.x + (tx - a.x) * e; tf.y = a.y + (ty - a.y) * e; tf.k = a.k + (k - a.k) * e; draw(); if (p < 1) requestAnimationFrame(st); };
    requestAnimationFrame(st);
  }
  G.setArr = setArr; G.fit = () => fitAll(true);
  function resize() { const r = box.getBoundingClientRect(); if (!r.width) return; dpr = window.devicePixelRatio || 1; const first = !W; W = r.width; H = r.height; cv.width = W * dpr; cv.height = H * dpr; if (first) { tf.x = W / 2; tf.y = H / 2; tf.k = Math.min(W, H) / 820; } draw(); }
  new ResizeObserver(resize).observe(box);
  const toW = (x, y) => [(x - tf.x) / tf.k, (y - tf.y) / tf.k];
  function hit(x, y) { const [wx, wy] = toW(x, y); let best = null, bd = (17 / tf.k) ** 2; for (const n of nodes) { if (vis && !vis.has(n.i)) continue; const d = (n.x - wx) ** 2 + (n.y - wy) ** 2; if (d < bd) { bd = d; best = n; } } return best; }
  function focusSet() { const f = sel || hover; if (!f) return null; const s = new Set([f.i, ...adj[f.i]]); if (sel && showNext) for (const n of nextList) s.add(n.i); return s; }
  function draw() {
    if (!W) return;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
    let step = 26 * tf.k; while (step < 14) step *= 2; while (step > 52) step /= 2;
    ctx.fillStyle = C.dot; const ox = ((tf.x % step) + step) % step, oy = ((tf.y % step) + step) % step;
    for (let x = ox; x < W; x += step) for (let y = oy; y < H; y += step) ctx.fillRect(x - .7, y - .7, 1.4, 1.4);
    ctx.translate(tf.x, tf.y); ctx.scale(tf.k, tf.k);
    const f = sel || hover, fs = focusSet(), lw = 1 / tf.k;
    for (const l of links) {
      const a = l.source, b = l.target; if (vis && (!vis.has(a.i) || !vis.has(b.i))) continue;
      const hot = f && (a === f || b === f);
      ctx.globalAlpha = hot ? .95 : f ? .08 : groups.length ? (a.g === b.g ? (l.t === 'c' ? .45 : .25) : .06) : (l.t === 'c' ? .5 : .3); ctx.strokeStyle = hot ? C.acc : C.link; ctx.lineWidth = (hot ? 1.6 : 1) * lw;
      ctx.setLineDash(l.t === 'r' ? [4 * lw, 3.5 * lw] : []); ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.stroke();
    }
    ctx.setLineDash([]);
    if (sel && showNext) {
      ctx.strokeStyle = C.acc; ctx.fillStyle = C.acc; ctx.lineWidth = 1.4 * lw;
      for (const n of nextList) {
        if (vis && !vis.has(n.i)) continue; ctx.globalAlpha = .8;
        const dx = n.x - sel.x, dy = n.y - sel.y, d = Math.hypot(dx, dy) || 1, ux = dx / d, uy = dy / d;
        const mx = (sel.x + n.x) / 2 - uy * d * .18, my = (sel.y + n.y) / 2 + ux * d * .18;
        const ex = n.x - ux * (n.rad + 3), ey = n.y - uy * (n.rad + 3);
        ctx.beginPath(); ctx.moveTo(sel.x, sel.y); ctx.quadraticCurveTo(mx, my, ex, ey); ctx.stroke();
        const tx = ex - mx, ty = ey - my, tl = Math.hypot(tx, ty) || 1, hx = tx / tl, hy = ty / tl, hl = 7 * lw * 1.2;
        ctx.beginPath(); ctx.moveTo(ex, ey); ctx.lineTo(ex - hx * hl - hy * hl * .45, ey - hy * hl + hx * hl * .45); ctx.lineTo(ex - hx * hl * .6, ey - hy * hl * .6); ctx.lineTo(ex - hx * hl + hy * hl * .45, ey - hy * hl - hx * hl * .45); ctx.closePath(); ctx.fill();
      }
    }
    for (const n of nodes) {
      if (vis && !vis.has(n.i)) continue;
      ctx.globalAlpha = fs && !fs.has(n.i) ? .22 : 1;
      const k = n.cls.length;
      if (k === 1) { ctx.fillStyle = C.k[n.cls[0]] || C.mute; ctx.beginPath(); ctx.arc(n.x, n.y, n.rad, 0, 7); ctx.fill(); }
      else for (let j = 0; j < k; j++) { ctx.fillStyle = C.k[n.cls[j]] || C.mute; ctx.beginPath(); ctx.moveTo(n.x, n.y); ctx.arc(n.x, n.y, n.rad, -Math.PI / 2 + j * 2 * Math.PI / k, -Math.PI / 2 + (j + 1) * 2 * Math.PI / k); ctx.closePath(); ctx.fill(); }
      ctx.lineWidth = 1.4 * lw; ctx.strokeStyle = C.bg; ctx.beginPath(); ctx.arc(n.x, n.y, n.rad + .6 * lw, 0, 7); ctx.stroke();
      if (n === sel || n === hover) { ctx.lineWidth = 2 * lw; ctx.strokeStyle = n === sel ? C.acc : C.ink2; ctx.beginPath(); ctx.arc(n.x, n.y, n.rad + 3.2 * lw, 0, 7); ctx.stroke(); }
    }
    ctx.globalAlpha = 1; ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    if (groups.length && !vis) {
      ctx.font = '600 10.5px "Source Code Pro", Consolas, monospace'; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
      for (const g of groups) {
        let gx = 0, gy = 1e9, cx = 0; for (const n of g.list) { cx += n.x; gy = Math.min(gy, n.y - n.rad); } cx /= g.list.length;
        const sx = cx * tf.k + tf.x, sy = gy * tf.k + tf.y - 9; const hot = f && f.g === g;
        let t = g.label.toUpperCase(); if (t.length > 30) t = t.slice(0, 28) + '…';
        ctx.globalAlpha = f && !hot ? .35 : .9; ctx.lineWidth = 4; ctx.strokeStyle = C.bg; ctx.strokeText(t + '  ' + g.list.length, sx, sy); ctx.fillStyle = hot ? C.acc : C.mute; ctx.fillText(t + '  ' + g.list.length, sx, sy);
      }
      ctx.globalAlpha = 1; ctx.textAlign = 'left';
    }
    const lab = new Set(); if (f) { lab.add(f); if (adj[f.i].size <= 24) for (const j of adj[f.i]) lab.add(nodes[j]); } if (sel && showNext && nextList.length <= 16) for (const n of nextList) lab.add(n);
    if (tf.k > 1.9) for (const n of nodes) { const sx = n.x * tf.k + tf.x, sy = n.y * tf.k + tf.y; if (sx > -50 && sx < W + 50 && sy > -20 && sy < H + 20) lab.add(n); }
    ctx.font = '12.5px "Nimbus Sans", Helvetica, Arial, sans-serif'; ctx.textBaseline = 'middle'; ctx.lineJoin = 'round';
    for (const n of lab) {
      if (vis && !vis.has(n.i)) continue;
      const sx = n.x * tf.k + tf.x + n.rad * tf.k + 6, sy = n.y * tf.k + tf.y; let t = title(n.r); if (t.length > 34 && n !== f) t = t.slice(0, 32) + '…';
      t = t.replace(/[₀-₉]/g, c => '0123456789'['₀₁₂₃₄₅₆₇₈₉'.indexOf(c)]);
      ctx.lineWidth = 4; ctx.strokeStyle = C.bg; ctx.strokeText(t, sx, sy); ctx.fillStyle = n === f ? C.ink : C.ink2; ctx.fillText(t, sx, sy);
    }
  }
  function warm() { sim.alphaTarget(.035); if (sim.alpha() < .04) sim.alpha(.04); sim.restart(); clearTimeout(idleT); idleT = setTimeout(() => sim.alphaTarget(0), 500); }
  box.addEventListener('pointermove', e => {
    const r = box.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
    const dt = Math.max(1, e.timeStamp - P.lt), d = Math.hypot(x - P.x, y - P.y);
    P.dt = dt; P.speed = P.inside ? P.speed * .6 + (d / dt) * .4 : 0; P.x = x; P.y = y; P.lt = e.timeStamp; P.inside = true; [P.wx, P.wy] = toW(x, y);
    if (drag) { drag.fx = P.wx; drag.fy = P.wy; if (P.down) P.down.moved += d; sim.alphaTarget(.12).restart(); return; }
    if (P.down && P.down.pan) { tf.x += x - P.down.x; tf.y += y - P.down.y; P.down.x = x; P.down.y = y; P.down.moved += d; draw(); return; }
    const h = hit(x, y); if (h !== hover) { hover = h; box.style.cursor = h ? 'grab' : 'default'; draw(); }
    // Keep the graph settled while selecting. The pointer only gives a tiny local nudge when moving slowly.
  });
  box.addEventListener('pointerleave', () => { P.inside = false; if (!drag) { hover = null; draw(); } });
  box.addEventListener('pointerdown', e => {
    box.setPointerCapture(e.pointerId); const r = box.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
    const h = hit(x, y); P.down = { x, y, moved: 0, pan: !h, node: h };
    if (h) { drag = h; h.fx = h.x; h.fy = h.y; box.style.cursor = 'grabbing'; sim.alphaTarget(.04).restart(); }
  });
  box.addEventListener('pointerup', e => {
    const d = P.down; P.down = null;
    if (drag) { const n = drag; drag = null; n.fx = null; n.fy = null; sim.alphaTarget(0); box.style.cursor = 'grab'; if (d && d.moved < 5) select(n); }
    else if (d && d.pan && d.moved < 5 && sel && e.detail < 2) { select(null); }
  });
  box.addEventListener('dblclick', e => { const r = box.getBoundingClientRect(); const h = hit(e.clientX - r.left, e.clientY - r.top); if (h) { select(h); setDepth(1); } else { setDepth(0); } });
  box.addEventListener('wheel', e => {
    e.preventDefault(); const r = box.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
    const k2 = Math.max(.25, Math.min(6, tf.k * Math.exp(-e.deltaY * .0015))); tf.x = x - (x - tf.x) * k2 / tf.k; tf.y = y - (y - tf.y) * k2 / tf.k; tf.k = k2; draw();
  }, { passive: false });
  document.addEventListener('keydown', e => { if (A.S.view !== 'graph' || e.key !== 'Escape' || !$('#viewer').hidden || !$('#cmdp').hidden) return; if (depth) setDepth(0); else select(null); });
  function centerOn(n) {
    const tx = W / 2 - 180 * 0 - n.x * tf.k, ty = H / 2 - n.y * tf.k, x0 = tf.x, y0 = tf.y, t0 = performance.now();
    const st = t => { const p = Math.min(1, (t - t0) / 420), e = 1 - Math.pow(1 - p, 3); tf.x = x0 + (tx - x0) * e; tf.y = y0 + (ty - y0) * e; draw(); if (p < 1) requestAnimationFrame(st); };
    requestAnimationFrame(st);
  }
  function setDepth(dp) {
    depth = dp;
    if (!dp || !sel) { vis = null; depth = 0; }
    else { vis = new Set([sel.i]); let fr = [sel.i]; for (let d = 0; d < dp; d++) { const nx = []; for (const i of fr) for (const j of adj[i]) if (!vis.has(j)) { vis.add(j); nx.push(j); } fr = nx; } if (showNext) for (const n of nextList) vis.add(n.i); }
    sim.alpha(.3).restart(); draw(); if (sel) sideSel(); else sideIdle();
  }
  function select(n) { sel = n; nextList = n ? nextOf(n) : []; if (!n) { vis = null; depth = 0; sideIdle(); } else { if (depth) setDepth(depth); sideSel(); } draw(); }
  function sideIdle() {
    side.innerHTML = `<h2>Reaction graph</h2><p class="mute" style="margin:8px 0 0">Each dot is one reaction from the library, coloured by its class. Solid lines join reactions that do the same conversion. Dashed lines join reactions that share an uncommon reagent.</p>
      <h4>Classes</h4><div class="chips">${Object.keys(CLS).filter(k => R.some(r => clsOf(r)[0] === k)).map(k => `<span class="chip" style="--k:${A.kcol(k)};cursor:default"><i class="sw8"></i>${CLS[k][0]}</span>`).join('')}</div>
      <h4>Arrange</h4><div class="seg" id="garr">${[['free', 'Free'], ['chapter', 'By chapter'], ['class', 'By class']].map(([v, l]) => `<button data-a="${v}" class="${arr === v ? 'on' : ''}">${l}</button>`).join('')}</div>
      <h4>Find a reaction</h4><input class="in" id="gfind" placeholder="Type a name, then Enter" autocomplete="off" spellcheck="false">
      <button class="btn sm" id="gfit" style="margin-top:14px">Fit everything</button>`;
    $('#garr', side).onclick = e => { const bt = e.target.closest('[data-a]'); if (!bt || bt.dataset.a === arr) return; setArr(bt.dataset.a); for (const x of $$('#garr button', side)) x.classList.toggle('on', x === bt); };
    $('#gfit', side).onclick = () => fitAll(true);
    const f = $('#gfind', side); f.onkeydown = e => { if (e.key !== 'Enter') return; const q = A.norm(f.value).trim(); if (!q) return; const n = nodes.find(n => A.norm(title(n.r)).includes(q)); if (n) { select(n); centerOn(n); } else A.toast('No reaction by that name'); };
  }
  function sideSel() {
    const r = sel.r; const [a, b] = r.fg;
    side.innerHTML = `<h2>${fmt(title(r))}</h2><div class="meta">${esc(chLine(r))}</div><div class="art" id="gart">${ART.rx[r.id] || ''}</div>
      <div class="meta">${esc(FGL[a])}${a !== b ? ' → ' + esc(FGL[b]) : ''}</div>
      <div class="acts"><button class="btn acc sm" id="gopen">Open</button><button class="btn sm" id="gclr">Clear</button></div>
      <div style="display:flex;flex-direction:column;gap:10px">${A.bondSw('gnext', showNext, `Show next steps (${nextList.length})`)}
      <div><h4 style="margin-top:6px">Neighbourhood</h4><div class="seg" id="gdep">${[['0', 'All'], ['1', '1 link'], ['2', '2 links'], ['3', '3 links']].map(([v, l]) => `<button data-d="${v}" class="${+v === depth ? 'on' : ''}">${l}</button>`).join('')}</div></div></div>
      ${nextList.length ? `<h4>What can come next</h4><div class="rows">${nextList.slice(0, 14).map(n => `<button class="row" data-gi="${n.i}"><span class="t">${fmt(title(n.r))}</span></button>`).join('')}</div>` : ''}`;
    $('#gopen', side).onclick = () => { const ids = [r.id, ...nextList.map(n => n.r.id)]; A.openViewer(r.id, ids, { from: $('#gart', side).getBoundingClientRect() }); };
    $('#gclr', side).onclick = () => select(null);
    $('#gnext', side).onclick = () => { showNext = !showNext; if (depth) setDepth(depth); else { sideSel(); draw(); } };
    $('#gdep', side).onclick = e => { const bt = e.target.closest('[data-d]'); if (bt) setDepth(+bt.dataset.d); };
    side.querySelectorAll('[data-gi]').forEach(bt => bt.onclick = () => { const n = nodes[+bt.dataset.gi]; select(n); centerOn(n); });
  }
  sideIdle(); resize(); fitAfter(1600);
}
})();
