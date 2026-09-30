// State -> image renderers. All colour math is copied verbatim from
// lgca-lab.html so the package shows the tool's own visuals:
//   density  = the tool's default "density (dark -> bright)" field view
//   direction = the tool's "net motion (hue = direction)" field view
//   particles = direction-coloured particle dots (the tool's particle palette)
// Layout matches the tool's field renderer: site (i,j) -> 2x2 block at
// x=(2i+j) mod 2W, y=2j (the hexagonal shear baked into pixel position).
// Display aspect: true hex geometry is sqrt(3)/2 ~ 0.866 vertical scale,
// applied via CSS in the table page (images are square in pixel space).

export const BG = [0x0b, 0x0e, 0x13]; // the tool's canvas background
export const DIRCOL = ["#e5484d", "#e0a336", "#7fbf3f", "#2bbfa3", "#4f9de6", "#b57bd6"]
  .map(h => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)]);
const EH = Math.sqrt(3.0) / 2.0;
export const CX = [1.0, 0.5, -0.5, -1.0, -0.5, 0.5];
export const CY = [0.0, EH, EH, 0.0, -EH, -EH];

// The tool's DLUT: 5 colour stops interpolated to 256 entries.
export const DLUT = (() => {
  const stops = [[8, 12, 28], [22, 48, 94], [42, 120, 214], [125, 184, 240], [255, 255, 255]];
  const L = [];
  for (let q = 0; q < 256; q++) {
    const t = q / 255 * (stops.length - 1), a = Math.min(stops.length - 2, Math.floor(t)), f = t - a;
    L.push([
      (stops[a][0] + (stops[a + 1][0] - stops[a][0]) * f) | 0,
      (stops[a][1] + (stops[a + 1][1] - stops[a][1]) * f) | 0,
      (stops[a][2] + (stops[a + 1][2] - stops[a][2]) * f) | 0,
    ]);
  }
  return L;
})();

export function densityIndex(n, dens) {
  return Math.min(255, (Math.sqrt(n / (24 * dens)) * 255) | 0);
}

export function hsl2rgb(h, s, l) {
  h = (((h % 360) + 360) % 360) / 360;
  const q = l < 0.5 ? l * (1 + s) : l + s - l * s, p = 2 * l - q;
  const f = t => {
    t = ((t % 1) + 1) % 1;
    if (t < 1 / 6) return p + (q - p) * 6 * t;
    if (t < 0.5) return q;
    if (t < 2 / 3) return p + (q - p) * (2 / 3 - t) * 6;
    return p;
  };
  return [(f(h + 1 / 3) * 255) | 0, (f(h) * 255) | 0, (f(h - 1 / 3) * 255) | 0];
}

function blank(w, h) {
  const img = Buffer.alloc(w * h * 4);
  for (let p = 0; p < w * h; p++) {
    img[p * 4] = BG[0]; img[p * 4 + 1] = BG[1]; img[p * 4 + 2] = BG[2]; img[p * 4 + 3] = 255;
  }
  return img;
}
function put2x2(img, W2, x, y, rgb) {
  for (const dx of [0, 1]) {
    const xx = (x + dx) % W2;
    for (const dy of [0, 1]) {
      const p = ((y + dy) * W2 + xx) * 4;
      img[p] = rgb[0]; img[p + 1] = rgb[1]; img[p + 2] = rgb[2]; img[p + 3] = 255;
    }
  }
}

// occ: Int32Array/array W*H*6 in the engine's layout occ[(j*W+i)*6+k]
export function renderDensity(occ, W, H, dens) {
  const W2 = 2 * W, img = blank(W2, 2 * H);
  for (let j = 0; j < H; j++) for (let i = 0; i < W; i++) {
    const b = (j * W + i) * 6;
    const n = occ[b] + occ[b + 1] + occ[b + 2] + occ[b + 3] + occ[b + 4] + occ[b + 5];
    if (!n) continue;
    put2x2(img, W2, (2 * i + j) % W2, 2 * j, DLUT[densityIndex(n, dens)]);
  }
  return {w: W2, h: 2 * H, rgba: img};
}

export function renderDirection(occ, W, H) {
  const W2 = 2 * W, img = blank(W2, 2 * H);
  for (let j = 0; j < H; j++) for (let i = 0; i < W; i++) {
    const b = (j * W + i) * 6;
    const n0 = occ[b], n1 = occ[b + 1], n2 = occ[b + 2], n3 = occ[b + 3], n4 = occ[b + 4], n5 = occ[b + 5];
    const n = n0 + n1 + n2 + n3 + n4 + n5;
    if (!n) continue;
    const fx = n0 - n3 + 0.5 * (n1 - n2 - n4 + n5), fy = 0.8660254 * (n1 + n2 - n4 - n5);
    const m = Math.hypot(fx, fy) / n;
    put2x2(img, W2, (2 * i + j) % W2, 2 * j,
      hsl2rgb(Math.atan2(fy, fx) * 57.29578, 0.15 + 0.75 * m, 0.28 + 0.37 * m));
  }
  return {w: W2, h: 2 * H, rgba: img};
}

// 4 px/site; dots offset along the channel direction like the tool's raw view.
export function renderParticles(occ, W, H) {
  const W4 = 4 * W, img = blank(W4, 4 * H);
  const put = (x, y, rgb) => {
    const p = (y * W4 + ((x % W4 + W4) % W4)) * 4;
    img[p] = rgb[0]; img[p + 1] = rgb[1]; img[p + 2] = rgb[2]; img[p + 3] = 255;
  };
  for (let j = 0; j < H; j++) for (let i = 0; i < W; i++) {
    const b = (j * W + i) * 6;
    for (let k = 0; k < 6; k++) {
      const n = occ[b + k];
      if (!n) continue;
      const x0 = (4 * i + 2 * j) % W4;
      const px = x0 + 1 + Math.round(1.2 * CX[k] + 0.5);
      const py = 4 * j + 1 + Math.round(-1.2 * CY[k] + 0.5);
      const yy = Math.min(4 * H - 1, Math.max(0, py));
      put(px, yy, DIRCOL[k]);
      if (n >= 2) { // bigger dot for multiple particles in one channel
        put(px + 1, yy, DIRCOL[k]);
        if (yy + 1 < 4 * H) { put(px, yy + 1, DIRCOL[k]); put(px + 1, yy + 1, DIRCOL[k]); }
      }
    }
  }
  return {w: W4, h: 4 * H, rgba: img};
}

// High-resolution particle view for the lightbox: S px per site (default 8),
// a d×d dot per occupied channel (d = S/4) placed 0.375·S from the cell centre
// along the channel direction (neighbouring dots do not overlap), one pixel
// larger when the channel holds 2+ particles. Same hex shear as above.
export function renderParticlesHi(occ, W, H, S = 8) {
  const WS = S * W, HS = S * H, img = blank(WS, HS);
  const d = Math.max(1, S >> 2), r = 0.375 * S, c = S / 2;
  const dot = (x, y, size, rgb) => {
    for (let dy = 0; dy < size; dy++) {
      const yy = y + dy; if (yy < 0 || yy >= HS) continue;
      for (let dx = 0; dx < size; dx++) {
        const xx = (((x + dx) % WS) + WS) % WS;
        const p = (yy * WS + xx) * 4;
        img[p] = rgb[0]; img[p + 1] = rgb[1]; img[p + 2] = rgb[2]; img[p + 3] = 255;
      }
    }
  };
  for (let j = 0; j < H; j++) for (let i = 0; i < W; i++) {
    const b = (j * W + i) * 6;
    const x0 = (S * i + (S / 2) * j) % WS;
    for (let k = 0; k < 6; k++) {
      const n = occ[b + k];
      if (!n) continue;
      const size = n >= 2 ? d + 1 : d;
      dot(Math.round(x0 + c + r * CX[k] - size / 2), Math.round(S * j + c - r * CY[k] - size / 2), size, DIRCOL[k]);
    }
  }
  return {w: WS, h: HS, rgba: img};
}

// GIF palette for density clips: index 0 = background, 1..255 = DLUT (entries
// 1..255; index 0 of DLUT is never produced by an occupied site at dens 0.4).
export function gifDensityPalette() {
  const pal = DLUT.map(c => c.slice());
  pal[0] = BG.slice();
  return pal;
}
// counts: Uint8Array W*H (per-site totals) -> palette-index frame at 2 px/site
export function gifDensityFrame(counts, W, H, dens) {
  const W2 = 2 * W, out = new Uint8Array(W2 * 2 * H);
  for (let j = 0; j < H; j++) for (let i = 0; i < W; i++) {
    const n = counts[j * W + i];
    const idx = n ? Math.max(1, densityIndex(n, dens)) : 0;
    const x = (2 * i + j) % W2, y = 2 * j;
    out[y * W2 + x] = idx; out[y * W2 + (x + 1) % W2] = idx;
    out[(y + 1) * W2 + x] = idx; out[(y + 1) * W2 + (x + 1) % W2] = idx;
  }
  return out;
}
