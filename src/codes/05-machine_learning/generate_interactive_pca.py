#!/usr/bin/env python

import numpy as np

svg_template = '''
<script>
// 1) Data array
const data = /*__DATA_PLACEHOLDER__*/;

// 2) Frame extents in data coords (customizable)
const frameExt = { xMin: -2.5, xMax: 3.2, yMin: -2.9, yMax: 2.2 };
const xTicks = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
const yTicks = [-2.0, -1.0, 0.0, 1.0, 2.0];
const tickFont = 15;
const cornerRadius = 5;

// 3) Canvas and margins
const W = 600, H = 540;
const margin = { top: 25, right: 25, bottom: 50, left: 60 };
const innerW = W - margin.left - margin.right;
const innerH = H - margin.top - margin.bottom;
const svgNS = "http://www.w3.org/2000/svg";
const svg = document.createElementNS(svgNS, "svg");
svg.setAttribute("width", W);
svg.setAttribute("height", H);
svg.setAttribute("style", "display:block; margin:0 auto; background-color:white;");
document.getElementById("scatter-container").appendChild(svg);

// 4) Equal scaling calculation
const dataWidth  = frameExt.xMax - frameExt.xMin;
const dataHeight = frameExt.yMax - frameExt.yMin;
const scale      = Math.min(innerW / dataWidth, innerH / dataHeight);
const xOffset    = margin.left + (innerW - scale * dataWidth) / 2;
const yOffset    = margin.top  + (innerH - scale * dataHeight) / 2;
const midX = W / 2, midY = H / 2;

function xScale(x) { return xOffset + (x - frameExt.xMin) * scale; }
function yScale(y) { return yOffset + (frameExt.yMax - y) * scale; }

// 5) Draw frame rectangle
const frame = document.createElementNS(svgNS, "rect");
frame.setAttribute("x", xOffset);
frame.setAttribute("y", yOffset);
frame.setAttribute("width",  dataWidth * scale);
frame.setAttribute("height", dataHeight * scale);
frame.setAttribute("fill",   "none");
frame.setAttribute("stroke", "black");
svg.appendChild(frame);

// 6) Draw custom ticks and labels
const tickLen = 8;
// X-axis ticks
xTicks.forEach(xVal => {
  const xPos = xScale(xVal);
  const yBase = yOffset + dataHeight * scale;
  // tick line
  const tx = document.createElementNS(svgNS, "line");
  tx.setAttribute("x1", xPos);
  tx.setAttribute("y1", yBase);
  tx.setAttribute("x2", xPos);
  tx.setAttribute("y2", yBase + tickLen);
  tx.setAttribute("stroke", "black");
  svg.appendChild(tx);
  // label
  const lbl = document.createElementNS(svgNS, "text");
  lbl.setAttribute("x", xPos);
  lbl.setAttribute("y", yBase + tickLen + tickFont);
  lbl.setAttribute("font-size", tickFont);
  lbl.setAttribute("text-anchor", "middle");
  lbl.textContent = xVal.toString();
  svg.appendChild(lbl);
});
// Y-axis ticks
yTicks.forEach(yVal => {
  const yPos = yScale(yVal);
  const xBase = xOffset;
  const ty = document.createElementNS(svgNS, "line");
  ty.setAttribute("x1", xBase - tickLen);
  ty.setAttribute("y1", yPos);
  ty.setAttribute("x2", xBase);
  ty.setAttribute("y2", yPos);
  ty.setAttribute("stroke", "black");
  svg.appendChild(ty);
  // label
  const lbl = document.createElementNS(svgNS, "text");
  lbl.setAttribute("x", xBase - tickLen - 4);
  lbl.setAttribute("y", yPos + tickFont/2 - 2);
  lbl.setAttribute("font-size", tickFont);
  lbl.setAttribute("text-anchor", "end");
  lbl.textContent = yVal.toString();
  svg.appendChild(lbl);
});

// 7} Draw axis labels
// X-axis label
const xlabel = document.createElementNS(svgNS, "text");
xlabel.setAttribute("x", xOffset + (dataWidth*scale)/2);
xlabel.setAttribute("y", yOffset + dataHeight*scale + margin.bottom - 7);
xlabel.setAttribute("text-anchor", "middle");
xlabel.setAttribute("font-size", tickFont);
xlabel.textContent = "Principal Component 1";
svg.appendChild(xlabel);

// Y-axis label
const ylabel = document.createElementNS(svgNS, "text");
// move to left of plot and rotate
ylabel.setAttribute(
  "transform",
  `translate(${margin.left - 35}, ${yOffset + (dataHeight*scale)/2}) rotate(-90)`
);
ylabel.setAttribute("text-anchor", "middle");
ylabel.setAttribute("font-size", tickFont);
ylabel.textContent = "Principal Component 2";
svg.appendChild(ylabel);

// 8) Plot points & smart bubbles
const bubbleW = 100, bubbleH = 100, tailLen = 10;
data.forEach(d => {
  const cx = xScale(d.x), cy = yScale(d.y);
  const c = document.createElementNS(svgNS, "circle");
  c.setAttribute("cx", cx);
  c.setAttribute("cy", cy);
  c.setAttribute("r", 5);
  c.setAttribute("fill", "steelblue");
  svg.appendChild(c);

  c.addEventListener("mouseover", () => {
    const g = document.createElementNS(svgNS, "g");
    g.setAttribute("id", "bubble-" + d.id);

    // Determine vertical placement
    let bubbleY;
    if (cy < midY) {
      // point in top half → bubble below
      bubbleY = cy + tailLen;
    } else {
      // bubble above
      bubbleY = cy - bubbleH - tailLen;
    }
    // Horizontal placement
    let bubbleX;
    if (cx > midX) {
      bubbleX = cx - bubbleW - tailLen;
    } else {
      bubbleX = cx + tailLen;
    }

    // Draw white background with rounded corners
    const bg = document.createElementNS(svgNS, "rect");
    bg.setAttribute("x", bubbleX);
    bg.setAttribute("y", bubbleY);
    bg.setAttribute("width", bubbleW);
    bg.setAttribute("height", bubbleH);
    bg.setAttribute("fill", "white");
    bg.setAttribute("rx", cornerRadius);
    bg.setAttribute("ry", cornerRadius);
    g.appendChild(bg);

    // Image sits on top of the background
    const img = document.createElementNS(svgNS, "image");
    img.setAttributeNS("http://www.w3.org/1999/xlink", "href", d.img);
    img.setAttribute("x", bubbleX);
    img.setAttribute("y", bubbleY);
    img.setAttribute("width", bubbleW);
    img.setAttribute("height", bubbleH);
    g.appendChild(img);

    // Border rectangle with same rounding
    const br = document.createElementNS(svgNS, "rect");
    br.setAttribute("x", bubbleX);
    br.setAttribute("y", bubbleY);
    br.setAttribute("width", bubbleW);
    br.setAttribute("height", bubbleH);
    br.setAttribute("fill", "none");
    br.setAttribute("stroke", "black");
    br.setAttribute("rx", cornerRadius);
    br.setAttribute("ry", cornerRadius);
    g.appendChild(br);

    // Tail (vertical)
    const p1 = `${cx},${cy}`;
    let p2, p3;
    if (cy < midY) {
      // tail on top of bubble
      p2 = `${bubbleX + bubbleW/2 - 5},${bubbleY}`;
      p3 = `${bubbleX + bubbleW/2 + 5},${bubbleY}`;
    } else {
      // tail on bottom
      p2 = `${bubbleX + bubbleW/2 - 5},${bubbleY + bubbleH}`;
      p3 = `${bubbleX + bubbleW/2 + 5},${bubbleY + bubbleH}`;
    }
    const tri = document.createElementNS(svgNS, "polygon");
    tri.setAttribute("points", [p1,p2,p3].join(" "));
    tri.setAttribute("fill", "white");
    tri.setAttribute("stroke", "black");
    g.appendChild(tri);

    svg.appendChild(g);
  });
  c.addEventListener("mouseout", () => {
    const old = document.getElementById("bubble-" + d.id);
    if (old) old.remove();
  });
});

// 9) Draw legend
// ——— Legend in top-right of plot area ———
const legend = document.createElementNS(svgNS, "g");
// position it 10px left of the frame’s right edge, 20px down from the top
const legendX = xOffset + dataWidth*scale - 70;
const legendY = yOffset + 20;
legend.setAttribute("id", "legend");

// create a background box
const bgLegend = document.createElementNS(svgNS, "rect");
// position it slightly above/left of your handle
bgLegend.setAttribute("x", legendX - 15);
bgLegend.setAttribute("y", legendY - 12);
bgLegend.setAttribute("width",  80);       // enough to cover circle + label
bgLegend.setAttribute("height", 24);
bgLegend.setAttribute("fill", "none");
bgLegend.setAttribute("stroke", "#ccc");
bgLegend.setAttribute("rx", cornerRadius);
bgLegend.setAttribute("ry", cornerRadius);
legend.appendChild(bgLegend);

// handle: a little circle
const handle = document.createElementNS(svgNS, "circle");
handle.setAttribute("cx", legendX);
handle.setAttribute("cy", legendY);
handle.setAttribute("r", 5);
handle.setAttribute("fill", "steelblue");
legend.appendChild(handle);

// label: “Data”
const lbl = document.createElementNS(svgNS, "text");
lbl.setAttribute("x", legendX + 20);
lbl.setAttribute("y", legendY + 4);        // center vertically on the circle
lbl.setAttribute("font-size", tickFont);
lbl.textContent = "Data";
legend.appendChild(lbl);

svg.appendChild(legend);

</script>
'''

IMG_ROOT = '../assets/figures/05-machine_learning/molecules'

X_pca = np.load('X_pca.npy')
n_imgs, n_components = X_pca.shape
assert n_components == 2

coord_lines = []
coord_lines.append('[')
for i in range(n_imgs):
    x, y = X_pca[i]
    img = f'{IMG_ROOT}/mol_{i:04d}.png'
    coord_lines.append(f'  {{ x: {x:.4f}, y: {y:.4f}, img: "{img}", id: "pt{i+1}" }},')
coord_lines.append(']')

svg = svg_template.replace('/*__DATA_PLACEHOLDER__*/', '\n'.join(coord_lines))

with open('interactive_pca.html', 'w') as f:
    f.write(svg)

