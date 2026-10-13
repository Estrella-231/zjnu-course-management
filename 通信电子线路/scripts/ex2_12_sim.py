# -*- coding: utf-8 -*-
"""习题 2.12  AM接收机中频回路 —— 解析计算 + 数值仿真
教材口径: 2.4.3 式(2-24)(2-27)(2-29)(2-33)(2-35); 2.2.3 式(2-10)(2-11)
题图判读(经像素级连线扫描确认, 原图坐标):
  端子 1  : y=1492 引线 -> x=1476 竖线向下 -> C 上极板(y=1610) -> C -> 底部地线(y=1704)
            => C 跨接在 1 端与地(=3端)之间, 即 C 并联在整个 1-3 绕组上, 全接入
  端子 2  : y=1539 横线(1333~1568), 同时接 I_g、R_g 上端
            => 信号源及其内阻经抽头 2 自耦接入, p1 = N23/N13 = 0.2
  端子 3  : x=1547, y=1586..1705 -> 底部地线
  端子 4  : y=1492 引线 -> R_L 上端 ; R_L 下端 -> 地线
  端子 5  : x=1606, y=1587..1705 -> 地线
  注意: y=1539(端子2横线) 与 x=1476(C/端子1竖线) 为交叉不相连
"""
import math
from PIL import Image, ImageDraw, ImageFont

# ---------------- 已知条件 ----------------
f0 = 455e3
Q0 = 100
N13, N23, N45 = 200, 40, 10
C  = 200e-12
Rs = 20e3          # 图上标 R_g, 信号源内阻
RL = 2e3

w0 = 2*math.pi*f0

# ---------------- 接入系数 ----------------
p1 = N23/N13        # 0.2  信号源/内阻经抽头2自耦接入
p2 = N45/N13        # 0.05 负载经 4-5 变压器接入

# ---------------- 1) 回路电感 L13 ----------------
L13  = 1.0/((w0**2)*C)
wL   = w0*L13

# ---------------- 2) 各电导 ----------------
g0   = 1.0/(w0*L13*Q0)      # 回路损耗(空载)
Rg_p = p1**2/Rs              # R_g 自耦折算到 1-3 端
RL_p = p2**2/RL              # R_L 变压器折算到 1-3 端
g_sum = g0 + Rg_p + RL_p
C_sum = C                    # 无晶体管结电容

# ---------------- 3) 有载品质因数 / 通频带 ----------------
QL  = w0*C_sum/g_sum         # = 1/(w0 L13 g_sum)
Rp  = 1.0/g_sum
BW  = f0/QL                  # 式(2-33)
Kr  = math.sqrt(99.0)        # 式(2-35)

# 空载对照
BW0 = f0/Q0
# 四级(若级联, 仅作对照)
m = 4
BW4 = BW*math.sqrt(2**(1/m)-1)
Kr4 = math.sqrt(100**(1/m)-1)/math.sqrt(2**(1/m)-1)

# ---------------- 仿真: 完整复导纳扫频 ----------------
def Yloop(f):
    return g_sum + 1j*(2*math.pi*f*C_sum - 1.0/(2*math.pi*f*L13))

def Vnorm(f, gs=None):
    """归一化电压传输 |1/(1+jQ_L(f/f0-f0/f))|, 由 y_Sigma/g_sum 给出"""
    g = g_sum if gs is None else gs
    return abs(g/Yloop(f))

def find_edge(fn, target, lo, hi):
    if fn(lo) == target: return lo
    rising = fn(hi) > fn(lo)
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if (fn(mid) < target) == rising: lo = mid
        else: hi = mid
    return 0.5*(lo+hi)

def e07(f): return Vnorm(f)
def e01(f): return Vnorm(f)
Bw07_lo = find_edge(e07, 0.70710678, f0*0.7, f0)
Bw07_hi = find_edge(e07, 0.70710678, f0, f0*1.3)
Bw01_lo = find_edge(e01, 0.1, f0*0.3, f0)
Bw01_hi = find_edge(e01, 0.1, f0, f0*3.0)
simBW = Bw07_hi - Bw07_lo
simKr = (Bw01_hi - Bw01_lo)/simBW

# |Z| 峰
best, sim_f0 = -1, f0
for i in range(4001):
    ff = f0*(0.90 + 0.20*i/4000.0)
    v = abs(1.0/Yloop(ff))
    if v > best: best, sim_f0 = v, ff
simZmax = 1.0/abs(Yloop(sim_f0))
simZp = 1.0/g_sum
# |Z| 的 -3dB(幅值 0.7071) 点
def z7(f): return abs(1.0/Yloop(f))
Zl_lo = find_edge(z7, 0.70710678*simZmax, f0*0.7, f0)
Zl_hi = find_edge(z7, 0.70710678*simZmax, f0, f0*1.3)

# 空载仿真
def Yloop0(f):
    return g0 + 1j*(2*math.pi*f*C_sum - 1.0/(2*math.pi*f*L13))
def z0(f): return abs(1.0/Yloop0(f))
Z0lo = find_edge(z0, 0.70710678*(1.0/g0), f0*0.7, f0)
Z0hi = find_edge(z0, 0.70710678*(1.0/g0), f0, f0*1.3)
simBW0 = Z0hi - Z0lo

# 仿真: 由 |Z| 峰所在频率反算回路电感 (独立的仿真验证)
L13_sim = 1.0/((2*math.pi*sim_f0)**2 * C)

# ---------------- 打印 ----------------
print("=" * 68)
print("接入系数:  p1 = N23/N13 = %.3f (信号源经抽头2自耦接入)" % p1)
print("           p2 = N45/N13 = %.3f (R_L 经 4-5 变压器接入)" % p2)
print("L13 = 1/((2*pi*f0)^2 C) = %.4f uH   (w0*L13 = %.3f ohm)" % (L13*1e6, wL))
print("-" * 68)
print("g0 = 1/(2 pi f0 L13 Q0) = %.6e S   (R_p0 = %.2f kohm)" % (g0, 1/g0/1e3))
print("p1^2/Rs = %.6e S   (R_g 折算 = %.1f kohm)" % (Rg_p, 1/Rg_p/1e3))
print("p2^2/RL = %.6e S   (R_L 折算 = %.1f kohm)" % (RL_p, 1/RL_p/1e3))
print("g_Sigma = %.6e S      R_p = %.2f kohm" % (g_sum, Rp/1e3))
print("-" * 68)
print("Q_L = %.4f      2df0.7 = %.6f kHz   K_r0.1 = %.4f" % (QL, BW/1e3, Kr))
print("空载对照: Q0 = %.1f   2df0.7(空载) = %.6f kHz" % (Q0, BW0/1e3))
print("四级级联对照: 2df0.7 = %.4f kHz   K_r0.1 = %.4f" % (BW4/1e3, Kr4))
print("-" * 68)
print("仿真 f0      = %.6f kHz" % (sim_f0/1e3))
print("仿真反算 L13 = %.4f uH   (偏差 %.4f%%)" % (L13_sim*1e6, 100*(L13_sim-L13)/L13))
print("仿真 |Z|max  = %.3f ohm  (解析 1/g_Sigma = %.3f)" % (simZmax, simZp))
print("仿真 2df0.7  = %.6f kHz   偏差 %.4f%%" % (simBW/1e3, 100*(simBW-BW)/BW))
print("仿真 K_r0.1  = %.4f        偏差 %.4f%%" % (simKr, 100*(simKr-Kr)/Kr))
print("仿真 2df0.7|Z| = %.6f kHz 偏差 %.4f%%" % ((Zl_hi-Zl_lo)/1e3, 100*((Zl_hi-Zl_lo)-BW)/BW))
print("仿真 2df0.7(空载) = %.6f kHz 偏差 %.4f%%" % (simBW0/1e3, 100*(simBW0-BW0)/BW0))
print("=" * 68)

# ---------------- 绘图 ----------------
FN = r"C:\Windows\Fonts\msyh.ttc"
FB = r"C:\Windows\Fonts\msyhbd.ttc"
W, H = 1560, 1230
img = Image.new("RGB", (W, H), "white")
dr  = ImageDraw.Draw(img)
f18 = ImageFont.truetype(FB, 19)
f16 = ImageFont.truetype(FN, 16)
f14 = ImageFont.truetype(FN, 14)
f22 = ImageFont.truetype(FB, 22)
BLUE, RED, GRAY = (20, 70, 190), (200, 30, 30), (110, 110, 110)
ORNG, GREEN = (215, 130, 0), (0, 140, 60)

def panel(x0, y0, x1, y1, title):
    dr.rectangle([x0, y0, x1, y1], outline=(60, 60, 60), width=2)
    dr.text((x0 + 10, y0 - 26), title, font=f18, fill=(0, 0, 0))
    return (x0, y0, x1, y1)

def axes(p, xlab, ylab, xlim, ylim):
    x0, y0, x1, y1 = p
    for i in range(5):
        yy = y0 + (y1-y0)*i/4.0
        dr.line([x0, yy, x1, yy], fill=(228, 228, 228))
        dr.line([x0, yy-9, x0, yy], fill=(60, 60, 60))
        dr.text((x0 - 74, yy - 9), ("%.3g" % (ylim[1] - (ylim[1]-ylim[0])*i/4.0)), font=f14, fill=GRAY)
    for i in range(7):
        xx = x0 + (x1-x0)*i/6.0
        dr.line([xx, y0, xx, y1], fill=(238, 238, 238))
        dr.line([xx, y1, xx, y1+9], fill=(60, 60, 60))
        dr.text((xx - 24, y1 + 13), ("%.4g" % (xlim[0] + (xlim[1]-xlim[0])*i/6.0)), font=f14, fill=GRAY)
    dr.line([x0, y1, x1, y1], fill=(60, 60, 60), width=2)
    dr.line([x0, y0, x0, y1], fill=(60, 60, 60), width=2)
    dr.text((x1 - 118, y1 + 34), xlab, font=f14, fill=(0, 0, 0))
    dr.text((x0 + 12, y0 + 8), ylab, font=f14, fill=(130, 130, 130))

def plot(p, xs, ys, xlim, ylim, color, width=2):
    x0, y0, x1, y1 = p
    pts = []
    for xv, yv in zip(xs, ys):
        px = x0 + (x1-x0)*(xv-xlim[0])/(xlim[1]-xlim[0])
        py = y1 - (y1-y0)*(yv-ylim[0])/(ylim[1]-ylim[0])
        if not (x0 <= px <= x1 and y0 <= py <= y1): continue
        pts.append((px, py))
    if len(pts) > 1: dr.line(pts, fill=color, width=width)

def vline(p, fv, xlim, color, txt=None, dy=6, dxa=None):
    x0, y0, x1, y1 = p
    px = x0 + (x1-x0)*(fv-xlim[0])/(xlim[1]-xlim[0])
    if x0 <= px <= x1: dr.line([px, y0, px, y1], fill=color, width=1)
    if txt:
        if dxa == "right": dr.text((px - 8 - dr.textlength(txt, font=f14), y0+dy), txt, font=f14, fill=color)
        else: dr.text((px + 6, y0 + dy), txt, font=f14, fill=color)

def hline(p, v, ylim, color, txt=None, dx=96, dy=-20):
    x0, y0, x1, y1 = p
    py = y1 - (y1-y0)*(v-ylim[0])/(ylim[1]-ylim[0])
    if y0 <= py <= y1: dr.line([x0, py, x1, py], fill=color, width=1)
    if txt: dr.text((x1 - dx - 6, py + dy), txt, font=f14, fill=color)

def legend(p, items):
    x0, y0, x1, y1 = p
    yy = y0 + 10
    for col, txt in items:
        dr.line([x1 - 250, yy, x1 - 215, yy], fill=col, width=3)
        dr.text((x1 - 208, yy - 9), txt, font=f14, fill=(0, 0, 0))
        yy += 22

dr.text((24, 14), "习题 2.12  AM 接收机中频回路 —— 解析解与数值仿真对照", font=f22, fill=(0, 0, 0))

FL, FH = 425.0, 485.0   # kHz
xs = [FL + (FH-FL)*i/1200.0 for i in range(1201)]

# ===== 面板1: 归一化谐振曲线 =====
PA = panel(80, 92, 770, 500, "面板 1  谐振曲线  |u|/|u|max   (教材式 2-30)")
ylim1 = (0, 1.12)
axes(PA, "f / kHz", "归一化电压传输", (FL, FH), ylim1)
plot(PA, xs, [Vnorm(f*1e3) for f in xs], (FL, FH), ylim1, BLUE, 3)
hline(PA, 0.70710678, ylim1, GRAY, "0.707 线")
vline(PA, f0/1e3, (FL, FH), GREEN, "f0 = 455 kHz")
vline(PA, Bw07_lo/1e3, (FL, FH), RED)
vline(PA, Bw07_hi/1e3, (FL, FH), RED)
dr.text((PA[0]+8, PA[1]+42), "%.3f kHz" % (Bw07_lo/1e3), font=f14, fill=RED)
dr.text((PA[0]+8, PA[1]+64), "%.3f kHz" % (Bw07_hi/1e3), font=f14, fill=RED)
legend(PA, [(BLUE, "有载 Q_L = %.2f" % QL), (GRAY, "-3 dB 边")])
dr.text((PA[0]+10, PA[3]+64),
        "Q_L = %.4f      2Δf0.7 = %.4f kHz      K_r0.1 = %.4f ( = sqrt99 )"
        % (QL, BW/1e3, Kr), font=f16, fill=(0, 0, 0))

# ===== 面板2: |Z| 曲线, 有载 vs 空载 =====
PB = panel(840, 92, 1530, 500, "面板 2  回路阻抗幅频特性  |Z|  有载 vs 空载")
ymax2 = max(1/g_sum, 1/g0)*1.15
ylim2 = (0, ymax2)
axes(PB, "f / kHz", "|Z| / ohm", (FL, FH), ylim2)
plot(PB, xs, [abs(1.0/Yloop(f*1e3)) for f in xs], (FL, FH), ylim2, RED, 3)
plot(PB, xs, [abs(1.0/Yloop0(f*1e3)) for f in xs], (FL, FH), ylim2, GRAY, 2)
vline(PB, f0/1e3, (FL, FH), GREEN, "f0")
vline(PB, Zl_lo/1e3, (FL, FH), RED)
vline(PB, Zl_hi/1e3, (FL, FH), RED)
dr.text((PB[0]+8, PB[1]+42), "%.3f kHz" % (Zl_lo/1e3), font=f14, fill=RED)
dr.text((PB[0]+8, PB[1]+64), "%.3f kHz" % (Zl_hi/1e3), font=f14, fill=RED)
legend(PB, [(RED, "有载 |Z|max = %.0f ohm" % simZp), (GRAY, "空载 Q0 = %d" % Q0)])
dr.text((PB[0]+10, PB[3]+64),
        "有载 2Δf0.7 = %.4f kHz   空载 2Δf0.7 = %.4f kHz   (展宽 %.2f 倍)"
        % (BW/1e3, BW0/1e3, BW/BW0), font=f16, fill=(0, 0, 0))

# ===== 面板3: 电导构成对 Q_L 的影响 =====
PC = panel(80, 640, 770, 1020, "面板 3  g_Σ 各项贡献与有载品质因数")
items = [("g_0 (回路损耗)", g0), ("p1^2/R_s (信号源内阻)", Rg_p), ("p2^2/R_L (负载)", RL_p)]
bx = PC[0] + 40
bw = 560
by = PC[1] + 56
for i, (nm, v) in enumerate(items):
    yy = by + i*70
    frac = v/g_sum
    dr.text((PC[0]+16, yy-22), "%d.  %s = %.4e S   (%.2f%%)" % (i+1, nm, v, frac*100),
            font=f14, fill=(0, 0, 0))
    dr.rectangle([bx, yy, bx+bw, yy+26], outline=(120, 120, 120), width=1)
    dr.rectangle([bx, yy, bx+bw*frac, yy+26], fill=(60, 110, 200))
    dr.text((bx+bw*frac+8, yy+4), "%.1f%%" % (frac*100), font=f14, fill=(60, 60, 60))
yy = by + 3*70
dr.rectangle([bx, yy, bx+bw, yy+26], outline=(180, 30, 30), width=2)
dr.rectangle([bx, yy, bx+bw, yy+26], fill=(200, 60, 60))
dr.text((PC[0]+16, yy-22), "g_Σ = %.6e S   (R_p = %.2f kohm)" % (g_sum, Rp/1e3), font=f14, fill=(180, 30, 30))
dr.text((bx+8, yy+4), "100%", font=f14, fill=(255, 255, 255))
dr.text((PC[0]+16, yy+44),
        "Q_L = w0 C / g_Σ = %.4f ;  若只有 g_0 则 Q_L = Q0 = %d" % (QL, Q0),
        font=f16, fill=(0, 0, 0))

# ===== 面板4: 解析 vs 仿真 =====
PD = panel(840, 640, 1530, 1020, "面板 4  解析值与仿真值对照")
rows = [
    ("回路电感 L13",   "%.4f uH" % (L13*1e6), "%.4f uH" % (L13_sim*1e6), "%+.4f%%" % (100*(L13_sim-L13)/L13)),
    ("谐振频率 f0",    "455.0000 kHz",  "%.4f kHz" % (sim_f0/1e3), "%+.4f%%" % (100*(sim_f0-f0)/f0)),
    ("回路并联谐振阻 Rp","%.3f ohm" % Rp, "%.3f ohm" % simZp, "%+.4f%%" % (100*(simZp-Rp)/Rp)),
    ("2Δf0.7 (归一化)", "%.6f kHz" % (BW/1e3), "%.6f kHz" % (simBW/1e3), "%+.4f%%" % (100*(simBW-BW)/BW)),
    ("2Δf0.7 (|Z|)",   "%.6f kHz" % (BW/1e3), "%.6f kHz" % ((Zl_hi-Zl_lo)/1e3), "%+.4f%%" % (100*((Zl_hi-Zl_lo)-BW)/BW)),
    ("K_r0.1",         "%.4f" % Kr,     "%.4f" % simKr,   "%+.4f%%" % (100*(simKr-Kr)/Kr)),
    ("2Δf0.7 空载",    "%.6f kHz" % (BW0/1e3), "%.6f kHz" % (simBW0/1e3), "%+.4f%%" % (100*(simBW0-BW0)/BW0)),
]
yy = PD[1] + 16
dr.line([PD[0]+8, yy-6, PD[2]-8, yy-6], fill=(60, 60, 60), width=1)
for a, b, c, d in [("项目", "解析值", "仿真值", "偏差")] + rows:
    fnt = f16 if a == "项目" else f14
    dr.text((PD[0]+16,  yy), a, font=fnt, fill=(0, 0, 0))
    dr.text((PD[0]+215, yy), b, font=fnt, fill=(30, 30, 30))
    dr.text((PD[0]+400, yy), c, font=fnt, fill=(30, 30, 30))
    dr.text((PD[0]+570, yy), d, font=fnt, fill=(0, 110, 0))
    yy += 32
dr.line([PD[0]+8, yy-6, PD[2]-8, yy-6], fill=(60, 60, 60), width=1)
dr.multiline_text((PD[0]+16, yy + 18),
    "接入系数  p1 = N23/N13 = %.2f (自耦)    p2 = N45/N13 = %.2f (变压器)\n"
    "折算负载  R_s' = (N13/N23)^2 R_s = %.0f kohm    R_L' = (N13/N45)^2 R_L = %.0f kohm\n"
    "矩形系数  K_r0.1 = sqrt99 = %.4f, 与 Q_L 无关(单调谐固有值)" % (p1, p2, 1/Rg_p/1e3, 1/RL_p/1e3, Kr),
    font=f16, fill=(0, 0, 0), spacing=12)

txt = ("L13 = 1/[(2πf0)^2 C] = %.2f uH        ω0 L13 = %.2f Ω        R_p0 = Q0 ω0 L13 = %.2f kΩ\n"
       "g0 = 1/(2πf0 L13 Q0) = %.4e S   p1^2/R_s = %.4e S   p2^2/R_L = %.4e S\n"
       "gΣ = %.4e S   R_p = %.2f kΩ   Q_L = %.4f   2Δf0.7 = %.4f kHz   K_r0.1 = %.4f"
       % (L13*1e6, wL, 1/g0/1e3, g0, Rg_p, RL_p, g_sum, Rp/1e3, QL, BW/1e3, Kr))
dr.multiline_text((80, 1108), txt, font=f16, fill=(0, 0, 0), spacing=10)

img.save(r"E:\tmp\通信电路\ex2-12.png")
cnt = sum(1 for py in range(H) for pxx in range(W) if img.getpixel((pxx, py)) != (255, 255, 255))
print("曲线像素自检: 非白像素 %d 个 (应 > 20000)" % cnt)
print("已保存 E:\\tmp\\通信电路\\ex2-12.png  (%dx%d)" % (W, H))