# -*- coding: utf-8 -*-
"""习题 2.10  单调谐放大器 —— 解析计算 + 数值仿真
教材口径: 2.4.3 式(2-24)(2-28)(2-33)(2-34)(2-35); 2.4.4 式(2-38)(2-39)(2-40);
          2.5.2 式(2-41)(2-52)
题图判读(已按像素级连线扫描复核修正):
  一次绕组 1-3 全绕组 L13=4uH/Q0=100;
  C 跨接在 1-3 之间 (上极板->端子1节点(+V_CC), 下极板->y=484 横线->端子3);
  10 kΩ 也跨接在 1-3 之间 (上端->黑点即端子1/+V_CC, 下端 x=944,y=486->同一条 y=484 横线);
  VT1 集电极接 2 抽头 (x=966,422 -> x=967,498 -> x=919,497..524 -> VT1 集电极) -> p1=N12/N13=0.7;
  二次绕组 4-5 (N45=5): 4 -> 0.01uF -> VT2 基极, 5 -> 地母线, p2=N45/N13=0.25
注: 扫描时 x=847,y=389..487 的"竖线"实为 "10kΩ" 标注文字的笔画, 勿当作导线
"""
import cmath, math
from PIL import Image, ImageDraw, ImageFont

# ---------------- 已知条件 ----------------
f0   = 10.7e6
w0   = 2 * math.pi * f0
L13  = 4e-6
Q0   = 100
N13, N23, N45 = 20, 6, 5
R10k = 10e3

y_ie = complex(2.86e-3, 3.4e-3)
y_oe = complex(0.2e-3, 1.3e-3)
y_fe = complex(26.4e-3, -36.4e-3)
y_re = complex(0.08e-3, -0.3e-3)

g_ie, C_ie = y_ie.real, y_ie.imag / w0
g_oe, C_oe = y_oe.real, y_oe.imag / w0
F_abs = abs(y_fe)

# 基极偏置 15k // 6.2k
g_s = 1/15e3 + 1/6.2e3

# ---------------- 接入系数 ----------------
p1  = (N13 - N23) / N13       # 0.7   集电极接 2 抽头, p1 = N12/N13
p2  = N45 / N13               # 0.25  输出绕组 4-5

g0 = 1.0 / (2 * math.pi * f0 * L13 * Q0)
g_R = 1.0 / R10k               # 10 kΩ 并联在 1-3 之间(教材/官方口径不计入 g_Sigma, 见下方对照)

g_sum = p1**2*g_oe + g0 + p2**2*g_ie                    # 教材式(2-24) 原式
C_sum = 1.0 / ((2*math.pi*f0)**2 * L13)
C     = C_sum - p1**2*C_oe - p2**2*C_ie
Q_L   = 1.0 / (2*math.pi*f0*L13*g_sum)
Au0   = p1*p2*F_abs/g_sum
BW    = f0 / Q_L                                                 # 式(2-33)
Kr    = math.sqrt(99.0)                                          # 式(2-35)

# 四级
m = 4
Au4    = Au0**m                                                 # 式(2-37)
BW4    = BW * math.sqrt(2**(1/m) - 1)                           # 式(2-39)
Kr4    = math.sqrt(100**(1/m) - 1)/math.sqrt(2**(1/m) - 1)       # 式(2-40)

# ---------------- (2) 稳定性 ----------------
# 折算到 1-3 端的负载复导纳 Y'_L (教材口径: 不含 10 kΩ)
Yp = g0 + p2**2*y_ie
# 式(2-41) 原始复数形式 (精确, 含偏置电阻 g_s)
S_cplx = (y_ie + g_s) * (y_oe + Yp) / (y_re * y_fe)
S_abs  = abs(S_cplx)
# 教材 2.5.2 近似式(2-44) 在 f=f0 取 xi=0
S_apx = (g_ie + g_s)*(g_oe + Yp.real)/(abs(y_re)*F_abs)
# 稳定电压增益 式(2-52):  官方口径 g1≈g_ie, C_re 取 |Im y_re|/w0 (即 y_re≈-j w0 Cre)
Cre_o   = abs(y_re.imag)/w0
Au0S    = 2*g_ie / (5*abs(y_re.imag))            # = 官方 3.82
# 精确口径: 计入偏置电阻 g_s, C_re 取 |y_re|/w0
Cre_p   = abs(y_re)/w0
Au0S_p  = 2*(g_ie + g_s) / (5*abs(y_re))

# ---------------- 对照: 若把并联在 1-3 的 10 kΩ 也计入 g_Sigma ----------------
g_sum_R = g_sum + g_R
Au0_R   = p1*p2*F_abs/g_sum_R
QL_R    = 1.0/(2*math.pi*f0*L13*g_sum_R)
BW_R    = f0/QL_R
Au4_R   = Au0_R**m
BW4_R   = BW_R*math.sqrt(2**(1/m)-1)

# ---------------- 仿真: 用完整复导纳直接扫频 ----------------
def Yloop(f):
    """回路总复导纳 y_Sigma = g_Sigma + j(w C_Sigma - 1/(w L13))"""
    return g_sum + 1j*(2*math.pi*f*C_sum - 1.0/(2*math.pi*f*L13))

def A(f, n=1):
    """式(2-27) 精确式 A_u = -p1 p2 y_fe / y_Sigma"""
    a = -p1*p2*y_fe / Yloop(f)
    return a**n

def find_edge(fn, target, lo, hi):
    """在 [lo,hi] 内二分找 fn(f)==target, 自动判定单调方向"""
    if fn(lo) == target: return lo
    rising = fn(hi) > fn(lo)
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if (fn(mid) < target) == rising: lo = mid
        else: hi = mid
    return 0.5*(lo+hi)

def r1(f): return abs(A(f)) / abs(Au0)
def r4(f): return abs(A(f, 4)) / abs(Au4)

# 单级 -3dB / -20dB
Bw07_lo = find_edge(r1, 0.70710678, f0*0.5, f0)
Bw07_hi = find_edge(r1, 0.70710678, f0, f0*2)
Bw01_lo = find_edge(r1, 0.1, f0*0.5, f0)
Bw01_hi = find_edge(r1, 0.1, f0, f0*2)
simBW    = Bw07_hi - Bw07_lo
simKr    = (Bw01_hi - Bw01_lo)/simBW
# 四级
Bw4_lo   = find_edge(r4, 0.70710678, f0*0.5, f0)
Bw4_hi   = find_edge(r4, 0.70710678, f0, f0*2)
Bw4_01_lo = find_edge(r4, 0.1, f0*0.5, f0)
Bw4_01_hi = find_edge(r4, 0.1, f0, f0*2)
simBW4 = Bw4_hi - Bw4_lo
simKr4 = (Bw4_01_hi - Bw4_01_lo)/simBW4
# 谐振峰所在频率(扫描找最大)
best, sim_f0 = -1, f0
for i in range(2001):
    ff = f0*(0.97 + 0.06*i/2000.0)
    v = abs(A(ff))
    if v > best: best, sim_f0 = v, ff
# 峰值与解析 A_u0 比较
simAu0 = best

# ---------------- 打印 ----------------
print("=" * 70)
print("接入系数: p1 = N12/N13 = %.4f (集电极接 2 抽头)   p2 = N45/N13 = %.4f" % (p1, p2))
print("g_ie = %.4f mS  C_ie = %.3f pF   g_oe = %.4f mS  C_oe = %.3f pF"
      % (g_ie*1e3, C_ie*1e12, g_oe*1e3, C_oe*1e12))
print("|y_fe| = %.4f mS   |y_re| = %.5f mS   C_re(官方) = %.3f pF"
      % (F_abs*1e3, abs(y_re)*1e3, Cre_o*1e12))
print("g_0 = %.6f mS      g_s = %.5f mS" % (g0*1e3, g_s*1e3))
print("-" * 70)
print("(1) g_Sigma = %.6f mS   [官方 3.14e-4 S]" % (g_sum*1e3))
print("(1) C_Sigma = %.4f pF ->  C = %.4f pF   [官方 55.37 / 42.73 pF]" % (C_sum*1e12, C*1e12))
print("(1) A_u0 = %.4f   [官方 25.06]      Q_L = %.4f   [官方 11.85]" % (Au0, Q_L))
print("(1) 2df0.7 = %.6f MHz [官方 0.90]   K_r0.1 = %.4f" % (BW/1e6, Kr))
print("(4) A_u4 = %.4e   [官方 3.94e5]" % Au4)
print("(4) 2df0.7 = %.6f MHz [官方 0.39]   K_r0.1 = %.4f  [官方 3.38]" % (BW4/1e6, Kr4))
print("-" * 70)
print("仿真 f0 = %.6f MHz" % (sim_f0/1e6))
print("仿真 2df0.7(单) = %.6f MHz   偏差 %.4f%%" % (simBW/1e6, 100*(simBW-BW)/BW))
print("仿真 K_r0.1(单) = %.4f        偏差 %.4f%%" % (simKr, 100*(simKr-Kr)/Kr))
print("仿真 2df0.7(四级) = %.6f MHz 偏差 %.4f%%" % (simBW4/1e6, 100*(simBW4-BW4)/BW4))
print("仿真 K_r0.1(四级) = %.4f       偏差 %.4f%%" % (simKr4, 100*(simKr4-Kr4)/Kr4))
print("-" * 70)
print("(2) S (式2-41 精确) = %.4f      S (式2-44 近似) = %.4f" % (S_abs, S_apx))
print("(2) 稳定电压增益 A_u0S = %.4f [官方 3.82]   精确口径(含 g_s)= %.4f" % (Au0S, Au0S_p))
print("(2) 实际 A_u0 = %.4f  ->  %s" % (Au0, "稳定" if Au0 <= Au0S else "不稳定(会自激)"))
print("-" * 70)
print("对照(若把并联在 1-3 的 10kΩ 计入 g_Sigma):")
print("   g_Sigma = %.6f mS   A_u0 = %.4f   Q_L = %.4f   2df0.7 = %.4f MHz"
      % (g_sum_R*1e3, Au0_R, QL_R, BW_R/1e6))
print("   A_u4 = %.4e   2df0.7(4) = %.4f MHz" % (Au4_R, BW4_R/1e6))
print("=" * 70)

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
        v = ylim[1] - (ylim[1]-ylim[0])*i/4.0
        dr.text((x0 - 66, yy - 9), ("%.3g" % v), font=f14, fill=GRAY)
    for i in range(7):
        xx = x0 + (x1-x0)*i/6.0
        dr.line([xx, y0, xx, y1], fill=(238, 238, 238))
        dr.line([xx, y1, xx, y1+9], fill=(60, 60, 60))
        v = xlim[0] + (xlim[1]-xlim[0])*i/6.0
        dr.text((xx - 22, y1 + 13), ("%.4g" % v), font=f14, fill=GRAY)
    dr.line([x0, y1, x1, y1], fill=(60, 60, 60), width=2)
    dr.line([x0, y0, x0, y1], fill=(60, 60, 60), width=2)
    dr.text((x1 - 108, y1 + 34), xlab, font=f14, fill=(0, 0, 0))
    dr.text((x0 + 12, y0 + 8), ylab, font=f14, fill=(130, 130, 130))

def plot(p, xs, ys, xlim, ylim, color, width=2):
    x0, y0, x1, y1 = p
    pts = []
    for xv, yv in zip(xs, ys):
        px = x0 + (x1-x0)*(xv - xlim[0])/(xlim[1]-xlim[0])
        py = y1 - (y1-y0)*(yv - ylim[0])/(ylim[1]-ylim[0])
        if not (x0 <= px <= x1 and y0 <= py <= y1): continue
        if pts and abs(pts[-1][0]-px) > (x1-x0) and abs(pts[-1][1]-py) > (y1-y0):
            continue
        pts.append((px, py))
    if len(pts) > 1: dr.line(pts, fill=color, width=width)

def vline(p, fv, xlim, color, txt=None, dy=6, dxa=None):
    x0, y0, x1, y1 = p
    px = x0 + (x1-x0)*(fv - xlim[0])/(xlim[1]-xlim[0])
    if x0 <= px <= x1: dr.line([px, y0, px, y1], fill=color, width=1)
    if txt:
        if dxa == "right": dr.text((px - 8 - dr.textlength(txt, font=f14), y0 + dy), txt, font=f14, fill=color)
        else: dr.text((px + 6, y0 + dy), txt, font=f14, fill=color)

def hline(p, v, ylim, color, txt=None, dx=96, dy=-20):
    x0, y0, x1, y1 = p
    py = y1 - (y1-y0)*(v - ylim[0])/(ylim[1]-ylim[0])
    if y0 <= py <= y1: dr.line([x0, py, x1, py], fill=color, width=1)
    if txt:
        dr.text((x1 - dx - 6, py + dy), txt, font=f14, fill=color)

def legend(p, items, anchor="tr"):
    x0, y0, x1, y1 = p
    yy = y0 + 10
    for col, txt in items:
        dr.line([x1 - 210, yy, x1 - 175, yy], fill=col, width=3)
        dr.text((x1 - 168, yy - 9), txt, font=f14, fill=(0, 0, 0))
        yy += 22

# ===== 标题 =====
dr.text((24, 14), "习题 2.10  单调谐放大器 —— 解析解与数值仿真对照", font=f22, fill=(0, 0, 0))

FL, FH = 9.2, 12.2   # MHz 坐标范围
xs_mhz = [FL + (FH-FL)*i/1200.0 for i in range(1201)]
xs = [v*1e6 for v in xs_mhz]

# ===== 面板1: 单级谐振曲线 =====
PA = panel(80, 92, 770, 500, "面板 1  单级放大器谐振曲线  A_u/A_u0   (教材式 2-30)")
ylim1 = (0, 1.12)
axes(PA, "f / MHz", "|A_u| / |A_u0|", (FL, FH), ylim1)
ys1 = [abs(A(f))/abs(Au0) for f in xs]
plot(PA, xs_mhz, ys1, (FL, FH), ylim1, BLUE, 3)
hline(PA, 0.70710678, ylim1, GRAY, "0.707 线")
vline(PA, f0/1e6, (FL, FH), GREEN, "f0 = 10.7")
vline(PA, Bw07_lo/1e6, (FL, FH), RED)
vline(PA, Bw07_hi/1e6, (FL, FH), RED)
dr.text((PA[0]+8, PA[1]+40), "%.3f MHz" % Bw07_lo, font=f14, fill=RED)
dr.text((PA[0]+8, PA[1]+62), "%.3f MHz" % Bw07_hi, font=f14, fill=RED)
legend(PA, [(BLUE, "单级 |A_u|"), (RED, "-3 dB 边")])
dr.text((PA[0]+10, PA[3]+64),
        "Q_L = %.4f     2Δf0.7 = %.4f MHz     K_r0.1 = %.4f ( = sqrt99 )"
        % (Q_L, BW/1e6, Kr), font=f16, fill=(0, 0, 0))

# ===== 面板2: 四级级联 =====
PB = panel(840, 92, 1530, 500, "面板 2  单级与四级级联谐振曲线对比   (教材式 2-38)")
axes(PB, "f / MHz", "|A_m| / |A_m0|", (FL, FH), ylim1)
plot(PB, xs_mhz, [abs(A(f, 4))/abs(Au4) for f in xs], (FL, FH), ylim1, RED, 3)
plot(PB, xs_mhz, [abs(A(f, 1))/abs(Au0) for f in xs], (FL, FH), ylim1, GRAY, 2)
hline(PB, 0.70710678, ylim1, GRAY, "0.707 线")
vline(PB, f0/1e6, (FL, FH), GREEN, "f0")
vline(PB, Bw4_lo/1e6, (FL, FH), RED)
vline(PB, Bw4_hi/1e6, (FL, FH), RED)
dr.text((PB[0]+8, PB[1]+40), "%.3f MHz" % Bw4_lo, font=f14, fill=RED)
dr.text((PB[0]+8, PB[1]+62), "%.3f MHz" % Bw4_hi, font=f14, fill=RED)
legend(PB, [(RED, "四级 m = 4"), (GRAY, "单级 m = 1")])
dr.text((PB[0]+10, PB[3]+64),
        "A_u4 = %.4e     2Δf0.7(4) = %.4f MHz     K_r0.1(4) = %.4f"
        % (Au4, BW4/1e6, Kr4), font=f16, fill=(0, 0, 0))

# ===== 面板3: 稳定性判据 =====
PC = panel(80, 640, 770, 1020, "面板 3  稳定性判据:  实际 A_u0 与稳定电压增益 A_u0S 的关系")
g2r = [0.10 + i*0.02 for i in range(150)]
ys_A = [p1*p2*F_abs/(g*1e-3) for g in g2r]
xlim3, ylim3 = (0, 3.0), (0, 120)
axes(PC, "回路总电导 g_Σ = p1²g_oe + g_0 + p2²g_ie  /mS", "A_u0", xlim3, ylim3)
plot(PC, g2r, ys_A, xlim3, ylim3, BLUE, 3)
hline(PC, Au0S, ylim3, ORNG, "A_u0S = %.3f" % Au0S, dx=110, dy=-22)
# 稳定临界 g_Sigma
g2crit = p1*p2*F_abs/Au0S*1e3
vline(PC, g2crit, xlim3, RED, "g_Σ 临界 %.3f mS" % g2crit, dy=40, dxa="right")
# 实际工作点
g2act = g_sum*1e3
vline(PC, g2act, xlim3, GREEN, "实际 g2 = %.4f mS" % g2act, dy=62)
dr.ellipse([PC[0] + (PC[2]-PC[0])*g2act/xlim3[1] - 5,
            PC[3] - (PC[3]-PC[1])*Au0/ylim3[1] - 5,
            PC[0] + (PC[2]-PC[0])*g2act/xlim3[1] + 5,
            PC[3] - (PC[3]-PC[1])*Au0/ylim3[1] + 5], outline=GREEN, width=3)
dr.text((PC[0]+180, PC[1]+8), "实际点  g_Σ = %.4f mS,  A_u0 = %.2f" % (g2act, Au0),
        font=f14, fill=GREEN)
legend(PC, [(BLUE, "A_u0 = p1 p2 |y_fe| / g2")])
dr.text((PC[0]+10, PC[3]+64),
        "稳定电压增益 A_u0S = 2g1/(5 w0 Cre) = %.4f ；  实际 A_u0 = %.4f  >  A_u0S"
        % (Au0S, Au0), font=f16, fill=(180, 0, 0))

# ===== 面板4: 解析 vs 仿真 对照数值 =====
PD = panel(840, 640, 1530, 1020, "面板 4  解析值与仿真值对照")
rows = [
    ("谐振频率 f0",  "10.700000 MHz",  "%.6f MHz" % (sim_f0/1e6),  "%+.4f%%" % (100*(sim_f0-f0)/f0)),
    ("单级 2Δf0.7",  "%.6f MHz" % (BW/1e6),  "%.6f MHz" % (simBW/1e6),  "%+.4f%%" % (100*(simBW-BW)/BW)),
    ("单级 K_r0.1",  "%.4f" % Kr,       "%.4f" % simKr,             "%+.4f%%" % (100*(simKr-Kr)/Kr)),
    ("四级 2Δf0.7",  "%.6f MHz" % (BW4/1e6), "%.6f MHz" % (simBW4/1e6), "%+.4f%%" % (100*(simBW4-BW4)/BW4)),
    ("四级 K_r0.1",  "%.4f" % Kr4,      "%.4f" % simKr4,            "%+.4f%%" % (100*(simKr4-Kr4)/Kr4)),
    ("单级 A_u0",    "%.4f" % Au0,      "%.4f" % simAu0,            "%+.4f%%" % (100*(simAu0-Au0)/Au0)),
]
yy = PD[1] + 16
dr.line([PD[0]+8, yy-6, PD[2]-8, yy-6], fill=(60, 60, 60), width=1)
for a, b, c, d in [("项目", "解析值", "仿真值", "偏差")] + rows:
    fnt = f16 if a == "项目" else f14
    dr.text((PD[0]+16,  yy), a, font=fnt, fill=(0, 0, 0))
    dr.text((PD[0]+200, yy), b, font=fnt, fill=(30, 30, 30))
    dr.text((PD[0]+380, yy), c, font=fnt, fill=(30, 30, 30))
    dr.text((PD[0]+560, yy), d, font=fnt, fill=(0, 110, 0))
    yy += 32
dr.line([PD[0]+8, yy-6, PD[2]-8, yy-6], fill=(60, 60, 60), width=1)
dr.multiline_text((PD[0]+16, yy + 18),
    "稳定性:  S(式2-41 精确) = %.4f     S(式2-44 近似) = %.4f\n"
    "          稳定电压增益 A_u0S = %.4f  (官方 3.82)   实际 A_u0 = %.4f\n"
    "          A_u0 / A_u0S = %.2f  > 1   =>   不稳定"
    % (S_abs, S_apx, Au0S, Au0, Au0/Au0S), font=f16, fill=(170, 0, 0), spacing=12)

txt = ("关键中间量:   g_Σ = %.6f mS      C_Σ = %.4f pF      C = %.4f pF\n"
       "              g_0 = %.6f mS     Q_L = %.4f     |y_re| = %.5f mS\n"
       "              接入系数  p1 = N12/N13 = %.2f (集电极接 2 抽头)   p2 = N45/N13 = %.2f\n"
       "对照: 若把并联在 1-3 的 10kΩ 计入 g_Σ ->  g_Σ = %.4f mS, A_u0 = %.2f, 2Δf0.7 = %.3f MHz"
       % (g_sum*1e3, C_sum*1e12, C*1e12, g0*1e3, Q_L, abs(y_re)*1e3, p1, p2,
          g_sum_R*1e3, Au0_R, BW_R/1e6))
dr.multiline_text((80, 1108), txt, font=f16, fill=(0, 0, 0), spacing=10)

img.save(r"E:\tmp\通信电路\ex2-10.png")

# ---- 曲线自检 ----
cnt = sum(1 for py in range(H) for pxx in range(W)
          if img.getpixel((pxx, py)) != (255, 255, 255))
print("曲线像素自检: 非白像素 %d 个 (应 > 20000)" % cnt)
print("已保存 E:\\tmp\\通信电路\\ex2-10.png  (%dx%d)" % (W, H))
