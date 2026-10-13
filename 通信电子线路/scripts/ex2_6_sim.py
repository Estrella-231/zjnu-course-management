# -*- coding: utf-8 -*-
"""教材习题 2.6  复杂并联谐振回路

题2.6图拓扑（全部并在上下母线之间）：
    R_i  ∥  C_i  ∥  L(Q0)  ∥  [ C1 串联 (C2 ∥ R_L) ]
R_L 接在 C1 与 C2 的中点上。

求：f0、Rp、QL、通频带 2Δf0.7
"""
import cmath, math, os

OUT = r"E:\tmp\通信电路"

# ---------------- 已知 ----------------
L  = 1e-6        # H
Q0 = 100.0       # 电感空载品质因数
C1 = 25e-12      # F
C2 = 15e-12      # F
Ci = 5e-12       # F
Ri = 10e3        # Ω
RL = 5e3         # Ω

# 等效并联电容（C1 串 (C2∥RL) 在 R_L 远大于容抗时，C2∥RL 的容抗≈X_C2）
Cser = C1 * C2 / (C1 + C2)
Ceq  = Ci + Cser

# ---------------- 解析解 ----------------
w0_a = 1.0 / math.sqrt(L * Ceq)
f0_a = w0_a / (2 * math.pi)
Rp_a = 1.0 / (1.0 / Ri + 1.0 / (Q0 * w0_a * L) + 1.0 / RL_eff()) if False else None
def RL_eff():
    return RL
QL_a = w0_a * Ceq * Rp_a if Rp_a else None

# ---------------- 精确导纳模型（数值求解） ----------------
def Y(w):
    r = w * L / Q0                      # 电感串联等效损耗
    YL = 1.0 / (r + 1j * w * L)
    Z2 = 1.0 / (1.0 / RL + 1j * w * C2)  # C2 ∥ R_L
    ZD = 1.0 / (1j * w * C1) + Z2        # C1 串联 (C2∥RL)
    return 1.0 / Ri + 1j * w * Ci + YL + 1.0 / ZD

def Z(w):
    return 1.0 / Y(w)

# 求 Im(Y)=0
lo, hi = 1e6, 5e8
def imY(w):
    return Y(w).imag
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if imY(lo) * imY(mid) <= 0:
        hi = mid
    else:
        lo = mid
w0 = 0.5 * (lo + hi)
f0 = w0 / (2 * math.pi)
Y0 = Y(w0)
Rp = 1.0 / Y0.real
QL = w0 * Ceq * Rp
BW = f0 / QL

# ---------------- 频域扫描 ----------------
N = 40001
f1, f2 = 32e6, 52e6
fs = [f1 + (f2 - f1) * i / (N - 1) for i in range(N)]
ws = [2 * math.pi * f for f in fs]
zs = [abs(Z(w)) for w in ws]
imys = [Y(w).imag for w in ws]
reys = [Y(w).real for w in ws]
angs = [math.degrees(cmath.phase(Z(w))) for w in ws]

imax = max(range(N), key=lambda i: zs[i])
def cross(target, side):
    rng = range(imax, 0, -1) if side == "low" else range(imax, N - 1)
    for i in rng:
        j = i - 1 if side == "low" else i + 1
        if (zs[i] - target) * (zs[j] - target) <= 0 and zs[j] != zs[i]:
            t = (target - zs[i]) / (zs[j] - zs[i])
            return fs[i] + t * (fs[j] - fs[i])
    return None
Zmax = zs[imax]
t3 = Zmax / math.sqrt(2)
fL, fH = cross(t3, "low"), cross(t3, "high")
QL_bw = fs[imax] / (fH - fL)

# ---------------- 绘图 ----------------
from PIL import Image, ImageDraw, ImageFont

S = 2
W, H = 1400, 1265
img = Image.new("RGB", (W * S, H * S), "white")
dr = ImageDraw.Draw(img)
FT = ImageFont.truetype(r"C:\Windows\Fonts\msyhbd.ttc", 27 * S)
FL = ImageFont.truetype(r"C:\Windows\Fonts\msyh.ttc", 20 * S)
FM = ImageFont.truetype(r"C:\Windows\Fonts\msyh.ttc", 17 * S)

PL, PR = 125, 50
p1t, p1b = 200, 470
p2t, p2b = 610, 880
p3t, p3b = 1020, 1195

def fx(f):
    return (PL + (f - f1) / (f2 - f1) * (W - PR - PL)) * S

def fy(v, lo, hi, top, bot):
    return (bot - (v - lo) / (hi - lo) * (bot - top)) * S

def fylog(v, lo, hi, top, bot):
    return fy(math.log10(v), math.log10(lo), math.log10(hi), top, bot)

def xgrid(top, bot):
    for f in range(32, 53, 2):
        x = fx(f * 1e6)
        dr.line([x, top * S, x, bot * S], fill=(236, 236, 236), width=S)
        dr.text((x - 15 * S, (bot + 8) * S), str(f), font=FM, fill=(70, 70, 70))

BLU, RED, GRN = (20, 90, 200), (210, 40, 40), (0, 150, 60)

# ---- 面板 1：导纳电纳过零 ----
# y 轴范围按实际数据自适应，避免曲线被压成一条直线
_vmin, _vmax = min(imys), max(imys)
_span = max(_vmax - _vmin, 1e-4)
lo1, hi1 = _vmin - _span * 0.45, _vmax + _span * 0.45
_step = _span / 2.5
dr.rectangle([PL * S, p1t * S, (W - PR) * S, p1b * S], outline="black", width=2 * S)
_k = math.floor(lo1 / _step)
while _k * _step <= hi1:
    v = _k * _step
    y = fy(v, lo1, hi1, p1t, p1b)
    dr.line([PL * S, y, (W - PR) * S, y],
            fill=(150, 150, 150) if abs(v) < 1e-12 else (228, 228, 228), width=S)
    dr.text((34 * S, y - 11 * S), f"{v*1e3:+.1f}mS", font=FM, fill=(70, 70, 70))
    _k += 1
xgrid(p1t, p1b)
dr.text((PL * S, (p1t - 36) * S), "导纳虚部 Im Y（mS）—— 过零点即谐振频率 f₀",
        font=FL, fill="black")
dr.line([(fx(f), fy(v, lo1, hi1, p1t, p1b)) for f, v in zip(fs, imys)],
        fill=BLU, width=3 * S)
x0 = fx(f0)
dr.line([x0, p1t * S, x0, p1b * S], fill=RED, width=2 * S)
lab = f"f₀ = {f0/1e6:.4f} MHz"
dr.text((x0 + 8 * S, (p1t + 8) * S), lab, font=FM, fill=RED)
dr.text(((W - PR) * S / 2, (p1b + 38) * S), "频率 f / MHz", font=FL, fill="black")
dr.text((PL * S, (p1b + 72) * S),
        f"Y = 1/Ri + jωCi + 1/(ωL/Q₀ + jωL) + 1/[1/(jωC₁) + 1/(1/R_L + jωC₂)]"
        f"        Im Y(f₀) = {imys[imax]*1e3:.4f} mS ≈ 0",
        font=FM, fill=(110, 110, 110))

# ---- 面板 2：|Z| 谐振峰 ----
LO, HI = 1e2, 1e5
dr.rectangle([PL * S, p2t * S, (W - PR) * S, p2b * S], outline="black", width=2 * S)
for v, lab2 in [(100, "100"), (1000, "1k"), (10000, "10k"), (100000, "100k")]:
    y = fylog(v, LO, HI, p2t, p2b)
    dr.line([PL * S, y, (W - PR) * S, y], fill=(228, 228, 228), width=S)
    dr.text((30 * S, y - 11 * S), lab2, font=FM, fill=(70, 70, 70))
xgrid(p2t, p2b)
dr.text((PL * S, (p2t - 36) * S), "回路阻抗幅值 |Z|（Ω，对数轴）", font=FL, fill="black")
dr.line([(fx(f), fylog(z, LO, HI, p2t, p2b)) for f, z in zip(fs, zs)],
        fill=BLU, width=3 * S)
y3 = fylog(t3, LO, HI, p2t, p2b)
dr.line([PL * S, y3, (W - PR) * S, y3], fill=GRN, width=S)
dr.text(((W - PR) - 140 * S, y3 - 24 * S), "-3dB", font=FM, fill=GRN)
for f, col, lb, al in [(fL, GRN, f"fL = {fL/1e6:.3f}", "r"),
                       (fH, GRN, f"fH = {fH/1e6:.3f}", "l"),
                       (f0, RED, f"f₀ = {f0/1e6:.3f}", "c")]:
    x = fx(f)
    dr.line([x, p2t * S, x, p2b * S], fill=col, width=2 * S)
    w = dr.textlength(lb, font=FM)
    xx = {"r": x - w - 7 * S, "l": x + 7 * S, "c": x - w / 2}[al]
    dr.text((xx, (p2b - 120) * S), lb, font=FM, fill=col)
dr.text(((W - PR) * S / 2, (p2b + 38) * S), "频率 f / MHz", font=FL, fill="black")
dr.text((PL * S, (p2b + 72) * S),
        f"谐振点 Rp = {Rp:.0f} Ω        2Δf₀.₇ = {BW/1e6:.4f} MHz"
        f"        f₀/2Δf = {QL_bw:.2f}",
        font=FM, fill=(40, 40, 40))

# ---- 面板 3：相位 ----
dr.rectangle([PL * S, p3t * S, (W - PR) * S, p3b * S], outline="black", width=2 * S)
for a in (-90, 0, 90):
    y = fy(a, -90, 90, p3t, p3b)
    dr.line([PL * S, y, (W - PR) * S, y],
            fill=(150, 150, 150) if a == 0 else (228, 228, 228), width=S)
    dr.text((34 * S, y - 11 * S), f"{a}°", font=FM, fill=(70, 70, 70))
xgrid(p3t, p3b)
dr.text((PL * S, (p3t - 36) * S), "阻抗相位 ∠Z（°）", font=FL, fill="black")
dr.line([(fx(f), fy(a, -90, 90, p3t, p3b)) for f, a in zip(fs, angs)],
        fill=BLU, width=3 * S)
dr.line([x0, p3t * S, x0, p3b * S], fill=RED, width=2 * S)
dr.text(((W - PR) * S / 2, (p3b + 38) * S), "频率 f / MHz", font=FL, fill="black")

dr.text((PL * S, 45 * S), "习题 2.6  复杂并联谐振回路仿真", font=FT, fill="black")
dr.text((PL * S, 92 * S),
        f"L = 1μH   Q₀ = 100   C₁ = 25pF   C₂ = 15pF   Cᵢ = 5pF   "
        f"Rᵢ = 10kΩ   R_L = 5kΩ",
        font=FL, fill=(90, 90, 90))
dr.text((PL * S, 126 * S),
        f"→   f₀ = {f0/1e6:.3f} MHz   Rp = {Rp:.0f} Ω   "
        f"Q_L = {QL:.2f}   2Δf₀.₇ = {BW/1e6:.3f} MHz",
        font=FL, fill=(90, 90, 90))

img.resize((W, H), Image.LANCZOS).save(os.path.join(OUT, "ex2-6.png"), "PNG")

px = img.load()
for nm, (t, b) in [("ImY", (p1t, p1b)), ("|Z|", (p2t, p2b)), ("∠Z", (p3t, p3b))]:
    n = sum(1 for y in range(t * S, b * S, 7)
            for x in range(PL * S, (W - PR) * S, 7) if px[x, y][:3] == BLU)
    print(f"  自检 {nm:3s} 曲线像素 = {n:4d}  {'OK' if n > 50 else '未画出!'}")

print("=" * 64)
print("习题 2.6  复杂并联谐振回路")
print("=" * 64)
print(f"C1 串 C2 等效电容 Cser = {Cser*1e12:.3f} pF")
print(f"总并联电容 Ceq = Ci + Cser = {Ci*1e12:.1f} + {Cser*1e12:.3f} = {Ceq*1e12:.3f} pF")
print("-" * 64)
print(f"1) f0 = 1/(2π√(L·Ceq)) = {f0_a/1e6:.4f} MHz   (数值解 {f0/1e6:.4f} MHz)")
print("-" * 64)
print("2) 谐振点各支路电导：")
r0 = w0 * L / Q0
print(f"   R_i 支路      : 1/Ri          = {1/Ri*1e3:.4f} mS")
print(f"   L 支路        : 1/(Q0·ω0·L)  = {1/(Q0*w0*L)*1e3:.4f} mS   (r0={r0:.3f} Ω)")
YD = 1.0 / (1.0 / (1j * w0 * C1) + 1.0 / (1.0 / RL + 1j * w0 * C2))
print(f"   负载支路      : Re(1/Z_D)     = {YD.real*1e3:.4f} mS   "
      f"(等效 {1/YD.real:.0f} Ω)")
print(f"   合计 G        = {Y0.real*1e3:.4f} mS  →  Rp = {Rp:.1f} Ω = {Rp/1e3:.3f} kΩ")
print("-" * 64)
print(f"3) QL = ω0·Ceq·Rp = {w0:.5e}×{Ceq:.4e}×{Rp:.1f} = {QL:.3f}")
Qi_, Qri, Qd = Q0, w0 * Ceq * Ri, w0 * Ceq / YD.real
print(f"   校核 1/QL = 1/Q0 + 1/Q(Ri) + 1/Q(D)")
print(f"        Q0={Qi_:.1f}, Q(Ri)={Qri:.2f}, Q(D)={Qd:.2f}"
      f"  →  QL = {1/(1/Qi_+1/Qri+1/Qd):.3f}")
print(f"4) 2Δf0.7 = f0/QL = {BW/1e6:.4f} MHz")
print("-" * 64)
print("仿真验证：")
print(f"  数值 f0 = {f0/1e6:.6f} MHz,  Rp = {1/Y(fs[imax]*2*math.pi).real:.1f} Ω")
print(f"  -3dB 频率 = {fL/1e6:.4f} / {fH/1e6:.4f} MHz,  f0/2Δf = {QL_bw:.3f}")
print(f"  与解析 QL 偏差 = {abs(QL_bw-QL)/QL*100:.3f} %")