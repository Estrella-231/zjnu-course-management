# -*- coding: utf-8 -*-
"""教材习题 2.3 数值验证 + 仿真

已知 LCR 并联谐振回路，f0=12MHz，电感 L 在 12MHz 时测得 L=3uH, Q0=100，
并联电阻 R=10kohm。求 C、Rp、QL。
"""
import cmath, math, os

OUT = r"E:\tmp\通信电路"

# ---------------- 已知条件 ----------------
f0  = 12e6          # 谐振频率 Hz
L   = 3e-6          # 电感 H
Q0  = 100.0         # 电感空载品质因数
R   = 10e3          # 并联（负载）电阻 ohm

w0  = 2 * math.pi * f0
rho = w0 * L                # 特性阻抗 = w0*L

# ---------------- 1. 谐振电容 C ----------------
C = 1.0 / (w0 * w0 * L)

# ---------------- 2. 电感的等效并联损耗电阻 ----------------
# 串联等效：Q0 = w0*L/r0  ->  r0 = w0*L/Q0
r0 = w0 * L / Q0
# 变换到并联：Rp0 = Q0^2 * r0 = Q0 * w0*L
Rp0 = Q0 * w0 * L

# ---------------- 3. 回路谐振电阻 Rp ----------------
Rp = 1.0 / (1.0 / Rp0 + 1.0 / R)

# ---------------- 4. 有载品质因数 QL ----------------
# 并联回路：QL = Rp * w0*C = 2*pi*f0*C / G0
QL = Rp * w0 * C
# 独立校核：1/QL = 1/Q0 + 1/QR，其中 QR = R*w0*C
QR = R * w0 * C
QL_chk = 1.0 / (1.0 / Q0 + 1.0 / QR)

BW = f0 / QL

# ---------------- 频域仿真（并联回路导纳） ----------------
def Y(f):
    w = 2 * math.pi * f
    return (1.0 / Rp0 + 1.0 / R) + 1j * (w * C - 1.0 / (w * L))

def Z(f):
    return 1.0 / Y(f)

N = 20001
f1, f2 = 9e6, 15e6
fs = [f1 + (f2 - f1) * i / (N - 1) for i in range(N)]
zs = [abs(Z(f)) for f in fs]
angs = [math.degrees(cmath.phase(Z(f))) for f in fs]

# 数值找谐振点（|Z| 最大）
imax = max(range(N), key=lambda i: zs[i])
f_num = fs[imax]
Zmax = zs[imax]

def cross(target, side):
    """在谐振点某一侧线性插值找 |Z|=target 的频率。side: 'low' / 'high'。"""
    rng = range(imax, 0, -1) if side == "low" else range(imax, N - 1)
    for i in rng:
        j = i - 1 if side == "low" else i + 1
        if (zs[i] - target) * (zs[j] - target) <= 0 and zs[j] != zs[i]:
            t = (target - zs[i]) / (zs[j] - zs[i])
            return fs[i] + t * (fs[j] - fs[i])
    return None

target = Zmax / math.sqrt(2)
fL, fH = cross(target, "low"), cross(target, "high")
QL_bw = f_num / (fH - fL) if (fL and fH) else float("nan")

# ---------------- 绘图 ----------------
from PIL import Image, ImageDraw, ImageFont

S = 2
W, H = 1360, 1000
img = Image.new("RGB", (W * S, H * S), "white")
dr = ImageDraw.Draw(img)
FT = ImageFont.truetype(r"C:\Windows\Fonts\msyhbd.ttc", 27 * S)
FL = ImageFont.truetype(r"C:\Windows\Fonts\msyh.ttc", 20 * S)
FM = ImageFont.truetype(r"C:\Windows\Fonts\msyh.ttc", 17 * S)

PL, PR = 120, 50
p1t, p1b = 195, 495
p2t, p2b = 645, 870

def fx(f):
    return (PL + (f - f1) / (f2 - f1) * (W - PR - PL)) * S

def fy(v, lo, hi, top, bot):
    """线性映射，返回像素坐标(已乘 S)。"""
    return (bot - (v - lo) / (hi - lo) * (bot - top)) * S

def fylog(v, lo, hi, top, bot):
    return fy(math.log10(v), math.log10(lo), math.log10(hi), top, bot)

def xgrid(top, bot):
    for f in range(9, 16):
        x = fx(f * 1e6)
        dr.line([x, top * S, x, bot * S], fill=(236, 236, 236), width=S)
        dr.text((x - 14 * S, (bot + 9) * S), str(f), font=FM, fill=(70, 70, 70))

# ---------- 面板 1：|Z| 对数轴 ----------
LO, HI = 10.0, 3e4
dr.rectangle([PL * S, p1t * S, (W - PR) * S, p1b * S], outline="black", width=2 * S)
for v, lab in [(10, "10"), (100, "100"), (1000, "1k"), (10000, "10k")]:
    y = fylog(v, LO, HI, p1t, p1b)
    dr.line([PL * S, y, (W - PR) * S, y], fill=(228, 228, 228), width=S)
    dr.text((38 * S, y - 11 * S), lab, font=FM, fill=(70, 70, 70))
xgrid(p1t, p1b)
dr.text((PL * S, (p1t - 36) * S), "并联谐振回路阻抗幅值 |Z|（Ω，对数轴）",
        font=FL, fill="black")

curve = [(fx(f), fylog(z, LO, HI, p1t, p1b)) for f, z in zip(fs, zs)]
dr.line(curve, fill=(20, 90, 200), width=3 * S)

y3 = fylog(Zmax / math.sqrt(2), LO, HI, p1t, p1b)
dr.line([PL * S, y3, (W - PR) * S, y3], fill=(0, 160, 60), width=S)
dr.text(((W - PR) - 150 * S, y3 - 22 * S), "-3dB 电平", font=FM, fill=(0, 130, 50))

for f, col, lab, dy, align in [
        (fL, (0, 150, 60), f"fL = {fL/1e6:.4f}", 6, "r"),
        (fH, (0, 150, 60), f"fH = {fH/1e6:.4f}", 6, "l"),
        (f_num, (210, 40, 40), f"f0 = {f_num/1e6:.4f} MHz", 34, "c")]:
    if f:
        x = fx(f)
        dr.line([x, p1t * S, x, p1b * S], fill=col, width=2 * S)
        w = dr.textlength(lab, font=FM)
        x0 = {"r": x - w - 7 * S, "l": x + 7 * S, "c": x - w / 2}[align]
        dr.text((x0, (p1t + dy) * S), lab, font=FM, fill=col)

dr.text(((W - PR) * S / 2, (p1b + 40) * S), "频率 f / MHz", font=FL, fill="black")
dr.text((PL * S, (p1b + 74) * S),
        f"谐振点 |Z|max = Rp = {Rp:.0f} Ω        通频带 2Δf = {(fH - fL)/1e6:.4f} MHz"
        f"        f0 / 2Δf = {QL_bw:.2f}",
        font=FM, fill=(40, 40, 40))

# ---------- 面板 2：相位 ----------
dr.rectangle([PL * S, p2t * S, (W - PR) * S, p2b * S], outline="black", width=2 * S)
for a in (-90, -45, 0, 45, 90):
    y = fy(a, -90, 90, p2t, p2b)
    dr.line([PL * S, y, (W - PR) * S, y],
            fill=(150, 150, 150) if a == 0 else (228, 228, 228), width=S)
    dr.text((34 * S, y - 11 * S), f"{a}°", font=FM, fill=(70, 70, 70))
xgrid(p2t, p2b)
dr.text((PL * S, (p2t - 36) * S), "阻抗相位 ∠Z（°）", font=FL, fill="black")
dr.line([(fx(f), fy(a, -90, 90, p2t, p2b)) for f, a in zip(fs, angs)],
        fill=(20, 90, 200), width=3 * S)
x0 = fx(f_num)
dr.line([x0, p2t * S, x0, p2b * S], fill=(210, 40, 40), width=2 * S)
dr.text(((W - PR) * S / 2, (p2b + 40) * S), "频率 f / MHz", font=FL, fill="black")

dr.text((PL * S, 45 * S), "习题 2.3  LCR 并联谐振回路仿真验证", font=FT, fill="black")
dr.text((PL * S, 92 * S),
        f"L = 3μH    Q0 = 100    R = 10kΩ    f0 = 12MHz     →     "
        f"C = {C*1e12:.2f} pF     Rp = {Rp:.0f} Ω     QL = {QL:.2f}",
        font=FL, fill=(90, 90, 90))
dr.text((PL * S, 126 * S),
        f"导纳 Y(f) = 1/Rp0 + 1/R + j(ωC − 1/ωL)，Rp0 = Q0·ω0·L = {Rp0/1e3:.2f} kΩ"
        f"（电感损耗）；解析 QL = Rp·ω0·C 与仿真 f0/2Δf 偏差 {abs(QL_bw-QL)/QL*100:.3f}%",
        font=FM, fill=(120, 120, 120))

img.resize((W, H), Image.LANCZOS).save(os.path.join(OUT, "ex2-3.png"), "PNG")

# 自检：确认两条曲线确实画上了
px = img.load()
for nm, (t, b) in [("|Z|", (p1t, p1b)), ("∠Z", (p2t, p2b))]:
    n = sum(1 for y in range(t * S, b * S, 7)
            for x in range(PL * S, (W - PR) * S, 7)
            if px[x, y][:3] == (20, 90, 200))
    print(f"  自检 面板{nm} 曲线像素 = {n}  {'OK' if n > 50 else '未画出!'}")

# ---------------- 输出 ----------------
print("=" * 62)
print("习题 2.3  LCR 并联谐振回路")
print("=" * 62)
print(f"f0  = {f0/1e6:.2f} MHz      w0 = {w0:.4e} rad/s")
print(f"L   = {L*1e6:.2f} uH        w0*L = {rho:.2f} ohm")
print(f"Q0  = {Q0:.0f}   R = {R/1e3:.1f} kohm")
print("-" * 62)
print(f"1) C  = 1/(w0^2 L)   = {C*1e12:.3f} pF")
print(f"2) r0 = w0 L / Q0    = {r0:.3f} ohm  (电感等效串联损耗)")
print(f"   Rp0 = Q0*w0*L    = {Rp0/1e3:.3f} kohm (电感等效并联损耗)")
print(f"   Rp = Rp0 || R    = {Rp0/1e3:.2f}k || {R/1e3:.1f}k = {Rp:.1f} ohm = {Rp/1e3:.3f} kohm")
print(f"3) QL = Rp*w0*C      = {QL:.3f}")
print(f"   校核 QR = R*w0*C  = {QR:.2f},  1/QL = 1/Q0 + 1/QR -> QL = {QL_chk:.3f}")
print(f"   2df = f0/QL       = {BW/1e6:.4f} MHz")
print("-" * 62)
print("仿真:")
print(f"  数值谐振频率 = {f_num/1e6:.6f} MHz")
print(f"  数值 |Z|max  = {Zmax:.1f} ohm")
print(f"  -3dB 频率    = {fL/1e6:.4f} / {fH/1e6:.4f} MHz")
print(f"  f0/(2df)     = {QL_bw:.3f}   (与解析值偏差 {abs(QL_bw-QL)/QL*100:.3f} %)")