# -*- coding: utf-8 -*-
"""教材习题 2.5  变压器耦合谐振回路

题2.5图接法（已逐点核对原图）：
  节点T(上轨) : R, C, 端子1, 端子4, C_L, R_L
  节点M(中轨) : Ig, R_g, C_g, 端子2(抽头)
  节点B(下轨) : 端子3, 端子5 及各支路下端
  R、C 两支路穿过中轨不相连（交叉无结点）

绕组：1-2 (N12=6) 与 2-3 (N23=4) 同属一线圈 L13；4-5 (N45=3) 为同芯次级。
求：电感 L13 与有载品质因数 QL
"""
import cmath, math, os

OUT = r"E:\tmp\通信电路"

# ---------------- 已知 ----------------
f  = 20e6          # 工作频率 Hz
C  = 25e-12        # 回路电容
Q0 = 60.0          # 线圈空载品质因数
N12, N23, N45 = 6, 4, 3
R  = 10e3          # 回路并联电阻
Rg = 2.5e3         # 信号源内阻
RL = 900.0         # 负载电阻
Cg = 9e-12         # 源输出电容
CL = 12e-12        # 负载电容

N13 = N12 + N23    # 一次线圈总匝数 = 10

# ---------------- 1) 把各支路折算到 1-3 全线圈侧 ----------------
# 源侧：跨 2-3 绕组接入（自耦抽头），用式(2-11)
k_g = N13 / N23            # 匝比 = 10/4 = 2.5
Rg_p = k_g ** 2 * Rg       # 式(2-11)
Cg_p = Cg / k_g ** 2       # 电容按阻抗平方反比折算
# 负载侧：4-5 次级折算到 1-3 一次侧，用式(2-10)
k_L = N13 / N45            # 匝比 = 10/3
RL_p = k_L ** 2 * RL       # 式(2-10)
CL_p = CL / k_L ** 2

C_tot = C + Cg_p + CL_p

# ---------------- 2) 谐振条件求 L13 ----------------
w0 = 2 * math.pi * f
L13 = 1.0 / (w0 * w0 * C_tot)

# ---------------- 3) 有载品质因数 ----------------
# 线圈并联等效损耗电阻
X_L = w0 * L13
Rp = Q0 * X_L                       # 空载谐振电阻
# 外部总负载电导
G_e = 1.0 / R + 1.0 / Rg_p + 1.0 / RL_p
R_e = 1.0 / G_e
R_tot = 1.0 / (1.0 / Rp + G_e)
QL = R_tot * w0 * C_tot

# 损耗相加校核
Q_R, Q_gp, Q_Lp = w0 * C_tot * R, w0 * C_tot * Rg_p, w0 * C_tot * RL_p
QL_chk = 1.0 / (1.0 / Q0 + 1.0 / Q_R + 1.0 / Q_gp + 1.0 / Q_Lp)

# ---------------- 等效回路频域仿真 ----------------
def Y(w):
    return (1.0 / Rp + G_e) + 1j * (w * C_tot - 1.0 / (w * L13))

def Z(w):
    return 1.0 / Y(w)

N = 40001
f1, f2 = 17e6, 23e6
fs = [f1 + (f2 - f1) * i / (N - 1) for i in range(N)]
zs = [abs(Z(2 * math.pi * x)) for x in fs]
angs = [math.degrees(cmath.phase(Z(2 * math.pi * x))) for x in fs]

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
W, H = 1360, 940
img = Image.new("RGB", (W * S, H * S), "white")
dr = ImageDraw.Draw(img)
FT = ImageFont.truetype(r"C:\Windows\Fonts\msyhbd.ttc", 27 * S)
FL = ImageFont.truetype(r"C:\Windows\Fonts\msyh.ttc", 20 * S)
FM = ImageFont.truetype(r"C:\Windows\Fonts\msyh.ttc", 17 * S)

PL, PR = 120, 50
p1t, p1b = 200, 520
p2t, p2b = 660, 880

def fx(x):
    return (PL + (x - f1) / (f2 - f1) * (W - PR - PL)) * S

def fy(v, lo, hi, top, bot):
    return (bot - (v - lo) / (hi - lo) * (bot - top)) * S

def fylog(v, lo, hi, top, bot):
    return fy(math.log10(v), math.log10(lo), math.log10(hi), top, bot)

def xgrid(top, bot):
    for x in range(17, 24):
        px = fx(x * 1e6)
        dr.line([px, top * S, px, bot * S], fill=(236, 236, 236), width=S)
        dr.text((px - 26 * S, (bot + 8) * S), f"{x}", font=FM, fill=(70, 70, 70))

BLU, RED, GRN = (20, 90, 200), (210, 40, 40), (0, 150, 60)

# 面板1：|Z|
LO, HI = 200.0, 3e4
dr.rectangle([PL * S, p1t * S, (W - PR) * S, p1b * S], outline="black", width=2 * S)
for v, lb in [(300, "300"), (1000, "1k"), (3000, "3k"), (10000, "10k")]:
    y = fylog(v, LO, HI, p1t, p1b)
    dr.line([PL * S, y, (W - PR) * S, y], fill=(228, 228, 228), width=S)
    dr.text((40 * S, y - 11 * S), lb, font=FM, fill=(70, 70, 70))
xgrid(p1t, p1b)
dr.text((PL * S, (p1t - 36) * S), "等效并联谐振回路阻抗幅值 |Z|（Ω，对数轴）",
        font=FL, fill="black")
dr.line([(fx(x), fylog(z, LO, HI, p1t, p1b)) for x, z in zip(fs, zs)],
        fill=BLU, width=3 * S)
y3 = fylog(t3, LO, HI, p1t, p1b)
dr.line([PL * S, y3, (W - PR) * S, y3], fill=GRN, width=S)
dr.text(((W - PR) - 120 * S, y3 - 24 * S), "-3dB", font=FM, fill=GRN)
for x, col, lb, al in [(fL, GRN, f"fL = {fL/1e6:.3f}", "r"),
                       (fH, GRN, f"fH = {fH/1e6:.3f}", "l"),
                       (f, RED, f"f₀ = 20.000", "c")]:
    px = fx(x)
    dr.line([px, p1t * S, px, p1b * S], fill=col, width=2 * S)
    w = dr.textlength(lb, font=FM)
    xx = {"r": px - w - 7 * S, "l": px + 7 * S, "c": px - w / 2}[al]
    dr.text((xx, (p1b - 130) * S), lb, font=FM, fill=col)
dr.text(((W - PR) * S / 2, (p1b + 38) * S), "频率 f / MHz", font=FL, fill="black")
dr.text((PL * S, (p1b + 74) * S),
        f"谐振点 R_p(有载) = {R_tot:.0f} Ω        2Δf = {(fH-fL)/1e6:.4f} MHz"
        f"        f₀/2Δf = {QL_bw:.2f}",
        font=FM, fill=(40, 40, 40))

# 面板2：相位
dr.rectangle([PL * S, p2t * S, (W - PR) * S, p2b * S], outline="black", width=2 * S)
for a in (-90, 0, 90):
    y = fy(a, -90, 90, p2t, p2b)
    dr.line([PL * S, y, (W - PR) * S, y],
            fill=(150, 150, 150) if a == 0 else (228, 228, 228), width=S)
    dr.text((34 * S, y - 11 * S), f"{a}°", font=FM, fill=(70, 70, 70))
xgrid(p2t, p2b)
dr.text((PL * S, (p2t - 36) * S), "阻抗相位 ∠Z（°）", font=FL, fill="black")
dr.line([(fx(x), fy(a, -90, 90, p2t, p2b)) for x, a in zip(fs, angs)],
        fill=BLU, width=3 * S)
dr.line([fx(f), p2t * S, fx(f), p2b * S], fill=RED, width=2 * S)
dr.text(((W - PR) * S / 2, (p2b + 38) * S), "频率 f / MHz", font=FL, fill="black")

dr.text((PL * S, 45 * S), "习题 2.5  变压器耦合谐振回路仿真", font=FT, fill="black")
dr.text((PL * S, 92 * S),
        f"f = 20MHz   C = 25pF   Q0 = 60   N12 = 6, N23 = 4, N45 = 3   "
        f"R = 10kΩ   Rg = 2.5kΩ   RL = 900Ω   Cg = 9pF   CL = 12pF",
        font=FL, fill=(90, 90, 90))
dr.text((PL * S, 126 * S),
        f"→   L13 = {L13*1e6:.3f} μH      Q_L = {QL:.2f}      "
        f"通频带 2Δf = {f/QL/1e6:.3f} MHz",
        font=FL, fill=(90, 90, 90))

img.resize((W, H), Image.LANCZOS).save(os.path.join(OUT, "ex2-5.png"), "PNG")

px_ = img.load()
for nm, (t, b) in [("|Z|", (p1t, p1b)), ("∠Z", (p2t, p2b))]:
    n = sum(1 for y in range(t * S, b * S, 7)
            for x in range(PL * S, (W - PR) * S, 7) if px_[x, y][:3] == BLU)
    print(f"  自检 {nm:3s} 曲线像素 = {n:4d}  {'OK' if n > 50 else '未画出!'}")

print("=" * 66)
print("习题 2.5  变压器耦合谐振回路")
print("=" * 66)
print(f"一次线圈总匝数 N13 = N12 + N23 = {N12} + {N23} = {N13}")
print("-" * 66)
print("第一步  各支路折算到 1-3 全线圈侧")
print(f"  源侧(跨 2-3 自耦接入, 式2-11):  匝比 k = N13/N23 = {k_g:.3f}")
print(f"    Rg' = k²Rg   = {Rg_p:.0f} Ω")
print(f"    Cg' = Cg/k²  = {Cg_p*1e12:.3f} pF")
print(f"  负载侧(4-5 变压器耦合, 式2-10): 匝比 k = N13/N45 = {k_L:.3f}")
print(f"    RL' = k²RL   = {RL_p:.0f} Ω")
print(f"    CL' = CL/k²  = {CL_p*1e12:.3f} pF")
print("-" * 66)
print("第二步  谐振条件求 L13")
print(f"  C总 = C + Cg' + CL' = {C*1e12:.0f} + {Cg_p*1e12:.3f} + {CL_p*1e12:.3f} = {C_tot*1e12:.3f} pF")
print(f"  L13 = 1/(w0²·C总) = 1/({w0:.4e}² × {C_tot:.4e}) = {L13*1e6:.4f} μH")
print("-" * 66)
print("第三步  有载品质因数 QL")
print(f"  w0·L13 = {X_L:.2f} Ω")
print(f"  线圈空载损耗 Rp = Q0·w0L = {Rp:.0f} Ω")
print(f"  外部等效电阻 Re = R‖Rg'‖RL' = {R/1e3:g}k‖{Rg_p/1e3:.3f}k‖{RL_p/1e3:g}k"
      f" = {R_e:.0f} Ω")
print(f"  总谐振电阻 = Rp‖Re = {R_tot:.0f} Ω")
print(f"  QL = R总·w0·C总 = {R_tot:.1f} × {w0*C_tot:.5e} = {QL:.3f}")
print(f"  校核 1/QL = 1/Q0+1/QR+1/Qg'+1/QL'")
print(f"       Q0={Q0:.1f} QR={Q_R:.2f} Qg'={Q_gp:.2f} QL'={Q_Lp:.2f}"
      f"  →  QL = {QL_chk:.3f}")
print(f"  通频带 2Δf = f0/QL = {f/QL/1e6:.4f} MHz")
print("-" * 66)
print("仿真验证：")
print(f"  数值 f0 = {fs[imax]/1e6:.6f} MHz")
print(f"  数值 |Z|max = {Zmax:.1f} Ω")
print(f"  -3dB 频率 = {fL/1e6:.4f} / {fH/1e6:.4f} MHz")
print(f"  f0/(2Δf) = {QL_bw:.3f}   与解析 QL 偏差 = {abs(QL_bw-QL)/QL*100:.3f} %")