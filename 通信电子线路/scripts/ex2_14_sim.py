# -*- coding: utf-8 -*-
"""习题 2.14  有高频放大器的接收机方框图 —— 总噪声系数与变频噪声要求对比
教材口径: 2.4.3 节 式(2-80) 多级级联噪声系数, 式(2-81) 无源网络噪声系数 N_F = 1/A_PH
方框图(题 2.14 图)级联次序:
   天线 -> 传输线(A_PH1=0.8) -> 高频小信号放大器(N_F2=2, A_PH2=60)
        -> 变频器(N_F3=6, A_PH3=0.3) -> 中频放大器(N_F4=4, A_PH4=1e4)
仿真: 逐级等效输入噪声功率递推 + 复高斯蒙特卡洛(直接叠加各级噪声, 不调用闭式公式)
"""
import math
import random
from PIL import Image, ImageDraw, ImageFont

# ---------------- 已知条件 ----------------
G1 = 0.8                       # 传输线功率增益 (无源)
NF_HF, G2 = 2.0, 60.0          # 高频小信号放大器
NF_MX, G3 = 6.0, 0.3           # 变频器
NF_IF, G4 = 4.0, 1.0e4          # 中频放大器
NF_LINE = 1.0 / G1              # 式(2-81): 无源网络噪声系数 = 1/A_PH

T0 = 290.0                     # 参考噪声温度, K

# ---------------- 级联(有高放) ----------------
STAGE_HF = [(G1, NF_LINE), (G2, NF_HF), (G3, NF_MX), (G4, NF_IF)]
# 级联(无高放): 变频器直接接在传输线之后
STAGE_NO = [(G1, NF_LINE), (G3, NF_MX), (G4, NF_IF)]


def friis(stages):
    """式(2-80): N_F = F1 + (F2-1)/G1 + (F3-1)/(G1G2) + ...
    处理 stages[k] 时分母 pre = G1*G2*...*G_{k-1}, 用完再乘本级增益"""
    f = stages[0][1]
    pre = stages[0][0]
    for k in range(1, len(stages)):
        f += (stages[k][1] - 1.0) / pre
        pre *= stages[k][0]
    return f


def cascade_power(stages):
    """逐级等效输入噪声功率递推。
    输入端先注入参考噪声 kT0B(归一化为 1), Pn <- (Pn + (F-1)) * G , Ps <- Ps * G"""
    pn, ps, rec = 1.0, 1.0, []
    for i, (g, F) in enumerate(stages):
        pn = (pn + (F - 1.0)) * g
        ps *= g
        rec.append((i + 1, g, F, pn, ps, pn / ps))
    return pn / ps, rec


def contrib(stages):
    """各级对总噪声系数的贡献量 (式(2-80) 各项)"""
    out, pre = [stages[0][1]], stages[0][0]
    for (g, F) in stages[1:]:
        out.append((F - 1.0) / pre)
        pre *= g
    return out


def cascade_mc(stages, nsamp=200000, seed=20251014):
    """蒙特卡洛: 复高斯白噪声逐级叠加, 输出噪声功率 / 总功率增益 = 噪声系数。
    不调用任何闭式噪声系数公式。
    输入端注入功率归一化为 1 的复高斯噪声(即参考噪声 kT0B);
    每级: n <- sqrt(G) * (n + sqrt(F-1) * w), 于是 E|n|^2 递推 = G*(Pn + F - 1)。"""
    rnd = random.Random(seed)
    sq2 = math.sqrt(2.0)
    n = [complex(rnd.gauss(0.0, 1.0), rnd.gauss(0.0, 1.0)) / sq2
         for _ in range(nsamp)]                       # 输入参考噪声 kT0B
    s = 1.0
    for (g, F) in stages:
        ag = math.sqrt(g)
        an = math.sqrt(F - 1.0)                       # 本级等效输入噪声 (相对 kT0B)
        for i in range(nsamp):
            w = complex(rnd.gauss(0.0, 1.0), rnd.gauss(0.0, 1.0)) / sq2
            n[i] = ag * (n[i] + an * w)
        s *= g
    p = 0.0
    for x in n:
        p += abs(x) ** 2
    p /= nsamp
    return p / s, s


def db(v):
    return 10.0 * math.log10(v)


# ---------------- 解析主结果 ----------------
NF_hf, rec_hf = cascade_power(STAGE_HF)
NF_no, rec_no = cascade_power(STAGE_NO)
NF_hf_friis = friis(STAGE_HF)
NF_no_friis = friis(STAGE_NO)

c_hf = contrib(STAGE_HF)          # [1.25, 1.25, 0.104167, 0.208333]
c_no = contrib(STAGE_NO)          # [1.25, 6.25, 12.5]

S_hf = 1.0 / (G1 * G2)            # 变频器噪声系数灵敏度: 1/48
S_no = 1.0 / G1                   # 1.25
Te_hf, Te_no = (NF_hf - 1) * T0, (NF_no - 1) * T0
BASE_HF = NF_hf - c_hf[2]         # 变频器理想时 N_F 的下限 2.708333
BASE_NO = NF_no - c_no[1]         # 变频器理想时 N_F 的下限 13.75


def nf_hf(f3):
    return BASE_HF + (f3 - 1.0) / (G1 * G2)


def nf_no(f3):
    return BASE_NO + (f3 - 1.0) / G1


print("=" * 78)
print("【有高放】级联: 传输线 -> 高放 -> 变频器 -> 中放")
print("  各级贡献 (式2-80): ", ["%.6f" % v for v in c_hf], " 合计 %.6f" % sum(c_hf))
print("  逐级递推 N_F = %.6f   式(2-80) N_F = %.6f   互差 %.3e" % (NF_hf, NF_hf_friis, abs(NF_hf-NF_hf_friis)))
print("  N_F = %.6f = %.4f dB     等效噪声温度 T_e = %.1f K" % (NF_hf, db(NF_hf), Te_hf))
print("  变频器贡献 %.6f (%.2f%%)   灵敏度 dNF/dNF3 = 1/(G1*G2) = %.6f" % (c_hf[2], 100*c_hf[2]/NF_hf, S_hf))
print("-" * 78)
print("【无高放】级联: 传输线 -> 变频器 -> 中放")
print("  各级贡献 (式2-80): ", ["%.6f" % v for v in c_no], " 合计 %.6f" % sum(c_no))
print("  逐级递推 N_F = %.6f   式(2-80) N_F = %.6f   互差 %.3e" % (NF_no, NF_no_friis, abs(NF_no-NF_no_friis)))
print("  N_F = %.6f = %.4f dB     等效噪声温度 T_e = %.1f K" % (NF_no, db(NF_no), Te_no))
print("  变频器贡献 %.6f (%.2f%%)   灵敏度 dNF/dNF3 = 1/G1 = %.6f" % (c_no[1], 100*c_no[1]/NF_no, S_no))
print("-" * 78)
print("变频器噪声贡献之比  无高放/有高放 = %.2f  (= A_PH2 = %.0f)" % (c_no[1]/c_hf[2], G2))
print("灵敏度系数之比      无高放/有高放 = %.2f  (= A_PH2)" % (S_no/S_hf))
print("-" * 78)
print("允许 N_F3 上限 (由总噪声预算反解):")
for b in (3.0, 5.0, 10.0, 20.0, 40.0):
    a = 1 + (G1*G2)*(b - BASE_HF)
    n_ = 1 + (b - BASE_NO)/(1.0/G1)
    print("    预算 N_F = %5.1f : 有高放 N_F3 <= %10.2f   无高放 N_F3 <= %s"
          % (b, a, ("%8.2f" % n_) if n_ >= 1.0 else "  不可行"))
print("=" * 78)

# ---------------- 蒙特卡洛仿真 ----------------
MC = []
print("蒙特卡洛仿真 (复高斯逐级噪声叠加, nsamp = 200000):")
for tag, base in (("有高放", STAGE_HF), ("无高放", STAGE_NO)):
    for f3 in (3.0, 6.0, 12.0):
        st = [(g, (f3 if F == NF_MX else F)) for (g, F) in base]
        f_mc, gt = cascade_mc(st, seed=20251014 + int(f3*7) + len(MC))
        f_an = friis(st)
        MC.append((tag, f3, f_an, f_mc, gt))
        print("  %s N_F3 = %4.1f : 解析 %10.6f   蒙特卡洛 %10.6f   偏差 %+.3f%%"
              % (tag, f3, f_an, f_mc, 100*(f_mc-f_an)/f_an))
print("=" * 78)

# ---------------- 绘图 ----------------
FN = r"C:\Windows\Fonts\msyh.ttc"
FB = r"C:\Windows\Fonts\msyhbd.ttc"
W, H = 1620, 1500
img = Image.new("RGB", (W, H), "white")
dr = ImageDraw.Draw(img)
f22 = ImageFont.truetype(FB, 22)
f18 = ImageFont.truetype(FB, 18)
f16 = ImageFont.truetype(FN, 16)
f15 = ImageFont.truetype(FN, 15)
f14 = ImageFont.truetype(FN, 14)
f13 = ImageFont.truetype(FN, 13)
BLUE, RED, GRAY = (20, 70, 190), (200, 30, 30), (110, 110, 110)
ORNG, GREEN = (225, 140, 0), (0, 145, 70)


def panel(x0, y0, x1, y1, title):
    dr.rectangle([x0, y0, x1, y1], outline=(60, 60, 60), width=2)
    dr.text((x0 + 10, y0 - 27), title, font=f18, fill=(0, 0, 0))


def axes(p, xlab, ylab, xlim, ylim, logy=False):
    x0, y0, x1, y1 = p
    if logy:
        lo, hi = math.log10(ylim[0]), math.log10(ylim[1])
        k = 0
        while 10.0**k <= ylim[1]*1.001:
            yv = 10.0**k
            py = y1 - (y1-y0)*(math.log10(yv)-lo)/(hi-lo)
            if y0-1 <= py <= y1+1:
                dr.line([x0, py, x1, py], fill=(228, 228, 228))
                dr.line([x0, py, x0, py-9], fill=(60, 60, 60))
                dr.text((x0-58, py-9), ("%g" % yv), font=f14, fill=GRAY)
            k += 1
    else:
        for i in range(5):
            yy = y0 + (y1-y0)*i/4.0
            dr.line([x0, yy, x1, yy], fill=(228, 228, 228))
            dr.line([x0, yy-9, x0, yy], fill=(60, 60, 60))
            dr.text((x0-78, yy-9), ("%.3g" % (ylim[1]-(ylim[1]-ylim[0])*i/4.0)), font=f14, fill=GRAY)
    for i in range(7):
        xx = x0 + (x1-x0)*i/6.0
        dr.line([xx, y0, xx, y1], fill=(238, 238, 238))
        dr.line([xx, y1, xx, y1+9], fill=(60, 60, 60))
        dr.text((xx-24, y1+13), ("%.4g" % (xlim[0]+(xlim[1]-xlim[0])*i/6.0)), font=f14, fill=GRAY)
    dr.line([x0, y1, x1, y1], fill=(60, 60, 60), width=2)
    dr.line([x0, y0, x0, y1], fill=(60, 60, 60), width=2)
    dr.text((x1-150, y1+34), xlab, font=f14, fill=(0, 0, 0))
    dr.text((x0+14, y0+8), ylab, font=f14, fill=(130, 130, 130))


def plot(p, xs, ys, xlim, ylim, color, width=3, logy=False):
    x0, y0, x1, y1 = p
    lo, hi = (math.log10(ylim[0]), math.log10(ylim[1])) if logy else (ylim[0], ylim[1])
    pts = []
    for xv, yv in zip(xs, ys):
        px = x0 + (x1-x0)*(xv-xlim[0])/(xlim[1]-xlim[0])
        t = (math.log10(yv)-lo)/(hi-lo) if logy else (yv-lo)/(hi-lo)
        py = y1 - (y1-y0)*t
        if not (x0-1 <= px <= x1+1 and y0-1 <= py <= y1+1):
            continue
        pts.append((px, py))
    if len(pts) > 1:
        dr.line(pts, fill=color, width=width)


def vline(p, xv, xlim, color, txt=None, dy=8):
    x0, y0, x1, y1 = p
    px = x0 + (x1-x0)*(xv-xlim[0])/(xlim[1]-xlim[0])
    if x0 <= px <= x1:
        dr.line([px, y0, px, y1], fill=color, width=1)
    if txt:
        dr.text((px + 6, y0 + dy), txt, font=f14, fill=color)


dr.text((24, 12), "习题 2.14  有高频放大器的接收机方框图 —— 总噪声系数与变频噪声要求对比",
        font=f22, fill=(0, 0, 0))

# ========== 面板A: 各级噪声贡献分解 ==========
panel(80, 100, 810, 540, "面板 A  总噪声系数的级间分解  (教材式 2-80 / 2-81)")
BX, BW = 110, 640
rowsA = [("有高放", c_hf, [("传输线", BLUE), ("高频小信号放大器", GREEN),
                          ("变频器", ORNG), ("中频放大器", RED)], NF_hf),
         ("无高放", c_no, [("传输线", BLUE), ("变频器", ORNG),
                          ("中频放大器", RED)], NF_no)]
yy = 128
for tag, cs, cols, tot in rowsA:
    dr.text((BX, yy), "%s    N_F = %.4f = %.2f dB    T_e = %.0f K"
            % (tag, tot, db(tot), (tot-1)*T0), font=f16, fill=(0, 0, 0))
    x = BX
    for (nm, col), v in zip(cols, cs):
        w = BW * v / tot
        dr.rectangle([x, yy+28, x+w, yy+62], fill=col, outline=(255, 255, 255))
        if w > 56:
            dr.text((x + w/2 - 17, yy+36), "%.1f%%" % (100*v/tot), font=f14, fill=(255, 255, 255))
        x += w
    ly = yy + 72
    for i in range(0, len(cols), 2):
        lx = BX
        for (nm, col), v in zip(cols[i:i+2], cs[i:i+2]):
            dr.rectangle([lx, ly, lx+13, ly+13], fill=col)
            t = "%s %.4f (%.2f%%)" % (nm, v, 100*v/tot)
            dr.text((lx+19, ly-1), t, font=f13, fill=(40, 40, 40))
            lx += dr.textlength(t, font=f13) + 60
        ly += 24
    yy += 134
dr.line([BX, 400, BX+BW, 400], fill=(200, 200, 200), width=1)
dr.text((BX, 416),
        "变频器噪声贡献:  有高放 %.4f  vs  无高放 %.4f   相差 %.0f 倍 (= 高放功率增益 A_PH2 = %.0f)"
        % (c_hf[2], c_no[1], c_no[1]/c_hf[2], G2), font=f15, fill=(180, 30, 30))
dr.text((BX, 444),
        "有高放时变频器段(橙)在条形中细窄到几乎看不见, 说明其噪声已被 60 倍前置增益淹没;",
        font=f14, fill=(60, 60, 60))
dr.text((BX, 466),
        "无高放时它一举占到 %.2f%%, 与中频放大器(%.2f%%)共同构成噪声主体。"
        % (100*c_no[1]/NF_no, 100*c_no[2]/NF_no), font=f14, fill=(60, 60, 60))

# ========== 面板B: N_F3 扫描 ==========
panel(890, 100, 1570, 540, "面板 B  总噪声系数随变频器噪声系数 N_F3 的变化")
PB = (955, 168, 1520, 335)
xlimB, ylimB = (1.0, 40.0), (2.0, 100.0)
axes(PB, "变频器噪声系数 N_F3 (线性)", "总噪声系数 N_F", xlimB, ylimB, logy=True)
xsB = [1.0 + 39.0*i/800.0 for i in range(801)]
plot(PB, xsB, [nf_hf(f) for f in xsB], xlimB, ylimB, BLUE, 3, True)
plot(PB, xsB, [nf_no(f) for f in xsB], xlimB, ylimB, RED, 3, True)
vline(PB, 6.0, xlimB, ORNG, "N_F3 = 6 (题给)", dy=30)
ty = 392
for t, c in [("有高放: N_F = %.5f + (N_F3 - 1)/48" % BASE_HF, BLUE),
             ("无高放: N_F = %.5f + (N_F3 - 1)/0.8" % BASE_NO, RED),
             ("灵敏度 dN_F/dN_F3 :   有高放 %.5f      无高放 %.4f" % (S_hf, S_no), (0, 0, 0)),
             ("N_F3 从 6 劣化到 60 (10 倍) 时:", (0, 0, 0)),
             ("     有高放 %.4f -> %.4f   (%.2f -> %.2f dB)   总噪声只恶化 %.2f dB"
              % (nf_hf(6.0), nf_hf(60.0), db(nf_hf(6.0)), db(nf_hf(60.0)), db(nf_hf(60.0))-db(nf_hf(6.0))), BLUE),
             ("     无高放 %.4f -> %.4f  (%.2f -> %.2f dB)  总噪声恶化 %.2f dB"
              % (nf_no(6.0), nf_no(60.0), db(nf_no(6.0)), db(nf_no(60.0)), db(nf_no(60.0))-db(nf_no(6.0))), RED)]:
    dr.text((PB[0]+14, ty), t, font=f14, fill=c)
    ty += 23

# ========== 面板C: 噪声预算 -> 允许 N_F3 ==========
panel(80, 590, 810, 1090, "面板 C  给定总噪声预算时, 变频器噪声系数的允许上限")
PC = (150, 655, 765, 845)
xlimC, ylimC = (2.0, 40.0), (1.0, 2000.0)
axes(PC, "总噪声预算 N_Fbudget (线性)", "允许上限 N_F3max", xlimC, ylimC, logy=True)
xsC = [2.0 + 38.0*i/400.0 for i in range(401)]
plot(PC, xsC, [max(1.0, 1 + (G1*G2)*(b - BASE_HF)) for b in xsC], xlimC, ylimC, BLUE, 3, True)
plot(PC, xsC, [max(1.0, 1 + (b - BASE_NO)/(1.0/G1)) for b in xsC], xlimC, ylimC, RED, 3, True)
vline(PC, BASE_HF, xlimC, GRAY)
vline(PC, BASE_NO, xlimC, GRAY)
ty = 900
for t, c in [("有高放: N_F3max = 1 + 48 (N_Fbudget - %.4f)" % BASE_HF, BLUE),
             ("无高放: N_F3max = 1 + 0.8 (N_Fbudget - %.2f)" % BASE_NO, RED),
             ("同一预算下两条线的间距, 就是高放带来的噪声要求宽裕度:", (0, 0, 0)),
             ("   预算 = %.0f (无高放实达值):  有高放 N_F3 <= %.0f ,  无高放 N_F3 <= %.0f"
              % (NF_no, 1 + (G1*G2)*(NF_no - BASE_HF), NF_MX), (0, 0, 0)),
             ("   预算 < %.4f (有高放) 或 < %.2f (无高放) 时该方案不可行:"
              % (BASE_HF, BASE_NO), (0, 0, 0)),
             ("   即使变频器理想 (N_F3 = 1) 也达不到预算", (0, 0, 0))]:
    dr.text((PC[0]+16, ty), t, font=f14, fill=c)
    ty += 23

# ========== 面板D: 解析 vs 蒙特卡洛 ==========
panel(890, 590, 1570, 1090, "面板 D  解析解与蒙特卡洛噪声仿真对照 (nsamp = 200000)")
yy = 610
dr.line([PD0 := 898, yy-6, 1562, yy-6], fill=(60, 60, 60), width=1)
for a, b, c, d, e in [("方案", "N_F3", "解析 N_F", "蒙特卡洛 N_F", "偏差")] + \
        [(t, "%.0f" % f3, "%.6f" % an, "%.6f" % mc, "%+.3f%%" % (100*(mc-an)/an)) for (t, f3, an, mc, gt) in MC]:
    fnt = f15 if a == "方案" else f14
    dr.text((910, yy), a, font=fnt, fill=(0, 0, 0))
    dr.text((990, yy), b, font=fnt, fill=(30, 30, 30))
    dr.text((1066, yy), c, font=fnt, fill=(30, 30, 30))
    dr.text((1216, yy), d, font=fnt, fill=(30, 30, 30))
    dr.text((1394, yy), e, font=fnt, fill=(0, 110, 0))
    yy += 30
dr.line([898, yy-6, 1562, yy-6], fill=(60, 60, 60), width=1)
dr.multiline_text((910, yy + 18),
    "仿真方法: 第 1 级输入端注入功率归一化为 1 的复高斯噪声(即参考噪声 kT0B, T0 = 290 K),\n"
    "各级按 N_eq,i = (F_i - 1) kT0B 独立注入白噪声, 逐级乘 sqrt(G_i) 后叠加; 末级输出噪声功率\n"
    "除以总功率增益即为噪声系数。全程不调用式(2-80) 闭式公式。\n"
    "式(2-80) 闭式值与逐级功率递推值互差 %.1e, 三条计算路径完全一致;\n"
    "蒙特卡洛残余偏差 %.3f%% ~ %.3f%% 属有限样本随机涨落(理论量级 1/sqrt(N) = %.3f%%)。"
    % (max(abs(NF_hf-NF_hf_friis), abs(NF_no-NF_no_friis)),
       min(100*(mc-an)/an for (_, _, an, mc, _) in MC),
       max(100*(mc-an)/an for (_, _, an, mc, _) in MC),
       100.0/math.sqrt(200000)),
    font=f14, fill=(0, 0, 0), spacing=9)

# ---------------- 结论 ----------------
dr.line([60, 1122, W-60, 1122], fill=(60, 60, 60), width=2)
dr.text((70, 1130), "结论", font=f18, fill=(0, 0, 0))
dr.multiline_text((70, 1160),
    "1) 接收机总噪声系数\n"
    "   有高放: N_F = %.2f + (%.0f-1)/%.1f + (%.0f-1)/(%.1f x %g) + (%.0f-1)/(%.1f x %g x %.1f) = %.4f = %.2f dB   (T_e = %.0f K)\n"
    "   无高放: N_F = %.2f + (%.0f-1)/%.1f + (%.0f-1)/(%.1f x %.1f) = %.4f = %.2f dB   (T_e = %.0f K)\n"
    "   => 无高放使整机噪声系数劣化 %.2f 倍 (恶化 %.2f dB)\n"
    "2) 变频器的噪声贡献\n"
    "   有高放 (N_F3-1)/(A_PH1 A_PH2) = %.0f/%g = %.4f, 只占总噪声 %.2f%%;\n"
    "   无高放 (N_F3-1)/A_PH1        = %.0f/%.1f = %.4f, 占总噪声 %.2f%% —— 相差 %.0f 倍, 恰等于高放的功率增益 A_PH2。\n"
    "3) 对变频噪声系数的要求\n"
    "   有高放: 变频器前面已有 %g 倍功率增益把它的噪声压下去, N_F3 每增大 1 总噪声仅增 %.5f;"
    "即使 N_F3 从 6 劣化到 60, 总噪声也只从 %.2f dB 升到 %.2f dB —— 要求非常宽松, 可用普通变频器。\n"
    "   无高放: 变频器升为接收机第一级有源器件, N_F3 直接进入总噪声的主导项, 每增大 1 总噪声就增 %.2f;"
    "从 6 劣化到 60, 总噪声从 %.2f dB 升到 %.2f dB; 即使 N_F3 = 1 总噪声仍有 %.2f = %.2f dB —— 必须用低噪声变频器。"
    % (NF_LINE, NF_HF, G1, NF_MX, G1, G2, NF_IF, G1, G2, G3, NF_hf, db(NF_hf), Te_hf,
       NF_LINE, NF_MX, G1, NF_IF, G1, G3, NF_no, db(NF_no), Te_no,
       NF_no/NF_hf, db(NF_no)-db(NF_hf),
       NF_MX-1, G1*G2, c_hf[2], 100*c_hf[2]/NF_hf,
       NF_MX-1, G1, c_no[1], 100*c_no[1]/NF_no, c_no[1]/c_hf[2],
       G2, S_hf, db(nf_hf(6.0)), db(nf_hf(60.0)),
       S_no, db(nf_no(6.0)), db(nf_no(60.0)), BASE_NO, db(BASE_NO)),
    font=f15, fill=(0, 0, 0), spacing=10)

img.save(r"E:\tmp\通信电路\ex2-14.png")
cnt = sum(1 for py in range(H) for pxx in range(W) if img.getpixel((pxx, py)) != (255, 255, 255))
print("非白像素 %d 个" % cnt)
print("已保存 E:\\tmp\\通信电路\\ex2-14.png  (%dx%d)" % (W, H))