# -*- coding: utf-8 -*-
"""生成教材图 1-3「调幅无线广播发射机的基本组成方框图」：文本版 + PNG 版。
关键点：CJK 字符在等宽终端占 2 列，必须按显示宽度对齐，不能用 len()。
"""
import unicodedata
from PIL import Image, ImageDraw, ImageFont

OUT_DIR = r"E:\tmp\通信电路"

# ---------- 显示宽度工具 ----------
def dw(s):
    """字符串在等宽终端的显示宽度（东亚全角/宽字符 = 2 列）。"""
    return sum(2 if unicodedata.east_asian_width(c) in ("W", "F") else 1 for c in s)

def pad(s, w):
    return s + " " * max(0, w - dw(s))


# ================= 文本版 =================
def build_text():
    W, H = 136, 14
    g = [[" "] * W for _ in range(H)]

    def put(x, y, s):
        """写入文本。x,y 以「显示列」计；CJK 字符占 2 列，故按宽度推进。"""
        for ch in s:
            if 0 <= y < H and 0 <= x < W:
                g[y][x] = ch
            x += 2 if unicodedata.east_asian_width(ch) in ("W", "F") else 1

    def hline(y, x1, x2, ch="\u2500"):
        for x in range(min(x1, x2), max(x1, x2) + 1):
            put(x, y, ch)

    def vline(x, y1, y2, ch="\u2502"):
        for y in range(min(y1, y2), max(y1, y2) + 1):
            put(x, y, ch)

    def box(x, y, w, lines, style="round"):
        """画方框并居中写文本。x,y 为左上角；w 为内部宽度（不含左右竖边）。"""
        tl, mid, tr, bl, br = ("\u256d", "\u2502", "\u256e", "\u2570", "\u256f")
        hline(y, x + 1, x + w, "\u2500")
        put(x, y, tl); put(x + w, y, tr)
        for j, ln in enumerate(lines):
            yy = y + 1 + j
            put(x, yy, mid); put(x + w, yy, mid)
            put(x + 1 + (w - 1 - dw(ln)) // 2, yy, ln)
        yb = y + len(lines) + 1
        hline(yb, x + 1, x + w, "\u2500")
        put(x, yb, bl); put(x + w, yb, br)
        return yb

    # ---- 上排：载波链 ----
    y0 = 1
    boxes = [
        ("正弦波振荡器",), ("缓冲器",), ("高频电压放大器",),
        ("振幅", "调制器"), ("高频电压及功率放大器",), ("发射天线",),
    ]
    xs, ws = [], []
    x = 2
    for b in boxes:
        w = max(dw(t) for t in b) + 6
        xs.append(x); ws.append(w)
        box(x, y0, w, list(b))
        x += w
    y1 = y0 + 3                      # 上排方框底线
    for i in range(len(boxes) - 1):
        ax, bx = xs[i] + ws[i], xs[i + 1]
        # ax/bx 是竖边所在列，连线必须从 ax+1 画到 bx-1，否则会覆盖方框边框
        hline(y0 + 1, ax + 1, bx - 2)
        put(bx - 1, y0 + 1, "\u25b6")

    # 振幅调制器的引出点（底边中点）
    i4 = 3
    tap_x = xs[i4] + ws[i4] // 2

    # ---- 下排：音频支路 ----
    y2 = y1 + 5
    a1x, a1w = 8, dw("声/电") + 8
    a2x, a2w = a1x + a1w + 4, dw("低频电压放大器") + 6
    box(a1x, y2, a1w, ["声/电", "变换器"])
    box(a2x, y2, a2w, ["低频电压放大器"])
    y2b = y2 + 3
    hline(y2 + 1, a1x + a1w + 1, a2x - 2)
    put(a2x - 1, y2 + 1, "\u25b6")

    # 折线：低频放大器输出 -> 向右到 tap_x -> 向上进入振幅调制器
    y_hor = y2 + 1
    hline(y_hor, a2x + a2w + 1, tap_x - 1)   # 先画横线，止于拐角前一列
    vline(tap_x, y1 + 2, y_hor)              # 再画竖线
    put(tap_x, y_hor, "\u251c")              # 最后落拐角，避免被覆盖
    put(tap_x, y1 + 1, "\u25b2")

    # 渲染：网格按「显示列」存放；输出时遇到宽字符直接跳过其占用的续列，
    # 否则 CJK 后面那个占位空格会被当成真实空格，宽度就错开了。
    out = []
    for row in g:
        line, x = [], 0
        while x < W:
            ch = row[x]
            line.append(ch)
            x += 2 if ch != " " and unicodedata.east_asian_width(ch) in ("W", "F") else 1
        out.append("".join(line).rstrip())
    return "\n".join(out)


# ================= 图片版 =================
def build_png(path):
    """按字体实测宽度自动定框，再排版，避免文字溢出方框。"""
    S = 2                                        # 超采样
    FS = 28                                      # 正文字号（成品像素）
    PAD, GAP, MARGIN = 34, 72, 50
    BH = 118                                     # 方框高度

    probe = ImageDraw.Draw(Image.new("RGB", (10, 10)))
    F  = ImageFont.truetype(r"C:\Windows\Fonts\msyh.ttc", FS * S)
    FB = ImageFont.truetype(r"C:\Windows\Fonts\msyhbd.ttc", FS * S)
    FT = ImageFont.truetype(r"C:\Windows\Fonts\msyhbd.ttc", 36 * S)
    FSUB = ImageFont.truetype(r"C:\Windows\Fonts\msyh.ttc", 22 * S)

    top = [["正弦波振荡器"], ["缓冲器"], ["高频电压放大器"],
           ["振幅", "调制器"], ["高频电压及功率放大器"], ["发射天线"]]
    bot = [["声/电", "变换器"], ["低频电压放大器"]]

    def width_of(lines, f):
        return max(probe.textlength(l, font=f) for l in lines) / S + 2 * PAD

    w_top = [width_of(b, FB) for b in top]
    w_bot = [width_of(b, F) for b in bot]

    W = MARGIN * 2 + sum(w_top) + GAP * (len(top) - 1)
    W = int(max(W, MARGIN * 2 + sum(w_bot) + GAP)) 
    H = 660

    img = Image.new("RGB", (int(W) * S, H * S), "white")
    d = ImageDraw.Draw(img)

    def rect(x, y, w, h, lines, bold=False):
        f = FB if bold else F
        d.rectangle([x * S, y * S, (x + w) * S, (y + h) * S],
                    fill="white", outline="black", width=2 * S)
        lh = FS * 1.5
        y0 = y * S + (h * S - lh * S * len(lines)) / 2
        for j, ln in enumerate(lines):
            tw = d.textlength(ln, font=f)
            d.text((x * S + (w * S - tw) / 2, y0 + lh * S * j), ln, font=f, fill="black")

    def arrow_h(x1, x2, y, head=15):
        d.line([x1 * S, y * S, (x2 - head) * S, y * S], fill="black", width=2 * S)
        d.polygon([(x2 * S, y * S), ((x2 - head) * S, (y - head * .6) * S),
                   ((x2 - head) * S, (y + head * .6) * S)], fill="black")

    def arrow_v(x, y1, y2, head=15):
        d.line([x * S, y1 * S, x * S, (y2 + head) * S], fill="black", width=2 * S)
        d.polygon([(x * S, y2 * S), ((x - head * .6) * S, (y2 + head) * S),
                   ((x + head * .6) * S, (y2 + head) * S)], fill="black")

    # 上排：载波链
    ty = 200
    xs, x = [], MARGIN
    for b, w in zip(top, w_top):
        rect(x, ty, w, BH, b, bold=True)
        xs.append((x, w)); x += w + GAP
    for i in range(len(xs) - 1):
        (ax, aw), (bx, _) = xs[i], xs[i + 1]
        arrow_h(ax + aw + 10, bx - 10, ty + BH / 2)

    # 下排：音频支路
    by = ty + BH + 160
    a1x = xs[1][0]
    a2x = a1x + w_bot[0] + GAP
    rect(a1x, by, w_bot[0], BH, bot[0])
    rect(a2x, by, w_bot[1], BH, bot[1])
    arrow_h(a1x + w_bot[0] + 10, a2x - 10, by + BH / 2)

    # 折线接入振幅调制器
    tx = xs[3][0] + xs[3][1] / 2
    y_mid = by + BH / 2
    d.line([(a2x + w_bot[1] + 10) * S, y_mid * S, tx * S, y_mid * S],
           fill="black", width=2 * S)
    arrow_v(tx, y_mid, ty + BH + 10)

    d.text((MARGIN * S, 55 * S), "图 1-3  调幅无线广播发射机的基本组成方框图",
           font=FT, fill="black")
    d.text((MARGIN * S, 120 * S), "载波链（上排，左→右）/ 音频支路（下排，注入振幅调制器）",
           font=FSUB, fill=(95, 95, 95))

    img.resize((int(W), H), Image.LANCZOS).save(path, "PNG")
    return path

if __name__ == "__main__":
    import os, sys
    os.makedirs(OUT_DIR, exist_ok=True)
    txt = build_text()

    # ---- 结构校验：肉眼数不出「显示列」，只能程序化核对 ----
    def dcol(line, idx):
        """字符下标 -> 显示列。"""
        return dw(line[:idx])

    errs = []
    lines = txt.splitlines()
    col_up = col_corner = None
    for row in lines:
        if "▲" in row:
            col_up = dcol(row, row.index("▲"))
        if "├" in row:
            col_corner = dcol(row, row.index("├"))

    # 振幅调制器底边的中心列
    bottom = [r for r in lines if r.count("╰") == 1 and "─" in r and "╯" in r
              and r.index("╰") < r.index("╯") and "调制" not in r]
    if bottom:
        r = bottom[0]
        c0 = dcol(r, r.index("╰"))
        c1 = dcol(r, r.rindex("╯"))
        center = (c0 + c1) // 2
        if col_up != center:
            errs.append(f"入射箭头列 {col_up} != 振幅调制器中心列 {center}")
        if col_corner != center:
            errs.append(f"折线拐角列 {col_corner} != 振幅调制器中心列 {center}")

    # 每个方框顶边/底边的显示宽度必须与其中框线一致
    for n, row in enumerate(lines, 1):
        if row.strip().startswith("╭") or "╭" in row:
            for seg in row.split("╭")[1:]:
                seg = seg.split("╮")[0]
                if seg and set(seg) != {"─"}:
                    errs.append(f"L{n}: 顶边有缺口 {seg!r}")
                    break

    if errs:
        print("校验未通过:", *errs, sep="\n  ")
        sys.exit(1)

    with open(os.path.join(OUT_DIR, "fig1-3.txt"), "w", encoding="utf-8") as f:
        f.write(txt + "\n")
    p = build_png(os.path.join(OUT_DIR, "fig1-3.png"))
    print("校验通过：入射列 =", col_up, " 拐角列 =", col_corner)
    print("PNG ->", p)