from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
FIGURES.mkdir(exist_ok=True)

d = pd.read_csv(RESULTS / "phase3b_hc3_precision_results.csv")
d = d[d.model_status == "fit"].copy()

# Figure 2: candidate drugs and negative controls on a common scale.
c = d[d.role == "candidate"].sort_values("hc3_adjusted_ms", ascending=False)
n = d[d.role == "negative_control"].sort_values("hc3_adjusted_ms", ascending=False)
plot = pd.concat([c, n], ignore_index=True)
y = np.arange(len(plot))[::-1]
colors = np.where(plot.role.eq("candidate"), "#B51F1F", "#174F7A")
markers = np.where(plot.role.eq("candidate"), "o", "s")
fig, ax = plt.subplots(figsize=(8.4, 6.6))
ax.axvspan(-3.81, -0.29, color="#DCE5EF", zorder=0)
ax.axvline(0, color="#555555", ls="--", lw=1)
ax.axvline(-2.05, color="#174F7A", ls=":", lw=1.2)
for i, row in plot.iterrows():
    ax.errorbar(row.hc3_adjusted_ms, y[i],
                xerr=[[row.hc3_adjusted_ms-row.hc3_ci_low],
                      [row.hc3_ci_high-row.hc3_adjusted_ms]],
                fmt=markers[i], color=colors[i], capsize=0, ms=5, lw=1.2)
ax.set_yticks(y, [x.capitalize() for x in plot.ingredient])
ax.set_xlabel("Adjusted mean exposed-minus-unexposed QTcF difference (ms)")
ax.set_title("Candidate drugs and negative controls on a common scale", loc="left", weight="bold")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIGURES / "Figure_2_recreated.tiff", dpi=600, pil_kwargs={"compression":"tiff_lzw"})
plt.close(fig)

# Figure 3: empirical calibration diagnostic.
nc = n.copy()
se = nc.hc3_se.to_numpy(float)
theta = nc.hc3_adjusted_ms.to_numpy(float)
x = np.linspace(-10.5, 6.5, 500)
mu, sigma = -2.05, 0.90
upper_cal = mu + 1.96*np.sqrt(sigma*sigma + x*0 + np.linspace(0,3.1,500)**2)
lower_cal = mu - 1.96*np.sqrt(sigma*sigma + np.linspace(0,3.1,500)**2)
sy = np.linspace(0,3.1,500)
fig, ax = plt.subplots(figsize=(8.4, 6.2))
ax.fill_betweenx(sy, lower_cal, upper_cal, color="#DCE5EF", label="Calibrated 95% region")
ax.plot(-1.96*sy, sy, "--", color="#777777", label="Uncalibrated 95% region")
ax.plot(1.96*sy, sy, "--", color="#777777")
ax.axvline(mu, color="#174F7A", ls=":")
ax.axvline(0, color="#555555", ls="--")
ax.scatter(theta, se, marker="s", s=52, color="#174F7A", zorder=3)
for drug, xx, yy in zip(nc.ingredient, theta, se):
    ax.annotate(drug.capitalize(), (xx, yy), xytext=(5, 0), textcoords="offset points", va="center", fontsize=8)
ax.set_xlim(-10.5,6.5); ax.set_ylim(0,3.1)
ax.set_xlabel("Adjusted QTcF difference (ms)"); ax.set_ylabel("Standard error (ms)")
ax.set_title("Empirical calibration against the negative controls", loc="left", weight="bold")
ax.legend(frameon=False, loc="upper left")
ax.spines[["top","right"]].set_visible(False)
fig.tight_layout()
fig.savefig(FIGURES / "Figure_3_recreated.tiff", dpi=600, pil_kwargs={"compression":"tiff_lzw"})
plt.close(fig)

# Supplementary Figure S11: aggregate drug-selection flow.
fig, ax = plt.subplots(figsize=(7.2, 8.5))
ax.set_xlim(0, 10); ax.set_ylim(0, 12); ax.axis("off")
def box(x, y, w, h, text, color="#EAF2F8"):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.03",facecolor=color,edgecolor="#59728A"))
    ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=8)
def arrow(x1,y1,x2,y2):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="->",mutation_scale=10,color="#6B7280"))
box(2,10.5,6,1,"Prespecified drug map\n35 drugs: 27 candidates and 8 negative controls","#FFFFFF")
box(.7,8.7,3.8,1,"Candidate drugs\n27 mapped")
box(5.5,8.7,3.8,1,"Negative controls\n8 mapped","#F2F3F4")
box(.7,6.8,3.8,1,"Administration identified\n26 candidates; droperidol excluded")
box(5.5,6.8,3.8,1,"Administration identified\n8 controls","#F2F3F4")
box(.7,4.6,3.8,1.4,"Feasibility threshold met\n16 candidates\n10 below 200 discordant pairs")
box(5.5,4.8,3.8,1,"Feasibility threshold met\n8 controls","#F2F3F4")
box(.7,2.7,3.8,1,"Primary analysis\n16 candidate-drug contrasts")
box(5.5,2.7,3.8,1,"Primary analysis\n8 negative controls","#F2F3F4")
box(.7,.7,3.8,1.2,"13 estimable models\n3 below 80 complete cases","#FDEBD0")
for a in [(5,10.5,2.6,9.7),(5,10.5,7.4,9.7),(2.6,8.7,2.6,7.8),(7.4,8.7,7.4,7.8),(2.6,6.8,2.6,6.0),(7.4,6.8,7.4,5.8),(2.6,4.6,2.6,3.7),(7.4,4.8,7.4,3.7),(2.6,2.7,2.6,1.9)]: arrow(*a)
ax.set_title("Supplementary Figure S11. Drug selection and analytic flow", loc="left", weight="bold")
fig.tight_layout()
fig.savefig(FIGURES / "Supplementary_Figure_S11_recreated.tiff", dpi=600, pil_kwargs={"compression":"tiff_lzw"})
plt.close(fig)
