from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "phase3b_patient_covariate_pairs.csv"
ORIG = ROOT / "data" / "phase3b_primary_results.csv"
RESULTS = ROOT / "results"

PRED = [
    "delta_hours_from_icu_admit", "delta_heart_rate", "delta_potassium",
    "delta_magnesium", "delta_calcium_total", "delta_creatinine"
]

def fit_ols(y, z, weights=None, hc3=False):
    n, k = z.shape
    if weights is None:
        weights = np.ones(n)
    sw = np.sqrt(weights)
    zw = z * sw[:, None]
    yw = y * sw
    bread = np.linalg.inv(zw.T @ zw)
    beta = bread @ (zw.T @ yw)
    resid = y - z @ beta
    if hc3:
        h = weights * np.einsum("ij,jk,ik->i", z, bread, z)
        h = np.minimum(h, 1 - 1e-10)
        score_scale = weights * resid / (1 - h)
        meat = z.T @ ((score_scale ** 2)[:, None] * z)
        cov = bread @ meat @ bread
        se = np.sqrt(np.diag(cov))
        df = n - k
    else:
        rss = np.sum(weights * resid**2)
        df = n-k
        cov = bread * (rss / df)
        se = np.sqrt(np.diag(cov))
    p = 2 * stats.t.sf(np.abs(beta / se), df=df)
    return beta, se, p, resid

def bh(p):
    p=np.asarray(p,float); out=np.full(len(p),np.nan)
    ok=np.isfinite(p); vals=p[ok]; order=np.argsort(vals); ranked=vals[order]
    adj=np.minimum.accumulate((ranked*len(vals)/np.arange(1,len(vals)+1))[::-1])[::-1]
    tmp=np.empty_like(adj); tmp[order]=np.minimum(adj,1); out[np.where(ok)[0]]=tmp
    return out

usecols=["outcome","stratum","window_h","subject_id","ingredient","role",
         "n_informative_stays","n_exposed_ecg","n_unexposed_ecg","delta_ms",
         *PRED,"mean_ecg_count_stay"]
chunks=[]
for c in pd.read_csv(DATA,usecols=usecols,chunksize=150000):
    q=c[(c.outcome=="qtcf_ms")&(c.stratum=="narrow_qrs")&(c.window_h==24)]
    chunks.append(q)
d=pd.concat(chunks,ignore_index=True)
orig=pd.read_csv(ORIG)

rows=[]
for drug,g in d.groupby("ingredient",sort=True):
    role=g.role.iloc[0]
    g=g.copy()
    g["log_ecg_monitoring"]=np.log1p(g.mean_ecg_count_stay)
    cols=PRED+["log_ecg_monitoring"]
    cc=g.dropna(subset=["delta_ms",*cols]).copy()
    status="fit" if len(cc)>=80 else "insufficient_complete_cases"
    base=dict(ingredient=drug,role=role,n_all=len(g),complete_cases=len(cc),model_status=status,
              unadjusted_all_ms=g.delta_ms.mean(),unadjusted_cc_ms=cc.delta_ms.mean())
    if status!="fit":
        rows.append(base); continue
    x=cc[cols].to_numpy(float)
    x=x-x.mean(axis=0)
    z=np.column_stack([np.ones(len(cc)),x]); y=cc.delta_ms.to_numpy(float)
    b,se,p,_=fit_ols(y,z,hc3=False)
    bhc,seh,ph,_=fit_ols(y,z,hc3=True)
    ne=cc.n_exposed_ecg.to_numpy(float); nu=cc.n_unexposed_ecg.to_numpy(float)
    w=ne*nu/(ne+nu)
    # Normalize (does not change WLS estimates or inference) and prespecify a
    # 99th-percentile cap to avoid allowing a handful of densely monitored
    # patients to dominate the sensitivity analysis.
    wcap=np.minimum(w,np.quantile(w,.99)); wcap=wcap/np.mean(wcap)
    # Weighted centering makes the intercept the weighted adjusted mean contrast.
    xw=cc[cols].to_numpy(float); xw=xw-np.average(xw,axis=0,weights=wcap)
    zw=np.column_stack([np.ones(len(cc)),xw])
    bw,sew,pw,_=fit_ols(y,zw,weights=wcap,hc3=True)
    base.update(original_ols_ms=b[0],original_se=se[0],original_p=p[0],
                hc3_adjusted_ms=bhc[0],hc3_se=seh[0],hc3_p=ph[0],
                hc3_ci_low=bhc[0]-stats.t.ppf(.975,len(cc)-z.shape[1])*seh[0],
                hc3_ci_high=bhc[0]+stats.t.ppf(.975,len(cc)-z.shape[1])*seh[0],
                precision_weighted_ms=bw[0],precision_weighted_hc3_se=sew[0],precision_weighted_p=pw[0],
                precision_weighted_ci_low=bw[0]-stats.t.ppf(.975,len(cc)-zw.shape[1])*sew[0],
                precision_weighted_ci_high=bw[0]+stats.t.ppf(.975,len(cc)-zw.shape[1])*sew[0],
                weight_p99=np.quantile(w,.99),weight_max=w.max())
    rows.append(base)

res=pd.DataFrame(rows)
for role in ["candidate","negative_control"]:
    ix=(res.role==role)&(res.model_status=="fit")
    res.loc[ix,"hc3_fdr_q"]=bh(res.loc[ix,"hc3_p"])
    res.loc[ix,"precision_weighted_fdr_q"]=bh(res.loc[ix,"precision_weighted_p"])

# Refit the empirical null using HC3 standard errors for all estimable controls.
nc=res[(res.role=="negative_control")&(res.model_status=="fit")].copy()
theta=nc.hc3_adjusted_ms.to_numpy(float); tau=nc.hc3_se.to_numpy(float)
def nll(par):
    mu=par[0]; sigma=np.exp(par[1]); v=tau*tau+sigma*sigma
    return .5*np.sum(np.log(2*np.pi*v)+(theta-mu)**2/v)
opt=minimize(nll,[theta.mean(),np.log(max(theta.std(ddof=1),.01))],method="BFGS")
mu=float(opt.x[0]); sigma=float(np.exp(opt.x[1]))
res["hc3_calibrated_ms"]=res.hc3_adjusted_ms-mu
res["hc3_calibrated_se"]=np.sqrt(res.hc3_se**2+sigma**2)
res["hc3_calibrated_ci_low"]=res.hc3_calibrated_ms-1.96*res.hc3_calibrated_se
res["hc3_calibrated_ci_high"]=res.hc3_calibrated_ms+1.96*res.hc3_calibrated_se
res["hc3_calibrated_p"]=2*stats.norm.sf(np.abs(res.hc3_calibrated_ms/res.hc3_calibrated_se))
ix=(res.role=="candidate")&(res.model_status=="fit")
res.loc[ix,"hc3_calibrated_fdr_q"]=bh(res.loc[ix,"hc3_calibrated_p"])

cmp=orig.merge(res,on=["ingredient","role"],how="left")
cmp["ols_difference_vs_original"]=cmp.original_ols_ms-cmp.adjusted_delta_ms
RESULTS.mkdir(exist_ok=True)
res.to_csv(RESULTS / "phase3b_hc3_precision_results.csv",index=False)
cmp.to_csv(RESULTS / "phase3b_replication_check.csv",index=False)
print("Primary rows",len(d),"drugs",d.ingredient.nunique())
print("Max abs OLS replication difference",np.nanmax(np.abs(cmp.ols_difference_vs_original)))
print("HC3 empirical null mu",mu,"sigma",sigma,"optimizer_success",opt.success)
print(res[["ingredient","role","complete_cases","unadjusted_all_ms","unadjusted_cc_ms","hc3_adjusted_ms","hc3_ci_low","hc3_ci_high","hc3_p","hc3_fdr_q","precision_weighted_ms","precision_weighted_ci_low","precision_weighted_ci_high"]].to_string(index=False))
