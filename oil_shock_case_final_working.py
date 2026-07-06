import json
import os
import time
import urllib.error
import urllib.request

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from statsmodels.tsa.api import VAR


# --------------------
# Settings
# --------------------
START_DATE = "2003-01-01"
END_DATE = pd.Timestamp.today().strftime("%Y-%m-%d")
MAX_LAGS = 12
IRF_HORIZON = 18
MC_DRAWS = 1000
SHOCK_Z = 1.5

os.makedirs("outputs/tables", exist_ok=True)


# --------------------
# BCB helpers (robust)
# --------------------
def _fetch_bcb_raw(code, start_date, end_date, retries=3):
    start_br = pd.to_datetime(start_date).strftime("%d/%m/%Y")
    end_br = pd.to_datetime(end_date).strftime("%d/%m/%Y")
    url = (
        f"https://api.bcb.gov.br/dados/serie/bcdata.sgs.{code}/dados"
        f"?formato=json&dataInicial={start_br}&dataFinal={end_br}"
    )

    headers = {
        "User-Agent": "Mozilla/5.0",
        "Accept": "application/json",
    }

    last_error = None

    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=30) as resp:
                raw = resp.read().decode("utf-8")

            if not raw or not raw.strip():
                raise ValueError(f"empty response from BCB for code={code}")

            data = json.loads(raw)

            if isinstance(data, dict) and ("error" in data or "message" in data):
                raise RuntimeError(str(data))

            return data

        except urllib.error.HTTPError:
            raise
        except Exception as e:
            last_error = e
            time.sleep(1 + i)

    raise last_error if last_error else RuntimeError("unknown BCB fetch error")


def _fetch_bcb_chunked(code, col_name, start_date, end_date, chunk_years=9):
    start = pd.to_datetime(start_date)
    end = pd.to_datetime(end_date)
    parts = []

    while start <= end:
        stop = min(start + pd.DateOffset(years=chunk_years) - pd.Timedelta(days=1), end)
        data = _fetch_bcb_raw(code, start, stop)

        if data:
            part = pd.DataFrame(data)
            part["date"] = pd.to_datetime(part["data"], dayfirst=True, errors="coerce")
            part[col_name] = pd.to_numeric(part["valor"], errors="coerce")
            part = part[["date", col_name]].dropna()
            parts.append(part)

        start = stop + pd.Timedelta(days=1)

    if not parts:
        return pd.DataFrame(columns=["date", col_name])

    out = (
        pd.concat(parts, ignore_index=True)
        .drop_duplicates(subset="date")
        .sort_values("date")
        .reset_index(drop=True)
    )
    return out


def get_bcb_series(code, col_name, start_date, end_date):
    try:
        data = _fetch_bcb_raw(code, start_date, end_date)
        df = pd.DataFrame(data)

        if df.empty:
            return pd.DataFrame(columns=["date", col_name])

        df["date"] = pd.to_datetime(df["data"], dayfirst=True, errors="coerce")
        df[col_name] = pd.to_numeric(df["valor"], errors="coerce")
        return df[["date", col_name]].dropna().sort_values("date").reset_index(drop=True)

    except urllib.error.HTTPError as e:
        if e.code == 406:
            return _fetch_bcb_chunked(code, col_name, start_date, end_date, chunk_years=9)
        raise

    except RuntimeError as e:
        msg = str(e).lower()
        if "10 anos" in msg or "periodicidade diária" in msg:
            return _fetch_bcb_chunked(code, col_name, start_date, end_date, chunk_years=9)
        raise


# --------------------
# FRED + transforms
# --------------------
def get_fred_series(series_id, col_name, start_date, end_date):
    url = f"https://fred.stlouisfed.org/graph/fredgraph.csv?id={series_id}"
    raw = pd.read_csv(url)

    date_col = "DATE" if "DATE" in raw.columns else "observation_date"
    value_col = [c for c in raw.columns if c != date_col][0]

    df = raw.rename(columns={date_col: "date", value_col: col_name})
    df["date"] = pd.to_datetime(df["date"])
    df[col_name] = pd.to_numeric(df[col_name], errors="coerce")

    mask = (df["date"] >= pd.to_datetime(start_date)) & (df["date"] <= pd.to_datetime(end_date))
    return df.loc[mask, ["date", col_name]].dropna().sort_values("date")


def daily_to_monthly_mean(df, col_name):
    out = df.set_index("date").resample("MS").mean().dropna().reset_index()
    return out[["date", col_name]]


# --------------------
# VAR helpers
# --------------------
def choose_lag(model, nobs, neq, max_lags=12):
    # safety cap to avoid "largest model cannot be estimated"
    safe_max = max(1, min(max_lags, int((nobs - 1) / (neq + 1))))

    for p in range(safe_max, 0, -1):
        try:
            selected = model.select_order(maxlags=p).aic
            if selected is None or (isinstance(selected, float) and np.isnan(selected)):
                return min(4, p)
            return max(1, int(selected))
        except Exception:
            continue

    return 1


def run_var(df, variables, max_lags=12):
    model = VAR(df[variables])
    nobs = len(df)
    neq = len(variables)
    lag = choose_lag(model, nobs, neq, max_lags)
    result = model.fit(lag)
    return result, lag


def get_irf_df(result, impulse, response, horizon=18, mc_draws=1000):
    irf = result.irf(horizon)
    i = result.names.index(impulse)
    r = result.names.index(response)

    values = irf.orth_irfs[:, r, i]
    low, high = irf.errband_mc(orth=True, repl=mc_draws, signif=0.05, seed=42)

    return pd.DataFrame(
        {
            "h": np.arange(horizon + 1),
            "irf": values,
            "low": low[:, r, i],
            "high": high[:, r, i],
        }
    )


def get_fevd_df(result, response_var, horizon=18):
    fevd_obj = result.fevd(horizon)
    arr = fevd_obj.decomp
    r = result.names.index(response_var)

    if arr.shape[0] == len(result.names):
        mat = arr[r, :, :]
    else:
        mat = arr[:, r, :]

    out = pd.DataFrame(mat, columns=result.names)
    out.insert(0, "h", np.arange(1, len(out) + 1))
    return out


def irf_summary(irf_df, label):
    temp = irf_df[irf_df["h"] >= 1].copy()
    peak_idx = temp["irf"].idxmax()

    return {
        "model": label,
        "peak_month": int(temp.loc[peak_idx, "h"]),
        "peak_response_pp": float(temp.loc[peak_idx, "irf"]),
        "cum_12m_pp": float(temp[temp["h"] <= 12]["irf"].sum()),
    }


# --------------------
# Main workflow
# --------------------
print("[1/6] Downloading data...")
ipca = get_bcb_series(433, "ipca_mom", START_DATE, END_DATE)
usdbrl_daily = get_bcb_series(1, "usdbrl", START_DATE, END_DATE)
selic_daily = get_bcb_series(432, "selic", START_DATE, END_DATE)
ibcbr = get_bcb_series(24363, "ibcbr", START_DATE, END_DATE)
brent_daily = get_fred_series("DCOILBRENTEU", "brent", START_DATE, END_DATE)
wti_daily = get_fred_series("DCOILWTICO", "wti", START_DATE, END_DATE)

print("Rows:", {"ipca": len(ipca), "usdbrl_daily": len(usdbrl_daily), "selic_daily": len(selic_daily), "ibcbr": len(ibcbr), "brent_daily": len(brent_daily), "wti_daily": len(wti_daily)})

# convert daily series to monthly BEFORE merge
usdbrl = daily_to_monthly_mean(usdbrl_daily, "usdbrl")
selic = daily_to_monthly_mean(selic_daily, "selic")
brent = daily_to_monthly_mean(brent_daily, "brent")
wti = daily_to_monthly_mean(wti_daily, "wti")

print("[2/6] Building monthly panel...")
df = (
    ipca.merge(selic, on="date", how="inner")
    .merge(ibcbr, on="date", how="inner")
    .merge(usdbrl, on="date", how="inner")
    .merge(brent, on="date", how="inner")
    .merge(wti, on="date", how="inner")
    .sort_values("date")
    .reset_index(drop=True)
)

# Transformations on aligned monthly frame
df["ipca"] = df["ipca_mom"]
df["oil_ret_brent"] = 100 * np.log(df["brent"]).diff()
df["oil_ret_wti"] = 100 * np.log(df["wti"]).diff()
df["fx_ret"] = 100 * np.log(df["usdbrl"]).diff()
df["d_selic"] = df["selic"].diff()
df["ibc_growth"] = 100 * np.log(df["ibcbr"]).diff().shift(1)

z = (df["oil_ret_brent"] - df["oil_ret_brent"].mean()) / df["oil_ret_brent"].std(ddof=0)
df["oil_shock_dummy"] = (z >= SHOCK_Z).astype(int)

# Baseline sample
df_base = df[["date", "ipca", "oil_ret_brent", "fx_ret", "ibc_growth", "d_selic", "oil_shock_dummy"]].dropna().copy()
df_base = df_base.rename(columns={"oil_ret_brent": "oil_ret"})

# Robustness sample
df_rob = df[["date", "ipca", "oil_ret_wti", "fx_ret", "ibc_growth", "d_selic"]].dropna().copy()
df_rob = df_rob.rename(columns={"oil_ret_wti": "oil_ret"})

# common sample for fair comparison
common_dates = set(df_base["date"]).intersection(set(df_rob["date"]))
df_base = df_base[df_base["date"].isin(common_dates)].sort_values("date")
df_rob = df_rob[df_rob["date"].isin(common_dates)].sort_values("date")

base_model_df = df_base.drop(columns=["date", "oil_shock_dummy"])
rob_model_df = df_rob.drop(columns=["date"])

print("Samples:", {"base": len(base_model_df), "rob": len(rob_model_df)})

print("[3/6] Estimating VAR models...")
vars_with_fx = ["oil_ret", "fx_ret", "ibc_growth", "d_selic", "ipca"]
vars_no_fx = ["oil_ret", "ibc_growth", "d_selic", "ipca"]

res_base, lag_base = run_var(base_model_df, vars_with_fx, MAX_LAGS)
res_nofx, lag_nofx = run_var(base_model_df, vars_no_fx, MAX_LAGS)
res_rob, lag_rob = run_var(rob_model_df, vars_with_fx, MAX_LAGS)

print(f"Lags selected -> baseline: {lag_base}, no-FX: {lag_nofx}, robustness: {lag_rob}")

print("[4/6] IRF and FEVD...")
irf_base = get_irf_df(res_base, "oil_ret", "ipca", horizon=IRF_HORIZON, mc_draws=MC_DRAWS)
irf_nofx = get_irf_df(res_nofx, "oil_ret", "ipca", horizon=IRF_HORIZON, mc_draws=MC_DRAWS)
irf_rob = get_irf_df(res_rob, "oil_ret", "ipca", horizon=IRF_HORIZON, mc_draws=MC_DRAWS)
fevd_base = get_fevd_df(res_base, "ipca", horizon=IRF_HORIZON)

print("[5/6] Plotly charts...")
fig1 = go.Figure()
fig1.add_trace(go.Scatter(x=df_base["date"], y=df_base["oil_ret"], mode="lines", name="Oil return (% m/m)"))
shock_points = df_base[df_base["oil_shock_dummy"] == 1]
fig1.add_trace(go.Scatter(x=shock_points["date"], y=shock_points["oil_ret"], mode="markers", name=f"Shock months (z >= {SHOCK_Z})"))
fig1.update_layout(template="plotly_white", title="Oil shocks (Brent)", xaxis_title="Date", yaxis_title="% m/m")
fig1.show()

fig2 = go.Figure()
fig2.add_trace(go.Scatter(x=irf_base["h"], y=irf_base["irf"], mode="lines", name="Baseline (with FX)"))
fig2.add_trace(go.Scatter(x=irf_nofx["h"], y=irf_nofx["irf"], mode="lines", name="No-FX model"))
fig2.add_trace(go.Scatter(x=irf_rob["h"], y=irf_rob["irf"], mode="lines", name="Robustness (WTI)"))
fig2.add_trace(go.Scatter(x=pd.concat([irf_base["h"], irf_base["h"][::-1]]), y=pd.concat([irf_base["high"], irf_base["low"][::-1]]), fill="toself", name="95% CI (baseline)", line=dict(width=0), opacity=0.2))
fig2.add_hline(y=0, line_dash="dash", line_width=1)
fig2.update_layout(template="plotly_white", title="Inflation response to oil shock (IRF)", xaxis_title="Months after shock", yaxis_title="IPCA response (p.p.)")
fig2.show()

fig3 = go.Figure()
for col in [c for c in fevd_base.columns if c != "h"]:
    fig3.add_trace(go.Scatter(x=fevd_base["h"], y=100 * fevd_base[col], mode="lines", stackgroup="one", name=col))
fig3.update_layout(template="plotly_white", title="FEVD of inflation (baseline)", xaxis_title="Horizon (months)", yaxis_title="Variance share (%)")
fig3.show()

print("[6/6] Saving tables...")
summary = pd.DataFrame([
    irf_summary(irf_base, "Baseline (with FX)"),
    irf_summary(irf_nofx, "No-FX model"),
    irf_summary(irf_rob, "Robustness (WTI)"),
])

base_12 = summary.loc[summary["model"] == "Baseline (with FX)", "cum_12m_pp"].iloc[0]
nofx_12 = summary.loc[summary["model"] == "No-FX model", "cum_12m_pp"].iloc[0]
fx_channel_proxy = base_12 - nofx_12

summary.to_csv("outputs/tables/irf_summary.csv", index=False)
irf_base.to_csv("outputs/tables/irf_baseline.csv", index=False)
irf_nofx.to_csv("outputs/tables/irf_no_fx.csv", index=False)
irf_rob.to_csv("outputs/tables/irf_robust_wti.csv", index=False)
fevd_base.to_csv("outputs/tables/fevd_baseline.csv", index=False)

print("\nSummary:")
print(summary.to_string(index=False))
print(f"\nFX channel proxy (12m cumulative baseline - noFX): {fx_channel_proxy:.4f} p.p.")
print("\nDone.")
