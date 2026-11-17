#!/usr/bin/env python3
"""Launch-day Twitch CCU forecaster (make.fun gta6-twitch-ccu-forecast-01).

Reads one case (schema makefun/twitch-ccu-case/v1): the launch's known
hour-1 viewer count plus a library of past launches' hourly curves
("comparables"). For each of hours 2..24 it writes five quantiles of that
hour's viewer count.

Method: every comparable contributes its release-day log-ratio curve
ln(viewers_h / viewers_1). A comparable's vote is weighted by how close its
own hour-1 count is to this case's hour-1 count (Gaussian kernel on
log10 scale, bandwidth 0.5 dex): launches of a similar size trace the
shape this launch is likely to follow. The forecast for hour h is the
case's hour-1 count times exp(weighted quantile of the log-ratios).
Early-access-contaminated launches (their release day is measured after a
paid head start, so the shape is not a launch shape) are left out of the
fit. Hours with too little weighted evidence fall back first to the
unweighted pool, then to a persistence forecast with a widening spread.

Deterministic: the output is a pure function of the case file.
Standard library only; no network, no clock, no randomness.
"""

import json
import math
import os
import sys
from datetime import datetime

LEVELS_DEFAULT = [0.05, 0.25, 0.5, 0.75, 0.95]
BANDWIDTH_DEX = 0.5      # kernel width on log10(hour-1 viewers)
MIN_WEIGHT = 0.02        # drop comparables with a smaller kernel weight
MIN_HOUR1 = 300          # a comparable needs a real hour 1 to have a shape
MIN_POINTS = 5           # weighted points needed before trusting an hour
TAIL_CAP = 2.2           # no quantile sits further than this (log units)
                         # from the hour's median: trims the wildest
                         # ratio outliers without narrowing honest spreads


def parse_ts(text):
    return datetime.fromisoformat(text.replace("Z", "+00:00"))


def release_day_hours(comparable):
    """Map hour 1..24 of the comparable's release day to its viewers."""
    try:
        start = parse_ts(comparable["day_start"])
    except Exception:
        return {}
    hours = {}
    for point in comparable.get("hourly") or []:
        try:
            delta = (parse_ts(point["from"]) - start).total_seconds()
            viewers = point["viewers"]
        except Exception:
            continue
        if 0 <= delta < 24 * 3600 and viewers is not None:
            hour = int(delta // 3600) + 1
            if hour not in hours:
                hours[hour] = viewers
    return hours


def weighted_quantile(pairs, level):
    """Quantile of (value, weight) pairs by cumulative weight."""
    pairs = sorted(pairs)
    total = sum(w for _, w in pairs)
    if total <= 0:
        return pairs[len(pairs) // 2][0]
    acc = 0.0
    for value, weight in pairs:
        acc += weight
        if acc / total >= level:
            return value
    return pairs[-1][0]


def plain_quantile(values, level):
    values = sorted(values)
    if len(values) == 1:
        return values[0]
    pos = level * (len(values) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(values) - 1)
    return values[lo] + (values[hi] - values[lo]) * (pos - lo)


def normal_z(level):
    # Standard normal quantile, Acklam's approximation (|err| < 1.2e-9 in z
    # is not claimed; this is the classic rational fit, ample here).
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285984469534e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if level < plow:
        q = math.sqrt(-2 * math.log(level))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    if level > phigh:
        q = math.sqrt(-2 * math.log(1 - level))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
                ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    q = level - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)


def forecast_case(case):
    levels = case.get("levels") or LEVELS_DEFAULT
    known = case.get("known") or []
    p = 0
    for row in known:
        if row.get("hour") == 1 and row.get("viewers"):
            p = row["viewers"]
            break
    if not p and known and known[0].get("viewers"):
        p = known[0]["viewers"]
    if not p:
        hist = case.get("history_hourly") or []
        p = hist[-1]["viewers"] if hist and hist[-1].get("viewers") else 1
    p = max(float(p), 1.0)
    log_p = math.log(p)
    log10_p = math.log10(p)

    weighted = {h: [] for h in range(2, 25)}
    pooled = {h: [] for h in range(2, 25)}
    for comp in case.get("comparables") or []:
        if "early-access-contaminated" in (comp.get("unlock") or ""):
            continue
        hours = release_day_hours(comp)
        first = hours.get(1)
        if not first or first < MIN_HOUR1:
            continue
        log_first = math.log(first)
        dist = (math.log10(first) - log10_p) / BANDWIDTH_DEX
        weight = math.exp(-0.5 * dist * dist)
        for h in range(2, 25):
            viewers = hours.get(h)
            if viewers and viewers > 0:
                ratio = math.log(viewers) - log_first
                pooled[h].append(ratio)
                if weight >= MIN_WEIGHT:
                    weighted[h].append((ratio, weight))

    rows = []
    for h in range(2, 25):
        pairs = weighted[h]
        if len(pairs) >= MIN_POINTS:
            qs = [weighted_quantile(pairs, t) for t in levels]
        elif len(pooled[h]) >= MIN_POINTS:
            qs = [plain_quantile(pooled[h], t) for t in levels]
        else:
            spread = 0.15 * math.sqrt(h - 1)
            qs = [z * spread for z in (normal_z(t) for t in levels)]
        median_q = qs[len(qs) // 2]
        qs = [min(max(q, median_q - TAIL_CAP), median_q + TAIL_CAP) for q in qs]
        values = [math.exp(log_p + q) for q in qs]
        values = sorted(values)
        values = [min(max(v, 0.01), 1e9) for v in values]
        rows.append({"hour": h, "q": [round(v, 2) for v in values]})
    return {"hours": rows}


def main():
    in_dir = os.environ.get("MAKEFUN_INPUT_DIR") or "/input"
    out_dir = os.environ.get("MAKEFUN_OUTPUT_DIR") or "/output"
    os.makedirs(out_dir, exist_ok=True)
    names = sorted(n for n in os.listdir(in_dir) if n.endswith(".in"))
    for name in names:
        with open(os.path.join(in_dir, name), encoding="utf-8") as fh:
            case = json.load(fh)
        result = forecast_case(case)
        out_name = name[:-3] + ".out"
        with open(os.path.join(out_dir, out_name), "w", encoding="utf-8") as fh:
            json.dump(result, fh, separators=(",", ":"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
