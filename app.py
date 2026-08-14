#!/usr/bin/env python3
"""PD-L1 环肽抑制剂发现平台（本地、零第三方依赖演示版）。

此应用用于候选优先级排序，不提供实验效力、临床或安全性结论。
运行：python app.py，然后访问 http://127.0.0.1:8765
"""
from __future__ import annotations

import csv
import json
import math
import random
import re
import threading
import webbrowser
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse
from platform_core import TargetLiteratureKB, ExperimentalFeedbackLoop, run_pipeline

ROOT = Path(__file__).resolve().parent
RUNS = ROOT / "runs"
RUNS.mkdir(exist_ok=True)
AA = "ACDEFGHIKLMNPQRSTVWY"
HYDRO = set("AVILMFWY")
POLAR = set("STNQCH")
CHARGED = set("DEKR")
HELMPATTERN = re.compile(r"^PEPTIDE1\{([A-Za-z.]+)\}\$\$\$\$$")


def clean_sequence(value: str) -> str:
    value = value.strip()
    match = HELMPATTERN.fullmatch(value)
    if match:
        tokens = match.group(1).split(".")
    else:
        tokens = list(re.sub(r"[^A-Za-z]", "", value).upper())
    if not tokens or any(len(t) != 1 or t.upper() not in AA for t in tokens):
        raise ValueError("请输入仅由 20 种标准氨基酸组成的 HELM 或单字母序列。")
    if not 5 <= len(tokens) <= 40:
        raise ValueError("环肽长度需在 5–40 个残基之间。")
    return "".join(t.upper() for t in tokens)


def helm(seq: str) -> str:
    return "PEPTIDE1{" + ".".join(seq) + "}$$$$"


def _sigmoid(x: float) -> float:
    return 1 / (1 + math.exp(-max(-30, min(30, x))))


def score_sequence(seq: str, weights: dict[str, float]) -> dict:
    """透明的规则型 surrogate。数值只能用于相对排序。"""
    n = len(seq)
    hydro = sum(x in HYDRO for x in seq) / n
    polar = sum(x in POLAR for x in seq) / n
    charged = sum(x in CHARGED for x in seq) / n
    aromatic = sum(x in "FWY" for x in seq) / n
    proline = seq.count("P") / n
    cysteine = seq.count("C") / n
    unique = len(set(seq)) / min(n, len(AA))
    estimated_mw = round(n * 110.0, 1)
    # 这些是基于序列组成的可解释启发式，不是训练出的预测器。
    affinity = -5.2 - 2.8 * hydro - 1.0 * aromatic + 0.7 * charged + 0.35 * proline
    permeability = _sigmoid(2.2 * hydro + 0.5 * aromatic - 2.4 * charged - 0.25 * (n - 10))
    toxicity = _sigmoid(2.0 * aromatic + 1.0 * hydro + 0.8 * cysteine - 1.55)
    druglikeness = max(0.0, min(1.0, 0.88 - abs(n - 9) * 0.035 - charged * 0.20 + unique * 0.12))
    affinity_reward = max(0.0, 1.0 - abs(affinity + 8.0) / 4.0)
    reward = (
        weights["affinity"] * affinity_reward
        + weights["permeability"] * permeability
        + weights["toxicity"] * (1 - toxicity)
        + weights["druglikeness"] * druglikeness
    )
    flags = []
    if estimated_mw > 1200: flags.append("估算 MW 偏高")
    if charged > 0.30: flags.append("带电残基比例偏高")
    if toxicity > 0.65: flags.append("结构性风险提示")
    if cysteine >= 0.25: flags.append("Cys 比例偏高")
    return {
        "sequence": seq, "helm": helm(seq), "length": n, "estimated_mw": estimated_mw,
        "hydrophobic_fraction": round(hydro, 3), "charged_fraction": round(charged, 3),
        "diversity": round(unique, 3), "affinity": round(affinity, 2),
        "permeability": round(permeability, 3), "toxicity": round(toxicity, 3),
        "druglikeness": round(druglikeness, 3), "reward": round(reward, 4), "flags": flags,
    }


def mutate(parent: str, rng: random.Random) -> str:
    seq = list(parent)
    for _ in range(rng.choice((1, 1, 1, 2))):
        kind = rng.choice(("replace", "replace", "insert", "delete"))
        if kind == "replace": seq[rng.randrange(len(seq))] = rng.choice(AA)
        elif kind == "insert" and len(seq) < 20: seq.insert(rng.randrange(len(seq) + 1), rng.choice(AA))
        elif kind == "delete" and len(seq) > 6: seq.pop(rng.randrange(len(seq)))
    return "".join(seq)


def generate(payload: dict) -> dict:
    seed = clean_sequence(payload.get("seed", "ACLVIFWY"))
    count = max(10, min(500, int(payload.get("count", 80))))
    raw = payload.get("weights", {})
    weights = {k: max(0, float(raw.get(k, default))) for k, default in
               {"affinity": .40, "permeability": .30, "toxicity": .15, "druglikeness": .15}.items()}
    total = sum(weights.values()) or 1
    weights = {k: v / total for k, v in weights.items()}
    rng = random.Random(payload.get("random_seed") or datetime.now().timestamp())
    pool, seen = [seed], {seed}
    for _ in range(count * 8):
        child = mutate(rng.choice(pool[-min(len(pool), 40):]), rng)
        if child not in seen:
            pool.append(child); seen.add(child)
        if len(pool) >= count * 2: break
    candidates = [score_sequence(s, weights) for s in pool]
    candidates.sort(key=lambda x: x["reward"], reverse=True)
    for i, item in enumerate(candidates[:count], 1): item["rank"] = i
    results = candidates[:count]
    run_id = datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%SZ")
    report = {"run_id": run_id, "created_at": datetime.now(timezone.utc).isoformat(), "seed": helm(seed),
              "weights": weights, "candidate_count": len(results), "candidates": results,
              "disclaimer": "规则型探索排序结果；不能代替已验证的结合、ADMET、对接或实验数据。"}
    (RUNS / f"{run_id}.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return report


def to_csv(report: dict) -> str:
    fields = ["rank", "sequence", "helm", "length", "estimated_mw", "affinity", "permeability", "toxicity", "druglikeness", "reward", "flags"]
    rows = []
    for candidate in report["candidates"]:
        row = {k: candidate.get(k, "") for k in fields}; row["flags"] = "; ".join(row["flags"]); rows.append(row)
    from io import StringIO
    output = StringIO(); writer = csv.DictWriter(output, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    return output.getvalue()


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args): pass

    def send(self, body: str | bytes, content_type="application/json; charset=utf-8", status=200, headers=None):
        raw = body.encode("utf-8") if isinstance(body, str) else body
        self.send_response(status); self.send_header("Content-Type", content_type); self.send_header("Content-Length", str(len(raw)))
        for key, value in (headers or {}).items(): self.send_header(key, value)
        self.end_headers(); self.wfile.write(raw)

    def do_GET(self):
        path = urlparse(self.path).path
        if path == "/": return self.send((ROOT / "platform-dashboard.html").read_text(encoding="utf-8"), "text/html; charset=utf-8")
        if path == "/api/target": return self.send(json.dumps(TargetLiteratureKB().target(), ensure_ascii=False))
        if path == "/api/literature": return self.send(json.dumps(TargetLiteratureKB().search(), ensure_ascii=False))
        if path == "/api/runs":
            runs = [{"id": p.stem, "created": p.stat().st_mtime} for p in sorted(RUNS.glob("run_*.json"), reverse=True)[:20]]
            return self.send(json.dumps(runs))
        if path.startswith("/api/export/"):
            target = RUNS / (Path(path).name + ".json")
            if not target.exists(): return self.send(json.dumps({"error": "未找到该运行记录"}), status=404)
            report = json.loads(target.read_text(encoding="utf-8"))
            return self.send(to_csv(report), "text/csv; charset=utf-8", headers={"Content-Disposition": f"attachment; filename={target.stem}.csv"})
        return self.send("Not found", "text/plain; charset=utf-8", 404)

    def do_POST(self):
        path = urlparse(self.path).path
        if path not in ("/api/screen", "/api/pipeline", "/api/feedback"): return self.send(json.dumps({"error": "Not found"}), status=404)
        try:
            n = int(self.headers.get("Content-Length", "0")); payload = json.loads(self.rfile.read(n) or b"{}")
            if path == "/api/pipeline":
                return self.send(json.dumps(run_pipeline(payload.get("seed", "ACLVIFWY"), payload.get("count", 30), payload.get("weights", {})), ensure_ascii=False))
            if path == "/api/feedback":
                result = ExperimentalFeedbackLoop().record(payload.get("candidate_id", ""), payload.get("assay", ""), payload.get("value", 0), payload.get("unit", ""))
                return self.send(json.dumps(result, ensure_ascii=False))
            return self.send(json.dumps(generate(payload), ensure_ascii=False))
        except (ValueError, json.JSONDecodeError) as exc:
            return self.send(json.dumps({"error": str(exc)}, ensure_ascii=False), status=400)
        except Exception:
            return self.send(json.dumps({"error": "服务器处理失败，请检查输入。"}, ensure_ascii=False), status=500)


def main():
    server = ThreadingHTTPServer(("127.0.0.1", 8765), Handler)
    print("PD-L1 环肽发现平台已启动：http://127.0.0.1:8765")
    threading.Timer(.5, lambda: webbrowser.open("http://127.0.0.1:8765")).start()
    try: server.serve_forever()
    except KeyboardInterrupt: print("\n平台已停止")

if __name__ == "__main__": main()
