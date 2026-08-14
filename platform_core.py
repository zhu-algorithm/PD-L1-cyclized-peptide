"""PD-L1 环肽发现平台：七个可替换模块的本地演示实现。

所有数值为透明的规则型 surrogate，只适用于流程演示和候选相对排序。
"""
from __future__ import annotations
import math, random
from dataclasses import dataclass
from datetime import datetime, timezone

AA = "ACDEFGHIKLMNPQRSTVWY"; HYDRO = set("AVILMFWY"); CHARGED = set("DEKR")
LITERATURE = [
 {"id":"PDB:5J89","title":"PD-L1 dimer small-molecule binding pocket","type":"structure","evidence":"His66 / Gln116 / Tyr123 surface pocket"},
 {"id":"PDB:5N2F","title":"PD-1 / PD-L1 complex","type":"structure","evidence":"interface reference for epitope design"},
 {"id":"CURATED:PDL1-001","title":"PD-L1 cyclic-peptide discovery workflow","type":"workflow","evidence":"requires orthogonal experimental confirmation"},
]
def sigmoid(x): return 1/(1+math.exp(-max(-30,min(30,x))))
def helm(seq): return "PEPTIDE1{"+".".join(seq)+"}$$$$"
def parse(value):
 value=value.strip(); prefix="PEPTIDE1{"; seq=value[len(prefix):value.index("}")].replace(".","") if value.startswith(prefix) and "}" in value else "".join(c for c in value.upper() if c.isalpha())
 if not 5<=len(seq)<=40 or any(c not in AA for c in seq): raise ValueError("序列应为 5–40 个标准氨基酸，或合法 HELM。")
 return seq
@dataclass
class TargetLiteratureKB:
 def search(self, query="PD-L1"):
  q=query.lower(); return [x for x in LITERATURE if q in (x["title"]+x["evidence"]+x["id"]).lower() or q in ("pd-l1","pdl1")]
 def target(self): return {"target":"PD-L1 / CD274","focus":"protein-protein interaction surface", "key_residues":["His66","Gln116","Tyr123"],"note":"本地示例知识卡；生产环境应接入受许可的结构与文献数据库。"}
@dataclass
class CyclicPeptideGenerator:
 def generate(self, seed, count=30, random_seed=None):
  seed=parse(seed); rng=random.Random(random_seed); pool={seed}
  while len(pool)<min(max(int(count),5),200):
   s=list(rng.choice(tuple(pool))); action=rng.choice(("replace","replace","insert","delete"))
   if action=="replace": s[rng.randrange(len(s))]=rng.choice(AA)
   elif action=="insert" and len(s)<20: s.insert(rng.randrange(len(s)+1),rng.choice(AA))
   elif action=="delete" and len(s)>6: s.pop(rng.randrange(len(s)))
   pool.add("".join(s))
  return [{"sequence":s,"helm":helm(s),"cyclization":"head-to-tail (conceptual)"} for s in pool]
@dataclass
class BindingSelectivityPredictor:
 def predict(self, seq):
  n=len(seq); h=sum(c in HYDRO for c in seq)/n; a=sum(c in "FWY" for c in seq)/n; ch=sum(c in CHARGED for c in seq)/n
  affinity=-5.2-2.8*h-a+.7*ch
  # Normalized multi-target surrogate. Higher values mean stronger predicted binding.
  multi={"PD-L1":sigmoid(-affinity-5.0),
         "target_X":sigmoid(1.4*h+.9*ch+.15*abs(n-10)-.9),
         "target_Y":sigmoid(1.2*a+.7*h-.35*ch-.75)}
  profile=MultiTargetSelectivityOptimizer.calculate_profile("PD-L1", multi)
  return {"affinity_kcal_mol":round(affinity,2),"selectivity_score":round(profile["overall_selectivity"],3),
          "multi_target_affinities":{k:round(v,3) for k,v in multi.items()}, "selectivity_profile":profile,
          "model":"multi-target sequence-composition surrogate"}

@dataclass
class MultiTargetSelectivityOptimizer:
 """Adapted from the Colab workflow: optimize primary-vs-off-target selectivity."""
 primary_target_name: str = "PD-L1"
 target_names: tuple = ("PD-L1", "target_X", "target_Y")
 selectivity_threshold: float = .8
 primary_affinity_threshold: float = .6

 @staticmethod
 def calculate_profile(primary_target_name, affinities):
  primary=affinities.get(primary_target_name, 0.0)
  secondary={k:v for k,v in affinities.items() if k != primary_target_name}
  total=sum(secondary.values())
  return {"overall_selectivity":primary/(primary+total+1e-6), "primary_affinity":primary,
          "secondary_affinities_sum":total, "off_target_affinities":secondary}

 @staticmethod
 def _variant(sequence, rng):
  seq=list(sequence); seq[rng.randrange(len(seq))]=rng.choice(AA); return "".join(seq)

 def optimize_for_selectivity(self, sequence, optimization_steps=5):
  predictor=BindingSelectivityPredictor(); rng=random.Random(sequence)
  best=sequence; best_pred=predictor.predict(best); best_score=best_pred["selectivity_score"]; history=[]
  for step in range(optimization_steps):
   variant=self._variant(best, rng); pred=predictor.predict(variant); score=pred["selectivity_score"]
   if score > best_score: best, best_pred, best_score = variant, pred, score
   history.append({"step":step+1, "candidate_selectivity":round(score,4), "best_selectivity":round(best_score,4),
                   "primary_affinity":round(pred["selectivity_profile"]["primary_affinity"],4)})
  return {"sequence":best, "selectivity_optimization":{"final_score":round(best_score,4),
          "final_selectivity_profile":best_pred["selectivity_profile"], "history":history}}
@dataclass
class ADMETSynthesizability:
 def assess(self, seq):
  n=len(seq); h=sum(c in HYDRO for c in seq)/n; ch=sum(c in CHARGED for c in seq)/n; ar=sum(c in "FWY" for c in seq)/n
  permeability=sigmoid(2.2*h+.5*ar-2.4*ch-.25*(n-10)); toxicity=sigmoid(2*ar+h+.8*seq.count("C")/n-1.55)
  synth=max(.05,min(.99,.92-.035*max(0,n-8)-.10*(seq.count("W")+seq.count("C"))/n))
  return {"estimated_mw":round(n*110,1),"permeability_score":round(permeability,3),"toxicity_risk":round(toxicity,3),"synthesizability":round(synth,3),"model":"sequence-composition surrogate"}
@dataclass
class DockingModule:
 def dock(self, seq, affinity):
  # Deterministic pseudo-pose: an integration seam for Vina/MD services, not a docking calculation.
  pose=abs(hash(seq))%10000; score=affinity-(pose%31-15)/100
  return {"docking_score_kcal_mol":round(score,2),"pose_id":f"DEMO-POSE-{pose:04d}","pocket":"PD-L1 surface pocket (His66/Gln116/Tyr123)","model":"deterministic placeholder; replace with validated docking engine"}
@dataclass
class CandidateRanker:
 def rank(self, rows, weights):
  defaults={"affinity":.35,"selectivity":.2,"admet":.25,"synthesis":.1,"docking":.1}; w={k:max(0,float(weights.get(k,v))) for k,v in defaults.items()}; total=sum(w.values()) or 1; w={k:v/total for k,v in w.items()}
  for x in rows:
   b=x["binding"]; a=x["admet"]; d=x["docking"]; aff=max(0,1-abs(b["affinity_kcal_mol"]+8)/4); dock=max(0,1-abs(d["docking_score_kcal_mol"]+8)/4)
   x["priority_score"]=round(w["affinity"]*aff+w["selectivity"]*b["selectivity_score"]+w["admet"]*((a["permeability_score"]+1-a["toxicity_risk"])/2)+w["synthesis"]*a["synthesizability"]+w["docking"]*dock,4)
  return sorted(rows,key=lambda x:x["priority_score"],reverse=True),w
@dataclass
class ExperimentalFeedbackLoop:
 def record(self, candidate_id, assay, value, unit=""):
  return {"candidate_id":candidate_id,"assay":assay,"value":float(value),"unit":unit,"recorded_at":datetime.now(timezone.utc).isoformat(),"next_action":"加入校准队列；待积累足量已质控数据后重新训练。"}
def run_pipeline(seed,count,weights):
 gen=CyclicPeptideGenerator(); bind=BindingSelectivityPredictor(); optimizer=MultiTargetSelectivityOptimizer(); admet=ADMETSynthesizability(); dock=DockingModule(); rows=[]
 for i,item in enumerate(gen.generate(seed,count),1):
  optimized=optimizer.optimize_for_selectivity(item["sequence"]); seq=optimized["sequence"]; item.update({"sequence":seq,"helm":helm(seq)})
  b=bind.predict(seq); rows.append({"id":f"CP-{i:03d}",**item,"optimization":optimized["selectivity_optimization"],"binding":b,"admet":admet.assess(seq),"docking":dock.dock(seq,b["affinity_kcal_mol"])})
 ranked, norm=CandidateRanker().rank(rows,weights)
 for i,x in enumerate(ranked,1): x["rank"]=i
 return {"created_at":datetime.now(timezone.utc).isoformat(),"weights":norm,"candidates":ranked,"disclaimer":"所有预测和对接字段均为本地规则型演示 surrogate，不可替代已验证模型、专业对接或实验。"}
