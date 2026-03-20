import os
import re
import json
import time
import hashlib
import argparse
from pathlib import Path
from typing import List, Dict, Any

import fitz
import pandas as pd
from tqdm import tqdm
from openai import OpenAI

FIELD_SCHEMA = {
    "irrigation_target": [
        "Area equipped for irrigation",
        "Area actually irrigated",
        "Harvested irrigated area",
        "Mixed",
        "Unclear",
    ],
    "representation_unit": [
        "Binary presence",
        "Area magnitude",
        "Both",
        "Unclear",
    ],
    "temporal_concept": [
        "Annual",
        "Seasonal",
        "Reference-period snapshot",
        "Unclear",
    ],
    "paddy_treatment": [
        "Included with irrigation",
        "Separated from irrigation",
        "Unclear",
    ],
    "permanent_crops_treatment": [
        "Included with irrigation",
        "Separated from irrigation",
        "Unclear",
    ],
    "irrigation_basis": [
        "Infrastructure-based",
        "Water-use based",
        "Mixed",
        "Unclear",
    ],
}

KEYWORDS = [
    # Core irrigation definitions
    "irrigation", "irrigated", "area equipped for irrigation",
    "actually irrigated", "harvested irrigated area",

    # Representation
    "presence", "binary", "fraction", "area", "hectares",

    # Temporal concepts
    "annual", "seasonal", "monthly", "reference year",
    "snapshot", "time period", "multi-year",

    # Paddy / rice
    "paddy", "rice", "flooded", "flooding", "wetland rice",

    # Permanent crops
    "permanent crop", "orchard", "vineyard", "tree crops",

    # Irrigation basis
    "water use", "water withdrawal", "actual application",
    "infrastructure", "irrigation system", "equipped area",

    # Strong signals
    "definition", "defined as", "we define",
    "data source", "method", "methods", "dataset",
    "table", "appendix"
]

SYSTEM_PROMPT = """
You are a strict, adversarial extraction assistant for a scientific audit of irrigation-map definitions.

Your task is to classify how a paper defines irrigation across predefined dimensions.

-------------------------
EVIDENCE AND CONSERVATISM
-------------------------
1. Use ONLY explicit statements from the provided excerpts.
2. If the definition is not explicitly stated, return "Unclear".
3. Do NOT infer definitions from context, common practice, or prior knowledge.
4. If multiple definitions are used in the same paper, classify as "Mixed" where applicable.
5. When in doubt, choose "Unclear".

-------------------------
STRICT LABEL SELECTION
-------------------------
6. You MUST select labels ONLY from the provided allowed list.
7. Do NOT paraphrase, generalize, or invent labels.
8. Any deviation from allowed labels is invalid.

-------------------------
JUSTIFICATION AND CITATION
-------------------------
9. Each classification MUST include:
   - label (from allowed list)
   - confidence (0–1)
   - short reason (1–2 sentences)
   - citation with:
        • page number
        • verbatim snippet copied EXACTLY from the text

10. The snippet MUST directly support the label.
11. If no supporting snippet exists → label MUST be "Unclear".

-------------------------
CONFIDENCE CALIBRATION
-------------------------
12. Use high confidence (≥0.8) ONLY if:
    - definition is explicit and unambiguous
13. Use medium confidence (0.4–0.7) if:
    - partially supported or somewhat indirect
14. Use low confidence (≤0.3) if:
    - weak or ambiguous evidence

-------------------------
OUTPUT FORMAT
-------------------------
15. Return VALID JSON ONLY.
16. Do NOT include explanations outside JSON.
17. Do NOT omit any field.

-------------------------
PHILOSOPHY
-------------------------
You are auditing definitional consistency, not summarizing the paper.
Absence of evidence = "Unclear".
Overclassification is a critical error.
"""

USER_PROMPT_TEMPLATE = """
Classify this paper using ONLY the allowed labels below.

You MUST strictly select one label per category.

-------------------------
ALLOWED LABELS
-------------------------
{schema_json}

-------------------------
IMPORTANT INSTRUCTIONS
-------------------------
- Use ONLY explicit textual evidence
- If no explicit definition is found → return "Unclear"
- Do NOT infer or assume definitions
- If multiple definitions coexist → return "Mixed" where applicable; otherwise return "Unclear"
- Every classification MUST include a supporting quote

-------------------------
OUTPUT FORMAT (STRICT)
-------------------------
Return EXACTLY this JSON structure:

{{
  "paper_title": "...",
  "classifications": {{
    "irrigation_target": {{
      "label": "...",
      "confidence": 0.0,
      "reason": "...",
      "citation": {{"page": 0, "snippet": "..."}}
    }},
    "representation_unit": {{
      "label": "...",
      "confidence": 0.0,
      "reason": "...",
      "citation": {{"page": 0, "snippet": "..."}}
    }},
    "temporal_concept": {{
      "label": "...",
      "confidence": 0.0,
      "reason": "...",
      "citation": {{"page": 0, "snippet": "..."}}
    }},
    "paddy_treatment": {{
      "label": "...",
      "confidence": 0.0,
      "reason": "...",
      "citation": {{"page": 0, "snippet": "..."}}
    }},
    "permanent_crops_treatment": {{
      "label": "...",
      "confidence": 0.0,
      "reason": "...",
      "citation": {{"page": 0, "snippet": "..."}}
    }},
    "irrigation_basis": {{
      "label": "...",
      "confidence": 0.0,
      "reason": "...",
      "citation": {{"page": 0, "snippet": "..."}}
    }}
  }},
  "supports_adversarial_recode": true,
  "adversarial_recode_notes": "Short note describing how this paper should be treated in a definition audit"
}}

-------------------------
EXCERPTS
-------------------------
{excerpts}
"""

VERIFY_SYSTEM_PROMPT = """
You are a strict scientific verifier.

Your task is to audit whether each classification is supported by its cited snippet.

For each field:
- KEEP the label only if the snippet directly supports it
- Otherwise REPLACE the label with "Unclear"
- LOWER confidence if support is weak or indirect
- DO NOT introduce new interpretations
- DO NOT add new evidence

STRICT RULES:
- If the snippet does not explicitly justify the label → label = "Unclear"
- If snippet is vague or ambiguous → reduce confidence

Return VALID JSON ONLY.
"""

VERIFY_USER_TEMPLATE = """
Allowed labels:
{schema_json}

Review this classification and verify whether each label is supported by its cited snippet.

Input:
{draft_json}
"""


def normalize_whitespace(text: str) -> str:
    text = text.replace("\u00ad", "")
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def extract_pages(pdf_path: Path) -> List[Dict[str, Any]]:
    doc = fitz.open(pdf_path)
    pages = []
    for i, page in enumerate(doc):
        text = page.get_text("text")
        text = normalize_whitespace(text)
        if text:
            pages.append({"page": i + 1, "text": text})
    doc.close()
    return pages


def score_page(text: str) -> int:
    t = text.lower()
    score = sum(t.count(k.lower()) for k in KEYWORDS)
    for sec in ["abstract", "methods", "method", "data", "discussion", "conclusion", "table"]:
        score += 3 * t.count(sec)
    return score


def select_relevant_pages(pages: List[Dict[str, Any]], max_pages: int = 12) -> List[Dict[str, Any]]:
    scored = [(score_page(p["text"]), p) for p in pages]
    scored.sort(key=lambda x: x[0], reverse=True)
    selected = [p for s, p in scored if s > 0][:max_pages]
    if not selected:
        selected = pages[:min(8, len(pages))]
    selected.sort(key=lambda x: x["page"])
    return selected


def build_excerpt_bundle(pages: List[Dict[str, Any]], max_chars_per_page: int = 7000) -> str:
    parts = []
    for p in pages:
        txt = p["text"][:max_chars_per_page]
        parts.append(f"[PAGE {p['page']}]\n{txt}")
    return "\n\n".join(parts)


def detect_definition_signals(text: str) -> int:
    t = text.lower()

    phrases = [
        "defined as",
        "we define",
        "is defined as",
        "refers to",
        "for the purpose of this study",
    ]

    keyword_hits = sum(t.count(k.lower()) for k in KEYWORDS)
    phrase_hits = sum(t.count(p) for p in phrases)

    return keyword_hits + 3 * phrase_hits


def evaluate_page_selection(pages, selected_pages):
    total_scores = [detect_definition_signals(p["text"]) for p in pages]
    selected_scores = [detect_definition_signals(p["text"]) for p in selected_pages]

    total_signal = sum(total_scores)
    selected_signal = sum(selected_scores)

    coverage = any(s > 0 for s in selected_scores)
    recall_proxy = selected_signal / total_signal if total_signal > 0 else None

    return {
        "n_pages_total": len(pages),
        "n_pages_selected": len(selected_pages),
        "total_signal": total_signal,
        "selected_signal": selected_signal,
        "recall_proxy": recall_proxy,
        "coverage": coverage,
    }


def get_client() -> OpenAI:
    base_url = os.environ.get("OPENAI_BASE_URL", "http://localhost:11434/v1")
    api_key = os.environ.get("OPENAI_API_KEY", "ollama")
    return OpenAI(
        base_url=base_url,
        api_key=api_key,
        timeout=600
    )


def call_model(client: OpenAI, model: str, system_prompt: str, user_prompt: str, temperature: float = 0.0) -> str:
    response = client.chat.completions.create(
        model=model,
        temperature=temperature,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
    )
    return response.choices[0].message.content


def safe_json_loads(text: str) -> Dict[str, Any]:
    text = text.strip()
    text = re.sub(r"^```json\s*", "", text)
    text = re.sub(r"^```\s*", "", text)
    text = re.sub(r"\s*```$", "", text)
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        m = re.search(r"\{.*\}", text, flags=re.DOTALL)
        if m:
            return json.loads(m.group(0))
        raise


def validate_result(result: Dict[str, Any]) -> Dict[str, Any]:
    out = result.copy()
    cls = out.get("classifications", {})

    for field, allowed in FIELD_SCHEMA.items():
        if field not in cls:
            cls[field] = {
                "label": "Unclear",
                "confidence": 0.0,
                "reason": "Missing field in output.",
                "citation": {"page": None, "snippet": ""}
            }
        else:
            label = cls[field].get("label")
            if label not in allowed:
                cls[field]["label"] = "Unclear"
                cls[field]["confidence"] = 0.0
                cls[field]["reason"] = "Invalid label returned by model."
                cls[field]["citation"] = {"page": None, "snippet": ""}

            cls[field].setdefault("confidence", 0.0)
            cls[field].setdefault("reason", "")
            cls[field].setdefault("citation", {"page": None, "snippet": ""})

    out["classifications"] = cls
    out.setdefault("supports_adversarial_recode", False)
    out.setdefault("adversarial_recode_notes", "")
    return out


def file_fingerprint(pdf_path: Path) -> str:
    h = hashlib.sha256()
    stat = pdf_path.stat()
    h.update(str(pdf_path.name).encode())
    h.update(str(stat.st_size).encode())
    h.update(str(stat.st_mtime_ns).encode())
    return h.hexdigest()[:16]


def cache_path_for(pdf_path: Path, cache_dir: Path) -> Path:
    return cache_dir / f"{pdf_path.stem}__{file_fingerprint(pdf_path)}.json"


def process_pdf(
    pdf_path: Path,
    cache_dir: Path,
    client: OpenAI,
    model: str,
    verify: bool = True,
    sleep_s: float = 0.0
) -> Dict[str, Any]:
    cpath = cache_path_for(pdf_path, cache_dir)
    if cpath.exists():
        with open(cpath, "r", encoding="utf-8") as f:
            return json.load(f)

    pages = extract_pages(pdf_path)
    selected_pages = select_relevant_pages(pages, max_pages=8)
    selection_eval = evaluate_page_selection(pages, selected_pages)
    excerpts = build_excerpt_bundle(selected_pages)

    user_prompt = USER_PROMPT_TEMPLATE.format(
        schema_json=json.dumps(FIELD_SCHEMA, indent=2),
        excerpts=excerpts
    )

    raw = call_model(client, model, SYSTEM_PROMPT, user_prompt, temperature=0.0)
    draft = safe_json_loads(raw)
    draft = validate_result(draft)

    verified = draft
    verify_raw = None

    if verify:
        verify_prompt = VERIFY_USER_TEMPLATE.format(
            schema_json=json.dumps(FIELD_SCHEMA, indent=2),
            draft_json=json.dumps(draft, ensure_ascii=False, indent=2)
        )
        verify_raw = call_model(client, model, VERIFY_SYSTEM_PROMPT, verify_prompt, temperature=0.0)
        verified = safe_json_loads(verify_raw)
        verified = validate_result(verified)

    verified["_source_file"] = pdf_path.name
    verified["_source_path"] = str(pdf_path)
    verified["_pages_used"] = [p["page"] for p in selected_pages]
    verified["_selection_eval"] = selection_eval
    verified["_n_pages_total"] = len(pages)
    verified["_cache_file"] = str(cpath)
    verified["_raw_model_output"] = raw
    verified["_raw_verify_output"] = verify_raw

    with open(cpath, "w", encoding="utf-8") as f:
        json.dump(verified, f, ensure_ascii=False, indent=2)

    if sleep_s > 0:
        time.sleep(sleep_s)

    return verified


def flatten_result(result: Dict[str, Any]) -> Dict[str, Any]:
    row = {
        "source_file": result.get("_source_file"),
        "source_path": result.get("_source_path"),
        "paper_title": result.get("paper_title", ""),
        "supports_adversarial_recode": result.get("supports_adversarial_recode", False),
        "adversarial_recode_notes": result.get("adversarial_recode_notes", ""),
        "pages_used": ",".join(map(str, result.get("_pages_used", []))),
        "n_pages_total": result.get("_n_pages_total", None),
        "cache_file": result.get("_cache_file", "")
    }

    sel = result.get("_selection_eval", {})
    row["selection_recall_proxy"] = sel.get("recall_proxy")
    row["selection_total_signal"] = sel.get("total_signal")
    row["selection_selected_signal"] = sel.get("selected_signal")
    row["selection_coverage"] = sel.get("coverage")

    for field in FIELD_SCHEMA.keys():
        item = result["classifications"].get(field, {})
        row[field] = item.get("label", "")
        row[f"{field}_confidence"] = item.get("confidence", "")
        row[f"{field}_reason"] = item.get("reason", "")
        cit = item.get("citation", {})
        row[f"{field}_citation_page"] = cit.get("page", "")
        row[f"{field}_citation_snippet"] = cit.get("snippet", "")

    return row


def rebuild_outputs(cache_dir: Path, out_csv: Path, out_xlsx: Path, out_jsonl: Path):
    cache_files = sorted(cache_dir.glob("*.json"))
    results = []
    for fp in cache_files:
        with open(fp, "r", encoding="utf-8") as f:
            results.append(json.load(f))

    rows = [flatten_result(r) for r in results]
    df = pd.DataFrame(rows)

    if len(df) > 0:
        df = df.sort_values(by=["source_file"]).reset_index(drop=True)

    df.to_csv(out_csv, index=False)
    df.to_excel(out_xlsx, index=False)

    with open(out_jsonl, "w", encoding="utf-8") as f:
        for r in results:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdf_dir", required=True)
    parser.add_argument("--cache_dir", default="./cache")
    parser.add_argument("--out_dir", default="./outputs")
    parser.add_argument("--model", default=os.environ.get("OPENAI_MODEL", "deepseek-r1:32b"))
    parser.add_argument("--sleep_s", type=float, default=0.0)
    parser.add_argument("--no_verify", action="store_true")
    parser.add_argument("--rebuild_only", action="store_true")
    args = parser.parse_args()

    pdf_dir = Path(args.pdf_dir)
    cache_dir = Path(args.cache_dir)
    out_dir = Path(args.out_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    out_csv = out_dir / "irrigation_definition_table.csv"
    out_xlsx = out_dir / "irrigation_definition_table.xlsx"
    out_jsonl = out_dir / "irrigation_definition_details.jsonl"

    if args.rebuild_only:
        rebuild_outputs(cache_dir, out_csv, out_xlsx, out_jsonl)
        print(f"Rebuilt outputs from cache into {out_dir}")
        return

    pdf_files = sorted(pdf_dir.glob("*.pdf"))
    if not pdf_files:
        raise FileNotFoundError(f"No PDFs found in {pdf_dir}")

    client = get_client()
    failed = []

    for pdf_path in tqdm(pdf_files, desc="Processing PDFs"):
        try:
            process_pdf(
                pdf_path=pdf_path,
                cache_dir=cache_dir,
                client=client,
                model=args.model,
                verify=not args.no_verify,
                sleep_s=args.sleep_s
            )
            rebuild_outputs(cache_dir, out_csv, out_xlsx, out_jsonl)
        except Exception as e:
            failed.append({"file": pdf_path.name, "error": str(e)})
            fail_path = out_dir / "failures.json"
            with open(fail_path, "w", encoding="utf-8") as f:
                json.dump(failed, f, ensure_ascii=False, indent=2)

    print(f"Done. Outputs written to {out_dir}")
    if failed:
        print(f"{len(failed)} files failed. See {out_dir / 'failures.json'}")


if __name__ == "__main__":
    main()
