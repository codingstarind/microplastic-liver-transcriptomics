"""Export the prompts typed to Claude Code for this project from the local session transcripts.

Reads the Claude Code transcript files (JSONL) kept on the analysis machine under
~/.claude/projects/C--Users-Soham-Desktop-microplastic/ and writes one Markdown file per session
with every user-typed prompt (pasted text included) and its UTC timestamp. Assistant replies, tool
calls, tool outputs and automatic system context are not exported. The session about an unrelated
project (OlmoCR benchmark, 91b92e0b) is skipped. Messages shared by a forked session and its parent
are written only once. Email addresses are redacted.

Run from the repository root:  python supplementary/prompts/extract_prompts.py
"""
import json
import re
from collections import Counter
from pathlib import Path

SRC = Path.home() / ".claude/projects/C--Users-Soham-Desktop-microplastic"
OUT = Path(__file__).resolve().parent
UNRELATED = {"91b92e0b-a779-4fb0-878e-e4596dd08c4d"}
EMAIL = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")
AUTOMATIC = ("<command-name>", "<task-notification>", "<ci-monitor-event>", "<side-sessions-event>")
DROP_BLOCKS = re.compile(r"<(system-reminder|local-command-stdout|local-command-stderr|command-message|command-args)>.*?</\1>", re.S)


def prompt_text(rec):
    if rec.get("type") != "user" or rec.get("isMeta") or rec.get("isCompactSummary"):
        return None
    c = rec.get("message", {}).get("content")
    if isinstance(c, list):
        if any(isinstance(x, dict) and x.get("type") == "tool_result" for x in c):
            return None
        c = "\n".join(x.get("text", "") for x in c if isinstance(x, dict) and x.get("type") == "text")
    if not isinstance(c, str):
        return None
    c = DROP_BLOCKS.sub("", c).strip()
    if (not c or c.startswith("This session is being continued from a previous conversation")
            or c.startswith(AUTOMATIC) or c == "[Request interrupted by user]"):
        return None  # automatic app/agent notices, not typed prompts
    c = re.sub(r"</?pasted_content[^>]*>", "", c).strip()
    return EMAIL.sub("[email redacted]", c)


FORKS = {"094de3c2-76e3-4802-a7f2-f451c3958dff"}  # fork of 59dc8d15: export the parent first


def first_ts(p):
    for line in open(p, encoding="utf-8"):
        try:
            ts = json.loads(line).get("timestamp")
        except json.JSONDecodeError:
            continue
        if ts:
            return ts
    return ""



def main():
    seen, index = set(), []
    for f in sorted(SRC.glob("*.jsonl"), key=lambda p: (p.stem in FORKS, first_ts(p))):
        sid = f.stem
        if sid in UNRELATED:
            continue
        title, models, prompts, first, last = None, Counter(), [], None, None
        for line in open(f, encoding="utf-8"):
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            ts = rec.get("timestamp")
            if ts:
                first, last = first or ts, ts
            if rec.get("type") in ("custom-title", "ai-title"):
                title = rec.get("customTitle") or rec.get("aiTitle") or title
            if rec.get("type") == "assistant" and rec.get("message", {}).get("model") not in (None, "<synthetic>"):
                models[rec["message"]["model"]] += 1
            t = prompt_text(rec)
            if t is None:
                continue
            key = rec.get("uuid") or (ts, t)
            if key in seen:
                continue
            seen.add(key)
            prompts.append((ts, t))
        name = f"session_{(first or '')[:10]}_{sid[:8]}.md"
        with open(OUT / name, "w", encoding="utf-8") as o:
            o.write(f"# Prompts: {title or 'untitled session'}\n\n")
            o.write(f"- Claude Code session ID: `{sid}`\n- First / last activity (UTC): {first} / {last}\n")
            o.write(f"- Models recorded in the main conversation: {', '.join(f'{m} ({n} replies)' for m, n in models.most_common())}\n")
            o.write(f"- Prompts exported here: {len(prompts)} (prompts already exported from a parent session are not repeated)\n\n")
            for i, (ts, t) in enumerate(prompts, 1):
                o.write(f"## Prompt {i} ({ts})\n\n{t}\n\n")
        index.append(dict(file=name, session_id=sid, title=title, first=first, last=last, n_prompts=len(prompts),
                          models=dict(models)))
    json.dump(index, open(OUT / "prompt_index.json", "w"), indent=2)
    for r in index:
        print(r["file"], r["n_prompts"], r["title"])


if __name__ == "__main__":
    main()
